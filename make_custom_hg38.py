#!/usr/bin/env python3
"""Build a custom hg38 FASTA for telomere-aware mapping.

Default behavior:
  - read ./hg38.fa
  - keep only chr1-22, chrX, chrY, and chrM
  - mask chrY PAR1/PAR2 to N using GRCh38/UCSC coordinates
  - append one 100re kb telomeric-repeat decoy contig

The script streams the input FASTA, so it is safe for multi-GB references.
"""

from __future__ import annotations

import argparse
import gzip
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, TextIO


DEFAULT_INPUT = Path("hg38.fa")
DEFAULT_OUTPUT = Path("hg38.custom.primary.telomere_decoy.ypar_masked.fa")
DEFAULT_PAR_INTERVALS = (
    (10_001, 2_781_479),
    (56_887_903, 57_217_415),
)

CANONICAL_HG38_CONTIGS = {
    *(f"chr{i}" for i in range(1, 23)),
    "chrX",
    "chrY",
    "chrM",
}


class DefaultsFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def _get_help_string(self, action: argparse.Action) -> str:
        if isinstance(action, (argparse._StoreTrueAction, argparse._StoreFalseAction)):
            return action.help or ""
        return super()._get_help_string(action)


@dataclass
class Summary:
    total_records: int = 0
    kept_records: int = 0
    dropped_non_primary_records: int = 0
    y_par_bases_covered: int = 0
    y_par_bases_changed_to_n: int = 0
    telomere_decoy_bases: int = 0


def open_text(path: Path, mode: str) -> TextIO:
    if "b" in mode:
        raise ValueError("open_text only supports text modes")
    if path.suffix == ".gz":
        return gzip.open(path, mode, encoding="ascii")
    return path.open(mode, encoding="ascii")


def parse_interval(value: str) -> tuple[int, int]:
    try:
        start_text, end_text = value.replace(",", "").split("-", 1)
        start = int(start_text)
        end = int(end_text)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            f"expected START-END, got {value!r}"
        ) from exc
    if start < 1 or end < start:
        raise argparse.ArgumentTypeError(
            f"invalid 1-based closed interval: {value!r}"
        )
    return start, end


def keep_contig(name: str, include_all: bool, excluded_contigs: set[str]) -> bool:
    return (include_all or name in CANONICAL_HG38_CONTIGS) and name not in excluded_contigs


def mask_line(
    seq: str,
    line_start_0: int,
    intervals_1based: Iterable[tuple[int, int]],
) -> tuple[str, int, int]:
    """Mask sequence positions overlapping 1-based closed intervals."""
    if not seq:
        return seq, 0, 0

    line_end_0 = line_start_0 + len(seq)
    chars: list[str] | None = None
    covered = 0
    changed = 0

    for start_1, end_1 in intervals_1based:
        interval_start_0 = start_1 - 1
        interval_end_0 = end_1
        overlap_start = max(line_start_0, interval_start_0)
        overlap_end = min(line_end_0, interval_end_0)
        if overlap_start >= overlap_end:
            continue

        rel_start = overlap_start - line_start_0
        rel_end = overlap_end - line_start_0
        if chars is None:
            chars = list(seq)

        for idx in range(rel_start, rel_end):
            if chars[idx] != "N":
                changed += 1
            chars[idx] = "N"
        covered += overlap_end - overlap_start

    if chars is None:
        return seq, covered, changed
    return "".join(chars), covered, changed


def write_telomere_decoy(
    out_handle: TextIO,
    *,
    name: str,
    motif: str,
    length: int,
    line_width: int,
) -> int:
    if length <= 0:
        return 0
    if not motif:
        raise ValueError("telomere motif must not be empty")
    if line_width <= 0:
        raise ValueError("line width must be positive")

    repeat_count = (length + len(motif) - 1) // len(motif)
    seq = (motif.upper() * repeat_count)[:length]

    out_handle.write(f">{name}\n")
    for offset in range(0, length, line_width):
        out_handle.write(seq[offset : offset + line_width])
        out_handle.write("\n")
    return length


def process_fasta(args: argparse.Namespace) -> Summary:
    input_path = args.input_fasta
    output_path = args.output_fasta
    if input_path.resolve() == output_path.resolve():
        raise SystemExit("refusing to overwrite the input FASTA")
    if output_path.exists() and not args.force:
        raise SystemExit(f"{output_path} already exists; use --force to overwrite")

    summary = Summary()
    current_name: str | None = None
    skip_record = False
    mask_current_record = False
    seq_pos_0 = 0
    excluded_contigs = set(args.exclude_contig or [])

    with open_text(input_path, "rt") as in_handle, open_text(output_path, "wt") as out_handle:
        for raw_line in in_handle:
            if raw_line.startswith(">"):
                header = raw_line[1:].strip()
                current_name = header.split()[0] if header else ""
                keep_record = keep_contig(
                    current_name,
                    include_all=args.include_all,
                    excluded_contigs=excluded_contigs,
                )

                summary.total_records += 1
                skip_record = not keep_record
                mask_current_record = (
                    args.mask_y_par and current_name == args.y_contig and not skip_record
                )
                seq_pos_0 = 0

                if skip_record:
                    summary.dropped_non_primary_records += 1
                else:
                    summary.kept_records += 1
                    out_handle.write(raw_line if raw_line.endswith("\n") else raw_line + "\n")
                continue

            if current_name is None:
                if raw_line.strip():
                    raise SystemExit("input FASTA has sequence before the first header")
                continue
            if skip_record:
                continue

            seq = raw_line.rstrip("\r\n")
            if mask_current_record:
                masked_seq, covered, changed = mask_line(seq, seq_pos_0, args.par_interval)
                summary.y_par_bases_covered += covered
                summary.y_par_bases_changed_to_n += changed
                out_handle.write(masked_seq)
                out_handle.write("\n")
            else:
                out_handle.write(raw_line if raw_line.endswith("\n") else raw_line + "\n")
            seq_pos_0 += len(seq)

        if not args.no_telomere_decoy:
            summary.telomere_decoy_bases = write_telomere_decoy(
                out_handle,
                name=args.telomere_name,
                motif=args.telomere_motif,
                length=args.telomere_length,
                line_width=args.line_width,
            )
            summary.kept_records += 1

    return summary


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Create primary-contig hg38 with chrY PAR masked and a telomere decoy.",
        formatter_class=DefaultsFormatter,
    )
    parser.add_argument("input_fasta", nargs="?", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("output_fasta", nargs="?", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--force",
        action="store_true",
        help="overwrite output FASTA if it already exists",
    )
    parser.add_argument(
        "--include-all",
        action="store_true",
        help="keep all input contigs instead of only chr1-22, chrX, chrY, and chrM",
    )
    parser.add_argument(
        "--exclude-contig",
        action="append",
        default=[],
        metavar="NAME",
        help="drop a FASTA record by exact contig name; may be supplied multiple times",
    )
    parser.add_argument(
        "--no-mask-y-par",
        dest="mask_y_par",
        action="store_false",
        help="do not mask chrY pseudoautosomal regions",
    )
    parser.set_defaults(mask_y_par=True)
    parser.add_argument(
        "--y-contig",
        default="chrY",
        help="name of the Y chromosome contig to mask",
    )
    parser.add_argument(
        "--par-interval",
        action="append",
        type=parse_interval,
        default=[],
        metavar="START-END",
        help="1-based closed chrY interval to mask; may be supplied multiple times",
    )
    parser.add_argument(
        "--no-telomere-decoy",
        action="store_true",
        help="do not append the telomeric-repeat decoy contig",
    )
    parser.add_argument(
        "--telomere-name",
        default="telomere_decoy_TTAGGG_100k",
        help="FASTA record name for the appended telomere decoy",
    )
    parser.add_argument(
        "--telomere-motif",
        default="TTAGGG",
        help="repeat motif used for the telomere decoy",
    )
    parser.add_argument(
        "--telomere-length",
        type=int,
        default=100_000,
        help="length of the telomere decoy contig in bases",
    )
    parser.add_argument(
        "--line-width",
        type=int,
        default=50,
        help="line width for the appended decoy sequence",
    )
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    if not args.par_interval:
        args.par_interval = list(DEFAULT_PAR_INTERVALS)

    summary = process_fasta(args)
    print(
        "\n".join(
            [
                f"wrote: {args.output_fasta}",
                f"input records: {summary.total_records}",
                f"kept records: {summary.kept_records}",
                f"dropped non-primary records: {summary.dropped_non_primary_records}",
                f"chrY PAR bases covered: {summary.y_par_bases_covered}",
                f"chrY PAR bases changed to N: {summary.y_par_bases_changed_to_n}",
                f"telomere decoy bases appended: {summary.telomere_decoy_bases}",
            ]
        ),
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
