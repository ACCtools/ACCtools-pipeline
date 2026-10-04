"""Build immutable, content-bound hs1 alignment references from the full FASTA."""
import argparse
from contextlib import contextmanager
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import tempfile


RECIPE_SCHEMA = "ACCtools.hs1_reference_recipe.v1"
MANIFEST_SCHEMA = "ACCtools.hs1_reference_manifest.v1"
DEFAULT_RECIPE = Path(__file__).with_name("hs1_reference_recipe.json")
FASTA_NAME = "chm13v2.0.ypar_masked_noM.fa"
FAI_NAME = FASTA_NAME + ".fai"
MANIFEST_NAME = "reference_manifest.json"
LINE_BASES = 80


def canonical_json(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")


def file_binding(path):
    """Hash bytes and reject a file replaced or modified during the read."""
    path = Path(path).resolve(strict=True)
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        before = os.fstat(handle.fileno())
        for chunk in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(chunk)
        after = os.fstat(handle.fileno())
    current = path.stat()
    fields = ("st_dev", "st_ino", "st_size", "st_mtime_ns", "st_ctime_ns")
    if any(len({getattr(stat, field) for stat in (before, after, current)}) != 1
           for field in fields):
        raise RuntimeError(f"Input changed while hashing: {path}")
    return {"path": str(path), "size": before.st_size, "sha256": digest.hexdigest()}


def read_fai_lengths(path):
    lengths = {}
    with Path(path).open(encoding="ascii") as handle:
        for row, line in enumerate(handle, 1):
            values = line.rstrip("\r\n").split("\t")
            if (len(values) < 2 or not values[0] or
                    any(character.isspace() for character in values[0])):
                raise ValueError(f"Malformed FAI row {row}: {path}")
            name = values[0]
            if name in lengths:
                raise ValueError(f"Duplicate FAI contig: {name}")
            length = int(values[1])
            if length < 1:
                raise ValueError(f"Invalid FAI length for {name}")
            lengths[name] = length
    if not lengths:
        raise ValueError("The model FAI is empty")
    return lengths


def validate_recipe(recipe, lengths, source_sha256):
    if not isinstance(recipe, dict) or recipe.get("schema") != RECIPE_SCHEMA:
        raise ValueError("Unsupported reference recipe schema")
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", recipe.get("name", "")):
        raise ValueError("Invalid reference recipe name")
    if recipe.get("coordinate_system") != "0-based-half-open":
        raise ValueError("Mask coordinates must be 0-based-half-open")
    if recipe.get("source_fasta_sha256") != source_sha256:
        raise ValueError("Source FASTA content does not match the reference recipe")
    required = recipe.get("required_contig_lengths")
    if not isinstance(required, dict):
        raise ValueError("Required contig lengths must be explicit")
    for name, length in required.items():
        if lengths.get(name) != length:
            raise ValueError(f"Required contig length mismatch: {name}")
    excluded = recipe.get("exclude_contigs")
    if (not isinstance(excluded, list) or
            any(not isinstance(name, str) for name in excluded) or
            len(set(excluded)) != len(excluded)):
        raise ValueError("Excluded contigs must be a list without duplicates")
    if set(excluded) - set(lengths):
        raise ValueError("An excluded contig is absent from the model FAI")
    if len(excluded) == len(lengths):
        raise ValueError("The recipe excludes every contig")
    masks = recipe.get("masks")
    if not isinstance(masks, dict):
        raise ValueError("Masks must be a contig-to-interval mapping")
    for name, intervals in masks.items():
        if name not in lengths or name in excluded:
            raise ValueError(f"Mask contig is absent or excluded: {name}")
        if not isinstance(intervals, list):
            raise ValueError(f"Mask intervals must be a list: {name}")
        previous_end = 0
        for interval in intervals:
            if (not isinstance(interval, list) or len(interval) != 2 or
                    any(type(value) is not int for value in interval)):
                raise ValueError(f"Invalid mask interval for {name}")
            start, end = interval
            if start < previous_end or not 0 <= start < end <= lengths[name]:
                raise ValueError(f"Overlapping, unsorted, or out-of-bounds mask: {name}")
            previous_end = end


def iter_fasta(path, expected_lengths):
    """Hold one contig at a time; reject malformed records before certification."""
    name, parts, length = None, [], 0
    seen = set()
    with Path(path).open("rb") as handle:
        for row, raw in enumerate(handle, 1):
            line = raw.rstrip(b"\r\n")
            if not line:
                continue
            if line.startswith(b">"):
                if name is not None:
                    if length != expected_lengths[name]:
                        raise ValueError(f"Source/model length mismatch: {name}")
                    sequence = b"".join(parts)
                    parts = []
                    yield name, sequence
                fields = line[1:].split()
                if not fields:
                    raise ValueError(f"Empty FASTA header at row {row}")
                name = fields[0].decode("ascii")
                if name in seen or name not in expected_lengths:
                    raise ValueError(f"Duplicate or unexpected FASTA contig: {name}")
                seen.add(name)
                length = 0
            else:
                if name is None:
                    raise ValueError(f"Sequence precedes FASTA header at row {row}")
                length += len(line)
                if length > expected_lengths[name]:
                    raise ValueError(f"Source/model length mismatch: {name}")
                parts.append(line)
    if name is not None:
        if length != expected_lengths[name]:
            raise ValueError(f"Source/model length mismatch: {name}")
        sequence = b"".join(parts)
        parts = []
        yield name, sequence
    if seen != set(expected_lengths):
        raise ValueError("Source/model contig dictionary mismatch")


def write_fasta_record(handle, name, sequence):
    handle.write(b">" + name.encode("ascii") + b"\n")
    offset = handle.tell()
    for start in range(0, len(sequence), LINE_BASES * 16384):
        block = sequence[start:start + LINE_BASES * 16384]
        handle.write(b"\n".join(block[i:i + LINE_BASES]
                                 for i in range(0, len(block), LINE_BASES)) + b"\n")
    bases = min(LINE_BASES, len(sequence))
    return f"{name}\t{len(sequence)}\t{offset}\t{bases}\t{bases + 1}\n"


def output_binding(path):
    result = file_binding(path)
    result["file"] = Path(path).name
    del result["path"]
    return result


def input_bindings(source_fasta, model_fai, recipe_path):
    return {"source_fasta": file_binding(source_fasta),
            "model_fai": file_binding(model_fai),
            "recipe": file_binding(recipe_path),
            "builder": file_binding(__file__)}


@contextmanager
def build_lock(path):
    with path.open("a") as handle:
        fcntl.flock(handle, fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(handle, fcntl.LOCK_UN)


def verify_cache(directory, bindings, recipe, lengths, key, namespace):
    if directory.is_symlink():
        raise RuntimeError("Refusing a symlinked reference version directory")
    try:
        manifest = json.loads((directory / MANIFEST_NAME).read_text())
        expected_lengths = {name: length for name, length in lengths.items()
                            if name not in recipe["exclude_contigs"]}
        if (manifest.get("schema") != MANIFEST_SCHEMA or
                manifest.get("status") != "complete" or
                manifest.get("inputs") != bindings or
                manifest.get("recipe") != recipe or
                manifest.get("content_key") != key or
                manifest.get("cache_namespace") != namespace or
                manifest.get("retained_contig_lengths") != expected_lengths):
            raise ValueError("Reference manifest input/recipe mismatch")
        for label, filename in (("fasta", FASTA_NAME), ("fai", FAI_NAME)):
            path = directory / filename
            if path.is_symlink() or manifest["outputs"][label] != output_binding(path):
                raise ValueError(f"Reference {label} content mismatch")
        if read_fai_lengths(directory / FAI_NAME) != expected_lengths:
            raise ValueError("Reference FAI dictionary mismatch")
    except (OSError, ValueError, KeyError, TypeError) as exc:
        raise RuntimeError(
            f"Existing reference version is invalid and was preserved: {directory}"
        ) from exc
    return {"fasta": str(directory / FASTA_NAME), "fai": str(directory / FAI_NAME),
            "manifest": str(directory / MANIFEST_NAME), "cache_namespace": namespace}


def prepare_reference(source_fasta, model_fai, output_root, recipe_path=DEFAULT_RECIPE):
    """Publish a new version atomically or verify an existing version without writes.

    There is deliberately no force/overwrite option. Input bytes, recipe, builder,
    and model FAI bind the version and all downstream alignment cache namespaces.
    """
    bindings = input_bindings(source_fasta, model_fai, recipe_path)
    recipe = json.loads(Path(recipe_path).read_text())
    lengths = read_fai_lengths(model_fai)
    validate_recipe(recipe, lengths, bindings["source_fasta"]["sha256"])
    key = hashlib.sha256(canonical_json(bindings)).hexdigest()
    namespace = recipe["name"] + "." + key[:24]
    root = Path(output_root).resolve() / recipe["name"]
    root.mkdir(parents=True, exist_ok=True)
    directory = root / key
    with build_lock(root / (key + ".lock")):
        if not directory.exists() and not directory.is_symlink():
            temporary = Path(tempfile.mkdtemp(prefix=".building.", dir=root))
            try:
                contigs, observed = [], {}
                with (temporary / FASTA_NAME).open("wb") as fasta, \
                        (temporary / FAI_NAME).open("w", encoding="ascii") as fai:
                    for name, sequence in iter_fasta(source_fasta, lengths):
                        if re.search(rb"[^ACGTRYSWKMBDHVNacgtryswkmbdhvn]", sequence):
                            raise ValueError(f"Invalid DNA sequence in {name}")
                        observed[name] = len(sequence)
                        excluded = name in recipe["exclude_contigs"]
                        source_digest = hashlib.sha256(sequence).hexdigest()
                        intervals = recipe.get("masks", {}).get(name, [])
                        if intervals:
                            sequence = bytearray(sequence)
                            for start, end in intervals:
                                sequence[start:end] = b"N" * (end - start)
                        if not excluded:
                            fai.write(write_fasta_record(fasta, name, sequence))
                        contigs.append({
                            "name": name, "length": len(sequence), "excluded": excluded,
                            "source_sequence_sha256": source_digest,
                            "masked_base_count": sum(end - start for start, end in intervals),
                            "output_sequence_sha256": None if excluded else
                                hashlib.sha256(sequence).hexdigest(),
                        })
                if list(observed.items()) != list(lengths.items()):
                    raise ValueError("Source/model contig dictionary or order mismatch")
                if input_bindings(source_fasta, model_fai, recipe_path) != bindings:
                    raise RuntimeError("Reference inputs changed during construction")
                manifest = {
                    "schema": MANIFEST_SCHEMA, "status": "complete",
                    "created_utc": datetime.now(timezone.utc).isoformat(),
                    "content_key": key, "cache_namespace": namespace,
                    "inputs": bindings, "recipe": recipe, "contigs": contigs,
                    "retained_contig_lengths": {name: length for name, length in lengths.items()
                                                if name not in recipe["exclude_contigs"]},
                    "outputs": {"fasta": output_binding(temporary / FASTA_NAME),
                                "fai": output_binding(temporary / FAI_NAME)},
                }
                (temporary / MANIFEST_NAME).write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
                # Cooperating builders use the same lock. Never replace a published version.
                if directory.exists() or directory.is_symlink():
                    raise RuntimeError("Reference version appeared during construction")
                temporary.rename(directory)
            finally:
                if temporary.exists():
                    shutil.rmtree(temporary)
        result = verify_cache(directory, bindings, recipe, lengths, key, namespace)
        if input_bindings(source_fasta, model_fai, recipe_path) != bindings:
            raise RuntimeError("Reference inputs changed during verification")
        return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-fasta", required=True)
    parser.add_argument("--model-fai", required=True)
    parser.add_argument("--output-root", required=True)
    parser.add_argument("--recipe", default=str(DEFAULT_RECIPE))
    args = parser.parse_args()
    print(json.dumps(prepare_reference(args.source_fasta, args.model_fai,
                                       args.output_root, args.recipe), indent=2))


if __name__ == "__main__":
    main()
