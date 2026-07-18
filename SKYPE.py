import os
import filecmp
import shlex
import shutil
import psutil
import argparse
import subprocess
import re
import tempfile
from dataclasses import dataclass
from typing import Optional

from concurrent.futures import ThreadPoolExecutor

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
MEM_SAFE_RATIO = 0.8

REFERENCE_HS1 = "hs1"
REFERENCE_HG38 = "hg38"

HS1_DEPTH_DIR = "01_depth"
HG38_DEPTH_DIR = "02_depth_hg38"
HS1_ALIGNASM_DIR = "20_alignasm"
HG38_ALIGNASM_DIR = "21_alignasm_hg38"
HS1_SKYPE_DIR = "30_skype"
HG38_SKYPE_DIR = "31_skype_hg38"

# Hg38 const
HG38_URL = "https://hgdownload.soe.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz"
HG38_SOURCE_FA = "hg38.fa"
HG38_ALIGNASM_FA = "hg38.custom.primary.ypar_masked_noM.fa"
HG38_DEPTH_FA = "hg38.custom.primary.ypar_unmasked.fa"
HG38_PUBLIC_FAI = "hg38.fa.fai"
HG38_ACEN_BED = "hg38_acen.bed"
HG38_CYTOBAND_BED = "hg38_cytobands_allchrs.bed"
HG38_TELOMERE_BED = "hg38_telomere.bed"
HG38_ALIGNASM_EXCLUDE_CONTIGS = ("chrM",)

VCF_INS_MIN_SVLEN = 100_000
VCF_INS_SEQUENCE_INFO_KEYS = ("INSSEQ", "SVINSSEQ", "SEQ")


@dataclass
class ReferenceBundle:
    name: str
    alignasm_ref: str
    depth_ref: str
    chr_fai: str
    tel_bed: str
    rpt_bed: str
    rcs_bed: str
    cyt_bed: str
    ref_stat: Optional[str]


def get_reference_name(reference):
    return getattr(reference, "name", reference)


def reference_stage_dir(reference, hs1_dir, hg38_dir):
    return hg38_dir if get_reference_name(reference) == REFERENCE_HG38 else hs1_dir


def depth_dir_name(reference):
    return reference_stage_dir(reference, HS1_DEPTH_DIR, HG38_DEPTH_DIR)


def alignasm_dir_name(reference):
    return reference_stage_dir(reference, HS1_ALIGNASM_DIR, HG38_ALIGNASM_DIR)


def skype_dir_name(reference):
    return reference_stage_dir(reference, HS1_SKYPE_DIR, HG38_SKYPE_DIR)


def gfa_to_fa(gfa_file, out_fa):
    # Convert GFA format to FASTA format, extracting sequence lines.
    with open(out_fa, "w") as out_f:
        with open(gfa_file) as gfa:
            for gfa_line in gfa:
                if gfa_line.startswith("S"):
                    parts = gfa_line.strip().split("\t")
                    out_f.write(f">{parts[1]}\n{parts[2]}\n")

def normalize_extra_args(extra_args):
    if not extra_args:
        return []
    if isinstance(extra_args, str):
        return shlex.split(extra_args)
    return list(extra_args)

def sanitize_vcf_id(value):
    value = str(value) if value not in (None, "") else "unknown"
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value)

def vcf_ins_query_name(line_no, rec_id):
    return f"vcf_ins_{int(line_no)}_{sanitize_vcf_id(rec_id)}"

def parse_vcf_info(info_text):
    info = {}
    if info_text in ("", "."):
        return info
    for item in info_text.split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            info[key] = value
        else:
            info[item] = True
    return info

def parse_vcf_optional_int(value):
    if value in (None, True, ""):
        return None
    try:
        return int(float(str(value).split(",", 1)[0]))
    except ValueError:
        return None

def vcf_ins_effective_svlen(info, caller=None):
    svlen = parse_vcf_optional_int(info.get("SVLEN"))
    if svlen is None and caller == "nanomonsv":
        svlen = parse_vcf_optional_int(info.get("SVINSLEN"))
    return abs(svlen) if svlen is not None else None

def clean_vcf_sequence(value):
    if value in (None, True, ""):
        return None
    seq = str(value).split(",", 1)[0].strip().upper()
    if not seq or seq == ".":
        return None
    if any(base not in "ACGTN" for base in seq):
        return None
    return seq

def is_symbolic_or_breakend_alt(alt):
    if not alt or alt == ".":
        return True
    alt = alt.split(",", 1)[0]
    return (alt.startswith("<") and alt.endswith(">")) or "[" in alt or "]" in alt

def extract_vcf_ins_sequence(ref, alt, info):
    alt = str(alt).split(",", 1)[0] if alt is not None else ""
    if not is_symbolic_or_breakend_alt(alt):
        alt_seq = clean_vcf_sequence(alt)
        ref_seq = clean_vcf_sequence(ref)
        if alt_seq:
            if ref_seq and len(alt_seq) > len(ref_seq) and alt_seq.startswith(ref_seq):
                return alt_seq[len(ref_seq):]
            return alt_seq

    for key in VCF_INS_SEQUENCE_INFO_KEYS:
        seq = clean_vcf_sequence(info.get(key))
        if seq:
            return seq
    return None

def collect_vcf_ins_sequences(
    vcf_path,
    min_svlen=VCF_INS_MIN_SVLEN,
):
    records = []
    caller = None
    with open(vcf_path, "rt", encoding="ascii") as vcf:
        for line_no, line in enumerate(vcf, start=1):
            line = line.rstrip("\n")
            if line.startswith("##source="):
                source = line.split("=", 1)[1].strip().lower()
                if source.startswith("nanomonsv-") or source.startswith("nanomonsv"):
                    caller = "nanomonsv"
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 8:
                continue
            chrom, pos, rec_id, ref, alt, _, _, info_text = cols[:8]
            info = parse_vcf_info(info_text)
            if str(info.get("SVTYPE", "")).upper() != "INS":
                continue
            svlen = vcf_ins_effective_svlen(info, caller)
            if svlen is None or svlen < int(min_svlen):
                continue
            if rec_id in ("", "."):
                rec_id = f"VCF_RECORD_{line_no}"
            seq = extract_vcf_ins_sequence(ref, alt, info)
            if seq is None:
                continue
            records.append({
                "line_no": line_no,
                "id": rec_id,
                "name": vcf_ins_query_name(line_no, rec_id),
                "chrom": chrom,
                "pos": pos,
                "svlen": svlen,
                "sequence": seq,
            })
    return records

def write_fasta_record(out_f, name, seq, width=80):
    out_f.write(f">{name}\n")
    for st in range(0, len(seq), width):
        out_f.write(seq[st:st + width] + "\n")

def write_vcf_ins_fasta(
    vcf_path,
    out_fa,
    force=False,
    min_svlen=VCF_INS_MIN_SVLEN,
):
    records = collect_vcf_ins_sequences(vcf_path, min_svlen)
    output_dir = os.path.dirname(os.path.abspath(out_fa))
    fd, candidate_fa = tempfile.mkstemp(
        prefix=f".{os.path.basename(out_fa)}.",
        suffix=".tmp",
        dir=output_dir,
    )
    try:
        with os.fdopen(fd, "wt", encoding="ascii") as out:
            for record in records:
                write_fasta_record(out, record["name"], record["sequence"])

        changed = (
            force
            or not os.path.isfile(out_fa)
            or not filecmp.cmp(candidate_fa, out_fa, shallow=False)
        )
        if changed:
            os.replace(candidate_fa, out_fa)
        return len(records), changed
    finally:
        if os.path.exists(candidate_fa):
            os.unlink(candidate_fa)

def require_file(path, label):
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Missing {label}: {path}")

def ensure_fai(fasta_path, force):
    fai_path = f"{fasta_path}.fai"
    if force or not os.path.isfile(fai_path):
        subprocess.run(["samtools", "faidx", fasta_path], check=True)
    return fai_path

def ensure_hg38_source(dep_folder, force):
    source_fa = os.path.join(dep_folder, HG38_SOURCE_FA)
    source_gz = f"{source_fa}.gz"
    if os.path.isfile(source_fa) and not force:
        return source_fa

    if force or not os.path.isfile(source_gz):
        subprocess.run(["wget", "-nv", "-O", source_gz, HG38_URL], check=True)
    subprocess.run(["pigz", "-d", "-f", source_gz], check=True)
    return source_fa

def run_make_custom_hg38(input_fa, output_fa, target_contig_fai, force, extra_args=None):
    if os.path.isfile(output_fa) and not force:
        return output_fa

    cmd = [
        "python3", os.path.join(SCRIPT_DIR, "make_custom_hg38.py"),
        input_fa, output_fa, target_contig_fai, "--force"
    ]
    if extra_args:
        cmd.extend(extra_args)
    subprocess.run(cmd, check=True)
    return output_fa

def prepare_hg38_references(dep_folder, force):
    public_data = os.path.join(dep_folder, "SKYPE", "public_data")
    os.makedirs(public_data, exist_ok=True)
    public_fai = os.path.join(public_data, HG38_PUBLIC_FAI)
    tel_bed = os.path.join(public_data, HG38_TELOMERE_BED)
    acen_bed = os.path.join(public_data, HG38_ACEN_BED)
    cyt_bed = os.path.join(public_data, HG38_CYTOBAND_BED)
    require_file(public_fai, "hg38 FAI")
    require_file(tel_bed, "hg38 telomere BED")
    require_file(acen_bed, "hg38 acen BED")
    require_file(cyt_bed, "hg38 cytoband BED")

    alignasm_fa = os.path.join(dep_folder, HG38_ALIGNASM_FA)
    depth_fa = os.path.join(dep_folder, HG38_DEPTH_FA)

    needs_alignasm_fa = force or not os.path.isfile(alignasm_fa)
    needs_depth_fa = force or not os.path.isfile(depth_fa)

    if needs_alignasm_fa or needs_depth_fa:
        source_fa = ensure_hg38_source(dep_folder, force)
        if needs_alignasm_fa:
            exclude_args = [
                arg
                for contig in HG38_ALIGNASM_EXCLUDE_CONTIGS
                for arg in ("--exclude-contig", contig)
            ]
            run_make_custom_hg38(source_fa, alignasm_fa, public_fai, force, exclude_args)
        if needs_depth_fa:
            run_make_custom_hg38(
                source_fa, depth_fa, public_fai, force,
                ["--no-mask-y-par"]
            )

    ensure_fai(alignasm_fa, force)
    ensure_fai(depth_fa, force)

    return ReferenceBundle(
        name=REFERENCE_HG38,
        alignasm_ref=alignasm_fa,
        depth_ref=depth_fa,
        chr_fai=public_fai,
        tel_bed=tel_bed,
        rpt_bed="",
        rcs_bed=acen_bed,
        cyt_bed=cyt_bed,
        ref_stat=None,
    )

def resolve_reference_bundle(
    dep_folder,
    reference,
    force=False,
):
    dep_folder = os.path.abspath(dep_folder)
    skype_folder_loc = os.path.join(dep_folder, "SKYPE")
    public_data = os.path.join(skype_folder_loc, "public_data")

    if reference == REFERENCE_HG38:
        return prepare_hg38_references(dep_folder, force)
    if reference != REFERENCE_HS1:
        raise ValueError(f"Unsupported reference: {reference}")

    return ReferenceBundle(
        name=REFERENCE_HS1,
        alignasm_ref=os.path.join(dep_folder, 'chm13v2.0_maskedY_noM.fa'),
        depth_ref=os.path.join(dep_folder, 'chm13v2.0.fa'),
        chr_fai=os.path.join(public_data, "chm13v2.0.fa.fai"),
        tel_bed=os.path.join(public_data, "chm13v2.0_telomere.bed"),
        rpt_bed=os.path.join(public_data, "chm13v2.0_repeat.m.bed"),
        rcs_bed=os.path.join(public_data, "chm13v2.0_censat_v2.1.m.bed"),
        cyt_bed=os.path.join(public_data, "chm13v2.0_cytobands_allchrs.bed"),
        ref_stat=os.path.join(public_data, "CHM13.win.stat.gz"),
    )

def get_samtools_sort_memory_limit(thread):
    thread = max(int(thread), 1)
    memory_per_thread_mb = int(
        psutil.virtual_memory().available * MEM_SAFE_RATIO / thread / (1024 ** 2)
    )
    return f"{max(memory_per_thread_mb, 1)}M"

def index_bam_with_samtools(sorted_bam_file, thread, force):
    THREAD = str(thread)
    bai_file = f"{sorted_bam_file}.bai"

    if force or not os.path.isfile(bai_file):
        subprocess.run([
            "samtools", "index", "-@", THREAD, sorted_bam_file, bai_file
        ], check=True)

    return sorted_bam_file

def sort_sam_and_index_bam_with_samtools(sam_file, sorted_bam_file, thread, force):
    THREAD = str(thread)

    if force or not os.path.isfile(sorted_bam_file):
        subprocess.run([
            "samtools", "sort", "-@", THREAD,
            "-m", get_samtools_sort_memory_limit(THREAD),
            "-O", "BAM", "-o", sorted_bam_file, sam_file
        ], check=True)

    return index_bam_with_samtools(sorted_bam_file, THREAD, force)

def hifi_preprocess(
    CELL_LINE, PREFIX, hifi_fastq, thread, dep_folder, force, hifiasm_args,
    reference=REFERENCE_HS1,
):
    # Preprocessing pipeline for HiFi data using hifiasm.
    depth_window = 100 * 1000
    os.makedirs(PREFIX, exist_ok=True)

    if dep_folder is None:
        install_dependency(os.path.join(PREFIX, '99_dependency'), True)
        dep_folder = os.path.join(PREFIX, '99_dependency')
    reference_bundle = resolve_reference_bundle(dep_folder, reference, force)

    THREAD = str(thread)
    hifiasm_args = normalize_extra_args(hifiasm_args)
    minimap2_preset = "lr:hq" if "--ont" in hifiasm_args else "map-hifi"

    # De novo assembly
    hifiasm_folder = os.path.join(PREFIX, '00_hifiasm')
    os.makedirs(hifiasm_folder, exist_ok=True)

    gfa_loc_list = [os.path.join(hifiasm_folder, f"{CELL_LINE}.{gfa_suffix}") for gfa_suffix in ['p_ctg.gfa', 'r_utg.gfa']]
    is_file = all(map(os.path.isfile, gfa_loc_list))

    if not is_file or force:
        hifiasm_cmd = [
            'hifiasm', '-o', f'{hifiasm_folder}/{CELL_LINE}',
            '--telo-m', 'CCCTAA', '-t', THREAD] + hifi_fastq + ['--primary']
        if hifiasm_args:
            hifiasm_cmd.extend(hifiasm_args)
        subprocess.run(hifiasm_cmd, check=True)

    out_fa_list = []
    os.makedirs(os.path.join(PREFIX, "10_asm"), exist_ok=True)
    for gfa_suffix in ['p_ctg.gfa', 'r_utg.gfa']:
        gfa_file = os.path.join(hifiasm_folder, f"{CELL_LINE}.{gfa_suffix}")

        kind = gfa_suffix[0]
        out_fa = os.path.join(PREFIX, "10_asm", f"{CELL_LINE}.{kind}.fa")
        if not os.path.isfile(out_fa) or force:
            gfa_to_fa(gfa_file, out_fa)

        out_fa_list.append(out_fa)

    # Mapping read to reference
    refseq = reference_bundle.depth_ref

    depth_folder = os.path.join(PREFIX, depth_dir_name(reference_bundle))
    os.makedirs(depth_folder, exist_ok=True)
    sam_file = os.path.join(depth_folder, f'{CELL_LINE}.sam')
    sorted_bam_file = os.path.join(depth_folder, f'{CELL_LINE}.sorted.bam')

    if not os.path.isfile(os.path.join(depth_folder, f'{CELL_LINE}.win.stat.gz')) or force:
        if not os.path.isfile(sorted_bam_file) or force:
            subprocess.run([
                "minimap2", "-x", minimap2_preset, "-K", "10G", "-t", THREAD,
                "-a", refseq] + hifi_fastq + ["-o", sam_file
            ], check=True)

            depth_bam_file = sort_sam_and_index_bam_with_samtools(sam_file, sorted_bam_file, THREAD, force)
            os.remove(sam_file)
        else:
            depth_bam_file = index_bam_with_samtools(sorted_bam_file, THREAD, force)

        subprocess.run([
            os.path.join(dep_folder, 'PanDepth/bin/pandepth'), "-w", str(depth_window), "-t", THREAD,
            "-o", os.path.join(depth_folder, CELL_LINE),
            "-i", depth_bam_file
        ], check=True)

    return out_fa_list

def flye_preprocess(
    CELL_LINE, PREFIX, hifi_fastq, thread, dep_folder, force, flye_type,
    flye_args, minimap2_preset, reference=REFERENCE_HS1,
):
    # Preprocessing pipeline for long-read data using Flye.
    minimap2_preset_flye = get_minimap2_preset_from_flye(flye_type)
    if minimap2_preset_flye is None and minimap2_preset is None:
        raise AssertionError("Please use --minimap2_preset to get preset for long-read mapping at minimap2")
    
    minimap2_preset = minimap2_preset_flye if minimap2_preset is None else minimap2_preset

    depth_window = 100 * 1000
    os.makedirs(PREFIX, exist_ok=True)

    if dep_folder is None:
        install_dependency(os.path.join(PREFIX, '99_dependency'), True)
        dep_folder = os.path.join(PREFIX, '99_dependency')
    reference_bundle = resolve_reference_bundle(dep_folder, reference, force)

    THREAD = str(thread)
    flye_args = normalize_extra_args(flye_args)

    # De novo assembly
    flye_folder = os.path.join(PREFIX, '00_flye')
    os.makedirs(flye_folder, exist_ok=True)

    out_fa_list = [os.path.join(flye_folder, 'assembly.fasta'), os.path.join(flye_folder, 'assembly_graph.gfa')]
    is_file = all(map(os.path.isfile, out_fa_list))

    if not is_file or force:
        flye_cmd = [
            'flye', f'--{flye_type}'] + hifi_fastq + ['-o', f'{flye_folder}',
            '-t', THREAD]
        if flye_args:
            flye_cmd.extend(flye_args)
        subprocess.run(flye_cmd, check=True)
    
    # Mapping read to reference
    refseq = reference_bundle.depth_ref
    
    depth_folder = os.path.join(PREFIX, depth_dir_name(reference_bundle))
    os.makedirs(depth_folder, exist_ok=True)
    sam_file = os.path.join(depth_folder, f'{CELL_LINE}.sam')
    sorted_bam_file = os.path.join(depth_folder, f'{CELL_LINE}.sorted.bam')
        
    if not os.path.isfile(os.path.join(depth_folder, f'{CELL_LINE}.win.stat.gz')) or force:
        if not os.path.isfile(sorted_bam_file) or force:
            subprocess.run([
                "minimap2", "-x", minimap2_preset, "-K", "10G", "-t", THREAD,
                "-a", refseq] + hifi_fastq + ["-o", sam_file
            ], check=True)

            depth_bam_file = sort_sam_and_index_bam_with_samtools(sam_file, sorted_bam_file, THREAD, force)
            os.remove(sam_file)
        else:
            depth_bam_file = index_bam_with_samtools(sorted_bam_file, THREAD, force)

        subprocess.run([
            os.path.join(dep_folder, 'PanDepth/bin/pandepth'), "-w", str(depth_window), "-t", THREAD,
            "-o", os.path.join(depth_folder, CELL_LINE),
            "-i", depth_bam_file
        ], check=True)

    return out_fa_list

def install_dependency(dep_folder, force):
    # Install all necessary dependencies for the SKYPE pipeline.
    if force:
        if os.path.isdir(dep_folder):
            shutil.rmtree(dep_folder)
    
    os.makedirs(dep_folder, exist_ok=True)
    
    # download reference
    subprocess.run(["wget", "-nv", "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz"], cwd=dep_folder, check=True)
    subprocess.run(["wget", "-nv", "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0_maskedY.fa.gz"], cwd=dep_folder, check=True)

    subprocess.run(["pigz", "-d", "chm13v2.0.fa.gz"], cwd=dep_folder, check=True)
    subprocess.run(["pigz", "-d", "chm13v2.0_maskedY.fa.gz"], cwd=dep_folder, check=True)

    # Remove chrM
    subprocess.run("faidx chm13v2.0_maskedY.fa -g chrM --invert-match -o chm13v2.0_maskedY_noM.fa", cwd=dep_folder, shell=True, check=True)
    os.remove(os.path.join(dep_folder, 'chm13v2.0_maskedY.fa'))

    subprocess.run("git clone https://github.com/ACCtools/PanDepth && "\
                   "cd PanDepth && make", cwd=dep_folder, shell=True, check=True)
    
    subprocess.run("git clone https://github.com/ACCtools/alignasm && "\
                   "cd alignasm && "\
                   "git clone https://github.com/Microsoft/vcpkg.git && ./vcpkg/bootstrap-vcpkg.sh && "\
                   "mkdir build && cmake -B build -S . -DCMAKE_BUILD_TYPE=Release -DCMAKE_TOOLCHAIN_FILE=vcpkg/scripts/buildsystems/vcpkg.cmake && "\
                   "cmake --build build", cwd=dep_folder, shell=True, check=True)
    
    subprocess.run("git clone https://github.com/ACCtools/SKYPE", cwd=dep_folder, shell=True, check=True)

def update_dependency(dep_folder):
    # Update existing dependencies to the latest versions from their repositories.
    os.makedirs(dep_folder, exist_ok=True)

    shutil.rmtree(os.path.join(dep_folder, 'PanDepth'), ignore_errors=True)
    subprocess.run("git clone https://github.com/ACCtools/PanDepth && "\
                   "cd PanDepth && make", cwd=dep_folder, shell=True, check=True)
    
    shutil.rmtree(os.path.join(dep_folder, 'alignasm'), ignore_errors=True)
    subprocess.run("git clone https://github.com/ACCtools/alignasm && "\
                   "cd alignasm && "\
                   "git clone https://github.com/Microsoft/vcpkg.git && ./vcpkg/bootstrap-vcpkg.sh && "\
                   "mkdir build && cmake -B build -S . -DCMAKE_BUILD_TYPE=Release -DCMAKE_TOOLCHAIN_FILE=vcpkg/scripts/buildsystems/vcpkg.cmake && "\
                   "cmake --build build", cwd=dep_folder, shell=True, check=True)
    
    shutil.rmtree(os.path.join(dep_folder, 'SKYPE'), ignore_errors=True)
    subprocess.run("git clone https://github.com/ACCtools/SKYPE", cwd=dep_folder, shell=True, check=True)


def is_nonempty_file(path):
    return os.path.isfile(path) and os.path.getsize(path) > 0

def run_alignasm(PREFIX_PATH, thread, fa_loc, ref_loc, ALIGNASM_LOC, force):
    # Run the alignasm tool for sequence alignment.
    THREAD = str(thread)

    paf_file = f"{PREFIX_PATH}.paf"
    pat_fa_file = f"{PREFIX_PATH}.pat.fa"
    alt_paf_file = f"{PREFIX_PATH}.alt.paf"
    tar_paf_tuple = (paf_file, f"{PREFIX_PATH}.aln.paf")
    is_file = all(map(os.path.isfile, tar_paf_tuple))

    if not is_file or force:
        if force or not os.path.isfile(paf_file):
            subprocess.run([
                "minimap2", "--cs", "-t", THREAD, "-x", "asm20",
                "--no-long-join", "-r2k", "-K10G",
                ref_loc, fa_loc, "-o", paf_file
            ], check=True)

        if force or not os.path.isfile(alt_paf_file):
            if is_nonempty_file(paf_file):
                subprocess.run([
                    "python3", os.path.join(SCRIPT_DIR, "paf_gap_seq.py"),
                    fa_loc, paf_file, pat_fa_file
                ], check=True)

                subprocess.run([
                    "minimap2", "--cs", "-t", THREAD, "-x", "asm20",
                    "-r2k", "-K10G",
                    ref_loc, pat_fa_file, "-o", alt_paf_file
                ], check=True)
            else:
                open(pat_fa_file, "wt").close()
                open(alt_paf_file, "wt").close()

        alignasm_cmd = [
            ALIGNASM_LOC, paf_file,
            "-t", THREAD, "--non_skip_linkable"
        ]
        if is_nonempty_file(alt_paf_file):
            alignasm_cmd.extend(["-a", alt_paf_file])
        subprocess.run(alignasm_cmd, check=True)

    return tar_paf_tuple

def subprocess_print(args, **kwargs):
    print(*args)

def run_skype(CELL_LINE, PREFIX, ctg_paf, ctg_aln_paf, utg_paf, utg_aln_paf,
              depth_loc, thread, dep_folder, is_progress, skype_force, graph_depth,
              option_02="", skype_start_at=0, print_args=False,
              benchmark_vcf_loc=None, reference_bundle=None, vcf_ins_aln_paf=None):
    # Execute the core SKYPE analysis scripts.
    dep_folder = os.path.abspath(dep_folder)
    skype_folder_loc = os.path.join(dep_folder, "SKYPE")
    if reference_bundle is None:
        reference_bundle = resolve_reference_bundle(dep_folder, REFERENCE_HS1)

    TEL_BED = reference_bundle.tel_bed
    CHR_FAI = reference_bundle.chr_fai
    RPT_BED = reference_bundle.rpt_bed
    RCS_BED = reference_bundle.rcs_bed
    CYT_BED = reference_bundle.cyt_bed
    REF_STAT_LOC = reference_bundle.ref_stat


    # python to bash variable
    MAIN_STAT_LOC = depth_loc
    MAIN_STAT_NORM_LOC = depth_loc if REF_STAT_LOC is None else depth_loc.replace('.win.stat.gz', '_normalized.win.stat.gz')
    PAF_LOC = ctg_aln_paf
    PAF_UTG_LOC = utg_aln_paf
    graph_paf_loc = PAF_UTG_LOC if benchmark_vcf_loc else PAF_LOC
    PPC_PAF_LOC = os.path.join(PREFIX, f"{os.path.basename(graph_paf_loc)}.ppc.paf")
    READ_BAM_LOC = os.path.join(os.path.dirname(depth_loc), f'{CELL_LINE}.bam')
    SORTED_READ_BAM_LOC = os.path.join(os.path.dirname(depth_loc), f'{CELL_LINE}.sorted.bam')
    if os.path.isfile(SORTED_READ_BAM_LOC):
        READ_BAM_LOC = SORTED_READ_BAM_LOC

    THREAD = str(thread)
    PROGRESS = ["--progress"] if is_progress else []

    subprocess_run = subprocess_print if print_args else subprocess.run

    if not os.path.isfile(os.path.join(PREFIX, 'virtual_sky.png')) or skype_force or skype_start_at > 0:
        if skype_start_at <= 0 and REF_STAT_LOC is not None:
            subprocess_run([
                "python", os.path.join(skype_folder_loc, "00_depth_norm.py"),
                MAIN_STAT_LOC, REF_STAT_LOC, RCS_BED
            ], check=True)

        if skype_start_at <= 2:
            EXTRA_02 = shlex.split(option_02) if option_02 else []
            build_graph_cmd = [
                "python", os.path.join(skype_folder_loc, "02_Build_Breakend_Graph_Limited.py"),
                graph_paf_loc, CHR_FAI, TEL_BED, RPT_BED, RCS_BED, MAIN_STAT_NORM_LOC, PREFIX, READ_BAM_LOC,
                "-t", THREAD, "-d", str(graph_depth),
            ]
            if benchmark_vcf_loc:
                build_graph_cmd.extend(["--vcf_input", os.path.abspath(benchmark_vcf_loc)])
                if vcf_ins_aln_paf is not None:
                    build_graph_cmd.extend(["--alt", os.path.abspath(vcf_ins_aln_paf)])
            else:
                build_graph_cmd.extend(["--alt", PAF_UTG_LOC, "--original_paf_loc", ctg_paf, utg_paf])
            subprocess_run(build_graph_cmd + EXTRA_02 + PROGRESS, check=True)

        if skype_start_at <= 11:
            subprocess_run([
                "python", os.path.join(skype_folder_loc, "11_Ref_Outlier_Contig_Modify.py"),
                CHR_FAI, PPC_PAF_LOC, PREFIX,
            ], check=True)

        if skype_start_at <= 21:
            free_mem_gb = psutil.virtual_memory().available * MEM_SAFE_RATIO / (1024 ** 3)
            thread_lim = int(free_mem_gb / 6)

            subprocess_run([
                "python", os.path.join(skype_folder_loc, "21_run_depth.py"),
                PPC_PAF_LOC, PREFIX,
                "--pandepth_loc", os.path.join(dep_folder, 'PanDepth', 'bin', 'pandepth'),
                "-t", str(max(min(thread_lim, thread), 1))
            ] + PROGRESS, check=True)

        if skype_start_at <= 22:
            subprocess_run([
                "python", os.path.join(skype_folder_loc, "22_save_matrix.py"),
                RCS_BED, PPC_PAF_LOC, MAIN_STAT_NORM_LOC,
                TEL_BED, CHR_FAI, CYT_BED, PREFIX, "-t", THREAD
            ] + PROGRESS, check=True)

        if skype_start_at <= 23:
            subprocess_run([
                "python", "23_run_nnls.py", PPC_PAF_LOC,
                os.path.abspath(PREFIX), MAIN_STAT_NORM_LOC, RCS_BED, PAF_UTG_LOC, "-t", THREAD
            ], check=True, cwd=skype_folder_loc)

        if skype_start_at <= 24:
            core_num = psutil.cpu_count(logical=False)
            if core_num is None:
                JULIA_THREAD = THREAD
            else:
                JULIA_THREAD = min(int(THREAD), core_num)

            subprocess_run([
                "python", "-X", f"juliacall-threads={THREAD}", "-X", "juliacall-handle-signals=yes",
                "24_cluster_weight.py", PPC_PAF_LOC, MAIN_STAT_NORM_LOC,
                TEL_BED, CHR_FAI, os.path.abspath(PREFIX), "-t", str(JULIA_THREAD)
            ], check=True, cwd=skype_folder_loc)

        if skype_start_at <= 30:
            subprocess_run([
                "python", os.path.join(skype_folder_loc, "30_virtual_sky.py"),
                PPC_PAF_LOC, MAIN_STAT_NORM_LOC,
                TEL_BED, CHR_FAI, PREFIX, CELL_LINE
            ], check=True)

        if skype_start_at <= 31:
            subprocess_run([
                "python", os.path.join(skype_folder_loc, "31_depth_analysis.py"),
                RCS_BED, PPC_PAF_LOC, MAIN_STAT_NORM_LOC,
                TEL_BED, CHR_FAI, CYT_BED, PREFIX, "-t", THREAD
            ] + PROGRESS, check=True)
        

def analysis(CELL_LINE, PREFIX, contig_loc, unitig_loc, depth_loc, thread, dep_folder,
             is_progress, force, skype_force, run_skype_func, graph_depth,
             no_utg=False, skype_dir=None, option_02="", skype_start_at=0,
             print_args=False, reference=REFERENCE_HS1, benchmark_vcf_loc=None):
    # Main analysis function that orchestrates alignment and SKYPE execution.
    os.makedirs(PREFIX, exist_ok=True)

    if dep_folder is None:
        install_dependency(os.path.join(PREFIX, '99_dependency'), True)
        dep_folder = os.path.join(PREFIX, '99_dependency')
    reference_bundle = resolve_reference_bundle(dep_folder, reference, force)

    if skype_dir is None:
        skype_dir = os.path.join(PREFIX, skype_dir_name(reference_bundle))
    
    alignasm_ref_loc = reference_bundle.alignasm_ref
    alignasm_loc = os.path.join(dep_folder, 'alignasm', 'build', 'alignasm')
    alignasm_dir = os.path.join(PREFIX, alignasm_dir_name(reference_bundle))
    os.makedirs(alignasm_dir, exist_ok=True)

    alignasm_folder_ctg = os.path.join(alignasm_dir, f'{CELL_LINE}.ctg')
    alignasm_folder_utg = os.path.join(alignasm_dir, f'{CELL_LINE}.utg')
    alignasm_folder_vcf = os.path.join(alignasm_dir, f'{CELL_LINE}.vcf')

    if no_utg:
        alignasm_thread = thread
        alignasm_worker_thread = 1
    else:
        if thread > 1:
            alignasm_thread = thread // 2
            alignasm_worker_thread = 2
        else:
            alignasm_thread = 1
            alignasm_worker_thread = 1

    args_ctg = (
        alignasm_folder_ctg, alignasm_thread, contig_loc, alignasm_ref_loc,
        alignasm_loc, force
    )
    args_utg = (
        alignasm_folder_utg, alignasm_thread, unitig_loc, alignasm_ref_loc,
        alignasm_loc, force
    )

    if benchmark_vcf_loc:
        if no_utg:
            ctg_paf, ctg_aln_paf = run_alignasm(
                alignasm_folder_ctg, thread, contig_loc, alignasm_ref_loc,
                alignasm_loc, force
            )
            utg_paf, utg_aln_paf = ctg_paf, ctg_aln_paf
        else:
            utg_paf, utg_aln_paf = run_alignasm(
                alignasm_folder_utg, thread, unitig_loc, alignasm_ref_loc,
                alignasm_loc, force
            )
            ctg_paf, ctg_aln_paf = utg_paf, utg_aln_paf

        vcf_ins_aln_paf = None
        vcf_ins_fa = os.path.join(alignasm_dir, f"{CELL_LINE}.vcf.ins.fa")
        vcf_ins_count, vcf_ins_fasta_changed = write_vcf_ins_fasta(
            benchmark_vcf_loc, vcf_ins_fa, force=force
        )
        if vcf_ins_count > 0:
            vcf_ins_aln_paf = f"{alignasm_folder_vcf}.aln.paf"
            if vcf_ins_fasta_changed or not os.path.isfile(vcf_ins_aln_paf):
                _, vcf_ins_aln_paf = run_alignasm(
                    alignasm_folder_vcf, thread, vcf_ins_fa, alignasm_ref_loc,
                    alignasm_loc, force or vcf_ins_fasta_changed
                )

        effective_skype_force = skype_force or vcf_ins_fasta_changed
        depth_loc = os.path.abspath(depth_loc)
        return run_skype_func(
            CELL_LINE, os.path.abspath(skype_dir), ctg_paf, ctg_aln_paf,
            utg_paf, utg_aln_paf, depth_loc, thread, dep_folder, is_progress,
            effective_skype_force, graph_depth, option_02=option_02,
            skype_start_at=skype_start_at, print_args=print_args,
            benchmark_vcf_loc=benchmark_vcf_loc, reference_bundle=reference_bundle,
            vcf_ins_aln_paf=vcf_ins_aln_paf
        )
    
    if no_utg:
        with ThreadPoolExecutor(max_workers=alignasm_worker_thread) as executor:
            future_ctg = executor.submit(run_alignasm, *args_ctg)
            ctg_paf, ctg_aln_paf = future_ctg.result()

        depth_loc = os.path.abspath(depth_loc)
        return run_skype_func(
            CELL_LINE, os.path.abspath(skype_dir), ctg_paf, ctg_aln_paf,
            ctg_paf, ctg_aln_paf, depth_loc, thread, dep_folder, is_progress,
            skype_force, graph_depth, option_02=option_02,
            skype_start_at=skype_start_at, print_args=print_args,
            benchmark_vcf_loc=benchmark_vcf_loc, reference_bundle=reference_bundle
        )

    else:
        with ThreadPoolExecutor(max_workers=alignasm_worker_thread) as executor:
            future_ctg = executor.submit(run_alignasm, *args_ctg)
            future_utg = executor.submit(run_alignasm, *args_utg)

            ctg_paf, ctg_aln_paf = future_ctg.result()
            utg_paf, utg_aln_paf = future_utg.result()

        depth_loc = os.path.abspath(depth_loc)
        return run_skype_func(
            CELL_LINE, os.path.abspath(skype_dir), ctg_paf, ctg_aln_paf,
            utg_paf, utg_aln_paf, depth_loc, thread, dep_folder, is_progress,
            skype_force, graph_depth, option_02=option_02,
            skype_start_at=skype_start_at, print_args=print_args,
            benchmark_vcf_loc=benchmark_vcf_loc, reference_bundle=reference_bundle
        )

def get_skype_parser():
    # Set up the command-line argument parser.
    parser = argparse.ArgumentParser(description="SKYPE pipeline")

    # Define parser
    subparsers = parser.add_subparsers(dest="command", required=True)

    def add_work_dir_arg(subparser):
        subparser.add_argument("WORK_DIR", type=str, help="Working directory for pipeline")

    def add_hifi_fastq_args(subparser):
        subparser.add_argument(
            "HIFI_FASTQ", type=str, help="Location of HiFi fastq file",
            nargs=argparse.REMAINDER
        )

    def add_flye_input_args(subparser, fastq_help="Location of long-read fastq file(s)"):
        subparser.add_argument(
            "FLYE_TYPE", type=str,
            help="Flye assembly read type (ex: pacbio-raw, nano-raw). See flye CLI docs."
        )
        subparser.add_argument(
            "LONG_READ_FASTQ", type=str,
            help=fastq_help, nargs="+"
        )

    def add_common_pipeline_args(subparser, include_progress=False):
        subparser.add_argument("-t", "--thread", type=int, help="Number of thread", default=1)
        subparser.add_argument("-p", "--prefix", type=str, help="Prefix for pipeline", default="SKYPE")
        if include_progress:
            subparser.add_argument("--progress", help="Show progress bar", action='store_true')
        subparser.add_argument(
            "--dependency_loc", type=str,
            help="SKYPE dependency folder location. If no value is specified, it will be installed automatically"
        )
        subparser.add_argument(
            "--preprocess_force", help="Don't trust previous file for preprocess",
            action='store_true'
        )

    def add_hifiasm_args(subparser):
        subparser.add_argument(
            "--hifiasm_args", type=str, default="",
            help='Custom hifiasm args (single quoted string, e.g. --hifiasm_args="--ont --chem-c 0")'
        )

    def add_flye_option_args(subparser):
        subparser.add_argument(
            "--flye_args", type=str, default="",
            help='Custom flye args (single quoted string, e.g. --flye_args="--min-overlap 5000")'
        )
        subparser.add_argument("--minimap2_preset", type=str, help="Minimap2 preset for long-read mapping")

    def add_analysis_input_args(subparser):
        subparser.add_argument("CONTIG", type=str, help="Contig fasta location")
        subparser.add_argument("UNITIG", type=str, help="UNITIG fasta location")
        subparser.add_argument("DEPTH_LOC", type=str, help="Depth result for Pandepth location")

    def add_skype_args(subparser, include_skype_dir=False, include_start_controls=False):
        subparser.add_argument("--skype_force", help="Don't trust previous file for SKYPE pipeline", action='store_true')
        if include_skype_dir:
            subparser.add_argument("--skype_dir", type=str, help="Output directory for SKYPE analysis")
        subparser.add_argument("-d", "--graph_depth", help="Depth of breakend graph", type=int, default=4)
        subparser.add_argument(
            "--option_02", type=str, default="",
            help='Extra args forwarded to 02_Build_Breakend_Graph_Limited.py, e.g. "--variant_mode"'
        )
        if include_start_controls:
            subparser.add_argument(
                "--skype_start_at", type=int, default=0,
                help='Start SKYPE pipeline from the given stage number (e.g., 23 to start at 23_run_nnls.py). Stages with smaller numbers are skipped.'
            )
            subparser.add_argument(
                "--print_args", help="Print SKYPE subprocess commands instead of executing them",
                action='store_true'
            )

    def add_reference_args(subparser):
        subparser.add_argument(
            "--reference", choices=[REFERENCE_HS1, REFERENCE_HG38],
            default=REFERENCE_HS1,
            help="Reference build for alignment/depth/SKYPE resources"
        )
        subparser.add_argument(
            "--benchmark_vcf_loc", type=str,
            help="Benchmark VCF input; enables SKYPE VCF input mode"
        )
    

    parser_hifi_prepro = subparsers.add_parser(
        "preprocess_hifi",
        help="Hifi preprocessing for SKYPE pipeline",
        usage="%(prog)s [options] WORK_DIR HIFI_FASTQ [HIFI_FASTQ ...]",
    )

    add_work_dir_arg(parser_hifi_prepro)
    add_hifi_fastq_args(parser_hifi_prepro)
    add_common_pipeline_args(parser_hifi_prepro)
    add_hifiasm_args(parser_hifi_prepro)


    parser_flye_prepro = subparsers.add_parser("preprocess_flye", help="Flye preprocessing for SKYPE pipeline")

    add_work_dir_arg(parser_flye_prepro)
    add_flye_input_args(parser_flye_prepro, fastq_help="Location of long-read fastq file")
    add_common_pipeline_args(parser_flye_prepro)
    add_flye_option_args(parser_flye_prepro)
    

    parser_anl = subparsers.add_parser("analysis", help="Analysis cancer cell line karyotype using contig")

    add_work_dir_arg(parser_anl)
    add_analysis_input_args(parser_anl)
    add_common_pipeline_args(parser_anl, include_progress=True)
    add_skype_args(parser_anl, include_skype_dir=True, include_start_controls=True)
    add_reference_args(parser_anl)

    parser_dep = subparsers.add_parser("install_dependency", help="Hifi preprocessing for SKYPE pipeline")

    parser_dep.add_argument("dependency_loc", type=str, help="SKYPE dependency folder location")

    parser_dep.add_argument("--force", help="Reinstall dependency folder", action='store_true')

    
    parser_up_dep = subparsers.add_parser("update_dependency", help="Hifi preprocessing for SKYPE pipeline")

    parser_up_dep.add_argument("dependency_loc", type=str, help="SKYPE dependency folder location")
    
    
    parser_run = subparsers.add_parser(
        "run_hifi",
        help="Pipeline for cancer hifi sequencing data",
        usage="%(prog)s [options] WORK_DIR HIFI_FASTQ [HIFI_FASTQ ...]",
    )

    add_work_dir_arg(parser_run)
    add_hifi_fastq_args(parser_run)
    add_common_pipeline_args(parser_run, include_progress=True)
    add_skype_args(parser_run)
    add_hifiasm_args(parser_run)
    add_reference_args(parser_run)

    
    parser_run_flye = subparsers.add_parser("run_flye", help="Pipeline for cancer long-read sequencing data using Flye")

    add_work_dir_arg(parser_run_flye)
    add_flye_input_args(parser_run_flye)
    add_common_pipeline_args(parser_run_flye, include_progress=True)
    add_skype_args(parser_run_flye)
    add_flye_option_args(parser_run_flye)
    add_reference_args(parser_run_flye)

    return parser

def get_minimap2_preset_from_flye(FLYE_TYPE : str):
    # Determine the appropriate minimap2 preset based on the Flye read type.
    if FLYE_TYPE == 'pacbio-hifi':
        return 'map-hifi'
    elif FLYE_TYPE.startswith('pacbio'):
        return 'map-pb'
    elif FLYE_TYPE == 'nano-raw':
        return 'map-ont'
    elif FLYE_TYPE.startswith('nano'):
        return 'lr:hq'
    else:
        return None

def main():
    # Main function to parse arguments and execute the corresponding command.
    parser = get_skype_parser()
    args = parser.parse_args()

    if args.command == "preprocess_hifi":
        if not args.HIFI_FASTQ:
            parser.error("preprocess_hifi: at least one HIFI_FASTQ is required")
        hifi_preprocess(args.prefix, args.WORK_DIR, args.HIFI_FASTQ, args.thread, args.dependency_loc, args.preprocess_force, args.hifiasm_args)
    elif args.command == "install_dependency":
        install_dependency(args.dependency_loc, args.force)
    elif args.command == "update_dependency":
        update_dependency(args.dependency_loc)
    elif args.command == "analysis":
        analysis(
            args.prefix, args.WORK_DIR, args.CONTIG, args.UNITIG, args.DEPTH_LOC,
            args.thread, args.dependency_loc, args.progress, args.preprocess_force,
            args.skype_force, run_skype, args.graph_depth, skype_dir=args.skype_dir,
            option_02=args.option_02, skype_start_at=args.skype_start_at,
            print_args=args.print_args, reference=args.reference,
            benchmark_vcf_loc=args.benchmark_vcf_loc
        )
    elif args.command == 'run_hifi':
        if not args.HIFI_FASTQ:
            parser.error("run_hifi: at least one HIFI_FASTQ is required")
        ctg_loc, utg_loc = hifi_preprocess(
            args.prefix, args.WORK_DIR, args.HIFI_FASTQ, args.thread,
            args.dependency_loc, args.preprocess_force, args.hifiasm_args,
            reference=args.reference
        )
        
        depth_loc = os.path.join(args.WORK_DIR, depth_dir_name(args.reference), f'{args.prefix}.win.stat.gz')
        if args.dependency_loc is None:
            dep_folder = os.path.join(args.WORK_DIR, '99_dependency')
        else:
            dep_folder = args.dependency_loc
        analysis(
            args.prefix, args.WORK_DIR, ctg_loc, utg_loc, depth_loc, args.thread,
            dep_folder, args.progress, args.preprocess_force, args.skype_force,
            run_skype, args.graph_depth, option_02=args.option_02,
            reference=args.reference, benchmark_vcf_loc=args.benchmark_vcf_loc
        )
    elif args.command == 'preprocess_flye':
        flye_preprocess(args.prefix, args.WORK_DIR, args.LONG_READ_FASTQ, args.thread, args.dependency_loc, args.preprocess_force, args.FLYE_TYPE, args.flye_args, args.minimap2_preset)
    elif args.command == 'run_flye':
        ctg_loc, utg_loc = flye_preprocess(
            args.prefix, args.WORK_DIR, args.LONG_READ_FASTQ, args.thread,
            args.dependency_loc, args.preprocess_force, args.FLYE_TYPE,
            args.flye_args, args.minimap2_preset, reference=args.reference
        )
        depth_loc = os.path.join(args.WORK_DIR, depth_dir_name(args.reference), f'{args.prefix}.win.stat.gz')

        if args.dependency_loc is None:
            dep_folder = os.path.join(args.WORK_DIR, '99_dependency')
        else:
            dep_folder = args.dependency_loc

        analysis(
            args.prefix, args.WORK_DIR, ctg_loc, utg_loc, depth_loc, args.thread,
            dep_folder, args.progress, args.preprocess_force, args.skype_force,
            run_skype, args.graph_depth, no_utg=True, option_02=args.option_02,
            reference=args.reference, benchmark_vcf_loc=args.benchmark_vcf_loc
        )


if __name__ == "__main__":
    try:
        main()
    except FileNotFoundError as exc:
        raise SystemExit(f"{os.path.basename(__file__)}: error: {exc}") from None
