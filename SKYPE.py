import os
import shutil
import psutil
import argparse
import subprocess

from concurrent.futures import ThreadPoolExecutor

def gfa_to_fa(gfa_file, out_fa):
    with open(out_fa, "w") as out_f:
        with open(gfa_file) as gfa:
            for gfa_line in gfa:
                if gfa_line.startswith("S"):
                    parts = gfa_line.strip().split("\t")
                    out_f.write(f">{parts[1]}\n{parts[2]}\n")

def hifi_preprocess(CELL_LINE, PREFIX, hifi_fastq, thread, dep_folder):
    depth_window = 100 * 1000
    os.makedirs(PREFIX, exist_ok=True)

    if dep_folder is None:
        install_dependency(os.path.join(PREFIX, '99_dependency'), True)
        dep_folder = os.path.join(PREFIX, '99_dependency')

    THREAD = str(thread)
    hifiasm_folder = os.path.join(PREFIX, '00_hifiasm')

    
    # De novo assembly
    os.makedirs(hifiasm_folder, exist_ok=True)
    subprocess.run([
        'hifiasm', '-o', f'{hifiasm_folder}/{CELL_LINE}',
        '--telo-m', 'CCCTAA', '-t', THREAD] + hifi_fastq + ['--primary'
    ], check=True)
    
    # Mapping read to reference
    refseq = os.path.join(dep_folder, 'chm13v2.0.fa')

    depth_folder = os.path.join(PREFIX, '01_depth')
    os.makedirs(depth_folder, exist_ok=True)
    sam_file = os.path.join(depth_folder, f'{CELL_LINE}.sam')
    bam_file = os.path.join(depth_folder, f'{CELL_LINE}.bam')
    sorted_bam_file = os.path.join(depth_folder, f'{CELL_LINE}.sorted.bam')


    if not os.path.isfile(os.path.join(depth_folder, f'{CELL_LINE}.win.stat.gz')):
        subprocess.run([
            "minimap2", "-x", "map-hifi", "-K", "10G", "-t", THREAD,
            "-a", refseq] + hifi_fastq + ["-o", sam_file
        ], check=True)

        subprocess.run([
            "sambamba", "view", "-l", "5", "-t", THREAD, "-f", "bam", "-S",
            sam_file, "-o", bam_file
        ], check=True)
        os.remove(sam_file)

        subprocess.run([
            "sambamba", "sort", "-t", THREAD, bam_file, sorted_bam_file
        ], check=True)
        os.remove(bam_file)
        
        subprocess.run([
            os.path.join(dep_folder, 'PanDepth/bin/pandepth'), "-w", str(depth_window), "-t", THREAD,
            "-o", os.path.join(depth_folder, CELL_LINE),
            "-i", sorted_bam_file
        ], check=True)

    # Mapping to assembly
    out_fa_list = []
    os.makedirs(os.path.join(PREFIX, "10_asm"), exist_ok=True)
    for gfa_suffix in ['p_ctg.gfa', 'r_utg.gfa']:
        gfa_file = f"{hifiasm_folder}/{CELL_LINE}.{gfa_suffix}"
        kind = gfa_suffix[0]
        out_fa = os.path.join(PREFIX, f"10_asm/{CELL_LINE}.{kind}.fa")
        gfa_to_fa(gfa_file, out_fa)

        out_fa_list.append(out_fa)

    return out_fa_list

def install_dependency(dep_folder, force):
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

def run_alignasm(PREFIX_PATH, thread, fa_loc, ref_loc, ALIGNASM_LOC):
    THREAD = str(thread)

    subprocess.run([
        "minimap2", "--cs", "-t", THREAD, "-x", "asm20",
        "--no-long-join", "-r2k", "-K10G",
        ref_loc, fa_loc, "-o", f"{PREFIX_PATH}.paf"
    ], check=True)

    subprocess.run([
        "python3", "paf_gap_seq.py",
        fa_loc, f"{PREFIX_PATH}.paf", f"{PREFIX_PATH}.pat.fa"
    ], check=True)

    subprocess.run([
        "minimap2", "--cs", "-t", THREAD, "-x", "asm20",
        "-r2k", "-K10G",
        ref_loc, f"{PREFIX_PATH}.pat.fa", "-o", f"{PREFIX_PATH}.alt.paf"
    ], check=True)

    subprocess.run([
        ALIGNASM_LOC,
        f"{PREFIX_PATH}.paf", "-a", f"{PREFIX_PATH}.alt.paf",
        "-t", THREAD, "--non_skip_linkable"
    ], check=True)

    return f"{PREFIX_PATH}.paf", f"{PREFIX_PATH}.aln.paf"


def analysis(CELL_LINE, PREFIX, contig_loc, unitig_loc, depth_loc, thread, dep_folder, is_progress):
    os.makedirs(PREFIX, exist_ok=True)

    if dep_folder is None:
        install_dependency(os.path.join(PREFIX, '99_dependency'), True)
        dep_folder = os.path.join(PREFIX, '99_dependency')
    
    alignasm_ref_loc = os.path.join(dep_folder, 'chm13v2.0_maskedY_noM.fa')
    alignasm_loc = os.path.join(dep_folder, 'alignasm', 'build', 'alignasm')
    os.makedirs(os.path.join(PREFIX, '20_alignasm'), exist_ok=True)

    alignasm_folder_ctg = os.path.join(PREFIX, '20_alignasm', f'{CELL_LINE}.ctg')
    alignasm_folder_utg = os.path.join(PREFIX, '20_alignasm', f'{CELL_LINE}.utg')

    if thread > 1:
        alignasm_thread = thread // 2
        alignasm_worker_thread = 2
    else:
        alignasm_thread = 1
        alignasm_worker_thread = 1

    args_ctg = (alignasm_folder_ctg, alignasm_thread, contig_loc, alignasm_ref_loc, alignasm_loc)
    args_utg = (alignasm_folder_utg, alignasm_thread, unitig_loc, alignasm_ref_loc, alignasm_loc)

    with ThreadPoolExecutor(max_workers=alignasm_worker_thread) as executor:
        future_ctg = executor.submit(run_alignasm, *args_ctg)
        future_utg = executor.submit(run_alignasm, *args_utg)

        ctg_paf, ctg_aln_paf = future_ctg.result()
        utg_paf, utg_aln_paf = future_utg.result()

    depth_loc = os.path.abspath(depth_loc)
    run_skype(CELL_LINE, os.path.abspath(os.path.join(PREFIX, "30_skype")), ctg_paf, ctg_aln_paf, utg_paf, utg_aln_paf, depth_loc, thread, dep_folder, is_progress)

def run_skype(CELL_LINE, PREFIX, ctg_paf, ctg_aln_paf, utg_paf, utg_aln_paf, depth_loc, thread, dep_folder, is_progress):
    skype_folder_loc = os.path.join(dep_folder, "SKYPE")
    
    TEL_BED = os.path.join(skype_folder_loc, "public_data/chm13v2.0_telomere.bed")
    CHR_FAI = os.path.join(skype_folder_loc,"public_data/chm13v2.0.fa.fai")
    RPT_BED = os.path.join(skype_folder_loc,"public_data/chm13v2.0_repeat.m.bed")
    RCS_BED = os.path.join(skype_folder_loc,"public_data/chm13v2.0_censat_v2.1.m.bed")
    CYT_BED = os.path.join(skype_folder_loc,"public_data/chm13v2.0_cytobands_allchrs.bed")
    

    # python to bash variable
    MAIN_STAT_LOC = depth_loc
    PAF_LOC = ctg_aln_paf
    PAF_UTG_LOC = utg_aln_paf
    THREAD = str(thread)
    PROGRESS = ["--progress"] if is_progress else []

    subprocess.run([
        "python", os.path.join(skype_folder_loc, "00_Contig_Preprocessing.py"), 
        PAF_LOC, TEL_BED, CHR_FAI, RPT_BED, RCS_BED, MAIN_STAT_LOC, PREFIX,
        "--alt", PAF_UTG_LOC
    ], check=True)

    subprocess.run([
        "python", os.path.join(skype_folder_loc, "02_Build_Breakend_Graph_Limited.py"),
        f"{PAF_LOC}.ppc.paf", CHR_FAI, RCS_BED, PREFIX,
        "--orignal_paf_loc", ctg_paf, utg_paf,
        "-t", THREAD
    ] + PROGRESS, check=True)

    subprocess.run([
        "python", os.path.join(skype_folder_loc, "11_Ref_Outlier_Contig_Modify.py"),
        PAF_LOC, CHR_FAI, f"{PAF_LOC}.ppc.paf", PREFIX,
        "--alt", PAF_UTG_LOC
    ], check=True)

    thread_lim = int(psutil.virtual_memory().total / (1024**3) / 8)
    subprocess.run([
        "python", os.path.join(skype_folder_loc, "21_run_depth.py"),
        PAF_LOC, f"{PAF_LOC}.ppc.paf", PREFIX, "--alt", PAF_UTG_LOC,
        "--pandepth_loc", os.path.join(dep_folder, 'PanDepth', 'bin', 'pandepth'),
        "-t", str(max(min(thread_lim, thread), 1))
    ] + PROGRESS, check=True)

    subprocess.run([
        "python", os.path.join(skype_folder_loc, "22_save_matrix.py"),
        RCS_BED, f"{PAF_LOC}.ppc.paf", MAIN_STAT_LOC,
        TEL_BED, CHR_FAI, CYT_BED, PREFIX, "-t", THREAD
    ] + PROGRESS, check=True)

    subprocess.run([
        "python", "-X", f"juliacall-threads={THREAD}", "-X", "juliacall-handle-signals=yes",
        "23_run_nnls.py",
        os.path.abspath(PREFIX), "-t", THREAD
    ], check=True, cwd=skype_folder_loc)

    subprocess.run([
        "python", os.path.join(skype_folder_loc, "30_depth_analysis.py"),
        RCS_BED, f"{PAF_LOC}.ppc.paf", MAIN_STAT_LOC,
        TEL_BED, CHR_FAI, CYT_BED, PREFIX, "-t", THREAD
    ] + PROGRESS, check=True)

    subprocess.run([
        "python", os.path.join(skype_folder_loc, "31_virtual_sky.py"),
        f"{PAF_LOC}.ppc.paf", MAIN_STAT_LOC,
        TEL_BED, CHR_FAI, PREFIX, CELL_LINE
    ], check=True)


def main():
    parser = argparse.ArgumentParser(description="SKYPE pipeline")

    # Define parser
    subparsers = parser.add_subparsers(dest="command", required=True)
    
    parser_hifi_prepro = subparsers.add_parser("preprocess_hifi", help="Hifi preprocessing for SKYPE pipeline")

    parser_hifi_prepro.add_argument("WORK_DIR", type=str, help="Working directory for pipeline")
    
    parser_hifi_prepro.add_argument("HIFI_FASTQ", type=str, help="Location of HiFi fastq file", nargs="+")

    parser_hifi_prepro.add_argument("-t", "--thread", type=int, help="Number of thread", default=1)

    parser_hifi_prepro.add_argument("-p", "--prefix", type=str, help="Prefix for pipeline", default="SKYPE")

    parser_hifi_prepro.add_argument("--dependency_loc", type=str, help="SKYPE dependency folder location. If no value is specified, it will be installed automatically")

    
    parser_anl = subparsers.add_parser("analysis", help="Analysis cancer cell line karyotype using contig")

    parser_anl.add_argument("WORK_DIR", type=str, help="Working directory for pipeline")

    parser_anl.add_argument("CONTIG", type=str, help="Contig fasta location")

    parser_anl.add_argument("UNITIG", type=str, help="UNITIG fasta location")

    parser_anl.add_argument("DEPTH_LOC", type=str, help="Depth result for Pandepth location")

    parser_anl.add_argument("-t", "--thread", type=int, help="Number of thread", default=1)

    parser_anl.add_argument("-p", "--prefix", type=str, help="Prefix for pipeline", default="SKYPE")

    parser_anl.add_argument("--progress", help="Show progress bar", action='store_true')

    
    parser_dep = subparsers.add_parser("install_dependency", help="Hifi preprocessing for SKYPE pipeline")

    parser_dep.add_argument("dependency_loc", type=str, help="SKYPE dependency folder location")

    parser_dep.add_argument("--force", help="Reinstall dependency folder", action='store_true')

    
    parser_run = subparsers.add_parser("run_hifi", help="Pipeline for cancer hifi sequencing data")

    parser_run.add_argument("WORK_DIR", type=str, help="Working directory for pipeline")
    
    parser_run.add_argument("HIFI_FASTQ", type=str, help="Location of HiFi fastq file", nargs="+")

    parser_run.add_argument("-t", "--thread", type=int, help="Number of thread", default=1)

    parser_run.add_argument("-p", "--prefix", type=str, help="Prefix for pipeline", default="SKYPE")

    parser_run.add_argument("--progress", help="Show progress bar", action='store_true')

    parser_run.add_argument("--dependency_loc", type=str, help="SKYPE dependency folder location. If no value is specified, it will be installed automatically")
    

    args = parser.parse_args()

    if args.command == "hifi_preprocess":
        hifi_preprocess(args.prefix, args.WORK_DIR, args.HIFI_FASTQ, args.thread, args.dependency_loc)
    elif args.command == "install_dependency":
        install_dependency(args.dependency_loc, args.force)
    elif args.command == "analysis":
        analysis(args.prefix, args.WORK_DIR, args.CONTIG, args.UNITIG, args.DEPTH_LOC, args.thread, args.dependency_loc, args.progress)
    elif args.command == 'run_hifi':
        ctg_loc, utg_loc = hifi_preprocess(args.prefix, args.WORK_DIR, args.HIFI_FASTQ, args.thread, args.dependency_loc)
        
        depth_loc = os.path.join(args.WORK_DIR, '01_depth', f'{args.prefix}.win.stat.gz')
        if args.dependency_loc is None:
            dep_folder = os.path.join(args.WORK_DIR, '99_dependency')
        else:
            dep_folder = args.dependency_loc
        analysis(args.prefix, args.WORK_DIR, ctg_loc, utg_loc, depth_loc, args.thread, dep_folder, args.progress)

if __name__ == "__main__":
    main()
