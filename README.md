# ACCtools pipeline

Complete pipeline for analyzing a cancer genome using ACCtools.

## Dependency
Please make anaconda environment to use SKYPE pipeline

```bash
mamba create -n skype -c conda-forge -c bioconda \
python=3.12 gxx cmake=3 zip psutil aria2 pyfaidx hifiasm flye minimap2 samtools \
"numpy<2" scipy matplotlib tqdm pycirclize=1.9 pandas networkx graph-tool=2.98 \
seaborn h5py vcfpy scikit-learn=1.6 

mamba activate skype
pip install juliacall adelie
```

## Example
```bash
git clone https://github.com/ACCtools/ACCtools-pipeline
cd ACCtools-pipeline

mamba activate skype

# Pacbio HiFi
python SKYPE.py run_hifi <Working directory> <hifi.fastq(.gz) ...>

# ONT R10 (HQ)
python SKYPE.py run_hifi --hifiasm_args="--ont --chem-c 0" <Working directory> <ontr10.fastq(.gz) ...>

# PacBio CLR
python SKYPE.py run_flye <Working directory> pacbio-raw <clr.fastq(.gz) ...>

# ONT R9
python SKYPE.py run_flye <Working directory> nano-raw <ontr9.fastq(.gz) ...>
```

## Native stage 01/10 pipeline

The native pipeline runs `01_Preprocess_NClose.py` followed by
`10_Graph_Find_Paths.py`, then stages 11, 21, 22, 23, and 31.
Full-assembly mode keeps its separate entry point.

`--option_skype` (`--option-skype` also works) forwards a quoted option string
to stage 01. SKYPE owns option routing in `skype_options.py`: stage 01 applies
its preprocessing options and writes `skype_options.json` in the output
directory. Stage 10 reads its options from that file. Each fresh stage-01 run
replaces the settings, including when no extra options are supplied.

```bash
python SKYPE.py analysis \
  --option_skype="--add_indel_graph" \
  <Working directory> <contig.fa> <unitig.fa> <depth.win.stat.gz>
```

The following options are recognized in the native route:

| Argument | Default | Role and notes |
| --- | --- | --- |
| `--exclude_nclose_list_loc <PATH>` | None | Stage 01 user exclusion list. |
| `--check_nclose_count` / `--nclose_count_vaf_threshold <FLOAT>` | Disabled / `0.1` | Stage 01 raw-read junction VAF filter. |
| `--add_indel_graph` | Disabled | Adds selected depth-supported type-4 indel rescue edges to the graph without increasing its dimensions. In VCF mode, DEL, DUP, and indel-like BND events are eligible; INS events are excluded. |
| `--skip_bam_analysis` | Disabled | Skips raw-read BAM validation of translocation candidates and the removal of raw-read-supported virtual-inversion candidates. Do not combine it with `--check_nclose_count`, because the requested VAF filter cannot run when BAM analysis is skipped. |
| `--vcf_filter_pass <FILTER> [FILTER ...]` | `PASS .` | In VCF input mode, replaces the accepted `FILTER` values. Matching is exact and case-sensitive; outside VCF mode this argument is ignored. |
| `--verbose`, `--limit_combinations <PATH>` | Disabled / automatic | Stage 10 graph-search diagnostics or an exact limit pair. |

Inputs and execution controls are constructed by ACCtools and are rejected in
`--option_skype` in native mode:

| Argument | Use in ACCtools instead |
| --- | --- |
| `-t`, `--thread` | Use the top-level `-t` / `--thread` option. |
| `-d`, `--graph_depth` | Use the top-level `-d` / `--graph_depth` option. |
| `--progress` | Use the top-level `--progress` option. |
| `--vcf_input` | Use `--benchmark_vcf_loc`; ACCtools also prepares insertion-sequence alignments and selects VCF mode correctly. |
| `--alt`, `--original_paf_loc` | Do not set these manually. ACCtools derives them from the contig/unitig alignments or the VCF insertion-sequence alignment. |

For `run_hifi`, place `--option_skype` before `<Working directory>` because every
argument after the working directory is interpreted as an input read file.

Native restart stages are `0, 1, 10, 11, 21, 22, 23, 31`. A nonzero restart
validates the artifacts required from the skipped stages before launching
subprocesses. A stage-10 restart reuses saved options; an explicit nonempty
`--option_skype` replaces its graph options for that invocation. Preprocessing
options require restarting at stage 01. Standalone stage 10 also accepts
`--option_skype=""` to use default graph options, and its direct CLI options
take precedence over saved settings. Missing settings files use defaults.

From the workspace wrapper:

```bash
bash run.sh --option_skype="--skip_bam_analysis --add_indel_graph" HCC1937
```

`--simple` adds `--add_indel_graph` to the same option string.

## Analysis inputs

The default native workflow discovers rearrangements from assembly alignments,
constructs one observed-depth matrix, and fits every candidate column with one
raw NNLS solve. It produces `SV_call_result.vcf`, `SKYPE_result.bed`,
`nclose_report.tsv`, and `total_cov.*`. There are no karyotype/variant mode
flags, normal-chromosome prior, post-NNLS filtering, or `_filter`/`_cluster`
result sets.

### VCF input mode

VCF input mode uses structural variants from an existing VCF instead of discovering NClose junctions from assembly alignments. The assembly alignment is still used for telomere/neotelomere anchors, and the mapped reads are still used to estimate depth and copy-number support.

```bash
# VCF and reads must use the same reference build (hs1 is the default)
python SKYPE.py run_hifi \
  --benchmark_vcf_loc <input.vcf> \
  --reference hs1 \
  <Working directory> <hifi.fastq(.gz) ...>

# Use an existing assembly and PanDepth result
python SKYPE.py analysis \
  --benchmark_vcf_loc <input.vcf> \
  --reference hg38 \
  <Working directory> <contig.fa> <unitig.fa> <depth.win.stat.gz>
```

By default, only records whose `FILTER` value is exactly `PASS` or `.` are evaluated. Replace that set through the `--option_skype` option when necessary:

```bash
python SKYPE.py run_hifi \
  --benchmark_vcf_loc <input.vcf> \
  --option_skype="--vcf_filter_pass PASS . Candidate" \
  <Working directory> <hifi.fastq(.gz) ...>
```

The main VCF-mode result is `<SKYPE output directory>/SV_benchmark_result.vcf`. It preserves the input records and adds `SKYPE_CN` and `SKYPE_STATUS`; records with side-specific measurements also receive `SKYPE_CN_DETAIL` and `SKYPE_STATUS_DETAIL`. The default SKYPE output directory is `<Working directory>/30_skype` for `hs1` and `<Working directory>/31_skype_hg38` for `hg38`. Parsing diagnostics are written to `vcf_mode_summary.json`, `vcf_mode_summary.tsv`, `vcf_mode_skipped_records.tsv`, and `vcf_mode_orientation_mismatches.tsv` in the same directory.

### Full-assembly input mode

Full-assembly mode uses each record of a complete genome assembly FASTA as one
matrix path. It aligns the FASTA to the selected reference with minimap2 and
alignasm, then passes the resulting `*.aln.paf` to the standalone
`full_assembly_pipeline.py`. None of the normal numbered stage scripts are
entered; it retains its Virtual SKY/karyotype output through the reusable
plotting modules and is mutually exclusive with
`--benchmark_vcf_loc`.

```bash
python SKYPE.py analysis \
  --full_assembly /path/HG008T_v3.2.fasta \
  --reference hs1 \
  <Working directory> <ignored-contig.fa> <ignored-unitig.fa> <depth.win.stat.gz>
```

The option is also accepted by `run_hifi` and `run_flye`; those commands still
map the reads and calculate sample depth but skip hifiasm/Flye assembly. Place
the option before `<Working directory>` for `run_hifi`.

## Compatible VCF inputs

SKYPE has caller-aware handling for the following structural-variant VCFs:

| Caller | Compatibility notes |
| --- | --- |
| Sniffles2 | Supports standard paired or singleton BND records with an explicit remote locus. |
| Severus | The version must be identifiable from `##source`. For versions before 1.7, BND records must contain a valid `STRANDS` value so SKYPE can correct mixed-strand orientations. |
| nanomonsv | Supports standard BND records and uses `SVINSLEN` as a fallback insertion length. |
| SAVANA | Supports standard BND and symbolic SV records. |
| GRIPSS | Can be detected from `##gripssVersion`; paired BND records are supported, while single breakends without a remote locus are reported and skipped. |
| Manta | Uses `MATEID` to keep microhomology-shifted mate records as one junction. |
| SvABA | Supports standard paired BND records. |

VCFs from other callers are parsed using standard BND ALT notation and symbolic SV fields when possible, but compatibility is not guaranteed.

| SV type | Required representation |
| --- | --- |
| `BND` | A canonical bracketed ALT allele containing the remote chromosome and position. `MATEID` or `MATE_ID` is recommended; otherwise SKYPE attempts reciprocal-coordinate pairing. Paired ALT orientations must be reciprocal. |
| `INV` | `SVTYPE=INV` and a valid `END` value. |
| `DEL` / `DUP` | `SVTYPE`, a valid `END`, and a reference span or absolute `SVLEN` of at least 100 kbp. |
| `INS` | `SVTYPE=INS`, an absolute `SVLEN` of at least 100 kbp, and an insertion sequence in the sequence-resolved ALT allele or `INSSEQ`, `SVINSSEQ`, or `SEQ`. nanomonsv `SVINSLEN` is used when `SVLEN` is absent. The ALT alignment supplies the reference path but does not determine whether the event passes the 100 kbp threshold. |

The input must also satisfy these requirements:

- Use a plain-text, uncompressed VCF; `.vcf.gz` is not accepted.
- Include exactly one `#CHROM` header and at least the eight standard VCF columns. Sample columns may follow them.
- Use contig names and coordinates from the selected `--reference` build (`hs1` or `hg38`). A length mismatch for a contig shared by the VCF and reference is fatal; events on contigs absent from the selected primary reference are skipped.
- Only the first ALT allele of a multiallelic record is evaluated.
- Unsupported `SVTYPE` values and malformed records are reported and skipped rather than evaluated as graph events; input rows are preserved in `SV_benchmark_result.vcf`.
