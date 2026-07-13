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

## Forwarding SKYPE stage-02 options

`--option_02` forwards a quoted string of additional arguments to
`02_Build_Breakend_Graph_Limited.py`, the SKYPE breakend-graph construction
stage. Multiple arguments must be placed in the same quoted string:

```bash
python SKYPE.py analysis \
  --option_02="--add_indel_graph" \
  <Working directory> <contig.fa> <unitig.fa> <depth.win.stat.gz>
```

The following public stage-02 arguments can be forwarded through
`--option_02`:

| Argument | Default | Role and notes |
| --- | --- | --- |
| `--disable_alt_ctg_simple` | Not set (rescue enabled) | Disables the default rescue of rearrangement candidates from primary contigs that were not otherwise retained. The rescue trims telomere-like terminal chunks, ignores fragments up to 10 kbp, selects chromosomes covering 90% of the remaining span, and removes nearby duplicate same-direction chromosome changes. |
| `--add_indel_graph` | Disabled | Adds selected depth-supported type-4 indel rescue edges to the graph without increasing its dimensions. In VCF mode, DEL, DUP, and indel-like BND events are eligible; INS events are excluded. |
| `--skip_bam_analysis` | Disabled | Skips raw-read BAM validation of translocation candidates and the removal of raw-read-supported virtual-inversion candidates. Do not combine it with `--check_nclose_count`, because the requested VAF filter cannot run when BAM analysis is skipped. |
| `--karyotype_mode` | Enabled | Selects the default karyotype-oriented mode with aggressive filtering. It is mutually exclusive with `--variant_mode`. |
| `--variant_mode` | Disabled | Selects depth-preserving variant analysis and disables karyotype-specific filtering and clustering. VCF input mode enables it automatically. |
| `--vcf_filter_pass <FILTER> [FILTER ...]` | `PASS .` | In VCF input mode, replaces the accepted `FILTER` values. Matching is exact and case-sensitive; outside VCF mode this argument is ignored. |

Some stage-02 arguments are already constructed by ACCtools and should not be
overridden through `--option_02`:

| Stage-02 argument | Use in ACCtools instead |
| --- | --- |
| `-t`, `--thread` | Use the top-level `-t` / `--thread` option. |
| `-d`, `--graph_depth` | Use the top-level `-d` / `--graph_depth` option. |
| `--progress` | Use the top-level `--progress` option unless progress is needed only for stage 02. |
| `--vcf_input` | Use `--benchmark_vcf_loc`; ACCtools also prepares insertion-sequence alignments and selects VCF mode correctly. |
| `--alt`, `--original_paf_loc` | Do not set these manually. ACCtools derives them from the contig/unitig alignments or the VCF insertion-sequence alignment. |

For `run_hifi`, place `--option_02` before `<Working directory>` because every
argument after the working directory is interpreted as an input read file.

## Analysis modes

The default is **karyotype mode**, which applies karyotype-oriented filtering and produces the base, `_filter`, and `_cluster` result sets. The two modes below are intended for retaining and quantifying structural-variant candidates. Because `run_hifi` treats every argument after the working directory as a read file, place all options before `<Working directory>`.

### Variant mode

Variant mode discovers rearrangements from the assembly alignments, as in the default workflow, but preserves depth-supported variant candidates by disabling the normal-chromosome prior and karyotype-specific filtering and clustering.

```bash
# Assembly-based variant analysis with PacBio HiFi reads
python SKYPE.py run_hifi \
  --option_02="--variant_mode" \
  <Working directory> <hifi.fastq(.gz) ...>
```

Variant mode produces only the unsuffixed result set, including `virtual_sky.*`, `karyotype.txt`, `total_cov.*`, `SV_call_result.vcf`, and `SKYPE_result.bed`. It is also selected automatically for VCF input mode and when the graph contains more than 1,000 NClose nodes.

### VCF input mode

VCF input mode uses structural variants from an existing VCF instead of discovering NClose junctions from assembly alignments. The assembly alignment is still used for telomere/neotelomere anchors, and the mapped reads are still used to estimate depth and copy-number support. This mode automatically enables variant mode.

The same mode flags can be used with `run_hifi`, `run_flye`, or `analysis`.

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

By default, only records whose `FILTER` value is exactly `PASS` or `.` are evaluated. Replace that set through the stage-02 forwarding option when necessary:

```bash
python SKYPE.py run_hifi \
  --benchmark_vcf_loc <input.vcf> \
  --option_02="--vcf_filter_pass PASS . Candidate" \
  <Working directory> <hifi.fastq(.gz) ...>
```

The main VCF-mode result is `<SKYPE output directory>/SV_benchmark_result.vcf`. It preserves the input records and adds `SKYPE_CN` and `SKYPE_STATUS`; records with side-specific measurements also receive `SKYPE_CN_DETAIL` and `SKYPE_STATUS_DETAIL`. The default SKYPE output directory is `<Working directory>/30_skype` for `hs1` and `<Working directory>/31_skype_hg38` for `hg38`. Parsing diagnostics are written to `vcf_mode_summary.json`, `vcf_mode_summary.tsv`, `vcf_mode_skipped_records.tsv`, and `vcf_mode_orientation_mismatches.tsv` in the same directory.

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
