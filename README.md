# ACCtools pipeline

Complete pipeline for analyzing a de novo cancer genome with ACCtools.  
*ACCtools is currently under development, so it is given in bash, but will be developed into a Python CLI in the future.*

## Dependency
- minimap2
- pyfaidx
- alignasm

Use anaconda to install package `minimap2`, `pyfaidx`.  
See [README.md](https://github.com/ACCtools/alignasm) to install `alignasm`.

## Pipeline
```bash
PREFIX=<Cancer ID>

CANCER_GENOME_LOC=<cancer genome location>
REF_LOC=<reference location>

THREAD=<number of threads>
ALIGNASM_LOC=<alignasm executable location>

mkdir -p $(dirname $PREFIX)

ln -s $(realpath $CANCER_GENOME_LOC) $PREFIX.fa
minimap2 --cs -t $THREAD -x asm20 --no-long-join -r2k -K10G $REF_LOC $PREFIX.fa -o $PREFIX.paf
python3 paf_gap_seq.py $PREFIX.fa $PREFIX.paf $PREFIX.pat.fa
minimap2 --cs -t $THREAD -x asm20 -r2k -K10G -P $REF_LOC $PREFIX.pat.fa -o $PREFIX.alt.paf

$ALIGNASM_LOC $PREFIX.paf -a $PREFIX.alt.paf
sort -k6,6 -k8,8n $PREFIX.aln.paf | paftools.js call -l 0 -L 0 -q 0 -f $REF_LOC - > $PREFIX.vcf
```