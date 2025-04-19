# ACCtools pipeline

Complete pipeline for analyzing a cancer genome using ACCtools.

## Dependency
Please make anaconda environment to use SKYPE pipeline

```bash
mamba create -n skype -c conda-forge -c bioconda \
gxx cmake=3 zip psutil aria2 pyfaidx hifiasm flye minimap2 sambamba \
numpy scipy matplotlib tqdm pycirclize pandas networkx graph-tool seaborn h5py

mamba activate skype
pip install juliacall
```

## Example
```bash
git clone https://github.com/ACCtools/ACCtools-pipeline
cd ACCtools-pipeline

mamba activate skype
python SKYPE.py run_hifi <Working directory> <hifi.fastq(.gz) ...>
python SKYPE.py run_flye <Working directory> nano-raw <nanoR9.fastq(.gz) ...> # run as flye --nano-raw ...
```