# SCINKD3
development repo of SCINKD (v3)


SCINKD is a tool to identify sex chromosomes from a single individual with a haplotype-resolved reference genome (https://github.com/DrPintoThe2nd/SCINKD) with a number a caveats stemming from, although generally robust, a small amount of input data capable of running on a moderately powered desktop computer.

Increasing statistical power by including more individuals is essential in many cases where SCINKD has failed, e.g. a particular sex-limited region is too small or the homogametic sex was sequenced, or for validating initial SCINKD results. SCINKD3 is designed with this in mind. However, it is not possible to run an analysis at this scale on anything less than an HPC with significant amounts of RAM.

This is the first push development repo (pre-alpha) and is very much use at your own risk.

Assuming all read files are availabe in the current working directory one can run SCINKD3 for multiple individuals with 24 cores and 64Gb of RAM:
```
git clone https://github.com/DrPintoThe2nd/SCINKD3.git
mamba env create -f SCINKD3_environment.yml
mamba activate scinkd
snakemake --use-conda --rerun-incomplete --nolock --cores 24 -s SCINKD3/SCINKD.v3.1.0.py -np #dry-run
snakemake --use-conda --rerun-incomplete --nolock --cores 24 -s SCINKD3/SCINKD.v3.1.0.py 
```
