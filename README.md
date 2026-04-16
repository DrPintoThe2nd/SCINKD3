# SCINKD3
Initial development repo of SCINKD[v3]

SCINKD[v2] is a tool to identify sex chromosomes from a single individual with a haplotype-resolved reference genome (https://github.com/DrPintoThe2nd/SCINKD) with a number a caveats stemming from, although generally robust, a small amount of input data, but also capable of running on a moderately powered desktop computer.

Development work is ongoing and eventually a manuscript forthcoming, until then if you find this workflow useful please cite the principle source: Pinto BJ, Gable SM, Keating SE, Smith CH, Gamble T, Nielsen SV, Wilson MA. (2026). Sex chromosome identification and genome curation from a single individual with SCINKD. _Molecular Biology and Evolution_. _In press_.

Increasing statistical power by including more individuals is essential in many cases where SCINKD[v2] has failed, e.g. a particular sex-limited region is too small, the homogametic sex was sequenced, or for validating initial SCINKD results. SCINKD3 is designed with this in mind. However, it is not possible to run an analysis at this scale without large amounts of RAM. This pipeline was designed for WGS illumina data, but can theoretically be run on an data type (WGS, RADseq, WES, RNAseq) sequenced to reasonable coverage (probably >=10x for WGS), but has not been tested extensively on varying levels of coverage or data types (only WGS and RADseq).

This is the first push development repo (alpha) and is very much use at your own risk. However, please feel free to open issues with questions or suggestions. Most of the plotting scripts are somewhat bare at present, but will be updated in time.

Assuming all read files are availabe in the current working directory we can modify the config file (1) to reflect their assignments in the config file along with (2) the read suffixes, (3) per-job CPU/RAM limits, and (4) the reference genome. Then, we can run SCINKD3 for multiple individuals with 24 cores and 64Gb of RAM:

```
git clone https://github.com/DrPintoThe2nd/SCINKD3.git
mamba env create -f SCINKD3/SCINKD3_environment.yml
mamba activate scinkd3
snakemake --use-conda --rerun-incomplete --nolock --cores 24 -s SCINKD3/SCINKD.v3.1.5.snakemake -np #dry-run
snakemake --use-conda --rerun-incomplete --nolock --cores 24 -s SCINKD3/SCINKD.v3.1.5.snakemake
```

At present, the final output of this workflow will produce sex-specific reads and a manhattan plot of sex-specific coverage for each individual (in repaired_reads/ and [F/M]_cov/, respectively); additionally it will output three summary files (1) a summary manhattan plot for females, (2) a summary manhattan plot for males, and (3) a pseudostatistical dotplot comparison between the sexes akin to the primary output file from SCINKD[v2] with outlier sequence labels applied. Example plots for a gecko, Sphaerodactylus townsendi (6F/7M RADseq dataset), a species with a known XX/XY system on chromosome 3 (Pinto et al. 2022; https://doi.org/10.1093/jhered/esac016):

```
(1) Females; Fspec_total.regions.manhattan.png
```

<img width="1600" height="600" alt="image" src="https://github.com/user-attachments/assets/c48b1f53-4620-403b-b85f-3e60324cd010" />

```
(2) Males; Mspec_total.regions.manhattan.png
```

<img width="1600" height="600" alt="image" src="https://github.com/user-attachments/assets/e0bba67c-6c62-47f3-894e-7def58787da5" />

```
(3) SCINKD dotplot; {reference genome name}.dotplot.png
```

<img width="1800" height="750" alt="image" src="https://github.com/user-attachments/assets/0e95080a-aabd-48ac-820d-e0247624a49b" />


