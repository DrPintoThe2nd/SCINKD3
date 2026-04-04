#!/bin/bash

#SBATCH -c 24                           # number of cores
#SBATCH --mem=64GB                      # specify additional memory
#SBATCH -t 0-02:00:00               # wall time (D-HH:MM)
#SBATCH -p htc                  # Uses the
#SBATCH -q public                   # Submits job
#SBATCH -o test.%j.out                    # STDOUT (%j = JobId)
#SBATCH -e test.%j.err                    # STDERR (%j = JobId)
#SBATCH --mail-type=ALL             # Send a notification when a job starts, stops, or fails
#SBATCH --mail-user=<email@address> # send-to address
#SBATCH --export=NONE

module purge

source activate scinkd3

snakemake --use-conda --rerun-incomplete --nolock --cores 24 -s SCINKD3/SCINKD.v3.1.1.py
