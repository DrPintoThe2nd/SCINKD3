#!/bin/bash

#SBATCH -c 48                  		# number of cores
#SBATCH --mem=200GB              	# specify additional memory
#SBATCH -t 0-04:00:00          	    # wall time (D-HH:MM)
#SBATCH -p htc	      	        # Uses the
#SBATCH -q public     	            # Submits job
#SBATCH -o test.%j.out       	 	# STDOUT (%j = JobId)
#SBATCH -e test.%j.err       	 	# STDERR (%j = JobId)
#SBATCH --mail-type=ALL             # Send a notification when a job starts, stops, or fails
#SBATCH --mail-user=bjp004@morningside.edu # send-to address
#SBATCH --export=NONE

module purge

source activate scinkd

snakemake --use-conda --rerun-incomplete --nolock --cores 48 -j 1 -s SCINKD.v3.1.0_1-sample.py
