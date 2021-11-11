#!/bin/bash


#SBATCH -J tableKWAS
#SBATCH -N 1                      # Ensure that all cores are on one machine
#SBATCH -n 1                # Use n cores for one job
#SBATCH -t 0-03:00                # Runtime in D-HH:MM
#SBATCH --mem=10G            # Memory pool for all cores
#SBATCH -o tableKWAS.%A.out       # File to which STDOUT will be written
#SBATCH -e tableKWAS.%A.err       # File to which STDERR will be written
#SBATCH --account=def-coling

~/programs/kmerGWAS/bin/build_kmers_table -l final_kmers_clean_list_path.txt -k 31 -a kmers_to_use -o kmers_table
