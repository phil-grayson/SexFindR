#!/bin/bash


#SBATCH -J kmerGWAS 
#SBATCH -N 1                      # Ensure that all cores are on one machine
#SBATCH -n 1                # Use n cores for one job 
#SBATCH -t 0-03:00                # Runtime in D-HH:MM 
#SBATCH --mem=10G            # Memory pool for all cores 
#SBATCH -o combineKWAS.%A.out       # File to which STDOUT will be written 
#SBATCH -e combineKWAS.%A.err       # File to which STDERR will be written 
#SBATCH --account=def-coling

~/programs/kmerGWAS/bin/list_kmers_found_in_multiple_samples -l final_kmers_clean_list_path.txt -k 31 --mac 5 -p 0.2 -o kmers_to_use
