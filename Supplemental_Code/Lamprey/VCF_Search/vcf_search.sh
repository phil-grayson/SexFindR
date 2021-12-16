#!/bin/bash

#SBATCH -J vcfSearch   # Name for the job (keep it short and informative)
#SBATCH -N 1       # Number of nodes
#SBATCH -n 1       # Use n cores
#SBATCH -t 0-03:00     # Runtime in D-HH:MM 
#SBATCH --mem=1000 # Memory requested (mb default, or specify G for Gb) 
#SBATCH -o vcfSearch.%A.out       # File to which STDOUT will be written 
#SBATCH -e vcfSearch.%A.err       # File to which STDERR will be written 
#SBATCH --account=def-coling # Who are are going to charge it to?

python3 lamprey_vcf_search_population.py filtered_PASS_sea_lamprey_304_update.vcf filtered_PASS_sea_lamprey_304_update_candidates.txt $1
