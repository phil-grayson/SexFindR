#!/bin/bash


#SBATCH -J strand_KWAS 
#SBATCH -N 1                      # Ensure that all cores are on one machine
#SBATCH -n 1                # Use n cores for one job 
#SBATCH -t 0-11:59                # Runtime in D-HH:MM 
#SBATCH --mem=100G            # Memory pool for all cores 
#SBATCH -o strandKWAS.%A.out       # File to which STDOUT will be written 
#SBATCH -e strandKWAS.%A.err       # File to which STDERR will be written 
#SBATCH --account=def-coling

~/programs/kmerGWAS/bin/kmers_add_strand_information -c output_kmc_canon -n output_kmc_all -k 31 -o kmers_with_strand
#echo "done, beginning cleanup"
#rm *.kmc*
#echo "done cleanup"
