#!/bin/bash


#SBATCH -J canon_KWAS 
#SBATCH -N 1                      # Ensure that all cores are on one machine
#SBATCH -n 16                # Use n cores for one job 
#SBATCH -t 0-11:59                # Runtime in D-HH:MM 
#SBATCH --mem=64G            # Memory pool for all cores 
#SBATCH -o CanonKWAS.%A.out       # File to which STDOUT will be written 
#SBATCH -e CanonKWAS.%A.err       # File to which STDERR will be written 
#SBATCH --account=def-coling

~/programs/kmerGWAS/external_programs/kmc_v3 -t16 -k31 -ci2 @input_files.txt output_kmc_canon ./ 1> kmc_canon.1 2> kmc_canon.2
~/programs/kmerGWAS/external_programs/kmc_v3 -t16 -k31 -ci0 -b @input_files.txt output_kmc_all ./ 1> kmc_all.1 2> kmc_all.2
