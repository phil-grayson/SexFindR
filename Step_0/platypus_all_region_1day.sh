#!/bin/bash

#SBATCH -J platy
#SBATCH -N 1                      # Ensure that all cores are on one machine
#SBATCH -n 24                # Use n cores for one job 
#SBATCH -t 1-00:00                # Runtime in D-HH:MM 
#SBATCH --mem-per-cpu=3900 #--mem=10000            # Memory pool for all cores 
#SBATCH -o Platy.%A.out       # File to which STDOUT will be written 
#SBATCH -e Platy.%A.err       # File to which STDERR will be written 
#SBATCH --account=def-docker

module load samtools
module load bcftools
source ~/cython/bin/activate

python /home/pgrayson/programs/Platypus/bin/Platypus.py callVariants --verbosity=3 --bufferSize=50000 --nCPU=24 --bamFiles=${1} --refFile=${2} --output=fugu_15M_14F.vcf
