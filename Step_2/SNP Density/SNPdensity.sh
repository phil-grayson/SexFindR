#!/bin/bash

#SBATCH -J bcf   # Name for the job (keep it short and informative)
#SBATCH -N 1       # Number of nodes
#SBATCH -n 1       # Use n cores
#SBATCH -t 0-03:00     # Runtime in D-HH:MM 
#SBATCH --mem=3900 # Memory requested (mb default, or specify G for Gb) 
#SBATCH -o bcf.%A.out       # File to which STDOUT will be written 
#SBATCH -e bcf.%A.err       # File to which STDERR will be written 
#SBATCH --account=def-coling # Who are are going to charge it to?

module load bcftools/1.9

# $1 is name of file containing sample list
# $2 is output vcf name
# $3 is input vcf name

bcftools view -a -s $1 -o $2 $3

bgzip -c $2 > ${2}.gz

bcftools index ${2}.gz

module load vcftools/0.1.14

vcftools --vcf $2 --SNPdensity 10000 --out ${2}
