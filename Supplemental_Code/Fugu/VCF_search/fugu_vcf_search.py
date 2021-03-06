#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 16:52:49 2020

@author: phil
"""

# execute in directory with males/females file and with vcf (e.g., filtered_PASS_fugu_14M_13F.vcf) as $1, and outfile as $2

import csv
import sys

# get male names to index from the header of vcf
males=[]
with open("males.txt") as handle:
    reader=csv.reader(handle,delimiter='\t')
    for strLine in reader:
        if len(strLine)==1:
            males.append(strLine[0]+"_1.fastq")

# get females the same
females=[]
with open("females.txt") as handle:
    reader=csv.reader(handle,delimiter='\t')
    for strLine in reader:
        if len(strLine)==1:
            females.append(strLine[0]+"_1.fastq")

fileout=(sys.argv[2])

fout= open(fileout, 'w')
fout.write("chrom"+"\t"+"base"+"\t"+"male_geno"+"\t"+"female_geno""\n")

with open(sys.argv[1]) as handle:
    reader=csv.reader(handle,delimiter='\t')
    for strLine in reader:
        if len(strLine) > 1:
            if strLine[0] == "#CHROM":
                male_indexes = []
                female_indexes = []
                for male in males:
                    male_indexes.append(strLine.index(male))
                for female in females:
                    female_indexes.append(strLine.index(female))
            else:
                male_genotype = []
                for index in male_indexes:
                    male_genotype.append(strLine[index].split(":")[0])
                if len(set(male_genotype)) == 1:
                    female_genotype = []
                    for index in female_indexes:
                        female_genotype.append(strLine[index].split(":")[0])
                    if len(set(female_genotype)) == 1:
                        if set(male_genotype) != set(female_genotype):
                            fout.write(strLine[0]+"\t"+strLine[1]+"\t"+male_genotype[0]+"\t"+female_genotype[0]+"\n")
fout.close()
