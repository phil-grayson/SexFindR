#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 16:52:49 2020

@author: phil
"""

# execute with infile as $1, outfile as $2, number of permutations as $3

import csv
import sys

# get male names to index from the header of vcf
males=[]
with open("/Users/phil/Desktop/fugu/males.txt") as handle:
    reader=csv.reader(handle,delimiter='\t')
    for strLine in reader:
        if len(strLine)==1:
            males.append(strLine[0]+"_1.fastq")

# get females the same
females=[]
with open("/Users/phil/Desktop/fugu/females.txt") as handle:
    reader=csv.reader(handle,delimiter='\t')
    for strLine in reader:
        if len(strLine)==1:
            females.append(strLine[0]+"_1.fastq")

fileout=(sys.argv[2]) #"/Users/phil/Desktop/test_out.txt"#

fout= open(fileout, 'w') 
fout.write("chrom"+"\t"+"base"+"\t"+"male_geno"+"\t"+"female_geno""\n")

with open(sys.argv[1]) as handle:
#with open("/Users/phil/Desktop/test.txt") as handle:
    reader=csv.reader(handle,delimiter='\t')
    #print(next(reader))    
    for strLine in reader:
        if len(strLine) > 1:
            if strLine[0] == "#CHROM":
                #print(strLine)
                male_indexes = []
                female_indexes = []
                for male in males:
                    male_indexes.append(strLine.index(male))
                for female in females:
                    female_indexes.append(strLine.index(female))
                #print(male_indexes)
                #print(female_indexes)
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
