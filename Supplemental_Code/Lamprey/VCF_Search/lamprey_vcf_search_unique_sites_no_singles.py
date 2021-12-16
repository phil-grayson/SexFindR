#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 16:52:49 2020

@author: phil
"""

## execute with filtered_PASS_sea_lamprey_304_update.vcf as $1,

import csv
import sys

# get male names to index from the header of vcf
males=[]
females=[]
with open("sexed_update_for_python.txt") as handle:
    reader=csv.reader(handle,delimiter='\t')
    for strLine in reader:
        if strLine[1] == "M":
            males.append(strLine[0])
        if strLine[1] == "F":
            females.append(strLine[0])

fileout=(sys.argv[2])

fout= open(fileout, 'w')
fout.write("chrom"+"\t"+"base"+"\t"+"male_geno"+"\t"+"female_geno"+"\t"+"male_number"+"\t"+"female_number"+"\n")

with open(sys.argv[1]) as handle:
    reader=csv.reader(handle,delimiter='\t')
    #print(next(reader))
    for strLine in reader:
        if len(strLine) > 1:
            if strLine[0] == "#CHROM":
                male_indexes=[]
                for male in males:
                    x=-1
                    for meta in strLine:
                        x+=1
                        if male+"_R1" in meta:
                            male_indexes.append(x)
                female_indexes=[]
                for female in females:
                    x=-1
                    for meta in strLine:
                        x+=1
                        if female+"_R1" in meta:
                            female_indexes.append(x)
            else:
                male_genotype = []
                female_genotype = []
                for index in male_indexes:
                    male_genotype.append(strLine[index].split(":")[0])
                for index in female_indexes:
                    female_genotype.append(strLine[index].split(":")[0])
                if set(male_genotype) != set(female_genotype):
                    male_dict={}
                    female_dict={}
                    for genotype in set(male_genotype):
                        male_dict[genotype]=male_genotype.count(genotype)
                    for genotype in set(female_genotype):
                        female_dict[genotype]=female_genotype.count(genotype)
                    male_singles=[]
                    female_singles=[]
                    for key in male_dict.keys():
                        if male_dict[key] == 1:
                            male_singles.append(key)
                    for key in male_singles:
                        del male_dict[key]
                    for key in female_dict.keys():
                        if female_dict[key] == 1:
                            female_singles.append(key)
                    for key in female_singles:
                        del female_dict[key]
                    if male_dict.keys() != female_dict.keys():
                        fout.write(strLine[0]+"\t"+strLine[1]+"\t"+str(male_dict)+"\t"+str(female_dict)+"\t"+str(len(male_dict))+"\t"+str(len(female_dict))+"\n")
fout.close()
