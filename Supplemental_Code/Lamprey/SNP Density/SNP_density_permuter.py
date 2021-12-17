#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 16:52:49 2020

@author: phil
"""

# execute with infile as $1, outfile as $2, number of permutations as $3

import csv
import sys
import random

fileout=(sys.argv[2])
#fileout="/Users/phil/Desktop/test_permuttttter.txt"

fout= open(fileout, 'w') 

perms = sys.argv[3]


#header = []
##with open("/Users/phil/Desktop/SNPdensity_individual_306/SNPdensity_rows_location.txt") as handle:
#with open(sys.argv[1]) as handle:
#    reader=csv.reader(handle,delimiter='\t')    
#    for strLine in reader:
#        header.append(strLine[0])
#fout.write('\t'.join(header)+'\n')


#with open("/Users/phil/Desktop/SNPdensity_individual_306/SNPdensity_rows_location.txt") as handle:
with open(sys.argv[1]) as handle:
    reader=csv.reader(handle,delimiter='\t')
    scrambleLine = (next(reader))
    scrambleList = []
    males = []
    females = []
    for item in scrambleLine:
        index = scrambleLine.index(item)
        scrambleList.append(index)
        if "Female" in item:
            females.append(index)
        elif "Male" in item:
            males.append(index)           
all=males+females 

#with open("/Users/phil/Desktop/SNPdensity_individual_306/SNPdensity_rows_location.txt") as handle:
#with open(sys.argv[1]) as handle:
#    reader=csv.reader(handle,delimiter='\t')
#    next(reader)
#    every_m_v_f = []
#    for strLine in reader:
#        male_sum=0
#        for male in males:
#            male_sum+= float(strLine[male])
#        female_sum=0
#        for female in females:
#            female_sum+=float(strLine[female])
#        male_mean = male_sum/len(males)
#        female_mean = female_sum/len(females)
#        male_v_female = male_mean - female_mean
#        male_v_female_round = round(male_v_female, 5)
#        every_m_v_f.append(str(male_v_female_round))
#
#fout.write("true_m_v_f"+"\t"+'\t'.join(every_m_v_f)+'\n')

for i in range(0,int(perms)):
    with open(sys.argv[1]) as handle:
    #with open("/Users/phil/Desktop/SNPdensity_individual_306/SNPdensity_rows_location.txt") as handle:
        reader=csv.reader(handle,delimiter='\t')
        next(reader)
        every_m_v_f = []
        random.shuffle(all)
        perm_males = all[0:124]
        perm_females=all[124:]
        for strLine in reader:
            male_sum=0
            for male in perm_males:
                male_sum+= float(strLine[male])
            female_sum=0
            for female in perm_females:
                female_sum+=float(strLine[female])
            male_mean = male_sum/len(perm_males)
            female_mean = female_sum/len(perm_females)
            male_v_female = male_mean - female_mean
            male_v_female_round = round(male_v_female, 5)
            every_m_v_f.append(str(male_v_female_round))
    fout.write(sys.argv[4]+"_p"+str(i+1)+"\t"+'\t'.join(every_m_v_f)+'\n')
fout.close()