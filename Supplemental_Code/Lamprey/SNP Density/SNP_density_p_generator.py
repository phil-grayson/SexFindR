#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 16:52:49 2020

@author: phil
"""

# execute with infile as $1, outfile as $2, number of permutations as $3

import csv
import sys

fileout=(sys.argv[2])
#fileout="/Users/phil/Desktop/test_permuttttter.txt"

fout= open(fileout, 'w') 
fout.write("LOCATION"+"\t"+"p_value"+"\n")

with open(sys.argv[1]) as handle:
    reader=csv.reader(handle,delimiter='\t')
    next(reader)
    for strLine in reader:
        true_value = float(strLine[1])
        if true_value > 0:
            #print("greater")
            hits = [element for element in strLine[2:] if float(element) > true_value]
            #print(len(hits))
            p=(len(hits)+1)/(len(strLine)-1) # p is hits plus our true value out of everything (minus the location column)  
        elif true_value < 0:
            #print("less")
            hits = [element for element in strLine[2:] if float(element) < true_value]
            #print(len(hits))
            p=(len(hits)+1)/(len(strLine)-1)
        else: # if true_value = 0, no need to calculate
            p=1
        fout.write(strLine[0]+'\t'+str(p)+'\n')
fout.close()
