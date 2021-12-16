#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 16:52:49 2020

@author: phil
"""

# execute with filtered_PASS_sea_lamprey_304_update.vcf as $1, outfile as $2, site as $3

# options for $3 are Erie, Huron, Michigan, Ontario, Superior, Cayuga, Phil, Holyoke, Seneca (but we won't do Seneca since there is only one female)

import csv
import sys

waters = dict(Erie = ["BC","GO"],Huron = ["Augres","BL","CH","SM"],Michigan = ["GR","MAN"],Ontario = ["BLR","DF"],Superior = ["BR","RO","TQ"],Cayuga = ["NY"],Phil = ["d_"],Holyoke = ["CTR"],Seneca = ["CC"])


# get male names to index from the header of vcf
males=[]
females=[]
with open("sexed_update_for_python.txt") as handle:
    reader=csv.reader(handle,delimiter='\t')
    water = (sys.argv[1])
    population=waters[water]
    for strLine in reader:
        for code in population:
            if code in strLine[0]:
                if strLine[1] == "M":
                    males.append(strLine[0])
                if strLine[1] == "F":
                    females.append(strLine[0])

print("males")
print(len(males))

print("females")
print(len(females))
