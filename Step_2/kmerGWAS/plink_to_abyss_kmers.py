# -*- coding: utf-8 -*-
"""
Created on Fri Jan  3 09:58:32 2020

@author: Phil
"""

import csv
import sys

fileout=(sys.argv[2])
fout= open(fileout, 'w')

with open(sys.argv[1]) as handle:
    reader=csv.reader(handle,delimiter='\t')
    for strLine in reader:
        kmer = strLine[2]
        pval=strLine[9]
        fout.write(">"+pval+"\n"+kmer+"\n")
fout.close()

