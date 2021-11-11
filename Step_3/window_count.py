# -*- coding: utf-8 -*-
"""
Created on Wed Dec 18 09:49:52 2019

@author: Phil
"""

# first argv is the input file
# second argv is the outfile
# third argv is our window size

import csv
import sys

fileout=(sys.argv[2])
fout= open(fileout, 'w') 
fout.write("scaf"+'\t'+"base"+'\t'+"count"+'\n')

# get first scaffold
f = open(sys.argv[1])
scaf = f.readline().split('\t')[0]
f.close()

with open(sys.argv[1]) as handle: 
    reader=csv.reader(handle,delimiter='\t')
    window_size = int(sys.argv[3])  
    window_count = 0
    for strLine in reader:
        if strLine[0] == scaf: #if scaf is the same as previous line
            if int(strLine[1]) <= window_size: #if base position is within the current window
                window_count+=1 # add a count
            else: # if base is not within window
                if window_count > 0: # if we have counted any, write that out
                    fout.write(scaf+'\t'+str(window_size)+'\t'+str(window_count)+'\n')
                    window_count=0
                if int(strLine[1]) < window_size + int(sys.argv[3]): # if it's just in the next window over, add a window size to the window
                    window_size = window_size + int(sys.argv[3])
                else: # if it's more than a window over, we need to decide to round up or down to get the next window coordinates
                    if int(strLine[1]) > int(round(int(strLine[1]), -4)):
                        window_size = int(round(int(strLine[1]), -4)) + int(sys.argv[3])
                    else:
                        window_size = int(round(int(strLine[1]), -4))
                window_count+=1 # either way, add a count to the window
        else: # when we run out of scaffold, we write it out, go to the next one and 
           fout.write(scaf+'\t'+str(window_size)+'\t'+str(window_count)+'\n')
           if int(strLine[1]) > int(round(int(strLine[1]), -4)):
               window_size = int(round(int(strLine[1]), -4)) + int(sys.argv[3])
           else:
               window_size = int(round(int(strLine[1]), -4))
           window_count = 1
           scaf = strLine[0]
fout.close()
