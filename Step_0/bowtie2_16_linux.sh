#!/bin/bash

bowtie2 -x $3 -1 $1 -2 $2 -X 2000 -p 16 | samtools view -b -S - |samtools sort - -o $1.bam

samtools index $1.bam
