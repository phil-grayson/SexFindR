===============================
Step 2. Sequence-based analyses
===============================

Fst, SNP Density, GWAS, K-mer

kmerGWAS
--------

After testing a few other algorithms and other options within this algorithm, I settled on the following workflow given its consistent results across the sex chromosome divergence continuum.  Requires ``kmerGWAS`` (v0.2 beta; https://github.com/voichek/kmersGWAS), ``ABYSS`` (https://github.com/bcgsc/abyss), and ``plink`` (v 1.07-x86_64; https://www.cog-genomics.org/plink/1.9/).

``kmerGWAS`` requires a very specific directory structure and more information can be found her (https://github.com/voichek/kmersGWAS/blob/master/manual.pdf).

Generally, you want the following:

1. Files present in a top level directory
2. Individual directories created within the top level directory corresponding to the base name of each sample (e.g., top level directory contains ``SRR8585991_1.fastq.gz``, ``SRR8585991_2.fastq.gz``, and the ``SRR8585991`` directory)
3. Within the base name subdirectories, there are execution scripts and an ``input_files.txt`` file that contains the paths to the ``fastq`` files (e.g., ``../SRR8585991_1.fastq.gz`` and ``../SRR8585991_2.fastq.gz`` each on their own line).

The following is a set of example commands to set up the directory structure for the fugu samples, but these will need to be made specific to your file names, path, etc. This was executed on a SLURM system and the execution scripts will need to be modified to match your system as well.

.. code-block:: console

    cd /home/pgrayson/scratch/fugu_kmerGWAS
    for file in $(ls --color=none *_1.*gz); do baseName=$(echo $file | awk -F'[_]' '{print $1}'); mkdir $baseName; done
    ls --color=none -d */ > dirlist.txt
    sed 's/\///g' dirlist.txt > clean_dirlist.txt
    for dir in $(cat clean_dirlist.txt); do cd $dir; echo $dir; ls --color=none ../*${dir}*gz > input_files.txt; cd ..; done
    for dir in $(cat clean_dirlist.txt); do cp step1_kmerGWAS.sh $dir; done
    for dir in $(cat clean_dirlist.txt); do cd $dir; sbatch step1_kmerGWAS.sh; cd ..; done

Once those first jobs complete, you can run the following (again with modification to the execution script):

.. code-block:: console

    for dir in $(cat clean_dirlist.txt); do cp strand_kmerGWAS.sh $dir; done
    for dir in $(cat clean_dirlist.txt); do cd $dir; sbatch strand_kmerGWAS.sh; cd ..; done
