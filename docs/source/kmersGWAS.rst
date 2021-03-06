=========
kmersGWAS
=========

The ``kmersGWAS`` approach described below is possibly the most powerful single analysis within the ``SexFindR`` protocol. Apart from not requiring a reference genome, this single method has proved capable of detecting sex-linked sequences across the entire sex chromosome evolutionary continuum, from highly degenerate mammalian systems to takifugu, which has a single base that is heterozygous is males and homozygous in females.

After testing a few other k-mer association algorithms as well as a number of specific options within this algorithm, I settled on the following workflow given its consistent results across the sex chromosome divergence continuum. This protocol requires ``kmersGWAS`` (v0.2 beta; https://github.com/voichek/kmersGWAS), ``ABYSS`` (v2.2.5; https://github.com/bcgsc/abyss), ``plink`` (v 1.07-x86_64; https://www.cog-genomics.org/plink/1.9/), ``R`` and ``blast+`` ( v2.10.0; https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.10.0/). I have presented the necessary steps and configurations for the fugu analysis below, but it should be noted that run times and memory usage vary greatly with genome size, coverage, and sample number. For fewer than 30 samples of fugu, the config provided in the SLURM scripts should be sufficient, but some of the ``kmersGWAS`` steps required days to run on large machines (200 GB+ memory) for lamprey, and other species.

Directory structure
-------------------

``kmersGWAS`` requires a specific directory structure and has many separate steps. These are presented below, but more information can be found in the kmersGWAS manual here (https://github.com/voichek/kmersGWAS/blob/master/manual.pdf).

Generally, you want the following:

1. Files present in a top level directory
2. Individual directories created within the top level directory corresponding to the base name of each sample (e.g., top level directory contains ``SRR8585991_1.fastq.gz``, ``SRR8585991_2.fastq.gz``, and the ``SRR8585991`` directory)
3. Within the base name subdirectories (e.g., ``SRR8585991``), there are execution scripts and an ``input_files.txt`` file that contains the paths to the ``fastq`` files (e.g., ``../SRR8585991_1.fastq.gz`` and ``../SRR8585991_2.fastq.gz`` each on their own line).

Running kmersGWAS
-----------------

The following is a set of example commands to set up the directory structure for the fugu samples, but these will need to be made specific to your file names, path, etc. Many of these steps were executed on a SLURM system and the execution scripts will need to be modified to match your system as well. This process begins with the ``fastq`` files and one copy of each of the execution scripts in the top level (``fugu_kmerGWAS``) directory. All commands are executed from that top level directory. I have included an example of all necessary text files in the GitHub repository to check against as well.

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

We next need the individuals k-mers list file. As above, these should be generated following the ``kmersGWAS`` manual. An example set of commands used to generate this file for the fugu samples is provided below:

.. code-block:: console

    cd /home/pgrayson/
    ll scratch/fugu_kmerGWAS/ | tail -n +2 | awk '{printf "/scratch/fugu_kmerGWAS/%s/kmers_with_strand\t%s\n", $NF,$NF}' > kmers_list_paths.txt
    grep -v .gz kmers_list_paths.txt | grep SRR > kmers_clean_list_path.txt
    sed 's/\/scratch/\/home\/pgrayson\/scratch/g' kmers_clean_list_path.txt > final_kmers_clean_list_path.txt
    mv final_kmers_clean_list_path.txt ~/scratch/fugu_kmerGWAS


Example execution scripts are provided to combine the k-mers and create the k-mers table.

.. code-block:: console

    sbatch combine_kmersGWAS.sh
    sbatch create_table_kmersGWAS.sh

Once these steps are completed, we generate ``plink`` binary files using the ``kmers_table_to_bed`` function. The ``phenotype.pheno`` file here contains two tab-delimited columns, ``accession_id`` and ``phenotype_value`` with the sample ID in ``accession_id`` and ``1`` or ``2`` in the ``phenotype value`` column for the male or female phenotype. ``-b`` here is batch size. Depending on the power of your machine and the number of samples included in your analysis, you might have to produce smaller batches (I ran on a Linux machine with 256 Gb of memory and had to batch other species with more samples). Batches here are the number of k-mers within the analysis for ``plink``.

.. code-block:: console

    ~/programs/kmerGWAS/bin/kmers_table_to_bed -t kmers_table -k 31 -p phenotype.pheno --maf 0.05 --mac 5 -b 1000000000 -o fugu_kmerGWAS_plink

Running plink
-------------

Next, we run these files through ``plink`` to obtain ``p-values`` for each k-mer, representing the association between k-mer presence and phenotypic sex. If you had to batch your k-mers (above), you will also need to execute this command for batch the other batches (e.g., fugu_kmerGWAS_plink.1, fugu_kmerGWAS_plink.2, etc.). Make sure to also append something new to the ``--out`` name (e.g., ``fugu_kmers1``, ``fugu_kmers2``), so that previous results are not overwritten.

.. code-block:: console

    ~/programs/plink-1.07-x86_64/plink --noweb --bfile fugu_kmerGWAS_plink.0 --allow-no-sex --assoc --out fugu_kmers

The resulting outfile (e.g., ``fugu_kmers.assoc``) can explored and parsed based on the ``P`` column to identify and pull out the k-mers that have the highest association with sex. An example for fugu to obtain the k-mers with the most significant ``p-values`` is:

.. code-block:: console

    awk '$9 < 0.000000000001' fugu_kmers.assoc > most_significant_fugu_assoc.txt

Results from different batches can be concatenated with ``cat fugu_kmers*.assoc`` or something similar following analysis in ``plink`` (although a simple ``cat`` will also repeat the header line from each batch, so be aware that you will need to take this into account using ``tail`` or parsing the output file from ``cat`` to remove the extra header lines).

Running ABYSS
-------------

Once you are happy with the filtered k-mer set, you can use ``ABYSS`` to assemble those k-mers into small contigs for ``blastn`` analysis (if a reference genome exists).  A basic ``python`` script (``plink_to_abyss_kmers.py``) is included to parse the filtered ``plink`` output into ``ABYSS``-ready input. The ``ABYSS`` input should be in fasta format, with the ``p-value`` as the sample ID and the k-mer as the sequence. e.g.,

.. code-block:: console

    >2.005e-13
    AAAAAAAAAAAATCATTTCCCACCTCATCAA
    >2.005e-13
    AAAAAAAAAAATCATTTCCCACCTCATCAAT
    >2.005e-13
    AAAAAAAAAATCATTTCCCACCTCATCAATC
    ...

To generate this file, the following should work:

.. code-block:: console

    python plink_to_abyss_kmers.py most_significant_fugu_assoc.txt fugu_plink_abyss_input.txt

If your output does not match the example above, you might need to change the index positions in the python script to correctly grab the k-mer and the p-value columns (given different ``plink`` versions). This input file is then run through ABYSS to assemble small contigs.

.. code-block:: console

    ABYSS -k25 -c0 -e0 fugu_plink_abyss_input.txt -o fugu_plink_abyss_output.txt

Running blastn
--------------

If you have a reference genome, you might be interested in where the enriched k-mers map to. The easiest way to check this is through a ``blastn`` analysis.

You first need to create a blastdb:

.. code-block:: console

    makeblastdb -in GCF_901000725.2_fTakRub1.2_genomic.fna -dbtype nucl

Then you can run blastn:

.. code-block:: console

    blastn -query fugu_plink_abyss_output.txt -db GCF_901000725.2_fTakRub1.2_genomic.fna -outfmt 6

The output file from blastn will provide top candidate regions for each contig that was assembled from the k-mers. In the case of fugu, there are only 4 contigs that all blast to NC_042303.1, but in poplar and the golden monkey, these blast outputs were parsed with an ``R`` script (e.g., ``poplar_kmerGWAS_blast_results.R``) to visualize which genomic regions the sex-assocaited k-mers map to. This script and the necessary input files have been included in the GitHub repo for poplar.
