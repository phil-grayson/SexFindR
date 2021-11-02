===================================
Step 0. Mapping and variant calling
===================================

Prior to running Step 1 (Coverage-based analysis) or any of the analyses from Step 2 (Sequence-based analyses) except the reference-free k-mer analysis, you will need to map your reads to a reference genome.

In the SexFindR paper, we used ``Bowtie2`` (v 2.3.4.3; http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) for read mapping and we provide some scripts and configuration for this within the GitHub and below. If you already have mapped reads from ``bwa-mem`` or another widely-used algorithm, please feel free to use those.

For all the sequence-based analyses mapped to a reference genome, SNP calling is also required. In the SexFindR paper, we used ``Platypus`` (commit 3e72641; https://github.com/andyrimmer/Platypus) to jointly call SNPs across all samples, and we provide some scripts and configurations for this within the GitHub and below. Again, if you have already called SNPs for your samples using ``GATK`` or another widely-used algorithm, please feel free to use those.

``Platypus`` requires that the genome file be indexed with ``samtools`` (v1.10; https://github.com/samtools/samtools). We also filtered the raw ``vcf`` to only include sites that ``PASS`` the quality filters. This requires ``bcftools`` (v1.9; https://github.com/samtools/bcftools). ``vcftools`` is also used to filter for biallelic sites (v0.1.14; https://github.com/vcftools/vcftools).

Download raw reads from NCBI
----------------------------
Requires ``fasterq-dump`` from the ``SRA-tools`` (https://github.com/ncbi/sra-tools). In my experience, ``fasterq-dump`` can be fairly flakey depending on your system and connection, so this download might need to be repeated multiple times to get all the samples.  You might want to separate out the ``acc_list.txt`` files or add something like ``if [ ! -f ${file}_1.fastq ]`` to the code below to avoid downloading the same file repeatedly.

.. code-block:: console

    for file in $(cat acc_list.txt); do echo $file; date; fasterq-dump $file -e 36; done

Download the genome
-------------------
.. code-block:: console

    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/901/000/725/GCF_901000725.2_fTakRub1.2/GCF_901000725.2_fTakRub1.2_genomic.fna.gz
    gunzip GCF_901000725.2_fTakRub1.2_genomic.fna.gz

Create a bowtie2 index
----------------------
.. code-block:: console

    bash bowtie2_makeindex_linux.sh GCF_901000725.2_fTakRub1.2_genomic.fna fugu &> index_outerr.txt &

Map reads to the genome
-----------------------
The following must be executed for each sample:

.. code-block:: console

    bash bowtie2_16_linux.sh SRR8585991_* fugu &> bt2_SRR8585991_outerr.txt &

``bowtie2_16_long.sh`` is also included as a reference on SLURM systems.

Call variants using Platypus
----------------------------
Variants are jointly called (all at once) through the use of a bam list (e.g., 15M_14F_bams.txt). The following was run on a SLURM system:

.. code-block:: console

    samtools faidx GCF_901000725.2_fTakRub1.2_genomic.fna
    sbatch platypus_all_region_1day.sh 15M_14F_bams.txt GCF_901000725.2_fTakRub1.2_genomic.fna

Filter for calls that PASS quality filters
------------------------------------------
This script will keep only those variants that have ``PASS`` in the ``FILTER`` field, removing low quality calls.

.. code-block:: console

    sbatch bcf_filter.sh fugu_14M_13F.vcf

Filter for biallelic sites
--------------------------
Some downstream analyses (e.g., ``Fst``) require only biallelic sites to work properly (e.g., 0/0, 0/1, 1/1).

.. code-block:: console

    vcftools --vcf filtered_PASS_fugu_14M_13F.vcf --max-alleles 2 --stdout --recode --recode-INFO-all | gzip -c > biallelic_filtered_PASS_fugu_14M_13F.vcf.gz
