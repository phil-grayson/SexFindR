===================================
Step 0. Mapping and variant calling
===================================

Prior to running Step 1 (Coverage-based analysis) or any of the analyses from Step 2 (Sequence-based analyses) except the reference-free k-mer analysis, you will need to map your reads to a reference genome.

In the SexFindR paper, we used ``Bowtie2`` (v 2.3.4.3 - http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) for read mapping and we provide some scripts and configuration for this within the GitHub and below. If you already have mapped reads from ``bwa-mem`` or another widely-used algorithm, please feel free to use those.

For all the sequence-based analyses mapped to a reference genome, SNP calling is also required. In the SexFindR paper, we used ``Platypus`` (commit 3e72641; https://github.com/andyrimmer/Platypus) to jointly call SNPs across all samples, and we provide some scripts and configurations for this within the GitHub and below. Again, if you have already called SNPs for your samples using ``GATK`` or another widely-used algorithm, please feel free to use those.

Download raw reads from NCBI
----------------------------
Requires `fasterq-dump` from the `SRA-tools` (https://github.com/ncbi/sra-tools)

.. code-block:: console

    for file in $(cat acc_list.txt); do echo $file; date; fasterq-dump $file -e 36; done

Download the genome
-------------------
.. code-block:: console

    wget GCF_901000725.2_fTakRub1.2_genomic.fna.gz

Create a bowtie2 index
----------------------
.. code-block:: console

    bash bowtie2_makeindex_linux.sh GCF_901000725.2_fTakRub1.2_genomic.fna fugu &> index_outerr.txt &

Map reads to the genome
-----------------------
The following must be executed for each sample:

.. code-block:: console

    bash bowtie2_16_linux.sh SRR8585991_* fugu &> bt2_SRR8585991_outerr.txt &

Call variants using Platypus
----------------------------
The following was run on a SLURM system:

.. code-block:: console

    sbatch platypus_all_region_1day.sh 15M_14F_bams.txt GCF_901000725.2_fTakRub1.2_genomic.fna
