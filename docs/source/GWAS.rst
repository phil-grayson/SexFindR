====================================
Genome Wide Association Study (GWAS)
====================================

GWAS identifies associations between a phenotype and a genotype (Klein et al. 2005). Similar to ``Fst``, by carrying out a GWAS with the male and female populations as two different phenotypes, it is possible to identify SNPs that are strongly or weakly associated with sex. These SNPs can be fixed, or nearly fixed, in either males or females.

This analysis requires ``vcftools`` and a combination of ``plink`` (v 1.07-x86_64; https://www.cog-genomics.org/plink/1.9/) and ``GEMMA`` (v0.98.1; https://github.com/genetics-statistics/GEMMA). We tested numerous GWAS algorithms and this one performed best across the sex chromosome divergence continuum.

Filtering the vcf
-----------------

We use ``vcftools`` to convert the ``vcf`` to ``plink`` format, remove indels, remove sites with >50% missing data, and sites where the minor allele frequency is less than 5% or more than 95%.

.. code-block:: console

    vcftools --vcf filtered_PASS_fugu_14M_13F.vcf --plink --remove-indels --max-missing 0.5 --max-maf 0.95 --maf 0.05 --out gwas_fugu

Running the GWAS
----------------

As outlined in the ``GEMMA`` manual (https://github.com/genetics-statistics/GEMMA/blob/master/doc/manual.pdf), the program requires input files in the ``plink`` binary format. We first run ``plink`` (as recommended) to generate these files and then supply these files to ``GEMMA`` alongside the option ``-lm 2`` to specify a likelihood ratio test.

.. code-block:: console

    /home/pgrayson/programs/plink-1.07-x86_64/plink --file gwas_fugu --pheno sex_fugu_meta.txt --make-bed --out gwas_fugu_plink --noweb --allow-no-sex
    ~/programs/gemma-0.98.1-linux-static -bfile gwas_fugu_plink -lm 2 -o fugu_gemma_out

``sex_fugu_meta.txt`` is included in the GitHub repository as an example file. It contains 3 tab-delimited columns that contain the sample name (repeated in column 1 and column 2) and the phenotype (1 for male, 2 for female - or vice versa). An example is below:

.. code-block:: console

    SRR8585991_1.fastq      SRR8585991_1.fastq      1
    SRR8585992_1.fastq      SRR8585992_1.fastq      2
    SRR8585993_1.fastq      SRR8585993_1.fastq      2

Analysing the output
--------------------

The commands above will generate ``fugu_gemma_out.assoc`` and ``fugu_gemma_out.log.txt``. The ``.assoc`` file contains the ``GWAS`` results, including the genome position (``rs``) the number of individuals that had a call (``n_obs``) or missing data (``n_mis``) at that site, the allele frequency (``af``), and the ``p-value`` for the likelihood ratio test (``p_lrt``). This file can be searched and parsed to identify the top candidates from this analysis.
