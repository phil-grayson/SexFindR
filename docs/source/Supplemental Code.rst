=================
Supplemental Code
=================

This section is a work in progress. Code for the analysis of lamprey and the other species from the main manuscript on the GitHub repository and documentation (when necessary or helpful) will be added.

Running VCF Search
------------------

The VCF Search described in the supplements was carried out on fugu as a positive control and on lamprey to look for fixed or nearly-fixed differences between males and females across the entire dataset as well as within individual lakes. In-house python scripts were written for this purpose and can be found in the ``VCF Search`` folders for lamprey and fugu. The python scripts require an uncompressed vcf, which is too large to include on GitHub, but I have subsampled the fugu vcf to include only the chromosome that contains the SDR so that the fugu script can be used for demonstation purposes.

Subsetting the fugu vcf file:

.. code-block:: console

    grep "^#" ~/Desktop/fugu/filtered_PASS_fugu_14M_13F.vcf > fugu_14M_13F.vcf
    grep -v "^#" ~/Desktop/fugu/filtered_PASS_fugu_14M_13F.vcf | grep "NC_042303.1" >> fugu_14M_13F.vcf
    gzip fugu_14M_13F.vcf

Running VCF Search

.. code-block:: console

    gunzip fugu_14M_13F.vcf.gz
    python3 fugu_vcf_search.py fugu_14M_13F.vcf out_fugu_vcf_search_test.txt
