===
Fst
===

``Fst`` is an index of allelic fixation in populations, and so, if there are high levels of ``Fst`` within discrete genomic regions when comparing males to females, this would suggest that those regions differ between males and females and, therefore, are not recombining. For the ``SexFindR`` workflow, ``Fst`` is calculated from a biallelic ``vcf`` file using ``vcftools``.

Example command
---------------

.. code-block:: console

    vcftools --gzvcf biallelic_filtered_PASS_fugu_14M_13F.vcf.gz --weir-fst-pop snpDen_males.txt --weir-fst-pop snpDen_females.txt --out biallelic_fst

Here, we make use of the same ``snpDen_males.txt`` and ``snpDen_females.txt`` reference files from the SNP Density step to assign individuals to their correct "population".
