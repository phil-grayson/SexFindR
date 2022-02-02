=================
Supplemental Code
=================

This section is a work in progress. Code for the analysis of lamprey and the other species from the main manuscript on the GitHub repository and documentation (when necessary or helpful) will be added.

Running VCF Search
------------------

The VCF Search described in the supplements was carried out on fugu as a positive control and on lamprey to look for fixed or nearly-fixed differences between males and females across the entire dataset as well as within individual lakes. In-house ``Python`` scripts were written for this purpose and can be found in the ``VCF Search`` folders for lamprey and fugu. The ``Python`` scripts require an uncompressed vcf, which is too large to include on GitHub, but I have subsampled the fugu vcf to include only the chromosome that contains the SDR so that the fugu script can be used for demonstation purposes.

Subsetting the fugu ``vcf`` file:

.. code-block:: console

    grep "^#" ~/Desktop/fugu/filtered_PASS_fugu_14M_13F.vcf > fugu_14M_13F.vcf
    grep -v "^#" ~/Desktop/fugu/filtered_PASS_fugu_14M_13F.vcf | grep "NC_042303.1" >> fugu_14M_13F.vcf
    gzip fugu_14M_13F.vcf

Running VCF Search

.. code-block:: console

    gunzip fugu_14M_13F.vcf.gz
    python3 fugu_vcf_search.py fugu_14M_13F.vcf out_fugu_vcf_search_test.txt

Running 100k Permutations
-------------------------

For SNP density, the ``R`` code provided runs 1000 permutations, but I wanted to run more. For fugu and lamprey, custom ``Python`` and ``R`` code was written to carry out the permuations and pivot the resulting table to a useful format. Below is the example code for the lamprey run.

.. code-block:: console

    cd ~/SexFindR/Supplemental_Code/Lamprey/SNP\ Density/
    unzip SNPdensity_rows_location.txt.gz
    for run in $(cat runs.txt); do python SNP_density_permuter.py SNPdensity_rows_location.txt perm_${run}.txt 2500 $run & sleep 0.1; done

    cat out_true_SNPden.txt perm* > true_plus_100k_perms.txt

    # using R to pivot the dataframe
    R
    long <- read_tsv("true_plus_100k_perms.txt") 
    pivot <- t(long) 
    output <- as.data.frame(pivot) %>% rownames_to_column()
    write_tsv(output,"trans_true_plus_100k_perms.txt",col_names = F)

    python3 SNP_density_p_generator.py trans_true_plus_100k_perms.txt SNP_density_100k_permutations_p_values.txt &
