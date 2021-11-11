===========
SNP Density
===========

In nascent sex-linked sequences, we would expect to see differences in overall SNP density between males and females due to Y reads still mapping to the X (or W to the Z). This is due to the divergent evolutionary trajectories of the X and Y following recombination suppression on the Y. SNP density can identify regions where male or female-specific SNPs are still segregating within their respective populations.

This analysis requires ``vcftools`` and ``R``.

Running SNP Density
-------------------

We use ``vcftools SNPDensity`` to calculate SNP Density for each sample in 10kb windows across the genome.  This requires individual ``vcf`` files for each sample. The individual files are created and the ``SNPDensity`` calculations carried out using the commands below.


.. code-block:: console

    mkdir individual_SNP_density

    for file in $(cat snpDen_females.txt); do sbatch SNPdensity.sh $file individual_SNP_density/Female_${file%%.*}.vcf biallelic_filtered_PASS_fugu_14M_13F.vcf.gz; sleep 0.1; done

    for file in $(cat snpDen_males.txt); do sbatch SNPdensity.sh $file individual_SNP_density/Male_${file%%.*}.vcf biallelic_filtered_PASS_fugu_14M_13F.vcf.gz; sleep 0.1; done

Here, ``snpDen_males.txt`` and ``snpDen_females.txt`` reference files are simply lists of the sample IDs (as found in the ``vcf`` header) that contain the sample from those sexes. The submission script above appends ``Male`` or ``Female`` to the ``SNPDensity`` output file, which we then use to parse these output files in ``R``.


Analyzing SNP Density
---------------------

SNP density can be analyzed from the ``R`` script ``SNPdensity_permutations_fugu.R``. This file runs the user through reading in the individual SNP densities calculated above, then calculates the true male and female means (and differences) for each 10 kb window, before carrying out a small permutation test. This permutation test generates the p-values that allows us to determine the correct order for SNP Density windows in ``Step 3``. A larger permutation test was carried out for both the lamprey and fugu in the main manuscript to further analyze the SNP Density data. Code and documentation for this test is included in the ``Supplemental Code`` folder on the GitHub repo.
