===============================
Step 1. Coverage-based analysis
===============================

To determine if there are large sex-specific regions present in your species of interest, a coverage-based analysis is first carried out.

``DifCover`` (https://github.com/timnat/DifCover) requires only a single male and single female ``bam`` file aligned to a reference genome.

Example ``DifCover`` command with 1 male (``SRR8585998_1.fastq.bam``) versus 1 female (``SRR8585999_1.fastq.bam``).

.. code-block:: console

    bash run_difcover.sh SRR8585998_1.fastq.bam SRR8585999_1.fastq.bam 1 >> difCov_Male98_over_Female99_outerr.txt 2>&1

Note that you will need to modify ``run_difcover.sh`` to have a proper path for ``FOLDER_PATH`` on your system, and the ``1`` is the library-specific Adjustment Coefficient (``AC``) to account for differences in sequencing depth between samples. The correct ``AC`` value can be determined by taking the ratio of modal depths for the samples of interest through ``samtools stats`` as follows:

 .. code-block:: console

    chmod u+rwx samtools_Modaldepth.sh
    for file in $(ls *q.bam); do ./samtools_Modaldepth.sh $file; done

Once you have a modal depth for each sample you are interested in, ``AC`` is modal (coverage of sample 2 / modal coverage of sample 1). Of note, I have found that although this works well for some samples, others report a modal depth of 1, which is not helpful for determining ``AC``. In these cases, I have have mixed success simply using the ratio of the ``bam`` file sizes for ``AC``. This can be dialed in after the analysis as well by plotting and examining the center of the distribution (which should lie at 0 for autosomes if ``AC`` is properly set).
