===============================
Step 1. Coverage-based analysis
===============================

To determine if there are large sex-specific regions present in your species of interest, a coverage-based analysis is first carried out.

``DifCover`` (https://github.com/timnat/DifCover) requires only a single male and single female ``bam`` file aligned to a reference genome.

Example ``DifCover`` command with 1 male (``SRR8585998_1.fastq.bam``) versus 1 female (``SRR8585999_1.fastq.bam``).

.. code-block:: console

    bash run_difcover.sh SRR8585998_1.fastq.bam SRR8585999_1.fastq.bam 1 >> difCov_Male98_over_Female99_outerr.txt 2>&1

Note that you will need to modify run_difcover.sh to have a proper path for ``FOLDER_PATH`` on your system, and the ``1`` is the library-specific Adjustment Coefficient (``AC``) to account for differences in sequencing depth between samples.
