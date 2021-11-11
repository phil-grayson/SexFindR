========================================
Step 3. Combined sequence-based analysis
========================================

If you have not been able to identify a consistent signal within your genomic data using Steps 1 and 2, the scripts for Step 3, Combined sequence-based analysis, might be able to help.

In Step 3, non-overlapping signals from the reference-based methods are converted to 10kb windows and combined to increase the user’s power to focus on specific genomic regions. Once candidate windows are identified in Step 3, the user can return to GWAS, SNP density, and ``Fst`` results to determine if fixed or nearly-fixed differences exist between the sexes in the top candidate windows.

Once again, the example code ``Fugu_SexFindR.R`` and necessary input files are included in the GitHub repository in order to run the analysis for fugu and generate the plot in the SexFindR figure in the manuscript.

.. image:: images/FuguSexFindR.png

*Figure 3. SexFindR Step 3 combined results identified the fugu SDR (red dotted line) alongside 4 additional candidate windows (thin black dotted lines), all on NC_042303.1.*
