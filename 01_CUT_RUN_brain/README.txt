O1_CUT_RUN_brain contains scripts and code to start with RSEG calls from brain H3K9me3 data CUT&RUN (provided in the RSEG folder) to call FXS H3K9me3 brain BREACHes and produce plots in Figure 1. 

001_call_BREACHes.sh
002_upset_plot.ipynb produces the upset plot in 1B.

genes_brain_CV.txt contains the coefficient of variation per gene using the signal summations of H3K9me3 CUT&RUN signal across the FXS patient brain samples for each gene using its TSS - 2kb upstream and + 10kb downstream.

003_brainH3K9me3BREACHes_genes.ipynb produces the plots in 1 C,D,E.
