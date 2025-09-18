02_ChIP_RNA_seq contains scripts and code:
- to start with RSEG calls from iPSC-NPC H3K9me3 data (provided in the RSEG folder) to call FXS H3K9me3 NPC BREACHes, 
- and produce plots in Figure 2, Figure 3, Figure S2, and Figure S4. 

001_call_BREACHes.sh
002_upset_plot.ipynb produces the upset plot in 2C.

genes_NPC_CV.txt contains the coefficient of variation per gene using the signal summations of H3K9me3 ChIP-seq signal across the FXS patient iPSC-NPCs for each gene using its TSS - 2kb upstream and + 10kb downstream.

NPC_H3K9me3_gene_signal_summation.txt contains the signal summation of H3K9me3 ChIP-seq signal for genes using each gene's TSS - 2kb upstream and + 10kb downstream.

DESEQ_melted.txt contains the RNA-seq data from the 7 iPSC-NPCs.

003_NPC_H3K9me3BREACHes_genes.ipynb produces the plots in 2D,E,F,I,J,G,H; S2D,F; 3C,F,D,G; S4B,D,C,E.
