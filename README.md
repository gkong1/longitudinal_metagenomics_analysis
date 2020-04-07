# longitudinal_metagenomics_analysis
The code for the analysis of the fecal metagenomics data are provided here.

## Plasma metabolomics analysis
The plasma of the same mice cohort as the longitudinal shotgun sequencing experiment were extracted at 12 weeks of age with terminal cardiac puncture and LC-MS metabolomics were performed to determine metabolomic profile in the plasma, comparing WT and HD mice. The median-normalized data is provided in the file "Normalized_metabolites_table.txt", with PBQCs as internal control included in the file. The method of analysis is provided in the file "metabolites_analysis.R" which details the code for PCA, sPLS-DA and various plots for visualization.

## Short chain fatty acid (SCFA) assay analysis
For the independent SCFA assay analysis: the normalized concentrations of acetate, butyrate and propionate of WT and HD mice at 12 weeks of age were included in the file "matrix with FCS_conc.csv". The samples for the SCFA assay were the same ones used for the plasma metabolomics. The analysis file for the SCFA data ("scfa_analysis") includes the code for one-way ANOVA statistical test as well as for plotting barplots. 
