# Summary
The `FSRnSimBS_sva` function integrates variable selection (FSR) with surrogate variable analysis (SVA), incorporating two approaches: FSR_sva (combined FSR-SVA) and SVAall_FSR (combined SVA_FSR). The function also incorporates the SVA0 method (surrogate variable estimation using primary variable only). Performance is evaluated across 100 replications using two metric categories: (1) variable selection (average false selection rate, FSR; number of relevant covariates, R) and (2) differential expression analysis (empirical FDR; average true positives, NTP; partial AUC, PAUC) (3) Number of surrogate variables estimated, mean and the standard deviation of the of the $R^2$ values from linear regression models where the hidden unselected relevant covariates are the explanatory variable and the response is the estimated surrogate variables. 

# Data Analysis and Simulation Based on the RFI RNA-seq Dataset
- Codes for RFI RNA-seq data analysis are [here](https://github.com/Farzana001-Noorzahan/covariate-selection-hidden-factor-rnaseq/blob/main/Analysis/rcodes/Real_analysis_RFI.R).
- Codes for simulation study are [here](https://github.com/Farzana001-Noorzahan/covariate-selection-hidden-factor-rnaseq/blob/main/Analysis/rcodes/4-Simulation-sva.R) and [here](https://github.com/Farzana001-Noorzahan/covariate-selection-hidden-factor-rnaseq/blob/main/Analysis/rcodes/4-Simulation-sva.sh)
 - Codes to reproduce the tables and figures for the RFI RNA-Seq data analysis is [here](https://github.com/Farzana001-Noorzahan/covariate-selection-hidden-factor-rnaseq/blob/main/Analysis/rcodes/MDPI_2.R)

# Data Analysis and Simulation Based on the Zebrafish RNA-seq Dataset
- Zebrafish RNA-seq data is  [here](https://github.com/Farzana001-Noorzahan/covariate-selection-hidden-factor-rnaseq/blob/main/Analysis/extra-rna-seq-data).
- Codes for simulation study are [here](https://github.com/Farzana001-Noorzahan/covariate-selection-hidden-factor-rnaseq/blob/main/Analysis/rcodes/4-Simulation-sva-zebra.R).and [here](https://github.com/Farzana001-Noorzahan/covariate-selection-hidden-factor-rnaseq/blob/main/Analysis/rcodes/4-Simulation-sva-zebra.sh).
 - Codes to reproduce the tables and figures for the RFI RNA-Seq data analysis is [here](https://github.com/Farzana001-Noorzahan/covariate-selection-hidden-factor-rnaseq/blob/main/Analysis/rcodes/MPDI_2 _Zebra.R)
