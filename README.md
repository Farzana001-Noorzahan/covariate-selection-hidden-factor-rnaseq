# Summary
The `FSRnSimBS_sva` function integrates variable selection (FSR) with surrogate variable analysis (SVA), incorporating two approaches: FSR_sva (combined FSR-SVA) and SVAall_FSR (combined SVA_FSR). The function also incorporates the SVA0 method (surrogate variable estimation using primary variable only). Performance is evaluated across 100 replications using two metric categories: (1) variable selection (average false selection rate, FSR; number of relevant covariates, R) and (2) differential expression analysis (empirical FDR; average true positives, NTP; partial AUC, PAUC) (3) Number of surrogate variables estimated, mean and the standard deviation of the of the $R^2$ values from linear regression models where the hidden unselected relevant covariates are the explanatory variable and the response is the estimated surrogate variables. 

# Algorithm of the Integrated Methods
## FSR$\_$sva
\begin{enumerate}

\item[1.] Read the count data, which is a matrix containing genes as rows ($g=1,2,\ldots,G$) and observations as columns ($i=1,2,\ldots,n$), along with the FixCov (the primary variable of interest) and the VarCov (variables subject to selection) in data frames.

\item[2.] Keep only the genes for which at least two observations have at least five counts. Given the library size $R_i$ as the 75th percentile for each sample i, we transform the filtered counts as log counts using the equation below:
\begin{equation*}
y_{gi} = \log_2\left( \frac{c_{gi} + 0.5}{R_i + 1} \times 10^6\right).
\end{equation*}

\item[3.] The log counts, along with the FixCov and VarCov, are passed as arguments to the `FSRAnalysisBS` function in the `csrnaseq` package in R, which performs the backward variable selection method described in Section 2.2. This function returns the BestCovOut element containing the selected covariates using the FSR variable selection method. We choose the option as `OWN`, B = 100 pseudo variables and m = 3 generated pseudo variables as the default setting in the `FSRAnalysisBS` function.

\item[4.] The full model is generated using the `model.matrix` function with the FixCov (variable of interest) and the BestCovOut (adjustment variables), and the null model is created similarly using only the adjustment variables. The `sva` function in the `sva` package in R is then used to estimate the number of surrogate variables and compute the surrogate variables. The estimated surrogate variables along with the selected variables from FSR are combined and stored together as a data frame.

\item[5.] Differential expression analysis is then conducted using the `DEAeval` function in the `csrnaseq` package in R to calculate evaluation metrics including the number of true positives (NTP), false discovery rate (FDR), and the average partial area under the ROC curve (PAUC) for the combined set of variables from step 5.

\end{enumerate}
