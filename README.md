(Forked to test out some survival models using linked shrinkage)

# ThinkInteractions
R scripts and demo's for fitting a regression model with linked shrinkage.

Update 2-10-2024: Repaired a coding error in the Shapley calculation. This implies resulting Shapley plots may deviate from those in the manuscript. 

Scripts: Mark van de Wiel, mark.vdwiel@amsterdamumc.nl.
Software implements a novel linked shrinkage model for two-way interactions, as presented in manuscript: "Linked shrinkage to improve estimation of interaction effects in regression models."
In addition, it implements the Shapley values derived in the manuscript.

**Data**
1. datasynth_Chol.Rdata: synthethic copy of the Helius data used in the manuscript, systolic blood pressure (SBP) as outcome. Rdata file containing the data.frame datasynth with 21,570 rows (samples) and 100 columns; first is the response, 2-15 the main effects, 16-100 the two way interactions.
2. datasynth_Sbp.Rdata: As datasynth_Chol.Rdata, but with Cholesterol as outcome (first column)

**Scripts**
1. ThinkInteractionsDemo.R: Demo script illustrating use of our method, Bayint, to fit regression models with two-way interactions. Script also illustrates tools and plots to intepret the model, in particular Shapley values (feature importance). Moreover, it demonstrates how the alternative models were used as comparison (including OLS, two-step, lasso en ridge variations, Bayesian local regression and variations on Bayint). ThinkInteractionsDemo.html contains R markdown compiled report.
2. Synth_RMSE.R. Scripts that allows users to reproduce our computations on the synthetic data, in particular the root MSEs of the parameter estimates as produced by various methods. Time-consuming as it runs all methods on 25 subsamples of size n=1,000 of the entire set.
3. ThinkInteractionsLogistic: Example for logistic regression;  ThinkInteractionsLogistic.html contains R markdown compiled report.
4. Auxiliary files auxiliarycodeRstan.R and auxiliarycodeOther.R containing Rstan and other source functions, respectively. These are invoked by the source() command in the other scripts.  

