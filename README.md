# Causal-Gaussian-DAGs
This repository contains all the R codes implemented for my thesis "Causal inference from observational studies: a graphical model perspective".


- auxiliary_functions.R -> implements all the functions needed to estimate the multiset of causal effects associated with an estimated CPDAG followi the approach proposed in the thesis. 
- final_function.R -> implements the function to estimate the causal effect as function of the weight matrix L only (or B:= I-L).
- Frequentist_simulation_2.4.R -> contains the R code to run the simulation study in Section 2.4. Methods under comparisons are: IDA, optimal IDA, and the method developed in the thesis for which the causal effect is retrieved as a function of L.
- bayesian_comparison_3.2.R -> contains the R code to run the simulation study performed in Section 3.2. We empirically obtain an MSE distribution for: the average causal effect estimated via optimal adjustment, the average causal effect estimated with our new methodology f(L)-path, and a bayesian estimator based on the posterior mean of causal effect. 
