This repository covers a good portion of the coding done for my Statistics Simulation Class (STAT 624 at BYU). The main premise was to create a Monte Carlo study evaluating how three unique imputation methods handle missing data, and if the methods utilizes and values imputed influence the integrity of the dataset. 

The methods used are as follows:
- List-Wise Deletion (LWD) is the simplest method of the three. LWD removes any observation with at least one missing value, effectively reducing the dataset to complete cases. For our study, LWD will serve as the “control” algorithm. 
- Expectation Maximization is a computationally powerful method that iteratively 'guesses' on the missing data while maximizing the likelihood estimation of that guess, based on the observed data. EM alternates between estimating the expected values of the missing data given current parameter estimates (Expectation-step) and then updating the parameter estimates to maximize the likelihood of the complete data (Maximization-step). 
- Conditional Sampling (CS) involves generating multiple imputations for missing multivariate data, each conditioned on the observed mean and covariance structure. CS assumes the data follows a MVN distribution and imputes missing values by drawing samples from this distribution. Unlike EM, which iteratively estimates a single set of parameter estimates, CS generates multiple datasets through conditional stochastic sampling such that

For this repository, the sim.R and sim_functions.R carry the simulation study settings (dataset size, correlation parameter, missing data type, and imposed pattern). The graphics.R create the plots found in the images folder (except for the tables, as they are created using LaTeX). 

For a detailed description of the project, [view my blog post](https://jxnpass.github.io/2025/01/07/MissingValueImputation.html).

