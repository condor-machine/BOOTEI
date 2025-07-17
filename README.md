# BOOTEI

Code for **BOOTEI: Bootstrap Ensemble Inference**, a method to improve permutation tests for categorical and ordinal data by smoothing over the discreteness of the permutation distribution via a bootstrap ensemble.

This repository contains the core C++ implementation and Monte Carlo simulation scripts used in the paper:


---

##  Contents

- `BOOTEI.cpp`: core engine (Rcpp) implementing BOOTEI-enhanced versions of:
  - Pearson χ² test (for independence)
  - Mann–Whitney (Wilcoxon rank-sum) test
  - Spearman correlation test  
  All tests support `alternative = "two.sided"`, `"greater"`, or `"less"`.  
  Set `B = 1` to disable bootstrap smoothing (standard permutation).

- `sim_chi2.R`: simulation study for χ² tests on categorical × categorical data.

- `sim_mannwhitney.R`: simulation study for Mann–Whitney on ordinal two-sample data.

- `sim_spearman.R`: simulation study for Spearman correlation on bivariate Likert-type data.

---

##  Requirements

- R (≥ 4.0)
- `Rcpp`, `tidyverse`, `progress`, `Cairo`, `MASS`

---


##  How to run

Each simulation script sources the `BOOTEI.cpp` file and compares standard permutation vs BOOTEI (with bootstrap smoothing) under varying sample sizes, significance levels, and strength of association.


## Example Usage (from R)

```r
# Load the BOOTEI engine (compile once per session)
Rcpp::sourceCpp("BOOTEI.cpp")

# Simulated categorical variables
set.seed(910)
x <- sample(c("A", "B"), 10, replace = TRUE)
y <- sample(c("yes", "no"), 10, replace = TRUE)

# Run BOOTEI χ² test with 100 bootstrap replicates and 1000 permutations
bootei(x, y, test = "chisq", B = 100, R = 1000)

# Run standard permutation χ² test (no bootstrap averaging)
bootei(x, y, test = "chisq", B = 1, R = 1000)
```




