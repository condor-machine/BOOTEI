# BOOTEI

Code for **BOOTEI: Bootstrap Ensemble Inference**, a method to improve permutation tests for categorical and ordinal data by smoothing over the discreteness of the permutation distribution via a bootstrap ensemble.

This repository contains the core C++ implementation and Monte Carlo simulation scripts used in the paper:

> *[BOOTEI: Smoothing permutation tests with bootstrap ensembles]* (full citation coming soon)

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



