# Regularized entropy‑balancing code

This repository contains R source code for the manuscript:
> Regularized entropy balancing methods for high dimensional covariates with an application to reservoir impacts

## Contents
- `code/01_simulation_and_summary.R`: Monte‑Carlo simulation studies
- `code/02_empirical_and_summary.R`: Empirical analysis for 897 Brazilian catchments (CAMELS‑BR dataset)

## Dependencies
Required R packages:
`gbm, glmnet, Matching, CBPS, sbw, quadprog, MASS`

Install packages example:
```r
install.packages(c("gbm","glmnet","Matching","CBPS","sbw","quadprog","MASS"))
