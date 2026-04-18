# High-Dimensional Predictive Regression with Mixed Integration Orders and Endogeneity

This repository contains R code accompanying the paper:

> **"High-Dimensional Predictive Regression with Mixed Integration Orders and Endogeneity"**
> Alien Shrestha, Jiancheng Jiang, Daniel J. Henderson
> Department of Mathematics and Statistics, UNC Charlotte

---

## Overview

Standard predictive regression methods break down when predictors have uncertain orders of integration (stationary, near-integrated, or unit root) and when multiple predictors may be endogenous. This repository implements the **HD-WEE-SCAD** framework, which simultaneously addresses:

1. Mixed integration orders across predictors — no pre-testing for stationarity required
2. Endogeneity via a linear error projection approach
3. High-dimensional parameter spaces via SCAD and Adaptive LASSO regularization

Four method combinations are implemented and compared:
- **OLS + Adaptive LASSO**
- **OLS + SCAD**
- **WEE + Adaptive LASSO**
- **WEE + SCAD** *(recommended)*

---

## Repository Structure

```
├── simulation/
│   ├── ols_lasso_corrected.R       # OLS with Adaptive LASSO regularization
│   ├── ols_scad.R                  # OLS with SCAD regularization
│   ├── wee_lasso.R                 # WEE with Adaptive LASSO regularization
│   └── wee_scad.R                  # WEE with SCAD regularization (preferred method)
│
├── empirical/
│   ├── empirical_analysis.R        # Empirical application script
│   ├── period1_data.csv            # Period I: January 1970 -- April 1987
│   ├── period2_data.csv            # Period II: May 1988 -- December 2000
│   └── period3_data.csv            # Period III: January 2001 -- December 2021
│
└── README.md
```

---

## Simulation

The simulation scripts evaluate finite-sample performance in a **high-dimensional setting** with `p = 25` predictors (`2p + 1 = 51` parameters after endogeneity correction). The data generating process includes predictors spanning the full spectrum of integration orders:

- Stationary: `ρ = 0.4`, `ρ = 0.8`
- Unit root: `ρ = 1.0`
- Near-integrated: `ρ = 0.99`

Selective endogeneity is introduced for two of the four signal predictors. Performance is assessed via MSE, false positive rate, false negative rate, and bootstrap coverage probabilities (Universal/reverse-percentile CI method, `B = 500` bootstrap replications).

### Running the simulations

Each script is self-contained. To run, for example, the WEE-SCAD simulation:

```r
source("simulation/wee_scad.R")
```

Results are saved as output tables reporting point estimates, standard errors, MSE, variable selection counts, and coverage probabilities across 1,000 Monte Carlo replications.

### Key simulation findings

| Method | Unit root coverage | β FP | β FN |
|---|---|---|---|
| OLS-LASSO | ~0.91 | 0 | 0 |
| OLS-SCAD | ~0.91 | 0 | 0 |
| WEE-LASSO | ~0.94 | 37 | 0 |
| **WEE-SCAD** | **~0.94** | **0** | **0** |

Only WEE-SCAD achieves near-nominal coverage at the unit root boundary **and** perfect variable selection simultaneously.

---

## Empirical Application

The empirical script applies the WEE-SCAD method to U.S. stock return predictability using the [Welch-Goyal (2008)](https://academic.oup.com/rfs/article/21/4/1455/1567518) dataset, updated through 2021. Monthly excess S&P 500 returns are regressed on nine financial and macroeconomic predictors jointly across three economic periods.

### Data

Three period-specific data files are included:

| File | Period | Description |
|---|---|---|
| `period1_data.csv` | Jan 1970 -- Apr 1987 | Pre-crash era, high inflation, Volcker Fed |
| `period2_data.csv` | May 1988 -- Dec 2000 | Post-crash recovery, tech boom |
| `period3_data.csv` | Jan 2001 -- Dec 2021 | Post dot-com, GFC, COVID-19 |

### Running the empirical analysis

The script `empirical_analysis.R` is designed to run on one period at a time. To switch periods, update the data input path at the top of the script:

```r
# Set this to the desired period data file
data <- read.csv("empirical/period1_data.csv")  # Change to period2 or period3
```

Output includes point estimates, 95% Universal bootstrap confidence intervals, bootstrap zero-fractions, and endogeneity correction term estimates for all nine predictors.

---

## Dependencies

The following R packages are required:

```r
install.packages(c("glmnet", "ncvreg", "MASS"))
```

---

## Method Summary

The WEE estimator applies adaptive observation weights that downweight non-stationary predictors, restoring valid inference under mixed integration orders without requiring knowledge of which predictors are integrated. SCAD regularization is applied within the WEE framework via coordinate descent, achieving the oracle property: asymptotically zero false positives, zero false negatives, and nearly unbiased estimation for all nonzero coefficients. Bootstrap confidence intervals are constructed using the Universal (reverse-percentile) method with `N(1,1)` random weights, which remains valid when penalization induces non-standard bootstrap distributions.

For full methodological details see the accompanying paper.

---

## Citation

If you use this code, please cite:

> Shrestha, A., Jiang, J., & Henderson, D.J. (2024). *High-Dimensional Predictive Regression with Mixed Integration Orders and Endogeneity.* Working paper, UNC Charlotte.

---

## Contact

Alien Shrestha — ashrest8@charlotte.edu
