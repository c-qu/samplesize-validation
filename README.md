# Simulation code 

The code for the simulations in the paper "Estimation of required sample size for external validation of risk models for binary outcomes" can be found above.

There is also an R package available. Details on how to install 'sampsizeval' and use it are given below.


<!-- README.md is generated from README.Rmd. Please edit that file -->

# sampsizeval: An R package to calculate the sample size for validation of risk models for binary outcomes

<!-- badges: start -->
[![Travis build status](https://travis-ci.com/mpavlou/sampsizeval.svg?branch=master)](https://travis-ci.com/mpavlou/sampsizeval)
[![R-CMD-check](https://github.com/mpavlou/sampsizeval/workflows/R-CMD-check/badge.svg)](https://github.com/mpavlou/sampsizeval/actions)
<!-- badges: end -->

The purpose of **sampsizeval** is to perform sample size calculations for the
validation of risk models for binary outcomes.

## Installation

The released version of sampsizeval can be installed from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("sampsizeval")
```

The development version can be installed from [GitHub](https://github.com/mpavlou) with:

``` r
# install.packages("devtools")
devtools::install_github("mpavlou/sampsizeval")
```
## Example

This is an example of a sample size calculation to validate a risk model for a binary outcome. The anticipated values of the outcome prevalence and the C-statistic are p=0.1 and C=0.75, respectively.

```{r example}
library(sampsizeval)
```

The target is to calculate the size of the validation data so as to estimate the C-statistic, the Calibration Slope and the Calibration in the Large with sufficient precision. In this example, the required precision is reflected by a SE of the estimated C-statistic of at most 0.025, and SE of the estimated Calibration Slope and Calibration in the Large of at most 0.1. 

```{r}
sampsizeval(p=0.1, c=0.75, se_c=0.025, se_cs =0.1, se_cl = 0.1)
```

The recommended sample size is 1536 observations.
<!-- Sample size required to achieve a SE of the Calibration Slope of at most 0.15. -->

<!-- Simple formula: -->

<!-- ```{r} -->
<!-- size_cs_ni(0.057, 0.7, 0.15^2) -->
<!-- ``` -->

<!-- Numerical integration: -->

<!-- ```{r} -->
<!-- size_cs_ni(0.057, 0.7, 0.15^2) -->
<!-- ``` -->

<!-- Sample size required to achieve a SE of the Calibration in the Large  of at most 0.15: -->

<!-- Simple formula: -->

<!-- ```{r} -->
<!-- size_cil(0.057, 0.7, 0.15^2) -->
<!-- ``` -->

<!-- Numerical integration: -->

<!-- ```{r} -->
<!-- size_cil_ni(0.057, 0.7, 0.15^2) -->
<!-- ``` -->

<!-- For a given precision for the estimated C-statistic, calibration slope and calibration in the large, the required sample size varies depending on the anticipated values of the  C-statistic and outcome prevalence. For example, for required precisions SE(C)=0.025, SE(CS)=0.15 and SE(CiL)=0.15, the sample size varies as follows: -->


<!-- ![Paper image](images/Figure_2_events.png) -->


