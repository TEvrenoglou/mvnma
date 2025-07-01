# mvnma: Multivariate Network Meta-Analysis

Official Git repository of R package **mvnma**

[![License: GPL (>=2)](https://img.shields.io/badge/license-GPL-blue)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)
[![CRAN Version](https://www.r-pkg.org/badges/version/mvnma)](https://cran.r-project.org/package=mvnma)
[![GitHub develop](https://img.shields.io/badge/develop-0.1--0-purple)](https://img.shields.io/badge/develop-0.1--0-purple)
[![Monthly Downloads](https://cranlogs.r-pkg.org/badges/mvnma)](https://cranlogs.r-pkg.org/badges/mvnma)
[![Total Downloads](https://cranlogs.r-pkg.org/badges/grand-total/mvnma)](https://cranlogs.r-pkg.org/badges/grand-total/mvnma)


## Authors

[Theodoros Evrenoglou](https://orcid.org/0000-0003-3336-8058),
[Guido Schwarzer](https://orcid.org/0000-0001-6214-9087)


## Description

**mvnma** is an R package that provides R functions for Bayesian multivariate network meta-analysis (mvNMA). The mvNMA model supported by this package refers to the single correlation coefficient model, interpreted as an amalgam of within- and across-outcome correlations. In this way, the model does not depend on the extraction of within-study outcome correlations, which are seldom reported at the study level. The treatment effect estimates and confidence intervals can be summarized both in terms of per-outcome treatment hierarchies and in terms of an across-outcomes benefit-risk assessment. The former is possible using ranking methods such as SUCRA, probability of best value, and median (or mean) ranks, each accompanied by a credible interval. A benefit-risk assessment is possible through the VIKOR method. This approach, originally proposed in the field of multi-criteria decision analysis, uses a deterministic algorithm to provide an amalgamated treatment hierarchy across outcomes and explicitly identify the set of treatments that offer the best compromise between benefits and harms across all outcomes. Since the output of the method is related to Markov Chain Monte Carlo (MCMC), convergence can be checked using a series of options, including trace plots, density plots, and the R-hat statistic. Finally, this R package offers the option to visualize the results of the mvNMA model through forest plots, which display the treatment effect estimates, scatter plots, which show the per-outcome rankings for any pair of outcomes, and Hasse diagrams, which visualize the partial order of the treatments across all outcomes.


## Installation

<!--
### Current stable [![CRAN Version](https://www.r-pkg.org/badges/version/mvnma)](https://cran.r-project.org/package=mvnma) release:
```r
install.packages("mvnma")
```
-->

### Current [![GitHub develop](https://img.shields.io/badge/develop-0.1--0-purple)](https://img.shields.io/badge/develop-0.1--0-purple) release on GitHub:

Installation using R package
[**remotes**](https://cran.r-project.org/package=remotes):
```r
install.packages("remotes")
remotes::install_github("TEvrenoglou/mvnma", ref = "develop")
```

### Bug Reports:

```r
bug.report(package = "mvnma")
```

The bug.report function is not supported in RStudio. Please send an email to Theodoros Evrenoglou <<theodoros.evrenoglou@uniklinik-freiburg.de>> if you use RStudio.

You can also report bugs on GitHub under [Issues](https://github.com/TEvrenoglou/mvnma/issues/).


<!--
### Reference

[Evrenoglou T, Nikolakopoulou A, Schwarzer G, Rücker G, Chaimani A (2024): Producing treatment hierarchies in network meta-analysis using probabilistic models and treatment-choice criteria. Preprint on ArXiv.](https://doi.org/10.48550/arXiv.2406.10612)
-->
