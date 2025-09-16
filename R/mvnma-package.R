#' mvnma: Brief overview of methods
#'
#' @description
#' R package \bold{mvnma} provides R functions for Bayesian multivariate
#' network meta-analysis (mvNMA). The mvNMA model supported by this package 
#' refers to the single correlation coefficient model, interpreted as an amalgam 
#' of within- and across-outcome correlations. In this way, the model does not 
#' depend on the extraction of within-study outcome correlations, which are seldom 
#' reported at the study level. The treatment effect estimates and confidence intervals 
#' can be summarized both in terms of per-outcome treatment hierarchies and in terms of 
#' an across-outcomes benefit-risk assessment. The former is possible using ranking methods 
#' such as SUCRA, probability of best value, and median (or mean) ranks, each accompanied by a 
#' credible interval. A benefit-risk assessment is possible through the VIKOR method. 
#' This approach, originally proposed in the field of multi-criteria decision analysis, uses 
#' a deterministic algorithm to provide an amalgamated treatment hierarchy across outcomes and 
#' explicitly identify the set of treatments that offer the best compromise between benefits and 
#' harms across all outcomes. Since the output of the method is related to Markov Chain Monte Carlo (MCMC), 
#' convergence can be checked using a series of options, including trace plots, density plots, and the R-hat 
#' statistic. Finally, this R package offers the option to visualize the results of the mvNMA model through 
#' forest plots, which display the treatment effect estimates, scatter plots, which show the per-outcome 
#' rankings for any pair of outcomes, and Hasse diagrams, which visualize the partial order of the treatments 
#' across all outcomes.
#'
#' @details
#' The R package \bold{mvnma} provides the following functions:
#' \itemize{
#' \item Function \code{\link{mvnma}} to perform a Bayesian multivariate
#'   network meta-analysis.
#' \item Function \code{\link{mvrank}} to get outcome-specific treatment
#'  rankings.
#' \item Function \code{\link{vikor}} to rank treatments across all outcomes
#'   using the VIKOR multi-criteria decision analysis method. Additionally,
#'   the function evaluates the concrete conditions defined by the VIKOR method
#'   and identifies the set of treatments that offer the best compromise
#'   between benefits and harms across all outcomes
#' \item Function \code{\link{forest.mvnma}} to visualize the results of the
#'   mvNMA model in terms of treatment effect estimates
#' \item Function \code{\link{plot.mvrank}} to visualize per outcome ranking
#'   results for any pair of outcomes
#' \item Function \code{\link{hasse.mvrank}} to visualize the partial order of
#'   the treatment across all outcomes
#' \item Function \code{\link{heatplot.mvrank}} to visualize in a heatplot
#'   the results in terms of outcome specific rankings
#' \item Function \code{\link{linechart}} to visualize the results of the three metrics calculated by the VIKOR method
#' \item Function \code{\link{as.mcmc.mvnma}} an auxiliary function to extract
#'   an MCMC object. This makes any \bold{mvnnma} object compatible with the
#'   convergence checks performed by the R package \bold{coda}.
#' }
#' 
#' Type \code{help(package = "mvnma")} for a listing of R functions
#' available in \bold{mvnma}.
#'
#' Type \code{citation("mvnma")} on how to cite \bold{mvnma}
#' in publications.
#'
#' To report problems and bugs, please send an email to Theodoros
#' Evrenoglou <theodoros.evrenoglou@uniklinik-freiburg.de>.
#'
#' The development version of \bold{mvnma} is available on GitHub
#' \url{https://github.com/TEvrenoglou/mvnma}.
#'
#' @name mvnma-package
#'
#' @author Theodoros Evrenoglou <theodoros.evrenoglou@@uniklinik-freiburg.de>,
#'   Guido Schwarzer <guido.schwarzer@@uniklinik-freiburg.de>
#'
#' @keywords package
#'
#' @importFrom R2jags jags
#' @importFrom coda as.mcmc as.mcmc.list
#' @importFrom meta forest gs metagen pairwise
#' @importFrom netmeta hasse netposet rankogram heatplot
#' @importFrom matrixStats colSds
#' @importFrom dplyr %>% all_of any_of arrange bind_rows bind_cols desc
#'   distinct filter group_by mutate rename select pull
#' @importFrom magrittr %<>%
#' @importFrom rlist list.cbind list.rbind
#' @importFrom graphics text
#' @importFrom stats complete.cases quantile relevel
#' @importFrom utils combn packageVersion
#' @importFrom ggplot2 ggplot aes geom_tile geom_line geom_point geom_text scale_fill_gradient guides 
#'   guide_colourbar guide_legend labs xlab ylab ylim scale_y_discrete theme theme_void theme_minimal
#'   element_text element_blank
#' @importFrom forcats fct_rev
#' @export as.mcmc

"_PACKAGE"

NULL
