#' mvnma: Brief overview of methods
#'
#' @description
#' R package \bold{mvnma} provides R functions for Bayesian multivariate
#' network meta-analysis (mvNMA). The mvNMA model supported by this package
#' refers to the single correlation coefficient model, interpreted as an
#' amalgam of within- and across-outcome correlations (Efthimiou et al., 2015)
#' which is a generalisation of Riley et al. (2008).
#' In this way, the model does not depend on the extraction of within-study
#' outcome correlations, which are seldom reported at the study level.
#' 
#' DuMouchel priors assuming constant relative treatment effects across
#' outcomes and enabling information sharing can be used (DuMouchel & Harris,
#' 1983). This may improve precision but can introduce bias when outcomes from
#' different domains (e.g., efficacy and safety) are analyzed jointly.
#' 
#' The treatment effect estimates
#' and confidence intervals can be summarised both in terms of per-outcome
#' treatment hierarchies and in terms of an across-outcomes benefit-risk
#' assessment. The former is possible using ranking methods such as the surface
#' under the cumulative ranking curve (SUCRA) (Salanti et al., 2011), the
#' probability of best value, and median (or mean) ranks, each accompanied by a
#' credible interval.
#' 
#' A benefit-risk assessment is possible through the VišeKriterijumska
#' Optimizacija I Kompromisno Rešenje (VIKOR) method (Opricovic & Tzeng, 2004;
#' Opricovic, 2011). This approach, originally proposed in the field of
#' multi-criteria decision analysis, uses a deterministic algorithm to provide
#' an amalgamated treatment hierarchy across outcomes and explicitly identify
#' the set of treatments that offer the best compromise between benefits and
#' harms across all outcomes.
#' 
#' Since the output of the method relies on Markov Chain Monte Carlo (MCMC)
#' sampling, convergence can be checked using a series of diagnostics, including
#' trace plots, density plots, and the R-hat statistic. Finally, this package
#' offers the option to visualise the results of the mvNMA model through forest
#' plots, which display the treatment effect estimates; scatter plots, which
#' show the per-outcome rankings for any pair of outcomes; and Hasse diagrams
#' (Carlsen & Bruggemann, 2014), which visualise the partial order of treatments
#' across all outcomes (Rücker & Schwarzer, 2017), as well as alternative
#' methods to yield a between-outcomes hierarchy, such as the spie chart method
#' (Daly et al., 2020).
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
#'   between benefits and harms across all outcomes.
#' \item Function \code{\link{forest.mvnma}} to visualize the results of the
#'   mvNMA model in terms of treatment effect estimates.
#' \item Function \code{\link{plot.mvrank}} to visualize per outcome ranking
#'   results for any pair of outcomes.
#' \item Function \code{\link{hasse.mvrank}} to visualize the partial order of
#'   the treatment across all outcomes.
#' \item Function \code{\link{heatplot.mvrank}} to visualize in a heatplot
#'   the results in terms of outcome specific rankings.
#' \item Function \code{\link{linechart}} to visualize the results of the three
#'   metrics calculated by the VIKOR method.
#' \item Function \code{\link{as.mcmc.mvnma}} an auxiliary function to extract
#'   an MCMC object. This makes any \bold{mvnnma} object compatible with the
#'   convergence checks performed by the R package \bold{coda}.
#'  \item Function \code{\link{spiechart}} to calculate a between outcomes
#'    hierarchy using the spie charts method.
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
#' @references
#' Carlsen L, Bruggemann R (2014):
#' Partial order methodology: a valuable tool in chemometrics.
#' \emph{Journal of Chemometrics},
#' \bold{28}, 226--34
#' 
#' Daly CH, Mbuagbaw L, Thabane L, Straus SE, Hamid JS (2020):
#' Spie charts for quantifying treatment effectiveness and safety in multiple
#' outcome network meta-analysis: a proof-of-concept study.
#' \emph{BMC Med Res Methodol}, \bold{20}, 266
#' 
#' DuMouchel WH, Harris JE (1983):
#' Bayes methods for combining the results of cancer studies in humans and
#' other species.
#' \emph{Journal of the American Statistical Association},
#' \bold{78}, 293--308
#' 
#' Efthimiou O, Mavridis D, Riley RD, Cipriani A, Salanti G (2015):
#' Joint synthesis of multiple correlated outcomes in networks of interventions.
#' \emph{Biostatistics}, 
#' \bold{16}, 84--97
#'
#' Opricovic S, Tzeng GH (2004):
#' Compromise solution by MCDM methods: A comparative analysis of VIKOR and
#' TOPSIS.
#' \emph{European Journal of Operational Research},
#' \bold{156}, 445--55
#' 
#' Opricovic S (2011):
#' Fuzzy VIKOR with an application to water resources planning. 
#' \emph{Expert Systems with Applications}, 
#' \bold{38}, 12983--90
#' 
#' Riley RD, Thompson JR, Abrams KR (2008):
#' An alternative model for bivariate random-effects meta-analysis when the
#' within-study correlations are unknown.
#' \emph{Biostatistics},
#' \bold{9}, 172--86
#' 
#' Rücker G, Schwarzer G (2017):
#' Resolve conflicting rankings of outcomes in network meta-analysis:
#' Partial ordering of treatments.
#' \emph{Research Synthesis Methods},
#' \bold{8}, 526--36
#'                                             
#' Salanti G, Ades AE, Ioannidis JP (2011):
#' Graphical methods and numerical summaries for presenting results
#' from multiple-treatment meta-analysis: an overview and tutorial.
#' \emph{Journal of Clinical Epidemiology},
#' \bold{64}, 163--71
#' 
#' @keywords package
#'
#' @importFrom R2jags jags
#' @importFrom coda as.mcmc as.mcmc.list
#' @importFrom meta forest gs metagen pairwise
#' @importFrom netmeta hasse netposet rankogram heatplot
#' @importFrom matrixStats colSds
#' @importFrom dplyr %>% all_of any_of arrange bind_rows bind_cols desc distinct filter group_by mutate rename select pull n_distinct
#' @importFrom magrittr %<>%
#' @importFrom rlist list.cbind list.rbind
#' @importFrom graphics text
#' @importFrom stats complete.cases quantile relevel
#' @importFrom utils combn packageVersion
#' @importFrom ggplot2 ggplot aes geom_tile geom_line geom_point geom_text scale_fill_gradient guides guide_colourbar guide_legend labs xlab ylab ylim scale_y_discrete theme theme_void theme_minimal element_text element_blank
#' @importFrom forcats fct_rev
#' @export as.mcmc

"_PACKAGE"

NULL
