#' Node-splitting for multivariate network meta-analysis
#' 
#' @description This function performs local inconsistency checks using the node-splitting method.
#'
#' @param x An object of class \code{\link{mvnma}}.
#' @param method.direct A character string indicating how direct estimates
#'   are obtained, either \code{"pairwise"} (outcome-specific pairwise
#'   meta-analysis) or \code{"multivariate"} (re-fit the multivariate network
#'   meta-analysis model to the direct evidence). In both cases the
#'   between-study heterogeneity is fixed at the value estimated by the
#'   multivariate model.Can be abbreviated.
#' @param tol.direct A numeric defining the maximum deviation of the direct
#'   evidence proportion from 0 or 1 to classify a comparison as providing
#'   only indirect or direct evidence, respectively. No indirect estimate is
#'   reported for such comparisons.
#' @param seed An optional numeric value used to set the random number
#'   generator before fitting the node-splitting models, in order to obtain
#'   reproducible results. Only relevant for \code{method.direct =
#'   "multivariate"}, where direct estimates are obtained from an MCMC model
#'   fit. The state of the random number generator is restored on exit.
#' @param quiet A logical indicating whether to print information on the
#'   progress of the JAGS model fitting. Only relevant for
#'   \code{method.direct = "multivariate"}, where direct estimates are
#'   obtained from an MCMC model fit.
#' @param digits Minimal number of significant digits for treatment
#'   estimates and confidence intervals, see \code{print.default}.
#' @param digits.se Minimal number of significant digits for standard
#'   errors, see \code{print.default}.
#' @param digits.pval Minimal number of significant digits for p-values of
#'   the test for disagreement between direct and indirect evidence, see
#'   \code{print.default}.
#' @param digits.stat Minimal number of significant digits for z-values of
#'   the test for disagreement between direct and indirect evidence, see
#'   \code{print.default}.
#' @param digits.prop Minimal number of significant digits for direct
#'   evidence proportions, see \code{print.default}.
#' @param print.se A logical indicating whether standard errors should be
#'   printed.
#' @param ci A logical indicating whether confidence intervals should be
#'   printed.
#' @param backtransf A logical indicating whether results should be back
#'   transformed for printing. For example, if \code{backtransf = TRUE},
#'   results for \code{sm = "OR"} are printed as odds ratios rather than log
#'   odds ratios, and the difference between the direct and indirect
#'   estimate is printed as a ratio of ratios. Estimates are always stored on
#'   the original scale in the \code{nodesplit} object; only the printed
#'   output is affected. Standard errors and z-values are never back
#'   transformed.
#' @param k A logical indicating whether the number of studies providing
#'   direct evidence should be printed.
#' @param prop A logical indicating whether the direct evidence proportion
#'   should be printed.
#' @param overall A logical indicating whether estimates from the
#'   multivariate network meta-analysis should be printed.
#' @param direct A logical indicating whether estimates derived from direct
#'   evidence should be printed.
#' @param indirect A logical indicating whether estimates derived from
#'   indirect evidence should be printed.
#' @param diff A logical indicating whether the difference between direct
#'   and indirect treatment estimates should be printed.
#' @param z A logical indicating whether z-values of the test for
#'   disagreement between direct and indirect evidence should be printed.
#' @param test A logical indicating whether p-values of the test for
#'   disagreement between direct and indirect evidence should be printed.
#' @param scientific.pval A logical specifying whether p-values should be
#'   printed in scientific notation, e.g., 1.2345e-01 instead of 0.12345.
#' @param text.NA A character string specifying text printed for missing
#'   values.
#' @param legend A logical indicating whether a legend should be printed.
#' @param \dots Additional arguments (ignored).
#' 
#' @details
#' Node-splitting is applied to every treatment comparison that is informed
#' by both direct and indirect evidence. A comparison qualifies if at least
#' one study reports it for the outcome and if the two treatments remain
#' connected in the network once that comparison is removed. Only the
#' comparison itself is removed, not the studies providing it, so that the
#' remaining arms of a multi-arm study still contribute to the indirect
#' estimate (Dias et al., 2010). The set of comparisons that can be split
#' may differ between outcomes, as an outcome is typically not reported by
#' all studies in the network.
#'
#' For each qualifying comparison, the direct estimate is obtained from the
#' studies directly comparing the two treatments, either by an
#' outcome-specific pairwise meta-analysis or by re-fitting the multivariate
#' model to the direct evidence (see argument \code{method.direct}). The
#' indirect estimate is then derived by back-calculation from the
#' multivariate network meta-analysis estimate and the direct estimate,
#' using the direct evidence proportion, that is, the ratio of the network
#' variance to the direct variance.
#'
#' Back-calculation divides both the indirect estimate and its variance by
#' one minus the direct evidence proportion. Results therefore become
#' unstable as the proportion approaches 1, and estimates for such
#' comparisons should be interpreted with caution. No indirect estimate is
#' reported if the proportion is not between \code{tol.direct} and 1 -
#' \code{tol.direct}; this includes proportions above 1, which can occur
#' when the network estimate is less precise than the direct estimate.
#'
#' Estimates are stored on the original scale and back transformed only for
#' printing (see argument \code{backtransf}). For ratio measures, the
#' difference between the direct and indirect estimate is a ratio of ratios
#' and is printed as \code{RoR}. Standard errors, z-values and p-values
#' always refer to the original scale.
#' 
#' @return
#' An object of class \code{nodesplit}; a list with one data frame per
#' outcome, with one row per treatment comparison and the following columns:
#' \item{comparison}{Treatment comparison.}
#' \item{k}{Number of studies providing direct evidence.}
#' \item{prop}{Direct evidence proportion. Reported even when it falls
#'   outside the range for which an indirect estimate is calculated.}
#' \item{mvnma.TE, mvnma.seTE, mvnma.lb, mvnma.ub}{Estimated treatment
#'   effect, standard error and confidence limits in the multivariate
#'   network meta-analysis.}
#' \item{direct.TE, direct.seTE, direct.lb, direct.ub}{Estimated treatment
#'   effect, standard error and confidence limits derived from direct
#'   evidence.}
#' \item{indirect.TE, indirect.seTE, indirect.lb, indirect.ub}{Estimated
#'   treatment effect, standard error and confidence limits derived from
#'   indirect evidence.}
#' \item{diff, se.diff, diff.lb, diff.ub}{Difference between direct and
#'   indirect estimates, its standard error and confidence limits.}
#' \item{z, p.val}{z-value and p-value of the test for disagreement between
#'   direct and indirect evidence.}
#' \item{sign}{A logical indicating whether the disagreement is statistically
#'   significant.}
#'
#' All estimates are on the original scale, that is, not back transformed
#' (see argument \code{backtransf}). Outcomes without any treatment
#' comparison contributing both direct and indirect evidence are set to
#' \code{NULL}, and \code{NULL} is returned with a warning if no such
#' comparison exists for any outcome.
#' 
#' @references
#' Dias S, Welton NJ, Caldwell DM, Ades AE (2010):
#' Checking consistency in mixed treatment comparison meta-analysis.
#' \emph{Statistics in Medicine},
#' \bold{29}, 932--44
#' 
#' König J, Krahn U, Binder H (2013):
#' Visualizing the flow of evidence in network meta-analysis and characterizing mixed treatment comparisons
#' \emph{Statistics in Medicine},
#' \bold{32}, 5414--29
#' 
#' Rücker G, Nikolakopoulou A, Papakonstantinou T, Salanti G, Riley RD,
#' Schwarzer G (2020):
#' The statistical importance of a study for a network meta-analysis estimate.
#' \emph{BMC Medical Research Methodology},
#' \bold{20}, 190
#' 
#' @examples
#' # Locate file "mvnma_examples.rda" with mvnma() results
#' .fname <- system.file("extdata/mvnma_examples.rda", package = "mvnma")
#' load(.fname)
#' 
#' # Local checks for inconsistency 
#' print(nodesplit(mvnma12),backtransf = FALSE)
#'  
#' @export nodesplit    

nodesplit <- function(x,
                      method.direct = c("pairwise", "multivariate"),
                      tol.direct = 5e-04,
                      seed = NULL,quiet = TRUE, 
                      ...){
  
  if (!inherits(x, "mvnma")) {
    stop("x should be an object of class mvnma")
  }
  #
  chklogical(quiet)
  chknumeric(tol.direct, min = 0, max = 1, length = 1)
  #
  method.direct <- match.arg(method.direct)
  #
  if (quiet) {
    oldopts <- options(jags.pb = "none")
    on.exit(options(oldopts), add = TRUE)
  }
  
  if (!is.null(seed)) {
    if (!exists(".Random.seed", envir = .GlobalEnv))
      runif(1)
    old.seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(assign(".Random.seed", old.seed, envir = .GlobalEnv), add = TRUE)
    set.seed(seed)
  }
  
  pairs <- attr(x, "pair.objects")
  
  name.outcome <- attr(x, "outcomes")
  
  sm <- attr(x, "sm")
  
  split.outcome <- splittable.comparisons(pairs, "list")
  
  names(split.outcome) <- name.outcome
  
  split.any <- splittable.comparisons(pairs, "any")
  
  if (nrow(split.any) == 0) {
    warning("The node-splitting cannot be applied as there are no nodes ",
            "contributing both direct and indirect evidence across all ",
            "outcomes.", call. = FALSE)
    return(invisible(NULL))
  }
  
  split.any$comp <- paste0(split.any$treat1, ":", split.any$treat2)
  
  r <- vector("list", nrow(split.any))
  
  for (i in seq_len(nrow(split.any))) {
    
    if (quiet) {
      # capture.output() evaluates in the calling frame, so the assignment
      # to r[[i]] persists
      invisible(capture.output(
        suppressMessages(suppressWarnings(
          r[[i]] <- pair.nodesplit(x,
                                   treat1 = split.any$treat1[i],
                                   treat2 = split.any$treat2[i],
                                   method.direct = method.direct,
                                   tol.direct = tol.direct)
                                    ## ins
        ))
      ))
    }
    else {
      r[[i]] <- pair.nodesplit(x,
                               treat1 = split.any$treat1[i],
                               treat2 = split.any$treat2[i],
                               method.direct = method.direct,
                               tol.direct = tol.direct)
                               ## ins
    }
  }
  
  comp <- paste0(split.any$treat1, ":", split.any$treat2)
  names(r) <- comp
  
  res <- lapply(name.outcome, function(o) {
    d <- do.call(rbind, lapply(r, function(z) z[[o]]))
    d <- cbind(comparison = comp, d)
    rownames(d) <- NULL
    d
  })
  names(res) <- name.outcome
  
  # keep only splittable comparisons for each outcome
  for (i in seq_along(res)) {
    
    if (nrow(split.outcome[[i]]) == 0) {
      res[i] <- list(NULL)
      next
    }
    
    split.outcome[[i]]$comparison <-
      paste0(split.outcome[[i]]$treat1, ":", split.outcome[[i]]$treat2)
    
    res[[i]] <- res[[i]] %>%
      filter(comparison %in% split.outcome[[i]]$comparison)
    
  }
  
  attr(res, "level") <- attr(x, "level")
  attr(res, "sm") <- sm
  class(res) <- "nodesplit"
  
  res
  
}

#' @rdname nodesplit
#' @method print nodesplit
#' @export

print.nodesplit <- function(x,
                            backtransf = gs("backtransf"),
                            digits = gs("digits"),
                            digits.se = gs("digits.se"),
                            digits.pval = gs("digits.pval"),
                            digits.stat = gs("digits.stat"),
                            digits.prop = max(gs("digits.pval") - 2, 2),
                            print.se = FALSE,
                            ci = TRUE,
                            k = TRUE,
                            prop = TRUE,
                            overall = TRUE,
                            direct = TRUE,
                            indirect = TRUE,
                            diff = TRUE,
                            z = TRUE,
                            test = TRUE,
                            scientific.pval = gs("scientific.pval"),
                            text.NA = gs("lab.NA"),
                            legend = TRUE,
                            ...) {
  
  chkclass(x, "nodesplit")
  #
  chknumeric(digits, min = 0, length = 1)
  chknumeric(digits.se, min = 0, length = 1)
  chknumeric(digits.pval, min = 1, length = 1)
  chknumeric(digits.stat, min = 0, length = 1)
  chknumeric(digits.prop, min = 0, length = 1)
  #
  chklogical(backtransf)
  chklogical(print.se)
  chklogical(ci)
  chklogical(k)
  chklogical(prop)
  chklogical(diff)
  chklogical(z)
  chklogical(overall)
  chklogical(direct)
  chklogical(indirect)
  chklogical(test)
  chklogical(scientific.pval)
  chklogical(legend)
  #
  level <- attr(x, "level")
  sm <- attr(x, "sm")
  #
  if (is.null(level))
    level <- 0.95
  #
  ci.lab <- paste0(round(100 * level, 1), "%-CI")
  
  # all results in a single row - be on the safe side
  oldopts <- options(width = 200)
  on.exit(options(oldopts))
  
  nam <- names(x)
  
  for (i in seq_along(nam)) {
    
    cat(paste0(if (i > 1) "\n" else "", "Outcome: ", nam[i], "\n\n"))
    
    dat.i <- x[[i]]
    
    bt.i <- backtransf && !is.null(sm) && !is.na(sm[i])
    #
    tobt <- function(z) if (bt.i) backtransf(z, sm[i]) else z
    #
    # ratio measures: the difference is a ratio of ratios
    rel.i <- !is.null(sm) && !is.na(sm[i]) &&
      (is_relative_effect(sm[i]) | sm[i] == "VE")
    
    if (is.null(dat.i) || nrow(dat.i) == 0) {
      cat("No comparison with both direct and indirect evidence.\n")
      next
    }
    
    # build column by column so that repeated CI columns can share a label
    out <- list(comparison = dat.i$comparison)
    names.out <- "comparison"
    #
    if (k) {
      out$k <- formatN(dat.i$k, digits = 0, text.NA = text.NA)
      names.out <- c(names.out, "k")
    }
    
    if (prop) {
      out$prop <- formatPT(dat.i$prop, digits = digits.prop)
      out$prop[rmSpace(out$prop) == "--"] <- text.NA
      names.out <- c(names.out, "prop")
    }
    if (overall) {
      out$TE.nma <- formatN(tobt(dat.i$mvnma.TE), digits = digits,
                            text.NA = text.NA)
      names.out <- c(names.out, "mvnma")
      #
      if (print.se) {
        out$mvnma.seTE <- formatN(dat.i$mvnma.seTE, digits = digits.se,
                                  text.NA = text.NA)
        names.out <- c(names.out, "seTE")
      }
      #
      if (ci) {
        out$mvnma.ci <- formatCI(formatN(tobt(dat.i$mvnma.lb), digits = digits),
                                 formatN(tobt(dat.i$mvnma.ub), digits = digits))
        out$mvnma.ci[is.na(out$mvnma.ci)] <- text.NA
        names.out <- c(names.out, ci.lab)
      }
    }
    
    if (direct) {
      out$TE.dir <- formatN(tobt(dat.i$direct.TE), digits = digits,
                            text.NA = text.NA)
      names.out <- c(names.out, "direct")
      #
      if (print.se) {
        out$se.dir <- formatN(dat.i$direct.seTE, digits = digits.se,
                              text.NA = text.NA)
        names.out <- c(names.out, "seTE")
      }
      #
      if (ci) {
        out$ci.dir <- formatCI(formatN(tobt(dat.i$direct.lb), digits = digits),
                               formatN(tobt(dat.i$direct.ub), digits = digits))
        out$ci.dir[is.na(out$ci.dir)] <- text.NA
        names.out <- c(names.out, ci.lab)
      }
    }
    
    if (indirect) {
      out$TE.ind <- formatN(tobt(dat.i$indirect.TE), digits = digits,
                            text.NA = text.NA)
      names.out <- c(names.out, "indir.")
      #
      if (print.se) {
        out$se.ind <- formatN(dat.i$indirect.seTE, digits = digits.se,
                              text.NA = text.NA)
        names.out <- c(names.out, "seTE")
      }
      #
      if (ci) {
        out$ci.ind <- formatCI(formatN(tobt(dat.i$indirect.lb), digits = digits),
                               formatN(tobt(dat.i$indirect.ub), digits = digits))
        out$ci.ind[is.na(out$ci.ind)] <- text.NA
        names.out <- c(names.out, ci.lab)
      }
    }
    if (diff) {
      #
      diff.lab <- if (backtransf && rel.i) "RoR" else "Diff"
      
      out$diff <- formatN(tobt(dat.i$diff), digits = digits, text.NA = text.NA)
      names.out <- c(names.out, diff.lab)
      #
      if (ci) {
        out$ci.diff <- formatCI(formatN(tobt(dat.i$diff.lb), digits = digits),
                                formatN(tobt(dat.i$diff.ub), digits = digits))
        out$ci.diff[is.na(out$ci.diff)] <- text.NA
        names.out <- c(names.out, ci.lab)
      }
    }
    
    if (z) {
      out$z <- formatN(dat.i$z, digits = digits.stat, text.NA = text.NA)
      out$z[rmSpace(out$z) == "--"] <- text.NA
      names.out <- c(names.out, "z")
    }
    
    if (test) {
      out$p <- formatPT(dat.i$p.val, digits = digits.pval,
                        scientific = scientific.pval)
      out$p[rmSpace(out$p) == "--"] <- text.NA
      names.out <- c(names.out, "p-value")
    }
    
    out <- as.data.frame(out, stringsAsFactors = FALSE)
    names(out) <- names.out
    #
    out[is.na(out)] <- text.NA
    
    prmatrix(out, quote = FALSE, right = TRUE,
             rowlab = rep("", nrow(out)))
  }
  
  if (legend) {
    cat("\nLegend:\n")
    cat(" comparison - Treatment comparison\n")
    if (k)
      cat(" k          - Number of studies providing direct evidence\n")
    if (prop)
      cat(" prop       - Direct evidence proportion\n")
    if (overall)
      cat(" mvnma      - Estimated treatment effect in multivariate network meta-analysis\n")
    if (direct)
      cat(" direct     - Estimated treatment effect derived from direct evidence\n")
    if (indirect)
      cat(" indir.     - Estimated treatment effect derived from indirect evidence\n")
    if (diff) {
      labels.used <- character(0)
      #
      for (i in seq_along(nam)) {
        dat.i <- x[[i]]
        #
        if (is.null(dat.i) || nrow(dat.i) == 0)
          next
        rel.i <- !is.null(sm) && !is.na(sm[i]) &&
          (is_relative_effect(sm[i]) | sm[i] == "VE")
        #
        labels.used <- union(labels.used,
                             if (backtransf && rel.i) "RoR" else "Diff")
      }
      #
      if ("Diff" %in% labels.used)
        cat(" Diff       - Difference between direct and indirect treatment estimates\n")
      if ("RoR" %in% labels.used)
        cat(" RoR        - Ratio of ratios between direct and indirect treatment estimates\n")
    }
    if (z)
      cat(" z          - z-value of test for disagreement (direct versus indirect)\n")
    if (print.se)
      cat(" seTE       - Standard error of treatment estimate\n")
    if (test)
      cat(" p-value    - p-value of test for disagreement",
          "(direct versus indirect)\n")
  }
  
  invisible(NULL)
}
