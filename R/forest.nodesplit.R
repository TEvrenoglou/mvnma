#' Forest plot for direct and indirect evidence in multivariate network
#' meta-analysis
#' 
#' @description
#' Forest plot showing the multivariate network meta-analysis estimate, the
#' direct estimate and the indirect estimate for each treatment comparison
#' that has been node-split. A separate forest plot is produced for each
#' outcome.
#' 
#' @param x An object of class \code{nodesplit}.
#' @param outcome An optional character vector or numeric vector selecting
#'   the outcome(s) to plot. By default, a forest plot is produced for every
#'   outcome with at least one comparison providing both direct and indirect
#'   evidence.
#' @param show A character string indicating which comparisons should be
#'   plotted: \code{"both"} (comparisons with both direct and indirect
#'   evidence), \code{"all"}, \code{"with.direct"} (comparisons with direct
#'   evidence), \code{"direct.only"}, or \code{"indirect.only"}.
#' @param subgroup A character string indicating which layout should be used
#'   in the forest plot: subgroups by treatment comparison
#'   (\code{"comparison"}) or by type of estimate (\code{"estimate"}).
#' @param overall A logical indicating whether estimates from the
#'   multivariate network meta-analysis should be plotted.
#' @param direct A logical indicating whether estimates derived from direct
#'   evidence should be plotted.
#' @param indirect A logical indicating whether estimates derived from
#'   indirect evidence should be plotted.
#' @param sortvar An optional vector used to sort the treatment comparisons.
#' @param subset An optional logical vector selecting a subset of the
#'   treatment comparisons.
#' @param text.overall A character string used in the forest plot to label
#'   the multivariate network meta-analysis estimate.
#' @param text.direct A character string used in the forest plot to label
#'   the direct estimate.
#' @param text.indirect A character string used in the forest plot to label
#'   the indirect estimate.
#' @param type.overall A character string specifying how to plot the
#'   multivariate network meta-analysis estimate.
#' @param type.direct A character string specifying how to plot the direct
#'   estimate.
#' @param type.indirect A character string specifying how to plot the
#'   indirect estimate.
#' @param col.square The colour for squares.
#' @param col.square.lines The colour for the outer lines of squares.
#' @param col.diamond The colour for diamonds.
#' @param col.diamond.lines The colour for the outer lines of diamonds.
#' @param equal.size A logical indicating whether all squares should be of
#'   equal size.
#' @param leftcols A character vector specifying (additional) columns to be
#'   printed on the left side of the forest plot.
#' @param leftlabs A character vector specifying labels for the columns on
#'   the left side of the forest plot.
#' @param rightcols A character vector specifying (additional) columns to be
#'   printed on the right side of the forest plot.
#' @param rightlabs A character vector specifying labels for the columns on
#'   the right side of the forest plot.
#' @param digits Minimal number of significant digits for treatment
#'   estimates and confidence intervals, see \code{print.default}.
#' @param digits.prop Minimal number of significant digits for direct
#'   evidence proportions, see \code{print.default}.
#' @param backtransf A logical indicating whether results should be back
#'   transformed in the forest plot. For example, if \code{backtransf =
#'   TRUE}, results for \code{sm = "OR"} are shown as odds ratios rather
#'   than log odds ratios.
#' @param lab.NA A character string to label missing values.
#' @param smlab A label for the summary measure. By default, the name of the
#'   outcome is used.
#' @param file An optional character string with the name of a file to save
#'   the forest plot(s), see \code{\link[meta]{forest.meta}}. If more than
#'   one outcome is plotted, the outcome name is appended to the file name,
#'   e.g., \code{"forest_Early_Response.pdf"}.
#' @param \dots Additional arguments passed on to \code{\link[meta]{forest.meta}}.
#' 
#' @details
#' A separate forest plot is produced for each selected outcome, since the
#' set of comparisons providing both direct and indirect evidence can differ
#' between outcomes. An outcome is skipped, with a warning, if it holds no
#' such comparison or if argument \code{show} selects none of them.
#'
#' Estimates are stored in the \code{nodesplit} object on the original
#' scale, so argument \code{backtransf} controls the scale used in the
#' forest plot.
#'
#' If argument \code{file} is provided, one file is written per outcome, and
#' the outcome name is appended to the file name when several outcomes are
#' plotted so that the files do not overwrite each other. To collect all
#' plots in a single multi-page file, open a graphics device before calling
#' this function instead of using argument \code{file}.
#'
#' @return
#' A list with the data sets used to produce the forest plot(s), returned
#' invisibly, with one element per selected outcome. Each element holds the
#' estimates plotted for that outcome, on the original scale. Skipped
#' outcomes are \code{NULL}.
#' @seealso \code{\link{nodesplit}}
#' 
#' @examples
#' \dontrun{
#' .fname <- system.file("extdata/mvnma_examples.rda", package = "mvnma")
#' load(.fname)
#' 
#' ns <- nodesplit(mvnma_all)
#' forest(ns)
#' }
#' 
#' @method forest nodesplit
#' @export

forest.nodesplit <- function(x,
                             outcome = NULL,
                             show = "both",
                             subgroup = "comparison",
                             overall = TRUE,
                             direct = TRUE,
                             indirect = TRUE,
                             sortvar = NULL,
                             subset = NULL,
                             text.overall = "mvNMA estimate",
                             text.direct = "Direct estimate",
                             text.indirect = "Indirect estimate",
                             type.overall,
                             type.direct,
                             type.indirect,
                             col.square = "gray",
                             col.square.lines = col.square,
                             col.diamond = "gray",
                             col.diamond.lines = "black",
                             equal.size = TRUE,
                             leftcols,
                             leftlabs,
                             rightcols = c("effect", "ci"),
                             rightlabs = NULL,
                             digits = gs("digits.forest"),
                             digits.prop = max(gs("digits.pval") - 2, 2),
                             backtransf = gs("backtransf"),
                             lab.NA = "",
                             smlab,
                             file = NULL,
                             ...) {
  
  chkclass(x, "nodesplit")
  #
  chklogical(overall)
  chklogical(direct)
  chklogical(indirect)
  chklogical(equal.size)
  #
  chknumeric(digits, min = 0, length = 1)
  chknumeric(digits.prop, min = 0, length = 1)
  #
  chklogical(backtransf)
  chkchar(lab.NA)
  if (!is.null(file))
    chkchar(file, length = 1)
  #
  chkchar(text.overall)
  chkchar(text.direct)
  chkchar(text.indirect)
  #
  show <- setchar(show, c("all", "both", "with.direct",
                          "direct.only", "indirect.only"))
  subgroup <- setchar(subgroup, c("comparison", "estimate"))
  
  if (!any(c(overall, direct, indirect)))
    stop("At least, one of the following estimates ",
         "must be included in forest plot:\n",
         "- mvNMA estimates (argument 'overall')\n",
         "- direct estimates (argument 'direct')\n",
         "- indirect estimates (argument 'indirect')",
         call. = FALSE)
  
  # type of plotting symbol
  #
  missing.type.overall <- missing(type.overall)
  #
  if (missing.type.overall)
    type.overall <- "diamond"
  else
    type.overall <- setchar(type.overall, c("diamond", "square"))
  #
  if (missing(type.direct))
    type.direct <- "square"
  else
    type.direct <- setchar(type.direct, c("diamond", "square"))
  #
  if (missing(type.indirect))
    type.indirect <- "square"
  else
    type.indirect <- setchar(type.indirect, c("diamond", "square"))
  #
  n.subgroup <- overall + direct + indirect
  #
  if (n.subgroup == 1 & overall & missing.type.overall)
    type.overall <- "square"
  
  if (overall & n.subgroup > 1) {
    if (text.overall == text.direct)
      stop("Text must be different for arguments 'text.overall' and ",
           "'text.direct'.", call. = FALSE)
    if (text.overall == text.indirect)
      stop("Text must be different for arguments 'text.overall' and ",
           "'text.indirect'.", call. = FALSE)
  }
  #
  if (direct & indirect & text.direct == text.indirect)
    stop("Text must be different for arguments 'text.direct' and ",
         "'text.indirect'.", call. = FALSE)
  
  # columns on the left hand side
  #
  missing.leftcols <- missing(leftcols)
  #
  if (missing.leftcols) {
    if (direct)
      leftcols <- c("studlab", "k", "prop")
    else
      leftcols <- "studlab"
  }
  #
  missing.leftlabs <- missing(leftlabs)
  #
  if (missing.leftlabs) {
    leftlabs <- rep(NA, length(leftcols))
    leftlabs[leftcols == "studlab"] <- "Comparison"
    leftlabs[leftcols == "k"] <- "Number of\nStudies"
    leftlabs[leftcols == "prop"] <- "Direct\nEvidence"
  }
  
  missing.smlab <- missing(smlab)
  
  sm.all <- attr(x, "sm")
  nam <- names(x)
  
  # select outcomes
  #
  if (is.null(outcome))
    sel.out <- seq_along(nam)
  else if (is.character(outcome)) {
    sel.out <- match(outcome, nam)
    #
    if (anyNA(sel.out))
      stop("Argument 'outcome' must match the outcome name(s): ",
           paste(paste0("'", nam, "'"), collapse = ", "), ".", call. = FALSE)
  }
  else {
    chknumeric(outcome, min = 1, max = length(nam))
    sel.out <- outcome
  }
  
  res <- vector("list", length(sel.out))
  names(res) <- nam[sel.out]
  
  for (j in seq_along(sel.out)) {
    
    i <- sel.out[j]
    #
    dat.i <- x[[i]]
    
    if (is.null(dat.i) || nrow(dat.i) == 0) {
      warning("No comparison with both direct and indirect evidence for ",
              "outcome '", nam[i], "'.", call. = FALSE)
      next
    }
    
    sm.i <- if (!is.null(sm.all) && !is.na(sm.all[i])) sm.all[i] else ""
    #
    # one data set per type of estimate
    #
    dat.overall <-
      data.frame(comparison = dat.i$comparison,
                 TE = dat.i$mvnma.TE,
                 seTE = dat.i$mvnma.seTE,
                 lower = dat.i$mvnma.lb,
                 upper = dat.i$mvnma.ub,
                 k = NA, prop = NA,
                 evidence = text.overall,
                 type.study = type.overall,
                 stringsAsFactors = FALSE)
    #
    dat.direct <-
      data.frame(comparison = dat.i$comparison,
                 TE = dat.i$direct.TE,
                 seTE = dat.i$direct.seTE,
                 lower = dat.i$direct.lb,
                 upper = dat.i$direct.ub,
                 k = dat.i$k,
                 prop = formatPT(dat.i$prop, digits = digits.prop),
                 evidence = text.direct,
                 type.study = type.direct,
                 stringsAsFactors = FALSE)
    #
    dat.indirect <-
      data.frame(comparison = dat.i$comparison,
                 TE = dat.i$indirect.TE,
                 seTE = dat.i$indirect.seTE,
                 lower = dat.i$indirect.lb,
                 upper = dat.i$indirect.ub,
                 k = NA, prop = NA,
                 evidence = text.indirect,
                 type.study = type.indirect,
                 stringsAsFactors = FALSE)
    
    # colours
    #
    dat.overall$col.estimate <-
      if (type.overall == "square") col.square else col.diamond
    dat.direct$col.estimate <-
      if (type.direct == "square") col.square else col.diamond
    dat.indirect$col.estimate <-
      if (type.indirect == "square") col.square else col.diamond
    #
    dat.overall$col.lines <-
      if (type.overall == "square") col.square.lines else col.diamond.lines
    dat.direct$col.lines <-
      if (type.direct == "square") col.square.lines else col.diamond.lines
    dat.indirect$col.lines <-
      if (type.indirect == "square") col.square.lines else col.diamond.lines
    
    # select comparisons
    #
    if (show == "all")
      sel <- rep_len(TRUE, nrow(dat.i))
    else if (show == "with.direct")
      sel <- !is.na(dat.i$direct.TE)
    else if (show == "both")
      sel <- !is.na(dat.i$direct.TE) & !is.na(dat.i$indirect.TE)
    else if (show == "direct.only")
      sel <- !is.na(dat.i$direct.TE) & is.na(dat.i$indirect.TE)
    else if (show == "indirect.only")
      sel <- is.na(dat.i$direct.TE) & !is.na(dat.i$indirect.TE)
    #
    if (!any(sel)) {
      warning("No comparison selected for outcome '", nam[i],
              "'. Consider using argument 'show = \"all\"'.", call. = FALSE)
      next
    }
    #
    dat.overall <- dat.overall[sel, , drop = FALSE]
    dat.direct <- dat.direct[sel, , drop = FALSE]
    dat.indirect <- dat.indirect[sel, , drop = FALSE]
    
    # sort and subset
    #
    if (!is.null(sortvar)) {
      #
      if (length(sortvar) != nrow(dat.i))
        stop("Argument 'sortvar' must be of length ", nrow(dat.i), ".",
             call. = FALSE)
      #
      sortvar.i <- sortvar[sel]
      #
      if (!is.numeric(sortvar.i))
        sortvar.i <- match(sortvar.i, dat.overall$comparison)
      #
      o <- order(sortvar.i)
      #
      dat.overall <- dat.overall[o, , drop = FALSE]
      dat.direct <- dat.direct[o, , drop = FALSE]
      dat.indirect <- dat.indirect[o, , drop = FALSE]
    }
    #
    if (!is.null(subset)) {
      #
      if (!is.logical(subset))
        stop("Argument 'subset' must be a logical vector.", call. = FALSE)
      #
      if (length(subset) != nrow(dat.i))
        stop("Argument 'subset' must be of length ", nrow(dat.i), ".",
             call. = FALSE)
      #
      subset.i <- subset[sel]
      #
      dat.overall <- dat.overall[subset.i, , drop = FALSE]
      dat.direct <- dat.direct[subset.i, , drop = FALSE]
      dat.indirect <- dat.indirect[subset.i, , drop = FALSE]
    }
    
    dat <- rbind(if (direct) dat.direct,
                 if (indirect) dat.indirect,
                 if (overall) dat.overall)
    #
    if (nrow(dat) == 0) {
      warning("No comparison(s) selected for outcome '", nam[i], "'.",
              call. = FALSE)
      next
    }
    
    smlab.i <- if (missing.smlab) nam[i] else smlab
    
    # one file per outcome; the outcome name is appended to the file name
    # if more than one outcome is plotted, as forest() opens and closes the
    # graphics device on every call and would otherwise overwrite the file
    #
    file.i <- NULL
    #
    if (!is.null(file)) {
      if (length(sel.out) > 1)
        file.i <-
          sub("(\\.[[:alnum:]]+)$",
              paste0("_", gsub("[^[:alnum:]._-]", "_", nam[i]), "\\1"),
              file)
      else
        file.i <- file
    }
    
    # meta-analysis object holding the estimates to plot
    #
    if (subgroup == "comparison")
      m <- suppressWarnings(
        metagen(dat$TE, dat$seTE, lower = dat$lower, upper = dat$upper,
                studlab = if (n.subgroup > 1) dat$evidence else dat$comparison,
                data = dat, sm = sm.i,
                common = FALSE, random = FALSE,
                method.tau = "DL", method.tau.ci = "",
                subgroup = if (n.subgroup > 1) dat$comparison else NULL,
                print.subgroup.name = FALSE))
    else
      m <- suppressWarnings(
        metagen(dat$TE, dat$seTE, lower = dat$lower, upper = dat$upper,
                studlab = dat$comparison,
                data = dat, sm = sm.i,
                common = FALSE, random = FALSE,
                method.tau = "DL", method.tau.ci = "",
                subgroup = if (n.subgroup > 1) dat$evidence else NULL,
                print.subgroup.name = FALSE))
    #
    if (overall & n.subgroup > 1) {
      m$w.common[dat$evidence == text.overall] <-
        max(m$w.common, na.rm = TRUE)
      m$w.random[dat$evidence == text.overall] <-
        max(m$w.random, na.rm = TRUE)
    }
    
    if (subgroup == "comparison")
      forest(m,
             digits = digits, common = FALSE, random = FALSE,
             hetstat = FALSE, test.subgroup = FALSE,
             leftcols = leftcols, leftlabs = leftlabs,
             rightcols = rightcols, rightlabs = rightlabs,
             lab.NA = lab.NA, smlab = smlab.i, backtransf = backtransf,
             file = file.i,
             type.study = dat$type.study,
             col.square = dat$col.estimate,
             col.square.lines = dat$col.lines,
             weight.study = if (equal.size) "same" else "common",
             ...)
    else
      forest(m,
             digits = digits, overall = FALSE, common = FALSE, random = FALSE,
             hetstat = FALSE, test.subgroup = FALSE,
             subgroup.hetstat = FALSE, prediction.subgroup = FALSE,
             leftcols = leftcols, leftlabs = leftlabs,
             rightcols = rightcols, rightlabs = rightlabs,
             lab.NA = lab.NA, smlab = smlab.i, backtransf = backtransf,
             file = file.i,
             type.study = dat$type.study,
             col.square = dat$col.estimate,
             col.square.lines = dat$col.lines,
             weight.study = if (equal.size) "same" else "common",
             ...)
    
    res[[j]] <- dat
  }
  
  invisible(res)
}

