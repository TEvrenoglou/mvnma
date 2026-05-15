#' Forest plot for multivariate network meta-analysis results
#' 
#' @description
#' Draws a forest plot in the active graphics window (using grid graphics
#' system).
#' 
#' @param x An object of class \code{\link{mvnma}}.
#' @param backtransf A logical indicating whether results should be
#'   back transformed in forest plots. If \code{backtransf = TRUE},
#'   results for \code{sm = "OR"} are presented as odds ratios rather
#'   than log odds ratios, for example.
#' @param separate A logical indicating whether separate forest plots should be
#'   created for outcomes.
#' @param leftcols 	A character vector specifying (additional) columns to be
#'   printed on the left side of the forest plot or a logical value.
#' @param leftlabs 	A character vector specifying labels for (additional)
#'   columns on left side of the forest plot.
#' @param rightcols A character vector specifying (additional) columns to be 
#'   printed on the right side of the forest plot or a logical value.
#' @param rightlabs A character vector specifying labels for (additional)
#'   columns on right side of the forest plot. 
#' @param col.study The colour for treatment effect estimates and confidence
#'   interval limits. 
#' @param col.square The colour for squares reflecting the treatment effect
#'   estimates.
#' @param col.square.lines The colour for the outer lines of squares reflecting
#'   the treatment effect estimates.
#' @param col.subgroup The colour to print information on different outcomes.
#' @param squaresize The size of squares reflecting the treatment effect
#'   estimates (default: 0.7).
#' @param header.line A logical value indicating whether to print a header
#'   line (default: TRUE) or a character string ("both", "below", "").
#' @param text.outcome A text string printed before the name of the outcome.
#' @param \dots Additional arguments passed on to
#'   \code{\link[meta]{forest.meta}}.
#'
#' @return
#' A forest plot is shown in the active graphics window.
#' 
#' @keywords hplot
#' 
#' @examples
#' # Locate file "mvnma_examples.rda" with mvnma() results
#' .fname <- system.file("extdata/mvnma_examples.rda", package = "mvnma")
#' load(.fname)
#' 
#' # Print multivariate network meta-analysis results
#' mvnma12
#' 
#' # Generate a forest plot with the results
#' forest(mvnma12)
#' 
#' # Generate a separate forest plot for each outcome
#' forest(mvnma12, separate = TRUE)
#' 
#' @method forest mvnma 
#' @export

forest.mvnma <- function(x, backtransf = FALSE,
                         #
                         separate = n_distinct(attr(x, "sm")) > 1,
                         #
                         leftcols = "studlab", leftlabs,
                         rightcols = c("effect", "ci"), rightlabs,
                         #
                         col.study = "black",
                         col.square = "black",
                         col.square.lines = "black",
                         col.subgroup = "black",
                         #
                         squaresize = 0.7,
                         header.line = TRUE, 
                         #
                         text.outcome = "Outcome: ",
                         #
                         ...) {
  
  chkclass(x, "mvnma")
  #
  chklogical(backtransf)
  chklogical(separate)
  chkchar(text.outcome, length = 1)
  #
  method.model <- attr(x, "method.model")
  reference.group <- attr(x, "reference.group")
  n.domain <- attr(x,"n.domain")
  sm <- attr(x, "sm")
  #
  x <- x[names(x) != "cor"]
  #
  if (method.model == "DM") {
    if (is.null(n.domain)) {
      x <- x[names(x) != "sigma"]
    }
    else{
      x <- x[!(names(x) %in% c("sigma1", "sigma2"))]
    }
  }
  #
  n.out <- length(x)
  #
  if (missing(leftlabs)) {
    if (separate)
      leftlabs <- "Treatment"
    else
      leftlabs <- paste0("Comparison with '", reference.group, "'")
  }
  #
  if (missing(rightlabs))
    rightlabs <- rep(NA, length(rightcols))
  
  # Get rid of warning "no visible binding for global variable"
  treat <- mean <- sd <- lower <- upper <- studlab <- NULL
  
  # Get estimates for each outcome
  #
  ests <- vector("list")
  #
  for (i in seq_len(n.out)) {
    ests[[i]] <- x[[i]]$basic_estimates
    ests[[i]]$treat <- row.names(ests[[i]])
    row.names(ests[[i]]) <- NULL
    #
    ests[[i]] %<>% select(treat, mean, sd, lower, upper)
    ests[[i]]$outcome <- attr(x, "names")[i]
  }
  
  # Create forest plot(s)
  #
  if (!separate) {
    dat <- bind_rows(ests)
    row.names(dat) <- NULL
    names(dat) <- c("studlab", "mean", "seTE", "lower", "upper", "outcome")
    #
    # Drop rows for reference group 
    #
    dat %<>% filter(studlab != reference.group)
    
    if (length(unique(sm)) == 1)
      sm <- unique(sm)
    else
      sm <- ""
    #
    m <- metagen(dat$mean, dat$seTE, sm = sm,
                 subgroup = paste0(text.outcome, dat$outcome),
                 backtransf = backtransf,
                 print.subgroup.name = FALSE,
                 studlab = dat$studlab,
                 common = FALSE, random = FALSE, hetstat = FALSE,
                 method.tau = "DL", method.tau.ci = "")
    #
    res <- forest(m,
                  header.line = header.line,
                  col.subgroup = "black",
                  leftcols = leftcols, 
                  leftlabs = leftlabs,
                  rightcols = rightcols,
                  rightlabs = rightlabs,
                  weight.study = "same",
                  col.study = col.study,
                  col.square = col.square,
                  col.square.lines = col.square.lines,
                  squaresize = squaresize,
                  #
                  calcwidth.subgroup = TRUE,
                  ...)
  }
  else {
    res <- vector("list", n.out)
    #
    for (i in seq_len(n.out)) {
      dat.i <- ests[[i]]
      row.names(dat.i) <- NULL
      names(dat.i) <- c("studlab", "mean", "seTE", "lower", "upper", "outcome")
      #
      # Drop rows for reference group 
      #
      dat.i %<>% filter(studlab != reference.group)
      #
      m.i <- metagen(dat.i$mean, dat.i$seTE, sm = sm[i],
                     backtransf = backtransf,
                     studlab = dat.i$studlab,
                     common = FALSE, random = FALSE, hetstat = FALSE,
                     method.tau = "DL", method.tau.ci = "")
      m.i$.text.details.methods <-
        paste0("Details:\n- Comparison with '", reference.group, "'")
      #
      res[[i]] <- forest(m.i,
                         header.line = header.line,
                         col.subgroup = "black",
                         leftcols = leftcols, 
                         leftlabs = leftlabs,
                         rightcols = rightcols,
                         rightlabs = rightlabs,
                         weight.study = "same",
                         col.study = col.study,
                         col.square = col.square,
                         col.square.lines = col.square.lines,
                         squaresize = squaresize,
                         #
                         smlab = paste0(text.outcome, unique(dat.i$outcome)),
                         #
                         details = TRUE,
                         calcwidth.subgroup = TRUE,
                         ...)
    }
  }
  #
  invisible(res)
}
