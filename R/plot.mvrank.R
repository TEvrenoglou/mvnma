#' Scatter plot to visualize the ranking of two outcomes in a multivariate network meta-analysis
#' 
#' @description
#' Draw a scatter plot in the active graphics window.
#' 
#' @param x An object of class \code{\link{mvrank}}.
#' @param which A mandatory numeric vector of length 2 specifying which
#'   outcomes should be plotted. For example, setting "outcome = c(2, 3)"
#'   implies that a scatter plot will be generated plotting the rankings
#'   of outcomes 2 and 3.
#' @param pos Position of treatment labels.
#' @param cex.point a numeric value specifying the size of the points
#'   (default: 1)
#' @param cex.label a numeric value specifying the size of the point labels in
#'   the plot (default: 0.7)
#' @param pch a vector of plotting characters or symbols (default: 19)
#' @param xlim the x limits of the plot
#' @param ylim the y limits of the plot
#' @param \dots Additional arguments for \code{\link{plot}} function.
#'
#' @return
#' A scatter plot is shown in the active graphics window.
#' 
#' @keywords hplot
#' 
#' @examples
#' # Locate file "mvnma_example.rda" with mvnma() results
#' .fname <- system.file("extdata/mvnma_examples.rda", package = "mvnma")
#' load(.fname)
#' 
#' # Rank treatments using SUCRAs (default)
#' ranks12 <- mvrank(mvnma12, small.values = c("undes", "undes"))
#' ranks12
#' 
#' # Visualize SUCRAs in a scatter plot with outcome 1 on the x-axis and
#' # outcome 2 on the y-axis
#' plot(ranks12)
#' 
#' # Visualize SUCRAs in a scatter plot with outcome 2 on the x-axis and
#' # outcome 1 on the y-axis
#' plot(ranks12, which = 2:1)
#' 
#' @method plot mvrank
#' @export

plot.mvrank <- function(x, which = 1:2, pos = 1,
                        cex.point = 1, cex.label = 0.7, pch = 19,
                        xlim = c(0, 1), ylim = c(0, 1),
                        ...) {
  
  chkclass(x, "mvrank")
  #
  n.outcome <- length(names(x))
  common_trts <- attr(x, "common_trts")
  #
  chknumeric(which, min = 1, max = n.outcome, length = 2)
  chknumeric(cex.point, min = 0, zero = TRUE)
  chknumeric(cex.label, min = 0, zero = TRUE)
  chknumeric(pch, min = 1, zero = TRUE)
  chknumeric(xlim, length = 2)
  chknumeric(ylim, length = 2)
  #
  first <- which[1]
  second <- which[2]
  
  # Get rid of warning "no visible binding for global variable"
  treat <- NULL
  
  outcomes <- names(x)[c(first, second)]
  #
  dat1 <- x[[first]]
  names(dat1)[1:2] <- c("treat", "rank1")
  dat1$out1 <- outcomes[1]
  #
  dat1 %<>% filter(treat %in% common_trts)
  
  dat2 <- x[[second]]
  names(dat2)[1:2] <- c("treat", "rank2")
  dat2$out2 <- outcomes[2]
  #
  dat2 %<>% filter(treat %in% common_trts)
  #
  dat <- merge(dat1, dat2, by = "treat", all.x = TRUE, all.y = TRUE)
  #
  plot(dat$rank1, dat$rank2, main = "",
       cex = cex.point,
       xlab = outcomes[1], ylab = outcomes[2],
       pch = pch, xlim = xlim, ylim = ylim, ...)
  #
  text(dat$rank1, dat$rank2, labels = dat$treat,
       cex = cex.label, pos = pos, col = "black")
  #
  invisible(NULL)
}
