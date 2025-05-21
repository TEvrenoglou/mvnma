#' Scatter plot to visualize the ranking lists for two outcomes.
#' 
#' @description
#' Draw a scatter plot (using grid graphics system) in the active graphics
#' window or store the forest plot in a file.
#' 
#' @param x An object of class \code{\link{mvrank}}.
#' @param which A mandatory numeric vector of length 2 specifying which
#'   outcomes should be plotted. For example, setting "outcome = c(2, 3)"
#'   implies that a scatter plot will be generated plotting the rankings
#'   of outcomes 2 and 3.
#' @param pos Position of treatment labels.
#' @param cex.point ...
#' @param cex.label ...
#' @param pch ...
#' @param xlim ...
#' @param ylim ...
#' @param ... Additional arguments for \code{\link{plot}} function.
#' 
#' @examples
#' \donttest{
#' library("netmeta")
#' 
#' # Use 'pairwise' to obtain contrast based data for the first two outcomes
#' data("Linde2015")
#' # Early response
#' p1 <- pairwise(treat = list(treatment1, treatment2, treatment3),
#'   event = list(resp1, resp2, resp3), n = list(n1, n2, n3),
#'   studlab = id, data = dat.linde2015, sm = "OR")
#' # Early remissions
#' p2 <- pairwise(treat = list(treatment1, treatment2, treatment3),
#'   event = list(remi1, remi2, remi3), n = list(n1, n2, n3),
#'   studlab = id, data = dat.linde2015, sm = "OR")
#' 
#' # Define outcome labels
#' outcomes <- c("Early_Response", "Early_Remission")
#' 
#' # Fit the model combining the two efficacy outcomes
#' set.seed(1909)
#' mvnma12 <- mvnma(p1, p2,
#'   reference.group = "Placebo", outclab = outcomes[1:2],
#'   n.iter = 1000, n.burnin = 100)
#' mvnma12
#' 
#' # Rank treatments using SUCRAs
#' ranks12 <- mvrank(mvnma12, small.values = c("und", "und"), method = "sucra")
#' ranks12
#' 
#' # Visualize SUCRAs in a scatter plot with outcome 1
#' # (as specified in the mvdata() function) in the x-axis and outcome 2
#' # (as specified in the mvdata() function) in the y-axis
#' plot(ranks12)
#' 
#' # Visualize SUCRAs in a scatter plot with outcome 2
#' # (as specified in the mvdata() function) in the x-axis and outcome 1
#' # (as specified in the mvdata() function) in the y-axis
#' plot(ranks12, which = 2:1)
#' }
#' 
#' @method plot mvrank 
#' @export  

plot.mvrank <- function(x, which = 1:2,
                        pos = 1,
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
