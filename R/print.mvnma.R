#' Print results of multivariate network meta-analysis
#' 
#' @description
#' Print results of multivariate network meta-analysis
#' 
#' @param x An object of class \code{\link{mvnma}}.
#' @param digits Minimal number of significant digits, see
#'   \code{print.default}.
#' @param digits.sd Minimal number of significant digits for standard
#'   deviations
#' @param print.sd A logical specifying whether standard deviations should be
#'   printed.
#' @param \dots Additional arguments (ignored)
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
#' outcomes <- c("Early Response", "Early Remission")
#'  
#' # Fit the model combining the two efficacy outcomes
#' set.seed(1909)
#' mvnma12 <- mvnma(p1, p2,
#'   reference.group = "Placebo", outclab = outcomes,
#'   n.iter = 10, n.burnin = 1)
#' mvnma12
#' }
#'
#' @method print mvnma 
#' @export

print.mvnma <- function(x,
                        digits = gs("digits"),
                        digits.sd = gs("digits.sd"),
                        print.sd = FALSE,
                        ...) {
  
  chkclass(x, "mvnma")
  #
  chknumeric(digits, min = 0, length = 1)
  chknumeric(digits.sd, min = 0, length = 1)
  chklogical(print.sd)
  #
  level <- attr(x, "level")
  reference.group <- attr(x, "reference.group")
  #
  ci.lab <- paste0(round(100 * level, 1), "%-CI")
  #
  x <- x[names(x) != "cor"]
  nam <- names(x)
  
  # Get rid of warning "no visible binding for global variable"
  lower <- upper <- psi <- NULL
  #
  for (i in seq_along(nam)) {
    cat(paste0(if (i > 1) "\n" else "", "Outcome: ", nam[i], "\n\n"))
    #
    dat.i <- x[[i]]$basic_estimates
    dat.i <- dat.i[rownames(dat.i) != reference.group, ]
    #
    dat.i$mean <- formatN(dat.i$mean, digits = digits)
    #
    if (!print.sd)
      dat.i$sd <- NULL
    else
      dat.i$sd <- formatN(dat.i$sd, digits = digits.sd)
    #
    dat.i$lower <- formatCI(formatN(dat.i$lower, digits = digits),
                            formatN(dat.i$upper, digits = digits))
    dat.i %<>% select(-upper)
    names(dat.i)[names(dat.i) == "lower"] <- ci.lab
    #
    dat.i$Rhat <- formatN(dat.i$Rhat, digits = 4)
    #
    rownames(dat.i) <- paste0("d[", rownames(dat.i), "]")
    #
    prmatrix(dat.i, quote = FALSE, right = TRUE)
  }
  #
  invisible(NULL)
}
