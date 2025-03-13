#' Print results for treatment ranking in multivariate network meta-analysis
#' 
#' @description
#' Print results for treatment ranking in multivariate network meta-analysis
#' 
#' @param x An object of class \code{\link{mvrank}}.
#' @param digits Minimal number of significant digits, see
#'   \code{print.default}.
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
#' # Perform analysis considering the efficacy outcomes
#' p12 <- list(p1, p2)
#'
#' # Use 'mvdata()' to transform the data in suitable JAGS format
#' data12 <- mvdata(p12)
#'
#' # Define outcome labels
#' outcomes <- c("Early Response", "Early Remission")
#'  
#' # Fit the model combining only the two efficacy outcomes
#' set.seed(1909)
#' mvnma12 <- mvnma(data = data12, 
#'   reference.group = "Placebo", outclab = outcomes,
#'   n.iter = 10, n.burnin = 1)
#' 
#' mvrank(mvnma12, small = c("und", "und"))
#' }
#'
#' 
#' @method print mvrank 
#' @export

print.mvrank <- function(x,
                         digits = gs("digits"),
                         ...) {
  
  chkclass(x, "mvrank")
  #
  chknumeric(digits, min = 0, length = 1)
  #
  nam <- names(x)
  
  # Get rid of warning "no visible binding for global variable"
  treatment <- NULL
  #
  for (i in seq_along(nam)) {
    cat(paste0(if (i > 1) "\n" else "", "Outcome: ", nam[i], "\n\n"))
    #
    dat.i <- x[[i]]
    rownames(dat.i) <- dat.i$treatment
    dat.i %<>% select(-treatment)
    #
    nam.i <- names(dat.i)
    #
    for (j in nam.i)
      dat.i[[j]] <- formatN(dat.i[[j]], digits = digits)
    #
    prmatrix(dat.i, quote = FALSE, right = TRUE)
  }
  #
  invisible(NULL)
}
