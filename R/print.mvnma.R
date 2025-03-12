#' Print results of multivariate network meta-analysis
#' 
#' @description
#' Print results of multivariate network meta-analysis
#' 
#' @param x An object of class \code{\link{mvnma}}.
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
#' outcomes <- c("Early_Response", "Early_Remission")
#'  
#' # Fit the model combining only the two efficacy outcomes
#' set.seed(1909)
#' mvnma12 <- mvnma(data = data12, 
#'   reference.group = "Placebo", outlab = outcomes,
#'   n.iter = 1000, n.burnin = 100)
#' mvnma12
#' }
#'
#' 
#' @method print mvnma 
#' @export

print.mvnma <- function(x, ...) {
  
  chkclass(x, "mvnma")
  #
  x <- x[names(x) != "outcome_correlation"]
  nam <- names(x)
  
  for (i in seq_along(nam)) {
    cat(paste0(if (i > 1) "\n" else "", "Outcome: ", nam[i], "\n\n"))
    #
    dat.i <- x[[i]]$basic_estimates
    rownames(dat.i) <- paste0("d[", rownames(dat.i), "]")
    prmatrix(dat.i, quote = FALSE, right = TRUE)
  }
  #
  invisible(NULL)
}
