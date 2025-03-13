#' Calculate all the outcome specific treatment effect estimates
#' 
#' @description
#' This function utilizes the consistency assumption within each outcome.
#' Using the set of basic parameters, it calculates the treatment effect
#' estimates and standard deviations for all outcome-specific treatment
#' comparisons. 
#'
#' @param x An object of class \code{\link{mvnma}}.
#' 
#' @return A list of length equal to the total number of outcomes. Each
#' element of the list contains one matrix for all outcome-specific treatment
#' effect estimates and one matrix for all outcome-specific standard deviations
#' of the treatment effect estimates.
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
#'   n.iter = 1000, n.burnin = 100)
#' mvnma12
#' 
#' # Get all estimates
#' league12 <- league(mvnma12)
#' 
#' # Print results for the outcome Early Response 
#' league12$"Early Response"
#' }
#'
#' @export league 

league <- function(x) {
  
  chkclass(x, "mvnma")
  #
  x <- x[names(x) != "cor"]
  #
  outcomes <- names(x)
  #
  n.out <- length(x)
  
  res <- vector("list")
  #
  for (i in seq_len(n.out))
    res[[i]] <- list(TE.random = x[[i]]$TE.random,
                     seTE.random = x[[i]]$seTE.random)
  #
  names(res) <- outcomes
  
  res
}
