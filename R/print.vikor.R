#' Evaluate the conditions of the VIKOR method and return the set of compromise
#' solutions.
#' 
#' @description
#' This function uses the three ranking matrices Q, S and R obtained by the
#' VIKOR method and evaluates the conditions C1 and C2 the identify the set of
#' compromise solutions.
#' 
#' @param x An object of class \code{\link{vikor}}.
#' @param digits A numeric specifying the number of digits to print the
#'   ranking matrix Q.
#' @param \dots Additional arguments (ignored).
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
#' # Get the best compromise solution across efficacy outcomes
#' vikor(ranks12)
#' 
#' # Use larger weight for response than remission
#' vikor(ranks12, weights = c(0.6, 0.3))
#' }
#'
#' @method print vikor
#' @export

print.vikor <- function(x, digits = 4, ...) {
  
  chkclass(x, "vikor")
  #
  chknumeric(digits, min = 0, length = 1)
  
  Q <- x %>% select(Q)
  S <- x %>% select(S)
  R <- x %>% select(R)
  #
  trts <- row.names(Q)
  #
  DQ <- 1 / (length(trts) - 1)
  
  cond1 <- Q$Q[2] - Q$Q[1] >= DQ
  #
  cond2_1 <- isTRUE(row.names(Q)[1] == row.names(S)[1])
  cond2_2 <- isTRUE(row.names(Q)[1] == row.names(R)[1])
  #
  cond2 <- isTRUE(cond2_1 & cond2_2)
  #
  if (cond1 & cond2) {
    solution <- row.names(Q)[1]
    #
    txt <- paste("The compromise treatment across all outcomes is:", solution)
  }
  else if ((cond1) & (!cond2)) {
    solution <- paste(row.names(Q)[1:2], collapse = ", ")
    #
    txt <-paste("The compromise set of treatments across all outcomes are:",
                solution)
  }
  else if (!cond1) {
    compr <- Q$Q - Q$Q[1] < DQ
    #
    E <- which(compr)
    #
    solution <- paste(row.names(Q)[E], collapse = ", ")
    #
    txt <- paste("The compromise set of treatments across all outcomes are:",
                 solution)
  }
  else if (!cond1 & !cond2) {
    txt <- paste("No compromise solution was identified. Please consider",
                 "different outcome weights.")
  }
  
  prmatrix(round(Q, digits = digits), quote = FALSE, right = TRUE)
  #
  cat(paste0("\n", txt, "\n"))
  #
  invisible(NULL)
}
