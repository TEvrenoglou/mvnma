#' Hasse diagram
#' 
#' @description
#' This function generates a Hasse diagram for a partial order of treatment
#' ranks in a multivariate network meta-analysis.
#' 
#' @param x An object of class \code{\link{mvrank}}.
#' @param \dots Additional arguments passed on to
#'   \code{\link[netmeta]{hasse.netposet}}.
#' 
#' @details
#' Generate a Hasse diagram (Carlsen & Bruggemann, 2014) for a partial order of
#' treatment ranks in a network meta-analysis (Rücker & Schwarzer, 2017).
#' 
#' This R function is a wrapper function for
#' \code{\link[netmeta]{hasse.netposet}}.
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
#'   n.iter = 1000, n.burnin = 100)
#' mvnma12
#' 
#' # Get all estimates
#' league12 <- league(mvnma12)
#' league12
#' 
#' # Rank treatments using sucra
#' ranks12 <- mvrank(mvnma12, small.values = c("und", "und"), method = "sucra")
#' ranks12
#' 
#' # Get the Hasse diagram for the efficacy outcomes
#' hasse(ranks12)
#' }
#' 
#' @references
#' Carlsen L, Bruggemann R (2014):
#' Partial order methodology: a valuable tool in chemometrics.
#' \emph{Journal of Chemometrics}, \bold{28}, 226--34
#'
#' Rücker G, Schwarzer G (2017):
#' Resolve conflicting rankings of outcomes in network meta-analysis:
#' Partial ordering of treatments.
#' \emph{Research Synthesis Methods}, \bold{8}, 526--36
#' 
#' @method hasse mvrank
#' @export

hasse.mvrank <- function(x, ...) {
  
  if (!(attr(x, "method") %in% c("SUCRA", "pBV")))
    stop("Hasse diagram can only be produced for ",
         "'method=SUCRA' and 'method=pBV'.",
         call. = FALSE)
  
  # Get rid of warning "no visible binding for global variable"
  treatment <- NULL
  
  treats <- E <- new_treats <- vector("list")
  #
  for (i in seq_along(x))
    treats[[i]] <- x[[i]]$treatment  
  #
  all_treats <- unique(unlist(treats))
  
  for (i in seq_along(x)) {
    E[[i]] <- which(!(all_treats %in% x[[i]]$treatment))
    #
    if (length(E[[i]]) > 0) {
      new_treats[[i]] <- cbind.data.frame(all_treats[E[[i]]], NA)
      names(new_treats[[i]]) <- names(x[[i]])
      #
      x[[i]] <- rbind.data.frame(x[[i]], new_treats[[i]])
    }
    #
    x[[i]] %<>% arrange(treatment)
    row.names(x[[i]]) <- x[[i]]$treatment
    #
    x[[i]]$treatment <- NULL
  } 
  #
  ranking_mat <- list.cbind(x) 
  names(ranking_mat) <- names(x)
  #
  po <- netposet(ranking_mat, outcomes = names(ranking_mat))
  #
  hasse(po, ...)
  #
  invisible(NULL)
}
