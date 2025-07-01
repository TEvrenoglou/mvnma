#' Rank treatments across all outcomes using the VIKOR multi-criteria decision
#' analysis method.
#' 
#' @description
#' This function employs the VIKOR method to analyze all outcome-specific
#' ranking lists. It provides both an amalgamated ranking list and guidance on
#' which treatments correspond to the best compromise solutions.
#' 
#' @param x An object of class \code{\link{mvrank}} or a matrix.
#' @param weights Outcome weights. The weights should always sum to 1. If not
#'   then they are standardized. If NULL, the function will assume equal outcome
#'   weights. 
#' @param v A scalar from 0 to 1 interpreted as the weight of the decision
#'   making process. Following guidance from the multi-criteria decision
#'   analysis field it is set to 0.5.
#' @param \dots Additional arguments.
#'
#' @details
#' This function takes a single mandatory argument, which is either an object of class
#' \code{\link{mvrank}} or a matrix. It then uses the multi-criteria decision
#' analysis method VIKOR to produce an amalgamated ranking list across all
#' outcomes. The standard VIKOR approach is applied when the \code{method} argument
#' is set to \code{"sucra"} or \code{"pBV"} in \code{\link{mvrank}}. A fuzzy
#' VIKOR method is applied when outcome-specific rankings are expressed in terms of
#' median ranks and 95% credible intervals. The latter is possible when the
#' \code{\link{mvrank}} object is created with \code{method = "ranks"}.
#' In both cases, the final ranking list is calculated based on treatments
#' common across all outcomes. Treatments not present across all outcomes are
#' excluded internally.
#' 
#' Using the argument 'weights' the users can specify the weight that each
#' outcome should have in the decision making process. For each outcome this
#' argument should have a value from 0 to 1 while the sum of all outcome
#' weights should be 1. If the sum of all weights is not 1, then these are
#' internally standardize to achieve this. The standardized weight values
#' are returned as a message to the user. Finally, if NULL then equal weights
#' are assumed across all outcomes.
#'
#' The argument 'v' specifies the weight of the decision making process.
#' The VIKOR method is a compromise programming approach that aims to balance
#' between each treatments overall and worst performance across all outcomes.
#' The balance between these two criteria is achieved using the parameter 'v'
#' which takes values from 0 to 1. Values close to 1 will give more weight to
#' the treatment's overall performance while values close to 0 will give more
#' weight to penalize the treatment's worst performance. The most common
#' choice of 'v' is typically 0.5 (default also here), thereby allowing for a
#' balanced decision making between treatment's overall and worst performance.
#' 
#' @return
#' The function returns a 'vikor' object. This consists of three ranking lists
#' which are the following:
#' \itemize{
#' \item A ranking list Q referring to the ranking when balancing both each
#'   treatment's overall and worst performance. This is the main ranking list
#'   of the method. 
#' \item A ranking list S referring to the ranking in terms of each treatment's
#'   overall performance.
#' \item A ranking list R referring to the ranking in terms of penalizing each
#'   treatment's worst performance.
#' }
#' In addition to the ranking lists, the function also evaluates the necessary conditions
#' defined by the VIKOR method and returns a message indicating the set of compromise
#' solutions.
#'
#' @references 
#' Opricovic, S., Tzeng, G. H. Compromise solution by MCDM methods:
#' A comparative analysis of VIKOR and TOPSIS. \emph{European Journal of Operational
#' Research}. 2004; \bold{156} (2), 445–455 https://doi.org/10.1016/S0377-2217(03)00020-1
#' 
#' Opricovic, S. Fuzzy VIKOR with an application to water resources planning. 
#' \emph{Expert Systems with Applications} 2011; \bold{38}(10), 12983-12990. https://doi.org/10.1016/j.eswa.2011.04.097
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
#' # Fit the model combining only the two efficacy outcomes
#' set.seed(1909)
#' mvnma12 <- mvnma(p1, p2,
#'   reference.group = "Placebo", outclab = outcomes,
#'   n.iter = 1000, n.burnin = 100)
#' mvnma12
#' 
#' # Rank treatments using SUCRAs
#' ranks12 <- mvrank(mvnma12, small.values = c("und", "und"), method = "sucra")
#' ranks12
#' 
#' # Get the best compromise solution across the efficacy outcomes
#' vikor(ranks12)
#' 
#' # Use larger weight for response than remission
#' vikor(ranks12, weights = c(0.6, 0.3))
#' }
#'
#' @export vikor

vikor <- function(x, ...)
  UseMethod("vikor")


#' @rdname vikor
#' @method vikor mvrank
#' @export

vikor.mvrank <- function(x, weights = NULL, v = 0.5, ...) {
  
  chkclass(x, "mvrank")
  #
  trts <- sort(attr(x, "common_trts"))
  outcomes <- names(x)
  
  if (attr(x, "method") %in% c("SUCRA", "pBV")) {
    # Get rid of warning "no visible binding for global variable"
    treatment <- NULL
    
    s <- vector("list")
    #
    for (i in seq_len(length(x))) {
      dat.i <- x[[i]] %>% 
        filter(treatment %in% trts) %>% 
        arrange(treatment)
      #
      s[[i]] <- as.data.frame(dat.i[, 2])
      #
      names(s[[i]]) <- paste("ranks", i, sep = "_")
    }
    #
    rankings <- bind_cols(s)
    #
    row.names(rankings) <- dat.i$treatment
    names(rankings) <- outcomes
    #
    res <- vikor_internal(rankings, weights = weights, v = v)
  }
  else {
    decision <- performance_fuzzy(x,trts)
    res <- fuzzy_vikor_internal(decision, weights = weights,v = v )
  }
  #
  attr(res, "ranking.method") <- attr(x, "method")
  #
  res
}


#' @rdname vikor
#' @method vikor matrix
#' @export

vikor.matrix <- function(x, weights = NULL, v = 0.5, ...) {
  
  chkclass(x, "matrix")
  #
  res <- vikor_internal(x, weights = weights, v = v)
  #
  res
}
