#' Rank treatments across all outcomes using the VišeKriterijumska Optimizacija
#' I Kompromisno Rešenje (VIKOR) multi-criteria decision analysis method
#' 
#' @description
#' This function employs the VišeKriterijumska Optimizacija I Kompromisno
#' Rešenje (VIKOR) method to analyze all outcome-specific ranking lists.
#' It provides both an amalgamated ranking list and guidance on
#' which treatments correspond to the best compromise solutions.
#' 
#' @param x An object of class \code{\link{mvrank}} or a matrix.
#' @param weights Outcome weights. The weights should always sum to 1. If not
#'   then they are standardized. If NULL, the function will assume equal outcome
#'   weights. 
#' @param v A scalar from 0 to 1 interpreted as the weight of the decision
#'   making process. Following guidance from the multi-criteria decision
#'   analysis field it is set to 0.5.
#' @param digits A numeric specifying the number of digits to print the
#'   ranking matrix Q.
#' @param \dots Additional arguments (ignored).
#'
#' @details
#' This function takes a single mandatory argument, which is either an object
#' of class \code{\link{mvrank}} or a matrix. It then uses the multi-criteria
#' decision analysis method VišeKriterijumska Optimizacija I Kompromisno
#' Rešenje (VIKOR) to produce an amalgamated ranking list across
#' all outcomes (Opricovic & Tzeng, 2004).
#' 
#' The standard VIKOR approach is applied when the \code{method} argument is
#' set to \code{"sucra"} or \code{"pBV"} in \code{\link{mvrank}}.
#' 
#' A fuzzy VIKOR method (Opricovic, 2011) is applied when outcome-specific
#' rankings are expressed in terms of median ranks and 95\% credible intervals.
#' The latter is possible when the \code{\link{mvrank}} object is created with
#' \code{method = "ranks"}.
#' 
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
#' In addition to the ranking lists, the function also evaluates the necessary
#' conditions defined by the VIKOR method and returns a message indicating the
#' set of compromise solutions.
#'
#' @references
#' Opricovic S, Tzeng GH (2004):
#' Compromise solution by MCDM methods: A comparative analysis of VIKOR and
#' TOPSIS.
#' \emph{European Journal of Operational Research},
#' \bold{156}, 445--55
#' 
#' Opricovic S (2011):
#' Fuzzy VIKOR with an application to water resources planning. 
#' \emph{Expert Systems with Applications}, 
#' \bold{38}, 12983--90
#' 
#' @examples
#' # Locate file "mvnma_example.rda" with mvnma() results
#' .fname <- system.file("extdata/mvnma_examples.rda", package = "mvnma")
#' load(.fname)
#' 
#' # Print the results of a bivariate network meta-analysis
#' mvnma12
#' 
#' # Rank treatments using SUCRAs
#' ranks12 <- mvrank(mvnma12, method = "sucra",
#'   small.values = c("undes", "undes"))
#' ranks12
#' 
#' # Get the best compromise solution across the efficacy outcomes
#' vikor(ranks12)
#' 
#' # Use larger weight for response than remission
#' vikor(ranks12, weights = c(0.6, 0.3))
#'
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
    x.common <- attr(x,"ranks.common")
    #
    for (i in seq_len(length(x))) {
      dat.i <- x.common[[i]] %>% 
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
  class(res) <- c("vikor", class(res))
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


#' @rdname vikor
#' @export vikor

vikor <- function(x, ...)
  UseMethod("vikor")


#' @rdname vikor
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
    txt <- paste("The compromise set of treatments across all outcomes are:",
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
  else if (!cond1 & !cond2)
    txt <- paste("No compromise solution was identified. Please consider",
                 "different outcome weights.")
  
  res_mat <- cbind(Q, S, R)
  
  if (attr(x, "ranking.method") %in% c("SUCRA", "pBV"))
    cat("VIKOR results\n\n")
  else
    cat("Fuzzy VIKOR results\n\n")
  #
  prmatrix(round(res_mat, digits = digits), quote = FALSE, right = TRUE)
  #
  cat(paste0("\n", txt, "\n"))
  #
  invisible(NULL)
}
