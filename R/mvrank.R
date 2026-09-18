#' Outcome-specific treatment rankings in multivariate network meta-analysis.
#' 
#' @description
#' Produces outcome-specific treatment rankings in multivariate
#' network meta-analysis based on the output of the \code{\link{mvnma}}
#' function.
#' 
#' @param x An object of class \code{\link{mvnma}}.
#' @param small.values A character vector specifying for each outcome whether
#'   small treatment effects indicate a beneficial ("desirable") or harmful
#'   ("undesirable") effect, can be abbreviated.
#' @param method The ranking method to be used. Three methods are currently
#'   supported. The SUCRA method (specified as \code{method = "SUCRA"}) is the
#'   default approach. The probability of being best method (specified as
#'   \code{method = "pbest"}) and the mean and median ranks (specified as
#'   \code{method = "ranks"}) are also supported.
#' @param digits Minimal number of significant digits, see
#'   \code{print.default}.
#' @param \dots Additional arguments (ignored)
#' 
#' @description
#' This function produces outcome-specific treatment rankings in multivariate
#' network meta-analysis based on the output of the \code{\link{mvnma}}
#' function. Two ranking methods (argument \code{method}) are currently
#' supported (Salanti et al., 2011):
#' \itemize{
#' \item Surface under the cumulative ranking curve (SUCRA) method,
#' \item Probability of being best (pbest) method.
#' }
#'  
#' @return
#' The function returns an object of class \code{mvrank}. It is a list
#' containing the following components:
#' \item{ranks}{A named list containing one data frame for each outcome. The
#'   data frame contains the treatment rankings. For \code{method = "SUCRA"}
#'   or \code{method = "pbest"}, it contains the variables \code{treatment}
#'   and the corresponding ranking measure. For \code{method = "ranks"}, it
#'   contains \code{treatment}, \code{median_rank}, \code{mean_rank},
#'   \code{lower.CrI}, and \code{upper.CrI}.}
#' \item{trts}{All treatments occurring in the outcome-specific rankings.}
#' \item{ranks.shared}{A named list containing rankings recalculated using
#'   only treatments shared by all outcomes.}
#' \item{trts.shared}{Treatments occurring in every outcome.}
#' \item{outcomes}{Outcome labels.}
#' \item{method}{The ranking method used.}
#' \item{call, version}{The matched function call and package version used to
#'   create the object.}
#' 
#' @references
#' Salanti G, Ades AE, Ioannidis JP (2011):
#' Graphical methods and numerical summaries for presenting results
#' from multiple-treatment meta-analysis: an overview and tutorial.
#' \emph{Journal of Clinical Epidemiology},
#' \bold{64}, 163--71
#' 
#' @examples
#' # Locate file "mvnma_examples.rda" with mvnma() results
#' .fname <- system.file("extdata/mvnma_examples.rda", package = "mvnma")
#' load(.fname)
#' 
#' # Rank treatments using SUCRAs (default)
#' ranks_sucra <- mvrank(mvnma_all, 
#'   small.values = c("undes", "undes", "des", "des", "des"))
#' #
#' ranks_sucra
#' 
#' \donttest{
#' # Rank treatments using pbest
#' ranks_pbest <- mvrank(mvnma_all,
#'   small.values = c("undes", "undes", "des", "des", "des"),
#'   method = "pbest")
#' #
#' ranks_pbest         
#' 
#' # Rank treatments using mean and median ranks
#' ranks_mean_median <- mvrank(mvnma_all,
#'   small.values = c("undes", "undes", "des", "des", "des"),
#'   method = "ranks")
#' #
#' ranks_mean_median
#' }
#' 
#' @export mvrank                

mvrank <- function(x, small.values, method = "SUCRA") {
  
  chkclass(x, "mvnma")
  #
  small.values <- setchar(small.values, c("undesirable", "desirable"))
  #
  method <- setchar(method, c("SUCRA", "pbest", "ranks", "pBV"))
  if (any(method == "pBV")) {
    warning("Argument 'method = \"pBV\"' replaced by  'method = \"pbest\"'.",
            call. = FALSE)
    method[method == "pBV"] <- "pbest"
  }
  chkchar(method, length = 1)
  #
  method.model <- x$method.model
  n.domain <- x$n.domain
  #
  outcomes <- x$outcomes
  x <- x[outcomes]
  #
  colname_list <- lapply(x, function(k) colnames(k$samples))
  trts.common <- Reduce(intersect, colname_list)
  
  # Get rid of warning "no visible binding for global variable"
  treatment <- pbest <- SUCRA <- Freq <- median_rank <- mean_rank <-
    lower.CrI <- upper.CrI <- NULL
  
  # Extract samples and create rankograms for each outcome
  #
  n.out <- length(outcomes)
  #
  d <- d_common <- n.trts <- trts.list <- rank_out <- rank_out_common <-
    ranks <- ranks.common <- quant <- quant_common <-
    rnk <- rnk_common <- vector("list", n.out)
  #
  for (i in seq_len(n.out)) {
    d[[i]] <- x[[i]]$samples
    n.trts[[i]] <- ncol(d[[i]])
    trts.list[[i]] <- colnames(d[[i]])
    if(length(setdiff(trts.list[[i]], trts.common)) == 0){
      d_common[[i]] <- d[[i]]  
    }
    else{
      d_common[[i]] <- d[[i]] %>% select(trts.common)
    }
    
    rank_out[[i]] <- rankogram(d[[i]], small.values = small.values[i])
    rank_out_common[[i]] <- rankogram(d_common[[i]], small.values = small.values[i])
    #
    if (method == "pbest") {
      ranks.i <- data.frame(pbest = rank_out[[i]]$ranking.matrix.random[, 1])
      #
      ranks.i$treatment <- row.names(ranks.i)
      row.names(ranks.i) <- NULL
      #
      ranks.i %<>% select(treatment, pbest) %>% arrange(desc(pbest))
      
      # recalculate only for the common treatments
      ranks.i.common <- data.frame(pbest = rank_out_common[[i]]$ranking.matrix.random[, 1])
      #
      ranks.i.common$treatment <- row.names(ranks.i.common)
      row.names(ranks.i.common) <- NULL
      #
      ranks.i.common %<>% select(treatment, pbest) %>% arrange(desc(pbest))
    }
    else if (method == "SUCRA") {
      ranks.i <- rank_out[[i]]$ranking.random
      ranks.i <- cbind.data.frame(names(ranks.i), unname(ranks.i))
      names(ranks.i) <- c("treatment", "SUCRA")
      ranks.i %<>% arrange(desc(SUCRA))
      # Recalculate for the common treatments
      ranks.i.common <- rank_out_common[[i]]$ranking.random
      ranks.i.common <-
        cbind.data.frame(names(ranks.i.common), unname(ranks.i.common))
      names(ranks.i.common) <- c("treatment", "SUCRA")
      ranks.i.common %<>% arrange(desc(SUCRA))
    }
    else if (method == "ranks") {
      if (small.values[i] == "undesirable")
        rnk[[i]] <- apply(-d[[i]], 1, rank, ties.method = "random")
      else
        rnk[[i]] <- apply(d[[i]], 1, rank, ties.method = "random") 
      
      quant[[i]] <-
        as.data.frame(
          t(apply(rnk[[i]], 1,
                  function(row) {
                    quantile(row, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
                    }
                  )
            )
          )
      #
      quant[[i]]$treatment <- row.names(quant[[i]])
      quant[[i]]$mean_ranks <- rowMeans(rnk[[i]])
      #
      row.names(quant[[i]]) <- NULL
      names(quant[[i]]) <-
        c("lower.CrI", "median_rank", "upper.CrI", "treatment", "mean_rank")
      
      ranks.i <- quant[[i]] %<>%
        select(treatment, median_rank, mean_rank, lower.CrI, upper.CrI) %>%
        arrange(mean_rank)
      
      # Recalculate for the common treatments
      #
      if (small.values[i] == "undesirable")
        rnk_common[[i]] <-
        apply(-d_common[[i]], 1, rank, ties.method = "random")
      else
        rnk_common[[i]] <-
        apply(d_common[[i]], 1, rank, ties.method = "random") 
      
      quant_common[[i]] <-
        as.data.frame(
          t(apply(rnk_common[[i]], 1,
                  function(row) {
                    quantile(row, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
                    }
                  )
            )
          )
      #
      quant_common[[i]]$treatment <- row.names(quant_common[[i]])
      quant_common[[i]]$mean_ranks <- rowMeans(rnk_common[[i]])
      #
      row.names(quant_common[[i]]) <- NULL
      names(quant_common[[i]]) <-
        c("lower.CrI", "median_rank", "upper.CrI", "treatment", "mean_rank")
      
      ranks.i.common <- quant_common[[i]] %<>%
        select(treatment, median_rank, mean_rank, lower.CrI, upper.CrI) %>%
        arrange(mean_rank)
    }
    #
    ranks[[i]] <- ranks.i
    ranks.common[[i]] <- ranks.i.common
  }
  #
  names(ranks) <- outcomes
  class(ranks) <- "mvrank"
  #
  names(ranks.common) <- outcomes
  class(ranks.common) <- "mvrank"
  #
  ranks <- list(
    ranks = ranks,
    trts = sort(unique(unlist(trts.list))),
    ranks.shared = ranks.common,
    trts.shared = trts.common,
    outcomes = outcomes,
    method = method,
    call = match.call(),
    version = packageDescription("mvnma")$Version
  )
  #
  class(ranks) <- "mvrank"
  ranks
}


#' @rdname mvrank 
#' @method print mvrank
#' @export

print.mvrank <- function(x, digits = gs("digits"), ...) {
  
  chkclass(x, "mvrank")
  x <- updateversion(x)
  #
  chknumeric(digits, min = 0, length = 1)
  #
  nam <- x$outcomes
  
  # Get rid of warning "no visible binding for global variable"
  treatment <- upper.CrI <- NULL
  #
  x <- x$ranks
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
    if (all(c("lower.CrI", "upper.CrI") %in% names(dat.i))) {
      dat.i$lower.CrI <- formatCI(dat.i$lower.CrI, dat.i$upper.CrI)
      dat.i %<>% select(-upper.CrI)
      dat.i$spacer1 <- ""
      dat.i$spacer2 <- ""
      dat.i <- dat.i[c("median_rank", "spacer1", "mean_rank", "spacer2",
                       "lower.CrI")]
      names(dat.i) <- c("median_rank", "", "mean_rank", "", "lower.CrI")
      names(dat.i)[names(dat.i) == "lower.CrI"] <- "95% CrI for rank"
    }
    names(dat.i)[names(dat.i) == "pbest"] <- "P(best)"
    names(dat.i)[names(dat.i) == "median_rank"] <- "Median rank"
    names(dat.i)[names(dat.i) == "mean_rank"] <- "Mean rank"
    #
    prmatrix(dat.i, quote = FALSE, right = TRUE)
  }
  #
  invisible(NULL)
}
