#' Outcome-specific treatment rankings in multivariate network meta-analysis.
#' 
#' @description
#' This function produces outcome-specific treatment rankings in multivariate
#' network meta-analysis based on the output of the \code{\link{mvnma}}
#' function. Two ranking methods are currently supported, these are: (i)
#' the SUCRA method and (ii) the probability of best value (pBV) method.
#' 
#' @param x An object of class \code{\link{mvnma}}.
#' @param small.values A character vector specifying for each outcome whether
#'   small treatment effects indicate a beneficial ("desirable") or harmful
#'   ("undesirable") effect, can be abbreviated.
#' @param method The ranking method to be used. Two methods are supported,
#'   the SUCRA method (specified as method = "SUCRA", default) and the
#'   probability of best value method (specified as method = "pBV").
#'
#' @examples
#' \donttest{
#' library("netmeta")
#' 
#' data("Linde2015")
#' 
#' # Use 'pairwise' to obtain contrast based data for each one of the five
#' # available outcomes 
#'
#' # Early response
#' p1 <- pairwise(treat = list(treatment1, treatment2, treatment3),
#'   event = list(resp1, resp2, resp3), n = list(n1, n2, n3),
#'   studlab = id, data = dat.linde2015, sm = "OR")
#'
#' # Early remissions
#' p2 <- pairwise(treat = list(treatment1, treatment2, treatment3),
#'   event = list(remi1, remi2, remi3), n = list(n1, n2, n3),
#'   studlab = id, data = dat.linde2015, sm = "OR")
#'
#' # Adverse events
#' p3 <- pairwise(treat = list(treatment1, treatment2,treatment3),
#'   event = list(ae1, ae2, ae3),  n = list(n1, n2, n3),
#'   studlab = id, data = dat.linde2015, sm = "OR")
#'
#' # Loss to follow-up
#' p4 <- pairwise(treat = list(treatment1, treatment2, treatment3),
#'   event = list(loss1, loss2, loss3), n = list(n1, n2, n3),
#'   studlab = id, data = dat.linde2015, sm = "OR")
#'
#' # Loss_to_follow_up_(AE)
#' p5 <- pairwise(treat = list(treatment1, treatment2, treatment3),
#'   event = list(loss.ae1, loss.ae2, loss.ae3), n = list(n1, n2, n3),
#'   studlab = id, data = dat.linde2015, sm = "OR")
#'
#' # Define outcome labels
#' outcomes <- c("Early_Response", "Early_Remission",
#'   "Adverse_events", "Loss_to_follow_up", "Loss_to_follow_up_AE")
#' 
#' # Perform analysis considering all outcomes
#' p_all <- list(p1, p2, p3, p4, p5)
#' data_all <- mvdata(p_all)
#'  
#' # Fit the model combining only the two efficacy outcomes
#' set.seed(1909)
#' mvnma_all <- mvnma(data = data_all,
#'   reference.group = "Placebo", outclab = outcomes,
#'   n.iter = 1000, n.burnin = 100)
#'
#' # Rank treatments using pBV
#' ranks_pBV <- mvrank(mvnma_all,
#'   small.values = c("undes", "undes", "des", "des","des"),
#'   method = "pBV")
#' ranks_pBV                    
#'    
#' # Rank treatments using SUCRAs
#' ranks_sucra <- mvrank(mvnma_all, 
#'   small.values = c("undes", "undes", "des", "des","des"),
#'   method = "SUCRA")
#' ranks_sucra
#'                      
#' # Same results without method = "SUCRA" 
#' ranks_sucra1 <- mvrank(mvnma_all, 
#'   small.values = c("undes", "undes", "des", "des","des"))
#' ranks_sucra1
#' }
#' 
#' @export mvrank                

mvrank <- function(x, small.values, method = "SUCRA") {
  
  chkclass(x, "mvnma")
  #
  method <- setchar(method, c("SUCRA","pBV"))
  chkchar(method, length = 1)
  #
  x <- x[names(x) != "cor"]
  #
  outcomes <- attributes(x)$names
  
  # Get rid of warning "no visible binding for global variable"
  treatment <- pBV <- SUCRA <- Freq <- NULL
  
  # Extract samples and create rankograms for each outcome
  #
  n.out <- length(outcomes)
  d <- n.trts <- trts <- rank_out <- ranks <- vector("list")
  
  for (i in seq_len(n.out)) {
    d[[i]] <- x[[i]]$samples
    n.trts[[i]] <- ncol(d[[i]])
    trts[[i]] <- colnames(d[[i]])
    rank_out[[i]] <- rankogram(d[[i]], small.values = small.values[i])
    #
    if (method == "pBV") {
      ranks.i <- data.frame(pBV = rank_out[[i]]$ranking.matrix.random[, 1])
      #
      ranks.i$treatment <- row.names(ranks.i)
      row.names(ranks.i) <- NULL
      #
      ranks.i %<>% select(treatment, pBV) %>% arrange(desc(pBV))
    }
    else if (method == "SUCRA") {
      ranks.i <- rank_out[[i]]$ranking.random
      ranks.i <- cbind.data.frame(names(ranks.i), unname(ranks.i))
      names(ranks.i) <- c("treatment", "SUCRA")
      ranks.i %<>% arrange(desc(SUCRA))
    }
    #
    ranks[[i]] <- ranks.i
  }
  #
  names(ranks) <- outcomes
  class(ranks) <- "mvrank"
  #
  common_trts <- as.data.frame(table(unlist(trts))) %>% filter(Freq == n.out)
  common_trts <- common_trts$Var1
  #
  attr(ranks, "common_trts") <- common_trts
  attr(ranks, "method") <- method
  #
  ranks
}
