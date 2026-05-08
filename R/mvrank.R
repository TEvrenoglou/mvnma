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
#' @param method The ranking method to be used. Three methods are currently
#'   supported. The SUCRA method (specified as \code{method = "SUCRA"}) is the
#'   default approach. The probability of best value method (specified as
#'   \code{method = "pBV"}) and the mean and median ranks (specified as
#'   \code{method = "ranks"}) are also supported.
#'
#' @examples
#' \donttest{
#' # Use 'pairwise' to obtain contrast based data for each one of the five
#' # available outcomes 
#'
#' # Early response
#' pw1 <- pairwise(treat = list(treatment1, treatment2, treatment3),
#'   event = list(resp1, resp2, resp3), n = list(n1, n2, n3),
#'   studlab = id, data = Linde2015, sm = "OR")
#'
#' # Early remissions
#' pw2 <- pairwise(treat = list(treatment1, treatment2, treatment3),
#'   event = list(remi1, remi2, remi3), n = list(n1, n2, n3),
#'   studlab = id, data = Linde2015, sm = "OR")
#'
#' # Adverse events
#' pw3 <- pairwise(treat = list(treatment1, treatment2,treatment3),
#'   event = list(ae1, ae2, ae3),  n = list(n1, n2, n3),
#'   studlab = id, data = Linde2015, sm = "OR")
#'
#' # Loss to follow-up
#' pw4 <- pairwise(treat = list(treatment1, treatment2, treatment3),
#'   event = list(loss1, loss2, loss3), n = list(n1, n2, n3),
#'   studlab = id, data = Linde2015, sm = "OR")
#'
#' # Loss_to_follow_up_(AE)
#' pw5 <- pairwise(treat = list(treatment1, treatment2, treatment3),
#'   event = list(loss.ae1, loss.ae2, loss.ae3), n = list(n1, n2, n3),
#'   studlab = id, data = Linde2015, sm = "OR")
#'
#' # Define outcome labels
#' outcomes <- c("Early_Response", "Early_Remission",
#'   "Adverse_events", "Loss_to_follow_up", "Loss_to_follow_up_AE")
#' 
#' # Fit the model combining all five outcomes
#' set.seed(1909)
#' mvnma_all <- mvnma(pw1, pw2, pw3, pw4, pw5,
#'   reference.group = "Placebo", outclab = outcomes,
#'   n.iter = 1000, n.burnin = 100)
#'
#' # Rank treatments using pBV
#' ranks_pBV <- mvrank(mvnma_all,
#'   small.values = c("undes", "undes", "des", "des", "des"),
#'   method = "pBV")
#' ranks_pBV                    
#'    
#' # Rank treatments using mean and median ranks
#' ranks_mean_median <- mvrank(mvnma_all,
#'   small.values = c("undes", "undes", "des", "des", "des"),
#'   method = "ranks")
#' ranks_mean_median
#'    
#' # Rank treatments using SUCRAs
#' ranks_sucra <- mvrank(mvnma_all, 
#'   small.values = c("undes", "undes", "des", "des", "des"),
#'   method = "SUCRA")
#' ranks_sucra
#'                      
#' # Same results without method = "SUCRA" 
#' ranks_sucra1 <- mvrank(mvnma_all, 
#'   small.values = c("undes", "undes", "des", "des", "des"))
#' ranks_sucra1
#' }
#' 
#' @export mvrank                

mvrank <- function(x, small.values, method = "SUCRA") {
  
  chkclass(x, "mvnma")
  #
  small.values <- setchar(small.values,c("undesirable", "desirable"))
  method <- setchar(method, c("SUCRA", "pBV", "ranks"))
  chkchar(method, length = 1)
  method.model <- attr(x, "method.model")
  n.domain <- attr(x,"n.domain")
  #
  x <- x[names(x) != "cor"]
  if (method.model == "DM") {
    if(is.null(n.domain)){
    x <- x[names(x) != "sigma"]
    }else{
      x <- x[!(names(x) %in% c("sigma1","sigma2"))]
    }
  }
  #
  outcomes <- attr(x, "names")
  #
  colname_list <- lapply(x, function(k) colnames(k$samples))
  common_trts <- Reduce(intersect, colname_list)
  #
  
  # Get rid of warning "no visible binding for global variable"
  treatment <- pBV <- SUCRA <- Freq <-
    median_rank <- mean_rank <- lower.CrI <- upper.CrI <- NULL
  
  # Extract samples and create rankograms for each outcome
  #
  n.out <- length(outcomes)
  d <- d_common <- n.trts <- trts <- rank_out <- rank_out_common <-  ranks <- ranks.common <-
    quant <- quant_common <- rnk <- rnk_common <- vector("list")
  
  for (i in seq_len(n.out)) {
    d[[i]] <- x[[i]]$samples
    n.trts[[i]] <- ncol(d[[i]])
    trts[[i]] <- colnames(d[[i]])
    if(length(setdiff(trts[[i]],common_trts))==0){
    d_common[[i]] <- d[[i]]  
    }else{
      d_common[[i]] <- d[[i]] %>% select(common_trts)
    }
    
    rank_out[[i]] <- rankogram(d[[i]], small.values = small.values[i])
    rank_out_common[[i]] <- rankogram(d_common[[i]], small.values = small.values[i])
    #
    if (method == "pBV") {
      ranks.i <- data.frame(pBV = rank_out[[i]]$ranking.matrix.random[, 1])
      #
      ranks.i$treatment <- row.names(ranks.i)
      row.names(ranks.i) <- NULL
      #
      ranks.i %<>% select(treatment, pBV) %>% arrange(desc(pBV))
      
      # recalculate only for the common treatments
      ranks.i.common <- data.frame(pBV = rank_out_common[[i]]$ranking.matrix.random[, 1])
      #
      ranks.i.common$treatment <- row.names(ranks.i.common)
      row.names(ranks.i.common) <- NULL
      #
      ranks.i.common %<>% select(treatment, pBV) %>% arrange(desc(pBV))
    }
    else if (method == "SUCRA") {
      ranks.i <- rank_out[[i]]$ranking.random
      ranks.i <- cbind.data.frame(names(ranks.i), unname(ranks.i))
      names(ranks.i) <- c("treatment", "SUCRA")
      ranks.i %<>% arrange(desc(SUCRA))
      # recalculate for the common treatments
      ranks.i.common <- rank_out_common[[i]]$ranking.random
      ranks.i.common <- cbind.data.frame(names(ranks.i.common), unname(ranks.i.common))
      names(ranks.i.common) <- c("treatment", "SUCRA")
      ranks.i.common %<>% arrange(desc(SUCRA))
    }
    else if (method == "ranks") {
      if (small.values[i] == "undesirable")
        rnk[[i]] <- apply(-d[[i]], 1, rank, ties.method = "random")
      else
        rnk[[i]] <- apply(d[[i]], 1, rank, ties.method = "random") 
      
      quant[[i]] <-
        as.data.frame(t(apply(rnk[[i]], 1,
                              function(row) {
                                quantile(row, probs = c(0.025, 0.5, 0.975),
                                         na.rm = TRUE)
                              })))
      
      quant[[i]]$treatment <- row.names(quant[[i]])
      quant[[i]]$mean_ranks <- rowMeans(rnk[[i]])
      #
      row.names(quant[[i]]) <- NULL
      names(quant[[i]]) <-
        c("lower.CrI", "median_rank", "upper.CrI", "treatment", "mean_rank")
      
      ranks.i <- quant[[i]] %<>%
        select(treatment, median_rank, mean_rank, lower.CrI, upper.CrI) %>%
        arrange(mean_rank)
      
      # recalculate for the common treatments
      
      if (small.values[i] == "undesirable")
        rnk_common[[i]] <- apply(-d_common[[i]], 1, rank, ties.method = "random")
      else
        rnk_common[[i]] <- apply(d_common[[i]], 1, rank, ties.method = "random") 
      
      quant_common[[i]] <-
        as.data.frame(t(apply(rnk_common[[i]], 1,
                              function(row) {
                                quantile(row, probs = c(0.025, 0.5, 0.975),
                                         na.rm = TRUE)
                              })))
      
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
  # common_trts <- as.data.frame(table(unlist(trts))) %>% filter(Freq == n.out)
  # common_trts <- common_trts$Var1
  #
  attr(ranks, "common_trts") <- common_trts
  attr(ranks, "method") <- method
  attr(ranks,"ranks.common.trts") <- ranks.common
  #
  ranks
}
