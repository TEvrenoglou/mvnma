#' Rank treatments across all outcomes using the spie chart method 
#' 
#' @description
#' This function employs the spie chart approach to combine outcome-specific treatment hierarchies 
#' obtained in terms of the probabilistic ranking metrics: (i) Surface Under the Cumulative Ranking curve (SUCRA) 
#' and (ii) probability of best value (pBV).
#' 
#' @param x An object of class \code{\link{mvrank}}.
#' @param weights Outcome weights. The weights should always sum to 1. If not
#'   then they are standardized. If NULL, the function will assume equal outcome
#'   weights. 
#' 
#' @references 
#' Daly CH, Mbuagbaw L, Thabane L, Straus SE, Hamid JS. Spie charts for quantifying treatment effectiveness 
#' and safety in multiple outcome network meta-analysis: a proof-of-concept study. \emph{BMC Med Res Methodol} 
#' 2020;\bold{20}(1):266. doi: 10.1186/s12874-020-01128-2.
#' 
#' @details
#' The spie chart is a modified pie chart in which each sector corresponds to a treatment and 
#' its radius is proportional to a ranking metric, so that sector area encodes treatment performance. 
#' This function constructs one spie chart per outcome using SUCRA and probability of best values (pBV) metrics, 
#' where each metric determines the radius of the corresponding sector. 
#' An amalgamated treatment hierarchy across outcomes is then derived by averaging the sector areas 
#' across outcomes, providing a quantitative summary of overall treatment performance.
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
#' spie.chart(ranks12)
#' }
#'
#' @export spie.chart

spie.chart <- function(x, ...)
  UseMethod("spie.chart")

#' @rdname spie.chart
#' @method spie.chart mvrank
#' @export


spie.chart.mvrank <- function(x, weights = NULL,...) {
  
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
    res <- spie.chart_internal(rankings, weights = weights)
  }
  else {
  stop("Spie chart area can only be calculated for method SUCRA and pBV")
    
  }
  
  class(res) <- c("spie.chart", class(res))
  
  attr(res, "weights") <- weights
  
  res
  }
  
  