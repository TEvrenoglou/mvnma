#' Rank treatments across all outcomes using the spie chart method 
#' 
#' @description
#' This function employs the spie chart approach to combine
#' outcome-specific treatment hierarchies obtained in terms of the
#' probabilistic ranking metrics:
#' \itemize{
#' \item Surface Under the Cumulative Ranking curve (SUCRA), or
#' \item probability of best value (pBV).
#' }
#' 
#' @param x An object of class \code{\link{mvrank}}.
#' @param weights Outcome weights. The weights should always sum to 1. If not
#'   then they are standardized. If NULL, the function will assume equal outcome
#'   weights. 
#' @param \dots Additional arguments passed to \code{\link[ggplot2]{ggplot}}
#'   function.
#' 
#' @return
#' A data frame with the area under the spie chart as amalgamated treatment
#' hierarchy across outcomes.
#' 
#' @references 
#' Daly CH, Mbuagbaw L, Thabane L, Straus SE, Hamid JS (2020):
#' Spie charts for quantifying treatment effectiveness and safety in multiple
#' outcome network meta-analysis: a proof-of-concept study.
#' \emph{BMC Med Res Methodol}, \bold{20}, 266
#' 
#' @details
#' The spie chart is a modified pie chart in which each sector corresponds to
#' a treatment and its radius is proportional to a ranking metric, so that
#' sector area encodes treatment performance (Daly et al., 2020). This function
#' constructs one spie chart per outcome using SUCRA or probability of best
#' values (pBV) metrics, where each metric determines the radius of the
#' corresponding sector.
#' 
#' An amalgamated treatment hierarchy across outcomes is then derived by
#' averaging the sector areas across outcomes, providing a quantitative
#' summary of overall treatment performance.
#' 
#' @examples
#' # Locate file "mvnma_example.rda" with mvnma() results
#' .fname <- system.file("extdata/mvnma_examples.rda", package = "mvnma")
#' load(.fname)
#' 
#' # Rank treatments using SUCRAs
#' ranks12 <- mvrank(mvnma12, method = "sucra",
#'   small.values = c("undes", "undes"))
#' ranks12
#' 
#' # Get amalgamated treatment hierarchy
#' spiechart(ranks12)
#' 
#' @rdname spiechart
#' @method spiechart mvrank
#' @export

spiechart.mvrank <- function(x, weights = NULL, ...) {
  
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
    names(rankings) <- outcomes
    row.names(rankings) <- dat.i$treatment
    #
    res <- spiechart_internal(rankings, weights)
  }
  else
    stop("Spie chart area can only be calculated for method SUCRA and pBV")
  
  class(res) <- c("spiechart", class(res))
  attr(res, "weights") <- weights
  #
  res
}


#' @rdname spiechart
#' @export spiechart

spiechart <- function(x, ...)
  UseMethod("spiechart")


#' @method print spiechart
#' @export

print.spiechart <- function(x, digits = 2, ...) {
  chkclass(x, "spiechart")
  #
  # Get rid of warning 'Undefined global functions or variables'
  #
  area <- NULL
  #
  x %<>%
    rename(Area = area) %>%
    round(digits = digits)
  #
  prmatrix(x, quote = FALSE, right = TRUE)
  #
  invisible(NULL)
}
