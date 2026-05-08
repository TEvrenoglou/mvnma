spiechart_internal <- function(x, weights = NULL){
  
  # Get rid of warning 'Undefined global functions or variables'
  #
  area <- NULL
  
  if (is.null(weights))
    weights <- rep(1 /ncol(x), ncol(x))
  else if ((!is.null(weights) & length(weights) != ncol(x)))
    stop("Please provide weights for all outcomes.",
         call. = FALSE)
  #
  if (!is.null(weights) && sum(weights) != 1) {
    weights <- round(weights / sum(weights), digits = 2)
    #
    warning("Weights should always sum up to 1. To do so the given weights ",
            "are now standardized and the new weights are: ",
            paste(weights, collapse = ", "))
  }
  
  # Transform weights to sum up to 2 * pi
  #
  weights <- 2 * pi * weights
  
  n.trts <- nrow(x)
  res <- data.frame(area = rep(NA, n.trts))
  row.names(res) <- row.names(x)
  for (i in seq_len(n.trts))
    res$area[i] <- (1 / (2 * pi)) * sum(weights * x[i, ]^2)
  #
  res %<>% arrange(desc(area))
  #
  attr(res, "transformed.weights") <- weights
  #
  res
}
