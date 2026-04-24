spie.chart_internal <- function(x, weights = NULL){
  
  # spie.data <- data.frame(x = x, theta = theta)
  
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
  
  # transform weights to sum up to 2*pi
  
  weights <- 2*pi*weights
  
  area <- vector("list")
  
  for(i in 1:nrow(x)){
    
    for(j in 1:length(weights)){
    
  area[[i]] <- (1/(2*pi))*sum(weights[j]*x[i,]^2)
  
    }
  }
  
  area <- data.frame ("Area" = unlist(area))
  
  row.names(area) <- row.names(x)
  
  area %<>% arrange(desc(Area)) %<>% mutate("Area" = round(Area,digits = 2)) 
  
  attr(area, "transformed.weights") <- weights
  
  return(area)
}