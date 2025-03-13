mvdata <- function(p) {
  
  if (is.list.pairwise(p) != "pairwise"){
    
    stop("Argument 'p' must be a list of 'pairwise' objects.")  
    
  }
  else {
    data_format <- create_data(p) 
    jags_data <- make_jags_data(data_format) 
  }
  #
  class(jags_data) <- "mvdata"
  attr(jags_data, "structured_data") <- data_format
  
  sm <- vector("character")
  #
  for (i in seq_along(p))
    sm[i] <- attributes(p[[i]])$sm 
  #
  attr(jags_data, "sm") <- sm
  
  jags_data
}
