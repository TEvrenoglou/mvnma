#' Transform contrast based data into JAGS format.
#' 
#' @description
#' This function transforms contrast-based data obtained from the
#' \code{\link[netmeta]{netmeta}} or \code{\link[meta]{pairwise}} function
#' into a format suitable for the multivariate model in JAGS.
#' 
#' @param p A list of pairwise objects. 
#' 
#' @return 
#' The function returns a list of vectors and matrices to be used as in input to the \code{\link{mvnma}} function.
#'
#' @examples
#' library(netmeta)
#' 
#' data("Linde2015")
#' 
#' # use 'pairwise' to obtain contrast based data for each one of the five available outcomes 
#'
#' # Early response
#'
#' p1 <- pairwise(treat = list(treatment1, treatment2, treatment3),
#'              event = list(resp1, resp2, resp3), 
#'               n = list(n1, n2, n3),
#'               studlab = id,
#'               data = dat.linde2015,
#'               sm = "OR")
#'
#'
#' # Early remissions
#'
#' p2 <- pairwise(treat = list(treatment1, treatment2, treatment3),
#'               event = list(remi1, remi2, remi3),
#'               n = list(n1, n2, n3),
#'               studlab = id,
#'               data = dat.linde2015,
#'               sm = "OR")
#'
#' # Perform analysis in terms of the Efficacy outcomes
#'
#' p_effic <- list(p1,p2)
#'
#' # Use 'mvdata()' to transform the data in suitable JAGS format
#' 
#' data_effic <- mvdata(p_effic)
#'
#' # Define outcome labels
#' 
#' outlab <- c("Early_Response",
#'              "Early_Remission")
#'              
#'          
#' # Fit the model combining only the two efficacy outcomes
#' 
#' mvmodel_effic <- mvnma(data = data_effic,
#'                 reference.group = "Placebo",
#'                 outlab = outlab,
#'                 n.iter = 1000,
#'                 n.burnin = 100)       
#'                 
#' # Plot the results for efficacy outcomes
#' 
#' forest(mvmodel_effic)
#' 
#' @export mvdata  

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
