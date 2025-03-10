#' Calculate all the outcome specific treatment effect estimates
#' 
#' @description
#' This function utilizes the consistency assumption within each outcome. Using the set of basic parameters, it calculates the treatment effect estimates and standard deviations for all outcome-specific treatment comparisons. 
#'
#' @param x An object of class \code{\link{mvnma}}.
#' 
#' @return A list of length equal to the total number of outcomes. Each element of the list contains one matrix for all outcome-specific treatment effect estimates and one matrix for all outcome-specific standard deviations of the treatment effect estimates.

#' @examples
#' library(netmeta)
#' 
#' data("Linde2015")
#' 
#' # use 'pairwise' to obtain contrast based data for each one of the five available outcomes 
#'
#'   # Early response
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
#' outlab <- c("Early_Response","Early_Remission")
#'             
#' # Fit the model combining only the two efficacy outcomes
#' 
#' mvmodel_effic <- mvnma(data = data_effic,
#'                 reference.group = "Placebo",
#'                 outlab = outlab,
#'                 n.iter = 1000,
#'                 n.burnin = 100)
#'                 
#' # Extract treatment effect estimates and heterogeneity for Early_Response 
#' 
#' mvmodel_effic$Early_Response$basic_estimates
#' 
#' # Get all estimates                 
#' 
#' league.effic <- league(mvmodel_effic)
#' 
#' # print results in terms of the outcome Early Response 
#' 
#' league.effic$Early_Response
#'
#' @export league 

league <- function(x){

  if(inherits(x,"mvnma")){
    
all_ests <- attributes(x)$all_res 
  
x <- x[names(x) != "outcome_correlation"]

outlab <- names(x)

n.out <- length(x)

res <- list()

for(i in 1:n.out){
  
res[[i]] <- list("TE.random"=all_ests$TE.random[[i]],
                 "sd.random"=all_ests$sd.random[[i]]
                 )  
  
  
}

names(res) <- outlab
  
return(res)
  }else{
  
    stop("Argument 'x' should be a 'mvnma' object.")
}

}
