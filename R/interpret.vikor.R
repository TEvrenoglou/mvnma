#' Evaluate the conditions of the VIKOR method and return the set of compromise
#' solutions.
#' 
#' @description
#' This function uses the three ranking matrices Q, S and R obtained by the
#' VIKOR method and evaluates the conditions C1 and C2 the identify the set of
#' compromise solutions.
#' 
#' @param x An object of class \code{\link{vikor}}.
#' @param \dots Additional arguments (ignored).
#' 
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

#' # Rank treatments using sucra
#' 
#' ranks_sucra <- mvrank(mvmodel_effic,small.values = c("undesirable","undesirable"), method = "sucra")
#'                     
#' ranks_sucra
#' 
#' # Get the best compromise solution across all Efficacy outcomes
#' 
#' vikor(ranks_sucra)
#' 
#' # Add larger weight for Response than Remission
#' 
#' vikor(ranks_sucra,weights=c(0.6,0.3))
#'
#' @method print vikor
#' @export

print.vikor <- function(x, ...) {
  
  chkclass(x, "vikor")
  
  Q <- x$Q
  S <- x$S
  R <- x$R
  
  trts <- row.names(Q)
  
  DQ <- 1 / (length(trts) - 1)
  
  cond1 <- Q$Q[2] - Q$Q[1] >= DQ
  
  cond2_1 <- isTRUE(row.names(Q)[1] == row.names(S)[1])
  
  cond2_2 <- isTRUE(row.names(Q)[1] == row.names(R)[1])
  
  cond2 <- isTRUE(cond2_1 & cond2_2)
  
  if ((cond1) & (cond2)) {
    
    solution <- row.names(Q)[1]  
    
    res1 <- paste("The compromise treatment across all outcomes is: ", solution)
  }
  else if ((cond1) & (!cond2)) {
    solution <- paste(row.names(Q)[1:2], collapse = ", ")
    res1 <- paste("The compromise set of treatments across all outcomes are: ",
                  solution )
  }
  else if (!cond1) {
    
    compr <- ifelse(Q$Q-Q$Q[1]<DQ,TRUE,FALSE)
    
    E <- which(compr)
    
    solution <- paste(c(row.names(Q)[E]), collapse = ", ")
    
    res1 <- paste("The compromise set of treatments across all outcomes are: ", solution)
    
  }
  else if ((!cond1) & (!cond2)) {
    res1 <- paste("No compromise solution was identified. Please consider different outcome weights.")
  }
  
  cat(paste0(res1, "\n"))
  #
  invisible(NULL)
}
