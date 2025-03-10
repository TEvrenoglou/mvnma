#' Outcome-specific treatment rankings in multivariate network meta-analysis.
#' 
#' @description
#' This function produces outcome-specific treatment rankings in multivariate network meta-analysis based on the output of the \code{\link{mvnma}} function.
#' Two ranking methods are currently supported, these are: (i) the SUCRA method and (ii) the probability of best value (pBV) method.
#' 
#' @param x An object of class \code{\link{mvnma}}.
#' @param small.values A character vector specifying for each outcome whether small treatment effects indicate a beneficial ("desirable") or harmful ("undesirable") effect, can be abbreviated.
#' @param method The ranking method to be used. Two methods are supported, the SUCRA method (specified as method="sucra") and the probability of best value method (specified as method="pBV"). If NULL then method="sucra".
#' @param digits Minimal number of significant digits 
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
#'
#' # Adverse events
#'
#' p3 <- pairwise(treat = list(treatment1, treatment2,treatment3),
#'               event = list(ae1, ae2, ae3), 
#'               n = list(n1, n2, n3),
#'               studlab = id,
#'               data = dat.linde2015,
#'               sm = "OR")
#'
#' # Loss to follow-up
#'
#' p4 <- pairwise(treat = list(treatment1, treatment2, treatment3),
#'               event = list(loss1, loss2, loss3), 
#'               n = list(n1, n2, n3),
#'               studlab = id,
#'               data = dat.linde2015,
#'               sm = "OR")
#'
#' ## Loss_to_follow_up_(AE)
#'
#' p5 <- pairwise(treat = list(treatment1, treatment2, treatment3),
#'               event = list(loss.ae1, loss.ae2, loss.ae3),
#'               n = list(n1, n2, n3),
#'               studlab = id,
#'               data = dat.linde2015,
#'               sm = "OR")#
#'
#' # Perform analysis in terms of the Efficacy outcomes
#'
#' p_all <- list(p1,p2,p3,p4,p5)
#'
#' # Use 'mvdata()' to transform the data in suitable JAGS format
#'
#' data_all <- mvdata(p_all)
#'
#' # Define outcome labels
#' 
#' outlab <- c("Early_Response",
#'              "Early_Remission",
#'              "Adverse_events",
#'             "Loss_to_follow_up",
#'             "Loss_to_follow_up_AE")
#'             
#'  
#' # Perform analysis in terms of the all outcomes
#' 
#' 
#' # Fit the model combining all five outcomes
#' 
#' mvmodel_all <- mvnma(data = data_all,
#'                 reference.group = "Placebo",
#'                 outlab = outlab,
#'                 n.iter = 1000,
#'                 n.burnin = 100)
#'
#'
#' # Rank treatments using pBv
#' 
#' ranks_pBV <- mvrank(mvmodel_all,small.values = c("undesirable","undesirable","desirable","desirable","desirable"),
#'                      method = "pBV")
#'  
#' ranks_pBV                    
#'                      
#' # Rank treatments using sucra
#' 
#' ranks_sucra <- mvrank(mvmodel_all,small.values = c("undesirable","undesirable","desirable","desirable","desirable"),
#'                      method = "sucra")
#'                      
#' ranks_sucra
#'                      
#' # same results without method = "sucra" 
#' 
#' ranks_sucra1 <- mvrank(mvmodel_all,small.values = c("undesirable","undesirable","desirable","desirable","desirable"))
#'
#' ranks_sucra1
#' 
#' @export mvrank                              



mvrank <- function(x, small.values,method=NULL,digits=4){


  if(!inherits(x,"mvnma")){
    
   stop("Argument 'x' should be a 'mvnma' object.") 
  
    }
  
  if(is.null(method)){
    
    method = "sucra"  
    
  }
  
  if(method %!in% c("sucra","pBV")){
    
  stop("Argument method should be either 'sucra' or 'pBV'.")  
    
  }
  
x <- x[names(x) != "outcome_correlation"]

outlab <- attributes(x)$names

n.out <- length(outlab)

#### extract samples and create rankograms for each outcome

d <- list()

n.trts <- list()

rank_out <- list()

pBV <- list()

sucra <- list()

ranks <- list()

trts <- list()



for(i in 1:n.out){
  
d[[i]] <- x[[i]]$samples  

n.trts[[i]] <- ncol(d[[i]])

trts[[i]] <- names(d[[i]])

rank_out[[i]] <- rankogram(d[[i]],small.values=small.values[i])

if(method=="pBV"){
  
pBV[[i]] <- as.data.frame(rank_out[[i]]$ranking.matrix.random[,1])  
  
names(pBV[[i]]) <- "pBV"

pBV[[i]]$treatment <- row.names(pBV[[i]])

row.names(pBV[[i]]) <- NULL

pBV[[i]] <- pBV[[i]] %>% 
   dplyr::select(treatment,pBV) %>% 
   arrange(desc(pBV)) 

pBV[[i]]$pBV <- round(pBV[[i]]$pBV,digits = digits)

ranks[[i]] <- pBV[[i]]

}else if(method=="sucra"){
  
sucra[[i]] <- rank_out[[i]]$ranking.random
  
sucra[[i]] <- cbind.data.frame(names(sucra[[i]]),unname(sucra[[i]]))
  
names(sucra[[i]]) <- c("treatment","SUCRA")
  
sucra[[i]] <- sucra[[i]] %>% 
 arrange(desc(SUCRA))

sucra[[i]]$SUCRA <- round(sucra[[i]]$SUCRA,digits = digits)
 
ranks[[i]] <- sucra[[i]]
 
}
  
}

common_trts <- as.data.frame(table(unlist(trts)))

common_trts <- common_trts %>%
  filter(Freq==n.out)

common_trts <- common_trts$Var1

names(ranks) <- outlab

class(ranks) <- "mvrank"

attr(ranks,"common_trts") <- common_trts

return(ranks)

}
