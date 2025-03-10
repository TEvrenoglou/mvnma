#' Forest plot for multivariate network meta-analysis
#' 
#' @description
#' Draws a forest plot in the active graphics window (using grid graphics system).
#' 
#' @param x An object of class \code{\link{mvnma}}.
#' @param ... Additional arguments for \code{\link{forest.meta}} function.
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
#' # Generate a forest plot with the results 
#' 
#' forest.mvnma(mvmodel_effic)                 
#'                 
#' @export forest.mvnma                   

forest.mvnma <- function(x, treat = NULL, backtransf = FALSE,
                       #
                       leftcols = "studlab", leftlabs = "Comparison",
                       rightcols = c("effect", "ci"),
                       rightlabs = c("TE", NA),
                       col.study = "black",
                       col.square = "black",
                       col.square.lines = "black",
                       squaresize = 0.7,
                       header.line = TRUE, 
                       col.subgroup = "black",
                       ...) {
  
  if(!inherits(x,"mvnma")){
    stop("Argument x must be of class 'mvnma'.")
  }
  
  x <- x[names(x) != "outcome_correlation"]
  
  n.out <- length(x)
  
  sm <- attributes(x)$sm
  
  ests <- list()
  
  ## get estimates for each outcome
  for(i in 1:n.out){
  
    ests[[i]] <- x[[i]]$basic_estimates
    
    ests[[i]]$treat <- row.names(ests[[i]])
    
    row.names(ests[[i]]) <- NULL
    
    ests[[i]] <- ests[[i]] %>% 
      dplyr::select(treat,TE,sd,lb.ci,ub.ci)
    
    ests[[i]]$outcome <- attributes(x)$names[i]
  
    
  }
  
 
   # ests_1 <- x[[1]]$basic_estimates
   # 
   # ests_1 <- ests_1 %>% 
   #   mutate("treat" = row.names(ests_1)) %>% 
   #   dplyr::select(treat,mean,sd,`2.5%`,`97.5%`)
   #   
   # ests_1$outcome <- attributes(x)$outlab[1]
   # 
   # ests_2 <- x[[2]]$basic_estimates
   # 
   # ests_2 <- ests_2 %>% 
   #   mutate("treat" = row.names(ests_2)) %>% 
   #   dplyr::select(treat,mean,sd,`2.5%`,`97.5%`)
   # 
   # ests_2$outcome <- attributes(x)$outlab[2]
   
   
### construct final dataset    '
   
  dat <- bind_rows(ests)
  
  row.names(dat) <- NULL
   
  names(dat) <- c("studlab","TE","seTE","lb.ci","ub.ci","outcome")
  
  m <- metagen(dat$TE, dat$seTE, sm = NULL,
               subgroup = dat$outcome,
               backtransf = backtransf,
               print.subgroup.name = FALSE,
               studlab = dat$studlab,
               common = FALSE, random = FALSE, hetstat = FALSE,
               method.tau = "DL", method.tau.ci = "")
  #
  # if (is.null(leftlabs)) {
  #   if (is.null(m$subgroup))
  #     "Study"
  #   else
  #     "Comparison / Study"
  # }
  #
  forest(m,
         backtransf = FALSE,
         header.line = header.line,
         col.subgroup = "black",
         leftcols = leftcols, 
         leftlabs = leftlabs,
         rightcols = rightcols,
         rightlabs = rightlabs,
         weight.study = "same",
         col.study = col.study,
         col.square = col.square,
         col.square.lines = col.square.lines,
         squaresize = squaresize,
         #
         calcwidth.subgroup = TRUE,
         ...
         )
  
  invisible(NULL)
}
