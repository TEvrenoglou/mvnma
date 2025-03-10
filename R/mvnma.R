#' Perform a Bayesian multivariate network meta-analysis using a single-correlation coefficient model.
#' 
#' @description
#' This function fits a Bayesian multivariate network meta-analysis model. Currently, the function can 
#' simultaneously pool up to 5 outcomes. Additionally, the studies to be included should be of maximum 3 arms.
#' 
#' @param data An object of class \code{\link{mvdata}}.
#' @param reference.group A common reference treatment across all outcomes.
#' @param outlab An optional argument with labels for each outcome. If NULL, the each outcome is labelled as 'outcome_1', 'outcome_2' etc.
#' @param n.iter Number of iterations.
#' @param n.burnin Number of iterations for burnin.
#' @param lb.rho Lower bounds for the Uniform prior(s) used for the correlation coefficient. If NULL all bounds are set to -1.
#' @param ub.rho Upper bounds for the Uniform prior(s) used for the correlation coefficient. If NULL all bounds are set to 1.

#' @details
#' The function \code{\link{mvnma}} expects data in the format of a "mvdata" object. This function transforms the data into a suitable format for JAGS.
#' A common reference treatment across all outcomes is required. However, this requirement is only for enabling the calculations, as the function 
#' \code{\link{league}} generates all possible comparisons.
#' 
#' The Bayesian multivariate network meta-analysis model fitted in the \bold{mvnma}
#' package assumes uniform priors for the between-outcome correlation coefficients. The lower and upper bounds of these priors can be defined using
#' the arguments `lb.rho` and `ub.rho`. If not set, the model will assume a `Unif(-1, 1)` prior for all correlation coefficients. In cases with cases with two outcomes a single value should be provided for `lb.rho` and for `ub.rho`. For example when \(\rho_{12}~Unif(0.5,1)/) then `lb.rho`=0.5 and `ub.rho`=1.    
#' In cases with more than two outcomes, the order in which the bounds are provided matters. For example, when pooling 4 outcomes, the lower and upper bounds correspond to the following order of correlation coefficients: \(\rho_{12}, \rho_{13}, \rho_{14}, \rho_{23}, \rho_{24}, \rho_{34}\).

#' @return
#' The function return a 'mvnma' object. This consists of the results for each outcome and the correlation coefficient estimates between the combined outcomes. The outcome-specific estimates are expressed in the format of a list (one for each outcome) which contains:
#' \begin{itemize}
#' \item The basic estimates (i.e. treatment vs. reference.group) for each outcome.
#' \item The heterogeneity estimates for each outcome 
#' \item The posterior samples corresponding to the basic estimates.
#' \end{itemize}
#' 
#' @examples
#'
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
#'
#' ## Adverse events
#'
#' p3 <- pairwise(treat = list(treatment1, treatment2,treatment3),
#'               event = list(ae1, ae2, ae3), 
#'               n = list(n1, n2, n3),
#'               studlab = id,
#'               data = dat.linde2015,
#'               sm = "OR")

#' ## Loss to follow-up

#'p4 <- pairwise(treat = list(treatment1, treatment2, treatment3),
#'               event = list(loss1, loss2, loss3), 
#'               n = list(n1, n2, n3),
#'               studlab = id,
#'               data = dat.linde2015,
#'               sm = "OR")

#' ## Loss_to_follow_up_(AE)

#'p5 <- pairwise(treat = list(treatment1, treatment2, treatment3),
#'               event = list(loss.ae1, loss.ae2, loss.ae3),
#'               n = list(n1, n2, n3),
#'               studlab = id,
#'               data = dat.linde2015,
#'               sm = "OR")#
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
#'              "Early_Remission",
#'              "Adverse_events",
#'             "Loss_to_follow_up",
#'             "Loss_to_follow_up_AE")
#'             
#' # Fit the model combining only the two efficacy outcomes
#' 
#' mvmodel_effic <- mvnma(data = data_effic,
#'                 reference.group = "Placebo",
#'                 outlab = outlab[1:2],
#'                 n.iter = 1000,
#'                 n.burnin = 100)
#'                 
#' # Extract treatment effect estimates and heterogeneity for Early_Response 
#' 
#' mvmodel_effic$Early_Response$basic_estimates
#' 
#' mvmodel_effic$Early_Response$heterogeneity      
#'      
#' # Extract outcome correlation 
#' 
#' mvmodel_effic$outcome_correlation
#'              
#' # Plot the results for efficacy outcomes
#' 
#' forest.mvnma(mvmodel_effic)
#' 
#' # Get all estimates                 
#' 
#' league.effic <- league(mvmodel_effic)
#' 
#' # Perform analysis in terms of the all outcomes
#' 
#' p_all <- list(p1,p2,p3,p4,p5)
#' 
#' data_all <- mvdata(p_all)
#' 
#' # Fit the model combining all five outcomes
#' 
#' mvmodel_all <- mvnma(data = data_all,
#'                 reference.group = "Placebo",
#'                 outlab = outlab,
#'                 n.iter = 1000,
#'                 n.burnin = 100)
#'                 
#' # Extract treatment effect estimates and heterogeneity for Early_Response 
#' 
#' mvmodel_all$Early_Response$basic_estimates
#' 
#' mvmodel_all$Early_Response$heterogeneity      
#'      
#' # Extract outcome correlation 
#' 
#' mvmodel_all$outcome_correlation
#'
#'
#' # Plot the results for all outcomes
#' 
#' forest.mvnma(mvmodel_all)
#'
#' # Get all estimates                 
#' 
#' league.all <- league(mvmodel_all)
#'
#' @export mvnma

mvnma <- function(data = NULL,
                  reference.group=NULL,
                  outlab = NULL,
                  n.iter = 10000,
                  n.burnin = 2000,
                  lb.rho = NULL,
                  ub.rho = NULL
                  ) {
  
  
if(is.null(data)){
  
stop("Argument 'data' is mandatory.")  
  
}else if((!is.null(data)) & class(data)!="mvdata"){
  
 stop("Argument 'data' should be a 'mvdata' object.") 
  
}  
  
if(is.null(reference.group)){
  
stop("Argument 'reference.group' is mandatory.")  
  
}  
# extract number of outcomes  
n.out <- ncol(data$var)

## create bounds for correlation prior


if(n.out==2){
  
bounds <- c()  

if( (is.null(lb.rho)) && (is.null(ub.rho)) ){

# if not information is provided about the prior then set this to unif[-1,1]  

lb.rho = -1

ub.rho = 1 
  
}else if( (!is.null(lb.rho)) && (!is.null(ub.rho)) ){

# ensure proper set of lower and upper bounds for the uniform prior    
bounds <- c(lb.rho,ub.rho)

lb.rho <- min(bounds)

ub.rho <- max(bounds)
  
}else if( (!is.null(lb.rho)) && (is.null(ub.rho)) ){
  
ub.rho <- 1  
  
}else if( (is.null(lb.rho)) && (!is.null(ub.rho)) ){
  
lb.rho <- -1  
  
}

}else if(n.out==3){
  
bounds1 <- c()

bounds2 <- c()

bounds3 <- c()

if( (is.null(lb.rho)) && (is.null(ub.rho)) ){
 
lb.rho1 <- lb.rho2 <- lb.rho3 <- -1

ub.rho1 <- ub.rho2 <- ub.rho3 <- 1
   
}else if( !(is.null(lb.rho)) && !(is.null(ub.rho)) ){
  
bounds1 <- c(lb.rho[1],ub.rho[1])

bounds2 <- c(lb.rho[2],ub.rho[2])

bounds3 <- c(lb.rho[3],ub.rho[3])

lb.rho1 <- min(bounds1)

lb.rho2 <- min(bounds2)

lb.rho3 <- min(bounds3)

ub.rho1 <- max(bounds1)  

ub.rho2 <- max(bounds2)

ub.rho3 <- max(bounds3)

}else if( (is.null(lb.rho)) && !(is.null(ub.rho)) ){
  
lb.rho1 <- lb.rho2 <- lb.rho3 <- -1

  }else if( !(is.null(lb.rho)) && (is.null(ub.rho)) ){
  
    ub.rho1 <- ub.rho2 <- ub.rho3 <- 1
  }

}else if(n.out==4){
  
  bounds1 <- c()
  
  bounds2 <- c()
  
  bounds3 <- c()
  
  bounds4 <- c()
  
  bounds5 <- c()
  
  bounds6 <- c()
  
  if( (is.null(lb.rho)) && (is.null(ub.rho)) ){
    
    lb.rho1 <- lb.rho2 <- lb.rho3 <- lb.rho4 <- lb.rho5 <- lb.rho6 <- -1
    
    ub.rho1 <- ub.rho2 <- ub.rho3 <- ub.rho4 <- ub.rho5 <- ub.rho6 <- 1
    
  }else if( !(is.null(lb.rho)) && !(is.null(ub.rho)) ){
    
    bounds1 <- c(lb.rho[1],ub.rho[1])
    
    bounds2 <- c(lb.rho[2],ub.rho[2])
    
    bounds3 <- c(lb.rho[3],ub.rho[3])
    
    bounds4 <- c(lb.rho[4],ub.rho[4])
    
    bounds5 <- c(lb.rho[5],ub.rho[5])
    
    bounds6 <- c(lb.rho[6],ub.rho[6])
    
    lb.rho1 <- min(bounds1)
    
    lb.rho2 <- min(bounds2)
    
    lb.rho3 <- min(bounds3)
    
    lb.rho4 <- min(bounds4)
    
    lb.rho5 <- min(bounds5)
    
    lb.rho6 <- min(bounds6)
    
    
    ub.rho1 <- max(bounds1)  
    
    ub.rho2 <- max(bounds2)
    
    ub.rho3 <- max(bounds3)
    
    ub.rho4 <- max(bounds4)
    
    ub.rho5 <- max(bounds5)
    
    ub.rho6 <- max(bounds6)
    
  }else if( (is.null(lb.rho)) && !(is.null(ub.rho)) ){
    
    lb.rho1 <- lb.rho2 <- lb.rho3 <- lb.rho4 <- lb.rho5 <- lb.rho6 <- -1
    
  }else if( !(is.null(lb.rho)) && (is.null(ub.rho)) ){
    
    ub.rho1 <- ub.rho2 <- ub.rho3 <- ub.rho4 <- ub.rho5 <- ub.rho6 <- 1
  }
  
}else if(n.out==5){
  
  bounds1 <- c()
  
  bounds2 <- c()
  
  bounds3 <- c()
  
  bounds4 <- c()
  
  bounds5 <- c()
  
  bounds5 <- c()
  
  bounds7 <- c()
  
  bounds8 <- c()
  
  bounds9 <- c()
  
  bounds10 <- c()
  
  if( (is.null(lb.rho)) && (is.null(ub.rho)) ){
    
    lb.rho1 <- lb.rho2 <- lb.rho3 <- lb.rho4 <- lb.rho5 <- lb.rho6 <- lb.rho7 <- lb.rho8 <- lb.rho9 <- lb.rho10 <- -1
    
    ub.rho1 <- ub.rho2 <- ub.rho3 <- ub.rho4 <- ub.rho5 <- ub.rho6 <- ub.rho7 <- ub.rho8 <- ub.rho9 <- ub.rho10 <- 1
    
  }else if( !(is.null(lb.rho)) && !(is.null(ub.rho)) ){
    
    bounds1 <- c(lb.rho[1],ub.rho[1])
    
    bounds2 <- c(lb.rho[2],ub.rho[2])
    
    bounds3 <- c(lb.rho[3],ub.rho[3])
    
    bounds4 <- c(lb.rho[4],ub.rho[4])
    
    bounds5 <- c(lb.rho[5],ub.rho[5])
    
    bounds6 <- c(lb.rho[6],ub.rho[6])
    
    bounds7 <- c(lb.rho[7],ub.rho[7])
    
    bounds8 <- c(lb.rho[8],ub.rho[8])
    
    bounds9 <- c(lb.rho[9],ub.rho[9])
    
    bounds10 <- c(lb.rho[10],ub.rho[10])
    
    lb.rho1 <- min(bounds1)
    
    lb.rho2 <- min(bounds2)
    
    lb.rho3 <- min(bounds3)
    
    lb.rho4 <- min(bounds4)
    
    lb.rho5 <- min(bounds5)
    
    lb.rho6 <- min(bounds6)
    
    lb.rho7 <- min(bounds7)
    
    lb.rho8 <- min(bounds8)
    
    lb.rho9 <- min(bounds9)
    
    lb.rho10 <- min(bounds10)
    
    
    ub.rho1 <- max(bounds1)  
    
    ub.rho2 <- max(bounds2)
    
    ub.rho3 <- max(bounds3)
    
    ub.rho4 <- max(bounds4)
    
    ub.rho5 <- max(bounds5)
    
    ub.rho6 <- max(bounds6)
    
    ub.rho7 <- max(bounds7)
    
    ub.rho8 <- max(bounds8)
    
    ub.rho9 <- max(bounds9)
    
    ub.rho10 <- max(bounds10)
    
  }else if( (is.null(lb.rho)) && !(is.null(ub.rho)) ){
    
    lb.rho1 <- lb.rho2 <- lb.rho3 <- lb.rho4 <- lb.rho5 <- lb.rho6 <- lb.rho7 <- lb.rho8 <- lb.rho9 <- lb.rho10 <- -1
    
  }else if( !(is.null(lb.rho)) && (is.null(ub.rho)) ){
    
    ub.rho1 <- ub.rho2 <- ub.rho3 <- ub.rho4 <- ub.rho5 <- ub.rho6 <- ub.rho7 <- ub.rho8 <- ub.rho9 <- ub.rho10 <- 1
  }
  
}         


  
## create outcome labels if not provided
if(is.null(outlab)){
  
outlab1 <- rep("outcome",n.out)
  
outlab2 <- seq(1:n.out)  

outlab <- c()

for(i in 1:n.out){
  
outlab[i] <- paste(outlab1[i],outlab2[i],sep = "_")    

}

}else{

  outlab <- outlab  
  
  if(length(outlab)!=n.out){
    
  stop("Please provide a label for all outcomes.")
    
  }
  
}  
  
labtreat <- data$labtreat$treat

ref <- unname(which(labtreat==reference.group))  

if(n.out==2){
  
if(ncol(data$T)==2){
  
  run.data <- list(
    y = data$y,
    var1 = data$var$var1,
    var2 = data$var$var2,
    ref = ref,
    Ns = data$Ns,
    N2h = data$N2h,
    T1 = data$T[,1],
    T2 = data$T[,2],
    NT = data$NT,
    lb.rho = lb.rho,
    ub.rho = ub.rho
    
  ) 
  
}else if(ncol(data$T)==3){
   
   run.data <- list(
    y = data$y,
    var1 = data$var$var1,
    var2 = data$var$var2,
    ref = ref,
    Ns = data$Ns,
    N2h = data$N2h,
    T1 = data$T[,1],
    T2 = data$T[,2],
    T3 = data$T[,3],
    NT = data$NT,
    lb.rho = lb.rho,
    ub.rho = ub.rho
  ) 
   
   run = jags(
     data = run.data,
     inits = NULL,
     parameters.to.save = c(
       "res.ref1",
       "res.ref2",
       "rho1",
       #"res.1",
       #"res.2",
       "psi1",
       "psi2",
       "d1",
       "d2"
     ),
     n.chains = 2,
     n.iter = n.iter,
     n.burnin = n.burnin,
     DIC = F,
     model.file = modfile_2out
   )

}

  } else if(n.out==3){
  

  run.data <- list(
    y = data$y,
    var1 = data$var$var1,
    var2 = data$var$var2,
    var3 = data$var$var3,
    ref = ref,
    Ns = data$Ns,
    N2h = data$N2h,
    T1 = data$T[,1],
    T2 = data$T[,2],
    T3 = data$T[,3],
    NT = data$NT,
    lb.rho1 = lb.rho1,
    lb.rho2 = lb.rho2,
    lb.rho3 = lb.rho3,
    ub.rho1 = ub.rho1,
    ub.rho2 = ub.rho2,
    ub.rho3 = ub.rho3
  ) 
  
  
  run = jags(
    data = run.data,
    inits = NULL,
    parameters.to.save = c(
      "res.ref1",
      "res.ref2",
      "res.ref3",
      "rho12",
      "rho13",
      "rho23",
      "psi1",
      "psi2",
      "psi3",
      "d1",
      "d2",
      "d3"
    ),
    n.chains = 2,
    n.iter = n.iter,
    n.burnin = n.burnin,
    DIC = F,
    model.file = modfile_3out
  )  
    
  }else if(n.out==4){
    
    run.data <- list(
      y = data$y,
      var1 = data$var$var1,
      var2 = data$var$var2,
      var3 = data$var$var3,
      var4 = data$var$var4,
      ref = ref,
      Ns = data$Ns,
      N2h = data$N2h,
      T1 = data$T[,1],
      T2 = data$T[,2],
      T3 = data$T[,3],
      NT = data$NT,
      lb.rho1 = lb.rho1,
      lb.rho2 = lb.rho2,
      lb.rho3 = lb.rho3,
      lb.rho4 = lb.rho4,
      lb.rho5 = lb.rho5,
      lb.rho6 = lb.rho6,
      
      ub.rho1 = ub.rho1,
      ub.rho2 = ub.rho2,
      ub.rho3 = ub.rho3,
      ub.rho4 = ub.rho4,
      ub.rho5 = ub.rho5,
      ub.rho6 = ub.rho6
    ) 
    
    
    run = jags(
      data = run.data,
      inits = NULL,
      parameters.to.save = c(
        "res.ref1",
        "res.ref2",
        "res.ref3",
        "res.ref4",
        "rho12",
        "rho13",
        "rho14",
        "rho23",
        "rho24",
        "rho34",
        "psi1",
        "psi2",
        "psi3",
        "psi4",
        "d1",
        "d2",
        "d3",
        "d4"
        
      ),
      n.chains = 2,
      n.iter = n.iter,
      n.burnin = n.burnin,
      DIC = F,
      model.file = modfile_4out
    )  
  }else if(n.out==5){
    
    run.data <- list(
      y = data$y,
      var1 = data$var$var1,
      var2 = data$var$var2,
      var3 = data$var$var3,
      var4 = data$var$var4,
      var5 = data$var$var5,
      ref = ref,
      Ns = data$Ns,
      N2h = data$N2h,
      T1 = data$T[,1],
      T2 = data$T[,2],
      T3 = data$T[,3],
      NT = data$NT,
      lb.rho1 = lb.rho1,
      lb.rho2 = lb.rho2,
      lb.rho3 = lb.rho3,
      lb.rho4 = lb.rho4,
      lb.rho5 = lb.rho5,
      lb.rho6 = lb.rho6,
      lb.rho7 = lb.rho7,
      lb.rho8 = lb.rho8,
      lb.rho9 = lb.rho9,
      lb.rho10 = lb.rho10,
      
      ub.rho1 = ub.rho1,
      ub.rho2 = ub.rho2,
      ub.rho3 = ub.rho3,
      ub.rho4 = ub.rho4,
      ub.rho5 = ub.rho5,
      ub.rho6 = ub.rho6,
      ub.rho7 = ub.rho7,
      ub.rho8 = ub.rho8,
      ub.rho9 = ub.rho9,
      ub.rho10 = ub.rho10
    ) 
    
    
    run = jags(
      data = run.data,
      inits = NULL,
      parameters.to.save = c(
        "res.ref1",
        "res.ref2",
        "res.ref3",
        "res.ref4",
        "res.ref5",
        "rho12",
        "rho13",
        "rho14",
        "rho15",
        "rho23",
        "rho24",
        "rho25",
        "rho34",
        "rho35",
        "rho45",
        "psi1",
        "psi2",
        "psi3",
        "psi4",
        "psi5",
        "d1",
        "d2",
        "d3",
        "d4",
        "d5"
        
      ),
      n.chains = 2,
      n.iter = n.iter,
      n.burnin = n.burnin,
      DIC = F,
      model.file = modfile_5out
    )  
  }        


results = as.data.frame(run$BUGSoutput$summary)

results$ind = as.character(rownames(results))

sims_bugs_out <- as.data.frame(run$BUGSoutput$sims.matrix)

####### Manipulate the results and create suitable datasets

res <- gather_results(results,
                      n.out = n.out,
                      labtreat = labtreat,
                      reference.group = reference.group,
                      ref = ref,
                      data = data,
                      sims_bugs_out = sims_bugs_out,
                      outlab = outlab
                      )

all_res <- get_all_estimates(res)

attr(res,"outlab") <- outlab

sm <- attributes(data)$sm

attr(res, "sm") <- sm

attr(res, "all_res") <- all_res

class(res) <- "mvnma"

return(res)

}
