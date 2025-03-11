#' Perform a Bayesian multivariate network meta-analysis using a
#' single-correlation coefficient model
#' 
#' @description
#' This function fits a Bayesian multivariate network meta-analysis model.
#' Currently, the function can simultaneously pool up to five outcomes.
#' Additionally, the studies to be included should be of maximum three arms.
#' 
#' @param data An object of class \code{\link{mvdata}}.
#' @param reference.group A common reference treatment across all outcomes.
#' @param outlab An optional argument with labels for each outcome. If NULL,
#'   the each outcome is labelled as 'outcome_1', 'outcome_2' etc.
#' @param n.iter Number of iterations.
#' @param n.burnin Number of iterations for burnin.
#' @param lower.rho Lower bounds for the Uniform prior(s) used for the correlation
#'   coefficient. If NULL all bounds are set to -1.
#' @param upper.rho Upper bounds for the Uniform prior(s) used for the correlation
#'   coefficient. If NULL all bounds are set to 1.
#' 
#' @details
#' The function \code{\link{mvnma}} expects data in the format of a "mvdata"
#' object. This function transforms the data into a suitable format for JAGS.
#' A common reference treatment across all outcomes is required. However, this
#' requirement is only for enabling the calculations, as the function 
#' \code{\link{league}} generates all possible comparisons.
#' 
#' The Bayesian multivariate network meta-analysis model fitted in the
#' \bold{mvnma} package assumes uniform priors for the between-outcome
#' correlation coefficients. The lower and upper bounds of these priors can be
#' defined using the arguments `lower.rho` and `upper.rho`. If not set, the
#' model will assume a `Unif (-1, 1)` prior for all correlation coefficients.
#' In cases with cases with two outcomes a single value should be provided for
#' `lower.rho` and for `upper.rho`. For example when \(rho_{12} ~ Unif (0.5,1)/)
#' then `lower.rho` = 0.5 and `upper.rho` = 1.
#' In cases with more than two outcomes, the order in which the bounds are
#' provided matters. For example, when pooling four outcomes, the lower and
#' upper bounds correspond to the following order of correlation coefficients:
#' (rho_{12}, rho_{13}, rho_{14}, rho_{23}, rho_{24}, rho_{34}).

#' @return
#' The function return an 'mvnma' object. This consists of the results for each
#' outcome and the correlation coefficient estimates between the combined
#' outcomes. The outcome-specific estimates are expressed in the format of a
#' list (one for each outcome) which contains:
#' \itemize{
#' \item The basic estimates (i.e. treatment vs. reference.group) for each
#'   outcome.
#' \item The heterogeneity estimates for each outcome 
#' \item The posterior samples corresponding to the basic estimates.
#' }
#' 
#' @examples
#' \donttest{
#' library("netmeta")
#' 
#' data("Linde2015")
#' 
#' # Use 'pairwise' to obtain contrast based data for each one of the five
#' # available outcomes 
#'
#' # Early response
#'
#' p1 <- pairwise(treat = list(treatment1, treatment2, treatment3),
#'   event = list(resp1, resp2, resp3), n = list(n1, n2, n3),
#'   studlab = id, data = dat.linde2015, sm = "OR")
#'
#' # Early remissions
#'
#' p2 <- pairwise(treat = list(treatment1, treatment2, treatment3),
#'   event = list(remi1, remi2, remi3), n = list(n1, n2, n3),
#'   studlab = id, data = dat.linde2015, sm = "OR")
#'
#' # Adverse events
#'
#' p3 <- pairwise(treat = list(treatment1, treatment2,treatment3),
#'   event = list(ae1, ae2, ae3),  n = list(n1, n2, n3),
#'   studlab = id, data = dat.linde2015, sm = "OR")
#'
#' # Loss to follow-up
#' 
#' p4 <- pairwise(treat = list(treatment1, treatment2, treatment3),
#'   event = list(loss1, loss2, loss3), n = list(n1, n2, n3),
#'   studlab = id, data = dat.linde2015, sm = "OR")
#'
#' # Loss_to_follow_up_(AE)
#' 
#' p5 <- pairwise(treat = list(treatment1, treatment2, treatment3),
#'   event = list(loss.ae1, loss.ae2, loss.ae3), n = list(n1, n2, n3),
#'   studlab = id, data = dat.linde2015, sm = "OR")
#'
#' # Perform analysis in terms of the efficacy outcomes
#'
#' p_effic <- list(p1, p2)
#'
#' # Use 'mvdata()' to transform the data in suitable JAGS format
#'
#' data_effic <- mvdata(p_effic)
#'
#' # Define outcome labels
#' 
#' outlab <- c("Early_Response", "Early_Remission",
#'   "Adverse_events", "Loss_to_follow_up", "Loss_to_follow_up_AE")
#'  
#' # Fit the model combining only the two efficacy outcomes
#' set.seed(1909)
#' mvmodel_effic <- mvnma(data = data_effic,
#'   reference.group = "Placebo", outlab = outlab[1:2],
#'   n.iter = 1000, n.burnin = 100)
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
#' forest(mvmodel_effic)
#' 
#' # Get all estimates                 
#' 
#' league.effic <- league(mvmodel_effic)
#' 
#' # Perform analysis in terms of the all outcomes
#' 
#' p_all <- list(p1, p2, p3, p4, p5)
#' 
#' data_all <- mvdata(p_all)
#' 
#' # Fit the model combining all five outcomes
#' 
#' mvmodel_all <- mvnma(data = data_all,
#'   reference.group = "Placebo", outlab = outlab,
#'   n.iter = 1000, n.burnin = 100)
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
#' # Plot the results for all outcomes
#' 
#' forest(mvmodel_all)
#'
#' # Get all estimates                 
#' 
#' league.all <- league(mvmodel_all)
#' }
#'
#' @export mvnma

mvnma <- function(data,
                  reference.group = NULL,
                  outlab = NULL,
                  n.iter = 10000,
                  n.burnin = 2000,
                  lower.rho = NULL,
                  upper.rho = NULL
                  ) {
  
  
  chkclass(data, "mvdata")
  #

  if (is.null(reference.group)) {
    
    stop("Argument 'reference.group' is mandatory.")  
    
  }  
  # extract number of outcomes  
  n.out <- ncol(data$var)
  
  ## create bounds for correlation prior
  
  
  if (n.out==2) {
    
    bounds <- c()  
    
    if ( (is.null(lower.rho)) && (is.null(upper.rho)) ) {
      
      # if not information is provided about the prior then set this to unif[-1,1]  
      
      lower.rho = -1
      
      upper.rho = 1 
      
    }else if ( (!is.null(lower.rho)) && (!is.null(upper.rho)) ) {
      
      # ensure proper set of lower and upper bounds for the uniform prior    
      bounds <- c(lower.rho,upper.rho)
      
      lower.rho <- min(bounds)
      
      upper.rho <- max(bounds)
      
    }else if ( (!is.null(lower.rho)) && (is.null(upper.rho)) ) {
      
      upper.rho <- 1  
      
    }else if ( (is.null(lower.rho)) && (!is.null(upper.rho)) ) {
      
      lower.rho <- -1  
      
    }
    
  }else if (n.out==3) {
    
    bounds1 <- c()
    
    bounds2 <- c()
    
    bounds3 <- c()
    
    if ( (is.null(lower.rho)) && (is.null(upper.rho)) ) {
      
      lower.rho1 <- lower.rho2 <- lower.rho3 <- -1
      
      upper.rho1 <- upper.rho2 <- upper.rho3 <- 1
      
    }else if ( !(is.null(lower.rho)) && !(is.null(upper.rho)) ) {
      
      bounds1 <- c(lower.rho[1],upper.rho[1])
      
      bounds2 <- c(lower.rho[2],upper.rho[2])
      
      bounds3 <- c(lower.rho[3],upper.rho[3])
      
      lower.rho1 <- min(bounds1)
      
      lower.rho2 <- min(bounds2)
      
      lower.rho3 <- min(bounds3)
      
      upper.rho1 <- max(bounds1)  
      
      upper.rho2 <- max(bounds2)
      
      upper.rho3 <- max(bounds3)
      
    }else if ( (is.null(lower.rho)) && !(is.null(upper.rho)) ) {
      
      lower.rho1 <- lower.rho2 <- lower.rho3 <- -1
      
    }else if ( !(is.null(lower.rho)) && (is.null(upper.rho)) ) {
      
      upper.rho1 <- upper.rho2 <- upper.rho3 <- 1
    }
    
  }else if (n.out==4) {
    
    bounds1 <- c()
    
    bounds2 <- c()
    
    bounds3 <- c()
    
    bounds4 <- c()
    
    bounds5 <- c()
    
    bounds6 <- c()
    
    if ( (is.null(lower.rho)) && (is.null(upper.rho)) ) {
      
      lower.rho1 <- lower.rho2 <- lower.rho3 <- lower.rho4 <- lower.rho5 <- lower.rho6 <- -1
      
      upper.rho1 <- upper.rho2 <- upper.rho3 <- upper.rho4 <- upper.rho5 <- upper.rho6 <- 1
      
    }else if ( !(is.null(lower.rho)) && !(is.null(upper.rho)) ) {
      
      bounds1 <- c(lower.rho[1],upper.rho[1])
      
      bounds2 <- c(lower.rho[2],upper.rho[2])
      
      bounds3 <- c(lower.rho[3],upper.rho[3])
      
      bounds4 <- c(lower.rho[4],upper.rho[4])
      
      bounds5 <- c(lower.rho[5],upper.rho[5])
      
      bounds6 <- c(lower.rho[6],upper.rho[6])
      
      lower.rho1 <- min(bounds1)
      
      lower.rho2 <- min(bounds2)
      
      lower.rho3 <- min(bounds3)
      
      lower.rho4 <- min(bounds4)
      
      lower.rho5 <- min(bounds5)
      
      lower.rho6 <- min(bounds6)
      
      
      upper.rho1 <- max(bounds1)  
      
      upper.rho2 <- max(bounds2)
      
      upper.rho3 <- max(bounds3)
      
      upper.rho4 <- max(bounds4)
      
      upper.rho5 <- max(bounds5)
      
      upper.rho6 <- max(bounds6)
      
    }else if ( (is.null(lower.rho)) && !(is.null(upper.rho)) ) {
      
      lower.rho1 <- lower.rho2 <- lower.rho3 <- lower.rho4 <- lower.rho5 <- lower.rho6 <- -1
      
    }else if ( !(is.null(lower.rho)) && (is.null(upper.rho)) ) {
      
      upper.rho1 <- upper.rho2 <- upper.rho3 <- upper.rho4 <- upper.rho5 <- upper.rho6 <- 1
    }
    
  }else if (n.out==5) {
    
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
    
    if ( (is.null(lower.rho)) && (is.null(upper.rho)) ) {
      
      lower.rho1 <- lower.rho2 <- lower.rho3 <- lower.rho4 <- lower.rho5 <- lower.rho6 <- lower.rho7 <- lower.rho8 <- lower.rho9 <- lower.rho10 <- -1
      
      upper.rho1 <- upper.rho2 <- upper.rho3 <- upper.rho4 <- upper.rho5 <- upper.rho6 <- upper.rho7 <- upper.rho8 <- upper.rho9 <- upper.rho10 <- 1
      
    }else if ( !(is.null(lower.rho)) && !(is.null(upper.rho)) ) {
      
      bounds1 <- c(lower.rho[1],upper.rho[1])
      
      bounds2 <- c(lower.rho[2],upper.rho[2])
      
      bounds3 <- c(lower.rho[3],upper.rho[3])
      
      bounds4 <- c(lower.rho[4],upper.rho[4])
      
      bounds5 <- c(lower.rho[5],upper.rho[5])
      
      bounds6 <- c(lower.rho[6],upper.rho[6])
      
      bounds7 <- c(lower.rho[7],upper.rho[7])
      
      bounds8 <- c(lower.rho[8],upper.rho[8])
      
      bounds9 <- c(lower.rho[9],upper.rho[9])
      
      bounds10 <- c(lower.rho[10],upper.rho[10])
      
      lower.rho1 <- min(bounds1)
      
      lower.rho2 <- min(bounds2)
      
      lower.rho3 <- min(bounds3)
      
      lower.rho4 <- min(bounds4)
      
      lower.rho5 <- min(bounds5)
      
      lower.rho6 <- min(bounds6)
      
      lower.rho7 <- min(bounds7)
      
      lower.rho8 <- min(bounds8)
      
      lower.rho9 <- min(bounds9)
      
      lower.rho10 <- min(bounds10)
      
      
      upper.rho1 <- max(bounds1)  
      
      upper.rho2 <- max(bounds2)
      
      upper.rho3 <- max(bounds3)
      
      upper.rho4 <- max(bounds4)
      
      upper.rho5 <- max(bounds5)
      
      upper.rho6 <- max(bounds6)
      
      upper.rho7 <- max(bounds7)
      
      upper.rho8 <- max(bounds8)
      
      upper.rho9 <- max(bounds9)
      
      upper.rho10 <- max(bounds10)
      
    }else if ( (is.null(lower.rho)) && !(is.null(upper.rho)) ) {
      
      lower.rho1 <- lower.rho2 <- lower.rho3 <- lower.rho4 <- lower.rho5 <- lower.rho6 <- lower.rho7 <- lower.rho8 <- lower.rho9 <- lower.rho10 <- -1
      
    }else if ( !(is.null(lower.rho)) && (is.null(upper.rho)) ) {
      
      upper.rho1 <- upper.rho2 <- upper.rho3 <- upper.rho4 <- upper.rho5 <- upper.rho6 <- upper.rho7 <- upper.rho8 <- upper.rho9 <- upper.rho10 <- 1
    }
    
  }         
  
  
  
  ## create outcome labels if not provided
  if (is.null(outlab)) {
    
    outlab1 <- rep("outcome",n.out)
    
    outlab2 <- seq(1:n.out)  
    
    outlab <- c()
    
    for(i in 1:n.out) {
      
      outlab[i] <- paste(outlab1[i],outlab2[i],sep = "_")    
      
    }
    
  }else{
    
    outlab <- outlab  
    
    if (length(outlab)!=n.out) {
      
      stop("Please provide a label for all outcomes.")
      
    }
    
  }  
  
  labtreat <- data$labtreat$treat
  
  ref <- unname(which(labtreat==reference.group))  
  
  if (n.out==2) {
    
    if (ncol(data$T)==2) {
      
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
        lower.rho = lower.rho,
        upper.rho = upper.rho
        
      ) 
      
    }else if (ncol(data$T)==3) {
      
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
        lower.rho = lower.rho,
        upper.rho = upper.rho
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
        DIC = FALSE,
        model.file = system.file("model", "mvnma_2.txt", package = "mvnma")
      )
      
    }
    
  } else if (n.out==3) {
    
    
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
      lower.rho1 = lower.rho1,
      lower.rho2 = lower.rho2,
      lower.rho3 = lower.rho3,
      upper.rho1 = upper.rho1,
      upper.rho2 = upper.rho2,
      upper.rho3 = upper.rho3
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
      DIC = FALSE,
      model.file = system.file("model", "mvnma_3.txt", package = "mvnma")
    )  
    
  }else if (n.out==4) {
    
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
      lower.rho1 = lower.rho1,
      lower.rho2 = lower.rho2,
      lower.rho3 = lower.rho3,
      lower.rho4 = lower.rho4,
      lower.rho5 = lower.rho5,
      lower.rho6 = lower.rho6,
      
      upper.rho1 = upper.rho1,
      upper.rho2 = upper.rho2,
      upper.rho3 = upper.rho3,
      upper.rho4 = upper.rho4,
      upper.rho5 = upper.rho5,
      upper.rho6 = upper.rho6
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
      DIC = FALSE,
      model.file = system.file("model", "mvnma_4.txt", package = "mvnma")
    )  
  }else if (n.out==5) {
    
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
      lower.rho1 = lower.rho1,
      lower.rho2 = lower.rho2,
      lower.rho3 = lower.rho3,
      lower.rho4 = lower.rho4,
      lower.rho5 = lower.rho5,
      lower.rho6 = lower.rho6,
      lower.rho7 = lower.rho7,
      lower.rho8 = lower.rho8,
      lower.rho9 = lower.rho9,
      lower.rho10 = lower.rho10,
      
      upper.rho1 = upper.rho1,
      upper.rho2 = upper.rho2,
      upper.rho3 = upper.rho3,
      upper.rho4 = upper.rho4,
      upper.rho5 = upper.rho5,
      upper.rho6 = upper.rho6,
      upper.rho7 = upper.rho7,
      upper.rho8 = upper.rho8,
      upper.rho9 = upper.rho9,
      upper.rho10 = upper.rho10
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
      DIC = FALSE,
      model.file = system.file("model", "mvnma_5.txt", package = "mvnma")
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
