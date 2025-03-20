#' Perform a Bayesian multivariate network meta-analysis using a
#' single-correlation coefficient model
#' 
#' @description
#' This function fits a Bayesian multivariate network meta-analysis model.
#' Currently, the function can simultaneously pool up to five outcomes.
#' Additionally, the studies to be included should be of maximum three arms.
#' 
#' @param \dots Either two to five pairwise objects or a single list with
#'   two to five pairwise objects.
#' @param reference.group A common reference treatment across all outcomes.
#' @param outclab An optional argument with labels for each outcome. If NULL,
#'   the each outcome is labelled as 'outcome_1', 'outcome_2' etc.
#' @param n.iter Number of iterations.
#' @param n.burnin Number of iterations for burn-in.
#' @param level The level used to calculate confidence intervals
#'   for network estimates.
#' @param lower.rho Lower bounds for the Uniform prior(s) used for the
#'   correlation coefficient. If NULL all bounds are set to -1.
#' @param upper.rho Upper bounds for the Uniform prior(s) used for the
#'   correlation coefficient. If NULL all bounds are set to 1.
#' @param quiet A logical indicating whether to print information on the
#'   progress of the JAGS model fitting.
#' 
#' @details
#' The function \code{\link{mvnma}} expects two to five outcomes /
#' \code{\link[meta]{pairwise}} objects. A common reference treatment across
#' all outcomes is required. However, this requirement is only for enabling
#' the calculations, as the function  \code{\link{league}} extracts all
#' possible comparisons.
#' 
#' The Bayesian multivariate network meta-analysis model fitted in the
#' \bold{mvnma} package assumes uniform priors for the between-outcome
#' correlation coefficients. The lower and upper bounds of these priors can be
#' defined using the arguments `lower.rho` and `upper.rho`. If not set, the
#' model will assume a `Unif (-1, 1)` prior for all correlation coefficients.
#' For two outcomes, a single value can be provided for `lower.rho` and
#' `upper.rho`. For example, `lower.rho` = 0.5 and `upper.rho` = 1 for
#' rho12 ~ Unif (0.5,1)).
#' For more than two outcomes, the order in which the bounds are provided
#' matters. For example, when pooling four outcomes, the lower and
#' upper bounds correspond to the following order of correlation coefficients:
#' (rho12, rho13, rho14, rho23, rho24, rho34).

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
#' @seealso \code{\link[meta]{pairwise}}
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
#' p1 <- pairwise(treat = list(treatment1, treatment2, treatment3),
#'   event = list(resp1, resp2, resp3), n = list(n1, n2, n3),
#'   studlab = id, data = dat.linde2015, sm = "OR")
#'
#' # Early remissions
#' p2 <- pairwise(treat = list(treatment1, treatment2, treatment3),
#'   event = list(remi1, remi2, remi3), n = list(n1, n2, n3),
#'   studlab = id, data = dat.linde2015, sm = "OR")
#'
#' # Adverse events
#' p3 <- pairwise(treat = list(treatment1, treatment2,treatment3),
#'   event = list(ae1, ae2, ae3),  n = list(n1, n2, n3),
#'   studlab = id, data = dat.linde2015, sm = "OR")
#'
#' # Loss to follow-up
#' p4 <- pairwise(treat = list(treatment1, treatment2, treatment3),
#'   event = list(loss1, loss2, loss3), n = list(n1, n2, n3),
#'   studlab = id, data = dat.linde2015, sm = "OR")
#'
#' # Loss_to_follow_up_(AE)
#' p5 <- pairwise(treat = list(treatment1, treatment2, treatment3),
#'   event = list(loss.ae1, loss.ae2, loss.ae3), n = list(n1, n2, n3),
#'   studlab = id, data = dat.linde2015, sm = "OR")
#'
#' # Define outcome labels
#' outcomes <- c("Early_Response", "Early_Remission",
#'   "Adverse_events", "Loss_to_follow_up", "Loss_to_follow_up_AE")
#'  
#' # Fit the model combining only the two efficacy outcomes
#' set.seed(1909)
#' mvnma12 <- mvnma(p1, p2,
#'   reference.group = "Placebo", outclab = outcomes[1:2],
#'   n.iter = 1000, n.burnin = 100)
#' mvnma12
#'        
#' # Extract treatment effect estimates and heterogeneity for Early_Response 
#' mvnma12$Early_Response$basic_estimates
#' mvnma12$Early_Response$heterogeneity
#' 
#' # Extract outcome correlation
#' mvnma12$cor
#' 
#' # Plot the results for efficacy outcomes
#' forest(mvnma12)
#' 
#' # Get all estimates
#' league12 <- league(mvnma12)
#' 
#' # Fit the model combining all five outcomes
#' mvnma_all <- mvnma(p1, p2, p3, p4, p5,
#'   reference.group = "Placebo", outclab = outcomes,
#'   n.iter = 1000, n.burnin = 100)
#' 
#' # Extract treatment effect estimates and heterogeneity for Early_Response 
#' mvnma_all$Early_Response$basic_estimates
#' mvnma_all$Early_Response$heterogeneity      
#' 
#' # Extract outcome correlation 
#' mvnma_all$cor
#'
#' # Plot the results for all outcomes
#' forest(mvnma_all)
#'
#' # Get all estimates 
#' league_all <- league(mvnma_all)
#' }
#'
#' @export mvnma

mvnma <- function(...,
                  reference.group = NULL, outclab = NULL,
                  n.iter = 10000, n.burnin = 2000,
                  level = gs("level.ma"),
                  lower.rho, upper.rho,
                  quiet = FALSE) {
  
  is_pairwise <- function(x)
    inherits(x, "pairwise")
  #
  args <- list(...)
  #
  n.out <- length(args)
  n.i <- seq_len(n.out)
  #
  if (n.out == 1) {
    if (is_pairwise(args[[1]]))
      stop("Provide between two and five pairwise objects.",
           call. = FALSE)
    #
    if (!is.list(args[[1]]))
      stop("All elements of argument '...' must be of classes ",
           "'netmeta', 'netcomb', or 'discomb'.",
           call. = FALSE)
    #
    if (!is_pairwise(args[[1]])) {
      n.out <- length(args[[1]])
      n.i <- seq_len(n.out)
      #
      args2 <- list()
      for (i in n.i)
        args2[[i]] <- args[[1]][[i]]
    }
    args <- args2
  }
  #  
  for (i in n.i) {
    if (!is_pairwise(args[[i]]))
      stop("All elements of argument '...' must be of class ",
           "'pairwise'.",
           call. = FALSE)
  }
  #
  if (n.out < 2 | n.out > 5)
    stop("Provide between two and five pairwise objects.",
         call. = FALSE)
  #
  data <- mvdata(args)
  
  treat_out <- data$treat_out
  
  #
  chknull(reference.group)
  chklevel(level)
  chklogical(quiet)
  # extract number of outcomes  
  n.out <- ncol(data$var)
  n.cor <- c(0, 1, 3, 6, 10)[n.out]
  #
  miss.lower <- missing(lower.rho)
  miss.upper <- missing(upper.rho)
  #
  if (!miss.lower)
    chknumeric(lower.rho, min = -1, max = 1, length = n.cor, NA.ok = FALSE)
  #
  if (!miss.upper)
    chknumeric(upper.rho, min = -1, max = 1, length = n.cor, NA.ok = FALSE)
  #
  if (!miss.lower & !miss.upper) {
    if (any(lower.rho >= upper.rho))
      stop("Values for argument 'lower.rho' must be smaller than values for ",
           "argument 'upper.rho'.",
           call. = FALSE)
  }
    
  # Create bounds for correlation prior
  #
  if (miss.lower)
    lower.rho1 <- -1
  else
    lower.rho1 <- lower.rho[1]
  #
  if (miss.upper)
    upper.rho1 <- 1
  else
    upper.rho1 <- upper.rho[1]
  #
  if (n.out >= 3) {
    if (miss.lower) {
      lower.rho2 <- -1
      lower.rho3 <- -1
    }
    else {
      lower.rho2 <- lower.rho[2]
      lower.rho3 <- lower.rho[3]
    }
    #
    if (miss.upper) {
      upper.rho2 <- 1
      upper.rho3 <- 1
    }
    else {
      upper.rho2 <- upper.rho[2]
      upper.rho3 <- upper.rho[3]
    }
  }
  #
  if (n.out >= 4) {
    if (miss.lower) {
      lower.rho4 <- -1
      lower.rho5 <- -1
      lower.rho6 <- -1
    }
    else {
      lower.rho4 <- lower.rho[4]
      lower.rho5 <- lower.rho[5]
      lower.rho6 <- lower.rho[6]
    }
    #
    if (miss.upper) {
      upper.rho4 <- 1
      upper.rho5 <- 1
      upper.rho6 <- 1
    }
    else {
      upper.rho4 <- upper.rho[4]
      upper.rho5 <- upper.rho[5]
      upper.rho6 <- upper.rho[6]
    }
  }
  #
  if (n.out >= 5) {
    if (miss.lower) {
      lower.rho7  <- -1
      lower.rho8  <- -1
      lower.rho9  <- -1
      lower.rho10 <- -1
    }
    else {
      lower.rho7  <- lower.rho[7]
      lower.rho8  <- lower.rho[8]
      lower.rho9  <- lower.rho[9]
      lower.rho10 <- lower.rho[10]
    }
    #
    if (miss.upper) {
      upper.rho7  <- 1
      upper.rho8  <- 1
      upper.rho9  <- 1
      upper.rho10 <- 1
    }
    else {
      upper.rho7  <- upper.rho[7]
      upper.rho8  <- upper.rho[8]
      upper.rho9  <- upper.rho[9]
      upper.rho10 <- upper.rho[10]
    }
  }
  
  
  # Create outcome labels if not provided
  #
  if (is.null(outclab))
    outclab <- paste("outcome", seq_len(n.out), sep = "_")  
  else if (length(outclab) != n.out)
    stop("Please provide labels for all outcomes.")
  
  trts <- data$labtreat$treat
  #
  ref <- unname(which(trts == reference.group))  
  
  
  multiarm <- ncol(data$T) > 2
  #
  run.data <- list(
    y = data$y,
    #
    var1 = data$var$var1,
    var2 = data$var$var2,
    var3 = NA,
    var4 = NA,
    var5 = NA,
    #
    ref = ref,
    #
    k = data$Ns,
    k2 = data$N2h,
    n = data$NT,
    #
    treat1 = data$T[, 1], treat2 = data$T[, 2], treat3 = NA,
    #
    lower.rho1 = lower.rho1, upper.rho1 = upper.rho1,
    lower.rho2 = NA, upper.rho2 = NA,
    lower.rho3 = NA, upper.rho3 = NA,
    lower.rho4 = NA, upper.rho4 = NA,
    lower.rho5 = NA, upper.rho5 = NA,
    lower.rho6 = NA, upper.rho6 = NA,
    lower.rho7 = NA, upper.rho7 = NA,
    lower.rho8 = NA, upper.rho8 = NA,
    lower.rho9 = NA, upper.rho9 = NA,
    lower.rho10 = NA, upper.rho10 = NA)
  #
  if (n.out >= 3) {
    run.data$var3 <- data$var$var3
    #
    run.data$lower.rho2 <- lower.rho2
    run.data$lower.rho3 <- lower.rho3
    #
    run.data$upper.rho2 <- upper.rho2
    run.data$upper.rho3 <- upper.rho3
  }
  #
  if (n.out >= 4) {
    run.data$var4 <- data$var$var4
    #
    run.data$lower.rho4 <- lower.rho4
    run.data$lower.rho5 <- lower.rho5
    run.data$lower.rho6 <- lower.rho6
    #
    run.data$upper.rho4 <- upper.rho4
    run.data$upper.rho5 <- upper.rho5
    run.data$upper.rho6 <- upper.rho6
  }
  #
  if (n.out >= 5) {
    run.data$var5 <- data$var$var5
    #
    run.data$lower.rho7 <- lower.rho7
    run.data$lower.rho8 <- lower.rho8
    run.data$lower.rho9 <- lower.rho9
    run.data$lower.rho10 <- lower.rho10
    #
    run.data$upper.rho7 <- upper.rho7
    run.data$upper.rho8 <- upper.rho8
    run.data$upper.rho9 <- upper.rho9
    run.data$upper.rho10 <- upper.rho10
  }
  #
  if (multiarm)
    run.data$treat3 <- data$T[, 3]
  else
    run.data$treat3 <- NULL
  #
  if (n.out == 2) {
    run.data$var3 <- run.data$var4 <- run.data$var5 <- NULL
    #
    run.data$lower.rho2 <- run.data$lower.rho3 <- run.data$lower.rho4 <-
      run.data$lower.rho5 <- run.data$lower.rho6 <- run.data$lower.rho7 <-
      run.data$lower.rho8 <- run.data$lower.rho9 <- run.data$lower.rho10 <-
      NULL
    #
    run.data$upper.rho2 <- run.data$upper.rho3 <- run.data$upper.rho4 <-
      run.data$upper.rho5 <- run.data$upper.rho6 <- run.data$upper.rho7 <-
      run.data$upper.rho8 <- run.data$upper.rho9 <- run.data$upper.rho10 <-
      NULL
    #
    params <- c("d1", "d2", 
                "psi1", "psi2",
                "rho1")
    #
    if (multiarm)
      model.file <- system.file("model", "mvnma_2_3arm.txt", package = "mvnma")
    else {
      model.file <- system.file("model", "mvnma_2_2arm.txt", package = "mvnma")
      run.data$k2 <- NULL
    }
  }
  #
  else if (n.out == 3) {
    run.data$var4 <- run.data$var5 <- NULL
    #
    run.data$lower.rho4 <- run.data$lower.rho5 <- run.data$lower.rho6 <-
      run.data$lower.rho7 <- run.data$lower.rho8 <- run.data$lower.rho9 <-
      run.data$lower.rho10 <- NULL
    #
    run.data$upper.rho4 <- run.data$upper.rho5 <- run.data$upper.rho6 <-
      run.data$upper.rho7 <- run.data$upper.rho8 <- run.data$upper.rho9 <-
      run.data$upper.rho10 <- NULL
    #
    params <- c("d1", "d2", "d3", 
                "psi1", "psi2", "psi3",
                "rho1", "rho2", "rho3")
    #
    if (multiarm)
      model.file <- system.file("model", "mvnma_3_3arm.txt", package = "mvnma")
    else {
      model.file <- system.file("model", "mvnma_3_2arm.txt", package = "mvnma")
      run.data$k2 <- NULL
    }
  }
  #
  else if (n.out == 4) {
    run.data$var5 <- NULL
    #
    run.data$lower.rho7 <- run.data$lower.rho8 <- run.data$lower.rho9 <-
      run.data$lower.rho10 <- NULL
    #
    run.data$upper.rho7 <- run.data$upper.rho8 <- run.data$upper.rho9 <-
      run.data$upper.rho10 <- NULL
    #
    params <- c("d1", "d2", "d3", "d4", 
                "psi1", "psi2", "psi3", "psi4",
                "rho1", "rho2", "rho3", "rho4", "rho5", "rho6")
    #
    if (multiarm)
      model.file <- system.file("model", "mvnma_4_3arm.txt", package = "mvnma")
    else
      model.file <- system.file("model", "mvnma_4_2arm.txt", package = "mvnma")
  }
  #
  else if (n.out == 5) {
    params <- c("d1", "d2", "d3", "d4", "d5",
                "psi1", "psi2", "psi3", "psi4", "psi5",
                "rho1", "rho2", "rho3", "rho4", "rho5", "rho6",
                "rho7", "rho8", "rho9", "rho10")
    #
    if (multiarm)
      model.file <- system.file("model", "mvnma_5_3arm.txt", package = "mvnma")
    else
      model.file <- system.file("model", "mvnma_5_2arm.txt", package = "mvnma")
  }
  
  
  #
  # Run Bayesian analysis
  #
  
  fit <- jags(
    data = run.data,
    inits = NULL,
    #
    parameters.to.save = params,
    #
    n.chains = 2, n.iter = n.iter, n.burnin = n.burnin, n.thin = 1,
    #
    DIC = FALSE,
    #
    model.file = model.file,
    quiet = quiet)
  #
  samples <- fit$BUGSoutput$sims.list
  colnames(samples$d1) <- trts
  colnames(samples$d2) <- trts
  #
  # Manipulate the results and create suitable datasets
  #
  res <- gather_results(fit,
                        outcomes = outclab,
                        trts = trts,
                        treat_out = treat_out,
                        reference.group = reference.group,
                        level = level)
  #
  attr(res, "outcomes") <- outclab
  attr(res, "trts") <- trts
  attr(res, "reference.group") <- reference.group
  attr(res, "level") <- level
  attr(res, "sm") <- attributes(data)$sm
  #
  class(res) <- "mvnma"
  #
  res
}
