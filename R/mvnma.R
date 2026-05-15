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
#' @param n.chains Number of Markov chains (default=4). 
#' @param n.domain Integer indicating the position of the last outcome in the 
#' first outcome domain (based on the order of the supplied pairwise objects). 
#' Used with `method = "DM"` to restrict information sharing within outcome 
#' domains. Ignored when `method = "standard"`. Default is `NULL`.
#' @param n.thin Thinning rate. Default is equal to
#'   \code{max(1, floor((n.iter - n.burnin) / 1000))}.
#' @param n.iter Number of iterations (default: 10000).
#' @param n.burnin Number of iterations for burn-in (default: 2000).
#' @param level The level used to calculate confidence intervals
#'   for network estimates.
#' @param scale.psi Values for the scale parameter(s) of the Half-Normal prior
#'   used for the heterogeneity parameters within each outcome. If NULL, all
#'   values are set to 1. If specified, it should have a length equal to the
#'   number of outcomes.
#' @param lower.rho Lower bounds for the Uniform prior(s) used for the
#'   correlation coefficient. If NULL all bounds are set to -1.
#' @param upper.rho Upper bounds for the Uniform prior(s) used for the
#'   correlation coefficient. If NULL all bounds are set to 1.
#' @param method A character string specifying the method to be used for model
#'   fitting. This can be either "standard" (default), referring to the
#'   standard bivariate model, or "DM", referring to the bivariate model based
#'   on the DuMouchel method. The argument can be abbreviated.
#' @param quiet A logical indicating whether to print information on the
#'   progress of the JAGS model fitting.
#' @param x An object of class \code{\link{mvnma}}.
#' @param digits Minimal number of significant digits, see
#'   \code{print.default}.
#' @param digits.sd Minimal number of significant digits for standard
#'   deviations
#' @param print.sd A logical specifying whether standard deviations should be
#'   printed.
#' @param \dots Additional arguments (ignored)
#' 
#' @details
#' The function \code{\link{mvnma}} expects two to five outcomes /
#' \code{\link[meta]{pairwise}} objects. A common reference treatment across
#' all outcomes is required to only show comparisons with the reference in
#' forest plots.
#' 
#' The Bayesian multivariate network meta-analysis model fitted in the
#' \bold{mvnma} package assumes uniform priors for the between-outcome
#' correlation coefficients. The lower and upper bounds of these priors can be
#' defined using the arguments `lower.rho` and `upper.rho`. If not set, the
#' model will assume a `Unif (-1, 1)` prior for all correlation coefficients.
#' For two outcomes, a single value can be provided for `lower.rho` and
#' `upper.rho`. For example, `lower.rho` = 0.5 and `upper.rho` = 1 for
#' rho12 ~ Unif (0.5, 1)).
#' For more than two outcomes, the order in which the bounds are provided
#' matters. For example, when pooling four outcomes, the lower and
#' upper bounds correspond to the following order of correlation coefficients:
#' (rho12, rho13, rho14, rho23, rho24, rho34).
#' 
#' Two types of priors for the treatment effect parameters are supported via 
#' the argument `method`. Setting `method = "standard"` fits an mvNMA model 
#' using non-informative normal priors (e.g., `N(0, 10^3)`).
#' 
#' Alternatively, `method = "DM"` specifies the DuMouchel prior, which assumes 
#' constant relative treatment effects across outcomes and enables information 
#' sharing (DuMouchel & Harris, 1983). This may improve precision but can
#' introduce bias when outcomes from different domains (e.g., efficacy and
#' safety) are analyzed jointly.
#' 
#' The argument `n.domain` can be used to restrict information sharing to 
#' predefined outcome domains. It indicates the position (based on the order 
#' of the supplied pairwise objects) of the last outcome in the first domain. 
#' For example, with four outcomes, setting `n.domain = 2` assigns the first 
#' two outcomes to one domain and the remaining outcomes to a second domain. 
#' In this case, information is shared only within domains.
#' 
#' By default, `n.domain = NULL`, in which case information is shared across 
#' all outcomes when `method = "DM"`. This may be appropriate when all outcomes 
#' belong to the same domain or when cross-domain sharing is justified.
#' 
#' The argument `n.domain` is ignored when `method = "standard"`.
#'  
#' @return
#' The function returns an 'mvnma' object. This consists of the results for each
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
#'
#' @references
#' Achana FA, Cooper NJ, Bujkiewicz S, Hubbard SJ, Kendrick D, Jones DR,
#' Sutton AJ (2014):
#' Network meta-analysis of multiple outcome measures accounting for borrowing
#' of information across outcomes.
#' \emph{BMC Medical Research Methodology},
#' \bold{21}, 92
#' 
#' DuMouchel WH, Harris JE (1983):
#' Bayes methods for combining the results of cancer studies in humans and
#' other species.
#' \emph{Journal of the American Statistical Association},
#' \bold{78}, 293--308
#' 
#' Efthimiou O, Mavridis D, Riley RD, Cipriani A, Salanti G (2015):
#' Joint synthesis of multiple correlated outcomes in networks of interventions.
#' \emph{Biostatistics}, 
#' \bold{16}, 84--97
#' 
#' Riley RD, Thompson JR, Abrams KR (2008):
#' An alternative model for bivariate random-effects meta-analysis when the
#' within-study correlations are unknown.
#' \emph{Biostatistics},
#' \bold{9}, 172--86
#' 
#' @examples
#' # Use 'pairwise' to obtain contrast based data for each one of the five
#' # available outcomes 
#' 
#' # Early response
#' pw1 <- pairwise(treat = list(treatment1, treatment2, treatment3),
#'   event = list(resp1, resp2, resp3), n = list(n1, n2, n3),
#'   studlab = id, data = Linde2015, sm = "OR")
#'
#' # Early remissions
#' pw2 <- pairwise(treat = list(treatment1, treatment2, treatment3),
#'   event = list(remi1, remi2, remi3), n = list(n1, n2, n3),
#'   studlab = id, data = Linde2015, sm = "OR")
#'
#' # Adverse events
#' pw3 <- pairwise(treat = list(treatment1, treatment2,treatment3),
#'   event = list(ae1, ae2, ae3),  n = list(n1, n2, n3),
#'   studlab = id, data = Linde2015, sm = "OR")
#'
#' # Loss to follow-up
#' pw4 <- pairwise(treat = list(treatment1, treatment2, treatment3),
#'   event = list(loss1, loss2, loss3), n = list(n1, n2, n3),
#'   studlab = id, data = Linde2015, sm = "OR")
#'
#' # Loss_to_follow_up_(AE)
#' pw5 <- pairwise(treat = list(treatment1, treatment2, treatment3),
#'   event = list(loss.ae1, loss.ae2, loss.ae3), n = list(n1, n2, n3),
#'   studlab = id, data = Linde2015, sm = "OR")
#'
#' # Define outcome labels
#' outcomes <- c("Early_Response", "Early_Remission",
#'   "Adverse_events", "Loss_to_follow_up", "Loss_to_follow_up_AE")
#'  
#' # Fit the model combining only the two efficacy outcomes
#' set.seed(1909)
#' mvnma12 <- mvnma(pw1, pw2,
#'   reference.group = "Placebo", outclab = outcomes[1:2],
#'   n.iter = 100, n.burnin = 20)
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
#' # Print odds ratios for efficacy outcomes
#' outc <- names(mvnma12)[names(mvnma12) != "cor"]
#' #
#' for (i in outc) {
#'   cat(paste0("\nOutcome: ", i, "\n\n"))
#'   print(round(exp(mvnma12[[i]]$TE.random), 2))
#' }
#' 
#' \dontrun{
#' # Fit the model combining all five outcomes
#' set.seed(1904)
#' mvnma_all <- mvnma(pw1, pw2, pw3, pw4, pw5,
#'   reference.group = "Placebo", outclab = outcomes,
#'   n.iter = 100, n.burnin = 20)
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
#' # Print odds ratios for all outcomes
#' outc <- names(mvnma_all)[names(mvnma_all) != "cor"]
#' #
#' for (i in outc) {
#'   cat(paste0("\nOutcome: ", i, "\n\n"))
#'   print(round(exp(mvnma_all[[i]]$TE.random), 2))
#' }
#' }
#' 
#' @export mvnma

mvnma <- function(...,
                  reference.group = NULL, outclab = NULL,   
                  n.domain = NULL,
                  n.chains = 4, n.iter = 10000, 
                  n.burnin = 2000, 
                  n.thin = max(1, floor((n.iter - n.burnin) / 1000)), 
                  level = gs("level.ma"),
                  scale.psi,
                  lower.rho, upper.rho,
                  method = "standard",
                  quiet = FALSE) {
  
  is_pairwise <- function(x)
    inherits(x, "pairwise")
  #
  args <- list(...)
  #
  n.out <- length(args)
  #
  chknumeric(n.domain,min=1,max=n.out)
  n.dom <- n.domain
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
  #
  method <- setchar(method, c("standard", "DM"))
  #
  # extract number of outcomes  
  n.out <- ncol(data$var)
  n.cor <- choose(n.out, 2)
  #
  miss.lower <- missing(lower.rho)
  miss.upper <- missing(upper.rho)
  miss.scale.psi <- missing(scale.psi)
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
  
  #
  if (!miss.scale.psi)
    chknumeric(scale.psi,zero = TRUE, min = 0, length = n.out, NA.ok = FALSE)
  
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
  
  # Create values for the scale and precision of parameter psi
  
  if (miss.scale.psi) {
    scale.psi1 <- 1
    scale.psi2 <- 1
  }  
  else {
    scale.psi1 <- scale.psi[1]
    scale.psi2 <- scale.psi[2]
  }
  
  prec.psi1 <- 1 / scale.psi1^2
  
  prec.psi2 <- 1 / scale.psi2^2
  
  if (n.out >= 3) {
    if (miss.scale.psi)
      scale.psi3 <- 1
    else
      scale.psi3 <- scale.psi[3]
    #
    prec.psi3 <- 1 / scale.psi3^2
  }
  
  if (n.out >= 4) {
    if (miss.scale.psi)
      scale.psi4 <- 1
    else
      scale.psi4 <- scale.psi[4]
    #
    prec.psi4 <- 1 / scale.psi4^2
  }
  
  if (n.out >= 5) {
    if (miss.scale.psi)
      scale.psi5 <- 1
    else
      scale.psi5 <- scale.psi[5]
    #
    prec.psi5 <- 1 / scale.psi5^2
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
    prec.psi1 = prec.psi1, prec.psi2=prec.psi2,
    prec.psi3 = NA, prec.psi4 = NA, prec.psi5 = NA,
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
    run.data$prec.psi3 <- prec.psi3
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
    run.data$prec.psi4 <- prec.psi4
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
    run.data$prec.psi5 <- prec.psi5
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
    run.data$prec.psi3 <- run.data$prec.psi4 <- run.data$prec.psi5 <- NULL
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
    model.code <- mvnma_code(n.out, method, multiarm,n.dom)
  }
  #
  else if (n.out == 3) {
    run.data$var4 <- run.data$var5 <- NULL
    #
    run.data$prec.psi4 <- run.data$prec.psi5 <- NULL
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
    model.code <- mvnma_code(n.out, method, multiarm,n.dom)
  }
  #
  else if (n.out == 4) {
    run.data$var5 <- NULL
    #
    run.data$prec.psi5 <- NULL
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
    model.code <- mvnma_code(n.out, method, multiarm,n.dom)
  }
  #
  else if (n.out == 5) {
    params <- c("d1", "d2", "d3", "d4", "d5",
                "psi1", "psi2", "psi3", "psi4", "psi5",
                "rho1", "rho2", "rho3", "rho4", "rho5", "rho6",
                "rho7", "rho8", "rho9", "rho10")
    #
    model.code <- mvnma_code(n.out, method, multiarm,n.dom)
  }
  #
  if (method == "DM")
    if(is.null(n.domain)){
    params <- c(params, "sigma")  
    }else{
    params <- c(params, "sigma1","sigma2") 
    }
  #
  if (!multiarm)
    run.data$k <- NULL
  
  
  #
  # Run Bayesian analysis
  #
  
  fit <- jags(
    data = run.data,
    inits = NULL,
    #
    parameters.to.save = params,
    #
    n.chains = n.chains, n.iter = n.iter, 
    n.burnin = n.burnin, n.thin = n.thin,
    #
    DIC = FALSE,
    #
    model.file = textConnection(model.code),
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
                        level = level,
                        n.domain = n.domain,
                        method = method)
  #
  attr(res, "outcomes") <- outclab
  attr(res, "trts") <- trts
  attr(res,"n.domain") <- n.domain
  attr(res, "reference.group") <- reference.group
  attr(res, "level") <- level
  attr(res, "sm") <- attr(data, "sm")
  attr(res, "method.model") <- method
  attr(res, "model.code") <- model.code
  attr(res, "fit") <- fit
  attr(res, "params") <- params
  #
  class(res) <- "mvnma"
  #
  res
}


#' @rdname mvnma
#' @method print mvnma
#' @export

print.mvnma <- function(x,
                        digits = gs("digits"),
                        digits.sd = gs("digits.sd"),
                        print.sd = FALSE,
                        ...) {
  
  chkclass(x, "mvnma")
  #
  chknumeric(digits, min = 0, length = 1)
  chknumeric(digits.sd, min = 0, length = 1)
  chklogical(print.sd)
  #
  level <- attr(x, "level")
  reference.group <- attr(x, "reference.group")
  method <- attr(x, "method")
  n.domain <- attr(x,"n.domain")
  #
  ci.lab <- paste0(round(100 * level, 1), "%-CI")
  #
  x <- x[names(x) != "cor"]
  #
  if (method == "DM") {
    if (is.null(n.domain)) {
      x <- x[names(x) != "sigma"]
    }
    else {
      x <- x[!(names(x) %in% c("sigma1","sigma2"))]
    }
  }
  nam <- names(x)
  
  # Get rid of warning "no visible binding for global variable"
  lower <- upper <- psi <- NULL
  #
  for (i in seq_along(nam)) {
    cat(paste0(if (i > 1) "\n" else "", "Outcome: ", nam[i], "\n\n"))
    #
    dat.i <- x[[i]]$basic_estimates
    dat.i <- dat.i[rownames(dat.i) != reference.group, ]
    #
    dat.i$mean <- formatN(dat.i$mean, digits = digits)
    #
    if (!print.sd)
      dat.i$sd <- NULL
    else
      dat.i$sd <- formatN(dat.i$sd, digits = digits.sd)
    #
    dat.i$lower <- formatCI(formatN(dat.i$lower, digits = digits),
                            formatN(dat.i$upper, digits = digits))
    dat.i %<>% select(-upper)
    names(dat.i)[names(dat.i) == "lower"] <- ci.lab
    #
    dat.i$Rhat <- formatN(dat.i$Rhat, digits = 4)
    #
    rownames(dat.i) <- paste0("d[", rownames(dat.i), "]")
    #
    prmatrix(dat.i, quote = FALSE, right = TRUE)
  }
  #
  invisible(NULL)
}
