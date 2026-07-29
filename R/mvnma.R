#' Perform a Bayesian multivariate network meta-analysis using a
#' single-correlation coefficient model
#' 
#' @description
#' This function fits a Bayesian multivariate network meta-analysis model for
#' two or more outcomes. Additionally, the studies can have multiple arms.
#' 
#' @param \dots Either two or more pairwise objects or a single list with
#'   two or more pairwise objects.
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
#' @param varTE.missing Assumed (very large) variance for outcomes not reported
#'   in a study. By default, the largest variance times 1000000 is used. This is
#'   the same value used for argument \code{seTE.ignore} in
#'   \code{\link[netmeta]{netimpact}} to mimicking the removal of individual
#'   studies from the network meta-analysis.
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
#' The multivariate network meta-analysis (mvNMA) model supported by this
#' package refers to the single correlation coefficient model, interpreted as
#' an amalgam of within- and across-outcome correlations
#' (Efthimiou et al., 2015) which is a generalisation of Riley et al. (2008).
#' 
#' The function \code{\link{mvnma}} expects two or more outcomes /
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
#' # Use 'pairwise' to obtain contrast based data for the first two outcomes
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
#' # Define outcome labels
#' outcomes <- c("Early_Response", "Early_Remission",
#'   "Adverse_events", "Loss_to_follow_up", "Loss_to_follow_up_AE")
#' 
#' # Fit the model combining only the two efficacy outcomes
#' # (note, we are using only 10 iterations and 2 burnins to reduce the
#' #  runtime of the example; in real applications use larger numbers)
#' set.seed(1910)
#' mvnma(pw1, pw2,
#'   reference.group = "Placebo", outclab = outcomes[1:2],
#'   n.iter = 10, n.burnin = 2)
#' 
#' \donttest{
#' # Use 'pairwise' to obtain contrast based data for the third to fifth
#' # outcome
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
#' # Fit the model combining only the two efficacy outcomes
#' # (note, we are using only 100 iterations and 20 burnins to reduce the
#' #  runtime of the example; in real applications use larger numbers)
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
#' # Fit the model combining all five outcomes
#' # (note, we are using only 100 iterations and 20 burnins to reduce the
#' #  runtime of the example; in real applications use larger numbers)
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
                  #
                  method = "standard",
                  n.domain = NULL,
                  #
                  reference.group = NULL, outclab = NULL,   
                  #
                  n.chains = 4, n.iter = 10000, 
                  n.burnin = 2000, 
                  n.thin = max(1, floor((n.iter - n.burnin) / 1000)), 
                  #
                  level = gs("level.ma"),
                  #
                  scale.psi,
                  lower.rho, upper.rho,
                  #
                  varTE.missing = NULL,
                  quiet = FALSE) {
  
  # Get rid of warning "no visible binding for global variable"
  studlab <- NULL
  
  
  #
  #
  # (1) Extract pairwise() objects
  #
  #
  
  args <- list(...)
  #
  if (length(args) == 1) {
    if (inherits(args[[1]], "pairwise"))
      stop("Provide two or more pairwise objects.",
           call. = FALSE)
    #
    if (!is.list(args[[1]]))
      stop("All elements of argument '...' must be of class 'pairwise'.",
           call. = FALSE)
    #
    n.args <- length(args[[1]])
    #
    args2 <- vector("list", n.args)
    #
    for (i in seq_len(n.args))
      args2[[i]] <- args[[1]][[i]]
    #
    args <- args2
  }
  #
  n.out <- length(args)
  n.rho <- choose(n.out, 2)
  #
  if (n.out < 2)
    stop("Provide two or more pairwise objects.",
         call. = FALSE)
  #  
  for (i in seq_len(n.out)) {
    if (!inherits(args[[i]], "pairwise"))
      stop("All elements of argument '...' must be of class ",
           "'pairwise'.",
           call. = FALSE)
  }
  #
  sm <- vector("character", n.out)
  reference.groups <- vector("character", n.out)
  trts_list <- vector("list", n.out)
  #
  for (i in seq_len(n.out)) {
    sm[i] <- attr(args[[i]], "sm")
    reference.groups[i] <- attr(args[[i]], "reference.group")
    trts_list[[i]] <- sort(unique(c(args[[i]]$treat1, args[[i]]$treat2)))
  }
  
  
  #
  #
  # (2) Check and set additional arguments
  #
  #
  
  method <- setchar(method, c("standard", "DM"))
  #
  chknumeric(n.domain, min = 1, max = n.out)
  #
  if (is.null(reference.group)) {
    if (length(unique(reference.groups)) == 1)
      reference.group <- unique(reference.groups)
    else
      stop("Argument 'reference.group' must be specified as it differs in ",
           "pairwise() objects:\n  ",
           paste0("'", reference.groups, "'", collapse = ", "),
           call. = FALSE)
  }
  else {
    reference.group <-
      setchar(reference.group, sort(unique(unlist(trts_list))))
  }
  #
  if (is.null(outclab))
    outclab <- paste("outcome", seq_len(n.out), sep = "_")  
  else if (length(outclab) != n.out)
    stop("Please provide labels for all outcomes.")
  #
  chknumeric(n.chains, min = 1, length = 1)
  chknumeric(n.iter, min = 1, length = 1)
  chknumeric(n.burnin, min = 1, length = 1)
  chknumeric(n.thin, min = 1, length = 1)
  #
  chklevel(level)
  #
  if (missing(scale.psi))
    scale.psi <- rep_len(1, n.out)
  else
    chknumeric(scale.psi, min = 0, zero = TRUE, length = n.out, NA.ok = FALSE)
  #
  prec.psi <- 1 / scale.psi^2
  #
  miss.lower.rho <- missing(lower.rho)
  miss.upper.rho <- missing(upper.rho)
  #
  if (miss.lower.rho)
    lower.rho <- rep_len(-1, n.rho)
  else
    chknumeric(lower.rho, min = -1, max = 1, length = n.rho, NA.ok = FALSE)
  #
  if (miss.upper.rho)
    upper.rho <- rep_len(1, n.rho)
  else
    chknumeric(upper.rho, min = -1, max = 1, length = n.rho, NA.ok = FALSE)
  #
  if (!miss.lower.rho & !miss.upper.rho) {
    if (any(lower.rho >= upper.rho))
      stop("Values for argument 'lower.rho' must be smaller than values for ",
           "argument 'upper.rho'.",
           call. = FALSE)
  }
  #
  if (!is.null(varTE.missing))
    chknumeric(varTE.missing, min = 0, zero = TRUE, length = 1)
  #
  chklogical(quiet)
  
  
  #
  #
  # (3) Create list with JAGS input
  #
  #
  
  dat <- mvdata(args)
  #
  # Check number of extracted outcomes
  #
  if (n.out != ncol(dat$var %>% select(-studlab)))
    stop("Number of variances and outcomes differ.", call. = FALSE)
  #
  trts.list <- dat$trts.list
  trts <- dat$trts
  #
  id_reference.group <- unname(which(trts == reference.group))
  #
  dat_var <- dat$var %>% filter(!duplicated(studlab))
  rownames(dat_var) <- dat_var$studlab
  dat_var %<>% select(-studlab)
  #
  control_matrix <- 1L * !is.na(dat_var)
  #
  var_matrix <- as.matrix(dat$var %>% select(-studlab))
  #
  if (is.null(varTE.missing)) {
    varTE.missing <- 1000^2 * max(var_matrix, na.rm = TRUE)
  }
  else {
    if (varTE.missing < max(dat$var %>% select(-studlab), na.rm = TRUE))
      stop("The value provided for argument 'varTE.missing' must be larger ",
           "than the largest available variance in the dataset.",
           call. = FALSE)
  }
  #
  var_matrix[is.na(var_matrix)] <- varTE.missing
  #
  dat_jags <- list(
    y = dat$y,
    #
    varmat = var_matrix,
    contmat = control_matrix,
    trtmat = dat$treatments,
    ref = id_reference.group,
    #
    n.studies = dat$n.studies,
    n = dat$n,
    #
    prec.psi = prec.psi, lower.rho = lower.rho, upper.rho = upper.rho
  )
  
  
  #
  #
  # (3) Run Bayesian analysis
  #
  #
  
  params <- c(paste0("d", seq_len(n.out)), "psi", "rho")
  #
  if (method == "DM") {
    if (is.null(n.domain))
      params <- c(params, "sigma")
    else
      params <- c(params, c("sigma1", "sigma2"))
  }
  #
  model.code <- mvnma_code(n.out, dat$arms, method, n.domain)
  #
  text_conn <- textConnection(model.code)
  on.exit(close(text_conn), add = TRUE)
  #
  fit <- jags(
    data = dat_jags,
    inits = NULL,
    #
    parameters.to.save = params,
    #
    n.chains = n.chains, n.iter = n.iter, 
    n.burnin = n.burnin, n.thin = n.thin,
    #
    DIC = FALSE,
    #
    model.file = text_conn,
    quiet = quiet)
  #
  # Column names set to treatment names
  #
  for (i in seq_len(n.out))
    colnames(fit$BUGSoutput$sims.list[[paste0("d", i)]]) <- trts
  #
  # Manipulate the results and create suitable datasets
  #
  res <- gather_results(fit,
                        outcomes = outclab,
                        trts = trts,
                        trts.list = trts.list,
                        reference.group = reference.group,
                        level = level,
                        n.domain = n.domain,
                        method = method)
  #
  attr(res, "outcomes") <- outclab
  attr(res, "trts") <- trts
  attr(res, "n.domain") <- n.domain
  attr(res, "reference.group") <- reference.group
  attr(res, "level") <- level
  attr(res, "sm") <- attr(dat, "sm")
  attr(res, "method.model") <- method
  attr(res, "model.code") <- model.code
  attr(res, "fit") <- fit
  attr(res, "params") <- params
  attr(res, "varTE.missing") <- varTE.missing
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
  n.domain <- attr(x, "n.domain")
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
      x <- x[!(names(x) %in% c("sigma1", "sigma2"))]
    }
  }
  #
  nam <- names(x)
  
  # Get rid of warning "no visible binding for global variable"
  lower <- upper <- NULL
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
