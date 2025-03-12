#' Forest plot for multivariate network meta-analysis
#' 
#' @description
#' Draws a forest plot in the active graphics window (using grid graphics system).
#' 
#' @param x An object of class \code{\link{mvnma}}.
#' @param backtransf ...
#' @param leftcols ...
#' @param leftlabs ...
#' @param rightcols ...
#' @param rightlabs ...
#' @param col.study ...
#' @param col.square ...
#' @param col.square.lines ...
#' @param squaresize ...
#' @param header.line ...
#' @param col.subgroup ...
#' @param \dots Additional arguments for \code{\link[meta]{forest.meta}} function.
#'
#' @examples
#' \donttest{
#' library("netmeta")
#' 
#' # Use 'pairwise' to obtain contrast based data for the first two outcomes
#' data("Linde2015")
#' # Early response
#' p1 <- pairwise(treat = list(treatment1, treatment2, treatment3),
#'   event = list(resp1, resp2, resp3), n = list(n1, n2, n3),
#'   studlab = id, data = dat.linde2015, sm = "OR")
#' # Early remissions
#' p2 <- pairwise(treat = list(treatment1, treatment2, treatment3),
#'   event = list(remi1, remi2, remi3), n = list(n1, n2, n3),
#'   studlab = id, data = dat.linde2015, sm = "OR")
#'
#' # Perform analysis considering the efficacy outcomes
#' p12 <- list(p1, p2)
#'
#' # Use 'mvdata()' to transform the data in suitable JAGS format
#' data12 <- mvdata(p12)
#'
#' # Define outcome labels
#' outcomes <- c("Early_Response", "Early_Remission")
#'  
#' # Fit the model combining only the two efficacy outcomes
#' set.seed(1909)
#' mvnma12 <- mvnma(data = data12, 
#'   reference.group = "Placebo", outlab = outcomes,
#'   n.iter = 1000, n.burnin = 100)
#' mvnma12
#' 
#' # Generate a forest plot with the results 
#' forest(mvnma12)                 
#' }
#'
#' 
#' @method forest mvnma 
#' @export

forest.mvnma <- function(x, backtransf = FALSE,
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
  
  chkclass(x, "mvnma")
  #
  x <- x[names(x) != "outcome_correlation"]
  
  n.out <- length(x)
  
  sm <- attributes(x)$sm
  
  # Get rid of warning "no visible binding for global variable"
  treat <- mean <- sd <- lower <- upper <- NULL
  
  # Get estimates for each outcome
  #
  ests <- vector("list")
  #
  for (i in seq_len(n.out)) {
    ests[[i]] <- x[[i]]$basic_estimates
    ests[[i]]$treat <- row.names(ests[[i]])
    row.names(ests[[i]]) <- NULL
    #
    ests[[i]] %<>% select(treat, mean, sd, lower, upper)
    ests[[i]]$outcome <- attributes(x)$names[i]
  }
  
  
  # ests_1 <- x[[1]]$basic_estimates
  # 
  # ests_1 <- ests_1 %>% 
  #   mutate("treat" = row.names(ests_1)) %>% 
  #   select(treat,mean,sd,`2.5%`,`97.5%`)
  #   
  # ests_1$outcome <- attributes(x)$outlab[1]
  # 
  # ests_2 <- x[[2]]$basic_estimates
  # 
  # ests_2 <- ests_2 %>% 
  #   mutate("treat" = row.names(ests_2)) %>% 
  #   select(treat,mean,sd,`2.5%`,`97.5%`)
  # 
  # ests_2$outcome <- attributes(x)$outlab[2]
  
  
  ### construct final dataset    '
  
  dat <- bind_rows(ests)
  row.names(dat) <- NULL
  names(dat) <- c("studlab", "mean", "seTE", "lower", "upper", "outcome")
  
  m <- metagen(dat$mean, dat$seTE, sm = NULL,
               subgroup = dat$outcome,
               backtransf = backtransf,
               print.subgroup.name = FALSE,
               studlab = dat$studlab,
               common = FALSE, random = FALSE, hetstat = FALSE,
               method.tau = "DL", method.tau.ci = "")
  #
  ret <- forest(m,
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
                ...)
  #
  invisible(ret)
}
