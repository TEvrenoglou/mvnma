#' Forest plot for multivariate network meta-analysis
#' 
#' @description
#' Draws a forest plot in the active graphics window (using grid graphics
#' system).
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
#' @param \dots Additional arguments for \code{\link[meta]{forest.meta}}
#'   function.
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
#' # Define outcome labels
#' outcomes <- c("Early Response", "Early Remission")
#'  
#' # Fit the model combining the two efficacy outcomes
#' set.seed(1909)
#' mvnma12 <- mvnma(p1, p2,
#'   reference.group = "Placebo", outclab = outcomes,
#'   n.iter = 1000, n.burnin = 100)
#' mvnma12
#' 
#' # Generate a forest plot with the results
#' forest(mvnma12)
#' }
#' 
#' @method forest mvnma 
#' @export

forest.mvnma <- function(x, backtransf = FALSE,
                         #
                         leftcols = "studlab", leftlabs,
                         rightcols = c("effect", "ci"),
                         rightlabs,
                         col.study = "black",
                         col.square = "black",
                         col.square.lines = "black",
                         squaresize = 0.7,
                         header.line = TRUE, 
                         col.subgroup = "black",
                         ...) {
  
  chkclass(x, "mvnma")
  #
  method.model <- attr(x, "method.model")
  reference.group <- attr(x, "reference.group")
  sm <- attr(x, "sm")
  #
  x <- x[names(x) != "cor"]
  #
  if (method.model == "DM")
    x <- x[names(x) != "sigma"]
  #
  n.out <- length(x)
  #
  if (missing(leftlabs))
    leftlabs <- paste0("Comparison with '", reference.group, "'")
  #
  if (missing(rightlabs))
    rightlabs <- rep(NA, length(rightcols))
  
  # Get rid of warning "no visible binding for global variable"
  treat <- mean <- sd <- lower <- upper <- studlab <- NULL
  
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
    ests[[i]]$outcome <- attr(x, "names")[i]
  }
  
  # Construct final dataset
  #
  dat <- bind_rows(ests)
  row.names(dat) <- NULL
  names(dat) <- c("studlab", "mean", "seTE", "lower", "upper", "outcome")
  #
  # Drop rows for reference group 
  #
  dat %<>% filter(studlab != reference.group)
  
  if (length(unique(sm)) == 1)
    sm <- unique(sm)
  else
    sm <- ""
  #
  m <- metagen(dat$mean, dat$seTE, sm = sm,
               subgroup = dat$outcome,
               backtransf = backtransf,
               print.subgroup.name = FALSE,
               studlab = dat$studlab,
               common = FALSE, random = FALSE, hetstat = FALSE,
               method.tau = "DL", method.tau.ci = "")
  #
  ret <- forest(m,
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
