#' Scatter plot to visualize the ranking lists for two outcomes.
#' 
#' @description
#' Draw a scatter plot (using grid graphics system) in the active graphics
#' window or store the forest plot in a file.
#' 
#' @param x An object of class \code{\link{mvrank}}.
#' @param outcome A mandatory numeric vector of length 2 specifying which
#'   outcomes should be plotted. For example, setting "outcome=c(2,3)" implies
#'   that a scatter plot will be generated plotting the rankings of outcomes 2
#'   and 3.
#' @param pos ...
#' @param cex.point ...
#' @param cex.label ...
#' @param pch ...
#' @param ... Additional arguments for \code{\link{plot}} function.
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
#' outlab <- c("Early_Response", "Early_Remission")
#' 
#' # Fit the model combining only the two efficacy outcomes
#' set.seed(1909)
#' mvnma12 <- mvnma(data = data12,
#'   reference.group = "Placebo", outlab = outlab[1:2],
#'   n.iter = 1000, n.burnin = 100)
#'            
#' # Extract treatment effect estimates and heterogeneity for Early_Response 
#' mvnma12$Early_Response$basic_estimates
#'                 
#' # Extract treatment effect estimates and heterogeneity for Early_Response 
#' mvnma12$Early_Response$basic_estimates
#' 
#' # Get all estimates
#' league12 <- league(mvnma12)
#' league12
#' 
#' # Rank treatments using sucra
#' ranks12 <- mvrank(mvnma12, small.values = c("und","und"), method = "sucra")
#' ranks12
#' 
#' # Visualize SUCRAs in a scatter plot with outcome 1
#' # (as specified in the mvdata() function) in the x-axis and outcome 2
#' # (as specified in the mvdata() function) in the y-axis
#' 
#' plot(ranks12, outcome = c(1, 2))
#' 
#' # Visualize SUCRAs in a scatter plot with outcome 2
#' # (as specified in the mvdata() function) in the x-axis and outcome 1
#' # (as specified in the mvdata() function) in the y-axis
#' 
#' plot(ranks12, outcome = c(2, 1))
#' }
#' 
#' @method plot mvrank 
#' @export  

plot.mvrank <- function(x, outcome = NULL, pos = 1,
                        cex.point = 1, cex.label = 0.7, pch = 19,
                        ...) {
  
  chkclass(x, "mvrank")
  #
  # Get rid of warning "no visible binding for global variable"
  treat1 <- treat2 <- NULL
  
  chknull(outcome)
  #
  if (length(outcome) < 2)
    stop("Argument  'outcome' should be of length 2.")
  else if (length(outcome) > 2) {
    outcome <- c(1, 2)
    warning("Argument 'outcome' should be of length 2. The produced scatter ",
            "plot now refers to outcomes 1 and 2.")
  }
  #
  outlab <- names(x)[c(outcome[1], outcome[2])]
  
  common_trts <- attributes(x)$common_trts
  
  r1 <- x[[outcome[1]]]
  names(r1)[1:2] <- c("treat1", paste(names(r1)[2],"1",sep = ""))
  r1$out1 <- outlab[1]
  #
  r1 %<>% filter(treat1 %in% common_trts) %>% arrange(treat1)
  
  r2 <- x[[outcome[2]]]
  names(r2)[1:2] <- c("treat2", paste(names(r2)[2],"2",sep = ""))
  r2$out2 <- outlab[2]
  #
  r2 %<>% filter(treat2 %in% common_trts) %>% arrange(treat2)
  #
  ranking <- cbind.data.frame(r1, r2)  
  #
  plot(ranking[, 2], ranking[, 5], main = "",
       cex = cex.point,
       xlab = outlab[1], ylab = outlab[2],
       pch = pch, ...)
  #
  text(ranking[,2], ranking[,5],labels = ranking$treat1,
       cex = cex.label, pos = pos, col = "black")
  #
  invisible(NULL)  
}
