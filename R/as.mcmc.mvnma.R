#' Auxiliary function to extract the Markov Chain Monte Carlo object
#' 
#' @description
#' Extract the Markov Chain Monte Carlo object from an R object of class
#' \code{mvnma}.
#' 
#' @param x An object of class \code{mvnma}.
#' @param drop.reference.group A logical indicating whether to drop the MCMC
#'   samples for the reference group.
#' @param \dots Additional arguments passed on to
#'   \code{\link[coda]{as.mcmc.list}}.
#' 
#' @return
#' The function returns the Markov Chain Monte Carlo object which is a list
#' of 'mcmc' objects.
#' 
#' @aliases as.mcmc.mvnma as.mcmc
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{mvnma}}
#' 
#' @examples
#' # Locate file "mvnma_example.rda" with mvnma() results
#' .fname <- system.file("extdata/mvnma_examples.rda", package = "mvnma")
#' load(.fname)
#' 
#' mcmc12 <- as.mcmc(mvnma12)
#' 
#' library(coda)
#' oldpar <- par(mfrow = c(3, 3))
#' # Traceplot
#' traceplot(mcmc12)
#' # Density plot
#' densplot(mcmc12)
#' # Traceplot and density plot together
#' plot(mcmc12)
#' 
#' # Do not print the trace plot for the reference group, i.e., placebo
#' mcmc12.drop <- as.mcmc(mvnma12, drop.reference.group = TRUE)
#' par(mfrow = c(2, 4))
#' # traceplot
#' traceplot(mcmc12.drop)
#' # density plot
#' densplot(mcmc12.drop)
#' # traceplot and density plot together
#' plot(mcmc12.drop)
#' 
#' par(oldpar)
#' 
#' @method as.mcmc mvnma
#' @export

as.mcmc.mvnma <- function(x, drop.reference.group = FALSE, ...) {
  chkclass(x, "mvnma")
  #
  chklogical(drop.reference.group)
  #
  trts <- attr(x, "trts")
  reference.group <- attr(x, "reference.group")
  #
  res <- as.mcmc.list(attr(x, "fit")$BUGSoutput, ...)
  #
  for (i in seq_along(res)) {
    nam <- colnames(res[[i]])
    #
    for (j in seq_along(trts)) {
      patt <- paste0("\\b(d\\d+)\\[", j, "\\]")
      repl <-  paste0("\\1[", trts[j], "\\]")
      #
      nam <- gsub(pattern = patt, replacement = repl, x = nam)
    }
    #
    colnames(res[[i]]) <- nam
  }
  #
  if (drop.reference.group) {
    for (i in seq_along(res)) {
      nam <- colnames(res[[i]])
      #
      patt <- paste0("\\b(d\\d+)\\[", reference.group, "\\]")
      #
      keep <- !grepl(pattern = patt, x = nam)
      #
      res[[i]] <- res[[i]][, keep, drop = FALSE]
    }
  }
  #
  res
}
