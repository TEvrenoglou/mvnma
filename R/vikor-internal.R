vikor_internal <- function(x, weights, v, pos.sol, neg.sol) {
  
  minmax <- rep("max", ncol(x))
  
  if (is.null(pos.sol)) {
    pos.sol <- rep(1, ncol(x))
  }
  else if ((!is.null(pos.sol)) & (length(pos.sol)!=ncol(x))) {
    stop("Please specify a positive ideal solution for each outcome.",
         call. = FALSE)
  }
  
  if (is.null(neg.sol)) { 
    neg.sol <- rep(0, ncol(x))
  }
  else if ((!is.null(neg.sol)) & (length(neg.sol) != ncol(x))) {
    stop("Please specify a negative ideal solution for each outcome.",
         call. = FALSE)
  }
  
  # Assume equal outcome weights if argument 'weights' is NULL
  #
  if (is.null(weights)) {
    weights <- rep(1/ncol(x), ncol(x))
  }
  else if ((!is.null(weights) & length(weights) != ncol(x))) {
    stop("Please provide weights for all outcomes.",
         call. = FALSE)
  }
  
  if (!is.null(weights) && sum(weights) != 1) {
    weights <- round(weights / sum(weights), digits = 2)
    #
    warning("Weights should always sum up to 1. To do so the given weights ",
            "are now standardized and the new weights are: ",
            paste(weights, collapse = ", "))
  }
  #
  critno <- ncol(x)
  #
  dist <- matrix(ncol = critno, nrow = nrow(x),
                 dimnames = list(row.names(x), colnames(x)))
  #
  for (i in seq_len(critno))
    dist[, i] <- (pos.sol[i] - x[, i]) / (pos.sol[i] - neg.sol[i])
  
  altno <- nrow(x)
  wnm <- t(t(dist) * weights)
  #
  Q <- R <- S <- vector("numeric", altno)
  #
  for (i in seq_len(altno)) {
    R[i] <- max(wnm[i, ])
    S[i] <- sum(wnm[i, ])
    
  }
  
  for (i in seq_len(altno)) {
    Q[i] <- (v * (S[i] - min(S)) / (max(S) - min(S))) + 
      ((1 - v) * (R[i] - min(R)) / (max(R) - min(R)))
  }
  #
  res <- data.frame(Q, S, R, row.names = row.names(x)) %>% arrange(Q)
  
  #
  class(res) <- c("vikor", class(res))
  attr(res, "performance.table") <- x
  #
  res
}
