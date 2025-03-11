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
    
    neg.sol <- rep(0,ncol(x))
    
  }else if ((!is.null(neg.sol)) & (length(neg.sol)!=ncol(x))) {
    
    stop("Please specify a negative ideal solution for each outcome.",
         call. = FALSE)
    
  }
  
  ## assume equal outcome weights if "weights" argument is NULL
  if (is.null(weights)) {
    
    weights <- rep(1/ncol(x), ncol(x))
    
  }else if ((!is.null(weights) & (length(weights)!=ncol(x)))) {
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
  Q <- sj <- rj <- vector("numeric", altno)
  #
  for (i in seq_len(altno)) {
    sj[i] <- sum(wnm[i, ])
    rj[i] <- max(wnm[i, ])
  }
  #
  for (i in seq_len(altno)) {
    Q[i] <- (v * (sj[i] - min(sj)) / (max(sj) - min(sj))) + 
      ((1 - v) * (rj[i] - min(rj)) / (max(rj) - min(rj)))
  }
  #
  S <- as.data.frame(sj)
  R <- as.data.frame(rj)
  #
  names(Q) <- row.names(S) <- row.names(R) <- row.names(x)
  
  Q <- as.data.frame(Q)
  names(Q) <- "Q"
  Q %<>% arrange(Q)
  #
  S <- as.data.frame(S)
  names(S) <- "S"
  S <- S %<>% arrange(S)
  #  
  R <- as.data.frame(R)
  names(R) <- "R"
  R %<>% arrange(R)
  
  res <- list(Q = Q, S = S, R = R)
  class(res) <- c("vikor", class(res))
  attr(res, "performance.table") <- x
  #
  res
}