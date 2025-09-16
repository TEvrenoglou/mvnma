vikor_internal <- function(x, weights, v) {
  
  minmax <- rep("max", ncol(x))
  
  # Assume equal outcome weights if argument 'weights' is NULL
  #
  if (is.null(weights))
    weights <- rep(1 /ncol(x), ncol(x))
  else if ((!is.null(weights) & length(weights) != ncol(x)))
    stop("Please provide weights for all outcomes.",
         call. = FALSE)
  #
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
  pos.sol <- rep(1, critno)
  neg.sol <- rep(0, critno) 
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
  #class(res) <- c("vikor", class(res))
  attr(res, "performance.table") <- x
  #
  res
}


performance_fuzzy <- function(x, trts, ...) {
  
  # Get rid of warning "no visible binding for global variable"
  #
  treatment <- lower.CrI <- median_rank <- upper.CrI <- NULL
  
  class(x) <- "list"
  
  for (i in seq_along(x)) {
    x[[i]] %<>%
      filter(treatment %in% trts) %>% 
      arrange(treatment) %>% 
      select(treatment, lower.CrI, median_rank, upper.CrI)
    #
    row.names(x[[i]]) <- x[[i]]$treatment
    #
    x[[i]] %<>% select(-treatment)
    #
    names(x[[i]]) <- paste(names(x[[i]]), i, sep = ".")
  }
  
  y <- bind_cols(x)    
  names(y) <- seq(1:ncol(y))
  y <- as.matrix(y)
  #
  y
}

fuzzy_vikor_internal <- function(x, weights, v = 0.5) {
  
  # Adjust weights
  
  # Number of outcomes
  n.out <- ncol(x) / 3
  #
  if (is.null(weights))
    weights <- rep(1 / n.out, n.out)
  else if ((!is.null(weights) & length(weights) != n.out))
    stop("Please provide weights for all outcomes.",
         call. = FALSE)
  #
  if (!is.null(weights) && sum(weights) != 1) {
    weights <- round(weights / sum(weights), digits = 2)
    #
    warning("Weights should always sum up to 1. To do so the given weights ",
            "are now standardized and the new weights are: ",
            paste(weights, collapse = ", "))
  }
  #
  # Create weights fuzzy weights
  #
  weights_fuzzy <- rep(c(weights / 3), each = 3)
  #
  pos.sol <- rep(1, length(weights_fuzzy))
  neg.sol <- rep(nrow(x), length(weights_fuzzy))
  
  # 1. S and R index
  
  # Calculate differences
  d <- matrix(nrow = nrow(x), ncol = ncol(x))
  for (i in seq(1, ncol(x), 3)) {
    # denominator = neg.sol[i + 2] - posI[i]
    d[, i] <- (x[, i] -pos.sol[i + 2])
    d[, i + 1] <- (x[, i + 1] -pos.sol[i + 1])
    d[, i + 2] <- (x[, i + 2] -pos.sol[i])
  }
  
  W <- diag(weights_fuzzy)
  
  NW <- d %*% W
  
  S <- matrix(nrow = nrow(x), ncol = 3)
  R <- matrix(nrow = nrow(x), ncol = 3)
  
  S[, 1] <- apply(NW[,seq(1, ncol(x), 3)], 1, sum)
  S[, 2] <- apply(NW[,seq(2, ncol(x), 3)], 1, sum)
  S[, 3] <- apply(NW[,seq(3, ncol(x), 3)], 1, sum)
  
  R[, 1] <- apply(NW[,seq(1, ncol(x), 3)], 1, max)
  R[, 2] <- apply(NW[,seq(2, ncol(x), 3)], 1, max)
  R[, 3] <- apply(NW[,seq(3, ncol(x), 3)], 1, max)
  
  # 2. Q index
  
  Q1 <- matrix(nrow = nrow(x), ncol = 3)
  Q2 <- matrix(nrow = nrow(x), ncol = 3)
  
  denominatorS <- max(S[, 3]) - min(S[, 1])
  Q1[, 1] <- (S[, 1] - min(S[, 3])) / denominatorS
  Q1[, 2] <- (S[, 2] - min(S[, 2])) / denominatorS
  Q1[, 3] <- (S[, 3] - min(S[, 1])) / denominatorS
  
  denominatorR = max(R[, 3]) - min(R[, 1])
  Q2[, 1] <- (R[, 1] - min(R[, 3])) / denominatorR
  Q2[, 2] <- (R[, 2] - min(R[, 2])) / denominatorR
  Q2[, 3] <- (R[, 3] - min(R[, 1])) / denominatorR
  
  if (v == 1)
    Q <- Q1
  else if (v == 0)
    Q <- Q2
  else
    Q <- v * Q1 + (1 - v) * Q2
  
  # 3. De-fuzzification using the method by Opricovic
  
  Def_S <- (S[, 1] + S[, 2] * 2 + S[, 3]) / 4
  Def_R <- (R[, 1] + R[, 2] * 2 + R[, 3]) / 4
  Def_Q <- (Q[, 1] + Q[, 2] * 2 + Q[, 3]) / 4
  
  # 4. Return results
  
  res <- data.frame(Q = Def_Q, S = Def_S, R = Def_R,
                    row.names = row.names(x)) %>% arrange(Q)
  #
  #class(res) <- c("vikor", class(res))
  attr(res, "performance.table") <- x
  attr(res, "fuzzy.weights") <- weights_fuzzy
  #
  res
}
