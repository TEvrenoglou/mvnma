updateversion <- function(x, verbose = FALSE) {
  
  update.0.3.0 <- update_needed(x$version, 0, 3, verbose)
  update.0.4.0 <- update_needed(x$version, 0, 4, verbose)
  
  #
  #  (1) Update mvnma object
  #
  if (inherits(x, "mvnma")) {
    #
    if (update.0.4.0) {
      x$outcomes <- replaceNULL(x$outcomes, attr(x, "outcomes"))
      x$trts <- replaceNULL(x$trts, attr(x, "trts"))
      x$n.domain <- replaceNULL(x$n.domain, attr(x, "n.domain"))
      x$reference.group <-
        replaceNULL(x$reference.group, attr(x, "reference.group"))
      x$level <- replaceNULL(x$level, attr(x, "level"))
      x$sm <- replaceNULL(x$sm, attr(x, "sm"))
      x$method.model <-
        replaceNULL(x$method.model, attr(x, "method.model"))
      x$model.code <- replaceNULL(x$model.code, attr(x, "model.code"))
      x$fit <- replaceNULL(x$fit, attr(x, "fit"))
      x$params <- replaceNULL(x$params, attr(x, "params"))
      x$varTE.missing <-
        replaceNULL(x$varTE.missing, attr(x, "varTE.missing"))
      x$scale.psi <- replaceNULL(x$scale.psi, attr(x, "scale.psi"))
      x$pair.objects <- replaceNULL(x$pair.objects, attr(x, "pair.objects"))
      x$n.chains <- replaceNULL(x$n.chains, attr(x, "n.chains"))
      x$n.iter <- replaceNULL(x$n.iter, attr(x, "n.iter"))
      x$n.burnin <- replaceNULL(x$n.burnin, attr(x, "n.burnin"))
      x$lower.rho <- replaceNULL(x$lower.rho, attr(x, "lower.rho"))
      x$upper.rho <- replaceNULL(x$upper.rho, attr(x, "upper.rho"))
      #
      x$version <- "0.4-0"
    }
    #
    return(x)
  }
  
  
  #
  #  (2) Update mvrank object
  #
  if (inherits(x, "mvrank")) {
    #
    if (update.0.3.0) {
      if (attr(x, "method") == "pBV")
        attr(x, "method") <- "pbest"
    }
    #
    if (update.0.4.0) {
      x$outcomes <- replaceNULL(x$outcomes, names(x))
      x$ranks <- replaceNULL(x$ranks, x[x$outcomes])
      x$trts <- replaceNULL(
        x$trts,
        sort(unique(unlist(lapply(x$ranks, function(z) z$treatment))))
      )
      x$ranks.shared <-
        replaceNULL(x$ranks.shared, attr(x, "ranks.common.trts"))
      x$trts.shared <- replaceNULL(x$trts.shared, attr(x, "common_trts"))
      x$method <- replaceNULL(x$method, attr(x, "method"))
      #
      x <- x[c("ranks", "trts", "ranks.shared", "trts.shared",
               "outcomes", "method")]
      x$version <- "0.4-0"
      #
      class(x) <- "mvrank"
    }
    #
    return(x)
  }
  
  #
  #  (3) Update vikor object
  #
  if (inherits(x, "vikor")) {
    #
    if (update.0.3.0) {
      if (attr(x, "ranking.method") == "pBV")
        attr(x, "ranking.method") <- "pbest"
    }
    #
    return(x)
  }
  
  x
}


update_needed <- function(version, major = 0, minor = 0,
                          verbose = FALSE) {
  if (is.null(version)) {
    version <- 0.1
    major.cur <- 0
    minor.cur <- 1
  }
  else {
    version <- unlist(strsplit(version, "-")[1])
    major.cur <-
      as.numeric(unlist(strsplit(version, ".", fixed = TRUE))[1])
    minor.cur <-
      as.numeric(unlist(strsplit(version, ".", fixed = TRUE))[2])
  }
  #
  res <-
    ifelse(major.cur < major,
           TRUE, ifelse(major.cur > major,
                        FALSE, minor.cur < minor))
  if (res & verbose)
    message(paste0("Update to mvnma, version ", major, ".", minor))
  #
  res
}
