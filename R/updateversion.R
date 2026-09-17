updateversion <- function(x, verbose = FALSE) {
  
  update.0.3.0 <- update_needed(x$version, 0, 3, verbose)
  
  #
  #  (1) Update mvrank object
  #
  if (inherits(x, "mvrank")) {
    ##
    if (update.0.3.0) {
      if (attr(x, "method") == "pBV")
        attr(x, "method") <- "pbest"
      #
      return(x)
    }
  }
  
  #
  #  (1) Update vikor object
  #
  if (inherits(x, "vikor")) {
    ##
    if (update.0.3.0) {
      if (attr(x, "ranking.method") == "pBV")
        attr(x, "ranking.method") <- "pbest"
      #
      return(x)
    }
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
  ##
  res <-
    ifelse(major.cur < major,
           TRUE, ifelse(major.cur > major,
                        FALSE, minor.cur < minor))
  if (res & verbose)
    message(paste0("Update to mvnma, version ", major, ".", minor))
  ##
  res
}
