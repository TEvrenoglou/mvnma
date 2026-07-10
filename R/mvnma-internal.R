mvdata <- function(x)
  make_jags_data(create_data(x))

create_T <- function(data, max.arms) {
  # Get rid of warning "no visible binding for global variable"
  studlab <- NULL
  #
  studies <- unique(data$studlab)  
  #
  res <- matrix(NA, nrow = length(studies), ncol = max.arms)
  #
  for (i in seq_along(studies)) {
    dat.i <- data %>% filter(studlab == studies[i])
    #
    trts.i <- sort(unique(c(dat.i$id1, dat.i$id2)))
    #
    res[i, seq_along(trts.i)] <- trts.i
  }
  #
  rownames(res) <- studies
  colnames(res) <- paste0("id", seq_len(ncol(res)))
  #
  res
}

create_data <- function(p) {
  
  # Get rid of warning "no visible binding for global variable"
  studlab <- TE <- seTE <- treat1 <- treat2 <- outcome <- n.arms <- 
    ref_study <- NULL
  #
  dat1 <- dat2 <- studies <- vector("list")
  #
  n.out <- length(p)
  #
  for (i in seq_len(n.out)) {
    #
    # Rename variable names specified in argument 'varnames' of pairwise()
    #
    vn <- attr(p[[i]], "varnames")
    if (!identical(vn, c("TE", "seTE")))
      p[[i]] %<>% rename_with(~ c("TE", "seTE"), all_of(vn))
    #
    p[[i]] %<>% select(studlab, TE, seTE, treat1, treat2)
    #
    # Outcome number
    #
    p[[i]]$outcome <- i
    #
    # Add number of treatment arms
    #
    tab_narms <- table(p[[i]]$studlab)
    dat_narms <-
      data.frame(studlab = names(tab_narms),
                 n.arms = (1 + sqrt(8 * as.vector(tab_narms) + 1)) / 2)
    p[[i]] <- merge(p[[i]], dat_narms, by = "studlab",
                    all.x = TRUE, all.y = FALSE)
    #
    # For multi-arm studies only keep comparisons with
    # - the reference treatment or
    # - the first alpha-numeric treatment
    #
    p[[i]] %<>%
      group_by(studlab) %>%
      mutate(ref_study = names(which.max(table(treat2)))) %>%
      filter(treat2 == ref_study) %>%
      select(-ref_study) %>%
      as.data.frame()
  }
  #
  # Combine all pairwise objects
  #
  comb_p <- bind_rows(p)
  #
  # Add missing comparisons for three-arm studies
  #
  twoarm <- comb_p %>%
    filter(n.arms == 2) %>%
    select(studlab) %>% unique() %>% pull()
  #
  threearm <- comb_p %>%
    filter(n.arms == 3) %>%
    select(studlab) %>% unique() %>% pull()
  #
  both <- twoarm[twoarm %in% threearm]
  #
  for (i in both) {
    drugs.i <- comb_p %>% filter(studlab == i, n.arms == 3) %>%
      select(studlab, treat1, treat2) %>% distinct()
    #
    dat2.i <- comb_p %>% filter(studlab == i, n.arms == 2)
    #
    for (j in unique(dat2.i$outcome)) {
      drugs.ij <- drugs.i %>% mutate(outcome = j, n.arms = 3)
      dat2.ij <- dat2.i %>% filter(outcome == j)
      #
      drop.j <- vector()
      for (k in seq_len(nrow(drugs.ij))) {
        drop.j[k] <- any(dat2.ij$studlab == i & dat2.ij$outcome == j &
                           dat2.ij$treat1 == drugs.ij$treat1[k] &
                           dat2.ij$treat2 == drugs.ij$treat2[k])
      }
      #
      comb_p <- bind_rows(comb_p, drugs.ij %>% filter(!drop.j))
      comb_p$n.arms[comb_p$studlab == i & comb_p$outcome == j] <- 3
    }
  }
  
  for (i in seq_len(n.out)) {
    dat1[[i]] <- comb_p %>% filter(outcome == i)
    dat2[[i]] <- comb_p %>% filter(outcome != i) 
    #
    dat2[[i]]$TE <- dat2[[i]]$seTE <- NA
    dat2[[i]]$new_outcome <- i
    #
    studies[[i]] <- unique(which(dat2[[i]]$studlab %!in% dat1[[i]]$studlab))
    #
    if (length(studies[[i]]) > 0) {
      dat2[[i]] <- dat2[[i]][studies[[i]], ]
      dat2[[i]]$outcome <- dat2[[i]]$new_outcome
      dat2[[i]]$new_outcome <- NULL
      #
      row.names(dat2[[i]]) <- NULL
    }
    else
      dat2[[i]] <- list()
  }
  #
  if (length(dat2) != 0) {
    dat2 <- list.rbind(dat2)
    #
    comb_p <- rbind.data.frame(comb_p, dat2)
    comb_p <- comb_p %>% arrange(outcome)
    #
    row.names(comb_p) <- NULL
  }
  else
    comb_p %<>% arrange(outcome)
  #
  comb_p %<>% distinct()
  #
  res <- comb_p %>%
    distinct() %>%
    arrange(outcome) %>%
    group_by(studlab) %>%
    arrange(desc(treat1), .by_group = TRUE) %>%
    arrange(n.arms)
  #
  res <- as.data.frame(res)
  rownames(res) <- seq_len(nrow(res))
  #
  # Add IDs for treatment names
  #
  trts <- sort(unique(c(res$treat1, res$treat2)))
  res$id1 <- as.integer(factor(res$treat1, levels = trts))
  res$id2 <- as.integer(factor(res$treat2, levels = trts))
  #
  res
}

make_jags_data <- function(dat) {
  
  # Get rid of warning "no visible binding for global variable"
  outcome <- seTE <- studlab <- TE <- var <- NULL
  
  # Number of outcomes
  n.out <- length(unique(dat$outcome))
  
  # Number of studies  
  k <- length(unique(dat$studlab))    
  
  # Treatments
  trts <- sort(unique(c(dat$treat1, dat$treat2)))
  
  # Number of treatments
  n <- length(trts)
  
  # Arms per study
  arm_data <- dat[!duplicated(dat$studlab), ]
  
  # Number of two-arm studies
  k2 <- sum(arm_data$n.arms == 2, na.rm = TRUE)
  
  treat_data <- create_T(dat, max.arms = max(arm_data$n.arms))
  
  # Extract vector with treatment effects
  #
  y <- dat$TE
  y <- ifelse(is.na(y), 0, y)
  #
  trts.list <- vector("list", n.out)
  dat_var <- NULL
  #
  for (i in seq_len(n.out)) {
    dat_var_i <- subset(dat, outcome == i)
    #
    if (is.null(dat_var))
      dat_var <- data.frame(studlab = dat_var_i$studlab)
    #
    dat_var_i_complete <- subset(dat_var_i, complete.cases(TE))
    trts.list[[i]] <- with(dat_var_i_complete, sort(unique(c(treat1, treat2))))
    #
    dat_var_i %<>% mutate(var = seTE^2) %>% select(var)
    #
    names(dat_var_i) <- paste0("var", i)
    #
    dat_var <- cbind(dat_var, dat_var_i)
  }
  #
  res <- list(y = y, var = dat_var, treatments = treat_data,
              k = k, k2 = k2, n = n, 
              trts = trts, trts.list = trts.list)
  #
  res
}

gather_results <- function(x, outcomes, trts, reference.group,
                           level, trts.list, method,
                           n.domain, ...) {
  
  res <- as.data.frame(x$BUGSoutput$summary)
  samples <- x$BUGSoutput$sims.list
  #
  n.chain <- x$BUGSoutput$n.chain
  n.thin <- x$BUGSoutput$n.thin
  n.burnin <- x$BUGSoutput$n.burnin
  #
  ref <- which(trts == reference.group)
  rnames <- rownames(res)
  #
  n.out <- length(outcomes)
  #
  lower.level <- (1 - level) / 2
  upper.level <- 1 - (1 - level) / 2
  #
  # Make the elements of the list "trts.list" to be every possible treatment if
  # argument 'method == "DM"'
  if (method == "DM") {
    for (i in seq_along(trts.list))
      trts.list[[i]] <- sort(trts)
  }
  #
  basic <- dat_treat <- psi <- rho <- d <-
    TE.random <- seTE.random <- lower.random <- upper.random <-
    vector("list")
  # Get rid of warning "no visible binding for global variable"
  sd <- Rhat <- n.eff <- lower <- upper <- out_num <- NULL
  #
  for (i in seq_len(n.out)) {
    d.i <- paste0("d", i)
    #
    # Samples
    #
    d[[i]] <- samples[[d.i]]
    colnames(d[[i]]) <- trts
    rownames(d[[i]]) <- seq_len(nrow(d[[i]]))
    #
    d[[i]] <- as.data.frame(d[[i]]) %>% 
      select(any_of(trts.list[[i]]))
    #
    # Results for basic parameters
    #
    basic[[i]] <- res %>%
      filter(grepl(paste0(d.i, "["), rnames, fixed = TRUE)) %>%
      mutate(lower = NA, upper = NA)
    #
    row.names(basic[[i]]) <- trts
    basic[[i]] <- basic[[i]][which(row.names(basic[[i]]) %in% trts.list[[i]]), ]
    ##
    for (j in rownames(basic[[i]])) {
      basic[[i]][j, c("lower", "upper")] <-
        quantile(d[[i]][[j]], probs = c(lower.level, upper.level), na.rm = TRUE)
    }
    #
    basic[[i]] %<>% select(mean, sd, lower, upper, Rhat, n.eff)
    #
    basic[[i]][reference.group, ] <- NA
    #
    # psi's (similar to tau)
    #
    psi[[i]] <- res %>% filter(grepl("psi", rnames)) %>%
      select(mean, sd, "2.5%", "97.5%", Rhat, n.eff) %>%
      rename(psi = mean, lower = "2.5%", upper = "97.5%")
    #
    # rho
    #
    rho[[i]] <- res %>% filter(grepl("rho", rnames))
    # ensure the order of outcomes
    rho[[i]]$out_num <- as.numeric(gsub("rho", "", rownames(rho[[i]])))
    #
    rho[[i]] %<>% arrange(out_num)
    rho[[i]]$out_num <- NULL
    
    #
    # Matrices with all results
    #
    dmat.i <- d[[i]]
    #
    TE.random.i <- seTE.random.i <-
      lower.random.i <- upper.random.i <-
      matrix(NA, nrow = ncol(dmat.i), ncol = ncol(dmat.i),
             dimnames = list(trts.list[[i]], trts.list[[i]]))
    #
    for (j in seq_len(ncol(dmat.i)))
      for (k in seq_len(ncol(dmat.i)))
        if (j != k)
          TE.random.i[j, k] <- mean(dmat.i[, j] - dmat.i[, k])
    #
    for (j in seq_len(ncol(dmat.i)))
      for (k in seq_len(ncol(dmat.i)))
        if (j != k)
          seTE.random.i[j, k] <- sd(dmat.i[, j] - dmat.i[, k])
    #
    for (j in seq_len(ncol(dmat.i)))
      for (k in seq_len(ncol(dmat.i)))
        if (j != k)
          lower.random.i[j, k] <- quantile(dmat.i[, j] - dmat.i[, k],
                                           lower.level)
    #
    for (j in seq_len(ncol(dmat.i)))
      for (k in seq_len(ncol(dmat.i)))
        if (j != k)
          upper.random.i[j, k] <- quantile(dmat.i[, j] - dmat.i[, k],
                                           upper.level)
    #
    diag(TE.random.i) <- diag(seTE.random.i) <-
      diag(lower.random.i) <- diag(upper.random.i) <- 0
    #
    TE.random[[i]] <- TE.random.i
    seTE.random[[i]] <- seTE.random.i
    lower.random[[i]] <- lower.random.i
    upper.random[[i]] <- upper.random.i
  }
  #
  # Prepare output
  #
  cor <- rho[[1]]
  #
  cor %<>% select(mean, sd, "2.5%", "97.5%", Rhat, n.eff) %>%
    rename(lower = "2.5%", upper = "97.5%")
  #
  psi <- psi[[1]]
  row.names(psi) <- outcomes
  #
  if (method == "DM") {
    if (is.null(n.domain)) {
      sigma <- res %>% filter(grepl("sigma", rnames))
    }
    else {
      sigma1 <- res %>% filter(grepl("sigma1", rnames))
      sigma2 <- res %>% filter(grepl("sigma2", rnames))
    }
  }
  #
  # Create row.names for cor
  #
  r1 <- t(combn(seq_along(outcomes), 2))
  r.names <- vector("numeric", nrow(cor))
  #
  for (i in seq_along(r.names)) {
    for (j in seq_len(ncol(r1))) {
      r.names[i] <- paste(outcomes[r1[i, 1]], outcomes[r1[i, 2]], sep = "/")    
    }
  }
  #
  row.names(cor) <- r.names
  #
  res <- vector("list", n.out)
  #
  for (i in seq_len(n.out)) {
    res[[i]] <- list(basic_estimates = basic[[i]],
                     heterogeneity = psi[i, ],
                     TE.random = TE.random[[i]],
                     seTE.random = seTE.random[[i]],
                     lower.random = lower.random[[i]],
                     upper.random = upper.random[[i]],
                     samples = d[[i]])
    #
    names(res)[i] <- outcomes[i]
  }
  #
  res[[length(res) + 1]] <- cor
  names(res)[length(res)] <- "cor"
  #
  if (method == "DM") {
    if (is.null(n.domain)) {
      res[[length(res) + 1]] <- sigma
      names(res)[length(res)] <- "sigma"
    }
    else {
      res[[length(res) + 1]] <- sigma1
      names(res)[length(res)] <- "sigma1"
      #
      res[[length(res) + 1]] <- sigma2
      names(res)[length(res)] <- "sigma2"
    }
  }
  #
  res
}

setchar <- function(x, val, text, list = FALSE, name = NULL,
                    stop.at.error = TRUE, addtext = "",
                    return.NULL = TRUE, nchar.equal = FALSE,
                    setNA = FALSE, pre = "") {
  val <- unique(val)
  #
  if (is.null(name))
    name <- deparse(substitute(x))
  nval <- length(val)
  #
  if (is.numeric(x)) {
    numeric.x <- TRUE
    idx <- x
    idx[idx < 1] <- NA
    idx[idx >= nval + 1] <- NA
  }
  else {
    numeric.x <- FALSE
    #
    if (length(unique(tolower(x))) != length(unique(x)) |
        length(unique(tolower(val))) != length(unique(val)))
      idx <- charmatch(x, val, nomatch = NA)
    else
      idx <- charmatch(tolower(x), tolower(val), nomatch = NA)
  }
  #
  if ((anyNA(idx) || any(idx == 0)) && !setNA) {
    if (list)
      first <- "List element '"
    else
      first <- "Argument '"
    #
    if (missing(text)) {
      if (numeric.x) {
        if (nval == 1)
          vlist <- "1"
        else if (nval == 2)
          vlist <- "1 or 2"
        else
          vlist <- paste("between 1 and", nval)
      }
      else {
        if (nval == 1)
          vlist <- paste0('"', val, '"')
        else if (nval == 2)
          vlist <- paste0('"', val, '"', collapse = " or ")
        else
          vlist <- paste0(paste0('"', val[-nval], '"', collapse = ", "),
                          ', or ', '"', val[nval], '"')
      }
      #
      if (stop.at.error)
        stop(first, name, "' must be ", pre,
             vlist, addtext, ".", call. = FALSE)
      else {
        if (return.NULL)
          return(NULL)
        else
          return(x)
      }
    }
    else {
      if (stop.at.error)
        stop(first, name, "' ", text, ".", call. = FALSE)
      else {
        if (return.NULL)
          return(NULL)
        else
          return(x)
      }
    }
  }
  #
  if (is.null(x))
    return(NULL)
  else
    res <- val[idx]
  #
  if (nchar.equal && nchar(res) != nchar(x))
    res <- x
  #
  res
}

chkclass <- function(x, class, name = NULL) {
  #
  # Check class of R object
  #
  if (is.null(name))
    name <- deparse(substitute(x))
  #
  n.class <- length(class)
  if (n.class == 1)
    text.class <- paste0('"', class, '"')
  else if (n.class == 2)
    text.class <- paste0('"', class, '"', collapse = " or ")
  else
    text.class <- paste0(paste0('"', class[-n.class], '"', collapse = ", "),
                         ', or ', '"', class[n.class], '"')
  #
  if (!inherits(x, class))
    stop("Argument '", name, "' must be an object of class ", text.class, ".",
         call. = FALSE)
  #
  invisible(NULL)
}

chkchar <- function(x, length = 0, name = NULL, nchar = NULL, single = FALSE,
                    NULL.ok = FALSE) {
  if (is.null(x) & NULL.ok)
    return(invisible(NULL))
  #
  if (!missing(single) && single)
    length <- 1
  if (is.null(name))
    name <- deparse(substitute(x))
  #
  if (length && length(x) != length) {
    if (length == 1)
      stop("Argument '", name, "' must be a character string.",
           call. = FALSE)
    else
      stop("Argument '", name, "' must be a character vector of length ",
           length, ".",
           call. = FALSE)
  }
  #
  if (length == 1) {
    if (!is.null(nchar) && !(nchar(x) %in% nchar))
      if (length(nchar) == 1 && nchar == 1)
        stop("Argument '", name, "' must be a single character.",
             call. = FALSE)
    else
      stop("Argument '", name, "' must be a character string of length ",
           if (length(nchar) == 2)
             paste0(nchar, collapse = " or ")
           else
             paste0(nchar, collapse = ", "),
           ".",
           call. = FALSE)
  }
  #
  if (!is.character(x) & !is.numeric(x))
    stop("Argument '", name, "' must be a character vector.")
  else {
    if (!is.null(nchar) & any(!(nchar(x) %in% nchar)))
      if (length(nchar) == 1 && nchar == 1)
        stop("Argument '", name, "' must be a vector of single characters.",
             call. = FALSE)
    else
      stop("Argument '", name, "' must be a character vector where ",
           "each element has ",
           if (length(nchar) == 2)
             paste0(nchar, collapse = " or ")
           else
             paste0(nchar, collapse = ", "),
           " characters.",
           call. = FALSE)
  }
  #
  invisible(NULL)
}

chklevel <- function(x, length = 0, ci = TRUE, name = NULL, single = FALSE) {
  if (!missing(single) && single)
    length <- 1
  #
  # Check for levels of confidence interval / contour level
  #
  if (is.null(name))
    name <- deparse(substitute(x))
  if (ci)
    "level for confidence interval (range: 0-1)"
  else
    "contour levels (range: 0-1)"
  #
  if (!is.numeric(x))
    if (length && length(x) != length)
      stop("Argument '", name, "' must be a numeric of length ", length, ".",
           call. = FALSE)
  else
    stop("Argument '", name, "' must be numeric.",
         call. = FALSE)
  #
  if (length && length(x) != length)
    stop("Argument '", name, "' must be a numeric of length ", length, ".",
         call. = FALSE)
  #
  if (any(x <= 0, na.rm = TRUE) | any(x >= 1, na.rm = TRUE))
    stop("Argument '", name, "' must be a numeric between 0 and 1.",
         call. = FALSE)
  #
  invisible(NULL)
}

chklogical <- function(x, name = NULL, text = "") {
  #
  # Check whether argument is logical
  #
  if (is.null(name))
    name <- deparse(substitute(x))
  #
  if (is.numeric(x))
    x <- as.logical(x)
  #
  if (length(x) !=  1 || !is.logical(x) || is.na(x))
    stop("Argument '", name, "' must be a logical",
         if (text != "") " ", text, ".", call. = FALSE)
  #
  invisible(NULL)
}

chknull <- function(x, name = NULL) {
  #
  # Check whether argument is NULL
  #
  if (is.null(name))
    name <- deparse(substitute(x))
  #
  if (is.null(x))
    stop("Argument '", name, "' is NULL.", call. = FALSE)
  #
  invisible(NULL)
}

chknumeric <- function(x, min, max, zero = FALSE, length = 0,
                       name = NULL, single = FALSE, integer = FALSE,
                       NA.ok = TRUE) {
  if (!missing(single) && single)
    length <- 1
  #
  # Check numeric variable
  #
  if (is.null(name))
    name <- deparse(substitute(x))
  #
  if (NA.ok)
    x <- x[!is.na(x)]
  else if (anyNA(x))
    stop("Missing values not allowed in argument '", name, "'.",
         call. = FALSE)
  #
  if (length(x) == 0)
    return(NULL)
  #
  if (!is.numeric(x))
    stop("Non-numeric value for argument '", name, "'.",
         call. = FALSE)
  #
  if (length && length(x) != length)
    stop("Argument '", name, "' must be a numeric of length ", length, ".",
         call. = FALSE)
  #
  if (!missing(min) & missing(max)) {
    if (zero & min == 0 & any(x <= min, na.rm = TRUE))
      stop("Argument '", name, "' must be positive.",
           call. = FALSE)
    else if (any(x < min, na.rm = TRUE))
      stop("Argument '", name, "' must be larger equal ",
           min, ".", call. = FALSE)
  }
  #
  if (missing(min) & !missing(max)) {
    if (zero & max == 0 & any(x >= max, na.rm = TRUE))
      stop("Argument '", name, "' must be negative.",
           call. = FALSE)
    else if (any(x > max, na.rm = TRUE))
      stop("Argument '", name, "' must be smaller equal ",
           min, ".", call. = FALSE)
  }
  #
  if ((!missing(min) & !missing(max)) &&
      (any(x < min, na.rm = TRUE) | any(x > max, na.rm = TRUE)))
    stop("Argument '", name, "' must be between ",
         min, " and ", max, ".", call. = FALSE)
  #
  if (integer && any(!is_wholenumber(x))) {
    if (length(x) == 1)
      stop("Argument '", name, "' must be an integer.",
           call. = FALSE)
    else
      stop("Argument '", name, "' may only contain integers.",
           call. = FALSE)
  }
  #
  invisible(NULL)
}

is_wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
  if (is.numeric(x))
    res <- abs(x - round(x)) < tol
  else
    res <- NA
  #
  res
}

formatN <- function(x, digits = 2, text.NA = "--", big.mark = "",
                    format.whole.numbers = TRUE,
                    monospaced = FALSE) {
  
  outdec <- options()$OutDec
  
  if (!monospaced) {
    if (format.whole.numbers) {
      res <- format(ifelse(is.na(x),
                           text.NA,
                           formatC(x, decimal.mark = outdec,
                                   format = "f", digits = digits,
                                   big.mark = big.mark)
      )
      )
    }
    else {
      res <- format(ifelse(is.na(x),
                           text.NA,
                           ifelse(is_wholenumber(x),
                                  x,
                                  formatC(x, decimal.mark = outdec,
                                          format = "f", digits = digits,
                                          big.mark = big.mark)
                           )
      )
      )
    }
  }
  else {
    x <- round(x, digits)
    res <- ifelse(is.na(x),
                  text.NA,
                  format(x, decimal.mark = outdec, big.mark = big.mark))
  }
  #
  res <- rmSpace(res, end = TRUE)
  #
  res
}

rmSpace <- function(x, end = FALSE, pat = " ") {
  
  if (!end) {
    while (any(substring(x, 1, 1) == pat, na.rm = TRUE)) {
      sel <- substring(x, 1, 1) == pat
      x[sel] <- substring(x[sel], 2)
    }
  }
  else {
    last <- nchar(x)
    
    while (any(substring(x, last, last) == pat, na.rm = TRUE)) {
      sel <- substring(x, last, last) == pat
      x[sel] <- substring(x[sel], 1, last[sel] - 1)
      last <- nchar(x)
    }
  }
  
  x
}

formatCI <- function(lower, upper,
                     bracket.left = gs("CIbracket"),
                     separator = gs("CIseparator"),
                     bracket.right,
                     justify.lower = "right",
                     justify.upper = justify.lower,
                     lower.blank = gs("CIlower.blank"),
                     upper.blank = gs("CIupper.blank"),
                     ...) {
  
  # Change layout of CIs
  #
  chkchar(bracket.left, length = 1)
  chkchar(separator, length = 1)
  if (!missing(bracket.right))
    chkchar(bracket.right, length = 1)
  #
  if (missing(bracket.left)) {
    bracktype <- setchar(bracket.left, c("[", "(", "{", ""))
    #
    if (bracktype == "[") {
      bracketLeft <- "["
      bracketRight <- "]"
    }
    else if (bracktype == "(") {
      bracketLeft <- "("
      bracketRight <- ")"
    }
    else if (bracktype == "{") {
      bracketLeft <- "{"
      bracketRight <- "}"
    }
    else if (bracktype == "") {
      bracketLeft <- ""
      bracketRight <- ""
    }
    #
    bracket.left <- bracketLeft
  }
  #
  if (missing(bracket.right))
    bracket.right <- bracketRight
  
  format.lower <- format(lower, justify = justify.lower)
  format.upper <- format(upper, justify = justify.upper)
  #
  if (!lower.blank)
    format.lower <- rmSpace(format.lower)
  if (!upper.blank)
    format.upper <- rmSpace(format.upper)
  #
  if (separator == "-")
    format.upper <-
    paste0(ifelse(substring(format.upper, 1, 1) == "-", " ", ""),
           format.upper)
  #
  res <- ifelse(lower != "NA" & upper != "NA",
                paste0(bracket.left,
                       format.lower,
                       separator,
                       format.upper,
                       bracket.right),
                "")
  #
  res
}

catch <- function(argname, matchcall, data, encl)
  eval(matchcall[[match(argname, names(matchcall))]], data, enclos = encl)

'%!in%' <- function(x, y)
  !('%in%'(x, y))
