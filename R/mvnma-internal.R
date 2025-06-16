catch <- function(argname, matchcall, data, encl)
  eval(matchcall[[match(argname, names(matchcall))]], data, enclos = encl)

multi_arm <- function(data) {
  
  # Get rid of warning "no visible binding for global variable"
  studlab <- treat2 <- NULL
  
  studies <- unique(data$studlab)    
  #
  r <- t <- E <- vector("list")
  
  for (i in seq_along(studies)) {
    r[[i]] <- data %>% filter(studlab == studies[i])
    #
    if (nrow(r[[i]]) > 2) {
      t[[i]] <- as.data.frame(table(r[[i]]$treat2))
      E[[i]] <- t[[i]]$Var1[which(t[[i]]$Freq == max(t[[i]]$Freq))]
      r[[i]] %<>% filter(treat2 == E[[i]])
    }
  }
  #
  # df harmonized in terms of treat2 
  #
  res <- list.rbind(r)
  #
  res
}

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
    trts.i <- c(unique(dat.i$label2), unique(dat.i$label1))
    #
    res[i, seq_along(trts.i)] <- trts.i
  }
  #
  res
}

'%!in%' <- function(x, y)
  !('%in%'(x, y))


#
# Helpers for mvdata()
#

create_data <- function(p, ...) {
  
  if (!(any(class(p) == "list")))
    stop("Argument 'p' must be a list of pairwise objects.", 
         call. = FALSE)
  
  # Get rid of warning "no visible binding for global variable"
  studlab <- TE <- seTE <- treat1 <- treat2 <- outcome <- n.arms <- NULL
  #
  dat1 <- dat2 <- studies <- vector("list")
  #
  n_outcomes <- length(p)
  #
  for (i in 1:n_outcomes) {
    p[[i]] %<>% select(studlab, TE, seTE, treat1, treat2) 
    p[[i]] <- multi_arm(p[[i]])
    p[[i]]$outcome <- i
    p[[i]] <- add_arms(p[[i]])
  }
  
  # Combine all p's
  comb_p <- list.rbind(p)
  
  for (i in 1:n_outcomes) {
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
    comb_p %<>%  arrange(outcome)
    #
    row.names(comb_p) <- NULL
  }
  else
    comb_p %<>% arrange(outcome)
  #
  res <- comb_p %>%
    group_by(studlab) %>%
    arrange(desc(treat1), .by_group = TRUE) %>% 
    arrange(n.arms)
  #
  res <- add_labels(res)
  #
  res
}

add_arms <- function(x, ...) {
  # Get rid of warning "no visible binding for global variable"
  studlab <- NULL
  #
  studies <- unique(x$studlab)
  #
  res <- vector("list")
  #
  for (i in seq_along(studies)) {
    res[[i]] <- x %>% filter(studlab == studies[i])
    #
    res[[i]]$n.arms <- length(unique(c(res[[i]]$treat1, res[[i]]$treat2)))
  }
  #
  res <- list.rbind(res)  
  #
  res
}

add_labels <- function(data) {
  all_treats <- unique(c(data$treat1, data$treat2))  
  #
  levels_treats <- as.data.frame(levels(as.factor(all_treats)))
  names(levels_treats) <- c("treat")
  levels_treats$level <- 1:nrow(levels_treats)  
  #
  data$label1 <- NA
  data$label2 <- NA
  #
  for (i in 1:nrow(data)) {
    for (j in 1:nrow(levels_treats)) {
      if (data$treat1[i] == levels_treats$treat[j]) {
        data$label1[i] = levels_treats$level[j]
      }
      #
      if (data$treat2[i] == levels_treats$treat[j]) {
        data$label2[i] = levels_treats$level[j]
      }
    }
  }  
  #
  data
}

make_jags_data <- function(dat) {
  
  # Get rid of warning "no visible binding for global variable"
  label <- NULL
  
  # Number of studies  
  Ns <- length(unique(dat$studlab))    
  
  # labtreat 
  
  labtreat1 <- cbind.data.frame(dat$treat1, dat$label1)
  labtreat2 <- cbind.data.frame(dat$treat2, dat$label2)
  #
  names(labtreat1) <- names(labtreat2) <- c("treat", "label")
  #
  labtreat <- rbind.data.frame(labtreat1, labtreat2) %>%
    distinct() %>% arrange(label)
  
  # NT
  NT <- length(labtreat$treat)
  
  # Number of outcomes
  n_outcomes <- length(unique(dat$outcome))
  
  # Arms per study
  arm_data <- dat[!duplicated(dat$studlab), ]
  
  # Number of two-arm studies
  two_arm <- length(which(arm_data$n.arms == 2))
  
  treat_data <- create_T(dat, max.arms = max(arm_data$n.arms))
  
  
  # Extract vector with treatment effects
  y <- dat$TE
  y <- ifelse(is.na(y), 0, y)
  #
  dat_out <- var <- treat_out <- vector("list")
  names_vec <- vector("character")
  
  for (i in seq_len(n_outcomes)) {
    dat_out[[i]] <- dat[dat$outcome == i, ]
    #
    var[[i]] <- dat_out[[i]]$seTE^2
    var[[i]] <- ifelse(is.na(var[[i]]), 10000, var[[i]])
    #
    names_vec[i] <- paste0("var", i)
    #
    dat_out[[i]] <- dat_out[[i]][complete.cases(dat_out[[i]]$TE), ]
    treat_out[[i]] <- unique(c(dat_out[[i]]$treat1, dat_out[[i]]$treat2))
  }
  #
  var_f <- list.cbind(var)
  var_f <- as.data.frame(var_f)
  names(var_f) <- names_vec
  #
  res <- list(y = y, var = var_f, T = treat_data,
              Ns = Ns, N2h = two_arm, NT = NT, 
              labtreat = labtreat, treat_out = treat_out)
  
  res
}

is.list.pairwise <- function(p, ...) {
  all_class <- sapply(p, class)
  #
  check <- vector()
  #
  for (i in 1:ncol(all_class))
    check[i] <- isTRUE("pairwise" %in% all_class[, i])  
  #
  pair <- ifelse(sum(check) >= 2, "pairwise", NA)
  #
  pair
}

gather_results <- function(x, outcomes, trts, reference.group,
                           level,treat_out, method, ...) {
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
  # Make the elements of the list "treat_out" to be every possible treatment if
  # argument 'method == "DM"'
  if (method == "DM") {
    for (i in 1:length(treat_out)) {
      treat_out[[i]] <- sort(trts)
    }
  }
  #
  basic <- dat_treat <- psi <- rho <- d <-
    TE.random <- seTE.random <- lower.random <- upper.random <-
    vector("list")
  # Get rid of warning "no visible binding for global variable"
  sd <- Rhat <- n.eff <- lower <- upper <- NULL
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
      select(any_of(treat_out[[i]]))
    #
    # Results for basic parameters
    #
    basic[[i]] <- res %>%
      filter(grepl(paste0(d.i, "["), rnames, fixed = TRUE))
    #
    row.names(basic[[i]]) <- trts
    basic[[i]] <- basic[[i]][which(row.names(basic[[i]]) %in% treat_out[[i]]), ]
    #
    basic[[i]] %<>%
      mutate(lower = basic[[i]]$`2.5%`,
               # apply(d[[i]], 2, quantile, probs = lower.level,
               #       na.rm = TRUE),
             upper = basic[[i]]$`97.5%`) %>% 
               # apply(d[[i]], 2, quantile, probs = upper.level,
               #       na.rm = TRUE)) %>%
      select(mean, sd, lower, upper, Rhat, n.eff)
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
    #
    # Matrices with all results
    #
    dmat.i <- d[[i]]
    #
    TE.random.i <- seTE.random.i <-
      lower.random.i <- upper.random.i <-
      matrix(NA, nrow = ncol(dmat.i), ncol = ncol(dmat.i),
             dimnames = list(treat_out[[i]], treat_out[[i]]))
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
  if (method == "DM")
    sigma <- res %>% filter(grepl("sigma", rnames))
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
  out1 <- list(basic_estimates = basic[[1]],
               heterogeneity = psi[1, ],
               TE.random = TE.random[[1]],
               seTE.random = seTE.random[[1]],
               lower.random = lower.random[[1]],
               upper.random = upper.random[[1]],
               samples = d[[1]])
  #
  out2 <- list(basic_estimates = basic[[2]],
               heterogeneity = psi[2, ],
               TE.random = TE.random[[2]],
               seTE.random = seTE.random[[2]],
               lower.random = lower.random[[2]],
               upper.random = upper.random[[2]],
               samples = d[[2]])
  #
  res <- list(out1, out2)
  names(res) <- c(outcomes[1], outcomes[2])
  #
  if (n.out >= 3) {
    out3 <- list(basic_estimates = basic[[3]],
                 heterogeneity = psi[3, ],
                 TE.random = TE.random[[3]],
                 seTE.random = seTE.random[[3]],
                 lower.random = lower.random[[3]],
                 upper.random = upper.random[[3]],
                 samples = d[[3]])
    #
    res[[3]] <- out3
    names(res)[3] <- outcomes[3]
  }
  #
  if (n.out >= 4) {
    out4 <- list(basic_estimates = basic[[4]],
                 heterogeneity = psi[4, ],
                 TE.random = TE.random[[4]],
                 seTE.random = seTE.random[[4]],
                 lower.random = lower.random[[4]],
                 upper.random = upper.random[[4]],
                 samples = d[[4]])
    #
    res[[4]] <- out4
    names(res)[4] <- outcomes[4]
  }
  #
  if (n.out >= 5) {
    out5 <- list(basic_estimates = basic[[5]],
                 heterogeneity = psi[5, ],
                 TE.random = TE.random[[5]],
                 seTE.random = seTE.random[[5]],
                 lower.random = lower.random[[5]],
                 upper.random = upper.random[[5]],
                 samples = d[[5]])
    #
    res[[5]] <- out5
    names(res)[5] <- outcomes[5]
  }
  #
  res[[length(res) + 1]] <- cor
  names(res)[length(res)] <- "cor"
  #
  if (method == "DM") {
    res[[length(res) + 1]] <- sigma
    names(res)[length(res)] <- "sigma"  
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
