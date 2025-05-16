mvnma_code <- function(n.outcomes, method, multiarm) {
  txt <-
    paste0("model {\n",
           "  # k = total number of studies\n",
           "  # k2 = number of two-arm studies\n",
           "  # n = number of treatments\n\n")
  #
  txt <- paste0(txt, code_control(n.outcomes, multiarm))
  #
  txt <- paste0(txt, "  \n")
  #
  # Files with variances, covariances and estimates
  #
  if (n.outcomes <= 5) {
    file <- paste0("mvnma_", n.outcomes, "_2arm.txt")
    path <- system.file("model", file, package = "mvnma")
    txt <- paste0(txt, paste(readLines(path), collapse = "\n"))
    #
    if (multiarm) {
      file <- paste0("mvnma_", n.outcomes, "_3arm.txt")
      path <- system.file("model", file, package = "mvnma")
      txt <- paste0(txt, "\n", paste(readLines(path), collapse = "\n"))
    }
    #
    txt <- paste0(txt, "\n")
  }
  else
    warning("R function mvnma_code() not implemented for more than five ",
            "outcomes.",
            call. = FALSE)
  #
  txt <- paste0(txt, code_means(n.outcomes, multiarm))
  #
  txt <- paste0(txt, "\n")
  #
  txt <- paste0(txt, code_priors(n.outcomes, method))
  #
  txt <- paste0(txt, "}\n")
  #
  txt
}

code_control <- function(n.outcomes, multiarm) {
  
  txt <-
    paste0(
      "  #\n",
      "  #\n",
      "  # (1) Set correlation to zero for studies not reporting one outcome\n",
      "  #\n",
      "  #\n\n")
  #
  txt <-
    paste0(txt,
           "  #\n",
           "  # Two-arm studies\n",
           "  #\n")
  #
  txt <- paste0(txt, "  for (i in 1:k2) {\n")
  #
  for (i in seq_len(n.outcomes - 1)) {
    for (j in (i + 1):n.outcomes) {
      txt <-
        paste0(txt,
               "    control", i, j, "[i] <- step(9999 - var", i,
               "[i]) * step(9999 - var", j, "[i])\n")
    }
  }
  #
  txt <- paste0(txt, "  }\n")
  #
  if (multiarm) {
    #
    txt <-
      paste0(txt,
             "  #\n",
             "  # Three-arm studies\n",
             "  #\n")
    #
    txt <- paste0(txt, "  for (i in ((k2 + 1):(2 * k - k2))) {\n")
    #
    for (i in seq_len(n.outcomes - 1)) {
      for (j in (i + 1):n.outcomes) {
        txt <-
          paste0(txt,
                 "    control", i, j, "[i] <- step(9999 - var", i,
                 "[i]) * step(9999 - var", j, "[i])\n")
      }
    }
    #
    txt <- paste0(txt, "  }\n")
  }
  #
  txt
}

code_means <- function(n.outcomes, multiarm) {
  
  txt <-
    paste0(
      "  #\n",
      "  #\n",
      "  # (3) Parameterization of the means\n",
      "  #\n",
      "  #\n\n")
  #
  txt <-
    paste0(txt,
           "  #\n",
           "  # Two-arm studies\n",
           "  #\n")
  #
  txt <- paste0(txt, "  for (i in 1:k2) {\n")
  #
  for (i in seq_len(n.outcomes)) {
    txt <-
      paste0(txt, "    mean[", n.outcomes, " * i",
             if (i != n.outcomes)
               paste0(" - ", n.outcomes - i)
             else
               strrep(" ", nchar(n.outcomes - 1) + 3),
             "] <- d", i, "[treat2[i]] - d", i, "[treat1[i]]\n")
  }
  #
  txt <- paste0(txt, "  }\n")
  #
  if (multiarm) {
    #
    txt <-
      paste0(txt,
             "  #\n",
             "  # Three-arm studies\n",
             "  #\n")
    #
    txt <- paste0(txt, "  for (i in 1:(k - k2)) {\n")
    #
    idx <- rep(seq_len(n.outcomes), 2)
    #
    for (i in seq_len(2 * n.outcomes)) {
      txt <-
        paste0(txt, "    mean[", n.outcomes, " * k2 + ", 2 * n.outcomes, " * i",
               if (i != 2 * n.outcomes)
                 paste0(" - ", 2 * n.outcomes - i)
               else
                 strrep(" ", nchar(2 * n.outcomes - 1) + 3),
               "] <- d",
               idx[i], "[treat", 2 + (i > n.outcomes), "[k2 + i]] - d",
               idx[i], "[treat1[k2 + i]]\n")
      #
      if (i != 2 * n.outcomes & idx[i] == n.outcomes)
        txt <- paste0(txt, "    #\n")
    }
    #
    txt <- paste0(txt, "  }\n")
  }
  txt
}

code_priors <- function(n.outcomes, method) {
  if (method == "standard")
    txt <- code_priors_standard(n.outcomes)
  else
    txt <- code_priors_dumouchel(n.outcomes)
  #
  txt <- paste0(txt, "  #\n")
  #
  txt <- paste0(txt, code_priors_psi(n.outcomes))
  #
  txt <- paste0(txt, "  #\n")
  #
  txt <- paste0(txt, code_priors_rho(n.outcomes))
  #
  txt
}

code_priors_standard <- function(n.outcomes) {
  
  txt <-
    paste0(
      "  #\n",
      "  #\n",
      "  # (4) Priors (standard model)\n",
      "  #\n",
      "  #\n\n")
  #
  txt <-
    paste0(txt,
           "  for (i in 1:(ref - 1)) {\n")
  #
  for (i in seq_len(n.outcomes))
    txt <- paste0(txt, "    d", i, "[i] ~ dnorm(0, 1e-03)\n")
  #
  txt <- paste0(txt, "  }\n")
  #
  txt <- paste0(txt, "  #\n")
  for (i in seq_len(n.outcomes))
    txt <- paste0(txt, "  d", i, "[ref] <- 0\n")
  txt <- paste0(txt, "  #\n")
  #
  txt <-
    paste0(txt,
           "  for (i in (ref + 1):n) {\n")
  #
  for (i in seq_len(n.outcomes))
    txt <- paste0(txt, "    d", i, "[i] ~ dnorm(0, 1e-03)\n")
  #
  txt <- paste0(txt, "  }\n")
  #
  txt
}

code_priors_dumouchel <- function(n.outcomes) {
  
  txt <-
    paste0(
      "  #\n",
      "  #\n",
      "  # (5) Priors (DuMouchel model)\n",
      "  #\n",
      "  #\n\n")
  #
  txt <-
    paste0(txt,
           "  for (i in 1:(ref - 1)) {\n",
           "    for (m in 1:", n.outcomes, ") {\n",
           "      meand[m, i] <- alpha[i] + gamma[m]\n",
           "      d[m, i] ~ dnorm(meand[m, i], prec.exp)\n",
           "    }\n",
           "    #\n")
  #
  for (i in seq_len(n.outcomes))
    txt <- paste0(txt, "    d", i, "[i] <- d[", i, ", i]\n")
  #
  txt <- paste0(txt, "  }\n")
  #
  txt <- paste0(txt, "  #\n")
  for (i in seq_len(n.outcomes))
    txt <- paste0(txt, "  d", i, "[ref] <- 0\n")
  txt <- paste0(txt, "  #\n")
  #
  txt <-
    paste0(txt,
           "  for (i in (ref + 1):n) {\n",
           "    for (m in 1:", n.outcomes, ") {\n",
           "      meand[m, i] <- alpha[i] + gamma[m]\n",
           "      d[m, i] ~ dnorm(meand[m, i], prec.exp)\n",
           "    }\n",
           "    #\n")
  #
  for (i in seq_len(n.outcomes))
    txt <- paste0(txt, "    d", i, "[i] <- d[", i, ", i]\n")
  #
  txt <- paste0(txt, "  }\n")
  #
  txt <- paste0(txt, "  #\n")
  #
  txt <-
    paste0(txt,
           "  for (m in 1:", n.outcomes, ") {\n",
           "  gamma[m] ~ dnorm(0, 1e-03)\n",
           "  }\n",
           "  #\n")
  #
  txt <-
    paste0(txt,
           "  for (i in 1:(ref - 1)) {\n",
           "    alpha[i] ~ dnorm(0, 1e-03)\n",
           "  }\n",
           "  #\n")
  #
  txt <-
    paste0(txt,
           "  for (i in (ref + 1):n) {\n",
           "    alpha[i] ~ dnorm(0, 1e-03)\n",
           "  }\n",
           "  #\n")
  #
  txt <-
    paste0(txt,
           "  prec.exp <- 1 / sigma.sq\n",
           "  sigma.sq <- sigma * sigma\n",
           "  sigma ~ dnorm(0, 1)T(0, )\n")
  #
  txt
}

code_priors_psi <- function(n.outcomes) {
  txt <- ""
  #
  for (i in seq_len(n.outcomes))
    txt <- paste0(txt, "  psi", i, ".sq <- psi", i, " * psi", i, "\n")
  #
  txt <- paste0(txt, "  #\n")
  #
  for (i in seq_len(n.outcomes))
    txt <- paste0(txt, "  psi", i, "  ~ dnorm(0, prec.psi", i, ")T(0, )\n")
  #
  txt
}

code_priors_rho <- function(n.outcomes) {
  txt <- ""
  #
  for (i in seq_len(choose(n.outcomes, 2)))
    txt <-
      paste0(txt, "  rho", i, " ~ dunif(lower.rho", i,
             ", upper.rho", i, ")\n")
  #
  txt
}
