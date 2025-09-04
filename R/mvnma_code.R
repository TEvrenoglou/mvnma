mvnma_code <- function(n.out, method, multiarm) {
  chknumeric(n.out, min = 2, length = 1)
  method <- setchar(method, c("standard", "DM"))
  chklogical(multiarm)
  #
  txt <-
    paste0("model {\n",
           "  # k = total number of studies\n",
           "  # k2 = number of two-arm studies\n",
           "  # n = number of treatments\n\n")
  #
  txt <- paste0(txt, code_control(n.out, multiarm))
  #
  txt <- paste0(txt, "  \n")
  #
  # Files with variances, covariances and estimates
  #
  txt <- paste0(txt, code_covar_ests(n.out, multiarm))
  #
  txt <- paste0(txt, code_means(n.out, multiarm))
  #
  txt <- paste0(txt, "\n")
  #
  txt <- paste0(txt, code_priors(n.out, method))
  #
  txt <- paste0(txt, "}\n")
  #
  txt
}

code_control <- function(n.out, multiarm) {
  
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
  for (i in seq_len(n.out - 1)) {
    for (j in (i + 1):n.out) {
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
    for (i in seq_len(n.out - 1)) {
      for (j in (i + 1):n.out) {
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

code_covar_ests <- function(n.out, multiarm) {
  
  txt <-
    paste0(
      "  #\n",
      "  #\n",
      "  # (2) Variances, covariances and estimates\n",
      "  #\n",
      "  #\n  \n")
  #
  txt <-
    paste0(txt,
           "  #\n",
           "  # Two-arm studies\n",
           "  #\n")
  #
  txt <- paste0(txt, "  for (i in 1:k2) {\n")
  #
  txt <-
    paste0(txt,
           "    #\n",
           "    # Variances\n",
           "    #\n")
  #
  for (i in seq_len(n.out)) {
    txt <-
      paste0(txt,
             "    S2[i, ", i, ", ", i, "] <- var", i, "[i] + ",
             "psi", i, ".sq\n")
  }
  #
  txt <-
    paste0(txt,
           "    #\n",
           "    # Covariances\n",
           "    #\n")
  #
  r <- 0
  #
  for (i in seq_len(n.out - 1)) {
    for (j in (i + 1):n.out) {
      r <- r + 1
      #
      txt <-
        paste0(txt,
               "    S2[i, ", i, ", ", j, "] <- ",
               "sqrt(S2[i, ", i, ", ", i, "]) * ",
               "sqrt(S2[i, ", j, ", ", j, "]) * ",
               "control", i, j, "[i] * rho", r, "\n")
    }
  }
  #
  txt <- paste0(txt, "    #\n")
  #
  for (i in seq_len(n.out - 1)) {
    for (j in (i + 1):n.out) {
      txt <-
        paste0(txt,
               "    S2[i, ", j, ", ", i, "] <- S2[i, ", i, ", ", j, "]\n")
    }
  }
  #
  txt <-
    paste0(txt,
           "    #\n",
           "    # Estimates\n",
           "    #\n")
  #
  txt <-
    paste0(txt,
           "    y[(", n.out, " * i - ", n.out - 1,
           "):(", n.out, " * i)] ~ dmnorm.vcov(mean[(",
           n.out, " * i - ", n.out - 1, "):(", n.out,
           " * i)], S2[i, , ])\n")
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
    idx <- rep(seq_len(n.out), 2)
    #
    txt <-
      paste0(txt,
             "    #\n",
             "    # Variances\n",
             "    #\n")
    #
    for (i in seq_len(2 * n.out)) {
      txt <-
        paste0(txt, "    S3[i, ", i, ", ", i, "] <- var", idx[i],
               "[k2 + 2 * i",
               if (i <= n.out) " - 1",
               "] + psi", idx[i], ".sq\n")
      #
      if (idx[i] == n.out)
        txt <- paste0(txt, "    #\n")
    }
    txt <-
      paste0(txt,
             "    # Covariances\n",
             "    #\n")
    #
    r <- 0
    #
    for (i in seq_len(n.out - 1)) {
      for (j in (i + 1):n.out) {
        r <- r + 1
        #
        txt <-
          paste0(txt,
                 "    S3[i, ", i, ", ", j, "] <- ",
                 "sqrt(S3[i, ", i, ", ", i, "]) * ",
                 "sqrt(S3[i, ", j, ", ", j, "]) * ",
                 "control", i, j, "[i] * rho", r, "\n")
      }
    }
    #
    txt <- paste0(txt, "    #\n")
    #
    r <- 0
    #
    for (i in n.out + seq_len(n.out - 1)) {
      for (j in (i + 1):(2 * n.out)) {
        r <- r + 1
        #
        txt <-
          paste0(txt,
                 "    S3[i, ", i, ", ", j, "] <- ",
                 "sqrt(S3[i, ", i, ", ", i, "]) * ",
                 "sqrt(S3[i, ", j, ", ", j, "]) * ",
                 "control", i - n.out, j - n.out,
                 "[i] * rho", r, "\n")
      }
    }
    #
    txt <- paste0(txt, "    #\n")
    #
    for (i in seq_len(n.out)) {
      txt <-
        paste0(txt,
               "    S3[i, ", i, ", ", i + n.out, "] <- 0.5 * ",
               "sqrt(S3[i, ", i, ", ", i, "]) * ",
               "sqrt(S3[i, ", i + n.out, ", ", i + n.out, "])\n")
    }
    #
    txt <- paste0(txt, "    #\n")
    #
    r <- 0
    #
    for (i in seq_len(n.out - 1)) {
      for (j in (i + 1):n.out) {
        l <- j + n.out
        r <- r + 1
        #
        txt <-
          paste0(txt,
                 "    S3[i, ", i, ", ", l, "] <- 0.5 * ",
                 "sqrt(S3[i, ", i, ", ", i, "]) * ",
                 "sqrt(S3[i, ", l, ", ", l, "]) * ",
                 "control", i, j, "[i] * rho", r, "\n")
      }
    }
    #
    txt <- paste0(txt, "    #\n")
    #
    r <- 0
    #
    for (j in seq_len(n.out - 1)) {
      for (i in (j + 1):n.out){
        l <- n.out + j
        #
        r <- r + 1
        #
        txt <-
          paste0(
            txt,
            "    S3[i, ", i, ", ", l, "] <- 0.5 * ",
            "sqrt(S3[i, ", i, ", ", i, "]) * ",
            "sqrt(S3[i, ", l, ", ", l, "]) * control",
            j, i, "[i] * rho", r, "\n")
      }
    }
    #
    txt <- paste0(txt, "    #\n")
    #
    for (i in seq_len(2 * n.out - 1)) {
      for (j in (i + 1):(2 * n.out)) {
        txt <-
          paste0(txt,
                 "    S3[i, ", j, ", ", i, "] <- S3[i, ", i, ", ", j, "]\n")
      }
    }
    #
    txt <-
      paste0(txt,
             "    #\n",
             "    # Estimates\n",
             "    #\n")
    #
    txt <-
      paste0(txt,
             "    y[(", n.out, " * k2 + ", 2 * n.out, " * i - ",
             2 * n.out - 1, "):(", n.out, " * k2 + ",
             2 * n.out, " * i)] ~\n      dmnorm.vcov(mean[(",
             n.out, " * k2 + ", 2 * n.out, " * i - ",
             2 * n.out - 1,
             "):(",
             n.out, " * k2 + ", 2 * n.out,
             " * i)], S3[i, , ])\n")
    #
    txt <- paste0(txt, "  }\n")
  }
  #
  txt <- paste0(txt, "  \n")
  #
  txt
}

code_means <- function(n.out, multiarm) {
  
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
  for (i in seq_len(n.out)) {
    txt <-
      paste0(txt, "    mean[", n.out, " * i",
             if (i != n.out)
               paste0(" - ", n.out - i)
             else
               strrep(" ", nchar(n.out - 1) + 3),
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
    idx <- rep(seq_len(n.out), 2)
    #
    for (i in seq_len(2 * n.out)) {
      txt <-
        paste0(txt, "    mean[", n.out, " * k2 + ", 2 * n.out, " * i",
               if (i != 2 * n.out)
                 paste0(" - ", 2 * n.out - i)
               else
                 strrep(" ", nchar(2 * n.out - 1) + 3),
               "] <- d",
               idx[i], "[treat", 2 + (i > n.out), "[k2 + i]] - d",
               idx[i], "[treat1[k2 + i]]\n")
      #
      if (i != 2 * n.out & idx[i] == n.out)
        txt <- paste0(txt, "    #\n")
    }
    #
    txt <- paste0(txt, "  }\n")
  }
  txt
}

code_priors <- function(n.out, method) {
  if (method == "standard")
    txt <- code_priors_standard(n.out)
  else
    txt <- code_priors_dumouchel(n.out)
  #
  txt <- paste0(txt, "  #\n")
  #
  txt <- paste0(txt, code_priors_psi(n.out))
  #
  txt <- paste0(txt, "  #\n")
  #
  txt <- paste0(txt, code_priors_rho(n.out))
  #
  txt
}

code_priors_standard <- function(n.out) {
  
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
  for (i in seq_len(n.out))
    txt <- paste0(txt, "    d", i, "[i] ~ dnorm(0, 1e-03)\n")
  #
  txt <- paste0(txt, "  }\n")
  #
  txt <- paste0(txt, "  #\n")
  for (i in seq_len(n.out))
    txt <- paste0(txt, "  d", i, "[ref] <- 0\n")
  txt <- paste0(txt, "  #\n")
  #
  txt <-
    paste0(txt,
           "  for (i in (ref + 1):n) {\n")
  #
  for (i in seq_len(n.out))
    txt <- paste0(txt, "    d", i, "[i] ~ dnorm(0, 1e-03)\n")
  #
  txt <- paste0(txt, "  }\n")
  #
  txt
}

code_priors_dumouchel <- function(n.out) {
  
  txt <-
    paste0(
      "  #\n",
      "  #\n",
      "  # (4) Priors (DuMouchel model)\n",
      "  #\n",
      "  #\n\n")
  #
  txt <-
    paste0(txt,
           "  for (i in 1:(ref - 1)) {\n",
           "    for (m in 1:", n.out, ") {\n",
           "      meand[m, i] <- alpha[i] + gamma[m]\n",
           "      d[m, i] ~ dnorm(meand[m, i], prec.exp)\n",
           "    }\n",
           "    #\n")
  #
  for (i in seq_len(n.out))
    txt <- paste0(txt, "    d", i, "[i] <- d[", i, ", i]\n")
  #
  txt <- paste0(txt, "  }\n")
  #
  txt <- paste0(txt, "  #\n")
  for (i in seq_len(n.out))
    txt <- paste0(txt, "  d", i, "[ref] <- 0\n")
  txt <- paste0(txt, "  #\n")
  #
  txt <-
    paste0(txt,
           "  for (i in (ref + 1):n) {\n",
           "    for (m in 1:", n.out, ") {\n",
           "      meand[m, i] <- alpha[i] + gamma[m]\n",
           "      d[m, i] ~ dnorm(meand[m, i], prec.exp)\n",
           "    }\n",
           "    #\n")
  #
  for (i in seq_len(n.out))
    txt <- paste0(txt, "    d", i, "[i] <- d[", i, ", i]\n")
  #
  txt <- paste0(txt, "  }\n")
  #
  txt <- paste0(txt, "  #\n")
  #
  txt <-
    paste0(txt,
           "  for (m in 1:", n.out, ") {\n",
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
           "  sigma ~ dnorm(0, 1e-02)T(0, )\n")
  #
  txt
}

code_priors_psi <- function(n.out) {
  txt <- ""
  #
  for (i in seq_len(n.out))
    txt <- paste0(txt, "  psi", i, ".sq <- psi", i, " * psi", i, "\n")
  #
  txt <- paste0(txt, "  #\n")
  #
  for (i in seq_len(n.out))
    txt <- paste0(txt, "  psi", i, "  ~ dnorm(0, prec.psi", i, ")T(0, )\n")
  #
  txt
}

code_priors_rho <- function(n.out) {
  txt <- ""
  #
  for (i in seq_len(choose(n.out, 2)))
    txt <-
      paste0(txt, "  rho", i, " ~ dunif(lower.rho", i,
             ", upper.rho", i, ")\n")
  #
  txt
}
