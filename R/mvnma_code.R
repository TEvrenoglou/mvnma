mvnma_code <- function(n.out, method, multiarm, n.dom) {
  
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
  # Files with variances, covariances and estimates
  #
  txt <- paste0(txt, code_covar_ests(n.out, multiarm))
  #
  txt <- paste0(txt, code_means(n.out, multiarm))
  #
  txt <- paste0(txt, "\n")
  #
  txt <- paste0(txt, code_priors(n.out, method,n.dom))
  #
  txt <- paste0(txt, "}\n")
  #
  txt
}

code_covar_ests <- function(n.out, multiarm) {
  
  txt <-
    paste0(
      "  #\n",
      "  #\n",
      "  # (1) Variances, covariances and estimates\n",
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
             "    S2[i, ", i, ", ", i, "] <- varmat[i, ", i, "] + ",
             "psi.sq[", i, "]\n")
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
               "contmat", "[i, ", i, "] * ",
               "contmat", "[i, ", j, "] * ",
               "rho[", r, "]\n")
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
        paste0(txt, "    S3[i, ", i, ", ", i, "] <- ",
               "varmat[k2 + 2 * i", if (i <= n.out) " - 1", ", ",
               idx[i], "] + psi.sq[", idx[i], "]\n")
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
                 "contmat", "[k2 + i, ", i, "] * ",
                 "contmat", "[k2 + i, ", j, "] * ",
                 "rho[", r, "]\n")
      }
    }
    #
    txt <- paste0(txt, "    #\n")
    #
    r <- 0
    #
    for (i in seq_len(n.out - 1)) {
      for (j in (i + 1):n.out) {
        r <- r + 1
        #
        txt <-
          paste0(txt,
                 "    S3[i, ", n.out + i, ", ", n.out + j, "] <- ",
                 "sqrt(S3[i, ", n.out + i, ", ", n.out + i, "]) * ",
                 "sqrt(S3[i, ", n.out + j, ", ", n.out + j, "]) * ",
                 "contmat", "[k2 + i, ", i, "] * ",
                 "contmat", "[k2 + i, ", j, "] * ",
                 "rho[", r, "]\n")
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
        r <- r + 1
        #
        txt <-
          paste0(txt,
                 "    S3[i, ", i, ", ", j + n.out, "] <- 0.5 * ",
                 "sqrt(S3[i, ", i, ", ", i, "]) * ",
                 "sqrt(S3[i, ", j + n.out, ", ", j + n.out, "]) * ",
                 "contmat", "[k2 + i, ", i, "] * ",
                 "contmat", "[k2 + i, ", j, "] * ",
                 "rho[", r, "]\n")
      }
    }
    #
    txt <- paste0(txt, "    #\n")
    #
    r <- 0
    #
    for (j in seq_len(n.out - 1)) {
      for (i in (j + 1):n.out){
        r <- r + 1
        #
        txt <-
          paste0(
            txt,
            "    S3[i, ", i, ", ", n.out + j, "] <- 0.5 * ",
            "sqrt(S3[i, ", i, ", ", i, "]) * ",
            "sqrt(S3[i, ", n.out + j, ", ", n.out + j, "]) * ",
            "contmat", "[k2 + i, ", j, "] * ",
            "contmat", "[k2 + i, ", i, "] * ",
            "rho[", r, "]\n")
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
      "  # (2) Parameterization of the means\n",
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
             "] <- d", i, "[trtmat[i, 2]] - d", i, "[trtmat[i, 1]]\n")
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
               idx[i], "[trtmat[k2 + i, ", 2 + (i > n.out), "]] - ",
               "d", idx[i], "[trtmat[k2 + i, 1]]\n")
      #
      if (i != 2 * n.out & idx[i] == n.out)
        txt <- paste0(txt, "    #\n")
    }
    #
    txt <- paste0(txt, "  }\n")
  }
  txt
}

code_priors <- function(n.out, method,n.dom) {
  if (method == "standard")
    txt <- code_priors_standard(n.out)
  else
    txt <- code_priors_dumouchel(n.out,n.dom)
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
      "  # (3) Priors (standard model)\n",
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

code_priors_dumouchel <- function(n.out,n.dom) {
  
  txt <-
    paste0(
      "  #\n",
      "  #\n",
      "  # (4) Priors (DuMouchel model)\n",
      "  #\n",
      "  #\n\n")
  #
  if (is.null(n.dom)) {
    txt <-
      paste0(txt,
             "  for (i in 1:(ref - 1)) {\n",
             "    for (m in 1:", n.out, ") {\n",
             "      meand[m, i] <- alpha[i] + gamma[m]\n",
             "      d[m, i] ~ dnorm(meand[m, i], prec.exp)\n",
             "    }\n",
             "    #\n")
  }
  else {
    txt <-
      paste0(txt,
             "  for (i in 1:(ref - 1)) {\n",
             "    for (m in 1:", n.dom, ") {\n",
             "      meand[m, i] <- alpha1[i] + gamma1[m]\n",
             "      d[m, i] ~ dnorm(meand[m, i], prec.exp1)\n",
             "    }\n",
             "    for (l in ", n.dom+1, " : ", n.out, ") {\n",
             "      meand[l, i] <- alpha2[i] + gamma2[l]\n",
             "      d[l, i] ~ dnorm(meand[l, i], prec.exp2)\n",
             "    }\n",
             "    #\n")
  }
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
  if (is.null(n.dom)) {
    txt <-
      paste0(txt,
             "  for (i in (ref + 1):n) {\n",
             "    for (m in 1:", n.out, ") {\n",
             "      meand[m, i] <- alpha[i] + gamma[m]\n",
             "      d[m, i] ~ dnorm(meand[m, i], prec.exp)\n",
             "    }\n",
             "    #\n")
  }
  else {
    txt <-
      paste0(txt,
             "  for (i in (ref + 1):n) {\n",
             "    for (m in 1:", n.dom, ") {\n",
             "      meand[m, i] <- alpha1[i] + gamma1[m]\n",
             "      d[m, i] ~ dnorm(meand[m, i], prec.exp1)\n",
             "    }\n",
             "    for (l in ", n.dom+1, " : ", n.out, ") {\n",
             "      meand[l, i] <- alpha2[i] + gamma2[l]\n",
             "      d[l, i] ~ dnorm(meand[l, i], prec.exp2)\n",
             "    }\n",
             "    #\n")
  }
  #
  for (i in seq_len(n.out))
    txt <- paste0(txt, "    d", i, "[i] <- d[", i, ", i]\n")
  #
  txt <- paste0(txt, "  }\n")
  #
  txt <- paste0(txt, "  #\n")
  #
  if (is.null(n.dom)) {
    txt <-
      paste0(txt,
             "  for (m in 1:", n.out, ") {\n",
             "  gamma[m] ~ dnorm(0, 1e-03)\n",
             "  }\n",
             "  #\n")
  }
  else{
    txt <-
      paste0(txt,
             "  for (m in 1:", n.dom, ") {\n",
             "  gamma1[m] ~ dnorm(0, 1e-03)\n",
             "  }\n",
             "  for (l in ", n.dom+1, " : ", n.out, ") {\n",
             "  gamma2[l] ~ dnorm(0, 1e-03)\n",
             "  }\n",
             "  #\n")
  }
  #
  if (is.null(n.dom)) {
    txt <-
      paste0(txt,
             "  for (i in 1:(ref - 1)) {\n",
             "    alpha[i] ~ dnorm(0, 1e-03)\n",
             "  }\n",
             "  #\n")
  }
  else {
    txt <-
      paste0(txt,
             "  for (i in 1:(ref - 1)) {\n",
             "    alpha1[i] ~ dnorm(0, 1e-03)\n",
             "    alpha2[i] ~ dnorm(0, 1e-03)\n",
             "  }\n",
             "  #\n")
  }
  #
  if (is.null(n.dom)) {
    txt <-
      paste0(txt,
             "  for (i in (ref + 1):n) {\n",
             "    alpha[i] ~ dnorm(0, 1e-03)\n",
             "  }\n",
             "  #\n")
  }
  else {
    txt <-
      paste0(txt,
             "  for (i in (ref + 1):n) {\n",
             "    alpha1[i] ~ dnorm(0, 1e-03)\n",
             "    alpha2[i] ~ dnorm(0, 1e-03)\n",
             "  }\n",
             "  #\n")
    
  }
  #
  if (is.null(n.dom)) {
    txt <-
      paste0(txt,
             "  prec.exp <- 1 / sigma.sq\n",
             "  sigma.sq <- sigma * sigma\n",
             "  sigma ~ dnorm(0, 1e-02)T(0, )\n")
  }
  else {
    txt <-
      paste0(txt,
             "  prec.exp1 <- 1 / sigma.sq1\n",
             "  sigma.sq1 <- sigma1 * sigma1\n",
             "  sigma1 ~ dnorm(0, 1e-02)T(0, )\n",
             "  prec.exp2 <- 1 / sigma.sq2\n",
             "  sigma.sq2 <- sigma2 * sigma2\n",
             "  sigma2 ~ dnorm(0, 1e-02)T(0, )\n"
      )
  }
  #
  txt
}

code_priors_psi <- function(n.out) {
  txt <- ""
  #
  for (i in seq_len(n.out))
    txt <- paste0(txt, "  psi.sq[", i, "] <- psi[", i, "] * psi[", i, "]\n")
  #
  txt <- paste0(txt, "  #\n")
  #
  for (i in seq_len(n.out))
    txt <- paste0(txt, "  psi[", i, "]  ~ dnorm(0, prec.psi[", i, "])T(0, )\n")
  #
  txt
}

code_priors_rho <- function(n.out) {
  txt <- ""
  #
  for (i in seq_len(choose(n.out, 2)))
    txt <-
      paste0(txt, "  rho[", i, "] ~ dunif(lower.rho[", i,
             "], upper.rho[", i, "])\n")
  #
  txt
}
