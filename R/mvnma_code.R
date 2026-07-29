mvnma_code <- function(n.out, arms = 2, method, n.dom) {
  
  chknumeric(n.out, min = 2, length = 1)
  method <- setchar(method, c("standard", "DM"))
  chknumeric(arms, min = 2, NA.ok = FALSE)
  arms <- as.integer(arms)
  #
  txt <-
    paste0("model {\n",
           "  # n.studies = number of studies by number of arms\n",
           "  # n = number of treatments\n\n")
  #
  # Files with variances, covariances and estimates
  #
  txt <- paste0(txt, code_covar_ests(n.out, arms))
  #
  txt <- paste0(txt, code_means(n.out, arms))
  #
  txt <- paste0(txt, "\n")
  #
  txt <- paste0(txt, code_priors(n.out, method,n.dom))
  #
  txt <- paste0(txt, "}\n")
  #
  txt
}

code_covar_ests <- function(n.out, arms = 2) {
  code_covar_ests_arm(n.out, arms)
}

code_covar_ests_arm <- function(n.out, arms) {
  
  txt <-
    paste0(
      "  #\n",
      "  #\n",
      "  # (1) Variances, covariances and estimates\n",
      "  #\n",
      "  #\n  \n")
  #
  if (2 %in% arms) {
    txt <-
      paste0(txt,
             "  #\n",
             "  # 2-arm studies\n",
             "  #\n")
    #
    txt <- paste0(txt, "  for (i in 1:n.studies[", match(2, arms), "]) {\n")
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
  }
  #
  for (a in arms[arms > 2])
    txt <- paste0(txt, code_covar_ests_multiarm(n.out, a, arms))
  #
  txt <- paste0(txt, "  \n")
  #
  txt
}

code_covar_ests_multiarm <- function(n.out, a, arms) {
  n.contr <- a - 1
  n.par <- n.out * n.contr
  arm.index <- match(a, arms)
  var_offset <- code_var_offset(a, arms)
  y_offset <- code_y_offset(a, n.out, arms)
  txt <-
    paste0(
      "  #\n",
      "  # ", a, "-arm studies\n",
      "  #\n",
      "  for (i in 1:n.studies[", arm.index, "]) {\n",
      "    #\n",
      "    # Variances\n",
      "    #\n")
  #
  pos <- expand.grid(outcome = seq_len(n.out), contrast = seq_len(n.contr))
  pos <- pos[order(pos$contrast, pos$outcome), ]
  #
  for (p in seq_len(n.par)) {
    row <- code_multiarm_row(var_offset, n.contr, "i", pos$contrast[p])
    txt <-
      paste0(txt,
             "    S", a, "[i, ", p, ", ", p, "] <- ",
             "varmat[", row, ", ", pos$outcome[p], "] + ",
             "psi.sq[", pos$outcome[p], "]\n")
    #
    if (pos$outcome[p] == n.out)
      txt <- paste0(txt, "    #\n")
  }
  #
  txt <-
    paste0(txt,
           "    # Covariances\n",
           "    #\n")
  #
  for (p1 in seq_len(n.par - 1)) {
    for (p2 in (p1 + 1):n.par) {
      txt <- paste0(txt, code_multiarm_covariance(a, n.contr, var_offset,
                                                  pos, p1, p2))
    }
  }
  #
  txt <- paste0(txt, "    #\n")
  #
  for (p1 in seq_len(n.par - 1)) {
    for (p2 in (p1 + 1):n.par) {
      #
      txt <-
        paste0(txt,
               "    S", a, "[i, ", p2, ", ", p1, "] <- ",
               "S", a, "[i, ", p1, ", ", p2, "]\n")
    }
  }
  #
  txt <-
    paste0(txt,
           "    #\n",
           "    # Estimates\n",
           "    #\n",
           "    y[(", y_offset, " + ", n.par, " * i - ", n.par - 1,
           "):(", y_offset, " + ", n.par, " * i)] ~\n",
           "      dmnorm.vcov(mean[(", y_offset, " + ", n.par, " * i - ",
           n.par - 1, "):(", y_offset, " + ", n.par, " * i)], ",
           "S", a, "[i, , ])\n",
           "  }\n")
  #
  txt
}

code_multiarm_covariance <- function(a, n.contr, var_offset, pos, p1, p2) {
  o1 <- pos$outcome[p1]
  o2 <- pos$outcome[p2]
  c1 <- pos$contrast[p1]
  c2 <- pos$contrast[p2]
  row1 <- code_multiarm_row(var_offset, n.contr, "i", c1)
  row2 <- code_multiarm_row(var_offset, n.contr, "i", c2)
  txt <-
    paste0("    S", a, "[i, ", p1, ", ", p2, "] <- ",
           if (c1 == c2) "" else "0.5 * ",
           "sqrt(S", a, "[i, ", p1, ", ", p1, "]) * ",
           "sqrt(S", a, "[i, ", p2, ", ", p2, "])")
  #
  if (o1 != o2) {
    txt <-
      paste0(txt,
             " * contmat[", row1, ", ", o1, "] * ",
             "contmat[", row2, ", ", o2, "] * ",
             "rho[", code_rho_index(o1, o2, max(pos$outcome)), "]")
  }
  #
  paste0(txt, "\n")
}

code_rho_index <- function(i, j, n.out) {
  if (i > j) {
    tmp <- i
    i <- j
    j <- tmp
  }
  #
  sum(n.out - seq_len(i - 1)) + j - i
}

code_multiarm_row <- function(offset, n.contr, study, contrast) {
  code_index(offset, paste0(n.contr, " * ", study, " - ",
                            n.contr - contrast))
}

code_index <- function(offset, term) {
  if (offset == "0")
    term
  else
    paste0(offset, " + ", term)
}

code_var_offset <- function(a, arms) {
  prev <- seq_len(match(a, arms) - 1)
  if (!length(prev))
    return("0")
  #
  paste0((arms[prev] - 1), " * n.studies[", prev, "]", collapse = " + ")
}

code_y_offset <- function(a, n.out, arms) {
  prev <- seq_len(match(a, arms) - 1)
  if (!length(prev))
    return("0")
  #
  paste0(n.out * (arms[prev] - 1), " * n.studies[", prev, "]",
         collapse = " + ")
}

code_means <- function(n.out, arms = 2) {
  code_means_arm(n.out, arms)
}

code_means_arm <- function(n.out, arms) {

  txt <-
    paste0(
      "  #\n",
      "  #\n",
      "  # (2) Parameterization of the means\n",
      "  #\n",
      "  #\n\n")
  #
  if (2 %in% arms) {
    txt <-
      paste0(txt,
             "  #\n",
             "  # 2-arm studies\n",
             "  #\n")
    #
    txt <- paste0(txt, "  for (i in 1:n.studies[", match(2, arms), "]) {\n")
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
  }
  #
  for (a in arms[arms > 2])
    txt <- paste0(txt, code_means_multiarm(n.out, a, arms))
  #
  txt
}

code_means_multiarm <- function(n.out, a, arms) {
  n.contr <- a - 1
  n.par <- n.out * n.contr
  arm.index <- match(a, arms)
  y_offset <- code_y_offset(a, n.out, arms)
  study_offset <- code_study_offset(a, arms)
  txt <-
    paste0(
      "  #\n",
      "  # ", a, "-arm studies\n",
      "  #\n",
      "  for (i in 1:n.studies[", arm.index, "]) {\n")
  #
  p <- 0
  for (c in seq_len(n.contr)) {
    for (o in seq_len(n.out)) {
      p <- p + 1
      txt <-
        paste0(txt, "    mean[", y_offset, " + ", n.par, " * i",
               if (p != n.par)
                 paste0(" - ", n.par - p)
               else
                 strrep(" ", nchar(n.par - 1) + 3),
               "] <- d", o, "[trtmat[", code_index(study_offset, "i"),
               ", ", c + 1, "]] - ",
               "d", o, "[trtmat[", code_index(study_offset, "i"),
               ", 1]]\n")
    }
    #
    if (c != n.contr)
      txt <- paste0(txt, "    #\n")
  }
  #
  paste0(txt, "  }\n")
}

code_study_offset <- function(a, arms) {
  prev <- seq_len(match(a, arms) - 1)
  if (!length(prev))
    return("0")
  #
  paste0("n.studies[", prev, "]", collapse = " + ")
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
