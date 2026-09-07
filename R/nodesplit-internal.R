# ---------------------------------------------------------------------------
# splittable.comparisons()
#
# Identify the pairwise treatment comparisons for which BOTH a direct and an
# indirect estimate can be obtained, and which can therefore be node-split.
#
# A comparison (t1, t2) qualifies when:
#   (a) direct evidence exists   -> at least one study reports a non-missing
#       estimate (TE / seTE) for t1 vs t2 in this outcome;
#   (b) indirect evidence exists -> after removing every STUDY that contains
#       both t1 and t2, t1 and t2 are still connected in the network formed
#       by the remaining studies.
#
# ---------------------------------------------------------------------------

# --- main function ---------------------------------------------------------

# x       a pairwise data set (data.frame with studlab / treat1 / treat2), or
#         a list of such data sets - one per outcome
# combine how to combine results across outcomes when x is a list:
#         "any"  - comparison is splittable in at least one outcome (default)
#         "all"  - comparison is splittable in every outcome
#         "list" - no combining; return one data.frame per outcome

splittable.comparisons <- function(x, combine = c("any", "all", "list")) {
  
  combine <- match.arg(combine)
  
  # single data.frame -> nothing to combine
  if (is.data.frame(x)) {
    res <- .splittable.one(x)
    # sort first, then reset row names - resetting before sorting would
    # simply reattach the pre-sort numbering
    res <- res[order(res$treat1, res$treat2), , drop = FALSE]
    row.names(res) <- NULL
    return(res)
  }
  
  if (!is.list(x))
    stop("x should be a data.frame or a list of data.frames")
  
  if (length(x) == 0)
    stop("x is empty")
  
  if (!all(vapply(x, is.data.frame, logical(1))))
    stop("all elements of x should be data.frames")
  
  res.list <- lapply(x, function(d) {
    z <- .splittable.one(d)
    z <- z[order(z$treat1, z$treat2), , drop = FALSE]
    row.names(z) <- NULL
    z
  })
  
  if (combine == "list") {
    if (is.null(names(res.list)))
      names(res.list) <- paste0("outcome", seq_along(res.list))
    return(res.list)
  }
  
  # combine across outcomes on the "treat1 vs treat2" key
  keys <- lapply(res.list,
                 function(z) paste(z$treat1, z$treat2, sep = " \001 "))
  
  sel <- if (combine == "any")
    Reduce(union, keys)
  else
    Reduce(intersect, keys)
  
  if (length(sel) == 0) {
    return(data.frame(treat1 = character(0),
                      treat2 = character(0),
                      stringsAsFactors = FALSE))
  }
  
  parts <- do.call(rbind, strsplit(sel, " \001 ", fixed = TRUE))
  
  res <- data.frame(treat1 = parts[, 1],
                    treat2 = parts[, 2],
                    stringsAsFactors = FALSE)
  
  res <- res[order(res$treat1, res$treat2), , drop = FALSE]
  row.names(res) <- NULL
  
  res
}

# --- helpers ---------------------------------------------------------------

# study-level long format, with treatments sorted within each row so that
# "A vs B" and "B vs A" are recognised as the same comparison
.contrasts <- function(d) {
  
  if (!all(c("studlab", "treat1", "treat2") %in% names(d)))
    stop("pairwise data must contain columns 'studlab', 'treat1' and 'treat2'")
  
  studlab <- as.character(d$studlab)
  t1 <- as.character(d$treat1)
  t2 <- as.character(d$treat2)
  
  keep <- !is.na(studlab) & !is.na(t1) & !is.na(t2) & t1 != t2
  
  # A row only counts as evidence if it actually carries an estimate for this
  # outcome. In the multivariate setting the per-outcome pairwise objects keep
  # rows for contrasts a study did not report, with TE / seTE set to NA;
  # treating those as direct evidence would overstate the network.
  if ("TE" %in% names(d))
    keep <- keep & !is.na(d$TE)
  
  if ("seTE" %in% names(d))
    keep <- keep & !is.na(d$seTE)
  
  data.frame(studlab = studlab[keep],
             treat1  = pmin(t1, t2)[keep],
             treat2  = pmax(t1, t2)[keep],
             stringsAsFactors = FALSE)
}


# breadth-first search: are `from` and `to` connected via `edges`?
.connected <- function(edges, from, to) {
  
  if (nrow(edges) == 0)
    return(FALSE)
  
  nodes <- unique(c(edges$treat1, edges$treat2))
  
  if (!(from %in% nodes) || !(to %in% nodes))
    return(FALSE)
  
  # undirected adjacency list
  adj <- split(c(edges$treat2, edges$treat1),
               c(edges$treat1, edges$treat2))
  
  seen  <- from
  queue <- from
  
  while (length(queue) > 0) {
    
    cur   <- queue[1]
    queue <- queue[-1]
    
    nb <- setdiff(adj[[cur]], seen)
    
    if (to %in% nb)
      return(TRUE)
    
    seen  <- c(seen, nb)
    queue <- c(queue, nb)
  }
  
  FALSE
}


# splittable comparisons for a single pairwise data set
.splittable.one <- function(d) {
  
  x <- .contrasts(d)
  
  comps <- unique(x[c("treat1", "treat2")])
  
  if (nrow(comps) == 0)
    return(comps)
  
  ok <- logical(nrow(comps))
  
  for (k in seq_len(nrow(comps))) {
    
    t1 <- comps$treat1[k]
    t2 <- comps$treat2[k]
    
    # # every study providing direct evidence on t1 vs t2
    # drop.studies <- unique(x$studlab[x$treat1 == t1 & x$treat2 == t2])
    # 
    # # network formed by the remaining studies
    # rest <- x[!(x$studlab %in% drop.studies), , drop = FALSE]
    
    # remove only the t1 vs t2 contrast; remaining arms of a multi-arm study
    # still contribute to the indirect estimate (Dias et al., 2010)
    rest <- x[!(x$treat1 == t1 & x$treat2 == t2), , drop = FALSE]
    
    ok[k] <- .connected(unique(rest[c("treat1", "treat2")]), t1, t2)
  }
  
  comps[ok, , drop = FALSE]
}

# ---------------------------------------------------------------------------
# Node-splitting for a single treatment comparison
# ---------------------------------------------------------------------------
pair.nodesplit <- function(x, treat1, treat2, 
                           method.direct = "pairwise", 
                           tol.direct = 5e-04, ...){
  
  #
  method.model <- attr(x, "method.model")
  n.domain <- attr(x,"n.domain")
  outcomes <- attr(x, "outcomes")
  psi <- extract_het(x)
  #
  level <- attr(x, "level")
  #
  method.direct <- setchar(method.direct, c("pairwise", "multivariate"))
  #
  # z-critical value for the chosen confidence level - computed once and
  # reused below instead of being recomputed on every use
  z.crit <- qnorm(1 - (1 - level) / 2)
  
  # calculate the direct estimate, either by re-fitting the multivariate
  # model to the direct evidence or by an outcome-specific pairwise
  # meta-analysis with heterogeneity fixed at the multivariate estimate
  dir <- if (method.direct == "pairwise")
    direct.metagen(x, treat1 = treat1, treat2 = treat2, ...)
  else
    direct.mvnma(x, treat1 = treat1, treat2 = treat2, ...)
  
  # which outcomes were actually fitted, and where each sits in `dir`
  keep.out <- attr(dir, "keep")
  #
  if (is.null(keep.out))
    keep.out <- rep(TRUE, length(outcomes))
  #
  # direct.mvnma() drops outcomes without direct evidence from the model fit,
  # so results must be mapped back by position; direct.metagen() returns one
  # element per outcome and needs no mapping
  pos <- if (method.direct == "pairwise")
    seq_along(keep.out)
  else
    cumsum(keep.out)
  #
  # per-outcome direct-comparison data and study counts, pulled off `dir`
  # *before* it gets subsetted below (subsetting drops non-standard
  # attributes, same reason `outcomes`/`level` are captured from `x` early)
  direct.data <- attr(dir, "direct.data")
  n.studies   <- attr(dir, "n.studies")
  
  # keep only outcome treatment effect estimates for mvnma
  x <- x[names(x) != "cor"]
  #
  if (method.model == "DM") {
    if (is.null(n.domain)) {
      x <- x[names(x) != "sigma"]
    }
    else{
      x <- x[!(names(x) %in% c("sigma1", "sigma2"))]
    }
  }
  # keep only outcome treatment effect estimates for the direct comparison 
  dir <- dir[names(dir) != "cor"]
  #
  if (method.model == "DM") {
    if (is.null(n.domain)) {
      dir <- dir[names(dir) != "sigma"]
    }
    else{
      dir <- dir[!(names(dir) %in% c("sigma1", "sigma2"))]
    }
  }
  
  n.out <- length(x)
  
  direct <- indirect <- overall <- difference <- res.final <- vector("list", n.out)
  
  for(i in 1:n.out){
    
    if (!is.null(n.studies) && n.studies[i] == 0) {
      
      # no direct evidence for this comparison in this outcome
      direct[[i]] <- data.frame("mean"=NA,"sd"=NA,"lower"=NA,"upper"=NA)
      
    } else if (!is.null(direct.data) && !is.null(n.studies) && n.studies[i] == 1) {
      
      # direct.data keeps rows with missing estimates; under
      # method.direct = "pairwise" a single study can be one non-missing row
      # among several, so select it explicitly
      sel <- !is.na(direct.data[[i]]$TE) & !is.na(direct.data[[i]]$seTE)
      #
      #
      if (!any(sel)) {
        # the single row carries no estimate for this outcome
        direct[[i]] <- data.frame("mean"=NA,"sd"=NA,"lower"=NA,"upper"=NA)
      }
      else {
        study.TE   <- direct.data[[i]]$TE[sel]
        study.seTE <- sqrt(direct.data[[i]]$seTE[sel]^2 + psi[i]^2)
        #
        direct[[i]] <- data.frame("mean" = study.TE,
                                  "sd" = study.seTE,
                                  "lower" = study.TE - z.crit*study.seTE,
                                  "upper" = study.TE + z.crit*study.seTE)
      }
    } else if (!keep.out[i] || length(dir) == 0 ||
               all(is.na(dir[[pos[i]]]$basic_estimates))) {
      
      direct[[i]] <- data.frame("mean"=NA,"sd"=NA,"lower"=NA,"upper"=NA)  
      
    } else {
      
      direct[[i]] <- dir[[pos[i]]]$basic_estimates[complete.cases(dir[[pos[i]]]$basic_estimates$mean),]  
      direct[[i]] <- direct[[i]] %>% 
        select(mean,sd,lower,upper)
      
    }
    ## get mvnma estimate for full analysis
    overall.TE <- get_cell(x[[i]]$TE.random,treat1 = treat1,treat2 = treat2)
    overall.seTE <- get_cell(x[[i]]$seTE.random,treat1 = treat1,treat2 = treat2)
    
    overall[[i]] <- cbind.data.frame("mean"= overall.TE,
                                     "sd" = overall.seTE,
                                     "lower" = overall.TE-z.crit*overall.seTE,
                                     "upper" = overall.TE+z.crit*overall.seTE)
    
    # direct evidence proportion: ratio of the network variance to the
    # direct variance (netmeta:::netmeasures)
    prop.dir <- overall[[i]]$sd^2 / direct[[i]]$sd^2
    
    indirect[[i]] <- back_calc_indirect(mean.overall = overall[[i]]$mean,
                                        var.overall = overall[[i]]$sd^2,
                                        mean.dir = direct[[i]]$mean,
                                        prop = prop.dir,
                                        tol.direct = tol.direct)
    
    indirect[[i]]$lower <- indirect[[i]]$mean - z.crit*indirect[[i]]$sd
    indirect[[i]]$upper <- indirect[[i]]$mean + z.crit*indirect[[i]]$sd
    #
    difference[[i]] <- z_test_diff(direct.TE = direct[[i]]$mean,
                                   direct.seTE = direct[[i]]$sd,
                                   indirect.TE = indirect[[i]]$mean,
                                   indirect.seTE = indirect[[i]]$sd,
                                   level = level)
    #
    res.final[[i]] <- data.frame("k" = if (is.null(n.studies)) NA else n.studies[i],
                                 "prop" = prop.dir,
                                 "mvnma.TE" = overall[[i]]$mean,
                                 "mvnma.seTE" = overall[[i]]$sd,
                                 "mvnma.lb" = overall[[i]]$lower,
                                 "mvnma.ub"= overall[[i]]$upper,
                                 "direct.TE" = direct[[i]]$mean,
                                 "direct.seTE" = direct[[i]]$sd,
                                 "direct.lb" = direct[[i]]$lower,
                                 "direct.ub" = direct[[i]]$upper,
                                 "indirect.TE" = indirect[[i]]$mean,
                                 "indirect.seTE" = indirect[[i]]$sd,
                                 "indirect.lb" = indirect[[i]]$lower,
                                 "indirect.ub" = indirect[[i]]$upper,
                                 "diff" = difference[[i]]$diff,
                                 "se.diff" = difference[[i]]$se.diff,
                                 "z" = difference[[i]]$z,
                                 "diff.lb" = difference[[i]]$lower.diff,
                                 "diff.ub" = difference[[i]]$upper.diff,
                                 "p.val"= difference[[i]]$p.val,
                                 "sign" = difference[[i]]$significant)
    #
    row.names(res.final[[i]]) <- paste0(treat1,":",treat2)
  }
  #
  names(res.final) <- outcomes
  return(res.final)
}

# ---------------------------------------------------------------------------
# Direct estimates for a single comparison: multivariate refit
# (direct.mvnma) or outcome-specific pairwise meta-analysis (direct.metagen)
# ---------------------------------------------------------------------------

# direct estimate from a re-fit of the multivariate model
direct.mvnma <- function(x, treat1, treat2, ...){
  
  data.dir <- direct.pair(x, t1 = treat1, t2 = treat2)
  
  psi.preset <- extract_het(x)
  
  outcomes <- attr(x, "outcomes")
  n.out    <- length(data.dir)
  
  # number of studies contributing direct evidence, per outcome
  n.studies <- vapply(data.dir, nrow, integer(1))
  
  # Outcomes for which this comparison is directly informed. mvnma() cannot be
  # given an empty pairwise object, so outcomes without any direct evidence are
  # dropped from the model fit; their direct estimate is set to NA in
  # pair.nodesplit().
  keep <- n.studies > 0
  
  names(data.dir) <- NULL
  
  # heterogeneity estimates are outcome-specific and must be subset alongside
  if (length(psi.preset) == n.out)
    psi.preset <- psi.preset[keep]
  
  scale.psi <- attr(x, "scale.psi")
  #
  if (length(scale.psi) == n.out)
    scale.psi <- scale.psi[keep]
  #
  n.rho <- sum(keep) * (sum(keep) - 1) / 2
  
  lower.rho <- attr(x, "lower.rho")
  upper.rho <- attr(x, "upper.rho")
  #
  # correlations are indexed by pairs of outcomes, not by outcome, so
  # outcome-specific values cannot be carried over to a reduced model
  if (length(lower.rho) != 1 && length(lower.rho) != n.rho)
    lower.rho <- NULL
  #
  if (length(upper.rho) != 1 && length(upper.rho) != n.rho)
    upper.rho <- NULL
  # carry over the original model-fit settings
  fit.args <- list(
    scale.psi = scale.psi,
    lower.rho = lower.rho,
    upper.rho = upper.rho,
    n.iter    = attr(x, "n.iter"),
    n.burnin  = attr(x, "n.burnin"),
    n.chains  = attr(x, "n.chains"),
    method    = "standard", # run always mvnma with standard model
    outclab = outcomes[keep],
    reference.group = treat2,
    psi.preset = psi.preset
  )
  
  fit.args <- modifyList(fit.args, list(...))
  
  if (is.null(fit.args$lower.rho)) fit.args$lower.rho <- -1
  if (is.null(fit.args$upper.rho)) fit.args$upper.rho <- 1
  
  # mvnma() requires at least two pairwise objects
  if (sum(keep) >= 2) {
    fit <- do.call(mvnma, c(data.dir[keep], fit.args))
  }
  else {
    if (any(n.studies > 1))
      warning("Comparison '", treat1, ":", treat2, "' is directly informed ",
              "in fewer than two outcomes; no model-based direct estimate ",
              "available.", call. = FALSE)
    fit <- list()
  }
  
  # attach the raw per-outcome direct-comparison data, the study counts and the
  # outcomes actually fitted, so that pair.nodesplit() can map results back
  attr(fit, "direct.data") <- data.dir
  attr(fit, "n.studies")   <- n.studies
  attr(fit, "keep")        <- keep
  
  fit
}


# direct estimate from an outcome-specific pairwise meta-analysis (default)

direct.metagen <- function(x, treat1, treat2, ...) {
  
  data.dir <- direct.pair(x, t1 = treat1, t2 = treat2)
  
  psi <- extract_het(x)
  
  outcomes <- attr(x, "outcomes")
  n.out    <- length(data.dir)
  level    <- attr(x, "level")
  
  # number of rows contributing direct evidence, per outcome. NA's are now removed.
  n.studies <- vapply(data.dir,
                      function(d) sum(!is.na(d$TE) & !is.na(d$seTE)),
                      integer(1))
  
  keep <- n.studies > 0
  
  fit <- vector("list", n.out)
  names(fit) <- outcomes
  
  for (i in seq_len(n.out)) {
    
    d.i <- data.dir[[i]]
    
    # rows without an estimate for this outcome carry no direct evidence
    sel <- !is.na(d.i$TE) & !is.na(d.i$seTE)
    
    if (!keep[i] || !any(sel)) {
      fit[[i]] <- list(basic_estimates =
                         data.frame(mean = NA_real_, sd = NA_real_,
                                    lower = NA_real_, upper = NA_real_))
      next
    }
    
    psi.i <- if (length(psi) == n.out && !is.na(psi[i])) psi[i] else 0
    
    m.i <- suppressWarnings(
      metagen(TE = d.i$TE[sel], seTE = d.i$seTE[sel],
              studlab = d.i$studlab[sel],
              sm = if (!is.null(attr(x, "sm"))) attr(x, "sm")[i] else "",
              level.ma = level,
              tau.preset = psi.i,
              method.tau.ci = "",
              common = FALSE, random = TRUE))
    
    fit[[i]] <- list(basic_estimates =
                       data.frame(mean = m.i$TE.random,
                                  sd = m.i$seTE.random,
                                  lower = m.i$lower.random,
                                  upper = m.i$upper.random))
  }
  
  attr(fit, "direct.data") <- data.dir
  attr(fit, "n.studies")   <- n.studies
  attr(fit, "keep")        <- keep
  
  fit
}

# create data to calculate direct estimate for a given pair (t1,t2)
direct.pair <- function(x, t1, t2){
  
  if (!inherits(x, "mvnma")) {
    stop("x should be an object of class mvnma")
  }
  
  pair.objects <- attr(x, "pair.objects")
  outcomes     <- attr(x, "outcomes")
  n.out        <- length(outcomes)
  
  data.dir <- vector("list", n.out)
  
  for (i in seq_len(n.out)) {
    
    pair.direct1 <- pair.objects[[i]] %>%
      filter(treat1 == t1, treat2 == t2)
    
    pair.direct2 <- pair.objects[[i]] %>%
      filter(treat1 == t2, treat2 == t1)
    
    if (nrow(pair.direct2) > 0) {
      pair.direct2$TE <- -pair.direct2$TE
      pos <- which(names(pair.direct2) %in% c("treat1", "treat2"))
      names(pair.direct2)[pos] <- c("treat2", "treat1")
    }
    
    # bind_rows() aligns by column name (like rbind.data.frame()), which is
    # what makes the treat1/treat2 relabelling above line the reversed rows
    # up correctly with the direct rows
    data.dir[[i]] <- bind_rows(pair.direct1, pair.direct2) %>%
      select(studlab, treat1, treat2, TE, seTE)
  }
  
  names(data.dir) <- outcomes
  
  data.dir
}

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

# Back-calculate the indirect estimate from the network estimate, the direct
# estimate and the direct evidence proportion, following the SIDE
# back-calculation method (Dias et al., 2010; Koenig et al., 2013) as
# implemented in netmeta. Returns NA when the proportion is within
# tol.direct of 0 or 1, i.e. when the comparison provides essentially only
# indirect or only direct evidence, and when it lies outside (0, 1).
back_calc_indirect <- function(mean.overall, var.overall,
                               mean.dir, prop, tol.direct = 5e-04) {
  
  if (length(prop) != 1 || length(mean.dir) != 1 ||
      is.na(prop) || is.na(mean.dir) ||
      prop <= tol.direct || prop >= 1 - tol.direct) {
    return(data.frame(mean = NA_real_, sd = NA_real_))
  }
  
  mean.indir <- (mean.overall - prop * mean.dir) / (1 - prop)
  
  var.indir <- var.overall / (1 - prop)
  
  data.frame(mean = mean.indir, sd = sqrt(var.indir))
}

z_test_diff <- function(direct.TE, direct.seTE, indirect.TE, indirect.seTE, level = level) {
  
  diff <- direct.TE - indirect.TE
  se.diff <- sqrt(direct.seTE^2 + indirect.seTE^2)
  
  z <- diff / se.diff
  p <- 2 * pnorm(abs(z), lower.tail = FALSE)
  
  z.crit <- qnorm(1 - (1 - level) / 2)
  lower <- diff - z.crit * se.diff
  upper <- diff + z.crit * se.diff
  significant <- p < (1 - level)
  
  res <- data.frame("diff"=diff,"se.diff" = se.diff,
                    "lower.diff" = lower,  "upper.diff" = upper, 
                    "z"=z,"p.val"=p,"significant"=significant)
  
  return(res)
  
}

get_cell <- function(df, treat1, treat2) {
  if (!(treat1 %in% rownames(df)) || !(treat2 %in% colnames(df))) {
    return(NA)
  }
  df[treat1, treat2]
}

# function to extract outcome-specific heterogeneity estimates
extract_het <- function(x,...){
  
  method.model <- attr(x, "method.model")
  
  x <- x[names(x) != "cor"]
  
  if (method.model == "DM") {
    x <- x[names(x) != "sigma"]
  }
  
  het <- vector("list")
  
  for(i in 1:length(x)){
    
    het[[i]] <- x[[i]]$heterogeneity$psi  
  }
  het <- unlist(het)
  
  return(het)
}

# helpers from R package meta used in print.nodesplit
# Copied from meta (unexported); see meta:::formatPT

formatPT <- function(x, lab = FALSE, labval = "p", noblanks = FALSE,
                     digits = 4, zero = TRUE, scientific = FALSE,
                     lab.NA = "--", big.mark = "",
                     JAMA = FALSE) {
  
  if (is.null(x))
    return("")
  
  outdec <- options()$OutDec
  
  n.zeros <- digits - 1
  n.zeros[n.zeros < 0] <- 0
  
  if (!scientific) {
    if (lab) {
      if (!JAMA)
        res <- format(ifelse(is.na(x) | is.nan(x),
                             paste(labval, "=", lab.NA),
                             ifelse(x == 0,
                                    paste(labval, "= 0"),
                                    ifelse(x < 1 / 10^digits,
                                           paste0(labval, " < 0", outdec,
                                                  paste(rep("0",
                                                            n.zeros), collapse = ""),
                                                  "1"),
                                           paste(paste(labval, "="),
                                                 formatC(round(x, digits),
                                                         decimal.mark = outdec,
                                                         big.mark = big.mark,
                                                         format = "f", digits = digits)
                                           )
                                    )
                             )
        )
        )
      else
        res <- format(ifelse(is.na(x) | is.nan(x),
                             paste(labval, "=", lab.NA),
                             ifelse(x < 0.001,
                                    paste0(labval, " < 0", outdec,
                                           paste(rep("0", 2), collapse = ""), "1"),
                                    ifelse(x >= 0.001 & x < 0.01,
                                           paste(paste(labval, "="),
                                                 formatC(x,
                                                         decimal.mark = outdec,
                                                         big.mark = big.mark,
                                                         format = "f", digits = 3)),
                                           ifelse(x >= 0.01 & x <= 0.99,
                                                  paste(paste(labval, "="),
                                                        formatC(x,
                                                                decimal.mark = outdec,
                                                                big.mark = big.mark,
                                                                format = "f", digits = 2)),
                                                  paste(paste(labval, ">"),
                                                        formatC(0.99,
                                                                decimal.mark = outdec,
                                                                big.mark = big.mark,
                                                                format = "f", digits = 2)))
                                    )
                             )
        )
        )
      
    }
    else {
      if (!JAMA)
        res <- format(ifelse(is.na(x) | is.nan(x),
                             lab.NA,
                             ifelse(x == 0,
                                    0,
                                    ifelse(x < 1 / 10^digits,
                                           paste0("< 0", outdec,
                                                  paste(rep("0", n.zeros), collapse = ""),
                                                  "1"),
                                           formatC(round(x, digits),
                                                   decimal.mark = outdec,
                                                   big.mark = big.mark,
                                                   format = "f", digits = digits)
                                    )
                             )
        ),
        justify = "right")
      else
        res <- format(ifelse(is.na(x) | is.nan(x),
                             lab.NA,
                             ifelse(x < 0.001,
                                    paste0("< 0", outdec,
                                           paste(rep("0", 2), collapse = ""), "1"),
                                    ifelse(x >= 0.001 & x < 0.01,
                                           formatC(x,
                                                   decimal.mark = outdec,
                                                   big.mark = big.mark,
                                                   format = "f", digits = 3),
                                           ifelse(x >= 0.01 & x <= 0.99,
                                                  formatC(x,
                                                          decimal.mark = outdec,
                                                          big.mark = big.mark,
                                                          format = "f", digits = 2),
                                                  paste(">",
                                                        formatC(0.99,
                                                                decimal.mark = outdec,
                                                                big.mark = big.mark,
                                                                format = "f", digits = 2)))
                                    )
                             )
        ),
        justify = "right")
    }
  }
  else {
    if (lab)
      res <- format(ifelse(is.na(x) | is.nan(x),
                           paste(labval, "=", lab.NA),
                           paste(labval, "=",
                                 formatC(x, decimal.mark = outdec,
                                         big.mark = big.mark,
                                         format = "e", digits = digits)
                           )
      )
      )
    else
      res <- formatC(x, decimal.mark = outdec,
                     big.mark = big.mark, format = "e", digits = digits)
  }
  #
  if (noblanks)
    res <- gsub(" ", "", res)
  if (!zero)
    res <- gsub("0\\.", "\\.", res)
  #
  # Treat NaNs as NAs
  #
  res[grep("NaN", res)] <- lab.NA
  
  res
}

# Copied from meta (unexported); see meta:::is_relative_effect
is_relative_effect <- function(x)
  x %in% c("HR", "OR", "RR", "IRR", "ROM", "DOR")
