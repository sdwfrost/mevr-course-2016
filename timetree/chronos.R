## Takes a rooted tree as produced by RLRootToTip.jar and the
## coefficients produced by it, and pins the tips and the root to
## their known/estimated values.

## January 6, 2013: modified for processing Art's tree.

## January 6, 2013: modified from PinTips_art.r to quickly do an
## ad-hoc strict clock fitting.

require(ape)
require(adephylo)

RLchronos <- function (phy, lambda = 1, model = "correlated", quiet = FALSE, 
                       calibration = makeChronosCalib(phy),
                       control = chronos.control()) 
{
  model <- match.arg(tolower(model), c("correlated", "relaxed", 
                                       "discrete"))
  n <- Ntip(phy)
  ROOT <- n + 1L
  m <- phy$Nnode
  el <- phy$edge.length
  if (any(el < 0)) 
    stop("some branch lengths are negative")
  e1 <- phy$edge[, 1L]
  e2 <- phy$edge[, 2L]
  N <- length(e1)
  TIPS <- 1:n
  EDGES <- 1:N
  tol <- control$tol
  node <- calibration$node
  age.min <- calibration$age.min
  age.max <- calibration$age.max
  if (model == "correlated") {
    basal <- which(e1 == ROOT)
    Nbasal <- length(basal)
    ind1 <- EDGES[-basal]
    ind2 <- match(e1[EDGES[-basal]], e2)
  }
  age <- numeric(n + m)
  if (!quiet) 
    cat("\nSetting initial dates...\n")
  
  ## seq.nod is sorted by node index.
  seq.nod <- .Call(seq_root2tip, phy$edge, n, phy$Nnode)
  ii <- 1L
  repeat {
    ini.time <- age
    ini.time[ROOT:(n + m)] <- NA
    ini.time[node] <- if (is.null(age.max)) 
      age.min
    else runif(length(node), age.min, age.max)
    if (is.na(ini.time[ROOT])) 
      ini.time[ROOT] <- if (is.null(age.max)) 
        3 * max(age.min)
      else 3 * max(age.max)
    ISnotNA.ALL <- unlist(lapply(seq.nod, function(x) sum(!is.na(ini.time[x]))))
    
    ## o is a vector of tip indices, sorted by those which have the
    ## most non-NA nodes between them and the root.
    o <- order(ISnotNA.ALL, decreasing = TRUE)
    ## seq.nod[o] is the vector of nodes visited on the path from
    ## root to tip, starting with the root and ending with the
    ## tip.  We need to modify it so that it handles the more
    ## recent tips first.
    
    not.NA.tip.ages <- ini.time[o]
    
    tips.by.age <- order(not.NA.tip.ages, decreasing=TRUE)
    
    for (y in seq.nod[tips.by.age]) {
      ## Are any un-set still?
      ISNA <- is.na(ini.time[y])
      if (any(ISNA)) {
        ## i == 1 is the root, and we already set that.
        i <- 2L
        while (i <= length(y)) {
          if (ISNA[i]) {
            ## if the ith node on this path is unset, then find the
            ## next *set* node; j is its index.
            j <- i + 1L
            while (ISNA[j]) j <- j + 1L
            ## nb.val is now the number of unset nodes on this path.
            nb.val <- j - i
            ## by is the interval of time between the date-set nodes
            ## divided by the number of edges on the path between
            ## them (which is the number of nodes + 1).
            by <- (ini.time[y[i - 1L]] - ini.time[y[j]])/(nb.val + 1)
            ini.time[y[i:(j - 1L)]] <-
              ini.time[y[i - 1L]] - by * seq_len(nb.val)
            i <- j + 1L
          }
          else i <- i + 1L
        }
      }
    }
    if (all(ini.time[e1] - ini.time[e2] >= 0)) 
      break
    ii <- ii + 1L
    if (ii > 1000) 
      stop("cannot find reasonable starting dates after 1000 tries:\nmaybe you need to adjust the calibration dates")
  }
  cat("Initial dates set.\n")
  ini.rate <- el/(ini.time[e1] - ini.time[e2])
  if (model == "discrete") {
    Nb.rates <- control$nb.rate.cat
    minmax <- range(ini.rate)
    if (Nb.rates == 1) {
      ini.rate <- sum(minmax)/2
    }
    else {
      inc <- diff(minmax)/Nb.rates
      ini.rate <- seq(minmax[1] + inc/2, minmax[2] - inc/2, 
                      inc)
      ini.freq <- rep(1/Nb.rates, Nb.rates - 1)
      lower.freq <- rep(0, Nb.rates - 1)
      upper.freq <- rep(1, Nb.rates - 1)
    }
  } else Nb.rates <- N
  
  
  ## The internal nodes are numbered (n+1), ..., (n+m).
  unknown.ages <- 1:m + n
  ## These appear to only be for internal nodes.
  lower.age <- rep(tol, m)
  upper.age <- rep(1/tol, m)
  
  bdd.internal.node.indices <- which(node > n)
  bdd.internal.nodes <- node[bdd.internal.node.indices]
  lower.age[bdd.internal.nodes - n] <- age.min[bdd.internal.node.indices]
  upper.age[bdd.internal.nodes - n] <- age.max[bdd.internal.node.indices]
  
  ## Eliminate nodes where age.min != age.max; whatever's left,
  ## we can set the age from age.min.
  ii <- which(age.min != age.max)
  ## if (length(ii)) {
  ##     node <- node[-ii]
  ##     age.min <- age.min[-ii]
  ## }
  if (length(node[-ii])) {
    age[node[-ii]] <- age.min[-ii]
  } else {
    age[node] <- age.min
  }
  
  
  ## Get rid of the internal nodes that we have knowledge of.
  fixed.internal.nodes <-
    bdd.internal.nodes[age.min[bdd.internal.node.indices] ==
                       age.max[bdd.internal.node.indices]]
  if (length(fixed.internal.nodes)) {
    unknown.ages <- unknown.ages[n - fixed.internal.nodes]
    lower.age <- lower.age[n - fixed.internal.nodes]
    upper.age <- upper.age[n - fixed.internal.nodes]
  }
  
  known.ages <- c(TIPS, bdd.internal.nodes)
  lower.rate <- rep(tol, Nb.rates)
  upper.rate <- rep(100 - tol, Nb.rates)
  degree_node <- tabulate(phy$edge)
  eta_i <- degree_node[e1]
  eta_i[e2 <= n] <- 1L
  X <- vector("list", N)
  for (i in EDGES) {
    j <- integer()
    if (e1[i] != ROOT) 
      j <- c(j, which(e2 == e1[i]))
    if (e2[i] >= n) 
      j <- c(j, which(e1 == e2[i]))
    X[[i]] <- j
  }
  ## List of indices in e2 of nodes with unknown ages.
  D_ki <- match(unknown.ages, e2)
  ## List of indices in e1 of edges coming out of each of the
  ## unknown nodes.
  A_ki <- lapply(unknown.ages, function(x) which(x == e1))
  gradient.poisson <- function(rate, node.time) {
    age[unknown.ages] <- node.time
    real.edge.length <- age[e1] - age[e2]
    ## el is the list of edge lengths in the tree.
    gr <- el/rate - real.edge.length
    tmp <- el/real.edge.length - rate
    gr.dates <- sapply(A_ki, function(x) sum(tmp[x])) - tmp[D_ki]
    c(gr, gr.dates)
  }
  gradient <- switch(model, correlated = function(rate, node.time) {
    gr <- gradient.poisson(rate, node.time)
    gr[RATE] <- gr[RATE] - lambda * 2 * (eta_i * rate -
                                         sapply(X, 
                                                function(x) sum(rate[x])))
    if (Nbasal == 2) {
      i <- basal[1]
      j <- basal[2]
      gr[i] <- gr[i] - lambda * (rate[i] - rate[j])
      gr[j] <- gr[j] - lambda * (rate[j] - rate[i])
    } else {
      for (i in 1:Nbasal) j <- basal[i]
      gr[j] <- gr[j] - lambda * 2 * (rate[j] * (1 - 1/Nbasal) - 
                                     sum(rate[basal[-i]])/Nbasal)/(Nbasal - 1)
    }
    gr
  }, relaxed = function(rate, node.time) {
    gr <- gradient.poisson(rate, node.time)
    mean.rate <- mean(rate)
    gr[RATE] <- gr[RATE] + lambda * 2 * dgamma(rate, mean.rate) * 
      (rank(rate)/Nb.rates - pgamma(rate, mean.rate))
    gr
  }, discrete = NULL)
  log.lik.poisson <- function(rate, node.time) {
    age[unknown.ages] <- node.time
    real.edge.length <- age[e1] - age[e2]
    if (isTRUE(any(real.edge.length < 0))) 
      return(-1e+100)
    B <- rate * real.edge.length
    sum(el * log(B) - B - lfactorial(el))
  }
  penal.loglik <-
    switch(model,
           correlated = function(rate, node.time)
           {
             loglik <- log.lik.poisson(rate, node.time)
             if (!is.finite(loglik)) return(-1e+100)
             loglik - lambda * (sum((rate[ind1] - rate[ind2])^2) + 
                                var(rate[basal]))
           },
           relaxed = function(rate, node.time)
           {
             loglik <- log.lik.poisson(rate, node.time)
             if (!is.finite(loglik)) return(-1e+100)
             mu <- mean(rate)
             loglik - lambda * sum((1:N/N - pgamma(sort(rate), mean(rate)))^2)
           },
           discrete = if (Nb.rates == 1) function(rate, node.time) log.lik.poisson(rate, 
                            node.time) else function(rate, node.time, freq) {
                              if (isTRUE(sum(freq) > 1)) return(-1e+100)
                              rate.freq <- sum(c(freq, 1 - sum(freq)) * rate)
                              log.lik.poisson(rate.freq, node.time)
                            })
  opt.ctrl <- list(eval.max = control$eval.max, iter.max = control$iter.max)

  ## In the optimization, the rates come first and the ages come second.
  ## p is the vector of parameters.
  ## RATE and AGE are the corresponding indices in the parameter;
  ## LOW and UP are bounds on the parameters.
  RATE <- 1:Nb.rates
  AGE <- Nb.rates + 1:length(unknown.ages)
  if (model == "discrete") {
    if (Nb.rates == 1) {
      start.para <- c(ini.rate, ini.time[unknown.ages])
      f <- function(p) -penal.loglik(p[RATE], p[AGE])
      g <- NULL
      LOW <- c(lower.rate, lower.age)
      UP <- c(upper.rate, upper.age)
    }
    else {
      FREQ <- length(RATE) + length(AGE) + 1:(Nb.rates - 1)
      start.para <- c(ini.rate, ini.time[unknown.ages], 
                      ini.freq)
      f <- function(p) -penal.loglik(p[RATE], p[AGE], p[FREQ])
      g <- NULL
      LOW <- c(lower.rate, lower.age, lower.freq)
      UP <- c(upper.rate, upper.age, upper.freq)
    }
  } else {
    start.para <- c(ini.rate, ini.time[unknown.ages])
    f <- function(p) -penal.loglik(p[RATE], p[AGE])
    g <- function(p) -gradient(p[RATE], p[AGE])
    LOW <- c(lower.rate, lower.age)
    UP <- c(upper.rate, upper.age)
  }
  k <- length(LOW)
  if (!quiet) 
    cat("Fitting in progress... get a first set of estimates\n")
  out <- nlminb(start.para, f, g, control = opt.ctrl, lower = LOW, 
                upper = UP)
  if (model == "discrete") {
    if (Nb.rates == 1) {
      f.rates <- function(p) -penal.loglik(p, current.ages)
      f.ages <- function(p) -penal.loglik(current.rates, p)
    }
    else {
      f.rates <- function(p) -penal.loglik(p, current.ages, 
                                           current.freqs)
      f.ages <- function(p) -penal.loglik(current.rates, 
                                          p, current.freqs)
      f.freqs <- function(p) -penal.loglik(current.rates, 
                                           current.ages, p)
      g.freqs <- NULL
    }
    g.rates <- NULL
    g.ages <- NULL
  } else {
    f.rates <- function(p) -penal.loglik(p, current.ages)
    g.rates <- function(p) -gradient(p, current.ages)[RATE]
    f.ages <- function(p) -penal.loglik(current.rates, p)
    g.ages <- function(p) -gradient(current.rates, p)[AGE]
  }
  current.ploglik <- -out$objective
  current.rates <- out$par[RATE]
  current.ages <- out$par[AGE]
  if (model == "discrete" && Nb.rates > 1) 
    current.freqs <- out$par[FREQ]
  dual.iter.max <- control$dual.iter.max
  i <- 0L
  if (!quiet) 
    cat("         Penalised log-lik =", current.ploglik, 
        "\n")
  repeat {
    if (dual.iter.max < 1) 
      break
    if (!quiet) 
      cat("Optimising rates...")
    out.rates <- nlminb(current.rates, f.rates, g.rates, 
                        control = list(eval.max = 1000, iter.max = 1000, 
                          step.min = 1e-08, step.max = 0.1), lower = lower.rate, 
                        upper = upper.rate)
    new.rates <- out.rates$par
    if (-out.rates$objective > current.ploglik) 
      current.rates <- new.rates
    if (model == "discrete" && Nb.rates > 1) {
      if (!quiet) 
        cat(" frequencies...")
      out.freqs <- nlminb(current.freqs, f.freqs,
                          control = list(eval.max = 1000, 
                            iter.max = 1000, step.min = 0.001, step.max = 0.5), 
                          lower = lower.freq, upper = upper.freq)
      new.freqs <- out.freqs$par
    }
    if (!quiet) 
      cat(" dates...")
    out.ages <- nlminb(current.ages, f.ages, g.ages,
                       control = list(eval.max = 1000, 
                         iter.max = 1000, step.min = 0.001, step.max = 100), 
                       lower = lower.age, upper = upper.age)
    new.ploglik <- -out.ages$objective
    if (!quiet) 
      cat("", current.ploglik, "\n")
    if (new.ploglik - current.ploglik > 1e-06 && i <= dual.iter.max) {
      current.ploglik <- new.ploglik
      current.rates <- new.rates
      current.ages <- out.ages$par
      if (model == "discrete" && Nb.rates > 1) 
        current.freqs <- new.freqs
      out <- out.ages
      i <- i + 1L
    }
    else break
  }
  if (!quiet) 
    cat("\nDone.\n")
  if (model == "discrete") {
    rate.freq <- if (Nb.rates == 1) 
      current.rates
    else mean(c(current.freqs, 1 - sum(current.freqs)) * 
              current.rates)
    logLik <- log.lik.poisson(rate.freq, current.ages)
    PHIIC <- list(logLik = logLik, k = k, PHIIC = -2 * logLik + 
                  2 * k)
  } else {
    logLik <- log.lik.poisson(current.rates, current.ages)
    PHI <- switch(model,
                  correlated = (current.rates[ind1] - 
                                current.rates[ind2])^2 + var(current.rates[basal]), 
                  relaxed = (1:N/N - pgamma(sort(current.rates), mean(current.rates)))^2)
    PHIIC <- list(logLik = logLik, k = k, lambda = lambda, 
                  PHIIC = -2 * logLik + 2 * k + lambda * svd(PHI)$d)
  }
  ## DEBUGGING
  ## attr(phy, "call") <- "FOOBAR"
  attr(phy, "call") <- match.call()
  attr(phy, "ploglik") <- -out$objective
  attr(phy, "rates") <- current.rates
  if (model == "discrete" && Nb.rates > 1) 
    attr(phy, "frequencies") <- current.freqs
  attr(phy, "message") <- out$message
  attr(phy, "PHIIC") <- PHIIC
  age[unknown.ages] <- current.ages
  phy$edge.length <- age[e1] - age[e2]
  class(phy) <- c("chronos", class(phy))
  phy
}

