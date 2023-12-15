ttemapit <-
  function(fixed,
           data = parent.frame(),
           kins = NULL,
           id,
           tte,
           maxiter = 500,
           tol = 1e-5,
           taumin = 1e-5,
           taumax = 1e5,
           tauregion = 10,
           verbose = FALSE,
           ...) {
    call <- match.call()
    random.slope <- NULL
    groups <- NULL
    family <- binomial(link = "logit")
    method <- "REML"
    method.optim <- "AI"
    if (!is.null(kins) && !inherits(kins, c("matrix", "list"))) {
      if (is.null(attr(class(kins), "package")))
        stop("Error: \"kins\" must be a matrix or a list.")
      else if (attr(class(kins), "package") != "Matrix")
        stop("Error: if \"kins\" is a sparse matrix, it must be created using the Matrix package.")
    }
    
    if (!id %in% names(data))
      stop("Error: \"id\" must be one of the variables in the names of \"data\".")
    if (inherits(data, "data.frame") &&
        length(class(data)) > 1)
      data <- as.data.frame(data)
    mdl <- model.frame(formula = fixed,
                       data = data,
                       na.action = na.omit)
    time <- data[rownames(mdl), tte]
    idx <- match(rownames(mdl), rownames(model.frame(
      formula = fixed,
      data = data,
      na.action = na.pass
    )))
    multi.pheno <- FALSE
    rm(mdl)
    fit0 <- do.call("glm", list(
      formula = fixed,
      data = data,
      family = family,
      ...
    ))
    
    for (i in 1:length(kins)) {
      match.idx1 <- match(data[idx, id], rownames(kins[[i]]))
      match.idx2 <- match(data[idx, id], colnames(kins[[i]]))
      if (any(is.na(c(match.idx1, match.idx2))))
        stop("Error: kins matrix ",
             i,
             " does not include all individuals in the data.")
      kins[[i]] <- kins[[i]][match.idx1, match.idx2]
    }
    group.id <- rep(1, length(idx))
    time.var <- NULL
    fit <- ttemapit.fit(
      fit0,
      time,
      kins,
      group.id,
      maxiter = maxiter,
      tol = tol,
      taumin = taumin,
      taumax = taumax,
      tauregion = tauregion,
      verbose = verbose
    )
    fit$call <- call
    fit$id_include <- data[idx, id]
    class(fit) <- "tteMAPIT"
    return(fit)
  }

ttemapit.fit <- function(fit0,
                         time,
                         kins,
                         group.id,
                         maxiter = 500,
                         tol = 1e-5,
                         taumin = 1e-5,
                         taumax = 1e5,
                         tauregion = 10,
                         verbose = FALSE) {
  time.var <- NULL
  method = "REML"
  method.optim = "AI"
  group.unique <- unique(group.id)
  group.idx <- list()
  for (i in 1:length(group.unique))
    group.idx[[i]] <- which(group.id == group.unique[i])
  covariance.idx <- fixrho.old <- fixrho.new <- NULL
  n <- nrow(fit0$model[[1]])
  n.pheno <- ncol(fit0$model[[1]])
  q <- length(kins)
  ng <- length(group.idx)
  fixtau.old <- rep(0, length(kins) + ng)
  fit <- ttemapit.ai(
    fit0 = fit0,
    time = time,
    kins = kins,
    covariance.idx = covariance.idx,
    group.idx = group.idx,
    maxiter = maxiter,
    tol = tol,
    verbose = verbose
  )
  fixtau.new <- 1 * (fit$theta < 1.01 * tol)
  while (any(fixtau.new != fixtau.old) ||
         (!is.null(fixrho.new) &&
          any(fixrho.new != fixrho.old))) {
    break
    warning(
      "Variance estimate on the boundary of the parameter space observed, refitting model...",
      call. = FALSE
    )
    fixtau.old <- fixtau.new
    fit <-
      ttemapit.ai(
        fit0 = fit0,
        time = time,
        kins = kins,
        covariance.idx = covariance.idx,
        group.idx = group.idx,
        fixtau = fixtau.old,
        fixrho = fixrho.old,
        maxiter = maxiter,
        tol = tol,
        verbose = verbose
      )
    fixtau.new <- 1 * (fit$theta < 1.01 * tol)
  }
  names(fit$theta)[1] <- "dispersion"
  names(fit$theta)[ng + (1:q)] <- names(kins)
  if (!fit$converged) {
    warning("Average Information REML not converged",
            call. = FALSE)
  }
  return(fit)
}

ttemapit.ai <-
  function(fit0,
           time,
           kins,
           covariance.idx = NULL,
           group.idx,
           tau = rep(0, length(kins) + length(group.idx)),
           fixtau = rep(0, length(kins) + length(group.idx)),
           fixrho = NULL,
           maxiter = 500,
           tol = 1e-5,
           verbose = FALSE) {
    y <- fit0$y
    n <- length(y)
    offset <- fit0$offset
    if (is.null(offset))
      offset <- rep(0, n)
    family <- fit0$family
    eta <- fit0$linear.predictors
    mu <- fit0$fitted.values
    mu.eta <- family$mu.eta(eta)
    Y <- eta - offset + (y - mu) / mu.eta
    sqrtW <-
      mu.eta / sqrt(1 / as.vector(weights(fit0)) * family$variance(mu))
    X <- model.matrix(fit0)
    alpha <- fit0$coef
    if (verbose) {
      cat("Fixed-effect coefficients:\n")
      print(alpha)
    }
    if (family$family %in% c("poisson", "binomial")) {
      tau[1] <- 1
      fixtau[1] <- 1
    }
    q <- length(kins)
    ng <- length(group.idx)
    idxtau <- which(fixtau == 0)
    q2 <- sum(fixtau == 0)
    
    tau[idxtau] <- rep(var(Y) / (q + ng), q2)
    diagSigma <- rep(0, n)
    for (i in 1:ng)
      diagSigma[group.idx[[i]]] <-
      tau[i] / sqrtW[group.idx[[i]]] ^ 2
    Sigma <- diag(diagSigma)
    for (i in 1:q) {
      tau[i + ng] <- tau[i + ng] / mean(diag(kins[[i]]))
      Sigma <- Sigma + tau[i + ng] * kins[[i]]
    }
    Sigma <- Sigma + diag(rep(1e-10, ncol(Sigma)))
    Sigma_i <- chol2inv(chol(Sigma))
    rm(Sigma, diagSigma)
    gc()
    Sigma_iX <- crossprod(Sigma_i, X)
    XSigma_iX <- crossprod(X, Sigma_iX)
    cov <- chol2inv(chol(XSigma_iX))
    Sigma_iXcov <- tcrossprod(Sigma_iX, cov)
    PY <-
      crossprod(Sigma_i, Y) - tcrossprod(Sigma_iX, t(crossprod(Sigma_iXcov, Y)))
    tau0 <- tau
    for (i in 1:q2) {
      APY <- crossprod(kins[[idxtau[i] - ng]], PY)
      PAPY <-
        crossprod(Sigma_i, APY) - tcrossprod(Sigma_iX, t(crossprod(Sigma_iXcov, APY)))
      tau[idxtau[i]] <-
        max(0, tau0[idxtau[i]] + tau0[idxtau[i]] ^ 2 * (sum(Y * PAPY) - (
          sum(Sigma_i * kins[[idxtau[i] - ng]]) - sum(Sigma_iX * crossprod(kins[[idxtau[i] -
                                                                                   ng]], Sigma_iXcov))
        )) / n)
    }
    
    for (i in seq_len(maxiter)) {
      if (verbose) {
        cat("\nIteration ", i, ":\n")
      }
      alpha0 <- alpha
      tau0 <- tau
      fit <- .Call(C_fitglmm_ai,
                   Y,
                   X,
                   q,
                   kins,
                   ng,
                   group.idx,
                   sqrtW ^ 2,
                   tau,
                   fixtau)
      
      
      Dtau <- as.numeric(fit$Dtau)
      tau[idxtau] <- tau0[idxtau] + Dtau
      tau[tau < tol & tau0 < tol] <- 0
      while (any(tau < 0)) {
        Dtau <- Dtau / 2
        tau[idxtau] <- tau0[idxtau] + Dtau
        tau[tau < tol & tau0 < tol] <- 0
      }
      tau[tau < tol] <- 0
      
      cov <- as.matrix(fit$cov)
      alpha <- as.numeric(fit$alpha)
      eta <- as.numeric(fit$eta) + offset
      if (verbose) {
        cat("Variance component estimates:\n")
        print(tau)
        cat("Fixed-effect coefficients:\n")
        print(alpha)
      }
      # INSERT SURVIVAL
      Lambda0 <- GetLambda0(eta, y, time) + 1e-5 # TODO: Fix Lambda0
      mu <- Lambda0 * exp(eta)
      Y <- eta - offset + (y - mu) / mu
      sqrtW <- as.vector(sqrt(mu))
      # END INSERTION
      # mu <- family$linkinv(eta)
      # mu.eta <- family$mu.eta(eta)
      # Y <- eta - offset + (y - mu) / mu.eta
      # sqrtW <- mu.eta / sqrt(1 / as.vector(weights(fit0)) * family$variance(mu))
      if (2 * max(abs(alpha - alpha0) / (abs(alpha) + abs(alpha0) + tol),
                  abs(tau - tau0) / (abs(tau) + abs(tau0) + tol)) < tol) {
        break
      }
      if (max(abs(tau)) > tol ^ (-2)) {
        warning(
          "Large variance estimate observed in the iterations, model not converged...",
          call. = FALSE
        )
        i <- maxiter
        break
      }
    }
    converged <- ifelse(i < maxiter, TRUE, FALSE)
    res <- y - mu
    res.var <- rep(1, n)
    for (i in 1:ng)
      res.var[group.idx[[i]]] <- tau[i]
    fit$Sigma_i <- fit$Sigma_iX <- NULL
    names(alpha) <- names(fit0$coef)
    rownames(cov) <- rownames(vcov(fit0))
    colnames(cov) <- colnames(vcov(fit0))
    return(
      list(
        theta = tau,
        n.pheno = 1,
        n.groups = ng,
        coefficients = alpha,
        linear.predictors = eta,
        fitted.values = mu,
        Y = Y,
        X = X,
        P = fit$P,
        logLik = fit$logLik,
        residuals = res,
        scaled.residuals = res * as.vector(weights(fit0)) / res.var,
        cov = cov,
        Sigma_i = fit$Sigma_i,
        Sigma_iX = fit$Sigma_iX,
        converged = converged
      )
    )
  }

GetdenominLambda0 <- function(caseIndexwithTies,
                             lin.pred.new,
                             newIndexWithTies) {
  explin <- exp(lin.pred.new)
  getdenom <- function(i,
                      caseIndexwithTies,
                      explin,
                      newIndexWithTies) {
    nc <- length(explin)
    x <- 1 / sum(explin[which(newIndexWithTies >= caseIndexwithTies[i])])
    return(x)
  }
  demonVec <- sapply(seq(1, length(caseIndexwithTies)),
                    getdenom,
                    caseIndexwithTies,
                    explin,
                    newIndexWithTies)
  return(demonVec)
}

GetLambda0 <- function(lin.pred, status, time)
{
  inC <- GetIndexofCases(status, time)
  lin.pred.new <- lin.pred[inC$timedata$orgIndex]
  demonVec <- GetdenominLambda0(inC$caseIndexwithTies,
                               lin.pred.new,
                               inC$timedata$newIndexWithTies)
  Lambda0 <- rep(0, length(lin.pred))
  lambda0 <- 0
  for (i in 1:length(lin.pred))
  {
    j <- inC$timedata$newIndexWithTies[i]
    lambda0 <- sum(demonVec[which(inC$caseIndexwithTies <= j)])
    Lambda0[inC$timedata$orgIndex[i]] <- lambda0
  }
  return(Lambda0)
  
}

GetIndexofCases <- function(status, time) {
  timedata <- data.frame(time = time,
                        orgIndex = seq(1, length(time)),
                        status = status)
  timedata <- timedata[order(timedata$time), ]
  timedata$newIndex <- seq(1, length(time))
  caseIndex <- which(timedata$status == 1)
  caseIndexwithTies <- caseIndex
  for (i in 2:length(caseIndex)) {
    if (timedata$time[caseIndex[i]] == timedata$time[caseIndex[i - 1]]) {
      caseIndexwithTies[i] = caseIndexwithTies[i - 1]
    }
  }
  
  uniqTimeIndex <- unique(caseIndexwithTies)
  uniqTimeVec <- timedata$time[uniqTimeIndex]
  uniqTimeData <- data.frame(uniqTimeVec = uniqTimeVec, uniqTimeIndex = uniqTimeIndex)
  timedata$newIndexWithTies <- timedata$newIndex
  
  for (i in 1:nrow(timedata)) {
    if (timedata$time[i] %in% uniqTimeVec) {
      timedata$newIndexWithTies[i] <- uniqTimeIndex[which(uniqTimeVec == timedata$time[i])]
    }
  }
  
  return(
    list(
      timedata = timedata,
      caseIndex = caseIndex,
      caseIndexwithTies = caseIndexwithTies,
      uniqTimeIndex = uniqTimeIndex
    )
  )
}
