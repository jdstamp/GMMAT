test_that("tteMAPIT dev", {
  c <- 3
  p <- 300
  n <- 700
  maxT <- 100
  n_causal <- 20
  n_epistatic <- 4
  genotype_matrix <- matrix(
    sample(0:2, p * n, replace = TRUE),
    ncol = p
  )
  covariate_matrix <- matrix(
    runif(c * n),
    ncol = c
  )
  # when
  sims <- simulate_tte_data(
    genotype_matrix, covariate_matrix, n_causal = n_causal, n_epistatic = n_epistatic,
    heritability = 0.9, rho = 0.2, maxT = maxT, scale_effects = 1e-2,
    group_ratio = 1, maf_threshold = 0.01, seed = 67132,
    logLevel = "INFO", logFile = NULL
  )
  survival_time <-sims$survival_time
  event_indicator <- as.numeric(sims$event_indicator)
  pheno <- data.frame(status = event_indicator,
                      tte = survival_time,
                      age = covariate_matrix[, 1],
                      sex = covariate_matrix[, 2],
                      income = covariate_matrix[, 3])
  pheno$id <- 1:n
  X <- sims$genotypes[which(apply(sims$genotypes, 1, var) != 0), ]
  target_string <- "epi"
  matching_columns <- grep(target_string, colnames(X), ignore.case = TRUE)
  xk <- X[, matching_columns[1]]
  X <- X[, -matching_columns[1]]
  Xsd <- apply(X, 1, sd)
  Xmean <- apply(X, 1, mean)
  X <- (X - Xmean) / Xsd
  GRM <- X %*% t(X) / ncol(X)
  GRM <- pmax(GRM, t(GRM))
  EPI <- t( xk * t(GRM * xk))
  EPI <- pmax(EPI, t(EPI))
  ERR <- diag(nrow(EPI))
  colnames(GRM) <- pheno$id
  rownames(GRM) <- pheno$id
  colnames(EPI) <- pheno$id
  rownames(EPI) <- pheno$id
  colnames(ERR) <- pheno$id
  rownames(ERR) <- pheno$id
  kins <- list(GRM = GRM, EPI = EPI)
  kins2 <- list(GRM = GRM)
  pheno$xk <- xk
  obj <- ttemapit(status ~ age + sex + income + xk, data = pheno, kins = kins, id = "id", tte = "tte", maxiter = 500, tol = 1e-3)
  obj2 <- ttemapit(status ~ age + sex + income + xk, data = pheno, kins = kins2, id = "id", tte = "tte", maxiter = 500, tol = 1e-3)
  RSS1 <- sum((obj$res)^2)
  RSS0 <- sum((obj2$res)^2)
  p1 <- length(kins) + length(obj$coef)
  p0 <- length(kins2) + length(obj2$coef)
  print(str(sims$cov_effects))
  print("H1: ------")
  print(obj$theta)
  print(obj$coef)
  print(sprintf("RSS: %f", RSS1))
  print("H0: ------")
  print(obj2$theta)
  print(obj2$coef)
  print(sprintf("RSS: %f", RSS0))
  F <- ((RSS0 - RSS1) / (p1 - p0)) / (RSS1 / (n - p1))
  print(sprintf("F: %f", F))
  pv <- pf(F, p1 - p0, n - p1, lower.tail = FALSE)
  print(sprintf("p-value: %f", pv))
})
