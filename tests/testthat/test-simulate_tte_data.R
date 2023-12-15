test_that("simulate tte data", {
  # given
  c <- 3
  p <- 200
  n <- 50
  maxT <- 100
  n_causal <- 5
  n_epistatic <- 2
  genotype_matrix <- matrix(
    sample(1:2, p * n, replace = TRUE),
    ncol = p
  )
  covariate_matrix <- matrix(
    runif(c * n),
    ncol = c
  )
  # when
  sims <- simulate_tte_data(
    genotype_matrix, covariate_matrix, n_causal = n_causal, n_epistatic = n_epistatic,
    heritability = 0.6, rho = 0.8, maxT = maxT, scale_effects = 1e-2,
    group_ratio = 1, maf_threshold = 0.01, seed = 67132,
    logLevel = "INFO", logFile = NULL
  )
  # then
  print(sims$survival_time)
})
