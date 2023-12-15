#' Simulate time to event data
#'
#' @param genotype_matrix Genotype matrix with samples as rows, and SNPs as columns.
#' @param covariate_matrix Matrix with fixed effect covariates
#' @param n_causal Number of SNPs that are causal.
#' @param n_epistatic Number of SNPs that are epistatic
#' @param group_ratio Ratio of sizes of groups that interact, e.g. a ratio 1:3 would be value 3.
#' @param heritability Broad-sense heritability.
#' @param rho Proportion of heritability explained by additivity.
#' @param maxT Maximum time to be observed.
#' @param scale_effects Scale effects
#' @param seed Random seed for simulation.
#' @param logLevel is a string parameter defining the log level for the logging package.
#' @param logFile is a string parameter defining the name of the log file for the logging output.
#' @param maf_threshold is a float parameter defining the threshold for the minor allele frequency not included in causal SNPs.
#' @return A list object containing the trait data, the genotype data, as well as the causal SNPs and summary statistics.
#' @export
#' @import checkmate
#' @importFrom coxed sim.survdata
#' @import dplyr
#' @import foreach
#' @import parallel
#' @importFrom stats var cor sd complete.cases
#' @importFrom utils head
simulate_tte_data <- function(
    genotype_matrix, covariate_matrix, n_causal = 1000, n_epistatic = 20,
    heritability = 0.6, rho = 0.8, maxT = 100, scale_effects = 1,
    group_ratio = 1, maf_threshold = 0.01, seed = 67132,
    logLevel = "INFO", logFile = NULL
) {
  
  set.seed(seed)
  coll <- makeAssertCollection()
  assertInt(n_causal, lower = 0, add = coll)
  assertInt(n_epistatic, lower = 0, add = coll)
  assertDouble(group_ratio, lower = 1, add = coll)
  assertDouble(heritability, lower = 0, upper = 1, add = coll)
  assertDouble(rho, lower = 0, upper = 1, add = coll)
  assertDouble(maf_threshold, lower = 0, upper = 1, add = coll)
  assertInt(seed, lower = 1, add = coll)
  assertMatrix(genotype_matrix, all.missing = FALSE, add = coll)
  reportAssertions(coll)
  
  logging::logReset()
  logging::basicConfig(level = logLevel)
  log <- logging::getLogger("simulate_traits")
  if (!is.null(logFile)) {
    filePath <- file.path(getwd(), logFile)
    log$debug("Logging to file: %s", filePath)
    log$addHandler(logging::writeToFile, file = filePath)
  }
  
  snp.ids <- 1:ncol(genotype_matrix)
  maf <- colMeans(genotype_matrix) / 2
  X <- scale(genotype_matrix)
  maf_compliant <- (maf > maf_threshold) & (maf < 1 - maf_threshold)
  # scale produces NaN when the columns have zero variance
  snp.ids.filtered <- snp.ids[complete.cases(t(X)) & maf_compliant]
  
  n_samples <- nrow(X)  # number of genotype samples
  n_snp <- length(snp.ids.filtered)  # number of SNPs passing quality control
  log$debug("Scaled genotype matrix: %d x %d", n_samples, n_snp)
  log$debug(
    "Disregard %d variants due to zero variance or small minor allele frequency.",
    ncol(genotype_matrix) - length(snp.ids.filtered)
  )
  log$debug("Minor allele frequency threshold %f.", maf_threshold)
  
  # divide groups into ratios
  n_group1 <- ceiling(n_epistatic / (1 + group_ratio))
  
  coll <- makeAssertCollection()
  assertInt(n_causal, lower = 0, upper = n_snp, add = coll)
  assertInt(n_group1, lower = 0, upper = n_epistatic, add = coll)
  reportAssertions(coll)
  
  # factor vectors for splitting the groups
  f_trait <- get_factors(n_group1, n_epistatic)
  
  
  log$debug("Number of causal SNPs: %d", n_causal)
  log$debug("NA in raw genotype matrix: %d", sum(is.na(genotype_matrix)))
  log$debug("NA in scaled genotype matrix: %d", sum(is.na(X)))
  
  Y <- c()
  causal_snps <- list()
  
  colnames(genotype_matrix) <- seq_len(ncol(genotype_matrix)) %>%
    sprintf(fmt = "snp_%05d")  # column names names for SNPs
  
  
  ## select causal SNPs
  causal_snps <- sample(snp.ids.filtered, n_causal, replace = F)
  epistatic_snps <- sample(causal_snps, n_epistatic, replace = F)
  epistatic_grouped <- split(epistatic_snps, f_trait)
  group1 <- epistatic_grouped$group1
  group2 <- epistatic_grouped$group2
  log$debug("Length causal set: %d", length(causal_snps))
  log$debug("Length trait specific set 1: %d", length(group1))
  log$debug("Length trait specific set 2: %d", length(group2))
  
  log$debug("Head of causal SNPs: %s", head(causal_snps))
  log$debug("Head of trait specific SNPs group 1: %s", head(group2))
  log$debug("Head of trait specific SNPs group 2: %s", head(group2))
  
  # create trait_specific interaction matrix
  X_causal <- X[, causal_snps]  # all SNPs have additive effects
  X_group1 <- X[, group1, drop = FALSE]
  X_group2 <- X[, group2, drop = FALSE]
  
  start_interactions <- proc.time()
  log$debug("Computing interactions. This may take a while.")
  i <- NULL
  X_epi <- foreach(i = seq_len(length(group1)), .combine = cbind) %do% {
    X_group1[, i] * X_group2
  }
  
  time_interactions <- proc.time() - start_interactions
  log$debug("Interactions X_epi computed in %f", time_interactions[3])
  log$debug(
    "Dimension of interaction matrix X_epi: %d x %d", nrow(X_epi),
    ncol(X_epi)
  )
  
  # marginal effects
  X_marginal <- X_causal
  beta <- rnorm(ncol(X_marginal))
  y_marginal <- X_marginal %*% beta
  scaled_beta <- beta * sqrt(heritability * rho/c(var(y_marginal)))
  
  # pairwise epistatic effects
  alpha <- rnorm(ncol(X_epi))
  n_epistatic_effects <- n_group1 * (n_epistatic - n_group1)
  if (n_epistatic_effects > 0) {
    y_epi <- X_epi %*% alpha
    scaled_alpha <- alpha * sqrt(heritability * (1 - rho)/c(var(y_epi)))
  } else {
    y_epi <- 0 * y_marginal
    log$debug("y_epi vector of zeros size y_marginal: %s", y_epi)
    scaled_alpha <- rep(0, ncol(X_epi))
  }
  
  colnames(genotype_matrix)[causal_snps] <- paste0(colnames(genotype_matrix[, causal_snps]),
                                                    rep("_add",
                                                        length(causal_snps)))
  colnames(genotype_matrix)[epistatic_snps] <- paste0(colnames(genotype_matrix[, epistatic_snps]),
                                                     rep("_epi",
                                                         length(epistatic_snps)))
  
  # covariate variation
  gamma <- rnorm(ncol(covariate_matrix))
  y_fixed <- covariate_matrix %*% gamma
  scaled_gamma <- gamma * sqrt((1 - heritability)/c(var(y_fixed)))
  
  COVARIATES <- cbind(covariate_matrix, X_marginal, X_epi)
  EFFECTS <- c(scaled_gamma, scaled_beta, scaled_alpha) * scale_effects
  simdata <- sim.survdata(N=ncol(COVARIATES), T=maxT, num.data.frames=1, X=COVARIATES, beta=EFFECTS)
  survival_time <-simdata$data$y
  event_indicator <- as.numeric(simdata$data$failed)
  
  results <- list(
    survival_time = survival_time,
    event_indicator = event_indicator,
    covariates = covariate_matrix,
    genotypes = genotype_matrix,
    cov_effects = scaled_gamma,
    add_effects = scaled_beta,
    epi_effects = scaled_alpha,
    sim.data = simdata
  )
  
  return(results)
}

get_factors <- function(n1, n) {
  return(
    c(
      rep("group1", n1),
      rep("group2", n - n1)
    )
  )
}

