test_that("tteMAPIT dev", {
	data(example)
	set.seed(123)
	pheno <- example$pheno
	n_samples <- nrow(pheno)
	time <- rpois(n_samples, 30)
	n_missing <- 20
	pheno$id <- 1:n_samples
	pheno$tte <- time
	pheno$disease[sample(1:n_samples,n_missing)] <- NA
	pheno$age[sample(1:n_samples,n_missing)] <- NA
	pheno$sex[sample(1:n_samples,n_missing)] <- NA
	GRM <- example$GRM
	rand <- sample(nrow(GRM))
	EPI <- GRM[rand, rand]
	colnames(EPI) <- pheno$id
	rownames(EPI) <- pheno$id
	kins <- list(GRM = GRM, EPI = EPI)
	obj <- ttemapit(disease ~ age + sex, data = pheno, kins = kins, id = "id", tte = "tte")
  print(obj$theta)
  print(obj$coef)
})

