test_that("cross-sectional id le 400 binomial", {
	skip_on_cran()

	data(example)
	suppressWarnings(RNGversion("3.5.0"))
	set.seed(123)
	pheno <- rbind(example$pheno, example$pheno[1:100, ])
	pheno$id <- 1:500
	pheno$disease[sample(1:500,20)] <- NA
	pheno$age[sample(1:500,20)] <- NA
	pheno$sex[sample(1:500,20)] <- NA
	pheno <- pheno[sample(1:500,450), ]
	pheno <- pheno[pheno$id <= 400, ]
	kins <- example$GRM
	obj5 <- glmmkin(disease ~ age + sex, data = pheno, kins = kins, id = "id", family = binomial(link = "logit"), method = "REML", method.optim = "AI")
	expect_equal(signif(as.numeric(obj5$theta)), signif(c(1, 0.1925225)))
	expect_equal(signif(as.numeric(obj5$coef)), signif(c(1.01676230, -0.01506251, -0.33240659)))
	obj6 <- glmmkin(disease ~ age + sex, data = pheno, kins = NULL, id = "id", family = binomial(link = "logit"), method = "REML", method.optim = "AI")
	expect_equal(c(obj6$theta), 1)
	expect_equal(signif(as.numeric(obj6$coef)), signif(c(0.94253958, -0.01429532, -0.32823930)))
	obj <- glm(disease ~ age + sex, data = pheno, family = binomial(link = "logit"))
	expect_equal(obj6$coef, obj$coef)

	idx <- sample(nrow(pheno))
	pheno <- pheno[idx, ]
	obj <- glmmkin(disease ~ age + sex, data = pheno, kins = kins, id = "id", family = binomial(link = "logit"), method = "REML", method.optim = "AI")
	expect_equal(obj5$theta, obj$theta)
	expect_equal(obj5$coef, obj$coef)
	obj <- glmmkin(disease ~ age + sex, data = pheno, kins = NULL, id = "id", family = binomial(link = "logit"), method = "REML", method.optim = "AI")
	expect_equal(obj6$theta, obj$theta)
	expect_equal(obj6$coef, obj$coef)
	obj <- glm(disease ~ age + sex, data = pheno, family = binomial(link = "logit"))
	expect_equal(obj6$coef, obj$coef)

	idx <- sample(nrow(kins))
	kins <- kins[idx, idx]
	obj <- glmmkin(disease ~ age + sex, data = pheno, kins = kins, id = "id", family = binomial(link = "logit"), method = "REML", method.optim = "AI")
	expect_equal(obj5$theta, obj$theta)
	expect_equal(obj5$coef, obj$coef)
	obj <- glmmkin(disease ~ age + sex, data = pheno, kins = NULL, id = "id", family = binomial(link = "logit"), method = "REML", method.optim = "AI")
	expect_equal(obj6$theta, obj$theta)
	expect_equal(obj6$coef, obj$coef)
})

test_that("cross-sectional id gt 400 binomial", {
	skip_on_cran()

	data(example)
	suppressWarnings(RNGversion("3.5.0"))
	set.seed(123)
	pheno <- rbind(example$pheno, example$pheno[1:100, ])
	pheno$id <- 1:500
	pheno$disease[sample(1:500,20)] <- NA
	pheno$age[sample(1:500,20)] <- NA
	pheno$sex[sample(1:500,20)] <- NA
	pheno <- pheno[sample(1:500,450), ]
	kins <- diag(500)
	kins[1:400, 1:400] <- example$GRM
	rownames(kins) <- colnames(kins) <- 1:500
	obj5 <- glmmkin(disease ~ age + sex, data = pheno, kins = kins, id = "id", family = binomial(link = "logit"), method = "REML", method.optim = "AI")
	expect_equal(signif(as.numeric(obj5$theta)), signif(c(1, 0.1263611)))
	expect_equal(signif(as.numeric(obj5$coef)), signif(c(0.92300945, -0.01457307, -0.18165858)))
	obj6 <- glmmkin(disease ~ age + sex, data = pheno, kins = NULL, id = "id", family = binomial(link = "logit"), method = "REML", method.optim = "AI")
	expect_equal(c(obj6$theta), 1)
	expect_equal(signif(as.numeric(obj6$coef)), signif(c(0.86840325, -0.01402939, -0.17898775)))
	obj <- glm(disease ~ age + sex, data = pheno, family = binomial(link = "logit"))
	expect_equal(obj6$coef, obj$coef)

	idx <- sample(nrow(pheno))
	pheno <- pheno[idx, ]
	obj <- glmmkin(disease ~ age + sex, data = pheno, kins = kins, id = "id", family = binomial(link = "logit"), method = "REML", method.optim = "AI")
	expect_equal(obj5$theta, obj$theta)
	expect_equal(obj5$coef, obj$coef)
	obj <- glmmkin(disease ~ age + sex, data = pheno, kins = NULL, id = "id", family = binomial(link = "logit"), method = "REML", method.optim = "AI")
	expect_equal(obj6$theta, obj$theta)
	expect_equal(obj6$coef, obj$coef)
	obj <- glm(disease ~ age + sex, data = pheno, family = binomial(link = "logit"))
	expect_equal(obj6$coef, obj$coef)

	idx <- sample(nrow(kins))
	kins <- kins[idx, idx]
	obj <- glmmkin(disease ~ age + sex, data = pheno, kins = kins, id = "id", family = binomial(link = "logit"), method = "REML", method.optim = "AI")
	expect_equal(obj5$theta, obj$theta)
	expect_equal(obj5$coef, obj$coef)
	obj <- glmmkin(disease ~ age + sex, data = pheno, kins = NULL, id = "id", family = binomial(link = "logit"), method = "REML", method.optim = "AI")
	expect_equal(obj6$theta, obj$theta)
	expect_equal(obj6$coef, obj$coef)
})

test_that("cross-sectional id le 400 gaussian", {
	skip_on_cran()

	data(example)
	suppressWarnings(RNGversion("3.5.0"))
	set.seed(123)
	pheno <- rbind(example$pheno, example$pheno[1:100, ])
	pheno$id <- 1:500
	pheno$disease[sample(1:500,20)] <- NA
	pheno$age[sample(1:500,20)] <- NA
	pheno$sex[sample(1:500,20)] <- NA
	pheno <- pheno[sample(1:500,450), ]
	pheno <- pheno[pheno$id <= 400, ]
	kins <- example$GRM
	obj5 <- glmmkin(trait ~ age + sex, data = pheno, kins = kins, id = "id", family = gaussian(link = "identity"), method = "REML", method.optim = "AI")
	expect_equal(signif(as.numeric(obj5$theta)), signif(c(0.7682481, 1.3037467)))
	expect_equal(signif(as.numeric(obj5$coef)), signif(c(3.7634933, 0.0346562, 0.3062784)))
	obj6 <- glmmkin(trait ~ age + sex, data = pheno, kins = NULL, id = "id", family = gaussian(link = "identity"), method = "REML", method.optim = "AI")
	expect_equal(signif(as.numeric(obj6$theta)), signif(1.996857))
	expect_equal(signif(as.numeric(obj6$coef)), signif(c(3.89665633, 0.03156906, 0.27860778)))
	obj <- lm(trait ~ age + sex, data = pheno)
	expect_equal(c(obj6$theta), summary(obj)$sigma^2)
	expect_equal(obj6$coef, obj$coef)

	idx <- sample(nrow(pheno))
	pheno <- pheno[idx, ]
	obj <- glmmkin(trait ~ age + sex, data = pheno, kins = kins, id = "id", family = gaussian(link = "identity"), method = "REML", method.optim = "AI")
	expect_equal(obj5$theta, obj$theta)
	expect_equal(obj5$coef, obj$coef)
	obj <- glmmkin(trait ~ age + sex, data = pheno, kins = NULL, id = "id", family = gaussian(link = "identity"), method = "REML", method.optim = "AI")
	expect_equal(obj6$theta, obj$theta)
	expect_equal(obj6$coef, obj$coef)
	obj <- lm(trait ~ age + sex, data = pheno)
	expect_equal(c(obj6$theta), summary(obj)$sigma^2)
	expect_equal(obj6$coef, obj$coef)

	idx <- sample(nrow(kins))
	kins <- kins[idx, idx]
	obj <- glmmkin(trait ~ age + sex, data = pheno, kins = kins, id = "id", family = gaussian(link = "identity"), method = "REML", method.optim = "AI")
	expect_equal(obj5$theta, obj$theta)
	expect_equal(obj5$coef, obj$coef)
	obj <- glmmkin(trait ~ age + sex, data = pheno, kins = NULL, id = "id", family = gaussian(link = "identity"), method = "REML", method.optim = "AI")
	expect_equal(obj6$theta, obj$theta)
	expect_equal(obj6$coef, obj$coef)
})

test_that("cross-sectional id gt 400 gaussian", {
	skip_on_cran()

	data(example)
	suppressWarnings(RNGversion("3.5.0"))
	set.seed(123)
	pheno <- rbind(example$pheno, example$pheno[1:100, ])
	pheno$id <- 1:500
	pheno$disease[sample(1:500,20)] <- NA
	pheno$age[sample(1:500,20)] <- NA
	pheno$sex[sample(1:500,20)] <- NA
	pheno <- pheno[sample(1:500,450), ]
	kins <- diag(500)
	kins[1:400, 1:400] <- example$GRM
	rownames(kins) <- colnames(kins) <- 1:500
	obj5 <- glmmkin(trait ~ age + sex, data = pheno, kins = kins, id = "id", family = gaussian(link = "identity"), method = "REML", method.optim = "AI")
	expect_equal(signif(as.numeric(obj5$theta)), signif(c(0.7205609, 1.3795639)))
	expect_equal(signif(as.numeric(obj5$coef)), signif(c(3.90459568, 0.03650594, 0.37958367)))
	obj6 <- glmmkin(trait ~ age + sex, data = pheno, kins = NULL, id = "id", family = gaussian(link = "identity"), method = "REML", method.optim = "AI")
	expect_equal(signif(as.numeric(obj6$theta)), signif(2.051809))
	expect_equal(signif(as.numeric(obj6$coef)), signif(c(3.7836398, 0.0344577, 0.3602852)))
	obj <- lm(trait ~ age + sex, data = pheno)
	expect_equal(c(obj6$theta), summary(obj)$sigma^2)
	expect_equal(obj6$coef, obj$coef)

	idx <- sample(nrow(pheno))
	pheno <- pheno[idx, ]
	obj <- glmmkin(trait ~ age + sex, data = pheno, kins = kins, id = "id", family = gaussian(link = "identity"), method = "REML", method.optim = "AI")
	expect_equal(obj5$theta, obj$theta)
	expect_equal(obj5$coef, obj$coef)
	obj <- glmmkin(trait ~ age + sex, data = pheno, kins = NULL, id = "id", family = gaussian(link = "identity"), method = "REML", method.optim = "AI")
	expect_equal(obj6$theta, obj$theta)
	expect_equal(obj6$coef, obj$coef)
	obj <- lm(trait ~ age + sex, data = pheno)
	expect_equal(c(obj6$theta), summary(obj)$sigma^2)
	expect_equal(obj6$coef, obj$coef)

	idx <- sample(nrow(kins))
	kins <- kins[idx, idx]
	obj <- glmmkin(trait ~ age + sex, data = pheno, kins = kins, id = "id", family = gaussian(link = "identity"), method = "REML", method.optim = "AI")
	expect_equal(obj5$theta, obj$theta)
	expect_equal(obj5$coef, obj$coef)
	obj <- glmmkin(trait ~ age + sex, data = pheno, kins = NULL, id = "id", family = gaussian(link = "identity"), method = "REML", method.optim = "AI")
	expect_equal(obj6$theta, obj$theta)
	expect_equal(obj6$coef, obj$coef)
})

