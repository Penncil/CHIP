sim_glm <- function(n, family) {
  X <- matrix(rnorm(n * 5), n, 5)
  eta <- -0.85 + drop(X %*% c(0.5, -0.5, 0.3, 0, 0.2))
  y <- switch(family,
              binomial = rbinom(n, 1, plogis(eta)),
              poisson  = rpois(n, exp(eta)),
              gaussian = eta + rnorm(n))
  list(X = X, y = y)
}

for (fam in c("binomial", "gaussian", "poisson")) {
  test_that(paste("low-dimensional", fam, "fit is close to the full-data MLE"), {
    set.seed(1)
    d <- sim_glm(2e4, fam)
    mle <- coef(glm(d$y ~ d$X, family = fam))
    for (wm in c("uniform", "optimal")) {
      fit <- chip(d$X, d$y, family = fam, penalty = FALSE, weight_method = wm, n1 = 0.1)
      expect_s3_class(fit, "chip")
      expect_equal(unname(coef(fit)), unname(mle), tolerance = 0.05)
      expect_true(all(is.finite(fit$se)))
    }
  })
}

test_that("weights and probabilities are consistent", {
  set.seed(2)
  d <- sim_glm(1e4, "binomial")
  fit <- chip(d$X, d$y, penalty = FALSE, weight_method = "optimal", n1 = 500)
  expect_length(fit$pi_opt, 1e4)
  expect_true(all(fit$pi_opt > 0 & fit$pi_opt <= 1))
  expect_equal(sum(fit$pi_opt), 500, tolerance = 1e-6)
  expect_equal(fit$weights[fit$delta == 1], 1 / fit$pi_opt[fit$delta == 1])
  expect_true(all(fit$weights[fit$delta == 0] == 0))
})

test_that("auto chooses optimal when n1 < sqrt(n), uniform otherwise", {
  set.seed(3)
  d <- sim_glm(1e4, "binomial")
  expect_equal(chip(d$X, d$y, penalty = FALSE, n1 = 90)$weight_method, "optimal")
  expect_equal(chip(d$X, d$y, penalty = FALSE, n1 = 500)$weight_method, "uniform")
})

test_that("high-dimensional screening recovers the true support", {
  set.seed(4)
  d <- data.generate(cov.type = 2, n = 5e4, p = 101, supp = 6)
  fit <- chip(d$X, d$Y, penalty = TRUE, n1 = 0.05)
  expect_true(all(paste0("X", 1:6) %in% fit$s_hat))
  expect_true(all(is.na(fit$se[!names(fit$se) %in% c("(Intercept)", fit$s_hat)])))
})

test_that("subsample estimate is returned and predict works", {
  set.seed(5)
  d <- sim_glm(1e4, "poisson")
  fit <- chip(d$X, d$y, family = "poisson", penalty = FALSE, n1 = 0.1)
  expect_length(fit$beta_hat, 6)
  mu <- predict(fit, d$X[1:5, ], type = "response")
  expect_length(mu, 5)
  expect_true(all(mu > 0))
  expect_s3_class(summary(fit), "summary.chip")
})

test_that("bad inputs are rejected", {
  d <- sim_glm(100, "binomial")
  expect_error(chip(d$X, d$y[-1]), "different numbers")
  expect_error(chip(d$X, d$y, n1 = 1000), "exceeds")
  expect_error(chip(d$X, d$y, family = "gamma"))
})
