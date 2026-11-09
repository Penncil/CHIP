# Internal estimation steps. `X_design` is the design matrix with an intercept
# column, restricted to the selected covariates.

# Quasi-families give the same estimates as binomial/poisson but accept the
# non-integer inverse-probability weights without warnings.
chip_family <- function(family) {
  switch(family,
         binomial = stats::quasibinomial(),
         poisson  = stats::quasipoisson(),
         gaussian = stats::gaussian())
}

screen_adaptive_lasso <- function(X1, y1, family, gamma_power = 1) {
  cv1 <- glmnet::cv.glmnet(x = X1, y = y1, family = family, alpha = 1)
  beta_pilot <- as.numeric(stats::coef(cv1, s = "lambda.min"))
  eps_small <- 1e-4
  omega <- 1 / ((abs(beta_pilot[-1]) + eps_small)^gamma_power)
  cv_adapt <- glmnet::cv.glmnet(x = X1, y = y1, family = family, alpha = 1,
                                penalty.factor = omega)
  which(as.numeric(stats::coef(cv_adapt, s = "lambda.1se"))[-1] != 0)
}

# Inverse-probability-weighted GLM on J1 = {delta == 1}. Weights are rescaled
# to max 1 (estimates unchanged; avoids divergent glm starting values).
fit_ipw <- function(X_design, y, delta, pi_opt, fam) {
  J1_idx <- which(delta == 1)
  w_obs_J1 <- 1 / pi_opt[J1_idx]
  w_obs_J1 <- w_obs_J1 / max(w_obs_J1)
  fit <- stats::glm.fit(X_design[J1_idx, , drop = FALSE], y[J1_idx],
                        weights = w_obs_J1, family = fam, intercept = TRUE)
  if (anyNA(fit$coefficients))
    stop("Subsample fit has aliased coefficients; increase the subsample size.")
  if (!fit$converged) warning("Subsample fit did not converge.")
  unname(fit$coefficients)
}

compute_optimal_pi <- function(u_tilde_norm, n1) {
  n <- length(u_tilde_norm)
  m <- n1

  # Step 1: sort ||U_i|| in ascending order
  sorted_vals <- sort(u_tilde_norm)

  # Step 2: vectorized search for l
  l_seq <- (n - m + 1):n
  left_side <- m + l_seq - n
  right_side <- cumsum(sorted_vals)[l_seq] / sorted_vals[l_seq]
  valid_l <- which(left_side <= right_side)

  if (length(valid_l) == 0) {
    warning("No valid l found; using proportional sampling instead.")
    p <- n1 * u_tilde_norm / sum(u_tilde_norm)
    p[p > 1] <- 1
    return(p)
  }
  l_max <- max(l_seq[valid_l])

  # Step 3: determine A^k and compute final p_i
  threshold <- sorted_vals[l_max]
  A_idx <- which(u_tilde_norm > threshold)
  multiplier <- (m + l_max - n) / sum(sorted_vals[1:l_max])
  p <- multiplier * u_tilde_norm
  p[A_idx] <- 1
  p
}

# One-step update: Hessian estimated on J1, score on J2.
one_step <- function(beta_hat, X_design, y, delta, pi_opt, n1, fam) {
  n <- length(y)
  J1_idx <- which(delta == 1)
  J2_idx <- which(delta == 0)
  X1_design <- X_design[J1_idx, , drop = FALSE]
  X2_design <- X_design[J2_idx, , drop = FALSE]

  p1_hat <- fam$linkinv(drop(X1_design %*% beta_hat))
  p2_hat <- fam$linkinv(drop(X2_design %*% beta_hat))

  # score on J2
  mult_J2 <- -(y[J2_idx] - p2_hat) / (1 - pi_opt[J2_idx])
  score_mean_J2 <- colSums(X2_design * mult_J2) / n

  # Hessian estimate on J1
  W1 <- fam$variance(p1_hat) / pi_opt[J1_idx]
  H_hat <- crossprod(X1_design, X1_design * W1) / n

  delta_active <- solve(H_hat, score_mean_J2)
  beta_hat - ((n - n1) / n) * delta_active
}

# Sandwich covariance at beta: full-data model variance plus the first-order
# sampling variance of the Horvitz-Thompson score on J2.
sandwich_vcov <- function(beta, X_design, y, pi_opt, n1, fam) {
  n <- length(y)
  pi_hat <- fam$linkinv(drop(X_design %*% beta))
  H <- crossprod(X_design, X_design * fam$variance(pi_hat)) / n
  c2 <- ((n - n1) / n)^2
  infl <- ifelse(pi_opt < 1, 1 + c2 * pi_opt / (1 - pi_opt), 1)
  meat <- crossprod(X_design * ((y - pi_hat) * sqrt(infl))) / n
  H_inv <- solve(H)
  H_inv %*% meat %*% H_inv / n
}
