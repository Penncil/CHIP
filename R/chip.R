#' Subsampling with one-step correction for large-scale GLMs
#'
#' Fits a linear, logistic or Poisson regression by (1) optionally screening
#' covariates with an adaptive lasso on a uniform subsample of size `n1`,
#' (2) fitting a pilot model on a uniform subsample of size `n0`, (3) for
#' optimal sampling, drawing a subsample with optimal probabilities (expected
#' size `n1`) and refitting, and (4) applying a one-step Newton correction
#' that uses the score of the non-sampled observations.
#' The estimate before the correction is returned as `beta_hat`.
#'
#' @param X Numeric covariate matrix (no intercept column).
#' @param y Response vector: numeric for `"gaussian"`, 0/1 for `"binomial"`,
#'   non-negative counts for `"poisson"`.
#' @param family One of `"binomial"`, `"gaussian"`, `"poisson"`.
#' @param penalty Logical. `TRUE` (default, high-dimensional) screens covariates
#'   with an adaptive lasso first; `FALSE` (low-dimensional) keeps all
#'   covariates. Guidance: use `FALSE` when the number of covariates is small
#'   relative to the subsample (roughly `ncol(x) < n0 / 10`) and you want
#'   estimates for every covariate; use `TRUE` when there are many covariates
#'   or the true model is believed to be sparse.
#' @param weight_method `"uniform"`, `"optimal"`, or `"auto"` (optimal if
#'   `n1 < sqrt(n)`, otherwise uniform).
#' @param n1 Subsample budget for screening and for the optimal subsample. A
#'   value in (0, 1) is a proportion of `n`; otherwise a count. Default
#'   `ceiling(sqrt(n) * log(p))`, with `p` counting the intercept.
#' @param n0 Initial (pilot) uniform subsample size, as a proportion or a
#'   count. Default `n1`.
#' @param hessian_x Logical. If `TRUE`, optimal probabilities use
#'   \eqn{|y_i - \hat\pi_i| \, \|\hat M^{-1} x_i\|}{|y - pi| ||M^-1 x||}
#'   (includes the Hessian \eqn{\hat M}); if `FALSE` (default),
#'   \eqn{|y_i - \hat\pi_i| \, \|x_i\|}{|y - pi| ||x||}.
#'   Ignored for uniform sampling.
#'
#' @return An object of class `"chip"` with
#' \describe{
#'   \item{beta}{Named coefficients, intercept first; 0 for screened-out
#'     covariates.}
#'   \item{se}{Standard errors; `NA` for screened-out covariates.}
#'   \item{vcov}{Covariance matrix of the selected coefficients.}
#'   \item{s_hat}{Names of the selected covariates.}
#'   \item{pi_opt}{Each individual's probability of entering the final
#'     subsample.}
#'   \item{weights}{Each individual's inverse-probability weight: `1/pi_opt`
#'     if in the final subsample, 0 otherwise.}
#'   \item{delta}{0/1 indicator of being in the final subsample.}
#'   \item{beta_hat}{Subsample estimate before the one-step correction.}
#'   \item{family, weight_method, settings, call}{Fit details.}
#' }
#'
#' @details Standard errors are sandwich estimates computed in one pass over
#' the full data at the final estimate, plus the first-order sampling variance
#' of the one-step score. They assume `n1` is large enough that the one-step
#' estimator is close to the full-data estimator; with `n1` near
#' \eqn{\sqrt{n}} they can be too small. With `penalty = TRUE` they are
#' conditional on the selected covariates.
#'
#' @examples
#' set.seed(123)
#' dat <- data.generate(cov.type = 2, n = 1e5, p = 200, supp = 10, corr = 0.3)
#' fit <- chip(dat$X, dat$Y, family = "binomial", penalty = TRUE)
#' summary(fit)
#' @export
chip <- function(X, y,
                 family = c("binomial", "gaussian", "poisson"),
                 penalty = TRUE,
                 weight_method = c("auto", "optimal", "uniform"),
                 n1 = NULL, n0 = NULL,
                 hessian_x = FALSE) {
  family <- match.arg(family)
  weight_method <- match.arg(weight_method)
  X <- as.matrix(X)
  y <- as.numeric(y)
  n <- length(y)
  p <- ncol(X) + 1
  if (nrow(X) != n) stop("`x` and `y` have different numbers of rows.")
  if (anyNA(X) || anyNA(y)) stop("Missing values are not allowed.")
  if (is.null(colnames(X))) colnames(X) <- paste0("X", seq_len(p - 1))

  n1 <- as_count(if (is.null(n1)) ceiling(sqrt(n) * log(p)) else n1, n, "n1")
  n0 <- as_count(if (is.null(n0)) n1 else n0, n, "n0")
  if (weight_method == "auto") weight_method <- if (n1 < sqrt(n)) "optimal" else "uniform"
  fam <- chip_family(family)

  ##### Screening with adaptive lasso (high-dimensional only) #####
  if (penalty) {
    delta <- stats::rbinom(n, 1, n1 / n)
    J1_idx <- which(delta == 1)
    s_hat <- screen_adaptive_lasso(X[J1_idx, , drop = FALSE], y[J1_idx], family)
  } else {
    s_hat <- seq_len(p - 1)
  }
  if (!penalty && n0 < 10 * length(s_hat))
    warning("`n0` is small relative to the number of covariates; ",
            "consider a larger `n0` or `penalty = TRUE`.")
  X_design <- cbind(`(Intercept)` = 1, X[, s_hat, drop = FALSE])

  ##### Pilot estimator with uniform sampling #####
  pi_opt <- rep(n0 / n, n)
  delta <- stats::rbinom(n, 1, pi_opt)
  beta_hat_1 <- fit_ipw(X_design, y, delta, pi_opt, fam)
  beta_hat <- beta_hat_1

  ##### Optimal sampling #####
  if (weight_method == "optimal") {
    pi_mle <- fam$linkinv(drop(X_design %*% beta_hat_1))
    residuals <- abs(y - pi_mle)
    if (hessian_x) {
      W <- fam$variance(pi_mle)
      M_hat <- crossprod(X_design, X_design * W) / n
      X_tilde <- X_design %*% solve(M_hat)
      x_norms <- sqrt(rowSums(X_tilde^2))
    } else {
      x_norms <- sqrt(rowSums(X_design^2))
    }
    w_num <- residuals * x_norms
    pi_opt <- compute_optimal_pi(u_tilde_norm = w_num, n1 = n1)
    delta <- stats::rbinom(n, 1, pi_opt)
    beta_hat_2 <- fit_ipw(X_design, y, delta, pi_opt, fam)
    beta_hat <- beta_hat_2
  }

  ##### One-step correction and standard errors #####
  beta_final <- one_step(beta_hat, X_design, y, delta, pi_opt, n1, fam)
  V <- sandwich_vcov(beta_final, X_design, y, pi_opt, n1, fam)

  beta <- stats::setNames(numeric(p), c("(Intercept)", colnames(X)))
  se <- stats::setNames(rep(NA_real_, p), names(beta))
  beta[colnames(X_design)] <- beta_final
  se[colnames(X_design)] <- sqrt(diag(V))
  dimnames(V) <- list(colnames(X_design), colnames(X_design))

  structure(list(
    beta = beta,
    se = se,
    vcov = V,
    s_hat = colnames(X)[s_hat],
    pi_opt = pi_opt,
    weights = ifelse(delta == 1, 1 / pi_opt, 0),
    delta = delta,
    beta_hat = stats::setNames(beta_hat, colnames(X_design)),
    family = family,
    weight_method = weight_method,
    settings = list(n = n, p = p, n1 = n1, n0 = n0, penalty = penalty,
                    hessian_x = hessian_x),
    call = match.call()
  ), class = "chip")
}

as_count <- function(v, n, name) {
  if (!is.numeric(v) || length(v) != 1 || v <= 0) stop("`", name, "` must be a positive number.")
  v <- if (v < 1) ceiling(v * n) else ceiling(v)
  if (v > n) stop("`", name, "` exceeds the number of observations.")
  v
}
