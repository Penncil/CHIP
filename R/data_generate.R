#' Simulate logistic-regression data
#'
#' Data generator used in the CHIP simulation study.
#'
#' @param cov.type Covariance of X: 1 = identity, 2 = AR(1) `corr^|j-k|`,
#'   3 = exchangeable `corr`.
#' @param n Number of observations.
#' @param p Number of coefficients including the intercept (X has `p - 1`
#'   columns).
#' @param supp Number of strong signals (alternating `0.5, -0.5`).
#' @param gamma Intercept.
#' @param corr Correlation parameter for `cov.type` 2 and 3.
#' @param weak_s Magnitude of the remaining (weak) coefficients.
#' @param beta Optional coefficient vector (length `p - 1`) overriding the
#'   default.
#' @return List with `Y`, `X` and `beta.true` (intercept first).
#' @export
data.generate <- function(cov.type = 1, n = 200, p = 50, supp = 5, gamma = -1,
                          corr = 0.3, weak_s = 0, beta = NULL) {
  if (is.null(beta)) {
    beta <- c(rep(c(0.5, -0.5), length.out = supp),
              rep(c(weak_s, -weak_s), length.out = p - supp - 1))
  }
  if (cov.type == 1) {
    sigmax <- diag(rep(1, p - 1))
  } else if (cov.type == 2) {
    sigmax <- outer(1:(p - 1), 1:(p - 1), function(x, y) corr^(abs(x - y)))
  } else if (cov.type == 3) {
    sigmax <- matrix(corr, p - 1, p - 1) + diag(1 - corr, p - 1)
  }
  X <- mvtnorm::rmvnorm(n, rep(0, p - 1), sigmax)
  meanTreat <- 1 / (1 + exp(-(gamma + X %*% beta)))
  Y <- stats::rbinom(n, 1, meanTreat)
  list(Y = Y, X = X, beta.true = c(gamma, beta))
}
