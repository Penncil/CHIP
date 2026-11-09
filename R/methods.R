#' @export
print.chip <- function(x, ...) {
  st <- x$settings
  header(x)
  print(x$beta[c("(Intercept)", x$s_hat)], ...)
  invisible(x)
}

#' @export
coef.chip <- function(object, ...) object$beta

#' @export
vcov.chip <- function(object, ...) object$vcov

#' @export
summary.chip <- function(object, ...) {
  keep <- c("(Intercept)", object$s_hat)
  est <- object$beta[keep]
  se <- object$se[keep]
  z <- est / se
  object$coef_table <- cbind(Estimate = est, `Std. Error` = se, `z value` = z,
                             `Pr(>|z|)` = 2 * stats::pnorm(-abs(z)))
  class(object) <- c("summary.chip", class(object))
  object
}

#' @export
print.summary.chip <- function(x, ...) {
  st <- x$settings
  header(x)
  stats::printCoefmat(x$coef_table, ...)
  if (st$penalty) cat("\nSEs are conditional on the selected covariates.\n")
  invisible(x)
}

#' Predictions from a CHIP fit
#'
#' @param object A `"chip"` fit.
#' @param newx Covariate matrix with the same columns as `X`.
#' @param type `"link"` for the linear predictor, `"response"` for the mean.
#' @param ... Unused.
#' @export
predict.chip <- function(object, newx, type = c("link", "response"), ...) {
  type <- match.arg(type)
  eta <- drop(cbind(1, as.matrix(newx)) %*% object$beta)
  if (type == "link") eta else chip_family(object$family)$linkinv(eta)
}

header <- function(x) {
  st <- x$settings
  cat("CHIP fit:", x$family, "regression,", x$weight_method, "sampling + one-step\n")
  cat(sprintf("n = %d, p = %d (incl. intercept), selected = %d, n0 = %d, n1 = %d\n\n",
              st$n, st$p, length(x$s_hat), st$n0, st$n1))
}
