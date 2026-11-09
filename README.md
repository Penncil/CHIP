# CHIP: Efficient subsampling with one-step correction for large-scale data

CHIP provides a computationally efficient approach to fitting linear, logistic, and Poisson regression models on large datasets. It first fits a model using a small subsample, then applies a single one-step correction to incorporate information from the full dataset.

For high-dimensional settings, CHIP performs variable screening using the adaptive Lasso on a uniformly sampled subset before model estimation.


## Installation

```r
# install.packages("remotes")
remotes::install_github("penncil/CHIP")
```

## Examples

High-dimensional (199 covariates, 10 signals):

```r
library(CHIP)
set.seed(123)

dat <- data.generate(cov.type = 2, n = 1e5, p = 200, supp = 10, corr = 0.3)

fit <- chip(dat$X, dat$Y, family = "binomial",
            penalty = TRUE, weight_method = "optimal", hessian_x = FALSE)
summary(fit)   # estimates, SEs, z and p-values for selected covariates
fit$beta       # all coefficients (0 for screened-out covariates)
fit$se         # standard errors (NA for screened-out covariates)
fit$pi_opt     # each individual's sampling probability
fit$weights    # each individual's inverse-probability weight (0 if not sampled)
fit$beta_hat   # subsample estimate before the one-step correction
```

Low-dimensional (no screening):

```r
library(CHIP)
set.seed(123)

n <- 1e6
X <- matrix(rnorm(n * 5), nrow = n)
beta_true <- runif(5, -1, 1)
y <- rbinom(n, 1, 1 / (1 + exp(-(X %*% beta_true - 0.85))))

fit <- chip(X, y, family = "binomial", penalty = FALSE,
            weight_method = "optimal", n1 = 1000)
summary(fit)
```

## Options

| Argument        | Meaning                                                             | Default |
|-----------------|---------------------------------------------------------------------|---------|
| `family`        | `"binomial"`, `"gaussian"`, `"poisson"`                             | `"binomial"` |
| `penalty`       | Adaptive-lasso screening (high-D) or not (low-D)                    | `TRUE` |
| `weight_method` | `"uniform"`, `"optimal"`, or `"auto"` (optimal if `n1 < sqrt(n)`, else uniform) | `"auto"` |
| `n1`            | Subsample budget (screening and optimal subsample); proportion or count | `sqrt(n) * log(p)` |
| `n0`            | Initial (pilot) uniform subsample size; proportion or count         | `n1` |
| `hessian_x`     | Include the Hessian in the optimal sampling probabilities           | `FALSE` |
