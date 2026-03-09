# ==============================================================================
# comp_function.R - Compositional data transformations and regression
# ==============================================================================

library(mvtnorm)
library(dplyr)

# --- Compositional Transformations ---

Closure <- function(x) x / apply(x, 1, sum)

helm <- function(n) {
  mat <- matrix(0, n - 1, n)
  i <- 2:n
  r <- 1 / sqrt(i * (i - 1))
  for (j in 1:(n - 1)) {
    mat[n - j, 1:(j + 1)] <- c(rep(r[j], j), -j * r[j])
  }
  return(mat)
}

gmean <- function(x) exp(mean(log(x)))

clr <- function(x) {
  if (ncol(x) == 1) return(x)
  gm <- apply(x, 1, gmean)
  log(x / gm)
}

clr_i <- function(x) {
  Closure(exp(x))
}

ilr <- function(x) {
  D <- ncol(x)
  as.matrix(clr(x)) %*% t(helm(D))
}

ilr_i <- function(x) {
  D1 <- ncol(x)
  clr_i(x %*% helm(D1 + 1))
}

# --- Regression ---

reg <- function(x, y, inter = TRUE) {
  if (inter) x <- cbind(1, x)
  solve(crossprod(x, x), crossprod(x, y))
}

# x: numeric predictors, y: compositional response
# Returns compositional coefficients (beta_0, beta_1, ..., beta_r)^T
real_comp_reg <- function(x, y) {
  ilr(reg(x, ilr_i(y), inter = TRUE))
}

# x: compositional predictor (univariate simplex), y: numeric response
# Returns coefficients: b0 is real, b1 is compositional
comp_real_reg <- function(x, y) {
  ilr(reg(ilr_i(x), y, inter = TRUE))
}

# x: compositional predictors (n * pD), y: compositional response (same dimension)
# x and y should be centered before input; result needs inverse transformation
comp_comp_reg <- function(x, y, comp_dim) {
  D <- comp_dim
  n <- nrow(x)
  p <- ncol(x) / D
  if (ncol(x) %% D != 0) {
    stop("the columns of covariates x should be a multiple of D")
  }

  x_split <- split(t(data.frame(x)), list(rep(1:p, each = D)))

  ilr_x <- function(i, x_split) {
    x.i <- t(matrix(x_split[[i]], D, n))
    ilr(x.i)
  }
  x.ilr <- do.call(cbind, lapply(matrix(1:p), ilr_x, x_split = x_split))
  y.ilr <- ilr(y)

  xtx_comp <- crossprod(x.ilr, x.ilr)
  xty_comp <- crossprod(x.ilr, y.ilr)

  in_prod_v <- function(x, D) sum(diag(matrix(x, D, D)))

  xtx <- split(data.frame(xtx_comp), list(rep(1:p, each = D))) %>%
    lapply(t) %>%
    lapply(split, list(rep(1:p, each = D))) %>%
    lapply(lapply, in_prod_v, D) %>%
    unlist() %>%
    matrix(p, p) %>%
    t()

  xty <- split(data.frame(xty_comp), list(rep(1:p, each = D))) %>%
    lapply(as.matrix) %>%
    lapply(function(m) sum(diag(m))) %>%
    unlist() %>%
    matrix(ncol = 1)

  solve(xtx, xty)
}

# --- Inner Products for Compositional Data ---

innerprod_m <- function(x, y) {
  x.ilr <- ilr(x)
  y.ilr <- ilr(y)
  rowSums(x.ilr * y.ilr)
}

# Aitchison inner product (verifiable: equals inner product of ilr-transformed data)
innerprodcomp <- function(x, y) {
  D <- length(x)
  lx <- log(x)
  ly <- log(y)
  x_ <- outer(lx, lx, `-`)
  y_ <- outer(ly, ly, `-`)
  sum(x_ * y_) / (2 * D)
}

# Compositional RMSE between original and estimated compositions
Rmse <- function(x.ori, x.est) {
  sqrt(innerprodcomp(x.ori / x.est, x.ori / x.est))
}

# --- Vector Utilities ---

unit_v <- function(x) x / sqrt(sum(x^2))

inner_prod_v <- function(x, y) sum(unit_v(x) * unit_v(y))

inner_prod_m <- function(x, y) {
  rowSums(t(apply(x, 1, unit_v)) * t(apply(y, 1, unit_v)))
}
