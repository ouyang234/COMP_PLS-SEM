library(mvtnorm)
library(dplyr)
# ilr变换与逆变换函数
# 通过ilr变换生成成分数据
# 定义成分数据的运算
# 成分数据的回归拟合：成分+成分、数值+成分、成分+数值


Closure <- function(x) x / apply(x, 1, sum) # x是matrix不是向量

helm <- function(n) {
  mat <- matrix(0, n - 1, n)
  i <- 2:n
  r <- 1 / sqrt(i * (i - 1))
  for (j in 1:(n - 1)) mat[n - j, 1:c(j + 1)] <- c(rep(r[j], j), -j * r[j])
  return(mat)
}

gmean <- function(x) exp(mean(log(x)))

clr <- function(x) {
  if (dim(x)[2] == 1) {
    x.clr <- x
  } else {
    gm <- apply(x, 1, gmean)
    x.clr <- log(x / gm)
  }
  return(x.clr)
}

clr_i <- function(x) {
  dat <- Closure(exp(x))
  return(dat)
}

ilr <- function(x) {
  D <- dim(x)[2]
  x.ilr <- as.matrix(clr(x)) %*% as.matrix(t(helm(D)))
  return(x.ilr)
}

ilr_i <- function(x) {
  D1 <- dim(x)[2]
  dat <- clr_i(x %*% helm(D1 + 1))
  return(dat)
}

reg <- function(x, y, inter = TRUE) {
  if (inter) x <- cbind(1, x)
  coe <- solve(crossprod(x, x), crossprod(x, y))
  return(coe)
}

# x为数值型变量，y为成分数据
# 返回的系数为（beta_0, beta_1, ..., beta_r)T
# beta_i为成分数据
real_comp_reg <- function(x, y) {
  y <- ilr_i(y)
  beta <- reg(x, y, inter = TRUE)
  comp_coe <- ilr(beta)
  return(comp_coe)
}

# x为成分数据（一元SD），y为数值型数据
# 返回的系数为R*SD，b0为R，b1为成分数据
comp_real_reg <- function(x, y) {
  x <- ilr_i(x)
  beta <- reg(x, y, inter = TRUE)
  comp_coe <- ilr(beta)
  return(comp_coe)
}

# x为成分数据（n*pD)，y为成分数据（同分量）
# 回归系数为数值型数据
# x, y需经过中心化进入函数，将结果再逆变换
# 博士论文中说与ilr变换的结果保持一致，对此有疑问
comp_comp_reg <- function(x, y, comp_dim) {
  D <- comp_dim
  n <- dim(x)[1]
  p <- dim(x)[2] / D
  if (dim(x)[2] %% D != 0) {
    stop("the columns of covariates x should be multiple of D")
  }
  x_split <- split(t(data.frame(x)), list(rep(1:p, each = D)))
  ilr_x <- function(i, x_split) {
    x.i <- t(matrix(x_split[[i]], D, n))
    x.i.ilr <- ilr(x.i)
    return(x.i.ilr)
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
  in_prod_m <- function(x) sum(diag(x))
  xty <- split(data.frame(xty_comp), list(rep(1:p, each = D))) %>%
    lapply(as.matrix) %>%
    lapply(in_prod_m) %>%
    unlist() %>%
    matrix(ncol = 1)
  coe <- solve(xtx, xty)
}

# 成分数据的内积
innerprod_m <- function(x, y) {
  x.ilr <- ilr(x)
  y.ilr <- ilr(y)
  prod <- (x.ilr * y.ilr) %>% apply(1, sum)
  return(prod)
}

# 经计算可验证成分数据的内积等于其ilr变换后数据的内积
innerprodcomp <- function(x, y) {
  D <- length(x)
  x_ <- matrix(0, D, D)
  y_ <- matrix(0, D, D)
  lx <- log(x)
  ly <- log(y)
  for (i in 1:D) {
    x_[i, ] <- lx[i] - lx
    y_[i, ] <- ly[i] - ly
  }
  prod <- sum(x_ * y_) / (2 * D)
  return(prod)
}

Rmse <- function(x.ori, x.est) {
  residual <- sqrt(innerprodcomp(x.ori / x.est))
  return(residual)
}

unit_v <- function(x) {
  mod <- sqrt(sum(x^2))
  return(x / mod)
}

inner_prod_v <- function(x, y) sum(unit_v(x) * unit_v(y))
inner_prod_m <- function(x, y) t(apply(x, 1, unit_v) * apply(y, 1, unit_v)) %>% apply(1, sum)
