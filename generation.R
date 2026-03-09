# ==============================================================================
# generation.R - Simulation data generation for PLS-SEM Monte Carlo studies
# ==============================================================================

# --- Formative-Formative (given inter-composite correlations) ---

gscmcovff1 <- function(B, indicatorx, indicatory, wx, wy, Sxixi, R2 = NULL) {
  q  <- nrow(B)
  equals0 <- function(x) all(x == 0)
  q1 <- sum(apply(B, 1, equals0))
  q2 <- q - q1
  p1 <- length(indicatorx)
  p2 <- length(indicatory)

  # Exogenous indicator covariance
  Sxx <- matrix(rep(c(p1:0), p1), p1, p1) * 0.01
  Sxx[upper.tri(Sxx)] <- 0
  Sxx <- (Sxx + t(Sxx)) / p1
  diag(Sxx) <- 1

  # W1: exogenous weight matrix
  W1 <- 1 * (matrix(indicatorx, p1, q1) == matrix(1:q1, p1, q1, byrow = TRUE))
  P <- c(0, cumsum(colSums(W1)))
  P1 <- P
  W1 <- W1 * wx
  for (i in 1:q1) {
    idx <- (P[i] + 1):P[i + 1]
    f <- t(W1[idx, i, drop = FALSE]) %*% Sxx[idx, idx] %*% W1[idx, i, drop = FALSE]
    W1[idx, i] <- W1[idx, i] / sqrt(as.numeric(f))
  }

  if (q1 > 1) {
    for (j in 2:q1) {
      for (i in 1:(j - 1)) {
        idx_i <- (P[i] + 1):P[i + 1]
        idx_j <- (P[j] + 1):P[j + 1]
        subSxx <- Sxx[idx_i, idx_j]
        f <- t(W1[idx_i, i, drop = FALSE]) %*% subSxx %*% W1[idx_j, j, drop = FALSE]
        f <- as.numeric(Sxixi[i, j]) / as.numeric(f)
        Sxx[idx_i, idx_j] <- subSxx * f
        Sxx[idx_j, idx_i] <- t(subSxx * f)
      }
    }
  }

  # Composite covariance matrix and path scaling
  Sll <- diag(q)
  Sll[1:q1, 1:q1] <- Sxixi

  for (m in (q1 + 1):(q1 + q2)) {
    tau <- sqrt(R2[m - q1] / (B[m, 1:(m - 1), drop = FALSE] %*%
                                 Sll[1:(m - 1), 1:(m - 1)] %*% t(B[m, 1:(m - 1), drop = FALSE])))
    B[m, ] <- tau * B[m, ]
    for (j in 1:(m - 1)) {
      Sll[j, m] <- Sll[m, j] <- B[m, 1:(m - 1), drop = FALSE] %*% Sll[1:(m - 1), j, drop = FALSE]
    }
  }

  Setaeta <- Sll[(q1 + 1):(q1 + q2), (q1 + 1):(q1 + q2)]

  # Endogenous indicator covariance
  Syy <- matrix(rep(c(p2:0), p2), p2, p2) * 0.01
  Syy[upper.tri(Syy)] <- 0
  Syy <- (Syy + t(Syy)) / p2
  diag(Syy) <- 1

  # W2: endogenous weight matrix
  W2 <- 1 * (matrix(indicatory, p2, q2) == matrix(1:q2, p2, q2, byrow = TRUE))
  P <- c(0, cumsum(colSums(W2)))
  P2 <- P
  W2 <- W2 * wy
  for (i in 1:q2) {
    idx <- (P[i] + 1):P[i + 1]
    f <- t(W2[idx, i, drop = FALSE]) %*% Syy[idx, idx] %*% W2[idx, i, drop = FALSE]
    W2[idx, i] <- W2[idx, i] / sqrt(as.numeric(f))
  }

  if (q2 > 1) {
    for (j in 2:q2) {
      for (i in 1:(j - 1)) {
        idx_i <- (P[i] + 1):P[i + 1]
        idx_j <- (P[j] + 1):P[j + 1]
        subSyy <- Syy[idx_i, idx_j]
        f <- t(W2[idx_i, i, drop = FALSE]) %*% subSyy %*% W2[idx_j, j, drop = FALSE]
        f <- as.numeric(Setaeta[i, j]) / as.numeric(f)
        Syy[idx_i, idx_j] <- subSyy * f
        Syy[idx_j, idx_i] <- t(subSyy * f)
      }
    }
  }

  Br1 <- B[(q1 + 1):q, 1:q1, drop = FALSE]
  Br2 <- B[(q1 + 1):q, (q1 + 1):q, drop = FALSE]
  IB  <- solve(diag(q2) - t(Br2))

  Sxy <- Sxx %*% W1 %*% t(Br1) %*% IB %*% solve(Setaeta) %*% t(W2) %*% Syy
  S   <- rbind(cbind(Sxx, Sxy), cbind(t(Sxy), Syy))

  list(S = S, B = B, Scomp = Sll, wx = rowSums(W1), wy = rowSums(W2))
}

# --- Reflective-Reflective ---

gscmcovrr1 <- function(B, indicatorx, indicatory, lambdax, lambday, Sxixi, R2 = NULL) {
  q  <- nrow(B)
  equals0 <- function(x) all(x == 0)
  q1 <- sum(apply(B, 1, equals0))
  q2 <- q - q1
  p1 <- length(indicatorx)
  p2 <- length(indicatory)

  B1 <- B[(q1 + 1):q, 1:q1, drop = FALSE]
  B2 <- B[(q1 + 1):q, (q1 + 1):q, drop = FALSE]

  Lx <- 1 * (matrix(indicatorx, p1, q1) == matrix(1:q1, p1, q1, byrow = TRUE)) * lambdax
  Ly <- 1 * (matrix(indicatory, p2, q2) == matrix(1:q2, p2, q2, byrow = TRUE)) * lambday

  Sxx <- Lx %*% Sxixi %*% t(Lx)
  Sdd <- diag(1 - diag(Sxx))
  diag(Sxx) <- 1

  Sll <- diag(q)
  Sll[1:q1, 1:q1] <- Sxixi

  for (m in (q1 + 1):(q1 + q2)) {
    tau <- sqrt(R2[m - q1] / (B[m, 1:(m - 1), drop = FALSE] %*%
                                 Sll[1:(m - 1), 1:(m - 1)] %*% t(B[m, 1:(m - 1), drop = FALSE])))
    B[m, ] <- tau * B[m, ]
    for (j in 1:(m - 1)) {
      Sll[j, m] <- Sll[m, j] <- B[m, 1:(m - 1), drop = FALSE] %*% Sll[1:(m - 1), j, drop = FALSE]
    }
  }
  Setaeta <- Sll[(q1 + 1):(q1 + q2), (q1 + 1):(q1 + q2)]

  Br1 <- B[(q1 + 1):q, 1:q1, drop = FALSE]
  Br2 <- B[(q1 + 1):q, (q1 + 1):q, drop = FALSE]
  IB  <- solve(diag(q2) - t(Br2))

  Syy <- Ly %*% Setaeta %*% t(Ly)
  See <- diag(1 - diag(Syy))
  diag(Syy) <- 1

  Sxy <- Lx %*% Sxixi %*% t(Br1) %*% IB %*% t(Ly)
  S   <- rbind(cbind(Sxx, Sxy), cbind(t(Sxy), Syy))

  list(S = S, B = B, Scomp = Sll, Sdd = Sdd, See = See)
}

# --- Formative-Formative (derive path coefficients from MV covariance) ---

gscmcovff2 <- function(indicatorx, indicatory, wx, wy) {
  q1 <- indicatorx[length(indicatorx)]
  q2 <- indicatory[length(indicatory)]
  q  <- q1 + q2
  p1 <- length(indicatorx)
  p2 <- length(indicatory)
  p  <- p1 + p2

  Smv <- matrix(rep(c(p:0), p), p, p) * 0.01
  Smv[upper.tri(Smv)] <- 0
  Smv <- (Smv + t(Smv)) / p
  diag(Smv) <- 1
  Sxx <- Smv[1:p1, 1:p1]
  Syy <- Smv[(p1 + 1):p, (p1 + 1):p]

  W1 <- 1 * (matrix(indicatorx, p1, q1) == matrix(1:q1, p1, q1, byrow = TRUE))
  P <- c(0, cumsum(colSums(W1)))
  W1 <- W1 * wx
  for (i in 1:q1) {
    idx <- (P[i] + 1):P[i + 1]
    f <- t(W1[idx, i, drop = FALSE]) %*% Sxx[idx, idx] %*% W1[idx, i, drop = FALSE]
    W1[idx, i] <- W1[idx, i] / sqrt(as.numeric(f))
  }

  W2 <- 1 * (matrix(indicatory, p2, q2) == matrix(1:q2, p2, q2, byrow = TRUE))
  P <- c(0, cumsum(colSums(W2)))
  W2 <- W2 * wy
  for (i in 1:q2) {
    idx <- (P[i] + 1):P[i + 1]
    f <- t(W2[idx, i, drop = FALSE]) %*% Syy[idx, idx] %*% W2[idx, i, drop = FALSE]
    W2[idx, i] <- W2[idx, i] / sqrt(as.numeric(f))
  }

  W   <- rbind(cbind(W1, matrix(0, p1, q2)), cbind(matrix(0, p2, q1), W2))
  Sll <- t(W) %*% Smv %*% W

  B <- matrix(0, q, q)
  for (i in 1:q2) {
    B[q1 + i, 1:(q1 + i - 1)] <- t(solve(Sll[1:(q1 + i - 1), 1:(q1 + i - 1)]) %*%
                                       Sll[1:(q1 + i - 1), q1 + i])
  }

  list(S = Smv, B = B, Scomp = Sll, wx = rowSums(W1), wy = rowSums(W2))
}

# --- Formative-Formative (given path coefficients, derive MV covariance) ---

gscmcovff3 <- function(B, indicatorx, indicatory, wx, wy, R2 = NULL) {
  q  <- nrow(B)
  equals0 <- function(x) all(x == 0)
  q1 <- sum(apply(B, 1, equals0))
  q2 <- q - q1
  p1 <- length(indicatorx)
  p2 <- length(indicatory)

  Sxx <- matrix(rep(c(p1:0), p1), p1, p1) * 0.01
  Sxx[upper.tri(Sxx)] <- 0
  Sxx <- (Sxx + t(Sxx)) / p1
  diag(Sxx) <- 1

  W1 <- 1 * (matrix(indicatorx, p1, q1) == matrix(1:q1, p1, q1, byrow = TRUE))
  P <- c(0, cumsum(colSums(W1)))
  W1 <- W1 * wx
  for (i in 1:q1) {
    idx <- (P[i] + 1):P[i + 1]
    f <- t(W1[idx, i, drop = FALSE]) %*% Sxx[idx, idx] %*% W1[idx, i, drop = FALSE]
    W1[idx, i] <- W1[idx, i] / sqrt(as.numeric(f))
  }

  Sxixi <- t(W1) %*% Sxx %*% W1

  Sll <- diag(q)
  Sll[1:q1, 1:q1] <- Sxixi

  for (m in (q1 + 1):(q1 + q2)) {
    tau <- sqrt(R2[m - q1] / (B[m, 1:(m - 1), drop = FALSE] %*%
                                 Sll[1:(m - 1), 1:(m - 1)] %*% t(B[m, 1:(m - 1), drop = FALSE])))
    B[m, ] <- tau * B[m, ]
    for (j in 1:(m - 1)) {
      Sll[j, m] <- Sll[m, j] <- B[m, 1:(m - 1), drop = FALSE] %*% Sll[1:(m - 1), j, drop = FALSE]
    }
  }

  Setaeta <- Sll[(q1 + 1):(q1 + q2), (q1 + 1):(q1 + q2)]

  Syy <- matrix(rep(c(p2:0), p2), p2, p2) * 0.01
  Syy[upper.tri(Syy)] <- 0
  Syy <- (Syy + t(Syy)) / p2
  diag(Syy) <- 1

  W2 <- 1 * (matrix(indicatory, p2, q2) == matrix(1:q2, p2, q2, byrow = TRUE))
  P <- c(0, cumsum(colSums(W2)))
  W2 <- W2 * wy
  for (i in 1:q2) {
    idx <- (P[i] + 1):P[i + 1]
    f <- t(W2[idx, i, drop = FALSE]) %*% Syy[idx, idx] %*% W2[idx, i, drop = FALSE]
    W2[idx, i] <- W2[idx, i] / sqrt(as.numeric(f))
  }

  if (q2 > 1) {
    for (j in 2:q2) {
      for (i in 1:(j - 1)) {
        idx_i <- (P[i] + 1):P[i + 1]
        idx_j <- (P[j] + 1):P[j + 1]
        subSyy <- Syy[idx_i, idx_j]
        f <- t(W2[idx_i, i, drop = FALSE]) %*% subSyy %*% W2[idx_j, j, drop = FALSE]
        f <- as.numeric(Setaeta[i, j]) / as.numeric(f)
        Syy[idx_i, idx_j] <- subSyy * f
        Syy[idx_j, idx_i] <- t(subSyy * f)
      }
    }
  }

  Br1 <- B[(q1 + 1):q, 1:q1, drop = FALSE]
  Br2 <- B[(q1 + 1):q, (q1 + 1):q, drop = FALSE]
  IB  <- solve(diag(q2) - t(Br2))

  Sxy <- Sxx %*% W1 %*% t(Br1) %*% IB %*% solve(Setaeta) %*% t(W2) %*% Syy
  S   <- rbind(cbind(Sxx, Sxy), cbind(t(Sxy), Syy))

  list(S = S, B = B, Scomp = Sll, wx = rowSums(W1), wy = rowSums(W2))
}
