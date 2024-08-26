
# B <- matrix(c(rep(0, 3), .6, rep(0, 5), .5, .6, rep(0, 5), .5, rep(0, 6), .4, rep(0, 5), .4, rep(0, 6)), 6, 6)
# Sxixi = matrix(c(1, .4, .1, .4, 1, .3, .1, .3, 1), 3, 3)
# r2 = c(.8, .7, .6)

gscmcovff1 <- function(B,indicatorx,indicatory,wx,wy,Sxixi,R2=NULL){
  q <- nrow(B)                                         # number of latent variables
  b <- rep(1,q)
  equals0 <- function(x) all(x == 0)
  q1 <- length(b[apply(B,1,equals0)])         # number of exogenous composites
  q2 <- q - q1                                          # number of endogenous composites
  p1 <- length(indicatorx)                         # number of indicators of exogenous composites
  p2 <- length(indicatory)                         # number of indicators of endogenous composites
  p <- p1+p2
  
  
  Sxx <- matrix(rep(c(p1:0),p1),p1,p1)*0.01
  Sxx[upper.tri(Sxx)] <- 0
  Sxx <- (Sxx+t(Sxx))/p1
  diag(Sxx) <- 1
  
  W1 <- matrix(indicatorx,p1,q1)                # construction of W1
  W1 <- 1*(W1 == matrix(c(1:q1),p1,q1,byrow=T))
  P <- cumsum(colSums(W1))
  P <- c(0,P)
  P1 <- P
  
  W1 <- W1*wx
  for (i in c(1:q1)){                                      # scaling of weights
    subSxx <- Sxx[(P[i] + 1):P[i+1],(P[i] + 1):P[i+1]]
    f <-  t(W1[(P[i] + 1):P[i+1],i,drop=F])%*%subSxx%*%W1[(P[i] + 1):P[i+1],i,drop=F]
    W1[(P[i] + 1):P[i+1],i] <- W1[(P[i] + 1):P[i+1],i]/sqrt(as.numeric(f))
  }
  if (q1 > 1){
    for (j in c(2:q1)){                                     # scaling of off-diagonal of Sxx
      for (i in c(1:(j-1))){
        subSxx <-  Sxx[(P[i]+1):P[i+1],(P[j]+1):P[j+1]]
        f <-  t(W1[(P[i] + 1):P[i+1],i,drop=F])%*%subSxx%*%W1[(P[j] + 1):P[j+1],j,drop=F]
        f <- as.numeric(Sxixi[i,j])/as.numeric(f)
        Sxx[(P[i]+1):P[i+1],(P[j]+1):P[j+1]] <- subSxx*f
        Sxx[(P[j]+1):P[j+1],(P[i]+1):P[i+1]] <- t(subSxx*f)
      }
    }}
  
  #  computing covariance matrix of composites an scaling of path coefficients
  
  Sll <- diag(q)
  Sll[1:q1,1:q1] <- Sxixi
  
  # if (length(R2) > 0){                                      # scaling of path coefficients
  #   for (m in c((q1+1):(q1+q2))){
  #     tau <- sqrt(R2[m-q1]/(B[m,1:(m-1),drop=F]%*%Sll[1:(m-1),1:(m-1)]%*%t(B[m,1:(m-1),drop=F])))
  #     B[m,] <- tau*B[m,]
  #   }
  # }
  # for (m in c((q1+1):(q1+q2))){
  #   for (j in c(1:(m-1))){
  #     Sll[j,m] <- Sll[m,j] <- B[m,1:(m-1),drop=F]%*%Sll[1:(m-1),j,drop=F]
  #   }
  # }
  for (m in c((q1+1):(q1+q2))){
    tau <- sqrt(R2[m-q1]/(B[m,1:(m-1),drop=F]%*%Sll[1:(m-1),1:(m-1)]%*%t(B[m,1:(m-1),drop=F])))
    B[m,] <- tau*B[m,]
    for (j in c(1:(m-1))){
      Sll[j,m] <- Sll[m,j] <- B[m,1:(m-1),drop=F]%*%Sll[1:(m-1),j,drop=F]
    }
  }
  
  Setaeta <- Sll[(q1+1):(q1+q2),(q1+1):(q1+q2)]
  
  # construction of Syy, W2
  
  Syy <- matrix(rep(c(p2:0),p2),p2,p2)*0.01
  Syy[upper.tri(Syy)] <- 0
  Syy <- (Syy+t(Syy))/p2
  diag(Syy) <- 1
  
  
  W2 <- matrix(indicatory,p2,q2)                # construction of W2
  W2 <- 1*(W2 == matrix(c(1:q2),p2,q2,byrow=T))
  P <- cumsum(colSums(W2))
  P <- c(0,P)
  P2 <- P
  W2 <- W2*wy
  for (i in c(1:q2)){                                      # scaling of weights
    subSyy <- Syy[(P[i] + 1):P[i+1],(P[i] + 1):P[i+1]]
    f <-  t(W2[(P[i] + 1):P[i+1],i,drop=F])%*%subSyy%*%W2[(P[i] + 1):P[i+1],i,drop=F]
    W2[(P[i] + 1):P[i+1],i] <- W2[(P[i] + 1):P[i+1],i]/sqrt(as.numeric(f))
  }
  
  if (q2 > 1){
    for (j in c(2:q2)){                                     # scaling of off-diagonal of Syy
      for (i in c(1:(j-1))){
        subSyy <-  Syy[(P[i]+1):P[i+1],(P[j]+1):P[j+1]]
        f <-  t(W2[(P[i] + 1):P[i+1],i,drop=F])%*%subSyy%*%W2[(P[j] + 1):P[j+1],j,drop=F]
        f <- as.numeric(Setaeta[i,j])/as.numeric(f)
        Syy[(P[i]+1):P[i+1],(P[j]+1):P[j+1]] <- subSyy*f
        Syy[(P[j]+1):P[j+1],(P[i]+1):P[i+1]] <- t(subSyy*f)
      }
    }}
  
  Br1 <- B[(q1+1):q,1:q1,drop=FALSE]
  Br2 <- B[(q1+1):q,(q1+1):q,drop=FALSE]
  IB <- solve(diag(q2) - t(Br2))
  
  Sxy <- Sxx%*%W1%*%t(Br1)%*%IB%*%solve(Setaeta)%*%t(W2)%*%Syy
  S <- rbind(cbind(Sxx,Sxy),cbind(t(Sxy),Syy))
  
  out <- list(S=S,B=B,Scomp=Sll,wx=rowSums(W1),wy=rowSums(W2))
  return(out)
}


gscmcovrr1 <- function(B,indicatorx,indicatory,lambdax,lambday,Sxixi,R2=NULL){
  
  q <- nrow(B)                                     # number of latent variables
  b <- rep(1,q)
  equals0 <- function(x) all(x == 0)
  q1 <- length(b[apply(B,1,equals0)])    # number of exogenous composites
  q2 <- q - q1                                     # number of endogenous composites
  p1 <- length(indicatorx)                      # number of indicators of exogenous composites
  p2 <- length(indicatory)                       # number of indicators of endogenous composites
  p <- p1+p2                                # number of manifest variables
  
  B1 <- B[(q1+1):q,1:q1,drop=FALSE]
  B2 <- B[(q1+1):q,(q1+1):q,drop=FALSE]
  
  # construction of Lx
  Lx <- matrix(indicatorx,p1,q1)
  Lx <- 1*(Lx == matrix(c(1:q1),p1,q1,byrow=T))
  Lx <- Lx*lambdax
  # construction of Ly
  Ly <- matrix(indicatory,p2,q2)
  Ly <- 1*(Ly == matrix(c(1:q2),p2,q2,byrow=T))
  Ly <- Ly*lambday
  
  Sxx <- Lx%*%Sxixi%*%t(Lx)
  Sdd  <-  diag(1 - diag(Sxx))
  diag(Sxx) <- 1
  
  #  computing covariance matrix of composites an scaling of path coefficients
  
  Sll <- diag(q)
  Sll[1:q1,1:q1] <- Sxixi
  
  for (m in c((q1+1):(q1+q2))){
    tau <- sqrt(R2[m-q1]/(B[m,1:(m-1),drop=F]%*%Sll[1:(m-1),1:(m-1)]%*%t(B[m,1:(m-1),drop=F])))
    B[m,] <- tau*B[m,]
    for (j in c(1:(m-1))){
      Sll[j,m] <- Sll[m,j] <- B[m,1:(m-1),drop=F]%*%Sll[1:(m-1),j,drop=F]
    }
  }
  Setaeta <- Sll[(q1+1):(q1+q2),(q1+1):(q1+q2)]
  
  Br1 <- B[(q1+1):q,1:q1,drop=FALSE]
  Br2 <- B[(q1+1):q,(q1+1):q,drop=FALSE]
  IB <- solve(diag(q2) - t(Br2))
  
  Syy <- Ly%*%Setaeta%*%t(Ly)
  See <- diag(1-diag(Syy))
  diag(Syy) <-1
  
  Sxy <- Lx%*%Sxixi%*%t(Br1)%*%IB%*%t(Ly)
  S <- rbind(cbind(Sxx,Sxy),cbind(t(Sxy),Syy))
  
  out <- list(S=S,B=B,Scomp=Sll,Sdd=Sdd,See=See)
  return(out)
}

gscmcovff2 <- function(indicatorx,indicatory,wx,wy){
  # 给定显变量协方差矩阵、外部模型系数、显变量与隐变量的对应关系
  # Smv在需要调整时放在函数外面
  
  q1 <- indicatorx[length(indicatorx)]         # number of exogenous composites
  q2 <- indicatory[length(indicatory)]        # number of endogenous composites
  q <- q1+q2
  p1 <- length(indicatorx)                         # number of indicators of exogenous composites
  p2 <- length(indicatory)                         # number of indicators of endogenous composites
  p <- p1+p2
  
  Smv <- matrix(rep(c(p:0),p),p,p)*0.01
  Smv[upper.tri(Smv)] <- 0
  Smv <- (Smv+t(Smv))/p
  diag(Smv) <- 1
  Sxx <- Smv[1:p1, 1:p1]
  Syy <- Smv[(p1+1):p,(p1+1):p]
  
  W1 <- matrix(indicatorx,p1,q1)                # construction of W1
  W1 <- 1*(W1 == matrix(c(1:q1),p1,q1,byrow=T))
  P <- cumsum(colSums(W1))
  P <- c(0,P)
  P1 <- P
  
  W1 <- W1*wx
  for (i in c(1:q1)){                                      # scaling of weights
    subSxx <- Sxx[(P[i] + 1):P[i+1],(P[i] + 1):P[i+1]]
    f <-  t(W1[(P[i] + 1):P[i+1],i,drop=F])%*%subSxx%*%W1[(P[i] + 1):P[i+1],i,drop=F]
    W1[(P[i] + 1):P[i+1],i] <- W1[(P[i] + 1):P[i+1],i]/sqrt(as.numeric(f))
  }
  
  W2 <- matrix(indicatory,p2,q2)                # construction of W2
  W2 <- 1*(W2 == matrix(c(1:q2),p2,q2,byrow=T))
  P <- cumsum(colSums(W2))
  P <- c(0,P)
  P2 <- P
  W2 <- W2*wy
  for (i in c(1:q2)){                                      # scaling of weights
    subSyy <- Syy[(P[i] + 1):P[i+1],(P[i] + 1):P[i+1]]
    f <-  t(W2[(P[i] + 1):P[i+1],i,drop=F])%*%subSyy%*%W2[(P[i] + 1):P[i+1],i,drop=F]
    W2[(P[i] + 1):P[i+1],i] <- W2[(P[i] + 1):P[i+1],i]/sqrt(as.numeric(f))
  }
  
  W <- rbind(cbind(W1, matrix(0, p1, q2)), cbind(matrix(0, p2, q1), W2))
  Sll <- t(W) %*% Smv %*% W
  
  B <- matrix(0, q, q)
  for (i in 1:q2){
    B[q1+i,1:(q1+i-1)] = t(solve(Sll[1:(q1+i-1),1:(q1+i-1)])%*%Sll[1:(q1+i-1), q1+i])
  }
  
  out <- list(S=Smv,B=B,Scomp=Sll,wx=rowSums(W1),wy=rowSums(W2))
  return(out)
}

gscmcovff3 <- function(B,indicatorx,indicatory,wx,wy,R2=NULL){
  q <- nrow(B)                                         # number of latent variables
  b <- rep(1,q)
  equals0 <- function(x) all(x == 0)
  q1 <- length(b[apply(B,1,equals0)])         # number of exogenous composites
  q2 <- q - q1                                          # number of endogenous composites
  p1 <- length(indicatorx)                         # number of indicators of exogenous composites
  p2 <- length(indicatory)                         # number of indicators of endogenous composites
  p <- p1+p2
  
  
  Sxx <- matrix(rep(c(p1:0),p1),p1,p1)*0.01
  Sxx[upper.tri(Sxx)] <- 0
  Sxx <- (Sxx+t(Sxx))/p1
  diag(Sxx) <- 1
  
  W1 <- matrix(indicatorx,p1,q1)                # construction of W1
  W1 <- 1*(W1 == matrix(c(1:q1),p1,q1,byrow=T))
  P <- cumsum(colSums(W1))
  P <- c(0,P)
  P1 <- P
  
  W1 <- W1*wx
  for (i in c(1:q1)){                                      # scaling of weights
    subSxx <- Sxx[(P[i] + 1):P[i+1],(P[i] + 1):P[i+1]]
    f <-  t(W1[(P[i] + 1):P[i+1],i,drop=F])%*%subSxx%*%W1[(P[i] + 1):P[i+1],i,drop=F]
    W1[(P[i] + 1):P[i+1],i] <- W1[(P[i] + 1):P[i+1],i]/sqrt(as.numeric(f))
  }

  Sxixi <- t(W1) %*% Sxx %*% W1
  
  #  computing covariance matrix of composites an scaling of path coefficients
  
  Sll <- diag(q)
  Sll[1:q1,1:q1] <- Sxixi
  
  # if (length(R2) > 0){                                      # scaling of path coefficients
  #   for (m in c((q1+1):(q1+q2))){
  #     tau <- sqrt(R2[m-q1]/(B[m,1:(m-1),drop=F]%*%Sll[1:(m-1),1:(m-1)]%*%t(B[m,1:(m-1),drop=F])))
  #     B[m,] <- tau*B[m,]
  #   }
  # }
  # for (m in c((q1+1):(q1+q2))){
  #   for (j in c(1:(m-1))){
  #     Sll[j,m] <- Sll[m,j] <- B[m,1:(m-1),drop=F]%*%Sll[1:(m-1),j,drop=F]
  #   }
  # }
  for (m in c((q1+1):(q1+q2))){
    tau <- sqrt(R2[m-q1]/(B[m,1:(m-1),drop=F]%*%Sll[1:(m-1),1:(m-1)]%*%t(B[m,1:(m-1),drop=F])))
    B[m,] <- tau*B[m,]
    for (j in c(1:(m-1))){
      Sll[j,m] <- Sll[m,j] <- B[m,1:(m-1),drop=F]%*%Sll[1:(m-1),j,drop=F]
    }
  }
  
  Setaeta <- Sll[(q1+1):(q1+q2),(q1+1):(q1+q2)]
  
  # construction of Syy, W2
  
  Syy <- matrix(rep(c(p2:0),p2),p2,p2)*0.01
  Syy[upper.tri(Syy)] <- 0
  Syy <- (Syy+t(Syy))/p2
  diag(Syy) <- 1
  
  
  W2 <- matrix(indicatory,p2,q2)                # construction of W2
  W2 <- 1*(W2 == matrix(c(1:q2),p2,q2,byrow=T))
  P <- cumsum(colSums(W2))
  P <- c(0,P)
  P2 <- P
  W2 <- W2*wy
  for (i in c(1:q2)){                                      # scaling of weights
    subSyy <- Syy[(P[i] + 1):P[i+1],(P[i] + 1):P[i+1]]
    f <-  t(W2[(P[i] + 1):P[i+1],i,drop=F])%*%subSyy%*%W2[(P[i] + 1):P[i+1],i,drop=F]
    W2[(P[i] + 1):P[i+1],i] <- W2[(P[i] + 1):P[i+1],i]/sqrt(as.numeric(f))
  }
  
  if (q2 > 1){
    for (j in c(2:q2)){                                     # scaling of off-diagonal of Syy
      for (i in c(1:(j-1))){
        subSyy <-  Syy[(P[i]+1):P[i+1],(P[j]+1):P[j+1]]
        f <-  t(W2[(P[i] + 1):P[i+1],i,drop=F])%*%subSyy%*%W2[(P[j] + 1):P[j+1],j,drop=F]
        f <- as.numeric(Setaeta[i,j])/as.numeric(f)
        Syy[(P[i]+1):P[i+1],(P[j]+1):P[j+1]] <- subSyy*f
        Syy[(P[j]+1):P[j+1],(P[i]+1):P[i+1]] <- t(subSyy*f)
      }
    }}
  
  Br1 <- B[(q1+1):q,1:q1,drop=FALSE]
  Br2 <- B[(q1+1):q,(q1+1):q,drop=FALSE]
  IB <- solve(diag(q2) - t(Br2))
  
  Sxy <- Sxx%*%W1%*%t(Br1)%*%IB%*%solve(Setaeta)%*%t(W2)%*%Syy
  S <- rbind(cbind(Sxx,Sxy),cbind(t(Sxy),Syy))
  
  out <- list(S=S,B=B,Scomp=Sll,wx=rowSums(W1),wy=rowSums(W2))
  return(out)
}
