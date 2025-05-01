PCA_vb <-function(Y, L, iter=1000L, seed=1L){
  set.seed(seed)
  tau=1
  N <-nrow(Y)
  D <-ncol(Y)
  ND = N*D
  ahat = ND/2+1
  X <- matrix(rnorm(L*N),N,L)
  W <- matrix(rnorm(L*D),L,D)
  I_L <- diag(1, L)
  tau_hist <- numeric(iter)
  pb <- txtProgressBar(min = 1, max = iter, style = 3)
  for(i in 1:iter){
    Lam_w <- tau*t(X)%*%X+I_L
    W <- solve(Lam_w, tau*t(X)%*%Y)
    Lam_x <- tau*W%*%t(W)+I_L
    X <- t(solve(Lam_x, tau*W%*%t(Y)))
    SE=sum( (Y-X%*%W)^2 )
    tau = ahat/(SE+1)
    tau_hist[i] <- tau
    setTxtProgressBar(pb, i)
  }
  list(W=W, X=X, tau=tau, tau_hist=tau_hist)
}

Y <- as.matrix(iris[,-5])
out <- PCA_vb(Y, 2, 10)
1/out$tau
plot(out$tau_hist[-1], type = "l")

(mean(c(Y-out$X%*%out$W)^2))

plot(Y, out$X%*%out$W, cex=0.5, col=rgb(0,0,0,0.5))
abline(0,1,col="royalblue")
plot(out$X, col=as.integer(iris$Species)+1)

####
library(Matrix)
Y = as(Y, "TsparseMatrix")
L=2
PCA_vb <-function(Y, L, iter=1000L, seed=1L){
  set.seed(seed)
  tau=1
  N <-nrow(Y)
  D <-ncol(Y)
  ND = N*D
  ahat = ND/2+1
  V1 <- matrix(rnorm(L*N),N,L)
  V2 <- matrix(rnorm(L*D),D,L)
  I_L <- diag(1, L)
  tau_hist <- numeric(iter)
  pb <- txtProgressBar(min = 1, max = iter, style = 3)
  for(i in 1:iter){
    Lam_w <- tau*t(X)%*%X+I_L
    W <- solve(Lam_w, tau*t(X)%*%Y)
    Lam_x <- tau*W%*%t(W)+I_L
    X <- t(solve(Lam_x, tau*W%*%t(Y)))
    SE=sum( (Y-X%*%W)^2 )
    tau = ahat/(SE+1)
    tau_hist[i] <- tau
    setTxtProgressBar(pb, i)
  }
  list(W=W, X=X, tau=tau, tau_hist=tau_hist)
}