library(Matrix)
library(VBspPCA)
#library(onlinePCA)
###
#util
mat4img <- function(mat1){
  t(apply(mat1, 2, rev))
}
###

Vx <- t(W0)%*%W0
plot(solve(Vx, t(Y%*%W0)),VBspPCA:::ABsol(Vx,t(Y%*%W0)))
#Y%*%(W0)



doVB <-function(Y, M, X, W, maxit=10){
  #set.seed(seed)
  N <-nrow(Y)
  D <-ncol(Y)
  # X <- matrix(rnorm(M*N),N,M)
  # W <- matrix(rnorm(M*D),M,D)    
  #}
  S_inv <- diag(1,M)
  lp <- numeric(maxit)
  #pb <- txtProgressBar(min = 1, max = maxit, style = 3)
  for(i in 1:maxit){
    Vw <- t(X)%*%X+S_inv
    W <- solve(Vw, t(X)%*%Y)
    Vx <- W%*%t(W)+S_inv
    X <- t(solve(Vx, t(Y%*%t(W))))
    lp[i] <- sum((Y - X%*%W)^2)
    #setTxtProgressBar(pb, i)
  }
  list(W=W, X=X, Vw=Vw, Vx=Vx, lp=lp)
}
#
miris <- scale(as.matrix(iris[,-5]), center = TRUE, scale = FALSE)
miris <- as(miris, "TsparseMatrix")

M <- 2
N <- nrow(iris)
W0 <- matrix(rnorm(4*2), 4, 2)
Z0 <- matrix(rnorm(N*2), N, 2)

ca5 <- rgb(0,0,0,0.5)

out <- VBPCA(Z0, W0, miris, rank = 2, iter = 1, prior_prec = 1)
head.matrix(out$mean_z)
Y<- as.matrix(miris)
out0 <- doVB(Y, 2, Z0, t(W0), maxit = 1)

Y[1,]%*%t(out0$W)

head(out0$Vx)
t(out$prec_z)

plot(out$mean_z,out0$X,col=rep(1:2, each=150))

out$prec_w
out0$Vw
out$prec_z
out0$Vx

range(miris@i)

head(out0$X)
head.matrix(out$mean_z)

W0
Z0
W <- matrix(0, 4,2)
Z <- matrix(0, N,2)
for(i in 1:length(miris@x)){
  W[miris@j[i]+1,] <- W[miris@j[i]+1,] + miris@x[i]*Z0[miris@i[i]+1,]
  Z[miris@i[i]+1,] <- Z[miris@i[i]+1,] + miris@x[i]*W0[miris@j[i]+1,]
}

W
t(Z0)%*%miris

all(Z==miris%*%(W0))



#####

out <- doVB(Y,2)
plot(out$X%*%out$W, Y, col=rgb(0,0,0,0.5))
plot(out$X, col=iris$Species)
###

doVB <-function(Y, M, maxit=100, init="Nystrom"){
  #set.seed(seed)
  N <-ncol(Y)
  D <-nrow(Y)
  # if(init=="Nystrom"){
  #   O <- matrix(rnorm(M*D),M,D)
  #   X <- O%*%Y
  #   O <- matrix(rnorm(N*M),N,M)
  #   W <- Y%*%O    
  # }else{
  X <- matrix(rnorm(M*N),M,N)
  W <- matrix(rnorm(M*D),D,M)    
  #}
  I_D <- diag(1,D)
  S_Dinv <- diag(1,D)
  S_Winv <- diag(1,M)
  I_M <- diag(1,M)
  lp <- numeric(maxit)
  pb <- txtProgressBar(min = 1, max = maxit, style = 3)
  for(i in 1:maxit){
    R <- Y - W%*%X
    mu <-drop(rowSums(R)%*%solve(N*I_D+S_Dinv))
    W <- t(solve(X%*%t(X)+S_Winv,t((Y-mu)%*%t(X))))
    X <- solve(t(W)%*%W+I_M, t(t(Y-mu)%*%W))
    lp[i] <- sum((R - mu)^2) 
    setTxtProgressBar(pb, i)
  }
  list(W=W, X=X, mu=mu, lp=lp)
}

M <-2
Y <- t(as.matrix(iris[,-5]))
out <- doVB(Y,3)
plot(out$W%*%out$X+out$mu,Y)

####

solve(out$prec_z)
solve(out0$Vx)

solve(out$prec_z , t(out$mean_z)%*%out$mean_z)
out$mean_w
out$prec_w

#####
doVB <-function(Y, M,X,W, maxit=100){
  N <-nrow(Y)
  D <-ncol(Y)
  # X <- matrix(rnorm(M*N),N,M)
  # W <- matrix(rnorm(M*D),M,D)
  lp <- numeric(maxit)
  pb <- txtProgressBar(min = 1, max = maxit, style = 3)
  for(i in 1:maxit){
    W <- sweep(t(X)%*%Y, 2, colSums(X*X)+1, FUN = "/")
    X <- t(sweep(t(Y%*%t(W)), 2, rowSums(W*W)+1, FUN = "/"))
    lp[i] <- sum((Y - X%*%W)^2) 
    setTxtProgressBar(pb, i)
  }
  list(W=W, X=X, lp=lp)
}
