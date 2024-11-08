######
out0 <- doVB(miris, 2, Z0, t(W0), maxit = 100)

plot(out0$X%*%out0$W,Y)
abline(0,1)
plot(out0$lp,type = "l")

Y<-miris

head(out0$Vx)
t(out$prec_z)

plot(out$mean_z,out0$X,col=rep(1:2, each=150))

out$prec_w
out0$Vw
out$prec_z
out0$Vx
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

###

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

doVB <-function(Y, M, X, W, maxit=10){
  #set.seed(seed)
  N <-nrow(Y)
  D <-ncol(Y)
  # X <- matrix(rnorm(M*N),N,M)
  # W <- matrix(rnorm(M*D),M,D)
  S_inv <- diag(1,M)
  lp <- numeric(maxit)
  #pb <- txtProgressBar(min = 1, max = maxit, style = 3)
  for(it in 1:maxit){
    eta_w <- matrix(0, M, D)
    eta_x <- matrix(0, N, M)
    for(i in 1:length(Y@x)){
      eta_w[,Y@j[i]+1] <- eta_w[,Y@j[i]+1] + Y@x[i]*X[Y@i[i]+1,]
      eta_x[Y@i[i]+1,] <- eta_x[Y@i[i]+1,] + Y@x[i]*W[,Y@j[i]+1]
    }
    Vw <- t(X)%*%X+S_inv
    X <- t(solve(Vx, t(eta_x)))
    Vx <- W%*%t(W)+S_inv
    W <- solve(Vw, eta_w)
    lp[it] <- sum((Y - X%*%W)^2)
    #setTxtProgressBar(pb, i)
  }
  list(W=W, X=X, Vw=Vw, Vx=Vx, lp=lp)
}