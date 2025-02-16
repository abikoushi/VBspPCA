###
set_data <- function(L, nrow, ncol, center=0, scale=1){
  Z <- matrix(rnorm(L*nrow,0,scale), nrow, L)
  Z <- sweep(Z,1,rowMeans(Z)-center)
  W <- matrix(rnorm(L*ncol,0,scale), L, ncol)
  W <- sweep(W,1,rowMeans(W)-center)
  Y <- matrix(rpois(nrow*ncol, exp(Z)%*%exp(W)), nrow, ncol)
  list(Y=Y,Z=Z,W=t(W))
}

set_data <- function(L, nrow, ncol, center=0, scale=1){
  Z <- matrix(rnorm(L*nrow,0,scale), nrow, L)
  Z <- sweep(Z,1,rowMeans(Z)-center)
  W <- matrix(rnorm(L*ncol,0,scale), L, ncol)
  W <- sweep(W,1,rowMeans(W)-center)
  Y <- matrix(rnorm(nrow*ncol, Z%*%W), nrow, ncol)
  Y[Y<0] <- 0
  list(Y=Y,Z=Z,W=t(W))
}

set.seed(123); dat <- set_data(2, 199, 205)
mean(dat$Y==0)
length(dat$Y)
Y <- as(log1p(dat$Y), "TsparseMatrix")

writeMM(Y, "test.mtx")
hist(dat$Y, breaks = "FD")
system.time({
  out_s <- VBspPCA:::SVBPCA("test.mtx", rank = 2,
                  b_size = 10000,
                  subiter = 2,
                  n_epochs = 20,  forgetting=0.9, delay=1)
})

plot(out_s$logprob, type="l")

fit <- out_s$mean_row%*%t(out_s$mean_col)
plot(fit, as.matrix(Y), pch=".")
abline(0,1)

curve(cos(x),0,pi)
