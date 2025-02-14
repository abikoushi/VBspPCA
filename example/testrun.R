library(Matrix)
library(VBspPCA)
library(loon.data)
#library(oocPCA)
###
#M<-2
miris <- t(as.matrix(iris[,-5]))
miris <- as(miris, "TsparseMatrix")
ca5 <- rgb(0,0,0,0.5)
col3 <- hcl.colors(3, palette = "Set 2", alpha = 0.9)

out <- VBPCA(miris, rank = 2, iter = 50, prior_prec = 1, a=1, b=1)
out$obs_prec
plot(out$mean_col, col=col3[iris$Species], pch=16)
plot(sweep(out$mean_row%*%t(out$mean_col), 1, out$mean_bias,"+"), as.matrix(miris), col=ca5)
abline(0, 1, col=ca5)
plot(out$logprob, type="l")
rowMeans(miris)
out$mean_bias
sqrt(1/out$obs_prec)

file_path = "test.mtx"
rank = 2
b_size = 100
subiter = 5
n_epochs = 5
prior_prec = 1
prior_shape = 1
prior_rate = 1

writeMM(miris, "test.mtx")



size = VBspPCA:::size_mtx("test.mtx")

system.time({
out_s = VBspPCA:::doVB_norm_s_mtx(file_path="test.mtx",
                     Nr = size[1],
                     Nc = size[2],
                     N1 = size[3],
                     L = 2,
                     ns = 600,
                     iter = 2,
                     subiter = 1,
                     prior_prec = 1,
                     a = 1, b=1,
                     delay = 1,
                     forgetting=1)
})
#forgetting: (0.5, 1]
#delay: >0
out_s$obs_prec
plot(out_s$logprob, type="l")
plot(sweep(out_s$mean_row%*%t(out_s$mean_col), 1, out_s$mean_bias, "+"), as.matrix(miris), col=ca5)
abline(0, 1, col=ca5)
plot(out_s$mean_col, col=col3[iris$Species], pch=16)

# ###
# set_data <- function(L, nrow, ncol, center=0, scale=1){
#   Z <- matrix(rnorm(L*nrow,0,scale), nrow, L)
#   Z <- sweep(Z,1,rowMeans(Z)-center)
#   W <- matrix(rnorm(L*ncol,0,scale), L, ncol)
#   W <- sweep(W,1,rowMeans(W)-center)
#   Y <- matrix(rpois(nrow*ncol, exp(Z)%*%exp(W)), nrow, ncol)
#   list(Y=Y,Z=Z,W=t(W))
# }
# 
# set_data <- function(L, nrow, ncol, center=0, scale=1){
#   Z <- matrix(rnorm(L*nrow,0,scale), nrow, L)
#   Z <- sweep(Z,1,rowMeans(Z)-center)
#   W <- matrix(rnorm(L*ncol,0,scale), L, ncol)
#   W <- sweep(W,1,rowMeans(W)-center)
#   Y <- matrix(rnorm(nrow*ncol, Z%*%W), nrow, ncol)
#   Y[Y<0] <- 0
#   list(Y=Y,Z=Z,W=t(W))
# }
# 
# set.seed(123); dat <- set_data(2, 199, 205)
# mean(dat$Y==0)
# length(dat$Y)
# Y <- as(log1p(dat$Y), "TsparseMatrix")
# 
# writeMM(Y, "test.mtx")
# hist(dat$Y, breaks = "FD")
# system.time({
#   out_s <- VBspPCA:::SVBPCA("test.mtx", rank = 2,
#                   b_size = 10000,
#                   subiter = 2,
#                   n_epochs = 20,  forgetting=0.9, delay=1)
# })
# 
# plot(out_s$logprob, type="l")
# 
# fit <- out_s$mean_row%*%t(out_s$mean_col)
# plot(fit, as.matrix(Y), pch=".")
# abline(0,1)
