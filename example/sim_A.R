library(Matrix)
library(VBspPCA)
library(ggplot2)

L = 2
nr = 2000
nc = 1000
rv = seq(-4, 4, length.out=nr)
cv = seq(-4, 4, length.out=nc)

Y  = matrix(rnorm(nr*nc, plogis(rv)%o%plogis(cv), 1),nr,nc)
Y[Y<0] <- 0
image.default(Y)

res = prcomp(Y, rank. = L, center = FALSE)
matplot(res$x, type="l", lty=1)
matplot(res$rotation, type="l", lty=1)

res$scale
plot(res$x%*%t(res$rotation), Y)
plot(res$sdev)

system.time({
  out_an <- VBspPCA:::VBPCA_diag(Y, 
                                 constr_type = "AN", 
                                 lambda_ini = 1,
                                 rank = L,
                                 maxit = 100,
                                 tau = 500, a = 1, b = 1, tol = 0.1)
})
#   user  system elapsed 
# 22.030   0.630  22.655 
plot(out_an$logprob[-1], type = "l")
matplot(out_an$mean_row, type="l", lty=1)
matplot(out_an$mean_col, type="l", lty=1)

system.time({
  out_sn <- VBspPCA:::VBPCA_diag(Y, 
                                 constr_type = "SN", 
                                 lambda_ini = 1,
                                 rank = L,
                                 maxit = 100,
                                 tau = 500, a = 1, b = 1, tol = 0.1)
})

plot(out_sn$logprob[-1], type = "l")
matplot(out_sn$mean_row, type="l", lty=1)
matplot(out_sn$mean_col, type="l", lty=1)

system.time({
  out_nn <- VBspPCA:::VBPCA_diag(Y, 
                                 constr_type = "NN", 
                                 lambda_ini = 1,
                                 rank = L,
                                 maxit = 100,
                                 tau = 500, a = 1, b = 1, tol = 0.1)
})

plot(out_nn$logprob[-1], type = "l")
matplot(out_nn$mean_row, type="l", lty=1)
matplot(out_nn$mean_col, type="l", lty=1)
