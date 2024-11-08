library(Matrix)
library(VBspPCA)
#library(onlinePCA)
###
#util
mat4img <- function(mat1){
  t(apply(mat1, 2, rev))
}
###
#M<-2
miris <- scale(as.matrix(iris[,-5]), center = TRUE, scale = FALSE)
miris <- as(miris, "TsparseMatrix")
ca5 <- rgb(0,0,0,0.5)
col3 <- hcl.colors(3, palette = "Set 2", alpha = 0.9)

out <- VBPCA(miris, rank = 3, iter = 10, prior_prec = 1)
plot(out$mean_z, col=col3[iris$Species], pch=16)
plot(out$mean_z%*%t(out$mean_w),as.matrix(miris), col=ca5)
abline(0,1,col=ca5)
plot(out$logprob, type="l")

a=1; b=1



N1 <- length(miris@x)
ind <- sort(sample.int(N1, 100))
subx <- miris@x[ind]
subi <- miris@i[ind]
subj <- miris@j[ind]

N <- prod(miris@Dim)

out_t <- VBspPCA:::doVB_norm_s(y = subx, rowi = subi, coli = subj,
                               Nr = miris@Dim[1], Nc = miris@Dim[2],
                               L=2, iter = 5,
                               prior_prec = 1,
                               a, b,
                               N1 = N1)

plot(out_t$logprob, type = "l")
plot(out$mean_z%*%t(out$mean_w), as.matrix(miris), col=ca5)
abline(0,1,col=ca5)



###
set_data <- function(L, nrow, ncol, center=0, scale=1){
  Z <- matrix(rnorm(L*nrow,0,scale), nrow, L)
  Z <- sweep(Z,1,rowMeans(Z)-center)
  W <- matrix(rnorm(L*ncol,0,scale), L, ncol)
  W <- sweep(W,1,rowMeans(W)-center)
  Y <- matrix(rpois(nrow*ncol, exp(Z)%*%exp(W)), nrow, ncol)
  list(Y=Y,Z=Z,W=t(W))
}

set.seed(123); dat <- set_data(2, 199, 205)
mean(dat$Y==0)
length(dat$Y)
Y <- as(dat$Y, "TsparseMatrix")

writeMM(log1p(Y), "test.mtx")
hist(dat$Y, breaks = "FD")
system.time({
  out_s <- VBspPCA:::SVBPCA("test.mtx", rank = 2,
                  b_size = 10000,
                  subiter = 2,
                  n_epochs = 20,  forgetting=0.9, delay=1)
})

plot(out_s$logprob, type="l")

fit <- out_s$mean_row%*%t(out_s$mean_col)
plot(fit, as.matrix(log1p(Y)), pch=".")

