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
#
miris <- scale(as.matrix(iris[,-5]), center = TRUE, scale = FALSE)
miris <- as(miris, "TsparseMatrix")
M <- 2
N <- nrow(iris)
W0 <- matrix(rnorm(4*2), 4, 2)
Z0 <- matrix(rnorm(N*2), N, 2)
ca5 <- rgb(0,0,0,0.5)
out <- VBPCA(Z0, W0, miris, rank = 2, iter = 10, prior_prec = 1)
col3 <- hcl.colors(3, palette = "Set 2", alpha = 0.9)
plot(out$mean_z, col=col3[iris$Species], pch=16)
plot(out$mean_z%*%t(out$mean_w),as.matrix(Y))
abline(0,1,col=ca5)

