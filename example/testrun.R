library(Matrix)
library(VBspPCA)
#library(loon.data)
#library(oocPCA)
###
#M<-2
miris <- t(as.matrix(iris[,-5]))
miris <- as(miris, "TsparseMatrix")
ca5 <- rgb(0,0,0,0.5)
col3 <- hcl.colors(3, palette = "Set 2", alpha = 0.9)
system.time({
  out <- VBPCA(miris, rank = 2, iter = 50, prior_prec = 1, a=1, b=1,
               use_rowintercept = FALSE)
})

plot(out$logprob, type="l")
plot(out$mean_col, col=col3[iris$Species], pch=16)
plot(fit_pca(out), as.matrix(miris), col=ca5)
abline(0, 1, col=ca5)
rowMeans(miris)
sqrt(1/out$obs_prec)

