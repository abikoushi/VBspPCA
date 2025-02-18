library(Matrix)
library(VBspPCA)
miris <- t(as.matrix(iris[,-5]))
#writeMM(miris, "test.mtx")
size = VBspPCA:::size_mtx("test.mtx")
system.time({
  out_s = VBspPCA:::SVBPCA(file_path="test.mtx",
                           rank = 2, subiter = 1,
                           n_epochs = 500,
                           b_size = 50,
                           delay=1, forgetting=0.8,
                           use_rowintercept = FALSE)
})
#forgetting: (0.5, 1]
#delay: >0
out_s$mean_row
sqrt(1/out_s$obs_prec)
plot(out_s$logprob, type="l")
ca5 <- rgb(0,0,0,0.5)
dim(miris)
col4 <- hcl.colors(4, palette = "Set 2", alpha = 0.9)
plot(fit_pca(out_s), as.matrix(miris), col=col4)
abline(0, 1, col=ca5)
col3 <- hcl.colors(3, palette = "Set 2", alpha = 0.9)
plot(out_s$mean_col, col=col3[iris$Species], pch=16)

