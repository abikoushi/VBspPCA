library(Matrix)
library(VBspPCA)
miris <- t(as.matrix(iris[,-5]))
dim(miris)
###

system.time({
  out <- VBPCA(miris, rank = 2, iter = 50, prior_prec = 1, a=1, b=1,
               use_rowintercept = TRUE)
})


plot(out$logprob[-1], type="l")
plot(out$mean_col, col=col3[iris$Species], pch=16)
plot(fit_pca(out), as.matrix(miris), col=ca5)
abline(0, 1, col="royalblue")
rowMeans(miris)
sqrt(1/out$obs_prec)


#writeMM(miris, "test.mtx")
size = VBspPCA:::size_mtx("test.mtx")
system.time({
  out_s = VBspPCA:::SVBPCA(file_path="test.mtx",
                           rank = 2, 
                           n_epochs = 100,
                           b_size = 100,
                           lr_param = c(1,0.9),
                           prior_shape = 1, prior_rate = 1,
                           use_rowintercept = TRUE)
})

out_s$mean_col
out_s$mean_row

out_s$obs_prec
plot(out_s$logprob, type="l")
ca5 <- rgb(0,0,0,0.5)
dim(miris)
col4 <- hcl.colors(4, palette = "Set 2", alpha = 0.9)
plot(fit_pca(out_s), as.matrix(miris), col=col4)
abline(0, 1, col=ca5)
col3 <- hcl.colors(3, palette = "Set 2", alpha = 0.9)
plot(out_s$mean_col, col=col3[iris$Species], pch=16)

####

system.time({
  out_s = VBspPCA:::SVBPCA(file_path="test.mtx",
                           rank = 2, 
                           n_epochs = 100,
                           b_size = 600,
                           lr_param = c(15,1),
                           prior_shape = 1, prior_rate = 1,
                           use_rowintercept = FALSE)
})

hist(out_s$mean_row)
out_s$obs_prec
plot(out_s$logprob, type="l")
ca5 <- rgb(0,0,0,0.5)
dim(miris)
col4 <- hcl.colors(4, palette = "Set 2", alpha = 0.9)
plot(fit_pca(out_s), as.matrix(miris), col=col4)
abline(0, 1, col=ca5)
col3 <- hcl.colors(3, palette = "Set 2", alpha = 0.9)
plot(out_s$mean_col, col=col3[iris$Species], pch=16)

