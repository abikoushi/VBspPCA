file_path = "test.mtx"
rank = 2
b_size = 100
subiter = 5
n_epochs = 5
prior_prec = 1
prior_shape = 1
prior_rate = 1
library(Matrix)
library(VBspPCA)
writeMM(miris, "test.mtx")
size = VBspPCA:::size_mtx("test.mtx")
system.time({
  out_s = VBspPCA:::SVBPCA(file_path="test.mtx",
                           rank = 2,subiter = 1,
                           n_epochs = 100,
                           b_size = 100,
                           delay=1,forgetting=0.8)
})
#forgetting: (0.5, 1]
#delay: >0
1/out_s$obs_prec
plot(out_s$logprob, type="l")
plot(sweep(out_s$mean_row%*%t(out_s$mean_col), 1, out_s$mean_bias, "+"), as.matrix(miris), col=ca5)
abline(0, 1, col=ca5)
plot(out_s$mean_col, col=col3[iris$Species], pch=16)

