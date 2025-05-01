library(Matrix)
library(VBspPCA)

miris <- t(as.matrix(iris[,-5]))
miris <- as(log1p(miris), "TsparseMatrix")
ca5 <- rgb(0, 0, 0, 0.5)
col3 <- hcl.colors(3, palette = "Set 2", alpha = 0.9)
system.time({
  out_an <- VBPCA_diag(miris,  constr_type = "AN", 
                    rank = 2, iter = 100,
                    tau = 1, a = 1, b = 1)
})

system.time({
  out_sn <- VBPCA_diag(miris,  constr_type = "SN", 
                    rank = 2, iter = 100,
                    tau = 1, a = 1, b = 1)
})

system.time({
  out_nn <- VBPCA_diag(miris,  constr_type = "NN", 
                    rank = 2, iter = 100,
                    tau = 1, a = 1, b = 1)
})

matplot(cbind(out_an$logprob[-1], 
      out_sn$logprob[-1],
      out_nn$logprob[-1]), type="l", col = col3)

plot(out_an$mean_col, col=col3[iris$Species], pch=16)
plot(out_sn$mean_col, col=col3[iris$Species], pch=16)
plot(out_nn$mean_col, col=col3[iris$Species], pch=16)

col3 <- hcl.colors(3, palette = "Set 2", alpha = 0.5)
plot(fit_pca(out_an), as.matrix(miris), col=col3[1])
points(fit_pca(out_sn), as.matrix(miris), col=col3[2])
points(fit_pca(out_nn), as.matrix(miris), col=col3[3])
abline(0, 1, col="royalblue", lty=2)
