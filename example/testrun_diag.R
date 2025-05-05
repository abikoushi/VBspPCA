library(Matrix)
library(VBspPCA)

miris <- t(as.matrix(iris[,-5]))
miris <- as(log1p(miris), "TsparseMatrix")
ca5 <- rgb(0, 0, 0, 0.5)
col3 <- hcl.colors(3, palette = "Set 2", alpha = 0.9)

SVBPCA_diag <- function(Y, rank, iter, bsize,
                        lr_type = "exponential",
                        lr_param = c(1,0.9),
                        constr_type = "AN", 
                        tau=1, a=1, b=1, display_progress=TRUE){
  if(any(class(Y)=="dgTMatrix")){
    Y <- as(Y, "TsparseMatrix")    
  }
  dims = dim(Y)
  V = lapply(dims, VBspPCA:::initnorm, rank=rank)
  lambda = 1
  res = VBspPCA:::doSVB_norm_woi_diag(V, lambda = lambda, 
                                     y=Y@x, X = cbind(Y@i, Y@j),
                                     dims = dims,
                                     L = rank,
                                     constr_type = constr_type,
                                     lr_type = lr_type,
                                     lr_param = lr_param,
                                     bsize = bsize,
                                     iter=iter, tau=tau, a=a, b=b,
                                     display_progress = display_progress)
  return(res)
}

dim(miris)
system.time({
  out_nn <- SVBPCA_diag(miris,  bsize = 100,
                        constr_type = "NN", 
                        lr_param = c(20,0.9),
                        rank = 2, iter = 10,
                        tau = 1, a = 1, b = 1)
})
1/sqrt(out_nn$obs_prec)

plot(out_nn$mean_col, col=col3[iris$Species], pch=16)
plot(fit_pca(out_nn), as.matrix(miris))
out_nn$mean_row
out_nn$obs_prec
plot(out_nn$logprob, type = "l")

#####
system.time({
  out_an <- VBspPCA:::VBPCA_diag(miris, 
                                 constr_type = "AN", 
                                 rank = 2, iter = 100,
                                 tau = 1, a = 1, b = 1)
})

system.time({
  out_sn <- VBspPCA:::VBPCA_diag(miris,  constr_type = "SN", 
                    rank = 2, iter = 100,
                    tau = 1, a = 1, b = 1)
})

system.time({
  out_nn <- VBspPCA:::VBPCA_diag(miris,  constr_type = "NN", 
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
abline(0, 1, col="grey", lty=3)
