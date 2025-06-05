
####

SVBPCA_diag <- function(Y, rank, iter, bsize,
                        lr_type = "exponential",
                        lr_param = c(1,0.9),
                        constr_type = "NN", 
                        tau=1, a=1, b=1,
                        lambda = 1,
                        display_progress=TRUE){
  if(any(class(Y)=="dgTMatrix")){
    Y <- as(Y, "TsparseMatrix")    
  }
  dims = dim(Y)
  V = lapply(dims, VBspPCA:::initnorm, rank=rank)
  
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
                        lr_param = c(15,0.9),
                        rank = 2, iter = 200,
                        lambda=200,
                        tau = 1, a = 1, b = 1)
})
print(out_nn$obs_prec)
#[1] 301

plot(out_nn$mean_col, col=col3[iris$Species], pch=16)
plot(fit_pca(out_nn), as.matrix(miris))
abline(0,1,lty=3)
out_nn$mean_row

plot(out_nn$logprob, type = "l")

