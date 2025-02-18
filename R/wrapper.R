VBPCA <- function(Y, rank, iter=10, prior_prec=1, a = 1, b = 1,
                  use_rowintercept = TRUE){
  if(!any(class(Y)=="dgTMatrix")){
    Y <- as(Y, "TsparseMatrix")
  }
  if(use_rowintercept){
    out = doVB_norm(y = Y@x, rowi = Y@i, coli = Y@j, 
              Nr = Y@Dim[1], Nc = Y@Dim[2], 
              L = rank,
              iter = iter,
              prior_prec = prior_prec, a=a, b=b)    
  }else{
    out = doVB_norm_woi(y = Y@x, rowi = Y@i, coli = Y@j, 
                    Nr = Y@Dim[1], Nc = Y@Dim[2], 
                    L = rank,
                    iter = iter,
                    prior_prec = prior_prec, a=a, b=b)    
  }
  return(out)
}

SVBPCA <- function(file_path, rank,
                   subiter = 1,
                   n_epochs = 100,
                   b_size = 10000,
                   prior_prec = 1,
                   prior_shape = 1, prior_rate = 1,
                   delay=1,
                   forgetting=0.8,
                   use_rowintercept = TRUE){
  size = size_mtx(file_path)
  if(use_rowintercept){
    out = doVB_norm_s_mtx(file_path,
                  size[1], size[2], size[3],
                  L = rank,
                  ns = b_size,
                  iter = n_epochs,
                  subiter = subiter,
                  prior_prec = prior_prec,
                  a=prior_shape, b=prior_rate,
                  delay=delay,forgetting=forgetting)
  }else{
    out = doVB_norm_wo_s_mtx(file_path,
                    size[1], size[2], size[3],
                    L = rank,
                    ns = b_size,
                    iter = n_epochs,
                    subiter = subiter,
                    prior_prec = prior_prec,
                    a=prior_shape, b=prior_rate,
                    delay=delay,forgetting=forgetting)
  }
  return(out)
}

fit_pca <- function(out){
  if(is.null(out$mean_intercept)){
    res = out$mean_row%*%t(out$mean_col)
  }else{
    res = sweep(out$mean_row%*%t(out$mean_col), 1, out$mean_intercept, "+")      
  }
  return(res)
}
