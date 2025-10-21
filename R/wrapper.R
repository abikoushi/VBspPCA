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
                   lr_type="exponential",
                   lr_param = c(15,0.8),
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
                  lr_param=lr_param,
                  lr_type=lr_type)
  }else{
    out = doVB_norm_wo_s_mtx(file_path,
                    size[1], size[2], size[3],
                    L = rank,
                    ns = b_size,
                    iter = n_epochs,
                    prior_prec = prior_prec,
                    a = prior_shape, b = prior_rate,
                    lr_param=lr_param,
                    lr_type=lr_type)
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

initnorm <- function(D, rank, sd=1){
  matrix(abs(rnorm(D*rank, 0, sd)), D, rank)
}

VBPCA_diag <- function(Y, rank, maxit, constr_type = "AN", 
                       lambda_ini = 1,
                       tau=1, a=1, b=1,
                       tol=0.01,
                       init_mean=1){
  if(!any(class(Y)=="dgTMatrix")){
    Y <- as(Y, "TsparseMatrix")
  }
  dims = dim(Y)
  init_sd = (sqrt(init_mean) * sqrt(pi))/sqrt(2)
  V = lapply(dims, initnorm, rank=rank, sd = init_sd/sqrt(rank))
  res = doVB_norm_woi_diag_om(V, lambda = lambda_ini, 
                              y=Y@x, X = cbind(Y@i, Y@j),
                              dims = dims,
                              L = rank,
                              constr_type = constr_type,
                              maxit = maxit,
                              tau=tau, a=a, b=b,
                              tol=tol)
  return(res)
}


VBPCA_diag_mtx <- function(readtxt, rank, maxit,
                           constr_type="AN",
                           lambda_ini = 1,
                           tau=1, a=1, b=1,
                           tol=0.01){
  size = size_mtx(readtxt)
  dims = size[1:2]
  N = prod(dims)
  N1 = size[3]
  V = lapply(dims, initnorm, rank=rank)
  res = doVB_norm_woi_diag_mtx(V, lambda = lambda_ini, 
                               readtxt,
                               dims = dims,
                               L = rank,
                               constr_type=constr_type,
                               maxit=maxit, tau=tau, a=a, b=b,
                               tol=tol)
  return(res)
}

VBPCA_diag_bin <- function(read_x, read_y,
                           rank, maxit, 
                           size,
                           constr_type="AN",
                           lambda_ini = 1,
                           tau=1, a=1, b=1,
                           tol=0.01){
  dims = size[1:2]
  N1 = size[3]
  V = lapply(dims, initnorm, rank=rank)
  res = doVB_norm_woi_diag_bin(V, 
                               lambda = lambda_ini, 
                               N1 = N1,
                               read_x, read_y,
                               dims = dims,
                               L = rank,
                               constr_type=constr_type,
                               maxit=maxit, tau=tau, a=a, b=b,
                               tol=tol)
  return(res)
}
