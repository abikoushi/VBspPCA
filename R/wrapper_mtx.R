####
#Normal PCA
####
SVBPCA <- function(file_path, rank,
                   subiter = 5,
                   n_epochs = 100,
                   b_size = 10000,
                   learning_rate = lr_default,
                   dataloader = dataloader_mtx, 
                   transformfun = log1p,
                   prior_prec = 1,
                   prior_shape = 1, prior_rate = 1,
                   initZ=NULL, initW=NULL,
                   ...){
  #get matrix size
  con <- file(file_path, open = "r") #Open for reading in text mode
  rowsize <- scan(con, what=integer(), comment.char = "%", nmax=1,  quiet = TRUE)
  colsize <- scan(con, what=integer(), nmax=1, quiet = TRUE)
  N1 <- scan(con, what=integer(), nmax=1, quiet = TRUE)
  close(con)
  b_size <- min(b_size, N1)
  a <- prior_shape
  b <- prior_rate
  ahat <- prior_shape
  bhat <- prior_rate
  if(is.null(initZ)){
    mu_z = matrix(rnorm(rowsize * rank), rowsize, rank)    
  }else{
    mu_z = initZ
  }
  if(is.null(initW)){
    mu_w = matrix(rnorm(colsize * rank), colsize, rank)
  }else{
    mu_w = initW
  }
  Lambda_z = diag(rank)
  Lambda_w = diag(rank)
  rind <-sample.int(N1)
  if(N1%%b_size==0){
    sb <- floor(N1/b_size)
    subind <- vector("list",sb)
    for(i in 1:sb){
      subind[[i]] <- rind[1:b_size+b_size*(i-1)]  
    }
  }else{
    sb <- floor(N1/b_size)
    subind <- vector("list",sb+1)
    for(i in 1:sb){
      subind[[i]] <- rind[1:b_size+b_size*(i-1)]  
    }
    subind[[sb+1]] <- rind[(b_size*sb):length(rind)] 
  }
  #print(length(subind))
  lp <- numeric(n_epochs)
  pb <- txtProgressBar(0, n_epochs, style = 3)
  for(ep in 1:n_epochs){
    rho <- learning_rate(ep, ...)
    for(k in 1:length(subind)){
      Y <- dataloader(file_path, subind[[k]])
      out_t <- doVB_norm_s(y = transformfun(Y[,3]), rowi = Y[,1]-1L, coli = Y[,2]-1L,
                           Nr = rowsize, Nc = colsize,
                           L=rank, iter = subiter,
                           prior_prec = prior_prec,
                           a = a, b = b,
                           N1 = N1,
                           mu_z, mu_w, 
                           Lambda_z, Lambda_w, ahat/bhat)
      Ns <- length(subind[[k]])
      rho2 <- rho*(Ns/N1)
      mu_z <- (1-rho2)*mu_z + rho2*out_t$mean_row
      Lambda_z <- (1-rho2)*Lambda_z + rho2*out_t$prec_row
      mu_w <- (1-rho2)*mu_w + rho2*out_t$mean_col
      Lambda_w <- (1-rho2)*Lambda_z + rho2*out_t$prec_col
      ahat <- (1-rho2)*ahat + rho2*out_t$prec_shape
      bhat <-  (1-rho2)*ahat + rho2*out_t$prec_rate
      lp[ep] <- lp[ep] + out_t$logprob[subiter]
    }
    setTxtProgressBar(pb,ep)
  }
  return(list(mean_row = mu_z,
              mean_col = mu_w,
              prec_row = Lambda_z,
              prec_col = Lambda_w,
              logprob = lp, 
              family="normal"))
}
