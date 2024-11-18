lr_default <- function(t, delay=1, forgetting=0.9){
  (t+delay)^(-forgetting)
}

scan1_mtx <- function(con, skip = 0){
  base::scan(con, nmax = 1, quiet=TRUE, 
             what = list(i=integer(), j=integer(), v=numeric()),
             skip = skip)
}

dataloader_mtx <- function(file_path, bag){
  con <- file(file_path, open = "r") #Open for reading in text mode
  #get matrix size
  rowsize <- scan(con, what=integer(), comment.char = "%", nmax=1,  quiet = TRUE)
  colsize <- scan(con, what=integer(), nmax=1, quiet = TRUE)
  len <- scan(con, what=integer(), nmax=1, quiet = TRUE)
  bag <- sort(bag)
  newL <- scan1_mtx(con, skip = bag[1] - 1L) #initialize
  out <- matrix(0, length(bag), 3)
  out[1,] <- unlist(newL)
  bag <- diff(bag)
  for(i in 1:length(bag)){
    newL <- scan1_mtx(con, skip = bag[i] - 1L)
    out[i+1,] <- unlist(newL)
  }
  close(con)
  return(out)
}

SVBPCA <- function(file_path, rank,
                    subiter = 5,
                    n_epochs = 100,
                    b_size = 10000,
                    learning_rate = lr_default,
                    dataloader = dataloader_mtx, 
                    prior_prec = 1,
                    prior_shape = 1, prior_rate = 1, ...){
   #get matrix size
   con <- file(file_path, open = "r") #Open for reading in text mode
   rowsize <- scan(con, what=integer(), comment.char = "%", nmax=1,  quiet = TRUE)
   colsize <- scan(con, what=integer(), nmax=1, quiet = TRUE)
   N1 <- scan(con, what=integer(), nmax=1, quiet = TRUE)
   close(con)
   b_size <- min(b_size, N1)
   a <- prior_shape
   b <- prior_rate
   mu_z = matrix(rnorm(rowsize * rank), rowsize, rank)
   mu_w = matrix(rnorm(colsize * rank), colsize, rank)
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
       out_t <- doVB_norm_s(y = Y[,3], rowi = Y[,1]-1L, coli = Y[,2]-1L,
                            Nr = rowsize, Nc = colsize,
                            L=rank, iter = subiter,
                            prior_prec = prior_prec,
                            a = a, b = b,
                            N1 = N1,
                            mu_z, mu_w, 
                            Lambda_z, Lambda_w)
       Ns <- length(subind[[k]])
       rho2 <- rho*(Ns/N1)
       mu_z <- (1-rho2)*mu_z + rho2*out_t$mean_z
       Lambda_z <- (1-rho2)*Lambda_z + rho2*out_t$prec_z
       mu_w <- (1-rho2)*mu_w + rho2*out_t$mean_w
       Lambda_w <- (1-rho2)*Lambda_z + rho2*out_t$prec_w
       lp[ep] <- lp[ep] + out_t$logprob[subiter]
     }
     setTxtProgressBar(pb,ep)
   }
   return(list(mean_row = mu_z,
               mean_col = mu_w,
               prec_row = Lambda_z,
               prec_col = Lambda_w,
               logprob = lp))
}