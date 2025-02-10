####
#Normal PCA
####
# SVBPCA <- function(file_path, rank,
#                    n_epochs,
#                    b_size,
#                    learning_rate = lr_default,
#                    dataloader = dataloader_mtx, 
#                    transformfun = I,
#                    prior_prec = 1,
#                    prior_shape = 1, prior_rate = 1,
#                    initZ=NULL, initW=NULL,
#                    subiter = 1,
#                    ...){
#   Dims <- size_mtx(file_path)
#   rowsize <- Dims[1]
#   colsize <- Dims[2]
#   N1 <- Dims[3]
#   b_size <- min(b_size, N1)
#   a <- prior_shape
#   b <- prior_rate
#   ahat <- prior_shape
#   bhat <- prior_rate
#   if(is.null(initZ)){
#     mu_z = matrix(rnorm(rowsize * rank), rowsize, rank)
#   }else{
#     mu_z = initZ
#   }
#   if(is.null(initW)){
#     mu_w = matrix(rnorm(colsize * rank), colsize, rank)
#   }else{
#     mu_w = initW
#   }
#   mu_B <- rnorm(rowsize)
#   Lambda_z = diag(rank)
#   Lambda_w = diag(rank)
#   Lambda_B <- 1
#   rind <-sample.int(N1)
#   if(N1%%b_size==0){
#     sb <- floor(N1/b_size)
#     subind <- vector("list",sb)
#     for(i in 1:sb){
#       subind[[i]] <- rind[1:b_size+b_size*(i-1)]  
#     }
#   }else{
#     sb <- floor(N1/b_size)
#     subind <- vector("list",sb+1)
#     for(i in 1:sb){
#       subind[[i]] <- rind[1:b_size+b_size*(i-1)]  
#     }
#     subind[[sb+1]] <- rind[(b_size*sb):length(rind)] 
#   }
#   #print(length(subind))
#   lp <- numeric(n_epochs)
#   pb <- txtProgressBar(0, n_epochs, style = 3)
#   #rho <- lr_default(1)
#   for(ep in 1:n_epochs){
#     rho <- learning_rate(ep, ...)
#     for(k in 1:length(subind)){
#       Y <- dataloader(file_path, subind[[k]])
#       out_t <- doVB_norm_s(y = transformfun(Y[,3]),
#                            rowi = Y[,1]-1L, coli = Y[,2]-1L,
#                            Nr = rowsize, Nc = colsize,
#                            L=rank, iter = subiter,
#                            prior_prec = prior_prec,
#                            a = a, b = b,
#                            N1 = N1,
#                            Z = mu_z, W = mu_w,
#                            B = mu_B,
#                            cov_z = Lambda_z, cov_w = Lambda_w,
#                            cov_B = Lambda_B,
#                            obs_prec = ahat/bhat)
#       Ns <- length(subind[[k]])
#       rho2 <- rho #*(Ns/N1)
#       rho1 <- 1-rho2
#       mu_z <- rho1*mu_z + rho2*out_t$mean_row
#       Lambda_z <- rho1*Lambda_z + rho2*out_t$cov_row
#       mu_w <- rho1*mu_w + rho2*out_t$mean_col
#       Lambda_w <- rho1*Lambda_z + rho2*out_t$cov_col
#       mu_B <- rho1*mu_B + rho2*out_t$mean_bias
#       Lambda_B <- rho1*Lambda_B + rho2*out_t$cov_bias
#       ahat <- rho1*ahat + rho2*out_t$prec_shape
#       bhat <-  rho1*ahat + rho2*out_t$prec_rate
#       lp[ep] <- lp[ep] + out_t$logprob[subiter]
#     }
#     setTxtProgressBar(pb,ep)
#   }
#   return(list(mean_row = mu_z,
#               mean_col = mu_w,
#               mean_bias = mu_B,
#               cov_row = Lambda_z,
#               cov_col = Lambda_w,
#               cov_bias = Lambda_B,
#               logprob = lp, 
#               family="normal"))
# }
