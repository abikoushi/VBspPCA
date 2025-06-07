library(Matrix)
library(VBspPCA)
library(ggplot2)
library(irlba)

datdir = scan("datapath.txt", what = character())
datpath = dir(datdir, full.names = TRUE)
dat = readMM(datpath[2])
L=30
res = prcomp(dat,rank. = 30)


#image(dat)
dim(dat)
#[1] 10147  8772
system.time({
  out_an <- VBspPCA:::VBPCA_diag(dat, 
                                 constr_type = "AN", 
                                 rank = 20, iter = 100,
                                 tau = 1, a = 1, b = 1)
})
#    user   system  elapsed 
#2289.476  133.604 2428.359
saveRDS(out_an, file = "Song_out_an.rds")
matplot(t(out_an$mean_col), type = "l")

plot(out_an$logprob[-1],type = "l")

plot(fit_pca(out_an), as.matrix(dat), pch=".")
abline(0,1)