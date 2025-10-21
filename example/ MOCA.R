library(Matrix)
library(VBspPCA)
library(ggplot2)
library(irlba)
#devtools::install_github("KlugerLab/oocPCA")
datdir = scan("datapath.txt", what = character())
datpath = dir(datdir[2], full.names = TRUE)
#dat = readMM(datpath[3])

VBspPCA:::size_mtx(datpath[3])
#   row    column   nonzero 
# 26183   1331984 903942242 

datpath = dir(datdir[3], full.names = TRUE)
row_meanvar = readRDS(datpath[1])

cumvar = cumsum(sort(row_meanvar$var, decreasing = TRUE))/sum(row_meanvar$var)
plot(seq_len(length.out = length(row_meanvar$var)), 
     cumvar,type="l")
abline(v = 2000, lty=2)

size = VBspPCA:::size_mtx(datpath[2])
M = readMM(datpath[2])

writeBin(c(rbind(M@i,M@j)), "moca_x.bin")
writeBin(c(log1p(M@x)), "moca_y.bin")


L=20

system.time({
  out_an <- VBspPCA:::VBPCA_diag_bin(read_x="moca_x.bin",
                                     read_y="moca_y.bin",
                                     rank = L, maxit=2, 
                                     size = size,
                                     constr_type = "AN",
                                     tol = 0.1,
                                     tau = 1, a = 1, b = 1)
})


# system.time({
#   out_an <- VBspPCA:::VBPCA_diag_mtx(datpath[2],  constr_type = "AN",
#                                      rank = L, maxit = 2,
#                                      tol = 0.1,
#                                      tau = 1, a = 1, b = 1)
# })
