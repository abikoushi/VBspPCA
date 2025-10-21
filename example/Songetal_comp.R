library(Matrix)
library(VBspPCA)
library(ggplot2)
library(irlba)
library(rsvd)
library(bench)
library(readr)
library(bench)
library(tidyr)
library(dplyr)

datdir = scan("datapath.txt", what = character())
datpath = dir(datdir[1], full.names = TRUE)
dat = readMM(datpath[2])


wch = which(apply(dat==0,1,all))

dat = dat[-wch,]
saveRDS(wch, file = "wch.rdata")
# res = prcomp(dat,rank. = 30)

writeBin(c(rbind(dat@i,dat@j)), "song_x.bin")
writeBin(c(log1p(dat@x)), "song_y.bin")
writeMM(dat, "Songetal_filt.mtx")
L=20
size = VBspPCA:::size_mtx( "Songetal_filt.mtx")
#   row   column  nonzero 
# 10147     8772 11584999 

# 11584999 / (10147 * 8772)
# [1] 0.1301547

load(file = "irlba_10.RData")
load(file = "prcomp_10.RData")

bm_irlba

# L=10
# bm_pr = bench::mark({res_pr <- prcomp(dat, rank. = L)},
#                     iterations = 1)

#save(bm_pr, res_pr, file = "prcomp_10.RData")
load("prcomp_10.RData")

df_cell = read_csv(datpath[1])
dim(res_pr$rotation)

df_Cell = data.frame(res_pr$rotation, 
                     cell_type = factor(df_cell$cell_type),
                     source = factor(df_cell$source),
                     sample_id = factor(df_cell$sample_id)) %>% 
  mutate(id=row_number()) %>% 
  pivot_longer(1:10, names_to = "component") %>% 
  mutate(component=as.integer(gsub('PC',"",component))) %>% 
  dplyr::filter(!is.na(cell_type))
head(df_Cell)

ggplot(df_Cell, aes(x=component, y=value, colour = cell_type, group=id))+
  geom_line(linewidth = 0.05)+
  guides(colour=guide_legend(override.aes = list(linewidth=2)))+
  scale_color_brewer(palette = "Set1")+
  theme_bw(16)
ggsave("pr_10.pdf", width = 7, height = 5)


bm_irlba = bench::mark({res_irlba <- prcomp_irlba(dat, n = L)},
                       iterations = 1)
save(bm_irlba, res_irlba, file = "irlba_10.RData")

bm_pr$total_time
bm_pr$mem_alloc

bm_irlba$total_time
bm_irlba$mem_alloc

# bm_an = bench::mark({
#   out_an <- VBspPCA:::VBPCA_diag_bin(read_x="song_x.bin",
#                                      read_y="song_y.bin",
#                                      rank = L, maxit=2, 
#                                      size = size,
#                                      constr_type = "AN",
#                                      tol = 0.1,
#                                      tau = 1, a = 1, b = 1)
# }, iterations = 1)
# out_an$logprob

mean(dat@x)
L = 10
bm_an = bench::mark({
  out_an <- VBspPCA:::VBPCA_diag(dat,
                                 constr_type = "AN",
                                 rank =  L, maxit = 10,
                                 lambda_ini = 1,
                                 tau = 1, a = 1, b = 1, 
                                 init_mean =  4)
}, iterations = 1)
bm_an
out_an = VBspPCA:::reorder_cols(out_an)
out_an
plot(out_an$logprob, type = "l")
bm_an$mem_alloc
bm_an$total_time

names(out_an)

matplot(t(out_an$mean_col), type = "l", lty=1, col=rgb(0,0,0,0.1))

# #    user   system  elapsed 
# #2289.476  133.604 2428.359
# saveRDS(out_an, file = "Song_out_an.rds")
# matplot(t(out_an$mean_col), type = "l")
# 
# plot(out_an$logprob[-1],type = "l")

ggplot(data=NULL, aes(x=c(fit_pca(out_an)), y=c(as.matrix(log1p(dat)))))+
  geom_abline(color="lightgrey") +
  geom_bin2d(aes(fill = after_stat(log1p(count)))) +
  theme_bw()
# 
# L = 15
# mu = 1
# sigma = sqrt(mu)*sqrt(pi)/sqrt(2)
# dot = numeric(10000)
# for(i in 1:10000){
#   x = abs(rnorm(L, 0, sigma/sqrt(L)))
#   y = abs(rnorm(L, 0, sigma/sqrt(L)))
#   dot[i] <- sum(x*y)
# }
# 
# mean(dot)
# hist(dot, breaks = "FD")
