library(Matrix)
library(VBspPCA)
library(ggplot2)
library(irlba)
library(readr)
library(bench)
library(tidyr)
library(dplyr)

datdir = scan("datapath.txt", what = character())
datpath = dir(datdir[1], full.names = TRUE)
dat = readMM(datpath[2])

nnzero(dat)/length(dat)

png("heatmap.png")
image(dat)
dev.off()

wch = which(apply(dat==0,1,all))

dat = dat[-wch,]
saveRDS(wch, file = "wch.rdata")

mean(log1p(dat))
L = 10
bm_nn = bench::mark({
  out_nn <- VBspPCA:::VBPCA_diag(log1p(dat),
                                 constr_type = "NN",
                                 rank =  L, maxit = 10,
                                 lambda_ini = 10,
                                 tau = 1, a = 1, b = 1, 
                                 init_mean =  0.18)
}, iterations = 1)
bm_nn
plot(out_nn$logprob, type = "l")
#out_nn$obs_prec
#[1] 5.584225
out_nn = VBspPCA:::reorder_cols(out_nn)

df_cell = read_csv(datpath[1])
head(df_cell$source)

df_Cell = data.frame(out_nn$mean_col , 
                     cell_type = factor(df_cell$cell_type),
                     source = factor(df_cell$source),
                     sample_id = factor(df_cell$sample_id)) %>% 
  mutate(id=row_number()) %>% 
  pivot_longer(1:10, names_to = "component") %>% 
  mutate(component=as.integer(gsub('X',"",component))) %>% 
  dplyr::filter(!is.na(cell_type))
head(df_Cell)

ggplot(df_Cell, aes(x=component, y=value, colour = cell_type, group=id))+
  geom_line(linewidth = 0.05)+
  guides(colour=guide_legend(override.aes = list(linewidth=2)))+
  scale_color_brewer(palette = "Set1")+
  theme_bw(16)
ggsave("pcp_nn10.pdf", width = 7, height = 5)

ggplot(df_Cell, aes(x=component, y=value, colour = sample_id, group=id))+
  geom_line(linewidth = 0.025)+
  guides(colour=guide_legend(override.aes = list(linewidth=2)))+
  scale_color_viridis_d()+
  theme_bw(16)
ggsave("pcp_nn10_cc.pdf", width = 7, height = 5)

# matplot(t(out_nn$mean_col), type = "l", lty=1, col=factor(df_cell$cell_type), lwd=0.01)
#saveRDS(out_nn, file = "Song_out_nn.rds")



ggplot(data=NULL, aes(x=c(fit_pca(out_nn)), y=c(as.matrix(log1p(dat)))))+
  geom_abline(color="lightgrey") +
  geom_bin2d(aes(fill = after_stat(log1p(count)))) +
  theme_bw()
ggsave("fit.pdf")

bm_sn = bench::mark({
  out_sn <- VBspPCA:::VBPCA_diag(log1p(dat),
                                 constr_type = "SN",
                                 rank =  L, maxit = 10,
                                 lambda_ini = 5,
                                 tau = 1, a = 1, b = 1, 
                                 init_mean =  0.18)
}, iterations = 1)

save(bm_sn,out_sn, file = "out_sn10.Rdata")
out_sn$mean_col
out_sn$H
out_sn$obs_prec
out_sn$logprob

plot(out_sn$logprob)

L = 20
bm_nn = bench::mark({
  out_nn <- VBspPCA:::VBPCA_diag(log1p(dat),
                                 constr_type = "NN",
                                 rank =  L, maxit = 10,
                                 lambda_ini = 10,
                                 tau = 1, a = 1, b = 1, 
                                 init_mean =  0.18)
}, iterations = 1)

plot(out_nn$logprob)
save(bm_nn,out_nn, file = "out_nn20.Rdata")
###
out_nn = VBspPCA:::reorder_cols(out_nn)
df_cell = read_csv(datpath[1])
head(df_cell$source)

df_Cell = data.frame(out_nn$mean_col , 
                     cell_type = factor(df_cell$cell_type),
                     source = factor(df_cell$source),
                     sample_id = factor(df_cell$sample_id)) %>% 
  mutate(id=row_number()) %>% 
  pivot_longer(1:20, names_to = "component") %>% 
  mutate(component=as.integer(gsub('X',"",component))) %>% 
  dplyr::filter(!is.na(cell_type))
head(df_Cell)

ggplot(df_Cell, aes(x=component, y=value, colour = cell_type, group=id))+
  geom_line(linewidth = 0.025)+
  guides(colour=guide_legend(override.aes = list(linewidth=2)))+
  scale_color_brewer(palette = "Set1")+
  theme_bw(16)
ggsave("pcp_nn20.pdf", width = 7, height = 5)

####
out_sn = VBspPCA:::reorder_cols(out_sn)
df_cell = read_csv(datpath[1])
head(df_cell$source)

df_Cell = data.frame(out_sn$mean_col , 
                     cell_type = factor(df_cell$cell_type),
                     source = factor(df_cell$source),
                     sample_id = factor(df_cell$sample_id)) %>% 
  mutate(id=row_number()) %>% 
  pivot_longer(1:10, names_to = "component") %>% 
  mutate(component=as.integer(gsub('X',"",component))) %>% 
  dplyr::filter(!is.na(cell_type))
head(df_Cell)
plot(out_sn$logprob)

ggplot(data=NULL, aes(x=c(fit_pca(out_sn)), y=c(as.matrix(log1p(dat)))))+
  geom_abline(color="lightgrey") +
  geom_bin2d(aes(fill = after_stat(log1p(count)))) +
  theme_bw()
ggplot(df_Cell, aes(x=component, y=value, colour = cell_type, group=id))+
  geom_line(linewidth = 0.025)+
  guides(colour=guide_legend(override.aes = list(linewidth=2)))+
  scale_color_brewer(palette = "Set1")+
  theme_bw(16)
ggsave("pcp_sn10.pdf", width = 7, height = 5)
