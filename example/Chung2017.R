library(ggplot2)
lr_default <- function(t, delay=1, forgetting=0.9){
  (t+delay)^(-forgetting)
}

ggplot(data=NULL)+
  stat_function(geom="point",fun=lr_default, args = list(forgetting=0.8, delay=1), n=5)+
  stat_function(geom="line",fun=lr_default, args = list(forgetting=0.8, delay=1), n=5)+
  xlim(0,100)


library(Matrix)
library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(VBspPCA)
library(umap)

cells <- read_csv("/Users/abeko/project/data/Data_Chung2017_Breast/Cells.csv")
path <- "/Users/abeko/project/data/Data_Chung2017_Breast/Exp_data_TPM.mtx"
mat <- readMM(path)
rm_wch <- which(apply(mat, MARGIN = 1, FUN = function(x)all(x==0)))

mat <- mat[-rm_wch, ]

#ca005 <- rgb(0,0,0,0.05)
res20_vb <- VBPCA(Y = log1p(mat), rank = 20, iter = 20, prior_prec = 0)
plot(res20_vb$logprob, type="l")

res20_svd <- svd(log1p(as.matrix(mat)), nu = 20, nv=20)

v <- apply(res20_vb$mean_col,2,var)
ord <- order(v, decreasing = TRUE)

df_vb <- mutate(data.frame(res20_vb$mean_col[,ord]),
             cell_type=cells$cell_type, id=1:n()) %>% 
  pivot_longer(X1:X20) %>% 
  mutate(name = factor(name, levels=paste0("X",1:20)))

p1 <- ggplot(df_vb, aes(x=name,y=value,colour=cell_type, group=id))+
  geom_line(alpha=0.1, position = position_jitter(height = 0,width = 0.3, seed = 123))+
  scale_colour_brewer(palette = "Set2")+
  guides(colour=guide_legend(override.aes = list(alpha=1, linewidth=2)))+
  theme_classic()


df <- data.frame(res20_svd$v) %>% 
  mutate(cell_type=cells$cell_type, id=1:n()) %>% 
  pivot_longer(X1:X20) %>% 
  mutate(name = factor(name, levels=paste0("X",1:20)))

print(p1)
p1 %+% df

ures_vb <- umap(res20_vb$mean_col)
ures_svd <- umap(res20_svd$v)
h5 <- hcl.colors(5, palette = "Set 2", alpha = 0.5)
plot(ures_vb$layout, col=h5[factor(cells$cell_type)], pch=16)
plot(ures_svd$layout, col=h5[factor(cells$cell_type)], pch=16)
