library(Matrix)
library(VBspPCA)
library(ggplot2)

miris <- t(log1p(as.matrix(iris[,-5])))
miris <- as(miris, "TsparseMatrix")
writeMM(miris, "iris.mtx")
datpath = "iris.mtx"

system.time({
  out_an <- VBspPCA:::VBPCA_diag_mtx(datpath,  constr_type = "AN",
                                     lambda = 10,
                                     rank = 2, iter = 100,
                                     tau = 1, a = 1, b = 1)
})


system.time({
  out_sn <- VBspPCA:::VBPCA_diag_mtx(datpath,  constr_type = "SN",
                                 rank = 2, iter = 100,
                                 tau = 1, a = 1, b = 1)
})
 
system.time({
  out_nn <- VBspPCA:::VBPCA_diag_mtx(datpath,  constr_type = "NN",
                                     rank = 2, iter = 100,
                                     tau = 1, a = 1, b = 1)
})

ggplot(data = NULL)+
  geom_abline(slope = 1,intercept = 0, colour="lightgrey")+
  geom_point(aes(x=c(fit_pca(out_an)), y=c(as.matrix(miris)), colour="AN"), alpha=0.25)+
  geom_point(aes(x=c(fit_pca(out_sn)), y=c(as.matrix(miris)), colour="SN"), alpha=0.25)+
  geom_point(aes(x=c(fit_pca(out_nn)), y=c(as.matrix(miris)), colour="NN"), alpha=0.25)+
  guides(colour=guide_legend(override.aes = list(alpha=1)))+
  theme_bw()


matplot(cbind(out_an$logprob[-1],
              out_sn$logprob[-1],
              out_nn$logprob[-1]), type="l")


#####


system.time({
  out_an <- VBspPCA:::VBPCA_diag(miris, 
                                 constr_type = "AN", 
                                 rank = 2, iter = 100,
                                 tau = 1, a = 1, b = 1)
})

system.time({
  out_sn <- VBspPCA:::VBPCA_diag(miris,  constr_type = "SN", 
                                 rank = 2, iter = 100,
                                 tau = 1, a = 1, b = 1)
})

system.time({
  out_nn <- VBspPCA:::VBPCA_diag(miris,  constr_type = "NN", 
                                 rank = 2, iter = 100,
                                 tau = 1, a = 1, b = 1)
})

#round(c(out_sn$obs_prec,out_nn$obs_prec,out_an$obs_prec),2)
#[1] 96.26 94.63 92.58



ggplot(data=NULL)+
  geom_line(aes(x=2:100, y=out_an$logprob[-1], colour="AN", linetype="AN"))+
  geom_line(aes(x=2:100, y=out_sn$logprob[-1], colour="SN", linetype="SN"))+
  geom_line(aes(x=2:100, y=out_nn$logprob[-1], colour="NN", linetype="NN"))+
  labs(colour="method", linetype="method")+
  theme_bw()

col3 <- hcl.colors(3, palette = "Set 2", alpha = 0.9)
plot(out_an$mean_col, col=col3[iris$Species], pch=16)
plot(out_sn$mean_col, col=col3[iris$Species], pch=16)
plot(out_nn$mean_col, col=col3[iris$Species], pch=16)


ggplot(data = NULL)+
  geom_abline(slope = 1,intercept = 0, colour="lightgrey")+
  geom_point(aes(x=c(fit_pca(out_an)), y=c(as.matrix(miris)), colour="AN"), alpha=0.25)+
  geom_point(aes(x=c(fit_pca(out_sn)), y=c(as.matrix(miris)), colour="SN"), alpha=0.25)+
  geom_point(aes(x=c(fit_pca(out_nn)), y=c(as.matrix(miris)), colour="NN"), alpha=0.25)+
  guides(colour=guide_legend(override.aes = list(alpha=1)))+
  theme_bw()
