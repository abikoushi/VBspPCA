library(Matrix)

miris <- as(as.matrix(iris[,-5]),"TsparseMatrix")

iter <- 10
rank <- 2
out <- VBspPCA:::doVB_norm(y = miris@x, rowi = miris@i, coli = miris@j, 
                 Nr = miris@Dim[1], Nc = miris@Dim[2], 
                 L = rank,
                 iter = iter,
                 prior_prec = 1,
                 a=1, b=1)
plot(out$lp)
plot(out$mean_z, col=iris$Species)

