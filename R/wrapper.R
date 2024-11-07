VBPCA <- function(Zini, Wini, Y, rank, iter=100, prior_prec=1, a = 1, b = 1){
  if(!any(class(Y)=="dgTMatrix")){
    Y <- as(Y, "TsparseMatrix")
  }
  doVB_norm(Zini, Wini, 
            y = Y@x, rowi = Y@i, coli = Y@j, 
            Nr = Y@Dim[1], Nc = Y@Dim[2], 
            L = rank,
            iter = iter,
            prior_prec = prior_prec, a=a, b=b)
}
