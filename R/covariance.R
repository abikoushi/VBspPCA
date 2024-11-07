scan_n_mtx <- function(con, n, skip = 0){
  base::scan(con, nmax = n, quiet=TRUE, 
             what = list(i=integer(), j=integer(), v=numeric()),
             skip = skip)
}


dataloader_mtx2 <- function(file_path, bag){
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