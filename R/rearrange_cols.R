reorder_cols = function(out, statfun = var){
  ord = order(apply(out$mean_col, 2, statfun), decreasing = TRUE)
  out$mean_row = out$mean_row[,ord]
  out$mean_col = out$mean_col[,ord]
  for(i in 1:2){
    out$eta[[i]] <- out$eta[[i]][,ord]  
  }
  out$H <- out$H[,ord]
  return(out)  
}
