exchangable_cov <- function(main_var,re_var,block_sizes){
  accum <- c(0,cumsum(block_sizes))
  mat <- matrix(0,max(accum),max(accum))
  mat_inv <- matrix(0,max(accum),max(accum))
  diag(mat) <- main_var
  diag(mat_inv) <- 1/main_var
  for(i in seq_along(block_sizes)){
    if(block_sizes[i]==0) next
    mat[(accum[i]+1):(accum[i+1]),(accum[i]+1):(accum[i+1])] <-
      mat[(accum[i]+1):(accum[i+1]),(accum[i]+1):(accum[i+1])] + re_var
    mat_inv[(accum[i]+1):(accum[i+1]),(accum[i]+1):(accum[i+1])] <-
      mat_inv[(accum[i]+1):(accum[i+1]),(accum[i]+1):(accum[i+1])] -
      re_var/main_var/(main_var+re_var*block_sizes[i])
  }
  return(list(V=mat,Prec=mat_inv))
}