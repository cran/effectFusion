
correctRest = function(B0_unres,A_diag){
  
  if(is.null(A_diag)==FALSE){
    
    B0_w = A_diag%*%B0_unres%*%t(A_diag)
    B0_w2 = B0_unres%*%t(A_diag)
    B0_w3 = A_diag %*% B0_unres
    B0 = B0_unres -  B0_w2 %*% solve(B0_w) %*% B0_w3
    
  } else {
    B0 = B0_unres
  }
    
  return(B0)
  
  
}










