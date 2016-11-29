

getSNominal = function(categories,sel_mod){  
  
  K = categories 
  delta = sel_mod
  ncat = choose(K,2)
  
  #--- covariates with just 2 categories
  if(ncat==1){
    
    if(delta==0) S=matrix(0,1,1)
    if(delta==1) S=matrix(1,1,1)
    
  } else {  
    
    # definition accoring to pdf on reduced and full model (28.11.2014)
     A1 = getA(K-1)[[1]]
     A_abs=abs(A1[,1:(K-1)])
     delta1 = delta[1:(K-1)]
     delta2 = delta[-c(1:(K-1))]
    
    # --- if all delta=1, no categories are fused
    if(sum(delta)==length(delta)){
      S=diag(delta1)
    } else {
      
      # construct matrix for refit
      S = diag(delta1)
      eta = delta1
      H_nom = sum(A_abs%*%delta1==2 & delta2==0)
      H_ind = which(A_abs%*%delta1==2 & delta2==0)
      H = A_abs[H_ind,,drop=FALSE]
      
      if(H_nom>0){
        for (h in 1:H_nom){
          k_h=which(H[h,]==1)[1]
          l=which(H[h,]==1)[2]
          S[k_h,l]=1
          S[l,l]=0
          h_tild = which(H[1:H_nom,l]==1)
          h_tild = h_tild[h_tild>h]
          H[h_tild,k_h]=1
          H[h_tild,l]=0
          eta[l]=0
        }
      }
      
      S = t(S[-which(eta==0),,drop=FALSE])
    }
  }
  
  
  return(S)
  
}

