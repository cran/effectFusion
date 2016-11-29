
getReparmats = function(model){
  
  
  n_ord = model$n_ord
  n_nom = model$n_nom
  n_cont = model$n_cont
  diff = model$diff
  cat=model$categories
  
  
  if(n_ord!=0){ 
  
  d=diff[1]
  D=matrix(0,d,d)
  ind=lower.tri(matrix(1,d,d),diag=TRUE)
  D[ind]=1
  
  if(d==1) B=matrix(1,1,1) else B=Matrix::bandSparse(d,k=-1:0,diag=list(rep(-1,d), rep(1,d)))
   
  if(n_ord>1){
  
    for (i in 2:n_ord){
      d=diff[i]
      if(d==1){
        Bh=1
        Dh=1
      } else {
      Dh=matrix(0,d,d)
      ind=lower.tri(matrix(1,d,d),diag=TRUE)
      Dh[ind]=1
    
      Bh=Matrix::bandSparse(d,k=-1:0,diag=list(rep(-1,d), rep(1,d)))
      }
      B=Matrix::bdiag(B,Bh)
      D=Matrix::bdiag(D,Dh)
    }
  }
  
  
  } 
   
  if(n_nom!=0){
    
    for (i in 1:n_nom){
      hd=diff[i+model$lNom-1]    
    
      k=cat[n_ord+i]-1
      q=cat[n_ord+i]
      # note: first k elements of beta are also first elements of theta
      # q are elements for effect coding
    
      if(k>1){
      
      Bh=matrix(0,hd-k,k)
      istrt=1    
    
      for(j in 1:(k-1)){
        hi = k-j
        iend = istrt+hi-1
        indx = istrt:iend
      
        Bh[indx,j]=matrix(1,hi,1)
        Bh[indx,(j+1):k]=-diag(hi)
        istrt=iend+1
      }
      Bh=rbind(diag(k),Bh)
      
      } else {
        Bh = as.matrix(1)
      }
      
      if(i==1 && n_ord==0){
        B=Bh
      } else {
        B=Matrix::bdiag(B,Bh)    
      }
      
  
      if(i==1 && n_ord==0){
        D=cbind(diag(k),matrix(0,k,hd-k))
      } else {  
        D=Matrix::bdiag(D,cbind(diag(k),matrix(0,k,hd-k)))
      }
     
      
    }
  }
  
  
  ############################################
  
  # Reparameterization Matrix for prior version
  # incorporating linear restriction
  
  if(n_ord>0){
    
    lOrd=sum(diff[1:n_ord])
    D_comb = as.matrix(D)
    D_comb[1:lOrd,1:lOrd]=diag(rep(1,lOrd))
    
    B_comb = as.matrix(B)
    B_comb[1:lOrd,1:lOrd]=diag(rep(1,lOrd))
    
    
  } else {
    D_comb = as.matrix(D)
    B_comb = as.matrix(B)
  }
  
  
  # Reparameterization Matrix for sampling 
  # dummy-coded ordinal effects
  D=as.matrix(D)
  if(sum(colSums(D)==0)>0) D_dummy = D[,-which(colSums(D)==0)] else D_dummy = D
  
  
  #####################################
  
  # including continuous predictors
  
  if(n_cont>0){
    for(i in 1:n_cont){
      D = Matrix::bdiag(1,D)
      B = Matrix::bdiag(1,B)
      D_comb = Matrix::bdiag(1,D_comb)
      B_comb = Matrix::bdiag(1,B_comb)
      D_dummy = Matrix::bdiag(1,D_dummy)
    }
    
  }
  
  
  
  return(list(D_comb=D_comb,B_comb=B_comb,D=as.matrix(D),B=as.matrix(B),D_dummy=as.matrix(D_dummy)))
}

