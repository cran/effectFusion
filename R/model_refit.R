
modelRefit = function(model,sel_mod,data,prior_flat,mcmc_flat){
  
  lprint_true=TRUE 
  A=model$A_diag
  
  n_cont=model$n_cont
  nVar=sum(model$n_nom,model$n_ord,n_cont)

  if(model$n_cont>0){
    ind=c(0,cumsum(c(rep(1,n_cont),model$diff)))+1
    cat=c(rep(2,n_cont),model$categories)
  } else {
    ind=c(0,cumsum(model$diff))+1
    cat=model$categories
  }
  
  for(k in 1:nVar){
    
    if(data$ind_ord[k]) S=getSOrdinal(sel_mod[ind[k]:(ind[k+1]-1)])
    if(data$ind_nom[k])S=getSNominal(cat[k],sel_mod[ind[k]:(ind[k+1]-1)])
    if(data$ind_cont[k]) S=matrix(1,1,1)

    if(k==1){
      S_M=S
    } else {
      S_M = Matrix::bdiag(S_M,S)
    }
    
    if (lprint_true){
    show_model=list()
    if(sum(rowSums(S)==0)>0){
      show_model[[1]]=c(0,which(rowSums(S)==0))
    } else show_model[[1]]=0
    
    if(length(show_model[[1]])!=cat[k]){
      if(sum(S[,1])==0){
        for(j in 2:ncol(S)){
          show_model[[j]]=which(S[,j]==1)
        }
      } else {
        for(j in 1:ncol(S)){
          show_model[[j+1]] = which(S[,j]==1)
        }
      }
    }
    }
  }

  S_M = as.matrix(S_M)
  if(0%in%colSums(S_M)){
    S_M = S_M[,-which(colSums(S_M)==0),drop=FALSE]
  }

  orig_dummy=createModelvars(data)$X_dummy[,-1]
  X_dummy = cbind(rep(1,length(data$y)),orig_dummy %*% S_M)
  
  res=mcmcLinreg(data$y, X_dummy, prior_flat,mcmc_flat)
  betaM=as.matrix(res$beta[mcmc_flat$burnin1:mcmc_flat$M,])
  if(ncol(betaM)==1){
    beta=cbind(betaM[,1],matrix(0,2000,nrow(S_M)))
  } else {
    beta=cbind(betaM[,1],betaM[,-1]%*%t(S_M))
  }

  sgma2=res$sgma2[mcmc_flat$burnin1:mcmc_flat$M]
  
  res_parm=list(beta=as.matrix(beta),sgma2=sgma2)
  return(res_parm)
  
}
