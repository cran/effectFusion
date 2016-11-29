
selectModel = function(incl_prob, model, strategy, PIP){
  
  diff=model$diff_dummy
  nVar=length(model$categories)
  
if(strategy=="pip"){
  if(missing(PIP)) PIP=0.5
  sel_mod=1-fusionVector(incl_prob,model,pip=PIP)
}
 
if(strategy=="pam"){
  
  Sim_vector=incl_prob
  categories = as.numeric(model$categories)
  
  ncomp=0.5*(categories^2-categories)
  up_comp=cumsum(ncomp)
  low_comp=c(1,up_comp+1)
  
  dissim=list()
  for(i in 1:model$n_nom){
    dissim[[i]]=1-Sim_vector[low_comp[i]:up_comp[i]]
  }
  
  K0_pam=rep(0,model$n_nom)
  asw_pam=rep(0,model$n_nom)
  ass_pam=list()
  asw_list=list()
  classError_pam=rep(0,model$n_nom)
  adRI_pam=rep(0,model$n_nom)
  
  if(1 %in% model$diff){
    lNom=sum(model$diff==1)+1
    model_binary = as.numeric(incl_prob[1:(lNom-1)] > 0.5)
  } else {
    lNom=1
    model_binary=c()
  }
  
  for(i in lNom:model$n_nom){
    k_max=categories[i]-1
    asw <- numeric(k_max-1)
    clus_matrix=matrix(-9,k_max,categories[i])
    for (k in 2:k_max){
      pam_k=cluster::pam(x=dissim[[i]],diss=TRUE, k)
      asw[k] <- pam_k$silinfo$avg.width
      clus_matrix[k,]=pam_k$clustering
    }
    asw=asw[-1]
    clus_matrix=clus_matrix[-1,,drop=FALSE]
    k.best <- which.max(asw)
    K0_pam[i]=k.best
    asw_pam[i]=asw[which.max(asw)]
    ass_pam[[i]]=clus_matrix[k.best,]
    asw_list[[i]]=asw
  }
  
  model_pam=c()
  for(j in lNom:model$n_nom){
    Sv=ass_pam[[j]]
    x= outer(Sv, Sv, "==")
    model_pam=as.numeric(c(model_pam,t(x)[lower.tri(t(x), diag = FALSE)]))
  }
  
  sel_mod=c(model_binary,1-model_pam)
  
}  
   
return(sel_mod) 
  
}

 

