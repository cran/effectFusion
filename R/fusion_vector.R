
fusionVector = function(delta,model,pip=0.5){
  
  delta=1-delta
  type=c(rep("c",model$n_cont),rep("o",model$n_ord),rep("n",model$n_nom))
  nVar=length(type)
  ind=c(1,rep(1,model$n_cont),cumsum(model$diff)+1)  
  fusion=c()
  
  for(i in 1:nVar){
    
    indi=ind[i]:(ind[i+1]-1)
    d=delta[indi]
    
    if(type[i]!="c"){
      fusion=c(fusion,as.numeric(d>pip))
    } else {
      fusion=c(fusion,getFusion(d,model$categories[i],pip))
    }  
  }
  
  return(fusion)
  
}








