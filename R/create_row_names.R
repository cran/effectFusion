
createRowNames=function(model){
  var=rep(1:(model$n_nom+model$n_ord),model$categories-1)+model$n_cont
  cat=sequence(model$categories-1)+1
  names=c()
  for(i in 1:length(cat)){
    names[i]=paste("var",var[i],".cat",cat[i],sep="")
  }
  if(model$n_cont>0){
    cont_names=c()
    for(i in 1:model$n_cont){
      cont_names[i]=paste("var",i,sep="")
    }
    names=c(cont_names,names)
  }
  names=c("(Intercept)",names)
  return(names)
}








