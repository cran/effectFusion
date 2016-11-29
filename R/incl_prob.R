
inclProb = function(S,mvars,model){
  
  diff = as.numeric(model$diff)
  lowd=cumsum(diff)-diff+1
  upd=cumsum(diff)
  
  categories = as.numeric(model$categories)
  index=c(0,cumsum(c(categories-1)))
  low=index+1
  up=index[-1]
  
  jk=sum(diff)
  incl_prob=rep(0,jk)
  for (m in 1:length(S[,1])){
    delta=c()
    for(j in 1:model$n_nom){
      Sv=S[m,low[j]:up[j]]
      x= outer(Sv, Sv, "==")
      delta=c(delta,Sv==1,t(x)[lower.tri(t(x), diag = FALSE)])
    }
    incl_prob=incl_prob+delta/length(S[,1])
  } 
  
  return(incl_prob)
  
  
}













