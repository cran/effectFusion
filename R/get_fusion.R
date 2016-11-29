

getFusion = function(delta, cat, pip=0.5){
  
  Delta = diag(rep(NA,cat))
  Delta[lower.tri(Delta,diag=F)]=delta
  Delta=Delta+t(Delta)-diag(diag(Delta))  
  Delta = Delta > pip
  fusion=rep(0,length(delta))
  i=1
  
  for(r in 1:(cat-1)){
     for(c in (r+1):cat){
        if(Delta[r,c]|is.na(Delta[r,c])){
          if(!1 %in% apply(as.matrix(Delta[-c(r,c),c(r,c)]),1,sum)) fusion[i]=1
        }   
        i=i+1
      }
  }
  
  return(fusion)
  
}













