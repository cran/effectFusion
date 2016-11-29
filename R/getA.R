
getA = function(levels){

A_list=list()
lev = levels
  
for(j in 1:length(lev)){
  
  if(lev[j]>1){
  
  hd=choose(lev[j],2)
  k=lev[j]
  Aresh=matrix(0,hd,k)
  istrt=1
  
  for(i in 1:(k-1)){
    hi = k-i
    iend = istrt+hi-1
    indx = istrt:iend
    Aresh[indx,i]=matrix(-1,hi,1)
    Aresh[indx,i+(1:hi)]=diag(hi)
    istrt=iend+1
  }
  
  if(hd==1){
    A=cbind(Aresh,1)
  } else {
    A=cbind(Aresh,diag(rep(1,hd)))
  }
  
  } else {
    A=as.matrix(0)
  }
  
  A_list$A = A
  names(A_list) = 1:j
  
}

A = A_list

}










