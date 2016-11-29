
getGamma = function(model){

  cat = model$categories
  
  vec_h = rep(cat,cat-1)/2
  vec_h_long = rep(cat,model$diff)/2
  if(model$n_ord>0){
    ind_h = cumsum(model$diff)[model$lNom-1]
    vec_h[1:ind_h]=1
    vec_h_long[1:ind_h]=1
  }
  vec_h = c(rep(1,model$n_cont+1),vec_h)
  vec_h_long = c(rep(1,model$n_cont+1),vec_h_long)

  return(list(gamma=vec_h,gamma_long=vec_h_long))
  
}  
  