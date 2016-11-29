
mcmcLinreg = function(y, X, prior,mcmc){

   k=ncol(X)
   N=nrow(as.matrix(y))
   XX = crossprod(X)
   Xy = crossprod(X,y)

   sn = prior$s0 + N/2
   tau2= rep(prior$tau2_fix,k)

   if(length(tau2)>1){
    B0_inv = diag(1/tau2)
   } else {
     B0_inv = 1/tau2
   }
   cholx = chol(XX + B0_inv)
   BN =  backsolve(cholx, backsolve(cholx, diag(ncol(cholx)),transpose=TRUE))
   bN = BN %*% Xy
  sgma2=drop(stats::var(y-X%*%bN))

   result = list(beta=array(0,dim=c(mcmc$M,k)),
                  sgma2=rep(0,mcmc$M))

   result$mcmc=mcmc
   result$prior=prior

   #-------------------MCMC sampler-------------------------------------------#

   for (m in 1:mcmc$M){

       #------ step 1: sample the regression coefficients beta
       if (prior$conj){
        beta=MASS::mvrnorm(1,bN,BN*sgma2)
        }else{
          cholx = chol(XX/sgma2 + B0_inv)
          BN =  backsolve(cholx, backsolve(cholx, diag(ncol(cholx)),transpose=TRUE))
          bN = BN %*% Xy/sgma2
          beta=MASS::mvrnorm(1,bN,BN)
        }

        result$beta[m,]=beta

       #----- step 2: sample the error variance
       Sn = prior$S0 + 1/2 * t(y-X%*%beta)%*%(y-X%*%beta)
       sgma2= 1/stats::rgamma(1,sn,Sn)

       result$sgma2[m]=sgma2

    }

return(result)

}
