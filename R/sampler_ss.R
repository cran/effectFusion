
mcmcSs = function(y, X, model,prior=list(),mcmc,mats){

    defaultPrior = list(r=50000,g0=5,G0=25,s0=0,S0=0,tau2_fix=NULL)
    prior = utils::modifyList(defaultPrior, as.list(prior))

    N=nrow(as.matrix(y))
    yw = y
    Xw = X

    categories = model$categories
    jkcov= length(model$cov0)
    jk=jkcov+1
    jk_beta = 1+sum(categories-1)+model$n_cont
    nVar = sum(model$n_nom,model$n_ord,model$n_cont)

    ind=cumsum(c(1,rep(1,model$n_cont),categories-1,1))

    G=createSelmat(model$diff)
    G_beta=createSelmat(categories-1)

    if(model$n_cont>0){
      for(i in 1:model$n_cont){
        G=Matrix::bdiag(1,G)
        G_beta=Matrix::bdiag(1,G_beta)
      }
      G=as.matrix(G)
      G_beta=as.matrix(G_beta)
    }

    trG=t(G)
    trG_beta=t(G_beta)
    D = as.matrix(Matrix::bdiag(1,mats$D_comb))
    B = as.matrix(Matrix::bdiag(1,mats$B))
    E = as.matrix(Matrix::bdiag(1,mats$D_dummy))
    A_diag = model$A
    XX = crossprod(Xw)
    Xy = crossprod(Xw,yw)
    sn = prior$s0 + N/2
    qr=prior$r-1

    gamma = getGamma(model)$gamma
    gamma_long = getGamma(model)$gamma_long
    gn=prior$g0+(c(rep(1,model$n_cont),categories-1))/2

    r_delta = rep(1,jkcov)
    delta =rep(1,jkcov)

    if (is.null(prior$tau2_fix)){
      tau2 = c(10000, rep(1,nVar))
    } else {
      tau2= c(10000,prior$tau2_fix)
    }

    tau2_lfd=c(tau2[1],t(trG_beta%*%tau2[-1]))
    tau2_delta =  c(tau2[1],t(trG%*%tau2[-1]))

    B0_res = correctRest(diag(jk),A_diag)
    B0_beta = D%*%B0_res%*%t(D)
    B0 = B0_beta * gamma * tau2_lfd
    cholx = chol(XX + solve(B0))
    BN =  backsolve(cholx, backsolve(cholx, diag(ncol(cholx)),transpose=TRUE))
    beta_0 = BN %*% Xy
    sgma2=drop(stats::var(yw-Xw%*%beta_0))
    theta_0=B%*%beta_0

    burnin=mcmc$burnin
    M=mcmc$M

    result = list(beta=array(0,dim=c(M,jk_beta)),
                  delta=array(0,dim=c(M,jkcov)),
                  sgma2=rep(0,M),
                  tau2=matrix(0,M,nVar+1),  # incl intercept
                  prior=prior)

   #-------------------MCMC sampler-------------------------------------------#

    m=1
    while(m<=M|m<=mcmc$burnin){

      if(m==mcmc$burnin){m=1;mcmc$burnin=0}
      if(m%%500==0) cat("spike-slab: ",m,"\n")

       #------ step 1: sample the regression coefficients theta

        B0_res = correctRest(diag(c(1,r_delta)),A_diag)
        B0_beta = E %*% (D%*%B0_res%*%t(D)) %*% t(E)
        B0 = B0_beta * gamma * tau2_lfd

        cholx = chol(XX/sgma2 + solve(B0))
        BN =  backsolve(cholx, backsolve(cholx, diag(ncol(cholx)),transpose=TRUE))
        bN = BN %*% (Xy /sgma2)

        beta=MASS::mvrnorm(1,bN,BN)
        theta=B%*%beta

        result$beta[m,]=beta

       #----- step 2: sample the error variance

       Sn = prior$S0 + 1/2 * t(yw-Xw%*%beta)%*%(yw-Xw%*%beta)
       sgma2= 1/stats::rgamma(1,sn,Sn)

       result$sgma2[m]=sgma2


       #----- step 3: sample the scales tau

       if (is.null(prior$tau2_fix)){
         Qh=c()
         B0_work=B0_beta * gamma
         for(i in 2:(nVar+1)){
           b_i = ind[i]:(ind[i+1]-1)
           Qh = c(Qh,t(beta[b_i]) %*% solve(B0_work[b_i,b_i]) %*% beta[b_i])
         }
         t_cov = 1/stats::rgamma(nVar,gn,prior$G0+0.5*Qh)
         tau2=c(1/stats::rgamma(1,prior$g0+1/2,prior$G0+0.5*theta[1]^2),t_cov)
         tau2_lfd=c(tau2[1],t(trG_beta%*%tau2[-1]))
         tau2_delta =  c(tau2[1],t(trG%*%tau2[-1]))
       }
       result$tau2[m,]=tau2

      # ---- step 4: sample the indicator variable delta

       if(m > mcmc$startsel){
         vv = 2*tau2_delta[-1]*gamma_long[-1]
         L= sqrt(prior$r) * exp(-qr * theta[-1]^2 / vv)
         post_delta = 1 / (1 + L)
         delta = stats::rbinom(jkcov,1,post_delta)
         r_delta= delta + (1-delta)/prior$r
       }

       result$delta[m,]=delta

      m=m+1

    }

    return(result)

}
