
createPrior <- function(y, X, model, p, family) {
    
    if (family == "gaussian") {
        mcmc_flat <- list(M = 3000, burnin = 1000)
        res_flat <- mcmcLinreg(y, X, prior = list(s0 = 0, S0 = 0, tau2_fix = 1000, conj = FALSE), 
            M = mcmc_flat$M, burnin = mcmc_flat$burnin, returnBurnin = FALSE)
        sgma2_flat <- mean(res_flat$sgma2)
        
    }
    if (family == "binomial") {
        res_flat <- logit(y, X, samp = 3000, burn = 1000, P0 = diag(0.1, nrow = ncol(X), ncol = ncol(X)))
        sgma2_flat <- NULL
    }
    
    beta_nom_bay <- colMeans(res_flat$beta)[-(1:(model$n_con + 1))]
    beta_flat <- colMeans(res_flat$beta)
    
    categories <- as.numeric(model$categories)
    
    index <- c(0, cumsum(c(categories - 1)))
    low <- index + 1
    up <- index[-1]
    
    
    mj0 <- rep(0, sum(categories - 1))
    M0 <- rep(0, sum(categories - 1))
    cov_j <- rep(0, length(categories))
    
    for (i in 1:model$n_nom) {
        betaH <- beta_nom_bay[low[i]:up[i]]
        mj0[low[i]:up[i]] <- mean(betaH)
        if (length(low[i]:up[i]) > 1) {
            M0[low[i]:up[i]] <- (diff(range(betaH)))^2
            cov_j[i] <- stats::var(betaH)
        } else {
            M0[low[i]:up[i]] <- (betaH)^2
            cov_j[i] <- stats::var(c(betaH, 0))
        }
    }
    
    M0 <- 10 * M0
    psi_j <- cov_j/p
    
    c0 <- 100
    C_0j <- (c0 - 1) * psi_j
    
    
    return(list(psi_con_inv = 0.01, mu_con = 0, psi_inv = 1/psi_j, M0_inv = 1/M0, mj0 = mj0, 
        comp_means_ini = beta_nom_bay, c0 = c0, C_0j = C_0j, p = p, beta_flat = beta_flat, sgma2_flat = sgma2_flat))
    
}
