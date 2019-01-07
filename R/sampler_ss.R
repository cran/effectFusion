
mcmcSs <- function(y, X, model, prior = list(), mcmc, mats, returnBurnin) {
    
    defaultPrior <- list(r = 50000, g0 = 5, G0 = 25, tau2_fix = NULL, s0 = 0, S0 = 0)
    prior <- utils::modifyList(defaultPrior, as.list(prior))
    
    N <- nrow(as.matrix(y))
    yw <- y
    Xw <- X
    
    categories <- model$categories
    jkcov <- length(model$cov0)
    jk <- jkcov + 1
    jk_beta <- 1 + sum(categories - 1) + model$n_cont
    nVar <- sum(model$n_nom, model$n_ord, model$n_cont)
    
    if (!is.null(prior$tau2_fix)) {
      if (length(prior$tau2_fix) != nVar) {
        stop("'tau2_fix' must be specified for every covariate separately")
      }
    }
    
    ind <- cumsum(c(1, rep(1, model$n_cont), categories - 1, 1))
    
    G <- createSelmat(model$diff)
    G_beta <- createSelmat(categories - 1)
    
    if (model$n_cont > 0) {
        for (i in 1:model$n_cont) {
            G <- Matrix::bdiag(1, G)
            G_beta <- Matrix::bdiag(1, G_beta)
        }
        G <- as.matrix(G)
        G_beta <- as.matrix(G_beta)
    }
    
    trG <- t(G)
    trG_beta <- t(G_beta)
    
    #D <- as.matrix(Matrix::bdiag(1, mats$D_comb))
    B <- as.matrix(Matrix::bdiag(1, mats$B))
    #E <- as.matrix(Matrix::bdiag(1, mats$D_dummy))
    #A_diag <- model$A
    
    TM <- getTransMat(model)
    XX <- crossprod(Xw)
    Xy <- crossprod(Xw, yw)
    sn <- prior$s0 + N/2
    qr <- prior$r - 1
    
    gamma <- getGamma(model)$gamma
    gamma_long <- getGamma(model)$gamma_long
    gn <- prior$g0 + (c(rep(1, model$n_cont), categories - 1))/2
    
    r_delta <- rep(1, jkcov)
    delta <- rep(1, jkcov)
    
    if (is.null(prior$tau2_fix)) {
        # starting values for sampling
        tau2 <- c(10^9, rep(1000, nVar))
    } else {
        # assign fix values
        tau2 <- c(10^9, prior$tau2_fix)
    }
    
    tau2_lfd <- c(tau2[1], t(trG_beta %*% tau2[-1]))
    tau2_delta <- c(tau2[1], t(trG %*% tau2[-1]))
    
    beta <- solve(XX + diag(1/tau2_lfd)) %*% Xy
    sgma2 <- drop(stats::var(yw - Xw %*% beta))
    #theta_0 <- B %*% beta
    
    burnin <- mcmc$burnin
    M <- mcmc$M
    
    result <- list(beta = array(0, dim = c(M + burnin, jk_beta)), delta = array(0, dim = c(M + 
        burnin, jkcov)), tau2 = matrix(0, M + burnin, nVar + 1), sgma2 = rep(0, M + burnin))  # incl intercept
    
    
    #-------------------MCMC sampler-------------------------------------------#
    
    for (m in 1:(M + burnin)) {
        
        if (mcmc$burnin == 0 & m%%500 == 0) 
            cat("spike-slab: ", m - burnin, "\n")
        if (mcmc$burnin != 0 & m%%500 == 0) 
            cat("spike-slab: ", m, "\n")
        if (m == mcmc$burnin + 1) {
            mcmc$burnin <- 0
        }
        
        #------ step 1: sample the regression coefficients beta
        
        B0_invh <- TM %*% diag(c(1, 1/r_delta)) %*% t(TM)/gamma
        B0_inv <- B0_invh/tau2_lfd
        
        BN <- solve(B0_inv + XX/sgma2)
        bN <- BN %*% (Xy/sgma2)
        
        beta <- MASS::mvrnorm(1, bN, BN)
        theta <- B %*% beta
        
        result$beta[m, ] <- beta
        
        #----- step 2: sample the error variance
        
        Sn <- prior$S0 + 1/2 * t(yw - Xw %*% beta) %*% (yw - Xw %*% beta)
        sgma2 <- 1/stats::rgamma(1, sn, Sn)
        
        result$sgma2[m] <- sgma2
        
        
        #----- step 3: sample the scales tau
        
        if (is.null(prior$tau2_fix)) {
            Qh <- c()
            
            for (i in 1:nVar) {
                b_i <- (ind[i] + 1):(ind[i + 1])
                Qh <- c(Qh, t(beta[b_i]) %*% (B0_invh[b_i, b_i] %*% beta[b_i]))
            }
            
            t_cov <- 1/stats::rgamma(nVar, gn, prior$G0 + 0.5 * Qh)
            tau2 <- c(1/stats::rgamma(1, prior$g0 + 1/2, prior$G0 + 0.5 * theta[1]^2), t_cov)
            tau2_lfd <- c(tau2[1], t(trG_beta %*% tau2[-1]))
            tau2_delta <- c(tau2[1], t(trG %*% tau2[-1]))
        }
        result$tau2[m, ] <- tau2
        
        # ---- step 4: sample the indicator variable delta
        
        if (m > mcmc$startsel) {
            vv <- 2 * tau2_delta[-1] * gamma_long[-1]
            L <- sqrt(prior$r) * exp(-qr * theta[-1]^2/vv)
            post_delta <- 1/(1 + L)
            delta <- stats::rbinom(jkcov, 1, post_delta)
            r_delta <- delta + (1 - delta)/prior$r
        }
        
        result$delta[m, ] <- delta
        
    }
    
    if (!returnBurnin) {
        result <- lapply(result, function(x, burnin) {
            if (is.matrix(x)) 
                return(x[-(1:burnin), ])
            if (is.vector(x)) 
                return(x[-(1:burnin)])
        }, burnin = burnin)
    }
    result[["prior"]] <- prior
    
    return(result)
    
}
