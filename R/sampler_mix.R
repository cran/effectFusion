
mcmcMix <- function(y, X, model, prior = list(), mcmc, returnBurnin) {
    
    defaultPrior <- list(e0 = 0.01, p = 100, hyperprior = FALSE, s0 = 0, S0 = 0)
    prior <- utils::modifyList(defaultPrior, as.list(prior))
    
    prior_par <- prior
    prior <- createPrior(y, X, model, prior_par$p, family = "gaussian")
    prior$e0 <- prior_par$e0
    prior$s0 <- prior_par$s0
    prior$S0 <- prior_par$S0
    
    N <- length(y)
    categories <- as.numeric(model$categories)  #vector with number of categories for each variable
    levels <- categories - 1
    n_nom <- model$n_nom
    n_con <- model$n_con
    r <- ncol(X)  #number of  all regressors (with intercept,cont.Var, dummy coded Var.)
    jk <- sum(levels)  #number of  effects (beta_jk)
    
    L <- categories - 1
    Lall <- sum(L)  # =jk (number of total mixture components(wihout 0-comp)=number of all levels)
    
    s0 <- prior$s0
    S0 <- prior$S0
    c0 <- prior$c0
    e0 <- rep(prior$e0, Lall)
    e0_j0 <- rep(prior$e0, length(L))
    mu_con <- prior$mu_con
    
    c0 <- prior$c0
    C_0j <- prior$C_0j
    cl <- c0 + 0.5 * L
    
    #--------------------------------
    # Initialize ....
    #---------------------------------
    
    if (prior_par$hyperprior) {
        N_jl <- rep(1, Lall)
        N_j0 <- rep(0, n_nom)
        Nh_j <- c()
        for (j in 1:n_nom) {
            Nh_j <- c(Nh_j, 0, rep(1, levels[j]))
        }
    }
    
    #---indices to address the beta_jk belonging to a specific covariate j:
    index <- c(0, cumsum(levels))
    low <- index + 1
    up <- index[-1]
    #---indices to address the mixtures components 0,1,...,Lj,
    # with 0-components included,i.e.1-33 for SILC):
    upc <- cumsum(categories)
    lowc <- cumsum(categories) - L
    
    #------- constant parameters for each iteration
    XX <- crossprod(X)
    Xy <- crossprod(X, y)
    
    sn <- s0 + N/2
    
    psi_inv <- prior$psi_inv
    psi_con_inv <- prior$psi_con_inv
    psi_inv_itcpt <- 0.01
    
    M0_inv <- prior$M0_inv * rep(1, Lall)
    Mm <- M0_inv * prior$mj0
    
    
    #---- initialize MCMC
    
    ###### initial value for beta:
    beta <- prior$beta_flat
    sgma2 <- prior$sgma2_flat
    
    ##### initial value for allocation vector S:
    S <- rep(0, jk)  #indicator vector for beta_jk
    
    ##### initial value for the vector reporting the number of effects allocated to a specific
    ##### component
    N_jl <- rep(1, Lall)
    N_j0 <- rep(0, n_nom)
    
    comp_means <- (beta)[(r - jk + 1):r]
    mu <- comp_means
    mean_beta_jl <- (beta)[(r - jk + 1):r]
    
    ###### initial value for comp_prec and B0_inv
    comp_prec <- rep(psi_inv, categories - 1)
    psi_vector <- rep(1/psi_inv, L + 1)
    
    eta0 <- rep(0, n_nom)  #to store the weights of the 0-components
    eta <- rep(0, Lall)  #to store the weights of the non-0 components
    
    M <- mcmc$M
    burnin <- mcmc$burnin
    M_thin <- M
    thinning <- 1
    
    if (sum(model$diff) > 999 & M > 15000) {
        thinning <- ceiling(M/15000)
        M_thin <- M/thinning
    }
    
    
    ##--- generate matrices for storing the draws
    result <- list(beta = matrix(NA, M_thin + burnin, r), S = matrix(NA, 
        M_thin + burnin, jk), eta = matrix(NA, M_thin + burnin, Lall), eta0 = matrix(NA, M_thin + 
        burnin, n_nom), mu = matrix(NA, M_thin + burnin, Lall), N_jl_matrix = matrix(NA, M_thin + 
        burnin, Lall), N_j0_matrix = matrix(NA, M_thin + burnin, n_nom), sgma2 = rep(NA, M_thin + burnin))
    
    
    #-------------------MCMC sampler-------------------------------------------#
    m_thin <- burnin + 1
    
    for (m in 1:(M + burnin)) {
        
        if (mcmc$burnin == 0 & m%%500 == 0) 
            cat("finite mixture: ", m - burnin, "\n")
        if (mcmc$burnin != 0 & m%%500 == 0) 
            cat("finite mixture: ", m, "\n")
        if (m == mcmc$burnin + 1) {
            mcmc$burnin <- 0
        }
        
        
        
        #------ step 1: sample the regression coefficients beta
        
        B0_inv <- diag(c(psi_inv_itcpt, rep(psi_con_inv, n_con), comp_prec))
        
        cholx <- chol(XX + sgma2 * B0_inv)
        b0 <- c(prior$beta_flat[1], rep(mu_con, n_con), comp_means)  ## component means
        BN <- sgma2 * backsolve(cholx, backsolve(cholx, diag(ncol(cholx)), transpose = TRUE))
        bN <- BN %*% (Xy/sgma2 + B0_inv %*% b0)
        
        beta <- as.vector(MASS::mvrnorm(1, bN, BN))
        if (thinning > 1 & mcmc$burnin == 0) {
            if (m%%thinning == 0) 
                result$beta[m_thin, ] <- beta
        } else {
            result$beta[m, ] <- beta
        }
        beta_nom <- beta[(r - jk + 1):length(beta)]
        
        
        #----- step 2: sample the error variance
        
        Sn <- S0 + 1/2 * t(y - X %*% beta) %*% (y - X %*% beta)
        sgma2 <- 1/stats::rgamma(1, sn, Sn)
        if (thinning > 1 & mcmc$burnin == 0 & (m%%thinning == 0)) {
            if (m%%thinning == 0) 
                result$sgma2[m_thin] <- sgma2
        } else {
            result$sgma2[m] <- sgma2
        }
        
        
        if (m > mcmc$startsel) {
            #----- step 3: sample the component weight
            
            e_jl <- e0 + N_jl  #without 0-component
            e_j0 <- e0_j0 + N_j0
            
            for (j in (1:n_nom)) {
                ind_j <- low[j]:up[j]
                e_j <- c(e_j0[j], e_jl[ind_j])
                etaj <- bayesm::rdirichlet(e_j)
                eta0[j] <- etaj[1]
                eta[ind_j] <- etaj[-1]
            }
            
            if (thinning > 1 & mcmc$burnin == 0) {
                if (m%%thinning == 0) {
                  result$eta[m_thin, ] <- eta
                  result$eta0[m_thin, ] <- eta0
                }
            } else {
                result$eta[m, ] <- eta
                result$eta0[m, ] <- eta0
            }
            
            #----- step 4: sample the mixture component means mu
            
            MN <- 1/(N_jl * comp_prec + M0_inv)
            mN <- MN * (mean_beta_jl * N_jl * comp_prec + Mm)
            mu <- stats::rnorm(Lall, mN, MN)
            if (thinning > 1 & mcmc$burnin == 0) {
                if (m%%thinning == 0) 
                  result$mu[m_thin, ] <- mu
            } else {
                result$mu[m, ] <- mu
            }
            
            
            #---- step 5a: sample the mixture component variances psi
            
            if (prior_par$hyperprior) {
                
                for (j in (1:n_nom)) {
                  ind_j <- low[j]:up[j]
                  indc_j <- lowc[j]:upc[j]
                  
                  ## sample ONE sigma2 for all components:
                  mu_j <- c(0, mu[ind_j])
                  W_l <- rep(0, L[j] + 1)  #for storing the within-component variance
                  
                  for (l in 1:(L[j] + 1)) {
                    # for each mixture component l, with 0-component
                    if (Nh_j[indc_j][l] > 0) {
                      W_l[l] <- sum((beta[ind_j][S[ind_j] == l] - mu_j[l])^2)
                    }
                  }
                  Cl_j <- C_0j[j] + 0.5 * sum(W_l)
                  
                  psi_vector[indc_j] <- 1/stats::rgamma(1, shape = cl[j], rate = Cl_j)
                }
            }
            
            
            #----- step 5: sample the allocation indicators and update parameters
            
            for (j in (1:n_nom)) {
                
                # For each varable j:
                
                Lj1 <- L[j] + 1  #number of mixture components (with 0-component included)
                Lj1_vec <- 1:Lj1
                ind_j <- low[j]:up[j]
                indc_j <- lowc[j]:upc[j]
                
                # read the actual beta_jk,component means,weights and sd
                betaj <- beta_nom[ind_j]
                eta_j <- c(eta0[j], eta[ind_j])
                mu_j <- c(0, mu[ind_j])
                sq_psi_j <- sqrt(psi_vector[indc_j])
                
                mat <- sapply(Lj1_vec, function(l) log(eta_j[l]) + stats::dnorm(betaj, mu_j[l], 
                  sq_psi_j[l], log = T))
                
                # Max_vector because of numerical issues necessary
                if (is.matrix(mat)) {
                  Max_vector <- apply(mat, 1, max)
                  S_j <- apply(exp(mat - Max_vector), 1, function(x) sample(1:Lj1, 1, prob = x))  # das sind die S-Werte!
                } else {
                  Max_vector <- max(mat)
                  S_j <- sample(Lj1_vec, 1, prob = exp(mat - Max_vector))
                }
                
                #---- Updating parameters
                
                Nh_jl <- tabulate(S_j, nbins = Lj1)
                mean.betah_j <- sapply(Lj1_vec, function(l) mean(betaj[S_j == l, drop = FALSE]))
                if (sum(is.na(mean.betah_j)) > 0) {
                  ## to catch the case if a group is empty: NA values are substituted by zeros
                  mean.betah_j[is.na(mean.betah_j)] <- 0
                }
                N_jl[ind_j] <- Nh_jl[Lj1_vec[-1]]
                N_j0[j] <- Nh_jl[1]
                
                mean_beta_jl[ind_j] <- mean.betah_j[Lj1_vec[-1]]
                
                S[ind_j] <- S_j
                comp_means[ind_j] <- mu_j[S_j]
                comp_prec[ind_j] <- 1/(psi_vector[indc_j])[S_j]
            }
            
            if (thinning > 1 & mcmc$burnin == 0) {
                if (m%%thinning == 0) {
                  result$S[m_thin, ] <- S[]
                  result$N_jl_matrix[m_thin, ] <- N_jl
                  result$N_j0_matrix[m_thin, ] <- N_j0
                  m_thin <- m_thin + 1
                }
            } else {
                result$S[m, ] <- S[]
                result$N_jl_matrix[m, ] <- N_jl
                result$N_j0_matrix[m, ] <- N_j0
            }
        }
        
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
    
    return(result[!names(result) %in% c("N_jl_matrix", "N_j0_matrix")])
}
