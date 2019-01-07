#' @importFrom mcclust minbinder

selectModel <- function(incl_prob, model, strategy) {
    
    diff <- model$diff_dummy
    nVar <- model$n_nom + model$n_ord + model$n_cont
    
    if (strategy == "spikeslab_binder") {
        fuseProb <- 1 - incl_prob
        type <- c(rep("c", model$n_cont), rep("o", model$n_ord), rep("n", model$n_nom))
        
        ind <- c(rep(1, model$n_cont), model$diff)
        ind <- c(0, cumsum(ind))
        
        for (i in 1:nVar) {
            if (type[i] == "c") {
                indi <- ind[i + 1]
            } else {
                indi <- (ind[i] + 1):(ind[i + 1])
            }
            
            fpi <- fuseProb[indi]
            
            if (type[i] == "c") {
                selMod_i <- (fpi < 0.5)
                S <- matrix(selMod_i, 1, 1)
            } else {
                if (type[i] == "o") {
                  selMod_i <- (fpi < 0.5)
                  S <- getSOrdinal(selMod_i)
                } else {
                  ncat <- model$categories[i - model$n_cont]
                  
                  if (ncat == 2) {
                    selMod_i <- (fpi < 0.5)
                    S <- getSNominal(2, selMod_i)
                  } else {
                    
                    mat_data <- diag(ncat)
                    mat_data[lower.tri(mat_data, diag = FALSE)] <- fpi
                    mat_data <- mat_data + t(mat_data) - diag(ncat)
                    
                    clusterRes <- mcclust::minbinder(psm = mat_data, method = "laugreen")
                    ncl <- max(clusterRes$cl)
                    cl <- clusterRes$cl
                    
                    # compute model
                    diff <- t(matrix(cl, ncol = ncat, nrow = ncat)) - matrix(cl, ncol = ncat, 
                      nrow = ncat)
                    h <- abs(diff) > 0
                    selMod_i <- as.vector(h[lower.tri(h, diag = FALSE)])
                    
                    matC <- matrix(rep(clusterRes$cl, ncl), ncol = ncl)
                    hC <- t(matrix(rep((1:ncl), ncat), nrow = ncl))
                    S <- as.matrix(((matC == hC) * 1)[-1, -1])
                  }
                }
            }
            if (i == 1) {
                sel_mod <- selMod_i
                S_M <- S
            } else {
                S_M <- Matrix::bdiag(S_M, S)
                sel_mod <- c(sel_mod, selMod_i)
            }
            
        }
    }
    
    if (strategy == "finmix_pam") {
        
        Sim_vector <- incl_prob
        categories <- as.numeric(model$categories)
        
        ncomp <- 0.5 * (categories^2 - categories)
        up_comp <- cumsum(ncomp)
        low_comp <- c(1, up_comp + 1)
        
        dissim <- list()
        for (i in 1:model$n_nom) {
            dissim[[i]] <- 1 - Sim_vector[low_comp[i]:up_comp[i]]
        }
        
        K0_pam <- rep(0, model$n_nom)
        asw_pam <- rep(0, model$n_nom)
        ass_pam <- list()
        gs_pam <- list()
        asw_list <- list()
        classError_pam <- rep(0, model$n_nom)
        adRI_pam <- rep(0, model$n_nom)
        
        Var <- 1:nVar
        ind <- c(rep(1, model$n_cont), model$diff)
        ind <- c(0, cumsum(ind))
        pam_model <- rep(NA, length(incl_prob))
        
        bin_vars <- c()
        if (1 %in% model$diff) {
            bin_vars <- Var[model$diff == 1]
            
            for (i in bin_vars) {
                indi <- ind[i + 1]
                pam_model[indi] <- as.numeric(incl_prob[indi] < 0.5)
            }
        }
        
        if (length(bin_vars) != 0) {
            not_bin_vars <- Var[-which(Var %in% bin_vars)]
        } else {
            not_bin_vars <- Var
        }
        
        
        for (i in not_bin_vars) {
            k_max <- categories[i] - 1
            asw <- numeric(k_max - 1)
            clus_matrix <- matrix(-9, k_max, categories[i])
            for (k in 2:k_max) {
                pam_k <- cluster::pam(x = dissim[[i]], diss = TRUE, k)
                asw[k] <- pam_k$silinfo$avg.width
                clus_matrix[k, ] <- pam_k$clustering
            }
            asw <- asw[-1]
            clus_matrix <- clus_matrix[-1, , drop = FALSE]
            
            k.best <- which(asw == max(asw))
            if (length(k.best) > 1) {
                clus_number <- apply(clus_matrix[k.best, ], 1, function(z) length(unique(z)))
                k.best <- k.best[which(clus_number == max(clus_number))[1]]
            }
            
            K0_pam[i] <- k.best
            asw_pam[i] <- asw[k.best]
            ass_pam[[i]] <- clus_matrix[k.best, ]
            asw_list[[i]] <- asw
        }
        
        for (j in not_bin_vars) {
            Sv <- ass_pam[[j]]
            x <- outer(Sv, Sv, "==")
            indi <- (ind[j] + 1):(ind[j + 1])
            pam_model[indi] <- 1 - as.numeric(t(x)[lower.tri(t(x), diag = FALSE)])
        }
        
        sel_mod <- pam_model
        
        
    }
    
    
    if (strategy == "finmix_binder") {
      type <- c(rep("c", model$n_cont), rep("o", model$n_ord), rep("n", model$n_nom))
      
      ind <- c(rep(1, model$n_cont), model$categories - 1)
      ind <- c(0, cumsum(ind))
      binder_model <- c()
      
      for (i in 1:nVar) {
        if (type[i] == "c") {
          indi <- ind[i + 1]
        } else {
          indi <- (ind[i] + 1):(ind[i + 1])
        }
        Sv <- GreedyEPL::MinimiseEPL(cbind(rep(1, nrow(incl_prob)), incl_prob[, indi]), 
                               list(Kup = model$categories[i], loss_type = "B"))$decision
        
        x <- outer(Sv, Sv, "==")
        indi <- (ind[i] + 1):(ind[i + 1])
        binder_model <- c(binder_model, 1 - as.numeric(t(x)[lower.tri(t(x), diag = FALSE)]))
      }
      
      sel_mod <- binder_model
    }
    
    return(as.numeric(sel_mod))
    
}
