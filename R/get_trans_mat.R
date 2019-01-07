
getTransMat <- function(model) {
    TM <- diag(model$n_cont + 1)
    
    if (model$n_ord != 0) {
        
        Aord_list <- list()
        cat_ord <- model$categories[1:model$n_ord]
        
        for (j in 1:length(cat_ord)) {
            
            lj <- cat_ord[j] - 1
            A <- diag(lj)
            if (is.matrix(A[-lj, -1])) {
                diag(A[-lj, -1]) <- -1
            } else {
                A[-lj, -1] <- -1
            }
            
            Aord_list$A <- A
            names(Aord_list) <- 1:j
        }
        TM <- Matrix::bdiag(TM, as.matrix(Matrix::bdiag(Aord_list)))
    }
    
    
    if (model$n_nom != 0) {
        Anom_list <- list()
        cat_nom <- model$categories[(model$n_ord + 1):(model$n_ord + model$n_nom)]
        
        for (j in 1:length(cat_nom)) {
            
            k <- cat_nom[j] - 1
            if (k > 1) {
                
                hd <- choose(k, 2)
                Aresh <- matrix(0, hd, k)
                istrt <- 1
                
                for (i in 1:(k - 1)) {
                  hi <- k - i
                  iend <- istrt + hi - 1
                  indx <- istrt:iend
                  Aresh[indx, i] <- matrix(-1, hi, 1)
                  Aresh[indx, i + (1:hi)] <- diag(hi)
                  istrt <- iend + 1
                }
                A <- cbind(diag(rep(1, k)), t(Aresh))
            } else {
                A <- as.matrix(1)
            }
            
            Anom_list$A <- A
            names(Anom_list) <- 1:j
            
        }
        TM <- Matrix::bdiag(TM, as.matrix(Matrix::bdiag(Anom_list)))
    }
    
    
    TM <- as.matrix(TM)
    return(TM)
}
