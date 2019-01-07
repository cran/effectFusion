
getSOrdinal <- function(delta) {
    
    delta_gr <- cumsum(delta)
    ncat <- length(delta)
    
    #--- covariates with just 2 categories
    if (ncat == 1) {
        
        if (delta == 0) 
            S <- matrix(0, 1, 1)
        if (delta == 1) 
            S <- matrix(1, 1, 1)
        
    } else {
        
        #--- covariates with more than 2 categories
        
        # -- ALL categories should be fused (variable is excluded)
        if (sum(delta_gr) == 0) {
            S <- matrix(0, ncat, 1)
        } else {
            
            
            # starting with diagonal matrix
            S <- diag(ncat)
            # -- categories should be fused, if any delta=0
            if (0 %in% delta) {
                
                for (i in 1:ncat) {
                  S[i, i] <- 0
                  S[i, delta_gr[i]] <- 1
                }
            }
            
            
            # delete empty columns
            if (any(colSums(S) == 0)) 
                S <- S[, -which(colSums(S) == 0), drop = FALSE]
            
            
        }
    }
    
    S
    
}
