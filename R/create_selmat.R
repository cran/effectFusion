
createSelmat <- function(cats) {
    nvar <- length(cats)
    ncol <- sum(cats)
    
    G <- matrix(0, nvar, ncol)
    ind_last <- 0
    
    for (i in 1:nvar) {
        icat <- cats[i]
        ind <- ind_last + (1:icat)
        G[i, ind] <- rep(1, icat)
        ind_last <- ind[icat]
    }
    return(G)
}
