
createModelvars <- function(data) {
    
    N <- length(data$y)
    ind_ord <- data$ind_ord
    ind_nom <- data$ind_nom
    ind_cont <- data$ind_cont
    n_ord <- sum(ind_ord)
    n_nom <- sum(ind_nom)
    n_cont <- sum(ind_cont)
    
    if (n_ord != 0) {
        U <- as.matrix(data$X[, ind_ord])
        cat_ord <- apply(U, 2, max)
        diff <- cat_ord - 1
        col_ord <- sum(diff)
        cov0 <- rep(1, col_ord)
        
    } else {
        diff <- c()
        cat_ord <- c()
        col_ord <- 0
        cov0 <- c()
    }
    
    data$X <- as.matrix(data$X)
    if (n_nom != 0) {
        V <- as.matrix(data$X[, ind_nom])
        cat_nom <- apply(V, 2, max)
        for (i in 1:n_nom) {
            cov0 <- c(cov0, rep(1, cat_nom[i] - 1), rep(0, choose(cat_nom[i] - 1, 2)))
        }
        diff <- c(diff, choose(cat_nom, 2))
    } else {
        cat_nom <- c()
    }
    categories <- c(cat_ord, cat_nom)
    
    if (n_cont > 0) {
        X_dummy <- cbind(matrix(1, N, 1), data$X[, ind_cont], coding(data$X[, !ind_cont], constant = FALSE, 
            splitcod = FALSE))
        cov0 <- c(rep(1, n_cont), cov0)
    } else {
        X_dummy <- coding(data$X, constant = TRUE, splitcod = FALSE)
    }
    
    
    return(list(X_dummy = X_dummy, diff = diff, cov0 = cov0, categories = categories))
    
}
