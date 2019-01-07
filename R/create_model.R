
createModel <- function(data, mvars) {
    
    model <- list(categories = mvars$categories, diff = mvars$diff, n_cont = sum(data$ind_cont),  
        n_ord = sum(data$ind_ord),  n_nom = sum(data$ind_nom), cov0 = mvars$cov0, lNom = sum(data$ind_ord) + 
            1)
    
    if (model$n_nom != 0) {
        A_diag <- getADiag(model)
    } else A_diag <- c()
    model$A_diag <- A_diag
    
    return(model)
    
}
