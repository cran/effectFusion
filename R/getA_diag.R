
getADiag <- function(model) {
    
    cat_nom <- model$categories[(model$n_ord + 1):(model$n_ord + model$n_nom)]
    A <- getA(cat_nom - 1)
    A_diag <- as.matrix(Matrix::bdiag(A))
    
    if (2 %in% cat_nom) {
        A_diag <- A_diag[-which(rowSums(A_diag) == 0), , drop = FALSE]
    }
    
    n_row <- nrow(A_diag)
    
    if (model$n_ord != 0) {
        col_ord <- sum(model$diff[1:model$n_ord])
        A_diag <- cbind(matrix(0, nrow(A_diag), col_ord), A_diag)
    }
    
    # add column for intercept and continuous predictors
    if (model$n_cont > 0) {
        A_diag <- cbind(matrix(0, nrow(A_diag), model$n_cont + 1), A_diag)
    } else {
        A_diag <- cbind(matrix(0, nrow(A_diag), 1), A_diag)
    }
    
    if (n_row == 0) {
        A_diag <- c()
    }
    
    return(A_diag)
}
