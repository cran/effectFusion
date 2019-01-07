#' @title Selected model of a \code{fusion} object
#' @description The function displays for categorical covariates the selected model of an
#' object of class \code{fusion} as list.
#'
#' @param x an object of class \code{fusion}
#'
#' @details The selected model for each categorical predictor is displayed as a list of length equal
#' to the number of categories after fusion. Fused categories are shown with their original labelling 
#' in one list element. 
#' The function is only available if effect fusion (method in \code{effectFusion} is unequal to \code{NULL})
#' and final model selection (argument \code{modelSelection} in \code{effectFusion} is not NULL) is performed.
#'
#' See \code{summary.fusion} for more details.
#' 
#' @author Daniela Pauger, Magdalena Leitner <magdalena.leitner@jku.at>
#' @seealso \code{\link{effectFusion}}
#' @export
#' 
#' @examples 
#' ## see example for effectFusion
#' 


model <- function(x) {
    
    stopifnot(is(x, "fusion"))
    
    if (is.null(x$modelSelection)) {
        cat("No model selection performed")
    } else {
        
        varCat <- x$model$n_ord + x$model$n_nom
        sel_mod <- x$refit$model
        ind <- c(0, cumsum(c(rep(1, x$model$n_cont), x$model$diff))) + 1
        cat <- x$model$categories
        
        for (i in 1:varCat) {
            k <- i + x$model$n_cont
            if (x$data$types[k] == "o") 
                S <- getSOrdinal(sel_mod[ind[k]:(ind[k + 1] - 1)])
            if (x$data$types[k] == "n") 
                S <- getSNominal(cat[i], sel_mod[ind[k]:(ind[k + 1] - 1)])
            
            show_model <- list()
            if (sum(rowSums(S) == 0) > 0) {
                show_model[[1]] <- c(0, which(rowSums(S) == 0))
            } else show_model[[1]] <- 0
            
            if (length(show_model[[1]]) != cat[i]) {
                if (sum(S[, 1]) == 0) {
                  for (j in 2:ncol(S)) {
                    show_model[[j]] <- which(S[, j] == 1)
                  }
                } else {
                  for (j in 1:ncol(S)) {
                    show_model[[j + 1]] <- which(S[, j] == 1)
                  }
                }
            }
            
            show_model <- lapply(show_model, function(x) {x + 1})
            show_model <- lapply(show_model, function(z, x, i) {x$data$levelnames[[i]][z]}, 
                   x = x, i = i)
            
            cat("Covariate '", names(x$data$levelnames)[i], "'", sep = "", "\n")
            print(show_model)
            
        }
        
    }
}
