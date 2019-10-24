#' @title Print object of class \code{fusion}
#'
#' @description The default print method for a \code{fusion} object.
#'
#' @param x an object of class \code{fusion}
#' @param ... further arguments passed to or from other methods (not used)
#'
#' @details Returns basic information about the full, initial model: the number of 
#' observations, the family, the number and types of used covariates with their 
#' number of categories and MCMC options.
#' 
#' @author Daniela Pauger, Magdalena Leitner <effectfusion.jku@gmail.com>
#'
#' @seealso \code{\link{effectFusion}}
#' @method print fusion
#' @export
#' 
#' @examples
#' ## see example for effectFusion

print.fusion <- function(x, ...) {
    
    stopifnot(is(x, "fusion"))
    
    if (x$method == "SpikeSlab") 
        cat("Bayesian effect fusion with spike and slab prior:\n")
    if (x$method == "FinMix") 
        cat("Bayesian effect fusion with finite mixture prior:\n")
    if (x$method == "No effect fusion performed. Full model was estimated.") {
        cat("\nNo effect fusion performed. Full model was estimated.")
        cat("\nA flat, uninformative prior was used.\n")
    }
    
    cat("\nCall:\n")
    print(x$call)
    
    cat("\nModel:", length(x$data$y), "observations")
    cat("\n Family:", x$family)
    cat("\n Covariates:", ncol(x$data$X))
    cat("\n Types:")
    cat("\n\t", x$model$n_cont, "continuous")
    cat("\n\t", x$model$n_ord, "ordinal")
    cat("\n\t", x$model$n_nom, "nominal")
    cat("\n with number of categories:", x$model$categories)
    
    cat("\n\nMCMC:")
    cat("\nM =", x$mcmc$M, "draws after a burn-in of", x$mcmc$burnin)
    if (x$method == "SpikeSlab" | x$method == "FinMix") {
      cat("\nVariable selection started after", x$mcmc$startsel, "iterations")
    }
}
