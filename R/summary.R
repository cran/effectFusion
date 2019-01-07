#' @title Summary of object of class \code{fusion}
#' 
#' @description Returns basic information about the model and the priors, MCMC details and posterior
#' means from the refit of the selected model or model averaged results as well as 95\%-HPD intervals for the regression effects.
#'
#' @param object an object of class \code{fusion}
#' @param ... further arguments passed to or from other methods (not used)
#'
#' @details The model selected with function \code{effectFusion} is refitted with a flat uninformative 
#' prior to get estimates for the regression coefficients \code{beta}. The posterior means and 95\%-HPD 
#' intervals resulting from this refit are shown with this function. Fused categories have the same 
#' regression coefficient estimates and the same HPD intervals. 
#' 
#' If a full model is fitted (\code{method} in \code{effectFusion} is NULL) or no final model selection is performed 
#' (argument \code{modelSelection} in \code{effectFusion} is NULL), the coefficient estimates are model
#' averaged results.
#'
#' @author Daniela Pauger, Magdalena Leitner <magdalena.leitner@jku.at>
#'
#' @seealso \code{\link{effectFusion}}
#' 
#' @method summary fusion
#' @export
#'
#'
#' @examples
#' ## see example for effectFusion

summary.fusion <- function(object, ...) {
    
    stopifnot(is(object, "fusion"))
    x <- object
    if (x$method == "SpikeSlab" | x$method == "FinMix") {
      if (x$method == "SpikeSlab") 
        cat("\nBayesian effect fusion with spike and slab prior:")
      if (x$method == "FinMix") 
        cat("\nBayesian effect fusion with finite mixture prior:")
    } else {
      cat("\nNo effect fusion performed")
      cat("\nFull model was estimated:")
    }
    
    cat("\n\nCall:\n")
    print(x$call)
    
    cat("\nMCMC:")
    cat("\nM =", x$mcmc$M, "draws after a burn-in of", x$mcmc$burnin)
    if (x$method == "SpikeSlab" | x$method == "FinMix") {
      cat("\nVariable selection started after", x$mcmc$startsel, "iterations\n")
    }
    
    if (x$method == "SpikeSlab") {
      if (is.null(x$prior$tau2_fix))
        cat("\nSpike and slab prior with r = ", x$prior$r, ", g0 = ", x$prior$g0, " and G0 = ", x$prior$G0, sep = "")
      if (!is.null(x$prior$tau2_fix))
        cat("\nSpike and slab prior with r =", x$prior$r, "and tau2_fix =", toString(x$prior$tau2_fix))
    }
    if (x$method == "FinMix") {
      cat("\nFinite mixture prior with e0 =", x$prior$e0, "and p =", x$prior$p)
    }
    if (x$method == "No effect fusion performed. Full model was estimated.") {
      cat("\nA flat, uninformative prior was used")
    }
    
    cat("\n\nPosterior means and 95%-HPD intervals of model refit:\n\n")
    
    if (is.null(x$modelSelection)) {
        postMean <- colMeans(x$fit$beta)
        hpd <- t(apply(x$fit$beta, 2, hpdMCMC))
    } else {
        postMean <- colMeans(x$refit$beta)
        hpd <- t(apply(x$refit$beta, 2, hpdMCMC))
    }
    tab <- cbind(postMean, hpd)
    colnames(tab) <- c("Estimate", "95%-HPD[l]", "95%-HPD[u]")
    rownames(tab) <- createRowNames(x$model, x$data$levelnames, colnames(x$data$X)[which(x$data$types == "c")])
    
    print(round(tab, 3))
    
}