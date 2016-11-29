#' @title Summary of object of class \code{fusion}
#' 
#' @description Returns basic information about the model and the priors, MCMC details and posterior
#' means from refit of the selected model and 95\%-HPD intervals for the regression effects.
#'
#' @param object an object of class \code{fusion}
#' @param ... further arguments passed to or from other methods (not used)
#'
#' @details The model selected with function \code{effectFusion} is refitted with a flat uninformative 
#' prior to get estimates of the regression coefficients \code{beta}. The posterior means and 95\%-HPD 
#' intervals resulting from this refit are shown with this function \code{summary.fusion}. Fused 
#' categories have the same estimates of the regression coefficients and the same HPD intervals. 
#' 
#' 
#' If no model selection is performed (argument \code{averaged} in \code{effectFusion} is TRUE), the coefficient estimates are model
#' averaged.
#'
#' @author Daniela Pauger <daniela.pauger@jku.at>
#'
#' @seealso \code{\link{effectFusion}}
#' 
#' @method summary fusion
#' @export
#'
#'
#' @examples
#' ## see example for effectFusion

summary.fusion = function(object,...){

  stopifnot(class(object) == "fusion")
  x=object

  if(x$method=="SS") cat("Bayesian effect fusion with spike and slab prior:")
  if(x$method=="FM") cat("Bayesian effect fusion with finite mixture prior:")

  cat("\n\nCall:\n")
  print(x$call)

  cat("\nMCMC:")
  cat("\nM =", x$mcmc$M, "draws after a burn-in of", x$mcmc$burnin)
  cat("\nVariable selection started after", x$mcmc$startsel, "iterations\n")

  if(x$method=="SS") cat("\nSpike and slab prior with r =",x$prior$r, "and Q =",x$prior$Q)
  if(x$method=="FM") cat("\nFinite mixture prior with e0 =",x$prior$e0, "and p =",x$prior$p)

  cat("\n\nPosterior means and 95%-HPD intervals of model refit:\n\n")

  if(x$averaged){
    postMean=colMeans(x$fit$beta)
    hpd=t(apply(x$fit$beta,2,hpdMCMC))
  } else {
    postMean=colMeans(x$refit$beta)
    hpd=t(apply(x$refit$beta,2,hpdMCMC))
  }
  tab=cbind(postMean,hpd)
  colnames(tab)=c("Estimate","95%-HPD[l]","95%-HPD[u]")
  rownames(tab)=createRowNames(x$model)

  print(round(tab,3))

}







