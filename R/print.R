#' @title Print object of class \code{fusion}
#'
#' @description The default print method for a \code{fusion} object.
#'
#' @param x an object of class \code{fusion}
#' @param ... further arguments passed to or from other methods (not used)
#'
#' @details Returns basic information about the model, the number of observations, the number and
#' types of used covariates with their number of categories and MCMC options.
#'
#' @seealso \code{\link{effectFusion}}
#' @method print fusion
#' @export

print.fusion = function(x, ...){

  stopifnot(class(x) == "fusion")

  if(x$method=="SS") cat("Bayesian effect fusion with spike and slab prior:\n")
  if(x$method=="FM") cat("Bayesian effect fusion with finite mixture prior:\n")

  cat("\nCall:\n")
  print(x$call)

  cat("\nModel:", length(x$data$y), "observations")
  cat("\n Covariates:", ncol(x$data$X))
  cat("\n Types:", x$model$n_cont, "continuous,", x$model$n_ord, "ordinal and", x$model$n_nom, "nominal")
  cat("\n with number of categories:", x$model$categories)

  cat("\n\nMCMC:")
  cat("\nM =", x$mcmc$M, "draws after a burn-in of", x$mcmc$burnin)
  cat("\nVariable selection started after", x$mcmc$startsel, "iterations")
}







