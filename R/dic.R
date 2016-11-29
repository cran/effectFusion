#' DIC
#' @description
#' This function computes the DIC (deviance information criterion) for the estimated model in a \code{fusion} object.
#'
#' @param x an object of class \code{fusion}
#' 
#' @details The DIC can be easily computed from the MCMC output and is defined as 
#' \eqn{DIC = 2 \overline{D(\theta)} - D(\overline{\theta})}, where \eqn{\overline{D(\theta)} = 
#' \frac{1}{M} \sum \limits_{m=1}^{M} D(\theta^{(m)})} is the average posterior deviance and 
#' \eqn{D(\bar{\theta})} is the deviance evaluated at \eqn{\bar{\theta} = \frac{1}{M} \sum 
#' \limits_{m=1}^{M} \theta^{(m)}}. \eqn{\theta^{(m)}} are samples from the posterior of the model 
#' and M is the number of MCMC iterations. 
#'
#' @return The DIC for the estimated model in the \code{fusion} object.
#' @author Daniela Pauger <daniela.pauger@jku.at>
#' @export
#' 
#' @seealso \code{\link{effectFusion}}
#' 
#' @examples
#' ## see example for effectFusion
#' 
#' @references Spiegelhalter, D., Best, N., Carlin, B. and van der Linde, A. (2002). Bayesian measures of model 
#' complexity and fit. \emph{J. R. Statist. Soc. B}, \strong{64}(4), 583-639.
#' 


dic = function(x){

  stopifnot(class(x) == "fusion")

  y=x$data$y
  X=x$data$X_dummy
  if(x$averaged){
    beta=x$fit$beta
    s2=x$fit$sgma2
  } else {
    beta=x$refit$beta
    s2=x$refit$sgma2
  }


  it=nrow(beta)
  beta_hat=colMeans(beta)
  N=length(y)
  D=0

  for(i in 1:it){
    D = D + (-2*(-N/2 * log(2*pi*s2[i]) - 1/(2*s2[i]) * t(y-X%*%beta[i,]) %*% (y-X%*%beta[i,])))
  }
  D_quer=D/it
  DIC = 2*D_quer - (-2*(-N/2 * log(2*pi*mean(s2)) - 1/(2*mean(s2)) * t(y-X%*%beta_hat) %*% (y-X%*%beta_hat)))

  as.numeric(DIC)

}

