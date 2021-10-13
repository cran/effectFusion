#' Simulated data set 1
#' 
#' The simulated data set \code{sim1} illustrates a setting with 500 observations from a linear 
#' regression model with normal response, 4 ordinal and 4 nominal predictors. Two regressors 
#' have 8 and two have 4 categories for each type of covariate (ordinal and nominal). Regression 
#' effects are set to \eqn{\beta_1 = (0, 1, 1, 2, 2, 4, 4)} and \eqn{\beta_3 = (0, -2, -2)} for the 
#' ordinal and \eqn{\beta_5 = (0, 1, 1, 1, 1, -2, -2)} and \eqn{\beta_7 = (0, 2, 2)} for the nominal 
#' covariates, and \eqn{\beta_h = 0} for h = 2, 4, 6, 8. Levels of the predictors are generated with 
#' probabilities \eqn{\pi_h = (0.1, 0.1, 0.2, 0.05, 0.2, 0.1, 0.2, 0.05)} and \eqn{\pi_h = (0.1, 0.4, 
#' 0.2, 0.3)} for regressors with 8 and 4 levels, respectively. For more details on the 
#' simulation setting see Pauger and Wagner (2019).
#' 
#' @docType data
#' @usage data(sim1)
#' @format A named list containing the following four variables:
#' \describe{
#'  \item{\code{y}}{vector with 500 observations of a normal response variable}
#'  \item{\code{X}}{matrix with 8 categorical predictors}
#'  \item{\code{beta}}{vector with coefficients used for data generation}
#'  \item{\code{types}}{character vector with types of covariates, 'o' for ordinal and 'n' for 
#'  nominal covariates}
#' }
#' 
#' @references {Pauger, D., and Wagner, H. (2019). Bayesian Effect Fusion for Categorical Predictors.
#' \emph{Bayesian Analysis}, \strong{14(2)}, 341-369. \doi{10.1214/18-BA1096}
#' }
#'  
#' @seealso \code{\link{effectFusion}}
#' @name sim1
#' @keywords datasets
NULL
