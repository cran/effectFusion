#' Simulated data set 3
#' 
#' The simulated data set \code{sim3} considers a setting with 2000 observations from a logistic 
#' regression model. The number and types of predictors, the regression effects and the level probabilites 
#' of the predictors are the same as for \code{sim1}. The number of observations was increased as the uncertainty
#' is usually higher for logistic regression compared to linear regression with normal response.
#' 
#' @docType data
#' @usage data(sim3)
#' @format A named list containing the following four variables:
#' \describe{
#'  \item{\code{y}}{vector with 2000 observations of a binary response variable}
#'  \item{\code{X}}{matrix with 8 categorical predictors}
#'  \item{\code{beta}}{vector with coefficients used for data generation}
#'  \item{\code{types}}{character vector with types of covariates, 'o' for ordinal and 'n' for 
#'  nominal covariates}
#' }
#' 
#'  
#' @seealso \code{\link{effectFusion}, \link{sim1}}
#' @name sim3
#' @keywords datasets
NULL
