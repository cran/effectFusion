#' Bayesian effect fusion for categorical predictors
#'
#' @description This function performs Bayesian variable selection and effect fusion for
#' categorical predictors in a linear regression model. Effect fusion aims at the question which
#' categories of an ordinal or nominal predictor have a similar effect on the response and therefore
#' can be fused to obtain a sparser representation of the model. Effect fusion and variable selection
#' can be obtained either with a prior that has an interpretation as spike and slab prior on the
#' level effect differences or with a sparse finite mixture prior on the level effects. The
#' regression coefficients are estimated with a flat uninformative prior after model selection or
#' model averaged. For posterior inference, an MCMC sampling scheme is used that involves only
#' Gibbs sampling steps.
#'
#' @param y a vector of the response observations
#' @param  X a matrix with covariates, each column representing one covariate. The columns must be
#' ordered according to the type of covariate starting with metric followed by ordinal and nominal covariates.
#' @param types a character vector to specify the type of each covariate; 'c' indicates continuous or metric predictors, 'o'
#' ordinal predictors and 'n' nominal predictors.
#' @param method The method controls the main prior structure that is used for effect fusion. For the prior
#' that has an interpreation as spike and slab prior on level effect differences chosse \code{method = 'SS'}
#' and for the sparse finite mixture prior on the level effects choose \code{method = 'FM'}. See
#' details for a description of the two approaches and their advantages and drawbacks.
#' @param mcmc an (optional) list of MCMC sampling options (see details)
#' @param prior an (optional) list of prior settings and hyper-parameters for the prior (see details).
#' The specification of this list depends on the chosen \code{method}.
#' @param pfp value between 0 and 1 that specifies a threshold level for effect fusion when \code{method='SS'}.
#' Higher values increase the fusion probability, lower values decrease the fusion probability.
#' @param averaged If TRUE (default is FALSE) no model selection is performed and parameter estimates
#' are model averaged results.
#'
#' @details  This function provides identification of categories (of ordinal and nominal predictors) with the
#' same effect on the response and automatic fusion.
#'
#' Two different prior versions for effect fusion and variable selection are available. The first prior version
#' allows a priori for almost perfect as well as almost zero dependence between level effects. This prior has also
#' an interpretation as independent spike and slab prior on all pairwise differences of level effects and correction
#' for the linear dependence of the effect differences. Even though the prior is mainly designed for fusion of level
#' effects, excluding some categories from the model as well as the whole covariate (variable selection) can also be easily accomplished.
#' Excluding a category from the model corresponds to fusion of this category to the baseline and excluding the whole
#' covariate consequently to fusion of all categories to the baseline.
#'
#' The second prior is a modification of the usual spike and slab prior for the regression coefficients by combining
#' the spike at zero with a finite location mixture of normal components. It enables detection of categories with similar
#' effects on the response by clustering the regression effects. Categories with effects that are allocated to the same
#' cluster are fused. Due to the specification with one component located at zero also automatic exclusion of
#' categories without effect on the outcome is provided.
#'
#' In settings with large numbers of categories, we recommend to use the sparse finite mixture prior. The
#' sparse finite mixture prior on the level effects does not take into account the ordering information of ordinal
#' predictors and treats them like nominal predictors.
#'
#' Metric predictors can be included in the model as well and variable selection with spike and slab
#' prior will be performed for these predictors.
#'
#' Details for the model specification (see arguments):
#' \describe{
#'  \item{\code{mcmc}}{A list:}
#'  {\describe{
#'    \item{\code{M}}{number of MCMC iterations after the burn-in phase; defaults to 8000.}
#'    \item{\code{burnin}}{number of MCMC iterations discarded as burn-in; defaults to 2000.}
#'    \item{\code{startsel}}{number of MCMC iterations drawn from the model without performing variable selection;
#'    defaults to 1000.}
#'  }}
#'  \item{\code{prior}}{\describe{\item{}{A list (depending on used \code{method}):}
#'    \item{\code{r}}{variance ratio of slab to spike component; default to 50000. \code{r} has almost no influence on
#'    model selection for ordinal predictors but for nominal covariates low values of \code{r} encourage strong
#'    fusion. \code{r} should be chosen not too samll but still small enought to avoid stickiness of MCMC. We
#'    recommend a value of at least 20000.}
#'    \item{\code{g0}}{shape parameter of inverse gamma distribution when \code{method='SS'}; default
#'    to 5. The default value is a standard choice in variable selection where the tails of spike and slab are not too thin
#'    to cause mixing problems in MCMC.}
#'    \item{\code{G0}}{scale parameter of inverse gamma distribution when \code{method='SS'}; default
#'    to 25. \code{G0} controlls to some extend the sparsity of the model. Smaller values for \code{G0} help to detect also small
#'    level effect differences of nominal predictors, but result in less fusion of categories for ordinal predictors.}
#'    \item{\code{tau2_fix}}{fix value of slab variance when \code{method='SS'} and if no hyper-prior should be used;
#'    default to \code{NULL}. Similar to the scale parameter of its hyper-prior \code{G0}, the fixed variance of the slab
#'    component \code{tau2_fix} can control to some extend the sparsity. Samller values help to detect also small
#'    level effect differences but may result in less fusion of categories for ordinal predictors.}
#'    \item{\code{s0}}{hyper-parameter (shape) of inverse gamma distribution on error variance, used for both
#'    versions of \code{method}; default to 0.}
#'    \item{\code{S0}}{hyper-parameter (scale) of inverse gamma distribution on error variance, used for both
#'    versions of \code{method}; default to 0.}
#'    \item{\code{e0}}{parameter of Dirichlet hyper-prior on mixture weights when \code{method='FM'};
#'    default to 0.01. Due to asymptotic convergence \code{e0} should be chosen smaller than 1. Small values such as
#'    0.01 help to concentrate the model space on sparse solutions.}
#'    \item{\code{p}}{prior parameter to control mixture component variances when \code{method='FM'}; default
#'    to 1000. When \code{hyperprior=FALSE}, larger values of \code{p} lead to less sparsity and it should be chosen
#'    not smaller than 100. If a hyper-prior on the mixture component variances is used (\code{hyperprior=TRUE}), \code{p}
#'    has almost no effect on the sparsity of the model and it should again be not smaller than 100.}
#'    \item{\code{hyperprior}}{logical value if inverse gamma hyper-prior on component variance when \code{method='FM'}
#'    should be used; default to FALSE. The hyper-prior leads to robust results concerning the specification of \code{p}
#'    but also to very sparse solutions.}
#' }}}
#'
#' @return The function returns an object of class \code{fusion} with methods \code{\link{dic}},
#' \code{\link{model}}, \code{\link{print.fusion}}, \code{\link{summary.fusion}} and
#' \code{\link{plot.fusion}}.
#'
#' An object of class \code{fusion} is a named list containing the following elements:
#'
#' \item{\code{fit}}{a named list containing the samples from the posterior distributions of the parameters
#' depending on the used prior structure (\code{method='SS'} or \code{method='FM'}):
#' \describe{
#' \item{\code{beta}}{regression coefficients \eqn{\beta_0} (intercept) and \eqn{\beta}}
#' \item{\code{delta}}{indicator variable \eqn{\delta} for slab component when \code{method='SS'}. The differences of
#' the level effects are assigned either to the spike (\code{delta=0}) or the slab component (\code{delta=0}). If an
#' effect difference is assigned to the spike, the difference is almost zero and the corresponding level effect
#' should be fused.}
#' \item{\code{sgma2}}{error variance \eqn{\sigma^2} of the model}
#' \item{\code{tau2}}{variance \eqn{\tau^2} of slab component when \code{method='SS'}}
#' \item{\code{S}}{latent allocation variable \eqn{S} for mixture components when \code{method='FM'}}
#' \item{\code{eta}}{mixture component weights \eqn{\eta} when \code{method='FM'}}
#' \item{\code{eta0}}{weight of component located at zero \eqn{\eta_0} when \code{method='FM'}}
#' \item{\code{mu}}{mixture component means \eqn{\mu} when \code{method='FM'}}
#' \item{\code{N_jl_matrix}}{matrix with numbers of regression effects of each covariate \eqn{j} assigned
#' to mixture component \eqn{l} called \eqn{N_{jl}} when \code{method='FM'}}
#' \item{\code{N_j0_matrix}}{matrix with numbers of regression effects of each covariate \eqn{j} assigned
#' to component located at zero called \eqn{N_{j0}} when \code{method='FM'}}
#' }}
#' \item{\code{refit}}{a named list containing samples from the posterior distributions of the parameters of
#' the model refit:
#' \describe{
#' \item{\code{beta}}{regression coefficients \eqn{\beta_0} (intercept) and \eqn{\beta}}
#' \item{\code{sgma2}}{error variance \eqn{\sigma^2} of the model}
#' \item{\code{model}}{vector of 0 and 1 representing the selected model and based on pairs of categories}
#' }}
#' \item{\code{method}}{see arguments}
#' \item{\code{data}}{a list containing the data \code{y}, \code{X}, \code{types} and dummy coded design
#' matrix \code{X_dummy}}
#' \item{model}{a named list containing information on the model
#' \describe{
#' \item{\code{categories}}{number of categories for categorical predictors}
#' \item{\code{diff}}{number of pairwise level effect differences}
#' \item{\code{n_nom}}{number of nominal predictors}
#' \item{\code{n_ord}}{number of ordinal predictors}
#' \item{\code{n_cont}}{number of metric predictors}
#' \item{\code{lNom}}{an index for categorical predictors in the design matrix}
#' \item{\code{A_diag}}{ matrix for the linear restrictions of the pairwise level effect differences}
#' }}
#' \item{\code{prior}}{see details for prior}
#' \item{\code{mcmc}}{see details for mcmc}
#' \item{\code{call}}{function call}
#'
#' @note The function can be used for ordinal and/or nominal predictors and metric covariates
#' can additionally be included in the model. Binary covariates as a special case of nominal predictors
#' can be included as well.
#'
#' The sparse finite mixture prior approach does not take into account the ordering information of ordinal
#' predictors. Ordinal predictors are treated as nominal predictors.
#'
#' For large models and more than 15,000 MCMC iterations, some thinning of the MCMC when using the sparse finite
#' mixture prior is performed due to computational issues.
#'
#' @author Daniela Pauger <daniela.pauger@jku.at>, Helga Wagner, Gertraud Malsiner-Walli
#'
#' @references  Pauger, D. and Wagner, H. (2016). Bayesian effect fusion for categorical predictors.
#' Submitted manuscript.
#'
#' @exportClass fusion
#' @export
#'
#' @examples
#' \dontrun{
#' # ----------- Load simulated data set 'sim1'
#' data(sim1)
#' y=sim1$y
#' X=sim1$X
#' types=sim1$types
#'
#' # ----------- Bayesian effect fusion for simulated data set with spike and slab prior
#' m1 <- effectFusion(y,X,types,method="SS")
#'
#' # print, summarize and plot results
#' print(m1)
#' summary(m1)
#' plot(m1)
#'
#' # evaluate model and model criteria
#' model(m1)
#' dic(m1)
#'
#' # ----------- Use finite mixture prior for comparison
#' # change prior parameter specification and use hyper-prior on component variances
#' m2 <- effectFusion(y,X,types,prior=list(p=1000,hyperprior=TRUE),method="FM")
#'
#' # summarize and plot results
#' print(m2)
#' summary(m2)
#' plot(m2)
#' model(m2)
#' dic(m2)
#'
#' # ------------  Use model averaged coefficient estimates
#' m3 <- effectFusion(y,X,types,method="SS",averaged=TRUE)
#' summary(m3)
#'
#'}

effectFusion = function(y,X,types,method,mcmc=list(),prior=list(),pfp=0.5,averaged=FALSE){

  cl <- match.call()

  if (is.null(y) || is.null(X))
    stop("need 'y' and 'X' argument")
  if (!is.matrix(X))
    stop("'X' must be a matrix")
  if (any(is.na(X)))
    stop("NA values in 'X' not allowed")
  if (any(is.na(y)))
    stop("NA values in 'y' not allowed")
  if (length(y) != nrow(X))
    stop("'y' and 'nrow(X)' must have same length")
  if (length(types) != NCOL(X))
    stop("'types' and 'ncol(X)' must have same length")
  if (any(is.na(match(types,c("c","o","n")))))
    stop("invalid argument in 'types'")
  if(!"o" %in% types & !"n" %in% types)
    stop("No categorical predictors")

  if("o" %in% types & method=="FM") {
    cat("WARNING: Finite mixture prior treats ordinal predictors as nominal\n\n")
    types[types=="o"]="n"
  }

  if(any(apply(X[,types=="o" | types == "n"],2,min)!=1))
    stop("start category labels with '1'")

  defaultMCMC = list(M = 8000, burnin = 2000, startsel = 1000)
  mcmc = utils::modifyList(defaultMCMC, as.list(mcmc))

  nVar=ncol(X)
  ind_cont=ind_ord=ind_nom=rep(F,nVar)
  ind_cont[types=="c"]=T
  ind_nom[types=="n"]=T
  ind_ord[types=="o"]=T

  data=list(y=y,X=X,ind_cont=ind_cont,ind_nom=ind_nom,ind_ord=ind_ord)
  mvars=createModelvars(data)
  model = createModel(data,mvars)

  # for variable selection for continuous variables
  # treat them as nominal with 2 categories
  if(method=="FM" & model$n_cont>0){
    model$categories=c(rep(2,model$n_cont),model$categories)
    model$diff=c(rep(1,model$n_cont),model$diff)
    model$n_cont=0
    model$n_nom=model$n_nom+1
  }

  mats = getReparmats(model)

  if(method=="SS"){
    mcmc_res = mcmcSs(y,X=mvars$X_dummy,model,prior,mcmc,mats)
    if(!averaged){
      incl_prob=colMeans(mcmc_res$delta)
      model_sel=selectModel(incl_prob,model,strategy="pip",PIP=1-pfp)
      refit_res=modelRefit(model,model_sel,data,prior_flat=list(s0=0,S0=0,tau2_fix=1000,conj=FALSE),mcmc_flat=list(M=3000,burnin1=1001))
    }
  }
  if(method=="FM"){
    mcmc_res = mcmcMix(y,X=mvars$X_dummy,model,prior,mcmc)
    if(!averaged){
      incl_prob = inclProb(S=mcmc_res$S,mvars,model)
      model_sel=selectModel(incl_prob,model,strategy="pam")
      if("c" %in% types) model_sel=c(rep(1,sum(ind_cont)),model_sel)
      refit_res=modelRefit(model,model_sel,data,prior_flat=list(s0=0,S0=0,tau2_fix=1000,conj=FALSE),mcmc_flat=list(M=3000,burnin1=1001))
    }
  }

  if(method=="FM" & sum(types=="c")>0){
    cont=sum(types=="c")
    model$n_cont=cont
    model$n_nom=model$n_nom-cont
    model$categories=model$categories[-c(1:cont)]
    model$diff=model$diff[-c(1:cont)]
  }


  if(averaged){
    ret=list(fit=mcmc_res[names(mcmc_res)!="prior"],method=method,data=list(y=y,X=X,X_dummy=mvars$X_dummy,types=types),model=model,prior=mcmc_res$prior,mcmc=mcmc,call=cl,averaged=averaged)
  } else {
    refit_res$model=model_sel
    ret=list(fit=mcmc_res[names(mcmc_res)!="prior"],refit=refit_res,method=method,data=list(y=y,X=X,X_dummy=mvars$X_dummy,types=types),model=model,prior=mcmc_res$prior,mcmc=mcmc,call=cl,averaged=averaged)

  }

  class(ret) <- "fusion"
  return(ret)


}




