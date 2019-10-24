#' Bayesian effect fusion for categorical predictors
#'
#' @description This function performs Bayesian variable selection and effect fusion for
#' categorical predictors in linear and logistic regression models. Effect fusion aims at the question which
#' categories of an ordinal or nominal predictor have a similar effect on the response and therefore
#' can be fused to obtain a sparser representation of the model. Effect fusion and variable selection
#' can be obtained either with a prior that has an interpretation as spike and slab prior on the
#' level effect differences or with a sparse finite mixture prior on the level effects. The
#' regression coefficients are estimated with a flat uninformative prior after model selection or
#' by using model averaged results. For posterior inference, a MCMC sampling scheme is used that involves only
#' Gibbs sampling steps. The sampling schemes for linear and logistic regression are almost identical as
#' in the case of logistic regression a data augmentation strategy (Polson et al. (2013)) is used that requires
#' only one additional step to sample from the Polya-Gamma distribution.
#'
#' @param y a vector of the response observations (continuous if \code{family =} \code{'gaussian'} and 0/1 if \code{family =} \code{'binomial'})
#' @param X a data frame with covariates, each column representing one covariate. Ordinal and nominal covariates should
#' be of class \code{factor}
#' @param types a character vector to specify the type of each covariate; 'c' indicates continuous or metric predictors, 'o'
#' ordinal predictors and 'n' nominal predictors
#' @param method controls the main prior structure that is used for effect fusion. For the prior
#' that has an interpretation as spike and slab prior on the level effect differences choose \code{method =} \code{'SpikeSlab'}
#' and for the sparse finite mixture prior on the level effects choose \code{method =} \code{'FinMix'}. See
#' details for a description of the two approaches and their advantages and drawbacks. For comparison purposes
#' it is also possible to fit a full model without performing any effect fusion (\code{method =} \code{NULL})
#' @param prior an (optional) list of prior settings and hyper-parameters for the prior (see details).
#' The specification of this list depends on the chosen \code{method} and the selected \code{family}
#' @param mcmc an (optional) list of MCMC sampling options (see details)
#' @param mcmcRefit an (optional) list of MCMC sampling options for the refit of the selected model (see details)
#' @param family indicates whether linear (default, \code{family =} \code{'gaussian'}) or logistic regression (\code{family =} \code{'binomial'})
#' should be performed
#' @param modelSelection if \code{modelSelection =} \code{'binder'} the final model is selected by minimising the expected posterior binder's loss
#' using an algorithm of Lau and Green (2008) for the spike and slab model and an algorithm of Rastelli and Friel (2016)
#' for the finite mixture approach. Alternatively, \code{modelSelection =} \code{'pam'} can be specified for the sparse finite mixture
#' model. In that case, the final model is selected by using pam clustering and the silhouette coefficient (see Malsiner-Walli et al., 2018 for details). 
#' If \code{modelSelection} is \code{NULL} no final model selection is performed and parameter estimates
#' are model averaged results. \code{modelSelection =} \code{'binder'} is the default value. \code{modelSelection =} \code{'pam'}
#' is only available for \code{method =} \code{'FinMix'}. If \code{method =} \code{'SpikeSlab'} and \code{modelSelection =} \code{'pam'},
#' \code{modelSelection} is automatically set to \code{'binder'}. For the finite mixture approach
#' we recommend to use \code{modelSelection =} \code{'binder'}, as this algorithm provides - in contrast to pam clustering and the
#' the silhouette coefficient - the opportunity to exclude whole covariates.
#' @param returnBurnin if \code{TRUE} (default is \code{FALSE}) the burn-in iterations of the MCMC sampling process are returned as well. 
#' This can be for example used to check convergence. Returning the burn-in does not influence the results of \code{\link{dic}},
#' \code{\link{model}}, \code{\link{plot}}, \code{\link{print}} and \code{\link{summary}}.
#'
#' @details This function provides identification of categories (of ordinal and nominal predictors) with the
#' same effect on the response and their automatic fusion.
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
#' a spike at zero with a finite location mixture of normal components. It enables detection of categories with similar
#' effects on the response by clustering the regression effects. Categories with effects that are allocated to the same
#' cluster are fused. Due to the specification with one component located at zero also automatic exclusion of
#' levels and whole covariates without any effect on the outcome is provided. However, when using \code{modelSelection =} \code{'pam'},
#' it is not possible to exclude whole covariates as the Silhouette coefficient does not allow for one cluster solutions. Therefore,
#' we recommend to use \code{modelSelection =} \code{'binder'}.
#'
#' In settings with large numbers of categories, we recommend to use the sparse finite mixture prior for computational reasons.
#' It is important to note that the sparse finite mixture prior on the level effects does not take into account the ordering information of ordinal
#' predictors and treats them like nominal predictors, whereas in the spike and slab case fusion is restricted to adjacent categories for ordinal predictors.
#'
#' Metric predictors can be included in the model as well and variable selection will be performed also for these predictors.
#'
#' If \code{modelSelection =} \code{NULL}, no final model selection is performed and model averaged results are returned. When \code{modelSelection =} \code{'binder'} the final
#' model is selected by minimizing the expected posterior Binder loss for each covariate separately.
#' Additionally, there is a second option for the finite mixture prior (\code{modelSelection =} \code{'pam'}) which performs model selection
#' by identifying the optimal partition of the effects using PAM clustering and the silhouette coefficient.
#' 
#' For comparison purposes it is also possible to fit a full model instead of performing effect fusion (\code{method =} \code{NULL}).
#' All other functions provided in this package, such as \code{\link{dic}} or \code{\link{summary}}, do also work for the full model.
#'
#' Details for the model specification (see arguments):
#' \describe{
#'  
#'  \item{\code{prior}}{\describe{\item{}{a list (depending on used \code{method} and specified \code{family}). If \code{method =} \code{NULL}, all prior
#'  specifications are ignored and a flat, uninformative prior is assigned to the level effects}
#'    \item{\code{r}}{variance ratio of slab to spike component; default to 50000 if \code{family =} \code{'gaussian'} and 
#'    5000000 if \code{family =} \code{'binomial'}. \code{r} should be chosen not too small but still small enough to avoid stickiness of MCMC. We
#'    recommend a value of at least 20000.}
#'    \item{\code{g0}}{shape parameter of inverse gamma distribution on \eqn{\tau^2} when \code{tau2_fix =} \code{NULL} and \code{method =} \code{'SpikeSlab'}; default
#'    to 5. The default value is a standard choice in variable selection where the tails of spike and slab component are not too thin
#'    to cause mixing problems in MCMC.}
#'    \item{\code{G0}}{scale parameter of inverse gamma distribution on \eqn{\tau^2} when \code{tau2_fix =} \code{NULL} and \code{method =} \code{'SpikeSlab'}; default
#'    to 25. \code{G0} controls to some extend the sparsity of the model. Smaller values for \code{G0} help to detect also small
#'    level effect differences, but result in less fusion of categories.}
#'    \item{\code{tau2_fix}}{If \code{tau2_fix =} \code{NULL}, an inverse gamma hyper-prior is specified on \eqn{\tau^2}. 
#'    However, the value of the slab variance can also be fixed for each covariate instead of using a hyperprior.
#'    \code{tau2_fix} is only of interest if \code{method =} \code{'SpikeSlab'}.
#'    Default to \code{NULL}. Similar to the scale parameter \code{G0}, the fixed variance of the slab
#'    component \code{tau2_fix} can control to some extend the sparsity.}
#'    \item{\code{e0}}{parameter of Dirichlet hyper-prior on mixture weights when \code{method =} \code{'FinMix'};
#'    default to 0.01. \code{e0} should be chosen smaller than 1 in order to encourage empty components. Small values such as
#'    0.01 help to concentrate the model space on sparse solutions.}
#'    \item{\code{p}}{prior parameter to control mixture component variances when \code{method =} \code{'FinMix'}; values between 100 and 100000
#'    led to good results in simulation studies; default to 100 if \code{family =} \code{'gaussian'} and 1000 if \code{family =} \code{'binomial'}. We recommend to try different values for this
#'    prior parameter and compare the models using \code{\link{dic}}. When \code{hyperprior =} \code{FALSE}, larger values of \code{p} lead to less sparsity and it should be chosen
#'    not smaller than 100. If a hyper-prior on the mixture component variances is used (\code{hyperprior =} \code{TRUE}), \code{p}
#'    has almost no effect on the sparsity of the model and it should again be not smaller than 100.}
#'    \item{\code{hyperprior}}{logical value if inverse gamma hyper-prior on component variance should be specified when \code{method =} \code{'FinMix'};
#'    default to FALSE. The hyper-prior leads to robust results concerning the specification of \code{p}
#'    but also to very sparse solutions.}
#'    \item{\code{s0}}{hyper-parameter (shape) of inverse gamma distribution on error variance, used for both
#'    versions of \code{method}, but only for \code{family =} \code{'gaussian'}; default to 0.}
#'    \item{\code{S0}}{hyper-parameter (scale) of inverse gamma distribution on error variance, used for both
#'    versions of \code{method}, but only for \code{family =} \code{'gaussian'}; default to 0.}
#' }}
#' 
#' \item{\code{mcmc}}{\describe{\item{}{a list:}
#'    \item{\code{M}}{number of MCMC iterations after the burn-in phase; default to 20000 for effect fusion models and 3000 for full models.}
#'    \item{\code{burnin}}{number of MCMC iterations discarded as burn-in; default to 5000 for effect fusion models and 1000 for full models.}
#'    \item{\code{startsel}}{number of MCMC iterations drawn from the model without performing effect fusion;
#'    default to 1000 for effect fusion models and 0 for full models.}
#' }}
#' 
#' \item{\code{mcmcRefit}}{\describe{\item{ }{a list (not necessary if \code{modelSelection =} \code{NULL} or \code{method =} \code{NULL}):}
#'    \item{\code{M_refit}}{number of MCMC iterations after the burn-in phase for the refit of the selected model; default to 3000.}
#'    \item{\code{burnin_refit}}{number of MCMC iterations discarded as burn-in for the refit of the selected model; default to 1000.}
#' }}}
#' 
#'
#' @return The function returns an object of class \code{fusion} with methods \code{\link{dic}},
#' \code{\link{model}}, \code{\link{print}}, \code{\link{summary}} and \code{\link{plot}}.
#'
#' An object of class \code{fusion} is a named list containing the following elements:
#'
#' \describe{
#' \item{\code{fit}}{\describe{\item{}{a named list containing the samples from the posterior distributions of the parameters
#' depending on the used prior structure (\code{method =} \code{'SpikeSlab'}, \code{method =} \code{'FinMix'} or \code{method =} \code{NULL}):}
#' \item{\code{beta}}{regression coefficients \eqn{\beta_0} (intercept) and \eqn{\beta}}
#' \item{\code{delta}}{indicator variable \eqn{\delta} for slab component when \code{method =} \code{'SpikeSlab'}. The differences of
#' the level effects are assigned either to the spike (\code{delta = 0}) or the slab component (\code{delta = 1}). If an
#' effect difference is assigned to the spike component, the difference is almost zero and the corresponding level effects
#' should be fused.}
#' \item{\code{tau2}}{variance \eqn{\tau^2} of slab component when \code{method =} \code{'SpikeSlab'}. If no hyperprior on \eqn{\tau^2} is specified, \code{tau2} contains the fixed values for \eqn{\tau^2}.}
#' \item{\code{S}}{latent allocation variable \eqn{S} for mixture components when \code{method =} \code{'FinMix'}}
#' \item{\code{eta}}{mixture component weights \eqn{\eta} when \code{method =} \code{'FinMix'}}
#' \item{\code{eta0}}{weights of components located at zero \eqn{\eta_0} when \code{method =} \code{'FinMix'}}
#' \item{\code{mu}}{mixture component means \eqn{\mu} when \code{method =} \code{'FinMix'}}
#' \item{\code{sgma2}}{error variance \eqn{\sigma^2} of the model (only for \code{family =} \code{'gaussian'})}
#' }}
#' \item{\code{fit_burnin}}{\describe{\item{}{a named list containing the same elements as \code{fit} including the burnin-phase, if \code{returnBurnin =} \code{TRUE},
#' \code{NULL} otherwise. The elements that correspond to the model selection procedure, e.g. \code{delta} or \code{S}, are \code{NA} for the first \code{startsel} iterations.}}}
#' \item{\code{refit}}{\describe{\item{}{a named list containing samples from the posterior distributions of the parameters of
#' the model refit (only if \code{method} and \code{modelSelection} are unequal to \code{NULL}):}
#' \item{\code{beta}}{\describe{regression coefficients including the intercept in the model with fused levels}}
#' \item{\code{sgma2}}{\describe{error variance of the model with fused levels (only for \code{family =} \code{'gaussian'})}}
#' \item{\code{X_dummy_fused}}{\describe{the dummy coded design matrix with fused levels}}
#' \item{\code{model}}{\describe{vector of zeros and ones representing the selected model based on pairs of categories}}
#' }}
#' \item{\code{method}}{see arguments}
#' \item{\code{family}}{see arguments}
#' \item{\code{data}}{\describe{\item{}{a named list containing the data \code{y}, \code{X}, the dummy coded design matrix \code{X_dummy},
#'  \code{types} and \code{levelnames} of ordinal and nominal predictors}
#' \item{\code{model}}{a named list containing information on the full, initial model}
#' \item{\code{categories}}{number of categories for categorical predictors}
#' \item{\code{diff}}{number of pairwise level effect differences}
#' \item{\code{n_cont}}{number of metric predictors}
#' \item{\code{n_ord}}{number of ordinal predictors}
#' \item{\code{n_nom}}{number of nominal predictors}
#' }}
#' \item{\code{prior}}{see details for prior}
#' \item{\code{mcmc}}{see details for mcmc}
#' \item{\code{mcmcRefit}}{see details for mcmcRefit}
#' \item{\code{modelSelection}}{see arguments}
#' \item{\code{returnBurnin}}{see arguments}
#' \item{\code{numbCoef}}{number of estimated regression coefficients (based on the refitted model if effect fusion and final model selection is performed, otherwise based on model averaged results or the full model, respectively)}
#' \item{\code{call}}{function call}
#' }
#'
#' @note The function can be used for ordinal and/or nominal predictors and metric covariates
#' can additionally be included in the model. Binary covariates as a special case of nominal predictors
#' can be included as well.
#'
#' The sparse finite mixture prior approach does not take into account the ordering information of ordinal
#' predictors. Ordinal predictors are treated as nominal predictors, whereas in the spike and slab case fusion is restricted to adjacent categories for ordinal predictors.
#'
#' For large models and more than 15,000 MCMC iterations, some thinning of the MCMC when using the sparse finite
#' mixture prior is performed due to computational issues.
#'
#' @author Daniela Pauger, Magdalena Leitner <effectfusion.jku@gmail.com>, Helga Wagner, Gertraud Malsiner-Walli
#'
#' @references {Pauger, D., and Wagner, H. (2018). Bayesian Effect Fusion for Categorical Predictors.
#' \emph{Bayesian Analysis}, in print.
#'
#' Malsiner-Walli, G., Pauger, D., and Wagner, H. (2018). Effect Fusion Using Model-Based Clustering.
#' \emph{Statistical Modelling}, \strong{18(2)}, 175-196.
#'
#' Polson, N.G., Scott, J.G., and Windle, J. (2013). Bayesian Inference for Logistic Models Using
#' Polya-Gamma Latent Variables. \emph{Journal of the American Statistical Association}, \strong{108(504)}, 1339-1349.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # ----------- Load simulated data set 'sim1' for linear regression
#' data(sim1)
#' y = sim1$y
#' X = sim1$X
#' types = sim1$types
#'
#' # ----------- Bayesian effect fusion for simulated data set with spike and slab prior
#' m1 <- effectFusion(y, X, types, method = 'SpikeSlab')
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
#' m2 <- effectFusion(y, X, types, method = 'FinMix')
#'
#' # summarize and plot results
#' print(m2)
#' summary(m2)
#' plot(m2)
#' model(m2)
#' dic(m2)
#' 
#' # change prior parameter specification
#' m3 <- effectFusion(y, X, types, prior= list(p = 10^3), method = 'FinMix')
#' plot(m3)
#'
#' # ------------  Use model averaged coefficient estimates
#' m4 <- effectFusion(y, X, types, method = 'SpikeSlab', modelSelection = NULL)
#' summary(m4)
#' 
#' # ------------  Estimate full model for comparison purposes
#' m5 <- effectFusion(y, X, types, method = NULL)
#' summary(m5)
#' plot(m5)
#' dic(m5)
#'
#' # ----------- Load simulated data set 'sim3' for logistic regression
#' data(sim3)
#' y = sim3$y
#' X = sim3$X
#' types = sim3$types
#'
#' # ----------- Bayesian effect fusion for simulated data set with finite mixture prior
#' m6 <- effectFusion(y, X, types, method = 'FinMix', prior = list(p = 10^4), family = 'binomial')
#'
#' # look at the results
#' print(m6)
#' summary(m6)
#' plot(m6)
#' model(m6)
#' dic(m6)
#'
#' # ----------- Use spike and slab prior for comparison
#' m7 <- effectFusion(y, X, types, method = 'SpikeSlab', family = 'binomial', returnBurnin = TRUE)
#'
#' # summarize and evaluate results
#' print(m7)
#' summary(m7)
#' plot(m7)
#' model(m7)
#' dic(m7)
#'}

effectFusion <- function(y, X, types, method, prior = list(), mcmc = list(), mcmcRefit = list(), family = "gaussian", 
    modelSelection = "binder", returnBurnin = FALSE) {
    
    cl <- match.call()
    if (is.null(y) || is.null(X))
        stop("need 'y' and 'X' argument")
    if (!is.vector(y) & !is.matrix(y))
        stop("'y' must be a vector or a matrix")
    if (is.matrix(y)) {
        if(ncol (y) != 1) {
          stop("'y' must be a matrix with one column")
        }
    }
    if (!is.data.frame(X)) 
        stop("'X' must be a data.frame")
    if (any(is.na(X))) 
        stop("NA values in 'X' not allowed")
    if (any(is.na(y))) 
        stop("NA values in 'y' not allowed")
    if (length(y) != nrow(X)) 
        stop("'y' and 'nrow(X)' must have same length")
    if (!is.vector(types))
        stop("'types' must be a vector")
    if (length(types) != NCOL(X)) 
        stop("'types' and 'ncol(X)' must have same length")
    if (any(is.na(match(types, c("c", "o", "n"))))) 
        stop("invalid argument in 'types'")
    if (!"o" %in% types & !"n" %in% types) 
        stop("No categorical predictors")
    ordering <- c(which(types == "c"), which(types == "o"), which(types == "n"))
    X <- X[, ordering, drop = F]
    types <- types[ordering]
    if (any(names(prior) == "tau2_fix")) {
      if (!is.null(prior$tau2_fix)) {
        prior$tau2_fix <- prior$tau2_fix[ordering]
      }
    }
    if (!is.null(method)) {
      if (method != "SpikeSlab" & method != "FinMix") {
        stop("no valid 'method' specified")
      }
      if ("o" %in% types & method == "FinMix") {
        warning("Finite mixture prior treats ordinal predictors as nominal.")
        types[types == "o"] <- "n"
      }
    }

    X_out <- X
    levelnames <- lapply(X, levels)
    levelnames[sapply(levelnames, is.null)] <- NULL
    X <- sapply(X, as.numeric)
    if (family != "gaussian" & family != "binomial")
        stop("'family' can either be 'gaussian' or 'binomial'")
    if (family == "binomial" & any(y != 1 & y != 0))
        stop("'y' has to be binary when 'family' is of type binomial")
    if (length(prior) > 0) {
      if (!is.null(method)) {
        if (method == "SpikeSlab") {
          if (family == "gaussian") {
            if (any(!names(prior) %in% c("r", "g0", "G0", "tau2_fix", "s0", "S0"))) {
              stop("Invalid prior parameters specified.")
            }
          }
          if (family == "binomial") {
            if (any(!names(prior) %in% c("r", "g0", "G0", "tau2_fix"))) {
              stop("Invalid prior parameters specified.")
            }
          }
        }
        if (method == "FinMix") {
          if (family == "gaussian") {
            if (any(!names(prior) %in% c("e0", "p", "hyperprior", "s0", "S0"))) {
              stop("Invalid prior parameters specified.")
            }
          }
          if (family == "binomial") {
            if (any(!names(prior) %in% c("e0", "p", "hyperprior"))) {
              stop("Invalid prior parameters specified.")
            }
          }
        }
      } else {
        stop("Full model is estimated. Invalid prior parameters specified.")
      }
    }
    if (length(mcmc) > 0) {
      if (any(!names(mcmc) %in% c("M", "burnin", "startsel"))) {
        stop("Invalid mcmc parameters specified.")
      }
    }
    if (length(mcmcRefit) > 0) {
      if (any(!names(mcmcRefit) %in% c("M_refit", "burnin_refit"))) {
        stop("Invalid mcmc parameters for the refit specified.")
      }
    }
    if (!is.null(modelSelection)) {
      if (modelSelection != "binder" & modelSelection != "pam") {
        stop("'modelSelection' has to be either 'binder' or 'pam' or 'NULL'")
      }
    }
    if (!isFALSE(returnBurnin) & !isTRUE(returnBurnin))
      stop("'returnBurnin' has to be either 'TRUE' or 'FALSE'")
    
    if (!is.null(method)) {
      defaultMCMC <- list(M = 20000, burnin = 5000, startsel = 1000)
    } else {
      defaultMCMC <- list(M = 3000, burnin = 1000, startsel = 0)
    }
    mcmc <- utils::modifyList(defaultMCMC, as.list(mcmc))
    if (mcmc$startsel > mcmc$burnin) {
      stop("Increase 'burnin' or decrease 'startsel'. Model selection has to start within the burnin-phase.")
    }
    defaultMCMCrefit <- list(M_refit = 3000, burnin_refit = 1000)
    mcmcRefit <- utils::modifyList(defaultMCMCrefit, as.list(mcmcRefit))
    
    nVar <- ncol(X)
    ind_cont <- ind_ord <- ind_nom <- rep(F, nVar)
    ind_cont[types == "c"] <- T
    ind_nom[types == "n"] <- T
    ind_ord[types == "o"] <- T
    
    data <- list(y = y, X = X, ind_cont = ind_cont, ind_nom = ind_nom, ind_ord = ind_ord)
    mvars <- createModelvars(data)
    model <- createModel(data, mvars)
    
    # for variable selection for continuous variables treat them as nominal with 2 categories
    if (!is.null(method)) {
      if (method == "FinMix" & model$n_cont > 0) {
        model$categories <- c(rep(2, model$n_cont), model$categories)
        model$diff <- c(rep(1, model$n_cont), model$diff)
        model$n_nom <- model$n_nom + model$n_cont
        model$n_cont <- 0
      }
    
      mats <- getReparmats(model)
    
      if (method == "SpikeSlab") {
        if (family == "gaussian") {
          if (!returnBurnin) {
            mcmc_res <- mcmcSs(y, X = mvars$X_dummy, model, prior, mcmc, mats, returnBurnin)
            mcmc_res_burnin <- NULL
          } else {
            mcmc_res_burnin <- mcmcSs(y, X = mvars$X_dummy, model, prior, mcmc, mats, returnBurnin)
          }
          
        }
        if (family == "binomial") {
          if (!returnBurnin) {
            mcmc_res <- mcmcSsLogit(y, X = mvars$X_dummy, model, prior, mcmc, mats, returnBurnin)
            mcmc_res_burnin <- NULL
          } else {
            mcmc_res_burnin <- mcmcSsLogit(y, X = mvars$X_dummy, model, prior, mcmc, mats, 
                                           returnBurnin)
          }
        }
        if (returnBurnin) {
          prior <- mcmc_res_burnin$prior
          mcmc_res_burnin$prior <- NULL
          mcmc_res <- lapply(mcmc_res_burnin, function(x, burnin) {
            if (is.matrix(x)) 
              return(x[-(1:burnin), ])
            if (is.vector(x)) 
              return(x[-(1:burnin)])
          }, burnin = mcmc$burnin)
          mcmc_res[["prior"]] <- prior
        }
        if (!is.null(modelSelection)) {
          if (modelSelection == "pam") {
            modelSelection <- "binder"
            warning("'method' = 'SpikeSlab' supports only 'modelSelection' = 'binder' or 'modelSelection' = 'NULL'. Default value 'modelSelection' = 'binder' is used.")
          }
          incl_prob <- colMeans(mcmc_res$delta)
          model_sel <- selectModel(incl_prob, model, strategy = "spikeslab_binder")
          refit_res <- modelRefit(model, model_sel, data, mcmcRefit, family)
        }
      }
      if (method == "FinMix") {
        if (family == "gaussian") {
          if (!returnBurnin) {
            mcmc_res <- mcmcMix(y, X = mvars$X_dummy, model, prior, mcmc, returnBurnin)
            mcmc_res_burnin <- NULL
          } else {
            mcmc_res_burnin <- mcmcMix(y, X = mvars$X_dummy, model, prior, mcmc, returnBurnin)
          }
        }
        if (family == "binomial") {
          if (!returnBurnin) {
            mcmc_res <- mcmcMixLogit(y, X = mvars$X_dummy, model, prior, mcmc, returnBurnin)
            mcmc_res_burnin <- NULL
          } else {
            mcmc_res_burnin <- mcmcMixLogit(y, X = mvars$X_dummy, model, prior, mcmc, returnBurnin)
          }
        }
        if (returnBurnin) {
          prior <- mcmc_res_burnin$prior
          mcmc_res_burnin$prior <- NULL
          mcmc_res <- lapply(mcmc_res_burnin, function(x, burnin) {
            if (is.matrix(x)) 
              return(x[-(1:burnin), ])
            if (is.vector(x)) 
              return(x[-(1:burnin)])
          }, burnin = mcmc$burnin)
          mcmc_res[["prior"]] <- prior
        }
        if (!is.null(modelSelection)) {
          if (modelSelection == "binder") {
            model_sel <- selectModel(mcmc_res$S, model, strategy = "finmix_binder")
          }
          if (modelSelection == "pam") {
            incl_prob <- inclProb(S = mcmc_res$S, mvars, model)
            model_sel <- selectModel(incl_prob, model, strategy = "finmix_pam")
          }
          
          refit_res <- modelRefit(model, model_sel, data, mcmcRefit, 
                                  family)
        }
      }
      
      if (method == "FinMix" & sum(types == "c") > 0) {
        cont <- sum(types == "c")
        model$n_cont <- cont
        model$n_nom <- model$n_nom - cont
        model$categories <- model$categories[-c(1:cont)]
        model$diff <- model$diff[-c(1:cont)]
      }
      
      if (is.null(modelSelection)) {
        ret <- list(fit = mcmc_res[names(mcmc_res) != "prior"], fit_burnin = mcmc_res_burnin[names(mcmc_res_burnin) != 
                    "prior"], method = method, family = family, data = list(y = y, X = X_out, X_dummy = mvars$X_dummy, 
                    types = types, levelnames = levelnames), model = model[!names(model) %in% c("lNom", "A_diag", "cov0")], prior = mcmc_res$prior, mcmc = mcmc, mcmcRefit = NULL,
                    modelSelection = modelSelection, returnBurnin = returnBurnin, numbCoef = sum(unique(colMeans(mcmc_res$beta)) != 0), call = cl)
      } else {
        refit_res$model <- model_sel
        ret <- list(fit = mcmc_res[names(mcmc_res) != "prior"], fit_burnin = mcmc_res_burnin[names(mcmc_res_burnin) != 
                    "prior"], refit = refit_res, method = method, family = family, data = list(y = y, 
                    X = X_out, X_dummy = mvars$X_dummy, types = types, levelnames = levelnames), model = model[!names(model) %in% c("lNom", "A_diag", "cov0")], prior = mcmc_res$prior, 
                    mcmc = mcmc, mcmcRefit = mcmcRefit, modelSelection = modelSelection, returnBurnin = returnBurnin, numbCoef = sum(unique(colMeans(refit_res$beta)) != 0), call = cl)
      }
      
    } else {
      if (family == "gaussian") {
        if (!returnBurnin) {
          mcmc_res <- mcmcLinreg(y, mvars$X_dummy, 
                                 prior = list(s0 = 0, S0 = 0, tau2_fix = 1000, conj = FALSE), 
                                 M = mcmc$M, burnin = mcmc$burnin, returnBurnin)
          fit <- mcmc_res[names(mcmc_res) == "beta" | names(mcmc_res) == "sgma2"]
          fit_burnin <- NULL
        } else {
          mcmc_res_burnin <- mcmcLinreg(y, mvars$X_dummy, 
                                 prior = list(s0 = 0, S0 = 0, tau2_fix = 1000, conj = FALSE), 
                                 M = mcmc$M, burnin = mcmc$burnin, returnBurnin)
          fit_burnin <- mcmc_res_burnin[names(mcmc_res_burnin) == "beta" | names(mcmc_res_burnin) == "sgma2"]
        }
      }
      if (family == "binomial") {
        if (!returnBurnin) {
          mcmc_res <- logit(y, mvars$X_dummy, samp = mcmc$M, burn = mcmc$burnin, 
                            P0 = diag(0.1, nrow = ncol(mvars$X_dummy), ncol = ncol(mvars$X_dummy)))
          fit <- mcmc_res[names(mcmc_res) == "beta"]
          fit_burnin <- NULL
        } else {
          mcmc_res_burnin <- logit(y, mvars$X_dummy, samp = mcmc$burnin + mcmc$M, burn = 0,
                                   P0 = diag(0.1, nrow = ncol(mvars$X_dummy), ncol = ncol(mvars$X_dummy)))
          fit_burnin <- mcmc_res_burnin[names(mcmc_res_burnin) == "beta"]
        }
      }
      if (returnBurnin) {
        mcmc_res <- lapply(mcmc_res_burnin, function(x, burnin) {
          if (is.matrix(x)) 
            return(x[-(1:burnin), ])
          if (is.vector(x)) 
            return(x[-(1:burnin)])
        }, burnin = mcmc$burnin)
        if (family == "gaussian") {
          fit <- mcmc_res[names(mcmc_res) == "beta" | names(mcmc_res) == "sgma2"]
        }
        if (family == "binomial") {
          fit <- mcmc_res[names(mcmc_res) == "beta"]
        }
      }
      
      ret <- list(fit = fit, fit_burnin = fit_burnin, method = "No effect fusion performed. Full model was estimated.", 
                  family = family, data = list(y = y, X = X_out, X_dummy = mvars$X_dummy, types = types, levelnames = levelnames), 
                  model = model[!names(model) %in% c("lNom", "A_diag", "cov0")], prior = "A flat, uninformative prior was used for model fitting.",
                  mcmc = mcmc[names(mcmc) != "startsel"], mcmcRefit = NULL, modelSelection = NULL, returnBurnin = returnBurnin, numbCoef = sum(unique(colMeans(fit$beta)) != 0), call = cl)
    }
    
    class(ret) <- "fusion"
    return(ret)
    
}
