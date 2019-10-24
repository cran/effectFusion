#' Plot an object of class \code{fusion}
#'
#' @description This function provides plots of posterior means and 95\%-HPD intervals for the
#' regression effects. Plots are based on the refitted MCMC samples of the
#' selected model in an object of class \code{fusion} or on model averaged results if no
#' final model selection was performed.
#'
#' @param x an object of class \code{fusion}
#' @param maxPlots maximum number of plots on a single page, default argument to 4
#' @param ... further arguments passed to or from other methods (not used)
#'
#' @details If no effect fusion or no final model selection is performed, posterior means and HPD intervals are model averaged
#' results. Otherwise, the parameters of the selected model are reestimated
#' with a flat uninformative prior. Thus, fused categories have the same posterior mean and HPD interval. Single categories that are excluded
#' from the model are fused to the reference category and therefore only a posterior mean at zero and no
#' interval is plotted.
#'
#' @author Daniela Pauger, Magdalena Leitner <effectfusion.jku@gmail.com>
#'
#' @method plot fusion
#' @export
#'
#' @seealso \code{\link{effectFusion}}
#'
#' @examples
#' ## see example for effectFusion

#' @importFrom gridExtra grid.arrange

plot.fusion <- function(x, maxPlots = 4, ...) {
    
    stopifnot(is(x, "fusion"))
    
    nVar <- sum(x$model$n_nom, x$model$n_ord, x$model$n_cont)
    nPlots <- nVar - x$model$n_cont
    upper <- cumsum(c(1 + x$model$n_cont, x$model$categories - 1))[-1]
    lower <- cumsum(c(2 + x$model$n_cont, x$model$categories - 1))
    levelnames <- lapply(x$data$levelnames, function(x) {x[-1]})
    
    
    if (nPlots <= maxPlots) {
        p <- list()
        
        for (i in 1:nPlots) {
            j <- i + x$model$n_cont
            plot_name <- paste("Covariate '", names(x$data$levelnames)[i], "'", sep = "")
            if (is.null(x$modelSelection)) {
                pH <- plotHPD(x$fit$beta[, lower[i]:upper[i]], title = plot_name, labels = levelnames[[i]])
            } else {
                pH <- plotHPD(x$refit$beta[, lower[i]:upper[i]], title = plot_name, labels = levelnames[[i]])
            }
            p[[i]] <- pH
        }
        
        nCol <- floor(sqrt(length(p)))
        do.call("grid.arrange", c(p, ncol = nCol))
        
    } else {
        
        nPages <- ceiling(nPlots/maxPlots)
        j <- 1 + x$model$n_cont
        i <- 1
        
        for (n in 1:nPages) {
            p <- list()
            if (n == nPages && (nPlots%%maxPlots > 0)) 
                maxPlots <- nPlots%%maxPlots
            for (k in 1:maxPlots) {
                plot_name <- paste("Covariate '", names(x$data$levelnames)[i], "'", sep = "")
                j <- j + 1
                if (is.null(x$modelSelection)) {
                  pH <- plotHPD(x$fit$beta[, lower[i]:upper[i]], title = plot_name, labels = levelnames[[i]])
                } else {
                  pH <- plotHPD(x$refit$beta[, lower[i]:upper[i]], title = plot_name, labels = levelnames[[i]])
                }
                i <- i + 1
                p[[k]] <- pH
            }
            nCol <- floor(sqrt(length(p)))
            if (n > 1) {
                oask <- grDevices::devAskNewPage(ask = TRUE)
                on.exit(grDevices::devAskNewPage(oask))
            }
            do.call("grid.arrange", c(p, ncol = nCol))
        }
    }
    
}