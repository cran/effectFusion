
createRowNames <- function(model, catlabels, contnames) {
    var <- rep(1:(model$n_nom + model$n_ord), model$categories - 1) + model$n_cont
    cat <- sequence(model$categories - 1) + 1
    catnames <- unlist(lapply(catlabels, function(x) x[-1]))
    names <- c()
    for (i in 1:length(cat)) {
        names[i] <- paste(names(catlabels)[var[i] - model$n_cont], ".", catnames[i], sep = "")
    }
    if (model$n_cont > 0) {
        cont_names <- contnames
        names <- c(cont_names, names)
    }
    names <- c("(Intercept)", names)
    return(names)
}
