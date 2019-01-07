## Copyright 2016 Nick Polson, James Scott, and Jesse Windle.
## Extracted from the R package 'BayesLogit' (Version: 0.6).
#' @useDynLib 'effectFusion', .registration = TRUE, .fixes = 'C_'

## Posterior by Gibbs
logit <- function(y, X, n = rep(1, length(y)), m0 = rep(0, ncol(X)), P0 = matrix(0, nrow = ncol(X), 
    ncol = ncol(X)), samp = 1000, burn = 500) {
    ## In the event X is one dimensional.
    X <- as.matrix(X)
    
    ## Combine data.  We do this so that the auxiliary variable matches the data.
    new.data <- logit.combine(y, X, n)
    y <- new.data$y
    X <- new.data$X
    n <- new.data$n
    
    ## Check that the data and priors are okay.
    N <- dim(X)[1]
    P <- dim(X)[2]
    
    ok <- check.parameters(y, n, m0, P0, N, P, samp, burn)
    if (!ok) 
        return(-1)
    
    ## Initialize output.
    output <- list()
    
    ## w = array(known.w, dim=c(N, samp)); beta = array(known.beta, dim=c(P , samp));
    w <- array(0, dim = c(N, samp))
    beta <- array(0, dim = c(P, samp))
    
    ## our Logit function, written in C, uses t(X).
    tX <- t(X)
    
    OUT <- .C(C_gibbs, w, beta, as.double(y), as.double(tX), as.double(n), as.double(m0), as.double(P0), 
        as.integer(N), as.integer(P), as.integer(samp), as.integer(burn))
    
    N <- OUT[[8]]
    
    tempw <- array(as.numeric(OUT[[1]]), dim = c(N, samp))
    
    output <- list(w = t(tempw), beta = t(OUT[[2]]), y = y, X = X, n = n)
    
    output
}
