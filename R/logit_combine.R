## Copyright 2016 Nick Polson, James Scott, and Jesse Windle.
## Extracted from the R package 'BayesLogit' (Version: 0.6).
#' @useDynLib 'effectFusion', .registration = TRUE, .fixes = 'C_'

## Combine
logit.combine <- function(y, X, n = rep(1, length(y))) {
    X <- as.matrix(X)
    
    N <- dim(X)[1]
    P <- dim(X)[2]
    
    m0 <- matrix(0, nrow = P)
    P0 <- matrix(0, nrow = P, ncol = P)
    
    ok <- check.parameters(y, n, m0, P0, N, P, 1, 0)
    if (!ok) 
        return(-1)
    
    ## Our combine_data function, written in C, uses t(X).
    tX <- t(X)
    
    OUT <- .C(C_combine, as.double(y), as.double(tX), as.double(n), as.integer(N), as.integer(P))
    
    N <- OUT[[4]]
    
    y <- array(as.numeric(OUT[[1]]), dim = c(N))
    tX <- array(as.numeric(OUT[[2]]), dim = c(P, N))
    n <- array(as.numeric(OUT[[3]]), dim = c(N))
    
    list(y = as.numeric(y), X = t(tX), n = as.numeric(n))
}
