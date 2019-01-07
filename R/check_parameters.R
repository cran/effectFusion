## Copyright 2016 Nick Polson, James Scott, and Jesse Windle.
## Extracted from the R package 'BayesLogit' (Version: 0.6).

## Check parameters to prevent an obvious error.
check.parameters <- function(y, n, m0, P0, R.X, C.X, samp, burn) {
    ok <- rep(TRUE, 8)
    ok[1] <- all(y >= 0)
    ok[2] <- all(n > 0)
    ok[3] <- C.X == nrow(P0)
    ok[4] <- C.X == ncol(P0)
    ok[5] <- (length(y) == length(n) && length(y) == R.X)
    ok[6] <- C.X == length(m0)
    ok[7] <- (samp > 0)
    ok[8] <- (burn >= 0)
    ok[9] <- all(y <= 1)
    
    if (!ok[1]) 
        print("y must be >= 0.")
    if (!ok[9]) 
        print("y is a proportion; it must be <= 1.")
    if (!ok[2]) 
        print("n must be > 0.")
    if (!ok[3]) 
        print(paste("col(X) != row(P0)", C.X, nrow(P0)))
    if (!ok[4]) 
        print(paste("col(X) != col(P0)", C.X, ncol(P0)))
    if (!ok[5]) 
        print(paste("Dimensions do not conform for y, X, and n.", "len(y) =", length(y), "dim(x) =", 
            R.X, C.X, "len(n) =", length(n)))
    if (!ok[6]) 
        print(paste("col(X) != length(m0)", C.X, length(m0)))
    if (!ok[7]) 
        print("samp must be > 0.")
    if (!ok[8]) 
        print("burn must be >=0.")
    
    ok <- all(ok)
}
