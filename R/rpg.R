## Copyright 2016 Nick Polson, James Scott, and Jesse Windle.
## Extracted from the R package 'BayesLogit' (Version: 0.6).
#' @useDynLib 'effectFusion', .registration = TRUE, .fixes = 'C_'

## Draw PG(n, z)
rpg <- function(num = 1, h = 1, z = 0) {
    ## Check Parameters.
    if (any(h <= 0)) {
        print("h must be > 0.")
        return(NA)
    }
    
    x <- rep(0, num)
    
    if (length(h) != num) {
        h <- array(h, num)
    }
    if (length(z) != num) {
        z <- array(z, num)
    }
    
    ## Faster if we do not track iter.
    OUT <- .C(C_rpg_hybrid, x, h, z, as.integer(num))
    
    OUT[[1]]
}
