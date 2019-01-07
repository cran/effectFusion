
coding <- function(x, constant = TRUE, splitcod = TRUE) {
    x <- as.matrix(x)
    n <- nrow(x)
    kx <- apply(x, 2, max) - 1
    xds <- matrix(0, n, sum(kx))
    for (j in 1:ncol(x)) {
        j1col <- ifelse(j > 1, sum(kx[1:(j - 1)]), 0)
        for (i in 1:n) {
            if (x[i, j] > 1) {
                if (splitcod) 
                  xds[i, j1col + 1:(x[i, j] - 1)] <- 1 else xds[i, j1col + (x[i, j] - 1)] <- 1
            }
        }
    }
    if (constant) 
        return(cbind(1, xds)) else return(xds)
}
