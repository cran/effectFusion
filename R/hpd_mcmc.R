
hpdMCMC = function (x, alpha = 0.05, ...) 
{
  if (!is.vector(x)) 
    x <- as.vector(x)
  nr <- NROW(x)
  l <- ceiling(alpha * nr)
  lb <- sort(x)[1:l]
  ub <- sort(x)[(nr - l + 1):nr]
  minl <- min(ub - lb)
  lower <- lb[which((ub - lb) == minl)][1]
  upper <- ub[which((ub - lb) == minl)][1]
  return(c(lower = lower, upper = upper))
}


