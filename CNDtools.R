# The code in this file allows for evaluations of the cdf and quantile functions
# for Canonical Noise Distributions (CNDs), as well as sampling functions for 
# the CNDs

# Code corresponds to the paper:
# Optimizing Noise for f-Differential Privacy
# via Anti-Concentration and Stochastic Dominance
# by Jordan Awan and Aishwarya Ramasethu

# Note that the tradeoff function f is defined as the type II error as a 
# function of the 1 - type I error

# Numerically calculate c_f if there is no closed form
# f need not be vectorized, but is faster if it is vectorized.
c_f <- function(f){
  obj = function(c){
    return((f(1-c)-c)^2)
  }
  opt = optim(par=.5,fn=obj,lower=0,upper=.5,method="Brent")
  return(opt$par)
}
# quantile function for the CND for f
qCND_base <- function(u, f, c) {
  
  if (u < c) {
    return(qCND(1 - f(1-u), f, c) - 1)
  } else if (u >= c && u <= 1 - c) {
    return((u - 0.5) / (1 - 2 * c))
  } else {
    return(qCND(f(u), f, c) + 1)
  }
}
# vectorized quantile function
qCND <- function(u,f,c){
  return(sapply(X=u,FUN=qCND_base,f=f,c=c))
}

# CDF of the CND for f
pCND_base <- function(x, c, f) {
  if (x < -1/2) {
    return(f(pCND_base(x + 1, c, f)))
  } else if (x >= -1/2 && x <= 1/2) {
    return(c * (1/2 - x) + (1 - c) * (x + 1/2))
  } else {  # x > 1/2
    return(1 - f(1-pCND_base(x - 1, c, f)))
  }
}
# vectorized cdf
pCND <- function(x,c,f){
  return(sapply(X=x,FUN=pCND_base,c=c,f=f))
}
# obtain n i.i.d. samples from the CND for f
rCND <- function(n, f, c) {
  u <- runif(n)
  return(qCND(u,f,c))
}
# obtain n i.i.d. samples from a discrete CND for f at sensitivity Delta
rDiscCND <- function(n, f, c, delta) {
  random_samples <- rCND(n, f, c)
  rounded_samples <- round(random_samples * delta)
  return(rounded_samples)
}