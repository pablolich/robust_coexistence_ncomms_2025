library(sn) # skewed normal distribution
library(mvtnorm) # multinormal distribution

# This is a discrete distribution having
# support only at a (with prob p) and b 
# a > 1 makes it right-skewed
# a < 1 makes it left-skewed
# we set a = sk, b = -1/sk
# setting a = 1 returns a matrix 
# with zero skewness and entries +/-1 
build_discrete <- function(n, rho, sk){
  a <- sk
  b <- -1/sk
  p <- 1 / (1 + a^2)
  M <- matrix(runif(n * n), n, n)
  M[M <= p] <- a
  M[M != a] <- b
  diag(M) <- 0
  if (rho == 1){
    M[lower.tri(M)] <- 0
    M <- M + t(M)
  }
  if (rho == -1){
    M[lower.tri(M)] <- 0
    M <- M - t(M)
  }
  return(M)
}

# This is a continuous distribution having
# finite support
# sk > 1 makes it right-skewed
# sk < 1 makes it left-skewed
# setting sk = 1 returns entries 
# uniformly distributed with mean 0 and var 1
build_beta <- function(n, rho, sk){
  alpha <- 1 / sk
  beta <- sk
  coeffs <- rbeta(n * n, alpha, beta)
  # subtract mean
  coeffs <- coeffs - 1 / (1 + sk^2)
  # standardize
  coeffs <- coeffs / sqrt(
    sk^3 / ((1 + sk^2)^2  * (1 + sk + sk^2))
  )
  M <- matrix(coeffs, n, n)
  diag(M) <- 0
  if (rho == 1){
    M[lower.tri(M)] <- 0
    M <- M + t(M)
  }
  if (rho == -1){
    M[lower.tri(M)] <- 0
    M <- M - t(M)
  }
  return(M)
}

# This is a continuous distribution having
# infinite support
# sk > 0 makes it right-skewed
# sk < 0 makes it left-skewed
# to have two comparable distributions
# call the function with sk = x and sk = -x
# we handle the case of sk = 0 and arbitrary rho
# separately (multinormal distribution)
build_normal <- function(n, rho, sk){
  if (sk == 0) return(build_multinormal(n, rho))
  # this controls the mean
  my_xi <- -((sqrt(2) * sk)/sqrt(pi + (-2 + pi) * sk^2))
  # this controls the standard deviation
  my_omega <- sqrt(pi) * sqrt(1 + sk^2) / sqrt(pi + (-2 + pi) * sk^2)
  coeffs <- rsn(n * n, xi = my_xi, omega = my_omega, alpha = sk)
  M <- matrix(coeffs, n, n)
  diag(M) <- 0
  if (rho == 1){
    M[lower.tri(M)] <- 0
    M <- M + t(M)
  }
  if (rho == -1){
    M[lower.tri(M)] <- 0
    M <- M - t(M)
  }
  return(M)
}

# build a matrix with interactions sampled in pairs
# from bivariate normal
build_multinormal <- function(n, rho){
  M <- matrix(rnorm(n * n, 0, 1), n, n)
  MijMji <-
    mvtnorm::rmvnorm(n * (n - 1) / 2, c(0, 0), matrix(c(1, rho, rho, 1), 2, 2, byrow = TRUE))
  M[upper.tri(M)] <- MijMji[, 1]
  M <- t(M)
  M[upper.tri(M)] <- MijMji[, 2]
  diag(M) <- 0
  return(M)
}

build_matrix <- function(n, dist, rho, sk){
  # build the matrix
  if (dist == "normal") M <- build_normal(n, rho, sk)
  if (dist == "beta") M <- build_beta(n, rho, sk)
  if (dist == "discrete") M <- build_discrete(n, rho, sk)
  return(M)
}
