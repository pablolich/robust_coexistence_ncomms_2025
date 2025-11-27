library(tidyverse)
# input needed:
# alpha_values, n, rho, skewness par
# expectation for polynomials of the random variables
# mus[k] = E[Bij^k] (needed up to k = 6)
# it is assumed that E[Bij] = mus[1] = 0, E[Bij^2] = mus[2] 1

# these are the values for (Bij, Bji) ~ N((0,0), ((1,rho), (rho,1)))
get_expectations_multinormal <- function(rho){
  mus <- c(0, 1, 0, 3, 0, 15)  
  rhos <- c(rho, 1 + 2 * rho^2, 3 * rho * (3 + 6 * rho^2))
  return(list(mus = mus, rhos = rhos))
}

# moments for normal and skewed normal
get_expectations_normal <- function(rho, sk){
  if (sk == 0) return(get_expectations_multinormal(rho))
  mus <- c(
    0, 
    1, 
    -((sqrt(2) * (-4 + pi) * sk^3)/(pi + (-2 + pi) * sk^2)^(3/2)), 
    (-12 * sk^4 + 3 * pi^2 * (1 + sk^2)^2 - 4 * pi * sk^2 * (3 + sk^2))/(pi + (-2 + pi) * sk^2)^2, 
    (sqrt(2) * sk^3 * (-10 * (-4 + pi) * pi + (16 + (20 - 7 * pi) * pi) * sk^2))/(pi + (-2 + pi) * sk^2)^(5/2), 
    (-40 * sk^6 + 15 * pi^3 * (1 + sk^2)^3 - 20 * pi * sk^4 * (9 + 5 * sk^2) - 6 * pi^2 * sk^2 * (15 + 10 * sk^2 + sk^4))/(pi + (-2 + pi) * sk^2)^3
  )
  return(list(mus = mus, rhos = NULL))
}

get_expectations_beta <- function(rho, sk){
  mus <- c(
    0, 
    1, 
    (2 * (-1 + sk) * sqrt(1 + 1/sk + sk))/(1 + sk), 
    -21 + 6/sk + sk *(6 + 54/(1 + sk *(3 + sk))), 
    (4 *(-1 + sk) * (1 + 1/sk + sk)^(3/2) * (6 + 5 * sk + 5 * sk^3 + 6 * sk^4))/((1 + sk) * (1 + sk * (3 + sk)) * (1 + sk * (4 + sk))), 
    (5 * (1 + sk + sk^2)^2 * (24 + sk * (-22 + sk * (-1 + sk * (16 + sk * (-1 - 22 * sk + 24 * sk^2))))))/(sk^2 * (1 + sk * (3 + sk)) * (1 + sk * (4 + sk)) * (1 + sk * (5 + sk)))
  )
  return(list(mus = mus, rhos = NULL))
}

get_expectations_twovals <- function(rho, sk){
  mus <- c(
    0, 
    1, 
    -(1/sk) + sk, 
    -1 + 1/sk^2 + sk^2, 
    ((-1 + sk) * (1 + sk) * (1 + sk^4))/sk^3, 
    1 + 1/sk^4 - 1/sk^2 - sk^2 + sk^4 
  )
  return(list(mus = mus, rhos = NULL))
}

get_variance_symmetric <- function(alpha, n, m){
  Mu2 <- m$mus[2]; Mu3 <- m$mus[3]; Mu4 <- m$mus[4]; Mu5 <- m$mus[5]; Mu6 <- m$mus[6];
  Alpha <- alpha
  Omega <- alpha^2 / (alpha^2 -2 * n)
  variances <-(
    ((-2 + n) * (-1 + n) * Mu2) / (n * Alpha^4) +
      ((-2 + n) * (-1 + n) * 
         (-432 + n * (724 + n * (-392 + n * (105 + n * (-17 + 2 * n))))) * Mu2^3 * Omega^2) / 
      (n^4 * Alpha^8) +
      (2 * (-2 + n) * (-1 + n) * 
         (64 + n * (-92 + n * (56 + (-14 + n) * n))) * Mu3^2 * Omega^2) / 
      (n^4 * Alpha^8) -
      (2 * (-2 + n)^4 * (-1 + n) * Mu5 * Omega^2) / 
      (n^4 * Alpha^7) +
      ((-2 + n)^4 * (-1 + n) * Mu6 * Omega^2) / 
      (n^4 * Alpha^8) +
      ((-2 + n) * (-1 + n) * 
         (-12 + n * (18 + (-5 + n) * n)) * Mu2^2 * Omega * (2 + Omega)) / 
      (n^3 * Alpha^6) +
      Mu3 * (
        -((4 * (-2 + n)^2 * (-1 + n) * Omega) / (n^2 * Alpha^5)) - 
          (4 * (-2 + n) * (-1 + n) * 
             (40 + n * (-88 + n * (49 + (-12 + n) * n))) * Mu2 * Omega^2) / 
          (n^4 * Alpha^7)
      ) +
      Mu4 * (
        (2 * (-2 + n) * (-1 + n) * 
           (92 + n * (-140 + n * (76 + n * (-17 + 2 * n)))) * Mu2 * Omega^2) / 
          (n^4 * Alpha^8) +
          ((-2 + n)^3 * (-1 + n) * Omega * (2 + Omega)) / 
          (n^3 * Alpha^6)
      ))
  return(variances)
}

get_variance_skew_symmetric <- function(alpha, n, m){
  Mu2 <- m$mus[2]; Mu3 <- m$mus[3]; Mu4 <- m$mus[4]; Mu5 <- m$mus[5]; Mu6 <- m$mus[6];
  Alpha <- alpha
  Omega <- alpha^2 / (alpha^2 + 2 * n)
  variances <-  (1 / (3 * n^2 * Alpha^8)) * (-1 + n) * (
    3 * n^2 * Alpha^4 * Mu2 -
      6 * (-2 + n) * n * Alpha^2 * ((-3 + n) * Mu2^2 + Mu4) * Omega +
      (-2 + n) * (
        6 * n^3 * Mu2^3 +
          28 * Mu3^2 +
          90 * Mu2 * (-2 * Mu2^2 + Mu4) +
          n^2 * (2 * Mu3^2 + 3 * Mu2 * ((Alpha^2 - 15 * Mu2) * Mu2 + 4 * Mu4)) -
          6 * Mu6 +
          n * (135 * Mu2^3 - 16 * Mu3^2 - 60 * Mu2 * Mu4 + 
                 3 * Alpha^2 * (-3 * Mu2^2 + Mu4) + 3 * Mu6)
      ) * Omega^2
  )
  return(variances)
}

get_variance_independent <- function(alpha, n, m){
  Mu2 <- m$mus[2]; Mu3 <- m$mus[3]; Mu4 <- m$mus[4]; Mu5 <- m$mus[5]; Mu6 <- m$mus[6];
  Alpha <- alpha
  Omega <- 1
  variances <-(
    ((-1 + n) * (n^3 * (Alpha^2 - 3 * Mu2) * Mu2 +
                   3 * Mu2^2 + 
                   n^4 * Mu2^2 + 
                   n^2 * (-Alpha^2 * Mu2 + 4 * Mu2^2 + 2 * Alpha * Mu3) - 
                   Mu4 + 
                   n * (-6 * Mu2^2 - 2 * Alpha * Mu3 + Mu4))) /
      (n^3 * Alpha^6)
  )
  return(variances)
}

# compute the variance approximations 
variance_approximations <- function(alpha_values = c(1,2,3), n = 10, rho = 0, dist = "normal", sk = 2){
  # get the expectations for the polynomials of random variables
  if (dist == "normal") m <- get_expectations_normal(rho, sk)
  if (dist == "beta") m <- get_expectations_beta(rho, sk)
  if (dist == "discrete") m <- get_expectations_twovals(rho, sk)
  if (rho == 0) variances <- get_variance_independent(alpha = alpha_values, n = n, m = m)
  if (rho == 1) variances <- get_variance_symmetric(alpha = alpha_values, n = n, m = m)
  if (rho == -1) variances <- get_variance_skew_symmetric(alpha = alpha_values, n = n, m = m)
  tb <- tibble(alpha = alpha_values, approx = "Richardson", values = variances)
  return(tb %>% arrange(alpha))
}

# compute the variance approximations 
pf_approximations <- function(alpha_values = c(1,2,3), n = 10, rho = 0, dist = "normal", sk = 2){
  variances <- variance_approximations(alpha_values = alpha_values, n = n, rho = rho, dist = dist, sk = sk)$values
  critvals <- 1 / (alpha_values * sqrt(variances))
  tb <- tibble(alpha = alpha_values, approx = "Richardson", values = pnorm(critvals)^n)
  return(tb %>% arrange(alpha))
}