library(tidyverse)

variance_approximations <- function(alpha_values = c(1,2,3), n = 10, rho = 0, dist = "normal", sk = 2){
  # get the expectations for the polynomials of random variables
  a <- alpha_values 
  inner_value <- a^2 - 4 * rho * (n - 1)
  inner_value[inner_value < 0]  <- 0
  vars <- 1/(a^2 * ((a^2 + a * sqrt(inner_value)) / (2 * (n - 1)) - (1 + rho)))
  tb <- tibble(alpha = alpha_values, approx = "resolvent", values = vars)
  return(tb %>% arrange(alpha))
}

# compute the variance approximations 
pf_approximations <- function(alpha_values = c(1,2,3), n = 10, rho = 0, dist = "normal", sk = 2){
  variances <- variance_approximations(alpha_values = alpha_values, n = n, rho = rho, dist = dist, sk = sk)$values
  critvals <- 1 / (alpha_values * sqrt(variances))
  tb <- tibble(alpha = alpha_values, approx = "resolvent", values = pnorm(critvals)^n)
  return(tb %>% arrange(alpha))
}