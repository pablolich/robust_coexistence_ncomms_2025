# compute the variance approximations 
pf_approximations <- function(alpha_values = c(1,2,3), n = 10, rho = 0, dist = "normal", sk = 2){
  atilde <- alpha_values / sqrt(n-1)
  bn <- sqrt(pracma::lambertWp(n^2 / (2 * pi)))
  crit_vals <- ((atilde + sqrt(atilde^2 - 4 * rho))/2 - bn) * bn
  tb <- tibble(alpha = alpha_values, 
               approx = "small-rank", 
               values = exp(-exp(-crit_vals)))
  return(tb %>% arrange(alpha))
}