library(tidyverse)
# requires to have loaded both
# build_matrices.R
# find_alpha_F.R
output_dir <- "results" # put here output directory


run_simulation <- function(n, dist, rho, sk, nsim){
  results <- tibble(n = rep(n, nsim), dist = dist, rho = rho, sk = sk, alpha_F = 0)
  for (i in 1:nsim){
    print(i)
    B <- build_matrix(n, dist, rho, sk)
    results$alpha_F[i] <- get_alpha_F_large_random_matrix(B, n, rho)
  }
  write_csv(results, file = paste0(output_dir, "/n_", n, "_rho_", rho, "_dist_", dist, "_sk_", sk, ".csv"))
}