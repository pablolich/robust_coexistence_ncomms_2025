library(tidyverse)
plot_approximations <- function(filename) {
  dt <- read_csv(filename)
  n <- dt$n[1]
  rho <- dt$rho[1]
  dist <- dt$dist[1]
  sk <- dt$sk[1]
  # observed values of alpha_F
  range_alpha <- range(dt$alpha_F)
  alpha_values <-
    seq(range_alpha[1], range_alpha[2], length.out = 500)
  # plot observed values
  pl <-
    ggplot(dt, aes(x = alpha_F)) + stat_ecdf() + theme_bw() + ylab("p_F")
  approximations <- tibble()
  # Richardson approximation
  if (rho %in% c(0, 1,-1)) {
    source("code/approximation_richardson.R")
    approximations <- bind_rows(approximations,
                                pf_approximations(alpha_values, n, rho, dist, sk))
  }
  source("code/approximation_small_rank.R")
  approximations <- bind_rows(approximations,
                              pf_approximations(alpha_values, n, rho, dist, sk))
  
  source("code/approximation_resolvent.R")
  approximations <- bind_rows(approximations,
                              pf_approximations(alpha_values, n, rho, dist, sk))
  
  
  
  pl <-
    pl + geom_line(data = approximations, aes(x = alpha, y = values, colour = approx))
  return(pl)
}