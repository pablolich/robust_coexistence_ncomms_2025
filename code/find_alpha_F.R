THRESH_IMAGINARY <- 10^-12 # call imaginary part 0 if it falls below this threshold

is_feasible <- function(alpha, B, ones){
  # solve the linear system
  y <- solve(diag(ones * alpha) + B, ones)
  min_y <- min(y)
  if (min_y < 0) {
    k <- which.min(y)
    return(k)
  }
  return(0)
}

update_alpha <- function(k, B, ones, old_alpha){
  b <- B[k, -k]
  Btilde <- B[-k, -k]
  tmp <- rARPACK::eigs(A = ones[-1] %o% b - Btilde, 
                       k = 1, which = "LR", opts = list(retvec = FALSE))
  new_alpha <- Re(tmp$values)
  if (new_alpha > old_alpha) return(new_alpha)
  return(1.01 * old_alpha)
}

get_alpha_F_large_random_matrix <- function(B, n, rho){
  C <- B - rep(1, n) %o% colMeans(B)
  # get initial guess
  L1 <- log(n^2 / (2 * pi))
  L2 <- log(L1)
  bn <- sqrt(L1 - L2 + L2/L1)
  # starting guess
  x_guess <- 0
  gamma_guess <- x_guess / bn + bn
  alpha_guess <- (gamma_guess + rho / gamma_guess) * sqrt(n-1)
  ones <- rep(1, n)
  # make sure it's not already feasible
  while(is_feasible(alpha_guess, C, ones) == 0) alpha_guess <- 0.95 * alpha_guess
  new_alpha <- alpha_guess
  k <- is_feasible(new_alpha, C, ones)
  while(k != 0){
    new_alpha <- update_alpha(k, C, ones, new_alpha)
    k <- is_feasible(new_alpha, C, ones)
  }
  return(new_alpha)
}

get_alpha_F <- function(B, n, rho){
  # form matrix C by removing col means
  C <- B - rep(1, n) %o% colMeans(B)
  # compute alpha_0
  eC <- eigen(-C, only.values = TRUE)$values
  eC <- Re(eC[abs(Im(eC)) < THRESH_IMAGINARY])
  alpha_0 <- max(eC)
  # compute alpha_i
  alpha_i <- rep(0, n)
  for (i in 1:n){
    Qi <- rep(1, n-1) %o% C[i,-i] - C[-i,-i]
    eQi <- eigen(Qi, only.values = TRUE)$values
    eQi <- Re(eQi[abs(Im(eQi)) < THRESH_IMAGINARY])
    alpha_i[i] <- max(eQi)
  }
  return(max(c(alpha_0, alpha_i)))
}
