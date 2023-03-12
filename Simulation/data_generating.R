######################################################################################################
# In this file we simulate our data,
# 
# Observed data:
# T0: censored event time
# Delta: event indicator, event = 0 censored, event time is less than or equal to censoring time
# A: treatment indicator
# X: observed covariates, type 1 proxy
# W: type 2 proxy, not caused by A
# Z: type 3 proxy, not causing Y, W
# 
# Unobserved data:
# U: unmeasured confounders, correlated with Z and W, causing A and Y
# true event time and censoring time are buried
#######################################################################################################

para_set <- list(mu_X = 1.1,
                 sigma_X = 0.75,
                 mu_U = 1.1,
                 sigma_U = 0.75,
                 alpha_A = c(0.3, 0.4, -0.6),
                 mu_Z = c(-0.2, -0.3, 0.65),
                 sigma_Z = 0.5,
                 mu_W = c(-0.6, 0.4, 0.65),
                 sigma_W = 0.5,
                 mu_T0 = c(0.1, 0.6, 0.25, 0.5),
                 mu_C = 0.2,
                 admin_C = 2
)



true_b <- function(para_set) {
  bw <- para_set$mu_T0[4]/para_set$mu_W[3]
  bx <- para_set$mu_T0[3] - bw * para_set$mu_W[2]
  b1 <- bw^2 * para_set$sigma_W^2/2
  ba <- para_set$mu_T0[2] - bw * 0
  b0 <- para_set$mu_T0[1] - bw * para_set$mu_W[1]
  print("b0, ba, bw, bx, b1")
  b <- c(b0, ba, bw, bx, b1)
  return(b)
}

true_t <- function(para_set) {
  tz <- -para_set$alpha_A[3]/para_set$mu_Z[3]
  tx <- -para_set$alpha_A[2] - tz * para_set$mu_Z[2]
  ta <- 0 - tz * 0 - tz^2 * para_set$sigma_Z^2
  t0 <- -para_set$alpha_A[1] - tz * para_set$mu_Z[1] + tz^2 * para_set$sigma_Z^2/2
  t <- c(t0, ta, tz,tx)
  print("t0, ta, tz, tx")
  return(t)
}


data_gen <- function(N, para_set, a = NULL) {
  # generate X, U
  X <- para_set$mu_X + rnorm(N, 0, para_set$sigma_X)
  U <- para_set$mu_U + rnorm(N, 0, para_set$sigma_U)
  X <- pmax(X, 0)
  U <- pmax(U, 0)
  
  if (is.null(a)) {
    # generate A
    prop_score_0 <- 1/(1 + exp(-cbind(1, X, U) %*% para_set$alpha_A))
    A <- rbinom(N, 1, prop_score_0)
  } else {
    A <- rep(a, N)
  }
  
  
  # generate Z
  Z <- cbind(1, X, U) %*% para_set$mu_Z + rnorm(N, 0, para_set$sigma_Z)
  
  # generate W
  W <- cbind(1, X, U) %*% para_set$mu_W + rnorm(N, 0, para_set$sigma_W)
  
  
  #generate Y
  T0 <- rexp(N, rate = cbind(1, A, X, U) %*% para_set$mu_T0)
  
  C <- rexp(N, rate = para_set$mu_C)
  C <- pmin(C, para_set$admin_C)
  if (is.null(a)) {
    df <- data.frame(X, U, A, Z, W, T0 = pmin(T0, C), Delta = (T0 <= C))
  } else {
    df <- data.frame(X, U, A, Z, W, T0 = T0, Delta = rep(1, N))
  }
  return(df)
}


