censor_cumhaz_est <- function(censor_sorted_time, T0, Delta) {
  K_C <- length(censor_sorted_time)
  N <- length(T0)
  ## q value, N * 1
  IF <- matrix(0, nrow = K_C, ncol = N)
  for(j in 1:K_C) {
    IF[j, ] <- (1 - Delta) * (T0 == censor_sorted_time[j])/sum(T0 >= censor_sorted_time[j])
  }
  IF <- apply(IF, 2, cumsum)
  return(IF)
}


noncensor_cumhaz_compute <- function(sorted_time, censor_sorted_time, cum_haz) {
  noncensor_cumhaz <- c()
  K = length(sorted_time)
  K_C = length(censor_sorted_time)
  censor_sorted_time <- c(0, censor_sorted_time)
  cum_haz <- c(0, cum_haz)
  l = 1
  j = 1
  while(j <= K) {
    if(l + 1 <= K_C) {
      if(sorted_time[j] < censor_sorted_time[l + 1])
        noncensor_cumhaz <- c(noncensor_cumhaz, cum_haz[l])
      else {
        l = l + 1
        next
      }
    } else {
      noncensor_cumhaz <- c(noncensor_cumhaz, cum_haz[l])
    }
    j = j + 1
  }
  return(noncensor_cumhaz)
}

noncensor_cumhaz_IF_compute <- function(sorted_time, censor_sorted_time, cum_haz_IF) {
  noncensor_cumhaz_IF <- c()
  K = length(sorted_time)
  K_C = length(censor_sorted_time)
  censor_sorted_time <- c(0, censor_sorted_time)
  cum_haz_IF <- rbind(0, cum_haz_IF)
  l = 1
  j = 1
  while(j <= K) {
    if(l + 1 <= K_C) {
      if(sorted_time[j] < censor_sorted_time[l + 1])
        noncensor_cumhaz_IF <- rbind(noncensor_cumhaz_IF, cum_haz_IF[l, ])
      else {
        l = l + 1
        next
      }
    } else {
      noncensor_cumhaz_IF <- rbind(noncensor_cumhaz_IF, cum_haz_IF[l, ])
    }
    j = j + 1
  }
  return(noncensor_cumhaz_IF)
}



hbridge <- function(b, h_W, h_Z, T0, sorted_time, noncensor_cumhaz) { 
  K <- length(sorted_time)
  g <- 0
  
  for(j in 1:K) {
    res <- cbind(h_Z, sorted_time[j]) / K
    integrand <- (T0 > sorted_time[j]) / exp(-noncensor_cumhaz[j]) - exp(-cbind(h_W, sorted_time[j]) %*% b * sorted_time[j])
    g <- g + c(integrand) * res
  }
  g <- colMeans(g)
  gmmf <- sum(g^2)
  return(gmmf)
}


h_IF <- function(b, h_W, h_Z, T0, sorted_time, noncensor_cumhaz, noncensor_cumhaz_IF) {
  hbridge_tmp <- function(b, h_W, h_Z, T0, sorted_time, noncensor_cumhaz) { 
    K = length(sorted_time)
    g <- 0
    
    for(j in 1:K) {
      res <- cbind(h_Z, sorted_time[j]) / K
      integrand <- (T0 > sorted_time[j]) / exp(-noncensor_cumhaz[j]) - exp(-cbind(h_W, sorted_time[j]) %*% b * sorted_time[j])
      g <- g + c(integrand) * res
    }
    g <- colMeans(g)
    return(g)
  }
  design_mat <- jacobian(func = function(...) {return(hbridge_tmp(...))}, 
                         x = b, h_W = h_W, h_Z = h_Z, T0 = T0, sorted_time = sorted_time, noncensor_cumhaz = noncensor_cumhaz)
  censor_mat <- h_censor(b, h_Z, T0, sorted_time, noncensor_cumhaz)
  K <- length(sorted_time)
  N <- nrow(h_W)
  eps <- 0
  
  for(j in 1:K) {
    res <- cbind(h_Z, sorted_time[j]) / K
    integrand <- (T0 > sorted_time[j]) / exp(-noncensor_cumhaz[j]) - exp(-cbind(h_W, sorted_time[j]) %*% b * sorted_time[j])
    eps <- eps + c(integrand) * res
  }
  eps <- eps / N
  return(-solve(design_mat) %*% (t(eps) + censor_mat %*% noncensor_cumhaz_IF))
}

h_censor <- function(b, h_Z, T0, sorted_time, noncensor_cumhaz) { 
  K <- length(sorted_time)
  G <- matrix(0, nrow = length(b), ncol = K)
  for(j in 1:K) {
    res <- cbind(h_Z, sorted_time[j]) / K
    integrand <- (T0 > sorted_time[j]) / exp(-noncensor_cumhaz[j])
    G[, j] <- colMeans(c(integrand) * res)
  }
  return(G)
}



POR_surv <- function(sorted_time, b, h_W) {
  K <- length(sorted_time)
  N <- nrow(h_W)
  ## q value, N * 1
  IF <- matrix(0, nrow = K, ncol = N)
  
  for(j in 1:K) {
    dh <- exp(-cbind(h_W, sorted_time[j]) %*% b * sorted_time[j])
    IF[j, ] <- c(dh)/N
  }
  return(IF)
}



POR_h_dot <- function(sorted_time, b, h_W) {
  G <- jacobian(func = function(...) {return(rowSums(POR_surv(...)))}, 
                x = b, sorted_time = sorted_time, h_W = h_W)
  return(G)
}


MSE_func_q <- function(bridge_func, para, q_Y, q_W, q_Z){
  g0 <- bridge_func(para = para, q_Y = q_Y, q_W = q_W, q_Z = q_Z)
  g <- apply(g0, 2, mean)
  gmmf <- sum(g^2)
  return(gmmf)
}

qbridge <- function(para, q_Y, q_W, q_Z) { 
  tlink <- 1 + exp(q_Z %*% para)
  g0 <- q_W
  g <- c(tlink) * g0 - q_Y
  return(g)
}


q_IF <- function(t, q_Y, q_W, q_Z) {
  return(-solve(t(q_W) %*% (c(exp(q_Z %*% t)) * q_Z)) %*% t(c(1 + exp(q_Z %*% t)) * q_W - q_Y))
}




PIPW_surv <- function(sorted_time, time, noncensor_cumhaz, A, q_Z, t, a = 1) {
  K <- length(sorted_time)
  N <- length(time)
  ## q value, N * 1
  q <- c(1 + exp(q_Z %*% t))
  IF <- matrix(0, nrow = K, ncol = N)
  for(j in 1:K) {
    dN <- (A == a) * q * (time > sorted_time[j]) / exp(-noncensor_cumhaz[j])
    IF[j, ] <- dN / N
  }
  return(IF)
}

PIPW_q_dot <- function(sorted_time, time, noncensor_cumhaz, A, q_Z, t, a) {
  G <- jacobian(func = function(...) {return(rowSums(PIPW_surv(...)))}, 
                x = t, sorted_time = sorted_time, time = time, 
                noncensor_cumhaz = noncensor_cumhaz, A = A, q_Z = q_Z, a = a)
  return(G)
}

PIPW_censor_dot <- function(sorted_time, time, noncensor_cumhaz, A, q_Z, t, a) {
  K <- length(sorted_time)
  N <- length(time)
  ## q value, N * 1
  q <- c(1 + exp(q_Z %*% t))
  G <- matrix(0, nrow = K, ncol = K)
  for(j in 1:K) {
    dN <- (A == a) * q * (time > sorted_time[j]) / exp(-noncensor_cumhaz[j])
    G[j, j] <- mean(dN)
  }
  return(G)
}


PDR_surv <- function(b, t, time, sorted_time, noncensor_cumhaz, h_W, q_Z, A, a = 1) {
  K <- length(sorted_time)
  N <- nrow(h_W)
  IF <- matrix(0, nrow = K, ncol = N)
  q <- c(1 + exp(q_Z %*% t))
  
  for(j in 1:K) {
    dh <- exp(-cbind(h_W, sorted_time[j]) %*% b * sorted_time[j])
    dN <- (A == a) * q * ((time > sorted_time[j]) / exp(-noncensor_cumhaz[j]) - dh)
    IF[j, ] <- c(dN + dh)/N
  }
  return(IF)
}


PDR_h_dot <- function(b, t, time, sorted_time, noncensor_cumhaz, h_W, q_Z, A, a) {
  G <- jacobian(func = function(...) {return(rowSums(PDR_surv(...)))}, 
                x = b, sorted_time = sorted_time, time = time, 
                noncensor_cumhaz = noncensor_cumhaz, A = A, h_W = h_W, t = t, q_Z = q_Z, a = a)
  return(G)
}

PDR_q_dot <- function(b, t, time, sorted_time, noncensor_cumhaz, h_W, q_Z, A, a) {
  G <- jacobian(func = function(...) {return(rowSums(PDR_surv(...)))}, 
                x = t, sorted_time = sorted_time, time = time, 
                noncensor_cumhaz = noncensor_cumhaz, A = A, h_W = h_W, b = b, q_Z = q_Z, a = a)
  return(G)
}

PDR_censor_dot <- function(b, t, time, sorted_time, noncensor_cumhaz, h_W, q_Z, A, a = 1) {
  K <- length(sorted_time)
  N <- nrow(h_W)
  IF <- matrix(0, nrow = N, ncol = K)
  q <- c(1 + exp(q_Z %*% t))
  G <- matrix(0, nrow = K, ncol = K)
  for(j in 1:K) {
    dh <- pmin(exp(-cbind(h_W, sorted_time[j]) %*% b * sorted_time[j]), 1)
    dN <- (A == a) * q * ((time > sorted_time[j]) / exp(-noncensor_cumhaz[j]) - dh)
    G[j, j] <- mean(dN + dh)
  }
  return(G)
}

trimmed_mean <- function(x, alpha = 0.1) {
  return(mean(x[x <= quantile(x, 1 - alpha/2) & x >= quantile(x, alpha/2)]))
}

trimmed_sd <- function(x) {
  x <- x[x <=  quantile(x, 0.975) & x >= quantile(x, 0.025)]
  return(sd(x))
}


supreme_test <- function(est, IF, n_sim) {
  N <- nrow(IF)
  null_simu <- c()
  for(rep in 1:n_sim) {
    Q_n <- rnorm(N, 0, 1)
    supreme <- max(abs(colSums(Q_n * apply(IF, 1, cumsum))))
    null_simu <- c(null_simu, supreme)
  }
  return(mean(null_simu >= max(abs(cumsum(est)))))
}



reg_DR <- function(T0, sorted_time, noncensor_cumhaz, or_pred, or_pred_A1, or_pred_A0, IPW_weights) {
  N <- length(T0)
  sorted_time <- c(0, sorted_time)
  
  sorted_time_mat <- matrix(rep(sorted_time, each = N), nrow = N)
  reg_DR_est <- IPW_weights * (t(t(T0 > sorted_time_mat) / exp(-noncensor_cumhaz)) - or_pred) + or_pred_A1 - or_pred_A0
  return(colMeans(reg_DR_est))
}


boot_reg_DR <- function(T0, sorted_time, noncensor_cumhaz, or_pred, or_pred_A1, or_pred_A0, IPW_weights) {
  N <- length(T0)
  #sorted_time <- c(0, sorted_time)
  
  sorted_time_mat <- matrix(rep(sorted_time, each = N), nrow = N)
  reg_DR_est <- IPW_weights * (t(t(T0 > sorted_time_mat) / exp(-noncensor_cumhaz)) - or_pred) + or_pred_A1 - or_pred_A0
  return(colMeans(reg_DR_est))
}





normalize <- function(x) {
  return((x - min(x))/mean((x - min(x))))
}




