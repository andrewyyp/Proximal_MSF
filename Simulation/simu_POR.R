POR_df <- df
## first step: fitting h confounding bridge function
inioptim_b <- c(0.1, 0.1, 0.1, 0.1, 0.1)
h_W <- with(POR_df, cbind(1, A, W, X))
h_Z <- with(POR_df, cbind(1, A, Z, X))


hlink <- optim(par = inioptim_b, fn = hbridge, 
        h_W = h_W, h_Z = h_Z, T0 = POR_df$T0, sorted_time = sorted_time, 
        noncensor_cumhaz = noncensor_cumhaz, method = "BFGS", hessian = FALSE)

b <- hlink$par
## extract the IF of q function, dim N * 5
b_IF <- h_IF(b, h_W, h_Z, POR_df$T0, sorted_time, noncensor_cumhaz, noncensor_cumhaz_IF)


h_W_1 <- with(POR_df, cbind(1, 1, W, X))
treated_surv_sub <- POR_surv(sorted_time, b, h_W_1)
treated_surv <- rowSums(treated_surv_sub)
treated_surv_IF <- treated_surv_sub - treated_surv / N
treated_surv_func <- stepfun(sorted_time, c(1, treated_surv))
h_dot_treated <- POR_h_dot(sorted_time, b, h_W_1)

h_W_0 <- with(POR_df, cbind(1, 0, W, X))
control_surv_sub <- POR_surv(sorted_time, b, h_W_0)
control_surv <- rowSums(control_surv_sub)
control_surv_IF <- control_surv_sub - control_surv / N
control_surv_func <- stepfun(sorted_time, c(1, control_surv))
h_dot_control <- POR_h_dot(sorted_time, b, h_W_0)


IF <- matrix(0, nrow = K, ncol = N)
for (i in 1:N) {
  for (j in 1:K) {
    IF[j, i] <- treated_surv_IF[j, i] - control_surv_IF[j, i] + 
      sum(b_IF[, i] * (h_dot_treated[j, ] - h_dot_control[j, ]))
  }
}

POR_sd_func <- stepfun(sorted_time, c(0, sqrt(rowSums(IF^2))))

POR_est <- c(treated_surv_func(time_gap * 1) - control_surv_func(time_gap * 1), 
             treated_surv_func(time_gap * 2) - control_surv_func(time_gap * 2), 
             treated_surv_func(time_gap * 3) - control_surv_func(time_gap * 3))
POR_sd_est <- c(POR_sd_func(time_gap * 1), POR_sd_func(time_gap * 2), POR_sd_func(time_gap * 3))


