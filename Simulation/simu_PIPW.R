PIPW_df <- df
## first step: fitting q confounding bridge function
## dimension of parameter t in q, e.g., 4 in this case
inioptim_t <- c(0, 0, 0, 0)
q_Y <- with(PIPW_df, matrix(rep(c(0, 1, 0, 0), each = length(A)), nrow = length(A), ncol = 4))
q_Z <- with(PIPW_df, (-1)^(1 - A) * cbind(1, A, Z, X))
q_W <- with(PIPW_df, (-1)^(1 - A) * cbind(1, A, W, X))
qlink <- optim(par = inioptim_t, fn = MSE_func_q,
                bridge_func = qbridge, q_Y = q_Y, q_W = q_W, q_Z = q_Z,
                method = "BFGS", hessian = FALSE)
## dim 4
t <- qlink$par
## extract the IF of q function, dim N * 4
t_IF <- q_IF(t, q_Y = q_Y, q_W = q_W, q_Z = q_Z)


## second step: the PIPW estimator, we will do this by each arm
## treated arm
treated_surv_sub <- PIPW_surv(sorted_time, PIPW_df$T0, noncensor_cumhaz = noncensor_cumhaz, PIPW_df$A, q_Z = q_Z, t, a = 1)
treated_surv <- rowSums(treated_surv_sub)
treated_surv_IF <- treated_surv_sub - treated_surv / N
treated_surv_func <- stepfun(sorted_time, c(1, treated_surv))
# dim 4 * K
q_dot_treated <- PIPW_q_dot(sorted_time, PIPW_df$T0, noncensor_cumhaz = noncensor_cumhaz, PIPW_df$A, q_Z = q_Z, t, a = 1)
q_censor_treated <- PIPW_censor_dot(sorted_time, PIPW_df$T0, noncensor_cumhaz = noncensor_cumhaz, PIPW_df$A, q_Z = q_Z, t, a = 1)


## control arm
control_surv_sub <- PIPW_surv(sorted_time, PIPW_df$T0, noncensor_cumhaz = noncensor_cumhaz, PIPW_df$A, q_Z = q_Z, t, a = 0)
control_surv <- rowSums(control_surv_sub)
control_surv_IF <- control_surv_sub - control_surv / N
control_surv_func <- stepfun(sorted_time, c(1, control_surv))
# dim 4 * K
q_dot_control <- PIPW_q_dot(sorted_time, PIPW_df$T0, noncensor_cumhaz = noncensor_cumhaz, PIPW_df$A, q_Z = q_Z, t, a = 0)
q_censor_control <- PIPW_censor_dot(sorted_time, PIPW_df$T0, noncensor_cumhaz = noncensor_cumhaz, PIPW_df$A, q_Z = q_Z, t, a = 0)


IF <- matrix(0, nrow = K, ncol = N)
for (i in 1:N) {
  for (j in 1:K) {
    IF[j, i] <- treated_surv_IF[j, i] - control_surv_IF[j, i] + 
                sum(t_IF[, i] * (q_dot_treated[j, ] - q_dot_control[j, ])) + 
                sum(noncensor_cumhaz_IF[, i] * (q_censor_treated[j, ] - q_censor_control[j, ]))
  }
}

PIPW_sd_func <- stepfun(sorted_time, c(0, sqrt(rowSums(IF^2))))

## point estimates and SD at each time point
PIPW_est <- c(treated_surv_func(time_gap * 1) - control_surv_func(time_gap * 1), 
              treated_surv_func(time_gap * 2) - control_surv_func(time_gap * 2), 
              treated_surv_func(time_gap * 3) - control_surv_func(time_gap * 3))
PIPW_sd_est <- c(PIPW_sd_func(time_gap * 1), PIPW_sd_func(time_gap * 2), PIPW_sd_func(time_gap * 3))

