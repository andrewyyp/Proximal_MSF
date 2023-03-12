PDR_BW_df <- df
PDR_BW_df$W <- abs(PDR_BW_df$W)^(1/2) 
PDR_BW_df$Z <- abs(PDR_BW_df$Z)^(1/2) 


b <- wrong_b
b_IF <- wrong_b_IF

t <- wrong_t
t_IF <- wrong_t_IF

q_Y <- with(PDR_BW_df, matrix(rep(c(0, 1, 0, 0), each = length(A)), nrow = length(A), ncol = 4))
q_Z <- with(PDR_BW_df, (-1)^(1 - A) * cbind(1, A, Z, X))
q_W <- with(PDR_BW_df, (-1)^(1 - A) * cbind(1, A, W, X))

h_W_1 <- with(PDR_BW_df, cbind(1, 1, W, X))
treated_surv_sub <- PDR_surv(b, t, PDR_BW_df$T0, sorted_time, noncensor_cumhaz, h_W = h_W_1, q_Z = q_Z, PDR_BW_df$A, a = 1)
treated_surv <- rowSums(treated_surv_sub)
treated_surv_IF <- treated_surv_sub - treated_surv / N
treated_surv_func <- stepfun(sorted_time, c(1, treated_surv))
h_dot_treated <- PDR_h_dot(b, t, PDR_BW_df$T0, sorted_time, noncensor_cumhaz, h_W = h_W_1, q_Z = q_Z, PDR_BW_df$A, a = 1)
q_dot_treated <- PDR_q_dot(b, t, PDR_BW_df$T0, sorted_time, noncensor_cumhaz, h_W = h_W_1, q_Z = q_Z, PDR_BW_df$A, a = 1)
q_censor_treated <- PDR_censor_dot(b, t, PDR_BW_df$T0, sorted_time, noncensor_cumhaz = noncensor_cumhaz, h_W_1, q_Z = q_Z, PDR_BW_df$A, a = 1)


h_W_0 <- with(PDR_BW_df, cbind(1, 0, W, X))
control_surv_sub <- PDR_surv(b, t, PDR_BW_df$T0, sorted_time, noncensor_cumhaz, h_W = h_W_0, q_Z = q_Z, PDR_BW_df$A, a = 0)
control_surv <- rowSums(control_surv_sub)
control_surv_IF <- control_surv_sub - control_surv / N
control_surv_func <- stepfun(sorted_time, c(1, control_surv))
h_dot_control <- PDR_h_dot(b, t, PDR_BW_df$T0, sorted_time, noncensor_cumhaz, h_W = h_W_0, q_Z = q_Z, PDR_BW_df$A, a = 0)
q_dot_control <- PDR_q_dot(b, t, PDR_BW_df$T0, sorted_time, noncensor_cumhaz, h_W = h_W_0, q_Z = q_Z, PDR_BW_df$A, a = 0)
q_censor_control <- PDR_censor_dot(b, t, PDR_BW_df$T0, sorted_time, noncensor_cumhaz = noncensor_cumhaz, h_W_0, q_Z = q_Z, PDR_BW_df$A, a = 0)



IF <- matrix(0, nrow = K, ncol = N)
for (i in 1:N) {
  for (j in 1:K) {
    IF[j, i] <- treated_surv_IF[j, i] - control_surv_IF[j, i] + 
      sum(b_IF[, i] * (h_dot_treated[j, ] - h_dot_control[j, ])) + 
      sum(t_IF[, i] * (q_dot_treated[j, ] - q_dot_control[j, ])) +
      sum(noncensor_cumhaz_IF[, i] * (q_censor_treated[j, ] - q_censor_control[j, ]))
  }
}

PDR_BW_sd_func <- stepfun(sorted_time, c(0, sqrt(rowSums(IF^2))))

## point estimates and SD at each time point
PDR_BW_est <- c(treated_surv_func(time_gap * 1) - control_surv_func(time_gap * 1), 
                  treated_surv_func(time_gap * 2) - control_surv_func(time_gap * 2), 
                  treated_surv_func(time_gap * 3) - control_surv_func(time_gap * 3))
PDR_BW_sd_est <- c(PDR_BW_sd_func(time_gap * 1), PDR_BW_sd_func(time_gap * 2), PDR_BW_sd_func(time_gap * 3))

