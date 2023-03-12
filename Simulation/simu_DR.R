DR_df <- df
DR_df_A1 <- DR_df
DR_df_A1$A <- 1
DR_df_A0 <- DR_df
DR_df_A0$A <- 0

ipw_stage <- glm(A ~ X + W + Z, family = "binomial", data = DR_df)

p_of_A <- fitted(ipw_stage)

IPW_weights <- with(DR_df, A/p_of_A - (1- A)/(1 - p_of_A))
IPW_weights_A1 <- with(DR_df, A/p_of_A)

IPW_weights_A0 <- with(DR_df, (1- A)/(1 - p_of_A))

or_stage <- aalen(Surv(time = T0, event = Delta) ~ const(A) + const(X) + const(W) + const(Z), data = DR_df, n.sim = 0)

or_pred <- predict(or_stage, newdata = DR_df, uniform = FALSE, se = FALSE, n.sim = 0)$S0

or_pred_A1 <- predict(or_stage, newdata = DR_df_A1, uniform = FALSE, se = FALSE, n.sim = 0)$S0

or_pred_A0 <- predict(or_stage, newdata = DR_df_A0, uniform = FALSE, se = FALSE, n.sim = 0)$S0

sorted_time <- with(df, sort(unique(T0[Delta == 1])))
K = length(sorted_time)
noncensor_cumhaz <- noncensor_cumhaz_compute(sorted_time, censor_sorted_time, censor_cumhaz)
noncensor_cumhaz <- c(0, noncensor_cumhaz)

DR_surv <- reg_DR(DR_df$T0, sorted_time, noncensor_cumhaz, or_pred, or_pred_A1, or_pred_A0, IPW_weights)
DR_surv <- stepfun(sorted_time, DR_surv)

DR_est <- c(DR_surv(time_gap * 1), DR_surv(time_gap * 2), DR_surv(time_gap * 3))



# bootstrap
boot_total <- 1000
boot_est_result <- c()
for(boot_rep in 1:boot_total) {
  index <- sample(1:N, N, replace = TRUE)
  tmp_df <- df[index, ]
  DR_df <- tmp_df
  DR_df_A1 <- DR_df
  DR_df_A1$A <- 1
  DR_df_A0 <- DR_df
  DR_df_A0$A <- 0
  
  ipw_stage <- glm(A ~ X + W + Z, family = "binomial", data = DR_df)
  
  p_of_A <- fitted(ipw_stage)
  
  IPW_weights <- with(DR_df, A/p_of_A - (1- A)/(1 - p_of_A))
  
  or_stage <- aalen(Surv(time = T0, event = Delta) ~ const(A) + const(X) + const(W) + const(Z), data = DR_df, n.sim = 0)
  
  or_pred <- predict(or_stage, newdata = DR_df, uniform = FALSE, se = FALSE, n.sim = 0)
  
  sorted_time <- or_pred$time
  or_pred <- or_pred$S0
  
  
  or_pred_A1 <- predict(or_stage, newdata = DR_df_A1, uniform = FALSE, se = FALSE, n.sim = 0)$S0
  
  or_pred_A0 <- predict(or_stage, newdata = DR_df_A0, uniform = FALSE, se = FALSE, n.sim = 0)$S0
  

  K = length(sorted_time)
  censor_sorted_time <- with(DR_df, sort(unique(T0[Delta == 0])))
  censor_cumhaz_sub <- with(DR_df, censor_cumhaz_est(censor_sorted_time, T0, Delta))
  censor_cumhaz <- rowSums(censor_cumhaz_sub)
  noncensor_cumhaz <- noncensor_cumhaz_compute(sorted_time, censor_sorted_time, censor_cumhaz)

  DR_surv <- boot_reg_DR(DR_df$T0, sorted_time, noncensor_cumhaz, or_pred, or_pred_A1, or_pred_A0, IPW_weights)
  DR_surv <- stepfun(sorted_time, c(0, DR_surv))
  boot_est_result <- rbind(boot_est_result, c(DR_surv(time_gap * 1), DR_surv(time_gap * 2), DR_surv(time_gap * 3)))
}

DR_sd_est <- apply(boot_est_result, 2, sd)
