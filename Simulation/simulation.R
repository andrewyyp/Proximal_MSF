rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
#set.seed(123232)
source("data_generating.R")
require(lava)
library(nleqslv)
library(MASS)
library(KRLS)
library(numDeriv)
library(gmm)
library(survival)
library(timereg)

N = 100000

T0_A0 <- data_gen(N, para_set, a = c(0))$T0
control <- ecdf(T0_A0)
T0_A1 <- data_gen(N, para_set, a = c(1))$T0
treated <- ecdf(T0_A1)
time_gap <- 0.25


truth <- c(treated(time_gap * 1) - control(time_gap * 1), 
            treated(time_gap * 2) - control(time_gap * 2), 
            treated(time_gap * 3) - control(time_gap * 3))

truth

true_b(para_set)
true_t(para_set)
    
source("proxicausal.R")

N = 1000

rep_num = 1000

POR_est_result <- c()
PIPW_est_result <- c()
POR_WOR_est_result <- c()
PIPW_WIPW_est_result <- c()
PDR_est_result <- c()
PDR_WOR_est_result <- c()
PDR_WIPW_est_result <- c()
PDR_BW_est_result <- c()
DR_est_result <- c()

POR_sd_result <- c()
PIPW_sd_result <- c()
POR_WOR_sd_result <- c()
PIPW_WIPW_sd_result <- c()
PDR_sd_result <- c()
PDR_WOR_sd_result <- c()
PDR_WIPW_sd_result <- c()
PDR_BW_sd_result <- c()
DR_sd_result <- c()

POR_covering <- c()
PIPW_covering <- c()
POR_WOR_covering <- c()
PIPW_WIPW_covering <- c()
PDR_covering <- c()
PDR_WOR_covering <- c()
PDR_WIPW_covering <- c()
PDR_BW_covering <- c()
DR_covering <- c()
b_est <- c()
t_est <- c()
for(rep in 1:rep_num) {
  df <- data_gen(N, para_set)
  sorted_time <- with(df, sort(unique(T0[Delta == 1])))
  K = length(sorted_time)
  censor_sorted_time <- with(df, sort(unique(T0[Delta == 0])))
  censor_cumhaz_sub <- with(df, censor_cumhaz_est(censor_sorted_time, T0, Delta))
  censor_cumhaz <- colSums(censor_cumhaz_sub)
  censor_cumhaz_IF <- t(t(censor_cumhaz_sub) - censor_cumhaz)
  noncensor_cumhaz <- noncensor_cumhaz_compute(sorted_time, censor_sorted_time, censor_cumhaz)
  noncensor_cumhaz_IF <- noncensor_cumhaz_IF_compute(sorted_time, censor_sorted_time, censor_cumhaz_IF)
  

  source("simu_POR.R")
  POR_est_result <- rbind(POR_est_result, POR_est)
  POR_sd_result <- rbind(POR_sd_result, POR_sd_est)
  POR_covering <- rbind(POR_covering, truth <= POR_est + 1.96 * POR_sd_est & truth >= POR_est - 1.96 * POR_sd_est)
  b_est <- rbind(b_est, b)

  source("simu_PIPW.R")
  PIPW_est_result <- rbind(PIPW_est_result, PIPW_est)
  PIPW_sd_result <- rbind(PIPW_sd_result, PIPW_sd_est)
  PIPW_covering <- rbind(PIPW_covering, truth <= PIPW_est + 1.96 * PIPW_sd_est & truth >= PIPW_est - 1.96 * PIPW_sd_est)
  t_est <- rbind(t_est, t)

  source("simu_PDR.R")
  PDR_est_result <- rbind(PDR_est_result, PDR_est)
  PDR_sd_result <- rbind(PDR_sd_result, PDR_sd_est)
  PDR_covering <- rbind(PDR_covering, truth <= PDR_est + 1.96 * PDR_sd_est &
                          truth >= PDR_est - 1.96 * PDR_sd_est)


  source("simu_POR_WOR.R")
  POR_WOR_est_result <- rbind(POR_WOR_est_result, POR_WOR_est)
  POR_WOR_sd_result <- rbind(POR_WOR_sd_result, POR_WOR_sd_est)
  POR_WOR_covering <- rbind(POR_WOR_covering, truth <= POR_WOR_est + 1.96 * POR_WOR_sd_est &
                              truth >= POR_WOR_est - 1.96 * POR_WOR_sd_est)


  source("simu_PDR_WOR.R")
  PDR_WOR_est_result <- rbind(PDR_WOR_est_result, PDR_WOR_est)
  PDR_WOR_sd_result <- rbind(PDR_WOR_sd_result, PDR_WOR_sd_est)
  PDR_WOR_covering <- rbind(PDR_WOR_covering, truth <= PDR_WOR_est + 1.96 * PDR_WOR_sd_est &
                              truth >= PDR_WOR_est - 1.96 * PDR_WOR_sd_est)


  source("simu_PIPW_WIPW.R")
  PIPW_WIPW_est_result <- rbind(PIPW_WIPW_est_result, PIPW_WIPW_est)
  PIPW_WIPW_sd_result <- rbind(PIPW_WIPW_sd_result, PIPW_WIPW_sd_est)
  PIPW_WIPW_covering <- rbind(PIPW_WIPW_covering, truth <= PIPW_WIPW_est + 1.96 * PIPW_WIPW_sd_est &
                                truth >= PIPW_WIPW_est - 1.96 * PIPW_WIPW_sd_est)


  source("simu_PDR_WIPW.R")
  PDR_WIPW_est_result <- rbind(PDR_WIPW_est_result, PDR_WIPW_est)
  PDR_WIPW_sd_result <- rbind(PDR_WIPW_sd_result, PDR_WIPW_sd_est)
  PDR_WIPW_covering <- rbind(PDR_WIPW_covering, truth <= PDR_WIPW_est + 1.96 * PDR_WIPW_sd_est &
                               truth >= PDR_WIPW_est - 1.96 * PDR_WIPW_sd_est)


  source("simu_PDR_BW.R")
  PDR_BW_est_result <- rbind(PDR_BW_est_result, PDR_BW_est)
  PDR_BW_sd_result <- rbind(PDR_BW_sd_result, PDR_BW_sd_est)
  PDR_BW_covering <- rbind(PDR_BW_covering, truth <= PDR_BW_est + 1.96 * PDR_BW_sd_est &
                               truth >= PDR_BW_est - 1.96 * PDR_BW_sd_est)


  source("simu_DR.R")
  DR_est_result <- rbind(DR_est_result, DR_est)
  DR_sd_result <- rbind(DR_sd_result, DR_sd_est)
  DR_covering <- rbind(DR_covering, truth <= DR_est + 1.96 * DR_sd_est &
                              truth >= DR_est - 1.96 * DR_sd_est)

  print(rep)
}

apply(POR_est_result, 2, mean) - truth
apply(POR_est_result, 2, sd)
apply(POR_sd_result, 2, mean)
apply(POR_covering, 2, mean)
true_b(para_set)
colMeans(b_est)
apply(PIPW_est_result, 2, mean) - truth
apply(PIPW_est_result, 2, sd)
apply(PIPW_sd_result, 2, mean)
apply(PIPW_covering, 2, mean)
true_t(para_set)
colMeans(t_est)

