library(aghq)
library(TMB)
library(Matrix)
library(xml2)
library(tidyverse)
library(mvQuad)
library(fda)
library(LaplacesDemon)
library(mvtnorm)
library(OSplines)
library(sGPfit)

source(paste0(source_path, "06_functions_seasonal.R"))
source(paste0(source_path, "00_functions_predictGP.R"))
source(paste0(source_path, "00_functions_simHyper.R"))


refined_x <- seq(min(x), max(x), by = 0.01)
refined_years <- refined_x + 1960
pred_upper <- as.numeric(as.Date("2030-01-01") - timeOrigin)/365.25
pred_x <- seq(max(refined_x), pred_upper, by = 0.01)
pred_years <- pred_x + 1960


period1 <- 1 ## 1-year
a1 = 2*period1*pi
B1 <- Compute_B_sB(x = refined_x, a = a1, k = k2, region = region, boundary = T)
B1 <- as(B1,"dgTMatrix")
Q1 <- Compute_Q_sB(a = a1, k = k2, region = region)
Q1 <- as(as(Q1, "matrix"),"dgTMatrix")
X1 <- as(cbind(cos(a1*refined_x), sin(a1*refined_x)), "dgTMatrix")
period2 <- 2 ## half-year
a2 = 2*period2*pi
B2 <- Compute_B_sB(x = refined_x, a = a2, k = k2, region = region, boundary = T)
B2 <- as(B2,"dgTMatrix")
Q2 <- Compute_Q_sB(a = a2, k = k2, region = region)
Q2 <- as(as(Q2, "matrix"),"dgTMatrix")
X2 <- as(cbind(cos(a2*refined_x), sin(a2*refined_x)), "dgTMatrix")

Q6 <- as(compute_weights_precision(x = seq(0, max(x), length.out = k1)), "dgTMatrix")
B6 <- as(local_poly(knots = seq(0, max(x), length.out = k1), refined_x = refined_x, p = 3), "dgTMatrix")


set.seed(123)
nsamps <- 3000

### Load the corresponding samples:
## Obtain g1,g2,g3,g4,g5
samps_g1 <- data.frame(x = refined_x)
samps_g2 <- data.frame(x = refined_x)
samps_g3 <- data.frame(x = refined_x)
samps_g4 <- data.frame(x = refined_x)
samps_g5 <- data.frame(x = refined_x)

samps_g1_1st <- data.frame(x = refined_x)
samps_g2_1st <- data.frame(x = refined_x)
samps_g3_1st <- data.frame(x = refined_x)
samps_g4_1st <- data.frame(x = refined_x)
samps_g5_1st <- data.frame(x = refined_x)

samps_g6 <- data.frame(x = refined_x)
samps_IWP_1st <- data.frame(x = refined_x)
samps_IWP_2nd <- data.frame(x = refined_x)

pred_g1 <- data.frame(x = pred_x)
pred_g2 <- data.frame(x = pred_x)
pred_g3 <- data.frame(x = pred_x)
pred_g4 <- data.frame(x = pred_x)
pred_g5 <- data.frame(x = pred_x)
pred_g6 <- data.frame(x = pred_x)

pred_IWP_1st <- data.frame(x = pred_x)


sd1 <- c()
sd2 <- c()
sd3 <- c()
sd4 <- c()
sd5 <- c()
sd6 <- c()
sd7 <- c()

years_cyl3 <- 44/12
years_cyl4 <- 9.1
years_cyl5 <- 10.4

period3 <- 1/years_cyl3
period4 <- 1/years_cyl4
period5 <- 1/years_cyl5

a3 = 2*period3*pi
B3 <- Compute_B_sB(x = refined_x, a = a3, k = k2, region = region, boundary = T)
B3 <- as(B3,"dgTMatrix")
X3 <- as(cbind(cos(a3*refined_x), sin(a3*refined_x)), "dgTMatrix")
a4 = 2*period4*pi
B4 <- Compute_B_sB(x = refined_x, a = a4, k = k2, region = region, boundary = T)
B4 <- as(B4,"dgTMatrix")
X4 <- as(cbind(cos(a4*refined_x), sin(a4*refined_x)), "dgTMatrix")
a5 = 2*period5*pi
B5 <- Compute_B_sB(x = refined_x, a = a5, k = k2, region = region, boundary = T)
B5 <- as(B5,"dgTMatrix")
X5 <- as(cbind(cos(a5*refined_x), sin(a5*refined_x)), "dgTMatrix")

X <- as(cbind(X1,X2,X3,X4,X5, 1, refined_x, refined_x^2), "dgTMatrix") 
B <- cbind(B1,B2,B3,B4,B5,B6)

load(file = paste0(result_path, "samps/", "1", "_sample.rda"))
sGP_samps$samps <- sGP_samps$samps[,(1:nsamps), drop = F]
samps_w <- sGP_samps$samps[1:ncol(B), , drop = F]
samps_beta <- sGP_samps$samps[(ncol(B) + 1):(ncol(B) + ncol(X)),, drop = F]

samps_g1_new <- as.matrix(B1 %*% samps_w[1:ncol(B1),, drop = F]) + as.matrix(X[,1:2] %*% samps_beta[1:2,, drop = F])
samps_g1_new <- as.data.frame(samps_g1_new)
samps_g1 <- cbind(samps_g1,samps_g1_new)

samps_g2_new <- as.matrix(B2 %*% samps_w[(ncol(B1) + 1):(ncol(B1) + ncol(B2)),, drop = F]) + as.matrix(X[,3:4] %*% samps_beta[3:4,, drop = F])
samps_g2_new <- as.data.frame(samps_g2_new)
samps_g2 <- cbind(samps_g2,samps_g2_new)

samps_g3_new <- as.matrix(B3 %*% samps_w[(ncol(B1) + ncol(B2) + 1):(ncol(B1) + ncol(B2) + ncol(B3)),, drop = F]) + as.matrix(X[,5:6] %*% samps_beta[5:6,, drop = F])
samps_g3_new <- as.data.frame(samps_g3_new)
samps_g3 <- cbind(samps_g3,samps_g3_new)

samps_g4_new <- as.matrix(B4 %*% samps_w[(ncol(B1) + ncol(B2) + ncol(B3) + 1):(ncol(B1) + ncol(B2) + ncol(B3) + ncol(B4)),, drop = F]) + as.matrix(X[,7:8] %*% samps_beta[7:8,, drop = F])
samps_g4_new <- as.data.frame(samps_g4_new)
samps_g4 <- cbind(samps_g4,samps_g4_new)

samps_g5_new <- as.matrix(B5 %*% samps_w[(ncol(B1) + ncol(B2) + ncol(B3) + ncol(B4) + 1):(ncol(B1) + ncol(B2) + ncol(B3) + ncol(B4) + ncol(B5)),, drop = F]) + as.matrix(X[,9:10] %*% samps_beta[9:10,, drop = F])
samps_g5_new <- as.data.frame(samps_g5_new)
samps_g5 <- cbind(samps_g5,samps_g5_new)

samps_g6_new <- as.matrix(B6 %*% samps_w[(ncol(B1) + ncol(B2) + ncol(B3) + ncol(B4) + ncol(B5) + 1):(ncol(B1) + ncol(B2) + ncol(B3) + ncol(B4) + ncol(B5) + ncol(B6)),, drop = F]) + as.matrix(X[,11:13] %*% samps_beta[11:13,, drop = F])
samps_g6_new <- as.data.frame(samps_g6_new)
samps_g6 <- cbind(samps_g6,samps_g6_new)

B6_new <- as(local_poly(knots = seq(0, max(x), length.out = k1), refined_x = refined_x, p = 2), "dgTMatrix")
samps_IWP_1st_new <- as.matrix(B6_new %*% samps_w[(ncol(B1) + ncol(B2) + ncol(B3) + ncol(B4) + ncol(B5) + 1):(ncol(B1) + ncol(B2) + ncol(B3) + ncol(B4) + ncol(B5) + ncol(B6)),]) + as.matrix(cbind(rep(1,length(refined_x)),2*refined_x) %*% samps_beta[12:13,, drop = F])
samps_IWP_1st <- cbind(samps_IWP_1st, samps_IWP_1st_new)

B6_new <- as(local_poly(knots = seq(0, max(x), length.out = k1), refined_x = refined_x, p = 1), "dgTMatrix")
samps_IWP_2nd_new <- as.matrix(B6_new %*% samps_w[(ncol(B1) + ncol(B2) + ncol(B3) + ncol(B4) + ncol(B5) + 1):(ncol(B1) + ncol(B2) + ncol(B3) + ncol(B4) + ncol(B5) + ncol(B6)),]) + as.matrix(cbind(rep(2,length(refined_x))) %*% samps_beta[13,, drop = F])
samps_IWP_2nd <- cbind(samps_IWP_2nd, samps_IWP_2nd_new)

sd1 <- c(sd1, exp(-0.5*sGP_samps$theta[(1:nsamps),1]))
sd2 <- c(sd2, exp(-0.5*sGP_samps$theta[(1:nsamps),2]))
sd3 <- c(sd3, exp(-0.5*sGP_samps$theta[(1:nsamps),3]))
sd4 <- c(sd4, exp(-0.5*sGP_samps$theta[(1:nsamps),4]))
sd5 <- c(sd5, exp(-0.5*sGP_samps$theta[(1:nsamps),5]))
sd6 <- c(sd6, exp(-0.5*sGP_samps$theta[(1:nsamps),6]))
sd7 <- c(sd7, exp(-0.5*sGP_samps$theta[(1:nsamps),7]))

samps_g1_1st <- extract_deriv_samples(t = samps_g1$x, fitted_samps = samps_g1[,-1, drop = F], degree = 1)
samps_g2_1st <- extract_deriv_samples(t = samps_g2$x, fitted_samps = samps_g2[,-1, drop = F], degree = 1)
samps_g3_1st <- extract_deriv_samples(t = samps_g3$x, fitted_samps = samps_g3[,-1, drop = F], degree = 1)
samps_g4_1st <- extract_deriv_samples(t = samps_g4$x, fitted_samps = samps_g4[,-1, drop = F], degree = 1)
samps_g5_1st <- extract_deriv_samples(t = samps_g5$x, fitted_samps = samps_g5[,-1, drop = F], degree = 1)

### Hyper-parameter summary:
scale1 <- compute_d_step_sGPsd(d = d_step,a = a1)
c1 <- 1/scale1
scale2 <- compute_d_step_sGPsd(d = d_step,a = a2)
c2 <- 1/scale2
scale3 <- compute_d_step_sGPsd(d = d_step,a = a3)
c3 <- 1/scale3
scale4 <- compute_d_step_sGPsd(d = d_step,a = a4)
c4 <- 1/scale4
scale5 <- compute_d_step_sGPsd(d = d_step,a = a5)
c5 <- 1/scale5

prior_set_IWP_used <- prior_conversion(d = d_step, prior = prior_set_IWP, p = 3)
c6 <- prior_set_IWP$u/prior_set_IWP_used$u
# summary(c1*sd1)
# summary(c2*sd2)
# summary(c3*sd3)
# summary(c4*sd4)
# summary(c5*sd5)
# summary(c6*sd6)
# summary(sd7)
for (i in 1:nsamps) {
  a3 = 2*period3*pi
  a4 = 2*period4*pi
  a5 = 2*period5*pi
  
  pred_g1_new <- predict_sGP_Var(observed_f_vec = c(samps_g1[nrow(samps_g1),(1+i)], samps_g1_1st[nrow(samps_g1_1st),(1+i)]),
                                 last_t = min(pred_g1$x), max_t = max(pred_g1$x), a = a1, sd = sd1[i])
  pred_g1 <- cbind(pred_g1, pred_g1_new[,2])
  
  
  pred_g2_new <- predict_sGP_Var(observed_f_vec = c(samps_g2[nrow(samps_g2),(1+i)], samps_g2_1st[nrow(samps_g2_1st),(1+i)]),
                                 last_t = min(pred_g2$x), max_t = max(pred_g2$x), a = a2, sd = sd2[i])
  pred_g2 <- cbind(pred_g2, pred_g2_new[,2])
  
  
  pred_g3_new <- predict_sGP_Var(observed_f_vec = c(samps_g3[nrow(samps_g3),(1+i)], samps_g3_1st[nrow(samps_g3_1st),(1+i)]),
                                 last_t = min(pred_g3$x), max_t = max(pred_g3$x), a = a3, sd = sd3[i])
  pred_g3 <- cbind(pred_g3, pred_g3_new[,2])
  
  
  pred_g4_new <- predict_sGP_Var(observed_f_vec = c(samps_g4[nrow(samps_g4),(1+i)], samps_g4_1st[nrow(samps_g4_1st),(1+i)]),
                                 last_t = min(pred_g4$x), max_t = max(pred_g4$x), a = a4, sd = sd4[i])
  pred_g4 <- cbind(pred_g4, pred_g4_new[,2])
  
  pred_g5_new <- predict_sGP_Var(observed_f_vec = c(samps_g5[nrow(samps_g5),(1+i)], samps_g4_1st[nrow(samps_g5_1st),(1+i)]),
                                 last_t = min(pred_g5$x), max_t = max(pred_g5$x), a = a5, sd = sd5[i])
  pred_g5 <- cbind(pred_g5, pred_g5_new[,2])
  
  pred_g6_all_new <- predict_IWP_Var(observed_f_vec = c(samps_g6[nrow(samps_g6),(1+i)], samps_IWP_1st[nrow(samps_IWP_1st),(1+i)], samps_IWP_2nd[nrow(samps_IWP_2nd),(1+i)]),
                                     last_t = min(pred_g6$x), max_t = max(pred_g6$x), p = 3, sd = sd6[i])
  pred_g6 <- cbind(pred_g6, pred_g6_all_new[,2])
  pred_IWP_1st <- cbind(pred_IWP_1st, pred_g6_all_new[,3])
  
}


### Plotting:
### Overall:
samps_overall<- samps_g1 + samps_g2 + samps_g3 + samps_g4 + samps_g5 + samps_g6
samps_overall$x <- refined_x
upper_int <- samps_overall[,-1] %>% apply(MARGIN = 1, FUN = quantile, p = 0.975)
lower_int <- samps_overall[,-1] %>% apply(MARGIN = 1, FUN = quantile, p = 0.025)
mean <- samps_overall[,-1] %>% apply(MARGIN = 1, FUN = mean)

pred_overall<- pred_g1 + pred_g2 + pred_g3 + pred_g4 + pred_g5 + pred_g6
pred_overall$x <- pred_x
pred_upper_int <- pred_overall[,-1] %>% apply(MARGIN = 1, FUN = quantile, p = 0.975)
pred_lower_int <- pred_overall[,-1] %>% apply(MARGIN = 1, FUN = quantile, p = 0.025)
pred_mean <- pred_overall[,-1] %>% apply(MARGIN = 1, FUN = mean)

pdf(file = paste0(figure_path, "sGP_overall.pdf"), height = 5, width = 5)
plot(mean~refined_years, type = 'l', col = 'blue', ylab = 'CO2', xlab = 'years', ylim = c(310,450), xlim = c(1960,2030))
lines(upper_int~refined_years, lty = 'dashed', col = 'red')
lines(lower_int~refined_years, lty = 'dashed', col = 'red')
lines(pred_mean~pred_years, lty = 'solid', col = "blue")
lines(pred_upper_int~pred_years, lty = 'dashed', col = 'red')
lines(pred_lower_int~pred_years, lty = 'dashed', col = 'red')
dev.off()

pdf(file = paste0(figure_path, "sGP_overall_pred.pdf"), height = 5, width = 5)
plot(mean~refined_years, type = 'l', col = 'blue', ylab = 'CO2', xlab = 'years', ylim = c(380,450), xlim = c(2010,2030))
points(observed_dataset$co2 ~ I(observed_dataset$timeYears + 1960), lwd = 1, cex = 0.1)
lines(upper_int~refined_years, lty = 'dashed', col = 'red')
lines(lower_int~refined_years, lty = 'dashed', col = 'red')
lines(pred_mean~pred_years, lty = 'solid', col = "blue")
lines(pred_upper_int~pred_years, lty = 'dashed', col = 'red')
lines(pred_lower_int~pred_years, lty = 'dashed', col = 'red')
dev.off()

### Trend
samps_trend<- samps_g6
samps_trend$x <- refined_x
upper_int <- samps_trend[,-1] %>% apply(MARGIN = 1, FUN = quantile, p = 0.975)
lower_int <- samps_trend[,-1] %>% apply(MARGIN = 1, FUN = quantile, p = 0.025)
mean <- samps_trend[,-1] %>% apply(MARGIN = 1, FUN = mean)
pred_trend<- pred_g6
pred_trend$x <- pred_x
pred_upper_int <- pred_trend[,-1] %>% apply(MARGIN = 1, FUN = quantile, p = 0.975)
pred_lower_int <- pred_trend[,-1] %>% apply(MARGIN = 1, FUN = quantile, p = 0.025)
pred_mean <- pred_trend[,-1] %>% apply(MARGIN = 1, FUN = mean)

pdf(file = paste0(figure_path, "sGP_trend.pdf"), height = 5, width = 5)
plot(mean~refined_years, type = 'l', col = 'blue', ylab = 'CO2', xlab = 'years', ylim = c(310,450), xlim = c(1960,2030))
lines(upper_int~refined_years, lty = 'dashed', col = 'red')
lines(lower_int~refined_years, lty = 'dashed', col = 'red')
lines(pred_mean~pred_years, lty = 'solid', col = "blue")
lines(pred_upper_int~pred_years, lty = 'dashed', col = 'red')
lines(pred_lower_int~pred_years, lty = 'dashed', col = 'red')
dev.off()

### Trend Derivative
samps_1st<- samps_IWP_1st
samps_1st$x <- refined_x
upper_int <- samps_1st[,-1] %>% apply(MARGIN = 1, FUN = quantile, p = 0.975)
lower_int <- samps_1st[,-1] %>% apply(MARGIN = 1, FUN = quantile, p = 0.025)
mean <- samps_1st[,-1] %>% apply(MARGIN = 1, FUN = mean)
pred_1st<- pred_IWP_1st
pred_1st$x <- pred_x
pred_upper_int <- pred_1st[,-1] %>% apply(MARGIN = 1, FUN = quantile, p = 0.975)
pred_lower_int <- pred_1st[,-1] %>% apply(MARGIN = 1, FUN = quantile, p = 0.025)
pred_mean <- pred_1st[,-1] %>% apply(MARGIN = 1, FUN = mean)
pdf(file = paste0(figure_path, "sGP_trend_1st.pdf"), height = 5, width = 5)
plot(mean~refined_years, type = 'l', col = 'blue', ylab = 'CO2', xlab = 'years', ylim = c(-2,5), xlim = c(1960,2030))
lines(upper_int~refined_years, lty = 'dashed', col = 'red')
lines(lower_int~refined_years, lty = 'dashed', col = 'red')
lines(pred_mean~pred_years, lty = 'solid', col = "blue")
lines(pred_upper_int~pred_years, lty = 'dashed', col = 'red')
lines(pred_lower_int~pred_years, lty = 'dashed', col = 'red')
dev.off()

### seasonal
samps_season<- samps_g1 + samps_g2 + samps_g3 + samps_g4 + samps_g5
samps_season$x <- refined_x
upper_int <- samps_season[,-1] %>% apply(MARGIN = 1, FUN = quantile, p = 0.975)
lower_int <- samps_season[,-1] %>% apply(MARGIN = 1, FUN = quantile, p = 0.025)
mean <- samps_season[,-1] %>% apply(MARGIN = 1, FUN = mean)
pred_season<- pred_g1 + pred_g2 + pred_g3 + pred_g4 + pred_g5
pred_season$x <- pred_x
pred_upper_int <- pred_season[,-1] %>% apply(MARGIN = 1, FUN = quantile, p = 0.975)
pred_lower_int <- pred_season[,-1] %>% apply(MARGIN = 1, FUN = quantile, p = 0.025)
pred_mean <- pred_season[,-1] %>% apply(MARGIN = 1, FUN = mean)

pdf(file = paste0(figure_path, "sGP_seasonal.pdf"), height = 5, width = 5)
plot(mean~refined_years, type = 'l', col = 'blue', ylab = 'CO2', xlab = 'years', ylim = c(-10,10), xlim = c(1960,2030))
lines(upper_int~refined_years, lty = 'dashed', col = 'red')
lines(lower_int~refined_years, lty = 'dashed', col = 'red')
lines(pred_mean~pred_years, lty = 'solid', col = "blue")
lines(pred_upper_int~pred_years, lty = 'dashed', col = 'red')
lines(pred_lower_int~pred_years, lty = 'dashed', col = 'red')
dev.off()

pdf(file = paste0(figure_path, "sGP_seasonal_pred.pdf"), height = 5, width = 5)
plot(mean~refined_years, type = 'l', col = 'blue', ylab = 'CO2', xlab = 'years', ylim = c(-10,10), xlim = c(2010,2030))
lines(upper_int~refined_years, lty = 'dashed', col = 'red')
lines(lower_int~refined_years, lty = 'dashed', col = 'red')
lines(pred_mean~pred_years, lty = 'solid', col = "blue")
lines(pred_upper_int~pred_years, lty = 'dashed', col = 'red')
lines(pred_lower_int~pred_years, lty = 'dashed', col = 'red')
dev.off()


### seasonal part 1:
samps_season<- samps_g1 + samps_g2
samps_season$x <- refined_x
upper_int <- samps_season[,-1] %>% apply(MARGIN = 1, FUN = quantile, p = 0.975)
lower_int <- samps_season[,-1] %>% apply(MARGIN = 1, FUN = quantile, p = 0.025)
mean <- samps_season[,-1] %>% apply(MARGIN = 1, FUN = mean)
pred_season<- pred_g1 + pred_g2
pred_season$x <- pred_x
pred_upper_int <- pred_season[,-1] %>% apply(MARGIN = 1, FUN = quantile, p = 0.975)
pred_lower_int <- pred_season[,-1] %>% apply(MARGIN = 1, FUN = quantile, p = 0.025)
pred_mean <- pred_season[,-1] %>% apply(MARGIN = 1, FUN = mean)

pdf(file = paste0(figure_path, "sGP_seasonal_annual.pdf"), height = 5, width = 5)
plot(mean~refined_years, type = 'l', col = 'blue', ylab = 'CO2', xlab = 'years', ylim = c(-10,10))
lines(upper_int~refined_years, lty = 'dashed', col = 'red')
lines(lower_int~refined_years, lty = 'dashed', col = 'red')
lines(pred_mean~pred_years, lty = 'solid', col = "blue")
lines(pred_upper_int~pred_years, lty = 'dashed', col = 'red')
lines(pred_lower_int~pred_years, lty = 'dashed', col = 'red')
abline(v = seq(min(refined_years),max(pred_years), by = (1/period1)), col = "green", lty = "dashed")
dev.off()

pdf(file = paste0(figure_path, "sGP_seasonal_annual_pred.pdf"), height = 5, width = 5)
plot(mean~refined_years, type = 'l', col = 'blue', ylab = 'CO2', xlab = 'years', ylim = c(-10,10), xlim = c(2010,2030))
lines(upper_int~refined_years, lty = 'dashed', col = 'red')
lines(lower_int~refined_years, lty = 'dashed', col = 'red')
lines(pred_mean~pred_years, lty = 'solid', col = "blue")
lines(pred_upper_int~pred_years, lty = 'dashed', col = 'red')
lines(pred_lower_int~pred_years, lty = 'dashed', col = 'red')
abline(v = seq(min(refined_years),max(pred_years), by = (1/period1)), col = "green", lty = "dashed")
dev.off()

### seasonal part 2:
samps_season<- samps_g3
samps_season$x <- refined_x
upper_int <- samps_season[,-1] %>% apply(MARGIN = 1, FUN = quantile, p = 0.975)
lower_int <- samps_season[,-1] %>% apply(MARGIN = 1, FUN = quantile, p = 0.025)
mean <- samps_season[,-1] %>% apply(MARGIN = 1, FUN = mean)
pred_season<-  pred_g3
pred_season$x <- pred_x
pred_upper_int <- pred_season[,-1] %>% apply(MARGIN = 1, FUN = quantile, p = 0.975)
pred_lower_int <- pred_season[,-1] %>% apply(MARGIN = 1, FUN = quantile, p = 0.025)
pred_mean <- pred_season[,-1] %>% apply(MARGIN = 1, FUN = mean)

pdf(file = paste0(figure_path, "sGP_seasonal_ENSO.pdf"), height = 5, width = 5)
plot(mean~refined_years, type = 'l', col = 'blue', ylab = 'CO2', xlab = 'years', ylim = c(-10,10))
lines(upper_int~refined_years, lty = 'dashed', col = 'red')
lines(lower_int~refined_years, lty = 'dashed', col = 'red')
lines(pred_mean~pred_years, lty = 'solid', col = "blue")
lines(pred_upper_int~pred_years, lty = 'dashed', col = 'red')
lines(pred_lower_int~pred_years, lty = 'dashed', col = 'red')
abline(v = seq(min(refined_years),max(pred_years), by = (1/period3)), col = "green", lty = "dashed")
dev.off()

pdf(file = paste0(figure_path, "sGP_seasonal_ENSO_pred.pdf"), height = 5, width = 5)
plot(mean~refined_years, type = 'l', col = 'blue', ylab = 'CO2', xlab = 'years', ylim = c(-10,10), xlim = c(2010,2030))
lines(upper_int~refined_years, lty = 'dashed', col = 'red')
lines(lower_int~refined_years, lty = 'dashed', col = 'red')
lines(pred_mean~pred_years, lty = 'solid', col = "blue")
lines(pred_upper_int~pred_years, lty = 'dashed', col = 'red')
lines(pred_lower_int~pred_years, lty = 'dashed', col = 'red')
abline(v = seq(min(refined_years),max(pred_years), by = (1/period3)), col = "green", lty = "dashed")
dev.off()


### seasonal part 3:
samps_season<- samps_g4
samps_season$x <- refined_x
upper_int <- samps_season[,-1] %>% apply(MARGIN = 1, FUN = quantile, p = 0.975)
lower_int <- samps_season[,-1] %>% apply(MARGIN = 1, FUN = quantile, p = 0.025)
mean <- samps_season[,-1] %>% apply(MARGIN = 1, FUN = mean)
pred_season<- pred_g4
pred_season$x <- pred_x
pred_upper_int <- pred_season[,-1] %>% apply(MARGIN = 1, FUN = quantile, p = 0.975)
pred_lower_int <- pred_season[,-1] %>% apply(MARGIN = 1, FUN = quantile, p = 0.025)
pred_mean <- pred_season[,-1] %>% apply(MARGIN = 1, FUN = mean)

pdf(file = paste0(figure_path, "sGP_seasonal_lunar.pdf"), height = 5, width = 5)
plot(mean~refined_years, type = 'l', col = 'blue', ylab = 'CO2', xlab = 'years', ylim = c(-5,5))
lines(upper_int~refined_years, lty = 'dashed', col = 'red')
lines(lower_int~refined_years, lty = 'dashed', col = 'red')
lines(pred_mean~pred_years, lty = 'solid', col = "blue")
lines(pred_upper_int~pred_years, lty = 'dashed', col = 'red')
lines(pred_lower_int~pred_years, lty = 'dashed', col = 'red')
abline(v = seq(min(refined_years),max(pred_years), by = (1/period4)), col = "green", lty = "dashed")
dev.off()

pdf(file = paste0(figure_path, "sGP_seasonal_lunar_pred.pdf"), height = 5, width = 5)
plot(mean~refined_years, type = 'l', col = 'blue', ylab = 'CO2', xlab = 'years', ylim = c(-10,10), xlim = c(2010,2030))
lines(upper_int~refined_years, lty = 'dashed', col = 'red')
lines(lower_int~refined_years, lty = 'dashed', col = 'red')
lines(pred_mean~pred_years, lty = 'solid', col = "blue")
lines(pred_upper_int~pred_years, lty = 'dashed', col = 'red')
lines(pred_lower_int~pred_years, lty = 'dashed', col = 'red')
abline(v = seq(min(refined_years),max(pred_years), by = (1/period4)), col = "green", lty = "dashed")
dev.off()

### seasonal part 4:
samps_season<- samps_g5
samps_season$x <- refined_x
upper_int <- samps_season[,-1] %>% apply(MARGIN = 1, FUN = quantile, p = 0.975)
lower_int <- samps_season[,-1] %>% apply(MARGIN = 1, FUN = quantile, p = 0.025)
mean <- samps_season[,-1] %>% apply(MARGIN = 1, FUN = mean)
pred_season<- pred_g5
pred_season$x <- pred_x
pred_upper_int <- pred_season[,-1] %>% apply(MARGIN = 1, FUN = quantile, p = 0.975)
pred_lower_int <- pred_season[,-1] %>% apply(MARGIN = 1, FUN = quantile, p = 0.025)
pred_mean <- pred_season[,-1] %>% apply(MARGIN = 1, FUN = mean)

pdf(file = paste0(figure_path, "sGP_seasonal_solar.pdf"), height = 5, width = 5)
plot(mean~refined_years, type = 'l', col = 'blue', ylab = 'CO2', xlab = 'years', ylim = c(-5,5))
lines(upper_int~refined_years, lty = 'dashed', col = 'red')
lines(lower_int~refined_years, lty = 'dashed', col = 'red')
lines(pred_mean~pred_years, lty = 'solid', col = "blue")
lines(pred_upper_int~pred_years, lty = 'dashed', col = 'red')
lines(pred_lower_int~pred_years, lty = 'dashed', col = 'red')
abline(v = seq(min(refined_years),max(pred_years), by = (1/period4)), col = "green", lty = "dashed")
dev.off()

pdf(file = paste0(figure_path, "sGP_seasonal_solar_pred.pdf"), height = 5, width = 5)
plot(mean~refined_years, type = 'l', col = 'blue', ylab = 'CO2', xlab = 'years', ylim = c(-10,10), xlim = c(2010,2030))
lines(upper_int~refined_years, lty = 'dashed', col = 'red')
lines(lower_int~refined_years, lty = 'dashed', col = 'red')
lines(pred_mean~pred_years, lty = 'solid', col = "blue")
lines(pred_upper_int~pred_years, lty = 'dashed', col = 'red')
lines(pred_lower_int~pred_years, lty = 'dashed', col = 'red')
abline(v = seq(min(refined_years),max(pred_years), by = (1/period4)), col = "green", lty = "dashed")
dev.off()


### seasonal parts lunar_solar together:
samps_season<- samps_g4 + samps_g5
samps_season$x <- refined_x
upper_int <- samps_season[,-1] %>% apply(MARGIN = 1, FUN = quantile, p = 0.975)
lower_int <- samps_season[,-1] %>% apply(MARGIN = 1, FUN = quantile, p = 0.025)
mean <- samps_season[,-1] %>% apply(MARGIN = 1, FUN = mean)
pred_season<- pred_g4 + pred_g5
pred_season$x <- pred_x
pred_upper_int <- pred_season[,-1] %>% apply(MARGIN = 1, FUN = quantile, p = 0.975)
pred_lower_int <- pred_season[,-1] %>% apply(MARGIN = 1, FUN = quantile, p = 0.025)
pred_mean <- pred_season[,-1] %>% apply(MARGIN = 1, FUN = mean)

pdf(file = paste0(figure_path, "sGP_seasonal_lunar_solar.pdf"), height = 5, width = 5)
plot(mean~refined_years, type = 'l', col = 'blue', ylab = 'CO2', xlab = 'years', ylim = c(-5,5))
lines(upper_int~refined_years, lty = 'dashed', col = 'red')
lines(lower_int~refined_years, lty = 'dashed', col = 'red')
lines(pred_mean~pred_years, lty = 'solid', col = "blue")
lines(pred_upper_int~pred_years, lty = 'dashed', col = 'red')
lines(pred_lower_int~pred_years, lty = 'dashed', col = 'red')
abline(v = seq(min(refined_years),max(pred_years), by = (1/period4)), col = "green", lty = "dashed")
dev.off()

pdf(file = paste0(figure_path, "sGP_seasonal_lunar_solar_pred.pdf"), height = 5, width = 5)
plot(mean~refined_years, type = 'l', col = 'blue', ylab = 'CO2', xlab = 'years', ylim = c(-10,10), xlim = c(2010,2030))
lines(upper_int~refined_years, lty = 'dashed', col = 'red')
lines(lower_int~refined_years, lty = 'dashed', col = 'red')
lines(pred_mean~pred_years, lty = 'solid', col = "blue")
lines(pred_upper_int~pred_years, lty = 'dashed', col = 'red')
lines(pred_lower_int~pred_years, lty = 'dashed', col = 'red')
dev.off()




### Post Summary of Hyper:
d_step <- 10
scale1 <- compute_d_step_sGPsd(d = d_step,a = a1)
c1 <- 1/scale1
scale2 <- compute_d_step_sGPsd(d = d_step,a = a2)
c2 <- 1/scale2
scale3 <- compute_d_step_sGPsd(d = d_step,a = a3)
c3 <- 1/scale3
scale4 <- compute_d_step_sGPsd(d = d_step,a = a4)
c4 <- 1/scale4
scale5 <- compute_d_step_sGPsd(d = d_step,a = a5)
c5 <- 1/scale5

prior_set_IWP_used <- prior_conversion(d = d_step, prior = prior_set_IWP, p = 3)
c6 <- prior_set_IWP$u/prior_set_IWP_used$u

load(file = paste0(result_path, "samps/", "1", "_hyper_sample.rda"))
sampsHyper1 <- simulate_hyper(hyper_result = hyper_result, nsamps = 3000, which_index = 1)
sampsHyper1 <- c1 * sampsHyper1
summary(sampsHyper1)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.5031  1.8179  2.4542  2.5480  3.1709  7.4120 

sampsHyper2 <- simulate_hyper(hyper_result = hyper_result, nsamps = 3000, which_index = 2)
sampsHyper2 <- c2 * sampsHyper2
summary(sampsHyper2)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.5031  1.8179  2.4542  2.5480  3.1709  7.4120 

sampsHyper3 <- simulate_hyper(hyper_result = hyper_result, nsamps = 3000, which_index = 3)
sampsHyper3 <- c3 * sampsHyper3
summary(sampsHyper3)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.2236  0.4004  0.4514  0.4562  0.5019  0.9460 

sampsHyper4 <- simulate_hyper(hyper_result = hyper_result, nsamps = 3000, which_index = 4)
sampsHyper4 <- c4 * sampsHyper4
summary(sampsHyper4)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.001699 0.004248 0.010436 0.016632 0.024445 0.151977 

sampsHyper5 <- simulate_hyper(hyper_result = hyper_result, nsamps = 3000, which_index = 5)
sampsHyper5 <- c5 * sampsHyper5
summary(sampsHyper5)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.004215 0.010895 0.020350 0.024036 0.033345 0.204786 

sampsHyper6 <- simulate_hyper(hyper_result = hyper_result, nsamps = 3000, which_index = 6)
sampsHyper6 <- c6 * sampsHyper6
summary(sampsHyper6)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.07827 0.43707 0.64460 0.75794 0.94617 7.32720 

sampsHyper7 <- simulate_hyper(hyper_result = hyper_result, nsamps = 3000, which_index = 7)
summary(sampsHyper7)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.5543  0.5870  0.5928  0.5930  0.5988  0.6352 





