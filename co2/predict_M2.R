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

years_cyl1 <- 1
years_cyl2 <- 0.5
years_cyl3 <- 44/12
years_cyl4 <- 9.1
years_cyl5 <- 10.4
period1 <- 1/years_cyl1
period2 <- 1/years_cyl2
period3 <- 1/years_cyl3
period4 <- 1/years_cyl4
period5 <- 1/years_cyl5
a1 = 2*period1*pi
a2 = 2*period2*pi
a3 = 2*period3*pi
a4 = 2*period4*pi
a5 = 2*period5*pi

X1 <- as(cbind(cos(a1*refined_x), sin(a1*refined_x)), "dgTMatrix")
X1new <- as(cbind(cos(a1*pred_x), sin(a1*pred_x)), "dgTMatrix")

X2 <- as(cbind(cos(a2*refined_x), sin(a2*refined_x)), "dgTMatrix")
X2new <- as(cbind(cos(a2*pred_x), sin(a2*pred_x)), "dgTMatrix")

X3 <- as(cbind(cos(a3*refined_x), sin(a3*refined_x)), "dgTMatrix")
X3new <- as(cbind(cos(a3*pred_x), sin(a3*pred_x)), "dgTMatrix")

X4 <- as(cbind(cos(a4*refined_x), sin(a4*refined_x)), "dgTMatrix")
X4new <- as(cbind(cos(a4*pred_x), sin(a4*pred_x)), "dgTMatrix")

X5 <- as(cbind(cos(a5*refined_x), sin(a5*refined_x)), "dgTMatrix")
X5new <- as(cbind(cos(a5*pred_x), sin(a5*pred_x)), "dgTMatrix")

Q6 <- as(compute_weights_precision(x = seq(0, max(x), length.out = k)), "dgTMatrix")
B6 <- as(local_poly(knots = seq(0, max(x), length.out = k), refined_x = refined_x, p = 3), "dgTMatrix")

set.seed(123)
nsamps <- 3000

### Load the corresponding samples:
## Obtain g1,g2,g3,g4,g5
samps_g1 <- data.frame(x = refined_x)
samps_g2 <- data.frame(x = refined_x)
samps_g3 <- data.frame(x = refined_x)
samps_g4 <- data.frame(x = refined_x)
samps_g5 <- data.frame(x = refined_x)
samps_g6 <- data.frame(x = refined_x)
samps_IWP_1st <- data.frame(x = refined_x)
samps_IWP_2nd <- data.frame(x = refined_x)

sd1 <- c()
sd2 <- c()

X3 <- as(cbind(cos(a3*refined_x), sin(a3*refined_x)), "dgTMatrix")
X3new <- as(cbind(cos(a3*pred_x), sin(a3*pred_x)), "dgTMatrix")
X4 <- as(cbind(cos(a4*refined_x), sin(a4*refined_x)), "dgTMatrix")
X4new <- as(cbind(cos(a4*pred_x), sin(a4*pred_x)), "dgTMatrix")
X5 <- as(cbind(cos(a5*refined_x), sin(a5*refined_x)), "dgTMatrix")
X5new <- as(cbind(cos(a5*pred_x), sin(a5*pred_x)), "dgTMatrix")

X <- as(cbind(X1,X2,X3,X4,X5, 1, refined_x, refined_x^2), "dgTMatrix") 
Xnew <- as(cbind(X1new, X2new, X3new, X4new,X5new, 1, pred_x, pred_x^2), "dgTMatrix") 
B <- cbind(B6)

load(file = paste0(result_path, "fixed_samps/", "1", "_sample.rda"))
sGP_samps$samps <- sGP_samps$samps[,(1:nsamps), drop = F]
samps_w <- sGP_samps$samps[1:ncol(B), , drop = F]
samps_beta <- sGP_samps$samps[(ncol(B) + 1):(ncol(B) + ncol(X)),, drop = F]

samps_g1_new <- as.matrix(as.matrix(X[,1:2] %*% samps_beta[1:2,, drop = F]))
samps_g1_new <- as.data.frame(samps_g1_new)
samps_g1 <- cbind(samps_g1,samps_g1_new)

samps_g2_new <- as.matrix(as.matrix(X[,3:4] %*% samps_beta[3:4,, drop = F]))
samps_g2_new <- as.data.frame(samps_g2_new)
samps_g2 <- cbind(samps_g2,samps_g2_new)

samps_g3_new <- as.matrix(as.matrix(X[,5:6] %*% samps_beta[5:6,, drop = F]))
samps_g3_new <- as.data.frame(samps_g3_new)
samps_g3 <- cbind(samps_g3,samps_g3_new)

samps_g4_new <- as.matrix(X[,7:8] %*% samps_beta[7:8,, drop = F])
samps_g4_new <- as.data.frame(samps_g4_new)
samps_g4 <- cbind(samps_g4,samps_g4_new)

samps_g5_new <- as.matrix(X[,9:10] %*% samps_beta[9:10,, drop = F])
samps_g5_new <- as.data.frame(samps_g5_new)
samps_g5 <- cbind(samps_g5,samps_g5_new)

samps_g6_new <- as.matrix(B6 %*% samps_w) + as.matrix(X[,11:13] %*% samps_beta[11:13,, drop = F])
samps_g6_new <- as.data.frame(samps_g6_new)
samps_g6 <- cbind(samps_g6,samps_g6_new)

B6_new <- as(local_poly(knots = seq(0, max(x), length.out = k), refined_x = refined_x, p = 2), "dgTMatrix")
samps_IWP_1st_new <- as.matrix(B6_new %*% samps_w) + as.matrix(cbind(rep(1,length(refined_x)),2*refined_x) %*% samps_beta[12:13,, drop = F])
samps_IWP_1st <- cbind(samps_IWP_1st, samps_IWP_1st_new)

B6_new <- as(local_poly(knots = seq(0, max(x), length.out = k), refined_x = refined_x, p = 1), "dgTMatrix")
samps_IWP_2nd_new <- as.matrix(B6_new %*% samps_w) + as.matrix(cbind(rep(2,length(refined_x))) %*% samps_beta[13,, drop = F])
samps_IWP_2nd <- cbind(samps_IWP_2nd, samps_IWP_2nd_new)

sd1 <- c(sd1, exp(-0.5*sGP_samps$theta[(1:nsamps),1]))
sd2 <- c(sd2, exp(-0.5*sGP_samps$theta[(1:nsamps),2]))

### Hyper-parameter summary:
prior_set_IWP_used <- prior_conversion(d = d_step, prior = prior_set_IWP, p = 3)
c <- prior_set_IWP$u/prior_set_IWP_used$u
# summary(c1*sd1)
# summary(sd2)

pred_g6 <- data.frame(x = pred_x)
pred_IWP_1st <- data.frame(x = pred_x)

for (i in 1:nsamps) {
  pred_g6_all_new <- predict_IWP_Var(observed_f_vec = c(samps_g6[nrow(samps_g6),(1+i)], samps_IWP_1st[nrow(samps_IWP_1st),(1+i)], samps_IWP_2nd[nrow(samps_IWP_2nd),(1+i)]),
                                     last_t = min(pred_g6$x), max_t = max(pred_g6$x), p = 3, sd = sd1[i])
  pred_g6 <- cbind(pred_g6, pred_g6_all_new[,2])
  pred_IWP_1st <- cbind(pred_IWP_1st, pred_g6_all_new[,3])
}

pred_g1 <- data.frame(x = pred_x)
pred_g2 <- data.frame(x = pred_x)
pred_g3 <- data.frame(x = pred_x)
pred_g4 <- data.frame(x = pred_x)
pred_g5 <- data.frame(x = pred_x)


pred_g1_new <- as.matrix(as.matrix(Xnew[,1:2] %*% samps_beta[1:2,, drop = F]))
pred_g1_new <- as.data.frame(pred_g1_new)
pred_g1 <- cbind(pred_g1, pred_g1_new)

pred_g2_new <- as.matrix(as.matrix(Xnew[,3:4] %*% samps_beta[3:4,, drop = F]))
pred_g2_new <- as.data.frame(pred_g2_new)
pred_g2 <- cbind(pred_g2, pred_g2_new)

pred_g3_new <- as.matrix(as.matrix(Xnew[,5:6] %*% samps_beta[5:6,, drop = F]))
pred_g3_new <- as.data.frame(pred_g3_new)
pred_g3 <- cbind(pred_g3, pred_g3_new)

pred_g4_new <- as.matrix(Xnew[,7:8] %*% samps_beta[7:8,, drop = F])
pred_g4_new <- as.data.frame(pred_g4_new)
pred_g4 <- cbind(pred_g4, pred_g4_new)

pred_g5_new <- as.matrix(Xnew[,9:10] %*% samps_beta[9:10,, drop = F])
pred_g5_new <- as.data.frame(pred_g5_new)
pred_g5 <- cbind(pred_g5, pred_g5_new)


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

pdf(file = paste0(figure_path, "fixed_overall.pdf"), height = 5, width = 5)
plot(mean~refined_years, type = 'l', col = 'blue', ylab = 'CO2', xlab = 'years', ylim = c(310,450), xlim = c(1960,2030))
lines(upper_int~refined_years, lty = 'dashed', col = 'red')
lines(lower_int~refined_years, lty = 'dashed', col = 'red')
lines(pred_mean~pred_years, lty = 'solid', col = "blue")
lines(pred_upper_int~pred_years, lty = 'dashed', col = 'red')
lines(pred_lower_int~pred_years, lty = 'dashed', col = 'red')
dev.off()
pdf(file = paste0(figure_path, "fixed_overall_pred.pdf"), height = 5, width = 5)
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

pdf(file = paste0(figure_path, "fixed_trend.pdf"), height = 5, width = 5)
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
pdf(file = paste0(figure_path, "fixed_trend_1st.pdf"), height = 5, width = 5)
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

pdf(file = paste0(figure_path, "fixed_seasonal.pdf"), height = 5, width = 5)
plot(mean~refined_years, type = 'l', col = 'blue', ylab = 'CO2', xlab = 'years', ylim = c(-10,10), xlim = c(1960,2030))
lines(upper_int~refined_years, lty = 'dashed', col = 'red')
lines(lower_int~refined_years, lty = 'dashed', col = 'red')
lines(pred_mean~pred_years, lty = 'solid', col = "blue")
lines(pred_upper_int~pred_years, lty = 'dashed', col = 'red')
lines(pred_lower_int~pred_years, lty = 'dashed', col = 'red')
dev.off()

pdf(file = paste0(figure_path, "fixed_seasonal_pred.pdf"), height = 5, width = 5)
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

pdf(file = paste0(figure_path, "fixed_seasonal_part1.pdf"), height = 5, width = 5)
plot(mean~refined_years, type = 'l', col = 'blue', ylab = 'CO2', xlab = 'years', ylim = c(-10,10))
lines(upper_int~refined_years, lty = 'dashed', col = 'red')
lines(lower_int~refined_years, lty = 'dashed', col = 'red')
lines(pred_mean~pred_years, lty = 'solid', col = "blue")
lines(pred_upper_int~pred_years, lty = 'dashed', col = 'red')
lines(pred_lower_int~pred_years, lty = 'dashed', col = 'red')
abline(v = seq(min(refined_years),max(pred_years), by = (1/period1)), col = "green", lty = "dashed")
dev.off()

pdf(file = paste0(figure_path, "fixed_seasonal_part1_pred.pdf"), height = 5, width = 5)
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

pdf(file = paste0(figure_path, "fixed_seasonal_part2.pdf"), height = 5, width = 5)
plot(mean~refined_years, type = 'l', col = 'blue', ylab = 'CO2', xlab = 'years', ylim = c(-10,10))
lines(upper_int~refined_years, lty = 'dashed', col = 'red')
lines(lower_int~refined_years, lty = 'dashed', col = 'red')
lines(pred_mean~pred_years, lty = 'solid', col = "blue")
lines(pred_upper_int~pred_years, lty = 'dashed', col = 'red')
lines(pred_lower_int~pred_years, lty = 'dashed', col = 'red')
abline(v = seq(min(refined_years),max(pred_years), by = (1/period3)), col = "green", lty = "dashed")
dev.off()

pdf(file = paste0(figure_path, "fixed_seasonal_part2_pred.pdf"), height = 5, width = 5)
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

pdf(file = paste0(figure_path, "fixed_seasonal_part3.pdf"), height = 5, width = 5)
plot(mean~refined_years, type = 'l', col = 'blue', ylab = 'CO2', xlab = 'years', ylim = c(-5,5))
lines(upper_int~refined_years, lty = 'dashed', col = 'red')
lines(lower_int~refined_years, lty = 'dashed', col = 'red')
lines(pred_mean~pred_years, lty = 'solid', col = "blue")
lines(pred_upper_int~pred_years, lty = 'dashed', col = 'red')
lines(pred_lower_int~pred_years, lty = 'dashed', col = 'red')
abline(v = seq(min(refined_years),max(pred_years), by = (1/period4)), col = "green", lty = "dashed")
dev.off()

pdf(file = paste0(figure_path, "fixed_seasonal_part3_pred.pdf"), height = 5, width = 5)
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

pdf(file = paste0(figure_path, "fixed_seasonal_part4.pdf"), height = 5, width = 5)
plot(mean~refined_years, type = 'l', col = 'blue', ylab = 'CO2', xlab = 'years', ylim = c(-5,5))
lines(upper_int~refined_years, lty = 'dashed', col = 'red')
lines(lower_int~refined_years, lty = 'dashed', col = 'red')
lines(pred_mean~pred_years, lty = 'solid', col = "blue")
lines(pred_upper_int~pred_years, lty = 'dashed', col = 'red')
lines(pred_lower_int~pred_years, lty = 'dashed', col = 'red')
abline(v = seq(min(refined_years),max(pred_years), by = (1/period5)), col = "green", lty = "dashed")
dev.off()

pdf(file = paste0(figure_path, "fixed_seasonal_part4_pred.pdf"), height = 5, width = 5)
plot(mean~refined_years, type = 'l', col = 'blue', ylab = 'CO2', xlab = 'years', ylim = c(-10,10), xlim = c(2010,2030))
lines(upper_int~refined_years, lty = 'dashed', col = 'red')
lines(lower_int~refined_years, lty = 'dashed', col = 'red')
lines(pred_mean~pred_years, lty = 'solid', col = "blue")
lines(pred_upper_int~pred_years, lty = 'dashed', col = 'red')
lines(pred_lower_int~pred_years, lty = 'dashed', col = 'red')
abline(v = seq(min(refined_years),max(pred_years), by = (1/period4)), col = "green", lty = "dashed")
dev.off()





### Post Summary of Hyper:
d_step <- 10
prior_set_IWP_used <- prior_conversion(d = d_step, prior = prior_set_IWP, p = 3)
c <- prior_set_IWP$u/prior_set_IWP_used$u

load(file = paste0(result_path, "fixed_samps/", "1", "_hyper_sample.rda"))

sampsHyperIWP <- simulate_hyper(hyper_result = hyper_result, nsamps = 3000, which_index = 1)
sampsHyperIWP <- c * sampsHyperIWP
summary(sampsHyperIWP)

sampsHyperResidual <- simulate_hyper(hyper_result = hyper_result, nsamps = 3000, which_index = 2)
summary(sampsHyperResidual)




