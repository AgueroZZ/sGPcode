.libPaths( c( .libPaths(), "~/./lib") )
library(tidyverse)
library(BayesGP)
set.seed(123)

### Read in the full data:
cUrl = paste0("http://scrippsco2.ucsd.edu/assets/data/atmospheric/",
              "stations/flask_co2/daily/daily_flask_co2_mlo.csv")
cFile = basename(cUrl)
if (!file.exists(cFile)) download.file(cUrl, cFile)
co2s = read.table(cFile, header = FALSE, sep = ",",
                  skip = 69, stringsAsFactors = FALSE, col.names = c("day",
                                                                     "time", "junk1", "junk2", "Nflasks", "quality",
                                                                     "co2"))
co2s$date = strptime(paste(co2s$day, co2s$time), format = "%Y-%m-%d %H:%M",
                     tz = "UTC")
co2s$day = as.Date(co2s$date)
timeOrigin = as.Date("1960/03/30")
co2s$timeYears = round(as.numeric(co2s$day - timeOrigin)/365.25,
                       3)
co2s$dayInt = as.integer(co2s$day)
allDays = seq(from = min(co2s$day), to = max(co2s$day),
              by = "7 day")
observed_dataset <- co2s %>% filter(!is.na(co2s$co2)) %>% dplyr::select(c("co2", "timeYears", "day"))
observed_dataset$quality <- ifelse(co2s$quality > 0, 1, 0)
# remove low-quality measurements
observed_dataset <- observed_dataset %>% filter(quality == 0)
train_dataset <- observed_dataset %>% filter(timeYears < 30)

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

d_step <- 10
prior_set_IWP <- list(alpha = 0.5, u = 30)
prior_set_sGP <- list(alpha = 0.5, u = 1)


x_full <- observed_dataset$timeYears
x <- train_dataset$timeYears
region <- c(0,max(x))


# Full data
all_data <- data.frame(y = train_dataset$co2, timeYears = train_dataset$timeYears, x1 = x, x2 = x, x3 = x, x4 = x, x5 = x)

mod_full <- model_fit(formula = y ~ f(x1, model = "sgp", a = a1, sd.prior = list(param = prior_set_sGP, h = d_step), region = range(x_full)) +
                        f(x2, model = "sgp", a = a2, k = 30, sd.prior = list(param = prior_set_sGP, h = d_step), region = range(x_full)) +
                        f(x3, model = "sgp", a = a3, k = 30, sd.prior = list(param = prior_set_sGP, h = d_step), region = range(x_full)) +
                        f(x4, model = "sgp", a = a4, k = 30, sd.prior = list(param = prior_set_sGP, h = d_step), region = range(x_full)) +
                        f(x5, model = "sgp", a = a5, k = 30, sd.prior = list(param = prior_set_sGP, h = d_step), region = range(x_full)) +
                        f(timeYears, model = "iwp", order = 3, 
                          initial_location = "left",
                          sd.prior = list(param = prior_set_IWP, h = d_step),
                          knots = seq(min(x_full), max(x_full), length.out = 50)), 
                      data = all_data, aghq_k = 4, family = "gaussian", M = 6000,
                      control.family = list(sd.prior = list(param = list(alpha = 0.5, u = 1))))
summary(mod_full)
# plot(mod_full)
save(mod_full, file = "mod_M1.rda")

set.seed(123)

mod_pred_iwp_samps <- predict(mod_full, variable = "timeYears",
                              newdata = data.frame(timeYears = x_full), only.samples = T)

X1_samps <- predict(mod_full, variable = "x1",
                    newdata = data.frame(x1 = x_full), only.samples = T, include.intercept = F)

X2_samps <- predict(mod_full, variable = "x2",
                    newdata = data.frame(x2 = x_full), only.samples = T, include.intercept = F)

X3_samps <- predict(mod_full, variable = "x3",
                    newdata = data.frame(x3 = x_full), only.samples = T, include.intercept = F)

X4_samps <- predict(mod_full, variable = "x4",
                    newdata = data.frame(x4 = x_full), only.samples = T, include.intercept = F)

X5_samps <- predict(mod_full, variable = "x5",
                    newdata = data.frame(x5 = x_full), only.samples = T, include.intercept = F)

sd_noise <- 1/exp(0.5*mod_full$samps$theta[,7])
noise_pred <- matrix(nrow = nrow(mod_pred_iwp_samps[,-1]), ncol = ncol(mod_pred_iwp_samps[,-1]))
for (ii in 1:length(sd_noise)) {
  noise_pred[,ii] <- rnorm(n = 2323, sd = sd_noise[ii])
}

M1_pred <- mod_pred_iwp_samps[,-1] + X1_samps[,-1] + X2_samps[,-1] + X3_samps[,-1] + X4_samps[,-1] + X5_samps[,-1] + noise_pred

mod_pred <- data.frame(timeYears = x_full, lower = apply(M1_pred, 1, quantile, probs = 0.1), 
                       mean = apply(M1_pred, 1, mean),
                       med = apply(M1_pred, 1, quantile, probs = 0.5), 
                       upper = apply(M1_pred, 1, quantile, probs = 0.9))

mod_pred$timeYears = observed_dataset$day

png(filename = "fit_M1.png", width = 800, height = 800)
plot(mod_pred$mean ~ mod_pred$timeYears, 
     type = "l", col = "blue", lwd = 1, 
     xlab = "Years", 
     ylab = "CO2 concentration", ylim = c(200, 620),
     cex.lab = 2.0, cex.axis = 2.0, cex.main = 1.5, cex.sub = 1.5)
lines(mod_pred$lower ~ mod_pred$timeYears, col = "red", lwd = 1, lty = 2)
lines(mod_pred$upper ~ mod_pred$timeYears, col = "red", lwd = 1, lty = 2)
points(observed_dataset$day, observed_dataset$co2, 
       pch = 16, col = "black", cex = 0.2)
abline(v = max(train_dataset$day), col = "purple", lty = "dashed")
dev.off()

pdf(file = "fit_M1.pdf", width = 8, height = 8)
plot(mod_pred$mean ~ mod_pred$timeYears, 
     type = "l", col = "blue", lwd = 1, 
     xlab = "Years", 
     ylab = "CO2 concentration", ylim = c(200, 620),
     cex.lab = 2.0, cex.axis = 2.0, cex.main = 1.5, cex.sub = 1.5)
lines(mod_pred$lower ~ mod_pred$timeYears, col = "red", lwd = 1, lty = 2)
lines(mod_pred$upper ~ mod_pred$timeYears, col = "red", lwd = 1, lty = 2)
points(observed_dataset$day, observed_dataset$co2, 
       pch = 16, col = "black", cex = 0.2)
abline(v = max(train_dataset$day), col = "purple", lty = "dashed")
dev.off()



### Compute training rMSE and CPR
mod_pred_train <- mod_pred %>% filter(timeYears <= max(train_dataset$day))
sqrt(mean((mod_pred_train$mean - train_dataset$co2)^2))
# 0.5653827
mean((mod_pred_train$lower <= train_dataset$co2) & (train_dataset$co2 <= mod_pred_train$upper))
# 0.8385714

### Compute test rMSE and CPR
mod_pred_test <- mod_pred %>% filter(timeYears > max(train_dataset$day))
test_data <- observed_dataset %>% filter(timeYears > max(train_dataset$timeYears))
sqrt(mean((mod_pred_test$mean - test_data$co2)^2))
# 12.15266
mean((mod_pred_test$lower <= test_data$co2) & (test_data$co2 <= mod_pred_test$upper))
# 0.8626001




#### Each component:
x_refined <- seq(from = min(observed_dataset$timeYears), to = max(observed_dataset$timeYears), by = 1/365.25)

X1_samps_refined <- predict(mod_full, variable = "x1",
                    newdata = data.frame(x1 = x_refined), only.samples = T, include.intercept = F)

X2_samps_refined <- predict(mod_full, variable = "x2",
                    newdata = data.frame(x2 = x_refined), only.samples = T, include.intercept = F)

X3_samps_refined <- predict(mod_full, variable = "x3",
                    newdata = data.frame(x3 = x_refined), only.samples = T, include.intercept = F)

X4_samps_refined <- predict(mod_full, variable = "x4",
                    newdata = data.frame(x4 = x_refined), only.samples = T, include.intercept = F)

X5_samps_refined <- predict(mod_full, variable = "x5",
                    newdata = data.frame(x5 = x_refined), only.samples = T, include.intercept = F)



M1_pred_seasonal <- X1_samps_refined[,-1] + X2_samps_refined[,-1] + X3_samps_refined[,-1] + X4_samps_refined[,-1] + X5_samps_refined[,-1]
mod_pred_seasonal <- data.frame(timeYears = x_refined, 
                                lower = apply(M1_pred_seasonal, 1, quantile, probs = 0.1), 
                                mean = apply(M1_pred_seasonal, 1, mean),
                                med = apply(M1_pred_seasonal, 1, quantile, probs = 0.5), 
                                upper = apply(M1_pred_seasonal, 1, quantile, probs = 0.9))

mod_pred_seasonal$timeYears = mod_pred_seasonal$timeYears * 365.25 + timeOrigin

png(filename = "seasonal_M1.png", width = 800, height = 800)
plot(mod_pred_seasonal$mean ~ mod_pred_seasonal$timeYears, 
     type = "l", col = "blue", lwd = 1, 
     xlab = "Years", 
     ylab = "",
     ylim = c(-10,10),
     cex.lab = 2.0, cex.axis = 2.0, cex.main = 1.5, cex.sub = 1.5)
lines(mod_pred_seasonal$lower ~ mod_pred_seasonal$timeYears, col = "red", lwd = 1, lty = 2)
lines(mod_pred_seasonal$upper ~ mod_pred_seasonal$timeYears, col = "red", lwd = 1, lty = 2)
abline(v = max(train_dataset$day), col = "purple", lty = "dashed")
dev.off()

png(filename = "seasonal_M1_zoomed.png", width = 800, height = 800)
plot(mod_pred_seasonal$mean ~ mod_pred_seasonal$timeYears, 
     type = "l", col = "blue", lwd = 1, 
     xlab = "Years", 
     ylab = "",
     ylim = c(-10,10),
     xlim = as.Date(c("1980-01-01","2010-01-01")),
     cex.lab = 2.0, cex.axis = 2.0, cex.main = 1.5, cex.sub = 1.5)
lines(mod_pred_seasonal$lower ~ mod_pred_seasonal$timeYears, col = "red", lwd = 1, lty = 2)
lines(mod_pred_seasonal$upper ~ mod_pred_seasonal$timeYears, col = "red", lwd = 1, lty = 2)
abline(v = max(train_dataset$day), col = "purple", lty = "dashed")
dev.off()


mod_pred_refined <- predict(mod_full, variable = "timeYears",
                            newdata = data.frame(timeYears = x_refined),
                            quantiles = c(0.1, 0.5, 0.9))
mod_pred_refined$timeYears = mod_pred_refined$timeYears * 365.25 + timeOrigin


png(filename = "trend_M1.png", width = 800, height = 800)
plot(mod_pred_refined$mean ~ mod_pred_refined$timeYears, 
     type = "l", col = "blue", lwd = 1, 
     xlab = "Years", 
     ylab = "",
     ylim = c(280,480),
     cex.lab = 2.0, cex.axis = 2.0, cex.main = 1.5, cex.sub = 1.5)
lines(mod_pred_refined$q0.1 ~ mod_pred_refined$timeYears, col = "red", lwd = 1, lty = 2)
lines(mod_pred_refined$q0.9 ~ mod_pred_refined$timeYears, col = "red", lwd = 1, lty = 2)
abline(v = max(train_dataset$day), col = "purple", lty = "dashed")
dev.off()


mod_pred_refined_deriv <- predict(mod_full, variable = "timeYears",
                            newdata = data.frame(timeYears = x_refined),
                            deriv = 1,
                            quantiles = c(0.1, 0.5, 0.9))
mod_pred_refined_deriv$timeYears = mod_pred_refined_deriv$timeYears * 365.25 + timeOrigin


png(filename = "trend_deriv_M1.png", width = 800, height = 800)
plot(mod_pred_refined_deriv$mean ~ mod_pred_refined_deriv$timeYears, 
     type = "l", col = "blue", lwd = 1, 
     xlab = "Years", 
     ylab = "",
     ylim = c(-5,8),
     cex.lab = 2.0, cex.axis = 2.0, cex.main = 1.5, cex.sub = 1.5)
lines(mod_pred_refined_deriv$q0.1 ~ mod_pred_refined_deriv$timeYears, col = "red", lwd = 1, lty = 2)
lines(mod_pred_refined_deriv$q0.9 ~ mod_pred_refined_deriv$timeYears, col = "red", lwd = 1, lty = 2)
abline(v = max(train_dataset$day), col = "purple", lty = "dashed")
dev.off()



