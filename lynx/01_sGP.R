set.seed(123)

#### PATH:
working_path <- getwd()
source_path <- paste0(working_path,"/source/")
figure_path <- paste0(working_path,"/figures/")
result_path <- paste0(working_path,"/results/")

library(tidyverse)
library(BayesGP)

data <- data.frame(year = seq(1821, 1934, by = 1), logy = log(as.numeric(lynx)), y = as.numeric(lynx))
data$x <- data$year - min(data$year)
x <- data$x
y <- data$y
data_reduced <- data[1:80,]
### Region of prediction
region_lynx <- c(1821,1960)

### Define a prior on 50-years predictive SD:
pred_SD <- list(u = 0.5, alpha = 0.01)

results_sGP <- BayesGP::model_fit(
  formula = y ~ f(x = year, model = "sgp", k = 100,
                  period = 10,
                  sd.prior = list(param = pred_SD, h = 50), 
                  initial_location = "left", region = region_lynx) +
    f(x = x, model = "IID", sd.prior = list(param = list(u = 1, alpha = 0.5))),
  data = data_reduced,
  family = "poisson")

### Forecast
pred_g1 <- predict(results_sGP, newdata = data.frame(x = x, year = data$year), variable = "year", include.intercept = T)
plot(pred_g1$mean ~ pred_g1$year, type = 'l', col = "blue", ylim = c(0,12))
# plot(pred_g1$q0.5 ~ pred_g1$year, type = 'l', col = "blue", ylim = c(0,12))
lines(pred_g1$q0.975 ~ pred_g1$year, col = "red", lty = "dashed")
lines(pred_g1$q0.025 ~ pred_g1$year, col = "red", lty = "dashed")
pred_g1_samps <- predict(results_sGP, newdata = data.frame(x = x, year = data$year), variable = "year", include.intercept = T, only.samples = T)


### Overall (Log):
pdf(file = paste0(figure_path, "lynx_log_predict.pdf"), height = 5, width = 5)
plot(logy~year, data = data, type = 'p', xlab = "year", ylab = "lynx (log)", xlim = c(1820,1940), lwd = 0.5, cex = 1, ylim = c(0,15))
lines(mean~year, data = pred_g1, col = 'blue', xlab = "year", ylab = "lynx (log)", xlim = c(1820,1940))
# lines(q0.5~year, data = pred_g1, col = 'blue', xlab = "year", ylab = "lynx (log)", xlim = c(1820,1940))
lines(q0.975~year, data = pred_g1, col = 'red', lty = "dashed")
lines(q0.025~year, data = pred_g1, col = 'red', lty = "dashed")
abline(v = 1900, col = "purple", lty = "dashed")
dev.off()

pdf(file = paste0(figure_path, "lynx_predict_samps.pdf"), height = 5, width = 5)
plot(y~year, log='y' ,data = data, type = 'p', xlab = "year", ylab = "lynx", xlim = c(1820,1940), lwd = 0.5, cex = 1, ylim = c(5,50000) ) #ylim = c(0,15))
matlines(y = exp(pred_g1_samps[,2:102]), x = data$year, col = "#FF000010", lty = 1)
Nshow = 4
matlines(y = exp(pred_g1_samps[,seq(1,Nshow)]), x = data$year, 
         col = paste0(RColorBrewer::brewer.pal(Nshow,'Dark2'),"99"),
         lty = 1)
dev.off()



### Overall:
pred_g1_samps <- predict(results_sGP, newdata = data.frame(year = data$year), variable = "year", include.intercept = T, only.samples = T)
pred_g1 <- data.frame(year = pred_g1_samps[,1], med = apply(exp(pred_g1_samps[,-1]), 1, quantile, probs = 0.5), q0.025 = apply(exp(pred_g1_samps[,-1]), 1, quantile, probs = 0.025), q0.975 = apply(exp(pred_g1_samps[,-1]), 1, quantile, probs = 0.975))
plot(y~year, data = data, type = 'p',
     log = "y", ylim = c(1, 160000),
     xlab = "year", ylab = "lynx", xlim = c(1820,1971),
     lwd = 0.5, cex = 0.2)
lines(pred_g1$med ~ pred_g1$year, type = 'l', col = "blue")
lines(pred_g1$q0.975 ~ pred_g1$year, col = "red", lty = "dashed")
lines(pred_g1$q0.025 ~ pred_g1$year, col = "red", lty = "dashed")


pred_g1_samps <- predict(results_sGP, newdata = data.frame(year = data$year), variable = "year", include.intercept = T, only.samples = T)
pred_g1 <- data.frame(year = pred_g1_samps[,1], 
                      mean = apply((pred_g1_samps[,-1]), 1, mean),
                      med = apply((pred_g1_samps[,-1]), 1, quantile, probs = 0.5), 
                      q0.025 = apply((pred_g1_samps[,-1]), 1, quantile, probs = 0.025), 
                      q0.975 = apply((pred_g1_samps[,-1]), 1, quantile, probs = 0.975))
### MSE rMSE MAE:
test_data <- data[-c(1:80),]
mean((test_data$logy -pred_g1$med[-c(1:80)])^2)
# 0.9280063
sqrt(mean((test_data$logy -pred_g1$med[-c(1:80)])^2))
# 0.9633308
mean(abs(test_data$logy -pred_g1$med[-c(1:80)]))
# 0.7984736
