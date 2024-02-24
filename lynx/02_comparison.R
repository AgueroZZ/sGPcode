set.seed(123)

#### PATH:
working_path <- getwd()
source_path <- paste0(working_path,"/source/")
figure_path <- paste0(working_path,"/figures/")
result_path <- paste0(working_path,"/results/")

library(tidyverse)
library(BayesGP)
library(sGPfit)

data <- data.frame(year = seq(1821, 1934, by = 1), logy = log(as.numeric(lynx)), y = as.numeric(lynx))
data$x <- data$year - min(data$year)
x <- data$x
y <- data$y
data_reduced <- data[1:80,]
### Region of prediction
region_lynx <- c(1821,1960)

### Define a prior on the marginal SD:
pred_SD <- list(u = 0.5, alpha = 0.01)
p1 <- 10

## Forecast using the comparison
Customized_RE1 <- list()
Customized_RE1$compute_B <- function(x){
  B <- Compute_B_B(x = x, k = 300, region = region_lynx, boundary = F, order = 4)
  return(B)
}
Customized_RE1$compute_P <- function(x, const){
  kappa1 = 2*pi/p1 + const
  theta1 = 2*asin(2*pi/(kappa1*p1))/pi
  scaling = sqrt(4 * cos(pi*theta1/2)*(kappa1^3))
  Prec <- (1/(scaling^2))*sGPfit:::Compute_Q_B_compare(k = 300, kappa = kappa1, theta = theta1, order = 4, simp = F, region = region_lynx, boundary = F)
  return(as(Prec, "matrix"))
}

results_const <- BayesGP::model_fit_loop(
  loop_holder = "const",
  loop_values = seq(0.001, 0.03, length.out = 20),
  formula = y ~ f(x = year, model = "customized", sd.prior = list(param = pred_SD)) +
    f(x = x, model = "IID", sd.prior = list(param = list(u = 1, alpha = 0.5))) ,
  data = data_reduced,
  family = "poisson",
  Customized_RE = Customized_RE1)
results_const
save(results_const, file = paste0("results_const.rda"))
plot(results_const$var, results_const$post, type = "o", xlab = "const", ylab = "post")
const_max = results_const$var[which.max(results_const$post)]


### Fit model:
fit_comparison_once <- function(const){
  data_reduced2 <- data_reduced
  
  kappa1 = 2*pi/p1 + const
  theta1 = 2*asin(2*pi/(kappa1*p1))/pi
  scaling = sqrt(4 * cos(pi*theta1/2)*(kappa1^3))
  
  Customized_RE1 <- list()
  Customized_RE1$compute_B <- function(x){
    B <- Compute_B_B(x = x, k = 300, region = region_lynx, boundary = F, order = 4)
    return(B)
  }
  Customized_RE1$compute_P <- function(x, kappa = kappa1, theta = theta1){
    Prec <- (1/(scaling^2)) * sGPfit:::Compute_Q_B_compare(k = 300, kappa = kappa, theta = theta, order = 4, simp = F, region = region_lynx, boundary = F)
    return(as(Prec, "matrix"))
  }
  results <- BayesGP::model_fit(
    formula = y ~ f(x = year, model = "customized", sd.prior = list(param = pred_SD)) +
      f(x = x, model = "IID", sd.prior = list(param = list(u = 1, alpha = 0.5))) ,
    data = data_reduced2,
    family = "poisson",
    Customized_RE = Customized_RE1)
  results
}
results_compare <- fit_comparison_once(const_max)
summary(results_compare)

year_vec <- data$year

full_design <- cbind(Customized_RE1$compute_B(x = year_vec), 1)
full_samps <- (full_design %*% results_compare$samps$samps[c(results_compare$random_samp_indexes$year, results_compare$fixed_samp_indexes$intercept), ])
full_summary <- list()
full_summary$year = year_vec
full_summary$mean <- (full_samps) %>% apply(1, mean)
full_summary$q0.5 <- (full_samps) %>% apply(1, median)
full_summary$q0.975 <- (full_samps) %>% apply(1, quantile, probs = 0.975)
full_summary$q0.025 <- (full_samps) %>% apply(1, quantile, probs = 0.025)


### Overall (Log):
pdf(file = paste0(figure_path, "lynx_log_predict_compare.pdf"), height = 5, width = 5)
plot(logy~year, data = data, type = 'p', xlab = "year", ylab = "lynx (log)", xlim = c(1820,1940), lwd = 0.5, cex = 1, ylim = c(0,15))
# lines(q0.5~year, data = full_summary, col = 'blue', xlab = "year", ylab = "lynx (log)", xlim = c(1820,1940))
lines(mean~year, data = full_summary, col = 'blue', xlab = "year", ylab = "lynx (log)", xlim = c(1820,1940))
lines(q0.975~year, data = full_summary, col = 'red', lty = "dashed")
lines(q0.025~year, data = full_summary, col = 'red', lty = "dashed")
abline(v = 1900, col = "purple", lty = "dashed")
dev.off()

pdf(file = paste0(figure_path, "lynx_predict_compare_samps.pdf"), height = 5, width = 5)
plot(y~year, log='y' ,data = data, type = 'p', xlab = "year", ylab = "lynx", xlim = c(1820,1940), lwd = 0.5, cex = 1, ylim = c(5,50000) ) #ylim = c(0,15))
# lines(q0.5~year, data = full_summary, col = 'blue', xlab = "year", ylab = "lynx (log)", xlim = c(1820,1940))
# lines(q0.975~year, data = full_summary, col = 'red', lty = "dashed")
# lines(q0.025~year, data = full_summary, col = 'red', lty = "dashed")
matlines(y = exp(full_samps[,2:102]), x = year_vec, col = "#FF000010", lty = 1)
Nshow = 4
matlines(y = exp(full_samps[,seq(1,Nshow)]), x = year_vec, 
         col = paste0(RColorBrewer::brewer.pal(Nshow,'Dark2'),"99"),
         lty = 1)
dev.off()


matplot(y = full_samps[,2:10], x = year_vec, type = "l", ylim = c(0,15))



### MSE and rMSE:
test_data <- data[-c(1:80),]
mean((test_data$logy - full_summary$mean[-c(1:80)])^2)
# 1.297235
sqrt(mean((test_data$logy - full_summary$mean[-c(1:80)])^2))
# 1.138962
mean(abs((test_data$logy - full_summary$mean[-c(1:80)])))
# 0.9534734




