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
test_data <- data[-c(1:80),]
### Region of prediction
region_lynx <- c(1821,1960)
fit_once <- function(u, alpha){
  pred_SD <- list(u = u, alpha = alpha)
  results_sGP <- BayesGP::model_fit(
    formula = y ~ f(x = year, model = "sgp", k = 100,
                    period = 10,
                    sd.prior = list(param = pred_SD, h = 50), 
                    initial_location = "left", region = region_lynx) +
      f(x = x, model = "IID", sd.prior = list(param = list(u = 1, alpha = 0.5))),
    data = data_reduced,
    family = "poisson")
  pred_g1 <- predict(results_sGP, newdata = data.frame(x = x, year = data$year), 
                     only.samples = T,
                     variable = "year", 
                     include.intercept = F, quantiles = NULL)
  return(pred_g1)
}

alpha = 0.5
u_vec = c(0.05,0.1,0.15,0.2,0.25)
pred_summary <- lapply(u_vec, fit_once, alpha = alpha)
pred_first_moment <- do.call(rbind, lapply(pred_summary, function(x) x[,-1] %>% apply(1, mean)))
pred_second_moment <- do.call(rbind, lapply(pred_summary, function(x) x[,-1]^2 %>% apply(1, mean)))

pdf(paste0(figure_path,"lynx_sensitivity_first.pdf"), width = 5, height = 5)
matplot(data$year, t(pred_first_moment), col = 1:length(u_vec), lty = 1:length(u_vec), type = 'l',
        xlab = "Year", ylab = "Post First Moment", cex = 1.5)
dev.off()

pdf(paste0(figure_path,"lynx_sensitivity_second.pdf"), width = 5, height = 5)
matplot(data$year, t(pred_second_moment), col = 1:length(u_vec), lty = 1:length(u_vec), type = 'l',
        xlab = "Year", ylab = "Post Second Moment", cex = 1.5)
legend("topleft", legend = paste(u_vec), col = 1:length(u_vec), lty = 1:length(u_vec), cex = 1,
       title = "Prior Median")
dev.off()
