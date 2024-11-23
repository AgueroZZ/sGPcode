library(BayesGP)
library(tidyverse)
library(parallel)
source("code/functions.R")

n <- 500
location_of_interest <- seq(0, 10, length.out = 500)
true_f <- function(x) sin(0.5* 2 * pi * x) * log(1 + x)
set.seed(123)
data <- simulate_data_poisson(func = true_f, n = n, sigma = 0.5, region = c(0,10), offset = 2, spacing = "equal")
mean(data$y == 0)

par(mfrow = c(1,2))

plot(data$x, data$y, type = "p", col = "black", 
     pch = 20, cex = 0.5,
     ylab = "y", xlab = "x")
lines(location_of_interest, exp(true_f(location_of_interest)), col = "red", lwd = 2)

plot(location_of_interest, true_f(location_of_interest),
     type = "l", col = "black",
     pch = 20, cex = 0.5,
     ylab = "y", xlab = "x")

par(mfrow = c(1,1))


u_vec = c(2, 1, 0.5, 0.1, 0.05, 0.03, 0.01, 0.005)
model_list <- list()

for (i in 1:length(u_vec)) {
  mod <- BayesGP::model_fit(
    y ~ f(
      x,
      model = "sgp",
      region = c(0,10),
      freq = 1/2,
      k = 12,
      sd.prior = list(param = list(u = u_vec[i], alpha = 0.5), h = 1)
    ) +
      f(index, model = "iid", sd.prior = 1),
    data = data,
    family = "Poisson"
  )
  model_list[[i]] <- mod
}

first_moments <- matrix(nrow = length(location_of_interest), ncol = length(model_list))
second_moments <- matrix(nrow = length(location_of_interest), ncol = length(model_list))
for (i in 1:length(model_list)) {
  mod <- model_list[[i]]
  post_samples <- predict(mod, 
                          newdata = data.frame(x = location_of_interest),
                          variable = "x",
                          only.samples = T,
                          include.intercept = F)
  first_moments[,i] <- apply(post_samples, 1, mean)
  second_moments[,i] <- apply(post_samples^2, 1, mean)
}

matplot(location_of_interest, first_moments, 
        type = "l", col = 1:length(u_vec), lwd = 1, lty = 1:length(u_vec), 
        ylim = c(min(first_moments), max(first_moments) + 3), 
        ylab = "first moment", xlab = "x")
legend("topleft", legend = u_vec, col = 1:length(u_vec), lty = 1:length(u_vec))

matplot(location_of_interest, second_moments,
        type = "l", col = 1:length(u_vec), lwd = 1, lty = 1:length(u_vec), 
        ylim = c(min(second_moments), max(second_moments) + 0), 
        ylab = "second moment", xlab = "x")
legend("topleft", legend = u_vec, col = 1:length(u_vec), lty = 1:length(u_vec))


assess_sensitivity_u_given_alpha <- function(u_vec, alpha, level = 0.95, data, location_of_interest, true_f, freq){
  # return mse_vec, coverage_vec and u_vec, in a data frame
  mse_vec <- c()
  coverage_vec <- c()
  
  for (i in 1:length(u_vec)) {
    mod <- BayesGP::model_fit(
      y ~ f(
        x,
        model = "sgp",
        region = c(0,10),
        freq = freq,
        k = 12,
        sd.prior = list(param = list(u = u_vec[i], alpha = alpha), h = 1)
      ) +
        f(index, model = "iid", sd.prior = 1),
      data = data,
      family = "Poisson"
    )
    res <- compute_mse_coverage(mod, location_of_interest, true_f, level = level)
    mse_vec <- c(mse_vec, res$mse)
    coverage_vec <- c(coverage_vec, res$coverage)
  }
  
  return(data.frame(u = u_vec, mse = mse_vec, coverage = coverage_vec))
}

alpha_0.05 <- assess_sensitivity_u_given_alpha(u_vec, 0.05, data = data, location_of_interest = location_of_interest, true_f = true_f, freq = 1/2)
alpha_0.1 <- assess_sensitivity_u_given_alpha(u_vec, 0.1, data = data, location_of_interest = location_of_interest, true_f = true_f, freq = 1/2)
alpha_0.3 <- assess_sensitivity_u_given_alpha(u_vec, 0.3, data = data, location_of_interest = location_of_interest, true_f = true_f, freq = 1/2)
alpha_0.5 <- assess_sensitivity_u_given_alpha(u_vec, 0.5, data = data, location_of_interest = location_of_interest, true_f = true_f, freq = 1/2)

save(alpha_0.05, alpha_0.1, alpha_0.3, alpha_0.5, file = "output/sensitivity_u_alpha_easy.RData")
```

Plot the rMSE and coverage rate as a function of `u` for different `alpha`

```{r}
load("output/sensitivity_u_alpha_easy.RData")
true_f_vec <- true_f(location_of_interest)
pdf(file = "output/figures/sensitivity_u_alpha_easy_mse.pdf", width = 5, height = 5)

plot(alpha_0.05$u, sqrt(alpha_0.05$mse), type = "l", 
     log = "x",
     col = 1, lwd = 2, xlab = "u", 
     lty = 1, 
     ylim = c(0.05,0.2),
     xaxt = "n",
     ylab = "rMSE",
     cex.lab = 1.2, 
     font.lab = 2,
     font = 2, 
     cex.axis = 1.2)
axis(1, at = c(0.005, 0.01, 0.05, 0.1, 0.5, 1, 2),
     cex.axis = 1.2, font.axis = 2,
     labels = c("0.005", "0.01", "0.05", "0.1", "0.5", "1", "2"))
lines(alpha_0.1$u, sqrt(alpha_0.1$mse), col = 2, lwd = 2, lty = 2)
lines(alpha_0.3$u, sqrt(alpha_0.3$mse), col = 3, lwd = 2, lty = 3)
lines(alpha_0.5$u, sqrt(alpha_0.5$mse), col = 4, lwd = 2, lty = 4)
legend("topright", title = "Prob", legend = c(0.05, 0.1, 0.3, 0.5), col = 1:4, lty = 1:4, lwd = 2, text.font = 2, title.cex = 1.2, cex = 1.2)

dev.off()



pdf(file = "output/figures/sensitivity_u_alpha_easy_cov.pdf", width = 5, height = 5)
plot(alpha_0.05$u, alpha_0.05$coverage, type = "l", 
     log = "x",
     ylim = c(0.2, 1),
     col = 1, lwd = 2, xlab = "u", 
     lty = 1,
     xaxt = "n",
     ylab = "Coverage Rate",
     cex.lab = 1.2,  # Adjust axis label size
     font.lab = 2,   # Bold axis labels
     font = 2,       # Bold text in plot
     cex.axis = 1.2) # Adjust axis tick size

# Customize x-axis ticks
axis(1, at = c(0.005, 0.01, 0.05, 0.1, 0.5, 1, 2), 
     labels = c("0.005", "0.01", "0.05", "0.1", "0.5", "1", "2"),
     cex.axis = 1.2,  # Adjust tick label size
     font.axis = 2)   # Bold tick labels

# Add lines
lines(alpha_0.1$u, alpha_0.1$coverage, col = 2, lwd = 2, lty = 2)
lines(alpha_0.3$u, alpha_0.3$coverage, col = 3, lwd = 2, lty = 3)
lines(alpha_0.5$u, alpha_0.5$coverage, col = 4, lwd = 2, lty = 4)

# Add horizontal line at 0.95
abline(h = 0.95, col = "purple", lty = 5)

# Add legend
legend("bottomright", 
       legend = c(0.05, 0.1, 0.3, 0.5), 
       col = 1:4, 
       lty = 1:4, 
       lwd = 2,
       text.font = 2,    # Bold legend text
       title = "Prob",   # Add legend title
       title.cex = 1.2,  # Adjust legend title size
       cex = 1.2)        # Adjust legend text size

dev.off()



n <- 500
location_of_interest <- seq(0, 10, length.out = 500)
true_f <- function(x){
  if(x < 2){
    return(2*sin(2 * 2 * pi * x) * (3-x))
  } else if (x > 2 && x < 4){
    return(2*sin(2 * 2 * pi * x))
  } else{
    return(2*sin(2 * 2 * pi * x) * (log(x-3) + 1))
  }
}
# vectorize the function
true_f <- Vectorize(true_f)


set.seed(123)
data <- simulate_data_poisson(func = true_f, n = n, sigma = 0.5, region = c(0,10), offset = 0)
mean(data$y == 0)

par(mfrow = c(1,2))

plot(data$x, data$y, type = "p", col = "black", 
     pch = 20, cex = 0.5,
     ylab = "y", xlab = "x")
lines(location_of_interest, exp(true_f(location_of_interest)), col = "red", lwd = 2)

plot(location_of_interest, true_f(location_of_interest),
     type = "l", col = "black",
     pch = 20, cex = 0.5,
     ylab = "y", xlab = "x")

par(mfrow = c(1,1))


u_vec = c(2, 1, 0.5, 0.1, 0.05, 0.03, 0.01, 0.005)
model_list <- list()

for (i in 1:length(u_vec)) {
  mod <- BayesGP::model_fit(
    y ~ f(
      x,
      model = "sgp",
      region = c(0,10),
      freq = 2,
      k = 12,
      sd.prior = list(param = list(u = u_vec[i], alpha = 0.5), h = 1)
    ) +
      f(index, model = "iid", sd.prior = 1),
    data = data,
    family = "Poisson"
  )
  model_list[[i]] <- mod
}

first_moments <- matrix(nrow = length(location_of_interest), ncol = length(model_list))
second_moments <- matrix(nrow = length(location_of_interest), ncol = length(model_list))
for (i in 1:length(model_list)) {
  mod <- model_list[[i]]
  post_samples <- predict(mod, 
                          newdata = data.frame(x = location_of_interest),
                          variable = "x",
                          only.samples = T,
                          include.intercept = F)
  first_moments[,i] <- apply(post_samples, 1, mean)
  second_moments[,i] <- apply(post_samples^2, 1, mean)
}

matplot(location_of_interest, first_moments, 
        type = "l", col = 1:length(u_vec), lwd = 1, lty = 1:length(u_vec), 
        ylim = c(min(first_moments), max(first_moments) + 5), 
        ylab = "first moment", xlab = "x")
legend("topright", legend = u_vec, col = 1:length(u_vec), lty = 1:length(u_vec))

matplot(location_of_interest, second_moments,
        type = "l", col = 1:length(u_vec), lwd = 1, lty = 1:length(u_vec), 
        ylim = c(min(second_moments), max(second_moments) + 0), 
        ylab = "second moment", xlab = "x")
legend("topright", legend = u_vec, col = 1:length(u_vec), lty = 1:length(u_vec))


alpha_0.05 <- assess_sensitivity_u_given_alpha(u_vec, 0.05, data = data, location_of_interest = location_of_interest, true_f = true_f, freq = 2)
alpha_0.1 <- assess_sensitivity_u_given_alpha(u_vec, 0.1, data = data, location_of_interest = location_of_interest, true_f = true_f, freq = 2)
alpha_0.3 <- assess_sensitivity_u_given_alpha(u_vec, 0.3, data = data, location_of_interest = location_of_interest, true_f = true_f, freq = 2)
alpha_0.5 <- assess_sensitivity_u_given_alpha(u_vec, 0.5, data = data, location_of_interest = location_of_interest, true_f = true_f, freq = 2)

save(alpha_0.05, alpha_0.1, alpha_0.3, alpha_0.5, file = "output/sensitivity_u_alpha_hard.RData")



load("output/sensitivity_u_alpha_hard.RData")
true_f_vec <- true_f(location_of_interest)
pdf(file = "output/figures/sensitivity_u_alpha_hard_mse.pdf", width = 5, height = 5)

plot(alpha_0.05$u, sqrt(alpha_0.05$mse), type = "l", 
     log = "x",
     col = 1, lwd = 2, xlab = "u", 
     lty = 1, 
     ylim = c(0.1,1),
     xaxt = "n",
     ylab = "rMSE",
     cex.lab = 1.2, 
     font.lab = 2,
     font = 2, 
     cex.axis = 1.2)
axis(1, at = c(0.005, 0.01, 0.05, 0.1, 0.5, 1, 2),
     cex.axis = 1.2, font.axis = 2,
     labels = c("0.005", "0.01", "0.05", "0.1", "0.5", "1", "2"))
lines(alpha_0.1$u, sqrt(alpha_0.1$mse), col = 2, lwd = 2, lty = 2)
lines(alpha_0.3$u, sqrt(alpha_0.3$mse), col = 3, lwd = 2, lty = 3)
lines(alpha_0.5$u, sqrt(alpha_0.5$mse), col = 4, lwd = 2, lty = 4)
legend("topright", title = "Prob", legend = c(0.05, 0.1, 0.3, 0.5), col = 1:4, lty = 1:4, lwd = 2, text.font = 2, title.cex = 1.2, cex = 1.2)

dev.off()



pdf(file = "output/figures/sensitivity_u_alpha_hard_cov.pdf", width = 5, height = 5)

plot(alpha_0.05$u, alpha_0.05$coverage, type = "l", 
     log = "x",
     ylim = c(0.2, 1),
     col = 1, lwd = 2, xlab = "u", 
     lty = 1,
     xaxt = "n",
     ylab = "Coverage Rate",
     cex.lab = 1.2,  # Adjust axis label size
     font.lab = 2,   # Bold axis labels
     font = 2,       # Bold text in plot
     cex.axis = 1.2) # Adjust axis tick size

# Customize x-axis ticks
axis(1, at = c(0.005, 0.01, 0.05, 0.1, 0.5, 1, 2), 
     labels = c("0.005", "0.01", "0.05", "0.1", "0.5", "1", "2"),
     cex.axis = 1.2,  # Adjust tick label size
     font.axis = 2)   # Bold tick labels

# Add lines
lines(alpha_0.1$u, alpha_0.1$coverage, col = 2, lwd = 2, lty = 2)
lines(alpha_0.3$u, alpha_0.3$coverage, col = 3, lwd = 2, lty = 3)
lines(alpha_0.5$u, alpha_0.5$coverage, col = 4, lwd = 2, lty = 4)

# Add horizontal line at 0.95
abline(h = 0.95, col = "purple", lty = 5)

# Add legend
legend("bottomright", 
       legend = c(0.05, 0.1, 0.3, 0.5), 
       col = 1:4, 
       lty = 1:4, 
       lwd = 2,
       text.font = 2,    # Bold legend text
       title = "Prob",   # Add legend title
       title.cex = 1.2,  # Adjust legend title size
       cex = 1.2)        # Adjust legend text size

dev.off()