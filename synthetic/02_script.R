library(BayesGP)
library(sGPfit)
library(tidyverse)
library(parallel)
source("code/functions.R")



location_of_interest <- seq(0, 10, length.out = 500)
true_f <- function(x) sin(0.5 * 2 * pi * x) * log(1 + x)
n_vec <- c(100, 200, 300, 500, 800, 1000, 2000, 3000, 5000)
time_vec <- numeric(length(n_vec))
num_replicate <- 10


runtime <- matrix(NA, nrow = length(n_vec), ncol = num_replicate)
for(i in 1:length(n_vec)){
  data <- simulate_data_poisson(func = true_f, n = n_vec[i], sigma = 0.5, region = c(0,10), offset = 2, spacing = "equal")
  result <- microbenchmark::microbenchmark(
    mod <- BayesGP::model_fit(y ~ f(x, model = "sgp", freq = 1/2, 
                                    k = (30/3), # 10 knots to 30 basis!
                                    range = c(0, 10),
                                    sd.prior = list(param = 1, h = 1)) +
                                f(index, model = "iid", sd.prior = 1),
                              data = data,
                              family = "Poisson"),
    times = num_replicate
  )
  runtime[i, ] <- result$time/1e9
}
save(runtime, file = "output/runtime_FEM.RData")


runtime <- matrix(NA, nrow = length(n_vec), ncol = num_replicate)
for(i in 1:length(n_vec)){
  data <- simulate_data_poisson(func = true_f, n = n_vec[i], sigma = 0.5, region = c(0,10), offset = 2, spacing = "equal")
  result <- microbenchmark::microbenchmark(
    mod <- BayesGP::model_fit(y ~ f(x, model = "sgp", freq = 1/2, 
                                    k = (60/3), # 20 knots to 60 basis!
                                    range = c(0, 10),
                                    sd.prior = list(param = 1, h = 1)) +
                                f(index, model = "iid", sd.prior = 1),
                              data = data,
                              family = "Poisson"),
    times = num_replicate
  )
  runtime[i, ] <- result$time/1e9
}
save(runtime, file = "output/runtime_FEM2.RData")





true_f_hard <- function(x){
  if(x < 2){
    return(2*sin(2 * 2 * pi * x) * (3-x))
  } else if (x > 2 && x < 4){
    return(2*sin(2 * 2 * pi * x))
  } else{
    return(2*sin(2 * 2 * pi * x) * (log(x-3) + 1))
  }
}
true_f_hard <- Vectorize(true_f_hard)
runtime <- matrix(NA, nrow = length(n_vec), ncol = num_replicate)
for(i in 1:length(n_vec)){
  data <- simulate_data_poisson(func = true_f_hard, n = n_vec[i], sigma = 0.8, region = c(0,10), offset = 2, spacing = "equal")
  result <- microbenchmark::microbenchmark(
    mod <- BayesGP::model_fit(y ~ f(x, model = "sgp", freq = 1/2, 
                                    k = (30/3) # 10 knots to 30 basis!
                                    range = c(0, 10),
                                    sd.prior = list(param = 1, h = 1)) +
                                f(index, model = "iid", sd.prior = 1),
                              data = data,
                              family = "Poisson"),
    times = num_replicate
  )
  runtime[i, ] <- result$time/1e9
}
save(runtime, file = "output/runtime_FEM_hard.RData")

runtime <- matrix(NA, nrow = length(n_vec), ncol = num_replicate)
for(i in 1:length(n_vec)){
  data <- simulate_data_poisson(func = true_f_hard, n = n_vec[i], sigma = 0.8, region = c(0,10), offset = 2, spacing = "equal")
  result <- microbenchmark::microbenchmark(
    mod <- BayesGP::model_fit(y ~ f(x, model = "sgp", freq = 1/2, 
                                    k = (60/3), # 20 knots to 60 basis!
                                    range = c(0, 10),
                                    sd.prior = list(param = 1, h = 1)) +
                                f(index, model = "iid", sd.prior = 1),
                              data = data,
                              family = "Poisson"),
    times = num_replicate
  )
  runtime[i, ] <- result$time/1e9
}
save(runtime, file = "output/runtime_FEM2_hard.RData")



load("output/runtime_FEM.RData")
runtime_df <- data.frame(n = n_vec, time = rowMeans(runtime), k = 30)
ggplot(runtime_df, aes(x = n, y = time)) +
  geom_point() +
  geom_line() +
  # scale_x_log10() +
  # scale_y_log10() +
  labs(x = "Number of Observations", y = "Runtime (s)") +
  theme_minimal()

# put into the same figure
load("output/runtime_FEM2.RData")
runtime_df2 <- data.frame(n = n_vec, time = rowMeans(runtime), k = 60)
runtime_df_combined <- rbind(runtime_df, runtime_df2)

ggplot(runtime_df_combined, aes(x = n, y = time, color = factor(k))) +
  geom_point() +
  geom_line() +
  labs(x = "Number of Observations", y = "Runtime (s)") +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.7, 0.7), 
    axis.title = element_text(size = 14), # Axis labels text size
    axis.text = element_text(size = 14, face = "bold"),  # Axis breaks text size
    legend.text = element_text(size = 15), # Legend labels text size
  ) +
  scale_color_discrete(labels = c("k = 30", "k = 60"))
ggsave("output/figures/runtime_FEM.pdf", width = 5, height = 5)



load("output/runtime_FEM_hard.RData")
runtime_df <- data.frame(n = n_vec, time = rowMeans(runtime), k = 30)
ggplot(runtime_df, aes(x = n, y = time)) +
  geom_point() +
  geom_line() +
  labs(x = "Number of Observations", y = "Runtime (s)") +
  theme_minimal()

# put into the same figure
load("output/runtime_FEM2_hard.RData")
runtime_df2 <- data.frame(n = n_vec, time = rowMeans(runtime), k = 60)
runtime_df_combined <- rbind(runtime_df, runtime_df2)

ggplot(runtime_df_combined, aes(x = n, y = time, color = factor(k))) +
  geom_point() +
  geom_line() +
  labs(x = "Number of Observations", y = "Runtime (s)") +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.5, 0.8), 
    axis.title = element_text(size = 14), # Axis labels text size
    axis.text = element_text(size = 14, face = "bold"),  # Axis breaks text size
    legend.text = element_text(size = 15), # Legend labels text size
  ) +
  scale_color_discrete(labels = c("k = 30", "k = 60"))
ggsave("output/figures/runtime_FEM_hard.pdf", width = 5, height = 5)



n <- 500
set.seed(123)
data <- simulate_data_poisson(func = true_f, n = n, sigma = 0.5, region = c(0,10), offset = 2, spacing = "equal")
sorted_data <- data %>% arrange(x)
prior_PSD <- list(u = 1, alpha = 0.5)
prior_SD <- list(u = prior_PSD$u/compute_d_step_sGPsd(d = 1, a = pi), alpha = prior_PSD$alpha)
mod_exact <- fit_exact(data, prior_SD)
samps_exact_overall <- sample_exact(mod_exact, M = 3000)
samps_exact_overall_df <- data.frame(x = sorted_data$x)
samps_exact_overall_df$mean <- rowMeans(samps_exact_overall)
samps_exact_overall_df$second_moment <- rowMeans(samps_exact_overall^2)
samps_exact_overall_df$lower <- apply(samps_exact_overall, 1, quantile, probs = 0.025)
samps_exact_overall_df$upper <- apply(samps_exact_overall, 1, quantile, probs = 0.975)

ggplot(samps_exact_overall_df, aes(x = x)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey", alpha = 0.5) +
  geom_line(aes(y = mean), color = "black", linetype = "solid") +
  geom_line(aes(y = true_f(x)), color = "red", linetype = "dotted") +
  labs(x = "x", y = "f(x)") +
  theme_minimal() +
  guides(fill = FALSE)
samps_approx_df_list <- list()
k_vec <- c(4, 5, 6, 7, 8, 9, 10, 12, 17, 22, 27, 32, 37)
num_basis <- 3*(k_vec -2)
for (i in 1:length(k_vec)){
  mod_approx <- BayesGP::model_fit(y ~ f(x, model = "sgp", freq = 1/2, 
                                         k = k_vec[i],
                                         range = c(0, 10),
                                         sd.prior = list(param = 1, h = 1)) +
                                     f(index, model = "iid", sd.prior = 1),
                                   data = data,
                                   family = "Poisson")
  
  samps_approx <- predict(mod_approx, newdata = data.frame(x = sort(data$x)), only.samples = TRUE, variable = "x", quantiles = NULL, include.intercept = F)
  
  samps_approx_df_list[[i]] <- data.frame(x = data$x, 
                                          mean = rowMeans(samps_approx), 
                                          second_moment = rowMeans(samps_approx^2),
                                          lower = apply(samps_approx, 1, quantile, probs = 0.025),
                                          upper = apply(samps_approx, 1, quantile, probs = 0.975))
  
}
save(samps_approx_df_list, file = "output/samps_approx_df_list_FEM.RData")
moment_matrix <- do.call(rbind, lapply(samps_approx_df_list, function(x) x$mean))
second_moment_matrix <- do.call(rbind, lapply(samps_approx_df_list, function(x) x$second_moment))

par(mfrow = c(2, 1))

matplot(data$x, t(moment_matrix), type = "l", lty = 1, col = 1:ncol(moment_matrix), xlab = "x", ylab = "First moment")
lines(samps_exact_overall_df$x, samps_exact_overall_df$mean, col = "red", lty = 2)

matplot(data$x, t(second_moment_matrix), type = "l", lty = 1, col = 1:ncol(second_moment_matrix), xlab = "x", ylab = "Second moment")
lines(samps_exact_overall_df$x, samps_exact_overall_df$second_moment, col = "red", lty = 2)

par(mfrow = c(1, 1))

error_to_exact_first_moment <- moment_matrix %>% apply(1, function(x) max((x - samps_exact_overall_df$mean)^2))
error_to_exact_second_moment <- second_moment_matrix %>% apply(1, function(x) max((x - samps_exact_overall_df$second_moment)^2))
error_df <- data.frame(num_basis = num_basis, error_to_exact_first_moment = error_to_exact_first_moment, error_to_exact_second_moment = error_to_exact_second_moment)

# Put them into one plot
error_df_long <- error_df %>% pivot_longer(cols = c(error_to_exact_first_moment, error_to_exact_second_moment), names_to = "moment", values_to = "error")
ggplot(error_df_long, aes(x = num_basis, y = error, color = moment)) +
  geom_point() +
  geom_line() +
  labs(x = "Number of Basis Functions", y = "Maximal Error") +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.6, 0.6), 
    axis.title = element_text(size = 14), # Axis labels text size
    axis.text = element_text(size = 14, face = "bold"),  # Axis breaks text size
    legend.text = element_text(size = 15), # Legend labels text size
  ) +
  scale_color_discrete(labels = c("First Moment", "Second Moment"))
ggsave("output/figures/error_FEM_convergence.pdf", width = 5, height = 5)


n <- 500
set.seed(123)
data <- simulate_data_poisson(func = true_f, n = n, sigma = 0.5, region = c(0,10), offset = 2, spacing = "equal")
## Runtime of the exact method:
result_exact <- microbenchmark::microbenchmark(
  fit_exact(data, prior_SD),
  times = num_replicate
)
runtime_exact <- result_exact$time/1e9

## Runtime of the approximate method:
runtime_approx <- matrix(NA, nrow = length(k_vec), ncol = num_replicate)
for (i in 1:length(k_vec)){
  result_approx <- microbenchmark::microbenchmark(
    BayesGP::model_fit(y ~ f(x, model = "sgp", freq = 1/2, 
                             k = k_vec[i],
                             range = c(0, 10),
                             sd.prior = list(param = 1, h = 1)) +
                         f(index, model = "iid", sd.prior = 1),
                       data = data,
                       family = "Poisson"),
    times = num_replicate
  )
  runtime_approx[i, ] <- result_approx$time/1e9
}
save(runtime_exact, runtime_approx, file = "output/runtime_FEM_convergence_n500.RData")
load("output/runtime_FEM_convergence_n500.RData")
runtime_df <- data.frame(num_basis = num_basis, time = rowMeans(runtime_approx)/mean(runtime_exact))
ggplot(runtime_df, aes(x = num_basis, y = time)) +
  geom_point() +
  geom_line() +
  labs(x = "Number of Basis Functions", y = "Relative Runtime") +
  theme_minimal()



n <- 1000
set.seed(123)
data <- simulate_data_poisson(func = true_f, n = n, sigma = 0.5, region = c(0,10), offset = 2, spacing = "equal")
```

```{r eval=FALSE}
## Runtime of the exact method:
result_exact <- microbenchmark::microbenchmark(
  fit_exact(data, prior_SD),
  times = num_replicate
)
runtime_exact <- result_exact$time/1e9

## Runtime of the approximate method:
runtime_approx <- matrix(NA, nrow = length(k_vec), ncol = num_replicate)
for (i in 1:length(k_vec)){
  result_approx <- microbenchmark::microbenchmark(
    BayesGP::model_fit(y ~ f(x, model = "sgp", freq = 1/2, 
                             k = k_vec[i],
                             range = c(0, 10),
                             sd.prior = list(param = 1, h = 1)) +
                         f(index, model = "iid", sd.prior = 1),
                       data = data,
                       family = "Poisson"),
    times = num_replicate
  )
  runtime_approx[i, ] <- result_approx$time/1e9
}
save(runtime_exact, runtime_approx, file = "output/runtime_FEM_convergence_n1000.RData")

load("output/runtime_FEM_convergence_n1000.RData")
runtime_df <- data.frame(num_basis = num_basis, time = rowMeans(runtime_approx)/mean(runtime_exact))
ggplot(runtime_df, aes(x = num_basis, y = time)) +
  geom_point() +
  geom_line() +
  labs(x = "Number of Basis Functions", y = "Relative Runtime") +
  theme_minimal()


load("output/runtime_FEM_convergence_n500.RData")
runtime_df500 <- data.frame(num_basis = num_basis, time = rowMeans(runtime_approx)/mean(runtime_exact))
runtime_df500$n <- rep(500, length(num_basis))
load("output/runtime_FEM_convergence_n1000.RData")
runtime_df1000 <- data.frame(num_basis = num_basis, time = rowMeans(runtime_approx)/mean(runtime_exact))
runtime_df1000$n <- rep(1000, length(num_basis))
runtime_df_combined <- rbind(runtime_df500, runtime_df1000)

ggplot(runtime_df_combined, aes(x = num_basis, y = time, color = factor(n))) +
  geom_point() +
  geom_line() +
  labs(x = "Number of Basis Functions", y = "Relative Runtime %") +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.7, 0.7), 
    axis.title = element_text(size = 14), # Axis labels text size
    axis.text = element_text(size = 14, face = "bold"),  # Axis breaks text size
    legend.text = element_text(size = 15), # Legend labels text size
  ) +
  scale_color_discrete(labels = c("n = 500", "n = 1000"))
ggsave("output/figures/runtime_FEM_convergence.pdf", width = 5, height = 5)





