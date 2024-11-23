simulate_data_poisson <- function(func, n = 500, sigma = 1, region = c(0,10), offset = 0, spacing = "random"){
  # Simulate the data
  if(spacing == "random"){
    x <- runif(n, min = region[1], max = region[2])
  } else{
    x <- seq(region[1], region[2], length.out = n)
  }

  truef <- func(x)
  eta <- truef + offset

  if(sigma > 0){
    eta <- eta + rnorm(n, sd = sigma)
  }

  lambda <- exp(eta)

  y <- rpois(n, lambda)

  # Return the data
  return(data.frame(x = x, y = y, offset = offset, index = 1:n, truef = truef))
}

compute_mse <- function(mod, locations, func){
  pred <- predict(mod, newdata = data.frame(x = locations), variable = "x", quantiles = NULL)
  mse <- mean((func(locations) - pred$mean)^2)
  return(mse)
}

compute_mse_coverage <- function(mod, locations, func, level = 0.95){
  pred <- predict(mod, newdata = data.frame(x = locations), variable = "x", quantiles = c((1 - level)/2, level + (1 - level)/2), include.intercept = F)
  mse <- mean((func(locations) - pred$mean)^2)
  coverage <- mean((func(locations) >= pred[,2]) & (func(locations) <= pred[,3]))
  return(list(mse = mse, coverage = coverage))
}

fit_exact <- function(data, prior_SD){
  sorted_data <- data %>% arrange(x)

  customized_list <- list()
  customized_list$compute_B <- function(x){
    n <- length(x)
    B <- Matrix::Diagonal(n = 2*n)[,1:(2*n)]
    B <- as(B[seq(1,2*n,by = 2),][, -c(1:2)], "matrix")
    return(B)
  }
  customized_list$compute_P <- function(x){
    Q <- sGPfit::joint_prec_construct(a = pi, t_vec = x[-1], sd = 1)
    Q <- as(Q, "matrix")
    return(Q)
  }
  sorted_data$cos_x <- cos(pi * sorted_data$x); sorted_data$sin_x <- sin(pi * sorted_data$x)

  mod_exact <- model_fit(y ~ sin_x + cos_x +
                           f(x, model = "customized",
                            sd.prior = list(param = prior_SD)) +
                          f(index, model = "iid", sd.prior = 1),
                          data = sorted_data,
                          family = "Poisson",
                          Customized_RE = customized_list)

  return(mod_exact)
}


sample_exact <- function(mod_exact, M = 3000){
  samps_exact <- aghq::sample_marginal(mod_exact$mod, M = M)
  samps_exact_coef <- samps_exact$samps[mod_exact$random_samp_indexes$x, ]
  samps_exact_f <- mod_exact$instances[[1]]@B %*% samps_exact_coef
  samps_exact_global <- mod_exact$design_mat_fixed[[3]] %*% samps_exact$samps[mod_exact$fixed_samp_indexes$cos_x,,drop = F ] +
    mod_exact$design_mat_fixed[[2]] %*% samps_exact$samps[mod_exact$fixed_samp_indexes$sin_x,,drop = F ]
  samps_exact_overall <- samps_exact_global + samps_exact_f
  samps_exact_overall
}
