### Making prediction for sGP:
predict_sGP_Var <- function(observed_f_vec, a, sd, t = NULL, mesh_size = 0.01, max_t = 10, last_t){
  if(is.null(t)){
    t <- seq(last_t, max_t, by = mesh_size)
  }
  ### Construct transition matrix:
  M_construct <- function(t0, t1, a){
    M <- matrix(nrow = 2, ncol = 2)
    d <- t1 - t0
    M[1,1] <- cos(a*d)
    M[1,2] <- (1/a)*sin(a*d)
    M[2,1] <- (-a)*sin(a*d)
    M[2,2] <- cos(a*d)
    M
  }
  ### Construct Sig matrix:
  Sig_construct <- function(t0, t1, a, sd){
    Sig <- matrix(nrow = 2, ncol = 2)
    d <- t1 - t0
    Sig[1,1] <- (d/(2*(a^2))) - sin(2*a*d)/(4*(a^3))
    Sig[1,2] <- (sin(a*d)^2)/(2*(a^2))
    Sig[2,1] <- Sig[1,2]
    Sig[2,2] <- (d/2) + (sin(2*a*d)/(4*a))
    (sd^2)*Sig
  }
  n <- length(t) - 1
  matT <- M_construct(t0 = last_t, t1 = t[2], a = a)
  Sig <- Sig_construct(t0 = last_t, t1 = t[2], a = a, sd = sd)
  sample_path <- tsDyn::VAR.sim(B = matT, lag = 1, n = n, starting = matrix(c(observed_f_vec),nrow = 1), varcov = Sig, include = "none", returnStarting = T)
  result <- data.frame(t = t, sample_path)
  names(result)[2] <- "function"
  names(result)[3] <- "1stDeriv"
  result
}

# a = 2*pi; sd = 1
# my_samps <- sim_sGP_Var(a = a, sd = sd)
# predict_samps <- predict_sGP_Var(a = a, sd = sd, observed_f_vec = as.numeric(my_samps[nrow(my_samps),2:3]), last_t = 10, max_t = 15)
# plot(my_samps[,2]~my_samps$t, type = 'l', xlim = c(0,15), ylim = c(-1,1))
# lines(predict_samps[,2]~predict_samps$t, col = 'red')



### Making prediction for IWP:
Compute_Ti <- function(svec,p = 2,i){
  Ti <- matrix(0,nrow = p, ncol = p)
  delta <- diff(c(0,svec))
  denom <- factorial(c(0:(p-1)))
  numr <- delta[i+1]^(0:(p-1))
  Ti[1,] <- numr/denom
  for (i in 2:p) {
    Ti[i,] <- c(rep(0,(i-1)),Ti[(i-1),((i-1):(p-1))])
  }
  Ti
}
Compute_Ci <- function(svec, p = 2, i, is.cov = FALSE){
  delta <- diff(c(0,svec))
  Result <- matrix(0,nrow = p, ncol = p)
  index <- i+1
  for (i in 1:p) {
    for (j in i:p) {
      Result[i,j] <- (delta[index]^(2*p + 1 - i - j))/((2*p + 1 - i - j)*factorial(p-i)*factorial(p-j))
    }
  }
  Result <- Matrix::forceSymmetric(Result)
  if(is.cov == T){
    return(Result)
  }
  else{
    round(solve(Result),digits = 5)
  }
}
predict_IWP_Var <- function(observed_f_vec, p, sd, t = NULL, mesh_size = 0.01, max_t = 10, last_t){
  if(is.null(t)){
    t <- seq(last_t, max_t, by = mesh_size)
  }
  n <- length(t) - 1
  matT <- Compute_Ti(svec = t, p = p, i = 1)
  Sig <- (sd^2) * Compute_Ci(svec = t, p = p, i = 1, is.cov = T)
  sample_path <- tsDyn::VAR.sim(B = matT, lag = 1, n = n, starting = matrix(observed_f_vec,nrow = 1), varcov = Sig, include = "none", returnStarting = T)
  result <- data.frame(t = t, sample_path)
  names(result)[2] <- "function"
  for (i in 2:ncol(result)) {
    names(result)[i] <- paste0("IW",(p-i+2))
  }
  result
}
# my_samps <- predict_IWP_Var(p = 3, sd = sd, observed_f_vec = as.numeric(c(0,0,0)), last_t = 0, max_t = 10)
# predict_samps <- predict_IWP_Var(p = 3, sd = sd, observed_f_vec = as.numeric(my_samps[nrow(my_samps),2:4]), last_t = 10, max_t = 15)
# plot(my_samps[,2]~my_samps$t, type = 'l', xlim = c(0,15), ylim = c(-200,200))
# lines(predict_samps[,2]~predict_samps$t, col = 'red')

# sd = 0.08
# predict_samps <- predict_IWP_Var(p = 3, sd = sd, observed_f_vec = as.numeric(my_samps[nrow(my_samps),2:4]), last_t = 10, max_t = 15)[,1:2]
# for (i in 1:1000) {
#   predict_one_samps <- predict_IWP_Var(p = 3, sd = sd, observed_f_vec = as.numeric(my_samps[nrow(my_samps),2:4]), last_t = 10, max_t = 15)
#   predict_samps <- cbind(predict_samps, predict_one_samps[,2])
# }
# plot(my_samps[,2]~my_samps$t, type = 'l', xlim = c(0,15), ylim = c(-200,200))
# matplot(y = predict_samps[,2:1000], x = predict_samps[,1], add = T, type = 'l')



### Compute derivatives based on function samples:
compute_numeric_deriv <- function(gx, h, degree = 1){
  diff(gx, differences = degree)/(h^degree)
}
extract_deriv_samples <- function(fitted_samps, t, degree){
  h <- mean(diff(t))
  fitted_samps_deriv <- apply(fitted_samps,2, compute_numeric_deriv, h = h, degree = degree)
  result <- cbind(t = t[-c(1:degree)], data.frame(as.matrix(fitted_samps_deriv)))
  result
}
# deriv_samples <- extract_deriv_samples(fitted_samps = predict_samps[,-1], t = predict_samps[,1], degree = 1)
# plot(deriv_samples$X1002 ~ deriv_samples$t, type = 'l')
# lines(predict_one_samps[,3]~predict_one_samps[,1], type = 'l', col = "red")

# deriv_samples <- extract_deriv_samples(fitted_samps = predict_samps[,-1], t = predict_samps[,1], degree = 2)
# plot(deriv_samples$X1002 ~ deriv_samples$t, type = 'l')
# lines(predict_one_samps[,4]~predict_one_samps[,1], type = 'l', col = "red")


