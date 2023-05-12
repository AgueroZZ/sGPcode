library(fda)
library(mvtnorm)
library(LaplacesDemon)


## Functions for seasonal smoothings:

## Construct weights precision with sB spline basis:
# a: number of period, cos(ax) + sin(ax)
# k: number of knots
# region: the region of interest
# accuracy: the grid size of numerical integration
# boundary: whether the boundary condition is applied
# Q <- Compute_Q_sB(a = 2, k = 10, region = c(0,10), boundary = T)
Compute_Q_sB <- function(a,k,region, accuracy = 0.01, boundary = TRUE){
  ss <- function(M) {forceSymmetric(M + t(M))}
  x <- seq(min(region),max(region),by = accuracy)
  if(boundary == TRUE){
    B_basis <- suppressWarnings(create.bspline.basis(rangeval = c(min(region),max(region)),
                                    nbasis = k,
                                    norder = 4,
                                    dropind = c(1,2)))
  }
  else{
    B_basis <- suppressWarnings(create.bspline.basis(rangeval = c(min(region),max(region)),
                                    nbasis = k,
                                    norder = 4))
  }
  Bmatrix <- eval.basis(x, B_basis, Lfdobj=0, returnMatrix=T)
  B1matrix <-  eval.basis(x, B_basis, Lfdobj=1, returnMatrix=T)
  B2matrix <-  eval.basis(x, B_basis, Lfdobj=2, returnMatrix=T)
  cos_matrix <- cos(a*x)
  sin_matrix <- sin(a*x)
  Bcos <- as(apply(Bmatrix, 2, function(x) x*cos_matrix), "dgCMatrix")
  B1cos <- as(apply(B1matrix, 2, function(x) x*cos_matrix), "dgCMatrix")
  B2cos <- as(apply(B2matrix, 2, function(x) x*cos_matrix), "dgCMatrix")
  Bsin <- as(apply(Bmatrix, 2, function(x) x*sin_matrix), "dgCMatrix")
  B1sin <- as(apply(B1matrix, 2, function(x) x*sin_matrix), "dgCMatrix")
  B2sin <- as(apply(B2matrix, 2, function(x) x*sin_matrix), "dgCMatrix")
  
  ### Compute I, L, T:
  Numerical_I <- as(diag(c(diff(c(0,x)))), "dgCMatrix")
  
  ### T
  T00 <- t(Bcos) %*% Numerical_I %*% Bcos
  T10 <- t(B1cos) %*% Numerical_I %*% Bcos
  T11 <- t(B1cos) %*% Numerical_I %*% B1cos
  T20 <- t(B2cos) %*% Numerical_I %*% Bcos
  T21 <- t(B2cos) %*% Numerical_I %*% B1cos
  T22 <- t(B2cos) %*% Numerical_I %*% B2cos
  
  ### L
  L00 <- t(Bsin) %*% Numerical_I %*% Bsin
  L10 <- t(B1sin) %*% Numerical_I %*% Bsin
  L11 <- t(B1sin) %*% Numerical_I %*% B1sin
  L20 <- t(B2sin) %*% Numerical_I %*% Bsin
  L21 <- t(B2sin) %*% Numerical_I %*% B1sin
  L22 <- t(B2sin) %*% Numerical_I %*% B2sin
  
  ### I
  I00 <- t(Bsin) %*% Numerical_I %*% Bcos
  I10 <- t(B1sin) %*% Numerical_I %*% Bcos
  I11 <- t(B1sin) %*% Numerical_I %*% B1cos
  I20 <- t(B2sin) %*% Numerical_I %*% Bcos
  I21 <- t(B2sin) %*% Numerical_I %*% B1cos
  I22 <- t(B2sin) %*% Numerical_I %*% B2cos
  
  ## G = <phi,phj>
  G <- rbind(cbind(T00, t(I00)), cbind(I00,L00))
  
  ## C = <D^2phi,D^2phj>
  C11 <- T22 - 2*a*ss(I21) - (a^2)*ss(T20) + 2*(a^3)*ss(I10) + 4 * (a^2) * L11 + (a^4)*T00
  C22 <- L22 + 2*a*ss(I21) - (a^2)*ss(L20) - 2*(a^3)*ss(I10) + 4 * (a^2) * T11 + (a^4)*L00
  C12 <- I22 + 2*a*T21 - (a^2)* ss(I20) - 2*a*t(L21) - 4*(a^2)*I11 + 2*(a^3)*L10 - 2*(a^3)*t(T10) + (a^4)*I00
  C <- rbind(cbind(C11,C12), cbind(t(C12), C22))
  
  ## M = <phi,D^2phj>
  M11 <- t(T20) - (2*a)*t(I10) - (a^2)*T00
  M12 <- t(I20) + (2*a)*t(T10) - (a^2)*I00
  M21 <- t(I20) - (2*a)*t(L10) - (a^2)*I00
  M22 <- t(L20) + (2*a)*t(I10) - (a^2)*L00
  M <- rbind(cbind(M11,M12), cbind(M21, M22))
  
  
  ### Compute the final precision matrix: Q
  Q <- (a^4)*G + C + (a^2)*ss(M)
  forceSymmetric(Q)
}




## Construct Design matrix with sB spline basis:
Compute_B_sB <- function(x, a, k, region, boundary = T){
  if(boundary){
    B_basis <- suppressWarnings(create.bspline.basis(rangeval = c(min(region),max(region)),
                                    nbasis = k,
                                    norder = 4,
                                    dropind = c(1,2)))
  }
  else{
    B_basis <- suppressWarnings(create.bspline.basis(rangeval = c(min(region),max(region)),
                                    nbasis = k,
                                    norder = 4))
  }
  Bmatrix <- eval.basis(x, B_basis, Lfdobj=0, returnMatrix=T)
  cos_matrix <- cos(a*x)
  sin_matrix <- sin(a*x)
  Bcos <- apply(Bmatrix, 2, function(x) x*cos_matrix)
  Bsin <- apply(Bmatrix, 2, function(x) x*sin_matrix)
  cbind(Bcos, Bsin)
}
# B <- Compute_B_sB(x = x, a = 2, k = 10, region = c(0,10), boundary = T)


## Compute marginal variance of sB spline approximation:
compute_marginal_var <- function(region, k, a, Prec, boundary = T){
  rng <- region
  period <- a
  nbasis <- k
  norder <- 4
  Cov <- solve(Prec)
  basis <- create.bspline.basis(rng, nbasis, norder)
  design <- Compute_B_sB(x = basis$params, a = a, region = rng, boundary = boundary, k = k)
  Cov_f <- design %*% Cov %*% t(design)
  ### Compute the average:
  exp(mean(log(diag(Cov_f))))
}
# compute_marginal_var(a = 2, k = 10, region = c(0,10), boundary = T, Prec = Q)

## Sampling sample path from the prior distribution
sampling_from_weights <- function(x, a, k, region, boundary = T, Prec, n = 1){
  B <- Compute_B_sB(x, a, k, region, boundary = boundary)
  coefs_samps <- rmvnp(n = n, mu = rep(0,(2*(k-2))), Omega = as.matrix(Prec))
  splfd <- B %*% t(coefs_samps)
  (splfd)
}
# samps <- sampling_from_weights(x, a = 1, k = 10, region = c(0,10), boundary = T, Prec = Q, n = 2)


## Computing the correponding sigmas value:
compute_sigma <- function(a, k, region, boundary = T, reference = 1){
  Prec <- Compute_Q_sB(a = a, k = k, region = region, boundary = boundary)
  marg_var <- compute_marginal_var(region = region, k = k, a = a, Prec = Prec, boundary = boundary)
  sigma <- sqrt(reference/marg_var)
  sigma
}
# compute_sigma <- compute_sigma(a = 1, k = 10, region = c(0,10), boundary = T)



## Construct a covariance matrix given a covariance function:
compute_matrix_given_cov <- function(from, to, m, K){
  t <- seq(from = from, to = to, length.out = m)
  Sigma <- sapply(t, function(s1) {
    sapply(t, function(s2) {
      K(s1, s2)
    })
  })
  Sigma
}


## Simulate (any) gaussian process given (any) SPD covariance function 
gaussprocess <- function(from = 0, to = 1, K = function(s, t) {min(s, t)},
                         start = 0, m = 1000) {
  # Simulates a Gaussian process with a given kernel
  #
  # args:
  #   from: numeric for the starting location of the sequence
  #   to: numeric for the ending location of the sequence
  #   K: a function that corresponds to the kernel (covariance function) of
  #      the process; must give numeric outputs, and if this won't produce a
  #      positive semi-definite matrix, it could fail; default is a Wiener
  #      process
  #   start: numeric for the starting position of the process
  #   m: positive integer for the number of points in the process to simulate
  #
  # return:
  #   A data.frame with variables "t" for the time index and "xt" for the value
  #   of the process
  
  t <- seq(from = from, to = to, length.out = m)
  Sigma <- sapply(t, function(s1) {
    sapply(t, function(s2) {
      K(s1, s2)
    })
  })
  
  path <- mvrnorm(mu = rep(0, times = m), Sigma = Sigma)
  path <- path - path[1] + start  # Must always start at "start"
  
  return(data.frame("t" = t, "xt" = path))
}


### Generate true covariance function:
generate_K_true <- function(sigma,alpha){
  K_true <- function(s,t){
    ((sigma/alpha)^2)*((min(s,t)/2)*cos(alpha*(abs(s-t))) - (cos(alpha*max(c(s,t)))*sin(alpha*min(c(s,t))))/(2*alpha))
  }
  K_true
}



## Measure relative error of an approximation:
compute_rel_Error <- function(approx, truth, method = "F"){
  norm((truth - approx), type = method)/norm(truth, type = method)
}


## State Space Representation:
joint_prec_construct <- function(t_vec, a, sd){
  n <- length(t_vec)
  Blist <- list()
  AClist <- list()
  Clist <- list()
  
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
  Compute_Ci <- function(t_vec, i){
    Ci <- forceSymmetric(solve(Sig_construct(t0 = t_vec[i], t1 = t_vec[i+1], a, sd)))
    Ci
  }
  
  Compute_Ai <- function(t_vec, i, Ci) {
    Ti <- M_construct(t0 = t_vec[i], t1 = t_vec[i+1], a)
    t(Ti) %*% Ci %*% Ti
  }
  
  Compute_Bi <- function(t_vec, i, Ci) {
    Ti <- M_construct(t0 = t_vec[i], t1 = t_vec[i+1], a)
    -t(Ti) %*% Ci
  }
  
  for (i in 1:(n - 1)) {
    Clist[[i]] <- Compute_Ci(t_vec = t_vec, i = i)
  }
  
  for (i in 2:(n - 1)) {
    AClist[[i]] <- Compute_Ai(t_vec = t_vec, i = i, Ci = Clist[[i]]) + Clist[[i-1]]
  }
  
  AClist[[1]] <- Compute_Ai(t_vec = t_vec, i = 1, Ci = Clist[[1]]) + Compute_Ci(t_vec = c(0,t_vec[1]), i = 1) 
  AClist[[n]] <- Compute_Ci(t_vec = t_vec, i = (n - 1))
  
  for (i in 1:(n - 1)) {
    Blist[[i]] <- Compute_Bi(t_vec = t_vec, i = i, Ci = Clist[[i]])
  }
  
  Qlist <- list()
  Q <- matrix(0, nrow = 0, ncol = n*2)
  for (i in 1:(n-1)) {
    Qlist[[i]] <- cbind(matrix(0,nrow = 2, ncol = 2 * (i-1)), AClist[[i]], Blist[[i]], matrix(0,nrow = 2, ncol = (2 * (n-i-1))) )
    Q <- rbind(Q,Qlist[[i]])
  }
  Q <- rbind(Q,cbind(matrix(0,nrow = 2, ncol = 2 * (n-1)), AClist[[n]]))
  Q <- Matrix::forceSymmetric(Q)
  as(as.matrix(Q), "dgTMatrix")
}



### Simulate sGP as VAR(1) model
sim_sGP_Var <- function(t = NULL, mesh_size = 0.01, max_t = 10, a, sd){
  if(is.null(t)){
    t <- seq(0, max_t, by = mesh_size)
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
  matT <- M_construct(t0 = t[1], t1 = t[2], a = a)
  Sig <- Sig_construct(t0 = t[1], t1 = t[2], a = a, sd = sd)
  sample_path <- tsDyn::VAR.sim(B = matT, lag = 1, n = n, starting = NULL, varcov = Sig, include = "none", returnStarting = T)
  result <- data.frame(t = t, sample_path)
  names(result)[2] <- "function"
  names(result)[3] <- "1stDeriv"
  result
}



### Scaling coefficient:
compute_d_step_sGPsd <- function(d,a){
  sqrt((1/(a^2))*((d/2) - (sin(2*a*d)/(4*a))))
}

