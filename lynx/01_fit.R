#### PATH:
working_path <- getwd()
source_path <- paste0(working_path,"/source/")
figure_path <- paste0(working_path,"/figures/")
result_path <- paste0(working_path,"/results/")

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
library(lubridate)
source(paste0(source_path, "05_functions_seasonal.R"))
source(paste0(source_path, "00_functions_predictGP.R"))

data <- data.frame(year = seq(1821, 1934, by = 1), logy = log10(as.numeric(lynx)), y = as.numeric(lynx))
data$x <- data$year - min(data$year)
x <- data$x
y <- data$y

### Define a prior on 50-years predictive SD:
pred_SD <- list(u = 1, alpha = 0.01)
compile(file = paste0(source_path, "01_lynx_two_sGP.cpp"))


### Put a discrete prior on a:
period_state <- seq(from = 6, to = 12, by = 0.1)
period <- 1/period_state
a_vec = 2*period*pi
log_ML_vec <- c()
for (i in 1:length(a_vec)) {
  period <- 0.1
  a1 = a_vec[i]
  scale1 <- compute_d_step_sGPsd(d = 50, a = a1)
  a2 = 2* a_vec[i]
  scale2 <- compute_d_step_sGPsd(d = 50, a = a2)
  n = length(x)
  B1 <- Diagonal(n = 2*n)[,1:(2*n)]
  B1 <- B1[seq(1,2*n,by = 2),][, -c(1:2)]
  Q1 <- joint_prec_construct(a = a1, t_vec = x[-1], sd = 1)
  Q1 <- as(as(Q1, "matrix"),"dgTMatrix")
  B1 <- as(as(B1, "matrix"),"dgTMatrix")
  B2 <- Diagonal(n = 2*n)[,1:(2*n)]
  B2 <- B2[seq(1,2*n,by = 2),][, -c(1:2)]
  Q2 <- joint_prec_construct(a = a2, t_vec = x[-1], sd = 1)
  Q2 <- as(as(Q2, "matrix"),"dgTMatrix")
  B2 <- as(as(B2, "matrix"),"dgTMatrix")
  B3 <- as(Matrix::diag(n), "dgTMatrix")
  Q3 = B3
  X <- as(cbind(cos(a1*x),sin(a1*x),cos(a2*x),sin(a2*x),1), "dgTMatrix")
  # X <- as(matrix(rep(1,n),ncol = 1), "dgTMatrix")
  dyn.load(dynlib(paste0(source_path,"01_lynx_two_sGP")))
  tmbdat <- list(
    # Design matrix
    B1 = B1,
    B2 = B2,
    B3 = B3,
    X = X,
    P1 = Q1,
    P2 = Q2,
    P3 = Q3,
    logP1det = as.numeric(determinant(Q1, logarithm = T)$modulus),
    logP2det = as.numeric(determinant(Q2, logarithm = T)$modulus),
    logP3det = as.numeric(determinant(Q3, logarithm = T)$modulus),
    # Response
    y = y,
    # PC Prior params
    u1 = pred_SD$u/scale1,
    alpha1 = pred_SD$alpha,
    u2 = pred_SD$u/scale2,
    alpha2 = pred_SD$alpha,
    u3 = 1,
    alpha3 = 0.01,
    betaprec = 0.001
  )
  tmbparams <- list(
    W = c(rep(0, (ncol(B1) + ncol(B2) + ncol(B3) + ncol(X)))), # W = c(U,beta); U = B-Spline coefficients
    theta1 = 0,
    theta2 = 0,
    theta3 = 0
  )
  ff <- TMB::MakeADFun(
    data = tmbdat,
    parameters = tmbparams,
    random = "W",
    DLL = "01_lynx_two_sGP",
    silent = TRUE
  )
  ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
  fitted_mod <- aghq::marginal_laplace_tmb(ff,5,c(0,0,0))
  hyper_result <- fitted_mod$marginals
  save(file = paste0(result_path, "samps/",i,"_hyper_sample.rda"), hyper_result)
  sGP_samps <- sample_marginal(fitted_mod, 3000)
  save(file = paste0(result_path, "samps/",i,"_sample.rda"), sGP_samps)
  log_ML_vec[i] <- fitted_mod$normalized_posterior$lognormconst
}
period_state[which.max(log_ML_vec)]

### Posterior of alpha:
log_ML_vec
ratio_ML <- exp(log_ML_vec - log_ML_vec[1])
Poster_alpha <- ratio_ML/sum(ratio_ML)
plot(Poster_alpha~a_vec, type = 'o')

pdf(file = paste0(figure_path, "lynx_poster_alpha.pdf"))
plot(Poster_alpha~period_state, type = 'p', xlim = c(8,12), xaxs='i',xlab = "periodicity (years)", ylab = "Posterior", ylim = c(0,0.28))
lines(predict(smooth.spline(x = period_state, y = Poster_alpha), x = seq(8,12,by = 0.01)), col = 'red')
abline(h = 1/length(period_state), col = "blue", lty = "dashed")
dev.off()

pos_alpha_data <- data.frame(alpha = a_vec, Prob = Poster_alpha)
save(file = paste0(result_path, "pos_alpha_data.rda"), pos_alpha_data)















