.libPaths( c( .libPaths(), "~/./lib") )

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
library(parallel)
library(foreach)
library(doMC)
library(sGPfit)

source(paste0(source_path, "06_functions_seasonal.R"))
compile(paste0(source_path,"04_Gaussian_Seasonal_CO2.cpp")) ## adding a new random component!
dyn.load(dynlib(paste0(source_path,"04_Gaussian_Seasonal_CO2")))


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
observed_dataset <- co2s %>% filter(!is.na(co2s$co2)) %>% dplyr::select(c("co2", "timeYears"))
observed_dataset$quality <- ifelse(co2s$quality > 0, 1, 0)
# remove low-quality measurements
observed_dataset <- observed_dataset %>% filter(quality == 0)



### Assume the d-year prediction SD has PC prior:
pred_SD_prior <- list(u = 1, alpha = 0.5)

### Implementing the seasonal sB basis:
d_step <- 10
period1 <- 1 ## 1-year
a1 = 2*period1*pi
scale1 <- compute_d_step_sGPsd(d = d_step,a = a1)
SD_prior1 <- list(u = pred_SD_prior$u/scale1, alpha = pred_SD_prior$alpha)

period2 <- 2 ## half-year
a2 = 2*period2*pi
scale2 <- compute_d_step_sGPsd(d = d_step,a = a2)
SD_prior2 <- list(u = pred_SD_prior$u/scale2, alpha = pred_SD_prior$alpha)

prior_set_IWP <- list(a = 0.5, u = 30)
prior_set_IWP_used <- prior_conversion(d = d_step, prior = prior_set_IWP, p = 3)


x <- observed_dataset$timeYears
region <- c(0,max(x))
k1 <- 50
k2 <- 30 # 30*3*5

B1 <- Compute_B_sB(x = x, a = a1, k = k2, region = region, boundary = T)
B1 <- as(B1,"dgTMatrix")
Q1 <- Compute_Q_sB(a = a1, k = k2, region = region, accuracy = 0.005)
Q1 <- as(as(Q1, "matrix"),"dgTMatrix")
X1 <- as(cbind(cos(a1*x), sin(a1*x)), "dgTMatrix")

B2 <- Compute_B_sB(x = x, a = a2, k = k2, region = region, boundary = T)
B2 <- as(B2,"dgTMatrix")
Q2 <- Compute_Q_sB(a = a2, k = k2, region = region, accuracy = 0.005)
Q2 <- as(as(Q2, "matrix"),"dgTMatrix")
X2 <- as(cbind(cos(a2*x), sin(a2*x)), "dgTMatrix")

Q6 <- as(compute_weights_precision(x = seq(0, max(x), length.out = k1)), "dgTMatrix")
B6 <- as(local_poly(knots = seq(0, max(x), length.out = k1), refined_x = x, p = 3), "dgTMatrix")

#### Posterior of alpha:
years_cyl3 <- 44/12
years_cyl4 <- 9.1
years_cyl5 <- 10.4

period3_vec <- 1/years_cyl3
period4_vec <- 1/years_cyl4
period5_vec <- 1/years_cyl5

period_vec <- expand.grid(period3_vec = period3_vec, period4_vec = period4_vec, period5_vec = period5_vec)
period_vec$years_cyl3 <- 1/period_vec$period3_vec
period_vec$years_cyl4 <- 1/period_vec$period4_vec
period_vec$years_cyl5 <- 1/period_vec$period5_vec

log_like_vec <- c()
do_once <- function(i){
  period3 <- period_vec$period3_vec[i] ## the third component, 44-month
  a3 = 2*period3*pi
  scale3 <- compute_d_step_sGPsd(d = d_step,a = a3)
  SD_prior3 <- list(u = pred_SD_prior$u/scale3, alpha = pred_SD_prior$alpha)
  B3 <- Compute_B_sB(x = x, a = a3, k = k2, region = region, boundary = T)
  B3 <- as(B3,"dgTMatrix")
  Q3 <- Compute_Q_sB(a = a3, k = k2, region = region, accuracy = 0.005)
  Q3 <- as(as(Q3, "matrix"),"dgTMatrix")
  X3 <- as(cbind(cos(a3*x), sin(a3*x)), "dgTMatrix")
  
  period4 <- period_vec$period4_vec[i] ## the fourth component, 9.1-year
  a4 = 2*period4*pi
  scale4 <- compute_d_step_sGPsd(d = d_step,a = a4)
  SD_prior4 <- list(u = pred_SD_prior$u/scale4, alpha = pred_SD_prior$alpha)
  B4 <- Compute_B_sB(x = x, a = a4, k = k2, region = region, boundary = T)
  B4 <- as(B4,"dgTMatrix")
  Q4 <- Compute_Q_sB(a = a4, k = k2, region = region, accuracy = 0.005)
  Q4 <- as(as(Q4, "matrix"),"dgTMatrix")
  X4 <- as(cbind(cos(a4*x), sin(a4*x)), "dgTMatrix")
  
  period5 <- period_vec$period5_vec[i] ## the fifth component, 10.4-year
  a5 = 2*period5*pi
  scale5 <- compute_d_step_sGPsd(d = d_step,a = a5)
  SD_prior5 <- list(u = pred_SD_prior$u/scale5, alpha = pred_SD_prior$alpha)
  B5 <- Compute_B_sB(x = x, a = a5, k = k2, region = region, boundary = T)
  B5 <- as(B5,"dgTMatrix")
  Q5 <- Compute_Q_sB(a = a5, k = k2, region = region, accuracy = 0.005)
  Q5 <- as(as(Q5, "matrix"),"dgTMatrix")
  X5 <- as(cbind(cos(a5*x), sin(a5*x)), "dgTMatrix")
  
  X <- as(cbind(X1,X2,X3,X4,X5, 1, x, x^2), "dgTMatrix")
  tmbdat <- list(
    # Design matrix
    B1 = B1,
    B2 = B2,
    B3 = B3,
    B4 = B4,
    B5 = B5,
    B6 = B6,
    X = X,
    P1 = Q1,
    P2 = Q2,
    P3 = Q3,
    P4 = Q4,
    P5 = Q5,
    P6 = Q6,
    logP1det = as.numeric(determinant(Q1, logarithm = T)$modulus),
    logP2det = as.numeric(determinant(Q2, logarithm = T)$modulus),
    logP3det = as.numeric(determinant(Q3, logarithm = T)$modulus),
    logP4det = as.numeric(determinant(Q4, logarithm = T)$modulus),
    logP5det = as.numeric(determinant(Q5, logarithm = T)$modulus),
    logP6det = as.numeric(determinant(Q6, logarithm = T)$modulus),
    # Response
    y = observed_dataset$co2,
    # PC Prior params
    u1 = SD_prior1$u,
    alpha1 = SD_prior1$alpha,
    u2 = SD_prior2$u,
    alpha2 = SD_prior2$alpha,
    u3 = SD_prior3$u,
    alpha3 = SD_prior3$alpha,
    u4 = SD_prior4$u,
    alpha4 = SD_prior4$alpha,
    u5 = SD_prior5$u,
    alpha5 = SD_prior5$alpha,
    u6 = prior_set_IWP_used$u,
    alpha6 = prior_set_IWP_used$a,
    u7 = 1,
    alpha7 = 0.5,
    betaprec = 0.001
  )
  tmbparams <- list(
    W = c(rep(0, (ncol(B1) + ncol(B2) + ncol(B3) + ncol(B4) + ncol(B5) + ncol(B6) + ncol(X)))), # W = c(U,beta); U = B-Spline coefficients
    theta1 = 0, # -2log(sigma)
    theta2 = 0,
    theta3 = 0,
    theta4 = 0,
    theta5 = 0,
    theta6 = 0,
    theta7 = 0
  )
  ff <- TMB::MakeADFun(
    data = tmbdat,
    parameters = tmbparams,
    random = "W",
    DLL = "04_Gaussian_Seasonal_CO2",
    silent = TRUE
  )
  ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
  # ff$fn(c(0,0,0,0,0,0,0)) # sanity check; if vmin = inf
  ## Try two choices of aghq_k, in case numerical failure..
  fitted_mod_all <- tryCatch(
    {aghq::marginal_laplace_tmb(ff,4,c(0,0,0,0,0,0,0))}, 
    error=function(e) {
      message(paste('An Error Occurred at task number', i, "...", "Trying a lower aghq_k..."))
      aghq::marginal_laplace_tmb(ff,3,c(0,0,0,0,0,0,0))
    }
  ) 
  hyper_result <- fitted_mod_all$marginals
  save(file = paste0(result_path, "samps/",i,"_hyper_sample.rda"), hyper_result)
  sGP_samps <- sample_marginal(fitted_mod_all, 3000)
  save(file = paste0(result_path, "samps/",i,"_sample.rda"), sGP_samps)
  cat(paste0("task ", i," just finished with success: \n"))
  cat(paste0("task ", i," has lognormconst being ",fitted_mod_all$normalized_posterior$lognormconst, ".\n"))
  fitted_mod_all$normalized_posterior$lognormconst
}
log_like_vec <- foreach(i = 1:nrow(period_vec), .combine='c', .packages = c('foreach', 'stats', 'OSplines', 'fda', 'aghq', 'LaplacesDemon')) %do% do_once(i)
ratio_ML <- exp(log_like_vec - log_like_vec[1])
Poster_alpha <- ratio_ML/sum(ratio_ML)


### Save the output:
pos_alpha_data_para <- cbind(data.frame(Prob = Poster_alpha), period_vec)
save(file = paste0(working_path,"/results/pos_alpha_data_para.rda"), pos_alpha_data_para)





