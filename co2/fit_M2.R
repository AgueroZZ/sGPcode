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
compile(paste0(source_path, "00_fixed_Seasonal_CO2.cpp"))
dyn.load(dynlib(paste0(source_path,"00_fixed_Seasonal_CO2")))


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
x_grid <- observed_dataset$co2
n <- length(x_grid)

years_cyl1 <- 1
years_cyl2 <- 0.5
years_cyl3 <- 44/12
years_cyl4 <- 9.1
years_cyl5 <- 10.4
period1 <- 1/years_cyl1
period2 <- 1/years_cyl2
period3 <- 1/years_cyl3
period4 <- 1/years_cyl4
period5 <- 1/years_cyl5
a1 = 2*period1*pi
a2 = 2*period2*pi
a3 = 2*period3*pi
a4 = 2*period4*pi
a5 = 2*period5*pi

d_step <- 10
prior_set_IWP <- list(a = 0.5, u = 30)
prior_set_IWP_used <- prior_conversion(d = d_step, prior = prior_set_IWP, p = 3)

x <- observed_dataset$timeYears
region <- c(0,max(x))
k <- 102
X1 <- as(cbind(cos(a1*x), sin(a1*x)), "dgTMatrix")
X2 <- as(cbind(cos(a2*x), sin(a2*x)), "dgTMatrix")
X3 <- as(cbind(cos(a3*x), sin(a3*x)), "dgTMatrix")
X4 <- as(cbind(cos(a4*x), sin(a4*x)), "dgTMatrix")
X5 <- as(cbind(cos(a5*x), sin(a5*x)), "dgTMatrix")


Q6 <- as(compute_weights_precision(x = seq(0, max(x), length.out = k)), "dgTMatrix")
B6 <- as(local_poly(knots = seq(0, max(x), length.out = k), refined_x = x, p = 3), "dgTMatrix")

#### Posterior of alpha:
### With parallel computation:
ncores <- detectCores() - 1
registerDoMC(cores = ncores)

period3_vec <- 1/years_cyl3
period4_vec <- 1/years_cyl4
period5_vec <- 1/years_cyl5

period_vec <- expand.grid(period3_vec = period3_vec, period4_vec = period4_vec, period5_vec = period5_vec)
period_vec$years_cyl3 <- 1/period_vec$period3_vec
period_vec$years_cyl4 <- 1/period_vec$period4_vec
period_vec$years_cyl5 <- 1/period_vec$period5_vec

log_like_vec <- c()
do_once <- function(i){
  X <- as(cbind(X1,X2,X3,X4,X5, 1, x, x^2), "dgTMatrix")
  tmbdat <- list(
    # Design matrix
    B = B6,
    X = X,
    P = Q6,
    logPdet = as.numeric(determinant(Q6, logarithm = T)$modulus),
    # Response
    y = observed_dataset$co2,
    # PC Prior params
    u1 = prior_set_IWP_used$u,
    alpha1 = prior_set_IWP_used$a,
    u2 = 1,
    alpha2 = 0.5,
    betaprec = 0.001
  )
  tmbparams <- list(
    W = c(rep(0, (ncol(B6) + ncol(X)))), # W = c(U,beta); U = B-Spline coefficients
    theta1 = 0, # -2log(sigma)
    theta2 = 0
  )
  ff <- TMB::MakeADFun(
    data = tmbdat,
    parameters = tmbparams,
    random = "W",
    DLL = "00_fixed_Seasonal_CO2",
    silent = TRUE
  )
  ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
  
  ## Try two choices of aghq_k, in case numerical failure..
  fitted_mod_all <- tryCatch(
    {aghq::marginal_laplace_tmb(ff,4,c(0,0))}, 
    error=function(e) {
      message(paste('An Error Occurred at task number', i, "...", "Trying a lower aghq_k..."))
      aghq::marginal_laplace_tmb(ff,3,c(0,0))
    }
  ) 
  hyper_result <- fitted_mod_all$marginals
  save(file = paste0(result_path, "fixed_samps/",i,"_hyper_sample.rda"), hyper_result)
  sGP_samps <- sample_marginal(fitted_mod_all, 3000)
  save(file = paste0(result_path, "fixed_samps/",i,"_sample.rda"), sGP_samps)
  cat(paste0("task ", i," just finished with success: \n"))
  cat(paste0("task ", i," has lognormconst being ",fitted_mod_all$normalized_posterior$lognormconst, ".\n"))
  fitted_mod_all$normalized_posterior$lognormconst
}
log_like_vec <- foreach(i = 1:nrow(period_vec), .combine='c', .packages = c('foreach', 'stats', 'OSplines', 'fda', 'aghq', 'LaplacesDemon')) %do% do_once(i)




