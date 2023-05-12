library(OSplines)
library(Matrix)
library(fda)
library(aghq)
library(tidyverse)
source("03_functions_seasonal.R")
### Comparing prediction method:


### Preparation: Read data and fit models
deadFile = Pmisc::downloadIfOld("https://www150.statcan.gc.ca/n1/tbl/csv/13100810-eng.zip",
                                path = "D:/seasonality/seasonal-spline/v5/influ")
(deadFileCsv = deadFile[which.max(file.info(deadFile)$size)])
x = read.csv(deadFileCsv)
x$date = as.Date(as.character(x[[grep("DATE", names(x))]]))
x$province = gsub("[,].*", "", x$GEO)
x = x[x$province ==
        "Ontario", ]
par(mfrow = c(2,2))
for (D in c("heart", "neoplasms", "respiratory", "Influenza")) {
  plot(x[grep(D, x$Cause), c("date", "VALUE")], ylab = D)
  abline(v = as.Date("2020/03/17"))
}
dateSeq = sort(unique(x$date))
dateSeqInt = as.integer(dateSeq)
x$dateInt = x$dateIid = as.integer(x$date)
x$cos12 = cos(2 * pi * x$dateInt/365.25)
x$cos6 = cos(2 * pi * x$dateInt * 2/365.25)
x$sin12 = sin(2 * pi * x$dateInt/365.25)
x$sin6 = sin(2 * pi * x$dateInt * 2/365.25)
xInfluenza = x[grepl("Influenza", x$Cause) & x$province == "Ontario", ]
xPreCovid = xInfluenza[xInfluenza$date < as.Date("2020/03/01"),]
par(mfrow = c(1,1))
# compile("00_Gaussian_Seasonal_Influ_reduced.cpp")
dyn.load(dynlib("00_Gaussian_Seasonal_Influ_reduced"))
t <- xPreCovid$dateInt/365.25 - min(xPreCovid$dateInt/365.25)
upper_region <- max(x$dateInt/365.25 - min(xPreCovid$dateInt/365.25))
period1 <- 2
a1 = period1*pi
region <- c(0,upper_region)
k1 <- 100
B1 <- Compute_B_sB(x = t, a = a1, k = k1, region = region, boundary = T)
B1 <- as(B1,"dgTMatrix")
Q1 <- Compute_Q_sB(a = a1, k = k1, region = region, accuracy = 0.001)
Q1 <- as(as(Q1, "matrix"),"dgTMatrix")
X1 <- as(cbind(cos(a1*t), sin(a1*t)), "dgTMatrix")
X <- as(cbind(cos(a1*t), sin(a1*t), 1, t), "dgTMatrix")
tmbdat <- list(
  # Design matrix
  B1 = B1,
  X = X,
  P1 = Q1,
  logP1det = as.numeric(determinant(Q1, logarithm = T)$modulus),
  # Response
  y = xPreCovid$VALUE,
  # PC Prior params
  u1 = 1,
  alpha1 = 0.5,
  betaprec = 0.001
)
tmbparams <- list(
  W = c(rep(0, (ncol(B1) + ncol(X)))), # W = c(U,beta); U = B-Spline coefficients
  theta1 = 0
)
ff <- TMB::MakeADFun(
  data = tmbdat,
  parameters = tmbparams,
  random = "W",
  DLL = "00_Gaussian_Seasonal_Influ_reduced",
  silent = TRUE
)

ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
fitted_FEM_mod <- aghq::marginal_laplace_tmb(ff,10,c(0))



### 1: Fit FEM approximation, use function extrapolation:
xPostCovid = xInfluenza
t_pred <- xPostCovid$dateInt/365.25 - min(xPreCovid$dateInt/365.25)
B_pred <- Compute_B_sB(x = t_pred, a = a1, k = k1, region = region, boundary = T)
X_pred <- as(cbind(cos(a1*t_pred), sin(a1*t_pred), 1, t_pred), "dgTMatrix")
samps <- sample_marginal(fitted_FEM_mod, 3000)
B <- as.matrix(B1)
samps_w <- samps$samps[1:ncol(B),]
samps_beta <- samps$samps[(ncol(B) + 1):(ncol(B) + 4),]
samps_g_pred <- as.matrix(B_pred %*% samps_w) + as.matrix(X_pred %*% samps_beta)
samps_g_pred <- as.data.frame(samps_g_pred)
samps_g_pred$t <- t_pred + 2010
samps_g_pred <- arrange(samps_g_pred, by = t)
matplot(y = exp(samps_g_pred[,1:30]), x = samps_g_pred$t, type = 'l', 
        add = F, lty = 'dashed', col = 'pink',
        ylim = c(0,200))
samps_g_pred$mean <- apply(exp(samps_g_pred)[,1:3000],1, mean)
samps_g_pred$upper <- apply(exp(samps_g_pred)[,1:3000],1, quantile, 0.975)
samps_g_pred$lower <- apply(exp(samps_g_pred)[,1:3000],1, quantile, 0.025)
plot(mean~t, data = samps_g_pred, type = 'l')
lines(upper~t, data = samps_g_pred, lty = 'dashed')
lines(lower~t, data = samps_g_pred, lty = 'dashed')
points(xInfluenza$VALUE~I(sort(t_pred) + 2010))



### 2: Fit exact method, use augmented space method
source("00_sGP_exact.R")
Q_exact <- Compute_precision_Aug_SGP(t_pred[-1], a = a1)
B_exact <- rbind(0,Compute_design_Aug_sGP(t_pred[-1])[1:(nrow(xPreCovid)-1), ])
tmbdat <- list(
  # Design matrix
  B1 = B_exact,
  X = X,
  P1 = Q_exact,
  logP1det = as.numeric(determinant(Q_exact, logarithm = T)$modulus),
  # Response
  y = xPreCovid$VALUE,
  # PC Prior params
  u1 = 1,
  alpha1 = 0.5,
  betaprec = 0.001
)
tmbparams <- list(
  W = c(rep(0, (ncol(B_exact) + ncol(X)))), # W = c(U,beta); U = B-Spline coefficients
  theta1 = 0
)
ff <- TMB::MakeADFun(
  data = tmbdat,
  parameters = tmbparams,
  random = "W",
  DLL = "00_Gaussian_Seasonal_Influ_reduced",
  silent = TRUE
)
ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
fitted_Exact_mod <- aghq::marginal_laplace_tmb(ff,10,c(0))
samps <- sample_marginal(fitted_Exact_mod, 3000)
samps_w <- samps$samps[1:ncol(B_exact),]
samps_beta <- samps$samps[(ncol(B_exact) + 1):(ncol(B_exact) + 4),]
samps_fw <- samps_w[seq(1,nrow(samps_w), by = 2), ]
samps_f1stw <- samps_w[seq(2,nrow(samps_w), by = 2), ]

X_full <- as(cbind(cos(a1*t_pred), sin(a1*t_pred), 1, t_pred), "dgTMatrix")
samps_fw_pred <- rbind(0,samps_fw) + as.matrix(X_full %*% samps_beta)
samps_fw_pred_result <- data.frame(t = t_pred + 2010)
samps_fw_pred_result$mean <- apply(exp(samps_fw_pred)[,1:3000],1, mean)
samps_fw_pred_result$upper <- apply(exp(samps_fw_pred)[,1:3000],1, quantile, 0.975)
samps_fw_pred_result$lower <- apply(exp(samps_fw_pred)[,1:3000],1, quantile, 0.025)
samps_fw_pred_result <- arrange(samps_fw_pred_result, by = t)

plot(mean~t, data = samps_fw_pred_result, type = 'l')
lines(upper~t, data = samps_fw_pred_result, lty = 'dashed')
lines(lower~t, data = samps_fw_pred_result, lty = 'dashed')
points(xInfluenza$VALUE~I(sort(t_pred) + 2010))


X_deriv <- as(cbind(-a1*sin(a1*t_pred), a1*cos(a1*t_pred), 1), "dgTMatrix")
samps_f1stw_pred <- rbind(0,samps_f1stw) + as.matrix(X_deriv %*% samps_beta[-3,])
samps_f1stw_pred_result <- data.frame(t = t_pred + 2010)
samps_f1stw_pred_result$mean <- apply((samps_f1stw_pred)[,1:3000],1, mean)
samps_f1stw_pred_result$upper <- apply((samps_f1stw_pred)[,1:3000],1, quantile, 0.975)
samps_f1stw_pred_result$lower <- apply((samps_f1stw_pred)[,1:3000],1, quantile, 0.025)
samps_f1stw_pred_result <- arrange(samps_f1stw_pred_result, by = t)

plot(mean~t, data = samps_f1stw_pred_result, type = 'l')
lines(upper~t, data = samps_f1stw_pred_result, lty = 'dashed')
lines(lower~t, data = samps_f1stw_pred_result, lty = 'dashed')


















