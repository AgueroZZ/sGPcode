.libPaths( c( .libPaths(), "~/./lib") )
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
library(R.utils)
library(OSplines)


#### PATH:
set.seed(123)
working_path <- getwd()
source_path <- paste0(working_path,"/source/")
figure_path <- paste0(working_path,"/figures/")
result_path <- paste0(working_path,"/results/")

### Functions and CPP files:
source(paste0(source_path,"05_functions_seasonal.R"))
source(paste0(source_path,"00_sGP_exact.R"))

TMB::compile(file = paste0(source_path, "01_Gaussian_Seasonal_Influ_reduced.cpp"))
dyn.load(dynlib(paste0(source_path, "01_Gaussian_Seasonal_Influ_reduced")))

#### Data:
deadFile = Pmisc::downloadIfOld("https://www150.statcan.gc.ca/n1/tbl/csv/13100810-eng.zip",
                                path = paste0(source_path))
(deadFileCsv = deadFile[which.max(file.info(deadFile)$size)])
x = read.csv(deadFileCsv)

x$date = as.Date(as.character(x[[grep("DATE", names(x))]]))
x$province = gsub("[,].*", "", x$GEO)
x = x[x$province ==
        "Ontario", ]
x = x[x$date <= "2022-01-01", ]
for (D in c("heart","Influenza")) {
  plot(x[grep(D, x$Cause), c("date", "VALUE")], ylab = D)
  abline(v = as.Date("2020/03/17"))
}

dateSeq = sort(unique(x$date))
z_seq = (as.integer(dateSeq) - min(as.integer(dateSeq)))/365.25
x$z = (as.integer(x$date)-min(as.integer(x$date)))/365.25
xInfluenza = x[grepl("Influenza", x$Cause) & x$province == "Ontario", ]
xPreCovid = xInfluenza[xInfluenza$date < as.Date("2020/03/01"),]
xPostCovid = xInfluenza[xInfluenza$date >= as.Date("2020/03/01"),]



### Fit Model: sGP1 + overdispersion
t <- xInfluenza$z
t_obs <- xPreCovid$z
t_pred <- xPostCovid$z
upper_region <- max(t)
period1 <- 2
a1 = period1*pi # 1-year
# a2 <- 2*a1 # half-year
# a3 <- 4*a1 # three-month
region <- c(0,upper_region)

n <- nrow(xPreCovid)
Q1_exact <- Compute_precision_Aug_SGP(t[-1], a = a1)
B1_exact <- rbind(0,Compute_design_Aug_sGP(t[-1])[1:(nrow(xPreCovid)-1), ])

# Q2_exact <- Compute_precision_Aug_SGP(t[-1], a = a2)
# B2_exact <- rbind(0,Compute_design_Aug_sGP(t[-1])[1:(nrow(xPreCovid)-1), ])
B2_exact <- as(Matrix::diag(n), "dgTMatrix")
Q2_exact = B2_exact

# Q3_exact <- Compute_precision_Aug_SGP(t[-1], a = a3)
# B3_exact <- rbind(0,Compute_design_Aug_sGP(t[-1])[1:(nrow(xPreCovid)-1), ])

X <- as(cbind(cos(a1*t_obs), sin(a1*t_obs), 1,t_obs), "dgTMatrix")
B <- cbind(B1_exact,B2_exact)

### Assume the d-year prediction SD has PC prior:
pred_SD_prior <- list(u = 0.01, alpha = 0.5)

### Implementing the seasonal sGP:
d_step <- 1
scale1 <- compute_d_step_sGPsd(d = d_step,a = a1)
SD_prior1 <- list(u = pred_SD_prior$u/scale1, alpha = pred_SD_prior$alpha)
# SD_prior1 <- list(u = 0.1, alpha = 0.5)
SD_prior2 <- list(u = 1, alpha = 0.5)

# scale2 <- compute_d_step_sGPsd(d = d_step,a = a2)
# SD_prior2 <- list(u = pred_SD_prior$u/scale2, alpha = pred_SD_prior$alpha)
# SD_prior2 <- list(u = pred_SD_prior$u, alpha = pred_SD_prior$alpha)

# scale3 <- compute_d_step_sGPsd(d = d_step,a = a3)
# SD_prior3 <- list(u = pred_SD_prior$u/scale3, alpha = pred_SD_prior$alpha)
# SD_prior3 <- list(u = pred_SD_prior$u, alpha = pred_SD_prior$alpha)

### Fitting with AGHQ:
tmbdat <- list(
  # Design matrix
  B1 = B1_exact,
  B2 = B2_exact,
  X = X,
  P1 = Q1_exact,
  P2 = Q2_exact,
  logP1det = as.numeric(determinant(Q1_exact, logarithm = T)$modulus),
  logP2det = as.numeric(determinant(Q2_exact, logarithm = T)$modulus),
  # Response
  y = xPreCovid$VALUE,
  # PC Prior params
  u1 = SD_prior1$u,
  alpha1 = SD_prior1$alpha,
  u2 = SD_prior2$u,
  alpha2 = SD_prior2$alpha,
  betaprec = 0.0001
)
tmbparams <- list(
  W = c(rep(0, (ncol(B) + ncol(X)))), # W = c(U,beta); U = B-Spline coefficients
  theta1 = 0,
  theta2 = 0
)
ff <- TMB::MakeADFun(
  data = tmbdat,
  parameters = tmbparams,
  random = "W",
  DLL = "01_Gaussian_Seasonal_Influ_reduced",
  silent = TRUE
)

ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
fitted_mod_all <- aghq::marginal_laplace_tmb(ff,4,c(0,0))


### Inference of hyper-parameter:
prec_marg <- fitted_mod_all$marginals[[1]]
logpostsigma <- compute_pdf_and_cdf(prec_marg,list(totheta = function(x) -2*log(x),fromtheta = function(x) exp(-x/2)),interpolation = 'spline')
with(logpostsigma,plot(transparam,pdf_transparam,type='l', xlab = "seasonal SD (annual period)", ylab = "density"))

prec_marg <- fitted_mod_all$marginals[[2]]
logpostsigma <- compute_pdf_and_cdf(prec_marg,list(totheta = function(x) -2*log(x),fromtheta = function(x) exp(-x/2)),interpolation = 'spline')
with(logpostsigma,plot(transparam,pdf_transparam,type='l', xlab = "overdispersion", ylab = "density"))


# Plot the overall estimate:
samps <- sample_marginal(fitted_mod_all, 3000)
samps_w <- samps$samps[1:ncol(B1_exact),]
samps_beta <- samps$samps[(ncol(B) + 1):(ncol(B) + ncol(X)),]

samps_g <- as.matrix(B1_exact %*% samps_w) + as.matrix(X %*% samps_beta)
samps_g <- as.data.frame(samps_g)
samps_g$t <- t_obs + 2010
samps_g <- arrange(samps_g, by = t)

plot((xPreCovid$VALUE)~I(sort(t_obs) + 2010), type = 'p', col = 'black', ylab = 'Death', xlab = 'x')
matplot(y = exp(samps_g[,1:30]), x = samps_g$t, type = 'l', add = T, lty = 'dashed', col = 'pink')


#### Predict into the COVID period!
B1_full_exact <- rbind(0,Compute_design_Aug_sGP(t[-1])[1:(nrow(xInfluenza)-1), ])
X_full <- as(cbind(cos(a1*t), sin(a1*t), 1, t), "dgTMatrix")
B_full <- cbind(B1_full_exact)

xPostCovid = xInfluenza
samps_g_pred <- as.matrix(B1_full_exact %*% samps_w) + as.matrix(X_full %*% samps_beta)
samps_g_pred <- as.data.frame(samps_g_pred)
samps_g_pred$t <- t + 2010
samps_g_pred <- arrange(samps_g_pred, by = t)
matplot(y = (samps_g_pred[,1:30]), x = samps_g_pred$t, type = 'l', 
        add = F, lty = 'dashed', col = 'pink',
        ylim = c(2,6))

### Interval:
samps_g_pred$mean <- apply(exp(samps_g_pred)[,1:3000],1, mean)
samps_g_pred$upper <- apply(exp(samps_g_pred)[,1:3000],1, quantile, 0.975)
samps_g_pred$lower <- apply(exp(samps_g_pred)[,1:3000],1, quantile, 0.025)


prec_marg <- fitted_mod_all$marginals[[2]]
logpostsigma <- compute_pdf_and_cdf(prec_marg,list(totheta = function(x) -2*log(x),fromtheta = function(x) exp(-x/2)),interpolation = 'spline')
with(logpostsigma,plot(transparam,pdf_transparam,type='l', xlab = "overdispersion", ylab = "density"))
samps_overdis <- sample(size = 3000, x = logpostsigma$transparam, 
                        prob = logpostsigma$pdf_transparam, replace = T)
samps_to_pred <- samps_g_pred %>% filter(t > (2010 + max(t_obs)))
samps_to_pred_over <- samps_to_pred[,1:3000]
samps_over_noise <- rnorm((3000*nrow(samps_to_pred_over)), sd = rep(samps_overdis, each = nrow(samps_to_pred_over)))
samps_g_pred_all <- samps_to_pred_over + matrix(samps_over_noise, nrow = nrow(samps_to_pred_over), ncol = ncol(samps_to_pred_over),
                                                byrow = F)
samps_g_pred_all$t <- samps_to_pred$t
samps_g_pred_all$upper <- apply(exp(samps_g_pred_all)[,1:3000],1, quantile, 0.975)
samps_g_pred_all$lower <- apply(exp(samps_g_pred_all)[,1:3000],1, quantile, 0.025)


pdf(file = paste0(figure_path, "influ_log.pdf"), height = 5, width = 5)
plot((xInfluenza$VALUE)~I(sort(t) + 2010), cex = 0.5, lwd = 0.3, type = 'p', log = "y",
     col = rgb(0,0,0,alpha = 0.3), xlab = "year", ylab = "Influenza Death", ylim = c(10,250))
lines(mean~t, data = samps_g_pred, col = "blue", lwd = 1.3)
lines(upper~t, data = samps_g_pred, lty = 'dashed', col = "red", lwd = 1.3)
lines(lower~t, data = samps_g_pred, lty = 'dashed', col = "red", lwd = 1.3)
dev.off()

pdf(file = paste0(figure_path, "influ.pdf"), height = 5, width = 5)
plot((xInfluenza$VALUE)~I(sort(t) + 2010), cex = 0.5, lwd = 0.3,
     type = 'p', ylim = c(0,250), 
     col = rgb(0,0,0,alpha = 0.3), xlab = "year", ylab = "Influenza Death")
lines(mean~t, data = samps_g_pred, col = "blue", lwd = 1.3)
lines(upper~t, data = samps_g_pred, lty = 'dashed', col = "red", lwd = 1.3)
lines(lower~t, data = samps_g_pred, lty = 'dashed', col = "red", lwd = 1.3)
polygon(c(samps_g_pred_all$t, rev(samps_g_pred_all$t)), 
        c(samps_g_pred_all$upper, rev(samps_g_pred_all$lower)),
        col = "#6BD7AF",
        density = 10, angle = 45)
dev.off()


plot((xInfluenza$VALUE)~I(sort(t) + 2010), type = 'p', col = 'black', ylab = 'Death', xlab = 'x')
matplot(y = exp(samps_g_pred[,1:30]), x = samps_g_pred$t, type = 'l', add = T, lty = 'dashed', col = 'pink')



### Compute Excess Death:

### Convert linear predictors to observations
samps_g_pred_ld <- samps_g_pred_all[,1:3000]
samps_g_pred_all_obs <- matrix(nrow = nrow(samps_g_pred_all), ncol = 0)
for (i in 1:3000) {
  samps_g_pred_all_obs <- cbind(samps_g_pred_all_obs, rpois(lambda = exp(samps_g_pred_ld[,i]), n = length(samps_g_pred_ld[,i])))
}


sample_excess_log <- log(xInfluenza$VALUE[-(1:nrow(xPreCovid))]) - log(samps_g_pred_all_obs[,1:3000])
sample_excess_log_total <- apply(na.omit(sample_excess_log), 2, sum)
pdf(file = paste0(figure_path, "influenza_excess_log.pdf"), height = 5, width = 5)
hist((sample_excess_log_total), xlab = "Influenza Total Excess Mortality (Log)", main = "Posterior Samples", breaks = 30)
dev.off()
sample_excess <- (xInfluenza$VALUE[-(1:nrow(xPreCovid))]) - (samps_g_pred_all_obs[,1:3000])
sample_excess_total <- apply(na.omit(sample_excess), 2, sum)
pdf(file = paste0(figure_path, "influenza_excess.pdf"), height = 5, width = 5)
hist((sample_excess_total), xlab = "Influenza Total Excess Mortality", main = "Posterior Samples", breaks = 30)
dev.off()

sample_excess_log <- cbind(t = (samps_g_pred_all$t), sample_excess_log)
pdf(file = paste0(figure_path, "influenza_log_excess_path.pdf"), height = 5, width = 5)
matplot(y = (sample_excess_log[,2:31]), x = (sample_excess_log[,1]), type = 'l',
        add = F, lty = 'dashed', col = 'pink', ylab = "Influenza Excess Mortality (Log)",
        xlab = "year", main = "Posterior Samples", xaxt="none")
axis(1, seq(2019,2022,1))

dev.off()
sample_excess <- cbind(t = (samps_g_pred_all$t), sample_excess)
pdf(file = paste0(figure_path, "influenza_excess_path.pdf"), height = 5, width = 5)
matplot(y = (sample_excess[,2:31]), x = sample_excess[,1], type = 'l', 
        add = F, lty = 'dashed', col = 'pink', ylab = "Influenza Excess Mortality",
        xlab = "year", main = "Posterior Samples", xaxt="none")
axis(1, seq(2019,2022,1))

dev.off()














