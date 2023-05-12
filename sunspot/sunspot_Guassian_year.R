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
library(sGPfit)

source(paste0(source_path, "06_functions_seasonal.R"))

### year Data:
year_data <- read.table(file = paste0(working_path, "/SN_y_tot_V2.0.txt"), sep = "")
names(year_data) <- c("time", "sunspot", "sd", "number")
year_data <- year_data %>% filter(sunspot > 0)
plot(year_data$sunspot~year_data$time, type = 'p')
plot(log(year_data$sunspot)~year_data$time, type = 'p')

### approximately a 11-year cycle
step_head <- 50
pred_SD <- list(u = 1, alpha = 0.5)
pred_SD2 <- list(u = 5, a = 0.5)

x <- year_data$time - min(year_data$time)
region <- c(0,max(x))
n = length(x)
compile(paste0(source_path,"02_sunspot_guassian.cpp"))
dyn.load(dynlib(paste0(source_path,"02_sunspot_guassian")))



#### Assign discrete prior on alpha: 
period_vec <- 1/seq(from = 8,to = 13, by = 0.1)
a_vec <- 2*period_vec*pi
log_ML <- c()
for (i in 1:length(period_vec)) {
  period <- period_vec[i]
  a = 2*period*pi
  scale <- compute_d_step_sGPsd(d = step_head, a = a)
  pred_SD2_used <- OSplines::prior_conversion(d = step_head, prior = pred_SD2, p = 2)
  x <- year_data$time - min(year_data$time)
  region <- c(0,max(x))
  n = length(x)
  k <- 102
  B1 <- Compute_B_sB(x = x, a = a, k = k, region = region, boundary = T)
  B1 <- as(B1,"dgTMatrix")
  B2 <- local_poly(knots = seq(min(x), max(x), length.out = k), p = 2, refined_x = x)
  B2 <- as(B2,"dgTMatrix")
  Q1 <- Compute_Q_sB(a = a, k = k, region = region)
  Q1 <- as(as(Q1, "matrix"),"dgTMatrix")
  Q2 <- as(OSplines::compute_weights_precision(x = seq(min(x), max(x), length.out = k)),"dgTMatrix")
  X1 <- as(cbind(sin(a*x), cos(a*x)), "dgTMatrix")
  X2 <- as(global_poly(x = x, p = 2), "dgTMatrix")
  X <- cbind(X1,X2)
  tmbdat <- list(
    # Design matrix
    B1 = B1,
    B2 = B2,
    X = X,
    P1 = Q1,
    P2 = Q2,
    logP1det = as.numeric(determinant(Q1, logarithm = T)$modulus),
    logP2det = as.numeric(determinant(Q2, logarithm = T)$modulus),
    # Response
    y = log(year_data$sunspot),
    # PC Prior params
    u1 = pred_SD$u/scale,
    alpha1 = pred_SD$alpha,
    u2 = pred_SD2_used$u,
    alpha2 = pred_SD2_used$a,
    u3 = 1,
    alpha3 = 0.1,
    betaprec = .001
  )
  tmbparams <- list(
    W = c(rep(0, (ncol(B1) + ncol(B2) + ncol(X)))), # W = c(U,beta); U = B-Spline coefficients
    theta1 = 0,
    theta2 = 0,
    theta3 = 0
  )
  ff <- TMB::MakeADFun(
    data = tmbdat,
    parameters = tmbparams,
    random = "W",
    DLL = "02_sunspot_guassian",
    silent = TRUE
  )
  ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
  fitted_mod <- aghq::marginal_laplace_tmb(ff,5,c(0,0,0))
  hyper_result <- fitted_mod$marginals
  save(file = paste0(result_path, "samps/",i,"_hyper_sample.rda"), hyper_result)
  sGP_samps <- sample_marginal(fitted_mod, 3000)
  save(file = paste0(result_path, "samps/",i,"_sample.rda"), sGP_samps)
  log_ML[i] <- fitted_mod$normalized_posterior$lognormconst
}
log_ML
ratio_ML <- exp(log_ML - max(log_ML))
Poster_alpha <- ratio_ML/sum(ratio_ML)
pdf(file = paste0(figure_path, "poster_alpha.pdf"))
plot(Poster_alpha~I(1/period_vec), type = 'p', xaxs='i', xlab = "periodicity", ylab = "Posterior")
lines(predict(smooth.spline(x = I(1/period_vec), y = Poster_alpha), x = seq(8,13,by = 0.01)), col = 'red')
abline(h = 1/length(period_vec), col = "blue", lty = "dashed")
dev.off()
pos_alpha_data <- data.frame(alpha = a_vec, Prob = Poster_alpha)
save(file = paste0(result_path, "pos_alpha_data.rda"), pos_alpha_data)





### Plot the posterior:
load(file = paste0(result_path, "pos_alpha_data.rda"))
pdf(file = paste0(figure_path, "poster_alpha.pdf"))
plot(pos_alpha_data$Prob~I(1/period_vec), type = 'p', xaxs='i', xlab = "periodicity (years)", ylab = "Posterior")
lines(predict(smooth.spline(x = I(1/period_vec), y = pos_alpha_data$Prob), x = seq(8,13,by = 0.01)), col = 'red')
abline(h = 1/length(period_vec), col = "blue", lty = "dashed")
dev.off()


### Fit the final model:
nsamps = 3000
pos_samples_alpha_vec <- sample(x = pos_alpha_data$alpha, size = nsamps, prob = pos_alpha_data$Prob, replace = T)
pos_samples_alpha <- table(pos_samples_alpha_vec) 
original_count <- rep(0, nrow(pos_alpha_data))
names(original_count) <- c(pos_alpha_data$alpha)
pos_samples_alpha <- rev(tapply(c(original_count, pos_samples_alpha), names(c(original_count, pos_samples_alpha)), "sum"))
sd1 <- c()
sd2 <- c()
sd3 <- c()
samps_g1 <- data.frame(x = x)
samps_g2 <- data.frame(x = x)
samps_g <- data.frame(x = x)


for (i in 1:length(pos_samples_alpha)) {
  nsampsi <- pos_samples_alpha[i]
  if(nsampsi != 0){
    a <- pos_alpha_data$alpha[i] ## the third component
    k <- 102
    B1 <- Compute_B_sB(x = x, a = a, k = k, region = region, boundary = T)
    B1 <- as(B1,"dgTMatrix")
    B <- cbind(B1,B2)
    X1 <- as(cbind(sin(a*x), cos(a*x)), "dgTMatrix")
    X <- cbind(X1,X2)
    
    load(file = paste0(result_path, "samps/", i, "_sample.rda"))
    sGP_samps$samps <- sGP_samps$samps[,(1:nsampsi), drop = F]
    samps_w <- sGP_samps$samps[1:ncol(B), , drop = F]
    samps_beta <- sGP_samps$samps[(ncol(B) + 1):(ncol(B) + ncol(X)),, drop = F]
    samps_g1_new <- as.matrix(B1 %*% samps_w[1:ncol(B1),, drop = F]) + as.matrix(X[,1:2] %*% samps_beta[1:2,, drop = F])
    samps_g1_new <- as.data.frame(samps_g1_new)
    samps_g1 <- cbind(samps_g1,samps_g1_new)
    
    samps_g2_new <- as.matrix(B2 %*% samps_w[(ncol(B1) + 1):(ncol(B)),, drop = F]) + as.matrix(X[,3:4] %*% samps_beta[3:4,, drop = F])
    samps_g2_new <- as.data.frame(samps_g2_new)
    samps_g2 <- cbind(samps_g2,samps_g2_new)
    
    samps_g_new <- samps_g1_new + samps_g2_new
    samps_g <- cbind(samps_g,samps_g_new)
    

    sd1 <- c(sd1, exp(-0.5*sGP_samps$theta[(1:nsampsi),1]))
    sd2 <- c(sd2, exp(-0.5*sGP_samps$theta[(1:nsampsi),2]))
    sd3 <- c(sd2, exp(-0.5*sGP_samps$theta[(1:nsampsi),3]))
  }
}

plot(year_data$sunspot~year_data$time, type = 'p', col = 'black', ylab = 'sunspot', xlab = 'year')
matplot(y = exp(samps_g[,2:30]), x = year_data$time, type = 'l', add = T, lty = 'dashed', col = 'pink')

pos_mean <- apply((samps_g[,2:3000]),1, mean)
pos_upper <- apply((samps_g[,2:3000]),1, quantile, p = 0.975)
pos_lower <- apply((samps_g[,2:3000]),1, quantile, p = 0.025)

pdf(file = paste0(figure_path, "log_sunspot.pdf"))
plot(log(sunspot)~time, data = year_data, lwd = 0.5, cex = 1, type = 'p', ylim = c(-2,8))
lines(pos_upper~year_data$time, col = "red", lty = "dashed")
lines(pos_lower~year_data$time, col = "red", lty = "dashed")
lines(pos_mean~year_data$time, lty = 'solid', col = "blue")
dev.off()

pos_mean <- apply(exp(samps_g[,2:3000]),1, mean)
pos_median <- apply(exp(samps_g[,2:3000]),1, quantile, p = 0.5)
pos_upper <- apply(exp(samps_g[,2:3000]),1, quantile, p = 0.975)
pos_lower <- apply(exp(samps_g[,2:3000]),1, quantile, p = 0.025)

pdf(file = paste0(figure_path, "sunspot.pdf"))
plot(pos_mean~year_data$time, type = 'n', col = "blue", ylim = c(0,400), ylab = "sunspot", xlab = "time")
points((sunspot)~time, data = year_data, lwd = 0.5, cex = 1)
lines(pos_upper~year_data$time, col = "red", lty = "dashed")
lines(pos_lower~year_data$time, col = "red", lty = "dashed")
lines(pos_mean~year_data$time, col = "blue")
dev.off()


pos_mean <- apply((samps_g1[,2:3000]),1, mean)
pos_median <- apply((samps_g1[,2:3000]),1, quantile, p = 0.5)
pos_upper <- apply((samps_g1[,2:3000]),1, quantile, p = 0.975)
pos_lower <- apply((samps_g1[,2:3000]),1, quantile, p = 0.025)

pdf(file = paste0(figure_path, "sunspot_seasonal.pdf"))
plot(pos_mean~year_data$time, type = 'l', col = "blue", ylim = c(-3,3), ylab = "seasonal", xlab = "time")
lines(pos_upper~year_data$time, col = "red", lty = "dashed")
lines(pos_lower~year_data$time, col = "red", lty = "dashed")
dev.off()


pos_mean <- apply((samps_g2[,2:3000]),1, mean)
pos_median <- apply((samps_g2[,2:3000]),1, quantile, p = 0.5)
pos_upper <- apply((samps_g2[,2:3000]),1, quantile, p = 0.975)
pos_lower <- apply((samps_g2[,2:3000]),1, quantile, p = 0.025)

pdf(file = paste0(figure_path, "sunspot_trend.pdf"))
plot(pos_mean~year_data$time, type = 'l', col = "blue", ylab = "trend", xlab = "time", ylim = c(0,8))
lines(pos_upper~year_data$time, col = "red", lty = "dashed")
lines(pos_lower~year_data$time, col = "red", lty = "dashed")
dev.off()
