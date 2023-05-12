source("01_fit.R")

n = length(x)
B1 <- Diagonal(n = 2*n)[,1:(2*n)]
B1 <- B1[seq(1,2*n,by = 2),][, -c(1:2)]
B1 <- as(as(B1, "matrix"),"dgTMatrix")
B2 <- Diagonal(n = 2*n)[,1:(2*n)]
B2 <- B2[seq(1,2*n,by = 2),][, -c(1:2)]
B2 <- as(as(B2, "matrix"),"dgTMatrix")
B3 <- as(Matrix::diag(n), "dgTMatrix")
Q3 = B3
B <- cbind(B1,B2,B3)
BB1 <- Diagonal(n = 2*n)[,1:(2*n)]
BB1 <- BB1[seq(2,2*n,by = 2),][, -c(1:2)]
BB2 <- BB1

load(file = paste0(result_path, "pos_alpha_data.rda"))
set.seed(123)
nsamps <- 3000
pos_samples_alpha_vec <- sample(x = pos_alpha_data$alpha, size = nsamps, prob = pos_alpha_data$Prob, replace = T)
pos_samples_alpha <- table(pos_samples_alpha_vec) 
original_count <- rep(0, nrow(pos_alpha_data))
names(original_count) <- c(pos_alpha_data$alpha)
pos_samples_alpha <- rev(tapply(c(original_count, pos_samples_alpha), names(c(original_count, pos_samples_alpha)), "sum"))

samps_const <- data.frame(x = x)
samps_g1 <- data.frame(x = x)
samps_g2 <- data.frame(x = x)
samps_g1_1st <- data.frame(x = x)
samps_g2_1st <- data.frame(x = x)
sd1 <- c()
sd2 <- c()

for (i in 1:length(pos_samples_alpha)) {
  nsampsi <- pos_samples_alpha[i]
  if(nsampsi != 0){
    a1 <- pos_alpha_data$alpha[i] ## the third component
    a2 <- 2*a1
    Q1 <- joint_prec_construct(a = a1, t_vec = x[-1], sd = 1)
    Q1 <- as(as(Q1, "matrix"),"dgTMatrix")
    Q2 <- joint_prec_construct(a = a2, t_vec = x[-1], sd = 1)
    Q2 <- as(as(Q2, "matrix"),"dgTMatrix")
    X <- as(cbind(cos(a1*x),sin(a1*x),cos(a2*x),sin(a2*x),1), "dgTMatrix")
    XX1 <- as(cbind(-a1*sin(a1*x), a1*cos(a1*x)), "dgTMatrix")
    XX2 <- as(cbind(-a2*sin(a2*x), a2*cos(a2*x)), "dgTMatrix")
    load(file = paste0(result_path, "samps/", i, "_sample.rda"))
    sGP_samps$samps <- sGP_samps$samps[,(1:nsampsi), drop = F]
    samps_w <- sGP_samps$samps[1:ncol(B), , drop = F]
    samps_beta <- sGP_samps$samps[(ncol(B) + 1):(ncol(B) + ncol(X)),, drop = F]
    samps_const_new <- samps_beta[5,,drop = F]
    samps_const <- suppressWarnings(cbind(samps_const, samps_const_new))
    samps_g1_new <- as.matrix(B1 %*% samps_w[1:ncol(B1),, drop = F]) + as.matrix(X[,1:2] %*% samps_beta[1:2,, drop = F])
    samps_g1_new <- as.data.frame(samps_g1_new)
    samps_g1 <- cbind(samps_g1,samps_g1_new)
    samps_g1_1st_new <- as.matrix(XX1 %*% samps_beta[1:2,, drop = F]) + as.matrix(BB1 %*% samps_w[1:ncol(B1),, drop = F])
    samps_g1_1st_new <- as.data.frame(samps_g1_1st_new)
    samps_g1_1st <- cbind(samps_g1_1st,samps_g1_1st_new)
    samps_g2_new <- as.matrix(B2 %*% samps_w[(ncol(B1) + 1):(ncol(B1) + ncol(B2)),, drop = F]) + as.matrix(X[,3:4] %*% samps_beta[3:4,, drop = F])
    samps_g2_new <- as.data.frame(samps_g2_new)
    samps_g2 <- cbind(samps_g2,samps_g2_new)
    samps_g2_1st_new <- as.matrix(XX2 %*% samps_beta[3:4,, drop = F]) + as.matrix(BB2 %*% samps_w[(ncol(B1) + 1):(ncol(B1) + ncol(B2)),, drop = F])
    samps_g2_1st_new <- as.data.frame(samps_g2_1st_new)
    samps_g2_1st <- cbind(samps_g1_1st,samps_g2_1st_new)
    sd1 <- c(sd1, exp(-0.5*sGP_samps$theta[(1:nsampsi),1]))
    sd2 <- c(sd2, exp(-0.5*sGP_samps$theta[(1:nsampsi),2]))
  }
}

### Making prediction:
max_t <- 150 ### predict up to time 150
mesh_size <- 0.5
pred_x <- seq(max(x), max_t, by = mesh_size)
pred_g1 <- data.frame(x = pred_x)
pred_g2 <- data.frame(x = pred_x)

for (i in 1:length(pos_samples_alpha_vec)) {
  a1 <- pos_samples_alpha_vec[i] ## the third component
  a2 <- 2*a1
  pred_g1_new <- predict_sGP_Var(observed_f_vec = c(samps_g1[nrow(samps_g1),(1+i)], samps_g1_1st[nrow(samps_g1_1st),(1+i)]),
                                 last_t = min(pred_g1$x), max_t = max(pred_g1$x), a = a1, sd = sd1[i], mesh_size = mesh_size)
  pred_g1 <- cbind(pred_g1, pred_g1_new[,2])
  pred_g2_new <- predict_sGP_Var(observed_f_vec = c(samps_g2[nrow(samps_g2),(1+i)], samps_g2_1st[nrow(samps_g2_1st),(1+i)]),
                                 last_t = min(pred_g2$x), max_t = max(pred_g2$x), a = a2, sd = sd2[i], mesh_size = mesh_size)
  pred_g2 <- cbind(pred_g2, pred_g2_new[,2])
}


### Plotting:
### Overall:
pred_beta_const <- suppressWarnings(cbind(x = pred_g1$x, samps_const[1,-1]))
samps_overall<- samps_g1 + samps_g2 + samps_const
samps_overall$x <- x
upper_int <- exp(samps_overall[,-1]) %>% apply(MARGIN = 1, FUN = quantile, p = 0.975)
lower_int <- exp(samps_overall[,-1]) %>% apply(MARGIN = 1, FUN = quantile, p = 0.025)
mean <- exp(samps_overall[,-1]) %>% apply(MARGIN = 1, FUN = mean)
pred_overall<- pred_g1 + pred_g2 + pred_beta_const
pred_overall$x <- pred_x
pred_upper_int <- exp(pred_overall[,-1]) %>% apply(MARGIN = 1, FUN = quantile, p = 0.975)
pred_lower_int <- exp(pred_overall[,-1]) %>% apply(MARGIN = 1, FUN = quantile, p = 0.025)
pred_mean <- exp(pred_overall[,-1]) %>% apply(MARGIN = 1, FUN = mean)
pdf(file = paste0(figure_path, "lynx_predict.pdf"), height = 5, width = 5)
plot(y~I(x + min(data$year)), type = 'p', xlab = "year", ylab = "lynx (log10)", xlim = c(1820,1971), lwd = 0.5, cex = 1, ylim = c(0,10000))
lines(mean~I(x + min(data$year)), col = 'blue', xlab = "year", ylab = "lynx (log10)", xlim = c(1820,1971))
lines(upper_int~I(x + min(data$year)), lty = 'dashed', col = 'red')
lines(lower_int~I(x + min(data$year)), lty = 'dashed', col = 'red')
lines(pred_mean~I(pred_x + min(data$year)), lty = 'solid', col = "blue")
lines(pred_upper_int~I(pred_x + min(data$year)), lty = 'dashed', col = 'red')
lines(pred_lower_int~I(pred_x + min(data$year)), lty = 'dashed', col = 'red')
dev.off()


### Overall (Log):
pred_beta_const <- suppressWarnings(cbind(x = pred_g1$x, samps_const[1,-1]))
samps_overall<- samps_g1 + samps_g2 + samps_const
samps_overall$x <- x
upper_int <- (samps_overall[,-1]) %>% apply(MARGIN = 1, FUN = quantile, p = 0.975)
lower_int <- (samps_overall[,-1]) %>% apply(MARGIN = 1, FUN = quantile, p = 0.025)
mean <- (samps_overall[,-1]) %>% apply(MARGIN = 1, FUN = mean)
pred_overall<- pred_g1 + pred_g2 + pred_beta_const
pred_overall$x <- pred_x
pred_upper_int <- (pred_overall[,-1]) %>% apply(MARGIN = 1, FUN = quantile, p = 0.975)
pred_lower_int <- (pred_overall[,-1]) %>% apply(MARGIN = 1, FUN = quantile, p = 0.025)
pred_mean <- (pred_overall[,-1]) %>% apply(MARGIN = 1, FUN = mean)
pdf(file = paste0(figure_path, "lynx_log_predict.pdf"), height = 5, width = 5)
plot(log(y)~I(x + min(data$year)), type = 'p', xlab = "year", ylab = "lynx (log10)", xlim = c(1820,1971), lwd = 0.5, cex = 1, ylim = c(0,10))
lines(mean~I(x + min(data$year)), col = 'blue', xlab = "year", ylab = "lynx (log10)", xlim = c(1820,1971))
lines(upper_int~I(x + min(data$year)), lty = 'dashed', col = 'red')
lines(lower_int~I(x + min(data$year)), lty = 'dashed', col = 'red')
lines(pred_mean~I(pred_x + min(data$year)), lty = 'solid', col = "blue")
lines(pred_upper_int~I(pred_x + min(data$year)), lty = 'dashed', col = 'red')
lines(pred_lower_int~I(pred_x + min(data$year)), lty = 'dashed', col = 'red')
dev.off()



