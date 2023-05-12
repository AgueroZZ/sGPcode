### Simulate n samples of SD given hyper_result list:
simulate_hyper <- function(hyper_result, nsamps = 3000, which_index = 1){
  
  prec_marg <- hyper_result[[which_index]]
  logpostsigma <- aghq::compute_pdf_and_cdf(prec_marg,list(totheta = function(x) -2*log(x),fromtheta = function(x) exp(-x/2)),interpolation = 'spline')
  hyper_dist <- data.frame(SD = logpostsigma$transparam, 
                              density = logpostsigma$pdf_transparam)
  
  hyper_samps <- sample(x = hyper_dist$SD, size = nsamps, prob = hyper_dist$density, replace = T)
  hyper_samps
}