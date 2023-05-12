Compute_design_Aug_sGP <- function(svec){
  p = 2
  Design <- Diagonal((p * length(svec)), x = 0)
  diag(Design)[seq(1,nrow(Design), by = p)] <- 1
  as(as.matrix(Design[seq(1,nrow(Design), by = p),]), "dgTMatrix")
}

Compute_precision_Aug_SGP <- function(svec, a){
  joint_prec_construct(svec, a, sd = 1)
}