source("load_packages.R")
source("03_functions_seasonal.R")


### construct joint precision:
Q <- joint_prec_construct(t_vec = c(1,2,3), a = 1, sd = 1)

### Check the covariance matrix:
## Using the state space approach:
solve(Q)[c(1,3,5),c(1,3,5)]

## Using the original formula:
K_true <- generate_K_true(1,1)
compute_matrix_given_cov(from = 1, to = 3, m = 3, K = K_true)


## Simulate as VAR model:
samps <- sim_sGP_Var(t = NULL, mesh_size = 0.01, max_t = 30, a = 1, sd = 1)

par(mfrow = c(1,2))

## simulated sample path of f:
plot(samps[,2]~ samps[,1], type = 'l')

## simulated sample path of f':
plot(samps[,3]~ samps[,1], type = 'l')



