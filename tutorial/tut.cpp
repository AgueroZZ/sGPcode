#include <TMB.hpp>                                // Links in the TMB libraries
//#include <fenv.h>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(y); //response variable
  DATA_SPARSE_MATRIX(B); // Design matrix (for random effects 1)
  DATA_SPARSE_MATRIX(X); // Design matrix (for fixed effects)
  DATA_SPARSE_MATRIX(P); // Penalty matrix
  DATA_SCALAR(logPdet); // Determinant of (fixed) penalty matrix
  DATA_SCALAR(u); // pc prior, u param
  DATA_SCALAR(alpha); // pc prior, alpha param
  DATA_SCALAR(u_over); // pc prior, u param
  DATA_SCALAR(alpha_over); // pc prior, alpha param
  DATA_SCALAR(betaprec); // precision for beta

  // Parameter
  PARAMETER_VECTOR(W); 
  int d = P.cols(); // Number of random effect elements
  vector<Type> U(d);
  int betadim = X.cols();
  vector<Type> beta(betadim);
  int n = y.size(); // Number of data points
  vector<Type> O(n);
  for (int i=0;i<d;i++) U(i) = W(i);
  for (int i=0;i<betadim;i++) beta(i) = W(i + d);
  for (int i=0;i<n;i++) O(i) = W(i + d + betadim);

  PARAMETER(theta); // theta = -2log(sigma)
  PARAMETER(theta_over); // theta = -2log(sigma)

  // Transformations
  vector<Type> eta =  B * U + X * beta + O;

  // Log likelihood
  Type ll = 0;
  ll = sum(dpois(y, exp(eta), TRUE));
  
  // Log prior on W
  Type lpW = 0;

  // Cross product
  vector<Type> PU = P*U;
  Type UPU = (U * PU).sum();
  lpW += -0.5 * exp(theta) * UPU; // U part
  Type bb = (beta * beta).sum();
  lpW += -0.5 * betaprec * bb; // Beta part
  Type over_part = (O * O).sum();
  lpW += -0.5 * exp(theta_over) * over_part; // Beta part

  // Log determinant
  Type logdet = d * theta + logPdet;
  lpW += 0.5 * logdet; // P part
  Type logdet_over = n * theta_over;
  lpW += 0.5 * logdet_over; // P part
  
  // Log prior for theta
  Type lpT = 0;
  Type phi = -log(alpha) / u;
  lpT += log(0.5 * phi) - phi*exp(-0.5*theta) - 0.5*theta;
  Type phi_over = -log(alpha_over) / u_over;
  lpT += log(0.5 * phi_over) - phi_over*exp(-0.5*theta_over) - 0.5*theta_over;
  
  // Final result!
  Type logpost = -1 * (ll + lpW + lpT);
  
  return logpost;
}