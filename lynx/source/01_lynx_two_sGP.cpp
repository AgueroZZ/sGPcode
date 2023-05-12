#include <TMB.hpp>                                // Links in the TMB libraries
//#include <fenv.h>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(y); //response variable
  DATA_SPARSE_MATRIX(B1); // Design matrix (for random effects 1)
  DATA_SPARSE_MATRIX(B2); // Design matrix (for random effects 2)
  DATA_SPARSE_MATRIX(B3); // Design matrix (for random effects 3)


  DATA_SPARSE_MATRIX(X); // Design matrix (for fixed effects)
  DATA_SPARSE_MATRIX(P1); // Penalty matrix
  DATA_SPARSE_MATRIX(P2); // Penalty matrix
  DATA_SPARSE_MATRIX(P3); // Penalty matrix


  int d1 = P1.cols(); // Number of B-Spline coefficients
  int d2 = P2.cols(); // Number of B-Spline coefficients
  int d3 = P3.cols(); // Number of B-Spline coefficients


  DATA_SCALAR(logP1det); // Determinant of (fixed) penalty matrix
  DATA_SCALAR(logP2det); // Determinant of (fixed) penalty matrix
  DATA_SCALAR(logP3det); // Determinant of (fixed) penalty matrix


  DATA_SCALAR(u1); // pc prior, u1 param
  DATA_SCALAR(alpha1); // pc prior, alpha1 param
  DATA_SCALAR(u2); // pc prior, u2 param
  DATA_SCALAR(alpha2); // pc prior, alpha2 param
  DATA_SCALAR(u3); // pc prior, u3 param
  DATA_SCALAR(alpha3); // pc prior, alpha3 param
  DATA_SCALAR(betaprec); // precision for beta

  // Parameter
  PARAMETER_VECTOR(W); // W = c(U), eta = B * U
  int Wdim = W.size();
  vector<Type> U1(d1);
  vector<Type> U2(d2);
  vector<Type> U3(d3);
  int betadim = X.cols();
  vector<Type> beta(betadim);
  for (int i=0;i<d1;i++) U1(i) = W(i);
  for (int i=0;i<d2;i++) U2(i) = W(d1 + i);
  for (int i=0;i<d3;i++) U3(i) = W(d1 + d2 + i);
  for (int i=0;i<betadim;i++) beta(i) = W(i + d1 + d2 + d3);
  PARAMETER(theta1); // theta1 = -2log(sigma)
  PARAMETER(theta2); // theta2 = -2log(sigma)
  PARAMETER(theta3); // theta3 = -2log(sigma)



  // Transformations
  vector<Type> eta = X * beta + B1 * U1 + B2 * U2 + B3 * U3;

  // Log likelihood
  Type ll = 0;
  ll = sum(dpois(y, exp(eta), TRUE));
  REPORT(ll);
  
  // Log prior on W
  Type lpW = 0;

  // Cross product
  vector<Type> P1U1 = P1*U1;
  Type U1P1U1 = (U1 * P1U1).sum();
  lpW += -0.5 * exp(theta1) * U1P1U1; // U part
  vector<Type> P2U2 = P2*U2;
  Type U2P2U2 = (U2 * P2U2).sum();
  lpW += -0.5 * exp(theta2) * U2P2U2; // U part
  vector<Type> P3U3 = P3*U3;
  Type U3P3U3 = (U3 * P3U3).sum();
  lpW += -0.5 * exp(theta3) * U3P3U3; // U part
  Type bb = (beta * beta).sum();
  lpW += -0.5 * betaprec * bb; // Beta part


  // Log determinant
  Type logdet = d1 * theta1 + logP1det;
  logdet += d2 * theta2 + logP2det;
  logdet += d3 * theta3 + logP3det;
  lpW += 0.5 * logdet; // P part

  REPORT(logdet);
  REPORT(lpW);

  
  // Log prior for theta
  Type lpT = 0;
  Type phi1 = -log(alpha1) / u1;
  lpT += log(0.5 * phi1) - phi1*exp(-0.5*theta1) - 0.5*theta1;
  Type phi2 = -log(alpha2) / u2;
  lpT += log(0.5 * phi2) - phi2*exp(-0.5*theta2) - 0.5*theta2;
  Type phi3 = -log(alpha3) / u3;
  lpT += log(0.5 * phi3) - phi2*exp(-0.5*theta3) - 0.5*theta3;
  REPORT(lpT);
  
  // Final result!
  Type logpost = -1 * (ll + lpW + lpT);
  
  return logpost;
}