#include <TMB.hpp>                                // Links in the TMB libraries
//#include <fenv.h>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(y); //response variable
  DATA_SPARSE_MATRIX(B1); // Design matrix (for random effects 1; 1-year SGP)
  DATA_SPARSE_MATRIX(B2); // Design matrix (for random effects 2; 6-month SGP)
  DATA_SPARSE_MATRIX(B3); // Design matrix (for random effects 3; 44-month SGP)
  DATA_SPARSE_MATRIX(B4); // Design matrix (for random effects 4; 9.1-year lunar SGP)
  DATA_SPARSE_MATRIX(B5); // Design matrix (for random effects 5; 10.4-year solar SGP)
  DATA_SPARSE_MATRIX(B6); // Design matrix (for random effects 6; IWP trend)

  DATA_SPARSE_MATRIX(X); // Design matrix (for fixed effects)
  DATA_SPARSE_MATRIX(P1); // Penalty matrix
  DATA_SPARSE_MATRIX(P2); // Penalty matrix
  DATA_SPARSE_MATRIX(P3); // Penalty matrix
  DATA_SPARSE_MATRIX(P4); // Penalty matrix
  DATA_SPARSE_MATRIX(P5); // Penalty matrix
  DATA_SPARSE_MATRIX(P6); // Penalty matrix

  int d1 = P1.cols(); // Number of B-Spline coefficients
  int d2 = P2.cols(); // Number of B-Spline coefficients
  int d3 = P3.cols(); // Number of B-Spline coefficients
  int d4 = P4.cols(); // Number of B-Spline coefficients
  int d5 = P5.cols(); // Number of B-Spline coefficients
  int d6 = P6.cols(); // Number of B-Spline coefficients

  DATA_SCALAR(logP1det); // Determinant of (fixed) penalty matrix
  DATA_SCALAR(logP2det); // Determinant of (fixed) penalty matrix
  DATA_SCALAR(logP3det); // Determinant of (fixed) penalty matrix
  DATA_SCALAR(logP4det); // Determinant of (fixed) penalty matrix
  DATA_SCALAR(logP5det); // Determinant of (fixed) penalty matrix
  DATA_SCALAR(logP6det); // Determinant of (fixed) penalty matrix

  DATA_SCALAR(u1); // pc prior, u1 param
  DATA_SCALAR(alpha1); // pc prior, alpha1 param
  DATA_SCALAR(u2); // pc prior, u2 param
  DATA_SCALAR(alpha2); // pc prior, alpha2 param
  DATA_SCALAR(u3); // pc prior, u1 param
  DATA_SCALAR(alpha3); // pc prior, alpha1 param
  DATA_SCALAR(u4); // pc prior, u1 param
  DATA_SCALAR(alpha4); // pc prior, alpha1 param
  DATA_SCALAR(u5); // pc prior, u1 param
  DATA_SCALAR(alpha5); // pc prior, alpha1 param
  DATA_SCALAR(u6); // pc prior, u1 param
  DATA_SCALAR(alpha6); // pc prior, alpha1 param
  DATA_SCALAR(u7); // pc prior, u1 param
  DATA_SCALAR(alpha7); // pc prior, alpha1 param
  DATA_SCALAR(betaprec); // precision for beta

  // Parameter
  PARAMETER_VECTOR(W); //
  int Wdim = W.size();
  vector<Type> U1(d1);
  vector<Type> U2(d2);
  vector<Type> U3(d3);
  vector<Type> U4(d4);
  vector<Type> U5(d5);
  vector<Type> U6(d6);

  int betadim = X.cols();
  vector<Type> beta(betadim);
  for (int i=0;i<d1;i++) U1(i) = W(i);
  for (int i=0;i<d2;i++) U2(i) = W(i + d1);
  for (int i=0;i<d3;i++) U3(i) = W(i + d1 + d2);
  for (int i=0;i<d4;i++) U4(i) = W(i + d1 + d2 + d3);
  for (int i=0;i<d5;i++) U5(i) = W(i + d1 + d2 + d3 + d4);
  for (int i=0;i<d6;i++) U6(i) = W(i + d1 + d2 + d3 + d4 + d5);
  for (int i=0;i<betadim;i++) beta(i) = W(i + d1 + d2 + d3 + d4 + d5 + d6);

  PARAMETER(theta1); // theta = -2log(sigma1)
  PARAMETER(theta2); // theta = -2log(sigma2)
  PARAMETER(theta3); // theta = -2log(sigma3)
  PARAMETER(theta4); // theta = -2log(sigma4)
  PARAMETER(theta5); // theta = -2log(sigma5)
  PARAMETER(theta6); // theta = -2log(sigma6)
  PARAMETER(theta7); // theta = -2log(sigma7)

  // Transformations
  vector<Type> eta = X * beta + B1 * U1 + B2 * U2 + B3 * U3 + B4 * U4 + B5 * U5 + B6 * U6;
  Type sigma1 = exp(-0.5*theta1);
  Type sigma2 = exp(-0.5*theta2);
  Type sigma3 = exp(-0.5*theta3);
  Type sigma4 = exp(-0.5*theta4);
  Type sigma5 = exp(-0.5*theta5);
  Type sigma6 = exp(-0.5*theta6);
  Type sigma7 = exp(-0.5*theta7);

  // Log likelihood
  Type ll = 0;
  ll = sum(dnorm(y, eta, sigma7, TRUE));
  
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
  vector<Type> P4U4 = P4*U4;
  Type U4P4U4 = (U4 * P4U4).sum();
  lpW += -0.5 * exp(theta4) * U4P4U4; // U part
  vector<Type> P5U5 = P5*U5;
  Type U5P5U5 = (U5 * P5U5).sum();
  lpW += -0.5 * exp(theta5) * U5P5U5; // U part
  vector<Type> P6U6 = P6*U6;
  Type U6P6U6 = (U6 * P6U6).sum();
  lpW += -0.5 * exp(theta6) * U6P6U6; // U part
  Type bb = (beta * beta).sum();
  lpW += -0.5 * betaprec * bb; // Beta part


  // Log determinant
  Type logdet1 = d1 * theta1 + logP1det;
  Type logdet2 = d2 * theta2 + logP2det;
  Type logdet3 = d3 * theta3 + logP3det;
  Type logdet4 = d4 * theta4 + logP4det;
  Type logdet5 = d5 * theta5 + logP5det;
  Type logdet6 = d6 * theta6 + logP6det;


  Type logdetall = (logdet1 + logdet2 + logdet3 + logdet4 + logdet5 + logdet6);
  lpW += 0.5 * logdetall; // P part
  
  // Log prior for theta
  Type lpT = 0;
  Type phi1 = -log(alpha1) / u1;
  lpT += log(0.5 * phi1) - phi1*exp(-0.5*theta1) - 0.5*theta1;
  Type phi2 = -log(alpha2) / u2;
  lpT += log(0.5 * phi2) - phi2*exp(-0.5*theta2) - 0.5*theta2;
  Type phi3 = -log(alpha3) / u3;
  lpT += log(0.5 * phi3) - phi3*exp(-0.5*theta3) - 0.5*theta3;
  Type phi4 = -log(alpha4) / u4;
  lpT += log(0.5 * phi4) - phi4*exp(-0.5*theta4) - 0.5*theta4;
  Type phi5 = -log(alpha5) / u5;
  lpT += log(0.5 * phi5) - phi5*exp(-0.5*theta5) - 0.5*theta5;
  Type phi6 = -log(alpha6) / u6;
  lpT += log(0.5 * phi6) - phi6*exp(-0.5*theta6) - 0.5*theta6;
  Type phi7 = -log(alpha7) / u7;
  lpT += log(0.5 * phi7) - phi7*exp(-0.5*theta7) - 0.5*theta7;
  
  // Final result!
  Type logpost = -1 * (ll + lpW + lpT);
  
  return logpost;
}