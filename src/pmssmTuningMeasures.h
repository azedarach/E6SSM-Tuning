/*
  pmssmTuningMeasures.h provides the function headers and prototypes
  for the functions that calculate the values of the EWSB conditions
  and fine tunings in the pMSSM.
  Author: Dylan Harries
  Date created: 5/2/2014
 */

#ifndef PMSSMTUNINGMEASURES_H
#define PMSSMTUNINGMEASURES_H

#include <math.h>
#include <iostream>
#include "softsusy.h"

// pMSSM_EWSBConditioni_TreeLevel, i = 1, 2, calculates the value of the 
// EWSB condition for m_H_i^2.
// Inputs:
//     SoftParsMssm r = pMSSM model object
//     double tb = value of tan(beta)
//     int l = # loops used. Can be either 0 (tree level EWSB conditions)
//             or 1 (leading one-loop tadpole corrections included)
double pMSSM_EWSBCondition1(SoftParsMssm, double, int);
double pMSSM_EWSBCondition2(SoftParsMssm, double, int);

// pMSSM_ImplementEWSBConstraints fixes the soft masses m_H_1^2 and
// m_H_2^2 to satisfy the EWSB conditions at the requested order in
// perturbation theory.
// Inputs:
//     SoftParsMssm & r = pMSSM model object
//     double tb = value of tan(beta)
//     int l = # loops used. Can be either 0 (tree level EWSB conditions)
//             or 1 (leading one-loop tadpole corrections included)
void pMSSM_ImplementEWSBConstraints(SoftParsMssm &, double, int);

// pMSSM_FineTunings calculates the fine tuning sensitivities 
// \Delta_a for the parameters a = \mu, B, m_1^2, m_2^2, m_Q^2, m_u^2, A_t at
// the requested order in perturbation theory.
// Inputs:
//     SoftParsMssm const & r = pMSSM model object
//     double tb = the value of tan(beta) to use
//     int l = # loops used. Can be either 0 (tree level EWSB conditions)
//             or 1 (leading one-loop tadpole corrections included)
//     DoubleVector & deltas = a vector (of length 7) to store the sensitivies,
//                             in the order \mu, B, m_1^2, m_2^2, m_Q^2, m_u^2, 
//                             A_t.
//     int & sing = integer flag to indicate if QR solution is singular
//     double & detF = double to store the determinant of the matrix used in 
//                     calculating the tuning
void pMSSM_FineTunings(SoftParsMssm const &, double, int, DoubleVector &, int &, double &);

// Calculates the DRbar' masses of the stops in the pMSSM.
// Inputs:
//     SoftParsMssm r = pMSSM model
//     DoubleVector & mstop = masses of stops
//     DoubleVector & mstopsq = squared masses of stops
//     double tb = value of tan(beta)
void physical_pMSSM(SoftParsMssm, DoubleVector &, DoubleVector &, double);

// Function needed to evaluate tadpoles.
// Inputs:
//     double mSq = squared mass m^2
//     double Q = renormalisation scale Q
double a0pMSSM(double, double);

// Calculates the tadpole corrections to the EWSB conditions. Note
// returns (1/v_i)\partial \Delta V/\partial v_i (i.e. opposite sign
// to equivalent E_6SSM methods).
// Inputs:
//     SoftParsMssm r = pMSSM model
//     double tb = tan(beta)
double doCalcTadpolepMSSMH1(SoftParsMssm, double);
double doCalcTadpolepMSSMH2(SoftParsMssm, double);

// The following are helper functions useful for constructing the 
// tuning measures at one-loop order.

// pMSSMdmstop1sqdvi calculates the derivative m_stop_1^2 wrt v_i.
// Note our convention is m_stop_1 is the lighter stop.
// Inputs:
//     SoftParsMssm r = pMSSM model
//     double tb = tan(beta)
double doCalcpMSSMdmstop1sqdv1(SoftParsMssm, double);
double doCalcpMSSMdmstop1sqdv2(SoftParsMssm, double);

// As above but for m_stop_2^2.
double doCalcpMSSMdmstop2sqdv1(SoftParsMssm, double);
double doCalcpMSSMdmstop2sqdv2(SoftParsMssm, double);

// These functions calculate the second derivatives of the one-loop
// corrections to the effective potential. Returns
// \Delta_{ij}'=\partial^2\Delta V/\partial v_i\partial v_j.
// Inputs:
//     SoftParsMssm r = pMSSM model
//     double tb = tan(beta)
double doCalcpMSSMDeltaPrime11(SoftParsMssm, double);
double doCalcpMSSMDeltaPrime12(SoftParsMssm, double);
double doCalcpMSSMDeltaPrime22(SoftParsMssm, double);

// These helper functions calculate the second derivatives of 
// the one-loop corrections wrt the input parameters. They 
// return (1/v_i)\partial^2\Delta V/\partial p_j\partial v_i.
// Inputs:
//     SoftParsMssm r = pMSSM model
//     double tb = tan(beta)
double doCalcpMSSMd2DeltaVdMudv1(SoftParsMssm, double);
double doCalcpMSSMd2DeltaVdMudv2(SoftParsMssm, double);
double doCalcpMSSMd2DeltaVdAtdv1(SoftParsMssm, double);
double doCalcpMSSMd2DeltaVdAtdv2(SoftParsMssm, double);
double doCalcpMSSMd2DeltaVdmQlSqdv1(SoftParsMssm, double);
double doCalcpMSSMd2DeltaVdmQlSqdv2(SoftParsMssm, double);
double doCalcpMSSMd2DeltaVdmUrSqdv1(SoftParsMssm, double);
double doCalcpMSSMd2DeltaVdmUrSqdv2(SoftParsMssm, double);

/*
  --------------------------------------------------------------------------------------
  Following are methods to calculate the fine tuning numerically. They are used only
  for testing, as the analytic methods when properly checked should be reliable and
  fast enough.
  --------------------------------------------------------------------------------------
 */

// pMSSM_NumericalTadpoles calculates the tadpole corrections
// numerically to verify the analytic expressions.
// Inputs:
//     SoftParsMssm r = pMSSM model
//     double tb = tan(beta)
//     double epsilon = value of epsilon used to estimate partial derivatives
DoubleVector pMSSM_NumericalTadpoles(SoftParsMssm, double, double);

// pMSSM_NumericalDeltaPrime calculates the three second derivatives of the 
// loop corrections to the effective potential, using the analytic
// expressions for the tadpole corrections (i.e. the first derivatives).
// These analytic expressions have been checked against the numerical
// methods for calculating the first derivatives. Outputs the 4 second
// derivatives as a matrix with elements M_{ij}=\Delta_{ij}'.
// Inputs:
//     SoftParsMssm r = pMSSM model
//     double tb = tan(beta)
//     double epsilon = value of epsilon used in estimating the derivative
DoubleMatrix pMSSM_NumericalDeltaPrime(SoftParsMssm, double, double);

// pMSSM_EWSBNewtonSolver computes the solution to the 
// EWSB conditions at either tree level or one loop order using Newton's method.
// Returns the solution for the VEVs v_1, v_2 as a DoubleVector.
// Note that the provided object is taken to contain both the values of
// the input parameters and the initial guess for the VEVs 
// Inputs:
//     SoftParsMssm r = the object to calculate the solution for
//     double tb = the initial guess for the value of tan(beta)
//     int l = # loops to use (0 or 1)
//     double tol = the required tolerance to test for convergence
//     int maxIters = the maximum number of allowed iterations
//     double epsilon = value of epsilon used to numerically calculate Jacobian
//     int & sing = integer flag to indicate if singular (in which case sing is non-zero)
DoubleVector pMSSM_EWSBNewtonSolver(SoftParsMssm, double, int, double, int, double, int &);

// pMSSM_vevNumericalDerivatives computes estimates for the partial
// derivatives of the Higgs VEVs wrt the low energy input parameters
// \mu, B,... at either tree level or one loop order. This is done
// using Newton's method and a 5-pt stencil approximation to the derivative. Note that
// the result calculated is the directional derivative in the direction specified by
// epsVec. The parameter point at which the derivative is to be evaluated, together
// with the initial estimate for the VEVs, is assumed to be contained in the provided
// MSSM model object.
// Inputs:
//     SoftParsMssm r = the object to calculate the derivatives for
//     double tb = the value of tan(beta) to use
//     int l = # loops to use(0 or 1)
//     DoubleVector epsVec = a vector indicating the direction in the low-energy
//                           parameter space in which to calculate the directional
//                           derivative. For the tree level method, epsVec must be
//                           of length 5, corresponding to the vector in parameter
//                           space (\mu, B, m_1^2, m_2^2). For
//                           the one loop method epsVec must be of length 7, containing
//                           the above 4 values followed by (m_Q^2, m_U^2, A_t). For
//                           ordinary partial derivatives, epsVec should have one element
//                           equal to unity, and all others vanishing.
//     double epsilon = the step length to use in the finite difference approximation
//     double tol = the tolerance required for Newton's method to converge
//     int maxIters = the maximum number of allowed iterations in Newton's method
//     int & sing = an integer flag indicating if singular (in which case sing is non-zero)
DoubleMatrix pMSSM_vevNumericalDerivatives(SoftParsMssm,double,int,DoubleVector,double,double,int,int &);

// pMSSM_NumericalFineTunings calculates the fine tunings numerically using the above methods.
// Inputs:
//     SoftParsMssm r = pMSSM object to calculate fine tuning for
//     double tb = tan(beta)
//     int l = # loops to use (0 or 1)
//     DoubleVector & deltas = vector to store fine tunings
//     int & sing = flag if singular or poor accuracy (non-zero if problem)
//     double & detF = determinant of tuning matrix
//     double epsilon = the step length to use in the finite difference approximation
//     double tol = the tolerance required for Newton's method to converge
//     int maxIters = the maximum number of allowed iterations in Newton's method
void pMSSM_NumericalFineTunings(SoftParsMssm, double, int, DoubleVector &, int &, double &, double, double, int);

/*
  --------------------------------------------------------------------------------------
  End of numerical methods
  --------------------------------------------------------------------------------------
 */

// pMSSM_Ztunings calculates the sensitivities using the expressions given in
// arXiv:1206.5800. The tunings are returned in a double vector, which has elements
//     1 = Z_Ab^LL
//     2 = Z_Atau^LL
//     3 = Z_M1^LL
//     4 = Z_M2^LL
//     5 = Z_ML3^LL
//     6 = Z_Me3^LL
//     7 = Z_mu^LL
//     8 = Z_mAsq^TL
//     9 = Z_tb^TL
//    10 = Z_At^NLL
//    11 = Z_MQ3^NLL
//    12 = Z_Mu3^NLL
//    13 = Z_Md3^LL
//    14 = Z_M3^NLL
// Note that in the assumptions made in the pMSSM, 5 of the 19 parameters have
// zero or negligible contribution to fine tuning.
// Inputs:
//     SoftParsMssm r = pMSSM model
//     double tb = tan(beta)
//     double X = cut-off scale
DoubleVector pMSSM_Ztunings(SoftParsMssm, double, double);


#endif
