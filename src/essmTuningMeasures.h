/*
  higgsNumericalTuningMeasures.h provides the function headers
  and prototypes for the functions that calculate the values of the
  E_6SSM EWSB conditions, with or without tadpole corrections,
  implement the EWSB constraints, and numerically calculate the 
  tuning measures (two possible methods?).
  Author: Dylan Harries
  Date created: 30/9/2013
 */

#ifndef ESSMTUNINGMEASURES_H
#define ESSMTUNINGMEASURES_H

#include "Higgs_E6SSM.h"
#include <math.h>
#include <iostream>

/*
  NOTE: in all functions involving one-loop quantities
  there is the option to fix m_top = 165 GeV at
  the renormalisation scale (usually assumed to be M_t).
  This is currently turned ON.
 */

// EWSBConditioni_TreeLevel, i = 1, 2, 3, calculates the value
// of the EWSB condition for m_H_i^2 (i=3 is m_s^2) at tree level (no tadpole
// corrections). Provided to allow for comparison with previous codes done
// at tree level.
// Inputs:
//     SoftParsEssm essmSusy = object containing the ESSM parameters
//     double s = value of the VEV s
//     double tb = value of tan(beta)
double EWSBCondition1_TreeLevel(SoftParsEssm,double,double);
double EWSBCondition2_TreeLevel(SoftParsEssm,double,double);
double EWSBCondition3_TreeLevel(SoftParsEssm,double,double);


// EWSBConditioni_Loops, i = 1, 2, 3, calculates the value
// of the EWSB condition for m_H_i^2 (i=3 is m_s^2), including the 
// contributions from the tadpole corrections. Evaluates the tadpoles
// at the renormalisation scale chosen in the model, and includes U(1)
// D-terms.
// Inputs:
//     SoftParsEssm essmSusy = object containing the ESSM parameters
//     double s = value of the VEV s
//     double tb = value of tan(beta)
double EWSBCondition1_OneLoop(SoftParsEssm,double,double);
double EWSBCondition2_OneLoop(SoftParsEssm,double,double);
double EWSBCondition3_OneLoop(SoftParsEssm,double,double);


// As above, but evaluates the tadpoles at the scale of M_t, and includes U(1)
// D-terms.
// Inputs:
//     SoftParsEssm essmSusy = object containing the ESSM parameters
//     double s = value of the VEV s
//     double tb = value of tan(beta)
double EWSBCondition1_OneLoop_atMt(SoftParsEssm,double,double);
double EWSBCondition2_OneLoop_atMt(SoftParsEssm,double,double);
double EWSBCondition3_OneLoop_atMt(SoftParsEssm,double,double);

// As above, but evaluates the tadpoles at the scale of M_t, and neglects
// U(1) D-terms.
// Inputs:
//     SoftParsEssm essmSusy = object containing the ESSM parameters
//     double s = value of the VEV s
//     double tb = value of tan(beta)
double EWSBCondition1_OneLoop_Roman(SoftParsEssm,double,double);
double EWSBCondition2_OneLoop_Roman(SoftParsEssm,double,double);
double EWSBCondition3_OneLoop_Roman(SoftParsEssm,double,double);

// As above, but evaluates the tadpoles at the scale Q defined in the model, and 
// neglects U(1) D-terms.
// Inputs:
//     SoftParsEssm essmSusy = object containing the ESSM parameters
//     double s = value of the VEV s
//     double tb = value of tan(beta)
double EWSBCondition1_OneLoop_Roman_atQ(SoftParsEssm,double,double);
double EWSBCondition2_OneLoop_Roman_atQ(SoftParsEssm,double,double);
double EWSBCondition3_OneLoop_Roman_atQ(SoftParsEssm,double,double);

// ImplementEWSBConstraints_TreeLevel imposes the EWSB constraints
// on the given ESSM model, setting the value of the soft masses
// m_H_1^2, m_H_2^2 and m_s^2 to satisfy the three EWSB conditions
// at tree level. Allows for comparison with previous codes done at
// tree level.
// Inputs:
//     SoftParsEssm & essmSusy = the object to impose the EWSB constraints on.
void ImplementEWSBConstraints_TreeLevel(SoftParsEssm &,double,double);

// ImplementEWSBConstraints_Loops imposes the EWSB constraints
// on the given ESSM model, setting the value of the soft masses
// m_H_1^2, m_H_2^2 and m_s^2 to satisfy the three EWSB conditions
// after including contributions from tadpole diagrams.
// Inputs:
//     SoftParsEssm & essmSusy = the object to impose the EWSB constraints on.
void ImplementEWSBConstraints_OneLoop(SoftParsEssm &,double,double);

// As above, except the tadpole corrections used include contributions from the D fields 
// and are evaluated at the renormalisation scale M_t.
// Inputs:
//     SoftParsEssm & essmSusy = the object to impose the EWSB constraints on.
//     double s = the value of the VEV s
//     double tb = the value of tan(beta)
void ImplementEWSBConstraints_OneLoop_atMt(SoftParsEssm &,double,double);

// As above, except the tadpole corrections used don't include contributions from the D fields 
// and are evaluated at the renormalisation scale M_t.
// Inputs:
//     SoftParsEssm & essmSusy = the object to impose the EWSB constraints on.
//     double s = the value of the VEV s
//     double tb = the value of tan(beta)
void ImplementEWSBConstraints_OneLoop_Roman(SoftParsEssm &,double,double);

// As above, except the tadpole corrections used don't include contributions from the D fields 
// and are evaluated at the renormalisation scale Q.
// Inputs:
//     SoftParsEssm & essmSusy = the object to impose the EWSB constraints on.
//     double s = the value of the VEV s
//     double tb = the value of tan(beta)
void ImplementEWSBConstraints_OneLoop_Roman_atQ(SoftParsEssm &,double,double);

// FineTunings_TreeLevel calculates the fine-tuning sensitivities \Delta_a for each 
// a = \lambda, A_\lambda, m_1^2, m_2^2, m_3^2, m_Q^2, m_U^2, A_t at tree-level,
// that is, neglecting the tadpole corrections to the EWSB conditions,
// at the given point in parameter space. 
// Note that at tree-level m_Q^2, m_U^2 and A_t do not appear, so the
// associated sensitivities vanish. 
// Inputs:
//     SoftParsEssm const & essmSusy = the object to evaluate the fine tuning at
//     double s = the value of the VEV s
//     double tb = the value of tan(beta) to use
//     DoubleVector & deltas = a vector (of length 8) to store the calculated
//     sensitivities, in the order \lambda, A_\lambda, m_1^2, m_2^2, m_3^2,
//     m_Q^2, m_U^2, A_t.
//     int & sing = integer flag to indicate if QR solution is singular
void FineTunings_TreeLevel(SoftParsEssm const &, double, double, DoubleVector &, int &, double &); 

// As above, but use Cramer's rule to obtain the solution.
// Additional input double & detF = determinant of tunings coefficients matrix
void FineTunings_TreeLevel_CR(SoftParsEssm const &, double, double, DoubleVector &, int &, double &);

// Methods to generate the matrix of coefficients of the 
// derivatives necessary to evaluate the tuning measures,
// and vectors of values for the right hand side of the tuning
// measure equations.
// Inputs:
//     SoftParsEssm const & essmSusy = object to calculate coefficients matrix from.
//     double s = the value of the VEV s to use
//     double tb = the value of tan(beta) to use
DoubleMatrix TuningsCoefficientsMatrix_TreeLevel(SoftParsEssm const &,double,double);


// Note: have reworked tuning measures to be evaluated in terms 
// of derivatives wrt v_1, v_2 and s instead of M_Z. 
DoubleVector TuningsRHS_Lambda_TreeLevel(SoftParsEssm const &,double,double);
DoubleVector TuningsRHS_Alambda_TreeLevel(SoftParsEssm const &,double,double);
DoubleVector TuningsRHS_Mh1Sq_TreeLevel(SoftParsEssm const &,double,double);
DoubleVector TuningsRHS_Mh2Sq_TreeLevel(SoftParsEssm const &,double,double);
DoubleVector TuningsRHS_MsSq_TreeLevel(SoftParsEssm const &,double,double);
DoubleVector TuningsRHS_MQlSq_TreeLevel(SoftParsEssm const &,double,double);
DoubleVector TuningsRHS_MUrSq_TreeLevel(SoftParsEssm const &,double,double);
DoubleVector TuningsRHS_At_TreeLevel(SoftParsEssm const &,double,double);

// As above, but including the leading one-loop contributions from
// top and stop loops.
void FineTunings_OneLoop(SoftParsEssm const &, double, double, DoubleVector &, int &, double &); 

void FineTunings_OneLoop_CR(SoftParsEssm const &, double, double, DoubleVector &, int &, double &);

DoubleMatrix TuningsCoefficientsMatrix_OneLoop(SoftParsEssm const &,double,double);

DoubleVector TuningsRHS_Lambda_OneLoop(SoftParsEssm const &,double,double);
DoubleVector TuningsRHS_Alambda_OneLoop(SoftParsEssm const &,double,double);
DoubleVector TuningsRHS_Mh1Sq_OneLoop(SoftParsEssm const &,double,double);
DoubleVector TuningsRHS_Mh2Sq_OneLoop(SoftParsEssm const &,double,double);
DoubleVector TuningsRHS_MsSq_OneLoop(SoftParsEssm const &,double,double);
DoubleVector TuningsRHS_MQlSq_OneLoop(SoftParsEssm const &,double,double);
DoubleVector TuningsRHS_MUrSq_OneLoop(SoftParsEssm const &,double,double);
DoubleVector TuningsRHS_At_OneLoop(SoftParsEssm const &,double,double);

// The following is a final fine tuning method that should incorporate the
// functionality of all of the previously defined methods above. It uses
// the version of the EWSB conditions in which the VEVs are not divided out
// to start with (i.e. unrealistic vacua are permitted).
// Inputs:
//     SoftParsEssm const & essmSusy = E6SSM model object, assumed to contain correct
//     values of parameters appearing in the Higgs scalar potential (e.g. lambda, 
//     A_lambda, and the gauge couplings).
//     double s = the value of s to use
//     double tb = the value of tan(beta) to use
//     int l = # loops to use (0 or 1)
//     DoubleMatrix const & partials = a Nx4 matrix containing the partial derivatives of the
//     EWSB conditions wrt the N input parameters. The element (i, 1) contains the value of
//     the parameter p_i. The (i, j)th, j>1,  element is taken to be
//     df_(j-1)/dp_i. These should be calculated consistently with the order requested when
//     this function is called, i.e. tadpole corrections should be included if l = 1.
//     int & sing = flag to indicate if singular (non-zero if problem)
//     double & detF = determinant of CP-even Higgs mass matrix
//     
//     Returns a DoubleVector of length N containing the fine tuning sensitivities. The ith
//     element of the vector is the tuning corresponding to the parameter in the ith row
//     of the input matrix.
DoubleVector E6SSM_FineTunings(SoftParsEssm const &, double, double, int, DoubleMatrix  const &, 
			       int &, double &);

// This function produces the partials matrix for the specific parameter set 
// {lambda, A_lambda, m_H1^2, m_H2^2, m_S^2, m_Q^2, m_U^2, A_t} with NO RG RUNNING.
// Inputs:
//     SoftParsEssm const & essmSusy = E6SSM model object to calculate tuning for
//     double s = the value of s to use
//     double tb = the value of tan(beta) to use
//     int l = # loops to use (0 or 1)
//     DoubleMatrix & partials = a matrix to store the resulting derivatives
void calculatePartialsMatrix(SoftParsEssm const &, double, double, int, DoubleMatrix &);

// The following are helper functions that are useful for constructing
// the leading one-loop tuning measures.

// dmstop1dvi calculates the derivative of m_stop_1^2 with respect to
// v_i for the given ESSM model, i = 1, 2, 3 (i=3 is s). 
// NB: m_stop_1^2 is taken to be the lighter stop in this method 
// -i.e. minus sign applies.
// Inputs:
//     SoftParsEssm essmSusy = object to calculate derivatives on
//     double s = the value of the VEV s to use
//     double tb = the value of tan(beta) to use
double doCalcdmstop1sqdv1(SoftParsEssm,double,double);
double doCalcdmstop1sqdv2(SoftParsEssm,double,double);
double doCalcdmstop1sqdv3(SoftParsEssm,double,double);

// dmstop2dvi is as above but for m_stop_2^2.
double doCalcdmstop2sqdv1(SoftParsEssm,double,double);
double doCalcdmstop2sqdv2(SoftParsEssm,double,double);
double doCalcdmstop2sqdv3(SoftParsEssm,double,double);

// The three functions below are purely used for testing the first
// derivatives of the stop squared masses. They calculate the 
// stop/top tadpole correction to the effective potential, in the
// form (1/v_i)\frac{\partial \Delta V}{\partial v_i}. The results
// should be compared against those obtained using the functions
// doCalcTadpolesESSMH<1,2,S> which are actually used in practical
// calculations. In the notation used in my notes, these functions
// calculate \Delta^{tadpole}_i for i = 1, 2, 3.
// Inputs:
//     SoftParsEssm essmSusy = object to calculate derivatives on
//     double s = the value of the VEV s to use
//     double tb = the value of tan(beta) to use
double calculateDeltaTadpole1(SoftParsEssm,double,double);
double calculateDeltaTadpole2(SoftParsEssm,double,double);
double calculateDeltaTadpole3(SoftParsEssm,double,double);

// These functions calculate the leading one-loop corrections
// to the CP-even Higgs mass matrix in the basis (v_1, v_2, v_3=s),
// given by \Delta_{ij}'=\frac{\partial^2\Delta V}{\partial v_i\partial v_j}.
// Inputs:
//     SoftParsEssm essmSusy = object to calculate derivatives on
//     double s = the value of the VEV s to use
//     double tb = the value of tan(beta) to use.
// CHECK IF NEED TO SUBTRACT TADPOLES FROM THESE!!
double doCalcDeltaPrime11(SoftParsEssm,double,double);
double doCalcDeltaPrime12(SoftParsEssm,double,double);
double doCalcDeltaPrime13(SoftParsEssm,double,double);
double doCalcDeltaPrime22(SoftParsEssm,double,double);
double doCalcDeltaPrime23(SoftParsEssm,double,double);
double doCalcDeltaPrime33(SoftParsEssm,double,double);

// These helper functions calculate the explicit derivatives
// of the tadpole contributions wrt the input parameters, i.e. returns
// the partial derivatives (1/v_i)*\frac{\partial^2 \Delta V}{\partial a\partial v_i}
// where the VEVs v_i, i = 1, 2, 3, are treated as being fixed.
// Inputs:
//     SoftParsEssm essmSusy = the object to calculate the derivatives on
//     double s = the value of the VEV s to use
//     double tb = the value of tan(beta) to use
double doCalcd2DeltaVdLambdadv1(SoftParsEssm,double,double);
double doCalcd2DeltaVdLambdadv2(SoftParsEssm,double,double);
double doCalcd2DeltaVdLambdadv3(SoftParsEssm,double,double);
double doCalcd2DeltaVdAtdv1(SoftParsEssm,double,double);
double doCalcd2DeltaVdAtdv2(SoftParsEssm,double,double);
double doCalcd2DeltaVdAtdv3(SoftParsEssm,double,double);
double doCalcd2DeltaVdmQlsqdv1(SoftParsEssm,double,double);
double doCalcd2DeltaVdmQlsqdv2(SoftParsEssm,double,double);
double doCalcd2DeltaVdmQlsqdv3(SoftParsEssm,double,double);
double doCalcd2DeltaVdmUrsqdv1(SoftParsEssm,double,double);
double doCalcd2DeltaVdmUrsqdv2(SoftParsEssm,double,double);
double doCalcd2DeltaVdmUrsqdv3(SoftParsEssm,double,double);

/*
  --------------------------------------------------------------------------------------
  The following are methods for evaluating the derivatives of the VEVs and the tuning
  measures numerically. Wherever possible DO NOT use them, as to work they require
  turning off the use m_t(M_t) option in the Higgs codes, which may lead to bugs
  and inconsistencies. If they must be used, note also that the EWSB conditions should
  be changed to those in which the VEVs are divided out (assumes realisitic vacua in which
  all VEVs are non-vanishing). This improves the convergence of Newton's method.
  --------------------------------------------------------------------------------------
 */


// The following methods provide functionality to numerically differentiate
// the EWSB conditions at one loop order, using Newton's method.

// EWSBJacobian_TreeLevel and _OneLoop compute the Jacobian matrix of
// the system defined by the three EWSB conditions, with elements
// J_{ij}= \partial f_i/\partial v_j (v_3 = s). The functions f_i
// are the EWSB conditions, with f_i=\partial V/\partial v_i. The
// two versions differ in the effective potential V used; one loop
// includes the leading one loop contributions from stop and top loops.
// Inputs:
//     SoftParsEssm essmSusy = the object to calculate the Jacobian for
//     double s = the value of the VEV s to use
//     double tb = the value of tan(beta) to use 
DoubleMatrix EWSBJacobian_TreeLevel(SoftParsEssm,double,double);
DoubleMatrix EWSBJacobian_OneLoop(SoftParsEssm,double,double);

// qrUpdateEWSBSoln_TreeLevel and _OneLoop use a QR decomposition
// to compute the update step in Newton's method, to solve for
// the VEVs at the parameter point contained in the ESSM object.
// Returns the updated estimate for the solution to the EWSB 
// conditions. Note that the provided object is assumed to 
// contain the initial (current) guess for the solutions to
// the EWSB conditions.
// Inputs:
//     SoftParsEssm essmSusy = the object to calculate the solution for
//     double s = the value of the VEV s to use
//     double tb = the value of tan(beta) to use
//     int & sing = integer flag to indicate if singular (in which case sing is non-zero)
DoubleVector qrUpdateEWSBSoln_TreeLevel(SoftParsEssm, double, double, int &);
DoubleVector qrUpdateEWSBSoln_OneLoop(SoftParsEssm, double, double, int &);

// EWSBNewtonSolver_TreeLevel and _OneLoop computes the solution to the 
// EWSB conditions at either tree level or one loop order using Newton's method.
// Returns the solution for the VEVs v_1, v_2 and s as a DoubleVector.
// Note that the provided object is taken to contain both the values of
// the input parameters and the initial guess for the VEVs (with the exception
// of s, of course). 
// Inputs:
//     SoftParsEssm essmSusy = the object to calculate the solution for
//     double s = the initial guess for the value of the VEV s
//     double tb = the initial guess for the value of tan(beta)
//     double tol = the required tolerance to test for convergence
//     int maxIters = the maximum number of allowed iterations
//     int & sing = integer flag to indicate if singular (in which case sing is non-zero)
DoubleVector EWSBNewtonSolver_TreeLevel(SoftParsEssm, double, double, double, int, int &);
DoubleVector EWSBNewtonSolver_OneLoop(SoftParsEssm, double, double, double, int, int &);

// vevNumericalDerivatives_TreeLevel and _OneLoop compute estimates for the partial
// derivatives of the Higgs and singlet VEVs wrt the low energy input parameters
// \lambda, A_\lambda,... at either tree level or one loop order. This is done
// using Newton's method and a 5-pt stencil approximation to the derivative. Note that
// the result calculated is the directional derivative in the direction specified by
// epsVec. The parameter point at which the derivative is to be evaluated, together
// with the initial estimate for the VEVs, is assumed to be contained in the provided
// ESSM model object.
// Inputs:
//     SoftParsEssm essmSusy = the object to calculate the derivatives for
//     double s = the value of the VEV s to use
//     double tb = the value of tan(beta) to use
//     DoubleVector epsVec = a vector indicating the direction in the low-energy
//                           parameter space in which to calculate the directional
//                           derivative. For the tree level method, epsVec must be
//                           of length 5, corresponding to the vector in parameter
//                           space (\lambda, A_\lambda, m_1^2, m_2^2, m_S^2). For
//                           the one loop method epsVec must be of length 8, containing
//                           the above 5 values followed by (m_Q^2, m_U^2, A_t). For
//                           ordinary partial derivatives, epsVec should have one element
//                           equal to unity, and all others vanishing.
//     double epsilon = the step length to use in the finite difference approximation
//     double tol = the tolerance required for Newton's method to converge
//     int maxIters = the maximum number of allowed iterations in Newton's method
//     int & sing = an integer flag indicating if singular (in which case sing is non-zero)
DoubleVector vevNumericalDerivatives_TreeLevel(SoftParsEssm,double,double,DoubleVector,double,double,int,int &);
DoubleVector vevNumericalDerivatives_OneLoop(SoftParsEssm,double,double,DoubleVector,double,double,int,int &);

/*
  --------------------------------------------------------------------------------------
*/


// A method to generate the matrices Q and R in the QR decomposition
// A = QR of the n x n square matrix A. On output the original input matrix A is
// replaced by the upper triangular matrix R (the original matrix
// can always be recovered as A = QR if necessary). Not exactly optimised,
// but should be sufficient for now when only working with 3 x 3 matrices.
// Inputs:
//     DoubleMatrix & A = the matrix to get the QR decomposition of
//     DoubleMatrix & Q = matrix storing orthogonal matrix Q
//     int n = dimension of matrix A
//     int & sing = flag indicating if singularity occurred.
void QR_decomp(DoubleMatrix &,DoubleMatrix &,int,int &);


// A method to solve a linear system by using the QR decomposistion to
// do back-substitution. Requires the matrices Q and R to be supplied.
// Inputs:
//     DoubleMatrix & Q = the orthogonal matrix Q in the QR decomposition
//     DoubleMatrix & R = the upper triangular matrix R
//     int n = the dimension of the matrices Q and R (square)
//     DoubleVector const & b = the right hand side of the linear system
//     DoubleVector & x = a vector to store the solution in
//     int & sing = flag indicating if singularity occurred.
void QR_solve(DoubleMatrix const &,DoubleMatrix const &,int,DoubleVector const &,DoubleVector &,int &);

// A method to solve a 3x3 linear system using Cramer's rule.
// Inputs:
//     DoubleMatrix const & A = the matrix of coefficients
//     DoubleVector const & = the vector on the RHS
//     double & detA = the calculated value of the determinant of A
//     int & sing = integer flag to indicate if no unique soln (non-zero if singular)
//     DoubleVector & soln = a vector of length 3 containing the solution
void CramersRule(DoubleMatrix const &, DoubleVector const &, double &, int &, DoubleVector &);


#endif
