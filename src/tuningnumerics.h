/*
  tuningnumerics provides some useful numerical methods for calculating
  the fine tuning in a given model. Most are just easier to use
  wrapper functions for some of the numerical utilities in SOFTSUSY. The 
  exceptions are the functions for evaluating the fine tuning itself.
 */

#ifndef TUNINGNUMERICS_H
#define TUNINGNUMERICS_H

#include "E6SSM_Spectrum_Generators/src/linalg.h"
#include "E6SSM_Spectrum_Generators/src/utils.h"
#include "E6SSM_Spectrum_Generators/src/numerics.h"

// The various doCalcFineTuning functions provide a set of overloaded
// function for calculating the fine tuning in a model. Hopefully they
// are general enough to apply to either the MSSM or the E6SSM without
// any substantial modifications when going between the two.

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

// Helper function to find index of maximum element in a vector.
int findMaxIndex(DoubleVector const & vec);

/* // The following functions overload some of the numerical routines in numerics.h */
/* // to allow an extra vector of auxiliary data to be passed to define a function. */
/* // NB to make these work easily I have modified numerics.h, which means this will not */
/* // compile if used with out-of-the-box SOFTSUSY-3.3.10. The following small changes */
/* // are necessary: */
/* //   - a function prototype for ludcmp has to be added to numerics.h */
/* //   - the function prototype of fdjac has to have the return type changed from void */
/* //     to DoubleMatrix, and the argument 'DoubleMatrix & df' has to be deleted from it. */

/* bool lnsrch(const DoubleVector & xold, const DoubleVector & auxData, double fold, const DoubleVector & g,  */
/* 	    DoubleVector & p,  */
/* 	    DoubleVector & x, double & f, double stpmax,  */
/* 	    void (*vecfunc)(const DoubleVector &, const DoubleVector &, DoubleVector &),  */
/* 	    DoubleVector & fvec); */

/* bool newt(DoubleVector & x, const DoubleVector & auxData, */
/* 	  void (*vecfunc)(const DoubleVector &, const DoubleVector &, DoubleVector &)); */

/* DoubleMatrix fdjac(int n, DoubleVector x, const DoubleVector & auxData, const DoubleVector & fvec,  */
/* 	   void (*vecfunc)(const DoubleVector &, const DoubleVector &, DoubleVector &)); */


#endif
