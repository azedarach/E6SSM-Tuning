/*
  tuningutils provides what are hopefully general methods for
  calculating the fine tuning in a SUSY model.
 */

#ifndef TUNINGUTILS_H
#define TUNINGUTILS_H

#include "tuningnumerics.h"

// The various doCalcFineTuning functions provide a set of overloaded
// function for calculating the fine tuning in a model. Hopefully they
// are general enough to apply to either the MSSM or the E6SSM without
// any substantial modifications when going between the two.

// The first version of doCalcFineTuning takes as inputs a model object, 
// a vector of parameters, a vector of VEVs, the scale at which the input parameters are defined, 
// an integer flag indicating if singular (non-zero if problem), 
// and three function references. The first function
// is assumed to provide as output the matrix appearing on the LHS of the 
// calculation of the derivatives dv_i/dp_j. As input it takes the SUSY model and the vector of VEVs.
// The second function should return the vector appearing on the RHS of the same calculation. 
// As input it takes the SUSY model, the parameters, the VEVs, the index of the parameter to calculate
// the vector for, and the scale at which the parameters are defined. 
// The third function should be the function used to compute d log(M_Z^2)/d log(p_j) from the derivatives
// of the VEVs. As input it takes the SUSY model, the parameter of interest, and vectors containing the VEVs and the
// derivatives of the VEVs wrt that parameter.
// Returns a vector containing the fine tuning wrt each of the
// given parameters, calculated according to \Delta = \frac{\partial \log M_Z^2}{\partial \log p_i}.

template <class SoftPars>
DoubleVector doCalcFineTuning(SoftPars susyModel, void (*BCatMX)(SoftPars &, DoubleVector &), 
			      DoubleVector const & pars, DoubleVector const & vevs, double mx, 
			      int & sing, 
			      DoubleMatrix (*lhsMatrix)(SoftPars, DoubleVector const &),
			      DoubleVector (*rhsVector)(SoftPars, void (*ftBC)(SoftPars &, DoubleVector &), 
							DoubleVector , DoubleVector const &, int, double, bool &),
			      double (*dlogMzSqdlogParam)(SoftPars, double, DoubleVector const &, DoubleVector const &))
{
  // First get the LHS matrix; this is the same for all parameters.
  DoubleMatrix lhs_mat = lhsMatrix(susyModel, vevs);
  DoubleMatrix Q_mat = lhs_mat;

  int n = lhs_mat.displayRows();

  // Create a vector to store the fine tunings.
  DoubleVector tunings = pars;
  
  // Also create a vector to store the derivatives of the VEVs.
  DoubleVector dVevs = vevs;
     
  if (n != lhs_mat.displayCols()) // matrix should be square
    {
      cerr << "WARNING: non-square matrix passed to fine tuning calculation: skipping calculating tunings." << endl;
      sing = 1;
    }
  else
    {
      QR_decomp(lhs_mat, Q_mat, n, sing);
      bool hasError = false;
      // Now loop over the given parameters and calculate the fine tuning for each.
      for (int i = pars.displayStart(); i <= pars.displayEnd(); i++)
  	{
	  hasError = false;

  	  QR_solve(Q_mat, lhs_mat, n, rhsVector(susyModel, BCatMX, pars, vevs, i, mx, hasError), dVevs, sing);

	  if (hasError)
	    {
	      cerr << "WARNING: calculated fine tuning for parameter " << i << " is inaccurate." << endl;
	      sing = 1;
	    }

  	  tunings(i) = fabs(dlogMzSqdlogParam(susyModel, pars(i), vevs, dVevs));

  	}
    }

  return tunings;
}

// The second version of doCalcFineTuning differs from the first only in that it takes an additional
// argument containing auxiliary data as a vector, which is assumed to be passed to the function calculating
// the right-hand side vector as its final argument. Used e.g. to give the values of the gauge and Yukawa
// couplings at MX (since the model is assumed to be provided at M_{SUSY}).
template <class SoftPars>
DoubleVector doCalcFineTuning(SoftPars susyModel, void (*BCatMX)(SoftPars &, DoubleVector &), 
			      DoubleVector const & pars, DoubleVector const & vevs, 
			      DoubleVector const & auxData, double mx, int & sing, 
			      DoubleMatrix (*lhsMatrix)(SoftPars, DoubleVector const &), 
			      DoubleVector (*rhsVector)(SoftPars, void (*ftBC)(SoftPars &, DoubleVector &), DoubleVector , 
							DoubleVector const &, int, double, bool &, DoubleVector const & ),
			      double (*dlogMzSqdlogParam)(SoftPars, double, DoubleVector const &, DoubleVector const &))
{
  // First get the LHS matrix; this is the same for all parameters.
  DoubleMatrix lhs_mat = lhsMatrix(susyModel, vevs);
  DoubleMatrix Q_mat = lhs_mat;

  int n = lhs_mat.displayRows();

  // Create a vector to store the fine tunings.
  DoubleVector tunings = pars;
  
  // Also create a vector to store the derivatives of the VEVs.
  DoubleVector dVevs = vevs;
     
  if (n != lhs_mat.displayCols()) // matrix should be square
    {
      cerr << "WARNING: non-square matrix passed to fine tuning calculation: skipping calculating tunings." << endl;
      sing = 1;
    }
  else
    {
      QR_decomp(lhs_mat, Q_mat, n, sing);
      bool hasError = false;
      // Now loop over the given parameters and calculate the fine tuning for each.
      for (int i = pars.displayStart(); i <= pars.displayEnd(); i++)
	{
	  QR_solve(Q_mat, lhs_mat, n, rhsVector(susyModel, BCatMX, pars, vevs, i, mx, hasError, auxData), dVevs, sing);

	  if (hasError)
	    {
	      cerr << "WARNING: calculated fine tunings are inaccurate." << endl;
	      sing = 1;
	    }

	  tunings(i) = fabs(dlogMzSqdlogParam(susyModel, pars(i), vevs, dVevs));

	}
    }

  return tunings;
}


#endif
