/*
  essmtuningutils contains all of the functions needed
  for calculating the fine tuning in the (p)E6SSM.
 */

#ifndef ESSMTUNINGUTILS_H
#define ESSMTUNINGUTILS_H


#include "E6SSM_Spectrum_Generators/models/genericE6SSM/genericE6SSM_two_scale_soft_parameters.hpp"
#include "E6SSM_Spectrum_Generators/src/wrappers.hpp"
#include "flags.h"
#include "tuningnumerics.h"
#include "tuningutils.h"

// To make including more parameters easier later and to avoid name clashes
// (may be even better for encapsulation later on if we define a class
// with a model as a member variable):
namespace essm_tuning_utils {

  enum tuning_parameters : unsigned { lam3, M1, M2, M3, M1p, 
	       Au3, Alam3, mH13Sq, mH23Sq, mS3Sq, mqL3Sq, mtRSq, NUMESSMTUNINGPARS };

// A function to calculate the 2x2 matrix appearing on the LHS of the calculation
// of the derivatives of the VEVs wrt the parameters in the MSSM. Note that this
// assumes that the SoftParsMssm object is already set to have its renormalisation
// scale set at M_SUSY.
 Eigen::Matrix<double,3,3> doCalcLHSTuningMatrix(flexiblesusy::genericE6SSM_soft_parameters essmSusy, Eigen::VectorXd const & vevs);

//double rgeDerivCalc(double x);

// A function to calculate the vectors appearing on the RHS of the calculation of the
// derivatives d v_i/ dp_j. Takes as arguments a SoftParsMssm object, assumed to be evaluated
// at the scale M_SUSY, vectors of the parameters to calculate the derivatives wrt. and the vevs,
// an integer labelling which particular parameter to calculate the derivative for, and the scale
// at which that parameter is defined.
/* DoubleVector doCalcRHSTuningVector(genericE6SSM_soft_parameters r, void (*ftBCatMX)(genericE6SSM_soft_parameters &, DoubleVector &),  */
/* 				   DoubleVector pars, DoubleVector const & vevs,  */
/* 				   int i, double mx, bool & hasError); */

// A function to calculate d log(M_Z^2)/d log(p) using the value of the parameter p
// and the already calculated derivatives d v_i/ d p_j. 
 double doCalcdLogMzSqdLogParam(flexiblesusy::genericE6SSM_soft_parameters r, double p, Eigen::VectorXd const & vevs, 
				Eigen::VectorXd const & dVevsdp);

// ESSM_EWSBConditioni, i = 1, 2, calculates the value of the 
// EWSB condition for m_H_i^2 and m_S^2.
// Inputs:
//     SoftParsMssm const & r = pMSSM model object
//     DoubleVector const & vevs = the values of v1 and v2
double ESSM_EWSBCondition1(flexiblesusy::genericE6SSM_soft_parameters const & r);
double ESSM_EWSBCondition2(flexiblesusy::genericE6SSM_soft_parameters const & r);
double ESSM_EWSBCondition3(flexiblesusy::genericE6SSM_soft_parameters const & r);

// Implement the EWSB conditions by varying the values of m_Hu^2 and m_Hd^2 at MX.
bool ESSM_ImplementEWSBConstraints_SoftMasses(flexiblesusy::genericE6SSM_soft_parameters r, double mx, double ms, 
					      bool useMxEqualsMs,
					      DoubleVector & updatedSoln, double tol = TOLERANCE);


bool ESSM_EWSB_NewtonShooter(flexiblesusy::genericE6SSM_soft_parameters const & r, DoubleVector & estimate, 
			     double tol = TOLERANCE);

double ESSM_Msusy_Cond(flexiblesusy::genericE6SSM_soft_parameters r, double ms);

void ESSM_EWSB_Shooter_Functions(DoubleVector const & parVals, DoubleVector & f);

double ccbSqrt(double);

// Function a0 defined in paper. Implement appropriately.
double a0Peter(double, double);

// Calculates the masses of the stops and exotic quarks for tadpoles.  The former are always included 
// in our tadpole contributions, the latter could be, but provide small contributions in comparison to 
// errors due to threshold treatment.
// Inputs:
//    flexiblesusy::genericE6SSM_soft_parameters r = ESSM model
//    DoubleVector mstop = mass of stops
//    DoubleVector mstopsq = squared masses of stops
//    DoubleVector mD1sq = squared masses of exotic D quarks
//    DoubleVector mD2sq = squared masses of exotic D quarks
//    double s = singlet VEV
//    double tb = tan(beta)
void physical_ESSM(flexiblesusy::genericE6SSM_soft_parameters,DoubleVector &,DoubleVector &,DoubleVector &,DoubleVector &,double,double);

//This version neglects U(1) D-terms and was written for comparison with Roman's program
void physical_ESSM_Roman(flexiblesusy::genericE6SSM_soft_parameters,DoubleVector &,DoubleVector &,DoubleVector &,DoubleVector &,double,double);

// Determines the dominant H1 tadpoles. There are a number of different versions of this depending on choice of
// scale. The difference  should be higher order. Exotics can also be included in this version though by default
// they switched off at the moment.   This version  includes U(1) D-terms and exotics, and evaluates the 
// scale where RG evolution is halted.
// Inputs:
//    flexiblesusy::genericE6SSM_soft_parameters r = ESSM model
//    double s = singlet VEV
//    double tb = tan(beta)
double doCalcTadpoleESSMH1(flexiblesusy::genericE6SSM_soft_parameters,double,double);

// Calculates tadpoles using logarithms at m_t. The assumption is coefficients do not substantially change 
// between scale at which RG evolution is halted and m_t. 
double doCalcTadpoleESSMH1_atMt(flexiblesusy::genericE6SSM_soft_parameters,double,double);

// This version neglects U(1) D-terms and calculates at m_t. 
double doCalcTadpoleESSMH1_Roman(flexiblesusy::genericE6SSM_soft_parameters,double,double);

// Yet another option, performs calculation like Roman, but uses logarithms at scale at which RGE evolution is halted.
double doCalcTadpoleESSMH1_Roman_atQ(flexiblesusy::genericE6SSM_soft_parameters,double,double);

// Now tadpoles for H2. The different versions are labled as for H1 tadpoles and contributions are matched with those.
double doCalcTadpoleESSMH2(flexiblesusy::genericE6SSM_soft_parameters,double,double);

double doCalcTadpoleESSMH2_atMt(flexiblesusy::genericE6SSM_soft_parameters,double,double);

double doCalcTadpoleESSMH2_Roman(flexiblesusy::genericE6SSM_soft_parameters,double,double);

double doCalcTadpoleESSMH2_Roman_atQ(flexiblesusy::genericE6SSM_soft_parameters,double,double);

// Tadpoles for S. Same labeling as above.
double doCalcTadpolesESSMS(flexiblesusy::genericE6SSM_soft_parameters,double,double);

double doCalcTadpolesESSMS_atMt(flexiblesusy::genericE6SSM_soft_parameters,double,double);

double doCalcTadpolesESSMS_Roman(flexiblesusy::genericE6SSM_soft_parameters,double,double);

double doCalcTadpolesESSMS_Roman_atQ(flexiblesusy::genericE6SSM_soft_parameters,double,double);

// Calculate the Higgs masses in the model. Returns true if tachyons problem.
// Inputs:
//    flexiblesusy::genericE6SSM_soft_parameters & r = ESSM model
//    double s = singlet VEV
//    double tb = tan(beta)
//    DoubleVector & mstop = stop masses
//    DoubleVector & mstopsq = squared stop masses
//    int WhatCorrections = corrections to use in calculation (see top of file)
//    bool speak = Boolean variable for trace printing
//    bool Bugspeak = Boolean variable for trace printing
//    DoubleVector & bounds = a vector containing lower (first element) and upper bound for valid Higgs masses
//    int & ExpValid = flag if Higgs masses are experimentally valid or not (0 if valid)
//    DoubleVector & mhout = Higgs masses
//    DoubleMatrix & mhmix = mixing matrix for Higgs fields
//    int & sing = flag to indicate poor accuracy in diagonalisation (non-zero if tolerance not met)
bool HiggsMasses(flexiblesusy::genericE6SSM_soft_parameters &,double,double,DoubleVector &,DoubleVector &,int,bool,bool,DoubleVector &,int &,DoubleVector &, DoubleMatrix &, DoubleMatrix &, int &);

// To the above we will add methods for calculating m_A^2 at tree level and one-loop order, and
// methods for calculating chargino and neutralino masses (at tree level).

// Calculate m_A^2 at tree level.
// Inputs:
//    flexiblesusy::genericE6SSM_soft_parameters const & essmSusy = the ESSM model to calculate m_A^2 for
//    double s = the value of the singlet VEV to use
//    double tb = the value of tan(beta) to use
double mAsq_TreeLevel(flexiblesusy::genericE6SSM_soft_parameters const &, double, double);

// Calculate m_A^2 at one loop order.
// Inputs:
//    flexiblesusy::genericE6SSM_soft_parameters essmSusy = the ESSM model to calculate m_A^2 for
//    double s = the value of the singlet VEV to use
//    double tb = the value of tan(beta) to use
double mAsq_OneLoop(flexiblesusy::genericE6SSM_soft_parameters, double, double, int &);

// A helper function used in calculating the mass m_A^2 at one loop order.
// Inputs:
//    double m1Sq = first squared mass m_1^2 to use
//    double m2Sq = second squared mass m_2^2 to use
//    double q = renormalisation scale
double f(double,double,double);

// The following are helper functions that are useful for constructing
// the leading one-loop tuning measures.

// dmstop1dvi calculates the derivative of m_stop_1^2 with respect to
// v_i for the given ESSM model, i = 1, 2, 3 (i=3 is s). 
// NB: m_stop_1^2 is taken to be the lighter stop in this method 
// -i.e. minus sign applies.
// Inputs:
//     flexiblesusy::genericE6SSM_soft_parameters essmSusy = object to calculate derivatives on
//     double s = the value of the VEV s to use
//     double tb = the value of tan(beta) to use
double doCalcdmstop1sqdv1(flexiblesusy::genericE6SSM_soft_parameters,double,double);
double doCalcdmstop1sqdv2(flexiblesusy::genericE6SSM_soft_parameters,double,double);
double doCalcdmstop1sqdv3(flexiblesusy::genericE6SSM_soft_parameters,double,double);

// dmstop2dvi is as above but for m_stop_2^2.
double doCalcdmstop2sqdv1(flexiblesusy::genericE6SSM_soft_parameters,double,double);
double doCalcdmstop2sqdv2(flexiblesusy::genericE6SSM_soft_parameters,double,double);
double doCalcdmstop2sqdv3(flexiblesusy::genericE6SSM_soft_parameters,double,double);

// The three functions below are purely used for testing the first
// derivatives of the stop squared masses. They calculate the 
// stop/top tadpole correction to the effective potential, in the
// form (1/v_i)\frac{\partial \Delta V}{\partial v_i}. The results
// should be compared against those obtained using the functions
// doCalcTadpolesESSMH<1,2,S> which are actually used in practical
// calculations. In the notation used in my notes, these functions
// calculate \Delta^{tadpole}_i for i = 1, 2, 3.
// Inputs:
//     flexiblesusy::genericE6SSM_soft_parameters essmSusy = object to calculate derivatives on
//     double s = the value of the VEV s to use
//     double tb = the value of tan(beta) to use
double calculateDeltaTadpole1(flexiblesusy::genericE6SSM_soft_parameters,double,double);
double calculateDeltaTadpole2(flexiblesusy::genericE6SSM_soft_parameters,double,double);
double calculateDeltaTadpole3(flexiblesusy::genericE6SSM_soft_parameters,double,double);

// These functions calculate the leading one-loop corrections
// to the CP-even Higgs mass matrix in the basis (v_1, v_2, v_3=s),
// given by \Delta_{ij}'=\frac{\partial^2\Delta V}{\partial v_i\partial v_j}.
// Inputs:
//     flexiblesusy::genericE6SSM_soft_parameters essmSusy = object to calculate derivatives on
//     double s = the value of the VEV s to use
//     double tb = the value of tan(beta) to use.
// CHECK IF NEED TO SUBTRACT TADPOLES FROM THESE!!
double doCalcDeltaPrime11(flexiblesusy::genericE6SSM_soft_parameters,double,double);
double doCalcDeltaPrime12(flexiblesusy::genericE6SSM_soft_parameters,double,double);
double doCalcDeltaPrime13(flexiblesusy::genericE6SSM_soft_parameters,double,double);
double doCalcDeltaPrime22(flexiblesusy::genericE6SSM_soft_parameters,double,double);
double doCalcDeltaPrime23(flexiblesusy::genericE6SSM_soft_parameters,double,double);
double doCalcDeltaPrime33(flexiblesusy::genericE6SSM_soft_parameters,double,double);

// These helper functions calculate the explicit derivatives
// of the tadpole contributions wrt the input parameters, i.e. returns
// the partial derivatives (1/v_i)*\frac{\partial^2 \Delta V}{\partial a\partial v_i}
// where the VEVs v_i, i = 1, 2, 3, are treated as being fixed.
// Inputs:
//     flexiblesusy::genericE6SSM_soft_parameters essmSusy = the object to calculate the derivatives on
//     double s = the value of the VEV s to use
//     double tb = the value of tan(beta) to use
double doCalcd2DeltaVdLambdadv1(flexiblesusy::genericE6SSM_soft_parameters,double,double);
double doCalcd2DeltaVdLambdadv2(flexiblesusy::genericE6SSM_soft_parameters,double,double);
double doCalcd2DeltaVdLambdadv3(flexiblesusy::genericE6SSM_soft_parameters,double,double);
double doCalcd2DeltaVdAtdv1(flexiblesusy::genericE6SSM_soft_parameters,double,double);
double doCalcd2DeltaVdAtdv2(flexiblesusy::genericE6SSM_soft_parameters,double,double);
double doCalcd2DeltaVdAtdv3(flexiblesusy::genericE6SSM_soft_parameters,double,double);
double doCalcd2DeltaVdmQlsqdv1(flexiblesusy::genericE6SSM_soft_parameters,double,double);
double doCalcd2DeltaVdmQlsqdv2(flexiblesusy::genericE6SSM_soft_parameters,double,double);
double doCalcd2DeltaVdmQlsqdv3(flexiblesusy::genericE6SSM_soft_parameters,double,double);
double doCalcd2DeltaVdmUrsqdv1(flexiblesusy::genericE6SSM_soft_parameters,double,double);
double doCalcd2DeltaVdmUrsqdv2(flexiblesusy::genericE6SSM_soft_parameters,double,double);
double doCalcd2DeltaVdmUrsqdv3(flexiblesusy::genericE6SSM_soft_parameters,double,double);

double doCalcMh2SquaredLogCoeff(flexiblesusy::genericE6SSM_soft_parameters, int);
double doCalcMh2SquaredLogSqCoeff(flexiblesusy::genericE6SSM_soft_parameters, int);

double doCalcMh1SquaredLogCoeff(flexiblesusy::genericE6SSM_soft_parameters , int);
double doCalcMh1SquaredLogSqCoeff(flexiblesusy::genericE6SSM_soft_parameters, int);

double doCalcMsSquaredLogCoeff(flexiblesusy::genericE6SSM_soft_parameters , int);
double doCalcMsSquaredLogSqCoeff(flexiblesusy::genericE6SSM_soft_parameters, int);

double doCalcmtRSquaredLogCoeff(flexiblesusy::genericE6SSM_soft_parameters , int);
double doCalcmtRSquaredLogSqCoeff(flexiblesusy::genericE6SSM_soft_parameters, int);

double doCalcmqL3SquaredLogCoeff(flexiblesusy::genericE6SSM_soft_parameters , int);
double doCalcmqL3SquaredLogSqCoeff(flexiblesusy::genericE6SSM_soft_parameters, int);



Eigen::Matrix<double,8,1> doCalcMh1SquaredDerivs(flexiblesusy::genericE6SSM_soft_parameters r, double ms, int gen,
				      bool & hasError);
Eigen::Matrix<double,8,1> doCalcMh2SquaredDerivs(flexiblesusy::genericE6SSM_soft_parameters r, double ms, int gen,
				      bool & hasError);
Eigen::Matrix<double,8,1> doCalcMsSquaredDerivs(flexiblesusy::genericE6SSM_soft_parameters r, double ms, int gen,
				      bool & hasError);
Eigen::Matrix<double,8,1> doCalcLambdaDerivs(flexiblesusy::genericE6SSM_soft_parameters r, double ms, int gen,
				bool & hasError);
Eigen::Matrix<double,8,1> doCalcGauginoDerivs(flexiblesusy::genericE6SSM_soft_parameters r, double ms, int whichGaugino,
				 bool & hasError);
Eigen::Matrix<double,8,1> doCalcSoftAuDerivs(flexiblesusy::genericE6SSM_soft_parameters r, double ms, int m, int n,
				bool & hasError);
Eigen::Matrix<double,8,1> doCalcSoftAlambdaDerivs(flexiblesusy::genericE6SSM_soft_parameters r, double ms, int gen,
				     bool & hasError);
Eigen::Matrix<double,8,1> doCalcMq2Derivs(flexiblesusy::genericE6SSM_soft_parameters r, double ms, int m, int n,
			     bool & hasError);
Eigen::Matrix<double,8,1> doCalcMu2Derivs(flexiblesusy::genericE6SSM_soft_parameters r, double ms, int m, int n,
			     bool & hasError);

Eigen::Matrix<double,8,1> doCalcRHSTuningVector_ESSM_Approx(flexiblesusy::genericE6SSM_soft_parameters modelAtMsusy, 
					       void (*ftBCatMX)(flexiblesusy::genericE6SSM_soft_parameters &, const Eigen::ArrayXd &), 
					       Eigen::ArrayXd const & vevs, tuning_parameters i, double mx, 
					       bool & hasError, Eigen::ArrayXd const & modelParsAtMx);

void pE6SSMftBCs(flexiblesusy::genericE6SSM_soft_parameters & m, Eigen::ArrayXd & tuningPars);

// A function for getting the predicted running value of M_Z^2, including 1-loop top and
// stop tadpole contributions. Used for numerically calculating the fine tuning in the
// E6SSM. Passed to SOFTSUSY's calcDerivative routine as part of the fine tuning calculation.
 double predpE6SSMMzSqRun(double parVal);


 Eigen::VectorXd doCalcESSMTuningNumerically(flexiblesusy::genericE6SSM_soft_parameters r, double ms, double mx, 
					     Eigen::ArrayXd pars,
					     void (*BCatMX)(flexiblesusy::genericE6SSM_soft_parameters & , Eigen::ArrayXd &));
 

} // namespace essm_tuning_utils
#endif
