/*
  mssmtuningutils contains all of the function needed for calculating the
  fine tuning in the (p)MSSM.

  Remaining checks as of 6/5/2014:
      - check all 2-loop approximate RGE solutions
      - double check tuning calculation (model-->EWSB-->tuning)
        is behaving as expected
      - settle error conditions (esp. small derivative, small error potential bug)
 */

#ifndef MSSMTUNINGUTILS_H
#define MSSMTUNINGUTILS_H

#include "softsusy_softsusy.h"
#include "flags.h"
#include "tuning_numerics.h"
#include "tuning_utils.h"

// A function to calculate the 2x2 matrix appearing on the LHS of the calculation
// of the derivatives of the VEVs wrt the parameters in the MSSM. Note that this
// assumes that the SoftParsMssm object is already set to have its renormalisation
// scale set at M_SUSY.
DoubleMatrix doCalcLHSTuningMatrix(SoftParsMssm mssmSusy, DoubleVector const & vevs);

double rgeDerivCalc(double x);

// A function to calculate the vectors appearing on the RHS of the calculation of the
// derivatives d v_i/ dp_j. Takes as arguments a SoftParsMssm object, assumed to be evaluated
// at the scale M_SUSY, vectors of the parameters to calculate the derivatives wrt. and the vevs,
// an integer labelling which particular parameter to calculate the derivative for, and the scale
// at which that parameter is defined.
DoubleVector doCalcRHSTuningVector(SoftParsMssm r, void (*ftBCatMX)(SoftParsMssm &, DoubleVector &), 
				   DoubleVector pars, DoubleVector const & vevs, 
				   int i, double mx, bool & hasError);

// A function to calculate d log(M_Z^2)/d log(p) using the value of the parameter p
// and the already calculated derivatives d v_i/ d p_j. 
double doCalcdLogMzSqdLogParam(SoftParsMssm r, double p, DoubleVector const & vevs, DoubleVector const & dVevsdp);

// MSSM_EWSBConditioni, i = 1, 2, calculates the value of the 
// EWSB condition for m_H_i^2.
// Inputs:
//     SoftParsMssm const & r = pMSSM model object
//     DoubleVector const & vevs = the values of v1 and v2
double MSSM_EWSBCondition1(SoftParsMssm const & r);
double MSSM_EWSBCondition2(SoftParsMssm const & r);

// Implement the EWSB conditions by varying the values of m_Hu^2 and m_Hd^2 at MX.
bool MSSM_ImplementEWSBConstraints_SoftMasses(SoftParsMssm r, double mx, double ms, bool useMxEqualsMs,
					      DoubleVector & updatedSoln, double tol = TOLERANCE);
bool MSSM_ImplementEWSBConstraints_SUGRA(SoftParsMssm r, double mx, double ms,
					      DoubleVector & updatedSoln, double tol = TOLERANCE);


bool MSSM_EWSB_NewtonShooter(SoftParsMssm const & r, DoubleVector & estimate, double tol = TOLERANCE);
bool SUGRA_EWSB_NewtonShooter(SoftParsMssm const & r, DoubleVector & estimate, double tol = TOLERANCE);

double MSSM_Msusy_Cond(SoftParsMssm r, double ms);

void MSSM_EWSB_Shooter_Functions(DoubleVector const & parVals, DoubleVector & f);
void SUGRA_EWSB_Shooter_Functions(DoubleVector const & parVals, DoubleVector & f);

// The function needed in the Newton iteration used to implement the EWSB conditions.
// Things are made slightly awkward by the fact that it should only take as arguments
// the values of the parameters and a vector to store the function values.  To get around
// it I have overloaded the SOFTSUSY numerics.h to allow an extra vector of 
// auxiliary data to be passed to define the function (see new definitions in tuningnumerics.h).
// In this case, guess contains the values of m_Hd^2(MX) and m_Hu^2(MX), 
// auxData contains all of the other parameters needed to define a SUSY model object,
// and on output the length 2 vector fVals contains the function values.
void MSSM_EWSBCondition_NewtonFunc(DoubleVector const & guess, DoubleVector & fVals);

DoubleVector MSSM_EWSBNewtonSolver(SoftParsMssm, double, int, int &);

// Calculates the DRbar' masses of the stops in the MSSM.
// Inputs:
//     SoftParsMssm r = MSSM model
//     DoubleVector & mstop = masses of stops
//     DoubleVector & mstop sq = squared masses of stops
//     double tb = value of tan(beta)
// Checked against SOFTSUSY results for model 2403883, seem correct.
void physical_MSSM(SoftParsMssm r, DoubleVector & mstop, DoubleVector & mstopsq, double tb);

// Function needed to evaluate tadpoles.
// Inputs:
//     double mSq = squared mass m^2
//     double Q = renormalisation scale Q
double a0MSSM(double mSq, double Q);

// Calculates the tadpole corrections to the EWSB conditions. Note
// returns (1/v_i)\partial \Delta V/\partial v_i (i.e. opposite sign
// to equivalent E_6SSM methods).
// Inputs:
//     SoftParsMssm r = MSSM model
//     double tb = tan(beta)
double doCalcTadpoleMSSMH1(SoftParsMssm r, double tb);
double doCalcTadpoleMSSMH2(SoftParsMssm r, double tb);

// The following are helper functions useful for constructing the 
// tuning measures at one-loop order.

// MSSMdmstop1sqdvi calculates the derivative m_stop_1^2 wrt v_i.
// Note our convention is m_stop_1 is the lighter stop.
// Inputs:
//     SoftParsMssm r = MSSM model
//     double tb = tan(beta)
double doCalcMSSMdmstop1sqdv1(SoftParsMssm r, double tb);
double doCalcMSSMdmstop1sqdv2(SoftParsMssm r, double tb);

// As above but for m_stop_2^2.
double doCalcMSSMdmstop2sqdv1(SoftParsMssm r, double tb);
double doCalcMSSMdmstop2sqdv2(SoftParsMssm r, double tb);

// These functions calculate the second derivatives of the one-loop
// corrections to the effective potential. Returns
// \Delta_{ij}'=\partial^2\Delta V/\partial v_i\partial v_j.
// Inputs:
//     SoftParsMssm r = MSSM model
//     double tb = tan(beta)
double doCalcMSSMDeltaPrime11(SoftParsMssm r, double tb);
double doCalcMSSMDeltaPrime12(SoftParsMssm r, double tb);
double doCalcMSSMDeltaPrime22(SoftParsMssm r, double tb);

// These helper functions calculate the second derivatives of 
// the one-loop corrections wrt the input parameters. They 
// return (1/v_i)\partial^2\Delta V/\partial p_j\partial v_i.
// Inputs:
//     SoftParsMssm r = pMSSM model
//     double tb = tan(beta)
double doCalcMSSMd2DeltaVdMudv1(SoftParsMssm, double);
double doCalcMSSMd2DeltaVdMudv2(SoftParsMssm, double);
double doCalcMSSMd2DeltaVdAtdv1(SoftParsMssm, double);
double doCalcMSSMd2DeltaVdAtdv2(SoftParsMssm, double);
double doCalcMSSMd2DeltaVdmQlSqdv1(SoftParsMssm, double);
double doCalcMSSMd2DeltaVdmQlSqdv2(SoftParsMssm, double);
double doCalcMSSMd2DeltaVdmUrSqdv1(SoftParsMssm, double);
double doCalcMSSMd2DeltaVdmUrSqdv2(SoftParsMssm, double);

double predpMSSMMzSqRun(double parVal);

// Because the EWSB conditions in the MSSM are so simple, it is straightforward
// to calculate the fine tuning numerically. This will be the check on our results.
// The object r is assumed to be provided at MX. The fine tuning is calculated
// as p / M_Z^2 \frac{\partial M_Z^2}{\partial p} where p is a set of parameters. M_Z^2
// is calculated using the tree level EWSB conditions + one-loop stop and top tadpole contributions
// at the scale M_SUSY. The function BCatMX is used to set the parameters at the high input scale MX.
 DoubleVector doCalcMSSMTuningNumerically(MssmSoftsusy r, double ms, double mx, 
 					 DoubleVector pars,
					  void (*BCatMX)(SoftParsMssm & , DoubleVector &));

// Functions for calculating the Taylor series approximate 
// solution to the RGEs for m_Hu^2 and m_Hd^2. These
// have been checked using the method doTestRGEApproxSolns
// and appear correct.
double doCalcMh2SquaredBeta1(SoftParsMssm);
double doCalcMh2SquaredBeta2(SoftParsMssm);
double doCalcMh2SquaredLogCoeff(SoftParsMssm, int);
double doCalcMh2SquaredLogSqCoeff(SoftParsMssm, int);

double doCalcMh1SquaredBeta1(SoftParsMssm);
double doCalcMh1SquaredBeta2(SoftParsMssm);
double doCalcMh1SquaredLogCoeff(SoftParsMssm , int);
double doCalcMh1SquaredLogSqCoeff(SoftParsMssm, int);

double doCalcMuBeta1(SoftParsMssm);
double doCalcMuBeta2(SoftParsMssm);
double doCalcMuLogCoeff(SoftParsMssm , int);
double doCalcMuLogSqCoeff(SoftParsMssm, int);

double doCalcBBeta1(SoftParsMssm);
double doCalcBBeta2(SoftParsMssm);
double doCalcBLogCoeff(SoftParsMssm , int);
double doCalcBLogSqCoeff(SoftParsMssm, int);

double doCalcAtBeta1(SoftParsMssm);
double doCalcAtBeta2(SoftParsMssm);
double doCalcAtLogCoeff(SoftParsMssm , int);
double doCalcAtLogSqCoeff(SoftParsMssm, int);

double doCalcmqL3SquaredBeta1(SoftParsMssm);
double doCalcmqL3SquaredBeta2(SoftParsMssm);
double doCalcmqL3SquaredLogCoeff(SoftParsMssm , int);
double doCalcmqL3SquaredLogSqCoeff(SoftParsMssm, int);

double doCalcmtRSquaredBeta1(SoftParsMssm);
double doCalcmtRSquaredBeta2(SoftParsMssm);
double doCalcmtRSquaredLogCoeff(SoftParsMssm , int);
double doCalcmtRSquaredLogSqCoeff(SoftParsMssm, int);

void doTestRGEApproxSolns(SoftParsMssm);
double getMuBeta1(double);
double getBBeta1(double);
double getmHdSqBeta1(double);
double getmHuSqBeta1(double);
double getmQlSqBeta1(double);
double getmUrSqBeta1(double);
double getAtBeta1(double);

// Functions for getting approximate derivatives of low scale parameters w.r.t high
// scale parameters. Returns the vector
// [ dmu/dp dB/dp dm_Hd^2/dp dm_Hu^2/dp dm_Ql^2/dp dm_uR^2/dp dA_t/dp]^T. The vector auxPars
// is assumed to contain the values of the gauge and Yukawa couplings at MX, in the order
// [ g1 g2 g3 yt yb ytau]^T. All of these still need to be checked for correctness.
DoubleVector doCalcMuDerivs(SoftParsMssm r, DoubleVector pars, double mx, bool & hasError,
			    DoubleVector const & auxPars);
DoubleVector doCalcBDerivs(SoftParsMssm r, DoubleVector pars, double mx, bool & hasError,
			    DoubleVector const & auxPars);
DoubleVector doCalcMh1SquaredDerivs(SoftParsMssm r, DoubleVector pars, double mx, 
				      bool & hasError, DoubleVector const & auxPars);
DoubleVector doCalcMh2SquaredDerivs(SoftParsMssm r, DoubleVector pars, double mx, 
				      bool & hasError, DoubleVector const & auxPars);
DoubleVector doCalcGauginoDerivs(SoftParsMssm r, DoubleVector pars, double mx, 
				 bool & hasError, DoubleVector const & auxPars, int whichGaugino);
DoubleVector doCalcSoftADerivs(SoftParsMssm r, DoubleVector pars, double mx, 
			       bool & hasError, DoubleVector const & auxPars, trilinears j, int m, int n);
DoubleVector doCalcSoftMassSquaredDerivs(SoftParsMssm r, DoubleVector pars, double mx, 
					 bool & hasError, DoubleVector const & auxPars, softMasses j, int m, int n);

DoubleVector doCalcRHSTuningVector_pMSSM_Approx(SoftParsMssm r, void (*ftBCatMX)(SoftParsMssm &, DoubleVector &), 
						DoubleVector pars, DoubleVector const & vevs, int i, double mx, 
						bool & hasError, DoubleVector const & auxPars);

void pMSSMftBCs(SoftParsMssm & m, DoubleVector & tuningPars);

void mSUGRAftBCs(SoftParsMssm & m, DoubleVector & tuningPars);

// This single function call will calculate the fine tuning in a pMSSM model. The given object
// is assumed to be provided at MX, and the parameters to calculate the fine tuning with
// are assumed to be already set in that object. The value ms is an initial guess for the scale
// M_{SUSY}, which will be recalculated when the EWSB conditions are imposed.
DoubleVector doCalcpMSSMFineTuning(SoftParsMssm r, double ms, bool & ewsbProblem, bool & hasError, 
				   bool useApproxSolns = false, double tol = TOLERANCE);

DoubleVector doCalcSUGRAFineTuning(SoftParsMssm r, double ms, bool & hasError, 
				   double tol = TOLERANCE);

#endif
