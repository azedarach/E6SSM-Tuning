/*
  essmtuningutils contains all of the functions needed
  for calculating the fine tuning in the (p)E6SSM.
 */

#ifndef ESSMTUNINGUTILS_H
#define ESSMTUNINGUTILS_H


#include "./E6SSM_Spectrum_Generators/models/genericE6SSM/genericE6SSM_two_scale_soft_parameters.hpp"
#include "./E6SSM_Spectrum_Generators/src/wrappers.hpp"
#include "flags.h"
#include "tuningnumerics.h"
#include "tuningutils.h"

// A function to calculate the 2x2 matrix appearing on the LHS of the calculation
// of the derivatives of the VEVs wrt the parameters in the MSSM. Note that this
// assumes that the SoftParsMssm object is already set to have its renormalisation
// scale set at M_SUSY.
DoubleMatrix doCalcLHSTuningMatrix(genericE6SSM_soft_parameters essmSusy, DoubleVector const & vevs);

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
double doCalcdLogMzSqdLogParam(genericE6SSM_soft_parameters r, double p, DoubleVector const & vevs, DoubleVector const & dVevsdp);

// ESSM_EWSBConditioni, i = 1, 2, calculates the value of the 
// EWSB condition for m_H_i^2 and m_S^2.
// Inputs:
//     SoftParsMssm const & r = pMSSM model object
//     DoubleVector const & vevs = the values of v1 and v2
double ESSM_EWSBCondition1(flexiblesusy::genericE6SSM_soft_parameters const & r);
double ESSM_EWSBCondition2(flexiblesusy::genericE6SSM_soft_parameters const & r);
double ESSM_EWSBCondition3(flexiblesusy::genericE6SSM_soft_parameters const & r);

// Implement the EWSB conditions by varying the values of m_Hu^2 and m_Hd^2 at MX.
/* bool ESSM_ImplementEWSBConstraints_SoftMasses(genericE6SSM_soft_parameters r, double mx, double ms, bool useMxEqualsMs, */
/* 					      DoubleVector & updatedSoln, double tol = TOLERANCE); */


/* bool ESSM_EWSB_NewtonShooter(SoftParsMssm const & r, DoubleVector & estimate, double tol = TOLERANCE); */

/* double ESSM_Msusy_Cond(SoftParsMssm r, double ms); */

/* void ESSM_EWSB_Shooter_Functions(DoubleVector const & parVals, DoubleVector & f); */

// The function needed in the Newton iteration used to implement the EWSB conditions.
// Things are made slightly awkward by the fact that it should only take as arguments
// the values of the parameters and a vector to store the function values.  To get around
// it I have overloaded the SOFTSUSY numerics.h to allow an extra vector of 
// auxiliary data to be passed to define the function (see new definitions in tuningnumerics.h).
// In this case, guess contains the values of m_Hd^2(MX) and m_Hu^2(MX), 
// auxData contains all of the other parameters needed to define a SUSY model object,
// and on output the length 2 vector fVals contains the function values.
/* void ESSM_EWSBCondition_NewtonFunc(DoubleVector const & guess, DoubleVector & fVals); */

/* DoubleVector ESSM_EWSBNewtonSolver(genericE6SSM_soft_parameters, double, int, int &); */

double ccbSqrt(double);

// Function a0 defined in paper. Implement appropriately.
double a0Peter(double, double);

// Calculates the masses of the stops and exotic quarks for tadpoles.  The former are always included 
// in our tadpole contributions, the latter could be, but provide small contributions in comparison to 
// errors due to threshold treatment.
// Inputs:
//    SoftParsEssm r = ESSM model
//    DoubleVector mstop = mass of stops
//    DoubleVector mstopsq = squared masses of stops
//    DoubleVector mD1sq = squared masses of exotic D quarks
//    DoubleVector mD2sq = squared masses of exotic D quarks
//    double s = singlet VEV
//    double tb = tan(beta)
void physical_ESSM(SoftParsEssm,DoubleVector &,DoubleVector &,DoubleVector &,DoubleVector &,double,double);

//This version neglects U(1) D-terms and was written for comparison with Roman's program
void physical_ESSM_Roman(SoftParsEssm,DoubleVector &,DoubleVector &,DoubleVector &,DoubleVector &,double,double);

// Determines the dominant H1 tadpoles. There are a number of different versions of this depending on choice of
// scale. The difference  should be higher order. Exotics can also be included in this version though by default
// they switched off at the moment.   This version  includes U(1) D-terms and exotics, and evaluates the 
// scale where RG evolution is halted.
// Inputs:
//    SoftParsEssm r = ESSM model
//    double s = singlet VEV
//    double tb = tan(beta)
double doCalcTadpoleESSMH1(SoftParsEssm,double,double);

// Calculates tadpoles using logarithms at m_t. The assumption is coefficients do not substantially change 
// between scale at which RG evolution is halted and m_t. 
double doCalcTadpoleESSMH1_atMt(SoftParsEssm,double,double);

// This version neglects U(1) D-terms and calculates at m_t. 
double doCalcTadpoleESSMH1_Roman(SoftParsEssm,double,double);

// Yet another option, performs calculation like Roman, but uses logarithms at scale at which RGE evolution is halted.
double doCalcTadpoleESSMH1_Roman_atQ(SoftParsEssm,double,double);

// Now tadpoles for H2. The different versions are labled as for H1 tadpoles and contributions are matched with those.
double doCalcTadpoleESSMH2(SoftParsEssm,double,double);

double doCalcTadpoleESSMH2_atMt(SoftParsEssm,double,double);

double doCalcTadpoleESSMH2_Roman(SoftParsEssm,double,double);

double doCalcTadpoleESSMH2_Roman_atQ(SoftParsEssm,double,double);

// Tadpoles for S. Same labeling as above.
double doCalcTadpolesESSMS(SoftParsEssm,double,double);

double doCalcTadpolesESSMS_atMt(SoftParsEssm,double,double);

double doCalcTadpolesESSMS_Roman(SoftParsEssm,double,double);

double doCalcTadpolesESSMS_Roman_atQ(SoftParsEssm,double,double);

// Calculate the Higgs masses in the model.
// Inputs:
//    SoftParsEssm & r = ESSM model
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
void HiggsMasses(SoftParsEssm &,double,double,DoubleVector &,DoubleVector &,int,bool,bool,DoubleVector &,int &,DoubleVector &, DoubleMatrix &, DoubleMatrix &, int &);

// To the above we will add methods for calculating m_A^2 at tree level and one-loop order, and
// methods for calculating chargino and neutralino masses (at tree level).

// Calculate m_A^2 at tree level.
// Inputs:
//    SoftParsEssm const & essmSusy = the ESSM model to calculate m_A^2 for
//    double s = the value of the singlet VEV to use
//    double tb = the value of tan(beta) to use
double mAsq_TreeLevel(SoftParsEssm const &, double, double);

// Calculate m_A^2 at one loop order.
// Inputs:
//    SoftParsEssm essmSusy = the ESSM model to calculate m_A^2 for
//    double s = the value of the singlet VEV to use
//    double tb = the value of tan(beta) to use
double mAsq_OneLoop(SoftParsEssm, double, double);

// A helper function used in calculating the mass m_A^2 at one loop order.
// Inputs:
//    double m1Sq = first squared mass m_1^2 to use
//    double m2Sq = second squared mass m_2^2 to use
//    double q = renormalisation scale
double f(double,double,double);


// Functions for getting approximate derivatives of low scale parameters w.r.t high
// scale parameters. Returns the vector
// [ dmu/dp dB/dp dm_Hd^2/dp dm_Hu^2/dp dm_Ql^2/dp dm_uR^2/dp dA_t/dp]^T. The vector auxPars
// is assumed to contain the values of the gauge and Yukawa couplings at MX, in the order
// [ g1 g2 g3 yt yb ytau]^T. All of these still need to be checked for correctness.
/* DoubleVector doCalcAlambdaDerivs(genericE6SSM_soft_parameters r, DoubleVector pars, double mx, int whichAlambda, */
/* 				 bool & hasError, DoubleVector const & auxPars); */
/* DoubleVector doCalcMh1SquaredDerivs(genericE6SSM_soft_parameters r, DoubleVector pars, double mx,  */
/* 				      bool & hasError, DoubleVector const & auxPars); */
/* DoubleVector doCalcMh2SquaredDerivs(genericE6SSM_soft_parameters r, DoubleVector pars, double mx,  */
/* 				      bool & hasError, DoubleVector const & auxPars); */
/* DoubleVector doCalcMsSquaredDerivs(genericE6SSM_soft_parameters r, DoubleVector pars, double mx,  */
/* 				      bool & hasError, DoubleVector const & auxPars); */
/* DoubleVector doCalcLambdaDerivs(genericE6SSM_soft_parameters r, DoubleVector pars, double mx, int whichLambda, */
/* 				bool & hasError, DoubleVector const & auxPars); */
/* DoubleVector doCalcKappaDerivs(genericE6SSM_soft_parameters r, DoubleVector pars, double mx, int whichKappa, */
/* 			       bool & hasError, DoubleVector const & auxPars); */
/* DoubleVector doCalcGauginoDerivs(genericE6SSM_soft_parameters r, DoubleVector pars, double mx,  */
/* 				 bool & hasError, DoubleVector const & auxPars, int whichGaugino); */
/* DoubleVector doCalcSoftADerivs(genericE6SSM_soft_parameters r, DoubleVector pars, double mx,  */
/* 			       bool & hasError, DoubleVector const & auxPars, trilinears j, int m, int n); */
/* DoubleVector doCalcSoftMassSquaredDerivs(genericE6SSM_soft_parameters r, DoubleVector pars, double mx,  */
/* 					 bool & hasError, DoubleVector const & auxPars, softMasses j, int m, int n); */

/* DoubleVector doCalcRHSTuningVector_pESSM_Approx(genericE6SSM_soft_parameters r, void (*ftBCatMX)(genericE6SSM_soft_parameters &, DoubleVector &),  */
/* 						DoubleVector pars, DoubleVector const & vevs, int i, double mx,  */
/* 						bool & hasError, DoubleVector const & auxPars); */

/* void pESSMftBCs(genericE6SSM_soft_parameters & m, DoubleVector & tuningPars); */



#endif
