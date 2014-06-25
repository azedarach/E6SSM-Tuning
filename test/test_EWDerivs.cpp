/*
  Test the expressions for the derivatives
  of the EWSB conditions w.r.t the VEVs and the 
  low scale parameters.
 */

#include "src/essmtuningutils.h"

using namespace flexiblesusy;
using namespace essm_tuning_utils;

double calcEWSBConditionVal(double vev);
double doCalcdEWSBConditiondVev(flexiblesusy::genericE6SSM_soft_parameters soft_model, int ewsb, int vev, bool & hasProblem);
double calcEWSBConditionValForParam(double param);
double doCalcdEWSBConditiondLambda(flexiblesusy::genericE6SSM_soft_parameters soft_model, int ewsb, bool & hasProblem);
double doCalcdEWSBConditiondAlambda(flexiblesusy::genericE6SSM_soft_parameters soft_model, int ewsb, bool & hasProblem);
double doCalcdEWSBConditiondMh1Squared(flexiblesusy::genericE6SSM_soft_parameters soft_model, int ewsb, bool & hasProblem);
double doCalcdEWSBConditiondMh2Squared(flexiblesusy::genericE6SSM_soft_parameters soft_model, int ewsb, bool & hasProblem);
double doCalcdEWSBConditiondMsSquared(flexiblesusy::genericE6SSM_soft_parameters soft_model, int ewsb, bool & hasProblem);
double doCalcdEWSBConditiondMqL3Squared(flexiblesusy::genericE6SSM_soft_parameters soft_model, int ewsb, bool & hasProblem);
double doCalcdEWSBConditiondMtRSquared(flexiblesusy::genericE6SSM_soft_parameters soft_model, int ewsb, bool & hasProblem);
double doCalcdEWSBConditiondAt(flexiblesusy::genericE6SSM_soft_parameters soft_model, int ewsb, bool & hasProblem);

/*
  --------------------------------------------------------------
  E6SSM U(1)_\psi and U(1)_\chi charges
  --------------------------------------------------------------
*/

const double QQPsi = 1.0/(2.0*flexiblesusy::Sqrt(6.0));
const double QuPsi = 1.0/(2.0*flexiblesusy::Sqrt(6.0));
const double QdPsi = 1.0/(2.0*flexiblesusy::Sqrt(6.0));
const double QLPsi = 1.0/(2.0*flexiblesusy::Sqrt(6.0));
const double QePsi = 1.0/(2.0*flexiblesusy::Sqrt(6.0));
const double QNPsi = 1.0/(2.0*flexiblesusy::Sqrt(6.0));
const double QSPsi = 4.0/(2.0*flexiblesusy::Sqrt(6.0));
const double QH1Psi = -1.0/flexiblesusy::Sqrt(6.0);
const double QH2Psi = -1.0/flexiblesusy::Sqrt(6.0);
const double QXPsi = -1.0/flexiblesusy::Sqrt(6.0);
const double QXbarPsi = -1.0/flexiblesusy::Sqrt(6.0);
const double QHPrPsi = 1.0/(2.0*flexiblesusy::Sqrt(6.0));
const double QHbarPrPsi = -1.0/(2.0*flexiblesusy::Sqrt(6.0));

const double QQChi = -1.0/(2.0*flexiblesusy::Sqrt(10.0));
const double QuChi = -1.0/(2.0*flexiblesusy::Sqrt(10.0));
const double QdChi = 3.0/(2.0*flexiblesusy::Sqrt(10.0));
const double QLChi = 3.0/(2.0*flexiblesusy::Sqrt(10.0));
const double QeChi = -1.0/(2.0*flexiblesusy::Sqrt(10.0));
const double QNChi = -5.0/(2.0*flexiblesusy::Sqrt(10.0));
const double QSChi = 0.0;
const double QH1Chi = -1.0/flexiblesusy::Sqrt(10.0);
const double QH2Chi = 1.0/flexiblesusy::Sqrt(10.0);
const double QXChi = 1.0/flexiblesusy::Sqrt(10.0);
const double QXbarChi = -1.0/flexiblesusy::Sqrt(10.0);
const double QHPrChi = 3.0/(2.0*flexiblesusy::Sqrt(10.0));
const double QHbarPrChi = -3.0/(2.0*flexiblesusy::Sqrt(10.0));

int main(int argc, char* argv[])
{

  /*
    Define maximum allowable deviation in %
   */
  double tol = 0.05;

  /*
    Pass or fail?
  */
  bool hasPassed = false;

  /*
    Set up model
   */

  // Values for model parameters
  Eigen::Matrix<double,3,3> yu, yd, ye, tyu, tyd, tye, kappa, tkappa;

  // Third generation couplings:
  yu(2,2) = 0.8;
  yd(2,2) = 0.1;
  ye(2,2) = 0.05;

  tyu(2,2) = yu(2,2) * 5000.0; // GeV
  tyd(2,2) = yd(2,2) * 5000.0; // GeV
  tye(2,2) = ye(2,2) * 5000.0; // GeV

  // We assume first and second generation Yukawas and A terms vanish,
  // so tests should be consistent with that assumption
  yu(0,0) = 0.0; yu(1,1) = 0.0;
  yd(0,0) = 0.0; yd(1,1) = 0.0;
  ye(0,0) = 0.0; ye(1,1) = 0.0;

  // Off-diagonal elements should vanish too
  yu(0,1) = 0.0; yu(0,2) = 0.0;
  yu(1,0) = 0.0; yu(1,2) = 0.0;
  yu(2,0) = 0.0; yu(2,1) = 0.0;

  yd(0,1) = 0.0; yd(0,2) = 0.0;
  yd(1,0) = 0.0; yd(1,2) = 0.0;
  yd(2,0) = 0.0; yd(2,1) = 0.0;

  ye(0,1) = 0.0; ye(0,2) = 0.0;
  ye(1,0) = 0.0; ye(1,2) = 0.0;
  ye(2,0) = 0.0; ye(2,1) = 0.0;

  tyu(0,1) = 0.0; tyu(0,2) = 0.0; // GeV
  tyu(1,0) = 0.0; tyu(1,2) = 0.0; // GeV
  tyu(2,0) = 0.0; tyu(2,1) = 0.0; // GeV

  tyd(0,1) = 0.0; tyd(0,2) = 0.0; // GeV
  tyd(1,0) = 0.0; tyd(1,2) = 0.0; // GeV
  tyd(2,0) = 0.0; tyd(2,1) = 0.0; // GeV

  tye(0,1) = 0.0; tye(0,2) = 0.0; // GeV
  tye(1,0) = 0.0; tye(1,2) = 0.0; // GeV
  tye(2,0) = 0.0; tye(2,1) = 0.0; // GeV

  // Kappa and T_kappa are assumed diagonal
  kappa(0,0) = 0.6; kappa(1,1) = 0.6; kappa(2,2) = 0.6;
  tkappa(0,0) = kappa(0,0) * 1000.0; tkappa(1,1) = kappa(1,1) * 1000.0; tkappa(2,2) = kappa(2,2) * 1000.0; // GeV

  kappa(0,1) = 0.0; kappa(0,2) = 0.0;
  kappa(1,0) = 0.0; kappa(1,2) = 0.0;
  kappa(2,0) = 0.0; kappa(2,1) = 0.0;

  tkappa(0,1) = 0.0; tkappa(0,2) = 0.0; // GeV
  tkappa(1,0) = 0.0; tkappa(1,2) = 0.0; // GeV
  tkappa(2,0) = 0.0; tkappa(2,1) = 0.0; // GeV

  Eigen::Matrix<double,2,2> lambda12, tlambda12, mH1ISq, mH2ISq, mSISq;

  // All soft parameters assumed diagonal
  lambda12(0,0) = 0.1; lambda12(0,1) = 0.0;
  lambda12(1,0) = 0.0; lambda12(1,1) = 0.2;

  tlambda12(0,0) = lambda12(0.0) * 500.0; tlambda12(0,1) = 0.0; // GeV
  tlambda12(1,0) = 0.0; tlambda12(1,1) = lambda12(1,1) * 500.0; // GeV

  mH1ISq(0,0) = Sqr(5000.0); mH1ISq(0,1) = 0.0; // GeV^2
  mH1ISq(1,0) = 0.0; mH1ISq(1,1) = Sqr(5000.0); // GeV^2

  mH2ISq(0,0) = Sqr(5000.0); mH2ISq(0,1) = 0.0; // GeV^2
  mH2ISq(1,0) = 0.0; mH2ISq(1,1) = Sqr(5000.0); // GeV^2

  mSISq(0,0) = Sqr(5000.0); mSISq(0,1) = 0.0; // GeV^2
  mSISq(1,0) = 0.0; mSISq(1,1) = Sqr(5000.0); // GeV^2

  double g1, g2, g3, gN, M1, M2, M3, M1p, mHpSq, mHpbarSq, lambda, tlambda, mHdSq, mHuSq, mSSq, mupr, Bmupr;

  g1 = 0.0;//0.45;
  gN = 0.0;//0.46;
  g2 = 0.0;//0.6;
  g3 = 0.0;//0.98;

  M1 = 300.0; // GeV
  M1p = 300.0; // GeV
  M2 = 500.0; // GeV
  M3 = 2000.0; // GeV

  mHpSq = Sqr(5000.0); // GeV^2
  mHpbarSq = Sqr(5000.0); // GeV^2

  lambda = 2.3;
  double Alambda = 5000.0;
  tlambda = lambda * Alambda; // GeV

  mHdSq = 1.0e8; // GeV^2
  mHuSq = 0.5e8; // GeV^2
  mSSq = -1.5e8; // GeV^2

  mupr = 5000.0; // GeV
  Bmupr = 5000.0; // GeV^2

  Eigen::Matrix<double,3,3> mQlSq, mLlSq, mUrSq, mDrSq, mErSq, mDxSq, mDxbarSq;

  mQlSq(2,2) = Sqr(2000.0); // GeV^2
  mUrSq(2,2) = Sqr(2000.0); // GeV^2

  // Soft masses are diagonal
  mQlSq(0,0) = Sqr(5000.0); mQlSq(1,1) = Sqr(5000.0); // GeV^2
  mUrSq(0,0) = Sqr(5000.0); mUrSq(1,1) = Sqr(5000.0); // GeV^2
  mDrSq(0,0) = Sqr(5000.0); mDrSq(1,1) = Sqr(5000.0); mDrSq(2,2) = Sqr(5000.0); // GeV^2
  mLlSq(0,0) = Sqr(5000.0); mLlSq(1,1) = Sqr(5000.0); mLlSq(2,2) = Sqr(5000.0); // GeV^2
  mErSq(0,0) = Sqr(5000.0); mErSq(1,1) = Sqr(5000.0); mErSq(2,2) = Sqr(5000.0); // GeV^2
  mDxSq(0,0) = Sqr(5000.0); mDxSq(1,1) = Sqr(5000.0); mDxSq(2,2) = Sqr(5000.0); // GeV^2
  mDxbarSq(0,0) = Sqr(5000.0); mDxbarSq(1,1) = Sqr(5000.0); mDxbarSq(2,2) = Sqr(5000.0); // GeV^2


  mQlSq(0,1) = 0.0; mQlSq(0,2) = 0.0; // GeV^2
  mQlSq(1,0) = 0.0; mQlSq(1,2) = 0.0; // GeV^2
  mQlSq(2,0) = 0.0; mQlSq(2,1) = 0.0; // GeV^2

  mUrSq(0,1) = 0.0; mUrSq(0,2) = 0.0; // GeV^2
  mUrSq(1,0) = 0.0; mUrSq(1,2) = 0.0; // GeV^2
  mUrSq(2,0) = 0.0; mUrSq(2,1) = 0.0; // GeV^2

  mDrSq(0,1) = 0.0; mDrSq(0,2) = 0.0; // GeV^2
  mDrSq(1,0) = 0.0; mDrSq(1,2) = 0.0; // GeV^2
  mDrSq(2,0) = 0.0; mDrSq(2,1) = 0.0; // GeV^2

  mLlSq(0,1) = 0.0; mLlSq(0,2) = 0.0; // GeV^2
  mLlSq(1,0) = 0.0; mLlSq(1,2) = 0.0; // GeV^2
  mLlSq(2,0) = 0.0; mLlSq(2,1) = 0.0; // GeV^2

  mErSq(0,1) = 0.0; mErSq(0,2) = 0.0; // GeV^2
  mErSq(1,0) = 0.0; mErSq(1,2) = 0.0; // GeV^2
  mErSq(2,0) = 0.0; mErSq(2,1) = 0.0; // GeV^2

  mDxSq(0,1) = 0.0; mDxSq(0,2) = 0.0; // GeV^2
  mDxSq(1,0) = 0.0; mDxSq(1,2) = 0.0; // GeV^2
  mDxSq(2,0) = 0.0; mDxSq(2,1) = 0.0; // GeV^2

  mDxbarSq(0,1) = 0.0; mDxbarSq(0,2) = 0.0; // GeV^2
  mDxbarSq(1,0) = 0.0; mDxbarSq(1,2) = 0.0; // GeV^2
  mDxbarSq(2,0) = 0.0; mDxbarSq(2,1) = 0.0; // GeV^2

  double v, s, tb, v1, v2;

  tb = 10.0;
  s = 6700.0; // GeV
  v = 246.0; // GeV
  v1 = v/Sqrt(1.0+tb*tb); // GeV
  v2 = v1*tb; // GeV

  double MX = 2000.0; // GeV

  int LOOPS = 2;
  int THRESH = 0;

  const double thetaE6 = ArcTan(Sqrt(15.0));

  // U(1)' charges
  double cE6 = cos(thetaE6);
  double sE6 = sin(thetaE6);
  
  // E6 input parameters
  genericE6SSM_input_parameters input = genericE6SSM_input_parameters();

  input.QQp = QQChi*cE6+QQPsi*sE6;
  input.Qup = QuChi*cE6+QuPsi*sE6;
  input.Qdp = QdChi*cE6+QdPsi*sE6;
  input.QLp = QLChi*cE6+QLPsi*sE6;
  input.Qep = QeChi*cE6+QePsi*sE6;
  input.QSp = QSChi*cE6+QSPsi*sE6;
  input.QH1p = QH1Chi*cE6+QH1Psi*sE6;
  input.QH2p = QH2Chi*cE6+QH2Psi*sE6;
  input.QDxp = QXChi*cE6+QXPsi*sE6;
  input.QDxbarp = QXbarChi*cE6+QXbarPsi*sE6;
  input.QHpp = QHPrChi*cE6+QHPrPsi*sE6;
  input.QHpbarp = QHbarPrChi*cE6+QHbarPrPsi*sE6;

  genericE6SSM_susy_parameters susy_model = genericE6SSM_susy_parameters(MX, LOOPS, THRESH, 
									 input, yd, ye, kappa, lambda12, 
									 lambda, yu, mupr, g1, g2, g3, gN, 
									 v1, v2, s);

  genericE6SSM_soft_parameters soft_model = genericE6SSM_soft_parameters(susy_model, tyd, tye, tkappa, tlambda12, 
									 tlambda, tyu, Bmupr, mQlSq, mLlSq, mHdSq, mHuSq,
									 mDrSq, mUrSq, mErSq, mSSq, mH1ISq, mH2ISq, mSISq,
									 mDxSq, mDxbarSq, mHpSq, mHpbarSq, M1, M2, M3, M1p);

  INCLUDE1LPTADPOLES = true;
  /*
    Do test
   */
try
  {
    Eigen::Matrix<double,3,1> vevs;
    vevs(0) = v1;
    vevs(1) = v2;
    vevs(2) = s;
    
    // Calculate derivatives analytically at the current scale
    Eigen::Matrix<double,3,3> dfdv = doCalcLHSTuningMatrix(soft_model, vevs);
    
    Eigen::Matrix<double,3,8> EWderivs;    

    EWderivs(0,0) = lambda*(v2*v2+s*s)-Alambda*s*tb/Sqrt(2.0);
    EWderivs(1,0) = lambda*(v1*v1+s*s)-Alambda*s/(Sqrt(2.0)*tb);
    EWderivs(2,0) = lambda*v*v-Alambda*v1*v2/(Sqrt(2.0)*s);
    
    EWderivs(0,1) = -lambda*s*tb/Sqrt(2.0);
    EWderivs(1,1) = -lambda*s/(Sqrt(2.0)*tb);
    EWderivs(2,1) = -lambda*v1*v2/(Sqrt(2.0)*s);
    
    EWderivs(0,2) = 1.0;
    EWderivs(1,2) = 0.0;
    EWderivs(2,2) = 0.0;
    
    EWderivs(0,3) = 0.0;
    EWderivs(1,3) = 1.0;
    EWderivs(2,3) = 0.0;
    
    EWderivs(0,4) = 0.0;
    EWderivs(1,4) = 0.0;
    EWderivs(2,4) = 1.0;
    
    EWderivs(0,5) = 0.0;
    EWderivs(1,5) = 0.0;
    EWderivs(2,5) = 0.0;
    
    EWderivs(0,6) = 0.0;
    EWderivs(1,6) = 0.0;
    EWderivs(2,6) = 0.0;
    
    EWderivs(0,7) = 0.0;
    EWderivs(1,7) = 0.0;
    EWderivs(2,7) = 0.0;
    
    
    if (INCLUDE1LPTADPOLES)
      {
	EWderivs(0,0) += doCalcd2DeltaVdLambdadv1(soft_model, s, tb);
	EWderivs(1,0) += doCalcd2DeltaVdLambdadv2(soft_model, s, tb);
	EWderivs(2,0) += doCalcd2DeltaVdLambdadv3(soft_model, s, tb);
	
	EWderivs(0,5) += doCalcd2DeltaVdmQlsqdv1(soft_model, s, tb);
	EWderivs(1,5) += doCalcd2DeltaVdmQlsqdv2(soft_model, s, tb);
	EWderivs(2,5) += doCalcd2DeltaVdmQlsqdv3(soft_model, s, tb);
	
	EWderivs(0,6) += doCalcd2DeltaVdmUrsqdv1(soft_model, s, tb);
	EWderivs(1,6) += doCalcd2DeltaVdmUrsqdv2(soft_model, s, tb);
	EWderivs(2,6) += doCalcd2DeltaVdmUrsqdv3(soft_model, s, tb);
	
	EWderivs(0,7) += doCalcd2DeltaVdAtdv1(soft_model, s, tb);
	EWderivs(1,7) += doCalcd2DeltaVdAtdv2(soft_model, s, tb);
	EWderivs(2,7) += doCalcd2DeltaVdAtdv3(soft_model, s, tb);
      }    
    
    
    // Calculate derivatives numerically
    bool hasProblem = false;
    double df1dv1 = doCalcdEWSBConditiondVev(soft_model, 1, 1, hasProblem);
    double df1dv2 = doCalcdEWSBConditiondVev(soft_model, 1, 2, hasProblem);
    double df1ds = doCalcdEWSBConditiondVev(soft_model, 1, 3, hasProblem);
    double df2dv1 = doCalcdEWSBConditiondVev(soft_model, 2, 1, hasProblem);
    double df2dv2 = doCalcdEWSBConditiondVev(soft_model, 2, 2, hasProblem);
    double df2ds = doCalcdEWSBConditiondVev(soft_model, 2, 3, hasProblem);
    double df3dv1 = doCalcdEWSBConditiondVev(soft_model, 3, 1, hasProblem);
    double df3dv2 = doCalcdEWSBConditiondVev(soft_model, 3, 2, hasProblem);
    double df3ds = doCalcdEWSBConditiondVev(soft_model, 3, 3, hasProblem);
    
    double df1dlambda = doCalcdEWSBConditiondLambda(soft_model, 1, hasProblem);
    double df2dlambda = doCalcdEWSBConditiondLambda(soft_model, 2, hasProblem);
    double df3dlambda = doCalcdEWSBConditiondLambda(soft_model, 3, hasProblem);
    
    double df1dAlambda = doCalcdEWSBConditiondAlambda(soft_model, 1, hasProblem);
    double df2dAlambda = doCalcdEWSBConditiondAlambda(soft_model, 2, hasProblem);
    double df3dAlambda = doCalcdEWSBConditiondAlambda(soft_model, 3, hasProblem);
    
    double df1dmHdSq = doCalcdEWSBConditiondMh1Squared(soft_model, 1, hasProblem);
    double df2dmHdSq = doCalcdEWSBConditiondMh1Squared(soft_model, 2, hasProblem);
    double df3dmHdSq = doCalcdEWSBConditiondMh1Squared(soft_model, 3, hasProblem);
    
    double df1dmHuSq = doCalcdEWSBConditiondMh2Squared(soft_model, 1, hasProblem);
    double df2dmHuSq = doCalcdEWSBConditiondMh2Squared(soft_model, 2, hasProblem);
    double df3dmHuSq = doCalcdEWSBConditiondMh2Squared(soft_model, 3, hasProblem);
    
    double df1dmSSq = doCalcdEWSBConditiondMsSquared(soft_model, 1, hasProblem);
    double df2dmSSq = doCalcdEWSBConditiondMsSquared(soft_model, 2, hasProblem);
    double df3dmSSq = doCalcdEWSBConditiondMsSquared(soft_model, 3, hasProblem);
    
    double df1dmqL3Sq = doCalcdEWSBConditiondMqL3Squared(soft_model, 1, hasProblem);
    double df2dmqL3Sq = doCalcdEWSBConditiondMqL3Squared(soft_model, 2, hasProblem);
    double df3dmqL3Sq = doCalcdEWSBConditiondMqL3Squared(soft_model, 3, hasProblem);
    
    double df1dmtRSq = doCalcdEWSBConditiondMtRSquared(soft_model, 1, hasProblem);
    double df2dmtRSq = doCalcdEWSBConditiondMtRSquared(soft_model, 2, hasProblem);
    double df3dmtRSq = doCalcdEWSBConditiondMtRSquared(soft_model, 3, hasProblem);
    
    double df1dAt = doCalcdEWSBConditiondAt(soft_model, 1, hasProblem);
    double df2dAt = doCalcdEWSBConditiondAt(soft_model, 2, hasProblem);
    double df3dAt = doCalcdEWSBConditiondAt(soft_model, 3, hasProblem);
    
    // Get relative errors
    double error_df1dv1, error_df1dv2, error_df1ds;
    double error_df2dv1, error_df2dv2, error_df2ds;
    double error_df3dv1, error_df3dv2, error_df3ds;

    error_df1dv1 = 100.0 * Abs((df1dv1-dfdv(0,0))/(0.5*(df1dv1+dfdv(0,0))));
    error_df1dv2 = 100.0 * Abs((df1dv2-dfdv(0,1))/(0.5*(df1dv2+dfdv(0,1))));
    error_df1ds = 100.0 * Abs((df1ds-dfdv(0,2))/(0.5*(df1ds+dfdv(0,2))));

    error_df2dv1 = 100.0 * Abs((df2dv1-dfdv(1,0))/(0.5*(df2dv1+dfdv(1,0))));
    error_df2dv2 = 100.0 * Abs((df2dv2-dfdv(1,1))/(0.5*(df2dv2+dfdv(1,1))));
    error_df2ds = 100.0 * Abs((df2ds-dfdv(1,2))/(0.5*(df2ds+dfdv(1,2))));

    error_df3dv1 = 100.0 * Abs((df3dv1-dfdv(2,0))/(0.5*(df3dv1+dfdv(2,0))));
    error_df3dv2 = 100.0 * Abs((df3dv2-dfdv(2,1))/(0.5*(df3dv2+dfdv(2,1))));
    error_df3ds = 100.0 * Abs((df3ds-dfdv(2,2))/(0.5*(df3ds+dfdv(2,2))));

    double error_df1dlambda, error_df2dlambda, error_df3dlambda;
    double error_df1dAlambda, error_df2dAlambda, error_df3dAlambda;
    double error_df1dmHdSq, error_df2dmHdSq, error_df3dmHdSq;
    double error_df1dmHuSq, error_df2dmHuSq, error_df3dmHuSq;
    double error_df1dmSSq, error_df2dmSSq, error_df3dmSSq;
    double error_df1dmqL3Sq, error_df2dmqL3Sq, error_df3dmqL3Sq;
    double error_df1dmtRSq, error_df2dmtRSq, error_df3dmtRSq;
    double error_df1dAt, error_df2dAt, error_df3dAt;

    error_df1dlambda = 100.0 * Abs((df1dlambda-EWderivs(0,0))/(0.5*(df1dlambda+EWderivs(0,0))));
    error_df2dlambda = 100.0 * Abs((df2dlambda-EWderivs(1,0))/(0.5*(df2dlambda+EWderivs(1,0))));
    error_df3dlambda = 100.0 * Abs((df3dlambda-EWderivs(2,0))/(0.5*(df3dlambda+EWderivs(2,0))));

    error_df1dAlambda = 100.0 * Abs((df1dAlambda-EWderivs(0,1))/(0.5*(df1dAlambda+EWderivs(0,1))));
    error_df2dAlambda = 100.0 * Abs((df2dAlambda-EWderivs(1,1))/(0.5*(df2dAlambda+EWderivs(1,1))));
    error_df3dAlambda = 100.0 * Abs((df3dAlambda-EWderivs(2,1))/(0.5*(df3dAlambda+EWderivs(2,1))));

    error_df1dmHdSq = 100.0 * Abs((df1dmHdSq-EWderivs(0,2))/(0.5*(df1dmHdSq+EWderivs(0,2))));
    error_df2dmHdSq = 100.0 * Abs((df2dmHdSq-EWderivs(1,2))/(0.5*(df2dmHdSq+EWderivs(1,2))));
    error_df3dmHdSq = 100.0 * Abs((df3dmHdSq-EWderivs(2,2))/(0.5*(df3dmHdSq+EWderivs(2,2))));

    error_df1dmHuSq = 100.0 * Abs((df1dmHuSq-EWderivs(0,3))/(0.5*(df1dmHuSq+EWderivs(0,3))));
    error_df2dmHuSq = 100.0 * Abs((df2dmHuSq-EWderivs(1,3))/(0.5*(df2dmHuSq+EWderivs(1,3))));
    error_df3dmHuSq = 100.0 * Abs((df3dmHuSq-EWderivs(2,3))/(0.5*(df3dmHuSq+EWderivs(2,3))));

    error_df1dmSSq = 100.0 * Abs((df1dmSSq-EWderivs(0,4))/(0.5*(df1dmSSq+EWderivs(0,4))));
    error_df2dmSSq = 100.0 * Abs((df2dmSSq-EWderivs(1,4))/(0.5*(df2dmSSq+EWderivs(1,4))));
    error_df3dmSSq = 100.0 * Abs((df3dmSSq-EWderivs(2,4))/(0.5*(df3dmSSq+EWderivs(2,4))));

    error_df1dmqL3Sq = 100.0 * Abs((df1dmqL3Sq-EWderivs(0,5))/(0.5*(df1dmqL3Sq+EWderivs(0,5))));
    error_df2dmqL3Sq = 100.0 * Abs((df2dmqL3Sq-EWderivs(1,5))/(0.5*(df2dmqL3Sq+EWderivs(1,5))));
    error_df3dmqL3Sq = 100.0 * Abs((df3dmqL3Sq-EWderivs(2,5))/(0.5*(df3dmqL3Sq+EWderivs(2,5))));

    error_df1dmtRSq = 100.0 * Abs((df1dmtRSq-EWderivs(0,6))/(0.5*(df1dmtRSq+EWderivs(0,6))));
    error_df2dmtRSq = 100.0 * Abs((df2dmtRSq-EWderivs(1,6))/(0.5*(df2dmtRSq+EWderivs(1,6))));
    error_df3dmtRSq = 100.0 * Abs((df3dmtRSq-EWderivs(2,6))/(0.5*(df3dmtRSq+EWderivs(2,6))));

    error_df1dAt = 100.0 * Abs((df1dAt-EWderivs(0,7))/(0.5*(df1dAt+EWderivs(0,7))));
    error_df2dAt = 100.0 * Abs((df2dAt-EWderivs(1,7))/(0.5*(df2dAt+EWderivs(1,7))));
    error_df3dAt = 100.0 * Abs((df3dAt-EWderivs(2,7))/(0.5*(df3dAt+EWderivs(2,7))));


    outputCharacteristics(8);
    
    // Print results
    cout << "**************************************************" << endl;
    cout << "* TEST RESULTS: e6ssm_EWSB_condition_derivs" << endl;
    cout << "**************************************************" << endl;
    cout << "* NUMERICAL DERIVATIVES: " << endl;
    cout << "*    df_1/dv_1 = " << df1dv1 << endl;
    cout << "*    df_1/dv_2 = " << df1dv2 << endl;
    cout << "*    df_1/ds = " << df1ds << endl;
    cout << "*    df_2/dv_1 = " << df2dv1 << endl;
    cout << "*    df_2/dv_2 = " << df2dv2 << endl;
    cout << "*    df_2/ds = " << df2ds << endl;
    cout << "*    df_3/dv_1 = " << df3dv1 << endl;
    cout << "*    df_3/dv_2 = " << df3dv2 << endl;
    cout << "*    df_3/ds = " << df3ds << endl;
    cout << "*    df_1/dlambda = " << df1dlambda << endl;
    cout << "*    df_2/dlambda = " << df2dlambda << endl;
    cout << "*    df_3/dlambda = " << df3dlambda << endl;
    cout << "*    df_1/dA_lambda = " << df1dAlambda << endl;
    cout << "*    df_2/dA_lambda = " << df2dAlambda << endl;
    cout << "*    df_3/dA_lambda = " << df3dAlambda << endl;
    cout << "*    df_1/dm_Hd^2 = " << df1dmHdSq << endl;
    cout << "*    df_2/dm_Hd^2 = " << df2dmHdSq << endl;
    cout << "*    df_3/dm_Hd^2 = " << df3dmHdSq << endl;
    cout << "*    df_1/dm_Hu^2 = " << df1dmHuSq << endl;
    cout << "*    df_2/dm_Hu^2 = " << df2dmHuSq << endl;
    cout << "*    df_3/dm_Hu^2 = " << df3dmHuSq << endl;
    cout << "*    df_1/dm_s^2 = " << df1dmSSq << endl;
    cout << "*    df_2/dm_s^2 = " << df2dmSSq << endl;
    cout << "*    df_3/dm_s^2 = " << df3dmSSq << endl;
    cout << "*    df_1/dm_q3L^2 = " << df1dmqL3Sq << endl;
    cout << "*    df_2/dm_q3L^2 = " << df2dmqL3Sq << endl;
    cout << "*    df_3/dm_q3L^2 = " << df3dmqL3Sq << endl;
    cout << "*    df_1/dm_tR^2 = " << df1dmtRSq << endl;
    cout << "*    df_2/dm_tR^2 = " << df2dmtRSq << endl;
    cout << "*    df_3/dm_tR^2 = " << df3dmtRSq << endl;
    cout << "*    df_1/dA_t = " << df1dAt << endl;
    cout << "*    df_2/dA_t = " << df2dAt << endl;
    cout << "*    df_3/dA_t = " << df3dAt << endl;
    cout << "**************************************************" << endl;
    cout << "* ANALYTIC DERIVATIVES: " << endl;
    cout << "*    df_1/dv_1 = " << dfdv(0,0) << endl;
    cout << "*    df_1/dv_2 = " << dfdv(0,1) << endl;
    cout << "*    df_1/ds = " << dfdv(0,2) << endl;
    cout << "*    df_2/dv_1 = " << dfdv(1,0) << endl;
    cout << "*    df_2/dv_2 = " << dfdv(1,1) << endl;
    cout << "*    df_2/ds = " << dfdv(1,2) << endl;
    cout << "*    df_3/dv_1 = " << dfdv(2,0) << endl;
    cout << "*    df_3/dv_2 = " << dfdv(2,1) << endl;
    cout << "*    df_3/ds = " << dfdv(2,2) << endl;
    cout << "*    df_1/dlambda = " << EWderivs(0,0) << endl;
    cout << "*    df_2/dlambda = " << EWderivs(1,0) << endl;
    cout << "*    df_3/dlambda = " << EWderivs(2,0) << endl;
    cout << "*    df_1/dA_lambda = " << EWderivs(0,1) << endl;
    cout << "*    df_2/dA_lambda = " << EWderivs(1,1) << endl;
    cout << "*    df_3/dA_lambda = " << EWderivs(2,1) << endl;
    cout << "*    df_1/dm_Hd^2 = " << EWderivs(0,2) << endl;
    cout << "*    df_2/dm_Hd^2 = " << EWderivs(1,2) << endl;
    cout << "*    df_3/dm_Hd^2 = " << EWderivs(2,2) << endl;
    cout << "*    df_1/dm_Hu^2 = " << EWderivs(0,3) << endl;
    cout << "*    df_2/dm_Hu^2 = " << EWderivs(1,3) << endl;
    cout << "*    df_3/dm_Hu^2 = " << EWderivs(2,3) << endl;
    cout << "*    df_1/dm_s^2 = " << EWderivs(0,4) << endl;
    cout << "*    df_2/dm_s^2 = " << EWderivs(1,4) << endl;
    cout << "*    df_3/dm_s^2 = " << EWderivs(2,4) << endl;
    cout << "*    df_1/dm_q3L^2 = " << EWderivs(0,5) << endl;
    cout << "*    df_2/dm_q3L^2 = " << EWderivs(1,5) << endl;
    cout << "*    df_3/dm_q3L^2 = " << EWderivs(2,5) << endl;
    cout << "*    df_1/dm_tR^2 = " << EWderivs(0,6) << endl;
    cout << "*    df_2/dm_tR^2 = " << EWderivs(1,6) << endl;
    cout << "*    df_3/dm_tR^2 = " << EWderivs(2,6) << endl;
    cout << "*    df_1/dA_t = " << EWderivs(0,7) << endl;
    cout << "*    df_2/dA_t = " << EWderivs(1,7) << endl;
    cout << "*    df_3/dA_t = " << EWderivs(2,7) << endl;
    cout << "**************************************************" << endl;
    cout << "* PERCENTAGE DIFFERENCES: " << endl;
    cout << "*    df_1/dv_1 = " << error_df1dv1 << "%" << endl;
    cout << "*    df_1/dv_2 = " << error_df1dv2 << "%" << endl;
    cout << "*    df_1/ds = " << error_df1ds << "%" << endl;
    cout << "*    df_2/dv_1 = " << error_df2dv1 << "%" << endl;
    cout << "*    df_2/dv_2 = " << error_df2dv2 << "%" << endl;
    cout << "*    df_2/ds = " << error_df2ds << "%" << endl;
    cout << "*    df_3/dv_1 = " << error_df3dv1 << "%" << endl;
    cout << "*    df_3/dv_2 = " << error_df3dv2 << "%" << endl;
    cout << "*    df_3/ds = " << error_df3ds << "%" << endl;
    cout << "*    df_1/dlambda = " << error_df1dlambda << "%" << endl;
    cout << "*    df_2/dlambda = " << error_df2dlambda << "%" << endl;
    cout << "*    df_3/dlambda = " << error_df3dlambda << "%" << endl;
    cout << "*    df_1/dA_lambda = " << error_df1dAlambda << "%" << endl;
    cout << "*    df_2/dA_lambda = " << error_df2dAlambda << "%" << endl;
    cout << "*    df_3/dA_lambda = " << error_df3dAlambda << "%" << endl;
    cout << "*    df_1/dm_Hd^2 = " << error_df1dmHdSq << "%" << endl;
    cout << "*    df_2/dm_Hd^2 = " << error_df2dmHdSq << "%" << endl;
    cout << "*    df_3/dm_Hd^2 = " << error_df3dmHdSq << "%" << endl;
    cout << "*    df_1/dm_Hu^2 = " << error_df1dmHuSq << "%" << endl;
    cout << "*    df_2/dm_Hu^2 = " << error_df2dmHuSq << "%" << endl;
    cout << "*    df_3/dm_Hu^2 = " << error_df3dmHuSq << "%" << endl;
    cout << "*    df_1/dm_s^2 = " << error_df1dmSSq << "%" << endl;
    cout << "*    df_2/dm_s^2 = " << error_df2dmSSq << "%" << endl;
    cout << "*    df_3/dm_s^2 = " << error_df3dmSSq << "%" << endl;
    cout << "*    df_1/dm_q3L^2 = " << error_df1dmqL3Sq << "%" << endl;
    cout << "*    df_2/dm_q3L^2 = " << error_df2dmqL3Sq << "%" << endl;
    cout << "*    df_3/dm_q3L^2 = " << error_df3dmqL3Sq << "%" << endl;
    cout << "*    df_1/dm_tR^2 = " << error_df1dmtRSq << "%" << endl;
    cout << "*    df_2/dm_tR^2 = " << error_df2dmtRSq << "%" << endl;
    cout << "*    df_3/dm_tR^2 = " << error_df3dmtRSq << "%" << endl;
    cout << "*    df_1/dA_t = " << error_df1dAt << "%" << endl;
    cout << "*    df_2/dA_t = " << error_df2dAt << "%" << endl;
    cout << "*    df_3/dA_t = " << error_df3dAt << "%" << endl;
    cout << "**************************************************" << endl;
    cout << "* RESULT: ";
    if (error_df1dv1 > tol || error_df1dv2 > tol || error_df1ds > tol ||
	error_df2dv1 > tol || error_df2dv2 > tol || error_df2ds > tol ||
	error_df3dv1 > tol || error_df3dv2 > tol || error_df3ds > tol || hasProblem || 
	error_df3dlambda > tol || error_df3dlambda > tol || error_df3dlambda > tol ||
	error_df3dAlambda > tol || error_df3dAlambda > tol || error_df3dAlambda > tol ||
	error_df3dmHdSq > tol || error_df3dmHdSq > tol || error_df3dmHdSq > tol ||
	error_df3dmHuSq > tol || error_df3dmHuSq > tol || error_df3dmHuSq > tol ||
	error_df3dmSSq > tol || error_df3dmSSq > tol || error_df3dmSSq > tol ||
	error_df3dmqL3Sq > tol || error_df3dmqL3Sq > tol || error_df3dmqL3Sq > tol ||
	error_df3dmtRSq > tol || error_df3dmtRSq > tol || error_df3dmtRSq > tol ||
	error_df3dAt > tol || error_df3dAt > tol || error_df3dAt > tol )
      {
	cout << "FAIL" << endl;
	if (error_df1dv1 > tol) cout << "*    maximum error exceeded for df_1/dv_1 derivative" << endl;
	if (error_df1dv2 > tol) cout << "*    maximum error exceeded for df_1/dv_2 derivative" << endl;
	if (error_df1ds > tol) cout << "*    maximum error exceeded for df_1/ds derivative" << endl;

	if (error_df2dv1 > tol) cout << "*    maximum error exceeded for df_2/dv_1 derivative" << endl;
	if (error_df2dv2 > tol) cout << "*    maximum error exceeded for df_2/dv_2 derivative" << endl;
	if (error_df2ds > tol) cout << "*    maximum error exceeded for df_2/ds derivative" << endl;

	if (error_df3dv1 > tol) cout << "*    maximum error exceeded for df_3/dv_1 derivative" << endl;
	if (error_df3dv2 > tol) cout << "*    maximum error exceeded for df_3/dv_2 derivative" << endl;
	if (error_df3ds > tol) cout << "*    maximum error exceeded for df_3/ds derivative" << endl;

	if (error_df1dlambda > tol) cout << "*    maximum error exceeded for df_1/dlambda derivative" << endl;
	if (error_df2dlambda > tol) cout << "*    maximum error exceeded for df_2/dlambda derivative" << endl;
	if (error_df3dlambda > tol) cout << "*    maximum error exceeded for df_3/dlambda derivative" << endl;

	if (error_df1dAlambda > tol) cout << "*    maximum error exceeded for df_1/dA_lambda derivative" << endl;
	if (error_df2dAlambda > tol) cout << "*    maximum error exceeded for df_2/dA_lambda derivative" << endl;
	if (error_df3dAlambda > tol) cout << "*    maximum error exceeded for df_3/dA_lambda derivative" << endl;

	if (error_df1dmHdSq > tol) cout << "*    maximum error exceeded for df_1/dm_Hd^2 derivative" << endl;
	if (error_df2dmHdSq > tol) cout << "*    maximum error exceeded for df_2/dm_Hd^2 derivative" << endl;
	if (error_df3dmHdSq > tol) cout << "*    maximum error exceeded for df_3/dm_Hd^2 derivative" << endl;

	if (error_df1dmHuSq > tol) cout << "*    maximum error exceeded for df_1/dm_Hu^2 derivative" << endl;
	if (error_df2dmHuSq > tol) cout << "*    maximum error exceeded for df_2/dm_Hu^2 derivative" << endl;
	if (error_df3dmHuSq > tol) cout << "*    maximum error exceeded for df_3/dm_Hu^2 derivative" << endl;

	if (error_df1dmSSq > tol) cout << "*    maximum error exceeded for df_1/dm_s^2 derivative" << endl;
	if (error_df2dmSSq > tol) cout << "*    maximum error exceeded for df_2/dm_s^2 derivative" << endl;
	if (error_df3dmSSq > tol) cout << "*    maximum error exceeded for df_3/dm_s^2 derivative" << endl;

	if (error_df1dmqL3Sq > tol) cout << "*    maximum error exceeded for df_1/dm_q3L^2 derivative" << endl;
	if (error_df2dmqL3Sq > tol) cout << "*    maximum error exceeded for df_2/dm_q3L^2 derivative" << endl;
	if (error_df3dmqL3Sq > tol) cout << "*    maximum error exceeded for df_3/dm_q3L^2 derivative" << endl;

	if (error_df1dmtRSq > tol) cout << "*    maximum error exceeded for df_1/dm_tR^2 derivative" << endl;
	if (error_df2dmtRSq > tol) cout << "*    maximum error exceeded for df_2/dm_tR^2 derivative" << endl;
	if (error_df3dmtRSq > tol) cout << "*    maximum error exceeded for df_3/dm_tR^2 derivative" << endl;

	if (error_df1dAt > tol) cout << "*    maximum error exceeded for df_1/dA_t derivative" << endl;
	if (error_df2dAt > tol) cout << "*    maximum error exceeded for df_2/dA_t derivative" << endl;
	if (error_df3dAt > tol) cout << "*    maximum error exceeded for df_3/dA_t derivative" << endl;

	hasPassed = false;
      }
    else
      {
	cout << "PASS" << endl; 
	hasPassed = true;
      }
    cout << "**************************************************" << endl;
    cout << "* END OF TEST" << endl;
    cout << "**************************************************" << endl;
    
  }
 catch(const string & a) 
   { 
     cout << "**************************************************" << endl;
     cout << "* " << a;
     cout << "\n* serious error during test: test aborted" << endl; 
     cout << "**************************************************" << endl;
     hasPassed = false;
   }
 catch(const char * a) 
   { 
     cout << "**************************************************" << endl;
     cout << "* " << a;
     cout << "\n* serious error during test: test aborted" << endl;
     cout << "**************************************************" << endl;
     hasPassed = false;
   }
 catch(...) 
   { 
     cout << "**************************************************" << endl;
     cout << "* WARNING: unknown error encountered during test: test aborted" << endl;
     cout << "**************************************************" << endl;
     hasPassed = false;
   }
 
 if (hasPassed)
   {
     return 0;
   }
 return 1;
 
}

// Useful variables to get data into the functions
static flexiblesusy::genericE6SSM_soft_parameters* tempmodel;
static int vevchoice;
static int ewsbchoice;

double calcEWSBConditionVal(double vev)
{
  genericE6SSM_soft_parameters model(*tempmodel);

  switch(vevchoice)
    {
    case 1: model.set_vd(vev); break;
    case 2: model.set_vu(vev); break;
    case 3: model.set_vs(vev); break;
    default:
      {
	ostringstream ii;
	ii << "ERROR: unrecognised choice of VEV in calcEWSBConditionVal" << endl;
	throw ii.str();
      }
    }

  switch(ewsbchoice)
    {
    case 1: return ESSM_EWSBCondition1(model); break;
    case 2: return ESSM_EWSBCondition2(model); break;
    case 3: return ESSM_EWSBCondition3(model); break;
    default:
      {
	ostringstream ii;
	ii << "ERROR: unrecognised choice of EWSB condition in calcEWSBConditionVal" << endl;
	throw ii.str();
      }
    }
  return 1; //< other errors
}

double doCalcdEWSBConditiondVev(flexiblesusy::genericE6SSM_soft_parameters soft_model, int ewsb, int vev, bool & hasProblem)
{
  tempmodel = &soft_model;
  ewsbchoice = ewsb;
  vevchoice = vev;

  double par;
  switch (vev)
    {
    case 1: par = soft_model.get_vd(); break;
    case 2: par = soft_model.get_vu(); break;
    case 3: par = soft_model.get_vs(); break;
    default:
      {
	ostringstream ii;
	ii << "ERROR: unrecognised choice of VEV in doCalcdEWSBConditiondVev" << endl;
	throw ii.str();
      }
    }

  double deriv, h, temp, err;

  double epsilon = 1.0e-5;

  // Initial estimate for step size h.
  h = epsilon*Abs(par);
  if (h == 0.0) h = epsilon;
  temp = par;
  par = temp + h;
  h = par - temp;

  deriv = calcDerivative(calcEWSBConditionVal, temp, h, &err);

  if (par > TOLERANCE && fabs(err / deriv) > 1.0) 
    {
      deriv = -numberOfTheBeast; // derivative is inaccurate, so flag using a large negative number
      hasProblem = true;
    }

  return deriv;


}

static int parchoice;
double calcEWSBConditionValForParam(double param)
{
  genericE6SSM_soft_parameters model(*tempmodel);

  switch(parchoice)
    {
    case tuning_parameters::lam3:
      {
	// Note that varying lambda ALSO
	// varies Tlambda, so we need to maintain
	// Alambda as constant
	double Tlambda = model.get_TLambdax();
	double lambda = model.get_Lambdax();
	double Alambda;
	
	if (Abs(Tlambda) < EPSTOL) 
	  {
	    Alambda = 0.0;
	  }
	else if (Abs(lambda) < 1.0e-100)
	  {
	    ostringstream ii;
	    ii << "WARNING: trying to calculate A_lambda where lambda3 coupling is " <<
	      Abs(lambda) << endl;
	    throw ii.str();
	  }
	else
	  {
	    Alambda = Tlambda/lambda;
	  }

	model.set_TLambdax(param*Alambda);
	model.set_Lambdax(param);
	break;
      }
    case tuning_parameters::Alam3:
      {
	model.set_TLambdax(model.get_Lambdax()*param);
	break;
      }
    case tuning_parameters::mH13Sq:
      {
	model.set_mHd2(param);
	break;
      }
    case tuning_parameters::mH23Sq:
      {
	model.set_mHu2(param);
	break;
      }
    case tuning_parameters::mS3Sq:
      {
	model.set_ms2(param);
	break;
      }
    case tuning_parameters::mqL3Sq:
      {
	model.set_mq2(2,2,param);
	break;
      }
    case tuning_parameters::mtRSq:
      {
	model.set_mu2(2,2,param);
	break;
      }
    case tuning_parameters::Au3:
      {
	model.set_TYu(2,2,model.get_Yu(2,2)*param);
	break;
      }
    default:
      {
	ostringstream ii;
	ii << "ERROR: unrecognised parameter choice in calcEWSBConditionValForParam" << endl;
	throw ii.str();
      }
    }

  switch(ewsbchoice)
    {
    case 1: return ESSM_EWSBCondition1(model); break;
    case 2: return ESSM_EWSBCondition2(model); break;
    case 3: return ESSM_EWSBCondition3(model); break;
    default:
      {
	ostringstream ii;
	ii << "ERROR: unrecognised EWSB condition choice in calcEWSBConditionValForParam" << endl;
	throw ii.str();
      }
    }

  return 1; //< other errors

}

double doCalcdEWSBConditiondLambda(flexiblesusy::genericE6SSM_soft_parameters soft_model, int ewsb, bool & hasProblem)
{
  tempmodel = &soft_model;
  ewsbchoice = ewsb;
  parchoice = tuning_parameters::lam3;

  double par = soft_model.get_Lambdax();

  double deriv, h, temp, err;

  double epsilon = 1.0e-5;

  // Initial estimate for step size h.
  h = epsilon*Abs(par);
  if (h == 0.0) h = epsilon;
  temp = par;
  par = temp + h;
  h = par - temp;

  deriv = calcDerivative(calcEWSBConditionValForParam, temp, h, &err);

  if (par > TOLERANCE && fabs(err / deriv) > 1.0) 
    {
      deriv = -numberOfTheBeast; // derivative is inaccurate, so flag using a large negative number
      hasProblem = true;
    }

  return deriv;

}

double doCalcdEWSBConditiondAlambda(flexiblesusy::genericE6SSM_soft_parameters soft_model, int ewsb, bool & hasProblem)
{
  tempmodel = &soft_model;
  ewsbchoice = ewsb;
  parchoice = tuning_parameters::Alam3;

  double par;

  double Tlambda = soft_model.get_TLambdax();
  double lambda = soft_model.get_Lambdax();
  double Alambda;
  
  if (Abs(Tlambda) < EPSTOL) 
    {
      Alambda = 0.0;
    }
  else if (Abs(lambda) < 1.0e-100)
    {
      ostringstream ii;
      ii << "WARNING: trying to calculate A_lambda where lambda3 coupling is " <<
	Abs(lambda) << endl;
      throw ii.str();
    }
  else
    {
      Alambda = Tlambda/lambda;
    }

  par = Alambda;  

  double deriv, h, temp, err;

  double epsilon = 1.0e-5;

  // Initial estimate for step size h.
  h = epsilon*Abs(par);
  if (h == 0.0) h = epsilon;
  temp = par;
  par = temp + h;
  h = par - temp;

  deriv = calcDerivative(calcEWSBConditionValForParam, temp, h, &err);

  if (par > TOLERANCE && fabs(err / deriv) > 1.0) 
    {
      deriv = -numberOfTheBeast; // derivative is inaccurate, so flag using a large negative number
      hasProblem = true;
    }

  return deriv;

}

double doCalcdEWSBConditiondMh1Squared(flexiblesusy::genericE6SSM_soft_parameters soft_model, int ewsb, bool & hasProblem)
{
  tempmodel = &soft_model;
  ewsbchoice = ewsb;
  parchoice = tuning_parameters::mH13Sq;

  double par = soft_model.get_mHd2();

  double deriv, h, temp, err;

  double epsilon = 1.0e-5;

  // Initial estimate for step size h.
  h = epsilon*Abs(par);
  if (h == 0.0) h = epsilon;
  temp = par;
  par = temp + h;
  h = par - temp;

  deriv = calcDerivative(calcEWSBConditionValForParam, temp, h, &err);

  if (par > TOLERANCE && fabs(err / deriv) > 1.0) 
    {
      deriv = -numberOfTheBeast; // derivative is inaccurate, so flag using a large negative number
      hasProblem = true;
    }

  return deriv;

}


double doCalcdEWSBConditiondMh2Squared(flexiblesusy::genericE6SSM_soft_parameters soft_model, int ewsb, bool & hasProblem)
{
  tempmodel = &soft_model;
  ewsbchoice = ewsb;
  parchoice = tuning_parameters::mH23Sq;

  double par = soft_model.get_mHu2();

  double deriv, h, temp, err;

  double epsilon = 1.0e-5;

  // Initial estimate for step size h.
  h = epsilon*Abs(par);
  if (h == 0.0) h = epsilon;
  temp = par;
  par = temp + h;
  h = par - temp;

  deriv = calcDerivative(calcEWSBConditionValForParam, temp, h, &err);

  if (par > TOLERANCE && fabs(err / deriv) > 1.0) 
    {
      deriv = -numberOfTheBeast; // derivative is inaccurate, so flag using a large negative number
      hasProblem = true;
    }

  return deriv;

}

double doCalcdEWSBConditiondMsSquared(flexiblesusy::genericE6SSM_soft_parameters soft_model, int ewsb, bool & hasProblem)
{
  tempmodel = &soft_model;
  ewsbchoice = ewsb;
  parchoice = tuning_parameters::mS3Sq;

  double par = soft_model.get_ms2();

  double deriv, h, temp, err;

  double epsilon = 1.0e-5;

  // Initial estimate for step size h.
  h = epsilon*Abs(par);
  if (h == 0.0) h = epsilon;
  temp = par;
  par = temp + h;
  h = par - temp;

  deriv = calcDerivative(calcEWSBConditionValForParam, temp, h, &err);

  if (par > TOLERANCE && fabs(err / deriv) > 1.0) 
    {
      deriv = -numberOfTheBeast; // derivative is inaccurate, so flag using a large negative number
      hasProblem = true;
    }

  return deriv;

}

double doCalcdEWSBConditiondMqL3Squared(flexiblesusy::genericE6SSM_soft_parameters soft_model, int ewsb, bool & hasProblem)
{
  tempmodel = &soft_model;
  ewsbchoice = ewsb;
  parchoice = tuning_parameters::mqL3Sq;

  double par = soft_model.get_mq2(2,2);

  double deriv, h, temp, err;

  double epsilon = 1.0e-5;

  // Initial estimate for step size h.
  h = epsilon*Abs(par);
  if (h == 0.0) h = epsilon;
  temp = par;
  par = temp + h;
  h = par - temp;

  deriv = calcDerivative(calcEWSBConditionValForParam, temp, h, &err);

  if (par > TOLERANCE && fabs(err / deriv) > 1.0) 
    {
      deriv = -numberOfTheBeast; // derivative is inaccurate, so flag using a large negative number
      hasProblem = true;
    }

  return deriv;

}

double doCalcdEWSBConditiondMtRSquared(flexiblesusy::genericE6SSM_soft_parameters soft_model, int ewsb, bool & hasProblem)
{
  tempmodel = &soft_model;
  ewsbchoice = ewsb;
  parchoice = tuning_parameters::mtRSq;

  double par = soft_model.get_mu2(2,2);

  double deriv, h, temp, err;

  double epsilon = 1.0e-5;

  // Initial estimate for step size h.
  h = epsilon*Abs(par);
  if (h == 0.0) h = epsilon;
  temp = par;
  par = temp + h;
  h = par - temp;

  deriv = calcDerivative(calcEWSBConditionValForParam, temp, h, &err);

  if (par > TOLERANCE && fabs(err / deriv) > 1.0) 
    {
      deriv = -numberOfTheBeast; // derivative is inaccurate, so flag using a large negative number
      hasProblem = true;
    }

  return deriv;

}

double doCalcdEWSBConditiondAt(flexiblesusy::genericE6SSM_soft_parameters soft_model, int ewsb, bool & hasProblem)
{
  tempmodel = &soft_model;
  ewsbchoice = ewsb;
  parchoice = tuning_parameters::Au3;

  double par;
  double yt = soft_model.get_Yu(2,2);
  double At;
  double TYt = soft_model.get_TYu(2,2);

  if (Abs(TYt) < EPSTOL)
    {
      At = 0.0;
    }
  else if (Abs(yt) < 1.0e-100)
    {
      ostringstream ii;
      ii << "WARNING: trying to calculate A_t where y_t coupling is " <<
	Abs(yt) << endl;
      throw ii.str();
    }
  else
    {
      At = TYt/yt;
    }

  par = At;

  double deriv, h, temp, err;

  double epsilon = 1.0e-5;

  // Initial estimate for step size h.
  h = epsilon*Abs(par);
  if (h == 0.0) h = epsilon;
  temp = par;
  par = temp + h;
  h = par - temp;

  deriv = calcDerivative(calcEWSBConditionValForParam, temp, h, &err);

  if (par > TOLERANCE && fabs(err / deriv) > 1.0) 
    {
      deriv = -numberOfTheBeast; // derivative is inaccurate, so flag using a large negative number
      hasProblem = true;
    }

  return deriv;

}
