/*
  Test the routines calculating the derivative w.r.t 
  t = log(Q/Q_0) of the 1-loop beta functions in the 
  E6SSM.
 */

#include "src/essmtuningutils.h"

using namespace flexiblesusy;
using namespace essm_tuning_utils;


double doCalcBeta1LoopMh1Squared(double);
double doCalcBeta1LoopMh2Squared(double);
double doCalcBeta1LoopMsSquared(double);
double doCalcBeta1LoopLambda(double);
double doCalcBeta1LoopAlambda(double);
double doCalcBeta1LoopAt(double);
double doCalcBeta1LoopMqL3Squared(double);
double doCalcBeta1LoopMtRSquared(double);
double calcBetaDeriv(genericE6SSM_soft_parameters, double, unsigned);


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

  g1 = 0.45;
  gN = 0.46;
  g2 = 0.6;
  g3 = 0.98;

  M1 = 300.0; // GeV
  M1p = 300.0; // GeV
  M2 = 500.0; // GeV
  M3 = 2000.0; // GeV

  mHpSq = Sqr(5000.0); // GeV^2
  mHpbarSq = Sqr(5000.0); // GeV^2

  lambda = -2.3;
  tlambda = lambda * (5000.0); // GeV

  mHdSq = 1.0e8; // GeV^2
  mHuSq = 05.e8; // GeV^2
  mSSq = -1.5e8; // GeV^2

  mupr = 5000.0; // GeV
  Bmupr = 5000.0; // GeV^2

  Eigen::Matrix<double,3,3> mQlSq, mLlSq, mUrSq, mDrSq, mErSq, mDxSq, mDxbarSq;

  mQlSq(2,2) = Sqr(1000.0); // GeV^2
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

  double MX = 20000.0; // GeV

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

  /*
    Do test
   */

  // Calculate derivatives of 1-loop beta functions numerically
  // and analytically and compare results.

  // Numerical calculation
  double numeric_dBeta1Mh1Sqdt, numeric_dBeta1Mh2Sqdt, numeric_dBeta1MsSqdt, numeric_dBeta1Atdt;
  double numeric_dBeta1Lambdadt, numeric_dBeta1MqSqdt, numeric_dBeta1MuSqdt, numeric_dBeta1Alambdadt;

  // Analytic calculation
  double analytic_dBeta1Mh1Sqdt, analytic_dBeta1Mh2Sqdt, analytic_dBeta1MsSqdt, analytic_dBeta1Atdt;
  double analytic_dBeta1Lambdadt, analytic_dBeta1MqSqdt, analytic_dBeta1MuSqdt, analytic_dBeta1Alambdadt;

  // % differences
  double error_dBeta1Mh1Sqdt, error_dBeta1Mh2Sqdt, error_dBeta1MsSqdt, error_dBeta1Atdt;
  double error_dBeta1Lambdadt, error_dBeta1MqSqdt, error_dBeta1MuSqdt, error_dBeta1Alambdadt;
  try
    {
      numeric_dBeta1Mh1Sqdt = calcBetaDeriv(soft_model, soft_model.get_scale(), tuning_parameters::mH13Sq);
      numeric_dBeta1Mh2Sqdt = calcBetaDeriv(soft_model, soft_model.get_scale(), tuning_parameters::mH23Sq);
      numeric_dBeta1MsSqdt = calcBetaDeriv(soft_model, soft_model.get_scale(), tuning_parameters::mS3Sq);
      numeric_dBeta1MqSqdt = calcBetaDeriv(soft_model, soft_model.get_scale(), tuning_parameters::mqL3Sq);
      numeric_dBeta1MuSqdt = calcBetaDeriv(soft_model, soft_model.get_scale(), tuning_parameters::mtRSq);
      numeric_dBeta1Atdt = calcBetaDeriv(soft_model, soft_model.get_scale(), tuning_parameters::Au3);
      numeric_dBeta1Lambdadt = calcBetaDeriv(soft_model, soft_model.get_scale(), tuning_parameters::lam3);
      numeric_dBeta1Alambdadt = calcBetaDeriv(soft_model, soft_model.get_scale(), tuning_parameters::Alam3);
      
      // Note that the analytic routines include an extra factor of 1/2 that needs to be accounted for
      analytic_dBeta1Mh1Sqdt = 2.0 * doCalcMh1SquaredLogSqCoeff(soft_model, 1);
      analytic_dBeta1Mh2Sqdt = 2.0 * doCalcMh2SquaredLogSqCoeff(soft_model, 1);
      analytic_dBeta1MsSqdt = 2.0 * doCalcMsSquaredLogSqCoeff(soft_model, 1);
      analytic_dBeta1MqSqdt = 2.0 * doCalcmqL3SquaredLogSqCoeff(soft_model, 1);
      analytic_dBeta1MuSqdt = 2.0 * doCalcmtRSquaredLogSqCoeff(soft_model, 1);
      analytic_dBeta1Atdt = 2.0 * doCalcAtLogSqCoeff(soft_model, 1);
      analytic_dBeta1Lambdadt = 2.0 * doCalcLambda3LogSqCoeff(soft_model, 1);
      analytic_dBeta1Alambdadt = 2.0 * doCalcAlambda3LogSqCoeff(soft_model,1);
      
      error_dBeta1Mh1Sqdt = 100.0*Abs((numeric_dBeta1Mh1Sqdt-analytic_dBeta1Mh1Sqdt)/(0.5*(numeric_dBeta1Mh1Sqdt+analytic_dBeta1Mh1Sqdt))); 
      error_dBeta1Mh2Sqdt = 100.0*Abs((numeric_dBeta1Mh2Sqdt-analytic_dBeta1Mh2Sqdt)/(0.5*(numeric_dBeta1Mh2Sqdt+analytic_dBeta1Mh2Sqdt))); 
      error_dBeta1MsSqdt = 100.0*Abs((numeric_dBeta1MsSqdt-analytic_dBeta1MsSqdt)/(0.5*(numeric_dBeta1MsSqdt+analytic_dBeta1MsSqdt))); 
      error_dBeta1MqSqdt = 100.0*Abs((numeric_dBeta1MqSqdt-analytic_dBeta1MqSqdt)/(0.5*(numeric_dBeta1MqSqdt+analytic_dBeta1MqSqdt))); 
      error_dBeta1MuSqdt = 100.0*Abs((numeric_dBeta1MuSqdt-analytic_dBeta1MuSqdt)/(0.5*(numeric_dBeta1MuSqdt+analytic_dBeta1MuSqdt))); 
      error_dBeta1Atdt = 100.0*Abs((numeric_dBeta1Atdt-analytic_dBeta1Atdt)/(0.5*(numeric_dBeta1Atdt+analytic_dBeta1Atdt))); 
      error_dBeta1Alambdadt = 100.0*Abs((numeric_dBeta1Alambdadt-analytic_dBeta1Alambdadt)/(0.5*(numeric_dBeta1Alambdadt+analytic_dBeta1Alambdadt))); 
      error_dBeta1Lambdadt = 100.0*Abs((numeric_dBeta1Lambdadt-analytic_dBeta1Lambdadt)/(0.5*(numeric_dBeta1Lambdadt+analytic_dBeta1Lambdadt))); 
      
      outputCharacteristics(8);
      
      // Print results
      cout << "**************************************************" << endl;
      cout << "* TEST RESULTS: e6ssm_1lp_beta_derivs" << endl;
      cout << "**************************************************" << endl;
      cout << "* NUMERICAL DERIVATIVES: " << endl;
      cout << "*    dbeta_{m_Hd^2}^{(1)}/dt = " << numeric_dBeta1Mh1Sqdt << endl;
      cout << "*    dbeta_{m_Hu^2}^{(1)}/dt = " << numeric_dBeta1Mh2Sqdt << endl;
      cout << "*    dbeta_{m_s^2}^{(1)}/dt = " << numeric_dBeta1MsSqdt << endl;
      cout << "*    dbeta_{m_q3L^2}^{(1)}/dt = " << numeric_dBeta1MqSqdt << endl;
      cout << "*    dbeta_{m_tR^2}^{(1)}/dt = " << numeric_dBeta1MuSqdt << endl;
      cout << "*    dbeta_{A_t}^{(1)}/dt = " << numeric_dBeta1Atdt << endl;
      cout << "*    dbeta_{A_lambda}^{(1)}/dt = " << numeric_dBeta1Alambdadt << endl;
      cout << "*    dbeta_{lambda}^{(1)}/dt = " << numeric_dBeta1Lambdadt << endl;
      cout << "**************************************************" << endl;
      cout << "* ANALYTIC DERIVATIVES: " << endl;
      cout << "*    dbeta_{m_Hd^2}^{(1)}/dt = " << analytic_dBeta1Mh1Sqdt << endl;
      cout << "*    dbeta_{m_Hu^2}^{(1)}/dt = " << analytic_dBeta1Mh2Sqdt << endl;
      cout << "*    dbeta_{m_s^2}^{(1)}/dt = " << analytic_dBeta1MsSqdt << endl;
      cout << "*    dbeta_{m_q3L^2}^{(1)}/dt = " << analytic_dBeta1MqSqdt << endl;
      cout << "*    dbeta_{m_tR^2}^{(1)}/dt = " << analytic_dBeta1MuSqdt << endl;
      cout << "*    dbeta_{A_t}^{(1)}/dt = " << analytic_dBeta1Atdt << endl;
      cout << "*    dbeta_{A_lambda}^{(1)}/dt = " << analytic_dBeta1Alambdadt << endl;
      cout << "*    dbeta_{lambda}^{(1)}/dt = " << analytic_dBeta1Lambdadt << endl;
      cout << "**************************************************" << endl;
      cout << "* PERCENTAGE DIFFERENCES: " << endl;
      cout << "*    dbeta_{m_Hd^2}^{(1)}/dt = " << error_dBeta1Mh1Sqdt << "%" << endl;
      cout << "*    dbeta_{m_Hu^2}^{(1)}/dt = " << error_dBeta1Mh2Sqdt << "%" << endl;
      cout << "*    dbeta_{m_s^2}^{(1)}/dt = " << error_dBeta1MsSqdt << "%" << endl;
      cout << "*    dbeta_{m_q3L^2}^{(1)}/dt = " << error_dBeta1MqSqdt << "%" << endl;
      cout << "*    dbeta_{m_tR^2}^{(1)}/dt = " << error_dBeta1MuSqdt << "%" << endl;
      cout << "*    dbeta_{A_t}^{(1)}/dt = " << error_dBeta1Atdt << "%" << endl;
      cout << "*    dbeta_{A_lambda}^{(1)}/dt = " << error_dBeta1Alambdadt << "%" << endl;
      cout << "*    dbeta_{lambda}^{(1)}/dt = " << error_dBeta1Lambdadt << "%" << endl;
      cout << "**************************************************" << endl;
      cout << "* RESULT: ";
      if (error_dBeta1Mh1Sqdt > tol || error_dBeta1Mh2Sqdt > tol || error_dBeta1MsSqdt > tol ||
	  error_dBeta1MqSqdt > tol || error_dBeta1MuSqdt > tol || error_dBeta1Atdt > tol ||
	  error_dBeta1Lambdadt > tol || error_dBeta1Alambdadt > tol)
	{
	  cout << "FAIL" << endl;
	  if (error_dBeta1Mh1Sqdt > tol) cout << "*    maximum error exceeded for beta m_Hd^2 derivative" << endl;
	  if (error_dBeta1Mh2Sqdt > tol) cout << "*    maximum error exceeded for beta m_Hu^2 derivative" << endl;
	  if (error_dBeta1MsSqdt > tol) cout << "*    maximum error exceeded for beta m_s^2 derivative" << endl;
	  if (error_dBeta1MqSqdt > tol) cout << "*    maximum error exceeded for beta m_qL3^2 derivative" << endl;
	  if (error_dBeta1MuSqdt > tol) cout << "*    maximum error exceeded for beta m_tR^2 derivative" << endl;
	  if (error_dBeta1Atdt > tol) cout << "*    maximum error exceeded for beta A_t derivative" << endl;
	  if (error_dBeta1Lambdadt > tol) cout << "*    maximum error exceeded for beta lambda derivative" << endl;
	  if (error_dBeta1Alambdadt > tol) cout << "*    maximum error exceeded for beta A_lambda derivative" << endl;
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
flexiblesusy::genericE6SSM_soft_parameters* tempmodel;

// Functions used for numerically calculating the derivative
// of the 1-loop beta functions. Each returns the value of
// the indicated beta function at the requested value of
// t = log(Q/Q_0).
double doCalcBeta1LoopMh1Squared(double t)
{
  
  // Get initial values
  genericE6SSM_soft_parameters model = *tempmodel; 

  model.set_loops(1); //< 1-loop betas only

  double q0 = model.get_scale();

  // Work out new scale to run to: t = log(q/q0)
  double q = q0 * exp(t);

  // Calculate beta functions at new scale
  model.run_to(q, PRECISION);

  genericE6SSM_soft_parameters beta = model.calc_beta();

  return beta.get_mHd2();
}

double doCalcBeta1LoopMh2Squared(double t)
{
  
  // Get initial values
  genericE6SSM_soft_parameters model = *tempmodel; 

  model.set_loops(1); //< 1-loop betas only

  double q0 = model.get_scale();

  // Work out new scale to run to: t = log(q/q0)
  double q = q0 * exp(t);

  // Calculate beta functions at new scale
  model.run_to(q, PRECISION);

  genericE6SSM_soft_parameters beta = model.calc_beta();

  return beta.get_mHu2();

}

double doCalcBeta1LoopMsSquared(double t)
{
  
  // Get initial values
  genericE6SSM_soft_parameters model = *tempmodel; 

  model.set_loops(1); //< 1-loop betas only

  double q0 = model.get_scale();

  // Work out new scale to run to: t = log(q/q0)
  double q = q0 * exp(t);

  // Calculate beta functions at new scale
  model.run_to(q, PRECISION);

  genericE6SSM_soft_parameters beta = model.calc_beta();

  return beta.get_ms2();

}

double doCalcBeta1LoopLambda(double t)
{
  
  // Get initial values
  genericE6SSM_soft_parameters model = *tempmodel; 

  model.set_loops(1); //< 1-loop betas only

  double q0 = model.get_scale();

  // Work out new scale to run to: t = log(q/q0)
  double q = q0 * exp(t);

  // Calculate beta functions at new scale
  model.run_to(q, PRECISION);

  genericE6SSM_soft_parameters beta = model.calc_beta();

  return beta.get_Lambdax();

}

double doCalcBeta1LoopAlambda(double t)
{
  
  // Get initial values
  genericE6SSM_soft_parameters model = *tempmodel; 

  model.set_loops(1); //< 1-loop betas only

  double q0 = model.get_scale();

  // Work out new scale to run to: t = log(q/q0)
  double q = q0 * exp(t);

  // Calculate beta functions at new scale
  model.run_to(q, PRECISION);

  genericE6SSM_soft_parameters beta = model.calc_beta();

  // Have to calculate beta_Alambda, since actual parameter
  // used in FlexibleSUSY is T_lambda = lambda * A_lambda

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
  
  double beta_Alambda = (beta.get_TLambdax() - Alambda * beta.get_Lambdax()) / lambda;

  return beta_Alambda;

}

double doCalcBeta1LoopMqL3Squared(double t)
{
  
  // Get initial values
  genericE6SSM_soft_parameters model = *tempmodel; 

  model.set_loops(1); //< 1-loop betas only

  double q0 = model.get_scale();

  // Work out new scale to run to: t = log(q/q0)
  double q = q0 * exp(t);

  // Calculate beta functions at new scale
  model.run_to(q, PRECISION);

  genericE6SSM_soft_parameters beta = model.calc_beta();

  return beta.get_mq2(2,2);
}

double doCalcBeta1LoopMtRSquared(double t)
{
  
  // Get initial values
  genericE6SSM_soft_parameters model = *tempmodel; 

  model.set_loops(1); //< 1-loop betas only

  double q0 = model.get_scale();

  // Work out new scale to run to: t = log(q/q0)
  double q = q0 * exp(t);

  // Calculate beta functions at new scale
  model.run_to(q, PRECISION);

  genericE6SSM_soft_parameters beta = model.calc_beta();

  return beta.get_mu2(2,2);
}

double doCalcBeta1LoopAt(double t)
{

  // Get initial values
  genericE6SSM_soft_parameters model = *tempmodel; 

  model.set_loops(1); //< 1-loop betas only

  double q0 = model.get_scale();

  // Work out new scale to run to: t = log(q/q0)
  double q = q0 * exp(t);

  // Calculate beta functions at new scale
  model.run_to(q, PRECISION);

  genericE6SSM_soft_parameters beta = model.calc_beta();

  double At;
  double TYt = model.get_TYu(2,2);
  double yt = model.get_Yu(2,2);

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

  double beta_At = (beta.get_TYu(2,2) - At * beta.get_Yu(2,2)) / yt;

  return beta_At;
}

// Calculates the derivative of the selected beta function at the 
// current scale of the provided object (i.e. the result of get_scale()
// is assumed to provide the value of Q_0 about which a Taylor
// series expansion would be constructed).
double calcBetaDeriv(genericE6SSM_soft_parameters r, double q, unsigned i)
{

  tempmodel = &r;

  double deriv, h, temp, err;

  double epsilon = 1.0e-5;

  double t0 = Log(q/r.get_scale()); //< estimate derivative at requested scale

  // Initial estimate for step size h.
  h = epsilon*Abs(t0);
  if (h == 0.0) h = epsilon;
  temp = t0;
  t0 = temp + h;
  h = t0 - temp;

  switch(i)
    {
    case tuning_parameters::lam3:
      {
	deriv = calcDerivative(doCalcBeta1LoopLambda, temp, h, &err);
	break;
      }
    case tuning_parameters::Alam3:
      {
	deriv = calcDerivative(doCalcBeta1LoopAlambda, temp, h, &err);
	break;
      }
    case tuning_parameters::Au3:
      {
	deriv = calcDerivative(doCalcBeta1LoopAt, temp, h, &err);
	break;
      }
    case tuning_parameters::mH13Sq:
      {
	deriv = calcDerivative(doCalcBeta1LoopMh1Squared, temp, h, &err);
	break;
      }
    case tuning_parameters::mH23Sq:
      {
	deriv = calcDerivative(doCalcBeta1LoopMh2Squared, temp, h, &err);
	break;
      }
    case tuning_parameters::mS3Sq:
      {
	deriv = calcDerivative(doCalcBeta1LoopMsSquared, temp, h, &err);
	break;
      }
    case tuning_parameters::mqL3Sq:
      {
	deriv = calcDerivative(doCalcBeta1LoopMqL3Squared, temp, h, &err);
	break;
      }
    case tuning_parameters::mtRSq:
      {
	deriv = calcDerivative(doCalcBeta1LoopMtRSquared, temp, h, &err);
	break;
      }
    default:
      {
	ostringstream ii;
	ii << "ERROR: unrecognised parameter requested: in calcBetaDeriv." << endl;
	throw ii.str();
      }
    }

  if (t0 > TOLERANCE && fabs(err / deriv) > 1.0) 
    {
      deriv = -numberOfTheBeast; // derivative is inaccurate, so flag using a large negative number
    }

  return deriv;

}
