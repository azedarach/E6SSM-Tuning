/*
  Test the routines that calculate the partial derivatives of
  the low scale parameters w.r.t the high scale parameters.
  The low scale parameters are calculated using a Taylor
  series approximation to the solutions of the E6SSM
  RGEs.
 */

#include "src/essmtuningutils.h"

using namespace flexiblesusy;
using namespace essm_tuning_utils;

Eigen::MatrixXd calcPercentageDifferences(const Eigen::MatrixXd & a, const Eigen::MatrixXd& b);

double getApproximateMh1Squared(double param);
double getApproximateMh2Squared(double param);
double getApproximateMsSquared(double param);
double getApproximateMqL3Squared(double param);
double getApproximateMtRSquared(double param);
double getApproximateLambda(double param);
double getApproximateAlambda(double param);
double getApproximateAt(double param);
Eigen::Matrix<double,8,1> doCalcNumericDerivs(genericE6SSM_soft_parameters r, double q, unsigned i, int lps, int logs);

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

  lambda = -1.9;
  tlambda = lambda * (5000.0); // GeV

  mHdSq = 1.0e8; // GeV^2
  mHuSq = -1.5e8; // GeV^2
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
  double q = 2300.0; // GeV

  int LOOPS = 2;
  int LEADINGLOGS = 1;
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

  // Calculate derivatives of Taylor series solutions
  // numerically and analytically

  int const NUMEWSBPARS = 8; //< number of parameters appearing in 1-loop EWSB equations

  // Numerical calculation
  Eigen::Matrix<double,NUMEWSBPARS,1> numeric_dParamsdMh2Sq, numeric_dParamsdMh1Sq, numeric_dParamsdMsSq;
  Eigen::Matrix<double, NUMEWSBPARS,1> numeric_dParamsdMqSq, numeric_dParamsdMuSq;
  Eigen::Matrix<double,NUMEWSBPARS,1> numeric_dParamsdLambda, numeric_dParamsdAlambda, numeric_dParamsdAt;
  Eigen::Matrix<double,NUMEWSBPARS,1> numeric_dParamsdM1, numeric_dParamsdM2, numeric_dParamsdM3, numeric_dParamsdM1p;

  // Analytic calculation
  Eigen::Matrix<double,NUMEWSBPARS,1> analytic_dParamsdMh2Sq, analytic_dParamsdMh1Sq, analytic_dParamsdMsSq;
  Eigen::Matrix<double, NUMEWSBPARS,1> analytic_dParamsdMqSq, analytic_dParamsdMuSq;
  Eigen::Matrix<double,NUMEWSBPARS,1> analytic_dParamsdLambda, analytic_dParamsdAlambda, analytic_dParamsdAt;
  Eigen::Matrix<double,NUMEWSBPARS,1> analytic_dParamsdM1, analytic_dParamsdM2, analytic_dParamsdM3, analytic_dParamsdM1p;

  // % differences
  Eigen::Matrix<double,NUMEWSBPARS,1> error_dParamsdMh2Sq, error_dParamsdMh1Sq, error_dParamsdMsSq;
  Eigen::Matrix<double, NUMEWSBPARS,1> error_dParamsdMqSq, error_dParamsdMuSq;
  Eigen::Matrix<double,NUMEWSBPARS,1> error_dParamsdLambda, error_dParamsdAlambda, error_dParamsdAt;
  Eigen::Matrix<double,NUMEWSBPARS,1> error_dParamsdM1, error_dParamsdM2, error_dParamsdM3, error_dParamsdM1p;

  bool hasError = false;

  bool hasPassedMh1Sq = false;
  bool hasPassedMh2Sq = false;
  bool hasPassedMsSq = false;
  bool hasPassedMqSq = false;
  bool hasPassedMuSq = false;
  bool hasPassedLambda = false;
  bool hasPassedAlambda = false;
  bool hasPassedAt = false;
  bool hasPassedM1 = false;
  bool hasPassedM2 = false;
  bool hasPassedM3 = false;
  bool hasPassedM1p = false;

  try
    {
      numeric_dParamsdMh2Sq = doCalcNumericDerivs(soft_model, q, tuning_parameters::mH23Sq, LOOPS, LEADINGLOGS);
      numeric_dParamsdMh1Sq = doCalcNumericDerivs(soft_model, q, tuning_parameters::mH13Sq, LOOPS, LEADINGLOGS);
      numeric_dParamsdMsSq = doCalcNumericDerivs(soft_model, q, tuning_parameters::mS3Sq, LOOPS, LEADINGLOGS);
      numeric_dParamsdMqSq = doCalcNumericDerivs(soft_model, q, tuning_parameters::mqL3Sq, LOOPS, LEADINGLOGS);
      numeric_dParamsdMuSq = doCalcNumericDerivs(soft_model, q, tuning_parameters::mtRSq, LOOPS, LEADINGLOGS);
      numeric_dParamsdLambda = doCalcNumericDerivs(soft_model, q, tuning_parameters::lam3, LOOPS, LEADINGLOGS);
      numeric_dParamsdAlambda = doCalcNumericDerivs(soft_model, q, tuning_parameters::Alam3, LOOPS, LEADINGLOGS);
      numeric_dParamsdAt = doCalcNumericDerivs(soft_model, q, tuning_parameters::Au3, LOOPS, LEADINGLOGS);
      numeric_dParamsdM1 = doCalcNumericDerivs(soft_model, q, tuning_parameters::M1, LOOPS, LEADINGLOGS);
      numeric_dParamsdM2 = doCalcNumericDerivs(soft_model, q, tuning_parameters::M2, LOOPS, LEADINGLOGS);
      numeric_dParamsdM3 = doCalcNumericDerivs(soft_model, q, tuning_parameters::M3, LOOPS, LEADINGLOGS);
      numeric_dParamsdM1p = doCalcNumericDerivs(soft_model, q, tuning_parameters::M1p, LOOPS, LEADINGLOGS);

      analytic_dParamsdMh2Sq = doCalcMh2SquaredDerivs(soft_model, q, 3, hasError);
      analytic_dParamsdMh1Sq = doCalcMh1SquaredDerivs(soft_model, q, 3, hasError);
      analytic_dParamsdMsSq = doCalcMsSquaredDerivs(soft_model, q, 3, hasError);
      analytic_dParamsdMqSq = doCalcMq2Derivs(soft_model, q, 3, 3, hasError);
      analytic_dParamsdMuSq = doCalcMu2Derivs(soft_model, q, 3, 3, hasError);
      analytic_dParamsdLambda = doCalcLambdaDerivs(soft_model, q, 3, hasError);
      analytic_dParamsdAlambda = doCalcSoftAlambdaDerivs(soft_model, q, 3, hasError);
      analytic_dParamsdAt = doCalcSoftAuDerivs(soft_model, q, 3, 3, hasError);
      analytic_dParamsdM1 = doCalcGauginoDerivs(soft_model, q, 1, hasError);
      analytic_dParamsdM2 = doCalcGauginoDerivs(soft_model, q, 2, hasError);
      analytic_dParamsdM3 = doCalcGauginoDerivs(soft_model, q, 3, hasError);
      analytic_dParamsdM1p = doCalcGauginoDerivs(soft_model, q, 4, hasError);

      error_dParamsdMh2Sq = calcPercentageDifferences(numeric_dParamsdMh2Sq, analytic_dParamsdMh2Sq);
      error_dParamsdMh1Sq = calcPercentageDifferences(numeric_dParamsdMh1Sq, analytic_dParamsdMh1Sq);
      error_dParamsdMsSq = calcPercentageDifferences(numeric_dParamsdMsSq, analytic_dParamsdMsSq);
      error_dParamsdMqSq = calcPercentageDifferences(numeric_dParamsdMqSq, analytic_dParamsdMqSq);
      error_dParamsdMuSq = calcPercentageDifferences(numeric_dParamsdMuSq, analytic_dParamsdMuSq);
      error_dParamsdLambda = calcPercentageDifferences(numeric_dParamsdLambda, analytic_dParamsdLambda);
      error_dParamsdAlambda = calcPercentageDifferences(numeric_dParamsdAlambda, analytic_dParamsdAlambda);
      error_dParamsdAt = calcPercentageDifferences(numeric_dParamsdAt, analytic_dParamsdAt);
      error_dParamsdM1 = calcPercentageDifferences(numeric_dParamsdM1, analytic_dParamsdM1);
      error_dParamsdM2 = calcPercentageDifferences(numeric_dParamsdM2, analytic_dParamsdM2);
      error_dParamsdM3 = calcPercentageDifferences(numeric_dParamsdM3, analytic_dParamsdM3);
      error_dParamsdM1p = calcPercentageDifferences(numeric_dParamsdM1p, analytic_dParamsdM1p);

      outputCharacteristics(8);
      
      // Print results
      cout << "**************************************************" << endl;
      cout << "* TEST RESULTS: e6ssm_low_scale_param_derivs" << endl;
      cout << "**************************************************" << endl;
      cout << "* SUBTEST 1: m_Hd^2 derivatives" << endl; 
      cout << "**************************************************" << endl;
      cout << "* NUMERICAL DERIVATIVES: " << endl;
      cout << "*    dlambda(M_SUSY)/dm_Hd^2(M_X) = " << numeric_dParamsdMh1Sq(0) << endl;
      cout << "*    dA_lambda(M_SUSY)/dm_Hd^2(M_X) = " << numeric_dParamsdMh1Sq(1) << endl;
      cout << "*    dm_Hd^2(M_SUSY)/dm_Hd^2(M_X) = " << numeric_dParamsdMh1Sq(2) << endl;
      cout << "*    dm_Hu^2(M_SUSY)/dm_Hd^2(M_X) = " << numeric_dParamsdMh1Sq(3) << endl;
      cout << "*    dm_s^2(M_SUSY)/dm_Hd^2(M_X) = " << numeric_dParamsdMh1Sq(4) << endl;
      cout << "*    dm_q3L^2(M_SUSY)/dm_Hd^2(M_X) = " << numeric_dParamsdMh1Sq(5) << endl;
      cout << "*    dm_tR^2(M_SUSY)/dm_Hd^2(M_X) = " << numeric_dParamsdMh1Sq(6) << endl;
      cout << "*    dA_t(M_SUSY)/dm_Hd^2(M_X) = " << numeric_dParamsdMh1Sq(7) << endl;
      cout << "**************************************************" << endl;
      cout << "* ANALYTIC DERIVATIVES: " << endl;
      cout << "*    dlambda(M_SUSY)/dm_Hd^2(M_X) = " << analytic_dParamsdMh1Sq(0) << endl;
      cout << "*    dA_lambda(M_SUSY)/dm_Hd^2(M_X) = " << analytic_dParamsdMh1Sq(1) << endl;
      cout << "*    dm_Hd^2(M_SUSY)/dm_Hd^2(M_X) = " << analytic_dParamsdMh1Sq(2) << endl;
      cout << "*    dm_Hu^2(M_SUSY)/dm_Hd^2(M_X) = " << analytic_dParamsdMh1Sq(3) << endl;
      cout << "*    dm_s^2(M_SUSY)/dm_Hd^2(M_X) = " << analytic_dParamsdMh1Sq(4) << endl;
      cout << "*    dm_q3L^2(M_SUSY)/dm_Hd^2(M_X) = " << analytic_dParamsdMh1Sq(5) << endl;
      cout << "*    dm_tR^2(M_SUSY)/dm_Hd^2(M_X) = " << analytic_dParamsdMh1Sq(6) << endl;
      cout << "*    dA_t(M_SUSY)/dm_Hd^2(M_X) = " << analytic_dParamsdMh1Sq(7) << endl;
      cout << "**************************************************" << endl;
      cout << "* PERCENTAGE DIFFERENCES: " << endl;
      cout << "*    dlambda(M_SUSY)/dm_Hd^2(M_X) = " << error_dParamsdMh1Sq(0) << endl;
      cout << "*    dA_lambda(M_SUSY)/dm_Hd^2(M_X) = " << error_dParamsdMh1Sq(1) << endl;
      cout << "*    dm_Hd^2(M_SUSY)/dm_Hd^2(M_X) = " << error_dParamsdMh1Sq(2) << endl;
      cout << "*    dm_Hu^2(M_SUSY)/dm_Hd^2(M_X) = " << error_dParamsdMh1Sq(3) << endl;
      cout << "*    dm_s^2(M_SUSY)/dm_Hd^2(M_X) = " << error_dParamsdMh1Sq(4) << endl;
      cout << "*    dm_q3L^2(M_SUSY)/dm_Hd^2(M_X) = " << error_dParamsdMh1Sq(5) << endl;
      cout << "*    dm_tR^2(M_SUSY)/dm_Hd^2(M_X) = " << error_dParamsdMh1Sq(6) << endl;
      cout << "*    dA_t(M_SUSY)/dm_Hd^2(M_X) = " << error_dParamsdMh1Sq(7) << endl;
      cout << "**************************************************" << endl;
      cout << "* SUBTEST 1 RESULT: ";
      if (error_dParamsdMh1Sq(0) > tol || error_dParamsdMh1Sq(1) > tol || error_dParamsdMh1Sq(2) > tol ||
	  error_dParamsdMh1Sq(3) > tol || error_dParamsdMh1Sq(4) > tol || error_dParamsdMh1Sq(5) > tol ||
	  error_dParamsdMh1Sq(6) > tol || error_dParamsdMh1Sq(7) > tol)
	{
	  cout << "FAIL" << endl;
	  if (error_dParamsdMh1Sq(0) > tol) cout << "*    maximum error exceeded for derivative of lambda(M_SUSY)" << endl;
	  if (error_dParamsdMh1Sq(1) > tol) cout << "*    maximum error exceeded for derivative of A_lambda(M_SUSY)" << endl;
	  if (error_dParamsdMh1Sq(2) > tol) cout << "*    maximum error exceeded for derivative of m_Hd^2(M_SUSY)" << endl;
	  if (error_dParamsdMh1Sq(3) > tol) cout << "*    maximum error exceeded for derivative of m_Hu^2(M_SUSY)" << endl;
	  if (error_dParamsdMh1Sq(4) > tol) cout << "*    maximum error exceeded for derivative of m_s^2(M_SUSY)" << endl;
	  if (error_dParamsdMh1Sq(5) > tol) cout << "*    maximum error exceeded for derivative of m_q3L^2(M_SUSY)" << endl;
	  if (error_dParamsdMh1Sq(6) > tol) cout << "*    maximum error exceeded for derivative of m_tR^2(M_SUSY)" << endl;
	  if (error_dParamsdMh1Sq(7) > tol) cout << "*    maximum error exceeded for derivative of A_t(M_SUSY)" << endl;
	  hasPassedMh1Sq = false;
	}
      else
	{
	  cout << "PASS" << endl; 
	  hasPassedMh1Sq = true;
	}
      cout << "**************************************************" << endl;
      cout << "* SUBTEST 2: m_Hu^2 derivatives" << endl; 
      cout << "**************************************************" << endl;
      cout << "* NUMERICAL DERIVATIVES: " << endl;
      cout << "*    dlambda(M_SUSY)/dm_Hu^2(M_X) = " << numeric_dParamsdMh2Sq(0) << endl;
      cout << "*    dA_lambda(M_SUSY)/dm_Hu^2(M_X) = " << numeric_dParamsdMh2Sq(1) << endl;
      cout << "*    dm_Hd^2(M_SUSY)/dm_Hu^2(M_X) = " << numeric_dParamsdMh2Sq(2) << endl;
      cout << "*    dm_Hu^2(M_SUSY)/dm_Hu^2(M_X) = " << numeric_dParamsdMh2Sq(3) << endl;
      cout << "*    dm_s^2(M_SUSY)/dm_Hu^2(M_X) = " << numeric_dParamsdMh2Sq(4) << endl;
      cout << "*    dm_q3L^2(M_SUSY)/dm_Hu^2(M_X) = " << numeric_dParamsdMh2Sq(5) << endl;
      cout << "*    dm_tR^2(M_SUSY)/dm_Hu^2(M_X) = " << numeric_dParamsdMh2Sq(6) << endl;
      cout << "*    dA_t(M_SUSY)/dm_Hu^2(M_X) = " << numeric_dParamsdMh2Sq(7) << endl;
      cout << "**************************************************" << endl;
      cout << "* ANALYTIC DERIVATIVES: " << endl;
      cout << "*    dlambda(M_SUSY)/dm_Hu^2(M_X) = " << analytic_dParamsdMh2Sq(0) << endl;
      cout << "*    dA_lambda(M_SUSY)/dm_Hu^2(M_X) = " << analytic_dParamsdMh2Sq(1) << endl;
      cout << "*    dm_Hd^2(M_SUSY)/dm_Hu^2(M_X) = " << analytic_dParamsdMh2Sq(2) << endl;
      cout << "*    dm_Hu^2(M_SUSY)/dm_Hu^2(M_X) = " << analytic_dParamsdMh2Sq(3) << endl;
      cout << "*    dm_s^2(M_SUSY)/dm_Hu^2(M_X) = " << analytic_dParamsdMh2Sq(4) << endl;
      cout << "*    dm_q3L^2(M_SUSY)/dm_Hu^2(M_X) = " << analytic_dParamsdMh2Sq(5) << endl;
      cout << "*    dm_tR^2(M_SUSY)/dm_Hu^2(M_X) = " << analytic_dParamsdMh2Sq(6) << endl;
      cout << "*    dA_t(M_SUSY)/dm_Hu^2(M_X) = " << analytic_dParamsdMh2Sq(7) << endl;
      cout << "**************************************************" << endl;
      cout << "* PERCENTAGE DIFFERENCES: " << endl;
      cout << "*    dlambda(M_SUSY)/dm_Hu^2(M_X) = " << error_dParamsdMh2Sq(0) << endl;
      cout << "*    dA_lambda(M_SUSY)/dm_Hu^2(M_X) = " << error_dParamsdMh2Sq(1) << endl;
      cout << "*    dm_Hd^2(M_SUSY)/dm_Hu^2(M_X) = " << error_dParamsdMh2Sq(2) << endl;
      cout << "*    dm_Hu^2(M_SUSY)/dm_Hu^2(M_X) = " << error_dParamsdMh2Sq(3) << endl;
      cout << "*    dm_s^2(M_SUSY)/dm_Hu^2(M_X) = " << error_dParamsdMh2Sq(4) << endl;
      cout << "*    dm_q3L^2(M_SUSY)/dm_Hu^2(M_X) = " << error_dParamsdMh2Sq(5) << endl;
      cout << "*    dm_tR^2(M_SUSY)/dm_Hu^2(M_X) = " << error_dParamsdMh2Sq(6) << endl;
      cout << "*    dA_t(M_SUSY)/dm_Hu^2(M_X) = " << error_dParamsdMh2Sq(7) << endl;
      cout << "**************************************************" << endl;
      cout << "* SUBTEST 2 RESULT: ";
      if (error_dParamsdMh2Sq(0) > tol || error_dParamsdMh2Sq(1) > tol || error_dParamsdMh2Sq(2) > tol ||
	  error_dParamsdMh2Sq(3) > tol || error_dParamsdMh2Sq(4) > tol || error_dParamsdMh2Sq(5) > tol ||
	  error_dParamsdMh2Sq(6) > tol || error_dParamsdMh2Sq(7) > tol)
	{
	  cout << "FAIL" << endl;
	  if (error_dParamsdMh2Sq(0) > tol) cout << "*    maximum error exceeded for derivative of lambda(M_SUSY)" << endl;
	  if (error_dParamsdMh2Sq(1) > tol) cout << "*    maximum error exceeded for derivative of A_lambda(M_SUSY)" << endl;
	  if (error_dParamsdMh2Sq(2) > tol) cout << "*    maximum error exceeded for derivative of m_Hd^2(M_SUSY)" << endl;
	  if (error_dParamsdMh2Sq(3) > tol) cout << "*    maximum error exceeded for derivative of m_Hu^2(M_SUSY)" << endl;
	  if (error_dParamsdMh2Sq(4) > tol) cout << "*    maximum error exceeded for derivative of m_s^2(M_SUSY)" << endl;
	  if (error_dParamsdMh2Sq(5) > tol) cout << "*    maximum error exceeded for derivative of m_q3L^2(M_SUSY)" << endl;
	  if (error_dParamsdMh2Sq(6) > tol) cout << "*    maximum error exceeded for derivative of m_tR^2(M_SUSY)" << endl;
	  if (error_dParamsdMh2Sq(7) > tol) cout << "*    maximum error exceeded for derivative of A_t(M_SUSY)" << endl;
	  hasPassedMh2Sq = false;
	}
      else
	{
	  cout << "PASS" << endl; 
	  hasPassedMh2Sq = true;
	}
      cout << "**************************************************" << endl;
      cout << "* SUBTEST 3: m_s^2 derivatives" << endl; 
      cout << "**************************************************" << endl;
      cout << "* NUMERICAL DERIVATIVES: " << endl;
      cout << "*    dlambda(M_SUSY)/dm_s^2(M_X) = " << numeric_dParamsdMsSq(0) << endl;
      cout << "*    dA_lambda(M_SUSY)/dm_s^2(M_X) = " << numeric_dParamsdMsSq(1) << endl;
      cout << "*    dm_Hd^2(M_SUSY)/dm_s^2(M_X) = " << numeric_dParamsdMsSq(2) << endl;
      cout << "*    dm_Hu^2(M_SUSY)/dm_s^2(M_X) = " << numeric_dParamsdMsSq(3) << endl;
      cout << "*    dm_s^2(M_SUSY)/dm_s^2(M_X) = " << numeric_dParamsdMsSq(4) << endl;
      cout << "*    dm_q3L^2(M_SUSY)/dm_s^2(M_X) = " << numeric_dParamsdMsSq(5) << endl;
      cout << "*    dm_tR^2(M_SUSY)/dm_s^2(M_X) = " << numeric_dParamsdMsSq(6) << endl;
      cout << "*    dA_t(M_SUSY)/dm_s^2(M_X) = " << numeric_dParamsdMsSq(7) << endl;
      cout << "**************************************************" << endl;
      cout << "* ANALYTIC DERIVATIVES: " << endl;
      cout << "*    dlambda(M_SUSY)/dm_s^2(M_X) = " << analytic_dParamsdMsSq(0) << endl;
      cout << "*    dA_lambda(M_SUSY)/dm_s^2(M_X) = " << analytic_dParamsdMsSq(1) << endl;
      cout << "*    dm_Hd^2(M_SUSY)/dm_s^2(M_X) = " << analytic_dParamsdMsSq(2) << endl;
      cout << "*    dm_Hu^2(M_SUSY)/dm_s^2(M_X) = " << analytic_dParamsdMsSq(3) << endl;
      cout << "*    dm_s^2(M_SUSY)/dm_s^2(M_X) = " << analytic_dParamsdMsSq(4) << endl;
      cout << "*    dm_q3L^2(M_SUSY)/dm_s^2(M_X) = " << analytic_dParamsdMsSq(5) << endl;
      cout << "*    dm_tR^2(M_SUSY)/dm_s^2(M_X) = " << analytic_dParamsdMsSq(6) << endl;
      cout << "*    dA_t(M_SUSY)/dm_s^2(M_X) = " << analytic_dParamsdMsSq(7) << endl;
      cout << "**************************************************" << endl;
      cout << "* PERCENTAGE DIFFERENCES: " << endl;
      cout << "*    dlambda(M_SUSY)/dm_s^2(M_X) = " << error_dParamsdMsSq(0) << endl;
      cout << "*    dA_lambda(M_SUSY)/dm_s^2(M_X) = " << error_dParamsdMsSq(1) << endl;
      cout << "*    dm_Hd^2(M_SUSY)/dm_s^2(M_X) = " << error_dParamsdMsSq(2) << endl;
      cout << "*    dm_Hu^2(M_SUSY)/dm_s^2(M_X) = " << error_dParamsdMsSq(3) << endl;
      cout << "*    dm_s^2(M_SUSY)/dm_s^2(M_X) = " << error_dParamsdMsSq(4) << endl;
      cout << "*    dm_q3L^2(M_SUSY)/dm_s^2(M_X) = " << error_dParamsdMsSq(5) << endl;
      cout << "*    dm_tR^2(M_SUSY)/dm_s^2(M_X) = " << error_dParamsdMsSq(6) << endl;
      cout << "*    dA_t(M_SUSY)/dm_s^2(M_X) = " << error_dParamsdMsSq(7) << endl;
      cout << "**************************************************" << endl;
      cout << "* SUBTEST 3 RESULT: ";
      if (error_dParamsdMsSq(0) > tol || error_dParamsdMsSq(1) > tol || error_dParamsdMsSq(2) > tol ||
	  error_dParamsdMsSq(3) > tol || error_dParamsdMsSq(4) > tol || error_dParamsdMsSq(5) > tol ||
	  error_dParamsdMsSq(6) > tol || error_dParamsdMsSq(7) > tol)
	{
	  cout << "FAIL" << endl;
	  if (error_dParamsdMsSq(0) > tol) cout << "*    maximum error exceeded for derivative of lambda(M_SUSY)" << endl;
	  if (error_dParamsdMsSq(1) > tol) cout << "*    maximum error exceeded for derivative of A_lambda(M_SUSY)" << endl;
	  if (error_dParamsdMsSq(2) > tol) cout << "*    maximum error exceeded for derivative of m_Hd^2(M_SUSY)" << endl;
	  if (error_dParamsdMsSq(3) > tol) cout << "*    maximum error exceeded for derivative of m_Hu^2(M_SUSY)" << endl;
	  if (error_dParamsdMsSq(4) > tol) cout << "*    maximum error exceeded for derivative of m_s^2(M_SUSY)" << endl;
	  if (error_dParamsdMsSq(5) > tol) cout << "*    maximum error exceeded for derivative of m_q3L^2(M_SUSY)" << endl;
	  if (error_dParamsdMsSq(6) > tol) cout << "*    maximum error exceeded for derivative of m_tR^2(M_SUSY)" << endl;
	  if (error_dParamsdMsSq(7) > tol) cout << "*    maximum error exceeded for derivative of A_t(M_SUSY)" << endl;
	  hasPassedMsSq = false;
	}
      else
	{
	  cout << "PASS" << endl; 
	  hasPassedMsSq = true;
	}
      cout << "**************************************************" << endl;
      cout << "* SUBTEST 4: m_q3L^2 derivatives" << endl; 
      cout << "**************************************************" << endl;
      cout << "* NUMERICAL DERIVATIVES: " << endl;
      cout << "*    dlambda(M_SUSY)/dm_q3L^2(M_X) = " << numeric_dParamsdMqSq(0) << endl;
      cout << "*    dA_lambda(M_SUSY)/dm_q3L^2(M_X) = " << numeric_dParamsdMqSq(1) << endl;
      cout << "*    dm_Hd^2(M_SUSY)/dm_q3L^2(M_X) = " << numeric_dParamsdMqSq(2) << endl;
      cout << "*    dm_Hu^2(M_SUSY)/dm_q3L^2(M_X) = " << numeric_dParamsdMqSq(3) << endl;
      cout << "*    dm_s^2(M_SUSY)/dm_q3L^2(M_X) = " << numeric_dParamsdMqSq(4) << endl;
      cout << "*    dm_q3L^2(M_SUSY)/dm_q3L^2(M_X) = " << numeric_dParamsdMqSq(5) << endl;
      cout << "*    dm_tR^2(M_SUSY)/dm_q3L^2(M_X) = " << numeric_dParamsdMqSq(6) << endl;
      cout << "*    dA_t(M_SUSY)/dm_q3L^2(M_X) = " << numeric_dParamsdMqSq(7) << endl;
      cout << "**************************************************" << endl;
      cout << "* ANALYTIC DERIVATIVES: " << endl;
      cout << "*    dlambda(M_SUSY)/dm_q3L^2(M_X) = " << analytic_dParamsdMqSq(0) << endl;
      cout << "*    dA_lambda(M_SUSY)/dm_q3L^2(M_X) = " << analytic_dParamsdMqSq(1) << endl;
      cout << "*    dm_Hd^2(M_SUSY)/dm_q3L^2(M_X) = " << analytic_dParamsdMqSq(2) << endl;
      cout << "*    dm_Hu^2(M_SUSY)/dm_q3L^2(M_X) = " << analytic_dParamsdMqSq(3) << endl;
      cout << "*    dm_s^2(M_SUSY)/dm_q3L^2(M_X) = " << analytic_dParamsdMqSq(4) << endl;
      cout << "*    dm_q3L^2(M_SUSY)/dm_q3L^2(M_X) = " << analytic_dParamsdMqSq(5) << endl;
      cout << "*    dm_tR^2(M_SUSY)/dm_q3L^2(M_X) = " << analytic_dParamsdMqSq(6) << endl;
      cout << "*    dA_t(M_SUSY)/dm_q3L^2(M_X) = " << analytic_dParamsdMqSq(7) << endl;
      cout << "**************************************************" << endl;
      cout << "* PERCENTAGE DIFFERENCES: " << endl;
      cout << "*    dlambda(M_SUSY)/dm_q3L^2(M_X) = " << error_dParamsdMqSq(0) << endl;
      cout << "*    dA_lambda(M_SUSY)/dm_q3L^2(M_X) = " << error_dParamsdMqSq(1) << endl;
      cout << "*    dm_Hd^2(M_SUSY)/dm_q3L^2(M_X) = " << error_dParamsdMqSq(2) << endl;
      cout << "*    dm_Hu^2(M_SUSY)/dm_q3L^2(M_X) = " << error_dParamsdMqSq(3) << endl;
      cout << "*    dm_s^2(M_SUSY)/dm_q3L^2(M_X) = " << error_dParamsdMqSq(4) << endl;
      cout << "*    dm_q3L^2(M_SUSY)/dm_q3L^2(M_X) = " << error_dParamsdMqSq(5) << endl;
      cout << "*    dm_tR^2(M_SUSY)/dm_q3L^2(M_X) = " << error_dParamsdMqSq(6) << endl;
      cout << "*    dA_t(M_SUSY)/dm_q3L^2(M_X) = " << error_dParamsdMqSq(7) << endl;
      cout << "**************************************************" << endl;
      cout << "* SUBTEST 4 RESULT: ";
      if (error_dParamsdMqSq(0) > tol || error_dParamsdMqSq(1) > tol || error_dParamsdMqSq(2) > tol ||
	  error_dParamsdMqSq(3) > tol || error_dParamsdMqSq(4) > tol || error_dParamsdMqSq(5) > tol ||
	  error_dParamsdMqSq(6) > tol || error_dParamsdMqSq(7) > tol)
	{
	  cout << "FAIL" << endl;
	  if (error_dParamsdMqSq(0) > tol) cout << "*    maximum error exceeded for derivative of lambda(M_SUSY)" << endl;
	  if (error_dParamsdMqSq(1) > tol) cout << "*    maximum error exceeded for derivative of A_lambda(M_SUSY)" << endl;
	  if (error_dParamsdMqSq(2) > tol) cout << "*    maximum error exceeded for derivative of m_Hd^2(M_SUSY)" << endl;
	  if (error_dParamsdMqSq(3) > tol) cout << "*    maximum error exceeded for derivative of m_Hu^2(M_SUSY)" << endl;
	  if (error_dParamsdMqSq(4) > tol) cout << "*    maximum error exceeded for derivative of m_s^2(M_SUSY)" << endl;
	  if (error_dParamsdMqSq(5) > tol) cout << "*    maximum error exceeded for derivative of m_q3L^2(M_SUSY)" << endl;
	  if (error_dParamsdMqSq(6) > tol) cout << "*    maximum error exceeded for derivative of m_tR^2(M_SUSY)" << endl;
	  if (error_dParamsdMqSq(7) > tol) cout << "*    maximum error exceeded for derivative of A_t(M_SUSY)" << endl;
	  hasPassedMqSq = false;
	}
      else
	{
	  cout << "PASS" << endl; 
	  hasPassedMqSq = true;
	}
      cout << "**************************************************" << endl;
      cout << "* SUBTEST 5: m_tR^2 derivatives" << endl; 
      cout << "**************************************************" << endl;
      cout << "* NUMERICAL DERIVATIVES: " << endl;
      cout << "*    dlambda(M_SUSY)/dm_tR^2(M_X) = " << numeric_dParamsdMuSq(0) << endl;
      cout << "*    dA_lambda(M_SUSY)/dm_tR^2(M_X) = " << numeric_dParamsdMuSq(1) << endl;
      cout << "*    dm_Hd^2(M_SUSY)/dm_tR^2(M_X) = " << numeric_dParamsdMuSq(2) << endl;
      cout << "*    dm_Hu^2(M_SUSY)/dm_tR^2(M_X) = " << numeric_dParamsdMuSq(3) << endl;
      cout << "*    dm_s^2(M_SUSY)/dm_tR^2(M_X) = " << numeric_dParamsdMuSq(4) << endl;
      cout << "*    dm_q3L^2(M_SUSY)/dm_tR^2(M_X) = " << numeric_dParamsdMuSq(5) << endl;
      cout << "*    dm_tR^2(M_SUSY)/dm_tR^2(M_X) = " << numeric_dParamsdMuSq(6) << endl;
      cout << "*    dA_t(M_SUSY)/dm_tR^2(M_X) = " << numeric_dParamsdMuSq(7) << endl;
      cout << "**************************************************" << endl;
      cout << "* ANALYTIC DERIVATIVES: " << endl;
      cout << "*    dlambda(M_SUSY)/dm_tR^2(M_X) = " << analytic_dParamsdMuSq(0) << endl;
      cout << "*    dA_lambda(M_SUSY)/dm_tR^2(M_X) = " << analytic_dParamsdMuSq(1) << endl;
      cout << "*    dm_Hd^2(M_SUSY)/dm_tR^2(M_X) = " << analytic_dParamsdMuSq(2) << endl;
      cout << "*    dm_Hu^2(M_SUSY)/dm_tR^2(M_X) = " << analytic_dParamsdMuSq(3) << endl;
      cout << "*    dm_s^2(M_SUSY)/dm_tR^2(M_X) = " << analytic_dParamsdMuSq(4) << endl;
      cout << "*    dm_q3L^2(M_SUSY)/dm_tR^2(M_X) = " << analytic_dParamsdMuSq(5) << endl;
      cout << "*    dm_tR^2(M_SUSY)/dm_tR^2(M_X) = " << analytic_dParamsdMuSq(6) << endl;
      cout << "*    dA_t(M_SUSY)/dm_tR^2(M_X) = " << analytic_dParamsdMuSq(7) << endl;
      cout << "**************************************************" << endl;
      cout << "* PERCENTAGE DIFFERENCES: " << endl;
      cout << "*    dlambda(M_SUSY)/dm_tR^2(M_X) = " << error_dParamsdMuSq(0) << endl;
      cout << "*    dA_lambda(M_SUSY)/dm_tR^2(M_X) = " << error_dParamsdMuSq(1) << endl;
      cout << "*    dm_Hd^2(M_SUSY)/dm_tR^2(M_X) = " << error_dParamsdMuSq(2) << endl;
      cout << "*    dm_Hu^2(M_SUSY)/dm_tR^2(M_X) = " << error_dParamsdMuSq(3) << endl;
      cout << "*    dm_s^2(M_SUSY)/dm_tR^2(M_X) = " << error_dParamsdMuSq(4) << endl;
      cout << "*    dm_q3L^2(M_SUSY)/dm_tR^2(M_X) = " << error_dParamsdMuSq(5) << endl;
      cout << "*    dm_tR^2(M_SUSY)/dm_tR^2(M_X) = " << error_dParamsdMuSq(6) << endl;
      cout << "*    dA_t(M_SUSY)/dm_tR^2(M_X) = " << error_dParamsdMuSq(7) << endl;
      cout << "**************************************************" << endl;
      cout << "* SUBTEST 5 RESULT: ";
      if (error_dParamsdMuSq(0) > tol || error_dParamsdMuSq(1) > tol || error_dParamsdMuSq(2) > tol ||
	  error_dParamsdMuSq(3) > tol || error_dParamsdMuSq(4) > tol || error_dParamsdMuSq(5) > tol ||
	  error_dParamsdMuSq(6) > tol || error_dParamsdMuSq(7) > tol)
	{
	  cout << "FAIL" << endl;
	  if (error_dParamsdMuSq(0) > tol) cout << "*    maximum error exceeded for derivative of lambda(M_SUSY)" << endl;
	  if (error_dParamsdMuSq(1) > tol) cout << "*    maximum error exceeded for derivative of A_lambda(M_SUSY)" << endl;
	  if (error_dParamsdMuSq(2) > tol) cout << "*    maximum error exceeded for derivative of m_Hd^2(M_SUSY)" << endl;
	  if (error_dParamsdMuSq(3) > tol) cout << "*    maximum error exceeded for derivative of m_Hu^2(M_SUSY)" << endl;
	  if (error_dParamsdMuSq(4) > tol) cout << "*    maximum error exceeded for derivative of m_s^2(M_SUSY)" << endl;
	  if (error_dParamsdMuSq(5) > tol) cout << "*    maximum error exceeded for derivative of m_q3L^2(M_SUSY)" << endl;
	  if (error_dParamsdMuSq(6) > tol) cout << "*    maximum error exceeded for derivative of m_tR^2(M_SUSY)" << endl;
	  if (error_dParamsdMuSq(7) > tol) cout << "*    maximum error exceeded for derivative of A_t(M_SUSY)" << endl;
	  hasPassedMuSq = false;
	}
      else
	{
	  cout << "PASS" << endl; 
	  hasPassedMuSq = true;
	}
      cout << "**************************************************" << endl;
      cout << "* SUBTEST 6: lambda derivatives" << endl; 
      cout << "**************************************************" << endl;
      cout << "* NUMERICAL DERIVATIVES: " << endl;
      cout << "*    dlambda(M_SUSY)/dlambda(M_X) = " << numeric_dParamsdLambda(0) << endl;
      cout << "*    dA_lambda(M_SUSY)/dlambda(M_X) = " << numeric_dParamsdLambda(1) << endl;
      cout << "*    dm_Hd^2(M_SUSY)/dlambda(M_X) = " << numeric_dParamsdLambda(2) << endl;
      cout << "*    dm_Hu^2(M_SUSY)/dlambda(M_X) = " << numeric_dParamsdLambda(3) << endl;
      cout << "*    dm_s^2(M_SUSY)/dlambda(M_X) = " << numeric_dParamsdLambda(4) << endl;
      cout << "*    dm_q3L^2(M_SUSY)/dlambda(M_X) = " << numeric_dParamsdLambda(5) << endl;
      cout << "*    dm_tR^2(M_SUSY)/dlambda(M_X) = " << numeric_dParamsdLambda(6) << endl;
      cout << "*    dA_t(M_SUSY)/dlambda(M_X) = " << numeric_dParamsdLambda(7) << endl;
      cout << "**************************************************" << endl;
      cout << "* ANALYTIC DERIVATIVES: " << endl;
      cout << "*    dlambda(M_SUSY)/dlambda(M_X) = " << analytic_dParamsdLambda(0) << endl;
      cout << "*    dA_lambda(M_SUSY)/dlambda(M_X) = " << analytic_dParamsdLambda(1) << endl;
      cout << "*    dm_Hd^2(M_SUSY)/dlambda(M_X) = " << analytic_dParamsdLambda(2) << endl;
      cout << "*    dm_Hu^2(M_SUSY)/dlambda(M_X) = " << analytic_dParamsdLambda(3) << endl;
      cout << "*    dm_s^2(M_SUSY)/dlambda(M_X) = " << analytic_dParamsdLambda(4) << endl;
      cout << "*    dm_q3L^2(M_SUSY)/dlambda(M_X) = " << analytic_dParamsdLambda(5) << endl;
      cout << "*    dm_tR^2(M_SUSY)/dlambda(M_X) = " << analytic_dParamsdLambda(6) << endl;
      cout << "*    dA_t(M_SUSY)/dlambda(M_X) = " << analytic_dParamsdLambda(7) << endl;
      cout << "**************************************************" << endl;
      cout << "* PERCENTAGE DIFFERENCES: " << endl;
      cout << "*    dlambda(M_SUSY)/dlambda(M_X) = " << error_dParamsdLambda(0) << endl;
      cout << "*    dA_lambda(M_SUSY)/dlambda(M_X) = " << error_dParamsdLambda(1) << endl;
      cout << "*    dm_Hd^2(M_SUSY)/dlambda(M_X) = " << error_dParamsdLambda(2) << endl;
      cout << "*    dm_Hu^2(M_SUSY)/dlambda(M_X) = " << error_dParamsdLambda(3) << endl;
      cout << "*    dm_s^2(M_SUSY)/dlambda(M_X) = " << error_dParamsdLambda(4) << endl;
      cout << "*    dm_q3L^2(M_SUSY)/dlambda(M_X) = " << error_dParamsdLambda(5) << endl;
      cout << "*    dm_tR^2(M_SUSY)/dlambda(M_X) = " << error_dParamsdLambda(6) << endl;
      cout << "*    dA_t(M_SUSY)/dlambda(M_X) = " << error_dParamsdLambda(7) << endl;
      cout << "**************************************************" << endl;
      cout << "* SUBTEST 6 RESULT: ";
      if (error_dParamsdLambda(0) > tol || error_dParamsdLambda(1) > tol || error_dParamsdLambda(2) > tol ||
	  error_dParamsdLambda(3) > tol || error_dParamsdLambda(4) > tol || error_dParamsdLambda(5) > tol ||
	  error_dParamsdLambda(6) > tol || error_dParamsdLambda(7) > tol)
	{
	  cout << "FAIL" << endl;
	  if (error_dParamsdLambda(0) > tol) cout << "*    maximum error exceeded for derivative of lambda(M_SUSY)" << endl;
	  if (error_dParamsdLambda(1) > tol) cout << "*    maximum error exceeded for derivative of A_lambda(M_SUSY)" << endl;
	  if (error_dParamsdLambda(2) > tol) cout << "*    maximum error exceeded for derivative of m_Hd^2(M_SUSY)" << endl;
	  if (error_dParamsdLambda(3) > tol) cout << "*    maximum error exceeded for derivative of m_Hu^2(M_SUSY)" << endl;
	  if (error_dParamsdLambda(4) > tol) cout << "*    maximum error exceeded for derivative of m_s^2(M_SUSY)" << endl;
	  if (error_dParamsdLambda(5) > tol) cout << "*    maximum error exceeded for derivative of m_q3L^2(M_SUSY)" << endl;
	  if (error_dParamsdLambda(6) > tol) cout << "*    maximum error exceeded for derivative of m_tR^2(M_SUSY)" << endl;
	  if (error_dParamsdLambda(7) > tol) cout << "*    maximum error exceeded for derivative of A_t(M_SUSY)" << endl;
	  hasPassedLambda = false;
	}
      else
	{
	  cout << "PASS" << endl; 
	  hasPassedLambda = true;
	}
      cout << "**************************************************" << endl;
      cout << "* SUBTEST 7: A_lambda derivatives" << endl; 
      cout << "**************************************************" << endl;
      cout << "* NUMERICAL DERIVATIVES: " << endl;
      cout << "*    dlambda(M_SUSY)/dA_lambda(M_X) = " << numeric_dParamsdAlambda(0) << endl;
      cout << "*    dA_lambda(M_SUSY)/dA_lambda(M_X) = " << numeric_dParamsdAlambda(1) << endl;
      cout << "*    dm_Hd^2(M_SUSY)/dA_lambda(M_X) = " << numeric_dParamsdAlambda(2) << endl;
      cout << "*    dm_Hu^2(M_SUSY)/dA_lambda(M_X) = " << numeric_dParamsdAlambda(3) << endl;
      cout << "*    dm_s^2(M_SUSY)/dA_lambda(M_X) = " << numeric_dParamsdAlambda(4) << endl;
      cout << "*    dm_q3L^2(M_SUSY)/dA_lambda(M_X) = " << numeric_dParamsdAlambda(5) << endl;
      cout << "*    dm_tR^2(M_SUSY)/dA_lambda(M_X) = " << numeric_dParamsdAlambda(6) << endl;
      cout << "*    dA_t(M_SUSY)/dA_lambda(M_X) = " << numeric_dParamsdAlambda(7) << endl;
      cout << "**************************************************" << endl;
      cout << "* ANALYTIC DERIVATIVES: " << endl;
      cout << "*    dlambda(M_SUSY)/dA_lambda(M_X) = " << analytic_dParamsdAlambda(0) << endl;
      cout << "*    dA_lambda(M_SUSY)/dA_lambda(M_X) = " << analytic_dParamsdAlambda(1) << endl;
      cout << "*    dm_Hd^2(M_SUSY)/dA_lambda(M_X) = " << analytic_dParamsdAlambda(2) << endl;
      cout << "*    dm_Hu^2(M_SUSY)/dA_lambda(M_X) = " << analytic_dParamsdAlambda(3) << endl;
      cout << "*    dm_s^2(M_SUSY)/dA_lambda(M_X) = " << analytic_dParamsdAlambda(4) << endl;
      cout << "*    dm_q3L^2(M_SUSY)/dA_lambda(M_X) = " << analytic_dParamsdAlambda(5) << endl;
      cout << "*    dm_tR^2(M_SUSY)/dA_lambda(M_X) = " << analytic_dParamsdAlambda(6) << endl;
      cout << "*    dA_t(M_SUSY)/dA_lambda(M_X) = " << analytic_dParamsdAlambda(7) << endl;
      cout << "**************************************************" << endl;
      cout << "* PERCENTAGE DIFFERENCES: " << endl;
      cout << "*    dlambda(M_SUSY)/dA_lambda(M_X) = " << error_dParamsdAlambda(0) << endl;
      cout << "*    dA_lambda(M_SUSY)/dA_lambda(M_X) = " << error_dParamsdAlambda(1) << endl;
      cout << "*    dm_Hd^2(M_SUSY)/dA_lambda(M_X) = " << error_dParamsdAlambda(2) << endl;
      cout << "*    dm_Hu^2(M_SUSY)/dA_lambda(M_X) = " << error_dParamsdAlambda(3) << endl;
      cout << "*    dm_s^2(M_SUSY)/dA_lambda(M_X) = " << error_dParamsdAlambda(4) << endl;
      cout << "*    dm_q3L^2(M_SUSY)/dA_lambda(M_X) = " << error_dParamsdAlambda(5) << endl;
      cout << "*    dm_tR^2(M_SUSY)/dA_lambda(M_X) = " << error_dParamsdAlambda(6) << endl;
      cout << "*    dA_t(M_SUSY)/dA_lambda(M_X) = " << error_dParamsdAlambda(7) << endl;
      cout << "**************************************************" << endl;
      cout << "* SUBTEST 7 RESULT: ";
      if (error_dParamsdAlambda(0) > tol || error_dParamsdAlambda(1) > tol || error_dParamsdAlambda(2) > tol ||
	  error_dParamsdAlambda(3) > tol || error_dParamsdAlambda(4) > tol || error_dParamsdAlambda(5) > tol ||
	  error_dParamsdAlambda(6) > tol || error_dParamsdAlambda(7) > tol)
	{
	  cout << "FAIL" << endl;
	  if (error_dParamsdAlambda(0) > tol) cout << "*    maximum error exceeded for derivative of lambda(M_SUSY)" << endl;
	  if (error_dParamsdAlambda(1) > tol) cout << "*    maximum error exceeded for derivative of A_lambda(M_SUSY)" << endl;
	  if (error_dParamsdAlambda(2) > tol) cout << "*    maximum error exceeded for derivative of m_Hd^2(M_SUSY)" << endl;
	  if (error_dParamsdAlambda(3) > tol) cout << "*    maximum error exceeded for derivative of m_Hu^2(M_SUSY)" << endl;
	  if (error_dParamsdAlambda(4) > tol) cout << "*    maximum error exceeded for derivative of m_s^2(M_SUSY)" << endl;
	  if (error_dParamsdAlambda(5) > tol) cout << "*    maximum error exceeded for derivative of m_q3L^2(M_SUSY)" << endl;
	  if (error_dParamsdAlambda(6) > tol) cout << "*    maximum error exceeded for derivative of m_tR^2(M_SUSY)" << endl;
	  if (error_dParamsdAlambda(7) > tol) cout << "*    maximum error exceeded for derivative of A_t(M_SUSY)" << endl;
	  hasPassedAlambda = false;
	}
      else
	{
	  cout << "PASS" << endl; 
	  hasPassedAlambda = true;
	}
      cout << "**************************************************" << endl;
      cout << "* SUBTEST 8: A_t derivatives" << endl; 
      cout << "**************************************************" << endl;
      cout << "* NUMERICAL DERIVATIVES: " << endl;
      cout << "*    dlambda(M_SUSY)/dA_t(M_X) = " << numeric_dParamsdAt(0) << endl;
      cout << "*    dA_lambda(M_SUSY)/dA_t(M_X) = " << numeric_dParamsdAt(1) << endl;
      cout << "*    dm_Hd^2(M_SUSY)/dA_t(M_X) = " << numeric_dParamsdAt(2) << endl;
      cout << "*    dm_Hu^2(M_SUSY)/dA_t(M_X) = " << numeric_dParamsdAt(3) << endl;
      cout << "*    dm_s^2(M_SUSY)/dA_t(M_X) = " << numeric_dParamsdAt(4) << endl;
      cout << "*    dm_q3L^2(M_SUSY)/dA_t(M_X) = " << numeric_dParamsdAt(5) << endl;
      cout << "*    dm_tR^2(M_SUSY)/dA_t(M_X) = " << numeric_dParamsdAt(6) << endl;
      cout << "*    dA_t(M_SUSY)/dA_t(M_X) = " << numeric_dParamsdAt(7) << endl;
      cout << "**************************************************" << endl;
      cout << "* ANALYTIC DERIVATIVES: " << endl;
      cout << "*    dlambda(M_SUSY)/dA_t(M_X) = " << analytic_dParamsdAt(0) << endl;
      cout << "*    dA_lambda(M_SUSY)/dA_t(M_X) = " << analytic_dParamsdAt(1) << endl;
      cout << "*    dm_Hd^2(M_SUSY)/dA_t(M_X) = " << analytic_dParamsdAt(2) << endl;
      cout << "*    dm_Hu^2(M_SUSY)/dA_t(M_X) = " << analytic_dParamsdAt(3) << endl;
      cout << "*    dm_s^2(M_SUSY)/dA_t(M_X) = " << analytic_dParamsdAt(4) << endl;
      cout << "*    dm_q3L^2(M_SUSY)/dA_t(M_X) = " << analytic_dParamsdAt(5) << endl;
      cout << "*    dm_tR^2(M_SUSY)/dA_t(M_X) = " << analytic_dParamsdAt(6) << endl;
      cout << "*    dA_t(M_SUSY)/dA_t(M_X) = " << analytic_dParamsdAt(7) << endl;
      cout << "**************************************************" << endl;
      cout << "* PERCENTAGE DIFFERENCES: " << endl;
      cout << "*    dlambda(M_SUSY)/dA_t(M_X) = " << error_dParamsdAt(0) << endl;
      cout << "*    dA_lambda(M_SUSY)/dA_t(M_X) = " << error_dParamsdAt(1) << endl;
      cout << "*    dm_Hd^2(M_SUSY)/dA_t(M_X) = " << error_dParamsdAt(2) << endl;
      cout << "*    dm_Hu^2(M_SUSY)/dA_t(M_X) = " << error_dParamsdAt(3) << endl;
      cout << "*    dm_s^2(M_SUSY)/dA_t(M_X) = " << error_dParamsdAt(4) << endl;
      cout << "*    dm_q3L^2(M_SUSY)/dA_t(M_X) = " << error_dParamsdAt(5) << endl;
      cout << "*    dm_tR^2(M_SUSY)/dA_t(M_X) = " << error_dParamsdAt(6) << endl;
      cout << "*    dA_t(M_SUSY)/dA_t(M_X) = " << error_dParamsdAt(7) << endl;
      cout << "**************************************************" << endl;
      cout << "* SUBTEST 8 RESULT: ";
      if (error_dParamsdAt(0) > tol || error_dParamsdAt(1) > tol || error_dParamsdAt(2) > tol ||
	  error_dParamsdAt(3) > tol || error_dParamsdAt(4) > tol || error_dParamsdAt(5) > tol ||
	  error_dParamsdAt(6) > tol || error_dParamsdAt(7) > tol)
	{
	  cout << "FAIL" << endl;
	  if (error_dParamsdAt(0) > tol) cout << "*    maximum error exceeded for derivative of lambda(M_SUSY)" << endl;
	  if (error_dParamsdAt(1) > tol) cout << "*    maximum error exceeded for derivative of A_lambda(M_SUSY)" << endl;
	  if (error_dParamsdAt(2) > tol) cout << "*    maximum error exceeded for derivative of m_Hd^2(M_SUSY)" << endl;
	  if (error_dParamsdAt(3) > tol) cout << "*    maximum error exceeded for derivative of m_Hu^2(M_SUSY)" << endl;
	  if (error_dParamsdAt(4) > tol) cout << "*    maximum error exceeded for derivative of m_s^2(M_SUSY)" << endl;
	  if (error_dParamsdAt(5) > tol) cout << "*    maximum error exceeded for derivative of m_q3L^2(M_SUSY)" << endl;
	  if (error_dParamsdAt(6) > tol) cout << "*    maximum error exceeded for derivative of m_tR^2(M_SUSY)" << endl;
	  if (error_dParamsdAt(7) > tol) cout << "*    maximum error exceeded for derivative of A_t(M_SUSY)" << endl;
	  hasPassedAt = false;
	}
      else
	{
	  cout << "PASS" << endl; 
	  hasPassedAt = true;
	}
      cout << "**************************************************" << endl;
      cout << "* SUBTEST 9: M_1 derivatives" << endl; 
      cout << "**************************************************" << endl;
      cout << "* NUMERICAL DERIVATIVES: " << endl;
      cout << "*    dlambda(M_SUSY)/dM_1(M_X) = " << numeric_dParamsdM1(0) << endl;
      cout << "*    dA_lambda(M_SUSY)/dM_1(M_X) = " << numeric_dParamsdM1(1) << endl;
      cout << "*    dm_Hd^2(M_SUSY)/dM_1(M_X) = " << numeric_dParamsdM1(2) << endl;
      cout << "*    dm_Hu^2(M_SUSY)/dM_1(M_X) = " << numeric_dParamsdM1(3) << endl;
      cout << "*    dm_s^2(M_SUSY)/dM_1(M_X) = " << numeric_dParamsdM1(4) << endl;
      cout << "*    dm_q3L^2(M_SUSY)/dM_1(M_X) = " << numeric_dParamsdM1(5) << endl;
      cout << "*    dm_tR^2(M_SUSY)/dM_1(M_X) = " << numeric_dParamsdM1(6) << endl;
      cout << "*    dA_t(M_SUSY)/dM_1(M_X) = " << numeric_dParamsdM1(7) << endl;
      cout << "**************************************************" << endl;
      cout << "* ANALYTIC DERIVATIVES: " << endl;
      cout << "*    dlambda(M_SUSY)/dM_1(M_X) = " << analytic_dParamsdM1(0) << endl;
      cout << "*    dA_lambda(M_SUSY)/dM_1(M_X) = " << analytic_dParamsdM1(1) << endl;
      cout << "*    dm_Hd^2(M_SUSY)/dM_1(M_X) = " << analytic_dParamsdM1(2) << endl;
      cout << "*    dm_Hu^2(M_SUSY)/dM_1(M_X) = " << analytic_dParamsdM1(3) << endl;
      cout << "*    dm_s^2(M_SUSY)/dM_1(M_X) = " << analytic_dParamsdM1(4) << endl;
      cout << "*    dm_q3L^2(M_SUSY)/dM_1(M_X) = " << analytic_dParamsdM1(5) << endl;
      cout << "*    dm_tR^2(M_SUSY)/dM_1(M_X) = " << analytic_dParamsdM1(6) << endl;
      cout << "*    dA_t(M_SUSY)/dM_1(M_X) = " << analytic_dParamsdM1(7) << endl;
      cout << "**************************************************" << endl;
      cout << "* PERCENTAGE DIFFERENCES: " << endl;
      cout << "*    dlambda(M_SUSY)/dM_1(M_X) = " << error_dParamsdM1(0) << endl;
      cout << "*    dA_lambda(M_SUSY)/dM_1(M_X) = " << error_dParamsdM1(1) << endl;
      cout << "*    dm_Hd^2(M_SUSY)/dM_1(M_X) = " << error_dParamsdM1(2) << endl;
      cout << "*    dm_Hu^2(M_SUSY)/dM_1(M_X) = " << error_dParamsdM1(3) << endl;
      cout << "*    dm_s^2(M_SUSY)/dM_1(M_X) = " << error_dParamsdM1(4) << endl;
      cout << "*    dm_q3L^2(M_SUSY)/dM_1(M_X) = " << error_dParamsdM1(5) << endl;
      cout << "*    dm_tR^2(M_SUSY)/dM_1(M_X) = " << error_dParamsdM1(6) << endl;
      cout << "*    dA_t(M_SUSY)/dM_1(M_X) = " << error_dParamsdM1(7) << endl;
      cout << "**************************************************" << endl;
      cout << "* SUBTEST 9 RESULT: ";
      if (error_dParamsdM1(0) > tol || error_dParamsdM1(1) > tol || error_dParamsdM1(2) > tol ||
	  error_dParamsdM1(3) > tol || error_dParamsdM1(4) > tol || error_dParamsdM1(5) > tol ||
	  error_dParamsdM1(6) > tol || error_dParamsdM1(7) > tol)
	{
	  cout << "FAIL" << endl;
	  if (error_dParamsdM1(0) > tol) cout << "*    maximum error exceeded for derivative of lambda(M_SUSY)" << endl;
	  if (error_dParamsdM1(1) > tol) cout << "*    maximum error exceeded for derivative of A_lambda(M_SUSY)" << endl;
	  if (error_dParamsdM1(2) > tol) cout << "*    maximum error exceeded for derivative of m_Hd^2(M_SUSY)" << endl;
	  if (error_dParamsdM1(3) > tol) cout << "*    maximum error exceeded for derivative of m_Hu^2(M_SUSY)" << endl;
	  if (error_dParamsdM1(4) > tol) cout << "*    maximum error exceeded for derivative of m_s^2(M_SUSY)" << endl;
	  if (error_dParamsdM1(5) > tol) cout << "*    maximum error exceeded for derivative of m_q3L^2(M_SUSY)" << endl;
	  if (error_dParamsdM1(6) > tol) cout << "*    maximum error exceeded for derivative of m_tR^2(M_SUSY)" << endl;
	  if (error_dParamsdM1(7) > tol) cout << "*    maximum error exceeded for derivative of A_t(M_SUSY)" << endl;
	  hasPassedM1 = false;
	}
      else
	{
	  cout << "PASS" << endl; 
	  hasPassedM1 = true;
	}
      cout << "**************************************************" << endl;
      cout << "* SUBTEST 10: M_2 derivatives" << endl; 
      cout << "**************************************************" << endl;
      cout << "* NUMERICAL DERIVATIVES: " << endl;
      cout << "*    dlambda(M_SUSY)/dM_2(M_X) = " << numeric_dParamsdM2(0) << endl;
      cout << "*    dA_lambda(M_SUSY)/dM_2(M_X) = " << numeric_dParamsdM2(1) << endl;
      cout << "*    dm_Hd^2(M_SUSY)/dM_2(M_X) = " << numeric_dParamsdM2(2) << endl;
      cout << "*    dm_Hu^2(M_SUSY)/dM_2(M_X) = " << numeric_dParamsdM2(3) << endl;
      cout << "*    dm_s^2(M_SUSY)/dM_2(M_X) = " << numeric_dParamsdM2(4) << endl;
      cout << "*    dm_q3L^2(M_SUSY)/dM_2(M_X) = " << numeric_dParamsdM2(5) << endl;
      cout << "*    dm_tR^2(M_SUSY)/dM_2(M_X) = " << numeric_dParamsdM2(6) << endl;
      cout << "*    dA_t(M_SUSY)/dM_2(M_X) = " << numeric_dParamsdM2(7) << endl;
      cout << "**************************************************" << endl;
      cout << "* ANALYTIC DERIVATIVES: " << endl;
      cout << "*    dlambda(M_SUSY)/dM_2(M_X) = " << analytic_dParamsdM2(0) << endl;
      cout << "*    dA_lambda(M_SUSY)/dM_2(M_X) = " << analytic_dParamsdM2(1) << endl;
      cout << "*    dm_Hd^2(M_SUSY)/dM_2(M_X) = " << analytic_dParamsdM2(2) << endl;
      cout << "*    dm_Hu^2(M_SUSY)/dM_2(M_X) = " << analytic_dParamsdM2(3) << endl;
      cout << "*    dm_s^2(M_SUSY)/dM_2(M_X) = " << analytic_dParamsdM2(4) << endl;
      cout << "*    dm_q3L^2(M_SUSY)/dM_2(M_X) = " << analytic_dParamsdM2(5) << endl;
      cout << "*    dm_tR^2(M_SUSY)/dM_2(M_X) = " << analytic_dParamsdM2(6) << endl;
      cout << "*    dA_t(M_SUSY)/dM_2(M_X) = " << analytic_dParamsdM2(7) << endl;
      cout << "**************************************************" << endl;
      cout << "* PERCENTAGE DIFFERENCES: " << endl;
      cout << "*    dlambda(M_SUSY)/dM_2(M_X) = " << error_dParamsdM2(0) << endl;
      cout << "*    dA_lambda(M_SUSY)/dM_2(M_X) = " << error_dParamsdM2(1) << endl;
      cout << "*    dm_Hd^2(M_SUSY)/dM_2(M_X) = " << error_dParamsdM2(2) << endl;
      cout << "*    dm_Hu^2(M_SUSY)/dM_2(M_X) = " << error_dParamsdM2(3) << endl;
      cout << "*    dm_s^2(M_SUSY)/dM_2(M_X) = " << error_dParamsdM2(4) << endl;
      cout << "*    dm_q3L^2(M_SUSY)/dM_2(M_X) = " << error_dParamsdM2(5) << endl;
      cout << "*    dm_tR^2(M_SUSY)/dM_2(M_X) = " << error_dParamsdM2(6) << endl;
      cout << "*    dA_t(M_SUSY)/dM_2(M_X) = " << error_dParamsdM2(7) << endl;
      cout << "**************************************************" << endl;
      cout << "* SUBTEST 10 RESULT: ";
      if (error_dParamsdM2(0) > tol || error_dParamsdM2(1) > tol || error_dParamsdM2(2) > tol ||
	  error_dParamsdM2(3) > tol || error_dParamsdM2(4) > tol || error_dParamsdM2(5) > tol ||
	  error_dParamsdM2(6) > tol || error_dParamsdM2(7) > tol)
	{
	  cout << "FAIL" << endl;
	  if (error_dParamsdM2(0) > tol) cout << "*    maximum error exceeded for derivative of lambda(M_SUSY)" << endl;
	  if (error_dParamsdM2(1) > tol) cout << "*    maximum error exceeded for derivative of A_lambda(M_SUSY)" << endl;
	  if (error_dParamsdM2(2) > tol) cout << "*    maximum error exceeded for derivative of m_Hd^2(M_SUSY)" << endl;
	  if (error_dParamsdM2(3) > tol) cout << "*    maximum error exceeded for derivative of m_Hu^2(M_SUSY)" << endl;
	  if (error_dParamsdM2(4) > tol) cout << "*    maximum error exceeded for derivative of m_s^2(M_SUSY)" << endl;
	  if (error_dParamsdM2(5) > tol) cout << "*    maximum error exceeded for derivative of m_q3L^2(M_SUSY)" << endl;
	  if (error_dParamsdM2(6) > tol) cout << "*    maximum error exceeded for derivative of m_tR^2(M_SUSY)" << endl;
	  if (error_dParamsdM2(7) > tol) cout << "*    maximum error exceeded for derivative of A_t(M_SUSY)" << endl;
	  hasPassedM2 = false;
	}
      else
	{
	  cout << "PASS" << endl; 
	  hasPassedM2 = true;
	}
      cout << "**************************************************" << endl;
      cout << "* SUBTEST 11: M_3 derivatives" << endl; 
      cout << "**************************************************" << endl;
      cout << "* NUMERICAL DERIVATIVES: " << endl;
      cout << "*    dlambda(M_SUSY)/dM_3(M_X) = " << numeric_dParamsdM3(0) << endl;
      cout << "*    dA_lambda(M_SUSY)/dM_3(M_X) = " << numeric_dParamsdM3(1) << endl;
      cout << "*    dm_Hd^2(M_SUSY)/dM_3(M_X) = " << numeric_dParamsdM3(2) << endl;
      cout << "*    dm_Hu^2(M_SUSY)/dM_3(M_X) = " << numeric_dParamsdM3(3) << endl;
      cout << "*    dm_s^2(M_SUSY)/dM_3(M_X) = " << numeric_dParamsdM3(4) << endl;
      cout << "*    dm_q3L^2(M_SUSY)/dM_3(M_X) = " << numeric_dParamsdM3(5) << endl;
      cout << "*    dm_tR^2(M_SUSY)/dM_3(M_X) = " << numeric_dParamsdM3(6) << endl;
      cout << "*    dA_t(M_SUSY)/dM_3(M_X) = " << numeric_dParamsdM3(7) << endl;
      cout << "**************************************************" << endl;
      cout << "* ANALYTIC DERIVATIVES: " << endl;
      cout << "*    dlambda(M_SUSY)/dM_3(M_X) = " << analytic_dParamsdM3(0) << endl;
      cout << "*    dA_lambda(M_SUSY)/dM_3(M_X) = " << analytic_dParamsdM3(1) << endl;
      cout << "*    dm_Hd^2(M_SUSY)/dM_3(M_X) = " << analytic_dParamsdM3(2) << endl;
      cout << "*    dm_Hu^2(M_SUSY)/dM_3(M_X) = " << analytic_dParamsdM3(3) << endl;
      cout << "*    dm_s^2(M_SUSY)/dM_3(M_X) = " << analytic_dParamsdM3(4) << endl;
      cout << "*    dm_q3L^2(M_SUSY)/dM_3(M_X) = " << analytic_dParamsdM3(5) << endl;
      cout << "*    dm_tR^2(M_SUSY)/dM_3(M_X) = " << analytic_dParamsdM3(6) << endl;
      cout << "*    dA_t(M_SUSY)/dM_3(M_X) = " << analytic_dParamsdM3(7) << endl;
      cout << "**************************************************" << endl;
      cout << "* PERCENTAGE DIFFERENCES: " << endl;
      cout << "*    dlambda(M_SUSY)/dM_3(M_X) = " << error_dParamsdM3(0) << endl;
      cout << "*    dA_lambda(M_SUSY)/dM_3(M_X) = " << error_dParamsdM3(1) << endl;
      cout << "*    dm_Hd^2(M_SUSY)/dM_3(M_X) = " << error_dParamsdM3(2) << endl;
      cout << "*    dm_Hu^2(M_SUSY)/dM_3(M_X) = " << error_dParamsdM3(3) << endl;
      cout << "*    dm_s^2(M_SUSY)/dM_3(M_X) = " << error_dParamsdM3(4) << endl;
      cout << "*    dm_q3L^2(M_SUSY)/dM_3(M_X) = " << error_dParamsdM3(5) << endl;
      cout << "*    dm_tR^2(M_SUSY)/dM_3(M_X) = " << error_dParamsdM3(6) << endl;
      cout << "*    dA_t(M_SUSY)/dM_3(M_X) = " << error_dParamsdM3(7) << endl;
      cout << "**************************************************" << endl;
      cout << "* SUBTEST 11 RESULT: ";
      if (error_dParamsdM3(0) > tol || error_dParamsdM3(1) > tol || error_dParamsdM3(2) > tol ||
	  error_dParamsdM3(3) > tol || error_dParamsdM3(4) > tol || error_dParamsdM3(5) > tol ||
	  error_dParamsdM3(6) > tol || error_dParamsdM3(7) > tol)
	{
	  cout << "FAIL" << endl;
	  if (error_dParamsdM3(0) > tol) cout << "*    maximum error exceeded for derivative of lambda(M_SUSY)" << endl;
	  if (error_dParamsdM3(1) > tol) cout << "*    maximum error exceeded for derivative of A_lambda(M_SUSY)" << endl;
	  if (error_dParamsdM3(2) > tol) cout << "*    maximum error exceeded for derivative of m_Hd^2(M_SUSY)" << endl;
	  if (error_dParamsdM3(3) > tol) cout << "*    maximum error exceeded for derivative of m_Hu^2(M_SUSY)" << endl;
	  if (error_dParamsdM3(4) > tol) cout << "*    maximum error exceeded for derivative of m_s^2(M_SUSY)" << endl;
	  if (error_dParamsdM3(5) > tol) cout << "*    maximum error exceeded for derivative of m_q3L^2(M_SUSY)" << endl;
	  if (error_dParamsdM3(6) > tol) cout << "*    maximum error exceeded for derivative of m_tR^2(M_SUSY)" << endl;
	  if (error_dParamsdM3(7) > tol) cout << "*    maximum error exceeded for derivative of A_t(M_SUSY)" << endl;
	  hasPassedM3 = false;
	}
      else
	{
	  cout << "PASS" << endl; 
	  hasPassedM3 = true;
	}
      cout << "**************************************************" << endl;
      cout << "* SUBTEST 12: M_1' derivatives" << endl; 
      cout << "**************************************************" << endl;
      cout << "* NUMERICAL DERIVATIVES: " << endl;
      cout << "*    dlambda(M_SUSY)/dM_1'(M_X) = " << numeric_dParamsdM1p(0) << endl;
      cout << "*    dA_lambda(M_SUSY)/dM_1'(M_X) = " << numeric_dParamsdM1p(1) << endl;
      cout << "*    dm_Hd^2(M_SUSY)/dM_1'(M_X) = " << numeric_dParamsdM1p(2) << endl;
      cout << "*    dm_Hu^2(M_SUSY)/dM_1'(M_X) = " << numeric_dParamsdM1p(3) << endl;
      cout << "*    dm_s^2(M_SUSY)/dM_1'(M_X) = " << numeric_dParamsdM1p(4) << endl;
      cout << "*    dm_q3L^2(M_SUSY)/dM_1'(M_X) = " << numeric_dParamsdM1p(5) << endl;
      cout << "*    dm_tR^2(M_SUSY)/dM_1'(M_X) = " << numeric_dParamsdM1p(6) << endl;
      cout << "*    dA_t(M_SUSY)/dM_1'(M_X) = " << numeric_dParamsdM1p(7) << endl;
      cout << "**************************************************" << endl;
      cout << "* ANALYTIC DERIVATIVES: " << endl;
      cout << "*    dlambda(M_SUSY)/dM_1'(M_X) = " << analytic_dParamsdM1p(0) << endl;
      cout << "*    dA_lambda(M_SUSY)/dM_1'(M_X) = " << analytic_dParamsdM1p(1) << endl;
      cout << "*    dm_Hd^2(M_SUSY)/dM_1'(M_X) = " << analytic_dParamsdM1p(2) << endl;
      cout << "*    dm_Hu^2(M_SUSY)/dM_1'(M_X) = " << analytic_dParamsdM1p(3) << endl;
      cout << "*    dm_s^2(M_SUSY)/dM_1'(M_X) = " << analytic_dParamsdM1p(4) << endl;
      cout << "*    dm_q3L^2(M_SUSY)/dM_1'(M_X) = " << analytic_dParamsdM1p(5) << endl;
      cout << "*    dm_tR^2(M_SUSY)/dM_1'(M_X) = " << analytic_dParamsdM1p(6) << endl;
      cout << "*    dA_t(M_SUSY)/dM_1'(M_X) = " << analytic_dParamsdM1p(7) << endl;
      cout << "**************************************************" << endl;
      cout << "* PERCENTAGE DIFFERENCES: " << endl;
      cout << "*    dlambda(M_SUSY)/dM_1'(M_X) = " << error_dParamsdM1p(0) << endl;
      cout << "*    dA_lambda(M_SUSY)/dM_1'(M_X) = " << error_dParamsdM1p(1) << endl;
      cout << "*    dm_Hd^2(M_SUSY)/dM_1'(M_X) = " << error_dParamsdM1p(2) << endl;
      cout << "*    dm_Hu^2(M_SUSY)/dM_1'(M_X) = " << error_dParamsdM1p(3) << endl;
      cout << "*    dm_s^2(M_SUSY)/dM_1'(M_X) = " << error_dParamsdM1p(4) << endl;
      cout << "*    dm_q3L^2(M_SUSY)/dM_1'(M_X) = " << error_dParamsdM1p(5) << endl;
      cout << "*    dm_tR^2(M_SUSY)/dM_1'(M_X) = " << error_dParamsdM1p(6) << endl;
      cout << "*    dA_t(M_SUSY)/dM_1'(M_X) = " << error_dParamsdM1p(7) << endl;
      cout << "**************************************************" << endl;
      cout << "* SUBTEST 12 RESULT: ";
      if (error_dParamsdM1p(0) > tol || error_dParamsdM1p(1) > tol || error_dParamsdM1p(2) > tol ||
	  error_dParamsdM1p(3) > tol || error_dParamsdM1p(4) > tol || error_dParamsdM1p(5) > tol ||
	  error_dParamsdM1p(6) > tol || error_dParamsdM1p(7) > tol)
	{
	  cout << "FAIL" << endl;
	  if (error_dParamsdM1p(0) > tol) cout << "*    maximum error exceeded for derivative of lambda(M_SUSY)" << endl;
	  if (error_dParamsdM1p(1) > tol) cout << "*    maximum error exceeded for derivative of A_lambda(M_SUSY)" << endl;
	  if (error_dParamsdM1p(2) > tol) cout << "*    maximum error exceeded for derivative of m_Hd^2(M_SUSY)" << endl;
	  if (error_dParamsdM1p(3) > tol) cout << "*    maximum error exceeded for derivative of m_Hu^2(M_SUSY)" << endl;
	  if (error_dParamsdM1p(4) > tol) cout << "*    maximum error exceeded for derivative of m_s^2(M_SUSY)" << endl;
	  if (error_dParamsdM1p(5) > tol) cout << "*    maximum error exceeded for derivative of m_q3L^2(M_SUSY)" << endl;
	  if (error_dParamsdM1p(6) > tol) cout << "*    maximum error exceeded for derivative of m_tR^2(M_SUSY)" << endl;
	  if (error_dParamsdM1p(7) > tol) cout << "*    maximum error exceeded for derivative of A_t(M_SUSY)" << endl;
	  hasPassedM1p = false;
	}
      else
	{
	  cout << "PASS" << endl; 
	  hasPassedM1p = true;
	}
      cout << "**************************************************" << endl;
      cout << "* TEST RESULT: ";
      if (!hasPassedMh1Sq || !hasPassedMh2Sq || !hasPassedMsSq || !hasPassedMqSq || !hasPassedMuSq
	  || !hasPassedLambda || !hasPassedAlambda || !hasPassedAt || !hasPassedM1 || !hasPassedM2
	  || !hasPassedM3 || !hasPassedM1p)
	{
	  cout << "FAIL" << endl;
	  if (!hasPassedMh1Sq) cout << "*    failed subtest 1: m_Hd^2 derivatives" << endl;
	  if (!hasPassedMh2Sq) cout << "*    failed subtest 2: m_Hu^2 derivatives" << endl;
	  if (!hasPassedMsSq) cout << "*    failed subtest 3: m_s^2 derivatives" << endl;
	  if (!hasPassedMqSq) cout << "*    failed subtest 4: m_q3L^2 derivatives" << endl;
	  if (!hasPassedMuSq) cout << "*    failed subtest 5: m_tR^2 derivatives" << endl;
	  if (!hasPassedLambda) cout << "*    failed subtest 6: lambda derivatives" << endl;
	  if (!hasPassedAlambda) cout << "*    failed subtest 7: A_lambda derivatives" << endl;
	  if (!hasPassedAt) cout << "*    failed subtest 8: A_t derivatives" << endl;
	  if (!hasPassedM1) cout << "*    failed subtest 9: M_1 derivatives" << endl;
	  if (!hasPassedM2) cout << "*    failed subtest 10: M_2 derivatives" << endl;
	  if (!hasPassedM3) cout << "*    failed subtest 11: M_3 derivatives" << endl;
	  if (!hasPassedM1p) cout << "*    failed subtest 12: M_1' derivatives" << endl;
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

// Calculates percentage difference between corresponding elements in the matrices a and b
Eigen::MatrixXd calcPercentageDifferences(const Eigen::MatrixXd & a, const Eigen::MatrixXd & b)
{
  // Matrix dimensions should agree
  if (a.rows() != b.rows() || a.cols() != b.cols())
    {
      ostringstream ii;
      ii << "ERROR: matrix dimensions do not agree: in calcPercentageDifferences" << endl;
      ii << "     : arg1 dimensions (" << a.rows() << " x " << a.cols() << ") do not match ";
      ii << " arg2 dimensions (" << b.rows() << " x " << b.cols() << ")" << endl;
      throw ii.str();
    } 

  // Convert to arrays for coefficient-wise division
  Eigen::ArrayXd array_a = a.array();
  Eigen::ArrayXd array_b = b.array();

  // Get percentage differences
  Eigen::ArrayXd numerator = array_a-array_b;
  Eigen::ArrayXd denominator = 0.5*(array_a+array_b);
  Eigen::ArrayXd errors = 100.0 * ((array_a-array_b)/(0.5*(array_a+array_b))).abs();

  // Check for elements where denominator vanishes: if the element is zero
  // in both a and b, return 0, else return large number
  for (int i = 0; i < errors.rows(); i++)
    {
      for (int j = 0; j < errors.cols(); j++)
	{
	  if (Abs(denominator(i,j)) < EPSTOL)
	    {
	      if (Abs(array_a(i,j)) < EPSTOL && Abs(array_b(i,j)) < EPSTOL)
		{
		  errors(i,j) = 0.0;
		}
	      else
		{
		  errors(i,j) = numberOfTheBeast;
		}
	    }
	}
    }

  return errors.matrix();
} 

genericE6SSM_soft_parameters* tempmodel;
double scale;
int parChoice;
Eigen::ArrayXd* pars;
int nLogs;
int nLps;

double getApproximateMh1Squared(double param)
{
  genericE6SSM_soft_parameters r(*tempmodel);
  Eigen::ArrayXd input_pars = *pars;

  input_pars(parChoice) = param;

  pE6SSMftBCs(r, input_pars);

  double logCoeff = doCalcMh1SquaredLogCoeff(r, nLps);
  double logSqCoeff = doCalcMh1SquaredLogSqCoeff(r, nLogs);

  double t = Log(scale/r.get_scale());

  double mH1Sq = r.get_mHd2() + t * logCoeff + Sqr(t) * logSqCoeff;

  return mH1Sq;

}

double getApproximateMh2Squared(double param)
{
  genericE6SSM_soft_parameters r(*tempmodel);
  Eigen::ArrayXd input_pars = *pars;

  input_pars(parChoice) = param;

  pE6SSMftBCs(r, input_pars);

  double logCoeff = doCalcMh2SquaredLogCoeff(r, nLps);
  double logSqCoeff = doCalcMh2SquaredLogSqCoeff(r, nLogs);

  double t = Log(scale/r.get_scale());

  double mH2Sq = r.get_mHu2() + t * logCoeff + Sqr(t) * logSqCoeff;

  return mH2Sq;
}

double getApproximateMsSquared(double param)
{
  genericE6SSM_soft_parameters r(*tempmodel);
  Eigen::ArrayXd input_pars = *pars;

  input_pars(parChoice) = param;

  pE6SSMftBCs(r, input_pars);

  double logCoeff = doCalcMsSquaredLogCoeff(r, nLps);
  double logSqCoeff = doCalcMsSquaredLogSqCoeff(r, nLogs);

  double t = Log(scale/r.get_scale());

  double mSSq = r.get_ms2() + t * logCoeff + Sqr(t) * logSqCoeff;

  return mSSq;
}

double getApproximateMqL3Squared(double param)
{
  genericE6SSM_soft_parameters r(*tempmodel);
  Eigen::ArrayXd input_pars = *pars;

  input_pars(parChoice) = param;

  pE6SSMftBCs(r, input_pars);

  double logCoeff = doCalcmqL3SquaredLogCoeff(r, nLps);
  double logSqCoeff = doCalcmqL3SquaredLogSqCoeff(r, nLogs);

  double t = Log(scale/r.get_scale());

  double mqL3Sq = r.get_mq2(2,2) + t * logCoeff + Sqr(t) * logSqCoeff;

  return mqL3Sq;
}

double getApproximateMtRSquared(double param)
{
  genericE6SSM_soft_parameters r(*tempmodel);
  Eigen::ArrayXd input_pars = *pars;

  input_pars(parChoice) = param;

  pE6SSMftBCs(r, input_pars);

  double logCoeff = doCalcmtRSquaredLogCoeff(r, nLps);
  double logSqCoeff = doCalcmtRSquaredLogSqCoeff(r, nLogs);

  double t = Log(scale/r.get_scale());

  double mtRSq = r.get_mu2(2,2) + t * logCoeff + Sqr(t) * logSqCoeff;

  return mtRSq;
}

double getApproximateLambda(double param)
{
  genericE6SSM_soft_parameters r(*tempmodel);
  Eigen::ArrayXd input_pars = *pars;

  input_pars(parChoice) = param;

  pE6SSMftBCs(r, input_pars);

  double logCoeff = doCalcLambda3LogCoeff(r, nLps);
  double logSqCoeff = doCalcLambda3LogSqCoeff(r, nLogs);

  double t = Log(scale/r.get_scale());

  double lambda = r.get_Lambdax() + t * logCoeff + Sqr(t) * logSqCoeff;

  return lambda;
}

double getApproximateAlambda(double param)
{
  genericE6SSM_soft_parameters r(*tempmodel);
  Eigen::ArrayXd input_pars = *pars;

  input_pars(parChoice) = param;

  pE6SSMftBCs(r, input_pars);

  double logCoeff = doCalcAlambda3LogCoeff(r, nLps);
  double logSqCoeff = doCalcAlambda3LogSqCoeff(r, nLogs);

  double t = Log(scale/r.get_scale());

  double Tlambda = r.get_TLambdax();
  double lambda = r.get_Lambdax();
  double Alambda_mx;


  if (Abs(Tlambda) < EPSTOL) 
    {
      Alambda_mx = 0.0;
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
      Alambda_mx = Tlambda/lambda;
    }

  double Alambda = Alambda_mx + t * logCoeff + Sqr(t) * logSqCoeff;

  return Alambda;
}

double getApproximateAt(double param)
{
  genericE6SSM_soft_parameters r(*tempmodel);
  Eigen::ArrayXd input_pars = *pars;

  input_pars(parChoice) = param;

  pE6SSMftBCs(r, input_pars);

  double logCoeff = doCalcAtLogCoeff(r, nLps);
  double logSqCoeff = doCalcAtLogSqCoeff(r, nLogs);

  double t = Log(scale/r.get_scale());

  double TYt = r.get_TYu(2,2);
  double yt = r.get_Yu(2,2);
  double At_mx;


  if (Abs(TYt) < EPSTOL) 
    {
      At_mx = 0.0;
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
      At_mx = TYt/yt;
    }

  double At = At_mx + t * logCoeff + Sqr(t) * logSqCoeff;

  return At;
}

Eigen::Matrix<double,8,1> doCalcNumericDerivs(genericE6SSM_soft_parameters r, double q, unsigned i, int lps, int logs)
{
  tempmodel = &r;
  parChoice = i;
  nLps = lps;
  nLogs = logs;

  Eigen::ArrayXd params;
  params.resize(tuning_parameters::NUMESSMTUNINGPARS,1);

  params(tuning_parameters::lam3) = r.get_Lambdax();
  params(tuning_parameters::mH13Sq) = r.get_mHd2();
  params(tuning_parameters::mH23Sq) = r.get_mHu2();
  params(tuning_parameters::mS3Sq) = r.get_ms2();
  params(tuning_parameters::mqL3Sq) = r.get_mq2(2,2);
  params(tuning_parameters::mtRSq) = r.get_mu2(2,2);
  params(tuning_parameters::M1) = r.get_MassB();
  params(tuning_parameters::M2) = r.get_MassWB();
  params(tuning_parameters::M3) = r.get_MassG();
  params(tuning_parameters::M1p) = r.get_MassBp();

  double TYt = r.get_TYu(2,2);
  double yt = r.get_Yu(2,2);
  double At;

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

  params(tuning_parameters::Au3) = At;

  double Tlambda = r.get_TLambdax();
  double lambda = r.get_Lambdax();
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

  params(tuning_parameters::Alam3) = Alambda;

  pars = &params;

  scale = q;

  double deriv, h, temp, err;

  double epsilon = 1.0e-5;

  // Initial estimate for step size h.
  h = epsilon*Abs(params(i));
  if (h == 0.0) h = epsilon;
  temp = params(i);
  params(i) = temp + h;
  h = params(i) - temp;

  Eigen::Matrix<double,8,1> derivs, errvec;

  derivs(0) = calcDerivative(getApproximateLambda, temp, h, &err); errvec(0) = err;
  derivs(1) = calcDerivative(getApproximateAlambda, temp, h, &err); errvec(1) = err;
  derivs(2) = calcDerivative(getApproximateMh1Squared, temp, h, &err); errvec(2) = err;
  derivs(3) = calcDerivative(getApproximateMh2Squared, temp, h, &err); errvec(3) = err;
  derivs(4) = calcDerivative(getApproximateMsSquared, temp, h, &err); errvec(4) = err;
  derivs(5) = calcDerivative(getApproximateMqL3Squared, temp, h, &err); errvec(5) = err;
  derivs(6) = calcDerivative(getApproximateMtRSquared, temp, h, &err); errvec(6) = err;
  derivs(7) = calcDerivative(getApproximateAt, temp, h, &err); errvec(7) = err;

  for (int j = 0; j < 8; j++)
    {
      if (Abs(errvec(j)) > 1.0e-8 && Abs(errvec(j)/derivs(j)) > 1.0)
	{
	  derivs(j) = numberOfTheBeast;
	}
    }

  return derivs;
}
