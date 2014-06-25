/*
  Does a small scan over the E6SSM parameter space
  and calculates the fine tuning numerically and approximately,
  printing the results for comparison. Also evaluates the time per point.
 */

#include "spectrum/models/genericE6SSM/genericE6SSM_input_parameters.hpp"
#include "spectrum/models/genericE6SSM/genericE6SSM_spectrum_generator.hpp"
#include "spectrum/models/genericE6SSM/genericE6SSM_info.hpp"
#include "spectrum/models/genericE6SSM/genericE6SSM_two_scale_model.hpp"
#include "spectrum/src/error.hpp"
#include "spectrum/src/spectrum_generator_settings.hpp"
#include "spectrum/src/lowe.h"
#include "spectrum/src/command_line_options.hpp"
#include "spectrum/src/wrappers.hpp"
#include "spectrum/src/problems.hpp"

#include <iostream>
#include <cstdlib>
#include <sys/time.h>


#include "src/flags.h"
#include "src/tuningnumerics.h"
#include "src/tuningutils.h"
#include "src/essmtuningutils.h"

using namespace flexiblesusy;
using namespace essm_tuning_utils;

double get_softAu(flexiblesusy::genericE6SSM<flexiblesusy::Two_scale> const & m, int i, int j);
double get_softAd(flexiblesusy::genericE6SSM<flexiblesusy::Two_scale> const & m, int i, int j);
double get_softAe(flexiblesusy::genericE6SSM<flexiblesusy::Two_scale> const & m, int i, int j);

flexiblesusy::genericE6SSM<flexiblesusy::Two_scale> doSimplifiedSpectrum(Eigen::Matrix<double,3,3> yuin, 
									 Eigen::Matrix<double,3,3> ydin, 
									 Eigen::Matrix<double,3,3> yein, 
									 double g1in, double g2in, double g3in, double gNin,
									 double lambda3in, Eigen::Matrix<double,2,2> lambda12in,
									 Eigen::Matrix<double,3,3> kappain, 
									 double mupr, double tanb, double hvev, double svev,
									 double MBin, double MWBin, double MGin, double MBpin,
									 Eigen::Matrix<double,3,3> tuin, 
									 Eigen::Matrix<double,3,3> tdin, 
									 Eigen::Matrix<double,3,3> tein, 
									 double tlambda3in, Eigen::Matrix<double,2,2> tlambdain, 
									 Eigen::Matrix<double,3,3> tkappain,
									 Eigen::Matrix<double,3,3> mQlSq, 
									 Eigen::Matrix<double,3,3> mUrSq, 
									 Eigen::Matrix<double,3,3> mDrSq, 
									 Eigen::Matrix<double,3,3> mLlSq, 
									 Eigen::Matrix<double,3,3> mErSq, 
									 Eigen::Matrix<double,3,3> mDxSq,
									 Eigen::Matrix<double,3,3> mDxbarSq, 
									 Eigen::Matrix<double,2,2> mHdISq, 
									 Eigen::Matrix<double,2,2> mHuISq, 
									 Eigen::Matrix<double,2,2> mSISq,
									 double Bmupr,  double mHpSq, double mHpbarSq, double mx, 
									 int l, int t, 
									 double tol, double & mHdSq, double & mHuSq, double & mSSq, 
									 double & ms, 
									 bool & hasEWSBProblem, bool & squarksTachyons, 
									 bool & vectorBosonsTachyons, bool & higgsTachyons, 
									 bool & tadpoleProblem, bool & poleHiggsTachyons, 
									 bool & inaccurateHiggsMass, bool & hasSeriousProblem);


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

  outputCharacteristics(8);

  /*
    Set up model
   */

  // Values for gauge and Yukawa couplings (not scanned over)
  Eigen::Matrix<double,3,3> yu, yd, ye, kappa;
  Eigen::Matrix<double,2,2> lambda12;
  double g1, g2, g3, gN;

  g1 = 0.45;
  gN = 0.46;
  g2 = 0.6;
  g3 = 0.98;

  yu(2,2) = 0.8;
  yd(2,2) = 0.1;
  ye(2,2) = 0.05;

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

  // Assume universal and diagonal kappa
  kappa(0,0) = 0.6; kappa(1,1) = 0.6; kappa(2,2) = 0.6;

  kappa(0,1) = 0.0; kappa(0,2) = 0.0;
  kappa(1,0) = 0.0; kappa(1,2) = 0.0;
  kappa(2,0) = 0.0; kappa(2,1) = 0.0;

  // Likewise for the inert lambda couplings
  lambda12(0,0) = 0.2; lambda12(0,1) = 0.0;
  lambda12(1,0) = 0.0; lambda12(1,1) = 0.2;

  // Corresponding soft trilinears: note we set A_b = A_tau = 0
  Eigen::Matrix<double,3,3> tyu, tyd, tye, tkappa;
  Eigen::Matrix<double,2,2> tlambda12;

  tyu(2,2) = yu(2,2) * 5000.0; // GeV
  tyd(2,2) = 0.0; // GeV
  tye(2,2) = 0.0; // GeV

  tyu(0,1) = 0.0; tyu(0,2) = 0.0; // GeV
  tyu(1,0) = 0.0; tyu(1,2) = 0.0; // GeV
  tyu(2,0) = 0.0; tyu(2,1) = 0.0; // GeV

  tyd(0,1) = 0.0; tyd(0,2) = 0.0; // GeV
  tyd(1,0) = 0.0; tyd(1,2) = 0.0; // GeV
  tyd(2,0) = 0.0; tyd(2,1) = 0.0; // GeV

  tye(0,1) = 0.0; tye(0,2) = 0.0; // GeV
  tye(1,0) = 0.0; tye(1,2) = 0.0; // GeV
  tye(2,0) = 0.0; tye(2,1) = 0.0; // GeV

  // Kappa and T_kappa are assumed diagonal, and we set A_kappa = 0
  tkappa(0,0) = 0.0; tkappa(1,1) = 0.0; tkappa(2,2) = 0.0; // GeV

  tkappa(0,1) = 0.0; tkappa(0,2) = 0.0; // GeV
  tkappa(1,0) = 0.0; tkappa(1,2) = 0.0; // GeV
  tkappa(2,0) = 0.0; tkappa(2,1) = 0.0; // GeV

  tlambda12(0,0) = 0.0; tlambda12(0,1) = 0.0; // GeV
  tlambda12(1,0) = 0.0; tlambda12(1,1) = 0.0; // GeV

  // Survival Higgs bilinears
  double mupr, Bmupr, mHpSq, mHpbarSq;

  mHpSq = Sqr(5000.0); // GeV^2
  mHpbarSq = Sqr(5000.0); // GeV^2

  mupr = 5000.0; // GeV
  Bmupr = 5000.0; // GeV^2

  // Soft Higgs and singlet masses (note 3rd
  // generation are output of EWSB conditions)
  Eigen::Matrix<double,2,2> mH1ISq, mH2ISq, mSISq;
  double mHdSq, mHuSq, mSSq;

  mH1ISq(0,0) = Sqr(5000.0); mH1ISq(0,1) = 0.0; // GeV^2
  mH1ISq(1,0) = 0.0; mH1ISq(1,1) = Sqr(5000.0); // GeV^2

  mH2ISq(0,0) = Sqr(5000.0); mH2ISq(0,1) = 0.0; // GeV^2
  mH2ISq(1,0) = 0.0; mH2ISq(1,1) = Sqr(5000.0); // GeV^2

  mSISq(0,0) = Sqr(5000.0); mSISq(0,1) = 0.0; // GeV^2
  mSISq(1,0) = 0.0; mSISq(1,1) = Sqr(5000.0); // GeV^2

  // Gaugino soft masses
  double M1, M2, M3, M1p;

  M1 = 300.0; // GeV
  M1p = 300.0; // GeV
  M2 = 500.0; // GeV
  M3 = 2000.0; // GeV

  mHdSq = 1.0e8; // GeV^2
  mHuSq = -1.5e8; // GeV^2
  mSSq = -1.5e8; // GeV^2

  // Soft scalar masses, all fixed to 5 TeV if not scanned over
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

  // VEVs, with s chosen so M_Z' > 2.5 TeV
  double v, s, tb, v1, v2;

  tb = 10.0;
  s = 6700.0; // GeV
  v = 246.0; // GeV
  v1 = v/Sqrt(1.0+tb*tb); // GeV
  v2 = v1*tb; // GeV

  // Input scale
  double MX = 20000.0; // GeV
  double MS = 2000.0; // GeV

  // Parameters for RGE running
  int LOOPS = 2;
  int LEADINGLOGS = 1;
  int THRESH = 0;

  // E6 mixing angle
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

  // Scanned variables
  double lambda, tlambda, Alambda, mqL3Sq, mtRSq, At;

  // Lower bounds
  double tb_low = 10.0;
  double lambda_low = 0.8;//-3.0;
  double Alambda_low = 1000.0; // GeV
  double mqL3Sq_low = Sqr(500.0); // GeV^2
  double mtRSq_low = Sqr(500.0); // GeV^2
  double At_low = 1000.0; // GeV
  double M2_low = 200.0;//200.0; // GeV
  double M3_low = 0.0;//200.0; // GeV

  // Upper bounds
  double tb_up = 50.0;
  double lambda_up = 3.0;
  double Alambda_up = 10000.0; // GeV
  double mqL3Sq_up = Sqr(2000.0); // GeV^2
  double mtRSq_up = Sqr(2000.0); // GeV^2
  double At_up = 10000.0; // GeV
  double M2_up = 2000.0; // GeV
  double M3_up = 5000.0; // GeV

  // Number of points (just a linear scan for this simple test)
  int tb_npts = 1;
  int lambda_npts = 1;
  int Alambda_npts = 1;
  int mqL3Sq_npts = 1;
  int mtRSq_npts = 1;
  int At_npts = 1;
  int M2_npts = 1;
  int M3_npts = 100;

  // Increments
  double tb_incr = 0.0;
  double lambda_incr = 0.0;
  double Alambda_incr = 0.0;
  double mqL3Sq_incr = 0.0;
  double mtRSq_incr = 0.0;
  double At_incr = 0.0;
  double M2_incr = 0.0;
  double M3_incr = 0.0;

  // Update increments if necessary
  if (tb_npts > 1)
    {
      tb_incr = (tb_up-tb_low)/(((double)tb_npts)-1.0);
    }
  if (lambda_npts > 1)
    {
      lambda_incr = (lambda_up-lambda_low)/(((double)lambda_npts)-1.0);
    }
  if (Alambda_npts > 1)
    {
      Alambda_incr = (Alambda_up-Alambda_low)/(((double)Alambda_npts)-1.0);
    }
  if (mqL3Sq_npts > 1)
    {
      mqL3Sq_incr = (mqL3Sq_up-mqL3Sq_low)/(((double)mqL3Sq_npts)-1.0);
    }
  if (mtRSq_npts > 1)
    {
      mtRSq_incr = (mtRSq_up-mtRSq_low)/(((double)mtRSq_npts)-1.0);
    }
  if (At_npts > 1)
    {
      At_incr = (At_up-At_low)/(((double)At_npts)-1.0);
    }
  if (M2_npts > 1)
    {
      M2_incr = (M2_up-M2_low)/(((double)M2_npts)-1.0);
    }
  if (M3_npts > 1)
    {
      M3_incr = (M3_up-M3_low)/(((double)M3_npts)-1.0);
    }

  // Variables for storing CPU and wall time
  struct timeval tvWall, tv2Wall, tvNumeric, tv2Numeric, tvApprox, tv2Approx;
  double cpuStart, cpuEnd, cpuStartNumeric, cpuEndNumeric, cpuStartApprox, cpuEndApprox;  
  int wall_time, wall_time_numeric, wall_time_approx;

  // Vectors for storing fine tunings
  Eigen::Matrix<double,tuning_parameters::NUMESSMTUNINGPARS,1> tuningNumeric, tuningApprox;

  // Error flags
  bool hasWrongVacuum = false; //< problem with getting the physical vacuum
  bool hasEWSBProblem = false; //< problem implementing EWSB conditions to requested tolerance
  bool hasSeriousProblem = false; //< serious numerical problem in spectrum calculation
  bool squarksTachyons = false; //< tachyonic DR bar stops or sbottoms
  bool vectorBosonsTachyons = false; //< tachyonic DR bar Z, W or Z'
  bool higgsTachyons = false; //< DR bar m_h0 or m_A0 tachyon
  bool poleHiggsTachyons = false; //< physical m_h0 or m_A0 tachyon
  bool inaccurateHiggsMass = false; //< inaccurate Higgs mass calculation
  bool tadpoleProblem = false; //< problem with calculating tadpoles

  bool hasCCBProblem = false; //< problem with CCB minima
  bool hasUFBProblem = false; //< UFB problem
  bool hasTuningError = false; //< numerical error in tuning calculation

  // Vector for storing fine tuning parameters
  Eigen::ArrayXd pars;
  pars.resize(tuning_parameters::NUMESSMTUNINGPARS);
  
  // Models that we will calculate the fine tuning for
  genericE6SSM<Two_scale> model;

  // Use approximate analytic expressions for tuning?
  // Turned off at the moment because I want to see the impact
  // of calculating the full 2-loop derivatives numerically
  bool useApproxSolns = false;

  // Suppress print out
  ENABLE_DEBUG = false;

  // Do scan and write output
  for (int i = 1; i <= tb_npts; i++)
  {
  for (int j = 1; j <= lambda_npts; j++)
  {
  for (int k = 1; k <= Alambda_npts; k++)
  {
  for (int l = 1; l <= mqL3Sq_npts; l++)
  {
  for (int m = 1; m <= mtRSq_npts; m++)
  {
  for (int n = 1; n <= At_npts; n++)
  {
  for (int o = 1; o <= M2_npts; o++)
  {
  for (int p = 1; p <= M3_npts; p++)
  {

    // Update variable values
    tb = tb_low+(((double)i)-1.0)*tb_incr;
    lambda = lambda_low+(((double)j)-1.0)*lambda_incr;
    Alambda = Alambda_low+(((double)k)-1.0)*Alambda_incr;
    mqL3Sq = mqL3Sq_low+(((double)l)-1.0)*mqL3Sq_incr;
    mtRSq = mtRSq_low+(((double)m)-1.0)*mtRSq_incr;
    At = At_low+(((double)n)-1.0)*At_incr;
    M2 = M2_low+(((double)o)-1.0)*M2_incr;
    M3 = M3_low+(((double)p)-1.0)*M3_incr;

    mQlSq(2,2) = mqL3Sq;
    mUrSq(2,2) = mtRSq;
    tlambda = lambda*Alambda;
    tyu(2,2) = yu(2,2) * At;

    // Start timing whole point
    gettimeofday(&tvWall, NULL);
    cpuStart = ((double) clock())/(CLOCKS_PER_SEC);
    
    hasEWSBProblem = false;
    squarksTachyons = false;
    vectorBosonsTachyons = false;
    higgsTachyons = false;
    tadpoleProblem = false;
    poleHiggsTachyons = false;
    inaccurateHiggsMass = false;
    hasSeriousProblem = false;

    bool hasTuningProblem = false;    

    // Model is returned at M_{SUSY}
    model = doSimplifiedSpectrum(yu, yd, ye, g1, g2, g3, gN, lambda, lambda12,
				   kappa, mupr, tb, v, s,
				   M1, M2, M3, M1p, tyu, tyd, tye, 
				   tlambda, tlambda12, tkappa, mQlSq, mUrSq, mDrSq, 
				   mLlSq, mErSq, mDxSq, mDxbarSq, mH1ISq, mH2ISq, mSISq,
				   Bmupr, mHpSq, mHpbarSq, MX, LOOPS, THRESH, 
				   TOLEWSB, mHdSq, mHuSq, mSSq, MS, 
				   hasEWSBProblem, squarksTachyons, vectorBosonsTachyons, higgsTachyons, 
				   tadpoleProblem, poleHiggsTachyons, inaccurateHiggsMass, hasSeriousProblem);

    model.run_to(MX);

    // Calculate the fine tuning approximately, and time it
    pars(tuning_parameters::lam3) = lambda;
    pars(tuning_parameters::Alam3) = Alambda;
    pars(tuning_parameters::M1) = M1;
    pars(tuning_parameters::M2) = M2;
    pars(tuning_parameters::M3) = M3;
    pars(tuning_parameters::M1p) = M1p;
    pars(tuning_parameters::Au3) = At;
    pars(tuning_parameters::mH13Sq) = model.get_mHd2();
    pars(tuning_parameters::mH23Sq) = model.get_mHu2();
    pars(tuning_parameters::mS3Sq) = model.get_ms2();
    pars(tuning_parameters::mqL3Sq) = mQlSq(2,2);
    pars(tuning_parameters::mtRSq) = mUrSq(2,2);

    gettimeofday(&tvApprox, NULL);
    cpuStartApprox = ((double) clock())/(CLOCKS_PER_SEC);
    tuningApprox = doCalcESSMTuningApprox(model, MS, MX, pars, hasTuningProblem, useApproxSolns);
    cpuEndApprox = ((double) clock())/(CLOCKS_PER_SEC);
    gettimeofday(&tv2Approx, NULL);
    wall_time_approx = (tv2Approx.tv_sec-tvApprox.tv_sec)*1000000+((int)tv2Approx.tv_usec-(int)tvApprox.tv_usec);

    // Finish timing whole point
    cpuEnd = ((double) clock())/(CLOCKS_PER_SEC);
    gettimeofday(&tv2Wall, NULL);
    wall_time = (tv2Wall.tv_sec-tvWall.tv_sec)*1000000+((int)tv2Wall.tv_usec-(int)tvWall.tv_usec);

    // Calculate the fine tuning numerically, and time it (note that this is not included
    // in the total time per point, as the tuning is not calculated numerically in the
    // actual scan)
    gettimeofday(&tvNumeric, NULL);
    cpuStartNumeric = ((double) clock())/(CLOCKS_PER_SEC);
    tuningNumeric = doCalcESSMTuningNumerically(model, MS, MX, pars, pE6SSMftBCs);
    cpuEndNumeric = ((double) clock())/(CLOCKS_PER_SEC);
    gettimeofday(&tv2Numeric, NULL);
    wall_time_numeric = (tv2Numeric.tv_sec-tvNumeric.tv_sec)*1000000+((int)tv2Numeric.tv_usec-(int)tvNumeric.tv_usec);

    // Calculate tunings semi-numerically
    Eigen::Matrix<double,tuning_parameters::NUMESSMTUNINGPARS,1> tuningSemianalytic;
    tuningSemianalytic = doCalcESSMTuningSemianalytic(model, MS, MX, hasTuningProblem);

    if (hasSeriousProblem)
      {
	cerr << "WARNING: serious problem encountered at point: declining to write data" << endl;
      }
    else
      {
	// Write results
	cout << wall_time << " ";
	cout << (int) ((cpuEnd-cpuStart)*1000000) << " ";
	cout << (int) ((cpuEndNumeric-cpuStartNumeric)*1000000) << " ";
	cout << M3 << " ";
	cout << tuningApprox.maxCoeff() << " ";
	cout << tuningNumeric.maxCoeff() << " ";
	cout << tuningSemianalytic.maxCoeff() << " ";
	cout << tuningApprox(tuning_parameters::lam3) << " ";
	cout << tuningNumeric(tuning_parameters::lam3) << " ";
	cout << tuningApprox(tuning_parameters::Alam3) << " ";
	cout << tuningNumeric(tuning_parameters::Alam3) << " ";
	cout << tuningApprox(tuning_parameters::mH13Sq) << " ";
	cout << tuningNumeric(tuning_parameters::mH13Sq) << " ";
	cout << tuningApprox(tuning_parameters::mH23Sq) << " ";
	cout << tuningNumeric(tuning_parameters::mH23Sq) << " ";
	cout << tuningApprox(tuning_parameters::mS3Sq) << " ";
	cout << tuningNumeric(tuning_parameters::mS3Sq) << " ";
	cout << tuningApprox(tuning_parameters::mqL3Sq) << " ";
	cout << tuningNumeric(tuning_parameters::mqL3Sq) << " ";
	cout << tuningApprox(tuning_parameters::mtRSq) << " ";
	cout << tuningNumeric(tuning_parameters::mtRSq) << " ";
	cout << tuningApprox(tuning_parameters::Au3) << " ";
	cout << tuningNumeric(tuning_parameters::Au3) << " ";
	cout << tuningApprox(tuning_parameters::M1) << " ";
	cout << tuningNumeric(tuning_parameters::M1) << " ";
	cout << tuningApprox(tuning_parameters::M2) << " ";
	cout << tuningNumeric(tuning_parameters::M2) << " ";
	cout << tuningApprox(tuning_parameters::M3) << " ";
	cout << tuningNumeric(tuning_parameters::M3) << " ";
	cout << tuningSemianalytic(tuning_parameters::M3) << " ";
	cout << tuningApprox(tuning_parameters::M1p) << " ";
	cout << tuningNumeric(tuning_parameters::M1p) << " ";
	if (hasTuningProblem || hasEWSBProblem)
	  {
	    cout << "1" << endl;
	  }
	else
	  {
	    cout << "0" << endl;
	  }
      }

  }//<M3
  }//<M2
  }//<At
  }//<mtRSq
  }//<mqL3Sq
  }//<Alambda
  }//<lambda
  }//<tb

  return 0;

}


flexiblesusy::genericE6SSM<flexiblesusy::Two_scale> doSimplifiedSpectrum(Eigen::Matrix<double,3,3> yuin, 
									 Eigen::Matrix<double,3,3> ydin, 
									 Eigen::Matrix<double,3,3> yein, 
									 double g1in, double g2in, double g3in, double gNin,
									 double lambda3in, Eigen::Matrix<double,2,2> lambda12in,
									 Eigen::Matrix<double,3,3> kappain, 
									 double mupr, double tanb, double hvev, double svev,
									 double MBin, double MWBin, double MGin, double MBpin,
									 Eigen::Matrix<double,3,3> tuin, 
									 Eigen::Matrix<double,3,3> tdin, 
									 Eigen::Matrix<double,3,3> tein, 
									 double tlambda3in, Eigen::Matrix<double,2,2> tlambda12in, 
									 Eigen::Matrix<double,3,3> tkappain,
									 Eigen::Matrix<double,3,3> mQlSq, 
									 Eigen::Matrix<double,3,3> mUrSq, 
									 Eigen::Matrix<double,3,3> mDrSq, 
									 Eigen::Matrix<double,3,3> mLlSq, 
									 Eigen::Matrix<double,3,3> mErSq, 
									 Eigen::Matrix<double,3,3> mDxSq,
									 Eigen::Matrix<double,3,3> mDxbarSq, 
									 Eigen::Matrix<double,2,2> mHdISq, 
									 Eigen::Matrix<double,2,2> mHuISq, 
									 Eigen::Matrix<double,2,2>  mSISq,
									 double Bmupr,  double mHpSq, double mHpbarSq, double mx, 
									 int l, int t, 
									 double tol, double & mHdSq, double & mHuSq, double & mSSq, 
									 double & ms, 
									 bool & hasEWSBProblem, bool & squarksTachyons, 
									 bool & vectorBosonsTachyons, bool & higgsTachyons, 
									 bool & tadpoleProblem, bool & poleHiggsTachyons, 
									 bool & inaccurateHiggsMass, bool & hasSeriousProblem)
{
  using namespace flexiblesusy;
  using namespace essm_tuning_utils;

  // While the expressions we have for the soft masses have hardcoded charges, we will
  // fix the E6 mixing angle to give us the U(1)_N E6 model.
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

  // In the pE6SSM, first and second generation Yukawas are assumed to vanish,
  // and there is no mixing.
  yuin(0,0) = 0.0;
  yuin(1,1) = 0.0;
  ydin(0,0) = 0.0;
  ydin(1,1) = 0.0;
  yein(0,0) = 0.0;
  yein(1,1) = 0.0;

  yuin(0,1) = 0.0; yuin(0,2) = 0.0;
  yuin(1,0) = 0.0; yuin(1,2) = 0.0;
  yuin(2,0) = 0.0; yuin(2,1) = 0.0;

  ydin(0,1) = 0.0; ydin(0,2) = 0.0;
  ydin(1,0) = 0.0; ydin(1,2) = 0.0;
  ydin(2,0) = 0.0; ydin(2,1) = 0.0;

  yein(0,1) = 0.0; yein(0,2) = 0.0;
  yein(1,0) = 0.0; yein(1,2) = 0.0;
  yein(2,0) = 0.0; yein(2,1) = 0.0;

  // Make sure off diagonals for kappa, T_kappa, and the inert lambdas should vanish
  kappain(0,1) = 0.0;
  kappain(0,2) = 0.0;
  kappain(1,0) = 0.0;
  kappain(1,2) = 0.0;
  kappain(2,0) = 0.0;
  kappain(2,1) = 0.0;

  lambda12in(0,1) = 0.0;
  lambda12in(1,0) = 0.0;

  double v1 = hvev/Sqrt(1.0+tanb*tanb);
  double v2 = v1*tanb;

  if (ENABLE_DEBUG)
    {
      cout << "# Initialising SUSY parameters..." << endl;
    }

  // SUSY parameters
  genericE6SSM_susy_parameters r_susy = genericE6SSM_susy_parameters(mx, l, t, input, ydin, yein, kappain, lambda12in, lambda3in, yuin, mupr, g1in, g2in, g3in, gNin, v1, v2, svev);

  // In the pE6SSM, the first and second generation A terms should vanish, as
  // should the off-diagonal mixings.
  tuin(0,0) = 0.0; tuin(1,1) = 0.0;
  tdin(0,0) = 0.0; tdin(1,1) = 0.0;
  tein(0,0) = 0.0; tein(1,1) = 0.0;

  tuin(0,1) = 0.0; tuin(0,2) = 0.0;
  tuin(1,0) = 0.0; tuin(1,2) = 0.0;
  tuin(2,0) = 0.0; tuin(2,1) = 0.0;

  tdin(0,1) = 0.0; tdin(0,2) = 0.0;
  tdin(1,0) = 0.0; tdin(1,2) = 0.0;
  tdin(2,0) = 0.0; tdin(2,1) = 0.0;

  tein(0,1) = 0.0; tein(0,2) = 0.0;
  tein(1,0) = 0.0; tein(1,2) = 0.0;
  tein(2,0) = 0.0; tein(2,1) = 0.0;

  tlambda12in(0,1) = 0.0;
  tlambda12in(1,0) = 0.0;

  tkappain(0,1) = 0.0;
  tkappain(0,2) = 0.0;
  tkappain(1,0) = 0.0;
  tkappain(1,2) = 0.0;
  tkappain(2,0) = 0.0;
  tkappain(2,1) = 0.0;

  // Make sure flavour non-diagonal masses are zero
  mHdISq(0,1) = 0.0; mHdISq(1,0) = 0.0;
  mHuISq(0,1) = 0.0; mHuISq(1,0) = 0.0;
  mSISq(0,1) = 0.0; mSISq(1,0) = 0.0;

  mQlSq(0,1) = 0.0; mQlSq(0,2) = 0.0;
  mQlSq(1,0) = 0.0; mQlSq(1,2) = 0.0;
  mQlSq(2,0) = 0.0; mQlSq(2,1) = 0.0;

  mLlSq(0,1) = 0.0; mLlSq(0,2) = 0.0;
  mLlSq(1,0) = 0.0; mLlSq(1,2) = 0.0;
  mLlSq(2,0) = 0.0; mLlSq(2,1) = 0.0;

  mUrSq(0,1) = 0.0; mUrSq(0,2) = 0.0;
  mUrSq(1,0) = 0.0; mUrSq(1,2) = 0.0;
  mUrSq(2,0) = 0.0; mUrSq(2,1) = 0.0;

  mDrSq(0,1) = 0.0; mDrSq(0,2) = 0.0;
  mDrSq(1,0) = 0.0; mDrSq(1,2) = 0.0;
  mDrSq(2,0) = 0.0; mDrSq(2,1) = 0.0;

  mErSq(0,1) = 0.0; mErSq(0,2) = 0.0;
  mErSq(1,0) = 0.0; mErSq(1,2) = 0.0;
  mErSq(2,0) = 0.0; mErSq(2,1) = 0.0;

  mDxSq(0,1) = 0.0; mDxSq(0,2) = 0.0;
  mDxSq(1,0) = 0.0; mDxSq(1,2) = 0.0;
  mDxSq(2,0) = 0.0; mDxSq(2,1) = 0.0;

  mDxbarSq(0,1) = 0.0; mDxbarSq(0,2) = 0.0;
  mDxbarSq(1,0) = 0.0; mDxbarSq(1,2) = 0.0;
  mDxbarSq(2,0) = 0.0; mDxbarSq(2,1) = 0.0;

  if (ENABLE_DEBUG)
    {
      cout << "# Initialising soft SUSY breaking parameters..." << endl;
    }

  // Soft SUSY breaking parameters. Note ordering of gauginos is 1 = M_1,
  // 2 = M_2, 3 = M_3 and 4 = M_1'
  genericE6SSM_soft_parameters r_soft = genericE6SSM_soft_parameters(r_susy, tdin, tein, tkappain, tlambda12in, tlambda3in, 
								     tuin, Bmupr, mQlSq, mLlSq, mHdSq, mHuSq, mDrSq, mUrSq, mErSq,
								     mSSq, mHdISq, mHuISq, mSISq, mDxSq, mDxbarSq, mHpSq, 
								     mHpbarSq, MBin, MWBin, MGin, 
								     MBpin);

  // Somehow need to construct the full model here
  genericE6SSM<Two_scale> r_approx = genericE6SSM<Two_scale>(r_soft);

  // M_{SUSY}, m_Hd^2(MX), m_s^2(MX) and m_Hu^2(MX) must be determined
  // simultaneously using an iterative procedure. 

  DoubleVector solnGuess(4);

  solnGuess(1) = mHdSq;
  solnGuess(2) = mHuSq;
  solnGuess(3) = mSSq;
  solnGuess(4) = ms;

  r_approx.set_mHd2(solnGuess(1));
  r_approx.set_mHu2(solnGuess(2));
  r_approx.set_ms2(solnGuess(3));


  // Try to determine the correct solution using Newton's method, 
  // catch any numerical errors and flag as a problem point if they occur.
  try
    {
      if (ENABLE_DEBUG)
	{
	  cout << "# Imposing EWSB conditions on soft masses..." << endl;
	}

      hasEWSBProblem = ESSM_ImplementEWSBConstraints_SoftMasses(r_approx, mx, ms, 
								false, solnGuess, tol);

      
      if (hasEWSBProblem)
	{
	  cerr << "WARNING: problem implementing EWSB conditions." << endl;
	  genericE6SSM<Two_scale> w = r_approx;
	  w.set_mHd2(solnGuess(1));
	  w.set_mHu2(solnGuess(2));
	  w.set_ms2(solnGuess(3));
	  cerr << "m_Hd^2 = " << w.get_mHd2() << endl;
	  cerr << "m_Hu^2 = " << w.get_mHu2() << endl;
	  cerr << "m_s^2 = " << w.get_ms2() << endl;
	  cerr << "M_{SUSY} = " << solnGuess(4) << endl;
	  w.run_to(solnGuess(4), PRECISION);
	  cerr << "f1 = " << ESSM_EWSBCondition1(w) << endl;
	  cerr << "f2 = " << ESSM_EWSBCondition2(w) << endl;
	  cerr << "f3 = " << ESSM_EWSBCondition3(w) << endl;
	}
      
      mHdSq = solnGuess(1);
      mHuSq = solnGuess(2);
      mSSq = solnGuess(3);
      ms = solnGuess(4);
      
      r_approx.set_mHd2(solnGuess(1));
      r_approx.set_mHu2(solnGuess(2));
      r_approx.set_ms2(solnGuess(3));
      
      // Try to get the spectrum. Catch any errors and flag as 
      // problem if they occur.
      // Get DR bar stop and sbottom masses.
      r_approx.run_to(solnGuess(4), PRECISION);

      if (ENABLE_DEBUG)
	{
	  cout << "# Calculating squark masses..." << endl;
	}      

      r_approx.calculate_MStop();
      r_approx.calculate_MSbot();

      if (ENABLE_DEBUG) 
	{
	  cout << "# Checking for tachyons..." << endl;
	}

      // Check for tachyonic squakrs      
      if (r_approx.get_problems().is_tachyon(genericE6SSM_info::Su) || r_approx.get_problems().is_tachyon(genericE6SSM_info::Sd))
	{
	  squarksTachyons = true;
	  if (ENABLE_DEBUG)
	    {
	      if (r_approx.get_problems().is_tachyon(genericE6SSM_info::Su)) cout << "#    tachyonic up-type squarks" << endl;
	      if (r_approx.get_problems().is_tachyon(genericE6SSM_info::Sd)) cout << "#    tachyonic down-type squarks" << endl;
	    }
	}

      // Values needed below, shouldn't be any problems here so long as tan(beta), v, and the gauge couplings
      // are sensible when input
      double v1 = r_approx.get_vd();
      double v2 = r_approx.get_vu();
      double s = r_approx.get_vs();

      double v = Sqrt(v1*v1+v2*v2);
      double tb = v2/v1;

      if (ENABLE_DEBUG)
	{
	  cout << "# Calculating quark masses..." << endl;
	}

      double sw = Sin(r_approx.ThetaW());
      double mt = r_approx.get_Yu(2, 2)*v*Sin(ArcTan(tb))/Sqrt(2.0);
      double mb = r_approx.get_Yd(2, 2)*v*Cos(ArcTan(tb))/Sqrt(2.0);

      if (ENABLE_DEBUG)
	{
	  cout << "# Calculating vector boson masses..." << endl;
	}

      r_approx.calculate_MVZ();
      r_approx.calculate_MVWm();
      r_approx.calculate_MVZp();

      if (ENABLE_DEBUG) 
	{
	  cout << "# Checking for tachyons..." << endl;
	}

      if (r_approx.get_problems().is_tachyon(genericE6SSM_info::VZ) || r_approx.get_problems().is_tachyon(genericE6SSM_info::VZp) 
	  || r_approx.get_problems().is_tachyon(genericE6SSM_info::VWm) )
	{
	  vectorBosonsTachyons = true;
	  if (ENABLE_DEBUG)
	    {
	      if (r_approx.get_problems().is_tachyon(genericE6SSM_info::VZ)) cout << "#    tachyonic Z" << endl;
	      if (r_approx.get_problems().is_tachyon(genericE6SSM_info::VZp)) cout << "#    tachyonic Z'" << endl;
	      if (r_approx.get_problems().is_tachyon(genericE6SSM_info::VWm)) cout << "#    tachyonic W" << endl;
	    }
	}
      

      // Calculate neutralino and chargino DR bar masses, there shouldn't be any problems here 
      // (or if there are they will be serious numerical ones in diagonalisation)
      if (ENABLE_DEBUG)
	{
	  cout << "# Calculating neutralino and chargino masses..." << endl;
	}

      r_approx.calculate_MChi();
      r_approx.calculate_MCha();  

      // Calculate DR bar Higgs masses using FlexibleSUSY routines...
      if (ENABLE_DEBUG)
	{
	  cout << "# Calculating DR bar Higgs masses..." << endl;
	}
      r_approx.calculate_Mhh();
      r_approx.calculate_MAh();
      r_approx.calculate_MHpm();
      
      if (ENABLE_DEBUG) 
	{
	  cout << "# Checking for tachyons..." << endl;
	}

      // Check if DR bar h0 or A0 are tachyons
      if (r_approx.get_problems().is_tachyon(genericE6SSM_info::hh) || r_approx.get_problems().is_tachyon(genericE6SSM_info::Ah) || 
	  r_approx.get_problems().is_tachyon(genericE6SSM_info::Hpm))
	{
	  higgsTachyons = true;
	  if (ENABLE_DEBUG)
	    {
	      if (r_approx.get_problems().is_tachyon(genericE6SSM_info::hh)) cout << "#    tachyonic hh" << endl;
	      if (r_approx.get_problems().is_tachyon(genericE6SSM_info::Ah)) cout << "#    tachyonic Ah" << endl;
	      if (r_approx.get_problems().is_tachyon(genericE6SSM_info::Hpm)) cout << "#    tachyonic Hpm" << endl;
	    }
	}

      if (ENABLE_DEBUG)
	{
	  cout << "# Calculating tadpole corrections..." << endl;
	}      

      // ... but use approximate expressions here for loop corrected Higgs masses
      DoubleVector mstop(2), mstopsq(2), mD1sq(3), mD2sq(3);
      physical_ESSM(r_approx, mstop, mstopsq, mD1sq, mD2sq, s, tb);      

      if (ENABLE_DEBUG)
	{
	  cout << "#    using m_~t1 = " << mstop(1) << ", m_~t2 = " << mstop(2) << endl;
	  cout << "#    using m_~t1^2 = " << mstopsq(1) << ", m_~t2^2 = " << mstopsq(2) << endl;
	}

      // In HiggsMasses, sing = 1 indicates inaccurate Higgs mass.
      // ExpValid = 15 indicates TADPOLEPROBLEM.
      // ExpValid = 30 indicates not experimentally valid, but otherwise no problems
      // WhatCorrections selects which version to use (use 1)
      int sing = 0;
      int ExpValid = 0;
      const int WhatCorrections = 1;
      DoubleVector expBounds(2);

      expBounds(1) = HIGGSCENT - HIGGSERROR;
      expBounds(2) = HIGGSCENT + HIGGSERROR;

      DoubleVector mh(3);
      DoubleMatrix mhmix(3,3), msq(3,3);

      if (ENABLE_DEBUG)
	{
	  cout << "# Calculating approximate CP-even Higgs pole masses..." << endl;
	}

      poleHiggsTachyons = HiggsMasses(r_approx, s, tb, mstop, mstopsq, WhatCorrections, 
				      false, false, expBounds, ExpValid, mh, mhmix, msq, sing);

      if (ENABLE_DEBUG)
	{
	  cout << "# Checking for problems with pole CP-even Higgs masses..." << endl;
	}

      // Other problems that can occur in physical Higgs mass calculation:
      // inaccurate result
      if (ExpValid == HIGGSPROBLEM || sing == 1)
	{
	  inaccurateHiggsMass = true;
	}
      if (ExpValid == TADPOLESPROBLEM)
	{
	  tadpoleProblem = true;
	  if (ENABLE_DEBUG)
	    {
	      cout << "#    problem with tadpoles" << endl;
	    }
	}
      if (ExpValid == POLEHIGGSTACHYON)
	{
	  poleHiggsTachyons = true;
	  if (ENABLE_DEBUG)
	    {
	      cout << "#    tachyonic m_h^2" << endl;
	    }
	}

      if (ENABLE_DEBUG)
	{
	  cout << "# Calculating approximate CP-odd Higgs pole mass..." << endl;
	}

      // Also calculate the 1-loop mass of the CP-odd Higgs, which seems
      // to be an important constraint on parameter space. Use approximate
      // expression for speed.
      int cpOddproblem = 0;
      double mAsq = mAsq_OneLoop(r_approx, s, tb, cpOddproblem);

      if (ENABLE_DEBUG)
	{
	  cout << "# Checking for problems with pole CP-odd Higgs mass..." << endl;
	}

      if (mAsq < 0.0)
	{
	  poleHiggsTachyons = true;
	  if (ENABLE_DEBUG)
	    {
	      cout << "#    tachyonic m_A^2" << endl;
	    }
	}
      if (cpOddproblem == TADPOLESPROBLEM)
	{
	  tadpoleProblem = true;
	  if (ENABLE_DEBUG)
	    {
	      cout << "#    problem with tadpoles" << endl;
	    }
	}

      if (ENABLE_DEBUG)
	{
	  cout << "# Setting approximate Higgs pole masses..." << endl;
	}

      // Update the object at M_{SUSY} with the new values for the pole masses.
      // (added in a method to do this directly). Later consider just using
      // the FlexibleSUSY routines to get the pole masses.
      r_approx.set_Mhh_pole(mh(1), mh(2), mh(3));
      r_approx.set_MAh_pole(ZeroSqrt(mAsq));

      if (ENABLE_DEBUG)
	{
	  cout << "#    Mhh(0) = " << r_approx.get_physical().Mhh(0)
	       << ", Mhh(1) = " << r_approx.get_physical().Mhh(1)
	       << ", Mhh(2) = " << r_approx.get_physical().Mhh(2) << endl;
	  cout << "#    MAh(0) = " << r_approx.get_physical().MAh(0) 
	       << ", MAh(1) = " << r_approx.get_physical().MAh(1) 
	       << ", MAh(2) = " << r_approx.get_physical().MAh(2) << endl;
	}

    }
  // DH::Check that there are no other possible failure conditions that may occur
  catch(const string & a) 
    { 
      cerr << "WARNING: serious numerical problem at point." << endl; 
      hasSeriousProblem = true;
    }
  catch(const char * a) 
    { 
      cerr << "WARNING: serious numerical problem at point." << endl; 
      hasSeriousProblem = true; 
    }
  catch(...) 
    { 
      cerr << "WARNING: unknown type of exception caught: in doSimplifiedSpectrum." << endl; 
      hasSeriousProblem = true; 
    }

  if (ENABLE_DEBUG)
    {
      cout << "# Finished spectrum calculation..." << endl;
    }  

  return r_approx;
}

double get_softAu(flexiblesusy::genericE6SSM<flexiblesusy::Two_scale> const & m, int i, int j)
{
  using namespace flexiblesusy;

  double tu = m.get_TYu(i,j);
  double yu = m.get_Yu(i,j);
  double au;
  
  if (Abs(tu) < EPSTOL) 
    {
      au = 0.0;
      return au;
    }
  
  if (Abs(yu) < 1.0e-100)
    {
      ostringstream ii;
      ii << "WARNING: trying to calculate A_u(" << i << ", " << j << ") where Y_u(" << i << ", " << j << ") coupling is " <<
	Abs(yu) << endl;
      throw ii.str();
    }
  else
    {
      au = tu/yu;
    }

  return au;
  
}

double get_softAd(flexiblesusy::genericE6SSM<flexiblesusy::Two_scale> const & m, int i, int j)
{
  using namespace flexiblesusy;

  double td = m.get_TYd(i,j);
  double yd = m.get_Yd(i,j);
  double ad;
  
  if (Abs(td) < EPSTOL) 
    {
      ad = 0.0;
      return ad;
    }
  
  if (Abs(yd) < 1.0e-100)
    {
      ostringstream ii;
      ii << "WARNING: trying to calculate A_d(" << i << ", " << j << ") where Y_d(" << i << ", " << j << ") coupling is " <<
	Abs(yd) << endl;
      throw ii.str();
    }
  else
    {
      ad = td/yd;
    }

  return ad;

}

double get_softAe(flexiblesusy::genericE6SSM<flexiblesusy::Two_scale> const & m, int i, int j)
{
  using namespace flexiblesusy;

  double te = m.get_TYe(i,j);
  double ye = m.get_Ye(i,j);
  double ae;
  
  if (Abs(te) < EPSTOL) 
    {
      ae = 0.0;
      return ae;
    }
  
  if (Abs(ye) < 1.0e-100)
    {
      ostringstream ii;
      ii << "WARNING: trying to calculate A_e(" << i << ", " << j << ") where Y_e(" << i << ", " << j << ") coupling is " <<
	Abs(ye) << endl;
      throw ii.str();
    }
  else
    {
      ae = te/ye;
    }

  return ae;
}
