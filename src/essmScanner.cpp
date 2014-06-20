/*
  Scan over the E6SSM parameter space and calculate the fine tuning 
  at each point.
  Modified to use the FlexibleSUSY E6SSM spectrum generator instead of
  SOFTSUSY routines.
 */

#include "./E6SSM_Spectrum_Generators/models/genericE6SSM/genericE6SSM_input_parameters.hpp"
#include "./E6SSM_Spectrum_Generators/models/genericE6SSM/genericE6SSM_spectrum_generator.hpp"
#include "./E6SSM_Spectrum_Generators/models/genericE6SSM/genericE6SSM_info.hpp"
#include "./E6SSM_Spectrum_Generators/models/genericE6SSM/genericE6SSM_two_scale_model.hpp"
#include "./E6SSM_Spectrum_Generators/src/error.hpp"
#include "./E6SSM_Spectrum_Generators/src/spectrum_generator_settings.hpp"
#include "./E6SSM_Spectrum_Generators/src/lowe.h"
#include "./E6SSM_Spectrum_Generators/src/command_line_options.hpp"
#include "./E6SSM_Spectrum_Generators/src/wrappers.hpp"
#include "./E6SSM_Spectrum_Generators/src/problems.hpp"

#include <iostream>
#include <cstdlib>
#include <sys/time.h>


#include "flags.h"
#include "tuningnumerics.h"
#include "tuningutils.h"
#include "essmtuningutils.h"

void errorCall();
string ToUpper(const string & s);

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

const double QQY = flexiblesusy::Sqrt(0.6)*(1.0/6.0);
const double QuY = flexiblesusy::Sqrt(0.6)*(-2.0/3.0);
const double QdY = flexiblesusy::Sqrt(0.6)*(1.0/3.0);
const double QLY = flexiblesusy::Sqrt(0.6)*(-1.0/2.0);
const double QeY = flexiblesusy::Sqrt(0.6);
const double QNY = 0.0;
const double QSY = 0.0;
const double QH1Y = flexiblesusy::Sqrt(0.6)*(-1.0/2.0);
const double QH2Y = flexiblesusy::Sqrt(0.6)*(1.0/2.0);
const double QXY = flexiblesusy::Sqrt(0.6)*(-1.0/3.0);
const double QXbarY = flexiblesusy::Sqrt(0.6)*(1.0/3.0);
const double QHPrY = flexiblesusy::Sqrt(0.6)*(-1.0/2.0);
const double QHbarPrY = flexiblesusy::Sqrt(0.6)*(1.0/2.0);
  
  /*
    --------------------------------------------------------------
  */

int main(int argc, const char* argv[])
{
  using namespace flexiblesusy;
  using namespace softsusy;
  using namespace essm_tuning_utils;
  typedef Two_scale algorithm_type;

  /*
    --------------------------------------------------------------
    SM Parameters
    --------------------------------------------------------------
   */
  double poleM_t = 173.2; // top quark pole mass in GeV
  double mbmb = 4.16; // bottom quark running mass m_b(m_b)^MSbar in GeV
  double poleM_tau = 1.777;// tau pole mass in GeV
  double poleM_Z = 91.1876; // Z boson pole mass in GeV
  double G_F = 1.16637900e-5; // G_F^MSbar in GeV^-2.
  double alphasmz = 0.1193; // alpha_s(mz)^MSbar (strong coupling constant = g_3^2/(4*PI))
  double alphaemmzinv = 127.9568; // 1 / alpha_em(mz)^MSbar (electromagnetic coupling = e^2/(4*PI))

  // DH:: neglect these and set them to zero, as the first and second generation
  // Yukawas are assumed to be zero.
  // double poleM_nu3 = 0.0;
  // double poleM_electron = 5.109989020e-04;
  // double poleM_muon = 1.056583570e-01;
  // double poleM_nu1 = 0.0;
  // double poleM_nu2 = 0.0;
  // double mdAt2GeV = 4.750000000e-03;
  // double muAt2GeV = 2.400000000e-03;
  // double msAt2GeV = 1.040000000e-01;
  // double mcAt2GeV = 1.270000000e+00;

  /*
    --------------------------------------------------------------
   */

  /*
    --------------------------------------------------------------
    E6 mixing angle (used to select model)
    --------------------------------------------------------------
   */
  // DH:: not varied at the moment
  //double thetaE6 = atan(Sqrt(15.0));

  /*
    --------------------------------------------------------------
   */

  /*
    --------------------------------------------------------------
    Default high scale input parameters (not currently used)
    --------------------------------------------------------------
   */

  double m0 = 100.0; // GeV
  double m12 = 500.0; // GeV
  double A0 = 500.0; // GeV
  double lambda0 = 1.0; 
  double kappa0 = 1.0;
  double mupr0 = 1000.0; // GeV
  double Bmupr0 = 1000.0; // GeV^2

  /*
    --------------------------------------------------------------
   */

  /*
    --------------------------------------------------------------
    Values for parameters that are not scanned over
    --------------------------------------------------------------
   */

  double MX = 20000.0; // GeV
  
  double M1 = 300.0; // GeV
  double M1p = 300.0; // GeV
  double M3 = 2000.0; // GeV
  
  double Ab = 0.0; // GeV
  double Atau = 0.0; // GeV
  
  double meL = 5000.0; // GeV
  double mmuL = 5000.0; // GeV
  double mtauL = 5000.0; // GeV
  
  double meR = 5000.0; // GeV
  double mmuR = 5000.0; // GeV
  double mtauR = 5000.0; // GeV
  
  double mqL1 = 5000.0; // GeV
  double mqL2 = 5000.0; // GeV

  double muR = 5000.0; // GeV
  double mcR = 5000.0; // GeV

  double mdR = 5000.0; // GeV
  double msR = 5000.0; // GeV
  double mbR = 5000.0; // GeV

  // Universal kappas
  double kappa1 = 0.6; // DH:: check
  double kappa2 = 0.6; // DH:: check
  double kappa3 = 0.6; // DH:: check

  double Akappa1 = 0.0; // GeV 
  double Akappa2 = 0.0; // GeV 
  double Akappa3 = 0.0; // GeV 

  // Off-diagonals zero as well
  // Choose diagonals so that inert Higgsinos
  // have mass >~ 1 TeV i.e. roughly 0.2
  double lambda1 = 0.2; 
  double lambda2 = 0.2; 

  double Alambda1 = 0.0; // GeV 
  double Alambda2 = 0.0; // GeV 

  double mupr = 5000.0; // GeV 
  double Bmupr = 5000.0; // GeV^2 

  double mH11Sq = Sqr(5000.0); // GeV^2 
  double mH12Sq = Sqr(5000.0); // GeV^2 

  double mH21Sq = Sqr(5000.0); // GeV^2 
  double mH22Sq = Sqr(5000.0); // GeV^2 

  double mS1Sq = Sqr(5000.0); // GeV^2 
  double mS2Sq = Sqr(5000.0); // GeV^2 

  double mHpSq = Sqr(5000.0); // GeV^2 
  double mHpbarSq = Sqr(5000.0); // GeV^2 

  double mDSq = Sqr(5000.0); // GeV^2
  double mDbarSq = Sqr(5000.0); // GeV^2

  /*
    --------------------------------------------------------------
   */

  /*
    --------------------------------------------------------------
    Values for scanned parameters
    --------------------------------------------------------------
   */

  double tanb = 10.0;

  double lambda3 = 3.0;
  double Alambda3 = 5000.0; // GeV
  double mqL3 = 500.0; // GeV
  double mtR = 500.0; // GeV
  double At = 1000.0; // GeV

  double M2 = 1000.0; // GeV

  // NB the scanned parameters are the squared masses
  // m_Q3^2 and m_u3^2, not m_Q3 and m_u3.  
  double mql3sq = sqr(mqL3);
  double mtRsq = sqr(mtR);  

  /*
    --------------------------------------------------------------
   */
  
  /*
    ----------------------------------------------------
    Scan parameters - upper and lower bounds, 
    number of points in each parameter direction,
    and increment step size.
    ----------------------------------------------------
   */

  // For the purposes of doing large scans, 
  // where we want to divide parameter space up into
  // many small subsections, it is easiest if these
  // are provided as inputs in a file.

  // Default number of points in each direction
  int tb_npts = 1;
  int lambda3_npts = 1;
  int Alambda3_npts = 1;
  int mqL3sq_npts = 1;
  int mtRsq_npts = 1;
  int At_npts = 1;
  int M2_npts = 1;
  
  // Default lower bounds
  double tb_low = 2.0;
  double lambda3_low = -1000.0; // GeV
  double Alambda3_low = -1000.0; // GeV
  double At_low = -10000.0; // GeV
  double mqL3_low = 200.0; // GeV
  double mtR_low = 200.0; // GeV
  
  double mqL3sq_low = sqr(mqL3_low);
  double mtRsq_low = sqr(mtR_low);
  
  double M2_low = 100.0; // GeV

  // Default upper bounds
  double tb_up = 50.0;
  double lambda3_up = 1000.0; // GeV
  double Alambda3_up = 1000.0; // GeV
  double At_up = 10000.0; // GeV
  double mqL3_up = 2000.0; // GeV
  double mtR_up = 2000.0; // GeV
  
  double mqL3sq_up = sqr(mqL3_up);
  double mtRsq_up = sqr(mtR_up);
  
  double M2_up = 2000.0; // GeV

  // Default step sizes
  double tb_incr = 0.0;
  double lambda3_incr = 0.0;
  double Alambda3_incr = 0.0;
  double mqL3sq_incr = 0.0;
  double mtRsq_incr = 0.0;
  double At_incr = 0.0;

  double M2_incr = 0.0;

  /*
    ----------------------------------------------------
   */

  /*
    ----------------------------------------------------
    Flags and constants
    ----------------------------------------------------
   */

  bool const USEAPPROXSOLNS = true; //< use approximate solutions to RGEs

  bool uselogscanAt = false; //< use a symmetric log-scan over A_t
  bool uselogscanAlambda3 = false; //< use a symmetric log-scan over A_lambda3

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


  /*
    ----------------------------------------------------
   */


  /*
    ----------------------------------------------------
    Actual calculations begin here. For ordinary scans, 
    you shouldn't have to touch any of the code below
    (except for debugging, of course). All of the 
    controls for driving the program are before this point.
    ----------------------------------------------------
   */

  /*
    ----------------------------------------------------
    Saved couplings and Higgs VEV
    ----------------------------------------------------
   */
  Eigen::Matrix<double,3,3> yuin, ydin, yein; //< Yukawas
  double g1in = 0.4, g2in = 0.45, g3in = 0.9, gNin = 0.4;
  double vin = 246.0; //< Higgs vev, GeV

  // Need M_Z' at least 2.5 TeV
  double vsin = 6700.0; //< singlet vev, GeV

  /*
    ----------------------------------------------------
   */

  /*
    ----------------------------------------------------
    Miscellaneous variable declarations. Any extra variables
    that are needed in the rest of the calculation should be
    here. With only a few exceptions, there shouldn't be any 
    new variables declared in the body of the scan 
    (a la FORTRAN rules, to hopefully avoid any local redefinitions 
    and therefore reduce bugs).
    ----------------------------------------------------
   */

  double lambda3in;
  Eigen::Matrix<double,2,2> lambda12in;
  lambda12in(0,0) = lambda1;
  lambda12in(1,1) = lambda2;
  lambda3in = lambda3;

  lambda12in(0,1) = 0.0; lambda12in(1,0) = 0.0;

  Eigen::Matrix<double,3,3> kappain;
  kappain(0,0) = kappa1;
  kappain(1,1) = kappa2;
  kappain(2,2) = kappa3;

  kappain(0,1) = 0.0; kappain(0,2) = 0.0;
  kappain(1,0) = 0.0; kappain(1,2) = 0.0;
  kappain(2,0) = 0.0; kappain(2,1) = 0.0;

  // Soft squared masses
  Eigen::Matrix<double,3,3> mQlSq, mUrSq, mDrSq, mLlSq, mErSq;

  mLlSq(0,0) = sqr(meL);
  mLlSq(1,1) = sqr(mmuL);
  mLlSq(2,2) = sqr(mtauL);

  mErSq(0,0) = sqr(meR);
  mErSq(1,1) = sqr(mmuR);
  mErSq(2,2) = sqr(mtauR);

  mQlSq(0,0) = sqr(mqL1);
  mQlSq(1,1) = sqr(mqL2);
  mQlSq(2,2) = sqr(mqL3);

  mUrSq(0,0) = sqr(muR);
  mUrSq(1,1) = sqr(mcR);
  mUrSq(2,2) = sqr(mtR);

  mDrSq(0,0) = sqr(mdR);
  mDrSq(1,1) = sqr(msR);
  mDrSq(2,2) = sqr(mbR);

  // Make sure off-diagonals are zero, as is assumed
  // in the cE6SSM

  mLlSq(0,1) = 0.0;
  mLlSq(0,2) = 0.0;
  mLlSq(1,0) = 0.0;
  mLlSq(1,2) = 0.0;
  mLlSq(2,0) = 0.0;
  mLlSq(2,1) = 0.0;

  mErSq(0,1) = 0.0;
  mErSq(0,2) = 0.0;
  mErSq(1,0) = 0.0;
  mErSq(1,2) = 0.0;
  mErSq(2,0) = 0.0;
  mErSq(2,1) = 0.0;

  mQlSq(0,1) = 0.0;
  mQlSq(0,2) = 0.0;
  mQlSq(1,0) = 0.0;
  mQlSq(1,2) = 0.0;
  mQlSq(2,0) = 0.0;
  mQlSq(2,1) = 0.0;

  mUrSq(0,1) = 0.0;
  mUrSq(0,2) = 0.0;
  mUrSq(1,0) = 0.0;
  mUrSq(1,2) = 0.0;
  mUrSq(2,0) = 0.0;
  mUrSq(2,1) = 0.0;

  mDrSq(0,1) = 0.0;
  mDrSq(0,2) = 0.0;
  mDrSq(1,0) = 0.0;
  mDrSq(1,2) = 0.0;
  mDrSq(2,0) = 0.0;
  mDrSq(2,1) = 0.0;

  // Soft trilinears = y*A
  Eigen::Matrix<double,3,3> tuin, tdin, tein, tkappain;
  Eigen::Matrix<double,2,2> tlambda12in;
  double tlambda3in;

  tuin(0,0) = 0.0; tuin(1,1) = 0.0; tuin(2,2) = yuin(2,2) * At;
  tdin(0,0) = 0.0; tdin(1,1) = 0.0; tdin(2,2) = ydin(2,2) * Ab;
  tein(0,0) = 0.0; tein(1,1) = 0.0; tuin(2,2) = yein(2,2) * Atau;

  tuin(0,1) = 0.0; tuin(0,2) = 0.0;
  tuin(1,0) = 0.0; tuin(1,2) = 0.0;
  tuin(2,0) = 0.0; tuin(2,1) = 0.0;

  tdin(0,1) = 0.0; tdin(0,2) = 0.0;
  tdin(1,0) = 0.0; tdin(1,2) = 0.0;
  tdin(2,0) = 0.0; tdin(2,1) = 0.0;

  tein(0,1) = 0.0; tein(0,2) = 0.0;
  tein(1,0) = 0.0; tein(1,2) = 0.0;
  tein(2,0) = 0.0; tein(2,1) = 0.0;

  tkappain(0,0) = kappa1 * Akappa1;
  tkappain(1,1) = kappa2 * Akappa2;
  tkappain(2,2) = kappa3 * Akappa3;

  tkappain(0,1) = 0.0; tkappain(0,2) = 0.0;
  tkappain(1,0) = 0.0; tkappain(1,2) = 0.0;
  tkappain(2,0) = 0.0; tkappain(2,1) = 0.0;

  tlambda12in(0,0) = lambda1 * Alambda1;
  tlambda12in(1,1) = lambda2 * Alambda2;
  tlambda12in(0,1) = 0.0; tlambda12in(1,0) = 0.0;

  tlambda3in = lambda3 * Alambda3;

  // Exotic soft scalar masses
  Eigen::Matrix<double,3,3> mDxSq, mDxbarSq;

  mDxSq(0,0) = mDSq;
  mDxSq(1,1) = mDSq;
  mDxSq(2,2) = mDSq;

  mDxbarSq(0,0) = mDbarSq;
  mDxbarSq(1,1) = mDbarSq;
  mDxbarSq(2,2) = mDbarSq;

  mDxSq(0,1) = 0.0;
  mDxSq(0,2) = 0.0;
  mDxSq(1,0) = 0.0;
  mDxSq(1,2) = 0.0;
  mDxSq(2,0) = 0.0;
  mDxSq(2,1) = 0.0;

  mDxbarSq(0,1) = 0.0;
  mDxbarSq(0,2) = 0.0;
  mDxbarSq(1,0) = 0.0;
  mDxbarSq(1,2) = 0.0;
  mDxbarSq(2,0) = 0.0;
  mDxbarSq(2,1) = 0.0;

  // Soft Higgs and singlet masses
  double mHdSq = 1.0e6;
  double mHuSq = 1.0e6;
  double mSSq = 1.0e6;

  Eigen::Matrix<double,2,2> mHdISq, mHuISq, mSISq;
  mHdISq(0,0) = mH11Sq;
  mHdISq(1,1) = mH12Sq;
  mHuISq(0,0) = mH21Sq;
  mHuISq(1,1) = mH22Sq;
  mSISq(0,0) = mS1Sq;
  mSISq(1,1) = mS2Sq;

  // No mixing...
  mHdISq(0,1) = 0.0;
  mHdISq(1,0) = 0.0;
  mHuISq(0,1) = 0.0;
  mHuISq(1,0) = 0.0;
  mSISq(0,1) = 0.0;
  mSISq(1,0) = 0.0;

  // M_{SUSY} = \sqrt{m_{~t_1}^DRbar(M_{SUSY})*m_{~t_2}^DRbar(M_{SUSY})}
  double MS = 1.0e3;

  // SM inputs...
  QedQcd oneset; 

  oneset.setAlpha(ALPHA, 1.0 / alphaemmzinv);			 
  oneset.setAlpha(ALPHAS, alphasmz);
  oneset.setMu(poleM_Z); MZ = poleM_Z;
  oneset.setMass(mBottom, mbmb);
  oneset.setMbMb(mbmb);
  oneset.setPoleMt(poleM_t);
  oneset.setMass(mTau, poleM_tau); 
  oneset.setPoleMtau(poleM_tau);

  QedQcd dataset = oneset;

  oneset.toMz();

  // Vectors for storing values of scanned variables
  Eigen::VectorXd tb_vals(1), lambda3_vals(1), Alambda3_vals(1), mqL3sq_vals(1), mtRsq_vals(1), At_vals(1), M2_vals(1);

  // Variables to store values of couplings at M_{SUSY}
  double lambda3AtMsusy, kappa3AtMsusy;

  // Vector for storing fine tunings
  Eigen::Matrix<double,tuning_parameters::NUMESSMTUNINGPARS,1> tunings;

  // Variables to store values of EWSB conditions
  double f1, f2, f3;

  // Models that we will calculate the fine tuning for
  genericE6SSM<Two_scale> m;

  /*
    ----------------------------------------------------
   */

  /*
    ----------------------------------------------------
    Flags for doing input file parsing, to check that
    the inputs are valid
    ----------------------------------------------------
   */

  // Flags to indicate whether the first option for
  // A_t requests a log or linear scan
  bool hasLogAtOption = false;
  bool hasNonLogAtOption = false;
  bool hasLogAlambda3Option = false;
  bool hasNonLogAlambda3Option = false;

  // Flags to indicate lower bounds are found
  bool hasTanbLL = false;
  bool hasLambda3LL = false;
  bool hasAlambda3LL = false;
  bool hasmqL3sqLL = false;
  bool hasmtRsqLL = false;
  bool hasAtLL = false;
  bool hasM2LL = false;

  // Flags to indicate upper bounds are found
  bool hasTanbUL = false;
  bool hasLambda3UL = false;
  bool hasAlambda3UL = false;
  bool hasmqL3sqUL = false;
  bool hasmtRsqUL = false;
  bool hasAtUL = false;
  bool hasM2UL = false;

  // Flags to indicate number of points are found
  bool hasTanbNpts = false;
  bool hasLambda3Npts = false;
  bool hasAlambda3Npts = false;
  bool hasmqL3sqNpts = false;
  bool hasmtRsqNpts = false;
  bool hasAtNpts = false;
  bool hasM2Npts = false;

  // Flags to indicate required couplings are found
  bool hasyu = false;
  bool hasyc = false;
  bool hasyt = false;
  bool hasyd = false;
  bool hasys = false;
  bool hasyb = false;
  bool hasye = false;
  bool hasymu = false;
  bool hasytau = false;
  bool hasg1 = false;
  bool hasgN = false;
  bool hasg2 = false;
  bool hasg3 = false;
  bool hasv = false;
  bool hass = false;

  /*
    ----------------------------------------------------
   */

  try 
    {

      /*
	----------------------------------------------------
	Read in scan parameters from file and calculate
	all values of the parameters (more convenient
	than doing this in the for loops)
	----------------------------------------------------
      */

     
      outputCharacteristics(8);

      string line; int lineNum = 0;
	
      while (getline(cin,line)) 
	{
	  lineNum++;

	  istringstream input(line); 
	  string word1, word2;
	  input >> word1; 

	  if (word1.find("#") != 0)
	    { 
	      // Remove anything after the first comment character
	      if (word1.find("#") != string::npos)
		{
		  word1 = word1.substr(0, word1.find("#"));
		}

	      // All valid options have an = in them, so check for that
	      if (word1.find("=") == string::npos)
		{
		  cerr << "WARNING: invalid option '" << word1 << "' at line " << lineNum;
		  cerr << ": ignoring it." << endl;
		}
	      else
		{
		  // Divide into two strings: option and value
		  word2 = word1.substr(word1.find("=")+1, string::npos); //< value
		  word1 = word1.substr(0, word1.find("=")); //< option
		  word1 = ToUpper(word1);

		  istringstream kk(word2);

		  // Check against list of valid options
		  if (word1 == "TBLL")
		    {
		      kk >> tb_low;
		      hasTanbLL = true;
		    }
		  else if (word1 == "TBUL")
		    {
		      kk >> tb_up;
		      hasTanbUL = true;
		    }
		  else if (word1 == "TBNPTS")
		    {
		      kk >> tb_npts;
		      hasTanbNpts = true;
		    }
		  else if (word1 == "LAMBDALL")
		    {
		      kk >> lambda3_low;
		      hasLambda3LL = true;
		    }
		  else if (word1 == "LAMBDAUL")
		    {
		      kk >> lambda3_up;
		      hasLambda3UL = true;
		    }
		  else if (word1 == "LAMBDANPTS")
		    {
		      kk >> lambda3_npts;
		      hasLambda3Npts = true;
		    }
		  else if (word1 == "MQLSQLL")
		    {
		      kk >> mqL3sq_low;
		      hasmqL3sqLL = true;
		    }
		  else if (word1 == "MQLSQUL")
		    {
		      kk >> mqL3sq_up;
		      hasmqL3sqUL = true;
		    }
		  else if (word1 == "MQLSQNPTS")
		    {
		      kk >> mqL3sq_npts;
		      hasmqL3sqNpts = true;
		    }
		  else if (word1 == "MURSQLL")
		    {
		      kk >> mtRsq_low;
		      hasmtRsqLL = true;
		    }
		  else if (word1 == "MURSQUL")
		    {
		      kk >> mtRsq_up;
		      hasmtRsqUL = true;
		    }
		  else if (word1 == "MURSQNPTS")
		    {
		      kk >> mtRsq_npts;
		      hasmtRsqNpts = true;
		    }
		  else if (word1 == "ATLL")
		    {
		      if (hasLogAtOption)
			{
			  cerr << "WARNING: log scan of A_t already requested:";
			  cerr << " ignoring option '" << word1 << "' at line " << lineNum << "." << endl;
			}
		      else
			{
			  hasNonLogAtOption = true;
			  kk >> At_low;
			  hasAtLL = true;
			}
		    }
		  else if (word1 == "ATUL")
		    {
		      if (hasLogAtOption)
			{
			  cerr << "WARNING: log scan of A_t already requested:";
			  cerr << " ignoring option '" << word1 << "' at line " << lineNum << "." << endl;
			}
		      else
			{
			  hasNonLogAtOption = true;
			  kk >> At_up;
			  hasAtUL = true;
			}
		    }
		  else if (word1 == "ATNPTS")
		    {
		      if (hasLogAtOption)
			{
			  cerr << "WARNING: log scan of A_t already requested:";
			  cerr << " ignoring option '" << word1 << "' at line " << lineNum << "." << endl;
			}
		      else
			{
			  hasNonLogAtOption = true;
			  kk >> At_npts;
			  hasAtNpts = true;
			}
		    }
		  else if (word1 == "LOGATLL")
		    {
		      if (hasNonLogAtOption)
			{
			  cerr << "WARNING: linear scan of A_t already requested:";
			  cerr << " ignoring option '" << word1 << "' at line " << lineNum << "." << endl;
			}
		      else
			{
			  hasLogAtOption = true;
			  kk >> At_low;
			  hasAtLL = true;
			}
		    }
		  else if (word1 == "LOGATUL")
		    {
		      if (hasNonLogAtOption)
			{
			  cerr << "WARNING: linear scan of A_t already requested:";
			  cerr << " ignoring option '" << word1 << "' at line " << lineNum << "." << endl;
			}
		      else
			{
			  hasLogAtOption = true;
			  kk >> At_up;
			  hasAtUL = true;
			}
		    }
		  else if (word1 == "LOGATNPTS")
		    {
		      if (hasNonLogAtOption)
			{
			  cerr << "WARNING: linear scan of A_t already requested:";
			  cerr << " ignoring option '" << word1 << "' at line " << lineNum << "." << endl;
			}
		      else
			{
			  hasLogAtOption = true;
			  kk >> At_npts;
			  hasAtNpts = true;
			}
		    }
		  else if (word1 == "ALAMBDALL")
		    {
		      if (hasLogAlambda3Option)
			{
			  cerr << "WARNING: log scan of A_lambda3 already requested:";
			  cerr << " ignoring option '" << word1 << "' at line " << lineNum << "." << endl;
			}
		      else
			{
			  hasNonLogAlambda3Option = true;
			  kk >> Alambda3_low;
			  hasAlambda3LL = true;
			}
		    }
		  else if (word1 == "ALAMBDAUL")
		    {
		      if (hasLogAlambda3Option)
			{
			  cerr << "WARNING: log scan of A_lambda3 already requested:";
			  cerr << " ignoring option '" << word1 << "' at line " << lineNum << "." << endl;
			}
		      else
			{
			  hasNonLogAlambda3Option = true;
			  kk >> Alambda3_up;
			  hasAlambda3UL = true;
			}
		    }
		  else if (word1 == "ALAMBDANPTS")
		    {
		      if (hasLogAlambda3Option)
			{
			  cerr << "WARNING: log scan of A_lambda3 already requested:";
			  cerr << " ignoring option '" << word1 << "' at line " << lineNum << "." << endl;
			}
		      else
			{
			  hasNonLogAlambda3Option = true;
			  kk >> Alambda3_npts;
			  hasAlambda3Npts = true;
			}
		    }
		  else if (word1 == "LOGALAMBDALL")
		    {
		      if (hasNonLogAlambda3Option)
			{
			  cerr << "WARNING: linear scan of A_lambda3 already requested:";
			  cerr << " ignoring option '" << word1 << "' at line " << lineNum << "." << endl;
			}
		      else
			{
			  hasLogAlambda3Option = true;
			  kk >> Alambda3_low;
			  hasAlambda3LL = true;
			}
		    }
		  else if (word1 == "LOGALAMBDAUL")
		    {
		      if (hasNonLogAlambda3Option)
			{
			  cerr << "WARNING: linear scan of A_lambda3 already requested:";
			  cerr << " ignoring option '" << word1 << "' at line " << lineNum << "." << endl;
			}
		      else
			{
			  hasLogAlambda3Option = true;
			  kk >> Alambda3_up;
			  hasAlambda3UL = true;
			}
		    }
		  else if (word1 == "LOGALAMBDANPTS")
		    {
		      if (hasNonLogAlambda3Option)
			{
			  cerr << "WARNING: linear scan of A_lambda3 already requested:";
			  cerr << " ignoring option '" << word1 << "' at line " << lineNum << "." << endl;
			}
		      else
			{
			  hasLogAlambda3Option = true;
			  kk >> Alambda3_npts;
			  hasAlambda3Npts = true;
			}
		    }
		  else if (word1 == "M2LL")
		    {
		      kk >> M2_low;
		      hasM2LL = true;
		    }
		  else if (word1 == "M2UL")
		    {
		      kk >> M2_up;
		      hasM2UL = true;
		    }
		  else if (word1 == "M2NPTS")
		    {
		      kk >> M2_npts;
		      hasM2Npts = true;
		    }
		  else if (word1 == "YU")
		    {
		      kk >> yuin(0,0);
		      hasyu = true;
		    }
		  else if (word1 == "YC")
		    {
		      kk >> yuin(1,1);
		      hasyc = true;
		    }
		  else if (word1 == "YT")
		    {
		      kk >> yuin(2,2);
		      hasyt = true;
		    }
		  else if (word1 == "YD")
		    {
		      kk >> ydin(0,0);
		      hasyd = true;
		    }
		  else if (word1 == "YS")
		    {
		      kk >> ydin(1,1);
		      hasys = true;
		    }
		  else if (word1 == "YB")
		    {
		      kk >> ydin(2,2);
		      hasyb = true;
		    }
		  else if (word1 == "YE")
		    {
		      kk >> yein(0,0);
		      hasye = true;
		    }
		  else if (word1 == "YMU")
		    {
		      kk >> yein(1,1);
		      hasymu = true;
		    }
		  else if (word1 == "YTAU")
		    {
		      kk >> yein(2,2);
		      hasytau = true;
		    }
		  else if (word1 == "G1")
		    {
		      kk >> g1in;
		      hasg1 = true;
		    }
		  else if (word1 == "GN")
		    {
		      kk >> gNin;
		      hasgN = true;
		    }
		  else if (word1 == "G2")
		    {
		      kk >> g2in;
		      hasg2 = true;
		    }
		  else if (word1 == "G3")
		    {
		      kk >> g3in;
		      hasg3 = true;
		    }
		  else if (word1 == "HVEV")
		    {
		      kk >> vin;
		      hasv = true;
		    }
		  else if (word1 == "SVEV")
		    {
		      kk >> vsin;
		      hass = true;
		    }
		  else
		    {
		      cerr << "WARNING: unrecognised option '" << word1 << "' requested:";
		      cerr << " ignoring it." << endl;
		    }
		} //< end of option check

	    } //< end of non-comment line

	} //< end of reading file

      // Check inputs are valid

      if (!hasyu)
	{
	  cerr << "WARNING: y_u value not found: default value is " << yuin(0,0) << "." << endl;
	}

      if (!hasyc)
	{
	  cerr << "WARNING: y_c value not found: default value is " << yuin(1,1) << "." << endl;
	}

      if (!hasyt)
	{
	  cerr << "WARNING: y_t value not found: default value is " << yuin(2,2) << "." << endl;
	}

      if (!hasyd)
	{
	  cerr << "WARNING: y_d value not found: default value is " << ydin(0,0) << "." << endl;
	}

      if (!hasys)
	{
	  cerr << "WARNING: y_s value not found: default value is " << ydin(1,1) << "." << endl;
	}

      if (!hasyb)
	{
	  cerr << "WARNING: y_b value not found: default value is " << ydin(2,2) << "." << endl;
	}

      if (!hasye)
	{
	  cerr << "WARNING: y_e value not found: default value is " << yein(0,0) << "." << endl;
	}

      if (!hasymu)
	{
	  cerr << "WARNING: y_mu value not found: default value is " << yein(1,1) << "." << endl;
	}

      if (!hasytau)
	{
	  cerr << "WARNING: y_tau value not found: default value is " << yein(2,2) << "." << endl;
	}

      if (!hasg1)
	{
	  cerr << "WARNING: g_1 value not found: default value is " << g1in << "." << endl;
	}

      if (!hasgN)
	{
	  cerr << "WARNING: g_N value not found: default value is " << gNin << "." << endl;
	}

      if (!hasg2)
	{
	  cerr << "WARNING: g_2 value not found: default value is " << g2in << "." << endl;
	}

      if (!hasg3)
	{
	  cerr << "WARNING: g_3 value not found: default value is " << g3in << "." << endl;
	}

      if (!hasv)
	{
	  cerr << "WARNING: Higgs vev v value not found: default value is " << vin << "." << endl;
	}

      if (!hass)
	{
	  cerr << "WARNING: Singlet vev s value not found: default value is " << vsin << "." << endl;
	}

      if (!hasTanbLL)
	{
	  cerr << "WARNING: lower tan(beta) limit not found: using default value " << tb_low << "." << endl;
	}

      if (!hasTanbUL)
	{
	  cerr << "WARNING: upper tan(beta) limit not found: using default value " << tb_up << "." << endl;
	}

      if (!hasTanbNpts)
	{
	  cerr << "WARNING: number of tan(beta) points not found: using default value " << tb_npts << "." << endl;
	}

      if (!hasLambda3LL)
	{
	  cerr << "WARNING: lower lambda3 limit not found: using default value " << lambda3_low << "." << endl;
	}

      if (!hasLambda3UL)
	{
	  cerr << "WARNING: upper lambda3 limit not found: using default value " << lambda3_up << "." << endl;
	}

      if (!hasLambda3Npts)
	{
	  cerr << "WARNING: number of lambda3 points not found: using default value " << lambda3_npts << "." << endl;
	}

      if (!hasmqL3sqLL)
	{
	  cerr << "WARNING: lower m_Q3^2 limit not found: using default value " << mqL3sq_low << "." << endl;
	}

      if (!hasmqL3sqUL)
	{
	  cerr << "WARNING: upper m_Q3^2 limit not found: using default value " << mqL3sq_up << "." << endl;
	}

      if (!hasmqL3sqNpts)
	{
	  cerr << "WARNING: number of m_Q3^2 points not found: using default value " << mqL3sq_npts << "." << endl;
	}

      if (!hasmtRsqLL)
	{
	  cerr << "WARNING: lower m_u3^2 limit not found: using default value " << mtRsq_low << "." << endl;
	}

      if (!hasmtRsqUL)
	{
	  cerr << "WARNING: upper m_u3^2 limit not found: using default value " << mtRsq_up << "." << endl;
	}

      if (!hasmtRsqNpts)
	{
	  cerr << "WARNING: number of m_u3^2 points not found: using default value " << mtRsq_npts << "." << endl;
	}

      if (!hasLogAtOption && !hasNonLogAtOption)
	{
	  cerr << "WARNING: A_t scan type not defined: using default ";
	  if (!uselogscanAt)
	    {
	      cerr << "linear scan." << endl;
	    }
	  else
	    {
	      cerr << "log scan." << endl;
	    }
	} 
      else if (hasLogAtOption)
	{
	  uselogscanAt = true;
	}
      else
	{
	  uselogscanAt = false;
	}

      if (!hasAtLL)
	{
	  cerr << "WARNING: lower A_t limit not found: using default value ";
	  if (uselogscanAt)
	    {
	      cerr << -At_up << "." << endl;
	    }
	  else
	    {
	      cerr << At_low << "." << endl;
	    }
	}

      if (!hasAtUL)
	{
	  cerr << "WARNING: upper A_t limit not found: using default value " << At_up << "." << endl;
	}

      if (!hasAtNpts)
	{
	  cerr << "WARNING: number of A_t points not found: using default value " << At_npts << "." << endl;
	}

      if (!hasLogAlambda3Option && !hasNonLogAlambda3Option)
	{
	  cerr << "WARNING: A_lambda3 scan type not defined: using default ";
	  if (!uselogscanAlambda3)
	    {
	      cerr << "linear scan." << endl;
	    }
	  else
	    {
	      cerr << "log scan." << endl;
	    }
	} 
      else if (hasLogAlambda3Option)
	{
	  uselogscanAlambda3 = true;
	}
      else
	{
	  uselogscanAlambda3 = false;
	}

      if (!hasAlambda3LL)
	{
	  cerr << "WARNING: lower A_lambda3 limit not found: using default value ";
	  if (uselogscanAlambda3)
	    {
	      cerr << -Alambda3_up << "." << endl;
	    }
	  else
	    {
	      cerr << Alambda3_low << "." << endl;
	    }
	}

      if (!hasAlambda3UL)
	{
	  cerr << "WARNING: upper A_lambda3 limit not found: using default value " << Alambda3_up << "." << endl;
	}

      if (!hasAlambda3Npts)
	{
	  cerr << "WARNING: number of A_lambda3 points not found: using default value " << Alambda3_npts << "." << endl;
	}

      if (!hasM2LL)
	{
	  cerr << "WARNING: lower M2 limit not found: using default value " << M2_low << "." << endl;
	}

      if (!hasM2UL)
	{
	  cerr << "WARNING: upper M2 limit not found: using default value " << M2_up << "." << endl;
	}

      if (!hasM2Npts)
	{
	  cerr << "WARNING: number of M2 points not found: using default value " << M2_npts << "." << endl;
	}

      // Check ordering of upper and lower bounds
      double temp;

      if (tb_low > tb_up)
	{
	  temp = tb_low;
	  tb_low = tb_up;
	  tb_up = temp;
	}

      if (lambda3_low > lambda3_up)
	{
	  temp = lambda3_low;
	  lambda3_low = lambda3_up;
	  lambda3_up = temp;
	}

      if (Alambda3_low > Alambda3_up)
	{
	  temp = Alambda3_low;
	  Alambda3_low = Alambda3_up;
	  Alambda3_up = temp;
	}

      if (mqL3sq_low > mqL3sq_up)
	{
	  temp = mqL3sq_low;
	  mqL3sq_low = mqL3sq_up;
	  mqL3sq_up = temp;
	}

      if (mtRsq_low > mtRsq_up)
	{
	  temp = mtRsq_low;
	  mtRsq_low = mtRsq_up;
	  mtRsq_up = temp;
	}

      if (At_low > At_up)
	{
	  temp = At_low;
	  At_low = At_up;
	  At_up = temp;
	}

      if (M2_low > M2_up)
	{
	  temp = M2_low;
	  M2_low = M2_up;
	  M2_up = temp;
	}

      // Check for valid values of tan(beta)
      if (tb_low < 1.0 || tb_low > 500.0)
	{
	  cerr << "WARNING: invalid lower tan(beta) limit of " << tb_low << ":";
	  cerr << " using default 2.0." << endl;
	  tb_low = 2.0;
	}

      if (tb_up < 1.0 || tb_up > 500.0)
	{
	  cerr << "WARNING: invalid upper tan(beta) limit of " << tb_up << ":";
	  cerr << " using default 50.0." << endl;
	  tb_low = 50.0;
	}

      // Check for valid numbers of points (i.e. at least 1 for each)
      if (tb_npts < 1)
	{
	  cerr << "WARNING: invalid number of tan(beta) points " << tb_npts << " requested:";
	  cerr << " using default 1." << endl;
	  tb_npts = 1;
	}

      if (lambda3_npts < 1)
	{
	  cerr << "WARNING: invalid number of mu points " << lambda3_npts << " requested:";
	  cerr << " using default 1." << endl;
	  lambda3_npts = 1;
	}

      if (Alambda3_npts < 1)
	{
	  cerr << "WARNING: invalid number of A_lambda3 points " << Alambda3_npts << " requested:";
	  cerr << " using default 1." << endl;
	  Alambda3_npts = 1;
	}

      if (mqL3sq_npts < 1)
	{
	  cerr << "WARNING: invalid number of m_Q3^2 points " << mqL3sq_npts << " requested:";
	  cerr << " using default 1." << endl;
	  mqL3sq_npts = 1;
	}

      if (mtRsq_npts < 1)
	{
	  cerr << "WARNING: invalid number of m_u3^2 points " << mtRsq_npts << " requested:";
	  cerr << " using default 1." << endl;
	  mtRsq_npts = 1;
	}

      if (At_npts < 1)
	{
	  cerr << "WARNING: invalid number of A_t points " << At_npts << " requested:";
	  cerr << " using default 1." << endl;
	  At_npts = 1;
	}

      if (M2_npts < 1)
	{
	  cerr << "WARNING: invalid number of M_2 points " << M2_npts << " requested:";
	  cerr << " using default 1." << endl;
	  M2_npts = 1;
	}

      // Get values for all of the scanned variables

      if (tb_npts > 1) tb_vals.resize(tb_npts);
      if (lambda3_npts > 1) lambda3_vals.resize(lambda3_npts);
      if (mqL3sq_npts > 1) mqL3sq_vals.resize(mqL3sq_npts);
      if (mtRsq_npts > 1) mtRsq_vals.resize(mtRsq_npts);

      if (uselogscanAt)
	{
	  At_vals.resize(2*At_npts); //< symmetric scan about the origin
	}
      else
	{
	  if (At_npts > 1) At_vals.resize(At_npts);
	}

      if (uselogscanAlambda3)
	{
	  Alambda3_vals.resize(2*Alambda3_npts); //< symmetric scan about the origin
	}
      else
	{
	  if (Alambda3_npts > 1) Alambda3_vals.resize(Alambda3_npts);
	}

      if (M2_npts > 1) M2_vals.resize(M2_npts);

      // Work out increments
      if (tb_npts != 1)
	{
	  tb_incr = (tb_up-tb_low)/(((double)tb_npts)-1.0);
	}

      if (lambda3_npts != 1)
	{
	  lambda3_incr = (lambda3_up-lambda3_low)/(((double)lambda3_npts)-1.0);
	}

      if (mqL3sq_npts != 1)
	{
	  mqL3sq_incr = (mqL3sq_up-mqL3sq_low)/(((double)mqL3sq_npts)-1.0);
	}

      if (mtRsq_npts != 1)
	{
	  mtRsq_incr = (mtRsq_up-mtRsq_low)/(((double)mtRsq_npts)-1.0);
	}

      if (M2_npts != 1)
	{
	  M2_incr = (M2_up-M2_low)/(((double)M2_npts)-1.0);
	}

      // Get values
      tuin(0,0) = 0.0; tuin(1,1) = 0.0; tuin(2,2) = yuin(2,2) * At;
      tdin(0,0) = 0.0; tdin(1,1) = 0.0; tdin(2,2) = ydin(2,2) * Ab;
      tein(0,0) = 0.0; tein(1,1) = 0.0; tuin(2,2) = yein(2,2) * Atau;

      for (int i = 1; i <= tb_npts; i++)
	{
	  tb_vals(i-1) = tb_low + (((double)i)-1.0)*tb_incr;
	}

      for (int i = 1; i <= lambda3_npts; i++)
	{
	  lambda3_vals(i-1) = lambda3_low + (((double)i)-1.0)*lambda3_incr;
	}

      for (int i = 1; i <= mqL3sq_npts; i++)
	{
	  mqL3sq_vals(i-1) = mqL3sq_low + (((double)i)-1.0)*mqL3sq_incr;
	}

      for (int i = 1; i <= mtRsq_npts; i++)
	{
	  mtRsq_vals(i-1) = mtRsq_low + (((double)i)-1.0)*mtRsq_incr;
	}

      for (int i = 1; i <= M2_npts; i++)
	{
	  M2_vals(i-1) = M2_low + (((double)i)-1.0)*M2_incr;
	}

      // Do the same for A_t and A_lambda3, taking into account whether log scan or not
      if (!uselogscanAt)
	{
	  if (At_npts != 1)
	    {
	      At_incr = (At_up-At_low)/(((double)At_npts)-1.0);
	    }

	  for (int i = 1; i <= At_npts; i++)
	    {
	      At_vals(i-1) = At_low + (((double)i)-1.0)*At_incr;
	    }

	}
      else
	{
	  // In this case the provided values are the log of the 
	  // limit.
	  if (At_npts != 1)
	    {
	      At_incr = (At_up-At_low)/(((double)At_npts)-1.0);
	    }

	  for (int i = 1; i <= At_npts; i++)
	    {
	      At_vals((At_npts+i)-1) = exp(At_low+(((double)i)-1.0)*At_incr);
	      At_vals((At_npts+1-i)-1) = -At_vals((At_npts+i)-1);
	    }

	}


      if (!uselogscanAlambda3)
	{
	  if (Alambda3_npts != 1)
	    {
	      Alambda3_incr = (Alambda3_up-Alambda3_low)/(((double)Alambda3_npts)-1.0);
	    }

	  for (int i = 1; i <= Alambda3_npts; i++)
	    {
	      Alambda3_vals(i-1) = Alambda3_low + (((double)i)-1.0)*Alambda3_incr;
	    }

	}
      else
	{
	  // In this case the provided values are the log of the 
	  // limit.
	  if (Alambda3_npts != 1)
	    {
	      Alambda3_incr = (Alambda3_up-Alambda3_low)/(((double)Alambda3_npts)-1.0);
	    }

	  for (int i = 1; i <= Alambda3_npts; i++)
	    {
	      Alambda3_vals((Alambda3_npts+i)-1) = exp(Alambda3_low+(((double)i)-1.0)*Alambda3_incr);
	      Alambda3_vals((Alambda3_npts+1-i)-1) = -Alambda3_vals((Alambda3_npts+i)-1);
	    }

	}

      if (ENABLE_DEBUG)
	{
	  cout << "# Scan summary: " << endl;
	  cout << "# Lower tan(beta) = " << tb_vals(0) << ",";
	  cout << " Upper tan(beta) = " << tb_vals(tb_vals.size()-1) << ",";
	  cout << " Number of tan(beta) points = " << tb_vals.size() << endl;
	  cout << "# Lower lambda_3 = " << lambda3_vals(0) << ",";
	  cout << " Upper lambda_3 = " << lambda3_vals(lambda3_vals.size()-1) << ",";
	  cout << " Number of lambda_3 points = " << lambda3_vals.size() << endl;
	  cout << "# Lower A_lambda_3 = " << Alambda3_vals(0) << ",";
	  cout << " Upper A_lambda_3 = " << Alambda3_vals(Alambda3_vals.size()-1) << ",";
	  cout << " Number of A_lambda_3 points = " << Alambda3_vals.size() << endl;
	  cout << "# Lower m_qL3^2 = " << mqL3sq_vals(0) << ",";
	  cout << " Upper m_qL3^2 = " << mqL3sq_vals(mqL3sq_vals.size()-1) << ",";
	  cout << " Number of m_qL3^2 points = " << mqL3sq_vals.size() << endl;
	  cout << "# Lower m_tR^2 = " << mtRsq_vals(0) << ",";
	  cout << " Upper m_tR^2 = " << mtRsq_vals(mtRsq_vals.size()-1) << ",";
	  cout << " Number of m_tR^2 points = " << mtRsq_vals.size() << endl;
	  cout << "# Lower A_t = " << At_vals(0) << ",";
	  cout << " Upper A_t = " << At_vals(At_vals.size()-1) << ",";
	  cout << " Number of A_t points = " << At_vals.size() << endl;
	  cout << "# Lower M_2 = " << M2_vals(0) << ",";
	  cout << " Upper M_2 = " << M2_vals(M2_vals.size()-1) << ",";
	  cout << " Number of M_2 points = " << M2_vals.size() << endl;
	}

      /*
	----------------------------------------------------
      */

      /*
	----------------------------------------------------
	Finally actually do the scan DH:: TODO here onwards
	----------------------------------------------------
      */
      
      for (int i = 0; i < tb_vals.size(); i++)
	{
      for (int j = 0; j < lambda3_vals.size(); j++)
	{
      for (int k = 0; k < Alambda3_vals.size(); k++)
	{
      for (int l = 0; l < mqL3sq_vals.size(); l++)
	{
      for (int p = 0; p < mtRsq_vals.size(); p++)
	{
      for (int q = 0; q < At_vals.size(); q++)
	{
      for (int s = 0; s < M2_vals.size(); s++)
	{

	  // Update variable values
	  mQlSq(2,2) = mqL3sq_vals(l);
	  mUrSq(2,2) = mtRsq_vals(p);
	  tuin(2,2) = yuin(2,2)*At_vals(q);
	  tlambda3in = lambda3_vals(j)*Alambda3_vals(k);
	  // Default guesses for values determined by iteration
	  mHdSq = 1.0e6;
	  mHuSq = 1.0e6;
	  mSSq = 1.0e6;
	  MS = 2.0e3;

	  if (ENABLE_DEBUG)
	    {
	      cout << "# SCAN: tan(beta) = " << tb_vals(i) << ", lambda_3 = " << lambda3_vals(j)
		   << ", A_lambda_3 = " << Alambda3_vals(k) << endl;
	      cout << "# SCAN: m_qL3^2 = " << mQlSq(2,2) << ", m_tR^2 = " << mUrSq(2,2)
		   << ", A_t = " << At_vals(q) << ", M_2 = " << M2_vals(s) << endl;
	      cout << "# SCAN: v1 = " << vin/Sqrt(1.0+tb_vals(i)*tb_vals(i))
		   << ", v2 = " << vin*tb_vals(i) << ", s = " << vsin << endl;
	    }

	  // Define model, first reset problem flags
	  hasEWSBProblem = false;
	  squarksTachyons = false;
	  vectorBosonsTachyons = false;
	  higgsTachyons = false;
	  tadpoleProblem = false;
	  poleHiggsTachyons = false;
	  inaccurateHiggsMass = false;
	  hasSeriousProblem = false;

	  struct timeval tvWall, tv2Wall;
	  
	  double cpuStart = ((double) clock())/(CLOCKS_PER_SEC);

	  gettimeofday(&tvWall, NULL);


	  // m is returned at M_{SUSY}
	  m = doSimplifiedSpectrum(yuin, ydin, yein, g1in, g2in, g3in, gNin, lambda3_vals(j), lambda12in,
				   kappain, mupr, tb_vals(i), vin, vsin,
				   M1, M2_vals(s), M3, M1p, tuin, tdin, tein, 
				   tlambda3in, tlambda12in, tkappain, mQlSq, mUrSq, mDrSq, 
				   mLlSq, mErSq, mDxSq, mDxbarSq, mHdISq, mHuISq, mSISq,
				   Bmupr, mHpSq, mHpbarSq, MX, LOOPS, THRESH, 
				   TOLEWSB, mHdSq, mHuSq, mSSq, MS, 
				   hasEWSBProblem, squarksTachyons, vectorBosonsTachyons, higgsTachyons, 
				   tadpoleProblem, poleHiggsTachyons, inaccurateHiggsMass, hasSeriousProblem);


	  lambda3AtMsusy = m.get_Lambdax(); kappa3AtMsusy = m.get_Kappa(2,2);

	  f1 = ESSM_EWSBCondition1(m);
	  f2 = ESSM_EWSBCondition2(m);
	  f3 = ESSM_EWSBCondition3(m);

	  // While at M_{SUSY}, check that we are not in a point where the
	  // unphysical vacuum v_1 = v_2 = 0 is stable.
	  hasWrongVacuum = false;

	  double v1AtMsusy = m.get_vd();
	  if (v1AtMsusy < 1.0e-100)
	    {
	      hasSeriousProblem = true;
	    }
	  double v2AtMsusy = m.get_vu();

	  double tbAtMsusy = v2AtMsusy/v1AtMsusy;

	  double wrongVacuumTest = (tbAtMsusy/(1.0-sqr(tbAtMsusy)))*
	    (m.get_mHu2()*tbAtMsusy-m.get_mHd2()/tbAtMsusy)-0.5*sqr(m.get_MVZ());

	  if (wrongVacuumTest < 0.0)
	    {
	      hasWrongVacuum = true;
	    }

	  // Check for UFB or CCB problems. Because I'm not sure where
	  // exactly is the best scale to do this, do it at both and
	  // if it fails at either one flag it (we can always go back
	  // and check the points later in more detail, anyway).

	  hasUFBProblem = false; hasCCBProblem = false;

	  // Test at M_{SUSY}

	  // According to hep-ph/9507294v3, to avoid UFB problems
	  // the constraint m_1^2+m_2^2 >= 2|m_3^2| should hold at all
	  // scales greater than M_{SUSY}, and supposedly in particular at MX.
	  double ufbTest = m.get_mHd2() + m.get_mHu2() +2.0*sqr(m.get_Lambdax()*m.get_vs()/Sqrt(2.0))
	    -2.0*fabs(m.get_Lambdax()*m.get_TLambdax()*m.get_vs()/Sqrt(2.0));

	  if (ufbTest < 0)
	    {
	      hasUFBProblem = true;
	    }

	  // For simplicity, apply the CCB criterion of hep-ph/9604417v3 Eq. (2a)
	  // and Eq. (5) of hep-ph/9507294v3 at M_{SUSY} (this might not be
	  // the optimum scale though).

	  double ccbTest = 3.0*(sqr(m.get_Lambdax()*m.get_vs()/Sqrt(2.0))+m.get_mHu2()
	  			+m.get_mq2( 2, 2)+m.get_mu2(2, 2))
	    -sqr(get_softAu(m, 2, 2));

	  if (ccbTest < 0.0)
	    {
	      hasCCBProblem = true;
	    }

	  // Check the other CCB conditions as well, even though the chance
	  // of them being a problem is small.
	  ccbTest = 3.0*(sqr(m.get_Lambdax()*m.get_vs()/Sqrt(2.0))+m.get_mHu2()
	  		 +m.get_mq2(1, 1)+m.get_mu2(1, 1))
	    -sqr(get_softAu(m, 1, 1));

	  if (ccbTest < 0.0)
	    {
	      hasCCBProblem = true;
	    }

	  ccbTest = 3.0*(sqr(m.get_Lambdax()*m.get_vs()/Sqrt(2.0))+m.get_mHu2()
	  		 +m.get_mq2(0, 0)+m.get_mu2(0, 0))
	    -sqr(get_softAu(m, 0, 0));

	  if (ccbTest < 0.0)
	    {
	      hasCCBProblem = true;
	    }

	  // A_b CCB tests
	  ccbTest = 3.0*(sqr(m.get_Lambdax()*m.get_vs()/Sqrt(2.0))+m.get_mHd2()
	  			+m.get_mq2( 2, 2)+m.get_md2(2, 2))
	    -sqr(get_softAd(m, 2, 2));

	  if (ccbTest < 0.0)
	    {
	      hasCCBProblem = true;
	    }

	  ccbTest = 3.0*(sqr(m.get_Lambdax()*m.get_vs()/Sqrt(2.0))+m.get_mHd2()
	  		 +m.get_mq2(1, 1)+m.get_md2(1, 1))
	    -sqr(get_softAd(m, 1, 1));

	  if (ccbTest < 0.0)
	    {
	      hasCCBProblem = true;
	    }

	  ccbTest = 3.0*(sqr(m.get_Lambdax()*m.get_vs()/Sqrt(2.0))+m.get_mHd2()
	  		 +m.get_mq2(0, 0)+m.get_md2(0, 0))
	    -sqr(get_softAd(m, 0, 0));

	  if (ccbTest < 0.0)
	    {
	      hasCCBProblem = true;
	    }

	  // A_tau test
	  ccbTest = 3.0*(sqr(m.get_Lambdax()*m.get_vs()/Sqrt(2.0))+m.get_mHd2()
	  			+m.get_ml2( 2, 2)+m.get_me2(2, 2))
	    -sqr(get_softAe(m, 2, 2));

	  if (ccbTest < 0.0)
	    {
	      hasCCBProblem = true;
	    }

	  ccbTest = 3.0*(sqr(m.get_Lambdax()*m.get_vs()/Sqrt(2.0))+m.get_mHd2()
	  		 +m.get_ml2(1, 1)+m.get_me2(1, 1))
	    -sqr(get_softAe(m, 1, 1));

	  if (ccbTest < 0.0)
	    {
	      hasCCBProblem = true;
	    }

	  ccbTest = 3.0*(sqr(m.get_Lambdax()*m.get_vs()/Sqrt(2.0))+m.get_mHd2()
	  		 +m.get_ml2(0, 0)+m.get_me2(0, 0))
	    -sqr(get_softAe(m, 0, 0));

	  if (ccbTest < 0.0)
	    {
	      hasCCBProblem = true;
	    }

	  m.run_to(MX, PRECISION);

	  // Test at MX
	  ccbTest = 3.0*(sqr(m.get_Lambdax()*m.get_vs()/Sqrt(2.0))+m.get_mHu2()
	  			+m.get_mq2( 2, 2)+m.get_mu2(2, 2))
	    -sqr(get_softAu(m, 2, 2));

	  if (ccbTest < 0.0)
	    {
	      hasCCBProblem = true;
	    }

	  // Check the other CCB conditions as well, even though the chance
	  // of them being a problem is small.
	  ccbTest = 3.0*(sqr(m.get_Lambdax()*m.get_vs()/Sqrt(2.0))+m.get_mHu2()
	  		 +m.get_mq2(1, 1)+m.get_mu2(1, 1))
	    -sqr(get_softAu(m, 1, 1));

	  if (ccbTest < 0.0)
	    {
	      hasCCBProblem = true;
	    }

	  ccbTest = 3.0*(sqr(m.get_Lambdax()*m.get_vs()/Sqrt(2.0))+m.get_mHu2()
	  		 +m.get_mq2(0, 0)+m.get_mu2(0, 0))
	    -sqr(get_softAu(m, 0, 0));

	  if (ccbTest < 0.0)
	    {
	      hasCCBProblem = true;
	    }

	  // A_b CCB tests
	  ccbTest = 3.0*(sqr(m.get_Lambdax()*m.get_vs()/Sqrt(2.0))+m.get_mHd2()
	  			+m.get_mq2( 2, 2)+m.get_md2(2, 2))
	    -sqr(get_softAd(m, 2, 2));

	  if (ccbTest < 0.0)
	    {
	      hasCCBProblem = true;
	    }

	  ccbTest = 3.0*(sqr(m.get_Lambdax()*m.get_vs()/Sqrt(2.0))+m.get_mHd2()
	  		 +m.get_mq2(1, 1)+m.get_md2(1, 1))
	    -sqr(get_softAd(m, 1, 1));

	  if (ccbTest < 0.0)
	    {
	      hasCCBProblem = true;
	    }

	  ccbTest = 3.0*(sqr(m.get_Lambdax()*m.get_vs()/Sqrt(2.0))+m.get_mHd2()
	  		 +m.get_mq2(0, 0)+m.get_md2(0, 0))
	    -sqr(get_softAd(m, 0, 0));

	  if (ccbTest < 0.0)
	    {
	      hasCCBProblem = true;
	    }

	  // A_tau test
	  ccbTest = 3.0*(sqr(m.get_Lambdax()*m.get_vs()/Sqrt(2.0))+m.get_mHd2()
	  			+m.get_ml2( 2, 2)+m.get_me2(2, 2))
	    -sqr(get_softAe(m, 2, 2));

	  if (ccbTest < 0.0)
	    {
	      hasCCBProblem = true;
	    }

	  ccbTest = 3.0*(sqr(m.get_Lambdax()*m.get_vs()/Sqrt(2.0))+m.get_mHd2()
	  		 +m.get_ml2(1, 1)+m.get_me2(1, 1))
	    -sqr(get_softAe(m, 1, 1));

	  if (ccbTest < 0.0)
	    {
	      hasCCBProblem = true;
	    }

	  ccbTest = 3.0*(sqr(m.get_Lambdax()*m.get_vs()/Sqrt(2.0))+m.get_mHd2()
	  		 +m.get_ml2(0, 0)+m.get_me2(0, 0))
	    -sqr(get_softAe(m, 0, 0));

	  if (ccbTest < 0.0)
	    {
	      hasCCBProblem = true;
	    }

	  
	  // Calculate fine tuning
	  hasTuningError = false;

	  bool tuningEWSBProblem = false;

	  // try
	  //   {
	  //     Eigen::ArrayXd pars;
	  //     pars.resize(tuning_parameters::NUMESSMTUNINGPARS);
	  //     pars(tuning_parameters::lam3) = lambda3_vals(j);
	  //     pars(tuning_parameters::Alam3) = Alambda3_vals(k);
	  //     pars(tuning_parameters::M1) = M1;
	  //     pars(tuning_parameters::M2) = M2_vals(s);
	  //     pars(tuning_parameters::M3) = M3;
	  //     pars(tuning_parameters::M1p) = M1p;
	  //     pars(tuning_parameters::Au3) = At_vals(q);
	  //     pars(tuning_parameters::mH13Sq) = m.get_mHd2();
	  //     pars(tuning_parameters::mH23Sq) = m.get_mHu2();
	  //     pars(tuning_parameters::mS3Sq) = m.get_ms2();
	  //     pars(tuning_parameters::mqL3Sq) = mQlSq(2,2);
	  //     pars(tuning_parameters::mtRSq) = mUrSq(2,2);
  	  //     //	      tunings = doCalcpMSSMFineTuning(m, MS, tuningEWSBProblem, hasTuningError, USEAPPROXSOLNS, TOLEWSB);
	  //     tunings = doCalcESSMTuningNumerically(m, MS, MX, pars, pE6SSMftBCs);
	  //     // Possible numerical errors here in running should be accounted for.
	  //     // if (!hasEWSBProblem && tuningEWSBProblem) //< tuning calculation hasn't found solution satisfies tolerance
	  //     // 	{
	  //     // 	  // Probably a numerical issue.
	  //     // 	  cerr << "WARNING: numerical issue in RG running for tuning calculation." << endl;
	  //     // 	  hasTuningError = true;
	  //     // 	}
	  //     cout << tunings << endl;
	  //   }
	  // catch(const string & a) 
	  //   { 
	  //     cerr << "WARNING: serious numerical problem encountered in fine tuning calculation." << endl; 
	  //     hasSeriousProblem = true; 
	  //   }
	  // catch(const char * a) 
	  //   { 
	  //     cerr << "WARNING: serious numerical problem encountered in fine tuning calculation." << endl; 
	  //     hasSeriousProblem = true; 
	  //   }
	  // catch(...) 
	  //   { 
	  //     cerr << "WARNING: unknown problem encountered in fine tuning calculation." << endl; 
	  //     hasSeriousProblem = true; 
	  //   }

	  double cpuEnd = ((double) clock())/(CLOCKS_PER_SEC);
	  gettimeofday(&tv2Wall, NULL);
	  int wall_time = (tv2Wall.tv_sec-tvWall.tv_sec)*1000000+((int)tv2Wall.tv_usec-(int)tvWall.tv_usec);


	  if (ENABLE_DEBUG)
	    {
	      cout << "# After spectrum calculation:" << endl;
	      cout << "# f1 = " << f1 << endl;
	      cout << "# f2 = " << f2 << endl;
	      cout << "# f3 = " << f3 << endl;
	      cout << "# m_Hd^2 = " << mHdSq << endl;
	      cout << "# m_Hu^2 = " << mHuSq << endl;
	      cout << "# m_s^2 = " << mSSq << endl;
	      cout << "# M_SUSY = " << MS << endl;
	      cout << "# Wall time = " << wall_time << endl;
	      cout << "# CPU time = " << (int) ((cpuEnd-cpuStart)*1000000) << endl;

	    }

	  if (!hasSeriousProblem)
	    {
	      // Display output
	      cout << tb_vals(i) << " ";
	      cout << lambda3_vals(j) << " ";
	      cout << Alambda3_vals(k) << " ";
	      cout << mqL3sq_vals(l) << " ";
	      cout << mtRsq_vals(p) << " ";
	      cout << At_vals(q) << " ";
	      cout << M2_vals(s) << " "; 
	      cout << lambda3AtMsusy << " ";
	      cout << kappa3AtMsusy << " ";
	      cout << m.get_physical().Mhh(0) << " ";
	      cout << m.get_physical().MAh(2) << " ";
	      cout << mHdSq << " ";
	      cout << mHuSq << " ";
	      cout << mSSq << " ";
	      cout << f1 << " ";
	      cout << f2 << " ";
	      cout << f3 << " ";
	      cout << MS << " ";
	      cout << m.get_MStop()(0) << " ";
	      cout << m.get_MStop()(1) << " ";
	      cout << m.get_MSbot()(0) << " ";
	      cout << m.get_MSbot()(1) << " ";
	      cout << m.get_MCha()(0) << " ";
	      cout << m.get_MCha()(1) << " ";
	      cout << m.get_MChi()(0)<< " ";
	      cout << m.get_MChi()(1)<< " ";
	      cout << m.get_MChi()(2)<< " ";
	      cout << m.get_MChi()(3)<< " ";
	      cout << m.get_MChi()(4)<< " ";
	      cout << m.get_MChi()(5)<< " ";
	      cout << m.get_Mhh()(0) << " ";
	      cout << m.get_Mhh()(1) << " ";
	      cout << m.get_Mhh()(2) << " ";
	      cout << m.get_MAh()(0) << " ";
	      cout << m.get_MAh()(1) << " ";
	      cout << m.get_MAh()(2) << " ";

	      // cout << tunings.max() << " ";
	      // for (int i = 1; i <= tunings.displayEnd(); i++)
	      // 	{
	      // 	  cout << tunings(i) << " ";
	      // 	}
	      if (hasUFBProblem)
		{
		  cout << UFBPROBLEM << " ";
		}
	      else
		{
		  cout << "0 ";
		}
	      if (hasCCBProblem)
		{
		  cout << CCBPROBLEM << " ";
		}
	      else
		{
		  cout << "0 ";
		}
	      
	      if (hasEWSBProblem)
		{
		  cout << EWSBPROBLEM << " ";
		}
	      else
		{
		  cout << "0 ";
		}
	      
	      if (hasWrongVacuum)
		{
		  cout << WRONGVACUUM << " ";
		}
	      else
		{
		  cout << "0 ";
		}
	      
	      if (squarksTachyons)
		{
		  cout << SQUARKTACHYON << " ";
		}
	      else
		{
		  cout << "0 ";
		}
	      if (vectorBosonsTachyons)
		{
		  cout << VECTORBOSONTACHYON << " ";
		}
	      else
		{
		  cout << "0 ";
		}
	      if (higgsTachyons)
		{
		  cout << HIGGSTACHYON << " ";
		}
	      else
		{
		  cout << "0 ";
		}
	      if (tadpoleProblem)
		{
		  cout << TADPOLESPROBLEM << " ";
		}
	      else
		{
		  cout << "0 ";
		}
	      if (m.get_physical().Mhh(0) < HIGGSCENT-HIGGSERROR || m.get_physical().Mhh(0) > HIGGSCENT+HIGGSERROR)
		{
		  cout << NOTEXPVALID << " ";
		}
	      else
		{
		  cout << "0 ";
		}
	      if (poleHiggsTachyons)
		{
		  cout << POLEHIGGSTACHYON << " ";
		}
	      else
		{
		  cout << "0 ";
		}
	      if (inaccurateHiggsMass)
		{
		  cout << HIGGSPROBLEM << " ";
		}
	      else
		{
		  cout << "0 ";
		}
	      
	      if (hasTuningError)
		{
		  cout << TUNINGERROR << " ";
		}
	      else
		{
		  cout << "0 ";
		}
	      
	      if (hasSeriousProblem)
		{
		  cout << NUMERICALPROBLEM << endl;
		}
	      else
		{
		  cout << "0" << endl;
		}
	    }
	} //< M_2 scan
	}//< A_t scan
	} //< m_u3^2 scan 
	} //< m_Q3^2 scan
	} //< Alambda scan
	} //< lambda scan
	} //< tan(beta) scan
      
      /*
	----------------------------------------------------
      */
      
      
      /*
	----------------------------------------------------
	End of calculations
	----------------------------------------------------
      */

    } 
  catch(const string & a) { cout << a; return -1; }
  catch(const char * a) { cout << a; return -1; }
  catch (const Error& error) 
    {
      ERROR(error.what());
      return EXIT_FAILURE;
     }
  catch(...) { cout << "Unknown type of exception caught.\n"; return -1; }

  return 0;
}

// Prints an error message
void errorCall()
{
  ostringstream ii;
  ii << "essmScanner.x called with incorrect arguments.\n";
  ii << "To run, you must provide as a file containing\n";
  ii << "the scan parameters you wish to use. Usage is:\n";
  ii << "./essmScanner.x < input_file.params\n";
  throw ii.str();
}

// Returns a string with all characters in upper case: very handy
string ToUpper(const string & s) {
        string result;
        unsigned int index;
        for (index = 0; index < s.length(); index++) {
	  char a = s[index];
	  a = toupper(a);
	  result = result + a;
        }
	
	return result;
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
								     mSSq, mHdISq, mHdISq, mSISq, mDxSq, mDxbarSq, mHpSq, 
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
