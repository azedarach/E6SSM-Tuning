/*
  Scan over the E6SSM parameter space and calculate the fine tuning 
  at each point.
  Modified to use the FlexibleSUSY E6SSM spectrum generator instead of
  SOFTSUSY routines.
 */

#include "./E6SSM_Spectrum_Generators/models/genericE6SSM/genericE6SSM_input_parameters.hpp"
#include "./E6SSM_Spectrum_Generators/models/genericE6SSM/genericE6SSM_spectrum_generator.hpp"

#include "./E6SSM_Spectrum_Generators/src/error.hpp"
#include "./E6SSM_Spectrum_Generators/src/spectrum_generator_settings.hpp"
#include "./E6SSM_Spectrum_Generators/src/lowe.h"
#include "./E6SSM_Spectrum_Generators/src/command_line_options.hpp"
#include "./E6SSM_Spectrum_Generators/src/wrappers.hpp"

#include <iostream>
#include <cstdlib>

#include "flags.h"
#include "tuningnumerics.h"
#include "tuningutils.h"
#include "essmtuningutils.h"

void errorCall();
string ToUpper(const string & s);

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

  double mupr = 5000.0; // GeV DH:: check
  double Bmupr = 5000.0; // GeV^2 DH:: check

  double mH11Sq = Sqr(5000.0); // GeV^2 DH:: check
  double mH12Sq = Sqr(5000.0); // GeV^2 DH:: check

  double mH21Sq = Sqr(5000.0); // GeV^2 DH:: check
  double mH22Sq = Sqr(5000.0); // GeV^2 DH:: check

  double mS1Sq = Sqr(5000.0); // GeV^2 DH:: check
  double mS2Sq = Sqr(5000.0); // GeV^2 DH:: check

  double mHpSq = Sqr(5000.0); // GeV^2 DH:: check
  double mHpbarSq = Sqr(5000.0); // GeV^2 DH:: check

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
  DoubleMatrix yuin(3,3), ydin(3,3), yein(3,3); //< Yukawas
  DoubleVector gin(3); //< gauge couplings
  double gNin = 0.0;
  double vin = 246.0; //< Higgs vev, GeV

  // Need M_Z' at least 2.5 TeV
  double vsin = 5000.0; //< singlet vev, GeV

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

  DoubleVector mGaugino(4);
  mGaugino(1) = M1;
  mGaugino(2) = M2;
  mGaugino(3) = M3;
  mGaugino(4) = M1p;

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
  Eigen::Matrix<double,3,3> tu, td, te, tkappa;
  Eigen::Matrix<double,2,2> tlambdaI;

  // Exotic soft scalar masses
  Eigen::Matrix<double,3,3> mXSq, mXbarSq;

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

  // Gravitino mass
  double mGrav = 0.0;

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
  DoubleVector tb_vals(1), lambda3_vals(1), Alambda3_vals(1), mqL3sq_vals(1), mtRsq_vals(1), At_vals(1), M2_vals(1);

  // Variables to store values of couplings at M_{SUSY}
  double lambda3AtMsusy, kappa3AtMsusy;

  // Vector for storing fine tunings
  DoubleVector tunings(NUMPMSSMPARS);

  // Variables to store values of EWSB conditions
  double f1, f2, f3;

  // Models that we will calculate the fine tuning for
  Spectrum_generator_settings spectrum_generator_settings;
  genericE6SSM_spectrum_generator<algorithm_type> spectrum_generator;

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
		  else if (word1 == "LOGLAMBDALL")
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
		      kk >> yuin(1,1);
		      hasyu = true;
		    }
		  else if (word1 == "YC")
		    {
		      kk >> yuin(2,2);
		      hasyc = true;
		    }
		  else if (word1 == "YT")
		    {
		      kk >> yuin(3,3);
		      hasyt = true;
		    }
		  else if (word1 == "YD")
		    {
		      kk >> ydin(1,1);
		      hasyd = true;
		    }
		  else if (word1 == "YS")
		    {
		      kk >> ydin(2,2);
		      hasys = true;
		    }
		  else if (word1 == "YB")
		    {
		      kk >> ydin(3,3);
		      hasyb = true;
		    }
		  else if (word1 == "YE")
		    {
		      kk >> yein(1,1);
		      hasye = true;
		    }
		  else if (word1 == "YMU")
		    {
		      kk >> yein(2,2);
		      hasymu = true;
		    }
		  else if (word1 == "YTAU")
		    {
		      kk >> yein(3,3);
		      hasytau = true;
		    }
		  else if (word1 == "G1")
		    {
		      kk >> gin(1);
		      hasg1 = true;
		    }
		  else if (word1 == "GN")
		    {
		      kk >> gNin;
		      hasgN = true;
		    }
		  else if (word1 == "G2")
		    {
		      kk >> gin(2);
		      hasg2 = true;
		    }
		  else if (word1 == "G3")
		    {
		      kk >> gin(3);
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
	  cerr << "WARNING: y_u value not found: default value is " << yuin(1,1) << "." << endl;
	}

      if (!hasyc)
	{
	  cerr << "WARNING: y_c value not found: default value is " << yuin(2,2) << "." << endl;
	}

      if (!hasyt)
	{
	  cerr << "WARNING: y_t value not found: default value is " << yuin(3,3) << "." << endl;
	}

      if (!hasyd)
	{
	  cerr << "WARNING: y_d value not found: default value is " << ydin(1,1) << "." << endl;
	}

      if (!hasys)
	{
	  cerr << "WARNING: y_s value not found: default value is " << ydin(2,2) << "." << endl;
	}

      if (!hasyb)
	{
	  cerr << "WARNING: y_b value not found: default value is " << ydin(3,3) << "." << endl;
	}

      if (!hasye)
	{
	  cerr << "WARNING: y_e value not found: default value is " << yein(1,1) << "." << endl;
	}

      if (!hasymu)
	{
	  cerr << "WARNING: y_mu value not found: default value is " << yein(2,2) << "." << endl;
	}

      if (!hasytau)
	{
	  cerr << "WARNING: y_tau value not found: default value is " << yein(3,3) << "." << endl;
	}

      if (!hasg1)
	{
	  cerr << "WARNING: g_1 value not found: default value is " << gin(1) << "." << endl;
	}

      if (!hasgN)
	{
	  cerr << "WARNING: g_N value not found: default value is " << gNin << "." << endl;
	}

      if (!hasg2)
	{
	  cerr << "WARNING: g_2 value not found: default value is " << gin(2) << "." << endl;
	}

      if (!hasg3)
	{
	  cerr << "WARNING: g_3 value not found: default value is " << gin(3) << "." << endl;
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

      if (tb_npts > 1) tb_vals.setEnd(tb_npts);
      if (lambda3_npts > 1) lambda3_vals.setEnd(lambda3_npts);
      if (mqL3sq_npts > 1) mqL3sq_vals.setEnd(mqL3sq_npts);
      if (mtRsq_npts > 1) mtRsq_vals.setEnd(mtRsq_npts);

      if (uselogscanAt)
	{
	  At_vals.setEnd(2*At_npts); //< symmetric scan about the origin
	}
      else
	{
	  if (At_npts > 1) At_vals.setEnd(At_npts);
	}

      if (uselogscanAlambda3)
	{
	  Alambda3_vals.setEnd(2*Alambda3_npts); //< symmetric scan about the origin
	}
      else
	{
	  if (Alambda3_npts > 1) Alambda3_vals.setEnd(Alambda3_npts);
	}

      if (M2_npts > 1) M2_vals.setEnd(M2_npts);

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
      for (int i = 1; i <= tb_npts; i++)
	{
	  tb_vals(i) = tb_low + (((double)i)-1.0)*tb_incr;
	}

      for (int i = 1; i <= lambda3_npts; i++)
	{
	  lambda3_vals(i) = lambda3_low + (((double)i)-1.0)*lambda3_incr;
	}

      for (int i = 1; i <= mqL3sq_npts; i++)
	{
	  mqL3sq_vals(i) = mqL3sq_low + (((double)i)-1.0)*mqL3sq_incr;
	}

      for (int i = 1; i <= mtRsq_npts; i++)
	{
	  mtRsq_vals(i) = mtRsq_low + (((double)i)-1.0)*mtRsq_incr;
	}

      for (int i = 1; i <= M2_npts; i++)
	{
	  M2_vals(i) = M2_low + (((double)i)-1.0)*M2_incr;
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
	      At_vals(i) = At_low + (((double)i)-1.0)*At_incr;
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
	      At_vals(At_npts+i) = exp(At_low+(((double)i)-1.0)*At_incr);
	      At_vals(At_npts+1-i) = -At_vals(At_npts+i);
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
	      Alambda3_vals(i) = Alambda3_low + (((double)i)-1.0)*Alambda3_incr;
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
	      Alambda3_vals(Alambda3_npts+i) = exp(Alambda3_low+(((double)i)-1.0)*Alambda3_incr);
	      Alambda3_vals(Alambda3_npts+1-i) = -Alambda3_vals(Alambda3_npts+i);
	    }

	}

      /*
	----------------------------------------------------
      */

      /*
	----------------------------------------------------
	Finally actually do the scan DH:: TODO here onwards
	----------------------------------------------------
      */
      
      for (int i = 1; i <= tb_vals.displayEnd(); i++)
	{
      for (int j = 1; j <= lambda3_vals.displayEnd(); j++)
	{
      for (int k = 1; k <= Alambda3_vals.displayEnd(); k++)
	{
      for (int l = 1; l <= mqL3sq_vals.displayEnd(); l++)
	{
      for (int p = 1; p <= mtRsq_vals.displayEnd(); p++)
	{
      for (int q = 1; q <= At_vals.displayEnd(); q++)
	{
      for (int s = 1; s <= M2_vals.displayEnd(); s++)
	{

	  // // Update variable values
	  // mQlSq(3,3) = mqL3sq_vals(l);
	  // mUrSq(3,3) = mtRsq_vals(p);
	  // tu(3,3) = yuin(3,3)*At_vals(q);
	  // mGaugino(2) = M2_vals(s);

	  // // Default guesses for values determined by iteration
	  // mHdSq = 1.0e6;
	  // mHuSq = 1.0e6;
	  // mSSq = 1.0e6;
	  // MS = 1.0e3;

	  // // Define model, first reset problem flags
	  // hasEWSBProblem = false;
	  // squarksTachyons = false;
	  // higgsTachyons = false;
	  // tadpoleProblem = false;
	  // poleHiggsTachyons = false;
	  // inaccurateHiggsMass = false;
	  // hasSeriousProblem = false;

	  // // m is returned at M_{SUSY}
	  // m = doSimplifiedSpectrum(yuin, ydin, yein, gin, gNin, lambda3_vals(j), tb_vals(i), vin, sin,
	  // 			   mGaugino, tu, td, te, mQlSq, mUrSq, mDrSq, mLlSq, mErSq, 
	  // 			   B_vals(k)*mu_vals(j), mGrav, MX, LOOPS, THRESH, oneset, 
	  // 			   physpars, TOLEWSB, mHdSq, mHuSq, mSSq, MS, 
	  // 			   hasEWSBProblem, squarksTachyons, higgsTachyons, 
	  // 			   tadpoleProblem, poleHiggsTachyons, 
	  // 			   inaccurateHiggsMass, hasSeriousProblem);

	  // f1 = ESSM_EWSBCondition1(m);
	  // f2 = ESSM_EWSBCondition2(m);
	  // f3 = ESSM_EWSBCondition3(m);

	  // // While at M_{SUSY}, check that we are not in a point where the
	  // // unphysical vacuum v_1 = v_2 = 0 is stable.
	  // hasWrongVacuum = false;

	  // double wrongVacuumTest = (m.displayTanb()/(1.0-sqr(m.displayTanb())))*
	  //   (m.displayMh2Squared()*m.displayTanb()-m.displayMh1Squared()/m.displayTanb())-0.5*sqr(m.displayMzRun());

	  // if (wrongVacuumTest < 0.0)
	  //   {
	  //     hasWrongVacuum = true;
	  //   }

	  // // Check for UFB or CCB problems. Because I'm not sure where
	  // // exactly is the best scale to do this, do it at both and
	  // // if it fails at either one flag it (we can always go back
	  // // and check the points later in more detail, anyway).

	  // hasUFBProblem = false; hasCCBProblem = false;

	  // // Test at M_{SUSY}

	  // // According to hep-ph/9507294v3, to avoid UFB problems
	  // // the constrain m_1^2+m_2^2 >= 2|m_3^2| should hold at all
	  // // scales greater than M_{SUSY}, and supposedly in particular at MX.
	  // double ufbTest = m.displayMh1Squared() + m.displayMh2Squared() +2.0*sqr(m.displaySusyMu())
	  //   -2.0*fabs(m.displayM3Squared());

	  // if (ufbTest < 0)
	  //   {
	  //     hasUFBProblem = true;
	  //   }

	  // // For simplicity, apply the CCB criterion of hep-ph/9604417v3 Eq. (2a)
	  // // and Eq. (5) of hep-ph/9507294v3 at M_{SUSY} (this might not be
	  // // the optimum scale though).

	  // double ccbTest = 3.0*(sqr(m.displaySusyMu())+m.displayMh2Squared()
	  // 			+m.displaySoftMassSquared(mQl, 3, 3)+m.displaySoftMassSquared(mUr, 3, 3))
	  //   -sqr(m.displaySoftA(UA, 3, 3));

	  // if (ccbTest < 0.0)
	  //   {
	  //     hasCCBProblem = true;
	  //   }

	  // // Check the other CCB conditions as well, even though the chance
	  // // of them being a problem is small.
	  // ccbTest = 3.0*(sqr(m.displaySusyMu())+m.displayMh2Squared()
	  // 		 +m.displaySoftMassSquared(mQl, 2, 2)+m.displaySoftMassSquared(mUr, 2, 2))
	  //   -sqr(m.displaySoftA(UA, 2, 2));

	  // if (ccbTest < 0.0)
	  //   {
	  //     hasCCBProblem = true;
	  //   }

	  // ccbTest = 3.0*(sqr(m.displaySusyMu())+m.displayMh2Squared()
	  // 		 +m.displaySoftMassSquared(mQl, 1, 1)+m.displaySoftMassSquared(mUr, 1, 1))
	  //   -sqr(m.displaySoftA(UA, 1, 1));

	  // if (ccbTest < 0.0)
	  //   {
	  //     hasCCBProblem = true;
	  //   }

	  // // A_b CCB tests
	  // ccbTest = 3.0*(sqr(m.displaySusyMu())+m.displayMh1Squared()
	  // 		 +m.displaySoftMassSquared(mQl, 3, 3)+m.displaySoftMassSquared(mDr, 3, 3))
	  //   -sqr(m.displaySoftA(DA, 3, 3));

	  // if (ccbTest < 0.0)
	  //   {
	  //     hasCCBProblem = true;
	  //   }

	  // ccbTest = 3.0*(sqr(m.displaySusyMu())+m.displayMh1Squared()
	  // 		 +m.displaySoftMassSquared(mQl, 2, 2)+m.displaySoftMassSquared(mDr, 2, 2))
	  //   -sqr(m.displaySoftA(DA, 2, 2));

	  // if (ccbTest < 0.0)
	  //   {
	  //     hasCCBProblem = true;
	  //   }

	  // ccbTest = 3.0*(sqr(m.displaySusyMu())+m.displayMh1Squared()
	  // 		 +m.displaySoftMassSquared(mQl, 1, 1)+m.displaySoftMassSquared(mDr, 1, 1))
	  //   -sqr(m.displaySoftA(DA, 1, 1));

	  // if (ccbTest < 0.0)
	  //   {
	  //     hasCCBProblem = true;
	  //   }

	  // // A_tau test
	  // ccbTest = 3.0*(sqr(m.displaySusyMu())+m.displayMh1Squared()
	  // 		 +m.displaySoftMassSquared(mLl, 3, 3)+m.displaySoftMassSquared(mEr, 3, 3))
	  //   -sqr(m.displaySoftA(EA, 3, 3));

	  // if (ccbTest < 0.0)
	  //   {
	  //     hasCCBProblem = true;
	  //   }

	  // ccbTest = 3.0*(sqr(m.displaySusyMu())+m.displayMh1Squared()
	  // 		 +m.displaySoftMassSquared(mLl, 2, 2)+m.displaySoftMassSquared(mEr, 2, 2))
	  //   -sqr(m.displaySoftA(EA, 2, 2));

	  // if (ccbTest < 0.0)
	  //   {
	  //     hasCCBProblem = true;
	  //   }

	  // ccbTest = 3.0*(sqr(m.displaySusyMu())+m.displayMh1Squared()
	  // 		 +m.displaySoftMassSquared(mLl, 1, 1)+m.displaySoftMassSquared(mEr, 1, 1))
	  //   -sqr(m.displaySoftA(EA, 1, 1));

	  // if (ccbTest < 0.0)
	  //   {
	  //     hasCCBProblem = true;
	  //   }

	  // m.runto(MX);

	  // // Test at MX

	  // ufbTest = m.displayMh1Squared() + m.displayMh2Squared() +2.0*sqr(m.displaySusyMu())
	  //   -2.0*fabs(m.displayM3Squared());

	  // if (ufbTest < 0)
	  //   {
	  //     hasUFBProblem = true;
	  //   }

	  // ccbTest = 3.0*(sqr(m.displaySusyMu())+m.displayMh2Squared()
	  // 		 +m.displaySoftMassSquared(mQl, 3, 3)+m.displaySoftMassSquared(mUr, 3, 3))
	  //   -sqr(m.displaySoftA(UA, 3, 3));

	  // if (ccbTest < 0.0)
	  //   {
	  //     hasCCBProblem = true;
	  //   }

	  // // Check the other CCB conditions as well, even though the chance
	  // // of them being a problem is small.
	  // ccbTest = 3.0*(sqr(m.displaySusyMu())+m.displayMh2Squared()
	  // 		 +m.displaySoftMassSquared(mQl, 2, 2)+m.displaySoftMassSquared(mUr, 2, 2))
	  //   -sqr(m.displaySoftA(UA, 2, 2));

	  // if (ccbTest < 0.0)
	  //   {
	  //     hasCCBProblem = true;
	  //   }

	  // ccbTest = 3.0*(sqr(m.displaySusyMu())+m.displayMh2Squared()
	  // 		 +m.displaySoftMassSquared(mQl, 1, 1)+m.displaySoftMassSquared(mUr, 1, 1))
	  //   -sqr(m.displaySoftA(UA, 1, 1));

	  // if (ccbTest < 0.0)
	  //   {
	  //     hasCCBProblem = true;
	  //   }

	  // // A_b CCB tests
	  // ccbTest = 3.0*(sqr(m.displaySusyMu())+m.displayMh1Squared()
	  // 		 +m.displaySoftMassSquared(mQl, 3, 3)+m.displaySoftMassSquared(mDr, 3, 3))
	  //   -sqr(m.displaySoftA(DA, 3, 3));

	  // if (ccbTest < 0.0)
	  //   {
	  //     hasCCBProblem = true;
	  //   }

	  // ccbTest = 3.0*(sqr(m.displaySusyMu())+m.displayMh1Squared()
	  // 		 +m.displaySoftMassSquared(mQl, 2, 2)+m.displaySoftMassSquared(mDr, 2, 2))
	  //   -sqr(m.displaySoftA(DA, 2, 2));

	  // if (ccbTest < 0.0)
	  //   {
	  //     hasCCBProblem = true;
	  //   }

	  // ccbTest = 3.0*(sqr(m.displaySusyMu())+m.displayMh1Squared()
	  // 		 +m.displaySoftMassSquared(mQl, 1, 1)+m.displaySoftMassSquared(mDr, 1, 1))
	  //   -sqr(m.displaySoftA(DA, 1, 1));

	  // if (ccbTest < 0.0)
	  //   {
	  //     hasCCBProblem = true;
	  //   }

	  // // A_tau test
	  // ccbTest = 3.0*(sqr(m.displaySusyMu())+m.displayMh1Squared()
	  // 		 +m.displaySoftMassSquared(mLl, 3, 3)+m.displaySoftMassSquared(mEr, 3, 3))
	  //   -sqr(m.displaySoftA(EA, 3, 3));

	  // if (ccbTest < 0.0)
	  //   {
	  //     hasCCBProblem = true;
	  //   }

	  // ccbTest = 3.0*(sqr(m.displaySusyMu())+m.displayMh1Squared()
	  // 		 +m.displaySoftMassSquared(mLl, 2, 2)+m.displaySoftMassSquared(mEr, 2, 2))
	  //   -sqr(m.displaySoftA(EA, 2, 2));

	  // if (ccbTest < 0.0)
	  //   {
	  //     hasCCBProblem = true;
	  //   }

	  // ccbTest = 3.0*(sqr(m.displaySusyMu())+m.displayMh1Squared()
	  // 		 +m.displaySoftMassSquared(mLl, 1, 1)+m.displaySoftMassSquared(mEr, 1, 1))
	  //   -sqr(m.displaySoftA(EA, 1, 1));

	  // if (ccbTest < 0.0)
	  //   {
	  //     hasCCBProblem = true;
	  //   }

	  // // Calculate fine tuning
	  // hasTuningError = false;

	  // bool tuningEWSBProblem = false;

	  // try
	  //   {
	  //     tunings = doCalcpMSSMFineTuning(m, MS, tuningEWSBProblem, hasTuningError, USEAPPROXSOLNS, TOLEWSB);

	  //     // Possible numerical errors here in running should be accounted for.
	  //     if (!hasEWSBProblem && tuningEWSBProblem) //< tuning calculation hasn't found solution satisfies tolerance
	  // 	{
	  // 	  // Probably a numerical issue.
	  // 	  cerr << "WARNING: numerical issue in RG running for tuning calculation." << endl;
	  // 	  hasTuningError = true;
	  // 	}
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


	  // // Display output
	  // cout << tb_vals(i) << " ";
	  // cout << mu_vals(j) << " ";
	  // cout << B_vals(k) << " ";
	  // cout << mqL3sq_vals(l) << " ";
	  // cout << mtRsq_vals(p) << " ";
	  // cout << At_vals(q) << " ";
	  // cout << M2_vals(s) << " "; 
	  // cout << m.displayPhys().mh0 << " ";
	  // cout << m.displayPhys().mA0 << " ";
	  // cout << mHdSq << " ";
	  // cout << mHuSq << " ";
	  // cout << f1 << " ";
	  // cout << f2 << " ";
	  // cout << m.displayMsusy() << " ";
	  // cout << m.displayDrBarPars().mu(1,3) << " ";
	  // cout << m.displayDrBarPars().mu(2,3) << " ";
	  // cout << m.displayDrBarPars().md(1,3) << " ";
	  // cout << m.displayDrBarPars().md(2,3) << " ";
	  // cout << m.displayDrBarPars().mchBpmz(1) << " ";
	  // cout << m.displayDrBarPars().mchBpmz(2) << " ";
	  // cout << m.displayDrBarPars().mnBpmz(1) << " ";
	  // cout << m.displayDrBarPars().mnBpmz(2) << " ";
	  // cout << m.displayDrBarPars().mnBpmz(3) << " ";
	  // cout << m.displayDrBarPars().mnBpmz(4) << " ";
	  // cout << m.displayDrBarPars().mh0 << " ";
	  // cout << m.displayDrBarPars().mA0 << " ";
	  // cout << tunings.max() << " ";
	  // for (int i = 1; i <= tunings.displayEnd(); i++)
	  //   {
	  //     cout << tunings(i) << " ";
	  //   }
	  // if (hasUFBProblem)
	  //   {
	  //     cout << UFBPROBLEM << " ";
	  //   }
	  // else
	  //   {
	  //     cout << "0 ";
	  //   }
	  // if (hasCCBProblem)
	  //   {
	  //     cout << CCBPROBLEM << " ";
	  //   }
	  // else
	  //   {
	  //     cout << "0 ";
	  //   }

	  // if (hasEWSBProblem)
	  //   {
	  //     cout << EWSBPROBLEM << " ";
	  //   }
	  // else
	  //   {
	  //     cout << "0 ";
	  //   }

	  // if (hasWrongVacuum)
	  //   {
	  //     cout << WRONGVACUUM << " ";
	  //   }
	  // else
	  //   {
	  //     cout << "0 ";
	  //   }

	  // if (squarksTachyons)
	  //   {
	  //     cout << SQUARKTACHYON << " ";
	  //   }
	  // else
	  //   {
	  //     cout << "0 ";
	  //   }
	  // if (higgsTachyons)
	  //   {
	  //     cout << HIGGSTACHYON << " ";
	  //   }
	  // else
	  //   {
	  //     cout << "0 ";
	  //   }
	  // if (tadpoleProblem)
	  //   {
	  //     cout << TADPOLESPROBLEM << " ";
	  //   }
	  // else
	  //   {
	  //     cout << "0 ";
	  //   }
	  // if (m.displayPhys().mh0 < HIGGSCENT-HIGGSERROR || m.displayPhys().mh0 > HIGGSCENT+HIGGSERROR)
	  //   {
	  //     cout << NOTEXPVALID << " ";
	  //   }
	  // else
	  //   {
	  //     cout << "0 ";
	  //   }
	  // if (poleHiggsTachyons)
	  //   {
	  //     cout << POLEHIGGSTACHYON << " ";
	  //   }
	  // else
	  //   {
	  //     cout << "0 ";
	  //   }
	  // if (inaccurateHiggsMass)
	  //   {
	  //     cout << HIGGSPROBLEM << " ";
	  //   }
	  // else
	  //   {
	  //     cout << "0 ";
	  //   }

	  // if (hasTuningError)
	  //   {
	  //     cout << TUNINGERROR << " ";
	  //   }
	  // else
	  //   {
	  //     cout << "0 ";
	  //   }

	  // if (hasSeriousProblem)
	  //   {
	  //     cout << NUMERICALPROBLEM << endl;
	  //   }
	  // else
	  //   {
	  //     cout << "0" << endl;
	  //   }

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
  ii << "pmssmScanner called with incorrect arguments.\n";
  ii << "To run, you must provide as a file containing\n";
  ii << "the scan parameters you wish to use. Usage is:\n";
  ii << "./pmssmScanner < input_file.params\n";
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

flexiblesusy::genericE6SSM<Two_scale> doSimplifiedSpectrum(DoubleMatrix const & yuin, DoubleMatrix const & ydin, 
							   DoubleMatrix const & yein, DoubleVector const & gin, 
							   DoubleVector const & lambdain, DoubleVector const & kappain, 
							   double mupr, double tanb, double hvev, double svev,
							   DoubleVector const & mGaugino, DoubleMatrix const & tuin, 
							   DoubleMatrix const & tdin, DoubleMatrix const & tein, 
							   DoubleVector const & tlambdain, DoubleVector const & tkappain,
							   DoubleMatrix const & mQlSq, DoubleMatrix const & mUrSq, 
							   DoubleMatrix const & mDrSq, DoubleMatrix const & mLlSq, 
							   DoubleMatrix const & mErSq, DoubleMatrix const & mDxSq,
							   DoubleMatrix const & mDxbarSq, DoubleVector const & mHdISq, 
							   DoubleVector const & mHuISq, DoubleVector const & mSISq,
							   double Bmupr,  double mHpSq, double mHpbarSq, double mx, 
							   int l, int t, QedQcd const & dataset, 
							   double tol, double & mHdSq, double & mHuSq, double & mSSq, 
							   double & ms, 
							   bool & hasEWSBProblem, bool & squarksTachyons, 
							   bool & vectorBosonsTachyons, bool & higgsTachyons, 
							   bool & tadpoleProblem, bool & poleHiggsTachyons, 
							   bool & inaccurateHiggsMass, bool & hasSeriousProblem)
{
  using namespace flexiblesusy;

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

  Eigen::Matrix<double,3,3> Yu, Yd, Ye, Kappa, TYu, TYd, TYe, TKappa;
  Eigen::Matrix<double,2,2> Lambda12, TLambda12, mH1I2, mH2I2; 

  for (int i = 1; i <= 3; i++)
    {
      Kappa(i-1,i-1) = kappain.display(i);
      TKappa(i-1,i-1) = tkappain.display();
      for (int j = 1; j <= 3; j++)
	{
	  Yu(i-1,j-1) = yuin.display(i,j);
	  Yd(i-1,j-1) = ydin.display(i,j);
	  Ye(i-1,j-1) = yein.display(i,j);
	  TYu(i-1,j-1) = tuin(i,j);
	  TYd(i-1,j-1) = tdin(i,j);
	  TYe(i-1,j-1) = tein(i,j);
	}
    }

  double Lambdax = lambdain(3);
  double TLambdax = tlambdain.display(3);

  for (int i = 1; i <= 2; i++)
    {
      Lambda12(i-1,i-1) = lambdain(i);
      TLambda12(i-1,i-1) = tlambdain(i);
    }

  // Off diagonals for kappa, T_kappa, and the inert lambdas should vanish
  Kappa(0,1) = 0.0;
  Kappa(0,2) = 0.0;
  Kappa(1,0) = 0.0;
  Kappa(1,2) = 0.0;
  Kappa(2,0) = 0.0;
  Kappa(2,1) = 0.0;

  TKappa(0,1) = 0.0;
  TKappa(0,2) = 0.0;
  TKappa(1,0) = 0.0;
  TKappa(1,2) = 0.0;
  TKappa(2,0) = 0.0;
  TKappa(2,1) = 0.0;

  double v1 = hvev/Sqrt(1.0+tanb*tanb);
  double v2 = v1*tanb;

  // SUSY parameters
  genericE6SSM_susy_parameters r_susy = genericE6SSM_susy_parameters(mx, l, t, input, Yd, Ye, Kappa, Lambda12, Lambdax, Yu, mupr, gin(1), gin(2), gin(3), gin(4), v1, v2, svev);

  // Soft masses
  Eigen::Matrix<double,3,3> mq2, mu2, md2, ml2, me2, mDx2, mDxbar2;
  for (int i = 1; i <= 3; i++)
    {
      for (int j = 1; j <= 3; j++)
	{
	  mq2(i-1,j-1) = mQlSq(i,j);
	  mu2(i-1,j-1) = mUrSq(i,j);
	  md2(i-1,j-1) = mDrSq(i,j);
	  ml2(i-1,j-1) = mLlSq(i,j);
	  me2(i-1,j-1) = mErSq(i,j);
	  mDx2(i-1,j-1) = mDxSq(i,j);
	  mDxbar2(i-1,j-1) = mDxbarSq(i,j);
	}
    }
  
  Eigen::Matrix<double,2,2> mH1I2, mH2I2, msI2;

  for (int i = 1; i <= 2; i++)
    {
      mH1I2(i-1,i-1) = mHdISq.display(i);
      mH2I2(i-1,i-1) = mHuISq.display(i);
      msI2(i-1,i-1) = mSISq.display(i);
    }

  mH1I2(0,1) = 0.0; mH1I2(1,0) = 0.0;
  mH2I2(0,1) = 0.0; mH2I2(1,0) = 0.0;
  msI2(0,1) = 0.0; msI2(1,0) = 0.0;

  // Soft SUSY breaking parameters. Note ordering of gauginos is 1 = M_1,
  // 2 = M_2, 3 = M_3 and 4 = M_1'
  genericE6SSM_soft_parameters r_soft = genericE6SSM_soft_parameters(r_susy, TYd, TYe, TKappa, TLambda12, TLambdax, 
								     TYu, Bmupr, mq2, ml2, mHdSq, mHuSq, md2, mu2, me2,
								     mSSq, mH1I2, mH2I2, msI2, mDx2, mDxbar2, mHpSq, 
								     mHpbarSq, mGaugino(1), mGaugino(2), mGaugino(3), 
								     mGaugino(4));

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
      hasEWSBProblem = ESSM_ImplementEWSBConstraints_SoftMasses(r_approx, mx, ms, 
								false, solnGuess, tol);
      
      if (hasEWSBProblem)
	{
	  cerr << "WARNING: problem implementing EWSB conditions." << endl;
	  genericE6SSM<Two_scale> w = r_approx;
	  w.set_mHd2(solnGuess(1));
	  w.set_mHu2(solnGuess(2));
	  w.set_ms2(solnGuess(3));
	  w.run_to(solnGuess(4));
	  cerr << "f1 = " << ESSM_EWSBCondition1(w) << endl;
	  cerr << "f2 = " << ESSM_EWSBCondition2(w) << endl;
	}
      
      mHdSq = solnGuess(1);
      mHuSq = solnGuess(2);
      mSSq = solnGuess(3)
      ms = solnGuess(4);
      
      r_approx.set_mHd2(solnGuess(1));
      r_approx.set_mHu2(solnGuess(2));
      r_approx.set_ms2(solnGuess(3));
      
      // Try to get the spectrum. Catch any errors and flag as 
      // problem if they occur.
      // Get DR bar stop and sbottom masses.
      r_approx.runto(solnGuess(4));
      
      r_approx.calculate_MSu();
      r_approx.calculate_MSd();

      // Check for tachyonic squakrs      
      if (r_approx.get_problems().is_tachyon(Su) || r_approx.get_problem().is_tachyon(Sd))
	{
	  squarksTachyons = true;
	}
      
      // Values needed below, shouldn't be any problems here so long as tan(beta), v, and the gauge couplings
      // are sensible when input
      double v1 = r_approx.get_vd();
      double v2 = r_approx.get_vu();
      double s = r_approx.get_vs();

      double v = Sqrt(v1*v1+v2*v2);
      double tb = v2/v1;

      double sw = Sin(r_approx.ThetaW());
      double mt = r_approx.get_Yu(3, 3)*v*Sin(ArcTan(tb))/Sqrt(2.0);
      double mb = r_approx.get_Yd(3, 3)*v*Cos(ArcTan(tb))/Sqrt(2.0);
      double mz2 = Sqr(r_approx.calculate_MVZ());
      double mw2 = Sqr(r_approx.calculate_MVWm());
      double mzp2 = Sqr(r_approx.calculate_MVZp());

      if (r_approx.get_problems().is_tachyon(VZ) || r_approx.get_problems().is_tachyon(VZp) 
	  || r_approx.get_problems().is_tachyon(VWm) )
	{
	  vectorBosonsTachyons = true;
	}
      
      drBarPars eg = r_approx.displayDrBarPars();
      
      // Calculate neutralino and chargino DR bar masses, there shouldn't be any problems here 
      // (or if there are they will be serious numerical ones in diagonalisation)
      r_approx.calculate_MChi();
      r_approx.calculate_MCha();  

      // Calculate DR bar Higgs masses using FlexibleSUSY routines...
      r_approx.calcDrBarHiggs(atan(r_approx.displayTanb()), mz2, mw2, sw, eg);
      
      // Check if DR bar h0 or A0 are tachyons
      if (r_approx.displayProblem().tachyon == h0 || r_approx.displayProblem().tachyon == A0)
	{
	  higgsTachyons = true;
	}
      
      // ... but use approximate expressions here for loop corrected Higgs masses
      DoubleVector mstop(2), mstopsq(2), mD1sq(3), mD2sq(3);
      physical_ESSM(r_approx, mstop, mstopsq, mD1sq, mD2sq, s, tb);      

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

      poleHiggsTachyons = HiggsMasses(r_approx, s, tb, mstop, mstopsq, WhatCorrections, 
				      false, false, expBounds, ExpValid, mh, mhmix, msq, sing);

      // Other problems that can occur in physical Higgs mass calculation:
      // inaccurate result
      if (ExpValid == HIGGSPROBLEM || sing == 1)
	{
	  inaccurateHiggsMass = true;
	}
      if (ExpValid == TADPOLESPROBLEM)
	{
	  tadpoleProblem = true;
	}
      if (ExpValid == POLEHIGGSTACHYON)
	{
	  poleHiggsTachyons = true;
	}


    }
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
  
  return r_approx;
}
