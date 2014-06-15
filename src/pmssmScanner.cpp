/*
  mssmScanner is the driver routine for doing scans over
  the five parameters of interest to us in the pMSSM
  parameters space. Data is written to standard output,
  so you should pipe the output to a file. Values saved are
    -- Column 1: tan(beta)
    -- Column 2: mu
    -- Column 3: B
    -- Column 4: m_Q3^2
    -- Column 5: m_u3^2
    -- Column 6: A_t
    -- Column 7: M_2
    -- Column 8: m_h0 (pole)
    -- Column 9: M_{SUSY}
    -- Column 10: m_{~t_1}^DRbar(M_{SUSY})
    -- Column 11: m_{~t_2}^DRbar(M_{SUSY})
    -- Column 12: m_{~b_1}^DRbar(M_{SUSY})
    -- Column 13: m_{~b_2}^DRbar(M_{SUSY})
    -- Column 14: m_{\chi_1^\pm}^DRbar(M_{SUSY})
    -- Column 15: m_{\chi_2^\pm}^Drbar(M_{SUSY})
    -- Column 16: m_{\chi_1^0}^DRbar(M_{SUSY})
    -- Column 17: m_{\chi_2^0}^DRbar(M_{SUSY})
    -- Column 18: m_{\chi_3^0}^DRbar(M_{SUSY})
    -- Column 19: m_{\chi_4^0}^DRbar(M_{SUSY})
    -- Column 20: Fine tuning \Delta
    -- Column 21: Fine tuning \Delta_{M_1}
    -- Column 22: Fine tuning \Delta_{M_2}
    -- Column 23: Fine tuning \Delta_{M_3}
    -- Column 24: Fine tuning \Delta_{A_t}
    -- Column 25: Fine tuning \Delta_{A_b}
    -- Column 26: Fine tuning \Delta_{A_\tau}
    -- Column 27: Fine tuning \Delta_{m_Hd^2}
    -- Column 28: Fine tuning \Delta_{m_Hu^2}
    -- Column 29: Fine tuning \Delta_{mu}
    -- Column 30: Fine tuning \Delta_{B}
    -- Column 31: Fine tuning \Delta_{m_{L12}^2}
    -- Column 32: Fine tuning \Delta_{m_{e12}^2}
    -- Column 33: Fine tuning \Delta_{m_{Q12}^2}
    -- Column 34: Fine tuning \Delta_{m_{u12}^2}
    -- Column 35: Fine tuning \Delta_{m_{d12}^2}
    -- Column 36: Fine tuning \Delta_{m_{L3}^2}
    -- Column 37: Fine tuning \Delta_{m_{e3}^2}
    -- Column 38: Fine tuning \Delta_{m_{Q3}^2}
    -- Column 39: Fine tuning \Delta_{m_{u3}^2}
    -- Column 40: Fine tuning \Delta_{m_{d3}^2}
    -- Column 41: UFB flag
    -- Column 42: CCB flag
    -- Column 43: stop/sbottom tachyon flag
    -- Column 44: Higgs mass calculation flag
    -- Column 45: Higgs tachyon flag 
    -- Column 46: General problem flag (returned value is most serious)
 */

#include <iostream>
#include <sstream>
#include <cstring>
#include <mycomplex.h>
#include <def.h>
#include <linalg.h>
#include <lowe.h>
#include <rge.h>
#include <softsusy.h>
#include <softpars.h>
#include <physpars.h>
#include <susy.h>
#include <utils.h>
#include <numerics.h>
#include <twoloophiggs.h>
#include "../src/tuningnumerics.h"
#include "../src/mssmtuningutils.h"
#include "../src/tuningutils.h"
#include <sys/time.h>

using namespace softsusy;


void errorCall();
string ToUpper(const string & s);

// This function takes in a set of input parameters at the scale
// MX, uses an iteration to determine M_{SUSY}, and calculates the
// DR bar stops, sbottoms, charginos, neutralinos, Higgses at M_{SUSY},
// and an approximate value for the physical m_h0. Object is returned 
// at M_{SUSY}. Error flag is true if problem with iteration. Note
// m_Hd^2 and m_Hu^2 are fixed by EWSB. Initially mHdSq and mHuSq
// should contain initial guesses for the values of the soft Higgs
// masses. On return they contain the values that satisfy the 1-loop
// tadpole equations, if a solution could be found. Likewise for ms.
MssmSoftsusy doSimplifiedSpectrum(DoubleMatrix const & yuin, DoubleMatrix const & ydin, DoubleMatrix const & yein,
				  DoubleVector const & gin, double susyMu, double tanb, double hvev, 
				  DoubleVector const & mGaugino, DoubleMatrix const & au, DoubleMatrix const & ad, 
				  DoubleMatrix const & ae, DoubleMatrix const & mQlSq, DoubleMatrix const & mUrSq, 
				  DoubleMatrix const & mDrSq, DoubleMatrix const & mLlSq, DoubleMatrix const & mErSq, 
				  double m3sq, double mGrav, double mx, int l, int t, QedQcd const & dataset, 
				  sPhysical const & physpars, double tol, double & mHdSq, double & mHuSq, double & ms, 
				  bool & hasEWSBProblem, bool & squarksTachyons, bool & higgsTachyons, 
				  bool & tadpoleProblem, bool & poleHiggsTachyons, 
				  bool & inaccurateHiggsMass, bool & hasSeriousProblem);

int main(int argc, char* argv[])
{
  tryToConvergeHard = true;
  // Sets up exception handling
  signal(SIGFPE, FPE_ExceptionHandler); 

  /*
    ----------------------------------------------------
    SM parameters
    ----------------------------------------------------
  */
  double poleM_t = 173.2; // top quark pole mass in GeV
  double mbmb = 4.16; // bottom quark running mass m_b(m_b)^MSbar in GeV
  double poleM_tau = 1.777;// tau pole mass in GeV
  double poleM_Z = 91.1876; // Z boson pole mass in GeV
  double G_F = 1.16637900e-5; // G_F^MSbar in GeV^-2.
  double alphasmz = 0.1193; // alpha_s(mz)^MSbar (strong coupling constant = g_3^2/(4*PI))
  double alphaemmz = 127.9568; // alpha_em(mz)^MSbar (electromagnetic coupling = e^2/(4*PI))

  /*
    ----------------------------------------------------
  */

  /*
    ----------------------------------------------------
    Values for SUSY breaking parameters that are not 
    scanned over
    ----------------------------------------------------
   */
  double MX = 20000.0; // GeV

  double M1 = 300.0; // GeV
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

  /*
    ----------------------------------------------------
   */

  /*
    ----------------------------------------------------
    Values for scanned parameters
    ----------------------------------------------------
   */

  // Scan variables
  double tanb = 10.0;

  double mu = 200.0; // GeV
  double B = 500.0; // GeV
  double mqL3 = 500.0; // GeV
  double mtR = 500.0; // GeV
  double At = 1000.0; // GeV

  double M2 = 1000.0; // GeV

  // NB the scanned parameters are the squared masses
  // m_Q3^2 and m_u3^2, not m_Q3 and m_u3.  
  double mql3sq = sqr(mqL3);
  double mtRsq = sqr(mtR);
  
  /*
    ----------------------------------------------------
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
  int mu_npts = 1;
  int B_npts = 1;
  int mqL3sq_npts = 1;
  int mtRsq_npts = 1;
  int At_npts = 1;
  int M2_npts = 1;
  
  // Default lower bounds
  double tb_low = 2.0;
  double mu_low = -1000.0; // GeV
  double B_low = -1000.0; // GeV
  double At_low = -10000.0; // GeV
  double mqL3_low = 200.0; // GeV
  double mtR_low = 200.0; // GeV
  
  double mqL3sq_low = sqr(mqL3_low);
  double mtRsq_low = sqr(mtR_low);
  
  double M2_low = 100.0; // GeV

  // Default upper bounds
  double tb_up = 50.0;
  double mu_up = 1000.0; // GeV
  double B_up = 1000.0; // GeV
  double At_up = 10000.0; // GeV
  double mqL3_up = 2000.0; // GeV
  double mtR_up = 2000.0; // GeV
  
  double mqL3sq_up = sqr(mqL3_up);
  double mtRsq_up = sqr(mtR_up);
  
  double M2_up = 2000.0; // GeV

  // Default step sizes
  double tb_incr = 0.0;
  double mu_incr = 0.0;
  double B_incr = 0.0;
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

  int const LOOPS = 2; //< number of loops to use in RG running
  int const THRESH = 0; //< threshold accuracy
  int const NUMLOOPSHIGGS = 2; //< number of loops in Higgs pole mass calculation
  int const NUMLOOPSREWSB = 2; //< number of loops in SOFTSUSY EWSB conditions

  double const TOLEWSB = 50.0; //< tolerance to impose on minimisation conditions, GeV^2

  int const UFBPROBLEM = 1; //< unbounded from below flag
  int const CCBPROBLEM = 2; //< charge/colour breaking minimum flag

  int const EWSBPROBLEM = 3; //< problem with iteration

  int const WRONGVACUUM = 4; //< problem with unphysical vacuum

  int const SQUARKTACHYON = 5; //< stop/sbottom tachyon

  int const TADPOLESPROBLEM = 15; //< problem calculating 1-loop tadpoles

  int const HIGGSPROBLEM = 10; //< problems with calculating the Higgs mass
  int const NOTEXPVALID = 30; //< point not experimentally valid
  int const HIGGSTACHYON = 31; //< ruled out because tachyonic
  int const POLEHIGGSTACHYON = 32; //< ruled out because physical Higgs is tachyon

  int const NUMERICALPROBLEM = 666; //< for serious numerical problems

  int const TUNINGERROR = 35;

  double const HIGGSCENT = 125.0; //< rough central value for Higgs mass, GeV
  double const HIGGSERROR = 7.0; //< theory error in Higgs calculation, GeV

  bool const USEAPPROXSOLNS = true; //< use approximate solutions to RGEs

  bool uselogscanAt = false; //< use a symmetric log-scan over A_t

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
    SOFTSUSY flags and inputs
    ----------------------------------------------------
   */

  // Vector for input parameters.
  DoubleVector pars(49);//< 49 = non-universal model

  double mgutGuess = 2.0e16;
  int sgnMu = 0.0;
  bool gaugeUnification = false; //< no gauge unification condition, since non-universal
  bool ewsbBCscale = false; //< MX != M_{SUSY}
  bool flavourViolation = false; //< no flavour violation
  bool RPVflag = false; //< no R-parity violation
  bool useAlternativeEWSB = true; //< m_Hd^2 and m_Hu^2 are outputs
  bool setTbAtMX = true; //< tan(beta) is set at the input scale

  double desiredMh = 0.0;

  void (*boundaryCondition)(MssmSoftsusy &, const DoubleVector &)=extendedSugraBcs;

  char * modelIdent = (char *)"nonUniversal"; 

  // Global SOFTSUSY flags
  numHiggsMassLoops = NUMLOOPSHIGGS;
  numRewsbLoops = NUMLOOPSREWSB;

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
  double vin = 246.0; //< Higgs vev, GeV



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

  // Gaugino masses
  DoubleVector mGaugino(3);
  mGaugino(1) = M1;
  mGaugino(2) = M2;
  mGaugino(3) = M3;

  // Soft squared masses
  DoubleMatrix mQlSq(3,3), mUrSq(3,3), mDrSq(3, 3), mLlSq(3,3), mErSq(3,3);

  mLlSq(1,1) = sqr(meL);
  mLlSq(2,2) = sqr(mmuL);
  mLlSq(3,3) = sqr(mtauL);

  mErSq(1,1) = sqr(meR);
  mErSq(2,2) = sqr(mmuR);
  mErSq(3,3) = sqr(mtauR);

  mQlSq(1,1) = sqr(mqL1);
  mQlSq(2,2) = sqr(mqL2);
  mQlSq(3,3) = sqr(mqL3);

  mUrSq(1,1) = sqr(muR);
  mUrSq(2,2) = sqr(mcR);
  mUrSq(3,3) = sqr(mtR);

  mDrSq(1,1) = sqr(mdR);
  mDrSq(2,2) = sqr(msR);
  mDrSq(3,3) = sqr(mbR);

  // Make sure off-diagonals are zero, as is required
  // in the pMSSM as defined by Rizzo et al.

  mLlSq(1,2) = 0.0;
  mLlSq(1,3) = 0.0;
  mLlSq(2,1) = 0.0;
  mLlSq(2,3) = 0.0;
  mLlSq(3,1) = 0.0;
  mLlSq(3,2) = 0.0;

  mErSq(1,2) = 0.0;
  mErSq(1,3) = 0.0;
  mErSq(2,1) = 0.0;
  mErSq(2,3) = 0.0;
  mErSq(3,1) = 0.0;
  mErSq(3,2) = 0.0;

  mQlSq(1,2) = 0.0;
  mQlSq(1,3) = 0.0;
  mQlSq(2,1) = 0.0;
  mQlSq(2,3) = 0.0;
  mQlSq(3,1) = 0.0;
  mQlSq(3,2) = 0.0;

  mUrSq(1,2) = 0.0;
  mUrSq(1,3) = 0.0;
  mUrSq(2,1) = 0.0;
  mUrSq(2,3) = 0.0;
  mUrSq(3,1) = 0.0;
  mUrSq(3,2) = 0.0;

  mDrSq(1,2) = 0.0;
  mDrSq(1,3) = 0.0;
  mDrSq(2,1) = 0.0;
  mDrSq(2,3) = 0.0;
  mDrSq(3,1) = 0.0;
  mDrSq(3,2) = 0.0;

  // Soft trilinears = y*A
  DoubleMatrix au(3,3), ad(3,3), ae(3,3);

  // Soft Higgs masses
  double mHdSq = 1.0e6;
  double mHuSq = 1.0e6;

  // Gravitino mass
  double mGrav = 0.0;

  // M_{SUSY} = \sqrt{m_{~t_1}^DRbar(M_{SUSY})*m_{~t_2}^DRbar(M_{SUSY})}
  double MS = 1.0e3;

  // SM inputs...
  QedQcd oneset; 

  oneset.setAlpha(ALPHA, 1.0 / alphaemmz);			 
  oneset.setAlpha(ALPHAS, alphasmz);
  oneset.setMu(poleM_Z); MZ = poleM_Z;
  oneset.setMass(mBottom, mbmb);
  oneset.setPoleMt(poleM_t);
  oneset.setMass(mTau, poleM_tau); oneset.setPoleMtau(poleM_tau);

  QedQcd dataset = oneset;

  oneset.toMz();

  // To store pole masses...
  sPhysical physpars;

  // Model that we will calculate fine tuning for
  MssmSoftsusy m; 

  // Model that will be used to calculate the couplings and VEV
  MssmSoftsusy m_couplings;

  // Vectors for storing values of scanned variables
  DoubleVector tb_vals(1), mu_vals(1), B_vals(1), mqL3sq_vals(1), mtRsq_vals(1), At_vals(1), M2_vals(1);

  // Vector for storing fine tunings
  DoubleVector tunings(NUMPMSSMPARS);

  // Variables to store values of EWSB conditions
  double f1, f2;

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

  // Flags to indicate lower bounds are found
  bool hasTanbLL = false;
  bool hasMuLL = false;
  bool hasBLL = false;
  bool hasmqL3sqLL = false;
  bool hasmtRsqLL = false;
  bool hasAtLL = false;
  bool hasM2LL = false;

  // Flags to indicate upper bounds are found
  bool hasTanbUL = false;
  bool hasMuUL = false;
  bool hasBUL = false;
  bool hasmqL3sqUL = false;
  bool hasmtRsqUL = false;
  bool hasAtUL = false;
  bool hasM2UL = false;

  // Flags to indicate number of points are found
  bool hasTanbNpts = false;
  bool hasMuNpts = false;
  bool hasBNpts = false;
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
  bool hasg2 = false;
  bool hasg3 = false;
  bool hasv = false;

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
		  else if (word1 == "MULL")
		    {
		      kk >> mu_low;
		      hasMuLL = true;
		    }
		  else if (word1 == "MUUL")
		    {
		      kk >> mu_up;
		      hasMuUL = true;
		    }
		  else if (word1 == "MUNPTS")
		    {
		      kk >> mu_npts;
		      hasMuNpts = true;
		    }
		  else if (word1 == "BLL")
		    {
		      kk >> B_low;
		      hasBLL = true;
		    }
		  else if (word1 == "BUL")
		    {
		      kk >> B_up;
		      hasBUL = true;
		    }
		  else if (word1 == "BNPTS")
		    {
		      kk >> B_npts;
		      hasBNpts = true;
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

      if (!hasMuLL)
	{
	  cerr << "WARNING: lower mu limit not found: using default value " << mu_low << "." << endl;
	}

      if (!hasMuUL)
	{
	  cerr << "WARNING: upper mu limit not found: using default value " << mu_up << "." << endl;
	}

      if (!hasMuNpts)
	{
	  cerr << "WARNING: number of mu points not found: using default value " << mu_npts << "." << endl;
	}

      if (!hasBLL)
	{
	  cerr << "WARNING: lower B limit not found: using default value " << B_low << "." << endl;
	}

      if (!hasBUL)
	{
	  cerr << "WARNING: upper B limit not found: using default value " << B_up << "." << endl;
	}

      if (!hasBNpts)
	{
	  cerr << "WARNING: number of B points not found: using default value " << B_npts << "." << endl;
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

      if (mu_low > mu_up)
	{
	  temp = mu_low;
	  mu_low = mu_up;
	  mu_up = temp;
	}

      if (B_low > B_up)
	{
	  temp = B_low;
	  B_low = B_up;
	  B_up = temp;
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

      if (mu_npts < 1)
	{
	  cerr << "WARNING: invalid number of mu points " << mu_npts << " requested:";
	  cerr << " using default 1." << endl;
	  mu_npts = 1;
	}

      if (B_npts < 1)
	{
	  cerr << "WARNING: invalid number of B points " << B_npts << " requested:";
	  cerr << " using default 1." << endl;
	  B_npts = 1;
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
      if (mu_npts > 1) mu_vals.setEnd(mu_npts);
      if (B_npts > 1) B_vals.setEnd(B_npts);
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

      if (M2_npts > 1) M2_vals.setEnd(M2_npts);

      // Work out increments
      if (tb_npts != 1)
	{
	  tb_incr = (tb_up-tb_low)/(((double)tb_npts)-1.0);
	}

      if (mu_npts != 1)
	{
	  mu_incr = (mu_up-mu_low)/(((double)mu_npts)-1.0);
	}

      if (B_npts != 1)
	{
	  B_incr = (B_up-B_low)/(((double)B_npts)-1.0);
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

      for (int i = 1; i <= mu_npts; i++)
	{
	  mu_vals(i) = mu_low + (((double)i)-1.0)*mu_incr;
	}

      for (int i = 1; i <= B_npts; i++)
	{
	  B_vals(i) = B_low + (((double)i)-1.0)*B_incr;
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

      // Do the same for A_t, taking into account whether log scan or not
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

      /*
	----------------------------------------------------
      */
      
      /*
	----------------------------------------------------
	Get saved values of couplings 
	----------------------------------------------------
      */
      

      /*
	----------------------------------------------------
      */
      
      /*
	----------------------------------------------------
	Finally actually do the scan
	----------------------------------------------------
      */
      
      for (int i = 1; i <= tb_vals.displayEnd(); i++)
	{
      for (int j = 1; j <= mu_vals.displayEnd(); j++)
	{
      for (int k = 1; k <= B_vals.displayEnd(); k++)
	{
      for (int l = 1; l <= mqL3sq_vals.displayEnd(); l++)
	{
      for (int p = 1; p <= mtRsq_vals.displayEnd(); p++)
	{
      for (int q = 1; q <= At_vals.displayEnd(); q++)
	{
      for (int s = 1; s <= M2_vals.displayEnd(); s++)
	{

	  // Update variable values
	  mQlSq(3,3) = mqL3sq_vals(l);
	  mUrSq(3,3) = mtRsq_vals(p);
	  au(3,3) = yuin(3,3)*At_vals(q);
	  mGaugino(2) = M2_vals(s);

	  // Default guesses for values determined by iteration
	  mHdSq = 1.0e6;
	  mHuSq = 1.0e6;
	  MS = 1.0e3;

	  // Define model, first reset problem flags
	  hasEWSBProblem = false;
	  squarksTachyons = false;
	  higgsTachyons = false;
	  tadpoleProblem = false;
	  poleHiggsTachyons = false;
	  inaccurateHiggsMass = false;
	  hasSeriousProblem = false;

	  // m is returned at M_{SUSY}
	  m = doSimplifiedSpectrum(yuin, ydin, yein, gin, mu_vals(j), tb_vals(i), vin, 
				   mGaugino, au, ad, ae, mQlSq, mUrSq, mDrSq, mLlSq, mErSq, 
				   B_vals(k)*mu_vals(j), mGrav, MX, LOOPS, THRESH, oneset, 
				   physpars, TOLEWSB, mHdSq, mHuSq, MS, 
				   hasEWSBProblem, squarksTachyons, higgsTachyons, 
				   tadpoleProblem, poleHiggsTachyons, 
				   inaccurateHiggsMass, hasSeriousProblem);

	  f1 = MSSM_EWSBCondition1(m);
	  f2 = MSSM_EWSBCondition2(m);

	  // While at M_{SUSY}, check that we are not in a point where the
	  // unphysical vacuum v_1 = v_2 = 0 is stable.
	  hasWrongVacuum = false;
	  //double wrongVacuumTest = sqr(m.displayM3Squared())-((sqr(m.displaySusyMu())+m.displayMh1Squared())
	  //						      *(sqr(m.displaySusyMu())+m.displayMh2Squared()));

	  double wrongVacuumTest = (m.displayTanb()/(1.0-sqr(m.displayTanb())))*
	    (m.displayMh2Squared()*m.displayTanb()-m.displayMh1Squared()/m.displayTanb())-0.5*sqr(m.displayMzRun());

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
	  // the constrain m_1^2+m_2^2 >= 2|m_3^2| should hold at all
	  // scales greater than M_{SUSY}, and supposedly in particular at MX.
	  double ufbTest = m.displayMh1Squared() + m.displayMh2Squared() +2.0*sqr(m.displaySusyMu())
	    -2.0*fabs(m.displayM3Squared());

	  if (ufbTest < 0)
	    {
	      hasUFBProblem = true;
	    }

	  // For simplicity, apply the CCB criterion of hep-ph/9604417v3 Eq. (2a)
	  // and Eq. (5) of hep-ph/9507294v3 at M_{SUSY} (this might not be
	  // the optimum scale though).

	  double ccbTest = 3.0*(sqr(m.displaySusyMu())+m.displayMh2Squared()
				+m.displaySoftMassSquared(mQl, 3, 3)+m.displaySoftMassSquared(mUr, 3, 3))
	    -sqr(m.displaySoftA(UA, 3, 3));

	  if (ccbTest < 0.0)
	    {
	      hasCCBProblem = true;
	    }

	  // Check the other CCB conditions as well, even though the chance
	  // of them being a problem is small.
	  ccbTest = 3.0*(sqr(m.displaySusyMu())+m.displayMh2Squared()
			 +m.displaySoftMassSquared(mQl, 2, 2)+m.displaySoftMassSquared(mUr, 2, 2))
	    -sqr(m.displaySoftA(UA, 2, 2));

	  if (ccbTest < 0.0)
	    {
	      hasCCBProblem = true;
	    }

	  ccbTest = 3.0*(sqr(m.displaySusyMu())+m.displayMh2Squared()
			 +m.displaySoftMassSquared(mQl, 1, 1)+m.displaySoftMassSquared(mUr, 1, 1))
	    -sqr(m.displaySoftA(UA, 1, 1));

	  if (ccbTest < 0.0)
	    {
	      hasCCBProblem = true;
	    }

	  // A_b CCB tests
	  ccbTest = 3.0*(sqr(m.displaySusyMu())+m.displayMh1Squared()
			 +m.displaySoftMassSquared(mQl, 3, 3)+m.displaySoftMassSquared(mDr, 3, 3))
	    -sqr(m.displaySoftA(DA, 3, 3));

	  if (ccbTest < 0.0)
	    {
	      hasCCBProblem = true;
	    }

	  ccbTest = 3.0*(sqr(m.displaySusyMu())+m.displayMh1Squared()
			 +m.displaySoftMassSquared(mQl, 2, 2)+m.displaySoftMassSquared(mDr, 2, 2))
	    -sqr(m.displaySoftA(DA, 2, 2));

	  if (ccbTest < 0.0)
	    {
	      hasCCBProblem = true;
	    }

	  ccbTest = 3.0*(sqr(m.displaySusyMu())+m.displayMh1Squared()
			 +m.displaySoftMassSquared(mQl, 1, 1)+m.displaySoftMassSquared(mDr, 1, 1))
	    -sqr(m.displaySoftA(DA, 1, 1));

	  if (ccbTest < 0.0)
	    {
	      hasCCBProblem = true;
	    }

	  // A_tau test
	  ccbTest = 3.0*(sqr(m.displaySusyMu())+m.displayMh1Squared()
			 +m.displaySoftMassSquared(mLl, 3, 3)+m.displaySoftMassSquared(mEr, 3, 3))
	    -sqr(m.displaySoftA(EA, 3, 3));

	  if (ccbTest < 0.0)
	    {
	      hasCCBProblem = true;
	    }

	  ccbTest = 3.0*(sqr(m.displaySusyMu())+m.displayMh1Squared()
			 +m.displaySoftMassSquared(mLl, 2, 2)+m.displaySoftMassSquared(mEr, 2, 2))
	    -sqr(m.displaySoftA(EA, 2, 2));

	  if (ccbTest < 0.0)
	    {
	      hasCCBProblem = true;
	    }

	  ccbTest = 3.0*(sqr(m.displaySusyMu())+m.displayMh1Squared()
			 +m.displaySoftMassSquared(mLl, 1, 1)+m.displaySoftMassSquared(mEr, 1, 1))
	    -sqr(m.displaySoftA(EA, 1, 1));

	  if (ccbTest < 0.0)
	    {
	      hasCCBProblem = true;
	    }

	  m.runto(MX);

	  // Test at MX

	  ufbTest = m.displayMh1Squared() + m.displayMh2Squared() +2.0*sqr(m.displaySusyMu())
	    -2.0*fabs(m.displayM3Squared());

	  if (ufbTest < 0)
	    {
	      hasUFBProblem = true;
	    }

	  ccbTest = 3.0*(sqr(m.displaySusyMu())+m.displayMh2Squared()
			 +m.displaySoftMassSquared(mQl, 3, 3)+m.displaySoftMassSquared(mUr, 3, 3))
	    -sqr(m.displaySoftA(UA, 3, 3));

	  if (ccbTest < 0.0)
	    {
	      hasCCBProblem = true;
	    }

	  // Check the other CCB conditions as well, even though the chance
	  // of them being a problem is small.
	  ccbTest = 3.0*(sqr(m.displaySusyMu())+m.displayMh2Squared()
			 +m.displaySoftMassSquared(mQl, 2, 2)+m.displaySoftMassSquared(mUr, 2, 2))
	    -sqr(m.displaySoftA(UA, 2, 2));

	  if (ccbTest < 0.0)
	    {
	      hasCCBProblem = true;
	    }

	  ccbTest = 3.0*(sqr(m.displaySusyMu())+m.displayMh2Squared()
			 +m.displaySoftMassSquared(mQl, 1, 1)+m.displaySoftMassSquared(mUr, 1, 1))
	    -sqr(m.displaySoftA(UA, 1, 1));

	  if (ccbTest < 0.0)
	    {
	      hasCCBProblem = true;
	    }

	  // A_b CCB tests
	  ccbTest = 3.0*(sqr(m.displaySusyMu())+m.displayMh1Squared()
			 +m.displaySoftMassSquared(mQl, 3, 3)+m.displaySoftMassSquared(mDr, 3, 3))
	    -sqr(m.displaySoftA(DA, 3, 3));

	  if (ccbTest < 0.0)
	    {
	      hasCCBProblem = true;
	    }

	  ccbTest = 3.0*(sqr(m.displaySusyMu())+m.displayMh1Squared()
			 +m.displaySoftMassSquared(mQl, 2, 2)+m.displaySoftMassSquared(mDr, 2, 2))
	    -sqr(m.displaySoftA(DA, 2, 2));

	  if (ccbTest < 0.0)
	    {
	      hasCCBProblem = true;
	    }

	  ccbTest = 3.0*(sqr(m.displaySusyMu())+m.displayMh1Squared()
			 +m.displaySoftMassSquared(mQl, 1, 1)+m.displaySoftMassSquared(mDr, 1, 1))
	    -sqr(m.displaySoftA(DA, 1, 1));

	  if (ccbTest < 0.0)
	    {
	      hasCCBProblem = true;
	    }

	  // A_tau test
	  ccbTest = 3.0*(sqr(m.displaySusyMu())+m.displayMh1Squared()
			 +m.displaySoftMassSquared(mLl, 3, 3)+m.displaySoftMassSquared(mEr, 3, 3))
	    -sqr(m.displaySoftA(EA, 3, 3));

	  if (ccbTest < 0.0)
	    {
	      hasCCBProblem = true;
	    }

	  ccbTest = 3.0*(sqr(m.displaySusyMu())+m.displayMh1Squared()
			 +m.displaySoftMassSquared(mLl, 2, 2)+m.displaySoftMassSquared(mEr, 2, 2))
	    -sqr(m.displaySoftA(EA, 2, 2));

	  if (ccbTest < 0.0)
	    {
	      hasCCBProblem = true;
	    }

	  ccbTest = 3.0*(sqr(m.displaySusyMu())+m.displayMh1Squared()
			 +m.displaySoftMassSquared(mLl, 1, 1)+m.displaySoftMassSquared(mEr, 1, 1))
	    -sqr(m.displaySoftA(EA, 1, 1));

	  if (ccbTest < 0.0)
	    {
	      hasCCBProblem = true;
	    }

	  // Calculate fine tuning
	  hasTuningError = false;

	  bool tuningEWSBProblem = false;

	  try
	    {
	      tunings = doCalcpMSSMFineTuning(m, MS, tuningEWSBProblem, hasTuningError, USEAPPROXSOLNS, TOLEWSB);

	      // Possible numerical errors here in running should be accounted for.
	      if (!hasEWSBProblem && tuningEWSBProblem) //< tuning calculation hasn't found solution satisfies tolerance
		{
		  // Probably a numerical issue.
		  cerr << "WARNING: numerical issue in RG running for tuning calculation." << endl;
		  hasTuningError = true;
		}
	    }
	  catch(const string & a) 
	    { 
	      cerr << "WARNING: serious numerical problem encountered in fine tuning calculation." << endl; 
	      hasSeriousProblem = true; 
	    }
	  catch(const char * a) 
	    { 
	      cerr << "WARNING: serious numerical problem encountered in fine tuning calculation." << endl; 
	      hasSeriousProblem = true; 
	    }
	  catch(...) 
	    { 
	      cerr << "WARNING: unknown problem encountered in fine tuning calculation." << endl; 
	      hasSeriousProblem = true; 
	    }


	  // Display output
	  cout << tb_vals(i) << " ";
	  cout << mu_vals(j) << " ";
	  cout << B_vals(k) << " ";
	  cout << mqL3sq_vals(l) << " ";
	  cout << mtRsq_vals(p) << " ";
	  cout << At_vals(q) << " ";
	  cout << M2_vals(s) << " "; 
	  cout << m.displayPhys().mh0 << " ";
	  cout << m.displayPhys().mA0 << " ";
	  cout << mHdSq << " ";
	  cout << mHuSq << " ";
	  cout << f1 << " ";
	  cout << f2 << " ";
	  cout << m.displayMsusy() << " ";
	  cout << m.displayDrBarPars().mu(1,3) << " ";
	  cout << m.displayDrBarPars().mu(2,3) << " ";
	  cout << m.displayDrBarPars().md(1,3) << " ";
	  cout << m.displayDrBarPars().md(2,3) << " ";
	  cout << m.displayDrBarPars().mchBpmz(1) << " ";
	  cout << m.displayDrBarPars().mchBpmz(2) << " ";
	  cout << m.displayDrBarPars().mnBpmz(1) << " ";
	  cout << m.displayDrBarPars().mnBpmz(2) << " ";
	  cout << m.displayDrBarPars().mnBpmz(3) << " ";
	  cout << m.displayDrBarPars().mnBpmz(4) << " ";
	  cout << m.displayDrBarPars().mh0 << " ";
	  cout << m.displayDrBarPars().mA0 << " ";
	  cout << tunings.max() << " ";
	  for (int i = 1; i <= tunings.displayEnd(); i++)
	    {
	      cout << tunings(i) << " ";
	    }
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
	  if (m.displayPhys().mh0 < HIGGSCENT-HIGGSERROR || m.displayPhys().mh0 > HIGGSCENT+HIGGSERROR)
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

	} //< M_2 scan
	}//< A_t scan
	} //< m_u3^2 scan 
	} //< m_Q3^2 scan
	} //< B scan
	} //< mu scan
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

// This function takes in a set of input parameters at the scale
// MX, uses an iteration to determine M_{SUSY}, and calculates the
// DR bar stops, sbottoms, charginos, neutralinos, Higgses at M_{SUSY},
// and an approximate value for the physical m_h0. Object is returned 
// at M_{SUSY}. Error flag is true if problem with iteration. Note
// m_Hd^2 and m_Hu^2 are fixed by EWSB. Initially mHdSq and mHuSq
// should contain initial guesses for the values of the soft Higgs
// masses. On return they contain the values that satisfy the 1-loop
// tadpole equations, if a solution could be found. Likewise for ms.
MssmSoftsusy doSimplifiedSpectrum(DoubleMatrix const & yuin, DoubleMatrix const & ydin, DoubleMatrix const & yein,
				  DoubleVector const & gin, double susyMu, double tanb, double hvev, 
				  DoubleVector const & mGaugino, DoubleMatrix const & au, DoubleMatrix const & ad, 
				  DoubleMatrix const & ae, DoubleMatrix const & mQlSq, DoubleMatrix const & mUrSq, 
				  DoubleMatrix const & mDrSq, DoubleMatrix const & mLlSq, DoubleMatrix const & mErSq, 
				  double m3sq, double mGrav, double mx, int l, int t, QedQcd const & dataset, 
				  sPhysical const & physpars, double tol, double & mHdSq, double & mHuSq, double & ms, 
				  bool & hasEWSBProblem, bool & squarksTachyons, bool & higgsTachyons, 
				  bool & tadpoleProblem, bool & poleHiggsTachyons, 
				  bool & inaccurateHiggsMass, bool & hasSeriousProblem)
{

  // SUSY parameters
  MssmSusy r_susy = MssmSusy(yuin, ydin, yein, gin, susyMu, tanb, mx, l, t, hvev);
  
  // Soft SUSY breaking parameters
  SoftParsMssm r_soft = SoftParsMssm(r_susy, mGaugino, au, ad, ae, mQlSq, mUrSq, mDrSq, mLlSq, mErSq,
				     m3sq, mHdSq, mHuSq, mGrav, mx, l, t);

  MssmSoftsusy r_approx = MssmSoftsusy(r_soft, physpars, mx, l, t, hvev);

  QedQcd oneset = dataset;

  oneset.toMz();

  r_approx.setData(oneset);

  // M_{SUSY}, m_Hd^2(MX) and m_Hu^2(MX) must be determined
  // simultaneously using an iterative procedure. 

  DoubleVector solnGuess(3);

  solnGuess(1) = mHdSq;
  solnGuess(2) = mHuSq;
  solnGuess(3) = ms;

  r_approx.setMh1Squared(solnGuess(1));
  r_approx.setMh2Squared(solnGuess(2));
  r_approx.setMsusy(ms);

  // Try to determine the correct solution using Newton's method, 
  // catch any numerical errors and flag as a problem point if they occur.
  try
    {
      hasEWSBProblem = MSSM_ImplementEWSBConstraints_SoftMasses(r_approx, mx, ms, 
  							false, solnGuess, tol);
      
      if (hasEWSBProblem)
	{
	  cerr << "WARNING: problem implementing EWSB conditions." << endl;
	  MssmSoftsusy s = r_approx;
	  s.setMh1Squared(solnGuess(1));
	  s.setMh2Squared(solnGuess(2));
	  s.setMsusy(solnGuess(3));
	  s.runto(s.displayMsusy());
	  cerr << "f1 = " << MSSM_EWSBCondition1(s) << endl;
	  cerr << "f2 = " << MSSM_EWSBCondition2(s) << endl;
	}
      
      mHdSq = solnGuess(1);
      mHuSq = solnGuess(2);
      ms = solnGuess(3);
      
      r_approx.setMsusy(solnGuess(3));
      r_approx.setMh1Squared(solnGuess(1));
      r_approx.setMh2Squared(solnGuess(2));
      
      // Try to get the spectrum. Catch any errors and flag as 
      // problem if they occur.
      // Get DR bar stop and sbottom masses.
      r_approx.runto(r_approx.displayMsusy());
      
      r_approx.stops();
      r_approx.sbottoms();

      // Check for tachyonic squakrs      
      if (r_approx.displayProblem().tachyon == stop || r_approx.displayProblem().tachyon == sbottom)
	{
	  squarksTachyons = true;
	}
      
      // Values needed below, shouldn't be any problems here so long as tan(beta), v, and the gauge couplings
      // are sensible when input
      double sw = r_approx.calcSinthdrbar();
      double mt = r_approx.displayYukawaElement(YU, 3, 3)*r_approx.displayHvev()*sin(atan(r_approx.displayTanb()))/sqrt(2.0);
      double mb = r_approx.displayYukawaElement(YD, 3, 3)*r_approx.displayHvev()*cos(atan(r_approx.displayTanb()))/sqrt(2.0);
      double mz2 = sqr(r_approx.displayMzRun());
      double mw2 = sqr(r_approx.displayMwRun());
      
      drBarPars eg = r_approx.displayDrBarPars();
      
      // Calculate DR bar Higgs masses
      r_approx.calcDrBarHiggs(atan(r_approx.displayTanb()), mz2, mw2, sw, eg);
      
      // Check if DR bar h0 or A0 are tachyons
      if (r_approx.displayProblem().tachyon == h0 || r_approx.displayProblem().tachyon == A0)
	{
	  higgsTachyons = true;
	}

      // Calculate neutralino and chargino DR bar masses, there shouldn't be any problems here 
      // (or if there are they will be serious numerical ones in diagonalisation)
      r_approx.calcDrBarNeutralinos(atan(r_approx.displayTanb()), sqrt(mz2), sqrt(mw2), sw, eg);
      r_approx.calcDrBarCharginos(atan(r_approx.displayTanb()), sqrt(mw2), eg);  
      
      // Set values
      eg.mu(1,3) = r_approx.displayDrBarPars().mu(1,3);
      eg.mu(2,3) = r_approx.displayDrBarPars().mu(2,3);
      eg.md(1,3) = r_approx.displayDrBarPars().md(1,3);
      eg.md(2,3) = r_approx.displayDrBarPars().md(2,3);
      eg.mt = mt;
      eg.mb = mb;
      eg.mz = sqrt(mz2);
      eg.mw = sqrt(mw2);
      
      r_approx.setDrBarPars(eg);
      
      // Tadpoles using only 1-loop stop and sbottom contributions
      r_approx.calcTadpole1Ms1loopSimple(mt, sw);
      r_approx.calcTadpole2Ms1loopSimple(mt, sw);
      
      // Check for problems with tadpole calculation. There shouldn't
      // be any, but just in case...
      tadpoleProblem = r_approx.displayProblem().noMuConvergence;

      double piwwt = 0.0, pizzt = 0.0;//< unused in simplified routine
      int accuracy = 1;//< number of loops = 0 or 1 as of SoftSUSY-3.3.10.
      
      // Calculate physical Higgs masses and flag if h0 or A0 tachyon
      poleHiggsTachyons = r_approx.higgsSimple(accuracy, piwwt, pizzt);
      
      // Other problems that can occur in physical Higgs mass calculation:
      // inaccurate result
      inaccurateHiggsMass = r_approx.displayProblem().inaccurateHiggsMass;

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
