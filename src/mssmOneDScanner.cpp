/*
  mssmOneDScanner is a program to scan over one MSSM parameter while 
  keeping all others fixed. Calculates the fine tuning at each 
  point and saves the result. Uses SOFTSUSY for MSSM calculations
  (currently uses SOFTSUSY-3.3.10). 
  To run for pMSSM points, the syntax is 
  "./mssmOneDScanner [block] [param] [lower_bound] [upper_bound] [n_pts] < [input_file]"
  [block] is the parameter block in which the parameter to be scanned over is found,
  so e.g. could be MINPAR or EXTPAR.
  [param] is an integer corresponding to the SLHA code for that parameter. Ones of 
  interest to us are 3 = M3, 11 = A_t, 43 = m_Q3, and 46 = m_u3.
  [lower_bound] is the lower bound to use in the scan.
  [upper_bound] is the upper bound to use in the scan.
  [n_pts] is an integer specifying the number of points in the scan.
  [input_file] should be an SLHA input file that contains the values for all of the
  other parameters that are not being scanned over. 
  The program outputs the values of the parameter scanned over, the fine tuning
  at that point calculated using numerical running of the RGEs, and the fine tuning
  calculated using the approximate solutions to the RGEs. For some parameters
  the program also outputs elements of the spectrum for comparison (e.g. if 
  M_3 is scanned over the gluino mass is output, if m_Q3, m_u3 or A_t are
  scanned over the DR bar and physical stop masses are output).

  The syntax for mSUGRA points is slightly different. For such a point, 
  [block] should always be MINPAR. The values of [param] are then slightly
  different, though, and do not correspond to the SLHA codes used for mSUGRA.
  The parameters that can be scanned over are m0 (param = 1), m12 (param = 2)
  and A0 (param = 3). The values of mu0 and B0 are obtained as outputs from
  the EWSB conditions.
  

  Notes:
    - I have only included SUGRA or general MSSM simulations, ignoring GMSB and AMSB for now.
    - for now I am excluding the possibility of flavour or R-parity violation
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
#include "tuningnumerics.h"
#include "mssmtuningutils.h"
#include "tuningutils.h"

using namespace softsusy;

string ToUpper(const string & );
void errorCall();

const int NUMPARSPMSSM = 25;

int main(int argc, char* argv[])
{
  tryToConvergeHard = true;
  /// Sets up exception handling
  signal(SIGFPE, FPE_ExceptionHandler); 

  // Necessary parameters, with some default values.
  double poleM_t = 173.2; // top quark pole mass in GeV (currently not actually used)
  double mbmb = 4.16; // bottom quark running mass m_b(m_b)^MSbar in GeV (not used)
  double poleM_tau = 1.777;// tau pole mass in GeV (also not used)
  double poleM_Z = 91.1876; // Z boson pole mass in GeV
  double G_F = 1.16637900e-5; // G_F^MSbar in GeV^-2.
  double alphasmz = 0.1193; // alpha_s(mz)^MSbar (strong coupling constant = g_3^2/(4*PI))
  double alphaemmz = 1.0/127.9568; // alpha_em(mz)^MSbar (electromagnetic coupling = e^2/(4*PI))

  QedQcd oneset; 

  oneset.setAlpha(ALPHA, 1.0 / alphaemmz);			 
  oneset.setAlpha(ALPHAS, alphasmz);
  oneset.setMu(poleM_Z); MZ = poleM_Z;
  oneset.setMass(mBottom, mbmb);
  oneset.setPoleMt(poleM_t);
  oneset.setMass(mTau, poleM_tau); oneset.setPoleMtau(poleM_tau);

  // Vector for input parameters - set to have length 3 for CMSSM/SUGRA
  // by default, extend if necessary.
  DoubleVector pars(3);

  double mgutGuess = 2.0e16;
  int sgnMu = 1;
  bool gaugeUnification = true, ewsbBCscale = false;
  bool flavourViolation = false;
  bool RPVflag = false;
  bool useAlternativeEWSB = false;
  bool setTbAtMX = false;

  double desiredMh = 0.0;

  void (*boundaryCondition)(MssmSoftsusy &, const DoubleVector &)=sugraBcs;

  char * modelIdent = (char *)""; 

  double tanb = 10.0;

  // SUSY and soft SUSY breaking parameters in the DR bar scheme, 
  // at the input scale MX.

  double mx = 20000.0; // GeV

  double M1 = 500.0; // GeV
  double M2 = 500.0; // GeV
  double M3 = 500.0; // GeV

  double At = 5000.0; // GeV
  double Ab = 500.0; // GeV
  double Atau = 500.0; // GeV

  double mu = 200.0; // GeV

  // While not very convenient for our scan, I'll just use the pole
  // mass until I can rewrite the BCs to use the DR bar m_A^2(MX).
  double mApole = 700.0; // GeV

  double meL = 500.0; // GeV
  double mmuL = 500.0; // GeV
  double mtauL = 500.0; // GeV

  double meR = 500.0; // GeV
  double mmuR = 500.0; // GeV
  double mtauR = 500.0; // GeV

  double mqL1 = 500.0; // GeV
  double mqL2 = 500.0; // GeV
  double mqL3 = 500.0; // GeV

  double muR = 500.0; // GeV
  double mcR = 500.0; // GeV
  double mtR = 500.0; // GeV

  double mdR = 500.0; // GeV
  double msR = 500.0; // GeV
  double mbR = 500.0; // GeV

  // Define scan defaults.
  // (default parameter is M3 = EXTPAR 3).
  string scan_block = "EXTPAR";
  int scan_par = 3;
  double lower_bound = -3000.0; // GeV for dimensionful parameters
  double upper_bound = 3000.0; // GeV for dimensionful parameters
  int n_pts = 1;

  try 
    {
      /// Sets format of output: 8 decimal places.
      outputCharacteristics(8);
      
      // Ideally read in SM data and SUSY parameters from benchmark files. 
      // Check arguments are correct first. 
      if (argc < 6)
	{
	  errorCall();
	} 

      // I don't want to be scanning over anything that is not in the MINPAR
      // or EXTPAR blocks, so show an error if neither of these blocks
      // are requested.

      scan_block = argv[1];
      scan_block = ToUpper(scan_block);

      if (strcmp(scan_block.c_str(), "EXTPAR") && strcmp(scan_block.c_str(), "MINPAR"))
	{
	  ostringstream ii;
	  ii << "WARNING: currently can only scan over parameters contained in the\n";
	  ii << "         blocks MINPAR or EXTPAR, you requested a scan over a parameter\n";
	  ii << "         in the block " << scan_block << ": exiting without doing scan.\n";
	  throw ii.str();
	}

      // Get parameter code to scan over. Note that if invalid input is given
      // this will set scan_par = 0, so we'd end up scanning over MX. For now
      // we'll allow that, even if it isn't really a useful scan.
      scan_par = strtol(argv[2], NULL, 10);

      // Get lower and upper bounds.
      lower_bound = strtod(argv[3], NULL);
      upper_bound = strtod(argv[4], NULL);

      // Get number of points to scan over.
      n_pts = strtol(argv[5], NULL, 10);

      // Check that the requested number of points is valid; if not, set n_pts = 1.
      if (n_pts < 1)
	{
	  cerr << "WARNING: requested " << n_pts << " points in scan, need at least 1: " 
	       << "using default 1 point." << endl;
	  n_pts = 1;
	}

      // Check ordering of upper and lower bounds.
      if (upper_bound < lower_bound)
	{
	  double temp = upper_bound;
	  upper_bound = lower_bound;
	  lower_bound = temp;
	}

      // So far everything is OK, so read from the file if possible.
      // Data that we need from the file is the type of MODEL (mSUGRA, pMSSM, ...),
      // the SM inputs, the MINPARs and the EXTPARs. Basically just the code
      // from softpoint with extraneous details removed...
    bool flag = false;
    
    string line, block;
    int model;
    
    while (getline(cin,line)) {
      //	  mgutGuess = mgutCheck("unified", gaugeUnification);
      
      //	cout << line << endl;
      istringstream input(line); 
      string word1, word2;
      input >> word1;
      
      if (word1.find("#") == string::npos) 
	{ 
	  // read in another word if there's no comment
	  input >> word2; 
	  
	  if (ToUpper(word1) == "BLOCK")  
	    { 
	      block = ToUpper(word2);
	  
	    } 
	  else 
	    { // ought to be data
	      istringstream kk(line);
	      if (block == "MODSEL") 
		{
		  int i; kk >> i; 
	    
		  switch(i) 
		    {
		    case 1: kk >> model; 
		      switch(model) 
			{
			case 0: boundaryCondition = &extendedSugraBcs;
			  modelIdent = (char *)"nonUniversal";
			  break;
			case 1: 
			  pars.setEnd(3); 
			  boundaryCondition = &sugraBcs; 
			  modelIdent = (char *)"sugra";
			  break;
			default: 
			  ostringstream ii;
			  ii << "Scanner cannot yet do model " 
			     << model << ": terminal error\n";
			  throw ii.str();
			}
		      break;
		    case 4: int i; kk >> i;
		      switch(i) 
			{
			case 0: RPVflag = false;
			  break;
			default:
			  ostringstream ii;
			  ii << "MODSEL 4 1: cannot yet do R-parity violation: exiting" << endl;
			  throw ii.str();
			}
		      break;
		    case 6: int j; kk >> j;
		      switch(j) 
			{
			case 0: flavourViolation = false; break;
			default:
			  ostringstream ii;
			  ii << "MODSEL 6: cannot yet do flavour violation: exiting." << endl;
			  throw ii.str();
			}
		      break;
		    default:
		      cout << "# WARNING: don't understand first integer " 
			   << word1 << " " << word2 << " in block " << block
			   << ": ignoring it\n";
		      break;
		    }
		}
	      else if (block == "MINPAR") 
		{
		  int i; double d; kk >> i >> d; 
		  switch (i) 
		    {
		    case 3: tanb = d; break;
		    case 4: sgnMu = int(d); break;
		    default: 
		      switch(model) {
		      case 0:
			// SUGRA inputs to fill out the pheno MSSM case
			switch(i) 
			  {
			  case 1: pars(1) = d; break;
			  case 2: pars(2) = d; break;
			  case 5: pars(3) = d; break;
			  default: 
			    ostringstream ii;
			    ii << "Didn't understand pheno MSSM input " << i << endl;
			    break;
			  } 
			break;
		      case 1: // SUGRA inputs
			switch(i)
			  {
			  case 1: 
			    pars(1) = d; 
			    break;
			  case 2: 
			    pars(2) = d; 
			    break;
			  case 5: 
			    pars(3) = d; 
			    break;
			  default: 
			    ostringstream ii;
			    ii << "Didn't understand SUGRA input " << i << endl;
			    break;
			  } 
			break;
		      default: 
			ostringstream ii;
			ii << "Didn't understand model input " << model << endl;
			break;
		      }
		      break;
		    }
		}
	      // Adding non-minimal options. 
	      else if (block == "EXTPAR") 
		{
		  /// First, we want to convert our input to EXTPAR if we have
		  /// mSUGRA already
		  if (!strcmp(modelIdent, "sugra")) 
		    {
		      modelIdent = (char *)"nonUniversal";
		      boundaryCondition = &extendedSugraBcs;
		      double m0 = pars(1), m12 = pars(2), a0 = pars(3);
		      pars.setEnd(49);
		      int i; for (i=1; i<=3; i++) pars(i) = m12;
		      for (i=11; i<=13; i++) pars(i) = a0;
		      pars(21) = m0*m0; pars(22) = m0*m0;
		      for (i=31; i<=36; i++) pars(i) = m0;		    
		      for (i=41; i<=49; i++) pars(i) = m0;		    	
		    }
		  
		  if (!strcmp(modelIdent, "nonUniversal")) 
		    {
		      int i; double d; kk >> i >> d;  
		      /// First, put parameters that depend not on
		      /// flavoured/unflavoured input
		      if (i == 0) 
			{ 
			  mgutGuess = d;
			  gaugeUnification = false;
			  // setting Minput=-1 should yield MSSM BCs at MSUSY
			  if (fabs(d + 1.0) < EPSTOL) 
			    {
			      mgutGuess = 1.0e3;
			      ewsbBCscale = true;
			      QEWSB = 1.0;
			      if (gaugeUnification) 
				cout << "# Gauge unification ignored since pheno MSSM"
				     << " assumes BC set at QEWSB\n"; 
			      gaugeUnification = false;
			    }
			}
		      else if (i == 25) 
			{
			  tanb = d;
			  if (pars.displayEnd() != 49) pars.setEnd(49);
			  pars(i) = d;
			  setTbAtMX = true;
			} 
		      else if (i == 23 || i == 26 ) 
			{
			  useAlternativeEWSB = true;
			  if (i == 23) { mu = d; 
			  }
			  if (i == 26) mApole = d; 
			}
		      else if (!flavourViolation) 
			{
			  if ((i > 0 && i <=  3) || (i >= 11 && i <= 13) || 
			      (i >= 21 && i <= 23) || (i == 26 || i == 25) 
			      || (i >= 31 && i <= 36) || 
			      (i >= 41 && i <= 49)) 
			    {
			      if (pars.displayEnd() != 49) pars.setEnd(49);
			      pars(i) = d;
			    }
			  else 
			    {
			      cout << "WARNING: did not understand parameter " 
				   << i << " in flavoured EXTPAR inputs\n";
			    }
			}
		    }
		}
	      else if (block == "SMINPUTS") 
		{
		  int i; double d; kk >> i >> d; 
		  switch (i) 
		    {
		    case 1: oneset.setAlpha(ALPHA, 1.0 / d); break;
		    case 2: GMU = d; break;
		    case 3: oneset.setAlpha(ALPHAS, d); break; 
		    case 4: oneset.setMu(d); MZ = d; break;
		    case 5: oneset.setMass(mBottom, d); flag = true; 
		      oneset.setMbMb(d); break;
		    case 6: oneset.setPoleMt(d); break;
		    case 7: oneset.setMass(mTau, d); 
		      oneset.setPoleMtau(d); break;
		    case 11: oneset.setMass(mElectron, d);  break;
		    case 13: oneset.setMass(mMuon, d); break;
		    case 21: oneset.setMass(mDown, d);
		      break;
		    case 22: oneset.setMass(mUp, d); break;
		    case 23: oneset.setMass(mStrange, d); break;
		    case 24: oneset.setMass(mCharm, d); break;
		    default: 
		      cout << "# WARNING: Don't understand data input " << i 
			   << " " << d << " in block "
			   << block << ": ignoring it\n"; break;
		    }
		} 
	      else if (block == "SOFTSUSY") 
		{
		  int i; double d; kk >> i >> d;
		  switch(i) {
		  case 1: TOLERANCE = d; break;
		  case 2: 
		    MIXING = int(d); 
		    //if (MIXING > 0) flavourViolation = true;
		    break;
		  case 3: PRINTOUT = int(d); break;
		  case 4: QEWSB = d; break;
		  case 5: INCLUDE_2_LOOP_SCALAR_CORRECTIONS = 
		      bool(int(d+EPSTOL)); break;
		  case 6: outputCharacteristics(int(d+EPSTOL)-1); break;  
		  case 14: 
		    {
		      int num = int(d + EPSTOL);
		      if (num == 1) tryToConvergeHard = true;		  
		    }

		    break;
		  default:
		    cout << "# WARNING: Don't understand data input " << i 
			 << " " << d << " in block "
			 << block << ": ignoring it\n"; break;
		  }
		}
	      else if (block == "TUNING") // < this block is specific to our code - it allows to turn on and off tadpole contributions
		{
		  int i; double d; kk >> i >> d;
		  switch(i)
		    {
		    case 1:
		      {
			INCLUDE1LPTADPOLES = bool(int(d+EPSTOL)); break;
		      }
		  default:
		    cout << "# WARNING: Don't understand data input " << i 
			 << " " << d << " in block "
			 << block << ": ignoring it\n"; break;
		    }
		}
	      else 
		{
		  cout << "# WARNING: don't recognise block " << block 
		       << ": ignoring all data in it" << endl;
		}
	      // end if blocks
	      
	    } // end of data
	} // end of no-comment
      
    } // end of file

    // intput error checking  
    if (sgnMu != 1 && sgnMu != -1 && sgnMu != 0) {
      ostringstream ii;
      ii << "Incorrect input for sign(mu)=" << sgnMu <<endl;
      throw ii.str();
    }
    if (tanb < 1.0 || tanb > 5.0e2)  {
      ostringstream ii;
      ii << "Incorrect input for tan beta=" << tanb <<endl;
      throw ii.str();
    }

    // Check that the requested block and SLHA code are consistent
    if (!strcmp(modelIdent, "nonUniversal"))
      {
	// For pMSSM models, we only allow scans over tan(beta) (MINPAR 3)
	if (!strcmp(scan_block.c_str(), "MINPAR"))
	  {
	    if (scan_par != 3)
	      {
		ostringstream ii;
		ii << "WARNING: invalid MINPAR " << scan_par << " request for non-universal model: exiting" << endl;
		throw ii.str();
	      }
	  }
	else if (!strcmp(scan_block.c_str(), "EXTPAR"))
	  {
	    // Note that this flags scans over MX as illegal, and also currently scans over m_A^2(MX).
	    if (!((scan_par > 0 && scan_par <=  3) || (scan_par >= 11 && scan_par <= 13) || 
		  (scan_par >= 21 && scan_par <= 23) || (scan_par == 26 || scan_par == 25) 
		  || (scan_par >= 31 && scan_par <= 36) || 
		  (scan_par >= 41 && scan_par <= 49)))
	      {
		ostringstream ii;
		ii << "WARNING: invalid EXTPAR " << scan_par << " request for non-universal model: exiting" << endl;
		throw ii.str();
	      }
	  }
      }
    else if (!strcmp(modelIdent, "sugra"))
      {
	if (!strcmp(scan_block.c_str(), "MINPAR"))
	  {
	    if (scan_par != 1 && scan_par != 2 && scan_par != 3 && scan_par != 5)
	      {
		ostringstream ii;
		ii << "WARNING: invalid MINPAR " << scan_par << " request for SUGRA model: exiting" << endl;
		throw ii.str();
	      }
	  }
	else if (!strcmp(scan_block.c_str(), "EXTPAR"))
	  {
	    ostringstream ii;
	    ii << "WARNING: SUGRA models should only scan over MINPAR parameters: exiting." << endl;
	    throw ii.str();
	  }
      }
    else
      {
	ostringstream ii;
	ii << "WARNING: have not yet implemented GMSB or AMSB models: exiting." << endl;
	throw ii.str();
      }
    
    QedQcd dataSet = oneset;

    oneset.toMz();

    // Check if we are scanning over tan(beta) in the MINPAR block or not.
    bool doScanOverTbAtMz = false;
    if (!strcmp(scan_block.c_str(), "MINPAR"))
      {
	if (scan_par == 3 && setTbAtMX == false)
	  {
	    doScanOverTbAtMz = true;
	  }
	else if (scan_par == 3 && setTbAtMX == true)
	  {
	    // In this case the value of tan(beta) at M_Z
	    // is ignored, so instead scan over EXTPAR 25.
	    scan_par = 25; 
	  }
      }

    // If we are in a SUGRA model and the user wants to scan over
    // A, need to change scan_par to 5 to match up with vectors below.
    if (!strcmp(modelIdent, "sugra"))
      {
    	if (scan_par == 5)
    	  {
    	    scan_par = 3;
    	  }
      }

    if (useAlternativeEWSB) 
      {
	sgnMu = 0; // Flags different BCs
      }
    
    // Do the scan.
    if (doScanOverTbAtMz)
      {
	tanb = lower_bound;
      }
    else
      {
	pars(scan_par) = lower_bound;
      }

    double par_incr = 0.0;

    if (n_pts != 1)
      {
	par_incr = (upper_bound-lower_bound)/(((double)n_pts)-1.0);
      }    

    double ms;

    // Vector for storing the values of the parameters we calculate the 
    // fine tunings for. Default is SUGRA, for which the parameters
    // in the vector are \{ m0^2, m12, mu, B0, A0 \} at the
    // scale MX.
    DoubleVector tuningPars(5);

    // If we want to calculate the fine tuning in the pMSSM, this
    // needs to be extended. In this case the parameters in the vector are
    // \{ M_1, M_2, M_3, A_t, A_b, A_tau, m_Hd^2, m_Hu^2, mu, B, 
    //    M_{L12}, M_{e12}, M_{Q12}, M_{u12}, M_{d12},
    //    M_{L3}, M_{e3}, M_{Q3}, M_{u3}, M_{d3}\}
    if (pars.displayEnd() == 49)
      {
	tuningPars.setEnd(NUMPMSSMPARS);
      }

    DoubleVector fineTunings = tuningPars;
    DoubleVector fineTuningsApprox = tuningPars;
    DoubleVector fineTuningsNumerical = tuningPars;

    // Vector for storing the VEVs at M_SUSY
    DoubleVector vevs(2);

    int sing;
    bool problemFlag;


    for (int j = 0; j < n_pts; j++)
      {
	sing = 0;
	problemFlag = false;

	if (doScanOverTbAtMz)
	  {
	    tanb = lower_bound+((double)j)*par_incr;
	  }
	else
	  {
	    pars(scan_par) = lower_bound + ((double)j)*par_incr;
	  }

	// Initialise model
	MssmSoftsusy r;
	if (setTbAtMX)
	  {
	    r.setSetTbAtMX(true);
	  }
	if (useAlternativeEWSB)
	  {
	    r.useAlternativeEwsb();
	    r.setMuCond(mu);
	    r.setSusyMu(mu);
	    r.setMaCond(mApole);
	  }
	if (scan_par == 23 && useAlternativeEWSB)
	  {
	    r.setMuCond(pars(scan_par));
	    r.setSusyMu(pars(scan_par));
	  }
	if (scan_par == 26 && useAlternativeEWSB)
	  {
	    r.setMaCond(pars(scan_par));
	  }

	r.setData(dataSet); // initial SM inputs


	// Calculate spectrum
	r.lowOrg(boundaryCondition, mgutGuess, pars, sgnMu, tanb, oneset, 
		 gaugeUnification, ewsbBCscale);
	

	// Find M_SUSY; this is the scale the fine tuning is calculated at.
	ms = r.displayMsusy();

	// Find the scale at which the BCs are applied.
	mx = r.displayMxBC();

	// cout << "ms = " << ms << endl;
	// cout << "mx = " << mx << endl;

	// Since I want to use the full evolution to calculate the spectrum, I will use
	// a separate object to calculate the fine tuning with just in case I have misunderstood
	// how the calculation works.
	MssmSoftsusy s = r;

	s.runto(mx);

	if (tuningPars.displayEnd() == NUMPMSSMPARS)
	  {
	    int maxIters = 200;
	    double tol = 1.0; // < this is reasonably low tolerance

	    bool hasEwsbProblem = false;
	    bool hasError = false;
	    DoubleVector fineTunings = doCalcpMSSMFineTuning(s, ms, hasEwsbProblem, hasError, false, tol);
	    DoubleVector fineTuningsApprox = doCalcpMSSMFineTuning(s, ms, hasEwsbProblem, hasError, true, tol);

	    DoubleVector updatedSoln(3);
	    bool problemFlag = MSSM_ImplementEWSBConstraints_SoftMasses(s, mx, ms, false,
									updatedSoln, tol);
	    s.setMh1Squared(updatedSoln(1));
	    s.setMh2Squared(updatedSoln(2));

	    tuningPars(1) = s.displayGaugino(1);
	    tuningPars(2) = s.displayGaugino(2);
	    tuningPars(3) = s.displayGaugino(3);
	    tuningPars(4) = s.displaySoftA(UA, 3, 3);
	    tuningPars(5) = s.displaySoftA(DA, 3, 3);
	    tuningPars(6) = s.displaySoftA(EA, 3, 3);
	    tuningPars(7) = s.displayMh1Squared();
	    tuningPars(8) = s.displayMh2Squared();
	    tuningPars(9) = s.displaySusyMu();
	    tuningPars(10) = s.displayM3Squared()/s.displaySusyMu();
	    tuningPars(11) = s.displaySoftMassSquared(mLl, 1, 1);
	    tuningPars(12) = s.displaySoftMassSquared(mEr, 1, 1);
	    tuningPars(13) = s.displaySoftMassSquared(mQl, 1, 1);
	    tuningPars(14) = s.displaySoftMassSquared(mUr, 1, 1);
	    tuningPars(15) = s.displaySoftMassSquared(mDr, 1, 1);
	    tuningPars(16) = s.displaySoftMassSquared(mLl, 3, 3);
	    tuningPars(17) = s.displaySoftMassSquared(mEr, 3, 3);
	    tuningPars(18) = s.displaySoftMassSquared(mQl, 3, 3);
	    tuningPars(19) = s.displaySoftMassSquared(mUr, 3, 3);
	    tuningPars(20) = s.displaySoftMassSquared(mDr, 3, 3);

	    DoubleVector fineTuningsNumerical = doCalcMSSMTuningNumerically(s, updatedSoln(3), mx, tuningPars, pMSSMftBCs);

	    s.runto(updatedSoln(3));
	    
	    if (problemFlag)
	      {
		cerr << "WARNING: failed to implement EWSB constraints: f1 = " <<  MSSM_EWSBCondition1(s) 
		     << ", f2 = " << MSSM_EWSBCondition2(s) << endl;
	      }

	    if (hasError)
	      {
		cerr << "WARNING: calculated fine tunings are inaccurate." << endl;
	      }

	DoubleVector softsusyTunings = s.fineTune(boundaryCondition, pars, mx, false);

	// If there are no problems, print the result
	if (!r.displayProblem().test() && problemFlag == false) 
	  {
	    // Results printed out depend on parameter scanner.
	    // For M_3, print value of M_3, value of physical gluino mass,
	    // fine tuning, lightest Higgs mass
	    
	    if (scan_par == 3)
	      {
		cout << pars(scan_par) << " " 
		     << r.displayPhys().mGluino << " "
		     << r.displayPhys().mh0 << " " 
		     << r.displayPhys().mch(1) << " " << r.displayPhys().mch(2) << " "
		     << r.displayPhys().mneut(1) << " " << r.displayPhys().mneut(2) << " "
		     << fineTunings.max() << " " << fineTuningsApprox.max() << " " << fineTuningsNumerical.max() << " "
		     << fineTunings(3) << " " << fineTuningsApprox(3) << " " << fineTuningsNumerical(3) << " " 
		     << softsusyTunings.max() << endl;
	      }
	    else if (scan_par == 43 || scan_par == 46)
	      {
		cout << pars(scan_par) << " " 
		     << r.displayPhys().mu(1,3) << " "
		     << r.displayPhys().mu(2,3) << " "
		     << r.displayPhys().mh0 << " "
		     << fineTunings.max() << " " << fineTuningsApprox.max() << " ";
		if (scan_par == 43)
		  {
		    cout << fineTunings(18) << " " << fineTuningsApprox(18) << " ";
		  }
		else
		  {
		    cout << fineTunings(19) << " " << fineTuningsApprox(19) << " ";
		  }
		cout << softsusyTunings.max() << endl;
	      }
	    else
	      {
		cout << pars(scan_par) << " " 
		     << r.displayPhys().mh0 << " " 
		     << r.displayPhys().mA0 << " " 
		     << r.displayPhys().mH0 << " " 
		     << r.displayPhys().mHpm << " "
		     << fineTunings.max() << " " << fineTuningsApprox.max() << " " << fineTuningsApprox.max() << " " << fineTuningsNumerical.max() << " "<< softsusyTunings.max() << endl;
	      }
	  }
	else if (r.displayProblem().test() && problemFlag == false) 
	  {
	    /// print out what the problem(s) is(are)
	    cerr << pars(scan_par) << " " << r.displayProblem() << endl;
	    /// in out stream might also print flag to indicate problem? 
	  }
	else if (problemFlag || hasError)
	  {
	    cerr << "WARNING: problem with fine tuning calculation." << endl;
	    cerr << pars(scan_par) << " "
		 << fineTunings.max() << " " << fineTuningsApprox.max() << " " << fineTuningsNumerical.max() << " "<< softsusyTunings.max() << endl;
	  }

	  }
	else
	  {
	    int maxIters = 200;
	    double tol = 1.0; // < this is reasonably low tolerance
	    bool hasError = false;
	    DoubleVector fineTunings = doCalcSUGRAFineTuning(s, ms, hasError, tol);

	    DoubleVector updatedSoln(3);

	    bool problemFlag = MSSM_ImplementEWSBConstraints_SUGRA(s, mx, ms,
								   updatedSoln, tol);
	    s.setSusyMu(updatedSoln(1));
	    s.setM3Squared(updatedSoln(2)*updatedSoln(1));

	    tuningPars(1) = s.displayMh1Squared();
	    tuningPars(2) = s.displayGaugino(1);
	    tuningPars(3) = s.displaySusyMu();
	    tuningPars(4) = s.displayM3Squared()/s.displaySusyMu();
	    tuningPars(5) = s.displaySoftA(UA, 3, 3);

	    DoubleVector fineTuningsNumerical = doCalcMSSMTuningNumerically(s, updatedSoln(3), mx, tuningPars, mSUGRAftBCs);

	    s.runto(updatedSoln(3));
	    
	    if (problemFlag)
	      {
		cerr << "WARNING: failed to implement EWSB constraints: f1 = " <<  MSSM_EWSBCondition1(s) 
		     << ", f2 = " << MSSM_EWSBCondition2(s) << endl;
	      }

	    if (hasError)
	      {
		cerr << "WARNING: calculated fine tunings are inaccurate." << endl;
	      }

	    DoubleVector softsusyTunings = s.fineTune(boundaryCondition, pars, mx, false);
	    
	    // If there are no problems, print the result
	    if (!r.displayProblem().test() && problemFlag == false) 
	      {
		cout << pars(scan_par)<< " " << r.displayPhys().mh0 << " " << fineTunings.max() << " " << fineTuningsNumerical.max() << " " << softsusyTunings.max() << endl;
	      }
	    else if (r.displayProblem().test() && problemFlag == false) 
	      {
		/// print out what the problem(s) is(are)
		cerr << pars(scan_par) << " " << r.displayProblem() << endl;
		/// in out stream might also print flag to indicate problem? 
	      }
	    else if (problemFlag || hasError)
	      {
		cerr << "WARNING: problem with fine tuning calculation." << endl;
		cerr << pars(scan_par) << " "
		     << fineTunings.max() << " " << fineTuningsNumerical.max() << " "<< softsusyTunings.max() << endl;
	      }
	    
	  }
      }
    }

  catch(const string & a) { cout << a; return -1; }
  catch(const char * a) { cout << a; return -1; }
  catch(...) { cout << "Unknown type of exception caught.\n"; return -1; }
  
  return 0;
}

// Returns a string with all characters in upper case: very handy
string ToUpper(const string & s) 
{
  string result;
  unsigned int index;
  for (index = 0; index < s.length(); index++) 
    {
      char a = s[index];
      a = toupper(a);
      result = result + a;
    }
  
  return result;
}

void errorCall() 
{
  ostringstream ii;
  ii << "mssmOneDScanner called with incorrect arguments. Calling syntax for pMSSM points is: \n";
  ii << "mssmOneDScanner [block] [param] [lower_bound] [upper_bound] [n_pts] < [input_file]\n";
  ii << "[block] should be a string describing the SLHA parameter block in which we can find\n";
  ii << "the parameter to be scanned over, e.g. MINPAR or EXTPAR.\n";
  ii << "[param] should be an integer corresponding to the SLHA code for the parameter to scan.\n";
  ii << "For example, EXTPAR 3 is M_3, EXTPAR 43 is m_Q3 etc.\n";
  ii << "[lower_bound] is the value of the lower bound for the scan.\n";
  ii << "[upper_bound] is the value of the upper bound.\n";
  ii << "[n_pts] is an integer corresponding to the number of points to scan over.\n";
  ii << "[input_file] should be an SLHA input file containing the values of the other \n";
  ii << "parameters that are not being scanned over.\n";
  throw ii.str();
}

	/* Commented out here 5/5/2014

	DoubleVector updatedSoln(3);
	int maxIters = 200;
	double tol = 1.0; // < this is reasonably low tolerance

	if (tuningPars.displayEnd() == NUMPMSSMPARS)
	  {
	    problemFlag = MSSM_ImplementEWSBConstraints_SoftMasses(s, mx, ms, false,
								   updatedSoln, maxIters, tol);
	    s.setMh1Squared(updatedSoln(1));
	    s.setMh2Squared(updatedSoln(2));
	  }
	else if (tuningPars.displayEnd() == 5)
	  {
	    // For SUGRA points it is mu and B0 that are fixed by EWSB, instead of 
	    // mHu^2 and mHd^2.
	    problemFlag = MSSM_ImplementEWSBConstraints_SUGRA(s, mx, ms, updatedSoln, maxIters, tol);

	    s.setSusyMu(updatedSoln(1));
	    s.setM3Squared(updatedSoln(1)*updatedSoln(2));

	  }

	// Get the values of the parameters at the input scale MX - this is horribly
	// messy, but it works...
	if (tuningPars.displayEnd() == 5) // sugra case
	  {
	    tuningPars(1) = s.displayMh1Squared();
	    tuningPars(2) = s.displayGaugino(1);
	    tuningPars(3) = s.displaySusyMu();
	    tuningPars(4) = s.displayM3Squared()/s.displaySusyMu();
	    tuningPars(5) = s.displaySoftA(UA, 3, 3);
	  }
	else if (tuningPars.displayEnd() == NUMPMSSMPARS) // pMSSM case
	  {
	    tuningPars(1) = s.displayGaugino(1);
	    tuningPars(2) = s.displayGaugino(2);
	    tuningPars(3) = s.displayGaugino(3);
	    tuningPars(4) = s.displaySoftA(UA, 3, 3);
	    tuningPars(5) = s.displaySoftA(DA, 3, 3);
	    tuningPars(6) = s.displaySoftA(EA, 3, 3);
	    tuningPars(7) = s.displayMh1Squared();
	    tuningPars(8) = s.displayMh2Squared();
	    tuningPars(9) = s.displaySusyMu();
	    tuningPars(10) = s.displayM3Squared()/s.displaySusyMu();
	    tuningPars(11) = s.displaySoftMassSquared(mLl, 1, 1);
	    tuningPars(12) = s.displaySoftMassSquared(mEr, 1, 1);
	    tuningPars(13) = s.displaySoftMassSquared(mQl, 1, 1);
	    tuningPars(14) = s.displaySoftMassSquared(mUr, 1, 1);
	    tuningPars(15) = s.displaySoftMassSquared(mDr, 1, 1);
	    tuningPars(16) = s.displaySoftMassSquared(mLl, 3, 3);
	    tuningPars(17) = s.displaySoftMassSquared(mEr, 3, 3);
	    tuningPars(18) = s.displaySoftMassSquared(mQl, 3, 3);
	    tuningPars(19) = s.displaySoftMassSquared(mUr, 3, 3);
	    tuningPars(20) = s.displaySoftMassSquared(mDr, 3, 3);
	  }
	else
	  {
	    ostringstream ii;
	    ii << "WARNING: unexpected number of fine tuning parameters encountered: exiting." << endl;
	    throw ii.str();
	  }

	// Also save gauge and Yukawa couplings at this scale
	DoubleVector highScaleCouplings(6);
	highScaleCouplings(1) = s.displayGaugeCoupling(1);
	highScaleCouplings(2) = s.displayGaugeCoupling(2);
	highScaleCouplings(3) = s.displayGaugeCoupling(3);
	highScaleCouplings(4) = s.displayYukawaElement(YU, 3, 3);
	highScaleCouplings(5) = s.displayYukawaElement(YD, 3, 3);
	highScaleCouplings(6) = s.displayYukawaElement(YE, 3, 3);


	// Also get the numerical fine tunings.
	if (tuningPars.displayEnd() == NUMPMSSMPARS)
	  {
	    fineTuningsNumerical = doCalcMSSMTuningNumerically(s, updatedSoln(3), mx, tuningPars, pMSSMftBCs);
	  }
	else
	  {
	    fineTuningsNumerical = doCalcMSSMTuningNumerically(s, updatedSoln(3), mx, tuningPars, mSUGRAftBCs);
	  }

	//cout << fineTuningsNumerical;
	//cout << "Fine tuning (numerical) = " << fineTuningsNumerical.max() << endl;
	s.runto(updatedSoln(3));

	if (problemFlag)
	  {
	    cerr << "WARNING: failed to implement EWSB constraints: f1 = " <<  MSSM_EWSBCondition1(s) 
		 << ", f2 = " << MSSM_EWSBCondition2(s) << endl;
	  }

	DoubleVector vevs(2);
	vevs(1) = s.displayHvev()/sqrt(1.0+sqr(s.displayTanb()));
	vevs(2) = vevs(1)*s.displayTanb();

	// Calculate fine tuning (using numerical running and also approximation) - also note that we should
	// calculate the fine tuning due to all parameters, so that
	// we can see where the scanned parameter is the dominant contribution to
	// fine tuning.

	// Next we want to combine all of these steps into a single function call.

	if (tuningPars.displayEnd() == NUMPMSSMPARS)
	  {
	    fineTunings = doCalcFineTuning<SoftParsMssm>(s, pMSSMftBCs, tuningPars, vevs, mx, sing, 
							 doCalcLHSTuningMatrix, doCalcRHSTuningVector, doCalcdLogMzSqdLogParam);
	    fineTuningsApprox = doCalcFineTuning<SoftParsMssm>(s, pMSSMftBCs, tuningPars, vevs, highScaleCouplings, mx, sing,
							       doCalcLHSTuningMatrix, 
							       doCalcRHSTuningVector_pMSSM_Approx, doCalcdLogMzSqdLogParam);
	  }
	else
	  {
	    fineTunings = doCalcFineTuning<SoftParsMssm>(s, mSUGRAftBCs, tuningPars, vevs, mx, sing, 
							 doCalcLHSTuningMatrix, doCalcRHSTuningVector, doCalcdLogMzSqdLogParam);
	  }


	//cout << fineTunings;
        //cout << "Fine tuning = " << fineTunings.max() << endl;

	//cout << fineTuningsApprox;
        //cout << "Fine tuning (approx) = " << fineTuningsApprox.max() << endl;


	// // if (tuningPars.displayEnd() == NUMPARSPMSSM) // pMSSM case
	// //   {
	// //     cout << "Fine tuning due to M3 = " << fineTunings(3) << endl;
	// //   }
	End commented out 5/5/2014 */

//	DoubleVector softsusyTunings = s.fineTune(boundaryCondition, pars, mx, false);

	//cout << softsusyTunings;
	/* Commented out here 5/5/2014
	   if (sing != 0)
	   {
	   cerr << "WARNING: problem with fine tuning calculation: tunings should not be trusted." << endl;
	   problemFlag = true;
	   }

	   // If there are no problems, print the result
	   if (!r.displayProblem().test() && problemFlag == false) 
	   {
	   if (tuningPars.displayEnd() == NUMPMSSMPARS)
	   {
	   // Results printed out depend on parameter scanner.
	   // For M_3, print value of M_3, value of physical gluino mass,
	   // fine tuning, lightest Higgs mass
	   
	   if (scan_par == 3)
	   {
		    cout << pars(scan_par) << " " 
			 << r.displayPhys().mGluino << " "
			 << r.displayPhys().mh0 << " " 
			 << fineTunings.max() << " " << fineTuningsApprox.max() << " " << fineTuningsNumerical.max() << " "
			 << fineTunings(3) << " " << fineTuningsApprox(3) << " " << fineTuningsNumerical(3) << " " 
			 << softsusyTunings.max() << endl;
		  }
		else if (scan_par == 43 || scan_par == 46)
		  {
		    cout << pars(scan_par) << " " 
			 << r.displayPhys().mu(1,3) << " "
			 << r.displayPhys().mu(2,3) << " "
			 << r.displayPhys().mh0 << " "
			 << fineTunings.max() << " " << fineTuningsApprox.max() << " ";
		    if (scan_par == 43)
		      {
			cout << fineTunings(19) << " " << fineTuningsApprox(19) << " ";
		      }
		    else
		      {
			cout << fineTunings(22) << " " << fineTuningsApprox(22) << " ";
		      }
		    cout << softsusyTunings.max() << endl;
		  }
		else
		  {
		    cout << pars(scan_par) << " " 
			 << r.displayPhys().mh0 << " " 
			 << r.displayPhys().mA0 << " " 
			 << r.displayPhys().mH0 << " " 
			 << r.displayPhys().mHpm << " "
			 << fineTunings.max() << " " << fineTuningsApprox.max() << " " << fineTuningsApprox.max() << " " << fineTuningsNumerical.max() << " "<< softsusyTunings.max() << endl;
		  }
	      }
	    else
	      {
		cout << pars(scan_par)<< " " << r.displayPhys().mh0 << " " << fineTunings.max() << " " << fineTuningsApprox.max() << " " << fineTuningsNumerical.max() << " " << softsusyTunings.max() << endl;
	      }
	  }
	else if (r.displayProblem().test() && problemFlag == false) 
	  {
	    /// print out what the problem(s) is(are)
	    cerr << pars(scan_par) << " " << r.displayProblem() << endl;
	    /// in out stream might also print flag to indicate problem? 
	    }
	  else if (problemFlag == true)
	{
	cerr << "WARNING: problem with fine tuning calculation." << endl;
	cerr << pars(scan_par) << " "
	<< fineTunings.max() << " " << fineTuningsApprox.max() << " " << fineTuningsNumerical.max() << " "<< softsusyTunings.max() << endl;
	}
	}
	End commented out here */
