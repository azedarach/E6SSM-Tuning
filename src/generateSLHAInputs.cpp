/*
  generateSLHAInputs produces the SLHA input files for the scan over the
  E6SSM and pMSSM parameter space. More convenient than the more general
  code generateE6SSMSLHAInputs, which is still a work in progress. Reads in
  file containing the location to save the input files to, a flag indicating
  if the scan is to be random or a grid scan, and the bounds for the parameters.
  Note that s is regarded as being fixed in each scan, and is provided as an
  input in the file.
  Currently scans can be done over the parameters
     - tan(beta)(MX) (EXTPAR 25) - this is preferable for matching the E6 and MSSM
       models at each scale
     - lambda
     - A_lambda
     - M_1
     - M_2
     - M_3
     - M_1'
     - m_Q3^2
     - m_u3^2
     - m_d3^2
     - m_L3^2
     - m_e3^2
     - A_t
     - A_b
     - A_tau
  which covers essentially all of the parameters found to have a non-negligible 
  contribution to fine tuning in hep-ph/1206.5800. Additional parameters
  may always be added later.
  The following parameters have now also been added to the scan:
     - m_Q1^2
     - m_u1^2
     - m_d1^2
     - m_L1^2
     - m_e1^2
     - m_Q2^2
     - m_u2^2
     - m_d2^2
     - m_L2^2
     - m_e2^2
     - kappa(3)
     - A_kappa(3)
  As in the pMSSM fine tuning study, the first and second generation masses
  can be set to be degenerate (so there are only five parameters in this sector).
  We continue to assume that all first and second generation Yukawa couplings 
  are zero, along with the associated A terms. 
  Author: Dylan Harries
  Date created: 26/3/2014
 */

#include "essmTuningMeasures.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <stdio.h>
#include <stdlib.h> 
#include <time.h>  
#include <math.h>

void errorCall();
string ToUpper(const string &) ;
void writeSLHAInputFile(ofstream &, const string &, QedQcd const &, double, 
			bool, double, double, DoubleVector const &, DoubleVector const &, double, int, bool);
void slhaInputFileHeader(ofstream &);
void slhaInputModselBlock(ofstream &, const string &);
void slhaInputSminputBlock(ofstream &, QedQcd const &, double);
void slhaInputMinparBlock(ofstream &, double);
void slhaInputExtparBlock(ofstream &, bool, double, DoubleVector const &, DoubleVector const &);
void slhaInputSoftSusyBlock(ofstream &, double, int, bool);
double uniformRan(long &);

using namespace std;

int main(int argc, char* argv[])
{

  // Default Standard Model inputs.
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

  bool flag = false;

  // Default SOFTSUSY options
  double tol = 1.0e-3;
  int mix = 0;
  bool include2lpScalars = true;

  // Default parameter lower bounds to be used if not supplied
  double tb_lower = 1.0;
  double lambda_lower = 0.0;
  double Alambda_lower = 0.0; // GeV
  double M1_lower = 0.0; // GeV
  double M2_lower = 0.0; // GeV
  double M3_lower = 0.0; // GeV
  double M1p_lower = 0.0; // GeV  
  double mQ3_lower = 0.0; // GeV
  double mu3_lower = 0.0; // GeV
  double md3_lower = 0.0; // GeV
  double mL3_lower = 0.0; // GeV
  double me3_lower = 0.0; // GeV
  double At_lower = 0.0; // GeV
  double Ab_lower = 0.0; // GeV
  double Atau_lower = 0.0; // GeV 
  // Add new parameter lower bounds here...
  double kappa_lower = 0.0;
  double Akappa_lower = 0.0; // GeV
  double mQ1_lower = 0.0; // GeV
  double mu1_lower = 0.0; // GeV
  double md1_lower = 0.0; // GeV
  double mL1_lower = 0.0; // GeV
  double me1_lower = 0.0; // GeV
  double mQ2_lower = 0.0; // GeV
  double mu2_lower = 0.0; // GeV
  double md2_lower = 0.0; // GeV
  double mL2_lower = 0.0; // GeV
  double me2_lower = 0.0; // GeV

  // Default parameter upper bounds to be used if not supplied
  double tb_upper = 1.0;
  double lambda_upper = 0.0;
  double Alambda_upper = 0.0; // GeV
  double M1_upper = 0.0; // GeV
  double M2_upper = 0.0; // GeV
  double M3_upper = 0.0; // GeV
  double M1p_upper = 0.0; // GeV  
  double mQ3_upper = 0.0; // GeV
  double mu3_upper = 0.0; // GeV
  double md3_upper = 0.0; // GeV
  double mL3_upper = 0.0; // GeV
  double me3_upper = 0.0; // GeV
  double At_upper = 0.0; // GeV
  double Ab_upper = 0.0; // GeV
  double Atau_upper = 0.0; // GeV 
  // Add new parameter upper bounds here...
  double kappa_upper = 0.0;
  double Akappa_upper = 0.0; // GeV
  double mQ1_upper = 0.0; // GeV
  double mu1_upper = 0.0; // GeV
  double md1_upper = 0.0; // GeV
  double mL1_upper = 0.0; // GeV
  double me1_upper = 0.0; // GeV
  double mQ2_upper = 0.0; // GeV
  double mu2_upper = 0.0; // GeV
  double md2_upper = 0.0; // GeV
  double mL2_upper = 0.0; // GeV
  double me2_upper = 0.0; // GeV

  // Default number of points to use in grid scan if not provided
  int tb_npts = 1;
  int lambda_npts = 1;
  int Alambda_npts = 1;
  int M1_npts = 1;
  int M2_npts = 1;
  int M3_npts = 1;
  int M1p_npts = 1;
  int mQ3_npts = 1;
  int mu3_npts = 1;
  int md3_npts = 1;
  int mL3_npts = 1;
  int me3_npts = 1;
  int At_npts = 1;
  int Ab_npts = 1;
  int Atau_npts = 1;
  // Add new default numbers of parameters here...
  int kappa_npts = 1;
  int Akappa_npts = 1;
  int mQ1_npts = 1;
  int mu1_npts = 1;
  int md1_npts = 1;
  int mL1_npts = 1;
  int me1_npts = 1;
  int mQ2_npts = 1;
  int mu2_npts = 1;
  int md2_npts = 1;
  int mL2_npts = 1;
  int me2_npts = 1;

  // Model types that can be produced
  string modelType = "E6SSM";
  string altModelType = "MSSM";

  // Default values for scan options. Note that E6 files are always produced;
  // the option you have is to also produce analogous MSSM points.
  bool isMSSM = false;
  bool hasModelRequest = false;
  bool isRandom = false;
  bool hasScanType = false;
  bool hasNumRandPoints = false;
  int randPoints = 1; // default number of models to produce if doing random points
  bool isFirstEqualsSecond = true; // degenerate first and second generation, default is true

  // Default value of the singlet VEV s. Together with lambda this fixes mu_eff.
  double s = 5500.0; // GeV
  bool hasVEV = false;

  // Default value for input scale MX values.
  bool hasMX = false;
  bool useMXequalsMSUSY = false;
  double MX = 20000.0; // GeV

  // Variable to store path to save files in
  string filePrefix;
  // Note that is no file path is provided, files produced will be saved
  // to the current folder. A warning should be produced indicating this.
  bool hasFilePath = false; 

  // Get seed for random number generator
  time_t timer = time(NULL);
  long crtTime = (long) timer;
  // Seed generator using current time.
  // WARNING: this could fail if the program
  // is called rapidly in quick succession (need
  // at least one second between successive calls)
  long idum = -1.0*crtTime; // random seed

  // Vectors for storing parameter numbers, bounds, scan types and number of points.
  // These will in general change length on each read as new things are added, which
  // could slow the code down. However as this is only a small part of the code and
  // the number of parameters that are expected to appear is small this is ok for now.
  DoubleVector extParNums(1); // selected parameters
  DoubleVector extLowBnds(1); // lower bounds for parameters
  DoubleVector extUpBnds(1); // upper bounds for parameters
  DoubleVector extNParVals(1); // number of values of each parameter
  
  int numExtPars = 0;
  bool isFirstExtPar = true;
  bool isReadingFirstExtPar = false;
  bool isNewExtParBlock = false;
  bool hasExtParBlock = false;

  bool hasSMInputs = false;

  // Attempt to read in specs. file
  try 
    {
      if (argc > 1)
  	{
  	  errorCall();
  	}
      else
  	{
  	  // # of decimal places to use in output. 8 is SLHA convention.
  	  outputCharacteristics(8);
	  
  	  // Read input parameters. Basically identical procedure to that used in
  	  // SOFTSUSY's softpoint.cpp. 
  	  string line, block;
  	  int model; 


  	  // NB as currently written this will read the last line of the input file
  	  // twice. Not exactly desirable, but also not too harmful for now.
  	  while (getline(cin, line))
  	    {
  	      istringstream instream(line);
  	      string word1, word2;
  	      instream >> word1;

  	      if (word1.find("#") == string::npos)
  		{
  		  // There's no comment, so read in the next word.
  		  instream >> word2;
  		  if (ToUpper(word1) == "BLOCK")
  		    {
  		      block = ToUpper(word2);
		      if (block == "EXTPARSEL")
			{
			  isNewExtParBlock = true;
			  if (isFirstExtPar == true && isReadingFirstExtPar == false)
			    {
			      isReadingFirstExtPar = true;
			      hasExtParBlock = true;
			    }
			  // You only get here once you've read the first EXTPAR block
			  else if (isFirstExtPar == true && isReadingFirstExtPar == true) 
			    {
			      isFirstExtPar = false;
			      isReadingFirstExtPar = false;
			    }
			}
  		    }
  		  else
  		    {
  		      // It's not a comment or a block line, so it should
  		      // be input data.
  		      istringstream datastream(line);

		      if (block == "FILENAMES")
			{
			  int i; string name; datastream >> i >> name;
			  if (i == 1)
			    {
			      filePrefix = name;
			      hasFilePath = true;
			    }
			  else
			    {
			      cout << "# WARNING: Don't understand data input " << i 
				   << " " << name << " in block "
				   << block << ": ignoring it\n"; break;
			    }
			}
		      else if (block == "OPTIONS")
			{
			  int i, d; datastream >> i >> d;
			  if (i == 3)
			    {
			      if (d >= 0)
				{
				  hasNumRandPoints = true;
				  randPoints = d;
				}
			      else
				{
				  cout << "# WARNING: requested number of random points should be a "
				       << "positive integer: using default " << randPoints << ".\n";
				}
			    }
			  else if (i == 1)
			    {
			      switch(d)
				{
				case 1: isMSSM = false; hasModelRequest = true; break;
				case 2: isMSSM = true; hasModelRequest = true; break;
				default : 
				  cout << "# WARNING: Don't understand data input " << i 
				       << " " << d << " in block "
				       << block << ": ignoring it\n"; break;
				}

			    }
			  else if (i == 2)
			    {
			      switch(d)
				{
				case 1: isRandom = false; hasScanType = true; break;
				case 2: isRandom = true; hasScanType = true; break;
				default :
				  cout << "# WARNING: Don't understand data input " << i
				       << " " << d << " in block "
				       << block << ": ignoring it\n"; break;
				}
			    }
			  else if (i == 4)
			    {
			      switch(d)
				{
				case 1: isFirstEqualsSecond = true; break;
				case 2: isFirstEqualsSecond = false; break;
				default :
				  cout << "# WARNING: Don't understand data input " << i
				       << " " << d << " in block "
				       << block << ": ignoring it\n"; break;
				}
			    }
			  else
			    {
			      cout << "# WARNING: Don't understand data input " << i 
				   << " " << d << " in block "
				   << block << ": ignoring it\n"; 
			    }
			}
		      else if (block == "SOFTSUSY")
			{
			  int i; double d; datastream >> i >> d;
			  switch(i)
			    {
			    case 1: tol = d; break;
			    case 2: mix = (int) d; break;
			    case 5: include2lpScalars = bool(int(d+EPSTOL)); break;
			    default :
			      cout << "# WARNING: Don't understand data input " << i
				   << " " << d << " in block "
				   << block << ": ignoring it\n"; break;
			    }
			}
		      else if (block == "MXVAL")
			{
			  int i; double d; datastream >> i >> d;
			  if (i == 1)
			    {
			      hasMX = true;
			      if (d < 0.0)
				{
				  useMXequalsMSUSY = true;
				  MX = -1;
				}
			      else
				{
				  MX = d;
				}
			    }
			  else
			    {
			      cout << "# WARNING: Don't understand data input " << i 
				   << " " << d << " in block "
				   << block << ": ignoring it\n"; 
			    }
			}
		      else if (block == "SVAL")
			{
			  int i; double d; datastream >> i >> d;

			  if (i == 1)
			    {
			      s = d;
			      hasVEV = true;
			    }
			  else
			    {
			      cout << "# WARNING: Don't understand data input " << i 
				   << " " << d << " in block "
				   << block << ": ignoring it\n"; 
			    }
			}
  		      else if (block == "EXTPARSEL")
  			{
  			  if (!strcmp(modelType.c_str(), "E6SSM"))
  			    {
  			      int i; double d; datastream >> i >> d;

			      if (isReadingFirstExtPar)
				{
				  
				  // Read in values
				  switch(i) 
				    {
				    case 1: extParNums(1) = d; break;
				    case 2: extLowBnds(1) = d; break;
				    case 3: extUpBnds(1) = d; break;
				    case 4: extNParVals(1) = d; break;
				    default:
				      cout << "# WARNING: Don't understand data input " << i 
					   << " " << d << " in block "
					   << block << ": ignoring it\n"; break;
				    }
				}
			      else
				{
				  if (isNewExtParBlock)
				    {
				      // Increment number of parameters
				      extParNums.setEnd(extParNums.displayEnd()+1);
				      extLowBnds.setEnd(extLowBnds.displayEnd()+1);
				      extUpBnds.setEnd(extUpBnds.displayEnd()+1);
				      extNParVals.setEnd(extNParVals.displayEnd()+1);

				      // Set to -1 to avoid confusion with valid MX entry
				      extParNums(extParNums.displayEnd()) = -1.0;

				      isNewExtParBlock = false;
				    }
				  // Read in values
				  switch(i) 
				    {
				    case 1: extParNums(extParNums.displayEnd()) = d; break;
				    case 2: extLowBnds(extLowBnds.displayEnd()) = d; break;
				    case 3: extUpBnds(extUpBnds.displayEnd()) = d; break;
				    case 4: extNParVals(extNParVals.displayEnd()) = d; break;
				    default:
				      cout << "# WARNING: Don't understand data input " << i 
					   << " " << d << " in block "
					   << block << ": ignoring it\n"; break;
				    }
				}
			    }
  			}
  		      else if (block == "SMINPUTS")
  			{
			  hasSMInputs = true;
  			  int i; double d; datastream >> i >> d;
  			  switch (i) {
  			  case 1: oneset.setAlpha(ALPHA, 1.0 / d); break;
  			  case 2: G_F = d; break;
  			  case 3: oneset.setAlpha(ALPHAS, d); break; 
  			  case 4: oneset.setMu(d); MZ = d; break;
  			  case 5: oneset.setMass(mBottom, d); flag = true; oneset.setMbMb(d); break;
  			  case 6: oneset.setPoleMt(d); break;
  			  case 7: oneset.setMass(mTau, d); oneset.setPoleMtau(d); break;
  			  case 11: oneset.setMass(mElectron, d);break;
  			  case 13: oneset.setMass(mMuon, d);break;
  			  case 21: oneset.setMass(mDown, d); break;
  			  case 22: oneset.setMass(mUp, d);  break;
  			  case 23: oneset.setMass(mStrange, d); break;
  			  case 24: oneset.setMass(mCharm, d); break;
  			  default: 
  			    cout << "# WARNING: Don't understand data input " << i 
  				 << " " << d << " in block "
  				 << block << ": ignoring it\n"; break;
  			  }
  			}
  		      else
  			{
  			  cout << "# WARNING: cannot recognise block " << block 
  			       << ": ignoring all data in it." << endl;
  			} // end of block checks
  		    } // end of data read in
  		} // end of not comment line
  	    } // end of file

	  // Error checking on whether a file path has been provided
	  if (hasFilePath == false)
	    {
	      cout << "# WARNING: No file path provided: model files will be saved to "
		   << "current directory.\n";
	      filePrefix = "";
	    }

	  // Error checking on options that are provided - mainly need to check if scan type
	  // is provided, and if a random scan is requested whether the number of points to
	  // produce has been provided.
	  if (hasModelRequest == false)
	    {
	      cout << "# WARNING: No model types requested. Only writing E6 models.\n";
	    }

	  if (hasScanType == false)
	    {
	      cout << "# WARNING: No scan type provided: doing default grid scan.\n";
	    }

	  if (isRandom)
	    {
	      if (hasNumRandPoints == false)
		{
		  cout << "# WARNING: Random scan requested but number of models not provided:\n"
		       << "#          producing default " << randPoints << " models.\n";
		}
	    }

	  // Error checking on s value that is provided
	  if (hasVEV == false)
	    {
	      cout << "# WARNING: No value for singlet VEV s provided: " 
		   << "using default s = " << s << " GeV.\n";
	    }

	  if (hasSMInputs == false)
	    {
	      cout << "# WARNING: No Standard Model inputs provided: using defaults." << endl;
	    }

	  // Error checking on whether an input scale is provided
	  if (hasMX == false)
	    {
	      cout << "# WARNING: No input scale provided: spectrum generator will use\n"
		   << "#          default MX = MGUT.\n";
	    }

	  // Error checking for lower and upper bounds for each parameter.
	  // For each parameter provided, we want to check that it is a valid
	  // E_6SSM parameter and that a sensible number of points has been
	  // requested, if that is the mode we are operating in. Lower bounds
	  // should also be less than upper bounds, and if this is not the case
	  // they are swapped.
	  numExtPars = extParNums.displayEnd()-extParNums.displayStart()+1;

	  DoubleVector validExtParNums(1), validExtNParVals(1), validExtLowBnds(1), validExtUpBnds(1);
	  int crtParNum, crtNParVals;
	  bool isValidPar;
	  bool isFirstValidPar = true;
	  bool hasValidExtPar = false;

	  if (hasExtParBlock)
	    {
	      for (int i = 1; i <= numExtPars; i++)
		{
		  isValidPar = true;
		  crtParNum = (int) extParNums(i);
		  // Compare against list of currently allowed SLHA parameter numbers to determine if valid
		  // parameter. Filter out invalid ones.
		  switch(crtParNum)
		    {
		      // This is an ugly way of checking, but it works.
		      // Note that \mu_MSSM (EXTPAR 23) is valid only when we are producing
		      // only MSSM model files. If E6 models are also to be produced, \mu is
		      // illegal and you must specify EXTPAR 65 mu_eff instead.
		    case 1: 
		    case 2: 
		    case 3: 
		    case 4: 
		    case 11: 
		    case 12: 
		    case 13: 
		    case 25: 
		    case 31:
		    case 32:
		    case 33:
		    case 34:
		    case 35: 
		    case 36:
		    case 41:
		    case 42: 
		    case 43:
		    case 44:
		    case 45: 
		    case 46:
		    case 47:
		    case 48: 
		    case 49: 
		    case 61:
		    case 62:
		    case 63: 
		    case 64:  break;
		    default:
		      cout << "# WARNING: Cannot currently recognise EXTPAR " << crtParNum
			   << ": ignoring it.\n"; isValidPar = false; break;
		    }

		  if (isValidPar)
		    { 
		      if (hasValidExtPar == false)
			{
			  hasValidExtPar = true;
			}
		      if (isRandom == false)
			{
			  // Check that number of points requested is greater than 0 if grid scanning. If not, print a
			  // warning and set to 1 point to evaluate.
			  crtNParVals = (int) extNParVals(i);
			  
			  if (crtNParVals <= 0)
			    {
			      cout << "# WARNING: Number of points for EXTPAR " << crtParNum<< " must be a positive integer:\n"
				   << "#          defaulting to 1 point.\n";
			      crtNParVals = 1;
			    }
			}
		      

		      // Check ordering of upper and lower bounds.
		      if (extUpBnds(i) < extLowBnds(i))
			{
			  double exttemp = extUpBnds(i);
			  extUpBnds(i) = extLowBnds(i);
			  extLowBnds(i) = exttemp;
			}
			
		      // Add to list of valid parameter specifications
		      if (isFirstValidPar)
			{
			  validExtParNums(1) = crtParNum;
			  validExtNParVals(1) = crtNParVals;
			  validExtLowBnds(1) = extLowBnds(i);
			  validExtUpBnds(1) = extUpBnds(i);
			  isFirstValidPar = false;
			}
		      else
			{
			  validExtParNums.setEnd(validExtParNums.displayEnd()+1);
			  validExtNParVals.setEnd(validExtNParVals.displayEnd()+1);
			  validExtLowBnds.setEnd(validExtLowBnds.displayEnd()+1);
			  validExtUpBnds.setEnd(validExtUpBnds.displayEnd()+1);

			  validExtParNums(validExtParNums.displayEnd()) = crtParNum;
			  validExtNParVals(validExtNParVals.displayEnd()) = crtNParVals;
			  validExtLowBnds(validExtLowBnds.displayEnd()) = extLowBnds(i);
			  validExtUpBnds(validExtUpBnds.displayEnd()) = extUpBnds(i);

			}
		    }
		  
		}
	    }
	  else
	    {
	      cout << "# WARNING: No EXTPAR blocks provided: default values will be used in spectrum generation.\n";
	    }


	  // Sort valid EXTPAR parameters into numerical order for neatness when writing
	  int numValidExtPars = validExtParNums.displayEnd()-validExtParNums.displayStart()+1;
	  DoubleVector parNums = validExtParNums.sort();
	  DoubleVector nParVals(numValidExtPars), lowBnds(numValidExtPars), upBnds(numValidExtPars);
	  int pos;
	  for (int k = 1; k <= numValidExtPars; k++)
	    {
	      pos  = validExtParNums.closest(parNums(k));
	      nParVals(k) = validExtNParVals(pos);
	      lowBnds(k) = validExtLowBnds(pos);
	      upBnds(k) = validExtUpBnds(pos);
	    }

	  // Then do summary print of requested EXTPAR parameters
	  cout << "# SUMMARY: " << endl;
	  cout << "# Requested singlet VEV s = " << s << " GeV" <<  endl;
	  if (hasMX == false)
	    {
	      cout << "# Assuming MX = M_GUT" << endl;
	    }
	  else
	    {
	      if (useMXequalsMSUSY)
		{
		  cout << "# Using MX = M_SUSY" << endl;
		}
	      else
		{
		  cout << "# Using MX = " << MX << endl;
		}
	    }
	  if (hasExtParBlock == true && hasValidExtPar == true)
	    {
	      for (int p = 1; p <= numValidExtPars; p++)
		{
		  cout << "# Requested EXTPAR " << (int) parNums(p) << ", with values: " << endl;
		  cout << "#           Lower bound = " << lowBnds(p) << endl;
		  cout << "#           Upper bound = " << upBnds(p) << endl;  
		  if (isRandom == false)
		    {
		      cout << "#           Number of values = " << (int) nParVals(p) << endl; 
		    }
		}
	    }
	  else if (hasExtParBlock == true && hasValidExtPar == false)
	    {
	      cout << "# WARNING: Could not find any valid EXTPAR requests:"
		   << " default values will be used in spectrum generation.\n";
	    }

	  // Either display the fixed number of models requested, or calculate the number that
	  // will be produced
	  if (isRandom)
	    {
	      cout << "# Requested uniformly distributed random points." << endl;
	      cout << "# Using " << randPoints << " random points for randomly scanned parameters." << endl;
	    }
	  else
	    {
	      cout << "# Requested grid scan of input parameters." << endl;
	    }

	  if (isFirstEqualsSecond)
	    {
	      cout << "# Requested degenerate first and second generation soft masses." << endl;
	    }
	  else
	    {
	      cout << "# Not requiring degenerate first and second generation soft masses." << endl;
	    }

	  cout << "# Writing model input files with file name structure: " << endl;
	  cout << "# " << filePrefix << "essm_[model #]" << endl;
	  if (isMSSM)
	    {
	      cout << "# " << filePrefix << "mssm_[model #]" << endl;
	    }


	  // Now update default values of allowed parameters with any values that
	  // were provided.
	  /*   - tan(beta)(MX) (EXTPAR 25) - this is preferable for matching the E6 and MSSM
	       models at each scale
	       - lambda = 61
	       - A_lambda = 63
	       - M_1 = 1
	       - M_2 = 2
	       - M_3 = 3
	       - M_1' = 4
	       - m_Q3^2 = 43 
	       - m_u3^2 = 46
	       - m_d3^2 = 49
	       - m_L3^2 = 33
	       - m_e3^2 = 36
	       - A_t = 11
	       - A_b = 12
	       - A_tau = 13
               - kappa = 62
               - A_kappa = 64
	       - m_Q1^2 = 41 
	       - m_u1^2 = 44
	       - m_d1^2 = 47
	       - m_L1^2 = 31
	       - m_e1^2 = 34
	       - m_Q2^2 = 42 
	       - m_u2^2 = 45
	       - m_d2^2 = 48
	       - m_L2^2 = 32
	       - m_e2^2 = 35
	  */

	  // We also need to check that, if degenerate soft masses are requested, values
	  // are used from only one request, by default the first generation.
	  bool hasmQ1Vals = false;
	  bool hasmQ2Vals = false;
	  bool hasmu1Vals = false;
	  bool hasmu2Vals = false;
	  bool hasmd1Vals = false;
	  bool hasmd2Vals = false;
	  bool hasmL1Vals = false;
	  bool hasmL2Vals = false;
	  bool hasme1Vals = false;
	  bool hasme2Vals = false;


	  for (int i = 1; i <= numValidExtPars; i++)
	    {
	      crtParNum = (int) parNums(i);
	      switch(crtParNum)
		{
		case 1:
		  {
		    M1_lower = lowBnds(i);
		    M1_upper = upBnds(i);
		    M1_npts = (int) nParVals(i);
		    break;
		  }
		case 2:
		  {
		    M2_lower = lowBnds(i);
		    M2_upper = upBnds(i);
		    M2_npts = (int) nParVals(i);
		    break;
		  }
		case 3:
		  {
		    M3_lower = lowBnds(i);
		    M3_upper = upBnds(i);
		    M3_npts = (int) nParVals(i);
		    break;
		  }
		case 4:
		  {
		    M1p_lower = lowBnds(i);
		    M1p_upper = upBnds(i);
		    M1p_npts = (int) nParVals(i);
		    break;
		  } 
		case 11:
		  {
		    At_lower = lowBnds(i);
		    At_upper = upBnds(i);
		    At_npts = (int) nParVals(i);
		    break;
		  } 
		case 12:
		  {
		    Ab_lower = lowBnds(i);
		    Ab_upper = upBnds(i);
		    Ab_npts = (int) nParVals(i);
		    break;
		  } 
		case 13: 
		  {
		    Atau_lower = lowBnds(i);
		    Atau_upper = upBnds(i);
		    Atau_npts = (int) nParVals(i);
		    break;
		  }
		case 25:
		  {
		    tb_lower = lowBnds(i);
		    tb_upper = upBnds(i);
		    tb_npts = (int) nParVals(i);
		    break;
		  }
		case 31: 
		  {
		    mL1_lower = lowBnds(i);
		    mL1_upper = upBnds(i);
		    mL1_npts = (int) nParVals(i);
		    hasmL1Vals = true;	  
		    break;
		  }
		case 32: 
		  {
		    mL2_lower = lowBnds(i);
		    mL2_upper = upBnds(i);
		    mL2_npts = (int) nParVals(i);
		    hasmL2Vals = true;	  
		    break;
		  }
		case 33: 
		  {
		    mL3_lower = lowBnds(i);
		    mL3_upper = upBnds(i);
		    mL3_npts = (int) nParVals(i);
		    break;
		  }
		case 34:
		  {
		    me1_lower = lowBnds(i);
		    me1_upper = upBnds(i);
		    me1_npts = (int) nParVals(i);
		    hasme1Vals = true;	  
		    break;
		  }
		case 35:
		  {
		    me2_lower = lowBnds(i);
		    me2_upper = upBnds(i);
		    me2_npts = (int) nParVals(i);
		    hasme2Vals = true;	  
		    break;
		  }
		case 36:
		  {
		    me3_lower = lowBnds(i);
		    me3_upper = upBnds(i);
		    me3_npts = (int) nParVals(i);
		    break;
		  }
		case 41:
		  {
		    mQ1_lower = lowBnds(i);
		    mQ1_upper = upBnds(i);
		    mQ1_npts = (int) nParVals(i);
		    hasmQ1Vals = true;	  
		    break;
		  }
		case 42:
		  {
		    mQ2_lower = lowBnds(i);
		    mQ2_upper = upBnds(i);
		    mQ2_npts = (int) nParVals(i);
		    hasmQ2Vals = true;	  
		    break;
		  }
		case 43:
		  {
		    mQ3_lower = lowBnds(i);
		    mQ3_upper = upBnds(i);
		    mQ3_npts = (int) nParVals(i);
		    break;
		  }
		case 44:
		  {
		    mu1_lower = lowBnds(i);
		    mu1_upper = upBnds(i);
		    mu1_npts = (int) nParVals(i);
		    hasmu1Vals = true;	  
		    break;
		  }
		case 45:
		  {
		    mu2_lower = lowBnds(i);
		    mu2_upper = upBnds(i);
		    mu2_npts = (int) nParVals(i);
		    hasmu2Vals = true;	  
		    break;
		  }
		case 46:
		  {
		    mu3_lower = lowBnds(i);
		    mu3_upper = upBnds(i);
		    mu3_npts = (int) nParVals(i);
		    break;
		  }
		case 47:
		  {
		    md1_lower = lowBnds(i);
		    md1_upper = upBnds(i);
		    md1_npts = (int) nParVals(i);
		    hasmd1Vals = true;	  
		    break;
		  }
		case 48:
		  {
		    md2_lower = lowBnds(i);
		    md2_upper = upBnds(i);
		    md2_npts = (int) nParVals(i);
		    hasmd2Vals = true;	  
		    break;
		  }
		case 49:
		  {
		    md3_lower = lowBnds(i);
		    md3_upper = upBnds(i);
		    md3_npts = (int) nParVals(i);
		    break;
		  }
		case 61:
		  {
		    lambda_lower = lowBnds(i);
		    lambda_upper = upBnds(i);
		    lambda_npts = (int) nParVals(i);
		    break;
		  }
		case 62:
		  {
		    kappa_lower = lowBnds(i);
		    kappa_upper = upBnds(i);
		    kappa_npts = (int) nParVals(i);
		    break;
		  }
		case 63:
		  {
		    Alambda_lower = lowBnds(i);
		    Alambda_upper = upBnds(i);
		    Alambda_npts = (int) nParVals(i);
		    break;
		  }
		case 64:
		  {
		    Akappa_lower = lowBnds(i);
		    Akappa_upper = upBnds(i);
		    Akappa_npts = (int) nParVals(i);
		    break;
		  }
		default :
		  cout << "# WARNING: encountered unexpected EXTPAR code " << crtParNum <<": ignoring it.\n"; break;
		}
	    }

	  // Check which parameters have been specified from the first and second generation, and
	  // select only one from each if we are keeping the two generations degenerate.
	  if (isFirstEqualsSecond)
	    {
	      if (hasmL1Vals == true && hasmL2Vals == false)
		{
		  mL2_lower = mL1_lower;
		  mL2_upper = mL1_upper;
		  mL2_npts = mL1_npts;
		}
	      else if (hasmL1Vals == false && hasmL2Vals == true)
		{
		  mL1_lower = mL2_lower;
		  mL1_upper = mL2_upper;
		  mL1_npts = mL2_npts;
		}	  
	      else if (hasmL1Vals == true && hasmL2Vals == true)
		{
		  mL2_lower = mL1_lower;
		  mL2_upper = mL1_upper;
		  mL2_npts = mL1_npts;
		}

	      if (hasme1Vals == true && hasme2Vals == false)
		{
		  me2_lower = me1_lower;
		  me2_upper = me1_upper;
		  me2_npts = me1_npts;
		}
	      else if (hasme1Vals == false && hasme2Vals == true)
		{
		  me1_lower = me2_lower;
		  me1_upper = me2_upper;
		  me1_npts = me2_npts;
		}	  
	      else if (hasme1Vals == true && hasme2Vals == true)
		{
		  me2_lower = me1_lower;
		  me2_upper = me1_upper;
		  me2_npts = me1_npts;
		}

	      if (hasmQ1Vals == true && hasmQ2Vals == false)
		{
		  mQ2_lower = mQ1_lower;
		  mQ2_upper = mQ1_upper;
		  mQ2_npts = mQ1_npts;
		}
	      else if (hasmQ1Vals == false && hasmQ2Vals == true)
		{
		  mQ1_lower = mQ2_lower;
		  mQ1_upper = mQ2_upper;
		  mQ1_npts = mQ2_npts;
		}	  
	      else if (hasmQ1Vals == true && hasmQ2Vals == true)
		{
		  mQ2_lower = mQ1_lower;
		  mQ2_upper = mQ1_upper;
		  mQ2_npts = mQ1_npts;
		}

	      if (hasmu1Vals == true && hasmu2Vals == false)
		{
		  mu2_lower = mu1_lower;
		  mu2_upper = mu1_upper;
		  mu2_npts = mu1_npts;
		}
	      else if (hasmu1Vals == false && hasmu2Vals == true)
		{
		  mu1_lower = mu2_lower;
		  mu1_upper = mu2_upper;
		  mu1_npts = mu2_npts;
		}	  
	      else if (hasmu1Vals == true && hasmu2Vals == true)
		{
		  mu2_lower = mu1_lower;
		  mu2_upper = mu1_upper;
		  mu2_npts = mu1_npts;
		}

	      if (hasmd1Vals == true && hasmd2Vals == false)
		{
		  md2_lower = md1_lower;
		  md2_upper = md1_upper;
		  md2_npts = md1_npts;
		}
	      else if (hasmd1Vals == false && hasmd2Vals == true)
		{
		  md1_lower = md2_lower;
		  md1_upper = md2_upper;
		  md1_npts = md2_npts;
		}	  
	      else if (hasmd1Vals == true && hasmd2Vals == true)
		{
		  md2_lower = md1_lower;
		  md2_upper = md1_upper;
		  md2_npts = md1_npts;
		}
	    }

	  // Define variable for MINPAR tan(beta), which we will just set equal
	  // to the EXTPAR value (seeing as it will be ignored anyway)
	  double minpar_tb;

	  // Define file stream for writing
	  int modelCount = 0; // current model number
	  string fileName = filePrefix;

	  // Vector to store current variable values and a vector
	  // to store their SLHA codes
	  DoubleVector crtParVals(28); // 28 = # currently implemented pars. + mu_eff
	  DoubleVector parCodes(28);
	  parCodes(1) = 1;
	  parCodes(2) = 2;
	  parCodes(3) = 3;
	  parCodes(4) = 4;
	  parCodes(5) = 11;
	  parCodes(6) = 12;
	  parCodes(7) = 13;
	  parCodes(8) = 25;
	  parCodes(9) = 31;
	  parCodes(10) = 32;
	  parCodes(11) = 33;
	  parCodes(12) = 34;
	  parCodes(13) = 35;
	  parCodes(14) = 36;
	  parCodes(15) = 41;
	  parCodes(16) = 42;
	  parCodes(17) = 43;
	  parCodes(18) = 44;
	  parCodes(19) = 45;
	  parCodes(20) = 46;
	  parCodes(21) = 47;
	  parCodes(22) = 48;
	  parCodes(23) = 49;
	  parCodes(24) = 61;
	  parCodes(25) = 62;
	  parCodes(26) = 63;
	  parCodes(27) = 64;
	  parCodes(28) = 65;

	  // Vectors for storing the MSSM versions of the parameters
	  DoubleVector mssmParCodes(24);
	  DoubleVector mssmParVals(24);

	  mssmParCodes(1) = 1;
	  mssmParCodes(2) = 2;
	  mssmParCodes(3) = 3;
	  mssmParCodes(4) = 11;
	  mssmParCodes(5) = 12;
	  mssmParCodes(6) = 13;
	  mssmParCodes(7) = 23;
	  mssmParCodes(8) = 24;
	  mssmParCodes(9) = 25;
	  mssmParCodes(10) = 31;
	  mssmParCodes(11) = 32;
	  mssmParCodes(12) = 33;
	  mssmParCodes(13) = 34;
	  mssmParCodes(14) = 35;
	  mssmParCodes(15) = 36;
	  mssmParCodes(16) = 41;
	  mssmParCodes(17) = 42;
	  mssmParCodes(18) = 43;
	  mssmParCodes(19) = 44;
	  mssmParCodes(20) = 45;
	  mssmParCodes(21) = 46;
	  mssmParCodes(22) = 47;
	  mssmParCodes(23) = 48;
	  mssmParCodes(24) = 49;


	  // 2 cases: either all parameters are randomly chosen, or we grid scan.
	  // All parameters are to be randomly scanned
	  if (isRandom)
	    {
	      // Generate a fixed number of random models and write each to file.
	      for (int x = 1; x <= randPoints; x++)
		{
		  // Increment model number
		  modelCount++;

		  // Generate random values for allowed parameters
		  crtParVals(1) = M1_lower + uniformRan(idum)*(M1_upper-M1_lower);
		  crtParVals(2) = M2_lower + uniformRan(idum)*(M2_upper-M2_lower);
		  crtParVals(3) = M3_lower + uniformRan(idum)*(M3_upper-M3_lower);
		  crtParVals(4) = M1p_lower + uniformRan(idum)*(M1p_upper-M1p_lower);
		  crtParVals(5) = At_lower + uniformRan(idum)*(At_upper-At_lower);
		  crtParVals(6) = Ab_lower + uniformRan(idum)*(Ab_upper-Ab_lower);
		  crtParVals(7) = Atau_lower + uniformRan(idum)*(Atau_upper-Atau_lower);
		  crtParVals(8) = tb_lower + uniformRan(idum)*(tb_upper-tb_lower);
		  // Update first and second generation appropriately
		  crtParVals(9) = mL1_lower + uniformRan(idum)*(mL1_upper-mL1_lower);
		  if (isFirstEqualsSecond)
		    {
		      crtParVals(10) = crtParVals(9);
		    }
		  else
		    {
		      crtParVals(10) = mL2_lower + uniformRan(idum)*(mL2_upper-mL2_lower);
		    }
		  crtParVals(11) = mL3_lower + uniformRan(idum)*(mL3_upper-mL3_lower);
		  // Update first and second generation appropriately
		  crtParVals(12) = me1_lower + uniformRan(idum)*(me1_upper-me1_lower);
		  if (isFirstEqualsSecond)
		    {
		      crtParVals(13) = crtParVals(12);
		    }
		  else
		    {
		      crtParVals(13) = me2_lower + uniformRan(idum)*(me2_upper-me2_lower);
		    }

		  crtParVals(14) = me3_lower + uniformRan(idum)*(me3_upper-me3_lower);
		  // Update first and second generation appropriately
		  crtParVals(15) = mQ1_lower + uniformRan(idum)*(mQ1_upper-mQ1_lower);
		  if (isFirstEqualsSecond)
		    {
		      crtParVals(16) = crtParVals(15);
		    }
		  else
		    {
		      crtParVals(16) = mQ2_lower + uniformRan(idum)*(mQ2_upper-mQ2_lower);
		    }
		  crtParVals(17) = mQ3_lower + uniformRan(idum)*(mQ3_upper-mQ3_lower);
		  // Update first and second generation appropriately
		  crtParVals(18) = mu1_lower + uniformRan(idum)*(mu1_upper-mu1_lower);
		  if (isFirstEqualsSecond)
		    {
		      crtParVals(19) = crtParVals(18);
		    }
		  else
		    {
		      crtParVals(19) = mu2_lower + uniformRan(idum)*(mu2_upper-mu2_lower);
		    }
		  crtParVals(20) = mu3_lower + uniformRan(idum)*(mu3_upper-mu3_lower);
		  // Update first and second generation appropriately
		  crtParVals(21) = md1_lower + uniformRan(idum)*(md1_upper-md1_lower);
		  if (isFirstEqualsSecond)
		    {
		      crtParVals(22) = crtParVals(21);
		    }
		  else
		    {
		      crtParVals(22) = md2_lower + uniformRan(idum)*(md2_upper-md2_lower);
		    }
		  crtParVals(23) = md3_lower + uniformRan(idum)*(md3_upper-md3_lower);
		  crtParVals(24) = lambda_lower + uniformRan(idum)*(lambda_upper-lambda_lower);
		  crtParVals(25) = kappa_lower + uniformRan(idum)*(kappa_upper-kappa_lower);
		  crtParVals(26) = Alambda_lower + uniformRan(idum)*(Alambda_upper-Alambda_lower);
		  crtParVals(27) = Akappa_lower + uniformRan(idum)*(Akappa_upper-Akappa_lower);
		  crtParVals(28) = s*crtParVals(24)/sqrt(2.0);

		  // Update MSSM parameter values
		  mssmParVals(1) = crtParVals(1); // M_1
		  mssmParVals(2) = crtParVals(2); // M_2
		  mssmParVals(3) = crtParVals(3); // M_3
		  mssmParVals(4) = crtParVals(5); // A_t
		  mssmParVals(5) = crtParVals(6); // A_b
		  mssmParVals(6) = crtParVals(7); // A_tau
		  mssmParVals(7) = crtParVals(28); // mu
		  mssmParVals(8) = 2.0*crtParVals(28)*crtParVals(26)/(2.0*crtParVals(8)/(1.0+sqr(crtParVals(8)))); // m_A^2
		  mssmParVals(9) = crtParVals(8);
		  mssmParVals(10) = crtParVals(9); // m_eL
		  mssmParVals(11) = crtParVals(10); // m_muL
		  mssmParVals(12) = crtParVals(11); // m_tauL
		  mssmParVals(13) = crtParVals(12); // m_eR
		  mssmParVals(14) = crtParVals(13); // m_muR
		  mssmParVals(15) = crtParVals(14); // m_tauR
		  mssmParVals(16) = crtParVals(15); // m_Q1
		  mssmParVals(17) = crtParVals(16); // m_Q2
		  mssmParVals(18) = crtParVals(17); // m_Q3
		  mssmParVals(19) = crtParVals(18); // m_uR
		  mssmParVals(20) = crtParVals(19); // m_cR
		  mssmParVals(21) = crtParVals(20); // m_tR
		  mssmParVals(22) = crtParVals(21); // m_dR
		  mssmParVals(23) = crtParVals(22); // m_sR
		  mssmParVals(24) = crtParVals(23); // m_bR

		  // Update value of MINPAR tan(beta)
		  minpar_tb = crtParVals(8);

		  // Write to file
		  ostringstream convert;

		  convert << modelCount;

		  if (isMSSM == false)
		    {
		      fileName = filePrefix + "essm_" + convert.str();
		      ofstream paramsFout(fileName.c_str());
		  
		      // Write SLHA style input file.
		      writeSLHAInputFile(paramsFout, modelType, oneset, G_F, hasMX, MX, minpar_tb, parCodes, crtParVals,
					 tol, mix, include2lpScalars);
		      
		      paramsFout.close();

		    }
		  else
		    {
		      fileName = filePrefix + "essm_" + convert.str();
		      ofstream paramsFout(fileName.c_str());
		  
		      // Write SLHA style input file for E6 model.
		      writeSLHAInputFile(paramsFout, modelType, oneset, G_F, hasMX, MX, minpar_tb, parCodes, crtParVals,
					 tol, mix, include2lpScalars);
		      
		      paramsFout.close();

		      fileName = filePrefix + "mssm_" + convert.str();

		      ofstream paramsFout2(fileName.c_str());

		      // Write SLHA style input file for MSSM model.
		      writeSLHAInputFile(paramsFout2, altModelType, oneset, G_F, hasMX, MX, minpar_tb, mssmParCodes, mssmParVals,
					 tol, mix, include2lpScalars);
		      
		      paramsFout2.close();
		    }

		}

	    }
	  // All parameters are to be grid scanned
	  else 
	    {
	      // Work out increments for each parameter
	      double tb_incr = 0.0;
	      double lambda_incr = 0.0;
	      double kappa_incr = 0.0;
	      double Alambda_incr = 0.0;
	      double Akappa_incr = 0.0;
	      double M1_incr = 0.0;
	      double M2_incr = 0.0;
	      double M3_incr = 0.0;
	      double M1p_incr = 0.0;
	      double mQ3_incr = 0.0;
	      double mu3_incr = 0.0;
	      double md3_incr = 0.0;
	      double mL3_incr = 0.0;
	      double me3_incr = 0.0;
	      double At_incr = 0.0;
	      double Ab_incr = 0.0;
	      double Atau_incr = 0.0;

	      double mQ1_incr = 0.0;
	      double mu1_incr = 0.0;
	      double md1_incr = 0.0;
	      double mL1_incr = 0.0;
	      double me1_incr = 0.0;
	      double mQ2_incr = 0.0;
	      double mu2_incr = 0.0;
	      double md2_incr = 0.0;
	      double mL2_incr = 0.0;
	      double me2_incr = 0.0;
	      
	      if (tb_npts != 1)
		{
		  tb_incr = (tb_upper-tb_lower)/(((double)tb_npts)-1.0);
		}
	      if (M1_npts != 1)
		{
		  M1_incr = (M1_upper-M1_lower)/(((double)M1_npts)-1.0);
		}
	      if (M2_npts != 1)
		{
		  M2_incr = (M2_upper-M2_lower)/(((double)M2_npts)-1.0);
		}
	      if (M3_npts != 1)
		{
		  M3_incr = (M3_upper-M3_lower)/(((double)M3_npts)-1.0);
		}
	      if (M1p_npts != 1)
		{
		  M1p_incr = (M1p_upper-M1p_lower)/(((double)M1p_npts)-1.0);
		}
	      // First generation soft masses
	      if (mQ1_npts != 1)
		{
		  mQ1_incr = (mQ1_upper-mQ1_lower)/(((double)mQ1_npts)-1.0);
		}
	      if (mu1_npts != 1)
		{
		  mu1_incr = (mu1_upper-mu1_lower)/(((double)mu1_npts)-1.0);
		}
	      if (md1_npts != 1)
		{
		  md1_incr = (md1_upper-md1_lower)/(((double)md1_npts)-1.0);
		}
	      if (mL1_npts != 1)
		{
		  mL1_incr = (mL1_upper-mL1_lower)/(((double)mL1_npts)-1.0);
		}
	      if (me1_npts != 1)
		{
		  me1_incr = (me1_upper-me1_lower)/(((double)me1_npts)-1.0);
		}
	      // Second generation soft masses
	      if (mQ2_npts != 1)
		{
		  mQ2_incr = (mQ2_upper-mQ2_lower)/(((double)mQ2_npts)-1.0);
		}
	      if (mu2_npts != 1)
		{
		  mu2_incr = (mu2_upper-mu2_lower)/(((double)mu2_npts)-1.0);
		}
	      if (md2_npts != 1)
		{
		  md2_incr = (md2_upper-md2_lower)/(((double)md2_npts)-1.0);
		}
	      if (mL2_npts != 1)
		{
		  mL2_incr = (mL2_upper-mL2_lower)/(((double)mL2_npts)-1.0);
		}
	      if (me2_npts != 1)
		{
		  me2_incr = (me2_upper-me2_lower)/(((double)me2_npts)-1.0);
		}
	      // Third generation soft masses
	      if (mQ3_npts != 1)
		{
		  mQ3_incr = (mQ3_upper-mQ3_lower)/(((double)mQ3_npts)-1.0);
		}
	      if (mu3_npts != 1)
		{
		  mu3_incr = (mu3_upper-mu3_lower)/(((double)mu3_npts)-1.0);
		}
	      if (md3_npts != 1)
		{
		  md3_incr = (md3_upper-md3_lower)/(((double)md3_npts)-1.0);
		}
	      if (mL3_npts != 1)
		{
		  mL3_incr = (mL3_upper-mL3_lower)/(((double)mL3_npts)-1.0);
		}
	      if (me3_npts != 1)
		{
		  me3_incr = (me3_upper-me3_lower)/(((double)me3_npts)-1.0);
		}
	      // Soft A trilinears
	      if (At_npts != 1)
		{
		  At_incr = (At_upper-At_lower)/(((double)At_npts)-1.0);
		}
	      if (Ab_npts != 1)
		{
		  Ab_incr = (Ab_upper-Ab_lower)/(((double)Ab_npts)-1.0);
		}
	      if (Atau_npts != 1)
		{
		  Atau_incr = (Atau_upper-Atau_lower)/(((double)Atau_npts)-1.0);
		}
	      // E6SSM specific parameters
	      if (lambda_npts != 1)
		{
		  lambda_incr = (lambda_upper-lambda_lower)/(((double)lambda_npts)-1.0);
		}
	      if (kappa_npts != 1)
		{
		  kappa_incr = (kappa_upper-kappa_lower)/(((double)kappa_npts)-1.0);
		}
	      if (Alambda_npts != 1)
		{
		  Alambda_incr = (Alambda_upper-Alambda_lower)/(((double)Alambda_npts)-1.0);
		}
	      if (Akappa_npts != 1)
		{
		  Akappa_incr = (Akappa_upper-Akappa_lower)/(((double)Akappa_npts)-1.0);
		}

	      if (isFirstEqualsSecond)
		{
	      // Loop over all parameter values and write to file
	      for (int i1 = 0; i1 < M1_npts; i1++)
		{
		  crtParVals(1) = M1_lower+((double)i1)*M1_incr;
		  mssmParVals(1) = crtParVals(1);

	      for (int i2 = 0; i2 < M2_npts; i2++)
		{
		  crtParVals(2) = M2_lower+((double)i2)*M2_incr;
		  mssmParVals(2) = crtParVals(2);

	      for (int i3 = 0; i3 < M3_npts; i3++)
		{
		  crtParVals(3) = M3_lower+((double)i3)*M3_incr;
		  mssmParVals(3) = crtParVals(3);

	      for (int i4 = 0; i4 < M1p_npts; i4++)
		{
		  crtParVals(4) = M1p_lower+((double)i4)*M1p_incr;

	      for (int i5 = 0; i5 < At_npts; i5++)
		{
		  crtParVals(5) = At_lower+((double)i5)*At_incr;
		  mssmParVals(4) = crtParVals(5);

	      for (int i6 = 0; i6 < Ab_npts; i6++)
		{
		  crtParVals(6) = Ab_lower+((double)i6)*Ab_incr;
		  mssmParVals(5) = crtParVals(6);

	      for (int i7 = 0; i7 < Atau_npts; i7++)
		{
		  crtParVals(7) = Atau_lower+((double)i7)*Atau_incr;
		  mssmParVals(6) = crtParVals(7);

	      for (int i8 = 0; i8 < tb_npts; i8++)
		{
		  crtParVals(8) = tb_lower+((double)i8)*tb_incr;
		  mssmParVals(9) = crtParVals(8);
		  minpar_tb = crtParVals(8);
	      for (int i9 = 0; i9 < mL1_npts; i9++)
		{
		  crtParVals(9) = mL1_lower+((double)i9)*mL1_incr;
		  mssmParVals(10) = crtParVals(9);
		  crtParVals(10) = crtParVals(9);
		  mssmParVals(11) = crtParVals(9);
	    
	      for (int i11 = 0; i11 < mL3_npts; i11++)
		{
		  crtParVals(11) = mL3_lower+((double)i11)*mL3_incr;
		  mssmParVals(12) = crtParVals(11);

	      for (int i12 = 0; i12 < me1_npts; i12++)
		{
		  crtParVals(12) = me1_lower+((double)i12)*me1_incr;
		  crtParVals(13) = crtParVals(12);
		  mssmParVals(13) = crtParVals(12);
		  mssmParVals(14) = crtParVals(12);

	      for (int i14 = 0; i14 < me3_npts; i14++)
		{
		  crtParVals(14) = me3_lower+((double)i14)*me3_incr;
		  mssmParVals(15) = crtParVals(14);

	      for (int i15 = 0; i15 < mQ1_npts; i15++)
		{
		  crtParVals(15) = mQ1_lower+((double)i15)*mQ1_incr;
		  crtParVals(16) = crtParVals(15);
		  mssmParVals(16) = crtParVals(15);
		  mssmParVals(17) = crtParVals(15);

	      for (int i17 = 0; i17 < mQ3_npts; i17++)
		{
		  crtParVals(17) = mQ3_lower+((double)i17)*mQ3_incr;
		  mssmParVals(18) = crtParVals(17);

	      for (int i18 = 0; i18 < mu1_npts; i18++)
		{
		  crtParVals(18) = mu1_lower+((double)i18)*mu1_incr;
		  crtParVals(19) = crtParVals(18);
		  mssmParVals(19) = crtParVals(18);
		  mssmParVals(20) = crtParVals(18);

	      for (int i20 = 0; i20 < mu3_npts; i20++)
		{
		  crtParVals(20) = mu3_lower+((double)i20)*mu3_incr;
		  mssmParVals(21) = crtParVals(20);

	      for (int i21 = 0; i21 < md1_npts; i21++)
		{
		  crtParVals(21) = md1_lower+((double)i21)*md1_incr;
		  crtParVals(22) = crtParVals(21);
		  mssmParVals(22) = crtParVals(21);
		  mssmParVals(23) = crtParVals(21);

	      for (int i23 = 0; i23 < md3_npts; i23++)
		{
		  crtParVals(23) = md3_lower+((double)i23)*md3_incr;
		  mssmParVals(24) = crtParVals(23);

	      for (int i24 = 0; i24 < lambda_npts; i24++)
		{
		  crtParVals(24) = lambda_lower+((double)i24)*lambda_incr;
		  crtParVals(28) = s*crtParVals(24)/sqrt(2.0);
		  mssmParVals(7) = crtParVals(28);

	      for (int i25 = 0; i25 < kappa_npts; i25++)
		{
		  crtParVals(25) = kappa_lower+((double)i25)*kappa_incr;

	      for (int i26 = 0; i26 < Alambda_npts; i26++)
		{
		  crtParVals(26) = Alambda_lower+((double)i26)*Alambda_incr;
		  mssmParVals(8) = 2.0*crtParVals(28)*crtParVals(26)/(2.0*crtParVals(8)/(1.0+sqr(crtParVals(8))));

	      for (int i27 = 0; i27 < Akappa_npts; i27++)
		{
		  crtParVals(27) = Akappa_lower+((double)i27)*Akappa_incr;

		  // Increment model count
		  modelCount++;

		  // Write files
		  // Write to file
		  ostringstream convert;

		  convert << modelCount;

		  if (isMSSM == false)
		    {
		      fileName = filePrefix + "essm_" + convert.str();
		      ofstream paramsFout(fileName.c_str());
		  
		      // Write SLHA style input file.
		      writeSLHAInputFile(paramsFout, modelType, oneset, G_F, hasMX, MX, minpar_tb, parCodes, crtParVals, 
					 tol, mix, include2lpScalars);
		      
		      paramsFout.close();

		    }
		  else
		    {
		      fileName = filePrefix + "essm_" + convert.str();
		      ofstream paramsFout(fileName.c_str());
		  
		      // Write SLHA style input file for E6 model.
		      writeSLHAInputFile(paramsFout, modelType, oneset, G_F, hasMX, MX, minpar_tb, parCodes, crtParVals,
					 tol, mix, include2lpScalars);
		      
		      paramsFout.close();

		      fileName = filePrefix + "mssm_" + convert.str();

		      ofstream paramsFout2(fileName.c_str());

		      // Write SLHA style input file for MSSM model.
		      writeSLHAInputFile(paramsFout2, altModelType, oneset, G_F, hasMX, MX, minpar_tb, mssmParCodes, mssmParVals,
					 tol, mix, include2lpScalars);
		      
		      paramsFout2.close();
		    }


		} // Akappa loop
		} // Alambda loop
	        } // kappa loop
		} // lambda loop
		} // md3 loop
	        } // md1 loop
		} // mu3 loop
	        } // mu1 loop
		} // mQ3 loop
	        } // mQ1 loop
		} // me3 loop
	        } // me1 loop
		} // mL3 loop
	        } // mL1 loop
		} // tb loop
		} // Atau loop
		} // Ab loop
		} // At loop
		} // M1p loop
		} // M3 loop
		} // M2 loop
		} // M1 loop
		} // End of degenerate first and second generation scan
	      else
		{
	      // Loop over all parameter values and write to file
	      for (int i1 = 0; i1 < M1_npts; i1++)
		{
		  crtParVals(1) = M1_lower+((double)i1)*M1_incr;
		  mssmParVals(1) = crtParVals(1);

	      for (int i2 = 0; i2 < M2_npts; i2++)
		{
		  crtParVals(2) = M2_lower+((double)i2)*M2_incr;
		  mssmParVals(2) = crtParVals(2);

	      for (int i3 = 0; i3 < M3_npts; i3++)
		{
		  crtParVals(3) = M3_lower+((double)i3)*M3_incr;
		  mssmParVals(3) = crtParVals(3);

	      for (int i4 = 0; i4 < M1p_npts; i4++)
		{
		  crtParVals(4) = M1p_lower+((double)i4)*M1p_incr;

	      for (int i5 = 0; i5 < At_npts; i5++)
		{
		  crtParVals(5) = At_lower+((double)i5)*At_incr;
		  mssmParVals(4) = crtParVals(5);

	      for (int i6 = 0; i6 < Ab_npts; i6++)
		{
		  crtParVals(6) = Ab_lower+((double)i6)*Ab_incr;
		  mssmParVals(5) = crtParVals(6);

	      for (int i7 = 0; i7 < Atau_npts; i7++)
		{
		  crtParVals(7) = Atau_lower+((double)i7)*Atau_incr;
		  mssmParVals(6) = crtParVals(7);

	      for (int i8 = 0; i8 < tb_npts; i8++)
		{
		  crtParVals(8) = tb_lower+((double)i8)*tb_incr;
		  mssmParVals(9) = crtParVals(8);
		  minpar_tb = crtParVals(8);

	      for (int i9 = 0; i9 < mL1_npts; i9++)
		{
		  crtParVals(9) = mL1_lower+((double)i9)*mL1_incr;
		  mssmParVals(10) = crtParVals(9);

	      for (int i10 = 0; i10 < mL2_npts; i10++)
		{
		  crtParVals(10) = mL2_lower+((double)i10)*mL2_incr;
		  mssmParVals(11) = crtParVals(10);
	    
	      for (int i11 = 0; i11 < mL3_npts; i11++)
		{
		  crtParVals(11) = mL3_lower+((double)i11)*mL3_incr;
		  mssmParVals(12) = crtParVals(11);

	      for (int i12 = 0; i12 < me1_npts; i12++)
		{
		  crtParVals(12) = me1_lower+((double)i12)*me1_incr;
		  mssmParVals(13) = crtParVals(12);

	      for (int i13 = 0; i13 < me2_npts; i13++)
		{
		  crtParVals(13) = me2_lower+((double)i13)*me2_incr;
		  mssmParVals(14) = crtParVals(13);

	      for (int i14 = 0; i14 < me3_npts; i14++)
		{
		  crtParVals(14) = me3_lower+((double)i14)*me3_incr;
		  mssmParVals(15) = crtParVals(14);

	      for (int i15 = 0; i15 < mQ1_npts; i15++)
		{
		  crtParVals(15) = mQ1_lower+((double)i15)*mQ1_incr;
		  mssmParVals(16) = crtParVals(15);

	      for (int i16 = 0; i16 < mQ2_npts; i16++)
		{
		  crtParVals(16) = mQ2_lower+((double)i16)*mQ2_incr;
		  mssmParVals(17) = crtParVals(16);

	      for (int i17 = 0; i17 < mQ3_npts; i17++)
		{
		  crtParVals(17) = mQ3_lower+((double)i17)*mQ3_incr;
		  mssmParVals(18) = crtParVals(17);

	      for (int i18 = 0; i18 < mu1_npts; i18++)
		{
		  crtParVals(18) = mu1_lower+((double)i18)*mu1_incr;
		  mssmParVals(19) = crtParVals(18);

	      for (int i19 = 0; i19 < mu2_npts; i19++)
		{
		  crtParVals(19) = mu2_lower+((double)i19)*mu2_incr;
		  mssmParVals(20) = crtParVals(19);

	      for (int i20 = 0; i20 < mu3_npts; i20++)
		{
		  crtParVals(20) = mu3_lower+((double)i20)*mu3_incr;
		  mssmParVals(21) = crtParVals(20);

	      for (int i21 = 0; i21 < md1_npts; i21++)
		{
		  crtParVals(21) = md1_lower+((double)i21)*md1_incr;
		  mssmParVals(22) = crtParVals(21);

	      for (int i22 = 0; i22 < md2_npts; i22++)
		{
		  crtParVals(22) = md2_lower+((double)i22)*md2_incr;
		  mssmParVals(23) = crtParVals(22);

	      for (int i23 = 0; i23 < md3_npts; i23++)
		{
		  crtParVals(23) = md3_lower+((double)i23)*md3_incr;
		  mssmParVals(24) = crtParVals(23);

	      for (int i24 = 0; i24 < lambda_npts; i24++)
		{
		  crtParVals(24) = lambda_lower+((double)i24)*lambda_incr;
		  crtParVals(28) = s*crtParVals(24)/sqrt(2.0);
		  mssmParVals(7) = crtParVals(28);

	      for (int i25 = 0; i25 < kappa_npts; i25++)
		{
		  crtParVals(25) = kappa_lower+((double)i25)*kappa_incr;

	      for (int i26 = 0; i26 < Alambda_npts; i26++)
		{
		  crtParVals(26) = Alambda_lower+((double)i26)*Alambda_incr;
		  mssmParVals(8) = 2.0*crtParVals(28)*crtParVals(26)/(2.0*crtParVals(8)/(1.0+sqr(crtParVals(8))));

	      for (int i27 = 0; i27 < Akappa_npts; i27++)
		{
		  crtParVals(27) = Akappa_lower+((double)i27)*Akappa_incr;

		  // Increment model count
		  modelCount++;

		  // Write files
		  // Write to file
		  ostringstream convert;

		  convert << modelCount;

		  if (isMSSM == false)
		    {
		      fileName = filePrefix + "essm_" + convert.str();
		      ofstream paramsFout(fileName.c_str());
		  
		      // Write SLHA style input file.
		      writeSLHAInputFile(paramsFout, modelType, oneset, G_F, hasMX, MX, minpar_tb, parCodes, crtParVals,
					 tol, mix, include2lpScalars);
		      
		      paramsFout.close();

		    }
		  else
		    {
		      fileName = filePrefix + "essm_" + convert.str();
		      ofstream paramsFout(fileName.c_str());
		  
		      // Write SLHA style input file for E6 model.
		      writeSLHAInputFile(paramsFout, modelType, oneset, G_F, hasMX, MX, minpar_tb, parCodes, crtParVals, 
					 tol, mix, include2lpScalars);
		      
		      paramsFout.close();

		      fileName = filePrefix + "mssm_" + convert.str();

		      ofstream paramsFout2(fileName.c_str());

		      // Write SLHA style input file for MSSM model.
		      writeSLHAInputFile(paramsFout2, altModelType, oneset, G_F, hasMX, MX, minpar_tb, mssmParCodes, mssmParVals, 
					 tol, mix, include2lpScalars);
		      
		      paramsFout2.close();
		    }


		} // Akappa loop
		} // Alambda loop
	        } // kappa loop
		} // lambda loop
		} // md3 loop
	        } // md2 loop
	        } // md1 loop
		} // mu3 loop
	        } // mu2 loop
	        } // mu1 loop
		} // mQ3 loop
	        } // mQ2 loop
	        } // mQ1 loop
		} // me3 loop
	        } // me2 loop
	        } // me1 loop
		} // mL3 loop
	        } // mL2 loop
	        } // mL1 loop
		} // tb loop
		} // Atau loop
		} // Ab loop
		} // At loop
		} // M1p loop
		} // M3 loop
		} // M2 loop
		} // M1 loop
		} // End of non-degenerate first and second generation scan
	    }
	}
    }
  // // Exception handling, in the unlikely event that something
  // // does go wrong.
  catch(const string & a) { cout << a; }
  catch(const char * a) { cout << a; }
  catch(...) { cout << "Unknown type of exception caught.\n"; }

  return 0;
}


void errorCall()
{
  ostringstream outstream;
  outstream << "Error: incorrect syntax.\n";
  outstream << "Usage: ./generateE6SSMSLHAInputs < specsFile\n";
  outstream << "The file specsFile should contain the required input\n";
  outstream << "file parameters.\n";
  throw outstream.str();
}
// ToUpper method from softpoint.cpp
string ToUpper(const string & s) 
{
  string result;
  unsigned int index;
  for (index = 0; index < s.length(); index++) {
    char a = s[index];
    a = toupper(a);
    result = result + a;
  }
  
  return result;
}

// Method to write an SLHA style input file for E6SSM models.
// Requires as inputs a file stream to write to, all necessary SM inputs, a boolean flag indicating
// whether MX is taken to be M_GUT, the value of MX to use if requested (-1 if MX = M_SUSY, > 0 otherwise)
// the value of MINPAR tan(beta), a vector of SLHA codes for the requested EXTPARS, and a vector
// of values for those EXTPARS.
void writeSLHAInputFile(ofstream & out, const string & modelType, QedQcd const & oneset, double G_F, bool hasMX, double MX, 
			double tbVal, DoubleVector const & parNums, DoubleVector const & parVals, double tolerance, 
			int mixing, bool include2lpScalars)
{

  int outprecis = out.precision();

  slhaInputFileHeader(out);
  slhaInputModselBlock(out, modelType);
  slhaInputSminputBlock(out, oneset, G_F);
  slhaInputMinparBlock(out, tbVal);  
  slhaInputExtparBlock(out, hasMX, MX, parNums, parVals);
  slhaInputSoftSusyBlock(out, tolerance, mixing, include2lpScalars);

  out.precision(outprecis);
}


void slhaInputFileHeader(ofstream & out)
{
  out.setf(ios::scientific, ios::floatfield);
  out.precision(8);

  out << "# Input in SLHA2-like format, and suitable for input to our" << endl;
  out << "# fine tuning code. " << endl;
  out << "# Note that for convenience E6SSM input is in an SLHA like format," << endl;
  out << "# with the parameter labels used corresponding to those for the analogous " << endl;
  out << "# parameters in the NMSSM (see hep-ph/0801.0045v3)." << endl;
}

void slhaInputModselBlock(ofstream & out, const string & model)
{
  out << "Block MODSEL                 # Model selection\n"; 
  int modsel = 2;
  if (!strcmp(model.c_str(), "E6SSM")) 
    {
      modsel = 2; // there is only one option for now
      out << "    3    " << modsel << "                   # " << model << "\n"; /// Les Houches
      /// accord codes
    }
  else if (!strcmp(model.c_str(), "MSSM"))
    {
      modsel = 0;
      out << "    1    " << modsel << "                   # " << model << "\n"; /// Les Houches
    }
}

void slhaInputSoftSusyBlock(ofstream & out, double tol, int mix, bool full2lp)
{
  out << "Block SOFTSUSY               # SOFTSUSY specific inputs\n"; 
  out << "    1  "; printRow(out, tol);
  out << "       # desired fractional accuracy in output\n";
  out << "    2  "; printRow(out, mix);
  out << "       # quark mixing option\n";
  if (full2lp)
    {
      out << "    5  "; printRow(out, 1.0);
      out << "       # full 2-loop running in RGEs\n";
    }
  else
    {
      out << "    5  "; printRow(out, 0.0);
      out << "       # full 2-loop running in RGEs\n";
    }
}
  
void slhaInputSminputBlock(ofstream & out, QedQcd const & oneset, double G_F)
{
  out << "Block SMINPUTS               # Standard Model inputs\n";
  out << "    1  "; printRow(out, 1.0/oneset.displayAlpha(ALPHA));
  out << "       # alpha_em^-1(M_Z)^MSbar\n";
  out << "    2  "; printRow(out, G_F);
  out << "       # G_F [GeV^-2]\n";
  out << "    3  "; printRow(out, oneset.displayAlpha(ALPHAS));
  out << "       # alpha_S(M_Z)^MSbar\n";
  out << "    4  "; printRow(out, oneset.displayMu());
  out << "       # M_Z pole mass\n";
  out << "    5  "; printRow(out, oneset.displayMass(mBottom));
  out << "       # mb(mb)^MSbar\n";
  out << "    6  "; printRow(out, oneset.displayPoleMt());
  out << "       # mt pole mass\n";
  out << "    7  "; printRow(out, oneset.displayPoleMtau());
  out << "       # mtau pole mass\n";
}

void slhaInputMinparBlock(ofstream & out, double tb)
{
  out << "Block MINPAR                 # Input parameters - minimal models\n";
  out << "    3  "; printRow(out, tb);
  out << "       # tanb\n";
}

void slhaInputExtparBlock(ofstream & out, bool hasMX, double MX, DoubleVector const & parNums, DoubleVector const & parVals)
{
  out << "Block EXTPAR                 # Input parameters - non-minimal models\n";
  int crtParNum;
  int numPars = parNums.displayEnd()-parNums.displayStart()+1;
  if (hasMX)
    {
      if (MX < 0.0)
	{
	  out << "    0  "; printRow(out, -1); out << "       # MX scale\n";
	}
      else
	{
	  out << "    0  "; printRow(out, MX); out << "       # MX scale\n"; 
	}
    }
  // Assume parameter list is already sorted
  for (int i = 1; i<= numPars; i++)
    {
      crtParNum = (int) parNums(i);
      switch(crtParNum)
	{
	case 1: out << "    1  "; printRow(out, parVals(i)); out << "       # M_1(MX)\n"; break;
	case 2: out << "    2  "; printRow(out, parVals(i)); out << "       # M_2(MX)\n"; break;
	case 3: out << "    3  "; printRow(out, parVals(i)); out << "       # M_3(MX)\n"; break;
	case 4: out << "    4  "; printRow(out, parVals(i)); out << "       # M_1'(MX)\n"; break;
	case 11: out << "   11  "; printRow(out, parVals(i)); out << "       # At(MX)\n"; break;
	case 12: out << "   12  "; printRow(out, parVals(i)); out << "       # Ab(MX)\n"; break;
	case 13: out << "   13  "; printRow(out, parVals(i)); out << "       # Atau(MX)\n"; break;
	case 21: out << "   21  "; printRow(out, parVals(i)); out << "       # mHd^2(MX)\n"; break;
	case 22: out << "   22  "; printRow(out, parVals(i)); out << "       # mHu^2(MX)\n"; break;
	case 23: out << "   23  "; printRow(out, parVals(i)); out << "       # mu(MX)\n"; break;
	case 24: out << "   24  "; printRow(out, parVals(i)); out << "       # mA^2(MX)\n"; break;
	case 25: out << "   25  "; printRow(out, parVals(i)); out << "       # tanb(MX)\n"; break;
	case 26: out << "   26  "; printRow(out, parVals(i)); out << "       # mA(pole)\n"; break;
	case 27: out << "   27  "; printRow(out, parVals(i)); out << "       # mHpm(pole)\n"; break;
	case 31: out << "   31  "; printRow(out, parVals(i)); out << "       # meL(MX)\n"; break;
	case 32: out << "   32  "; printRow(out, parVals(i)); out << "       # mmuL(MX)\n"; break;
	case 33: out << "   33  "; printRow(out, parVals(i)); out << "       # mtauL(MX)\n"; break;
	case 34: out << "   34  "; printRow(out, parVals(i)); out << "       # meR(MX)\n"; break;
	case 35: out << "   35  "; printRow(out, parVals(i)); out << "       # mmuR(MX)\n"; break;
	case 36: out << "   36  "; printRow(out, parVals(i)); out << "       # mtauR(MX)\n"; break;
	case 41: out << "   41  "; printRow(out, parVals(i)); out << "       # mqL1(MX)\n"; break;
	case 42: out << "   42  "; printRow(out, parVals(i)); out << "       # mqL2(MX)\n"; break;
	case 43: out << "   43  "; printRow(out, parVals(i)); out << "       # mqL3(MX)\n"; break;
	case 44: out << "   44  "; printRow(out, parVals(i)); out << "       # muR(MX)\n"; break;
	case 45: out << "   45  "; printRow(out, parVals(i)); out << "       # mcR(MX)\n"; break;
	case 46: out << "   46  "; printRow(out, parVals(i)); out << "       # mtR(MX)\n"; break;
	case 47: out << "   47  "; printRow(out, parVals(i)); out << "       # mdR(MX)\n"; break;
	case 48: out << "   48  "; printRow(out, parVals(i)); out << "       # msR(MX)\n"; break;
	case 49: out << "   49  "; printRow(out, parVals(i)); out << "       # mbR(MX)\n"; break;
	case 61: out << "   61  "; printRow(out, parVals(i)); out << "       # lambda_3(MX)\n"; break;
	case 62: out << "   62  "; printRow(out, parVals(i)); out << "       # kappa_3(MX)\n"; break;
	case 63: out << "   63  "; printRow(out, parVals(i)); out << "       # Alambda_3(MX)\n"; break;
	case 64: out << "   64  "; printRow(out, parVals(i)); out << "       # Akappa_3(MX)\n"; break;
	case 65: out << "   65  "; printRow(out, parVals(i)); out << "       # lambda_3<S_3>(MX)\n"; break;
	case 70: out << "   70  "; printRow(out, parVals(i)); out << "       # mS_3^2(MX)\n"; break;
	case 71: out << "   71  "; printRow(out, parVals(i)); out << "       # mS_2^2(MX)\n"; break;
	case 72: out << "   72  "; printRow(out, parVals(i)); out << "       # mS_1^2(MX)\n"; break;
	case 73: out << "   73  "; printRow(out, parVals(i)); out << "       # mHd2^2(MX)\n"; break;
	case 74: out << "   74  "; printRow(out, parVals(i)); out << "       # mHu2^2(MX)\n"; break;
	case 75: out << "   75  "; printRow(out, parVals(i)); out << "       # mHd1^2(MX)\n"; break;
	case 76: out << "   76  "; printRow(out, parVals(i)); out << "       # mHu1^2(MX)\n"; break;
	case 81: out << "   81  "; printRow(out, parVals(i)); out << "       # lambda_2(MX)\n"; break;
	case 82: out << "   82  "; printRow(out, parVals(i)); out << "       # kappa_2(MX)\n"; break;
	case 83: out << "   83  "; printRow(out, parVals(i)); out << "       # Alambda_2(MX)\n"; break;
	case 84: out << "   84  "; printRow(out, parVals(i)); out << "       # Akappa_2(MX)\n"; break;
	case 85: out << "   85  "; printRow(out, parVals(i)); out << "       # lambda_2<S_2>(MX)\n"; break;
	case 91: out << "   91  "; printRow(out, parVals(i)); out << "       # lambda_1(MX)\n"; break;
	case 92: out << "   92  "; printRow(out, parVals(i)); out << "       # kappa_1(MX)\n"; break;
	case 93: out << "   93  "; printRow(out, parVals(i)); out << "       # Alambda_1(MX)\n"; break;
	case 94: out << "   94  "; printRow(out, parVals(i)); out << "       # Akappa_1(MX)\n"; break;
	case 95: out << "   95  "; printRow(out, parVals(i)); out << "       # lambda_1<S_1>(MX)\n"; break;
	default: break;
	 	  
	}
    }
}

// We generate uniformly distributed random numbers using SOFTSUSY routines.
double uniformRan(long & idum)
{
  double n;

  // Get normally distributed random number
  n = gasdev(idum);

  // Convert to uniform distribution
  n = 0.5*(1.0+erf(n/sqrt(2.0)));

  return n;
}

