/*
  scanInputsGenerator produces a set of files,
  each of which contains the parameters for a subscan
  of the pE6SSM parameter space. Use in conjunction
  with emssmScanner to actually do the scans.
  Currently files are written as if scanning over
  tan(beta), lambda, A_lambda, m_Q3^2, m_u3^2, A_t and M_2. All 
  other variables are fixed. There is an optional
  flag to scan over A_t and A_lambda using logarithmically spaced points.
 */

#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <sstream>
#include <cstring>
#include "E6SSM_Spectrum_Generators/src/linalg.h"
#include "E6SSM_Spectrum_Generators/src/utils.h"
#include "E6SSM_Spectrum_Generators/src/numerics.h"

softsusy::DoubleVector generateVals(double lower, double upper, int npts, bool isLogSpaced = false);


void writeScanBoundsFile(ofstream & boundsfout, double tb_lb, double tb_ub, int tb_subnpts, 
			 double lambda_lb, double lambda_ub, int lambda_subnpts, 
			 double Alambda_lb, double Alambda_ub, int Alambda_subnpts, 
			 double mqL3sq_lb, double mqL3sq_ub, int mqL3sq_subnpts,
			 double mtRsq_lb, double mtRsq_ub, int mtRsq_subnpts, 
			 double At_lb, double At_ub, int At_subnpts,
			 double M2_lb, double M2_ub, int M2_subnpts, bool atlogscan, bool alambdalogscan,
			 double yu, double yc, double yt, double yd, double ys, double yb,
			 double ye, double ymu, double ytau, double g1, double g2, double g3, double g1p,
			 double v, double s);


int main(int argc, char* argv[])
{
using namespace softsusy;
  /*
    ----------------------------------------------------
    Parameter lower bounds over entire scan
    ----------------------------------------------------
  */
  
  double tb_low = 2.0;
  double lambda_low = -3.0;
  double Alambda_low = -10000.0; // GeV
  double mqL3sq_low = sqr(200.0); // GeV
  double mtRsq_low = sqr(200.0); // GeV
  double At_low = -10000.0; // GeV

  // Maybe we should do a coarse scan over M_2...
  double M2_low = 100.0;

  /*
    ----------------------------------------------------
   */

  /*
    ----------------------------------------------------
    Parameter upper bounds over entire scan
    ----------------------------------------------------
   */

  double tb_up = 50.0;
  double lambda_up = 3.0;
  double Alambda_up = 10000.0; // GeV
  double mqL3sq_up = sqr(2000.0); // GeV
  double mtRsq_up = sqr(2000.0); // GeV
  double At_up = 10000.0; // GeV

  // For scanning over M_2...
  double M2_up = 2000.0; // GeV

  /*
    ----------------------------------------------------
   */

  /*
    ----------------------------------------------------
    Number of points for each parameter in entire scan
    ----------------------------------------------------
   */

  int tb_npts = 1;
  int lambda_npts = 40;//6;
  int Alambda_npts = 40;//6;
  int mqL3sq_npts = 30;//6;
  int mtRsq_npts = 30;//6;
  int At_npts = 40;//6;

  int M2_npts = 3;

  /*
    ----------------------------------------------------
   */

  /*
    ----------------------------------------------------
    Number of subintervals for each parameter
    ----------------------------------------------------
   */

  int tb_nintvls = 1;
  int lambda_nintvls = 4;//2;
  int Alambda_nintvls = 4;//2;
  int mqL3sq_nintvls = 3;
  int mtRsq_nintvls = 3;
  int At_nintvls = 4;//3;

  int M2_nintvls = 1;

  /*
    ----------------------------------------------------
   */

  /*
    ----------------------------------------------------
    Values of the Higgs vev and dimensionless couplings 
    at MX, calculated using SOFTSUSY.
    ----------------------------------------------------
   */

  // pE6SSM assumption: 1st and 2nd generation Yukawas are zero
  double yu = 0.0;
  double yc = 0.0;
  double yt = 8.40236500e-01;

  double yd = 0.0;
  double ys = 0.0;
  double yb = 6.39364115e-01;

  double ye = 0.0;
  double ymu = 0.0;
  double ytau = 5.29502172e-01;

  double g1 = sqrt(5.0/3.0)*3.71453554e-01;
  double g2 = 6.43227110e-01;
  double g3 = 9.79054815e-01;
  // U(1)_N gauge coupling
  double g1p = g1;

  double v = 2.36246110e+02;
  // Singlet vev. Note current Z' mass limit
  // of ~2.5 TeV ==> s ~ 6.7 TeV
  double s = 6700.0; // GeV

  /*
    ----------------------------------------------------
   */

  /*
    ----------------------------------------------------
    Flags, switches and output file names
    ----------------------------------------------------
  */
  
  bool atlogscan = true;
  bool alambdalogscan = true;
  double const LOWERLOGAT = 1.0;
  double const LOWERLOGALAMBDA = 1.0;
  string outfileprefix = "/home/dylan/Documents/Postgraduate/E6SSM-Tuning/essmScans/essmScan1/essmScan1tb2/essmScan1tb2Inputs/essmScan1tb2";
  
  
  /*
    ----------------------------------------------------
  */
  
  /*
    ----------------------------------------------------
    Actual work starts here.
    ----------------------------------------------------
  */
  
  // Check at least one point is requested for each
  if (tb_npts < 1)
    {
      cerr << "WARNING: must have at least 1 point for tan(beta): using default 1." << endl;
      tb_npts = 1;
    }

  if (lambda_npts < 1)
    {
      cerr << "WARNING: must have at least 1 point for lambda: using default 1." << endl;
      lambda_npts = 1;
    }

  if (Alambda_npts < 1)
    {
      cerr << "WARNING: must have at least 1 point for Alambda: using default 1." << endl;
      Alambda_npts = 1;
    }

  if (mqL3sq_npts < 1)
    {
      cerr << "WARNING: must have at least 1 point for m_{qL_3}^2: using default 1." << endl;
      mqL3sq_npts = 1;
    }

  if (mtRsq_npts < 1)
    {
      cerr << "WARNING: must have at least 1 point for m_{tR}^2: using default 1." << endl;
      mtRsq_npts = 1;
    }

  if (At_npts < 1)
    {
      cerr << "WARNING: must have at least 1 point for A_t: using default 1." << endl;
      At_npts = 1;
    }

  if (M2_npts < 1)
    {
      cerr << "WARNING: must have at least 1 point for M_2: using default 1." << endl;
      M2_npts = 1;
    }

  // Most of the parameters are trivial, because we scan
  // linearly over them. Work out the relevant increments.
  
  DoubleVector tb_vals = generateVals(tb_low, tb_up, tb_npts);
  DoubleVector lambda_vals = generateVals(lambda_low, lambda_up, lambda_npts);

  DoubleVector M2_vals = generateVals(M2_low, M2_up, M2_npts);

  // Note that the parameters that we scan over are m_Q3^2 and
  // m_u3^2, not m_Q3 and m_u3.
  DoubleVector mqL3sq_vals = generateVals(mqL3sq_low, mqL3sq_up, mqL3sq_npts);
  DoubleVector mtRsq_vals = generateVals(mtRsq_low, mtRsq_up, mtRsq_npts);

  // Slightly trickier is A_t. First we need to know if we are
  // doing a linear scan (trivial), or log scan (less trivial).
  // If we do a log scan, ignore the given lower bound. Instead
  // we will do symmetric positive and negative values going from
  // a defined lower limit up to the bound with the largest absolute value.

  DoubleVector At_vals(At_npts);

  if (!atlogscan)
    {
      At_vals = generateVals(At_low, At_up, At_npts);
    }
  else
    {

      double tmplow = LOWERLOGAT;
      double tmpup;

      if (fabs(At_low) > fabs(At_up))
	{
	  tmpup = fabs(At_low);
	}
      else
	{
	  tmpup = fabs(At_up);
	}

      // Make sure result isn't zero. If it is, 
      // shift by lower limit.
      if (tmpup < EPSTOL)
	{
	  tmpup += LOWERLOGAT;
	}

      // Then everything should be safe. We only use half
      // of the requested points (the rest get used for
      // the negative values - this gets worked out by
      // pmssmScanner).
      At_npts = At_npts/2; //< integer division

      At_vals = generateVals(tmplow, tmpup, At_npts, atlogscan);

    }

  // Likewise for A_lambda
  DoubleVector Alambda_vals(Alambda_npts);

  if (!alambdalogscan)
    {
      Alambda_vals = generateVals(Alambda_low, Alambda_up, Alambda_npts);
    }
  else
    {

      double tmplow = LOWERLOGALAMBDA;
      double tmpup;

      if (fabs(Alambda_low) > fabs(Alambda_up))
	{
	  tmpup = fabs(Alambda_low);
	}
      else
	{
	  tmpup = fabs(Alambda_up);
	}

      // Make sure result isn't zero. If it is, 
      // shift by lower limit.
      if (tmpup < EPSTOL)
	{
	  tmpup += LOWERLOGALAMBDA;
	}

      // Then everything should be safe. We only use half
      // of the requested points (the rest get used for
      // the negative values - this gets worked out by
      // pmssmScanner).
      Alambda_npts = Alambda_npts/2; //< integer division

      Alambda_vals = generateVals(tmplow, tmpup, Alambda_npts, alambdalogscan);

    }


  // Check at least one interval is requested for each parameter
  if (tb_nintvls < 1)
    {
      cerr << "WARNING: must have at least 1 interval for tan(beta): using default 1." << endl;
      tb_nintvls = 1;
    }

  if (lambda_nintvls < 1)
    {
      cerr << "WARNING: must have at least 1 interval for lambda: using default 1." << endl;
      lambda_nintvls = 1;
    }

  if (Alambda_nintvls < 1)
    {
      cerr << "WARNING: must have at least 1 interval for A_lambda: using default 1." << endl;
      Alambda_nintvls = 1;
    }

  if (mqL3sq_nintvls < 1)
    {
      cerr << "WARNING: must have at least 1 interval for m_{qL_3}^2: using default 1." << endl;
      mqL3sq_nintvls = 1;
    }

  if (mtRsq_nintvls < 1)
    {
      cerr << "WARNING: must have at least 1 interval for m_{tR}^2: using default 1." << endl;
      mtRsq_nintvls = 1;
    }

  if (At_nintvls < 1)
    {
      cerr << "WARNING: must have at least 1 interval for A_t: using default 1." << endl;
      At_nintvls = 1;
    }

  if (M2_nintvls < 1)
    {
      cerr << "WARNING: must have at least 1 interval for M_2: using default 1." << endl;
      M2_nintvls = 1;
    }

  // Now break everything up into the requested number of 
  // subintervals to do the file writing, i.e. work out
  // the indices corresponding to the endpoint of each 
  // parameter interval. First and last entries in each
  // vector are the positions of the lower and upper bounds
  // (i.e. if no sub-division then vector of indices has length
  // 2).
  DoubleVector tb_indices(1+tb_nintvls);
  DoubleVector lambda_indices(1+lambda_nintvls);
  DoubleVector Alambda_indices(1+Alambda_nintvls);
  DoubleVector mqL3sq_indices(1+mqL3sq_nintvls);
  DoubleVector mtRsq_indices(1+mtRsq_nintvls);
  DoubleVector At_indices(1+At_nintvls);

  DoubleVector M2_indices(1+M2_nintvls);

  tb_indices(tb_indices.displayEnd()) = tb_npts+1;
  lambda_indices(lambda_indices.displayEnd()) = lambda_npts+1;
  Alambda_indices(Alambda_indices.displayEnd()) = Alambda_npts+1;
  mqL3sq_indices(mqL3sq_indices.displayEnd()) = mqL3sq_npts+1;
  mtRsq_indices(mtRsq_indices.displayEnd()) = mtRsq_npts+1;
  At_indices(At_indices.displayEnd()) = At_npts+1;
  M2_indices(M2_indices.displayEnd()) = M2_npts+1;

  int tb_step = tb_npts/tb_nintvls;
  int lambda_step = lambda_npts/lambda_nintvls;
  int Alambda_step = Alambda_npts/Alambda_nintvls;
  int mqL3sq_step = mqL3sq_npts/mqL3sq_nintvls;
  int mtRsq_step = mtRsq_npts/mtRsq_nintvls;
  int At_step = At_npts/At_nintvls;

  int M2_step = M2_npts/M2_nintvls;

  for (int i = 1; i <= tb_nintvls; i++)
    {
      tb_indices(i) = 1 + tb_step*(i-1);
    }

  for (int i = 1; i <= lambda_nintvls; i++)
    {
      lambda_indices(i) = 1 + lambda_step*(i-1);
    }

  for (int i = 1; i <= Alambda_nintvls; i++)
    {
      Alambda_indices(i) = 1 + Alambda_step*(i-1);
    }

  for (int i = 1; i <= mqL3sq_nintvls; i++)
    {
      mqL3sq_indices(i) = 1 + mqL3sq_step*(i-1);
    }

  for (int i = 1; i <= mtRsq_nintvls; i++)
    {
      mtRsq_indices(i) = 1 + mtRsq_step*(i-1);
    }

  for (int i = 1; i <= At_nintvls; i++)
    {
      At_indices(i) = 1 + At_step*(i-1);
    }

  for (int i = 1; i <= M2_nintvls; i++)
    {
      M2_indices(i) = 1 + M2_step*(i-1);
    }

  // Loop over every possible combination and write to file.

  double tb_lb, lambda_lb, Alambda_lb, mqL3sq_lb, mtRsq_lb, At_lb, M2_lb;
  double tb_ub, lambda_ub, Alambda_ub, mqL3sq_ub, mtRsq_ub, At_ub, M2_ub;
  int tb_subnpts, lambda_subnpts, Alambda_subnpts, mqL3sq_subnpts, mtRsq_subnpts, At_subnpts, M2_subnpts;

  int count = 0;

  // If doing a logarithmic scan in A_t, get the logs of
  // the corresponding points. Likewise for A_lambda.
  for (int i = 1; i < tb_indices.displayEnd(); i++)
    {
  for (int j = 1; j < lambda_indices.displayEnd(); j++)
    {
  for (int k = 1; k < Alambda_indices.displayEnd(); k++)
    {
  for (int l = 1; l < mqL3sq_indices.displayEnd(); l++)
    {
  for (int m = 1; m < mtRsq_indices.displayEnd(); m++)
    {
  for (int n = 1; n < At_indices.displayEnd(); n++)
    {
  for (int p = 1; p < M2_indices.displayEnd(); p++)
    {

      count++;

      tb_lb = tb_vals(tb_indices(i));
      tb_ub = tb_vals(tb_indices(i+1)-1);
      tb_subnpts = tb_indices(i+1)-tb_indices(i);

      lambda_lb = lambda_vals(lambda_indices(j));
      lambda_ub = lambda_vals(lambda_indices(j+1)-1);
      lambda_subnpts = lambda_indices(j+1)-lambda_indices(j);

      mqL3sq_lb = mqL3sq_vals(mqL3sq_indices(l));
      mqL3sq_ub = mqL3sq_vals(mqL3sq_indices(l+1)-1);
      mqL3sq_subnpts = mqL3sq_indices(l+1)-mqL3sq_indices(l);

      mtRsq_lb = mtRsq_vals(mtRsq_indices(m));
      mtRsq_ub = mtRsq_vals(mtRsq_indices(m+1)-1);
      mtRsq_subnpts = mtRsq_indices(m+1)-mtRsq_indices(m);


      At_lb = At_vals(At_indices(n));
      At_ub = At_vals(At_indices(n+1)-1);
      At_subnpts = At_indices(n+1)-At_indices(n);

      if (atlogscan)
	{
	  // Note that it is log(A_t) that gets written
	  // to file for log scans
	  At_lb = log(At_lb);
	  At_ub = log(At_ub);
	}

      Alambda_lb = Alambda_vals(Alambda_indices(n));
      Alambda_ub = Alambda_vals(Alambda_indices(n+1)-1);
      Alambda_subnpts = Alambda_indices(n+1)-Alambda_indices(n);

      if (alambdalogscan)
	{
	  // Note that it is log(A_lambda) that gets written
	  // to file for log scans
	  Alambda_lb = log(Alambda_lb);
	  Alambda_ub = log(Alambda_ub);
	}


      M2_lb = M2_vals(M2_indices(p));
      M2_ub = M2_vals(M2_indices(p+1)-1);
      M2_subnpts = M2_indices(p+1)-M2_indices(p);


      // Write bounds to file. File names are created by appending
      // the string "_{ijklmnp}" to the given file prefix
      ostringstream convert;

      // Construct linear index for file
      int number = (p-1)+M2_nintvls*((n-1)+At_nintvls*((m-1)+mtRsq_nintvls*
						       ((l-1)+mqL3sq_nintvls*
							((k-1)+Alambda_nintvls*((j-1)+lambda_nintvls*(i-1))))));

      convert << number;

      string outfilename = outfileprefix + "_" + convert.str() + ".params";

      ofstream boundsfout(outfilename.c_str());

      writeScanBoundsFile(boundsfout, tb_lb, tb_ub, tb_subnpts, lambda_lb, lambda_ub, lambda_subnpts, 
			  Alambda_lb, Alambda_ub, Alambda_subnpts, mqL3sq_lb, mqL3sq_ub, mqL3sq_subnpts,
			  mtRsq_lb, mtRsq_ub, mtRsq_subnpts, At_lb, At_ub, At_subnpts,
			  M2_lb, M2_ub, M2_subnpts, atlogscan, alambdalogscan, yu, yc, yt, yd, ys, yb,
			  ye, ymu, ytau, g1, g2, g3, g1p, v, s);

      boundsfout.close();

    } //< M_2 loop
    } //< A_t loop
    } //< m_u3^2 loop
    } //< m_Q3^2 loop
    } //< Alambda loop
    } //< lambda loop
    } //< tan(beta) loop


  cout << "Total files written = " << count << endl;
  /*
    ----------------------------------------------------
   */

  return 0;
}

// Give it a lower and upper bound and a number of points, and it will return
// a vector of values. Optionally specify isLog = true to get logarithmically
// spaced points.
DoubleVector generateVals(double lower, double upper, int npts, bool isLogSpaced)
{
  DoubleVector vals(npts);

  double incr = 0.0;

  if (npts < 1)
    {
      ostringstream ii;
      ii << "ERROR: in generateVals: number of values must be at least 1.\n";
      throw ii.str();
    }

  if (!isLogSpaced)
    {
      if (npts != 1)
	{
	  incr = (upper-lower)/(((double)npts)-1.0);
	}

      for (int i = 1; i <= npts; i++)
	{
	  vals(i) = lower + (((double)i)-1.0)*incr;
	}

    }
  else
    {
      // log scan = linear scan over log(vals),
      // but note the restriction on lower and upper
      // in this case: both should be greater than zero.

      if (lower <= 0.0 || upper <= 0.0)
	{
	  ostringstream ii;
	  ii << "ERROR: in generateVals: lower and upper bounds must be greater than zero for log scan.\n";
	  throw ii.str();
	}

      double loglwr = log(fabs(lower));
      double logupr = log(fabs(upper));

      if (npts != 1)
	{
	  incr = (logupr-loglwr)/(((double)npts)-1.0);
	}

      for (int i = 1; i <= npts; i++)
	{
	  vals(i) = exp(loglwr + (((double)i)-1.0)*incr);
	}

    }

  return vals;
}

void writeScanBoundsFile(ofstream & boundsfout, double tb_lb, double tb_ub, int tb_subnpts, 
			 double lambda_lb, double lambda_ub, int lambda_subnpts, 
			 double Alambda_lb, double Alambda_ub, int Alambda_subnpts, 
			 double mqL3sq_lb, double mqL3sq_ub, int mqL3sq_subnpts,
			 double mtRsq_lb, double mtRsq_ub, int mtRsq_subnpts, 
			 double At_lb, double At_ub, int At_subnpts,
			 double M2_lb, double M2_ub, int M2_subnpts, bool atlogscan, bool alambdalogscan,
			 double yu, double yc, double yt, double yd, double ys, double yb,
			 double ye, double ymu, double ytau, double g1, double g2, double g3, double g1p,
			 double v, double s)
{

  int outprecis = boundsfout.precision();

  boundsfout << "# Scan parameters file suitable for input into pmssmScanner" << endl;
  boundsfout << "# Summary: " << endl;
  boundsfout << "# lower tan(beta) = " << tb_lb << ", upper tan(beta) = " << tb_ub;
  boundsfout << ", points in sub-interval = " << tb_subnpts << endl;
  boundsfout << "# lower lambda = " << lambda_lb << ", upper lambda = " << lambda_ub;
  boundsfout << ", points in sub-interval = " << lambda_subnpts << endl;
  if (alambdalogscan)
    {
      boundsfout << "# A_lambda is log scanned" << endl;
      boundsfout << "# lower log(A_lambda) = " << Alambda_lb << ", upper log(A_lambda) = " << Alambda_ub;
      boundsfout << ", points in sub-interval = " << Alambda_subnpts << endl;
    }
  else
    {
      boundsfout << "# lower A_lambda = " << Alambda_lb << ", upper A_lambda = " << Alambda_ub;
      boundsfout << ", points in sub-interval = " << Alambda_subnpts << endl;
    }
  boundsfout << "# lower m_Q3^2 = " << mqL3sq_lb << ", upper m_Q3^2 = " << mqL3sq_ub;
  boundsfout << ", points in sub-interval = " << mqL3sq_subnpts << endl;
  boundsfout << "# lower m_u3^2 = " << mtRsq_lb << ", upper m_u3^2 = " << mtRsq_ub;
  boundsfout << ", points in sub-interval = " << mtRsq_subnpts << endl;
  if (atlogscan)
    {
      boundsfout << "# A_t is log scanned" << endl;
      boundsfout << "# lower log(A_t) = " << At_lb << ", upper log(A_t) = " << At_ub;
      boundsfout << ", points in sub-interval = " << At_subnpts << endl;
    }
  else
    {
      boundsfout << "# lower A_t = " << At_lb << ", upper A_t = " << At_ub;
      boundsfout << ", points in sub-interval = " << At_subnpts << endl;
    }
  boundsfout << "# lower M_2 = " << M2_lb << ", upper M_2 = " << M2_ub;
  boundsfout << ", points in sub-interval = " << M2_subnpts << endl;

  boundsfout.setf(ios::scientific, ios::floatfield);
  boundsfout.precision(8);

  boundsfout << "TBLL=" << tb_lb << endl;
  boundsfout << "TBUL=" << tb_ub << endl;
  boundsfout << "TBNPTS=" << tb_subnpts << endl;
  boundsfout << "LAMBDALL=" << lambda_lb << endl;
  boundsfout << "LAMBDAUL=" << lambda_ub << endl;
  boundsfout << "LAMBDANPTS=" << lambda_subnpts << endl;
  if (alambdalogscan)
    {
      boundsfout << "LOGALAMBDALL=" << Alambda_lb << endl;
      boundsfout << "LOGALAMBDAUL=" << Alambda_ub << endl;
      boundsfout << "LOGALAMBDANPTS=" << Alambda_subnpts << endl;
    }
  else
    {
      boundsfout << "ALAMBDALL=" << Alambda_lb << endl;
      boundsfout << "ALAMBDAUL=" << Alambda_ub << endl;
      boundsfout << "ALAMBDANPTS=" << Alambda_subnpts << endl;
    }
  boundsfout << "MQLSQLL=" << mqL3sq_lb << endl;
  boundsfout << "MQLSQUL=" << mqL3sq_ub << endl;
  boundsfout << "MQLSQNPTS=" << mqL3sq_subnpts << endl;
  boundsfout << "MURSQLL=" << mtRsq_lb << endl;
  boundsfout << "MURSQUL=" << mtRsq_ub << endl;
  boundsfout << "MURSQNPTS=" << mtRsq_subnpts << endl;
  if (atlogscan)
    {
      boundsfout << "LOGATLL=" << At_lb << endl;
      boundsfout << "LOGATUL=" << At_ub << endl;
      boundsfout << "LOGATNPTS=" << At_subnpts << endl;
    }
  else
    {
      boundsfout << "ATLL=" << At_lb << endl;
      boundsfout << "ATUL=" << At_ub << endl;
      boundsfout << "ATNPTS=" << At_subnpts << endl;
    }
  boundsfout << "M2LL=" << M2_lb << endl;
  boundsfout << "M2UL=" << M2_ub << endl;
  boundsfout << "M2NPTS=" << M2_subnpts << endl;

  boundsfout << "YU=" << yu << endl;
  boundsfout << "YC=" << yc << endl;
  boundsfout << "YT=" << yt << endl;

  boundsfout << "YD=" << yd << endl;
  boundsfout << "YS=" << ys << endl;
  boundsfout << "YB=" << yb << endl;

  boundsfout << "YE=" << ye << endl;
  boundsfout << "YMU=" << ymu << endl;
  boundsfout << "YTAU=" << ytau << endl;

  boundsfout << "G1=" << g1 << endl;
  boundsfout << "G2=" << g2 << endl;
  boundsfout << "G3=" << g3 << endl;
  boundsfout << "GN=" << g1p << endl;

  boundsfout << "HVEV=" << v << endl;
  boundsfout << "SVEV=" << s << endl;

  boundsfout.precision(outprecis);

}
