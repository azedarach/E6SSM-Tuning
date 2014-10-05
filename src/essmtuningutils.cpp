/*
  Implementations of the functions defined in essmtuningutils.h
 */

#include "essmtuningutils.h"

namespace essm_tuning_utils {

using namespace flexiblesusy;

// A function to calculate the 3x3 matrix appearing on the LHS of the calculation
// of the derivatives of the VEVs wrt the parameters in the E6SSM. Note that this
// assumes that the E6SSM object is already set to have its renormalisation
// scale set at M_SUSY.
  Eigen::Matrix<double,3,3> doCalcLHSTuningMatrix(genericE6SSM_soft_parameters essmSusy, Eigen::VectorXd const & vevs)
{
  
  Eigen::Matrix<double,3,3> lhsMat;
  
  if (vevs.size() != 3)
    {
      cerr << "WARNING: incorrect number of VEVs supplied to function: skipping calculating LHS matrix." << endl;
    }
  else
    {
      
      double v1 = vevs(0);
      double v2 = vevs(1);
      double s = vevs(2);
      double tb = v2/v1;
      double v = Sqrt(v1*v1+v2*v2);

      double g1 = essmSusy.get_g1();
      double g2 = essmSusy.get_g2();
      double gbar = Sqrt(g2*g2+0.6*g1*g1);
      double gdash_1 = essmSusy.get_gN();

      double cb = 1.0/Sqrt(1.0+tb*tb);
      double sb = tb*cb;
      double c2b = (1.0-tb*tb)/(1.0+tb*tb);

      genericE6SSM_input_parameters input = essmSusy.get_input();
      
      // Because we are neglecting U(1) mixing, 
      // we will for now approximate the effective charges
      // by the U(1)' charges
      double Qtilde_1 = input.QH1p;
      double Qtilde_2 = input.QH2p;
      double Qtilde_s = input.QSp;
      
      
      double Tlambda = essmSusy.get_TLambdax();
      double lambda = essmSusy.get_Lambdax();
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
      
      double df1dv1 = lambda*Alambda*s*v2/(Sqrt(2.0)*v1*v1)+gbar*gbar*v1/4.0+gdash_1*gdash_1*Qtilde_1*Qtilde_1*v1;
      double df1dv2 = lambda*lambda*v2-lambda*Alambda*s/(Sqrt(2.0)*v1)-gbar*gbar*v2/4.0+gdash_1*gdash_1*Qtilde_1*Qtilde_2*v2;
      double df1ds = lambda*lambda*s-lambda*Alambda*v2/(Sqrt(2.0)*v1)+gdash_1*gdash_1*Qtilde_1*Qtilde_s*s;
      double df2dv1 = lambda*lambda*v1-lambda*Alambda*s/(Sqrt(2.0)*v2)-gbar*gbar*v1/4.0+gdash_1*gdash_1*Qtilde_2*Qtilde_1*v1;
      double df2dv2 = lambda*Alambda*s*v1/(Sqrt(2.0)*v2*v2)+gbar*gbar*v2/4.0+gdash_1*gdash_1*Qtilde_2*Qtilde_2*v2;
      double df2ds = lambda*lambda*s-lambda*Alambda*v1/(Sqrt(2.0)*v2)+gdash_1*gdash_1*Qtilde_2*Qtilde_s*s;
      double df3dv1 = lambda*lambda*v1-lambda*Alambda*v2/(Sqrt(2.0)*s)+gdash_1*gdash_1*Qtilde_s*Qtilde_1*v1;
      double df3dv2 = lambda*lambda*v2-lambda*Alambda*v1/(Sqrt(2.0)*s)+gdash_1*gdash_1*Qtilde_s*Qtilde_2*v2;
      double df3ds = lambda*Alambda*v1*v2/(Sqrt(2.0)*s*s)+gdash_1*gdash_1*Qtilde_s*Qtilde_s*s;

      // Add in 1-loop corrections if requested
      if (INCLUDE1LPTADPOLES)
	{
	  // NB Peter's code defines the tadpoles with the opposite sign.
	  df1dv1 = df1dv1 + (1.0/v1)*(doCalcDeltaPrime11(essmSusy, s, tb) + doCalcTadpoleESSMH1(essmSusy, s, tb));
	  df1dv2 = df1dv2 + doCalcDeltaPrime12(essmSusy, s, tb)/v1;
	  df1ds = df1ds + doCalcDeltaPrime13(essmSusy, s, tb)/v1;

	  df2dv1 = df2dv1 + doCalcDeltaPrime12(essmSusy, s, tb)/v2;
	  df2dv2 = df2dv2 + (1.0/v2)*(doCalcDeltaPrime22(essmSusy, s, tb)+ doCalcTadpoleESSMH2(essmSusy, s, tb));
	  df2ds = df2ds + doCalcDeltaPrime23(essmSusy, s, tb)/v2;

	  df3dv1 = df3dv1 + doCalcDeltaPrime13(essmSusy, s, tb)/s;
	  df3dv2 = df3dv2 + doCalcDeltaPrime23(essmSusy, s, tb)/s;
	  df3ds = df3ds + (1.0/s)*(doCalcDeltaPrime33(essmSusy, s, tb) + doCalcTadpolesESSMS(essmSusy, s, tb));
	}      
      
            
      lhsMat(0,0) = df1dv1;
      lhsMat(0,1) = df1dv2;
      lhsMat(0,2) = df1ds;
      lhsMat(1,0) = df2dv1;
      lhsMat(1,1) = df2dv2;
      lhsMat(1,2) = df2ds;
      lhsMat(2,0) = df3dv1;
      lhsMat(2,1) = df3dv2;
      lhsMat(2,2) = df3ds;

    }
  
  return lhsMat;

}

// A function to calculate d log(M_Z^2)/d log(p) using the value of the parameter p
// and the already calculated derivatives d v_i/ d p_j. For the 
  double doCalcdLogMzSqdLogParam(genericE6SSM_soft_parameters r, double p, Eigen::VectorXd const & vevs, 
			       Eigen::VectorXd const & dVevsdp)
{

  double deriv = 0.0;

  if ((vevs.size() != 3) || (dVevsdp.size() != 3))
    {
      cerr << "WARNING: incorrect number of VEVs supplied to function: skipping calculating derivative." << endl;
    }
  else
    {
      double g1 = r.get_g1();
      double g2 = r.get_g2();
      double gbar = Sqrt(g2*g2+0.6*g1*g1);
      
      double v1 = vevs(0);
      double v2 = vevs(1);

      double v = Sqrt(v1*v1+v2*v2);

      // Neglect any neutral mixing for now, but later may want to include contribution from it.      
      deriv = (2.0*p/(v*v))*(v1*dVevsdp(0)+v2*dVevsdp(1));
    }

  return deriv;

}

// ESSM_EWSBConditioni, i = 1, 2, calculates the value of the 
// EWSB condition for m_H_i^2 and m_S^2.
// Inputs:
//     SoftParsMssm const & r = pMSSM model object
//     DoubleVector const & vevs = the values of v1 and v2
double ESSM_EWSBCondition1(genericE6SSM_soft_parameters const & r)
{

  double m1Sq = r.get_mHd2();

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

  double gdash = Sqrt(3.0/5.0)*r.get_g1();
  double g_2 = r.get_g2();
  double gbar = Sqrt(gdash*gdash+g_2*g_2);
  double gdash_1 = r.get_gN();

  genericE6SSM_input_parameters input = r.get_input();

  // Because we are neglecting U(1) mixing, 
  // we will for now approximate the effective charges
  // by the U(1)' charges
  double Qtilde_1 = input.QH1p;
  double Qtilde_2 = input.QH2p;
  double Qtilde_s = input.QSp;

  double v2 = r.get_vu();
  double v2Sq = Sqr(v2);
  double v1 = r.get_vd();
  double s = r.get_vs();

  double v = Sqrt(v2*v2+v1*v1);
  double tb = v2/v1;

  double cSqb = 1.0/(1.0+tb*tb);
  double sSqb = (tb*tb)/(1.0+tb*tb);
  double c2b = (1.0-tb*tb)/(1.0+tb*tb);

  double f1;

  if (!INCLUDE1LPTADPOLES)
    {
      f1 = m1Sq + 0.5*lambda*lambda*(v2Sq+s*s)-lambda*Alambda*s*tb/Sqrt(2.0)
	+gbar*gbar*v*v*c2b/8.0 + (0.5*Qtilde_1*gdash_1*gdash_1*
				  (Qtilde_1*v*v*cSqb + Qtilde_2*v*v*sSqb + Qtilde_s*s*s));
    }
  else
    {
      f1 = m1Sq + 0.5*lambda*lambda*(v2Sq+s*s)-lambda*Alambda*s*tb/Sqrt(2.0)
	+gbar*gbar*v*v*c2b/8.0 + (0.5*Qtilde_1*gdash_1*gdash_1*
				  (Qtilde_1*v*v*cSqb + Qtilde_2*v*v*sSqb + Qtilde_s*s*s))
	-doCalcTadpoleESSMH1(r, s, tb);
    }

  return f1;

}

double ESSM_EWSBCondition2(genericE6SSM_soft_parameters const & r)
{

  double m2Sq = r.get_mHu2();

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

  double gdash = Sqrt(3.0/5.0)*r.get_g1();
  double g_2 = r.get_g2();
  double gbar = Sqrt(gdash*gdash+g_2*g_2);
  double gdash_1 = r.get_gN();

  genericE6SSM_input_parameters input = r.get_input();

  // Because we are neglecting U(1) mixing, 
  // we will for now approximate the effective charges
  // by the U(1)' charges
  double Qtilde_1 = input.QH1p;
  double Qtilde_2 = input.QH2p;
  double Qtilde_s = input.QSp;

  double v1 = r.get_vd();
  double v2 = r.get_vu();
  double s = r.get_vs();

  double tb = v2/v1;
  double cSqb = 1.0/(1.0+tb*tb);
  double sSqb = (tb*tb)/(1.0+tb*tb);
  double c2b = (1.0-tb*tb)/(1.0+tb*tb);

  double v = Sqrt(v1*v1+v2*v2);
  double v1Sq = v1*v1;

  double f2;
  
  if (!INCLUDE1LPTADPOLES)
    {
      f2 = m2Sq + 0.5*lambda*lambda*(v1Sq+s*s)-lambda*Alambda*s/(Sqrt(2.0)*tb)
	-gbar*gbar*v*v*c2b/8.0 + (0.5*Qtilde_2*gdash_1*gdash_1*
				  (Qtilde_1*v*v*cSqb + Qtilde_2*v*v*sSqb + Qtilde_s*s*s));
    }
  else
    {
      f2 = m2Sq + 0.5*lambda*lambda*(v1Sq+s*s)-lambda*Alambda*s/(Sqrt(2.0)*tb)
	-gbar*gbar*v*v*c2b/8.0 + (0.5*Qtilde_2*gdash_1*gdash_1*
				  (Qtilde_1*v*v*cSqb + Qtilde_2*v*v*sSqb + Qtilde_s*s*s))
	-doCalcTadpoleESSMH2(r, s, tb);
    }


  return f2;
}

double ESSM_EWSBCondition3(genericE6SSM_soft_parameters const & r)
{

  double msSq = r.get_ms2();

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


  double gdash = Sqrt(3.0/5.0)*r.get_g1();
  double g_2 = r.get_g2();
  double gbar = Sqrt(gdash*gdash+g_2*g_2);
  double gdash_1 = r.get_gN();

  genericE6SSM_input_parameters input = r.get_input();

  // Because we are neglecting U(1) mixing, 
  // we will for now approximate the effective charges
  // by the U(1)' charges
  double Qtilde_1 = input.QH1p;
  double Qtilde_2 = input.QH2p;
  double Qtilde_s = input.QSp;

  double v1 = r.get_vd();
  double v2 = r.get_vu();
  double s = r.get_vs();

  double tb = v2/v1;
  double cSqb = 1.0/(1.0+tb*tb);
  double sSqb = (tb*tb)/(1.0+tb*tb);
  double c2b = (1.0-tb*tb)/(1.0+tb*tb);
  double s2b = 2.0*tb/(1.0+tb*tb);

  double v = Sqrt(v1*v1+v2*v2);

  double f3;

  if (!INCLUDE1LPTADPOLES)
    {
      f3 = msSq + 0.5*lambda*lambda*v*v-(lambda*Alambda*v1*v2/(Sqrt(2.0)*s))
	+0.5*Qtilde_s*gdash_1*gdash_1*(Qtilde_1*v*v*cSqb + Qtilde_2*v*v*sSqb + Qtilde_s*s*s);
    }
  else
    {
      f3 = msSq + 0.5*lambda*lambda*v*v-lambda*Alambda*v1*v2/(Sqrt(2.0)*s)
	+0.5*Qtilde_s*gdash_1*gdash_1*(Qtilde_1*v*v*cSqb + Qtilde_2*v*v*sSqb + Qtilde_s*s*s)
	-doCalcTadpolesESSMS(r, s, tb);
    }

  return f3;

}

// Implement the EWSB conditions by varying the values of m_Hu^2 and m_Hd^2 at MX.
// Scale factors needed for Newton's method implementation. 
static double mHuSqScaleFactor = 1.0e6;
static double mHdSqScaleFactor = 1.0e6;
static double mSSqScaleFactor = 1.0e6;
static double MsusyScaleFactor = 2.0e3;
bool ESSM_ImplementEWSBConstraints_SoftMasses(genericE6SSM_soft_parameters r, double mx, double ms, bool useMxEqualsMs,
					      DoubleVector & updatedSoln, double tol)
{
  // We have to use a shooting method to do this in general if MX is different from M_{SUSY}. 
  // We should start with an object given at the scale MX.
  double q = r.get_scale();

  bool hasProblem = false;

  

  if (fabs(q-mx) > TOLERANCE)
    {
  
      ostringstream ii;
      ii << "WARNING: trying to calculate required soft masses at inappropriate scale q = " << q << " GeV,\n"
	 << "         instead of at MX = " << mx << " GeV: exiting." << endl;
      throw ii.str();
    }
  else
    {
     
      // If MX = M_{SUSY} there is no need to use the shooting approach; just
      // rearrange the EWSB conditions to solve for m_Hu^2 and m_Hd^2.
      if (useMxEqualsMs)
	{
	  if (fabs(mx-ms) > TOLERANCE || fabs(q-ms) > TOLERANCE)
	    {
	      cerr << "WARNING: MX = M_{SUSY} requested but given values do not agree. Assuming M_{SUSY} = " << ms << "." << endl;
	      r.run_to(ms, PRECISION);
	    }

	  if (ENABLE_DEBUG)
	    {
	      cout << "#    setting MSUSY = MX" << endl;
	    }

	  // Solve for the updated soft Higgs masses

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

	  double mSSq, mH1Sq, mH2Sq;
	  double g1 = r.get_g1();
	  double g2 = r.get_g2();
	  double gbar = Sqrt(g2*g2+0.6*g1*g1);
	  double gdash_1 = r.get_gN();	  

	  double v1 = r.get_vd();
	  double v2 = r.get_vu();
	  double s = r.get_vs();
	  double v = Sqrt(v1*v1+v2*v2);
	  double tb = v2/v1;  

	  double cSqb = 1.0/(1.0+tb*tb);
	  double sSqb = (tb*tb)/(1.0+tb*tb);
	  double c2b = (1.0-tb*tb)/(1.0+tb*tb);
	  double s2b = 2.0*tb/(1.0+tb*tb);	  

	  genericE6SSM_input_parameters input = r.get_input();
	  
	  // Because we are neglecting U(1) mixing, 
	  // we will for now approximate the effective charges
	  // by the U(1)' charges
	  double Qtilde_1 = input.QH1p;
	  double Qtilde_2 = input.QH2p;
	  double Qtilde_s = input.QSp;



	  if (!INCLUDE1LPTADPOLES)
	    {

	      if (ENABLE_DEBUG)
		{
		  cout << "#    using tree level tadpole equations" << endl;
		}

	      mSSq = (1.0/s)*(lambda*Alambda*v*v*s2b/(2.0*Sqrt(2.0))-0.5*lambda*lambda*v*v*s
			      -0.5*gdash_1*gdash_1*Qtilde_s*s*(Qtilde_1*v*v*cSqb+Qtilde_2*v*v*sSqb+Qtilde_s*s*s));
	      mH1Sq = lambda*Alambda*s*tb/Sqrt(2.0)-0.5*lambda*lambda*(s*s+v*v*sSqb)-gbar*gbar*v*v*c2b/8.0
		-0.5*Qtilde_1*gdash_1*gdash_1*(Qtilde_1*v*v*cSqb+Qtilde_2*v*v*sSqb+Qtilde_s*s*s);
	      mH2Sq = lambda*Alambda*s/(Sqrt(2.0)*tb)-0.5*lambda*lambda*(s*s+v*v*cSqb)+gbar*gbar*v*v*c2b/8.0
		-0.5*Qtilde_2*gdash_1*gdash_1*(Qtilde_1*v*v*cSqb+Qtilde_2*v*v*sSqb+Qtilde_s*s*s);

	    }
	  else
	    {

	      if (ENABLE_DEBUG)
		{
		  cout << "#    using tree level tadpole equations + 1-loop tadpoles" << endl;
		}

	      // Calculate the values of the soft squared masses using the three EWSB
	      // conditions at tree level.
	      mSSq = (1.0/s)*(lambda*Alambda*v*v*s2b/(2.0*Sqrt(2.0))-0.5*lambda*lambda*v*v*s
			      -0.5*gdash_1*gdash_1*Qtilde_s*s*(Qtilde_1*v*v*cSqb+Qtilde_2*v*v*sSqb+Qtilde_s*s*s));
	      mH1Sq = lambda*Alambda*s*tb/Sqrt(2.0)-0.5*lambda*lambda*(s*s+v*v*sSqb)-gbar*gbar*v*v*c2b/8.0
		-0.5*Qtilde_1*gdash_1*gdash_1*(Qtilde_1*v*v*cSqb+Qtilde_2*v*v*sSqb+Qtilde_s*s*s);
	      mH2Sq = lambda*Alambda*s/(Sqrt(2.0)*tb)-0.5*lambda*lambda*(s*s+v*v*cSqb)+gbar*gbar*v*v*c2b/8.0
		-0.5*Qtilde_2*gdash_1*gdash_1*(Qtilde_1*v*v*cSqb+Qtilde_2*v*v*sSqb+Qtilde_s*s*s);

	      // Add in appropriate loop corrections (note sign!).
	      mSSq = mSSq + doCalcTadpolesESSMS(r, s, tb);
	      mH1Sq = mH1Sq + doCalcTadpoleESSMH1(r, s, tb);
	      mH2Sq = mH2Sq + doCalcTadpoleESSMH2(r, s, tb);

	    }
	  updatedSoln(1) = mH1Sq;
	  updatedSoln(2) = mH2Sq;
	  updatedSoln(3) = mSSq;
	  updatedSoln(4) = ms;
	}
      else
	{
	  if (ENABLE_DEBUG)
	    {
	      cout << "#    iterating to get soft masses" << endl;
	    }

	  // Otherwise we have to use a shooting method to determine the values
	  // of m_Hu^2, m_Hd^2 and m_S^2 at MX that satisfy the EWSB conditions at M_{SUSY}.

	  // Guess initial values of m_Hu^2, m_Hd^2 and m_S^2 at MX that solve the EWSB conditions,
	  // and guess the initial scale M_{SUSY}.

	  // The initial guess for M_{SUSY} is pretty obvious - it is the provided value ms.
	  double mSusy_guess = ms;

	  MsusyScaleFactor = ms;

	  // To guess the appropriate initial values of m_Hu^2(MX), m_Hd^2(MX) and m_S^2, we use 
	  // our approximate Taylor series solution to the RGEs. This should be reasonable
	  // for small values of log(MX/M_{SUSY}). Note that this assumes only non-zero
	  // couplings are third generation.
	  int nLps = r.get_loops();
	  if (nLps != 1 && nLps != 2)
	    {
	      cerr << "WARNING: requested RG running at " << nLps << " loop order: can only do 1 or 2 loop running:" << endl;
	      cerr << "         assuming 2 loop running." << endl;
	      nLps = 2;
	    }

	  if (ENABLE_DEBUG)
	    {
	      cout << "#    constructing initial guess";   
	    }

	  double t_run = log(mx/ms);

	  double mHuSq_guess, mHdSq_guess, mSSq_guess;

	  // If t is too large, this is an unreliable estimate, so just use a default value.
	  const int TLIMIT = 10;
	  if (t_run >= TLIMIT)
	    {
	      if (ENABLE_DEBUG)
		{
		  cout << " using default values" << endl;
		}
	      mHdSqScaleFactor = r.get_mHd2();
	      mHuSqScaleFactor = r.get_mHu2();
	      mSSqScaleFactor = r.get_ms2();

	      mHuSq_guess = mHuSqScaleFactor;
	      mHdSq_guess = mHdSqScaleFactor;
	      mSSq_guess = mSSqScaleFactor;
	    }
	  else
	    {
	      if (ENABLE_DEBUG)
		{
		  cout << " using approximate RGE solutions" << endl;
		}

	      // The Taylor series expansion is about the solution at M_{SUSY} initially.
	      genericE6SSM_soft_parameters w = r; // copy object to avoid numerical errors
	      w.run_to(ms, PRECISION);
	      
	      double mHuSqInit, mHdSqInit, mSSqInit;
	     
	      double Tlambda = w.get_TLambdax();
	      double lambda = w.get_Lambdax();
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
	      
	      double g1 = w.get_g1();
	      double g2 = w.get_g2();
	      double gbar = Sqrt(g2*g2+0.6*g1*g1);
	      double gdash_1 = w.get_gN();	      

	      double v1 = w.get_vd();
	      double v2 = w.get_vu();
	      double s = w.get_vs();
	      double v = Sqrt(v1*v1+v2*v2);
	      double tb = v2/v1;
  	      
	      double cSqb = 1.0/(1.0+tb*tb);
	      double sSqb = (tb*tb)/(1.0+tb*tb);
	      double c2b = (1.0-tb*tb)/(1.0+tb*tb);
	      double s2b = 2.0*tb/(1.0+tb*tb);	  
	      
	      genericE6SSM_input_parameters input = w.get_input();
	      
	      // Because we are neglecting U(1) mixing, 
	      // we will for now approximate the effective charges
	      // by the U(1)' charges
	      double Qtilde_1 = input.QH1p;
	      double Qtilde_2 = input.QH2p;
	      double Qtilde_s = input.QSp;

	      if (!INCLUDE1LPTADPOLES)
		{
		  mSSqInit = (1.0/s)*(lambda*Alambda*v*v*s2b/(2.0*Sqrt(2.0))-0.5*lambda*lambda*v*v*s
				      -0.5*gdash_1*gdash_1*Qtilde_s*s*(Qtilde_1*v*v*cSqb+Qtilde_2*v*v*sSqb+Qtilde_s*s*s));
		  mHdSqInit = lambda*Alambda*s*tb/Sqrt(2.0)-0.5*lambda*lambda*(s*s+v*v*sSqb)-gbar*gbar*v*v*c2b/8.0
		    -0.5*Qtilde_1*gdash_1*gdash_1*(Qtilde_1*v*v*cSqb+Qtilde_2*v*v*sSqb+Qtilde_s*s*s);
		  mHuSqInit = lambda*Alambda*s/(Sqrt(2.0)*tb)-0.5*lambda*lambda*(s*s+v*v*cSqb)+gbar*gbar*v*v*c2b/8.0
		    -0.5*Qtilde_2*gdash_1*gdash_1*(Qtilde_1*v*v*cSqb+Qtilde_2*v*v*sSqb+Qtilde_s*s*s);
		}
	      else
		{
		  // Calculate the values of the soft squared masses using the three EWSB
		  // conditions at tree level.
		  mSSqInit = (1.0/s)*(lambda*Alambda*v*v*s2b/(2.0*Sqrt(2.0))-0.5*lambda*lambda*v*v*s
				      -0.5*gdash_1*gdash_1*Qtilde_s*s*(Qtilde_1*v*v*cSqb+Qtilde_2*v*v*sSqb+Qtilde_s*s*s));
		  mHdSqInit = lambda*Alambda*s*tb/Sqrt(2.0)-0.5*lambda*lambda*(s*s+v*v*sSqb)-gbar*gbar*v*v*c2b/8.0
		    -0.5*Qtilde_1*gdash_1*gdash_1*(Qtilde_1*v*v*cSqb+Qtilde_2*v*v*sSqb+Qtilde_s*s*s);
		  mHuSqInit = lambda*Alambda*s/(Sqrt(2.0)*tb)-0.5*lambda*lambda*(s*s+v*v*cSqb)+gbar*gbar*v*v*c2b/8.0
		    -0.5*Qtilde_2*gdash_1*gdash_1*(Qtilde_1*v*v*cSqb+Qtilde_2*v*v*sSqb+Qtilde_s*s*s);
		  
		  // Add in appropriate loop corrections (note sign!).
		  mSSqInit = mSSqInit + doCalcTadpolesESSMS(w, s, tb);
		  mHdSqInit = mHdSqInit + doCalcTadpoleESSMH1(w, s, tb);
		  mHuSqInit = mHuSqInit + doCalcTadpoleESSMH2(w, s, tb);
		}

	      if (ENABLE_DEBUG)
		{
		  cout << "#    got tree approximation: m_Hd^2 = " << mHdSqInit
		       << ", m_Hu^2 = " << mHuSqInit << ", m_s^2 = " << mSSqInit << endl;
		  cout << "#    constructing log coefficients" << endl;
		}	      

	      // Construct required coefficients in Taylor series.
	      double mHdSqLogCoeff = doCalcMh1SquaredLogCoeff(w, nLps);
	      double mHdSqLogSqCoeff = doCalcMh1SquaredLogSqCoeff(w, 1); 
	      double mHuSqLogCoeff = doCalcMh2SquaredLogCoeff(w, nLps);
	      double mHuSqLogSqCoeff = doCalcMh2SquaredLogSqCoeff(w, 1); 
	      double mSSqLogCoeff = doCalcMsSquaredLogCoeff(w, nLps);
	      double mSSqLogSqCoeff = doCalcMsSquaredLogSqCoeff(w, 1); 

	      if (ENABLE_DEBUG)
		{
		  cout << "#    log coefficient m_Hd^2 = " << mHdSqLogCoeff << endl;
		  cout << "#    log^2 coefficient m_Hd^2 = " << mHdSqLogSqCoeff << endl;
		  cout << "#    log coefficient m_Hu^2 = " << mHuSqLogCoeff << endl;
		  cout << "#    log^2 coefficient m_Hu^2 = " << mHuSqLogSqCoeff << endl;
		  cout << "#    log coefficient m_s^2 = " << mSSqLogCoeff << endl;
		  cout << "#    log^2 coefficient m_s^2 = " << mSSqLogSqCoeff << endl;
		}
	      
	      mHuSq_guess = mHuSqInit + t_run*mHuSqLogCoeff + Sqr(t_run)*mHuSqLogSqCoeff;
	      mHdSq_guess = mHdSqInit + t_run*mHdSqLogCoeff + Sqr(t_run)*mHdSqLogSqCoeff;
	      mSSq_guess = mSSqInit + t_run*mSSqLogCoeff + Sqr(t_run)*mSSqLogSqCoeff;
	      
	      mHuSqScaleFactor = mHuSq_guess;
	      mHdSqScaleFactor = mHdSq_guess;
	      mSSqScaleFactor = mSSq_guess;

	    }

	  if (ENABLE_DEBUG)
	    {
	      cout << "#    initial guess: m_Hd^2 = " << mHdSq_guess
		   << ", m_Hu^2 = " << mHuSq_guess 
		   << ", m_s^2 = " << mSSq_guess
		   << ", MSUSY = " << mSusy_guess << endl;
	    }

	  // Then shoot using Newton's method to try to get the actual solution. 
	  DoubleVector solEstimate(4);

	  solEstimate(1) = mHdSq_guess/mHdSqScaleFactor;
	  solEstimate(2) = mHuSq_guess/mHuSqScaleFactor;
	  solEstimate(3) = mSSq_guess/mSSqScaleFactor;
	  solEstimate(4) = mSusy_guess/MsusyScaleFactor;

	  hasProblem = ESSM_EWSB_NewtonShooter(r, solEstimate, tol);

	  updatedSoln(1) = solEstimate(1)*mHdSqScaleFactor;
	  updatedSoln(2) = solEstimate(2)*mHuSqScaleFactor;
	  updatedSoln(3) = solEstimate(3)*mSSqScaleFactor;
	  updatedSoln(4) = solEstimate(4)*MsusyScaleFactor;

	}

    }

  if (ENABLE_DEBUG)
    {
      cout << "#    found solution: m_Hd^2 = " << updatedSoln(1)
	   << ", m_Hu^2 = " << updatedSoln(2) 
	   << ", m_s^2 = " << updatedSoln(3)
	   << ", MSUSY = " << updatedSoln(4) << endl;
    }

  return hasProblem;

}

double ESSM_Msusy_Cond(genericE6SSM_soft_parameters r, double ms)
{
  DoubleVector mstop(2), mstopsq(2), mD1sq(3), mD2sq(3);

  double v1 = r.get_vd();
  double v2 = r.get_vu();
  double s = r.get_vs();
  double tb = v2/v1;

  physical_ESSM(r, mstop, mstopsq, mD1sq, mD2sq, s, tb);

  double f = (ms/Sqrt(mstop(1)*mstop(2)))-1.0;

  return f;
}

static genericE6SSM_soft_parameters *tempsoftessm;
void ESSM_EWSB_Shooter_Functions(DoubleVector const & parVals, DoubleVector & f)
{
  genericE6SSM_soft_parameters s = *tempsoftessm;

  s.set_mHd2(parVals(1)*mHdSqScaleFactor);
  s.set_mHu2(parVals(2)*mHuSqScaleFactor);
  s.set_ms2(parVals(3)*mSSqScaleFactor);

  // DH:NOTE
  cout << "lambda = " << s.get_Lambdax() << endl;
  cout << "Alambda = " << s.get_TLambdax()/s.get_Lambdax() << endl;
  cout << "At = " << s.get_TYu(2,2)/s.get_Yu(2,2) << endl;
  cout << "mHd2 = " << s.get_mHd2() << endl;
  cout << "mHu2 = " << s.get_mHu2() << endl;
  cout << "ms2 = " << s.get_ms2() << endl;

  s.run_to(parVals(4)*MsusyScaleFactor, PRECISION);

  cout << "msusy = " << s.get_scale() << endl;

  f(1) = ESSM_EWSBCondition1(s)/mHdSqScaleFactor;
  f(2) = ESSM_EWSBCondition2(s)/mHuSqScaleFactor;
  f(3) = ESSM_EWSBCondition3(s)/mSSqScaleFactor;
  f(4) = ESSM_Msusy_Cond(s, parVals(4)*MsusyScaleFactor);

}

bool ESSM_EWSB_NewtonShooter(genericE6SSM_soft_parameters const & r, DoubleVector & estimate, double tol)
{
  // Copy object to actually do the running on
  genericE6SSM_soft_parameters w = r;

  // Set initial guess values
  double ms_guess = estimate(4)*MsusyScaleFactor;

  w.set_mHd2(estimate(1)*mHdSqScaleFactor);
  w.set_mHu2(estimate(2)*mHuSqScaleFactor);
  w.set_ms2(estimate(3)*mSSqScaleFactor);

  DoubleVector fVals(4);

  tempsoftessm = &w;

  // Now shoot to try to get solutions. Use globally convergent
  // Newton's method as the root finder.
  bool hasProblem = newt(estimate, ESSM_EWSB_Shooter_Functions);

  w.set_mHd2(estimate(1)*mHdSqScaleFactor);
  w.set_mHu2(estimate(2)*mHuSqScaleFactor);
  w.set_ms2(estimate(3)*mSSqScaleFactor);

  w.run_to(estimate(4)*MsusyScaleFactor, PRECISION);

  DoubleVector mstop(2), mstopsq(2), mD1sq(3), mD2sq(3);
  double v1 = w.get_vd();
  double v2 = w.get_vu();
  double vs = w.get_vs();
  double tb = v2/v1;
  physical_ESSM(w, mstop, mstopsq, mD1sq, mD2sq, vs, tb);

  // Have we actually converged or is newt just being ridiculous?
  double f1 = ESSM_EWSBCondition1(w);
  double f2 = ESSM_EWSBCondition2(w);
  double f3 = ESSM_EWSBCondition3(w);
  double f4 = estimate(4)*MsusyScaleFactor/Sqrt(mstop(1)*mstop(2))-1.0;

  if (true)//ENABLE_DEBUG)
    {
      cout << "#    Q                     = " << w.get_scale() << endl;
      cout << "#    m_Hd^2                = " << w.get_mHd2() << endl;
      cout << "#    m_Hu^2                = " << w.get_mHu2() << endl;
      cout << "#    m_s^2                 = " << w.get_ms2() << endl;
      cout << "#    f1                    = " << f1 << endl;
      cout << "#    f1/mHdSq Scale Factor = " << f1/mHdSqScaleFactor << endl;
      cout << "#    f2                    = " << f2 << endl;
      cout << "#    f2/mHdSq Scale Factor = " << f1/mHuSqScaleFactor << endl;
      cout << "#    f3                    = " << f3 << endl;
      cout << "#    f3/mSSq Scale Factor  = " << f3/mSSqScaleFactor << endl;
      cout << "#    f4                    = " << f4 << endl; 
    }

  if (fabs(f1) > tol || fabs(f2) > tol || fabs(f3) > tol || fabs(f4) > tol)
    {
      hasProblem = true;
    }

  return hasProblem;

}

//Should be part of softsusy but you may not want to link to the file which
// uses this.  If you do just comment this out.  
double ccbSqrt(double f){ return Sqrt(Abs(f)); }

// Function a0 defined in paper.
// Note sign difference to paper - a0 term
// appears with opposite sign below for some reason.
double a0Peter(double mSq, double Q)
{
  return mSq*(Log(mSq/(Q*Q))-1);
}

// Calculates the masses of the stops and exotic quarks for tadpoles.  the former are always included in our tadpole contributions, the latter could be, but provide small contributions in comparison to errors due to threshold treatment.  
void physical_ESSM(genericE6SSM_soft_parameters r,DoubleVector & mstop, DoubleVector & mstopsq, DoubleVector & mD1sq, DoubleVector & mD2sq, double s, double tb)
{
  bool speak = false;

  // cout << "mq2 = " << r.get_mq2(2,2) << endl;
  // cout << "mu2 = " << r.get_mu2(2,2) << endl;
  // cout << "stop mass matrix: " << endl;

  double yt = r.get_Yu(2, 2);

  double v1 = r.get_vd();
  double v2 = r.get_vu();

  double vev = Sqrt(v1*v1+v2*v2);

  double oneOrt2 = 1/(Sqrt(2.0));
  double mtop = yt*v2/(Sqrt(2.0));

  if(USEMTOFMT){
    mtop = 165;
    yt = mtop/(v2*oneOrt2);
    r.set_Yu(2,2, yt);
    
  }

  if(speak){
    cerr << " yt*v2/(sqrt(2.0)) = " <<  yt*v2/(Sqrt(2.0)) << endl;
    cerr << "mtop = " << mtop << endl;
  }
  
  double oneO40 = 1.0/(40.0);

  DoubleVector mDsq(3), mDbarsq(3), kappa(3), lambda(3);

  lambda(3) = r.get_Lambdax();

  for (int i = 1; i <= 2; i++)
    {
      lambda(i) = r.get_Lambda12(i-1,i-1);
    }

  for (int i = 1; i <= 3; i++)
    {
      kappa(i) = r.get_Kappa(i-1,i-1);
      mDsq(i) = r.get_mDx2(i-1,i-1);
      mDbarsq(i) = r.get_mDxbar2(i-1,i-1);
    }
 
  DoubleVector Akappa(3);
  DoubleVector Alambda(3);

  DoubleVector Tlambda(3), Tkappa(3);

  Tlambda(3) = r.get_TLambdax();

  for (int i = 1; i <= 2; i++)
    {
      Tlambda(i) = r.get_TLambda12(i-1,i-1);

      if (Abs(Tlambda(i)) < EPSTOL) 
	{
	  Alambda(i) = 0.0;
	}
      else if (Abs(lambda(i)) < 1.0e-100)
	{
	  ostringstream ii;
	  ii << "WARNING: trying to calculate A_lambda(" << i << ") where lambda(" << i << ") coupling is " <<
	    Abs(lambda(i)) << endl;
	  throw ii.str();
	}
      else
	{
	  Alambda(i) = Tlambda(i)/lambda(i);
	}

    }

  if (Abs(Tlambda(3)) < EPSTOL) 
    {
      Alambda(3) = 0.0;
    }
  else if (Abs(lambda(3)) < 1.0e-100)
    {
      ostringstream ii;
      ii << "WARNING: trying to calculate A_lambda(3) where lambda(3) coupling is " <<
	Abs(lambda(3)) << endl;
      throw ii.str();
    }
  else
    {
      Alambda(3) = Tlambda(3)/lambda(3);
    }

  double At;
  double TYt = r.get_TYu(2,2);

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

  for (int i = 1; i <= 3; i++)
    {
      Tkappa(i) = r.get_TKappa(i-1,i-1);
      
      if (Abs(Tkappa(i)) < EPSTOL)
	{
	  Akappa(i) = 0.0;
	}
      else if (Abs(kappa(i)) < 1.0e-100)
	{
	  ostringstream ii;
	  ii << "WARNING: trying to calculate A_kappa(" << i << ") where kappa(" << i << ") coupling is " <<
	    Abs(kappa(i)) << endl;
	  throw ii.str();
	}
      else
	{
	  Akappa(i) = Tkappa(i)/kappa(i);
	}
    }

  double mQlsq = r.get_mq2(2,2);
  double mUrsq =  r.get_mu2(2,2);

  double g1 = r.get_g1();
  double g2 = r.get_g2();
  double g1p = r.get_gN();

  //Note that the heavier stop is stop1 this matches Roman's notation.

  mstopsq(1) = 0.5*(mQlsq +  mUrsq + 2*mtop*mtop +0.125*Sqr(g2)*(Sqr(v1) - Sqr(v2)) +0.125*3.0*Sqr(g1)*(Sqr(v1) - Sqr(v2))/5.0 + oneO40*Sqr(g1p)*(-3.0*Sqr(v1) - 2.0*Sqr(v2) + 5.0*Sqr(s))  + Sqrt(Sqr(mQlsq- mUrsq+0.125*Sqr(g2)*(Sqr(v1) - Sqr(v2)) - 0.125*Sqr(g1)*(Sqr(v1) - Sqr(v2))) + 4.0*mtop*mtop*Sqr(At - lambda(3)*s*oneOrt2*v1/v2)));
  //      cout << "mstopsq(1) = "  << mstopsq(1) << endl;
  

  if (mstopsq(1) >= 0.0)
    {
  mstop(1) = Sqrt(0.5*(mQlsq + mUrsq + 2*mtop*mtop + 0.125*Sqr(g2)*(Sqr(v1) - Sqr(v2)) + 0.125*3.0*Sqr(g1)*(Sqr(v1) - Sqr(v2))/5.0 + oneO40*Sqr(g1p)*(-3.0*Sqr(v1) - 2.0*Sqr(v2) + 5.0*Sqr(s))  + Sqrt(Sqr(mQlsq- mUrsq+0.125*Sqr(g2)*(Sqr(v1) - Sqr(v2)) - 0.125*Sqr(g1)*(Sqr(v1) - Sqr(v2))) + 4.0*mtop*mtop*Sqr(At - lambda(3)*s*oneOrt2*v1/v2))));
    }
  else
    {
      mstop(1) = -Sqrt(Abs(mstopsq(1)));
    }

  mstopsq(2) = 0.5*(mQlsq + mUrsq + 2*mtop*mtop +0.125*Sqr(g2)*(Sqr(v1) - Sqr(v2)) +0.125*3.0*Sqr(g1)*(Sqr(v1) - Sqr(v2))/5.0 + oneO40*Sqr(g1p)*(-3.0*Sqr(v1) - 2.0*Sqr(v2) + 5.0*Sqr(s))  - Sqrt(Sqr(mQlsq- mUrsq+0.125*Sqr(g2)*(Sqr(v1) - Sqr(v2)) - 0.125*Sqr(g1)*(Sqr(v1) - Sqr(v2))) + 4.0*mtop*mtop*Sqr(At - lambda(3)*s*oneOrt2*v1/v2)));
  //cout << "mstopsq(2) = "  << mstopsq(2) << endl;
  

  if (mstopsq(2) >= 0.0)
    { 
  mstop(2) = Sqrt(0.5*(mQlsq + mUrsq + 2*mtop*mtop +0.125*Sqr(g2)*(Sqr(v1) - Sqr(v2)) +0.125*3.0*Sqr(g1)*(Sqr(v1) - Sqr(v2))/5.0 + oneO40*Sqr(g1p)*(-3.0*Sqr(v1) - 2.0*Sqr(v2) + 5.0*Sqr(s))  - Sqrt(Sqr(mQlsq- mUrsq+0.125*Sqr(g2)*(Sqr(v1) - Sqr(v2)) - 0.125*Sqr(g1)*(Sqr(v1) - Sqr(v2))) + 4.0*mtop*mtop*Sqr(At - lambda(3)*s*oneOrt2*v1/v2))));
    }
  else
    {
      mstop(2) = -Sqrt(Abs(mstopsq(2)));
    }
 
  if(speak){								   
    cerr << " mstop(1) = " << mstop(1) << endl;
    cerr << " mstop(2) = " << mstop(2) << endl;
  }
 
  // DH:: warn if stops are found to be tachyonic
  if (mstopsq(1) < 0.0 || mstopsq(2) < 0.0)
    {
      cerr << "Warning: tachyonic stop masses." << endl;
      cerr << "m_stop_1^2 = " << mstopsq(1) << " GeV^2." << endl;
      cerr << "m_stop_2^2 = " << mstopsq(2) << " GeV^2." << endl;
    }  

  int gen=1;
  for(gen=1; gen< 4; gen++){
    mD1sq(gen) = 0.5*(mDsq(gen) + mDbarsq(gen) + kappa(gen)*kappa(gen)*s*s - 2.5*oneO40*Sqr(g1p)*(-3.0*Sqr(v1) - 2.0*Sqr(v2) + 5.0*Sqr(s))  
		      - Sqrt(Sqr(mDsq(gen)- mDbarsq(gen) + 0.5*oneO40*Sqr(g1p)*(-3.0*Sqr(v1) - 2.0*Sqr(v2) + 5.0*Sqr(s))   + 0.1*Sqr(g1)*(Sqr(v1) - Sqr(v2))) + 4.0*Sqr(Akappa(gen)*kappa(gen)*s*oneOrt2 - 0.5*kappa(gen)*lambda(3)*v1*v2)));

    mD2sq(gen) = 0.5*(mDsq(gen) + mDbarsq(gen) + kappa(gen)*kappa(gen)*s*s - 2.5*oneO40*Sqr(g1p)*(-3.0*Sqr(v1) - 2.0*Sqr(v2) + 5.0*Sqr(s))  
		      + Sqrt(Sqr(mDsq(gen)- mDbarsq(gen) + 0.5*oneO40*Sqr(g1p)*(-3.0*Sqr(v1) - 2.0*Sqr(v2) + 5.0*Sqr(s))   + 0.1*Sqr(g1)*(Sqr(v1) - Sqr(v2))) + 4.0*Sqr(Akappa(gen)*kappa(gen)*s*oneOrt2 - 0.5*kappa(gen)*lambda(3)*v1*v2)));
}
return;
}

//This version neglects U(1) D-terms and was written for comparison with Romans program
void physical_ESSM_Roman(genericE6SSM_soft_parameters r,DoubleVector & mstop, DoubleVector & mstopsq, DoubleVector & mD1sq, DoubleVector & mD2sq, double s, double tb){

  bool speak = false;
  double yt = r.get_Yu(2, 2);

  double v1 = r.get_vd();
  double v2 = r.get_vu();

  double vev = Sqrt(v1*v1+v2*v2);

  double oneOrt2 = 1/(Sqrt(2.0));
  double mtop = yt*v2/(Sqrt(2.0));
  mtop = 165;
  yt = mtop/(v2*oneOrt2);
  r.set_Yu(2,2, yt);

  if(speak){
    cout << "mtop = " << mtop << endl;
    cout << "yt = " << yt << endl; 
  }
  
  double oneO40 = 1.0/(40.0);

  DoubleVector lambda(3), kappa(3), Akappa(3), Alambda(3), mDsq(3), mDbarsq(3);
  
  lambda(3) = r.get_Lambdax();
  
  for (int i = 1; i <= 2; i++)
    {
      lambda(i) = r.get_Lambda12(i-1,i-1);
    }
  
  for (int i = 1; i <= 3; i++)
    {
      kappa(i) = r.get_Kappa(i-1,i-1);
      mDsq(i) = r.get_mDx2(i-1,i-1);
      mDbarsq(i) = r.get_mDxbar2(i-1,i-1);
    }
  
  
  DoubleVector Tlambda(3), Tkappa(3);
  
  Tlambda(3) = r.get_TLambdax();
  
  for (int i = 1; i <= 2; i++)
    {
      Tlambda(i) = r.get_TLambda12(i-1,i-1);
      
      if (Abs(Tlambda(i)) < EPSTOL) 
	{
	  Alambda(i) = 0.0;
	}
      else if (Abs(lambda(i)) < 1.0e-100)
	{
	  ostringstream ii;
	  ii << "WARNING: trying to calculate A_lambda(" << i << ") where lambda(" << i << ") coupling is " <<
	    Abs(lambda(i)) << endl;
	  throw ii.str();
	}
      else
	{
	  Alambda(i) = Tlambda(i)/lambda(i);
	}

    }

  if (Abs(Tlambda(3)) < EPSTOL) 
    {
      Alambda(3) = 0.0;
    }
  else if (Abs(lambda(3)) < 1.0e-100)
    {
      ostringstream ii;
      ii << "WARNING: trying to calculate A_lambda(3) where lambda(3) coupling is " <<
	Abs(lambda(3)) << endl;
      throw ii.str();
    }
  else
    {
      Alambda(3) = Tlambda(3)/lambda(3);
    }

  double At;
  double TYt = r.get_TYu(2,2);

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

  for (int i = 1; i <= 3; i++)
    {
      Tkappa(i) = r.get_TKappa(i-1,i-1);
      
      if (Abs(Tkappa(i)) < EPSTOL)
	{
	  Akappa(i) = 0.0;
	}
      else if (Abs(kappa(i)) < 1.0e-100)
	{
	  ostringstream ii;
	  ii << "WARNING: trying to calculate A_kappa(" << i << ") where kappa(" << i << ") coupling is " <<
	    Abs(kappa(i)) << endl;
	  throw ii.str();
	}
      else
	{
	  Akappa(i) = Tkappa(i)/kappa(i);
	}
    }


  double mQlsq = r.get_mq2(2,2);
  double mUrsq =  r.get_mu2(2,2);

  double g1 = r.get_g1();
  double g2 = r.get_g2();
  double g1p = r.get_gN();
  //getting stop masses 
  //Note that the heavier stop is stop1, this matches Roman's notation.
  mstop(1) = Sqrt(0.5*(mQlsq +  mUrsq + 2*mtop*mtop  + Sqrt(Sqr(mQlsq- mUrsq) + 4.0*mtop*mtop*Sqr(At - lambda(3)*s*oneOrt2*v1/v2))));
  
  mstopsq(1) = 0.5*(mQlsq +  mUrsq + 2*mtop*mtop + Sqrt(Sqr(mQlsq- mUrsq) + 4.0*mtop*mtop*Sqr(At - lambda(3)*s*oneOrt2*v1/v2)));
  
  mstop(2) = Sqrt(0.5*(mQlsq +  mUrsq + 2*mtop*mtop   - Sqrt(Sqr(mQlsq- mUrsq) + 4.0*mtop*mtop*Sqr(At - lambda(3)*s*oneOrt2*v1/v2))));
  
  mstopsq(2) = 0.5*(mQlsq +  mUrsq + 2*mtop*mtop  - Sqrt(Sqr(mQlsq- mUrsq) + 4.0*mtop*mtop*Sqr(At - lambda(3)*s*oneOrt2*v1/v2)));									
  
  // DH:: warn if stops are found to be tachyonic
   if (mstopsq(1) < 0.0 || mstopsq(2) < 0.0)
    {
      cerr << "Warning: tachyonic stop masses." << endl;
      cerr << "m_stop_1^2 = " << mstopsq(1) << " GeV^2." << endl;
      cerr << "m_stop_2^2 = " << mstopsq(2) << " GeV^2." << endl;
    }

  int gen=1;
  for(gen=1; gen< 4; gen++){
    mD1sq(gen) = 0.5*(mDsq(gen) + mDbarsq(gen) + kappa(gen)*kappa(gen)*s*s - 2.5*oneO40*Sqr(g1p)*(-3.0*Sqr(v1) - 2.0*Sqr(v2) + 5.0*Sqr(s))  
		      - Sqrt(Sqr(mDsq(gen)- mDbarsq(gen) + 0.5*oneO40*Sqr(g1p)*(-3.0*Sqr(v1) - 2.0*Sqr(v2) + 5.0*Sqr(s))   + 0.1*Sqr(g1)*(Sqr(v1) - Sqr(v2))) + 4.0*Sqr(Akappa(gen)*kappa(gen)*s*oneOrt2 - 0.5*kappa(gen)*lambda(3)*v1*v2)));
    
    mD2sq(gen) = 0.5*(mDsq(gen) + mDbarsq(gen) + kappa(gen)*kappa(gen)*s*s - 2.5*oneO40*Sqr(g1p)*(-3.0*Sqr(v1) - 2.0*Sqr(v2) + 5.0*Sqr(s))  
		      + Sqrt(Sqr(mDsq(gen)- mDbarsq(gen) + 0.5*oneO40*Sqr(g1p)*(-3.0*Sqr(v1) - 2.0*Sqr(v2) + 5.0*Sqr(s))   + 0.1*Sqr(g1)*(Sqr(v1) - Sqr(v2))) + 4.0*Sqr(Akappa(gen)*kappa(gen)*s*oneOrt2 - 0.5*kappa(gen)*lambda(3)*v1*v2)));  
  }
  return;
}


//Determines the dominant H1 tadpoles. There are a number of different versions of this depending on choice of scale. The difference  should be higher order. Exotics can also be included in this version though by they switched off at the moment.   This version  includes U(1) D-terms and exotics, and evaluates the scale where RG evelotion is halted.
double doCalcTadpoleESSMH1(genericE6SSM_soft_parameters r,  double s , double tb)
{
  
  bool speak = false;

  double threeO32pisq = 3.0/(32.0*Sqr(PI));
  
  double oneO40 = 1.0/(40.0);
  double yt = r.get_Yu(2, 2);

  double v1 = r.get_vd();
  double v2 = r.get_vu();
  double vev = Sqrt(v1*v1+v2*v2);

  double oneOrt2 = 1/(Sqrt(2.0));

  double mtop = yt*v2*oneOrt2;

  double q = r.get_scale();

  if(USEMTOFMT){
    //    cout << "using mtop = " << mtop << endl; 
    mtop = 165;
    yt = mtop/(v2*oneOrt2);
    r.set_Yu(2,2, yt);
  }

  DoubleVector lambda(3), kappa(3), Akappa(3), Alambda(3), mDsq(3), mDbarsq(3);

  lambda(3) = r.get_Lambdax();

  for (int i = 1; i <= 2; i++)
    {
      lambda(i) = r.get_Lambda12(i-1,i-1);
    }

  for (int i = 1; i <= 3; i++)
    {
      kappa(i) = r.get_Kappa(i-1,i-1);
      mDsq(i) = r.get_mDx2(i-1,i-1);
      mDbarsq(i) = r.get_mDxbar2(i-1,i-1);
    }


  DoubleVector Tlambda(3), Tkappa(3);

  Tlambda(3) = r.get_TLambdax();

  for (int i = 1; i <= 2; i++)
    {
      Tlambda(i) = r.get_TLambda12(i-1,i-1);

      if (Abs(Tlambda(i)) < EPSTOL) 
	{
	  Alambda(i) = 0.0;
	}
      else if (Abs(lambda(i)) < 1.0e-100)
	{
	  ostringstream ii;
	  ii << "WARNING: trying to calculate A_lambda(" << i << ") where lambda(" << i << ") coupling is " <<
	    Abs(lambda(i)) << endl;
	  throw ii.str();
	}
      else
	{
	  Alambda(i) = Tlambda(i)/lambda(i);
	}

    }

  if (Abs(Tlambda(3)) < EPSTOL) 
    {
      Alambda(3) = 0.0;
    }
  else if (Abs(lambda(3)) < 1.0e-100)
    {
      ostringstream ii;
      ii << "WARNING: trying to calculate A_lambda(3) where lambda(3) coupling is " <<
	Abs(lambda(3)) << endl;
      throw ii.str();
    }
  else
    {
      Alambda(3) = Tlambda(3)/lambda(3);
    }

  double At;
  double TYt = r.get_TYu(2,2);

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

  for (int i = 1; i <= 3; i++)
    {
      Tkappa(i) = r.get_TKappa(i-1,i-1);
      
      if (Abs(Tkappa(i)) < EPSTOL)
	{
	  Akappa(i) = 0.0;
	}
      else if (Abs(kappa(i)) < 1.0e-100)
	{
	  ostringstream ii;
	  ii << "WARNING: trying to calculate A_kappa(" << i << ") where kappa(" << i << ") coupling is " <<
	    Abs(kappa(i)) << endl;
	  throw ii.str();
	}
      else
	{
	  Akappa(i) = Tkappa(i)/kappa(i);
	}
    }

  
  DoubleVector mstop(2), mstopsq(2);//stop1,stop2
  DoubleVector mD1sq(3), mD2sq(3);//3 gens

  double mQlsq = r.get_mq2(2,2);
  double mUrsq =  r.get_mu2(2,2);

  double g1 = r.get_g1();
  double g2 = r.get_g2();
  double g1p = r.get_gN();

  physical_ESSM(r,mstop,mstopsq, mD1sq,mD2sq, s, tb);
   
  if (mstopsq(1) < 0.0 || mstopsq(2) < 0.0)
    {
      mstopsq(1) = Abs(mstopsq(1));
      mstopsq(2) = Abs(mstopsq(2));
    }

  double  Delmtop = 0;
  double  Delmstop1 = 0.5*(0.6*0.25*Sqr(g1) + 0.25*Sqr(g2) - 6.0/(40.0) *Sqr(g1p) + 0.5*( 2.0*(mQlsq- mUrsq+0.125*Sqr(g2)*(Sqr(v1) - Sqr(v2)) - 0.125*Sqr(g1)*(Sqr(v1) - Sqr(v2)))*0.25*(Sqr(g2)- Sqr(g1)) + 8.0*Sqr(mtop)*(At  - lambda(3)*s*oneOrt2*v1/v2)*(-lambda(3)*s*oneOrt2/(v2*v1) ) ) /(   Sqrt(Sqr(mQlsq- mUrsq+0.125*Sqr(g2)*(Sqr(v1) - Sqr(v2)) - 0.125*Sqr(g1)*(Sqr(v1) - Sqr(v2))) + 4.0*mtop*mtop*Sqr(At - lambda(3)*s*oneOrt2*v1/v2))  )     );

  double Delmstop2 = 0.5*(0.6*0.25*Sqr(g1) + 0.25*Sqr(g2) - 6.0/(40.0) *Sqr(g1p)
			  -0.5*( 2.0*(mQlsq- mUrsq+0.125*Sqr(g2)*(Sqr(v1) - Sqr(v2)) - 0.125*Sqr(g1)*(Sqr(v1) - Sqr(v2)))*0.25*(Sqr(g2)- Sqr(g1)) + 8.0*Sqr(mtop)*(At  - lambda(3)*s*oneOrt2*v1/v2)*(-lambda(3)*s*oneOrt2/(v2*v1) )              ) /(   Sqrt(Sqr(mQlsq- mUrsq+0.125*Sqr(g2)*(Sqr(v1) - Sqr(v2)) - 0.125*Sqr(g1)*(Sqr(v1) - Sqr(v2))) + 4.0*mtop*mtop*Sqr(At - lambda(3)*s*oneOrt2*v1/v2))  )     );


  //tested and debugged up tothis point 27/06/07
  double delta =   threeO32pisq*(-2.0*a0Peter(mstopsq(1),q)*Delmstop1-2.0*a0Peter(mstopsq(2),q)*Delmstop2 + 4.0*a0Peter(Sqr(mtop),q)*Delmtop);		
		     
  return delta; 
}

//Calculates tadpoles using logarithms at m_t.  The assumption is coefficients do not substantially change between scale at which RG evolution is halted and m_t. 
double doCalcTadpoleESSMH1_atMt(genericE6SSM_soft_parameters r,  double s , double tb){
  bool speak = false;

  double threeO32pisq = 3.0/(32.0*Sqr(PI));
  //cout <<"threeO32pisq = "<< threeO32pisq << endl;
  double oneO40 = 1.0/(40.0);

  double yt = r.get_Yu(2, 2);

  double v1 = r.get_vd();
  double v2 = r.get_vu();

  double vev = Sqrt(v1*v1+v2*v2);

  double oneOrt2 = 1/(Sqrt(2.0));
  double mtop = yt*v2*oneOrt2;
  double q = r.get_scale();
  
  // mtop = 165;
  // q=mtop;
  q = 165;
  // cout << "q = " << q << endl;
  // cout << "mtop - " << mtop << endl;
 

  if(USEMTOFMT){
    mtop = 165;
    yt = mtop/(v2*oneOrt2);
    r.set_Yu(2,2, yt);    
  }

  DoubleVector lambda(3), kappa(3), Akappa(3), Alambda(3), mDsq(3), mDbarsq(3);

  lambda(3) = r.get_Lambdax();

  for (int i = 1; i <= 2; i++)
    {
      lambda(i) = r.get_Lambda12(i-1,i-1);
    }

  for (int i = 1; i <= 3; i++)
    {
      kappa(i) = r.get_Kappa(i-1,i-1);
      mDsq(i) = r.get_mDx2(i-1,i-1);
      mDbarsq(i) = r.get_mDxbar2(i-1,i-1);
    }


  DoubleVector Tlambda(3), Tkappa(3);

  Tlambda(3) = r.get_TLambdax();

  for (int i = 1; i <= 2; i++)
    {
      Tlambda(i) = r.get_TLambda12(i-1,i-1);

      if (Abs(Tlambda(i)) < EPSTOL) 
	{
	  Alambda(i) = 0.0;
	}
      else if (Abs(lambda(i)) < 1.0e-100)
	{
	  ostringstream ii;
	  ii << "WARNING: trying to calculate A_lambda(" << i << ") where lambda(" << i << ") coupling is " <<
	    Abs(lambda(i)) << endl;
	  throw ii.str();
	}
      else
	{
	  Alambda(i) = Tlambda(i)/lambda(i);
	}

    }

  if (Abs(Tlambda(3)) < EPSTOL) 
    {
      Alambda(3) = 0.0;
    }
  else if (Abs(lambda(3)) < 1.0e-100)
    {
      ostringstream ii;
      ii << "WARNING: trying to calculate A_lambda(3) where lambda(3) coupling is " <<
	Abs(lambda(3)) << endl;
      throw ii.str();
    }
  else
    {
      Alambda(3) = Tlambda(3)/lambda(3);
    }

  double At;
  double TYt = r.get_TYu(2,2);

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

  for (int i = 1; i <= 3; i++)
    {
      Tkappa(i) = r.get_TKappa(i-1,i-1);
      
      if (Abs(Tkappa(i)) < EPSTOL)
	{
	  Akappa(i) = 0.0;
	}
      else if (Abs(kappa(i)) < 1.0e-100)
	{
	  ostringstream ii;
	  ii << "WARNING: trying to calculate A_kappa(" << i << ") where kappa(" << i << ") coupling is " <<
	    Abs(kappa(i)) << endl;
	  throw ii.str();
	}
      else
	{
	  Akappa(i) = Tkappa(i)/kappa(i);
	}
    }

  DoubleVector mstop(2), mstopsq(2);//stop1,stop2
  DoubleVector mD1sq(3), mD2sq(3);//3 gens

  double mQlsq = r.get_mq2(2,2);
  double mUrsq =  r.get_mu2(2,2);

  double g1 = r.get_g1();
  double g2 = r.get_g2();
  double g1p = r.get_gN();
  physical_ESSM(r,mstop,mstopsq, mD1sq,mD2sq, s, tb);
  if (mstopsq(1) < 0.0 || mstopsq(2) < 0.0)
    {
      mstopsq(1) = Abs(mstopsq(1));
      mstopsq(2) = Abs(mstopsq(2));
    }
  double  Delmtop = 0;
  double  Delmstop1 = 0.5*(0.6*0.25*Sqr(g1) + 0.25*Sqr(g2) - 6.0/(40.0) *Sqr(g1p)
			  +0.5*( 2.0*(mQlsq- mUrsq+0.125*Sqr(g2)*(Sqr(v1) - Sqr(v2)) - 0.125*Sqr(g1)*(Sqr(v1) - Sqr(v2)))*0.25*(Sqr(g2)- Sqr(g1)) + 8.0*Sqr(mtop)*(At  - lambda(3)*s*oneOrt2*v1/v2)*(-lambda(3)*s*oneOrt2/(v2*v1) )              ) /(   Sqrt(Sqr(mQlsq- mUrsq+0.125*Sqr(g2)*(Sqr(v1) - Sqr(v2)) - 0.125*Sqr(g1)*(Sqr(v1) - Sqr(v2))) + 4.0*mtop*mtop*Sqr(At - lambda(3)*s*oneOrt2*v1/v2))  )     );
 
 double Delmstop2 = 0.5*(0.6*0.25*Sqr(g1) + 0.25*Sqr(g2) - 6.0/(40.0) *Sqr(g1p)
			 -0.5*( 2.0*(mQlsq- mUrsq+0.125*Sqr(g2)*(Sqr(v1) - Sqr(v2)) - 0.125*Sqr(g1)*(Sqr(v1) - Sqr(v2)))*0.25*(Sqr(g2)- Sqr(g1)) + 8.0*Sqr(mtop)*(At  - lambda(3)*s*oneOrt2*v1/v2)*(-lambda(3)*s*oneOrt2/(v2*v1) )              ) /(   Sqrt(Sqr(mQlsq- mUrsq+0.125*Sqr(g2)*(Sqr(v1) - Sqr(v2)) - 0.125*Sqr(g1)*(Sqr(v1) - Sqr(v2))) + 4.0*mtop*mtop*Sqr(At - lambda(3)*s*oneOrt2*v1/v2))  )     );
  
//tested and debugged up tothis point 27/06/07

 double delta =   threeO32pisq*(-2.0*a0Peter(mstopsq(1),q)*Delmstop1-2.0*a0Peter(mstopsq(2),q)*Delmstop2 + 4.0*a0Peter(Sqr(mtop),q)*Delmtop);

 return delta;
}

//This version neglects U(1) D-terms and calculates at mt. 
double doCalcTadpoleESSMH1_Roman(genericE6SSM_soft_parameters r,  double s , double tb){
  bool speak = false;
  double threeO32pisq = 3.0/(32.0*Sqr(PI));
  //cout <<"threeO32pisq = "<< threeO32pisq << endl;
  double oneO40 = 1.0/(40.0);
  double yt = r.get_Yu(2, 2);

  double v1 = r.get_vd();
  double v2 = r.get_vu();

  double vev = Sqrt(v1*v1+v2*v2);

  double oneOrt2 = 1/(Sqrt(2.0));
  double mtop = yt*v2*oneOrt2;
  mtop = 165; 
  yt = mtop/(v2*oneOrt2);
  r.set_Yu(2,2, yt);

  double q = r.get_scale();
  q = 165; //GeV. Fudgeing to match Romans code

  DoubleVector lambda(3), kappa(3), Akappa(3), Alambda(3), mDsq(3), mDbarsq(3);

  lambda(3) = r.get_Lambdax();

  for (int i = 1; i <= 2; i++)
    {
      lambda(i) = r.get_Lambda12(i-1,i-1);
    }

  for (int i = 1; i <= 3; i++)
    {
      kappa(i) = r.get_Kappa(i-1,i-1);
      mDsq(i) = r.get_mDx2(i-1,i-1);
      mDbarsq(i) = r.get_mDxbar2(i-1,i-1);
    }


  DoubleVector Tlambda(3), Tkappa(3);

  Tlambda(3) = r.get_TLambdax();

  for (int i = 1; i <= 2; i++)
    {
      Tlambda(i) = r.get_TLambda12(i-1,i-1);

      if (Abs(Tlambda(i)) < EPSTOL) 
	{
	  Alambda(i) = 0.0;
	}
      else if (Abs(lambda(i)) < 1.0e-100)
	{
	  ostringstream ii;
	  ii << "WARNING: trying to calculate A_lambda(" << i << ") where lambda(" << i << ") coupling is " <<
	    Abs(lambda(i)) << endl;
	  throw ii.str();
	}
      else
	{
	  Alambda(i) = Tlambda(i)/lambda(i);
	}

    }

  if (Abs(Tlambda(3)) < EPSTOL) 
    {
      Alambda(3) = 0.0;
    }
  else if (Abs(lambda(3)) < 1.0e-100)
    {
      ostringstream ii;
      ii << "WARNING: trying to calculate A_lambda(3) where lambda(3) coupling is " <<
	Abs(lambda(3)) << endl;
      throw ii.str();
    }
  else
    {
      Alambda(3) = Tlambda(3)/lambda(3);
    }

  double At;
  double TYt = r.get_TYu(2,2);

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

  for (int i = 1; i <= 3; i++)
    {
      Tkappa(i) = r.get_TKappa(i-1,i-1);
      
      if (Abs(Tkappa(i)) < EPSTOL)
	{
	  Akappa(i) = 0.0;
	}
      else if (Abs(kappa(i)) < 1.0e-100)
	{
	  ostringstream ii;
	  ii << "WARNING: trying to calculate A_kappa(" << i << ") where kappa(" << i << ") coupling is " <<
	    Abs(kappa(i)) << endl;
	  throw ii.str();
	}
      else
	{
	  Akappa(i) = Tkappa(i)/kappa(i);
	}
    }

  DoubleVector mstop(2), mstopsq(2);//stop1,stop2
  DoubleVector mD1sq(3), mD2sq(3);//3 gens
  double mQlsq = r.get_mq2(2,2);
  double mUrsq =  r.get_mu2(2,2);
  
  double g1 = r.get_g1();
  double g2 = r.get_g2();
  double g1p = r.get_gN();

  physical_ESSM_Roman(r,mstop,mstopsq, mD1sq,mD2sq, s, tb);
  if (mstopsq(1) < 0.0 || mstopsq(2) < 0.0)
    {
      mstopsq(1) = Abs(mstopsq(1));
      mstopsq(2) = Abs(mstopsq(2));
    }
  double  Delmtop = 0;
  double  Delmstop1 =  0.5*( 0.5* 8.0*Sqr(mtop)*(At  - lambda(3)*s*oneOrt2*v1/v2)*(-lambda(3)*s*oneOrt2/(v2*v1) ) ) /( Sqrt(Sqr(mQlsq- mUrsq) + 4.0*mtop*mtop*Sqr(At - lambda(3)*s*oneOrt2*v1/v2))  ) ;
  
  double Delmstop2 = -0.5*( 0.5*(  8.0*Sqr(mtop)*(At  - lambda(3)*s*oneOrt2*v1/v2)*(-lambda(3)*s*oneOrt2/(v2*v1) )              ) /(   Sqrt(Sqr(mQlsq- mUrsq) + 4.0*mtop*mtop*Sqr(At - lambda(3)*s*oneOrt2*v1/v2))  )     );

  double delta =   threeO32pisq*(-2.0*a0Peter(mstopsq(1),q)*Delmstop1-2.0*a0Peter(mstopsq(2),q)*Delmstop2 + 4.0*a0Peter(Sqr(mtop),q)*Delmtop);

  return delta;   
}
 
//Yet another option, performs Calculation like Roman, but uses logarthyms at scale at which RGE evolution is halted.
double doCalcTadpoleESSMH1_Roman_atQ(genericE6SSM_soft_parameters r,  double s , double tb){
  double threeO32pisq = 3.0/(32.0*Sqr(PI));
  //cout <<"threeO32pisq = "<< threeO32pisq << endl;
  double oneO40 = 1.0/(40.0);
  double yt = r.get_Yu(2, 2);

  double v1 = r.get_vd();
  double v2 = r.get_vu();

  double vev = Sqrt(v1*v1+v2*v2);

  double oneOrt2 = 1/(Sqrt(2.0));
  double mtop = yt*v2*oneOrt2;
  //mtop = 165; 
  yt = mtop/(v2*oneOrt2);
  r.set_Yu(2,2, yt);
  
  if(USEMTOFMT){
    mtop = 165;
    yt = mtop/(v2*oneOrt2);
    r.set_Yu(2,2, yt);
  }
 
  double q = r.get_scale();
  //q = 165; //GeV. Fudgeing to match Romans code

  DoubleVector lambda(3), kappa(3), Akappa(3), Alambda(3), mDsq(3), mDbarsq(3);

  lambda(3) = r.get_Lambdax();

  for (int i = 1; i <= 2; i++)
    {
      lambda(i) = r.get_Lambda12(i-1,i-1);
    }

  for (int i = 1; i <= 3; i++)
    {
      kappa(i) = r.get_Kappa(i-1,i-1);
      mDsq(i) = r.get_mDx2(i-1,i-1);
      mDbarsq(i) = r.get_mDxbar2(i-1,i-1);
    }


  DoubleVector Tlambda(3), Tkappa(3);

  Tlambda(3) = r.get_TLambdax();

  for (int i = 1; i <= 2; i++)
    {
      Tlambda(i) = r.get_TLambda12(i-1,i-1);

      if (Abs(Tlambda(i)) < EPSTOL) 
	{
	  Alambda(i) = 0.0;
	}
      else if (Abs(lambda(i)) < 1.0e-100)
	{
	  ostringstream ii;
	  ii << "WARNING: trying to calculate A_lambda(" << i << ") where lambda(" << i << ") coupling is " <<
	    Abs(lambda(i)) << endl;
	  throw ii.str();
	}
      else
	{
	  Alambda(i) = Tlambda(i)/lambda(i);
	}

    }

  if (Abs(Tlambda(3)) < EPSTOL) 
    {
      Alambda(3) = 0.0;
    }
  else if (Abs(lambda(3)) < 1.0e-100)
    {
      ostringstream ii;
      ii << "WARNING: trying to calculate A_lambda(3) where lambda(3) coupling is " <<
	Abs(lambda(3)) << endl;
      throw ii.str();
    }
  else
    {
      Alambda(3) = Tlambda(3)/lambda(3);
    }

  double At;
  double TYt = r.get_TYu(2,2);

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

  for (int i = 1; i <= 3; i++)
    {
      Tkappa(i) = r.get_TKappa(i-1,i-1);
      
      if (Abs(Tkappa(i)) < EPSTOL)
	{
	  Akappa(i) = 0.0;
	}
      else if (Abs(kappa(i)) < 1.0e-100)
	{
	  ostringstream ii;
	  ii << "WARNING: trying to calculate A_kappa(" << i << ") where kappa(" << i << ") coupling is " <<
	    Abs(kappa(i)) << endl;
	  throw ii.str();
	}
      else
	{
	  Akappa(i) = Tkappa(i)/kappa(i);
	}
    }

  DoubleVector mstop(2), mstopsq(2);//stop1,stop2
  DoubleVector mD1sq(3), mD2sq(3);//3 gens
  double mQlsq = r.get_mq2(2,2);
  double mUrsq =  r.get_mu2(2,2);
  
  double g1 = r.get_g1();
  double g2 = r.get_g2();
  double g1p = r.get_gN();

  physical_ESSM_Roman(r,mstop,mstopsq, mD1sq,mD2sq, s, tb);
  if (mstopsq(1) < 0.0 || mstopsq(2) < 0.0)
    {
      mstopsq(1) = Abs(mstopsq(1));
      mstopsq(2) = Abs(mstopsq(2));
    }
  double  Delmtop = 0;
  double  Delmstop1 =  0.5*( 0.5* 8.0*Sqr(mtop)*(At  - lambda(3)*s*oneOrt2*v1/v2)*(-lambda(3)*s*oneOrt2/(v2*v1) ) ) /( Sqrt(Sqr(mQlsq- mUrsq) + 4.0*mtop*mtop*Sqr(At - lambda(3)*s*oneOrt2*v1/v2))  ) ;
  
  double Delmstop2 = -0.5*( 0.5*(  8.0*Sqr(mtop)*(At  - lambda(3)*s*oneOrt2*v1/v2)*(-lambda(3)*s*oneOrt2/(v2*v1) )              ) /(   Sqrt(Sqr(mQlsq- mUrsq) + 4.0*mtop*mtop*Sqr(At - lambda(3)*s*oneOrt2*v1/v2))  )     );

//tested and debugged up tothis point 27/06/07

 double delta =   threeO32pisq*(-2.0*a0Peter(mstopsq(1),q)*Delmstop1-2.0*a0Peter(mstopsq(2),q)*Delmstop2 + 4.0*a0Peter(Sqr(mtop),q)*Delmtop);
 					    
 return delta;  
}


//Now tadpoles for H2. The different versions are labled as for H1 tadpoles and contributions are matched with those
double doCalcTadpoleESSMH2(genericE6SSM_soft_parameters r, double s , double tb  ){
   bool speak =false;
   double yt = r.get_Yu(2, 2);
   // cout << "yt = " << yt << endl;

   double v1 = r.get_vd();
   double v2 = r.get_vu();
 
   double vev = Sqrt(v1*v1+v2*v2);

   double oneOrt2 = 1/(Sqrt(2.0));
   double mtop = yt*v2/(Sqrt(2.0));
   double q = r.get_scale();
   
   if(USEMTOFMT){
     mtop = 165;
     yt = mtop/(v2*oneOrt2);
     r.set_Yu(2,2, yt);
   }
 
   double threeO32pisq = 3.0/(32.0*Sqr(PI));
   //double oneOrt2 = 1/(Sqrt(2));			      
   double top = - 6.0 * Sqr(yt) * a0Peter(mtop, q)/ (16.0 * Sqr(PI));
   double topPeter = - 6.0 * Sqr(yt) * a0Peter(Sqr(mtop), q)/ (16.0 * Sqr(PI));
   double oneO40 = 1.0/(40.0);

   DoubleVector lambda(3), kappa(3), Akappa(3), Alambda(3), mDsq(3), mDbarsq(3);
   
   lambda(3) = r.get_Lambdax();
   
   for (int i = 1; i <= 2; i++)
     {
       lambda(i) = r.get_Lambda12(i-1,i-1);
     }
   
   for (int i = 1; i <= 3; i++)
     {
       kappa(i) = r.get_Kappa(i-1,i-1);
       mDsq(i) = r.get_mDx2(i-1,i-1);
       mDbarsq(i) = r.get_mDxbar2(i-1,i-1);
     }
   
   
   DoubleVector Tlambda(3), Tkappa(3);
   
   Tlambda(3) = r.get_TLambdax();
   
   for (int i = 1; i <= 2; i++)
     {
       Tlambda(i) = r.get_TLambda12(i-1,i-1);
       
       if (Abs(Tlambda(i)) < EPSTOL) 
	 {
	   Alambda(i) = 0.0;
	 }
       else if (Abs(lambda(i)) < 1.0e-100)
	 {
	   ostringstream ii;
	   ii << "WARNING: trying to calculate A_lambda(" << i << ") where lambda(" << i << ") coupling is " <<
	     Abs(lambda(i)) << endl;
	   throw ii.str();
	 }
       else
	 {
	   Alambda(i) = Tlambda(i)/lambda(i);
	 }
       
     }
   
   if (Abs(Tlambda(3)) < EPSTOL) 
     {
       Alambda(3) = 0.0;
     }
   else if (Abs(lambda(3)) < 1.0e-100)
     {
       ostringstream ii;
       ii << "WARNING: trying to calculate A_lambda(3) where lambda(3) coupling is " <<
	 Abs(lambda(3)) << endl;
       throw ii.str();
     }
   else
     {
       Alambda(3) = Tlambda(3)/lambda(3);
     }
   
   double At;
   double TYt = r.get_TYu(2,2);
   
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
   
   for (int i = 1; i <= 3; i++)
     {
       Tkappa(i) = r.get_TKappa(i-1,i-1);
       
       if (Abs(Tkappa(i)) < EPSTOL)
	 {
	   Akappa(i) = 0.0;
	 }
       else if (Abs(kappa(i)) < 1.0e-100)
	 {
	   ostringstream ii;
	   ii << "WARNING: trying to calculate A_kappa(" << i << ") where kappa(" << i << ") coupling is " <<
	     Abs(kappa(i)) << endl;
	   throw ii.str();
	 }
       else
	 {
	   Akappa(i) = Tkappa(i)/kappa(i);
	 }
     }
   
   double mQlsq = r.get_mq2(2,2);
   double mUrsq =  r.get_mu2(2,2);
   
   double g1 = r.get_g1();
   double g2 = r.get_g2();
   double g1p = r.get_gN();
   double Delmtop, Delmstop1,Delmstop2; 

   
   DoubleVector mstop(2), mstopsq(2);//stop1,stop2
   DoubleVector mD1sq(3), mD2sq(3);//3 gens
   
   physical_ESSM(r,mstop, mstopsq, mD1sq,mD2sq, s, tb); 

   if (mstopsq(1) < 0.0 || mstopsq(2) < 0.0)
     {
       mstopsq(1) = Abs(mstopsq(1));
       mstopsq(2) = Abs(mstopsq(2));
     }

   Delmtop = Sqr(r.get_Yu(2,2));
 
   Delmstop1 = 0.5*(2.0*Sqr(r.get_Yu(2,2))- 0.6*0.25*Sqr(g1) -  0.25*Sqr(g2) - 0.1*Sqr(g1p)
		    +0.5*( 2.0*(mQlsq- mUrsq+0.125*Sqr(g2)*(Sqr(v1) - Sqr(v2)) - 0.125*Sqr(g1)*(Sqr(v1) - Sqr(v2)))*0.25*(Sqr(g1)- Sqr(g2)) + 4.0*Sqr( At - lambda(3)*s*oneOrt2*v1/v2)* Sqr(r.get_Yu(2,2)) + 8.0*Sqr(mtop)*(At - lambda(3)*s*oneOrt2*v1/v2)*(lambda(3)*s*v1*oneOrt2/(v2*v2*v2) )              ) /(   Sqrt(Sqr(mQlsq- mUrsq+0.125*Sqr(g2)*(Sqr(v1) - Sqr(v2)) - 0.125*Sqr(g1)*(Sqr(v1) - Sqr(v2))) + 4.0*mtop*mtop*Sqr(At - lambda(3)*s*oneOrt2*v1/v2))  ) );
 
 
   Delmstop2 = 0.5*(2.0*Sqr(r.get_Yu(2,2))- 0.6*0.25*Sqr(g1) -  0.25*Sqr(g2) - 0.1*Sqr(g1p)
		    -0.5*( 2.0*(mQlsq- mUrsq+0.125*Sqr(g2)*(Sqr(v1) - Sqr(v2)) - 0.125*Sqr(g1)*(Sqr(v1) - Sqr(v2)))*0.25*(Sqr(g1)- Sqr(g2)) + 4.0*Sqr( At - lambda(3)*s*oneOrt2*v1/v2)* Sqr(r.get_Yu(2,2)) + 8.0*Sqr(mtop)*(At - lambda(3)*s*oneOrt2*v1/v2)*(lambda(3)*s*v1*oneOrt2/(v2*v2*v2) )              ) /(   Sqrt(Sqr(mQlsq- mUrsq+0.125*Sqr(g2)*(Sqr(v1) - Sqr(v2)) - 0.125*Sqr(g1)*(Sqr(v1) - Sqr(v2))) + 4.0*mtop*mtop*Sqr(At - lambda(3)*s*oneOrt2*v1/v2))  ) );
  if(speak){
    cout << "Delmstop1 = " << Delmstop1 << endl;
    cout << "Delmstop2 = " << Delmstop2 << endl;
  } 
  
  
  //checked and debuugged to this point 28/6/07 
  double delta =   threeO32pisq*(-2.0*a0Peter(mstopsq(1),q)*Delmstop1-2.0*a0Peter(mstopsq(2),q)*Delmstop2 + 4.0*a0Peter(Sqr(mtop),q)*Delmtop);
  //cout <<"delta = " << delta << endl;
  return delta;
}
 
double doCalcTadpoleESSMH2_atMt(genericE6SSM_soft_parameters r, double s , double tb  ){
  bool speak = false;
  double yt = r.get_Yu(2,2);
  // cout << "yt = " << yt << endl; 

  double v1 = r.get_vd();
  double v2 = r.get_vu();

  double vev = Sqrt(v1*v1+v2*v2);
  double oneOrt2 = 1/(Sqrt(2.0));
  double mtop = yt*v2/(Sqrt(2.0));
  double q = r.get_scale();
  //fudging Mt and scale to agree with Roman's code 
  //mtop = 165;
  // q = mtop;
  q=165; 
  if(speak){
    cout << "q = " << q << endl;
    cout << "mtop = " << mtop << endl;
  }

  if(USEMTOFMT){
    if(speak){
      cout << "using mtop = " << mtop << endl; 
    }
    mtop = 165;
    yt = mtop/(v2*oneOrt2);
    r.set_Yu(2,2, yt);
  }
  
  double threeO32pisq = 3.0/(32.0*Sqr(PI));
  //double oneOrt2 = 1/(Sqrt(2));			      
  double oneO40 = 1.0/(40.0);
  
  DoubleVector lambda(3), kappa(3), Akappa(3), Alambda(3), mDsq(3), mDbarsq(3);
  
  lambda(3) = r.get_Lambdax();
  
  for (int i = 1; i <= 2; i++)
    {
      lambda(i) = r.get_Lambda12(i-1,i-1);
    }
  
  for (int i = 1; i <= 3; i++)
    {
      kappa(i) = r.get_Kappa(i-1,i-1);
      mDsq(i) = r.get_mDx2(i-1,i-1);
      mDbarsq(i) = r.get_mDxbar2(i-1,i-1);
    }
  
  
  DoubleVector Tlambda(3), Tkappa(3);
  
  Tlambda(3) = r.get_TLambdax();
  
  for (int i = 1; i <= 2; i++)
    {
      Tlambda(i) = r.get_TLambda12(i-1,i-1);
      
      if (Abs(Tlambda(i)) < EPSTOL) 
	{
	  Alambda(i) = 0.0;
	}
      else if (Abs(lambda(i)) < 1.0e-100)
	{
	  ostringstream ii;
	  ii << "WARNING: trying to calculate A_lambda(" << i << ") where lambda(" << i << ") coupling is " <<
	    Abs(lambda(i)) << endl;
	  throw ii.str();
	}
      else
	{
	  Alambda(i) = Tlambda(i)/lambda(i);
	}
      
    }
  
  if (Abs(Tlambda(3)) < EPSTOL) 
    {
      Alambda(3) = 0.0;
    }
  else if (Abs(lambda(3)) < 1.0e-100)
    {
      ostringstream ii;
      ii << "WARNING: trying to calculate A_lambda(3) where lambda(3) coupling is " <<
	Abs(lambda(3)) << endl;
      throw ii.str();
    }
  else
    {
      Alambda(3) = Tlambda(3)/lambda(3);
    }
  
  double At;
  double TYt = r.get_TYu(2,2);
  
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
  
  for (int i = 1; i <= 3; i++)
     {
       Tkappa(i) = r.get_TKappa(i-1,i-1);
       
       if (Abs(Tkappa(i)) < EPSTOL)
	 {
	   Akappa(i) = 0.0;
	 }
       else if (Abs(kappa(i)) < 1.0e-100)
	 {
	   ostringstream ii;
	   ii << "WARNING: trying to calculate A_kappa(" << i << ") where kappa(" << i << ") coupling is " <<
	     Abs(kappa(i)) << endl;
	   throw ii.str();
	 }
       else
	 {
	   Akappa(i) = Tkappa(i)/kappa(i);
	 }
     }

  double Delmtop, Delmstop1,Delmstop2; 

  double mQlsq = r.get_mq2(2,2);
  double mUrsq =  r.get_mu2(2,2);

  double g1 = r.get_g1();
  double g2 = r.get_g2();
  double g1p = r.get_gN();
  //getting stop masses
  DoubleVector mstop(2), mstopsq(2);//stop1,stop2
  DoubleVector mD1sq(3), mD2sq(3);//3 gens
  physical_ESSM(r,mstop, mstopsq, mD1sq,mD2sq, s, tb); 

  if (mstopsq(1) < 0.0 || mstopsq(2) < 0.0)
    {
      mstopsq(1) = Abs(mstopsq(1));
      mstopsq(2) = Abs(mstopsq(2));
    }

  Delmtop = Sqr(r.get_Yu(2,2));
  // cout << "Delmtop = " << Delmtop << endl; 
  Delmstop1 = 0.5*(2.0*Sqr(r.get_Yu(2,2))- 0.6*0.25*Sqr(g1) -  0.25*Sqr(g2) - 0.1*Sqr(g1p)
		   +0.5*( 2.0*(mQlsq- mUrsq+0.125*Sqr(g2)*(Sqr(v1) - Sqr(v2)) - 0.125*Sqr(g1)*(Sqr(v1) - Sqr(v2)))*0.25*(Sqr(g1)- Sqr(g2)) + 4.0*Sqr( At - lambda(3)*s*oneOrt2*v1/v2)* Sqr(r.get_Yu(2,2)) + 8.0*Sqr(mtop)*(At - lambda(3)*s*oneOrt2*v1/v2)*(lambda(3)*s*v1*oneOrt2/(v2*v2*v2) )              ) /(   Sqrt(Sqr(mQlsq- mUrsq+0.125*Sqr(g2)*(Sqr(v1) - Sqr(v2)) - 0.125*Sqr(g1)*(Sqr(v1) - Sqr(v2))) + 4.0*mtop*mtop*Sqr(At - lambda(3)*s*oneOrt2*v1/v2))  ) );
      
  Delmstop2 = 0.5*(2.0*Sqr(r.get_Yu(2,2))- 0.6*0.25*Sqr(g1) -  0.25*Sqr(g2) - 0.1*Sqr(g1p)
		   -0.5*( 2.0*(mQlsq- mUrsq+0.125*Sqr(g2)*(Sqr(v1) - Sqr(v2)) - 0.125*Sqr(g1)*(Sqr(v1) - Sqr(v2)))*0.25*(Sqr(g1)- Sqr(g2)) + 4.0*Sqr( At - lambda(3)*s*oneOrt2*v1/v2)* Sqr(r.get_Yu(2,2)) + 8.0*Sqr(mtop)*(At - lambda(3)*s*oneOrt2*v1/v2)*(lambda(3)*s*v1*oneOrt2/(v2*v2*v2) )              ) /(   Sqrt(Sqr(mQlsq- mUrsq+0.125*Sqr(g2)*(Sqr(v1) - Sqr(v2)) - 0.125*Sqr(g1)*(Sqr(v1) - Sqr(v2))) + 4.0*mtop*mtop*Sqr(At - lambda(3)*s*oneOrt2*v1/v2))  ) );
  if(speak){
    cout << "Delmstop1 = " << Delmstop1 << endl;
    cout << "Delmstop2 = " << Delmstop2 << endl;
  }
  

  //checked and debuugged to this point 28/6/07 
  
  double delta =   threeO32pisq*(-2.0*a0Peter(mstopsq(1),q)*Delmstop1-2.0*a0Peter(mstopsq(2),q)*Delmstop2 + 4.0*a0Peter(Sqr(mtop),q)*Delmtop);
  
  return delta;
}

double doCalcTadpoleESSMH2_Roman(genericE6SSM_soft_parameters r, double s , double tb  ){
  bool speak = false;  
  double yt = r.get_Yu(2,2);
  // cout << "yt = " << yt << endl;
  double v1 = r.get_vd();
  double v2 = r.get_vu();
 
  double vev = Sqrt(v1*v1+v2*v2);

  double oneOrt2 = 1/(Sqrt(2.0));
  double mtop = yt*v2/(Sqrt(2.0));
  mtop = 165;
  yt = mtop/(v2*oneOrt2);
  r.set_Yu(2,2,yt);
  double q = r.get_scale();
  q = 165; //GeV. Fudgeing to match Romans code
  double threeO32pisq = 3.0/(32.0*Sqr(PI));
  //double oneOrt2 = 1/(Sqrt(2));			      

  double oneO40 = 1.0/(40.0);

  DoubleVector lambda(3), kappa(3), Akappa(3), Alambda(3), mDsq(3), mDbarsq(3);
  
  lambda(3) = r.get_Lambdax();
  
  for (int i = 1; i <= 2; i++)
    {
      lambda(i) = r.get_Lambda12(i-1,i-1);
    }
  
  for (int i = 1; i <= 3; i++)
    {
      kappa(i) = r.get_Kappa(i-1,i-1);
      mDsq(i) = r.get_mDx2(i-1,i-1);
      mDbarsq(i) = r.get_mDxbar2(i-1,i-1);
    }
  
  
  DoubleVector Tlambda(3), Tkappa(3);
  
  Tlambda(3) = r.get_TLambdax();
  
  for (int i = 1; i <= 2; i++)
    {
      Tlambda(i) = r.get_TLambda12(i-1,i-1);
      
      if (Abs(Tlambda(i)) < EPSTOL) 
	{
	  Alambda(i) = 0.0;
	}
      else if (Abs(lambda(i)) < 1.0e-100)
	{
	  ostringstream ii;
	  ii << "WARNING: trying to calculate A_lambda(" << i << ") where lambda(" << i << ") coupling is " <<
	    Abs(lambda(i)) << endl;
	  throw ii.str();
	}
      else
	{
	  Alambda(i) = Tlambda(i)/lambda(i);
	}
      
    }
  
  if (Abs(Tlambda(3)) < EPSTOL) 
    {
      Alambda(3) = 0.0;
    }
  else if (Abs(lambda(3)) < 1.0e-100)
    {
      ostringstream ii;
      ii << "WARNING: trying to calculate A_lambda(3) where lambda(3) coupling is " <<
	Abs(lambda(3)) << endl;
      throw ii.str();
    }
  else
    {
      Alambda(3) = Tlambda(3)/lambda(3);
    }
  
  double At;
  double TYt = r.get_TYu(2,2);
  
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
  
  for (int i = 1; i <= 3; i++)
     {
       Tkappa(i) = r.get_TKappa(i-1,i-1);
       
       if (Abs(Tkappa(i)) < EPSTOL)
	 {
	   Akappa(i) = 0.0;
	 }
       else if (Abs(kappa(i)) < 1.0e-100)
	 {
	   ostringstream ii;
	   ii << "WARNING: trying to calculate A_kappa(" << i << ") where kappa(" << i << ") coupling is " <<
	     Abs(kappa(i)) << endl;
	   throw ii.str();
	 }
       else
	 {
	   Akappa(i) = Tkappa(i)/kappa(i);
	 }
     }

  double mQlsq = r.get_mq2(2,2);
  double mUrsq =  r.get_mu2(2,2);
  
  double g1 = r.get_g1();
  double g2 = r.get_g2();
  double g1p = r.get_gN();
  double Delmtop, Delmstop1,Delmstop2; 
  //   cout << "in H2 tads" << endl; 
  
  //getting stop masses
  DoubleVector mstop(2), mstopsq(2);//stop1,stop2
  DoubleVector mD1sq(3), mD2sq(3);//3 gens
  physical_ESSM_Roman(r,mstop, mstopsq, mD1sq,mD2sq, s, tb); 

  if (mstopsq(1) < 0.0 || mstopsq(2) < 0.0)
    {
      mstopsq(1) = Abs(mstopsq(1));
      mstopsq(2) = Abs(mstopsq(2));
    }

  Delmtop = Sqr(r.get_Yu(2,2));
  //cout << "Delmtop = " << Delmtop << endl; 
  
  Delmstop1 = 0.5*(2.0*Sqr(r.get_Yu(2,2))  +0.5*(  4.0*Sqr( At - lambda(3)*s*oneOrt2*v1/v2)* Sqr(r.get_Yu(2,2))+ 8.0*Sqr(mtop)*(At - lambda(3)*s*oneOrt2*v1/v2)*(lambda(3)*s*v1*oneOrt2/(v2*v2*v2) )             ) /(   Sqrt(Sqr(mQlsq- mUrsq) + 4.0*mtop*mtop*Sqr(At - lambda(3)*s*oneOrt2*v1/v2))  ) ) ;
  
  Delmstop2 = 0.5*(2.0*Sqr(r.get_Yu(2,2))
		   -0.5*(  4.0*Sqr( At - lambda(3)*s*oneOrt2*v1/v2)* Sqr(r.get_Yu(2,2)) + 8.0*Sqr(mtop)*(At - lambda(3)*s*oneOrt2*v1/v2)*(lambda(3)*s*v1*oneOrt2/(v2*v2*v2) )              ) /(   Sqrt(Sqr(mQlsq- mUrsq)+ 4.0*mtop*mtop*Sqr(At - lambda(3)*s*oneOrt2*v1/v2))       ));
  
  
  //checked and debuugged to this point 28/6/07 
  double delta =   threeO32pisq*(-2.0*a0Peter(mstopsq(1),q)*Delmstop1-2.0*a0Peter(mstopsq(2),q)*Delmstop2 + 4.0*a0Peter(Sqr(mtop),q)*Delmtop);
  //cout <<"delta = " << delta << endl;
  return delta;
}

double doCalcTadpoleESSMH2_Roman_atQ(genericE6SSM_soft_parameters r, double s , double tb  ){
  double yt = r.get_Yu(2,2);
  // cout << "yt = " << yt << endl; 
  double v1 = r.get_vd();
  double v2 = r.get_vu();

  double vev = Sqrt(v1*v1+v2*v2);//246; //NOTE CHANGE

  double oneOrt2 = 1/(Sqrt(2.0));
  double mtop = yt*v2/(Sqrt(2.0));
  //  mtop = 165;
  yt = mtop/(v2*oneOrt2);
  r.set_Yu(2,2,yt);
  double q = r.get_scale();
  
  if(USEMTOFMT){
    cout << "using mtop = " << mtop << endl; 
    mtop = 165;
    yt = mtop/(v2*oneOrt2);
    r.set_Yu(2,2,yt);
    
  }
  
  double threeO32pisq = 3.0/(32.0*Sqr(PI));
  //double oneOrt2 = 1/(Sqrt(2));			      

  double oneO40 = 1.0/(40.0);

  DoubleVector lambda(3), kappa(3), Akappa(3), Alambda(3), mDsq(3), mDbarsq(3);
  
  lambda(3) = r.get_Lambdax();
  
  for (int i = 1; i <= 2; i++)
    {
      lambda(i) = r.get_Lambda12(i-1,i-1);
    }
  
  for (int i = 1; i <= 3; i++)
    {
      kappa(i) = r.get_Kappa(i-1,i-1);
      mDsq(i) = r.get_mDx2(i-1,i-1);
      mDbarsq(i) = r.get_mDxbar2(i-1,i-1);
    }
  
  
  DoubleVector Tlambda(3), Tkappa(3);
  
  Tlambda(3) = r.get_TLambdax();
  
  for (int i = 1; i <= 2; i++)
    {
      Tlambda(i) = r.get_TLambda12(i-1,i-1);
      
      if (Abs(Tlambda(i)) < EPSTOL) 
	{
	  Alambda(i) = 0.0;
	}
      else if (Abs(lambda(i)) < 1.0e-100)
	{
	  ostringstream ii;
	  ii << "WARNING: trying to calculate A_lambda(" << i << ") where lambda(" << i << ") coupling is " <<
	    Abs(lambda(i)) << endl;
	  throw ii.str();
	}
      else
	{
	  Alambda(i) = Tlambda(i)/lambda(i);
	}
      
    }
  
  if (Abs(Tlambda(3)) < EPSTOL) 
    {
      Alambda(3) = 0.0;
    }
  else if (Abs(lambda(3)) < 1.0e-100)
    {
      ostringstream ii;
      ii << "WARNING: trying to calculate A_lambda(3) where lambda(3) coupling is " <<
	Abs(lambda(3)) << endl;
      throw ii.str();
    }
  else
    {
      Alambda(3) = Tlambda(3)/lambda(3);
    }
  
  double At;
  double TYt = r.get_TYu(2,2);
  
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
  
  for (int i = 1; i <= 3; i++)
     {
       Tkappa(i) = r.get_TKappa(i-1,i-1);
       
       if (Abs(Tkappa(i)) < EPSTOL)
	 {
	   Akappa(i) = 0.0;
	 }
       else if (Abs(kappa(i)) < 1.0e-100)
	 {
	   ostringstream ii;
	   ii << "WARNING: trying to calculate A_kappa(" << i << ") where kappa(" << i << ") coupling is " <<
	     Abs(kappa(i)) << endl;
	   throw ii.str();
	 }
       else
	 {
	   Akappa(i) = Tkappa(i)/kappa(i);
	 }
     }

  double Delmtop, Delmstop1,Delmstop2; 
  double mQlsq = r.get_mq2(2,2);
  double mUrsq =  r.get_mu2(2,2);
  
  double g1 = r.get_g1();
  double g2 = r.get_g2();
  double g1p = r.get_gN();

  
  DoubleVector mstop(2), mstopsq(2);//stop1,stop2
  DoubleVector mD1sq(3), mD2sq(3);//3 gens
  
  physical_ESSM_Roman(r,mstop, mstopsq, mD1sq,mD2sq, s, tb);
   
  if (mstopsq(1) < 0.0 || mstopsq(2) < 0.0)
    {
      mstopsq(1) = Abs(mstopsq(1));
      mstopsq(2) = Abs(mstopsq(2));
    }

  Delmtop = Sqr(r.get_Yu(2,2));
  //cout << "Delmtop = " << Delmtop << endl; 
  
  Delmstop1 = 0.5*(2.0*Sqr(r.get_Yu(2,2))  +0.5*(  4.0*Sqr( At - lambda(3)*s*oneOrt2*v1/v2)* Sqr(r.get_Yu(2,2))+ 8.0*Sqr(mtop)*(At - lambda(3)*s*oneOrt2*v1/v2)*(lambda(3)*s*v1*oneOrt2/(v2*v2*v2) )             ) /(   Sqrt(Sqr(mQlsq- mUrsq) + 4.0*mtop*mtop*Sqr(At - lambda(3)*s*oneOrt2*v1/v2))  ) ) ;
  
  Delmstop2 = 0.5*(2.0*Sqr(r.get_Yu(2,2))
		   -0.5*(  4.0*Sqr( At - lambda(3)*s*oneOrt2*v1/v2)* Sqr(r.get_Yu(2,2)) + 8.0*Sqr(mtop)*(At - lambda(3)*s*oneOrt2*v1/v2)*(lambda(3)*s*v1*oneOrt2/(v2*v2*v2) )              ) /(   Sqrt(Sqr(mQlsq- mUrsq)+ 4.0*mtop*mtop*Sqr(At - lambda(3)*s*oneOrt2*v1/v2))       ));
  
  
   //checked and debuugged to this point 28/6/07 
   double delta =   threeO32pisq*(-2.0*a0Peter(mstopsq(1),q)*Delmstop1-2.0*a0Peter(mstopsq(2),q)*Delmstop2 + 4.0*a0Peter(Sqr(mtop),q)*Delmtop);
   //cout <<"delta = " << delta << endl;
   return delta;   
}

//Tadpoles for S.  same labeling as above.
double doCalcTadpolesESSMS( genericE6SSM_soft_parameters r, double s , double tb  ){
  bool speak = false;  
  double threeO32pisq = 3.0/(32.0*Sqr(PI));
  double oneO40 = 1.0/(40.0);
  double oneOrt2 = 1/(Sqrt(2.0));
  double yt = r.get_Yu(2,2);

  double v1 = r.get_vd();
  double v2 = r.get_vu();

  double vev = Sqrt(v1*v1+v2*v2);
  //double oneOrt2 = 1/(Sqrt(2.0));
  double mtop = yt*v2*oneOrt2;
  // cout << "mtop = " << mtop << endl;
  double q = r.get_scale();
  
  if(USEMTOFMT){
    //cout << "using mtop = " << mtop << endl;
    mtop = 165;
    yt = mtop/(v2*oneOrt2);
    r.set_Yu(2,2,yt);
  }

  DoubleVector lambda(3), kappa(3), Akappa(3), Alambda(3), mDsq(3), mDbarsq(3);
  
  lambda(3) = r.get_Lambdax();
  
  for (int i = 1; i <= 2; i++)
    {
      lambda(i) = r.get_Lambda12(i-1,i-1);
    }
  
  for (int i = 1; i <= 3; i++)
    {
      kappa(i) = r.get_Kappa(i-1,i-1);
      mDsq(i) = r.get_mDx2(i-1,i-1);
      mDbarsq(i) = r.get_mDxbar2(i-1,i-1);
    }
  
  
  DoubleVector Tlambda(3), Tkappa(3);
  
  Tlambda(3) = r.get_TLambdax();
  
  for (int i = 1; i <= 2; i++)
    {
      Tlambda(i) = r.get_TLambda12(i-1,i-1);
      
      if (Abs(Tlambda(i)) < EPSTOL) 
	{
	  Alambda(i) = 0.0;
	}
      else if (Abs(lambda(i)) < 1.0e-100)
	{
	  ostringstream ii;
	  ii << "WARNING: trying to calculate A_lambda(" << i << ") where lambda(" << i << ") coupling is " <<
	    Abs(lambda(i)) << endl;
	  throw ii.str();
	}
      else
	{
	  Alambda(i) = Tlambda(i)/lambda(i);
	}
      
    }
  
  if (Abs(Tlambda(3)) < EPSTOL) 
    {
      Alambda(3) = 0.0;
    }
  else if (Abs(lambda(3)) < 1.0e-100)
    {
      ostringstream ii;
      ii << "WARNING: trying to calculate A_lambda(3) where lambda(3) coupling is " <<
	Abs(lambda(3)) << endl;
      throw ii.str();
    }
  else
    {
      Alambda(3) = Tlambda(3)/lambda(3);
    }
  
  double At;
  double TYt = r.get_TYu(2,2);
  
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
  
  for (int i = 1; i <= 3; i++)
     {
       Tkappa(i) = r.get_TKappa(i-1,i-1);
       
       if (Abs(Tkappa(i)) < EPSTOL)
	 {
	   Akappa(i) = 0.0;
	 }
       else if (Abs(kappa(i)) < 1.0e-100)
	 {
	   ostringstream ii;
	   ii << "WARNING: trying to calculate A_kappa(" << i << ") where kappa(" << i << ") coupling is " <<
	     Abs(kappa(i)) << endl;
	   throw ii.str();
	 }
       else
	 {
	   Akappa(i) = Tkappa(i)/kappa(i);
	 }
     }

  double mQlsq = r.get_mq2(2,2);
  double mUrsq =  r.get_mu2(2,2);
  
  double g1 = r.get_g1();
  double g2 = r.get_g2();
  double g1p = r.get_gN();
  
  DoubleVector mstop(2),  mstopsq(2);//stop1,stop2
  DoubleVector mD1sq(3),mD2sq(3);//3 gens

  physical_ESSM(r, mstop,  mstopsq, mD1sq, mD2sq, s, tb);

  if (mstopsq(1) < 0.0 || mstopsq(2) < 0.0)
    {
      mstopsq(1) = Abs(mstopsq(1));
      mstopsq(2) = Abs(mstopsq(2));
    }

  DoubleVector DelmuDsq(3), DelmD1sq(3),DelmD2sq(3);
  double Delmtop,Delmstop1,Delmstop2;  

  Delmtop = 0;
  
  Delmstop1 = 0.5*(0.25*Sqr(g1p)
		   +0.5*( 8.0*Sqr(mtop)*(At - lambda(3)*s*oneOrt2*v1/v2)*(-lambda(3)*v1*oneOrt2/(v2*s) ) ) /(   Sqrt(Sqr(mQlsq- mUrsq+0.125*Sqr(g2)*(Sqr(v1) - Sqr(v2)) - 0.125*Sqr(g1)*(Sqr(v1) - Sqr(v2))) + 4.0*mtop*mtop*Sqr(At - lambda(3)*s*oneOrt2*v1/v2))  ));
  
  Delmstop2 = 0.5*(0.25*Sqr(g1p)
		   -0.5*( 8.0*Sqr(mtop)*(At - lambda(3)*s*oneOrt2*v1/v2)*(-lambda(3)*v1*oneOrt2/(v2*s) ) ) /(   Sqrt(Sqr(mQlsq- mUrsq+0.125*Sqr(g2)*(Sqr(v1) - Sqr(v2)) - 0.125*Sqr(g1)*(Sqr(v1) - Sqr(v2))) + 4.0*mtop*mtop*Sqr(At - lambda(3)*s*oneOrt2*v1/v2))  ));

  if(speak){
    cout << "Delmstop1 = " << Delmstop1 << endl;
    cout << "Delmstop2 = " << Delmstop2 << endl;
  }
 
  double delta1 =   threeO32pisq*(-2.0*a0(mstop(1),q)*Delmstop1-2.0*a0(mstop(2),q)*Delmstop2 + 4.0*a0(mtop,q)*Delmtop);
  double delta =   threeO32pisq*(-2.0*a0Peter(mstopsq(1),q)*Delmstop1-2.0*a0Peter(mstopsq(2),q)*Delmstop2 + 4.0*a0Peter(Sqr(mtop),q)*Delmtop);
  //cout <<"delta = " << delta << endl;
  return delta;
}

double doCalcTadpolesESSMS_atMt( genericE6SSM_soft_parameters r, double s , double tb  )  {
  bool speak = false;
  double threeO32pisq = 3.0/(32.0*Sqr(PI));
  double oneO40 = 1.0/(40.0);
  double oneOrt2 = 1/(Sqrt(2.0));
  double yt = r.get_Yu(2,2);

  double v1 = r.get_vd();
  double v2 = r.get_vu();

  double vev = Sqrt(v1*v1+v2*v2);

  //double oneOrt2 = 1/(Sqrt(2.0));
  double mtop = yt*v2*oneOrt2;
  // cout << "mtop = " << mtop << endl;
  double q = r.get_scale();
  mtop = 165;
  q = mtop;

  
  if(USEMTOFMT){
    mtop = 165;
    yt = mtop/(v2*oneOrt2);
    r.set_Yu(2,2,yt);
    
  }

  if(speak){
   cout << "q = " << q << endl;
  cout << "mtop = " << mtop << endl;
  }


  DoubleVector lambda(3), kappa(3), Akappa(3), Alambda(3), mDsq(3), mDbarsq(3);
  
  lambda(3) = r.get_Lambdax();
  
  for (int i = 1; i <= 2; i++)
    {
      lambda(i) = r.get_Lambda12(i-1,i-1);
    }
  
  for (int i = 1; i <= 3; i++)
    {
      kappa(i) = r.get_Kappa(i-1,i-1);
      mDsq(i) = r.get_mDx2(i-1,i-1);
      mDbarsq(i) = r.get_mDxbar2(i-1,i-1);
    }
  
  
  DoubleVector Tlambda(3), Tkappa(3);
  
  Tlambda(3) = r.get_TLambdax();
  
  for (int i = 1; i <= 2; i++)
    {
      Tlambda(i) = r.get_TLambda12(i-1,i-1);
      
      if (Abs(Tlambda(i)) < EPSTOL) 
	{
	  Alambda(i) = 0.0;
	}
      else if (Abs(lambda(i)) < 1.0e-100)
	{
	  ostringstream ii;
	  ii << "WARNING: trying to calculate A_lambda(" << i << ") where lambda(" << i << ") coupling is " <<
	    Abs(lambda(i)) << endl;
	  throw ii.str();
	}
      else
	{
	  Alambda(i) = Tlambda(i)/lambda(i);
	}
      
    }
  
  if (Abs(Tlambda(3)) < EPSTOL) 
    {
      Alambda(3) = 0.0;
    }
  else if (Abs(lambda(3)) < 1.0e-100)
    {
      ostringstream ii;
      ii << "WARNING: trying to calculate A_lambda(3) where lambda(3) coupling is " <<
	Abs(lambda(3)) << endl;
      throw ii.str();
    }
  else
    {
      Alambda(3) = Tlambda(3)/lambda(3);
    }
  
  double At;
  double TYt = r.get_TYu(2,2);
  
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
  
  for (int i = 1; i <= 3; i++)
    {
      Tkappa(i) = r.get_TKappa(i-1,i-1);
      
      if (Abs(Tkappa(i)) < EPSTOL)
	{
	  Akappa(i) = 0.0;
	}
      else if (Abs(kappa(i)) < 1.0e-100)
	{
	  ostringstream ii;
	  ii << "WARNING: trying to calculate A_kappa(" << i << ") where kappa(" << i << ") coupling is " <<
	    Abs(kappa(i)) << endl;
	  throw ii.str();
	}
      else
	{
	  Akappa(i) = Tkappa(i)/kappa(i);
	}
    }
  
  
  double mQlsq = r.get_mq2(2,2);
  double mUrsq =  r.get_mu2(2,2);

  double g1 = r.get_g1();
  double g2 = r.get_g2();
  double g1p = r.get_gN();
  DoubleVector mstop(2),  mstopsq(2);//stop1,stop2
  DoubleVector mD1sq(3),mD2sq(3);//3 gens

  physical_ESSM(r, mstop,  mstopsq, mD1sq, mD2sq, s, tb);

  if (mstopsq(1) < 0.0 || mstopsq(2) < 0.0)
    {
      mstopsq(1) = Abs(mstopsq(1));
      mstopsq(2) = Abs(mstopsq(2));
    }

  double Delmtop,Delmstop1,Delmstop2;  
  Delmtop = 0;
  Delmstop1 = 0.5*(0.25*Sqr(g1p)
		   +0.5*( 8.0*Sqr(mtop)*(At - lambda(3)*s*oneOrt2*v1/v2)*(-lambda(3)*v1*oneOrt2/(v2*s) ) ) /(   Sqrt(Sqr(mQlsq- mUrsq+0.125*Sqr(g2)*(Sqr(v1) - Sqr(v2)) - 0.125*Sqr(g1)*(Sqr(v1) - Sqr(v2))) + 4.0*mtop*mtop*Sqr(At - lambda(3)*s*oneOrt2*v1/v2))  ));
 
  Delmstop2 = 0.5*(0.25*Sqr(g1p)
		   -0.5*( 8.0*Sqr(mtop)*(At - lambda(3)*s*oneOrt2*v1/v2)*(-lambda(3)*v1*oneOrt2/(v2*s) ) ) /(   Sqrt(Sqr(mQlsq- mUrsq+0.125*Sqr(g2)*(Sqr(v1) - Sqr(v2)) - 0.125*Sqr(g1)*(Sqr(v1) - Sqr(v2))) + 4.0*mtop*mtop*Sqr(At - lambda(3)*s*oneOrt2*v1/v2))  ));

  if(speak){
    cout << "Delmstop1 = " << Delmstop1 << endl;
    cout << "Delmstop2 = " << Delmstop2 << endl;
  }
  

  double delta =   threeO32pisq*(-2.0*a0Peter(mstopsq(1),q)*Delmstop1-2.0*a0Peter(mstopsq(2),q)*Delmstop2 + 4.0*a0Peter(Sqr(mtop),q)*Delmtop);
  
  //cout <<"delta = " << delta << endl;
  return delta;      
}
   

double doCalcTadpolesESSMS_Roman( genericE6SSM_soft_parameters r, double s , double tb  )  {
     
  double threeO32pisq = 3.0/(32.0*Sqr(PI));
  double oneO40 = 1.0/(40.0);
  double oneOrt2 = 1/(Sqrt(2.0));
  double yt = r.get_Yu(2,2);

  double v1 = r.get_vd();
  double v2 = r.get_vu();

  double vev = Sqrt(v1*v1+v2*v2);

  ////double oneOrt2 = 1/(Sqrt(2.0));
  double mtop = yt*v2*oneOrt2;
  mtop = 165;
  
  yt = mtop/(v2*oneOrt2);
  r.set_Yu(2,2,yt);
  //cout << "mtop = " << mtop << endl;
  double q = r.get_scale();
  q = 165; //GeV. Fudgeing to match Romans code
  //cout << "q = " << q << endl;
  
  DoubleVector lambda(3), kappa(3), Akappa(3), Alambda(3), mDsq(3), mDbarsq(3);
  
  lambda(3) = r.get_Lambdax();
  
  for (int i = 1; i <= 2; i++)
    {
      lambda(i) = r.get_Lambda12(i-1,i-1);
    }
  
  for (int i = 1; i <= 3; i++)
    {
      kappa(i) = r.get_Kappa(i-1,i-1);
      mDsq(i) = r.get_mDx2(i-1,i-1);
      mDbarsq(i) = r.get_mDxbar2(i-1,i-1);
    }
  
  
  DoubleVector Tlambda(3), Tkappa(3);
  
  Tlambda(3) = r.get_TLambdax();
  
  for (int i = 1; i <= 2; i++)
    {
      Tlambda(i) = r.get_TLambda12(i-1,i-1);
      
      if (Abs(Tlambda(i)) < EPSTOL) 
	{
	  Alambda(i) = 0.0;
	}
      else if (Abs(lambda(i)) < 1.0e-100)
	{
	  ostringstream ii;
	  ii << "WARNING: trying to calculate A_lambda(" << i << ") where lambda(" << i << ") coupling is " <<
	    Abs(lambda(i)) << endl;
	  throw ii.str();
	}
      else
	{
	  Alambda(i) = Tlambda(i)/lambda(i);
	}
      
    }
  
  if (Abs(Tlambda(3)) < EPSTOL) 
    {
      Alambda(3) = 0.0;
    }
  else if (Abs(lambda(3)) < 1.0e-100)
    {
      ostringstream ii;
      ii << "WARNING: trying to calculate A_lambda(3) where lambda(3) coupling is " <<
	Abs(lambda(3)) << endl;
      throw ii.str();
    }
  else
    {
      Alambda(3) = Tlambda(3)/lambda(3);
    }
  
  double At;
  double TYt = r.get_TYu(2,2);
  
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
  
  for (int i = 1; i <= 3; i++)
    {
      Tkappa(i) = r.get_TKappa(i-1,i-1);
      
      if (Abs(Tkappa(i)) < EPSTOL)
	{
	  Akappa(i) = 0.0;
	}
      else if (Abs(kappa(i)) < 1.0e-100)
	{
	  ostringstream ii;
	  ii << "WARNING: trying to calculate A_kappa(" << i << ") where kappa(" << i << ") coupling is " <<
	    Abs(kappa(i)) << endl;
	  throw ii.str();
	}
      else
	{
	  Akappa(i) = Tkappa(i)/kappa(i);
	}
    }

  double mQlsq = r.get_mq2(2,2);
  double mUrsq =  r.get_mu2(2,2);
  
  double g1 = r.get_g1();
  double g2 = r.get_g2();
  double g1p = r.get_gN(); 
  DoubleVector mstop(2),  mstopsq(2);//stop1,stop2
  DoubleVector mD1sq(3),mD2sq(3);//3 gens
  physical_ESSM_Roman(r, mstop,  mstopsq, mD1sq, mD2sq, s, tb);
  
  if (mstopsq(1) < 0.0 || mstopsq(2) < 0.0)
    {
      mstopsq(1) = Abs(mstopsq(1));
      mstopsq(2) = Abs(mstopsq(2));
    }

  DoubleVector DelmuDsq(3), DelmD1sq(3),DelmD2sq(3);
  double Delmtop,Delmstop1,Delmstop2;  
  
  Delmtop = 0;
  
  Delmstop1 = 0.5*(0.5*( 8.0*Sqr(mtop)*(At - lambda(3)*s*oneOrt2*v1/v2)*(-lambda(3)*v1*oneOrt2/(v2*s) ) )) /(   Sqrt(Sqr(mQlsq- mUrsq)  + 4.0*mtop*mtop*Sqr(At - lambda(3)*s*oneOrt2*v1/v2))  );
  
  Delmstop2 = 0.5*( -0.5*( 8.0*Sqr(mtop)*(At - lambda(3)*s*oneOrt2*v1/v2)*(-lambda(3)*v1*oneOrt2/(v2*s) ) )) /(   Sqrt(Sqr(mQlsq- mUrsq)+ 4.0*mtop*mtop*Sqr(At - lambda(3)*s*oneOrt2*v1/v2))  );
  

  double delta =   threeO32pisq*(-2.0*a0Peter(mstopsq(1),q)*Delmstop1-2.0*a0Peter(mstopsq(2),q)*Delmstop2 + 4.0*a0Peter(Sqr(mtop),q)*Delmtop);
  return delta;
}


double doCalcTadpolesESSMS_Roman_atQ( genericE6SSM_soft_parameters r, double s , double tb  ){
  double threeO32pisq = 3.0/(32.0*Sqr(PI));
  double oneO40 = 1.0/(40.0);
  double oneOrt2 = 1/(Sqrt(2.0));
  double yt = r.get_Yu(2,2);

  double v1 = r.get_vd();
  double v2 = r.get_vu();

  double vev = Sqrt(v1*v1+v2*v2);

  ////double oneOrt2 = 1/(Sqrt(2.0));
  double mtop = yt*v2*oneOrt2;
  //mtop = 165;
  
  yt = mtop/(v2*oneOrt2);
  r.set_Yu(2,2,yt);
  //cout << "mtop = " << mtop << endl;
  double q = r.get_scale();
  //q = 165; //GeV. Fudgeing to match Romans code
  cout << "q = " << q << endl;
  
  if(USEMTOFMT){
    cout << "using mtop = " << mtop << endl; 
    mtop = 165;
    yt = mtop/(v2*oneOrt2);
    r.set_Yu(2,2,yt);
  }
  
  DoubleVector lambda(3), kappa(3), Akappa(3), Alambda(3), mDsq(3), mDbarsq(3);
  
  lambda(3) = r.get_Lambdax();
  
  for (int i = 1; i <= 2; i++)
    {
      lambda(i) = r.get_Lambda12(i-1,i-1);
    }
  
  for (int i = 1; i <= 3; i++)
    {
      kappa(i) = r.get_Kappa(i-1,i-1);
      mDsq(i) = r.get_mDx2(i-1,i-1);
      mDbarsq(i) = r.get_mDxbar2(i-1,i-1);
    }
  
  
  DoubleVector Tlambda(3), Tkappa(3);
  
  Tlambda(3) = r.get_TLambdax();
  
  for (int i = 1; i <= 2; i++)
    {
      Tlambda(i) = r.get_TLambda12(i-1,i-1);
      
      if (Abs(Tlambda(i)) < EPSTOL) 
	{
	  Alambda(i) = 0.0;
	}
      else if (Abs(lambda(i)) < 1.0e-100)
	{
	  ostringstream ii;
	  ii << "WARNING: trying to calculate A_lambda(" << i << ") where lambda(" << i << ") coupling is " <<
	    Abs(lambda(i)) << endl;
	  throw ii.str();
	}
      else
	{
	  Alambda(i) = Tlambda(i)/lambda(i);
	}
      
    }
  
  if (Abs(Tlambda(3)) < EPSTOL) 
    {
      Alambda(3) = 0.0;
    }
  else if (Abs(lambda(3)) < 1.0e-100)
    {
      ostringstream ii;
      ii << "WARNING: trying to calculate A_lambda(3) where lambda(3) coupling is " <<
	Abs(lambda(3)) << endl;
      throw ii.str();
    }
  else
    {
      Alambda(3) = Tlambda(3)/lambda(3);
    }
  
  double At;
  double TYt = r.get_TYu(2,2);
  
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
  
  for (int i = 1; i <= 3; i++)
    {
      Tkappa(i) = r.get_TKappa(i-1,i-1);
      
      if (Abs(Tkappa(i)) < EPSTOL)
	{
	  Akappa(i) = 0.0;
	}
      else if (Abs(kappa(i)) < 1.0e-100)
	{
	  ostringstream ii;
	  ii << "WARNING: trying to calculate A_kappa(" << i << ") where kappa(" << i << ") coupling is " <<
	    Abs(kappa(i)) << endl;
	  throw ii.str();
	}
      else
	{
	  Akappa(i) = Tkappa(i)/kappa(i);
	}
    }

  double mQlsq = r.get_mq2(2,2);
  double mUrsq =  r.get_mu2(2,2);
  
  double g1 = r.get_g1();
  double g2 = r.get_g2();
  double g1p = r.get_gN();
  DoubleVector mstop(2),  mstopsq(2);//stop1,stop2
  DoubleVector mD1sq(3),mD2sq(3);//3 gens
  physical_ESSM_Roman(r, mstop,  mstopsq, mD1sq, mD2sq, s, tb);

  if (mstopsq(1) < 0.0 || mstopsq(2) < 0.0)
    {
      mstopsq(1) = Abs(mstopsq(1));
      mstopsq(2) = Abs(mstopsq(2));
    }

  double Delmtop,Delmstop1,Delmstop2;  
  Delmtop = 0;
  
  Delmstop1 = 0.5*(0.5*( 8.0*Sqr(mtop)*(At - lambda(3)*s*oneOrt2*v1/v2)*(-lambda(3)*v1*oneOrt2/(v2*s) ) )) /(   Sqrt(Sqr(mQlsq- mUrsq)  + 4.0*mtop*mtop*Sqr(At - lambda(3)*s*oneOrt2*v1/v2))  );
  
  Delmstop2 = 0.5*( -0.5*( 8.0*Sqr(mtop)*(At - lambda(3)*s*oneOrt2*v1/v2)*(-lambda(3)*v1*oneOrt2/(v2*s) ) )) /(   Sqrt(Sqr(mQlsq- mUrsq)+ 4.0*mtop*mtop*Sqr(At - lambda(3)*s*oneOrt2*v1/v2))  );
  
  
  double delta =   threeO32pisq*(-2.0*a0Peter(mstopsq(1),q)*Delmstop1-2.0*a0Peter(mstopsq(2),q)*Delmstop2 + 4.0*a0Peter(Sqr(mtop),q)*Delmtop);
  //cout <<"delta = " << delta << endl;
 
  return delta;
}

bool HiggsMasses(genericE6SSM_soft_parameters const & model, double s, double tb, DoubleVector & mstop, DoubleVector & mstopsq, int WhatCorrections, bool speak, bool Bugspeak, DoubleVector & bounds, int & ExpValid, DoubleVector & mhout, DoubleMatrix & mhmix, DoubleMatrix & msq, int & sing) {


  // Copy object for calculation
  genericE6SSM_soft_parameters r = model;

  bool higgsTachyon = false;
  ExpValid = 0; //< assumes no problems initially

  double mQLsq = r.get_mq2(2,2);
  double mURsq = r.get_mu2(2,2);
  double g2 = r.get_g2();
  double g1 = r.get_g1();

  double yt = r.get_Yu(2,2); //cout << "yt = " << yt << endl;
  double g1p = r.get_gN();

  double gbar = Sqrt(3.0*Sqr(g1)/5.0 + Sqr(g2)) ; 


  double v1 = r.get_vd();
  double v2 = r.get_vu();
 
  double vev = Sqrt(v1*v1+v2*v2);

  double oneOrt2 = 1.0/Sqrt(2.0);

  DoubleVector mylambda(3), mykappa(3), myAkappa(3), myAlambda(3), mDsq(3), mDbarsq(3);
  
  mylambda(3) = r.get_Lambdax();
  
  for (int i = 1; i <= 2; i++)
    {
      mylambda(i) = r.get_Lambda12(i-1,i-1);
    }
  
  for (int i = 1; i <= 3; i++)
    {
      mykappa(i) = r.get_Kappa(i-1,i-1);
      mDsq(i) = r.get_mDx2(i-1,i-1);
      mDbarsq(i) = r.get_mDxbar2(i-1,i-1);
    }
  
  
  DoubleVector Tlambda(3), Tkappa(3);
  
  Tlambda(3) = r.get_TLambdax();
  
  for (int i = 1; i <= 2; i++)
    {
      Tlambda(i) = r.get_TLambda12(i-1,i-1);
      
      if (Abs(Tlambda(i)) < EPSTOL) 
	{
	  myAlambda(i) = 0.0;
	}
      else if (Abs(mylambda(i)) < 1.0e-100)
	{
	  ostringstream ii;
	  ii << "WARNING: trying to calculate A_lambda(" << i << ") where lambda(" << i << ") coupling is " <<
	    Abs(mylambda(i)) << endl;
	  throw ii.str();
	}
      else
	{
	  myAlambda(i) = Tlambda(i)/mylambda(i);
	}
      
    }
  
  if (Abs(Tlambda(3)) < EPSTOL) 
    {
      myAlambda(3) = 0.0;
    }
  else if (Abs(mylambda(3)) < 1.0e-100)
    {
      ostringstream ii;
      ii << "WARNING: trying to calculate A_lambda(3) where lambda(3) coupling is " <<
	Abs(mylambda(3)) << endl;
      throw ii.str();
    }
  else
    {
      myAlambda(3) = Tlambda(3)/mylambda(3);
    }
  
  double At;
  double TYt = r.get_TYu(2,2);
  
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
  
  for (int i = 1; i <= 3; i++)
    {
      Tkappa(i) = r.get_TKappa(i-1,i-1);
      
      if (Abs(Tkappa(i)) < EPSTOL)
	{
	  myAkappa(i) = 0.0;
	}
      else if (Abs(mykappa(i)) < 1.0e-100)
	{
	  ostringstream ii;
	  ii << "WARNING: trying to calculate A_kappa(" << i << ") where kappa(" << i << ") coupling is " <<
	    Abs(mykappa(i)) << endl;
	  throw ii.str();
	}
      else
	{
	  myAkappa(i) = Tkappa(i)/mykappa(i);
	}
    }


  if( WhatCorrections ==2 || WhatCorrections == 4){
    r.set_scale(165);
  }

  double mu_eff = mylambda(3)*s/(Sqrt(2.0));

  genericE6SSM_input_parameters input = r.get_input();

  const bool UseEffectiveChargesIntreeEWSB = false;
  double Q1eff = input.QH1p, Q2eff = input.QH2p, QSeff = input.QSp;

  // DH: currently neglecting U(1) mixing and therefore off-diagonal
  // gauge coupling g_11 = 0.
  // if(UseEffectiveChargesIntreeEWSB){
  //   Q1eff = -3.0/(Sqrt(40.0)) - 0.5*Sqrt(3.0/5.0)*r.displayg_11()/(g1p);
  //   Q2eff = -2.0 / (Sqrt(40.0)) + 0.5*Sqrt(3.0/5.0)*r.displayg_11()/(g1p);
  //   QSeff = 5.0/(Sqrt(40.0));  
  // }
  // else{  
  //   Q1eff = input.QH1p;
  //   Q2eff = input.QH2p;
  //   QSeff = input.Qsp;
  // } 

  DoubleVector mD1sq(3), mD2sq(3);
  
 if(WhatCorrections == 1 || WhatCorrections == 2){
   physical_ESSM(r,mstop,mstopsq, mD1sq,mD2sq, s, tb);
 }
 else{   
   physical_ESSM_Roman(r,mstop,mstopsq, mD1sq,mD2sq, s, tb);
 }

 if((mstopsq(2) < 0)|| (mstopsq(1) < 0)){
   if(speak) cerr << "tachyonics stop masses so i'm not going to bother doing higgs mass corrections.  Point already ruled out. " << endl;

     // Dylan:: I have modified this so that in the event of tachyonic stops, the returned Higgs masses are zero
     // as a flag. Count it as a tadpole problem.
     mhout.set(1, 0.0);
     mhout.set(2, 0.0);
     mhout.set(3, 0.0);
     ExpValid = TADPOLESPROBLEM;
     return higgsTachyon;
   }  

 //The stops from physical have the opposite mass ordering
 // to that used in physicalo and the tadpoles
 //to the convention used here so now we mass order the stop
  
   if (mstop(1) > mstop(2)){
      double stoptemp = mstop(1);
      mstop(1) = mstop(2);
      mstop(2) = stoptemp;
    }
 
  if(mstopsq(1) > mstopsq(2)){
      double stoptempsq = mstopsq(1);
      mstopsq(1) = mstopsq(2);
      mstopsq(2) = stoptempsq;
    }
 double mtop = yt*v2*oneOrt2;
 mtop = 165;
 yt = mtop/(v2*oneOrt2);
 r.set_Yu(2,2, yt);
 //CP-Even Higgs 1 loop Corrections in the basis shown in Theory and penomenology paper.
 double DelMh11,  DelMh12,  DelMh22,  DelMh23,  DelMh13,  DelMh33;
 
 //CP-Even Higgs 1 loop Corrections in the v1,v2,s basis
 double DelMh11prime,  DelMh12prime,  DelMh22prime,  DelMh23prime,  DelMh13prime,  DelMh33prime;
 
 double DoubleStop1derriv11, DoubleStop1derriv12 ,  DoubleStop1derriv22 ,DoubleStop1derriv23, DoubleStop1derriv13 ,  DoubleStop1derriv33  ;

 double DoubleStop2derriv11, DoubleStop2derriv12 ,  DoubleStop2derriv22 ,DoubleStop2derriv23, DoubleStop2derriv13 ,  DoubleStop2derriv33;
 
 double  mstop1derriv_1, mstop2derriv_1, mstop1derriv_2, mstop2derriv_2, mstop1derriv_s, mstop2derriv_s;

 if(WhatCorrections == 1 || WhatCorrections == 2){  
 double rthalf = Sqrt(Sqr(mQLsq- mURsq+0.125*Sqr(g2)*(Sqr(v1) - Sqr(v2)) - 0.125*Sqr(g1)*(Sqr(v1) - Sqr(v2))) + 4.0*mtop*mtop*Sqr(At - mylambda(3)*s*oneOrt2*v1/v2));
 double rtminushalf = 1/rthalf;
 double rtminus3half = rtminushalf*rtminushalf*rtminushalf;
 double Mstop11MINUSMstop22 = mQLsq- mURsq + 0.125*Sqr(g2)*(Sqr(v1) - Sqr(v2)) - 0.125*Sqr(g1)*(Sqr(v1) - Sqr(v2));
 double Mstop11MINUSMstop22repeat = mQLsq- mURsq+0.125*Sqr(g2)*(Sqr(v1) - Sqr(v2)) - 0.125*Sqr(g1)*(Sqr(v1) - Sqr(v2));
 double mtoprun = yt*v2/root2;

 double Delrt1 = 0.5*(Mstop11MINUSMstop22)*(Sqr(g2) - Sqr(g1))*v1 + 8.0*Sqr(mtop)*(At  - mylambda(3)*s*oneOrt2*v1/v2)*(-mylambda(3)*s*oneOrt2/(v2));
 double Delrt2 = 0.5*( Mstop11MINUSMstop22)*(Sqr(g1)- Sqr(g2))*v2+ 4.0*Sqr( At - mylambda(3)*s*oneOrt2*v1/v2)* Sqr(yt)*v2+ 8.0*Sqr(mtop)*(At - mylambda(3)*s*oneOrt2*v1/v2)*(mylambda(3)*s*v1*oneOrt2/(v2*v2) );

 // DH:: commented out repeated declaration here
 // double At = r.displayA_top();

 double Delrts = 8.0*Sqr(mtoprun)*(At - mu_eff/tb)*(-mylambda(3)/(root2*tb));
 
 DoubleStop2derriv11 = 0.5*(0.25*(Sqr(g2) + 3.0/(5.0)*Sqr(g1))  -  0.15*Sqr(g1p) ) + 0.25*rtminushalf*(2*(Mstop11MINUSMstop22)*(0.25*(Sqr(g2) -  Sqr(g1)))  + 0.125*Sqr(Sqr(g2) -  Sqr(g1) )*v1*v1 + 8.0* Sqr(mtoprun)*Sqr((mu_eff/v2)))-0.125*rtminus3half*(Sqr(Delrt1));		     
 
 DoubleStop1derriv11 = 0.5*(0.25*(Sqr(g2) + 3.0/(5.0)*Sqr(g1))  -  0.15*Sqr(g1p) ) - 0.25*rtminushalf*(2*(Mstop11MINUSMstop22)*(0.25*(Sqr(g2) -  Sqr(g1)))  + 0.125*Sqr(Sqr(g2) -  Sqr(g1) )*v1*v1 + 8.0* Sqr(mtoprun)*Sqr(mu_eff/v2))   + 0.125*rtminus3half*(Sqr(Delrt1));
 
 DoubleStop2derriv12 =  0.25*rtminushalf*(- 0.125*Sqr(( Sqr(g2) -  Sqr(g1) ))*v1*v2+ 8.0* Sqr(mtoprun)*(At - mu_eff/tb)*mu_eff/(Sqr(v2)) -  8.0* Sqr(mtoprun)*(Sqr(mu_eff))/(tb*Sqr(v2))  - 8.0*v2*Sqr(yt)*(At - mu_eff/tb)*(mu_eff/v2))
   - 0.125*rtminus3half*Delrt1*Delrt2;
  
 DoubleStop1derriv12 = - DoubleStop2derriv12;  
 DoubleStop2derriv13 = 2.0*rtminushalf*Sqr(mtoprun)*( (At - mu_eff/tb)*(-mylambda(3)/(root2*v2)) + mu_eff*mylambda(3)/(root2*tb*v2) )  +  rtminus3half*Sqr(mtoprun)*(At - mu_eff/tb)*mylambda(3)/(root2*tb)*Delrt1;
 
 DoubleStop1derriv13 = - DoubleStop2derriv13;
 DoubleStop2derriv23 =  0.25*rtminushalf*( 8.0*v2*Sqr(yt)*(At - mu_eff/tb)*(- mylambda(3)/(root2*tb))  +   8.0* Sqr(mtoprun)*(At - mu_eff/tb)*(mylambda(3)/(root2*tb*v2)) - 8.0* Sqr(mtoprun)*( mylambda(3)/(root2*tb))*(mu_eff/(v2*tb))   ) - 0.125*rtminus3half*Delrt2*Delrts;
   
 DoubleStop1derriv23 = - DoubleStop2derriv23;
 //cout << "  DoubleStop2derriv23 = " <<  DoubleStop2derriv23 << endl;
 DoubleStop2derriv22 = Sqr(yt) - 0.125*Sqr(gbar) - .05*Sqr(g1p) + 0.25*rtminushalf*(0.5*Mstop11MINUSMstop22*(Sqr(g1)- Sqr(g2)) + 0.125*Sqr(Sqr(g1)- Sqr(g2))*v2*v2 + 4.0*Sqr(At - mu_eff/tb)* Sqr(yt)  +8.0*Sqr(yt)*(At - mu_eff/tb)*mu_eff/tb + 4.0*Sqr(yt)*Sqr(mu_eff)/(tb*tb)   )  - 0.125*rtminus3half*Delrt2*Delrt2;

 DoubleStop1derriv22 = Sqr(yt) - 0.125*Sqr(gbar) - 0.05*Sqr(g1p) - 0.25*rtminushalf*(0.5*Mstop11MINUSMstop22*(Sqr(g1)- Sqr(g2)) + 0.125*Sqr(Sqr(g1)- Sqr(g2))*v2*v2 + 4.0*Sqr(At - mu_eff/tb)* Sqr(yt)  +8.0*Sqr(yt)*(At - mu_eff/tb)*mu_eff/tb + 4.0*Sqr(yt)*Sqr(mu_eff)/(tb*tb)      )  + 0.125*rtminus3half*Delrt2*Delrt2;

 DoubleStop2derriv33 = 0.125*Sqr(g1p)    + 0.25*rtminushalf*(4.0*Sqr(mtoprun)*Sqr(mylambda(3))/(tb*tb)) - 0.125*rtminus3half*Delrts*Delrts;
  DoubleStop1derriv33 = 0.125*Sqr(g1p)    - 0.25*rtminushalf*(4.0*Sqr(mtoprun)*Sqr(mylambda(3))/(tb*tb)) + 0.125*rtminus3half*Delrts*Delrts;  
  mstop2derriv_1 = 0.5*(0.6*0.25*Sqr(g1) + 0.25*Sqr(g2) - 6.0/(40.0) *Sqr(g1p)
			+0.5*( 2.0*(mQLsq- mURsq+0.125*Sqr(g2)*(Sqr(v1) - Sqr(v2)) - 0.125*Sqr(g1)*(Sqr(v1) - Sqr(v2)))*0.25*(Sqr(g2)- Sqr(g1)) + 8.0*Sqr(mtop)*(At  - mylambda(3)*s*oneOrt2*v1/v2)*(-mylambda(3)*s*oneOrt2/(v2*v1) )              ) /(   Sqrt(Sqr(mQLsq- mURsq+0.125*Sqr(g2)*(Sqr(v1) - Sqr(v2)) - 0.125*Sqr(g1)*(Sqr(v1) - Sqr(v2))) + 4.0*mtop*mtop*Sqr(At - mylambda(3)*s*oneOrt2*v1/v2))  )     )*v1;
  
 mstop1derriv_1 = 0.5*(0.6*0.25*Sqr(g1) + 0.25*Sqr(g2) - 6.0/(40.0) *Sqr(g1p)
		       -0.5*( 2.0*(mQLsq- mURsq+0.125*Sqr(g2)*(Sqr(v1) - Sqr(v2)) - 0.125*Sqr(g1)*(Sqr(v1) - Sqr(v2)))*0.25*(Sqr(g2)- Sqr(g1)) + 8.0*Sqr(mtop)*(At  - mylambda(3)*s*oneOrt2*v1/v2)*(-mylambda(3)*s*oneOrt2/(v2*v1) )              ) /(   Sqrt(Sqr(mQLsq- mURsq+0.125*Sqr(g2)*(Sqr(v1) - Sqr(v2)) - 0.125*Sqr(g1)*(Sqr(v1) - Sqr(v2))) + 4.0*mtop*mtop*Sqr(At - mylambda(3)*s*oneOrt2*v1/v2))  )     )*v1;
 
 mstop2derriv_2 = 0.5*(2.0*Sqr(yt)- 0.6*0.25*Sqr(g1) -  0.25*Sqr(g2) - 0.1*Sqr(g1p)
		       +0.5*( 2.0*(mQLsq- mURsq+0.125*Sqr(g2)*(Sqr(v1) - Sqr(v2)) - 0.125*Sqr(g1)*(Sqr(v1) - Sqr(v2)))*0.25*(Sqr(g1)- Sqr(g2)) + 4.0*Sqr( At - mylambda(3)*s*oneOrt2*v1/v2)* Sqr(yt) + 8.0*Sqr(mtop)*(At - mylambda(3)*s*oneOrt2*v1/v2)*(mylambda(3)*s*v1*oneOrt2/(v2*v2*v2) )              ) /(   Sqrt(Sqr(mQLsq- mURsq+0.125*Sqr(g2)*(Sqr(v1) - Sqr(v2)) - 0.125*Sqr(g1)*(Sqr(v1) - Sqr(v2))) + 4.0*mtop*mtop*Sqr(At - mylambda(3)*s*oneOrt2*v1/v2))  )     )*v2;
 
 mstop1derriv_2  = 0.5*(2.0*Sqr(yt)- 0.6*0.25*Sqr(g1) -  0.25*Sqr(g2) - 0.1*Sqr(g1p)
			-0.5*( 2.0*(mQLsq- mURsq+0.125*Sqr(g2)*(Sqr(v1) - Sqr(v2)) - 0.125*Sqr(g1)*(Sqr(v1) - Sqr(v2)))*0.25*(Sqr(g1)- Sqr(g2)) + 4.0*Sqr( At - mylambda(3)*s*oneOrt2*v1/v2)* Sqr(yt) + 8.0*Sqr(mtop)*(At - mylambda(3)*s*oneOrt2*v1/v2)*(mylambda(3)*s*v1*oneOrt2/(v2*v2*v2) )              ) /(   Sqrt(Sqr(mQLsq- mURsq+0.125*Sqr(g2)*(Sqr(v1) - Sqr(v2)) - 0.125*Sqr(g1)*(Sqr(v1) - Sqr(v2))) + 4.0*mtop*mtop*Sqr(At - mylambda(3)*s*oneOrt2*v1/v2))  )     )*v2;
 
 mstop2derriv_s = 0.5*(0.25*Sqr(g1p)
		       +0.5*( 8.0*Sqr(mtop)*(At - mylambda(3)*s*oneOrt2*v1/v2)*(-mylambda(3)*v1*oneOrt2/(v2*s) ) ) /(   Sqrt(Sqr(mQLsq- mURsq+0.125*Sqr(g2)*(Sqr(v1) - Sqr(v2)) - 0.125*Sqr(g1)*(Sqr(v1) - Sqr(v2))) + 4.0*mtop*mtop*Sqr(At - mylambda(3)*s*oneOrt2*v1/v2))  ))*s;
 
 mstop1derriv_s = 0.5*(0.25*Sqr(g1p)
		       -0.5*( 8.0*Sqr(mtop)*(At - mylambda(3)*s*oneOrt2*v1/v2)*(-mylambda(3)*v1*oneOrt2/(v2*s) ) ) /(   Sqrt(Sqr(mQLsq- mURsq+0.125*Sqr(g2)*(Sqr(v1) - Sqr(v2)) - 0.125*Sqr(g1)*(Sqr(v1) - Sqr(v2))) + 4.0*mtop*mtop*Sqr(At - mylambda(3)*s*oneOrt2*v1/v2))  ))*s;
 }
 else{     
   double rthalf = Sqrt(Sqr(mQLsq- mURsq) + 4.0*mtop*mtop*Sqr(At - mylambda(3)*s*oneOrt2*v1/v2));
   double rtminushalf = 1/rthalf;
   double rtminus3half = rtminushalf*rtminushalf*rtminushalf;
   double Mstop11MINUSMstop22 = mQLsq- mURsq;
   double mtoprun = yt*v2/root2;
   double Delrt1 = 8.0*Sqr(mtop)*(At  - mylambda(3)*s*oneOrt2*v1/v2)*(-mylambda(3)*s*oneOrt2/(v2));
   double Delrt2 =  4.0*Sqr( At - mylambda(3)*s*oneOrt2*v1/v2)* Sqr(yt)*v2+ 8.0*Sqr(mtop)*(At - mylambda(3)*s*oneOrt2*v1/v2)*(mylambda(3)*s*v1*oneOrt2/(v2*v2) );

   // DH:: commented out repeated declaration here
   //double At = r.displayA_top();
   double Delrts = 8.0*Sqr(mtoprun)*(At - mu_eff/tb)*(-mylambda(3)/(root2*tb));
   
   DoubleStop2derriv11 =  0.25*rtminushalf*( 8.0* Sqr(mtoprun)*Sqr((mu_eff/v2))) -0.125*rtminus3half*(Sqr(Delrt1));		     
   
   DoubleStop1derriv11 = - 0.25*rtminushalf*(  8.0* Sqr(mtoprun)*Sqr(mu_eff/v2))   + 0.125*rtminus3half*(Sqr(Delrt1));
   
   DoubleStop2derriv12 =  0.25*rtminushalf*( 8.0* Sqr(mtoprun)*(At - mu_eff/tb)*mu_eff/(Sqr(v2)) -  8.0* Sqr(mtoprun)*(Sqr(mu_eff))/(tb*Sqr(v2))  - 8.0*v2*Sqr(yt)*(At - mu_eff/tb)*(mu_eff/v2))
     - 0.125*rtminus3half*Delrt1*Delrt2;
 DoubleStop1derriv12 = - DoubleStop2derriv12;
 
 DoubleStop2derriv13 = 2.0*rtminushalf*Sqr(mtoprun)*( (At - mu_eff/tb)*(-mylambda(3)/(root2*v2)) + mu_eff*mylambda(3)/(root2*tb*v2) )  +  rtminus3half*Sqr(mtoprun)*(At - mu_eff/tb)*mylambda(3)/(root2*tb)*Delrt1;
 
 DoubleStop1derriv13 = - DoubleStop2derriv13;
 DoubleStop2derriv23 =  0.25*rtminushalf*( 8.0*v2*Sqr(yt)*(At - mu_eff/tb)*(- mylambda(3)/(root2*tb))  +   8.0* Sqr(mtoprun)*(At - mu_eff/tb)*(mylambda(3)/(root2*tb*v2)) - 8.0* Sqr(mtoprun)*( mylambda(3)/(root2*tb))*(mu_eff/(v2*tb))   ) - 0.125*rtminus3half*Delrt2*Delrts;
 
 DoubleStop1derriv23 = - DoubleStop2derriv23;
 //cout << "  DoubleStop2derriv23 = " <<  DoubleStop2derriv23 << endl;
 DoubleStop2derriv22 = Sqr(yt) + 0.25*rtminushalf*(  4.0*Sqr(At - mu_eff/tb)* Sqr(yt)  +8.0*Sqr(yt)*(At - mu_eff/tb)*mu_eff/tb + 4.0*Sqr(yt)*Sqr(mu_eff)/(tb*tb)   )  - 0.125*rtminus3half*Delrt2*Delrt2;
 
 
 DoubleStop1derriv22 = Sqr(yt) - 0.25*rtminushalf*( 4.0*Sqr(At - mu_eff/tb)* Sqr(yt)  +8.0*Sqr(yt)*(At - mu_eff/tb)*mu_eff/tb + 4.0*Sqr(yt)*Sqr(mu_eff)/(tb*tb)      ) + 0.125*rtminus3half*Delrt2*Delrt2;
 
 DoubleStop2derriv33 =  + 0.25*rtminushalf*(4.0*Sqr(mtoprun)*Sqr(mylambda(3))/(tb*tb)) - 0.125*rtminus3half*Delrts*Delrts;
 
 DoubleStop1derriv33 =  - 0.25*rtminushalf*(4.0*Sqr(mtoprun)*Sqr(mylambda(3))/(tb*tb)) + 0.125*rtminus3half*Delrts*Delrts;
 
 mstop2derriv_1 = 0.5*(0.5*(  8.0*Sqr(mtop)*(At  - mylambda(3)*s*oneOrt2*v1/v2)*(-mylambda(3)*s*oneOrt2/(v2*v1) )              ) /(   Sqrt(Sqr(mQLsq- mURsq) + 4.0*mtop*mtop*Sqr(At - mylambda(3)*s*oneOrt2*v1/v2))  )     )*v1;

 mstop1derriv_1 = 0.5*(-0.5*( 8.0*Sqr(mtop)*(At  - mylambda(3)*s*oneOrt2*v1/v2)*(-mylambda(3)*s*oneOrt2/(v2*v1) )              ) /(   Sqrt(Sqr(mQLsq- mURsq) + 4.0*mtop*mtop*Sqr(At - mylambda(3)*s*oneOrt2*v1/v2))  )     )*v1;

 mstop2derriv_2 = 0.5*(2.0*Sqr(yt)
		       +0.5*( 4.0*Sqr( At - mylambda(3)*s*oneOrt2*v1/v2)* Sqr(yt) + 8.0*Sqr(mtop)*(At - mylambda(3)*s*oneOrt2*v1/v2)*(mylambda(3)*s*v1*oneOrt2/(v2*v2*v2) )              ) /(   Sqrt(Sqr(mQLsq- mURsq) + 4.0*mtop*mtop*Sqr(At - mylambda(3)*s*oneOrt2*v1/v2))  )     )*v2;
 
  //cout << "mstop1derriv_2 = " << mstop1derriv_2 << endl;
 mstop1derriv_2  = 0.5*(2.0*Sqr(yt)
			-0.5*( 4.0*Sqr( At - mylambda(3)*s*oneOrt2*v1/v2)* Sqr(yt) + 8.0*Sqr(mtop)*(At - mylambda(3)*s*oneOrt2*v1/v2)*(mylambda(3)*s*v1*oneOrt2/(v2*v2*v2) )              ) /(   Sqrt(Sqr(mQLsq- mURsq) + 4.0*mtop*mtop*Sqr(At - mylambda(3)*s*oneOrt2*v1/v2))  )     )*v2;
 
 mstop2derriv_s = 0.5*(0.5*( 8.0*Sqr(mtop)*(At - mylambda(3)*s*oneOrt2*v1/v2)*(-mylambda(3)*v1*oneOrt2/(v2*s) ) ) /(   Sqrt(Sqr(mQLsq- mURsq) + 4.0*mtop*mtop*Sqr(At - mylambda(3)*s*oneOrt2*v1/v2))  ))*s;
 
 mstop1derriv_s = 0.5*( -0.5*( 8.0*Sqr(mtop)*(At - mylambda(3)*s*oneOrt2*v1/v2)*(-mylambda(3)*v1*oneOrt2/(v2*s) ) ) /(   Sqrt(Sqr(mQLsq- mURsq) + 4.0*mtop*mtop*Sqr(At - mylambda(3)*s*oneOrt2*v1/v2))  ))*s;
 }
 
 
 DelMh11prime = 3.0/(32.0*Sqr(PI))*(2.0*(2.0*Log(mstop(1)/r.get_scale())*Sqr( mstop1derriv_1 )  + a0Peter(mstopsq(1),r.get_scale())* DoubleStop1derriv11)   +2.0*(2.0*Log(mstop(2)/r.get_scale())*Sqr( mstop2derriv_1 )  + a0Peter(mstopsq(2),r.get_scale())* DoubleStop2derriv11));
 
 DelMh12prime = 3.0/(32.0*Sqr(PI))*(2.0*(2.0*Log(mstop(1)/r.get_scale())*(mstop1derriv_1*mstop1derriv_2)+  a0Peter(mstopsq(1),r.get_scale())* DoubleStop1derriv12) + 2.0*( 2.0*Log(mstop(2)/r.get_scale())*mstop2derriv_1*mstop2derriv_2 +  a0Peter(mstopsq(2),r.get_scale())* DoubleStop2derriv12));
 
 DelMh13prime = 3.0/(32.0*Sqr(PI))*(2.0*(2.0*Log(mstop(1)/r.get_scale())*(mstop1derriv_1*mstop1derriv_s)+  a0Peter(mstopsq(1),r.get_scale())* DoubleStop1derriv13) + 2.0*( 2.0*Log(mstop(2)/r.get_scale())*mstop2derriv_1*mstop2derriv_s +  a0Peter(mstopsq(2),r.get_scale())* DoubleStop2derriv13));
 
 DelMh22prime = 3.0/(32.0*Sqr(PI))*(2.0*(2.0*Log(mstop(1)/r.get_scale())*Sqr( mstop1derriv_2 )  + a0Peter(mstopsq(1),r.get_scale())* DoubleStop1derriv22)   +2.0*(2.0*Log(mstop(2)/r.get_scale())*Sqr( mstop2derriv_2 )  + a0Peter(mstopsq(2),r.get_scale())* DoubleStop2derriv22) - 4.0*(2.0*Log(mtop/r.get_scale())*Sqr(Sqr(yt))*Sqr(v2) +a0Peter(Sqr(mtop),r.get_scale())*Sqr(yt)));
 
 DelMh23prime = 3.0/(32.0*Sqr(PI))*(2.0*(2.0*Log(mstop(1)/r.get_scale())*(mstop1derriv_2*mstop1derriv_s) +  a0Peter(mstopsq(1), r.get_scale())* DoubleStop1derriv23) + 2.0*( 2.0*Log(mstop(2)/r.get_scale())*mstop2derriv_2*mstop2derriv_s +  a0Peter(mstopsq(2),r.get_scale())* DoubleStop2derriv23));
 
 //cout << " DelMh23prime = " <<  DelMh23prime  << endl;
 DelMh33prime = 3.0/(32.0*Sqr(PI))*(2.0*(2.0*Log(mstop(1)/r.get_scale())*Sqr( mstop1derriv_s )  + a0Peter(mstopsq(1),r.get_scale())* DoubleStop1derriv33)   +2.0*(2.0*Log(mstop(2)/r.get_scale())*Sqr( mstop2derriv_s )  + a0Peter(mstopsq(2),r.get_scale())* DoubleStop2derriv33));
 //cout << " DelMh33prime = " <<  DelMh33prime  << endl;

 // cout << "Del11prime = " << DelMh11prime << endl;
 // cout << "Del12prime = " << DelMh12prime << endl;
 // cout << "Del13prime = " << DelMh13prime << endl;
 // cout << "Del22prime = " << DelMh22prime << endl;
 // cout << "Del23prime = " << DelMh23prime << endl;
 // cout << "Del33prime = " << DelMh33prime << endl;

//REMEMBER TO INCLUDE TADPOLE CORRECTION TERMS  BECAUSE I HAVE USED TREELEVEL EWSB CONDITIONS TO PRODUCE THE TRELEVEL MASS MATRIX

	

if(WhatCorrections ==1){
  DelMh11prime =  DelMh11prime  + doCalcTadpoleESSMH1(r,s,tb);
  DelMh22prime =  DelMh22prime  + doCalcTadpoleESSMH2( r, s ,tb  );
  DelMh33prime =  DelMh33prime  + doCalcTadpolesESSMS( r, s ,tb  );
 }
 if( WhatCorrections ==2){
   DelMh11prime =  DelMh11prime  + doCalcTadpoleESSMH1_atMt(r,s,tb);
   DelMh22prime =  DelMh22prime  + doCalcTadpoleESSMH2_atMt( r, s ,tb  );
   DelMh33prime =  DelMh33prime  + doCalcTadpolesESSMS_atMt( r, s ,tb  );  
 }
 
 if( WhatCorrections ==3){
   DelMh11prime =  DelMh11prime  + doCalcTadpoleESSMH1_Roman_atQ(r,s,tb);
   DelMh22prime =  DelMh22prime  + doCalcTadpoleESSMH2_Roman_atQ( r, s ,tb  );
   DelMh33prime =  DelMh33prime  + doCalcTadpolesESSMS_Roman_atQ( r, s ,tb  );
}
 
 if( WhatCorrections ==4){
   DelMh11prime =  DelMh11prime  + doCalcTadpoleESSMH1_Roman(r,s,tb);
   DelMh22prime =  DelMh22prime  + doCalcTadpoleESSMH2_Roman( r, s ,tb  );
   DelMh33prime =  DelMh33prime  + doCalcTadpolesESSMS_Roman( r, s ,tb  );
}


 // cout << "After adding tadpoles: " << endl;
 // cout << "Del11prime = " << DelMh11prime << endl;
 // cout << "Del12prime = " << DelMh12prime << endl;
 // cout << "Del13prime = " << DelMh13prime << endl;
 // cout << "Del22prime = " << DelMh22prime << endl;
 // cout << "Del23prime = " << DelMh23prime << endl;
 // cout << "Del33prime = " << DelMh33prime << endl;


 //Rotate by beta
 DelMh11  = Sqr(Cos(ArcTan(tb))) * DelMh11prime + 2*tb/(1+tb*tb)*DelMh12prime + Sqr(Sin(ArcTan(tb)))* DelMh22prime;
 DelMh22  = Sqr(Sin(ArcTan(tb))) * DelMh11prime - 2*tb/(1+tb*tb)*DelMh12prime + Sqr(Cos(ArcTan(tb)))* DelMh22prime;
 DelMh33 = DelMh33prime;
 DelMh12  = (Sqr(Cos(ArcTan(tb))) -Sqr(Sin(ArcTan(tb)))) * DelMh12prime  + tb/(1+tb*tb)*( DelMh22prime - DelMh11prime);
 DelMh13 = Cos(ArcTan(tb))* DelMh13prime + Sin(ArcTan(tb))* DelMh23prime; 
 DelMh23 = Cos(ArcTan(tb))* DelMh23prime - Sin(ArcTan(tb))* DelMh13prime; 

 //entries
 double HE11 = 2.0*Sqr(mylambda(3))*v1*v2*tb/(1 + tb*tb) + 0.25*Sqr(gbar)*(Sqr(v1) - Sqr(v2))*(1-tb*tb)/(1 + tb*tb) + Sqr(g1p)*Sqr(v1*Q1eff*Cos(ArcTan(tb)) + v2*Sin(ArcTan(tb))*Q2eff) + DelMh11;
 
 const bool Includeleadtwoloop = true; // True = use two loop. Should use two loop in final results?
 if(Includeleadtwoloop)
   {
     //Remove one loop because we will replace it with equivelent terms obtained from rge evolution
     HE11 = HE11 - DelMh11;
     //if(speak) cout << "Tree (hopefully!) HE11 = " << HE11 << endl;
     double l = Log(mstop(1)*mstop(2)/(Sqr(mtop)));
     //if(speak) cout << " l = " << l << endl;  
     double Xt = At - mu_eff/tb;
     //if(speak)cout << "Xt = " << Xt << endl;
     double Ut = 2.0*Sqr(Xt)/(mstop(1)*mstop(2))*(1 - Sqr(Xt)/(12.0*mstop(1)*mstop(2)));
     if(Bugspeak)
       { 
	 cerr << " Ut = " << Ut << endl;
	 cerr << "Should be like 1lp HE11 = " <<  HE11*(1 - 3*Sqr(r.get_Yu(2, 2))*l/(8.0*Sqr(PI)))+ 3.0 * Sqr(Sqr(r.get_Yu( 2, 2)) )*Sqr(vev)*Sqr(Sqr(Sin(ArcTan(tb))))/(8.0*Sqr(PI))*(0.5*Ut + l) << endl;
       }
  
     if(Bugspeak) cerr << "Should be likeDelMh11 = " <<  HE11*( - 3*Sqr(r.get_Yu( 2, 2))*l)/(8.0*Sqr(PI))+ 3.0 * Sqr(Sqr(r.get_Yu( 2, 2)) )*Sqr(vev)*Sqr(Sqr(Sin(ArcTan(tb))))/(8.0*Sqr(PI))*(0.5*Ut + l) << endl;
     
     HE11 = HE11*(1 - 3*Sqr(r.get_Yu( 2, 2))*l/(8.0*Sqr(PI)))+ 3.0 * Sqr(Sqr(r.get_Yu( 2, 2)) )*Sqr(vev)*Sqr(Sqr(Sin(ArcTan(tb))))/(8.0*Sqr(PI))*(0.5*Ut + l + 1/(16*Sqr(PI))*(1.5*Sqr(r.get_Yu( 2, 2)) - 8.0*Sqr(r.get_g3()))*(Ut + l)*l);
     
     //if(speak) cout << "After 2lp HE11 = " << HE11 << endl; 
     
   }
 double HE12 = (0.25*Sqr(mylambda(3)) - 0.125*Sqr(gbar))*Sqr(vev)* 4.0*tb*(1 - Sqr(tb))/(Sqr(1.0+tb*tb)) +   Sqr(g1p)*( Sqr(v1)*Q1eff + Sqr(v2)*Q2eff)*(Q2eff - Q1eff)*tb/(1+tb*tb)  + DelMh12;
 double HE13 = -mylambda(3)*myAlambda(3)*v2*Sqrt(2)*Cos(ArcTan(tb)) + Sqr( mylambda(3))*s*vev +  Sqr(g1p)*( Sqr(v1)*Q1eff + Sqr(v2)*Q2eff)*QSeff*s/(vev) + DelMh13;
 double HE22 = mylambda(3)*myAlambda(3)*Sqrt(2.0)*s*(1+tb*tb)/(2.0*tb) - (0.5*Sqr(mylambda(3)) - 0.25*Sqr(gbar))*Sqr(vev)* 4.0*tb*tb/(Sqr(1.0+tb*tb)) +Sqr(g1p)*Sqr(vev)*tb*tb/(Sqr(1.0+tb*tb))*Sqr(Q2eff - Q1eff) + DelMh22 ;
 double HE23 = - mylambda(3)*myAlambda(3)/(Sqrt(2))*vev*(1 - tb*tb)/(tb*tb +1.0) + Sqr(g1p)*(Q2eff - Q1eff)*QSeff*v2*s*Cos(ArcTan(tb)) + DelMh23;
 double HE33 = mylambda(3)*myAlambda(3)/(s*Sqrt(2))*v1*v2+ Sqr(g1p)*QSeff*QSeff*s*s + DelMh33;
 
 //Now do this with softsusy routines to replace NR routines
 DoubleVector mhphysq(3);
 DoubleMatrix MH_even(3,3);
 DoubleMatrix mixMH(3,3);
 MH_even(1,1) = HE11; 
 MH_even(1,2) = HE12; 
 MH_even(1,3) = HE13; 
 MH_even(2,1) = HE12;
 MH_even(2,2) = HE22;
 MH_even(2,3) = HE23;
 MH_even(3,1) = HE13;
 MH_even(3,2) = HE23;
 MH_even(3,3) = HE33;
 
 double tan_psi = vev/(2*s)*2*tb/(Sqr(tb) + 1); 
 double cos_psi = 1/Sqrt(1 + Sqr(tan_psi));
 double sin_psi = tan_psi/Sqrt(1 + Sqr(tan_psi));
  
 double Twolplight_higgs;
 if((mstopsq(2) < 0)|| (mstopsq(1) < 0))
   {
     if(speak)
       { 
	 cerr << "tachyonics stop masses so i'm not going to bother doing higgs mass corrections.  Point already ruled out. " << endl;
	 
	 // Dylan:: I have modified this so that in the event of tachyonic stops, the returned Higgs masses are zero
	 // as a flag.
	 mhout.set(1, 0.0);
	 mhout.set(2, 0.0);
	 mhout.set(3, 0.0);
	 ExpValid = TADPOLESPROBLEM;
       }  
   }
 else
   {

     if (MH_even.diagonaliseSym(mixMH, mhphysq) >  TOLERANCE * 1.0e-3) 
       { 
	 cerr << "Warning:  accuracy bad in CP-even Higgs diagonalisation" << endl;
	 // Flag to indicate poor accuracy
	 sing = 1;
	 ExpValid = HIGGSPROBLEM;
       }
     DoubleVector mhphy(3);
     
     if (mhphysq(1) < 0. || mhphysq(2) < 0. || mhphysq(3) < 0.) 
       {
	 
	 ExpValid = POLEHIGGSTACHYON;

	 higgsTachyon = true;	 

	 if (mhphysq(1) < 0.0)
	   {
	     mhphy.set(1, ZeroSqrt(mhphysq(1)));
	   }
	 else
	   {
	     mhphy.set(1, Sqrt(mhphysq(1)));
	   }
	 
	 if (mhphysq(2) < 0.0)
	   {
	     mhphy.set(2, ZeroSqrt(mhphysq(2)));
	   }
	 else
	   {
	     mhphy.set(2, Sqrt(mhphysq(2)));
	   }
	 
	 if (mhphysq(3) < 0.0)
	   {
	     mhphy.set(3, ZeroSqrt(mhphysq(3)));
	   }
	 else
	   {
	     mhphy.set(3, Sqrt(mhphysq(3)));
	   }
	 
       }
     else
       {
	 mhphy.set(1, Sqrt(mhphysq(1)));
	 mhphy.set(2, Sqrt(mhphysq(2)));
	 mhphy.set(3, Sqrt(mhphysq(3)));
	 
       }
     

     mhout = mhphy.sort();   
 
     mhmix = mixMH; // NOTE ADDITIONAL OUTPUT TO SAVE MIXING AS WELL
     msq = MH_even;
     
     if(speak)cerr <<  "2lp light_higgs = " <<mhphy(1) << endl;
     if(speak)cerr <<  "higgs = " <<mhphy << endl;
     
   }
 
 
 // Have been using 120 GeV and 130 GeV as limits. Note additional check to make sure 
 // not already flagged as tachyonic
 if (mhout.display(1) >= bounds.display(1) && mhout.display(1) <= bounds.display(2) && ExpValid != POLEHIGGSTACHYON)
   {
     ExpValid = 0;
   }
 
 return higgsTachyon;
 
}


// Calculate m_A^2 at tree level.
// Inputs:
//    SoftParsEssm const & essmSusy = the ESSM model to calculate m_A^2 for
//    double s = the value of the singlet VEV to use
//    double tb = the value of tan(beta) to use
double mAsq_TreeLevel(genericE6SSM_soft_parameters const & essmSusy, double s, double tb)
{

  double lambda = essmSusy.get_Lambdax();
  double Tlambda;
  double Alambda;
  
  Tlambda = essmSusy.get_TLambdax();
  
  if (Abs(Tlambda) < EPSTOL) 
    {
      Alambda = 0.0;
    }
  else if (Abs(lambda) < 1.0e-100)
    {
      ostringstream ii;
      ii << "WARNING: trying to calculate A_lambda(3) where lambda(3) coupling is " <<
	Abs(lambda) << endl;
      throw ii.str();
    }
  else
    {
      Alambda = Tlambda/lambda;
    }

  double v1 = essmSusy.get_vd();
  double v2 = essmSusy.get_vu();

  double v = Sqrt(v1*v1+v2*v2);

  double s2b = 2.0*tb/(1.0+tb*tb);

  double mAsq = Sqrt(2.0)*lambda*Alambda*s/s2b + Sqrt(2.0)*lambda*Alambda*v*v*s2b/(4.0*s);

  return mAsq;
}

// Calculate m_A^2 at one loop order.
// Inputs:
//    SoftParsEssm essmSusy = the ESSM model to calculate m_A^2 for
//    double s = the value of the singlet VEV to use
//    double tb = the value of tan(beta) to use
double mAsq_OneLoop(genericE6SSM_soft_parameters essmSusy, double s, double tb, int & problem)
{
  problem = 0;

  double v1 = essmSusy.get_vd();
  double v2 = essmSusy.get_vu();

  double v = Sqrt(v1*v1+v2*v2);

  double Q = essmSusy.get_scale(); // Renormalisation scale

  // Get stop mass (note option to fix m_t(M_t) = 165 GeV)
  double yt = essmSusy.get_Yu(2, 2);
  
  double lambda = essmSusy.get_Lambdax();
  double Tlambda;
  double Alambda;
  
  Tlambda = essmSusy.get_TLambdax();
  
  if (Abs(Tlambda) < EPSTOL) 
    {
      Alambda = 0.0;
    }
  else if (Abs(lambda) < 1.0e-100)
    {
      ostringstream ii;
      ii << "WARNING: trying to calculate A_lambda(3) where lambda(3) coupling is " <<
	Abs(lambda) << endl;
      throw ii.str();
    }
  else
    {
      Alambda = Tlambda/lambda;
    }
  
  double At;
  double TYt = essmSusy.get_TYu(2,2);
  
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

  double mtop = yt*v2/(Sqrt(2.0));
  double oneOrt2 = 1.0/Sqrt(2.0);

  if(USEMTOFMT){
    mtop = 165;
    yt = mtop/(v2*oneOrt2);
    essmSusy.set_Yu(2,2, yt);    
  }

  // Get stop masses
  DoubleVector mstop(2), mstopsq(2), mD1sq(3), mD2sq(3);

  physical_ESSM(essmSusy, mstop, mstopsq, mD1sq, mD2sq, s, tb);

  if (mstopsq(1) < 0.0 || mstopsq(2) < 0.0)
    {
      problem = TADPOLESPROBLEM;
    }

  double mstop1 = mstop.display(1); // mstop1 is heavier
  double mstop2 = mstop.display(2);

  double fVal;
  fVal = f(mstop2*mstop2, mstop1*mstop1, Q);

  double sb = tb/Sqrt(1.0+tb*tb);
  double cb = 1.0/Sqrt(1.0+tb*tb);
  double s2b = 2.0*tb/(1.0+tb*tb);

  double delta_A = (Sqrt(2.0)*lambda*s/s2b+Sqrt(2.0)*lambda*v*v*s2b/(4.0*s))*
    (3.0/((4*PI)*(4*PI)))*yt*yt*At*fVal;

  double eig1 = mAsq_TreeLevel(essmSusy, s, tb)+delta_A;

  return eig1;

}

// A helper function used in calculating the mass m_A^2 at one loop order.
// Inputs:
//    double m1Sq = first squared mass m_1^2 to use
//    double m2Sq = second squared mass m_2^2 to use
//    double q = renormalisation scale
double f(double m1Sq, double m2Sq, double Q)
{
  double fVal = (1.0/(m1Sq-m2Sq))*(m1Sq*Log(m1Sq/(Q*Q))-m2Sq*Log(m2Sq/(Q*Q))-m1Sq+m2Sq);

  return fVal;
}

// UP TO HERE

// The following are helper functions that are useful for constructing
// the the leading one-loop tuning measures.

// doCalcdmstop1sqdvi calculates the derivative of m_stop_1^2 with respect to
// v_i for the given ESSM model, i = 1, 2, 3 (i=3 is s). 
// NB: m_stop_1^2 is taken to be the lighter stop in this method 
// -i.e. minus sign applies.
// Inputs:
//     genericE6SSM_soft_parameters const & essmSusy = object to calculate derivatives on
//     double s = the value of the VEV s to use
//     double tb = the value of tan(beta) to use
double doCalcdmstop1sqdv1(genericE6SSM_soft_parameters essmSusy, double s, double tb)
{
  double v1 = essmSusy.get_vd();
  double v2 = essmSusy.get_vu(); 
  double v = Sqrt(v1*v1+v2*v2);

  double g1 = essmSusy.get_g1();
  double g2 = essmSusy.get_g2();
  double gd1 = essmSusy.get_gN();

  double yt = essmSusy.get_Yu( 2, 2);
  double mtop = yt*v2/Sqrt(2.0);

  // Option to fix m_t(M_t) = 165 GeV
  
  if(USEMTOFMT){
    mtop = 165;
    yt = Sqrt(2.0)*mtop/v2;
    essmSusy.set_Yu( 2, 2, yt);
  }

  
  double lambda = essmSusy.get_Lambdax();
  
  double At;
  double TYt = essmSusy.get_TYu(2,2);
  
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

  double mQlsq = essmSusy.get_mq2(2,2);
  double mUrsq =  essmSusy.get_mu2(2,2);

  double MQQsq = mQlsq-mUrsq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2);
  double Xt = At - lambda*s/(Sqrt(2.0)*tb);
  double rt = MQQsq*MQQsq+4.0*mtop*mtop*Xt*Xt;

  double deriv = 0.25*g2*g2+0.6*0.25*g1*g1-0.15*gd1*gd1;
  deriv = deriv-(1.0/Sqrt(rt))*(MQQsq*0.25*(g2*g2-g1*g1)
				+4.0*mtop*mtop*Xt*(-lambda*s/(Sqrt(2.0)*v1*v2)));
  deriv = deriv*(v1/2.0);

  return deriv;

}

double doCalcdmstop1sqdv2(genericE6SSM_soft_parameters essmSusy, double s, double tb)
{
  double v1 = essmSusy.get_vd(); 
  double v2 = essmSusy.get_vu(); 
  double v = Sqrt(v1*v1+v2*v2);

  double g1 = essmSusy.get_g1();
  double g2 = essmSusy.get_g2();
  double gd1 = essmSusy.get_gN();

  double yt = essmSusy.get_Yu( 2, 2);
  double mtop = yt*v2/Sqrt(2.0);
  // Option to fix m_t(M_t) = 165 GeV
  
  if(USEMTOFMT){
    mtop = 165;
    yt = Sqrt(2.0)*mtop/v2;
    essmSusy.set_Yu( 2, 2, yt);
  }

  
  double lambda = essmSusy.get_Lambdax();
  
  double At;
  double TYt = essmSusy.get_TYu(2,2);
  
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

  double mQlsq = essmSusy.get_mq2(2,2);
  double mUrsq =  essmSusy.get_mu2(2,2);

  double MQQsq = mQlsq-mUrsq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2);
  double Xt = At - lambda*s/(Sqrt(2.0)*tb);
  double rt = MQQsq*MQQsq+4.0*mtop*mtop*Xt*Xt;

  double deriv = 2.0*yt*yt-0.25*g2*g2-0.6*0.25*g1*g1-0.1*gd1*gd1;
  deriv = deriv-(1.0/Sqrt(rt))*(0.25*MQQsq*(g1*g1-g2*g2)+2.0*yt*yt*Xt*Xt
				+4.0*mtop*mtop*Xt*(lambda*s*v1/(Sqrt(2.0)*v2*v2*v2)));
  deriv = deriv*(v2/2.0);

  return deriv;
}

double doCalcdmstop1sqdv3(genericE6SSM_soft_parameters essmSusy, double s, double tb)
{
  double v1 = essmSusy.get_vd();
  double v2 = essmSusy.get_vu(); 
  double v = Sqrt(v1*v1+v2*v2);

  double g1 = essmSusy.get_g1();
  double g2 = essmSusy.get_g2();
  double gd1 = essmSusy.get_gN();

  double yt = essmSusy.get_Yu( 2, 2);
  double mtop = yt*v2/Sqrt(2.0);
  // Option to fix m_t(M_t) = 165 GeV
  
  if(USEMTOFMT){
    mtop = 165;
    yt = Sqrt(2.0)*mtop/v2;
    essmSusy.set_Yu(2,2, yt);
  }

  
  double lambda = essmSusy.get_Lambdax();
 
  double At;
  double TYt = essmSusy.get_TYu(2,2);
  
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

  double mQlsq = essmSusy.get_mq2(2,2);
  double mUrsq =  essmSusy.get_mu2(2,2);

  double MQQsq = mQlsq-mUrsq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2);
  double Xt = At - lambda*s/(Sqrt(2.0)*tb);
  double rt = MQQsq*MQQsq+4.0*mtop*mtop*Xt*Xt;

  double deriv = 0.25*gd1*gd1-(4.0*mtop*mtop/Sqrt(rt))*Xt*(-lambda*v1/(Sqrt(2.0)*v2*s));
  deriv = 0.5*s*deriv;

  return deriv;
}

// doCalcdmstop2sqdvi is as above but for m_stop_2^2 (+ sign applies).
double doCalcdmstop2sqdv1(genericE6SSM_soft_parameters essmSusy, double s, double tb)
{
  double v1 = essmSusy.get_vd(); 
  double v2 = essmSusy.get_vu(); 
  double v = Sqrt(v1*v1+v2*v2);

  double g1 = essmSusy.get_g1();
  double g2 = essmSusy.get_g2();
  double gd1 = essmSusy.get_gN();

  double yt = essmSusy.get_Yu( 2, 2);
  double mtop = yt*v2/Sqrt(2.0);
  // Option to fix m_t(M_t) = 165 GeV
  
  if(USEMTOFMT){
    mtop = 165;
    yt = Sqrt(2.0)*mtop/v2;
    essmSusy.set_Yu(2,2, yt);
  }

  
  double lambda = essmSusy.get_Lambdax();

  double At;
  double TYt = essmSusy.get_TYu(2,2);
  
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

  double mQlsq = essmSusy.get_mq2(2,2);
  double mUrsq =  essmSusy.get_mu2(2,2);

  double MQQsq = mQlsq-mUrsq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2);
  double Xt = At - lambda*s/(Sqrt(2.0)*tb);
  double rt = MQQsq*MQQsq+4.0*mtop*mtop*Xt*Xt;

  double deriv = 0.25*g2*g2+0.6*0.25*g1*g1-0.15*gd1*gd1;
  deriv = deriv+(1.0/Sqrt(rt))*(MQQsq*0.25*(g2*g2-g1*g1)
				+4.0*mtop*mtop*Xt*(-lambda*s/(Sqrt(2.0)*v1*v2)));
  deriv = deriv*(v1/2.0);

  return deriv;
}

double doCalcdmstop2sqdv2(genericE6SSM_soft_parameters essmSusy, double s, double tb)
{
  double v1 = essmSusy.get_vd();
  double v2 = essmSusy.get_vu(); 
  double v = Sqrt(v1*v1+v2*v2);

  double g1 = essmSusy.get_g1();
  double g2 = essmSusy.get_g2();
  double gd1 = essmSusy.get_gN();

  double yt = essmSusy.get_Yu( 2, 2);
  double mtop = yt*v2/Sqrt(2.0);
  // Option to fix m_t(M_t) = 165 GeV
  
  if(USEMTOFMT){
    mtop = 165;
    yt = Sqrt(2.0)*mtop/v2;
    essmSusy.set_Yu(2,2, yt);
  }

  
  double lambda = essmSusy.get_Lambdax();
 
  double At;
  double TYt = essmSusy.get_TYu(2,2);
  
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

  double mQlsq = essmSusy.get_mq2(2,2);
  double mUrsq =  essmSusy.get_mu2(2,2);

  double MQQsq = mQlsq-mUrsq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2);
  double Xt = At - lambda*s/(Sqrt(2.0)*tb);
  double rt = MQQsq*MQQsq+4.0*mtop*mtop*Xt*Xt;

  double deriv = 2.0*yt*yt-0.25*g2*g2-0.6*0.25*g1*g1-0.1*gd1*gd1;
  deriv = deriv+(1.0/Sqrt(rt))*(0.25*MQQsq*(g1*g1-g2*g2)+2.0*yt*yt*Xt*Xt
				+4.0*mtop*mtop*Xt*(lambda*s*v1/(Sqrt(2.0)*v2*v2*v2)));
  deriv = deriv*(v2/2.0);

  return deriv;
}

double doCalcdmstop2sqdv3(genericE6SSM_soft_parameters essmSusy, double s, double tb)
{
  double v1 = essmSusy.get_vd();
  double v2 = essmSusy.get_vu(); 
  double v = Sqrt(v1*v1+v2*v2);

  double g1 = essmSusy.get_g1();
  double g2 = essmSusy.get_g2();
  double gd1 = essmSusy.get_gN();

  double yt = essmSusy.get_Yu( 2, 2);
  double mtop = yt*v2/Sqrt(2.0);
  // Option to fix m_t(M_t) = 165 GeV
  
  if(USEMTOFMT){
    mtop = 165;
    yt = Sqrt(2.0)*mtop/v2;
    essmSusy.set_Yu(2,2, yt);
  }

  
  double lambda = essmSusy.get_Lambdax();
 
  double At;
  double TYt = essmSusy.get_TYu(2,2);
  
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

  double mQlsq = essmSusy.get_mq2(2,2);
  double mUrsq =  essmSusy.get_mu2(2,2);

  double MQQsq = mQlsq-mUrsq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2);
  double Xt = At - lambda*s/(Sqrt(2.0)*tb);
  double rt = MQQsq*MQQsq+4.0*mtop*mtop*Xt*Xt;

  double deriv = 0.25*gd1*gd1+(4.0*mtop*mtop/Sqrt(rt))*Xt*(-lambda*v1/(Sqrt(2.0)*v2*s));
  deriv = 0.5*s*deriv;

  return deriv;
}

// The three functions below are purely used for testing the first
// derivatives of the stop squared masses. They calculate the 
// stop/top tadpole correction to the effective potential, in the
// form (1/v_i)\frac{\partial \Delta V}{\partial v_i}. The results
// should be compared against those obtained using the functions
// doCalcTadpolesESSMH<1,2,S> which are actually used in practical
// calculations. In the notation used in my notes, these functions
// calculate \Delta^{tadpole}_i for i = 1, 2, 3.
// Inputs:
//     genericE6SSM_soft_parameters const & essmSusy = object to calculate derivatives on
//     double s = the value of the VEV s to use
//     double tb = the value of tan(beta) to use
double calculateDeltaTadpole1(genericE6SSM_soft_parameters essmSusy, double s, double tb)
{
  double v1 = essmSusy.get_vd();
  double v2 = essmSusy.get_vu(); 
  double v = Sqrt(v1*v1+v2*v2);

  double g1 = essmSusy.get_g1();
  double g2 = essmSusy.get_g2();
  double gd1 = essmSusy.get_gN();

  double q = essmSusy.get_scale();

  double yt = essmSusy.get_Yu( 2, 2);
  double mtop = yt*v2/Sqrt(2.0);
  // Option to fix m_t(M_t) = 165 GeV
  
  if(USEMTOFMT){
    mtop = 165;
    yt = Sqrt(2.0)*mtop/v2;
    essmSusy.set_Yu(2,2, yt);
  }

  
  double lambda = essmSusy.get_Lambdax();
  
  double At;
  double TYt = essmSusy.get_TYu(2,2);
  
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

  double mQlsq = essmSusy.get_mq2(2,2);
  double mUrsq =  essmSusy.get_mu2(2,2);

  double MQQsq = mQlsq-mUrsq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2);
  double Xt = At - lambda*s/(Sqrt(2.0)*tb);
  double rt = MQQsq*MQQsq+4.0*mtop*mtop*Xt*Xt;

  // Contributions from auxiliary D-terms. Might in future
  // want to use effective charges instead of U(1)_N charges.
  double DeltaQ3 = (gd1*gd1/80.0)*(-3.0*v1*v1-2.0*v2*v2+5.0*s*s);
  double Deltau3 = (gd1*gd1/80.0)*(-3.0*v1*v1-2.0*v2*v2+5.0*s*s);

  // The stop squared masses, in the notation m_stop_1^2 has -ve
  // sign contribution (is lighter).
  double mstop1sq = 0.5*(mQlsq+mUrsq+0.125*(g2*g2+0.6*g1*g1)*(v1*v1-v2*v2)
		     +DeltaQ3+Deltau3+2.0*mtop*mtop-Sqrt(rt));
  double mstop2sq = 0.5*(mQlsq+mUrsq+0.125*(g2*g2+0.6*g1*g1)*(v1*v1-v2*v2)
		     +DeltaQ3+Deltau3+2.0*mtop*mtop+Sqrt(rt));

  bool speak = false;
  if (speak)
    {
      cout << "Calculated stop masses: " << endl;
      cout << "m_stop_1^2 = " << mstop1sq << ", ";
      cout << "m_stop_1 = " << Sqrt(mstop1sq) << ", ";
      cout << "m_stop_2^2 = " << mstop2sq << ", ";
      cout << "m_stop_2 = " << Sqrt(mstop2sq) << endl;
    }


  // The derivatives of the stop masses wrt the VEVs.
  double dmstop1sqdv1 = doCalcdmstop1sqdv1(essmSusy, s, tb);
  double dmstop2sqdv1 = doCalcdmstop2sqdv1(essmSusy, s, tb);

  // The derivative of the top mass wrt the VEVs. This is non-zero
  // only for dm_top/dv2.
  double dmtopsqdv1 = 0.0;

  if (speak)
    {
      cout << "Calculated derivatives of stop/top masses: " << endl;
      cout << "dm_stop_1^2/dv1 = " << dmstop1sqdv1 << ", ";
      cout << "dm_stop_2^2/dv1 = " << dmstop2sqdv1 << ", ";
      cout << "dm_top/dv1 = " << dmtopsqdv1 << endl;
    }

  // The tadpole contribution in the form (1/v_i)\frac{\partial \Delta V}{\partial v_i}.
  double dDeltaVdv1 = (3.0/(32.0*PI*PI))*(2.0*a0Peter(mstop1sq, q)*dmstop1sqdv1
					  +2.0*a0Peter(mstop2sq, q)*dmstop2sqdv1
					  -4.0*a0Peter(mtop*mtop, q)*dmtopsqdv1);

  double deltaTadpole1 = (1.0/v1)*dDeltaVdv1;

  if (speak)
    {
      cout << "Calculated tadpole contributions: " << endl;
      cout << "dDeltaV/dv1 = " << dDeltaVdv1 << ", ";
      cout << "Delta^{tadpole}_1 = " << deltaTadpole1 << endl;
    }

  return deltaTadpole1;

}

double calculateDeltaTadpole2(genericE6SSM_soft_parameters essmSusy, double s, double tb)
{
  double v1 = essmSusy.get_vd(); 
  double v2 = essmSusy.get_vu(); 
  double v = Sqrt(v1*v1+v2*v2);

  double g1 = essmSusy.get_g1();
  double g2 = essmSusy.get_g2();
  double gd1 = essmSusy.get_gN();

  double q = essmSusy.get_scale();

  double yt = essmSusy.get_Yu( 2, 2);
  double mtop = yt*v2/Sqrt(2.0);
  // Option to fix m_t(M_t) = 165 GeV
  
  if(USEMTOFMT){
    mtop = 165;
    yt = Sqrt(2.0)*mtop/v2;
    essmSusy.set_Yu(2,2, yt);
  }

  
  double lambda = essmSusy.get_Lambdax();
 
  double At;
  double TYt = essmSusy.get_TYu(2,2);
  
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

  double mQlsq = essmSusy.get_mq2(2,2);
  double mUrsq =  essmSusy.get_mu2(2,2);

  double MQQsq = mQlsq-mUrsq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2);
  double Xt = At - lambda*s/(Sqrt(2.0)*tb);
  double rt = MQQsq*MQQsq+4.0*mtop*mtop*Xt*Xt;

  // Contributions from auxiliary D-terms. Might in future
  // want to use effective charges instead of U(1)_N charges.
  double DeltaQ3 = (gd1*gd1/80.0)*(-3.0*v1*v1-2.0*v2*v2+5.0*s*s);
  double Deltau3 = (gd1*gd1/80.0)*(-3.0*v1*v1-2.0*v2*v2+5.0*s*s);

  // The stop squared masses, in the notation m_stop_1^2 has -ve
  // sign contribution (is lighter).
  double mstop1sq = 0.5*(mQlsq+mUrsq+0.125*(g2*g2+0.6*g1*g1)*(v1*v1-v2*v2)
		     +DeltaQ3+Deltau3+2.0*mtop*mtop-Sqrt(rt));
  double mstop2sq = 0.5*(mQlsq+mUrsq+0.125*(g2*g2+0.6*g1*g1)*(v1*v1-v2*v2)
		     +DeltaQ3+Deltau3+2.0*mtop*mtop+Sqrt(rt));

  bool speak = false;
  if (speak)
    {
      cout << "Calculated stop masses: " << endl;
      cout << "m_stop_1^2 = " << mstop1sq << ", ";
      cout << "m_stop_1 = " << Sqrt(mstop1sq) << ", ";
      cout << "m_stop_2^2 = " << mstop2sq << ", ";
      cout << "m_stop_2 = " << Sqrt(mstop2sq) << endl;
    }

  // The derivatives of the stop masses wrt the VEVs.
  double dmstop1sqdv2 = doCalcdmstop1sqdv2(essmSusy, s, tb);
  double dmstop2sqdv2 = doCalcdmstop2sqdv2(essmSusy, s, tb);

  // The derivative of the top mass wrt the VEVs. This is non-zero
  // only for dm_top/dv2.
  double dmtopsqdv2 = yt*yt*v2;

  if (speak)
    {
      cout << "Calculated derivatives of stop/top masses: " << endl;
      cout << "dm_stop_1^2/dv2 = " << dmstop1sqdv2 << ", ";
      cout << "dm_stop_2^2/dv2 = " << dmstop2sqdv2 << ", ";
      cout << "dm_top/dv2 = " << dmtopsqdv2 << endl;
    }

  // The tadpole contribution in the form (1/v_i)\frac{\partial \Delta V}{\partial v_i}.
  double dDeltaVdv2 = (3.0/(32.0*PI*PI))*(2.0*a0Peter(mstop1sq, q)*dmstop1sqdv2
					  +2.0*a0Peter(mstop2sq, q)*dmstop2sqdv2
					  -4.0*a0Peter(mtop*mtop, q)*dmtopsqdv2);

  double deltaTadpole2 = (1.0/v2)*dDeltaVdv2;

  if (speak)
    {
      cout << "Calculated tadpole contributions: " << endl;
      cout << "dDeltaV/dv2 = " << dDeltaVdv2 << ", ";
      cout << "Delta^{tadpole}_2 = " << deltaTadpole2 << endl;
    }

  return deltaTadpole2;

}

double calculateDeltaTadpole3(genericE6SSM_soft_parameters essmSusy, double s, double tb)
{
  double v1 = essmSusy.get_vd(); 
  double v2 = essmSusy.get_vu(); 
  double v = Sqrt(v1*v1+v2*v2);

  double g1 = essmSusy.get_g1();
  double g2 = essmSusy.get_g2();
  double gd1 = essmSusy.get_gN();

  double q = essmSusy.get_scale();

  double yt = essmSusy.get_Yu( 2, 2);
  double mtop = yt*v2/Sqrt(2.0);
  // Option to fix m_t(M_t) = 165 GeV
  
  if(USEMTOFMT){
    mtop = 165;
    yt = Sqrt(2.0)*mtop/v2;
    essmSusy.set_Yu(2,2, yt);
  }

  
  double lambda = essmSusy.get_Lambdax();
  
  double At;
  double TYt = essmSusy.get_TYu(2,2);
  
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

  double mQlsq = essmSusy.get_mq2(2,2);
  double mUrsq =  essmSusy.get_mu2(2,2);

  double MQQsq = mQlsq-mUrsq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2);
  double Xt = At - lambda*s/(Sqrt(2.0)*tb);
  double rt = MQQsq*MQQsq+4.0*mtop*mtop*Xt*Xt;

  // Contributions from auxiliary D-terms. Might in future
  // want to use effective charges instead of U(1)_N charges.
  double DeltaQ3 = (gd1*gd1/80.0)*(-3.0*v1*v1-2.0*v2*v2+5.0*s*s);
  double Deltau3 = (gd1*gd1/80.0)*(-3.0*v1*v1-2.0*v2*v2+5.0*s*s);

  // The stop squared masses, in the notation m_stop_1^2 has -ve
  // sign contribution (is lighter).
  double mstop1sq = 0.5*(mQlsq+mUrsq+0.125*(g2*g2+0.6*g1*g1)*(v1*v1-v2*v2)
		     +DeltaQ3+Deltau3+2.0*mtop*mtop-Sqrt(rt));
  double mstop2sq = 0.5*(mQlsq+mUrsq+0.125*(g2*g2+0.6*g1*g1)*(v1*v1-v2*v2)
		     +DeltaQ3+Deltau3+2.0*mtop*mtop+Sqrt(rt));

  bool speak = false;
  if (speak)
    {
      cout << "Calculated stop masses: " << endl;
      cout << "m_stop_1^2 = " << mstop1sq << ", ";
      cout << "m_stop_1 = " << Sqrt(mstop1sq) << ", ";
      cout << "m_stop_2^2 = " << mstop2sq << ", ";
      cout << "m_stop_2 = " << Sqrt(mstop2sq) << endl;
    }

  // The derivatives of the stop masses wrt the VEVs.
  double dmstop1sqdv3 = doCalcdmstop1sqdv3(essmSusy, s, tb);
  double dmstop2sqdv3 = doCalcdmstop2sqdv3(essmSusy, s, tb);

  // The derivative of the top mass wrt the VEVs. This is non-zero
  // only for dm_top/dv2.
  double dmtopsqdv3 = 0.0;

  if (speak)
    {
      cout << "Calculated derivatives of stop/top masses: " << endl;
      cout << "dm_stop_1^2/ds = " << dmstop1sqdv3 << ", ";
      cout << "dm_stop_2^2/ds = " << dmstop2sqdv3 << ", ";
      cout << "dm_top/ds = " << dmtopsqdv3 << endl;
    }

  // The tadpole contribution in the form (1/v_i)\frac{\partial \Delta V}{\partial v_i}.
  double dDeltaVdv3 = (3.0/(32.0*PI*PI))*(2.0*a0Peter(mstop1sq, q)*dmstop1sqdv3
					  +2.0*a0Peter(mstop2sq, q)*dmstop2sqdv3
					  -4.0*a0Peter(mtop*mtop, q)*dmtopsqdv3);

  double deltaTadpole3 = (1.0/s)*dDeltaVdv3;

  if (speak)
    {
      cout << "Calculated tadpole contributions: " << endl;
      cout << "dDeltaV/ds = " << dDeltaVdv3 << ", ";
      cout << "Delta^{tadpole}_3 = " << deltaTadpole3 << endl;
    }

  return deltaTadpole3;
}

// These functions calculate the leading one-loop corrections
// to the CP-even Higgs mass matrix in the basis (v_1, v_2, v_3=s),
// given by \Delta_{ij}'=\frac{\partial^2\Delta V}{\partial v_i\partial v_j}.
// Inputs:
//     genericE6SSM_soft_parameters essmSusy = object to calculate derivatives on
//     double s = the value of the VEV s to use
//     double tb = the value of tan(beta) to use.
double doCalcDeltaPrime11(genericE6SSM_soft_parameters essmSusy, double s, double tb)
{
  double v1 = essmSusy.get_vd(); 
  double v2 = essmSusy.get_vu();
  double v = Sqrt(v1*v1+v2*v2);

  double g1 = essmSusy.get_g1();
  double g2 = essmSusy.get_g2();
  double gbar = Sqrt(g2*g2+0.6*g1*g1);
  double gd1 = essmSusy.get_gN();

  // Logarithms are calculated at the given Q
  double q = essmSusy.get_scale();

  double yt = essmSusy.get_Yu( 2, 2);
  double mtop = yt*v2/Sqrt(2.0);
  // Option to fix m_t(M_t) = 165 GeV
  
  if(USEMTOFMT){
    mtop = 165;
    yt = Sqrt(2.0)*mtop/v2;
    essmSusy.set_Yu(2,2, yt);
  }

  
  double lambda = essmSusy.get_Lambdax();
 
  double At;
  double TYt = essmSusy.get_TYu(2,2);
  
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

  double mQlsq = essmSusy.get_mq2(2,2);
  double mUrsq =  essmSusy.get_mu2(2,2);

  double MQQsq = mQlsq-mUrsq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2);
  double Xt = At - lambda*s/(Sqrt(2.0)*tb);
  double rt = MQQsq*MQQsq+4.0*mtop*mtop*Xt*Xt;
  double RQQ = MQQsq*(g2*g2-g1*g1);

  // Squared stop masses are calculated using physical_ESSM to make sure
  // auxiliary D-terms are included.
  DoubleVector mstop(2), mstopsq(2), mD1sq(3), mD2sq(3);
  physical_ESSM(essmSusy, mstop, mstopsq, mD1sq, mD2sq, s, tb);

  if (mstopsq(1) < 0.0 || mstopsq(2) < 0.0)
    {
      mstopsq(1) = Abs(mstopsq(1));
      mstopsq(2) = Abs(mstopsq(2));
    }

  // Have to swap the mass ordering of the stops, because
  // physical_ESSM uses the opposite sign.
  double mstop1sq, mstop2sq;
  mstop1sq = mstopsq(2);
  mstop2sq = mstopsq(1);

  double term1, term2, term3, term4;

  term1 = v1*v1*Sqr(0.125*gbar*gbar-0.075*gd1*gd1)
    +(1.0/rt)*Sqr(0.125*v1*RQQ-2.0*mtop*mtop*Xt*s*lambda/(Sqrt(2.0)*v2));
  term1 = term1*Log(mstop1sq*mstop2sq/(q*q*q*q));

  term2 = (v1/(32.0*Sqrt(rt)))*(gbar*gbar-0.6*gd1*gd1)
    *(v1*RQQ-16.0*mtop*mtop*Xt*s*lambda/(Sqrt(2.0)*v2));
  term2 = term2*Log(mstop2sq/mstop1sq);

  term3 = (0.125*gbar*gbar-0.075*gd1*gd1)*(a0Peter(mstop1sq, q)+a0Peter(mstop2sq,q));

  term4 = (1.0/Sqrt(rt))*(4.0*RQQ+Sqr(g2*g2-g1*g1)*v1*v1+16.0*yt*yt*s*s*lambda*lambda);
  term4 = term4-(1.0/(rt*Sqrt(rt)))*Sqr(v1*RQQ-16.0*mtop*mtop*Xt*s*lambda/(Sqrt(2.0)*v2));
  term4 = term4*(a0Peter(mstop2sq,q)-a0Peter(mstop1sq,q))/32.0;

  double Del11prime = (3.0/(16.0*PI*PI))*(term1+term2+term3+term4);

  return Del11prime;
}

double doCalcDeltaPrime12(genericE6SSM_soft_parameters essmSusy, double s, double tb)
{
  double v1 = essmSusy.get_vd(); 
  double v2 = essmSusy.get_vu(); 
  double v = Sqrt(v1*v1+v2*v2);

  double g1 = essmSusy.get_g1();
  double g2 = essmSusy.get_g2();
  double gbar = Sqrt(g2*g2+0.6*g1*g1);
  double gd1 = essmSusy.get_gN();

  // Logarithms are calculated at the given Q
  double q = essmSusy.get_scale();

  double yt = essmSusy.get_Yu( 2, 2);
  double mtop = yt*v2/Sqrt(2.0);
  // Option to fix m_t(M_t) = 165 GeV
  
  if(USEMTOFMT){
    mtop = 165;
    yt = Sqrt(2.0)*mtop/v2;
    essmSusy.set_Yu(2,2, yt);
  }

  
  double lambda = essmSusy.get_Lambdax();
  
  double At;
  double TYt = essmSusy.get_TYu(2,2);
  
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

  double mQlsq = essmSusy.get_mq2(2,2);
  double mUrsq =  essmSusy.get_mu2(2,2);
  double mueff = lambda*s/Sqrt(2.0);

  double MQQsq = mQlsq-mUrsq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2);
  double Xt = At - lambda*s/(Sqrt(2.0)*tb);
  double rt = MQQsq*MQQsq+4.0*mtop*mtop*Xt*Xt;
  double RQQ = MQQsq*(g2*g2-g1*g1);

  // Squared stop masses are calculated using physical_ESSM to make sure
  // auxiliary D-terms are included.
  DoubleVector mstop(2), mstopsq(2), mD1sq(3), mD2sq(3);
  physical_ESSM(essmSusy, mstop, mstopsq, mD1sq, mD2sq, s, tb);

  if (mstopsq(1) < 0.0 || mstopsq(2) < 0.0)
    {
      mstopsq(1) = Abs(mstopsq(1));
      mstopsq(2) = Abs(mstopsq(2));
    }

  // Have to swap the mass ordering of the stops, because
  // physical_ESSM uses the opposite sign.
  double mstop1sq, mstop2sq;
  mstop1sq = mstopsq(2);
  mstop2sq = mstopsq(1);

  double term1, term2, term3;

  term1 = (0.125*gbar*gbar-0.075*gd1*gd1)*(yt*yt-0.125*gbar*gbar-0.05*gd1*gd1)
    +(1.0/rt)*(0.125*RQQ-yt*yt*Xt*mueff*tb)*(yt*yt*Xt*At-0.125*RQQ);
  term1 = term1*v1*v2*Log(mstop1sq*mstop2sq/(q*q*q*q));

  term2 = (0.125*gbar*gbar-0.075*gd1*gd1)*(yt*yt*Xt*At-0.125*RQQ)
    +(0.125*RQQ-yt*yt*Xt*mueff*tb)*(yt*yt-0.125*gbar*gbar-0.05*gd1*gd1);
  term2 = term2*(v1*v2/Sqrt(rt))*Log(mstop2sq/mstop1sq);

  term3 = (Sqr(g2*g2-g1*g1)*v1*v2/32.0+At*yt*yt*mueff)
    +(0.125*RQQ-yt*yt*Xt*mueff*tb)*(2.0*yt*yt*Xt*At-0.25*RQQ)*(v1*v2/rt);
  term3= term3*(1.0/Sqrt(rt))*(a0Peter(mstop2sq,q)-a0Peter(mstop1sq,q));

  double Del12prime = (3.0/(16.0*PI*PI))*(term1+term2-term3);

  return Del12prime;
}

double doCalcDeltaPrime13(genericE6SSM_soft_parameters essmSusy, double s, double tb)
{
  double v1 = essmSusy.get_vd(); 
  double v2 = essmSusy.get_vu();
  double v = Sqrt(v1*v1+v2*v2);

  double g1 = essmSusy.get_g1();
  double g2 = essmSusy.get_g2();
  double gbar = Sqrt(g2*g2+0.6*g1*g1);
  double gd1 = essmSusy.get_gN();

  // Logarithms are calculated at the given Q
  double q = essmSusy.get_scale();

  double yt = essmSusy.get_Yu( 2, 2);
  double mtop = yt*v2/Sqrt(2.0);
  // Option to fix m_t(M_t) = 165 GeV
  
  if(USEMTOFMT){
    mtop = 165;
    yt = Sqrt(2.0)*mtop/v2;
    essmSusy.set_Yu(2,2, yt);
  }

  
  double lambda = essmSusy.get_Lambdax();
  
  double At;
  double TYt = essmSusy.get_TYu(2,2);
  
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

  double mQlsq = essmSusy.get_mq2(2,2);
  double mUrsq =  essmSusy.get_mu2(2,2);
  double mueff = lambda*s/Sqrt(2.0);

  double MQQsq = mQlsq-mUrsq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2);
  double Xt = At - lambda*s/(Sqrt(2.0)*tb);
  double rt = MQQsq*MQQsq+4.0*mtop*mtop*Xt*Xt;
  double RQQ = MQQsq*(g2*g2-g1*g1);

  // Squared stop masses are calculated using physical_ESSM to make sure
  // auxiliary D-terms are included.
  DoubleVector mstop(2), mstopsq(2), mD1sq(3), mD2sq(3);
  physical_ESSM(essmSusy, mstop, mstopsq, mD1sq, mD2sq, s, tb);

  if (mstopsq(1) < 0.0 || mstopsq(2) < 0.0)
    {
      mstopsq(1) = Abs(mstopsq(1));
      mstopsq(2) = Abs(mstopsq(2));
    }

  // Have to swap the mass ordering of the stops, because
  // physical_ESSM uses the opposite sign.
  double mstop1sq, mstop2sq;
  mstop1sq = mstopsq(2);
  mstop2sq = mstopsq(1);

  double term1, term2, term3;

  term1 = 0.125*s*v1*gd1*gd1*(0.125*gbar*gbar-0.075*gd1*gd1)
    -(0.125*RQQ*v1-yt*yt*Xt*mueff*v2)*(2.0*mtop*mtop*Xt*lambda)/(rt*Sqrt(2.0)*tb);
  term1 = term1*Log(mstop1sq*mstop2sq/(q*q*q*q));

  term2 = (0.125*gbar*gbar-0.075*gd1*gd1)*(2.0*v1*mtop*mtop*Xt*lambda)/(tb*Sqrt(2.0*rt))
    -(0.125*s*gd1*gd1/Sqrt(rt))*(0.125*RQQ*v1-2.0*mtop*mtop*Xt*s*lambda/(Sqrt(2.0)*v2));
  term2 = term2*Log(mstop2sq/mstop1sq);

  term3 = (yt*yt*v1*lambda*mueff/Sqrt(2.0*rt))
    *(1.0-Xt*tb/mueff-4.0*Xt*Xt*mtop*mtop/rt+v1*v2*RQQ*Xt/(4.0*mueff*rt));
  term3 = term3*(a0Peter(mstop2sq, q)-a0Peter(mstop1sq,q));

  double Del13prime = (3.0/(16.0*PI*PI))*(term1-term2+term3);

  return Del13prime;
}

double doCalcDeltaPrime22(genericE6SSM_soft_parameters essmSusy, double s, double tb)
{
  double v1 = essmSusy.get_vd(); 
  double v2 = essmSusy.get_vu(); 
  double v = Sqrt(v1*v1+v2*v2);

  double g1 = essmSusy.get_g1();
  double g2 = essmSusy.get_g2();
  double gbar = Sqrt(g2*g2+0.6*g1*g1);
  double gd1 = essmSusy.get_gN();

  // Logarithms are calculated at the given Q
  double q = essmSusy.get_scale();

  double yt = essmSusy.get_Yu( 2, 2);
  double mtop = yt*v2/Sqrt(2.0);
  // Option to fix m_t(M_t) = 165 GeV
  
  if(USEMTOFMT){
    mtop = 165;
    yt = Sqrt(2.0)*mtop/v2;
    essmSusy.set_Yu(2,2, yt);
  }

  
  double lambda = essmSusy.get_Lambdax();
  
  double At;
  double TYt = essmSusy.get_TYu(2,2);
  
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

  double mQlsq = essmSusy.get_mq2(2,2);
  double mUrsq =  essmSusy.get_mu2(2,2);
  double mueff = lambda*s/Sqrt(2.0);

  double MQQsq = mQlsq-mUrsq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2);
  double Xt = At - lambda*s/(Sqrt(2.0)*tb);
  double rt = MQQsq*MQQsq+4.0*mtop*mtop*Xt*Xt;
  double RQQ = MQQsq*(g2*g2-g1*g1);

  // Squared stop masses are calculated using physical_ESSM to make sure
  // auxiliary D-terms are included.
  DoubleVector mstop(2), mstopsq(2), mD1sq(3), mD2sq(3);
  physical_ESSM(essmSusy, mstop, mstopsq, mD1sq, mD2sq, s, tb);

  if (mstopsq(1) < 0.0 || mstopsq(2) < 0.0)
    {
      mstopsq(1) = Abs(mstopsq(1));
      mstopsq(2) = Abs(mstopsq(2));
    }

  // Have to swap the mass ordering of the stops, because
  // physical_ESSM uses the opposite sign.
  double mstop1sq, mstop2sq;
  mstop1sq = mstopsq(2);
  mstop2sq = mstopsq(1);

  double term1, term2, term3, term4, term5,term6;

  term1 = Sqr(yt*yt-0.125*gbar*gbar-0.05*gd1*gd1)
    +Sqr(8.0*Xt*At*yt*yt-RQQ)/(64.0*rt);
  term1 = term1*v2*v2*Log(mstop1sq*mstop2sq/(q*q*q*q));

  term2 = (v2*v2/(4.0*Sqrt(rt)))*(yt*yt-0.125*gbar*gbar-0.05*gd1*gd1)
    *(8.0*yt*yt*Xt*At-RQQ)*Log(mstop2sq/mstop1sq);

  term3 = (yt*yt-0.125*gbar*gbar-0.05*gd1*gd1)*(a0Peter(mstop1sq,q)+a0Peter(mstop2sq,q));

  term4 = (1.0/Sqrt(rt))*(Sqr(g2*g2-g1*g1)*v2*v2/32.0-0.125*RQQ+yt*yt*At*At
			  -(v2*v2*Sqr(8.0*Xt*At*yt*yt-RQQ))/(32.0*rt));
  term4 = term4*(a0Peter(mstop2sq,q)-a0Peter(mstop1sq,q));

  term5 = 2.0*yt*yt*yt*yt*v2*v2*Log(mtop*mtop/(q*q));

  term6 = 2.0*yt*yt*a0Peter(mtop*mtop, q);

  double Del22prime = (3.0/(16.0*PI*PI))*(term1+term2+term3+term4-term5-term6);

  return Del22prime;
}

double doCalcDeltaPrime23(genericE6SSM_soft_parameters essmSusy, double s, double tb)
{
  double v1 = essmSusy.get_vd(); 
  double v2 = essmSusy.get_vu(); 
  double v = Sqrt(v1*v1+v2*v2);

  double g1 = essmSusy.get_g1();
  double g2 = essmSusy.get_g2();
  double gbar = Sqrt(g2*g2+0.6*g1*g1);
  double gd1 = essmSusy.get_gN();

  // Logarithms are calculated at the given Q
  double q = essmSusy.get_scale();

  double yt = essmSusy.get_Yu( 2, 2);
  double mtop = yt*v2/Sqrt(2.0);
  // Option to fix m_t(M_t) = 165 GeV
  
  if(USEMTOFMT){
    mtop = 165;
    yt = Sqrt(2.0)*mtop/v2;
    essmSusy.set_Yu(2,2, yt);
  }

  
  double lambda = essmSusy.get_Lambdax();
  
  double At;
  double TYt = essmSusy.get_TYu(2,2);
  
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

  double mQlsq = essmSusy.get_mq2(2,2);
  double mUrsq =  essmSusy.get_mu2(2,2);
  double mueff = lambda*s/Sqrt(2.0);

  double MQQsq = mQlsq-mUrsq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2);
  double Xt = At - lambda*s/(Sqrt(2.0)*tb);
  double rt = MQQsq*MQQsq+4.0*mtop*mtop*Xt*Xt;
  double RQQ = MQQsq*(g2*g2-g1*g1);

  // Squared stop masses are calculated using physical_ESSM to make sure
  // auxiliary D-terms are included.
  DoubleVector mstop(2), mstopsq(2), mD1sq(3), mD2sq(3);
  physical_ESSM(essmSusy, mstop, mstopsq, mD1sq, mD2sq, s, tb);

  if (mstopsq(1) < 0.0 || mstopsq(2) < 0.0)
    {
      mstopsq(1) = Abs(mstopsq(1));
      mstopsq(2) = Abs(mstopsq(2));
    }

  // Have to swap the mass ordering of the stops, because
  // physical_ESSM uses the opposite sign.
  double mstop1sq, mstop2sq;
  mstop1sq = mstopsq(2);
  mstop2sq = mstopsq(1);

  double term1, term2, term3;

  term1 = 0.125*s*gd1*gd1*(yt*yt-0.125*gbar*gbar-0.05*gd1*gd1)
    -(2.0*mtop*mtop*Xt*lambda)*(yt*yt*Xt*At-0.125*RQQ)/(Sqrt(2.0)*rt*tb);
  term1 = term1*v2*Log(mstop1sq*mstop2sq/(q*q*q*q));
 
  term2 = (1.0/Sqrt(rt))*(0.125*s*gd1*gd1*(yt*yt*Xt*At-0.125*RQQ)
			  -(2.0*mtop*mtop*Xt*lambda)*(yt*yt-0.125*gbar*gbar-0.05*gd1*gd1)/(Sqrt(2.0)*tb));
  term2 = term2*v2*Log(mstop2sq/mstop1sq);

  term3 = ((yt*yt*lambda*v1*At)/Sqrt(2.0*rt))*(1.0-4.0*mtop*mtop*Xt*Xt/rt+v2*v2*Xt*RQQ/(4.0*At*rt));
  term3 = term3*(a0Peter(mstop2sq,q)-a0Peter(mstop1sq,q));
  
  double Del23prime = (3.0/(16.0*PI*PI))*(term1+term2-term3);

  return Del23prime;
}

double doCalcDeltaPrime33(genericE6SSM_soft_parameters essmSusy, double s, double tb)
{
  double v1 = essmSusy.get_vd(); 
  double v2 = essmSusy.get_vu(); 
  double v = Sqrt(v1*v1+v2*v2);

  double g1 = essmSusy.get_g1();
  double g2 = essmSusy.get_g2();
  double gbar = Sqrt(g2*g2+0.6*g1*g1);
  double gd1 = essmSusy.get_gN();

  // Logarithms are calculated at the given Q
  double q = essmSusy.get_scale();

  double yt = essmSusy.get_Yu( 2, 2);
  double mtop = yt*v2/Sqrt(2.0);
  // Option to fix m_t(M_t) = 165 GeV
  
  if(USEMTOFMT){
    mtop = 165;
    yt = Sqrt(2.0)*mtop/v2;
    essmSusy.set_Yu(2,2, yt);
  }

  
  double lambda = essmSusy.get_Lambdax();
  
  double At;
  double TYt = essmSusy.get_TYu(2,2);
  
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

  double mQlsq = essmSusy.get_mq2(2,2);
  double mUrsq =  essmSusy.get_mu2(2,2);
  double mueff = lambda*s/Sqrt(2.0);

  double MQQsq = mQlsq-mUrsq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2);
  double Xt = At - lambda*s/(Sqrt(2.0)*tb);
  double rt = MQQsq*MQQsq+4.0*mtop*mtop*Xt*Xt;
  double RQQ = MQQsq*(g2*g2-g1*g1);

  // Squared stop masses are calculated using physical_ESSM to make sure
  // auxiliary D-terms are included.
  DoubleVector mstop(2), mstopsq(2), mD1sq(3), mD2sq(3);
  physical_ESSM(essmSusy, mstop, mstopsq, mD1sq, mD2sq, s, tb);

  if (mstopsq(1) < 0.0 || mstopsq(2) < 0.0)
    {
      mstopsq(1) = Abs(mstopsq(1));
      mstopsq(2) = Abs(mstopsq(2));
    }

  // Have to swap the mass ordering of the stops, because
  // physical_ESSM uses the opposite sign.
  double mstop1sq, mstop2sq;
  mstop1sq = mstopsq(2);
  mstop2sq = mstopsq(1);

  double term1, term2, term3, term4;

  term1 = Sqr(gd1*gd1*s)/64.0+(2.0*Sqr(mtop*mtop)*Xt*Xt*lambda*lambda)/(rt*tb*tb);
  term1 = term1*Log(mstop1sq*mstop2sq/(q*q*q*q));

  term2 = (gd1*gd1*mtop*mtop*Xt*mueff)/(2.0*Sqrt(rt)*tb);
  term2 = term2*Log(mstop2sq/mstop1sq);

  term3 = 0.125*gd1*gd1*(a0Peter(mstop1sq,q)+a0Peter(mstop2sq,q));

  term4 = ((lambda*lambda*mtop*mtop)/(Sqrt(rt)*tb*tb))*(1.0-4.0*Xt*Xt*mtop*mtop/rt);
  term4 = term4*(a0Peter(mstop2sq,q)-a0Peter(mstop1sq,q));

  double Del33prime = (3.0/(16.0*PI*PI))*(term1-term2+term3+term4);

  return Del33prime;
}

// These helper functions calculate the explicit derivatives
// of the tadpole contributions wrt the input parameters, i.e. returns
// the partial derivatives (1/v_i)*\frac{\partial^2 \Delta V}{\partial a\partial v_i}
// where the VEVs v_i, i = 1, 2, 3, are treated as being fixed.
// Inputs:
//     genericE6SSM_soft_parameters essmSusy = the object to calculate the derivatives on
//     double s = the value of the VEV s to use
//     double tb = the value of tan(beta) to use
double doCalcd2DeltaVdLambdadv1(genericE6SSM_soft_parameters essmSusy, double s, double tb)
{
  double v1 = essmSusy.get_vd();
  double v2 = essmSusy.get_vu(); 
  double v = Sqrt(v1*v1+v2*v2);

  double g1 = essmSusy.get_g1();
  double g2 = essmSusy.get_g2();
  double gbar = Sqrt(g2*g2+0.6*g1*g1);
  double gd1 = essmSusy.get_gN();

  // Logarithms are calculated at the given Q
  double q = essmSusy.get_scale();

  double yt = essmSusy.get_Yu( 2, 2);
  double mtop = yt*v2/Sqrt(2.0);
  // Option to fix m_t(M_t) = 165 GeV
  
  if(USEMTOFMT){
    mtop = 165;
    yt = Sqrt(2.0)*mtop/v2;
    essmSusy.set_Yu(2,2, yt);
  }

  
  double lambda = essmSusy.get_Lambdax();
  
  double At;
  double TYt = essmSusy.get_TYu(2,2);
  
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

  double mQlsq = essmSusy.get_mq2(2,2);
  double mUrsq =  essmSusy.get_mu2(2,2);
  double mueff = lambda*s/Sqrt(2.0);

  double MQQsq = mQlsq-mUrsq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2);
  double Xt = At - lambda*s/(Sqrt(2.0)*tb);
  double rt = MQQsq*MQQsq+4.0*mtop*mtop*Xt*Xt;

  // Squared stop masses are calculated using physical_ESSM to make sure
  // auxiliary D-terms are included.
  DoubleVector mstop(2), mstopsq(2), mD1sq(3), mD2sq(3);
  physical_ESSM(essmSusy, mstop, mstopsq, mD1sq, mD2sq, s, tb);

  if (mstopsq(1) < 0.0 || mstopsq(2) < 0.0)
    {
      mstopsq(1) = Abs(mstopsq(1));
      mstopsq(2) = Abs(mstopsq(2));
    }

  // Have to swap the mass ordering of the stops, because
  // physical_ESSM uses the opposite labelling.
  double mstop1sq, mstop2sq;
  mstop1sq = mstopsq(2);
  mstop2sq = mstopsq(1);

  // Calculate the explicit input parameter dependence of the stop masses
  double dmstop1sqdlambda = 0.5*(2.0*Sqrt(2.0)*mtop*mtop*Xt*s)/(Sqrt(rt)*tb);
  double dmstop2sqdlambda = -0.5*(2.0*Sqrt(2.0)*mtop*mtop*Xt*s)/(Sqrt(rt)*tb);

  // Calculate the derivative of the stop masses wrt the appropriate VEV
  double dmstop1sqdv1 = doCalcdmstop1sqdv1(essmSusy, s, tb);
  double dmstop2sqdv1 = doCalcdmstop2sqdv1(essmSusy, s, tb);

  // Calculate the explicit input parameter dependence of the derivatives
  // of the stop masses wrt the appropriate VEV
  double coeff = ((4.0*mtop*mtop)/(Sqrt(rt)))*(lambda*s/(Sqrt(2.0)*v1*v2))*(s/(Sqrt(2.0)*tb));
  coeff = coeff-(4.0*mtop*mtop*Xt*s)/(Sqrt(2.0*rt)*v1*v2);
  coeff = coeff-(((4.0*mtop*mtop*Xt*s)/(Sqrt(2.0*rt)*rt*tb))*
		 (0.25*MQQsq*(g1*g1-g2*g2)+4.0*mtop*mtop*Xt*lambda*s/(Sqrt(2.0)*v1*v2)));

  double d2mstop1sqdldv1 = -0.5*v1*coeff;
  double d2mstop2sqdldv1 = 0.5*v1*coeff;

  // Calculate the full second derivative of the tadpole contribution
  double term1 = dmstop1sqdlambda*dmstop1sqdv1*Log(mstop1sq/(q*q))+a0Peter(mstop1sq,q)*d2mstop1sqdldv1;
  double term2 = dmstop2sqdlambda*dmstop2sqdv1*Log(mstop2sq/(q*q))+a0Peter(mstop2sq,q)*d2mstop2sqdldv1;

  double d2DeltaV = (3.0/(16.0*PI*PI))*(term1+term2);
  d2DeltaV = (1.0/v1)*d2DeltaV;

  return d2DeltaV;
}

double doCalcd2DeltaVdLambdadv2(genericE6SSM_soft_parameters essmSusy, double s, double tb)
{
  double v1 = essmSusy.get_vd(); 
  double v2 = essmSusy.get_vu(); 
  double v = Sqrt(v1*v1+v2*v2);

  double g1 = essmSusy.get_g1();
  double g2 = essmSusy.get_g2();
  double gbar = Sqrt(g2*g2+0.6*g1*g1);
  double gd1 = essmSusy.get_gN();

  // Logarithms are calculated at the given Q
  double q = essmSusy.get_scale();

  double yt = essmSusy.get_Yu( 2, 2);
  double mtop = yt*v2/Sqrt(2.0);
  // Option to fix m_t(M_t) = 165 GeV
  
  if(USEMTOFMT){
    mtop = 165;
    yt = Sqrt(2.0)*mtop/v2;
    essmSusy.set_Yu(2,2, yt);
  }

  
  double lambda = essmSusy.get_Lambdax();
  
  double At;
  double TYt = essmSusy.get_TYu(2,2);
  
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

  double mQlsq = essmSusy.get_mq2(2,2);
  double mUrsq =  essmSusy.get_mu2(2,2);
  double mueff = lambda*s/Sqrt(2.0);

  double MQQsq = mQlsq-mUrsq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2);
  double Xt = At - lambda*s/(Sqrt(2.0)*tb);
  double rt = MQQsq*MQQsq+4.0*mtop*mtop*Xt*Xt;

  // Squared stop masses are calculated using physical_ESSM to make sure
  // auxiliary D-terms are included.
  DoubleVector mstop(2), mstopsq(2), mD1sq(3), mD2sq(3);
  physical_ESSM(essmSusy, mstop, mstopsq, mD1sq, mD2sq, s, tb);

  if (mstopsq(1) < 0.0 || mstopsq(2) < 0.0)
    {
      mstopsq(1) = Abs(mstopsq(1));
      mstopsq(2) = Abs(mstopsq(2));
    }

  // Have to swap the mass ordering of the stops, because
  // physical_ESSM uses the opposite labelling.
  double mstop1sq, mstop2sq;
  mstop1sq = mstopsq(2);
  mstop2sq = mstopsq(1);

  // Calculate the explicit input parameter dependence of the stop masses
  double dmstop1sqdlambda = 0.5*(2.0*Sqrt(2.0)*mtop*mtop*Xt*s)/(Sqrt(rt)*tb);
  double dmstop2sqdlambda = -0.5*(2.0*Sqrt(2.0)*mtop*mtop*Xt*s)/(Sqrt(rt)*tb);

  // Calculate the derivative of the stop masses wrt the appropriate VEV
  double dmstop1sqdv2 = doCalcdmstop1sqdv2(essmSusy, s, tb);
  double dmstop2sqdv2 = doCalcdmstop2sqdv2(essmSusy, s, tb);

  // Calculate the explicit input parameter dependence of the derivatives
  // of the stop masses wrt the appropriate VEV
  double coeff = (1.0/Sqrt(rt))*(4.0*mtop*mtop*Xt*s*v1/(Sqrt(2.0)*v2*v2*v2)
				 -4.0*mtop*mtop*s*lambda*s*v1/(2.0*tb*v2*v2*v2)
				 -4.0*yt*yt*Xt*s/(Sqrt(2.0)*tb));
  coeff = coeff-((4.0*mtop*mtop*Xt*s)/(Sqrt(2.0*rt)*rt*tb))*
    (0.25*MQQsq*(g2*g2-g1*g1)-2.0*yt*yt*Xt*Xt-4.0*mtop*mtop*Xt*lambda*s*v1/(Sqrt(2.0)*v2*v2*v2));

  double d2mstop1sqdldv2 = -0.5*v2*coeff;
  double d2mstop2sqdldv2 = 0.5*v2*coeff;

  // Calculate the full second derivative of the tadpole contribution
  double term1 = dmstop1sqdlambda*dmstop1sqdv2*Log(mstop1sq/(q*q))+a0Peter(mstop1sq,q)*d2mstop1sqdldv2;
  double term2 = dmstop2sqdlambda*dmstop2sqdv2*Log(mstop2sq/(q*q))+a0Peter(mstop2sq,q)*d2mstop2sqdldv2;

  double d2DeltaV = (3.0/(16.0*PI*PI))*(term1+term2);
  d2DeltaV = (1.0/v2)*d2DeltaV;

  return d2DeltaV;
}

double doCalcd2DeltaVdLambdadv3(genericE6SSM_soft_parameters essmSusy, double s, double tb)
{
  double v1 = essmSusy.get_vd(); 
  double v2 = essmSusy.get_vu();
  double v = Sqrt(v1*v1+v2*v2);

  double g1 = essmSusy.get_g1();
  double g2 = essmSusy.get_g2();
  double gbar = Sqrt(g2*g2+0.6*g1*g1);
  double gd1 = essmSusy.get_gN();

  // Logarithms are calculated at the given Q
  double q = essmSusy.get_scale();

  double yt = essmSusy.get_Yu( 2, 2);
  double mtop = yt*v2/Sqrt(2.0);
  // Option to fix m_t(M_t) = 165 GeV
  
  if(USEMTOFMT){
    mtop = 165;
    yt = Sqrt(2.0)*mtop/v2;
    essmSusy.set_Yu(2,2, yt);
  }

  
  double lambda = essmSusy.get_Lambdax();
  
  double At;
  double TYt = essmSusy.get_TYu(2,2);
  
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

  double mQlsq = essmSusy.get_mq2(2,2);
  double mUrsq =  essmSusy.get_mu2(2,2);
  double mueff = lambda*s/Sqrt(2.0);

  double MQQsq = mQlsq-mUrsq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2);
  double Xt = At - lambda*s/(Sqrt(2.0)*tb);
  double rt = MQQsq*MQQsq+4.0*mtop*mtop*Xt*Xt;

  // Squared stop masses are calculated using physical_ESSM to make sure
  // auxiliary D-terms are included.
  DoubleVector mstop(2), mstopsq(2), mD1sq(3), mD2sq(3);
  physical_ESSM(essmSusy, mstop, mstopsq, mD1sq, mD2sq, s, tb);

  if (mstopsq(1) < 0.0 || mstopsq(2) < 0.0)
    {
      mstopsq(1) = Abs(mstopsq(1));
      mstopsq(2) = Abs(mstopsq(2));
    }

  // Have to swap the mass ordering of the stops, because
  // physical_ESSM uses the opposite labelling.
  double mstop1sq, mstop2sq;
  mstop1sq = mstopsq(2);
  mstop2sq = mstopsq(1);

  // Calculate the explicit input parameter dependence of the stop masses
  double dmstop1sqdlambda = 0.5*(2.0*Sqrt(2.0)*mtop*mtop*Xt*s)/(Sqrt(rt)*tb);
  double dmstop2sqdlambda = -0.5*(2.0*Sqrt(2.0)*mtop*mtop*Xt*s)/(Sqrt(rt)*tb);

  // Calculate the derivative of the stop masses wrt the appropriate VEV
  double dmstop1sqdv3 = doCalcdmstop1sqdv3(essmSusy, s, tb);
  double dmstop2sqdv3 = doCalcdmstop2sqdv3(essmSusy, s, tb);

  // Calculate the explicit input parameter dependence of the derivatives
  // of the stop masses wrt the appropriate VEV
  double coeff = s*lambda*v1*mtop*mtop/(2.0*Sqrt(rt)*tb*v2*s);
  coeff = coeff-mtop*mtop*Xt*v1/(Sqrt(2.0*rt)*v2*s);
  coeff = coeff-4.0*Sqr(mtop*mtop)*Xt*Xt*s*lambda*v1/(2.0*Sqrt(rt)*rt*tb*v2*s);

  double d2mstop1sqdldv3 = -2.0*s*coeff;
  double d2mstop2sqdldv3 = 2.0*s*coeff;

  // Calculate the full second derivative of the tadpole contribution
  double term1 = dmstop1sqdlambda*dmstop1sqdv3*Log(mstop1sq/(q*q))+a0Peter(mstop1sq,q)*d2mstop1sqdldv3;
  double term2 = dmstop2sqdlambda*dmstop2sqdv3*Log(mstop2sq/(q*q))+a0Peter(mstop2sq,q)*d2mstop2sqdldv3;

  double d2DeltaV = (3.0/(16.0*PI*PI))*(term1+term2);
  d2DeltaV = (1.0/s)*d2DeltaV;

  return d2DeltaV;
}

double doCalcd2DeltaVdAtdv1(genericE6SSM_soft_parameters essmSusy, double s, double tb)
{
  double v1 = essmSusy.get_vd(); 
  double v2 = essmSusy.get_vu(); 
  double v = Sqrt(v1*v1+v2*v2);

  double g1 = essmSusy.get_g1();
  double g2 = essmSusy.get_g2();
  double gbar = Sqrt(g2*g2+0.6*g1*g1);
  double gd1 = essmSusy.get_gN();

  // Logarithms are calculated at the given Q
  double q = essmSusy.get_scale();

  double yt = essmSusy.get_Yu( 2, 2);
  double mtop = yt*v2/Sqrt(2.0);
  // Option to fix m_t(M_t) = 165 GeV
  
  if(USEMTOFMT){
    mtop = 165;
    yt = Sqrt(2.0)*mtop/v2;
    essmSusy.set_Yu(2,2, yt);
  }

  
  double lambda = essmSusy.get_Lambdax();
  
  double At;
  double TYt = essmSusy.get_TYu(2,2);
  
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

  double mQlsq = essmSusy.get_mq2(2,2);
  double mUrsq =  essmSusy.get_mu2(2,2);
  double mueff = lambda*s/Sqrt(2.0);

  double MQQsq = mQlsq-mUrsq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2);
  double Xt = At - lambda*s/(Sqrt(2.0)*tb);
  double rt = MQQsq*MQQsq+4.0*mtop*mtop*Xt*Xt;

  // Squared stop masses are calculated using physical_ESSM to make sure
  // auxiliary D-terms are included.
  DoubleVector mstop(2), mstopsq(2), mD1sq(3), mD2sq(3);
  physical_ESSM(essmSusy, mstop, mstopsq, mD1sq, mD2sq, s, tb);

  if (mstopsq(1) < 0.0 || mstopsq(2) < 0.0)
    {
      mstopsq(1) = Abs(mstopsq(1));
      mstopsq(2) = Abs(mstopsq(2));
    }

  // Have to swap the mass ordering of the stops, because
  // physical_ESSM uses the opposite labelling.
  double mstop1sq, mstop2sq;
  mstop1sq = mstopsq(2);
  mstop2sq = mstopsq(1);

  // Calculate the explicit input parameter dependence of the stop masses
  double dmstop1sqdAt = -0.5*(4.0*mtop*mtop*Xt)/Sqrt(rt);
  double dmstop2sqdAt = 0.5*(4.0*mtop*mtop*Xt)/Sqrt(rt);

  // Calculate the derivative of the stop masses wrt the appropriate VEV
  double dmstop1sqdv1 = doCalcdmstop1sqdv1(essmSusy, s, tb);
  double dmstop2sqdv1 = doCalcdmstop2sqdv1(essmSusy, s, tb);

  // Calculate the explicit input parameter dependence of the derivatives
  // of the stop masses wrt the appropriate VEV
  double coeff = (4.0*mtop*mtop*Xt/(Sqrt(rt)*rt))*(0.25*MQQsq*(g1*g1-g2*g2)
						   +4.0*mtop*mtop*Xt*lambda*s/(Sqrt(2.0)*v1*v2));
  coeff = coeff-(4.0*mtop*mtop*lambda*s/(Sqrt(2.0*rt)*v1*v2));

  double d2mstop1sqdAtdv1 = -0.5*v1*coeff;
  double d2mstop2sqdAtdv1 = 0.5*v1*coeff;

  // Calculate the full second derivative of the tadpole contribution
  double term1 = dmstop1sqdAt*dmstop1sqdv1*Log(mstop1sq/(q*q))+a0Peter(mstop1sq,q)*d2mstop1sqdAtdv1;
  double term2 = dmstop2sqdAt*dmstop2sqdv1*Log(mstop2sq/(q*q))+a0Peter(mstop2sq,q)*d2mstop2sqdAtdv1;

  double d2DeltaV = (3.0/(16.0*PI*PI))*(term1+term2);
  d2DeltaV = (1.0/v1)*d2DeltaV;

  return d2DeltaV;
}

double doCalcd2DeltaVdAtdv2(genericE6SSM_soft_parameters essmSusy, double s, double tb)
{
  double v1 = essmSusy.get_vd(); 
  double v2 = essmSusy.get_vu(); 
  double v = Sqrt(v1*v1+v2*v2);

  double g1 = essmSusy.get_g1();
  double g2 = essmSusy.get_g2();
  double gbar = Sqrt(g2*g2+0.6*g1*g1);
  double gd1 = essmSusy.get_gN();

  // Logarithms are calculated at the given Q
  double q = essmSusy.get_scale();

  double yt = essmSusy.get_Yu( 2, 2);
  double mtop = yt*v2/Sqrt(2.0);
  // Option to fix m_t(M_t) = 165 GeV
  
  if(USEMTOFMT){
    mtop = 165;
    yt = Sqrt(2.0)*mtop/v2;
    essmSusy.set_Yu(2,2, yt);
  }

  
  double lambda = essmSusy.get_Lambdax();
  
  double At;
  double TYt = essmSusy.get_TYu(2,2);
  
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

  double mQlsq = essmSusy.get_mq2(2,2);
  double mUrsq =  essmSusy.get_mu2(2,2);
  double mueff = lambda*s/Sqrt(2.0);

  double MQQsq = mQlsq-mUrsq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2);
  double Xt = At - lambda*s/(Sqrt(2.0)*tb);
  double rt = MQQsq*MQQsq+4.0*mtop*mtop*Xt*Xt;

  // Squared stop masses are calculated using physical_ESSM to make sure
  // auxiliary D-terms are included.
  DoubleVector mstop(2), mstopsq(2), mD1sq(3), mD2sq(3);
  physical_ESSM(essmSusy, mstop, mstopsq, mD1sq, mD2sq, s, tb);

  if (mstopsq(1) < 0.0 || mstopsq(2) < 0.0)
    {
      mstopsq(1) = Abs(mstopsq(1));
      mstopsq(2) = Abs(mstopsq(2));
    }

  // Have to swap the mass ordering of the stops, because
  // physical_ESSM uses the opposite labelling.
  double mstop1sq, mstop2sq;
  mstop1sq = mstopsq(2);
  mstop2sq = mstopsq(1);

  // Calculate the explicit input parameter dependence of the stop masses
  double dmstop1sqdAt = -0.5*(4.0*mtop*mtop*Xt)/Sqrt(rt);
  double dmstop2sqdAt = 0.5*(4.0*mtop*mtop*Xt)/Sqrt(rt);

  // Calculate the derivative of the stop masses wrt the appropriate VEV
  double dmstop1sqdv2 = doCalcdmstop1sqdv2(essmSusy, s, tb);
  double dmstop2sqdv2 = doCalcdmstop2sqdv2(essmSusy, s, tb);

  // Calculate the explicit input parameter dependence of the derivatives
  // of the stop masses wrt the appropriate VEV
  double coeff = (1.0/Sqrt(rt))*(4.0*yt*yt*Xt+4.0*mtop*mtop*lambda*s*v1/(Sqrt(2.0)*v2*v2*v2));
  coeff = coeff+((4.0*mtop*mtop*Xt)/(rt*Sqrt(rt)))*(0.25*MQQsq*(g2*g2-g1*g1)-2.0*yt*yt*Xt*Xt
						    -4.0*mtop*mtop*Xt*lambda*s*v1/(Sqrt(2.0)*v2*v2*v2));

  double d2mstop1sqdAtdv2 = -0.5*v2*coeff;
  double d2mstop2sqdAtdv2 = 0.5*v2*coeff;

  // Calculate the full second derivative of the tadpole contribution
  double term1 = dmstop1sqdAt*dmstop1sqdv2*Log(mstop1sq/(q*q))+a0Peter(mstop1sq,q)*d2mstop1sqdAtdv2;
  double term2 = dmstop2sqdAt*dmstop2sqdv2*Log(mstop2sq/(q*q))+a0Peter(mstop2sq,q)*d2mstop2sqdAtdv2;

  double d2DeltaV = (3.0/(16.0*PI*PI))*(term1+term2);
  d2DeltaV = (1.0/v2)*d2DeltaV;

  return d2DeltaV;

}

double doCalcd2DeltaVdAtdv3(genericE6SSM_soft_parameters essmSusy, double s, double tb)
{
  double v1 = essmSusy.get_vd();
  double v2 = essmSusy.get_vu(); 
  double v = Sqrt(v1*v1+v2*v2);

  double g1 = essmSusy.get_g1();
  double g2 = essmSusy.get_g2();
  double gbar = Sqrt(g2*g2+0.6*g1*g1);
  double gd1 = essmSusy.get_gN();

  // Logarithms are calculated at the given Q
  double q = essmSusy.get_scale();

  double yt = essmSusy.get_Yu( 2, 2);
  double mtop = yt*v2/Sqrt(2.0);
  // Option to fix m_t(M_t) = 165 GeV
  
  if(USEMTOFMT){
    mtop = 165;
    yt = Sqrt(2.0)*mtop/v2;
    essmSusy.set_Yu(2,2, yt);
  }

  
  double lambda = essmSusy.get_Lambdax();
  
  double At;
  double TYt = essmSusy.get_TYu(2,2);
  
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

  double mQlsq = essmSusy.get_mq2(2,2);
  double mUrsq =  essmSusy.get_mu2(2,2);
  double mueff = lambda*s/Sqrt(2.0);

  double MQQsq = mQlsq-mUrsq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2);
  double Xt = At - lambda*s/(Sqrt(2.0)*tb);
  double rt = MQQsq*MQQsq+4.0*mtop*mtop*Xt*Xt;

  // Squared stop masses are calculated using physical_ESSM to make sure
  // auxiliary D-terms are included.
  DoubleVector mstop(2), mstopsq(2), mD1sq(3), mD2sq(3);
  physical_ESSM(essmSusy, mstop, mstopsq, mD1sq, mD2sq, s, tb);

  if (mstopsq(1) < 0.0 || mstopsq(2) < 0.0)
    {
      mstopsq(1) = Abs(mstopsq(1));
      mstopsq(2) = Abs(mstopsq(2));
    }

  // Have to swap the mass ordering of the stops, because
  // physical_ESSM uses the opposite labelling.
  double mstop1sq, mstop2sq;
  mstop1sq = mstopsq(2);
  mstop2sq = mstopsq(1);

  // Calculate the explicit input parameter dependence of the stop masses
  double dmstop1sqdAt = -0.5*(4.0*mtop*mtop*Xt)/Sqrt(rt);
  double dmstop2sqdAt = 0.5*(4.0*mtop*mtop*Xt)/Sqrt(rt);

  // Calculate the derivative of the stop masses wrt the appropriate VEV
  double dmstop1sqdv3 = doCalcdmstop1sqdv3(essmSusy, s, tb);
  double dmstop2sqdv3 = doCalcdmstop2sqdv3(essmSusy, s, tb);

  // Calculate the explicit input parameter dependence of the derivatives
  // of the stop masses wrt the appropriate VEV
  double coeff = 4.0*Sqr(mtop*mtop)*Xt*Xt*lambda*v1/(Sqrt(2.0*rt)*rt*v2*s);
  coeff = coeff-mtop*mtop*lambda*v1/(Sqrt(2.0*rt)*v2*s);

  double d2mstop1sqdAtdv3 = -2.0*s*coeff;
  double d2mstop2sqdAtdv3 = 2.0*s*coeff;

  // Calculate the full second derivative of the tadpole contribution
  double term1 = dmstop1sqdAt*dmstop1sqdv3*Log(mstop1sq/(q*q))+a0Peter(mstop1sq,q)*d2mstop1sqdAtdv3;
  double term2 = dmstop2sqdAt*dmstop2sqdv3*Log(mstop2sq/(q*q))+a0Peter(mstop2sq,q)*d2mstop2sqdAtdv3;

  double d2DeltaV = (3.0/(16.0*PI*PI))*(term1+term2);
  d2DeltaV = (1.0/s)*d2DeltaV;

  return d2DeltaV;

}

double doCalcd2DeltaVdmQlsqdv1(genericE6SSM_soft_parameters essmSusy, double s, double tb)
{
  double v1 = essmSusy.get_vd(); 
  double v2 = essmSusy.get_vu(); 
  double v = Sqrt(v1*v1+v2*v2);

  double g1 = essmSusy.get_g1();
  double g2 = essmSusy.get_g2();
  double gbar = Sqrt(g2*g2+0.6*g1*g1);
  double gd1 = essmSusy.get_gN();

  // Logarithms are calculated at the given Q
  double q = essmSusy.get_scale();

  double yt = essmSusy.get_Yu( 2, 2);
  double mtop = yt*v2/Sqrt(2.0);
  // Option to fix m_t(M_t) = 165 GeV
  
  if(USEMTOFMT){
    mtop = 165;
    yt = Sqrt(2.0)*mtop/v2;
    essmSusy.set_Yu(2,2, yt);
  }

  
  double lambda = essmSusy.get_Lambdax();
  
  double At;
  double TYt = essmSusy.get_TYu(2,2);
  
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

  double mQlsq = essmSusy.get_mq2(2,2);
  double mUrsq =  essmSusy.get_mu2(2,2);
  double mueff = lambda*s/Sqrt(2.0);

  double MQQsq = mQlsq-mUrsq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2);
  double Xt = At - lambda*s/(Sqrt(2.0)*tb);
  double rt = MQQsq*MQQsq+4.0*mtop*mtop*Xt*Xt;

  // Squared stop masses are calculated using physical_ESSM to make sure
  // auxiliary D-terms are included.
  DoubleVector mstop(2), mstopsq(2), mD1sq(3), mD2sq(3);
  physical_ESSM(essmSusy, mstop, mstopsq, mD1sq, mD2sq, s, tb);

  if (mstopsq(1) < 0.0 || mstopsq(2) < 0.0)
    {
      mstopsq(1) = Abs(mstopsq(1));
      mstopsq(2) = Abs(mstopsq(2));
    }

  // Have to swap the mass ordering of the stops, because
  // physical_ESSM uses the opposite labelling.
  double mstop1sq, mstop2sq;
  mstop1sq = mstopsq(2);
  mstop2sq = mstopsq(1);

  // Calculate the explicit input parameter dependence of the stop masses
  double dmstop1sqdmQlsq = 0.5*(1.0-MQQsq/Sqrt(rt));
  double dmstop2sqdmQlsq =0.5*(1.0+MQQsq/Sqrt(rt));

  // Calculate the derivative of the stop masses wrt the appropriate VEV
  double dmstop1sqdv1 = doCalcdmstop1sqdv1(essmSusy, s, tb);
  double dmstop2sqdv1 = doCalcdmstop2sqdv1(essmSusy, s, tb);

  // Calculate the explicit input parameter dependence of the derivatives
  // of the stop masses wrt the appropriate VEV
  double coeff = 0.25*(g2*g2-g1*g1)/Sqrt(rt);
  coeff = coeff+(MQQsq/(rt*Sqrt(rt)))*(0.25*MQQsq*(g1*g1-g2*g2)
				       +4.0*mtop*mtop*Xt*lambda*s/(Sqrt(2.0)*v1*v2));

  double d2mstop1sqdmQlsqdv1 = -0.5*v1*coeff;
  double d2mstop2sqdmQlsqdv1 = 0.5*v1*coeff;

  // Calculate the full second derivative of the tadpole contribution
  double term1 = dmstop1sqdmQlsq*dmstop1sqdv1*Log(mstop1sq/(q*q))+a0Peter(mstop1sq,q)*d2mstop1sqdmQlsqdv1;
  double term2 = dmstop2sqdmQlsq*dmstop2sqdv1*Log(mstop2sq/(q*q))+a0Peter(mstop2sq,q)*d2mstop2sqdmQlsqdv1;

  double d2DeltaV = (3.0/(16.0*PI*PI))*(term1+term2);
  d2DeltaV = (1.0/v1)*d2DeltaV;

  return d2DeltaV;

}

double doCalcd2DeltaVdmQlsqdv2(genericE6SSM_soft_parameters essmSusy, double s, double tb)
{
  double v1 = essmSusy.get_vd();
  double v2 = essmSusy.get_vu(); 
  double v = Sqrt(v1*v1+v2*v2);
  
  double g1 = essmSusy.get_g1();
  double g2 = essmSusy.get_g2();
  double gbar = Sqrt(g2*g2+0.6*g1*g1);
  double gd1 = essmSusy.get_gN();

  // Logarithms are calculated at the given Q
  double q = essmSusy.get_scale();

  double yt = essmSusy.get_Yu( 2, 2);
  double mtop = yt*v2/Sqrt(2.0);
  // Option to fix m_t(M_t) = 165 GeV
  
  if(USEMTOFMT){
    mtop = 165;
    yt = Sqrt(2.0)*mtop/v2;
    essmSusy.set_Yu(2,2, yt);
  }

  
  double lambda = essmSusy.get_Lambdax();
 
  double At;
  double TYt = essmSusy.get_TYu(2,2);
  
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

  double mQlsq = essmSusy.get_mq2(2,2);
  double mUrsq =  essmSusy.get_mu2(2,2);
  double mueff = lambda*s/Sqrt(2.0);

  double MQQsq = mQlsq-mUrsq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2);
  double Xt = At - lambda*s/(Sqrt(2.0)*tb);
  double rt = MQQsq*MQQsq+4.0*mtop*mtop*Xt*Xt;

  // Squared stop masses are calculated using physical_ESSM to make sure
  // auxiliary D-terms are included.
  DoubleVector mstop(2), mstopsq(2), mD1sq(3), mD2sq(3);
  physical_ESSM(essmSusy, mstop, mstopsq, mD1sq, mD2sq, s, tb);

  if (mstopsq(1) < 0.0 || mstopsq(2) < 0.0)
    {
      mstopsq(1) = Abs(mstopsq(1));
      mstopsq(2) = Abs(mstopsq(2));
    }

  // Have to swap the mass ordering of the stops, because
  // physical_ESSM uses the opposite labelling.
  double mstop1sq, mstop2sq;
  mstop1sq = mstopsq(2);
  mstop2sq = mstopsq(1);

  // Calculate the explicit input parameter dependence of the stop masses
  double dmstop1sqdmQlsq = 0.5*(1.0-MQQsq/Sqrt(rt));
  double dmstop2sqdmQlsq =0.5*(1.0+MQQsq/Sqrt(rt));

  // Calculate the derivative of the stop masses wrt the appropriate VEV
  double dmstop1sqdv2 = doCalcdmstop1sqdv2(essmSusy, s, tb);
  double dmstop2sqdv2 = doCalcdmstop2sqdv2(essmSusy, s, tb);

  // Calculate the explicit input parameter dependence of the derivatives
  // of the stop masses wrt the appropriate VEV
  double coeff = 0.25*(g1*g1-g2*g2)/Sqrt(rt);
  coeff = coeff+(MQQsq/(rt*Sqrt(rt)))*(0.25*MQQsq*(g2*g2-g1*g1)-2.0*yt*yt*Xt*Xt
				       -4.0*mtop*mtop*Xt*lambda*s*v1/(Sqrt(2.0)*v2*v2*v2));

  double d2mstop1sqdmQlsqdv2 = -0.5*v2*coeff;
  double d2mstop2sqdmQlsqdv2 = 0.5*v2*coeff;

  // Calculate the full second derivative of the tadpole contribution
  double term1 = dmstop1sqdmQlsq*dmstop1sqdv2*Log(mstop1sq/(q*q))+a0Peter(mstop1sq,q)*d2mstop1sqdmQlsqdv2;
  double term2 = dmstop2sqdmQlsq*dmstop2sqdv2*Log(mstop2sq/(q*q))+a0Peter(mstop2sq,q)*d2mstop2sqdmQlsqdv2;

  double d2DeltaV = (3.0/(16.0*PI*PI))*(term1+term2);
  d2DeltaV = (1.0/v2)*d2DeltaV;

  return d2DeltaV;

}

double doCalcd2DeltaVdmQlsqdv3(genericE6SSM_soft_parameters essmSusy, double s, double tb)
{
  double v1 = essmSusy.get_vd(); 
  double v2 = essmSusy.get_vu(); 
  double v = Sqrt(v1*v1+v2*v2);

  double g1 = essmSusy.get_g1();
  double g2 = essmSusy.get_g2();
  double gbar = Sqrt(g2*g2+0.6*g1*g1);
  double gd1 = essmSusy.get_gN();

  // Logarithms are calculated at the given Q
  double q = essmSusy.get_scale();

  double yt = essmSusy.get_Yu( 2, 2);
  double mtop = yt*v2/Sqrt(2.0);
  // Option to fix m_t(M_t) = 165 GeV
  
  if(USEMTOFMT){
    mtop = 165;
    yt = Sqrt(2.0)*mtop/v2;
    essmSusy.set_Yu(2,2, yt);
  }

  
  double lambda = essmSusy.get_Lambdax();
  
  double At;
  double TYt = essmSusy.get_TYu(2,2);
  
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

  double mQlsq = essmSusy.get_mq2(2,2);
  double mUrsq =  essmSusy.get_mu2(2,2);
  double mueff = lambda*s/Sqrt(2.0);

  double MQQsq = mQlsq-mUrsq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2);
  double Xt = At - lambda*s/(Sqrt(2.0)*tb);
  double rt = MQQsq*MQQsq+4.0*mtop*mtop*Xt*Xt;

  // Squared stop masses are calculated using physical_ESSM to make sure
  // auxiliary D-terms are included.
  DoubleVector mstop(2), mstopsq(2), mD1sq(3), mD2sq(3);
  physical_ESSM(essmSusy, mstop, mstopsq, mD1sq, mD2sq, s, tb);

  if (mstopsq(1) < 0.0 || mstopsq(2) < 0.0)
    {
      mstopsq(1) = Abs(mstopsq(1));
      mstopsq(2) = Abs(mstopsq(2));
    }

  // Have to swap the mass ordering of the stops, because
  // physical_ESSM uses the opposite labelling.
  double mstop1sq, mstop2sq;
  mstop1sq = mstopsq(2);
  mstop2sq = mstopsq(1);

  // Calculate the explicit input parameter dependence of the stop masses
  double dmstop1sqdmQlsq = 0.5*(1.0-MQQsq/Sqrt(rt));
  double dmstop2sqdmQlsq =0.5*(1.0+MQQsq/Sqrt(rt));

  // Calculate the derivative of the stop masses wrt the appropriate VEV
  double dmstop1sqdv3 = doCalcdmstop1sqdv3(essmSusy, s, tb);
  double dmstop2sqdv3 = doCalcdmstop2sqdv3(essmSusy, s, tb);

  // Calculate the explicit input parameter dependence of the derivatives
  // of the stop masses wrt the appropriate VEV
  double coeff = mtop*mtop*Xt*lambda*v1*MQQsq/(Sqrt(2.0*rt)*rt*v2*s);

  double d2mstop1sqdmQlsqdv3 = -2.0*s*coeff;
  double d2mstop2sqdmQlsqdv3 = 2.0*s*coeff;

  // Calculate the full second derivative of the tadpole contribution
  double term1 = dmstop1sqdmQlsq*dmstop1sqdv3*Log(mstop1sq/(q*q))+a0Peter(mstop1sq,q)*d2mstop1sqdmQlsqdv3;
  double term2 = dmstop2sqdmQlsq*dmstop2sqdv3*Log(mstop2sq/(q*q))+a0Peter(mstop2sq,q)*d2mstop2sqdmQlsqdv3;

  double d2DeltaV = (3.0/(16.0*PI*PI))*(term1+term2);
  d2DeltaV = (1.0/s)*d2DeltaV;

  return d2DeltaV;

}

double doCalcd2DeltaVdmUrsqdv1(genericE6SSM_soft_parameters essmSusy, double s, double tb)
{
  double v1 = essmSusy.get_vd(); 
  double v2 = essmSusy.get_vu(); 
  double v = Sqrt(v1*v1+v2*v2);

  double g1 = essmSusy.get_g1();
  double g2 = essmSusy.get_g2();
  double gbar = Sqrt(g2*g2+0.6*g1*g1);
  double gd1 = essmSusy.get_gN();

  // Logarithms are calculated at the given Q
  double q = essmSusy.get_scale();

  double yt = essmSusy.get_Yu( 2, 2);
  double mtop = yt*v2/Sqrt(2.0);
  // Option to fix m_t(M_t) = 165 GeV
  
  if(USEMTOFMT){
    mtop = 165;
    yt = Sqrt(2.0)*mtop/v2;
    essmSusy.set_Yu(2,2, yt);
  }

  
  double lambda = essmSusy.get_Lambdax();
  
  double At;
  double TYt = essmSusy.get_TYu(2,2);
  
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

  double mQlsq = essmSusy.get_mq2(2,2);
  double mUrsq =  essmSusy.get_mu2(2,2);
  double mueff = lambda*s/Sqrt(2.0);

  double MQQsq = mQlsq-mUrsq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2);
  double Xt = At - lambda*s/(Sqrt(2.0)*tb);
  double rt = MQQsq*MQQsq+4.0*mtop*mtop*Xt*Xt;

  // Squared stop masses are calculated using physical_ESSM to make sure
  // auxiliary D-terms are included.
  DoubleVector mstop(2), mstopsq(2), mD1sq(3), mD2sq(3);
  physical_ESSM(essmSusy, mstop, mstopsq, mD1sq, mD2sq, s, tb);

  if (mstopsq(1) < 0.0 || mstopsq(2) < 0.0)
    {
      mstopsq(1) = Abs(mstopsq(1));
      mstopsq(2) = Abs(mstopsq(2));
    }

  // Have to swap the mass ordering of the stops, because
  // physical_ESSM uses the opposite labelling.
  double mstop1sq, mstop2sq;
  mstop1sq = mstopsq(2);
  mstop2sq = mstopsq(1);

  // Calculate the explicit input parameter dependence of the stop masses
  double dmstop1sqdmUrsq = 0.5*(1.0+MQQsq/Sqrt(rt));
  double dmstop2sqdmUrsq =0.5*(1.0-MQQsq/Sqrt(rt));

  // Calculate the derivative of the stop masses wrt the appropriate VEV
  double dmstop1sqdv1 = doCalcdmstop1sqdv1(essmSusy, s, tb);
  double dmstop2sqdv1 = doCalcdmstop2sqdv1(essmSusy, s, tb);

  // Calculate the explicit input parameter dependence of the derivatives
  // of the stop masses wrt the appropriate VEV
  double coeff = 0.25*(g2*g2-g1*g1)/Sqrt(rt);
  coeff = coeff+(MQQsq/(rt*Sqrt(rt)))*(0.25*MQQsq*(g1*g1-g2*g2)
				       +4.0*mtop*mtop*Xt*lambda*s/(Sqrt(2.0)*v1*v2));

  double d2mstop1sqdmUrsqdv1 = 0.5*v1*coeff;
  double d2mstop2sqdmUrsqdv1 = -0.5*v1*coeff;

  // Calculate the full second derivative of the tadpole contribution
  double term1 = dmstop1sqdmUrsq*dmstop1sqdv1*Log(mstop1sq/(q*q))+a0Peter(mstop1sq,q)*d2mstop1sqdmUrsqdv1;
  double term2 = dmstop2sqdmUrsq*dmstop2sqdv1*Log(mstop2sq/(q*q))+a0Peter(mstop2sq,q)*d2mstop2sqdmUrsqdv1;

  double d2DeltaV = (3.0/(16.0*PI*PI))*(term1+term2);
  d2DeltaV = (1.0/v1)*d2DeltaV;

  return d2DeltaV;

}

double doCalcd2DeltaVdmUrsqdv2(genericE6SSM_soft_parameters essmSusy, double s, double tb)
{
  double v1 = essmSusy.get_vd(); 
  double v2 = essmSusy.get_vu(); 
  double v = Sqrt(v1*v1+v2*v2);

  double g1 = essmSusy.get_g1();
  double g2 = essmSusy.get_g2();
  double gbar = Sqrt(g2*g2+0.6*g1*g1);
  double gd1 = essmSusy.get_gN();

  // Logarithms are calculated at the given Q
  double q = essmSusy.get_scale();

  double yt = essmSusy.get_Yu( 2, 2);
  double mtop = yt*v2/Sqrt(2.0);
  // Option to fix m_t(M_t) = 165 GeV
  
  if(USEMTOFMT){
    mtop = 165;
    yt = Sqrt(2.0)*mtop/v2;
    essmSusy.set_Yu(2,2, yt);
  }

  
  double lambda = essmSusy.get_Lambdax();
  
  double At;
  double TYt = essmSusy.get_TYu(2,2);
  
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

  double mQlsq = essmSusy.get_mq2(2,2);
  double mUrsq =  essmSusy.get_mu2(2,2);
  double mueff = lambda*s/Sqrt(2.0);

  double MQQsq = mQlsq-mUrsq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2);
  double Xt = At - lambda*s/(Sqrt(2.0)*tb);
  double rt = MQQsq*MQQsq+4.0*mtop*mtop*Xt*Xt;

  // Squared stop masses are calculated using physical_ESSM to make sure
  // auxiliary D-terms are included.
  DoubleVector mstop(2), mstopsq(2), mD1sq(3), mD2sq(3);
  physical_ESSM(essmSusy, mstop, mstopsq, mD1sq, mD2sq, s, tb);

  if (mstopsq(1) < 0.0 || mstopsq(2) < 0.0)
    {
      mstopsq(1) = Abs(mstopsq(1));
      mstopsq(2) = Abs(mstopsq(2));
    }

  // Have to swap the mass ordering of the stops, because
  // physical_ESSM uses the opposite labelling.
  double mstop1sq, mstop2sq;
  mstop1sq = mstopsq(2);
  mstop2sq = mstopsq(1);

  // Calculate the explicit input parameter dependence of the stop masses
  double dmstop1sqdmUrsq = 0.5*(1.0+MQQsq/Sqrt(rt));
  double dmstop2sqdmUrsq =0.5*(1.0-MQQsq/Sqrt(rt));

  // Calculate the derivative of the stop masses wrt the appropriate VEV
  double dmstop1sqdv2 = doCalcdmstop1sqdv2(essmSusy, s, tb);
  double dmstop2sqdv2 = doCalcdmstop2sqdv2(essmSusy, s, tb);

  // Calculate the explicit input parameter dependence of the derivatives
  // of the stop masses wrt the appropriate VEV
  double coeff = 0.25*(g1*g1-g2*g2)/Sqrt(rt);
  coeff = coeff+(MQQsq/(rt*Sqrt(rt)))*(0.25*MQQsq*(g2*g2-g1*g1)-2.0*yt*yt*Xt*Xt
				       -4.0*mtop*mtop*Xt*lambda*s*v1/(Sqrt(2.0)*v2*v2*v2));

  double d2mstop1sqdmUrsqdv2 = 0.5*v2*coeff;
  double d2mstop2sqdmUrsqdv2 = -0.5*v2*coeff;

  // Calculate the full second derivative of the tadpole contribution
  double term1 = dmstop1sqdmUrsq*dmstop1sqdv2*Log(mstop1sq/(q*q))+a0Peter(mstop1sq,q)*d2mstop1sqdmUrsqdv2;
  double term2 = dmstop2sqdmUrsq*dmstop2sqdv2*Log(mstop2sq/(q*q))+a0Peter(mstop2sq,q)*d2mstop2sqdmUrsqdv2;

  double d2DeltaV = (3.0/(16.0*PI*PI))*(term1+term2);
  d2DeltaV = (1.0/v2)*d2DeltaV;

  return d2DeltaV;

}

double doCalcd2DeltaVdmUrsqdv3(genericE6SSM_soft_parameters essmSusy, double s, double tb)
{
  double v1 = essmSusy.get_vd(); 
  double v2 = essmSusy.get_vu(); 
  double v = Sqrt(v1*v1+v2*v2);

  double g1 = essmSusy.get_g1();
  double g2 = essmSusy.get_g2();
  double gbar = Sqrt(g2*g2+0.6*g1*g1);
  double gd1 = essmSusy.get_gN();

  // Logarithms are calculated at the given Q
  double q = essmSusy.get_scale();

  double yt = essmSusy.get_Yu(2, 2);
  double mtop = yt*v2/Sqrt(2.0);
  // Option to fix m_t(M_t) = 165 GeV
  
  if(USEMTOFMT){
    mtop = 165;
    yt = Sqrt(2.0)*mtop/v2;
    essmSusy.set_Yu(2,2, yt);
  }

  
  double lambda = essmSusy.get_Lambdax();
 
  double At;
  double TYt = essmSusy.get_TYu(2,2);
  
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

  double mQlsq = essmSusy.get_mq2(2,2);
  double mUrsq =  essmSusy.get_mu2(2,2);
  double mueff = lambda*s/Sqrt(2.0);

  double MQQsq = mQlsq-mUrsq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2);
  double Xt = At - lambda*s/(Sqrt(2.0)*tb);
  double rt = MQQsq*MQQsq+4.0*mtop*mtop*Xt*Xt;

  // Squared stop masses are calculated using physical_ESSM to make sure
  // auxiliary D-terms are included.
  DoubleVector mstop(2), mstopsq(2), mD1sq(3), mD2sq(3);
  physical_ESSM(essmSusy, mstop, mstopsq, mD1sq, mD2sq, s, tb);

  if (mstopsq(1) < 0.0 || mstopsq(2) < 0.0)
    {
      mstopsq(1) = Abs(mstopsq(1));
      mstopsq(2) = Abs(mstopsq(2));
    }

  // Have to swap the mass ordering of the stops, because
  // physical_ESSM uses the opposite labelling.
  double mstop1sq, mstop2sq;
  mstop1sq = mstopsq(2);
  mstop2sq = mstopsq(1);

  // Calculate the explicit input parameter dependence of the stop masses
  double dmstop1sqdmUrsq = 0.5*(1.0+MQQsq/Sqrt(rt));
  double dmstop2sqdmUrsq =0.5*(1.0-MQQsq/Sqrt(rt));

  // Calculate the derivative of the stop masses wrt the appropriate VEV
  double dmstop1sqdv3 = doCalcdmstop1sqdv3(essmSusy, s, tb);
  double dmstop2sqdv3 = doCalcdmstop2sqdv3(essmSusy, s, tb);

  // Calculate the explicit input parameter dependence of the derivatives
  // of the stop masses wrt the appropriate VEV
  double coeff = mtop*mtop*Xt*lambda*v1*MQQsq/(Sqrt(2.0*rt)*rt*v2*s);

  double d2mstop1sqdmUrsqdv3 = 2.0*s*coeff;
  double d2mstop2sqdmUrsqdv3 = -2.0*s*coeff;

  // Calculate the full second derivative of the tadpole contribution
  double term1 = dmstop1sqdmUrsq*dmstop1sqdv3*Log(mstop1sq/(q*q))+a0Peter(mstop1sq,q)*d2mstop1sqdmUrsqdv3;
  double term2 = dmstop2sqdmUrsq*dmstop2sqdv3*Log(mstop2sq/(q*q))+a0Peter(mstop2sq,q)*d2mstop2sqdmUrsqdv3;

  double d2DeltaV = (3.0/(16.0*PI*PI))*(term1+term2);
  d2DeltaV = (1.0/s)*d2DeltaV;

  return d2DeltaV;

}

double doCalcMh2SquaredLogCoeff(genericE6SSM_soft_parameters r, int nLps)
{

  if (nLps == 1 || nLps == 2)
    {
      r.set_loops(nLps);    
    }
  else
    {
      cerr << "WARNING: can only do 1 or 2 loops: returning nLps = 1 loop result." << endl;
      r.set_loops(1);      
    }

  genericE6SSM_soft_parameters w = r.calc_beta();

  return w.get_mHu2();

}

double doCalcMh2SquaredLogSqCoeff(genericE6SSM_soft_parameters r, int nLogs)
{
  // Assumes first and second generation Yukawas vanish, 
  // consistent with the previous pMSSM study
  r.set_Yu(0,0,0.0);
  r.set_Yu(1,1,0.0);

  r.set_Yd(0,0,0.0);
  r.set_Ye(1,1,0.0);

  r.set_Ye(0,0,0.0);
  r.set_Ye(1,1,0.0);

  // For consistency the correspond A terms should also vanish
  r.set_TYu(0,0,0.0);
  r.set_TYu(1,1,0.0);

  r.set_TYd(0,0,0.0);
  r.set_TYd(1,1,0.0);

  r.set_TYe(0,0,0.0);
  r.set_TYe(1,1,0.0);

  if (nLogs == 0)
    {
      return 0.0;
    }

  if (nLogs != 1)
    {
      cerr << "WARNING: cannot yet do more than LL calculation: using nLogs = 1." << endl;
    }

  r.set_loops(1);

  double QQp = r.get_input().QQp;
  double Qup = r.get_input().Qup;
  double Qdp = r.get_input().Qdp;
  double QLp = r.get_input().QLp;
  double Qep = r.get_input().Qep;
  double QH1p = r.get_input().QH1p;
  double QH2p = r.get_input().QH2p;
  double QSp = r.get_input().QSp;
  double QHpp = r.get_input().QHpp;
  double QHpbarp = r.get_input().QHpbarp;
  double QDxp = r.get_input().QDxp;
  double QDxbarp = r.get_input().QDxbarp;

  genericE6SSM_soft_parameters w = r.calc_beta();

  double mHu2, mHd2, ms2, mHp2, mHpbar2;
  Eigen::Matrix<double,3,3> mq2, ml2, mu2, md2, me2, mDx2, mDxbar2;
  Eigen::Matrix<double,2,2> mH1I2, mH2I2, msI2;

  double betamHu2, betamHd2, betams2, betamHp2, betamHpbar2;
  Eigen::Matrix<double,3,3> betamq2, betaml2, betamu2, betamd2, betame2, betamDx2, betamDxbar2;
  Eigen::Matrix<double,2,2> betamH1I2, betamH2I2, betamsI2;

  mHu2 = r.get_mHu2(); betamHu2 = w.get_mHu2();
  mHd2 = r.get_mHd2(); betamHd2 = w.get_mHd2();
  ms2 = r.get_ms2(); betams2 = w.get_ms2();
  mHp2 = r.get_mHp2(); betamHp2 = w.get_mHp2();
  mHpbar2 = r.get_mHpbar2(); betamHpbar2 = w.get_mHpbar2();

  mq2 = r.get_mq2(); betamq2 = w.get_mq2();
  ml2 = r.get_ml2(); betaml2 = w.get_ml2();
  mu2 = r.get_mu2(); betamu2 = w.get_mu2();
  md2 = r.get_md2(); betamd2 = w.get_md2();
  me2 = r.get_me2(); betame2 = w.get_me2();
  mDx2 = r.get_mDx2(); betamDx2 = w.get_mDx2();
  mDxbar2 = r.get_mDxbar2(); betamDxbar2 = w.get_mDxbar2();

  mH1I2 = r.get_mH1I2(); betamH1I2 = w.get_mH1I2();
  mH2I2 = r.get_mH2I2(); betamH2I2 = w.get_mH2I2();
  msI2 = r.get_msI2(); betamsI2 = w.get_msI2();

  double M1, M2, M1p;
  double betaM1, betaM2, betaM1p;

  M1 = r.get_MassB(); betaM1 = w.get_MassB();
  M2 = r.get_MassWB(); betaM2 = w.get_MassWB();
  M1p = r.get_MassBp(); betaM1p = w.get_MassBp();

  double g1, g2, g3, gN;
  double betag1, betag2, betag3, betagN;
  g1 = r.get_g1(); betag1 = w.get_g1();
  g2 = r.get_g2(); betag2 = w.get_g2();
  g3 = r.get_g3(); betag3 = w.get_g3();
  gN = r.get_gN(); betagN = w.get_gN();

  double yt;
  double betayt;
  yt = r.get_Yu(2,2); betayt = w.get_Yu(2,2);

  double Tyt;
  double betaTyt;
  Tyt = r.get_TYu(2,2); betaTyt = w.get_TYu(2,2);

  double lambda, betalambda, Tlambda, betaTlambda;
  lambda = r.get_Lambdax(); betalambda = w.get_Lambdax();
  Tlambda = r.get_TLambdax(); betaTlambda = w.get_TLambdax();

  double sigma11 = mHu2 - mHd2 + mHpbar2 - mHp2 + md2.trace()
    - mDx2.trace() + mDxbar2.trace() + me2.trace() - mH1I2.trace()
    + mH2I2.trace() - ml2.trace() + mq2.trace() - 2.0 * mu2.trace();

  double sigma14 = 2.0 * QH1p * mHd2 + 2.0 * QH2p * mHu2 + 2.0 * QHpbarp * mHpbar2
    + 2.0 * QHpp * mHp2 + QSp * ms2 + 3.0 * Qdp * md2.trace() + 3.0 * QDxp * mDx2.trace()
    + 3.0 * QDxbarp * mDxbar2.trace() + Qep * me2.trace() + 2.0 * QH1p * mH1I2.trace()
    + 2.0 * QH2p * mH2I2.trace() + 2.0 * QLp * ml2.trace() + 6.0 * QQp * mq2.trace() + QSp * msI2.trace()
    + 3.0 * Qup * mu2.trace();

  double betasigma11 = betamHu2 - betamHd2 + betamHpbar2 - betamHp2 + betamd2.trace()
    - betamDx2.trace() + betamDxbar2.trace() + betame2.trace() - betamH1I2.trace()
    + betamH2I2.trace() - betaml2.trace() + betamq2.trace() - 2.0 * betamu2.trace();

  double betasigma14 = 2.0 * QH1p * betamHd2 + 2.0 * QH2p * betamHu2 + 2.0 * QHpbarp * betamHpbar2
    + 2.0 * QHpp * betamHp2 + QSp * betams2 + 3.0 * Qdp * betamd2.trace() + 3.0 * QDxp * betamDx2.trace()
    + 3.0 * QDxbarp * betamDxbar2.trace() + Qep * betame2.trace() + 2.0 * QH1p * betamH1I2.trace()
    + 2.0 * QH2p * betamH2I2.trace() + 2.0 * QLp * betaml2.trace() + 6.0 * QQp * betamq2.trace() + QSp * betamsI2.trace()
    + 3.0 * Qup * betamu2.trace();

  double coeff = -2.4 * (g1 * betag1 * Sqr(M1) + Sqr(g1) * M1 * betaM1);
  coeff -= 16.0 * Sqr(QH2p) * (gN * betagN * Sqr(M1p) + Sqr(gN) * M1p * betaM1p);
  coeff -= 12.0 * (g2 * betag2 * Sqr(M2) + Sqr(g2) * M2 * betaM2);
  coeff += 4.0 * lambda * betalambda * (mHd2 + mHu2 + ms2);
  coeff += 2.0 * Sqr(lambda) * (betamHd2 + betamHu2 + betams2) + 4.0 * Tlambda * betaTlambda;
  coeff += 12.0 * yt * betayt * (mHu2 + mq2(2,2) + mu2(2,2));
  coeff += 6.0 * Sqr(yt) * (betamHu2 + betamq2(2,2) + betamu2(2,2)) + 12.0 * Tyt * betaTyt;
  coeff += 1.2 * g1 * betag1 * sigma11 + 0.6 * Sqr(g1) * betasigma11;
  coeff += 4.0 * QH2p * gN * betagN * sigma14 + 2.0 * QH2p * Sqr(gN) * betasigma14;

  coeff = 0.5*oneOver16PiSqr*coeff;

  return coeff;
}

double doCalcMh1SquaredLogCoeff(genericE6SSM_soft_parameters r, int nLps)
{
  if (nLps == 1 || nLps == 2)
    {
      r.set_loops(nLps);    
    }
  else
    {
      cerr << "WARNING: can only do 1 or 2 loops: returning nLps = 1 loop result." << endl;
      r.set_loops(1);      
    }

  genericE6SSM_soft_parameters w = r.calc_beta();

  return w.get_mHd2();
}

double doCalcMh1SquaredLogSqCoeff(genericE6SSM_soft_parameters r, int nLogs)
{
  // Assumes first and second generation Yukawas vanish, 
  // consistent with the previous pMSSM study
  r.set_Yu(0,0,0.0);
  r.set_Yu(1,1,0.0);

  r.set_Yd(0,0,0.0);
  r.set_Ye(1,1,0.0);

  r.set_Ye(0,0,0.0);
  r.set_Ye(1,1,0.0);

  // For consistency the correspond A terms should also vanish
  r.set_TYu(0,0,0.0);
  r.set_TYu(1,1,0.0);

  r.set_TYd(0,0,0.0);
  r.set_TYd(1,1,0.0);

  r.set_TYe(0,0,0.0);
  r.set_TYe(1,1,0.0);

  if (nLogs == 0)
    {
      return 0.0;
    }

  if (nLogs != 1)
    {
      cerr << "WARNING: cannot yet do more than LL calculation: using nLogs = 1." << endl;
    }

  r.set_loops(1);

  double QQp = r.get_input().QQp;
  double Qup = r.get_input().Qup;
  double Qdp = r.get_input().Qdp;
  double QLp = r.get_input().QLp;
  double Qep = r.get_input().Qep;
  double QH1p = r.get_input().QH1p;
  double QH2p = r.get_input().QH2p;
  double QSp = r.get_input().QSp;
  double QHpp = r.get_input().QHpp;
  double QHpbarp = r.get_input().QHpbarp;
  double QDxp = r.get_input().QDxp;
  double QDxbarp = r.get_input().QDxbarp;

  genericE6SSM_soft_parameters w = r.calc_beta();

  double mHu2, mHd2, ms2, mHp2, mHpbar2;
  Eigen::Matrix<double,3,3> mq2, ml2, mu2, md2, me2, mDx2, mDxbar2;
  Eigen::Matrix<double,2,2> mH1I2, mH2I2, msI2;

  double betamHu2, betamHd2, betams2, betamHp2, betamHpbar2;
  Eigen::Matrix<double,3,3> betamq2, betaml2, betamu2, betamd2, betame2, betamDx2, betamDxbar2;
  Eigen::Matrix<double,2,2> betamH1I2, betamH2I2, betamsI2;

  mHu2 = r.get_mHu2(); betamHu2 = w.get_mHu2();
  mHd2 = r.get_mHd2(); betamHd2 = w.get_mHd2();
  ms2 = r.get_ms2(); betams2 = w.get_ms2();
  mHp2 = r.get_mHp2(); betamHp2 = w.get_mHp2();
  mHpbar2 = r.get_mHpbar2(); betamHpbar2 = w.get_mHpbar2();

  mq2 = r.get_mq2(); betamq2 = w.get_mq2();
  ml2 = r.get_ml2(); betaml2 = w.get_ml2();
  mu2 = r.get_mu2(); betamu2 = w.get_mu2();
  md2 = r.get_md2(); betamd2 = w.get_md2();
  me2 = r.get_me2(); betame2 = w.get_me2();
  mDx2 = r.get_mDx2(); betamDx2 = w.get_mDx2();
  mDxbar2 = r.get_mDxbar2(); betamDxbar2 = w.get_mDxbar2();

  mH1I2 = r.get_mH1I2(); betamH1I2 = w.get_mH1I2();
  mH2I2 = r.get_mH2I2(); betamH2I2 = w.get_mH2I2();
  msI2 = r.get_msI2(); betamsI2 = w.get_msI2();

  double M1, M2, M1p;
  double betaM1, betaM2, betaM1p;

  M1 = r.get_MassB(); betaM1 = w.get_MassB();
  M2 = r.get_MassWB(); betaM2 = w.get_MassWB();
  M1p = r.get_MassBp(); betaM1p = w.get_MassBp();

  double g1, g2, g3, gN;
  double betag1, betag2, betag3, betagN;
  g1 = r.get_g1(); betag1 = w.get_g1();
  g2 = r.get_g2(); betag2 = w.get_g2();
  g3 = r.get_g3(); betag3 = w.get_g3();
  gN = r.get_gN(); betagN = w.get_gN();

  double yb, ytau;
  double betayb, betaytau;
  yb = r.get_Yd(2,2); betayb = w.get_Yd(2,2);
  ytau = r.get_Ye(2,2); betaytau = w.get_Ye(2,2);

  double Tyb, Tytau;
  double betaTyb, betaTytau;
  Tyb = r.get_TYd(2,2); betaTyb = w.get_TYd(2,2);
  Tytau = r.get_TYe(2,2); betaTytau = w.get_TYe(2,2);

  double lambda, betalambda, Tlambda, betaTlambda;
  lambda = r.get_Lambdax(); betalambda = w.get_Lambdax();
  Tlambda = r.get_TLambdax(); betaTlambda = w.get_TLambdax();

  double sigma11 = mHu2 - mHd2 + mHpbar2 - mHp2 + md2.trace()
    - mDx2.trace() + mDxbar2.trace() + me2.trace() - mH1I2.trace()
    + mH2I2.trace() - ml2.trace() + mq2.trace() - 2.0 * mu2.trace();

  double sigma14 = 2.0 * QH1p * mHd2 + 2.0 * QH2p * mHu2 + 2.0 * QHpbarp * mHpbar2
    + 2.0 * QHpp * mHp2 + QSp * ms2 + 3.0 * Qdp * md2.trace() + 3.0 * QDxp * mDx2.trace()
    + 3.0 * QDxbarp * mDxbar2.trace() + Qep * me2.trace() + 2.0 * QH1p * mH1I2.trace()
    + 2.0 * QH2p * mH2I2.trace() + 2.0 * QLp * ml2.trace() + 6.0 * QQp * mq2.trace() + QSp * msI2.trace()
    + 3.0 * Qup * mu2.trace();

  double betasigma11 = betamHu2 - betamHd2 + betamHpbar2 - betamHp2 + betamd2.trace()
    - betamDx2.trace() + betamDxbar2.trace() + betame2.trace() - betamH1I2.trace()
    + betamH2I2.trace() - betaml2.trace() + betamq2.trace() - 2.0 * betamu2.trace();

  double betasigma14 = 2.0 * QH1p * betamHd2 + 2.0 * QH2p * betamHu2 + 2.0 * QHpbarp * betamHpbar2
    + 2.0 * QHpp * betamHp2 + QSp * betams2 + 3.0 * Qdp * betamd2.trace() + 3.0 * QDxp * betamDx2.trace()
    + 3.0 * QDxbarp * betamDxbar2.trace() + Qep * betame2.trace() + 2.0 * QH1p * betamH1I2.trace()
    + 2.0 * QH2p * betamH2I2.trace() + 2.0 * QLp * betaml2.trace() + 6.0 * QQp * betamq2.trace() + QSp * betamsI2.trace()
    + 3.0 * Qup * betamu2.trace();

  double coeff = -2.4 * (g1 * betag1 * Sqr(M1) + Sqr(g1) * M1 * betaM1);
  coeff -= 16.0 * Sqr(QH1p) * (gN * betagN * Sqr(M1p) + Sqr(gN) * M1p * betaM1p);
  coeff -= 12.0 * (g2 * betag2 * Sqr(M2) + Sqr(g2) * M2 * betaM2);
  coeff += 4.0 * lambda * betalambda * (mHd2 + mHu2 + ms2);
  coeff += 2.0 * Sqr(lambda) * (betamHd2 + betamHu2 + betams2) + 4.0 * Tlambda * betaTlambda;
  coeff += 12.0 * yb * betayb * (mHd2 + mq2(2,2) + md2(2,2));
  coeff += 6.0 * Sqr(yb) * (betamHd2 + betamq2(2,2) + betamd2(2,2)) + 12.0 * Tyb * betaTyb;
  coeff += 4.0 * ytau * betaytau * (mHd2 + ml2(2,2) + me2(2,2));
  coeff += 2.0 * Sqr(ytau) * (betamHd2 + betaml2(2,2) + betame2(2,2)) + 4.0 * Tytau * betaTytau;
  coeff -= (1.2 * g1 * betag1 * sigma11 + 0.6 * Sqr(g1) * betasigma11);
  coeff += 4.0 * QH1p * gN * betagN * sigma14 + 2.0 * QH1p * Sqr(gN) * betasigma14;

  coeff = 0.5*oneOver16PiSqr*coeff;

  return coeff;
}

double doCalcMsSquaredLogCoeff(genericE6SSM_soft_parameters r, int nLps)
{
  if (nLps == 1 || nLps == 2)
    {
      r.set_loops(nLps);    
    }
  else
    {
      cerr << "WARNING: can only do 1 or 2 loops: returning nLps = 1 loop result." << endl;
      r.set_loops(1);      
    }

  genericE6SSM_soft_parameters w = r.calc_beta();

  return w.get_ms2();
}

double doCalcMsSquaredLogSqCoeff(genericE6SSM_soft_parameters r, int nLogs)
{
  // Assumes first and second generation Yukawas vanish, 
  // consistent with the previous pMSSM study
  r.set_Yu(0,0,0.0);
  r.set_Yu(1,1,0.0);

  r.set_Yd(0,0,0.0);
  r.set_Ye(1,1,0.0);

  r.set_Ye(0,0,0.0);
  r.set_Ye(1,1,0.0);

  // For consistency the corresponding A terms should also vanish
  r.set_TYu(0,0,0.0);
  r.set_TYu(1,1,0.0);

  r.set_TYd(0,0,0.0);
  r.set_TYd(1,1,0.0);

  r.set_TYe(0,0,0.0);
  r.set_TYe(1,1,0.0);

  if (nLogs == 0)
    {
      return 0.0;
    }

  if (nLogs != 1)
    {
      cerr << "WARNING: cannot yet do more than LL calculation: using nLogs = 1." << endl;
    }

  r.set_loops(1);

  double QQp = r.get_input().QQp;
  double Qup = r.get_input().Qup;
  double Qdp = r.get_input().Qdp;
  double QLp = r.get_input().QLp;
  double Qep = r.get_input().Qep;
  double QH1p = r.get_input().QH1p;
  double QH2p = r.get_input().QH2p;
  double QSp = r.get_input().QSp;
  double QHpp = r.get_input().QHpp;
  double QHpbarp = r.get_input().QHpbarp;
  double QDxp = r.get_input().QDxp;
  double QDxbarp = r.get_input().QDxbarp;

  genericE6SSM_soft_parameters w = r.calc_beta();

  double mHu2, mHd2, ms2, mHp2, mHpbar2;
  Eigen::Matrix<double,3,3> mq2, ml2, mu2, md2, me2, mDx2, mDxbar2;
  Eigen::Matrix<double,2,2> mH1I2, mH2I2, msI2;

  double betamHu2, betamHd2, betams2, betamHp2, betamHpbar2;
  Eigen::Matrix<double,3,3> betamq2, betaml2, betamu2, betamd2, betame2, betamDx2, betamDxbar2;
  Eigen::Matrix<double,2,2> betamH1I2, betamH2I2, betamsI2;

  mHu2 = r.get_mHu2(); betamHu2 = w.get_mHu2();
  mHd2 = r.get_mHd2(); betamHd2 = w.get_mHd2();
  ms2 = r.get_ms2(); betams2 = w.get_ms2();
  mHp2 = r.get_mHp2(); betamHp2 = w.get_mHp2();
  mHpbar2 = r.get_mHpbar2(); betamHpbar2 = w.get_mHpbar2();

  mq2 = r.get_mq2(); betamq2 = w.get_mq2();
  ml2 = r.get_ml2(); betaml2 = w.get_ml2();
  mu2 = r.get_mu2(); betamu2 = w.get_mu2();
  md2 = r.get_md2(); betamd2 = w.get_md2();
  me2 = r.get_me2(); betame2 = w.get_me2();
  mDx2 = r.get_mDx2(); betamDx2 = w.get_mDx2();
  mDxbar2 = r.get_mDxbar2(); betamDxbar2 = w.get_mDxbar2();

  mH1I2 = r.get_mH1I2(); betamH1I2 = w.get_mH1I2();
  mH2I2 = r.get_mH2I2(); betamH2I2 = w.get_mH2I2();
  msI2 = r.get_msI2(); betamsI2 = w.get_msI2();

  double M1, M1p;
  double betaM1, betaM1p;

  M1 = r.get_MassB(); betaM1 = w.get_MassB();
  M1p = r.get_MassBp(); betaM1p = w.get_MassBp();

  double g1, g2, g3, gN;
  double betag1, betag2, betag3, betagN;
  g1 = r.get_g1(); betag1 = w.get_g1();
  g2 = r.get_g2(); betag2 = w.get_g2();
  g3 = r.get_g3(); betag3 = w.get_g3();
  gN = r.get_gN(); betagN = w.get_gN();

  double lambda, betalambda, Tlambda, betaTlambda;

  lambda = r.get_Lambdax(); betalambda = w.get_Lambdax();
  Tlambda = r.get_TLambdax(); betaTlambda = w.get_TLambdax();

  Eigen::Matrix<double,3,3> Kappa, KappaAdj, betaKappa, betaKappaAdj;
  Eigen::Matrix<double,3,3> TKappa, TKappaAdj, betaTKappa, betaTKappaAdj;
  Eigen::Matrix<double,2,2> Lambda12, Lambda12Adj, betaLambda12, betaLambda12Adj;
  Eigen::Matrix<double,2,2> TLambda12, TLambda12Adj, betaTLambda12, betaTLambda12Adj;

  Kappa = r.get_Kappa(); KappaAdj = Kappa.adjoint();
  betaKappa = w.get_Kappa(); betaKappaAdj = betaKappa.adjoint();

  TKappa = r.get_TKappa(); TKappaAdj = TKappa.adjoint();
  betaTKappa = w.get_TKappa(); betaTKappaAdj = betaTKappa.adjoint();

  Lambda12 = r.get_Lambda12(); Lambda12Adj = Lambda12.adjoint();
  betaLambda12 = w.get_Lambda12(); betaLambda12Adj = betaLambda12.adjoint();

  TLambda12 = r.get_TLambda12(); TLambda12Adj = TLambda12.adjoint();
  betaTLambda12 = w.get_TLambda12(); betaTLambda12Adj = betaTLambda12.adjoint();

  double sigma14 = 2.0 * QH1p * mHd2 + 2.0 * QH2p * mHu2 + 2.0 * QHpbarp * mHpbar2
    + 2.0 * QHpp * mHp2 + QSp * ms2 + 3.0 * Qdp * md2.trace() + 3.0 * QDxp * mDx2.trace()
    + 3.0 * QDxbarp * mDxbar2.trace() + Qep * me2.trace() + 2.0 * QH1p * mH1I2.trace()
    + 2.0 * QH2p * mH2I2.trace() + 2.0 * QLp * ml2.trace() + 6.0 * QQp * mq2.trace() + QSp * msI2.trace()
    + 3.0 * Qup * mu2.trace();

  double betasigma14 = 2.0 * QH1p * betamHd2 + 2.0 * QH2p * betamHu2 + 2.0 * QHpbarp * betamHpbar2
    + 2.0 * QHpp * betamHp2 + QSp * betams2 + 3.0 * Qdp * betamd2.trace() + 3.0 * QDxp * betamDx2.trace()
    + 3.0 * QDxbarp * betamDxbar2.trace() + Qep * betame2.trace() + 2.0 * QH1p * betamH1I2.trace()
    + 2.0 * QH2p * betamH2I2.trace() + 2.0 * QLp * betaml2.trace() + 6.0 * QQp * betamq2.trace() + QSp * betamsI2.trace()
    + 3.0 * Qup * betamu2.trace();

  double coeff = -16.0 * Sqr(QSp) * (gN * betagN * Sqr(M1p) + Sqr(gN) * M1p * betaM1p);
  coeff += 8.0 * lambda * betalambda * (mHd2 + mHu2 + ms2);
  coeff += 4.0 * Sqr(lambda) * (betamHd2 + betamHu2 + betams2) + 8.0 * Tlambda * betaTlambda;
  coeff += 6.0 * betams2 * (Kappa * KappaAdj).trace() + 6.0 * ms2 * (betaKappa * KappaAdj).trace() + 6.0 * ms2 * (Kappa * betaKappaAdj).trace();
  coeff += 4.0 * betams2 * (Lambda12 * Lambda12Adj).trace() + 4.0 * ms2 * (betaLambda12 * Lambda12Adj).trace() + 4.0 * ms2 * (Lambda12 * betaLambda12Adj).trace();
  coeff += 6.0 * (betaTKappa * TKappaAdj).trace() + 6.0 * (TKappa * betaTKappaAdj).trace();
  coeff += 4.0 * (betaTLambda12 * TLambda12Adj).trace() + 4.0 * (TLambda12 * betaTLambda12Adj).trace();
  coeff += 4.0 * (betamH1I2 * Lambda12Adj * Lambda12 + mH1I2 * betaLambda12Adj * Lambda12 + mH1I2 * Lambda12Adj * betaLambda12).trace();
  coeff += 6.0 * (betaKappa * KappaAdj * mDx2 + Kappa * betaKappaAdj * mDx2 + Kappa * KappaAdj * betamDx2).trace();
  coeff += 6.0 * (betaKappa * mDxbar2 * KappaAdj + Kappa * betamDxbar2 * KappaAdj + Kappa * mDxbar2 * betaKappaAdj).trace();
  coeff += 4.0 * (betaLambda12 * Lambda12Adj * mH2I2 + Lambda12 * betaLambda12Adj * mH2I2 + Lambda12 * Lambda12Adj * betamH2I2).trace();
  coeff += 4.0 * QSp * gN * betagN * sigma14 + 2.0 * QSp * Sqr(gN) * betasigma14;

  coeff = 0.5*oneOver16PiSqr*coeff;

  return coeff;
}

double doCalcmtRSquaredLogCoeff(genericE6SSM_soft_parameters r, int nLps)
{
  if (nLps == 1 || nLps == 2)
    {
      r.set_loops(nLps);    
    }
  else
    {
      cerr << "WARNING: can only do 1 or 2 loops: returning nLps = 1 loop result." << endl;
      r.set_loops(1);      
    }

  genericE6SSM_soft_parameters w = r.calc_beta();

  return w.get_mu2(2,2);
}

double doCalcmtRSquaredLogSqCoeff(genericE6SSM_soft_parameters r, int nLogs)
{
  // Assumes first and second generation Yukawas vanish, 
  // consistent with the previous pMSSM study
  r.set_Yu(0,0,0.0);
  r.set_Yu(1,1,0.0);

  r.set_Yd(0,0,0.0);
  r.set_Ye(1,1,0.0);

  r.set_Ye(0,0,0.0);
  r.set_Ye(1,1,0.0);

  // For consistency the corresponding A terms should also vanish
  r.set_TYu(0,0,0.0);
  r.set_TYu(1,1,0.0);

  r.set_TYd(0,0,0.0);
  r.set_TYd(1,1,0.0);

  r.set_TYe(0,0,0.0);
  r.set_TYe(1,1,0.0);

  if (nLogs == 0)
    {
      return 0.0;
    }

  if (nLogs != 1)
    {
      cerr << "WARNING: cannot yet do more than LL calculation: using nLogs = 1." << endl;
    }

  r.set_loops(1);

  double QQp = r.get_input().QQp;
  double Qup = r.get_input().Qup;
  double Qdp = r.get_input().Qdp;
  double QLp = r.get_input().QLp;
  double Qep = r.get_input().Qep;
  double QH1p = r.get_input().QH1p;
  double QH2p = r.get_input().QH2p;
  double QSp = r.get_input().QSp;
  double QHpp = r.get_input().QHpp;
  double QHpbarp = r.get_input().QHpbarp;
  double QDxp = r.get_input().QDxp;
  double QDxbarp = r.get_input().QDxbarp;

  genericE6SSM_soft_parameters w = r.calc_beta();

  double mHu2, mHd2, ms2, mHp2, mHpbar2;
  Eigen::Matrix<double,3,3> mq2, ml2, mu2, md2, me2, mDx2, mDxbar2;
  Eigen::Matrix<double,2,2> mH1I2, mH2I2, msI2;

  double betamHu2, betamHd2, betams2, betamHp2, betamHpbar2;
  Eigen::Matrix<double,3,3> betamq2, betaml2, betamu2, betamd2, betame2, betamDx2, betamDxbar2;
  Eigen::Matrix<double,2,2> betamH1I2, betamH2I2, betamsI2;

  mHu2 = r.get_mHu2(); betamHu2 = w.get_mHu2();
  mHd2 = r.get_mHd2(); betamHd2 = w.get_mHd2();
  ms2 = r.get_ms2(); betams2 = w.get_ms2();
  mHp2 = r.get_mHp2(); betamHp2 = w.get_mHp2();
  mHpbar2 = r.get_mHpbar2(); betamHpbar2 = w.get_mHpbar2();

  mq2 = r.get_mq2(); betamq2 = w.get_mq2();
  ml2 = r.get_ml2(); betaml2 = w.get_ml2();
  mu2 = r.get_mu2(); betamu2 = w.get_mu2();
  md2 = r.get_md2(); betamd2 = w.get_md2();
  me2 = r.get_me2(); betame2 = w.get_me2();
  mDx2 = r.get_mDx2(); betamDx2 = w.get_mDx2();
  mDxbar2 = r.get_mDxbar2(); betamDxbar2 = w.get_mDxbar2();

  mH1I2 = r.get_mH1I2(); betamH1I2 = w.get_mH1I2();
  mH2I2 = r.get_mH2I2(); betamH2I2 = w.get_mH2I2();
  msI2 = r.get_msI2(); betamsI2 = w.get_msI2();

  double M1, M3, M1p;
  double betaM1, betaM3, betaM1p;

  M1 = r.get_MassB(); betaM1 = w.get_MassB();
  M3 = r.get_MassG(); betaM3 = w.get_MassG();
  M1p = r.get_MassBp(); betaM1p = w.get_MassBp();

  double g1, g2, g3, gN;
  double betag1, betag2, betag3, betagN;
  g1 = r.get_g1(); betag1 = w.get_g1();
  g3 = r.get_g3(); betag3 = w.get_g3();
  gN = r.get_gN(); betagN = w.get_gN();

  double yt;
  double betayt;
  yt = r.get_Yu(2,2); betayt = w.get_Yu(2,2);

  double Tyt;
  double betaTyt;
  Tyt = r.get_TYu(2,2); betaTyt = w.get_TYu(2,2);

  double sigma11 = mHu2 - mHd2 + mHpbar2 - mHp2 + md2.trace()
    - mDx2.trace() + mDxbar2.trace() + me2.trace() - mH1I2.trace()
    + mH2I2.trace() - ml2.trace() + mq2.trace() - 2.0 * mu2.trace();

  double sigma14 = 2.0 * QH1p * mHd2 + 2.0 * QH2p * mHu2 + 2.0 * QHpbarp * mHpbar2
    + 2.0 * QHpp * mHp2 + QSp * ms2 + 3.0 * Qdp * md2.trace() + 3.0 * QDxp * mDx2.trace()
    + 3.0 * QDxbarp * mDxbar2.trace() + Qep * me2.trace() + 2.0 * QH1p * mH1I2.trace()
    + 2.0 * QH2p * mH2I2.trace() + 2.0 * QLp * ml2.trace() + 6.0 * QQp * mq2.trace() + QSp * msI2.trace()
    + 3.0 * Qup * mu2.trace();

  double betasigma11 = betamHu2 - betamHd2 + betamHpbar2 - betamHp2 + betamd2.trace()
    - betamDx2.trace() + betamDxbar2.trace() + betame2.trace() - betamH1I2.trace()
    + betamH2I2.trace() - betaml2.trace() + betamq2.trace() - 2.0 * betamu2.trace();

  double betasigma14 = 2.0 * QH1p * betamHd2 + 2.0 * QH2p * betamHu2 + 2.0 * QHpbarp * betamHpbar2
    + 2.0 * QHpp * betamHp2 + QSp * betams2 + 3.0 * Qdp * betamd2.trace() + 3.0 * QDxp * betamDx2.trace()
    + 3.0 * QDxbarp * betamDxbar2.trace() + Qep * betame2.trace() + 2.0 * QH1p * betamH1I2.trace()
    + 2.0 * QH2p * betamH2I2.trace() + 2.0 * QLp * betaml2.trace() + 6.0 * QQp * betamq2.trace() + QSp * betamsI2.trace()
    + 3.0 * Qup * betamu2.trace();

  double coeff = -16.0 * Sqr(Qup) * (gN * betagN * Sqr(M1p) + Sqr(gN) * M1p * betaM1p);
  coeff -= (64.0/15.0) * (g1 * betag1 * Sqr(M1) + Sqr(g1) * M1 * betaM1);
  coeff -= (64.0/3.0) * (g3 * betag3 * Sqr(M3) + Sqr(g3) * M3 * betaM3);
  coeff += 8.0 * yt * betayt * (mHu2 + mu2(2,2) + mq2(2,2));
  coeff += 4.0 * Sqr(yt) *(betamHu2 + betamu2(2,2) + betamq2(2,2)) + 8.0 * Tyt * betaTyt;
  coeff -= (1.6 * g1 * betag1 * sigma11 + 0.8 * Sqr(g1) * betasigma11);
  coeff += 4.0 * Qup * gN * betagN * sigma14 + 2.0 * Qup * Sqr(gN) * betasigma14;

  coeff = 0.5*oneOver16PiSqr*coeff;

  return coeff;
}

double doCalcmqL3SquaredLogCoeff(genericE6SSM_soft_parameters r, int nLps)
{
  if (nLps == 1 || nLps == 2)
    {
      r.set_loops(nLps);    
    }
  else
    {
      cerr << "WARNING: can only do 1 or 2 loops: returning nLps = 1 loop result." << endl;
      r.set_loops(1);      
    }

  genericE6SSM_soft_parameters w = r.calc_beta();

  return w.get_mq2(2,2);
}

double doCalcmqL3SquaredLogSqCoeff(genericE6SSM_soft_parameters r, int nLogs)
{
  // Assumes first and second generation Yukawas vanish, 
  // consistent with the previous pMSSM study
  r.set_Yu(0,0,0.0);
  r.set_Yu(1,1,0.0);

  r.set_Yd(0,0,0.0);
  r.set_Ye(1,1,0.0);

  r.set_Ye(0,0,0.0);
  r.set_Ye(1,1,0.0);

  // For consistency the corresponding A terms should also vanish
  r.set_TYu(0,0,0.0);
  r.set_TYu(1,1,0.0);

  r.set_TYd(0,0,0.0);
  r.set_TYd(1,1,0.0);

  r.set_TYe(0,0,0.0);
  r.set_TYe(1,1,0.0);

  if (nLogs == 0)
    {
      return 0.0;
    }

  if (nLogs != 1)
    {
      cerr << "WARNING: cannot yet do more than LL calculation: using nLogs = 1." << endl;
    }

  r.set_loops(1);

  double QQp = r.get_input().QQp;
  double Qup = r.get_input().Qup;
  double Qdp = r.get_input().Qdp;
  double QLp = r.get_input().QLp;
  double Qep = r.get_input().Qep;
  double QH1p = r.get_input().QH1p;
  double QH2p = r.get_input().QH2p;
  double QSp = r.get_input().QSp;
  double QHpp = r.get_input().QHpp;
  double QHpbarp = r.get_input().QHpbarp;
  double QDxp = r.get_input().QDxp;
  double QDxbarp = r.get_input().QDxbarp;

  genericE6SSM_soft_parameters w = r.calc_beta();

  double mHu2, mHd2, ms2, mHp2, mHpbar2;
  Eigen::Matrix<double,3,3> mq2, ml2, mu2, md2, me2, mDx2, mDxbar2;
  Eigen::Matrix<double,2,2> mH1I2, mH2I2, msI2;

  double betamHu2, betamHd2, betams2, betamHp2, betamHpbar2;
  Eigen::Matrix<double,3,3> betamq2, betaml2, betamu2, betamd2, betame2, betamDx2, betamDxbar2;
  Eigen::Matrix<double,2,2> betamH1I2, betamH2I2, betamsI2;

  mHu2 = r.get_mHu2(); betamHu2 = w.get_mHu2();
  mHd2 = r.get_mHd2(); betamHd2 = w.get_mHd2();
  ms2 = r.get_ms2(); betams2 = w.get_ms2();
  mHp2 = r.get_mHp2(); betamHp2 = w.get_mHp2();
  mHpbar2 = r.get_mHpbar2(); betamHpbar2 = w.get_mHpbar2();

  mq2 = r.get_mq2(); betamq2 = w.get_mq2();
  ml2 = r.get_ml2(); betaml2 = w.get_ml2();
  mu2 = r.get_mu2(); betamu2 = w.get_mu2();
  md2 = r.get_md2(); betamd2 = w.get_md2();
  me2 = r.get_me2(); betame2 = w.get_me2();
  mDx2 = r.get_mDx2(); betamDx2 = w.get_mDx2();
  mDxbar2 = r.get_mDxbar2(); betamDxbar2 = w.get_mDxbar2();

  mH1I2 = r.get_mH1I2(); betamH1I2 = w.get_mH1I2();
  mH2I2 = r.get_mH2I2(); betamH2I2 = w.get_mH2I2();
  msI2 = r.get_msI2(); betamsI2 = w.get_msI2();

  double M1, M2, M3, M1p;
  double betaM1, betaM2, betaM3, betaM1p;

  M1 = r.get_MassB(); betaM1 = w.get_MassB();
  M2 = r.get_MassWB(); betaM2 = w.get_MassWB();
  M3 = r.get_MassG(); betaM3 = w.get_MassG();
  M1p = r.get_MassBp(); betaM1p = w.get_MassBp();

  double g1, g2, g3, gN;
  double betag1, betag2, betag3, betagN;
  g1 = r.get_g1(); betag1 = w.get_g1();
  g2 = r.get_g2(); betag2 = w.get_g2();
  g3 = r.get_g3(); betag3 = w.get_g3();
  gN = r.get_gN(); betagN = w.get_gN();

  double yt, yb;
  double betayt, betayb;
  yt = r.get_Yu(2,2); betayt = w.get_Yu(2,2);
  yb = r.get_Yd(2,2); betayb = w.get_Yd(2,2);

  double Tyt, Tyb;
  double betaTyt, betaTyb;
  Tyt = r.get_TYu(2,2); betaTyt = w.get_TYu(2,2);
  Tyb = r.get_TYd(2,2); betaTyb = w.get_TYd(2,2);

  double sigma11 = mHu2 - mHd2 + mHpbar2 - mHp2 + md2.trace()
    - mDx2.trace() + mDxbar2.trace() + me2.trace() - mH1I2.trace()
    + mH2I2.trace() - ml2.trace() + mq2.trace() - 2.0 * mu2.trace();

  double sigma14 = 2.0 * QH1p * mHd2 + 2.0 * QH2p * mHu2 + 2.0 * QHpbarp * mHpbar2
    + 2.0 * QHpp * mHp2 + QSp * ms2 + 3.0 * Qdp * md2.trace() + 3.0 * QDxp * mDx2.trace()
    + 3.0 * QDxbarp * mDxbar2.trace() + Qep * me2.trace() + 2.0 * QH1p * mH1I2.trace()
    + 2.0 * QH2p * mH2I2.trace() + 2.0 * QLp * ml2.trace() + 6.0 * QQp * mq2.trace() + QSp * msI2.trace()
    + 3.0 * Qup * mu2.trace();

  double betasigma11 = betamHu2 - betamHd2 + betamHpbar2 - betamHp2 + betamd2.trace()
    - betamDx2.trace() + betamDxbar2.trace() + betame2.trace() - betamH1I2.trace()
    + betamH2I2.trace() - betaml2.trace() + betamq2.trace() - 2.0 * betamu2.trace();

  double betasigma14 = 2.0 * QH1p * betamHd2 + 2.0 * QH2p * betamHu2 + 2.0 * QHpbarp * betamHpbar2
    + 2.0 * QHpp * betamHp2 + QSp * betams2 + 3.0 * Qdp * betamd2.trace() + 3.0 * QDxp * betamDx2.trace()
    + 3.0 * QDxbarp * betamDxbar2.trace() + Qep * betame2.trace() + 2.0 * QH1p * betamH1I2.trace()
    + 2.0 * QH2p * betamH2I2.trace() + 2.0 * QLp * betaml2.trace() + 6.0 * QQp * betamq2.trace() + QSp * betamsI2.trace()
    + 3.0 * Qup * betamu2.trace();

  double coeff = -16.0 * Sqr(QQp) * (gN * betagN * Sqr(M1p) + Sqr(gN) * M1p * betaM1p);
  coeff -= (4.0/15.0) * (g1 * betag1 * Sqr(M1) + Sqr(g1) * M1 * betaM1);
  coeff -= (64.0/3.0) * (g3 * betag3 * Sqr(M3) + Sqr(g3) * M3 * betaM3);
  coeff -= 12.0 * (g2 * betag2 * Sqr(M2) + Sqr(g2) * M2 * betaM2);
  coeff += 4.0 * yt * betayt * (mHu2 + mu2(2,2) + mq2(2,2));
  coeff += 2.0 * Sqr(yt) * (betamHu2 + betamu2(2,2) + betamq2(2,2)) + 4.0 * Tyt * betaTyt;
  coeff += 4.0 * yb * betayb * (mHd2 + md2(2,2) + mq2(2,2));
  coeff += 2.0 * Sqr(yb) * (betamHd2 + betamd2(2,2) + betamq2(2,2)) + 4.0 * Tyb * betaTyb;
  coeff += 0.4 * g1 * betag1 * sigma11 + 0.2 * Sqr(g1) * betasigma11;
  coeff += 4.0 * QQp * gN * betagN * sigma14 + 2.0 * QQp * Sqr(gN) * betasigma14;

  coeff = 0.5*oneOver16PiSqr*coeff;

  return coeff;
}

double doCalcLambda3LogCoeff(genericE6SSM_soft_parameters r, int nLps)
{
  if (nLps == 1 || nLps == 2)
    {
      r.set_loops(nLps);    
    }
  else
    {
      cerr << "WARNING: can only do 1 or 2 loops: returning nLps = 1 loop result." << endl;
      r.set_loops(1);      
    }

  genericE6SSM_soft_parameters w = r.calc_beta();

  return w.get_Lambdax();
}

double doCalcLambda3LogSqCoeff(genericE6SSM_soft_parameters r, int nLogs)
{
  // Assumes first and second generation Yukawas vanish, 
  // consistent with the previous pMSSM study
  r.set_Yu(0,0,0.0);
  r.set_Yu(1,1,0.0);

  r.set_Yd(0,0,0.0);
  r.set_Yd(1,1,0.0);

  r.set_Ye(0,0,0.0);
  r.set_Ye(1,1,0.0);

  // For consistency the corresponding A terms should also vanish
  r.set_TYu(0,0,0.0);
  r.set_TYu(1,1,0.0);

  r.set_TYd(0,0,0.0);
  r.set_TYd(1,1,0.0);

  r.set_TYe(0,0,0.0);
  r.set_TYe(1,1,0.0);

  if (nLogs == 0)
    {
      return 0.0;
    }

  if (nLogs != 1)
    {
      cerr << "WARNING: cannot yet do more than LL calculation: using nLogs = 1." << endl;
    }

  r.set_loops(1);

  double QQp = r.get_input().QQp;
  double Qup = r.get_input().Qup;
  double Qdp = r.get_input().Qdp;
  double QLp = r.get_input().QLp;
  double Qep = r.get_input().Qep;
  double QH1p = r.get_input().QH1p;
  double QH2p = r.get_input().QH2p;
  double QSp = r.get_input().QSp;
  double QHpp = r.get_input().QHpp;
  double QHpbarp = r.get_input().QHpbarp;
  double QDxp = r.get_input().QDxp;
  double QDxbarp = r.get_input().QDxbarp;

  genericE6SSM_soft_parameters w = r.calc_beta();

  // Save required matrices prior to calculation; this seems
  // to be necessary for successful running on the cloud
  double g1, g2, gN;
  double betag1, betag2, betagN;
  double lambda, betalambda;

  g1 = r.get_g1(); betag1 = w.get_g1();
  g2 = r.get_g2(); betag2 = w.get_g2();
  gN = r.get_gN(); betagN = w.get_gN();
  lambda = r.get_Lambdax(); betalambda = w.get_Lambdax();

  Eigen::Matrix<double,3,3> Yu, Yd, Ye, Kappa, BetaYu, BetaYd, BetaYe, BetaKappa;
  Eigen::Matrix<double,3,3> YuAdj, YdAdj, YeAdj, KappaAdj, BetaYuAdj, BetaYdAdj, BetaYeAdj, BetaKappaAdj;
  Eigen::Matrix<double,2,2> Lambda12, BetaLambda12, Lambda12Adj, BetaLambda12Adj;

  Yu = r.get_Yu(); YuAdj = Yu.adjoint();
  Yd = r.get_Yd(); YdAdj = Yd.adjoint();
  Ye = r.get_Ye(); YeAdj = Ye.adjoint();
  Kappa = r.get_Kappa(); KappaAdj = Kappa.adjoint();
  Lambda12 = r.get_Lambda12(); Lambda12Adj = Lambda12.adjoint();

  BetaYu = w.get_Yu(); BetaYuAdj = BetaYu.adjoint();
  BetaYd = w.get_Yd(); BetaYdAdj = BetaYd.adjoint();
  BetaYe = w.get_Ye(); BetaYeAdj = BetaYe.adjoint();
  BetaKappa = w.get_Kappa(); BetaKappaAdj = BetaKappa.adjoint();
  BetaLambda12 = w.get_Lambda12(); BetaLambda12Adj = BetaLambda12.adjoint();

  double coeff = betalambda*(-0.6 * Sqr(g1) - 3.0 * Sqr(g2) - 2.0 * Sqr(gN) * (Sqr(QH1p) + Sqr(QH2p) + Sqr(QSp))
  				  + 4.0 * Sqr(lambda) + 3.0 * (Yd*YdAdj).trace()
  				  +(Ye*YeAdj).trace() + 3.0 * (Yu*YuAdj).trace()
  				  + 3.0 * (Kappa*KappaAdj).trace() 
  				  + 2.0 * (Lambda12*Lambda12Adj).trace());
  coeff += lambda*(-1.2 * g1 * betag1 - 6.0 * g2 * betag2 
		   - 4.0 * gN * betagN * (Sqr(QH1p) + Sqr(QH2p) + Sqr(QSp)) + 8.0 * lambda * betalambda
		   + 3.0 * (BetaYd*YdAdj).trace() + 3.0 * (Yd*BetaYdAdj).trace()
		   + (BetaYe*YeAdj).trace() + (Ye*BetaYeAdj).trace()
		   + 3.0 * (BetaYu*YuAdj).trace() + 3.0 * (Yu*BetaYuAdj).trace()
		   + 3.0 * (BetaKappa*KappaAdj).trace() 
		   + 3.0 * (Kappa*BetaKappaAdj).trace()
		   + 2.0 * (BetaLambda12*Lambda12Adj).trace() 
		   + 2.0 * (Lambda12*BetaLambda12Adj).trace());

  coeff = 0.5*oneOver16PiSqr*coeff;

  return coeff;
}

double doCalcAlambda3LogCoeff(genericE6SSM_soft_parameters r, int nLps)
{
  if (nLps == 1 || nLps == 2)
    {
      r.set_loops(nLps);    
    }
  else
    {
      cerr << "WARNING: can only do 1 or 2 loops: returning nLps = 1 loop result." << endl;
      r.set_loops(1);      
    }

  genericE6SSM_soft_parameters w = r.calc_beta();

  // Note that A_lambda must be calculated from the beta functions 
  // for lambda and T_lambda
  double lambda, Tlambda, Alambda;
  
  lambda = r.get_Lambdax();
  Tlambda = r.get_TLambdax();
  
  if (Abs(Tlambda) < EPSTOL) 
    {
      Alambda = 0.0;
    }
  if (Abs(lambda) < 1.0e-100)
    {
      ostringstream ii;
      ii << "WARNING: trying to calculate A_lambda(3) where lambda(3) coupling is " <<
	Abs(lambda) << endl;
      throw ii.str();
    }
  else
    {
      Alambda = Tlambda/lambda;
    }

  double beta_Alambda = (w.get_TLambdax() - Alambda * w.get_Lambdax())/lambda; 

  return beta_Alambda;
}

double doCalcAlambda3LogSqCoeff(genericE6SSM_soft_parameters r, int nLogs)
{
  // Assumes first and second generation Yukawas vanish, 
  // consistent with the previous pMSSM study
  r.set_Yu(0,0,0.0);
  r.set_Yu(1,1,0.0);

  r.set_Yd(0,0,0.0);
  r.set_Ye(1,1,0.0);

  r.set_Ye(0,0,0.0);
  r.set_Ye(1,1,0.0);

  // For consistency the corresponding A terms should also vanish
  r.set_TYu(0,0,0.0);
  r.set_TYu(1,1,0.0);

  r.set_TYd(0,0,0.0);
  r.set_TYd(1,1,0.0);

  r.set_TYe(0,0,0.0);
  r.set_TYe(1,1,0.0);

  if (nLogs == 0)
    {
      return 0.0;
    }

  if (nLogs != 1)
    {      cerr << "WARNING: cannot yet do more than LL calculation: using nLogs = 1." << endl;
    }

  r.set_loops(1);

  double QQp = r.get_input().QQp;
  double Qup = r.get_input().Qup;
  double Qdp = r.get_input().Qdp;
  double QLp = r.get_input().QLp;
  double Qep = r.get_input().Qep;
  double QH1p = r.get_input().QH1p;
  double QH2p = r.get_input().QH2p;
  double QSp = r.get_input().QSp;
  double QHpp = r.get_input().QHpp;
  double QHpbarp = r.get_input().QHpbarp;
  double QDxp = r.get_input().QDxp;
  double QDxbarp = r.get_input().QDxbarp;

  // Note that A_lambda must be calculated from the beta functions 
  // for lambda and T_lambda
  double lambda, Tlambda, Alambda;
  
  lambda = r.get_Lambdax();
  Tlambda = r.get_TLambdax();
  
  if (Abs(Tlambda) < EPSTOL) 
    {
      Alambda = 0.0;
    }
  if (Abs(lambda) < 1.0e-100)
    {
      ostringstream ii;
      ii << "WARNING: trying to calculate A_lambda(3) where lambda(3) coupling is " <<
	Abs(lambda) << endl;
      throw ii.str();
    }
  else
    {
      Alambda = Tlambda/lambda;
    }

  genericE6SSM_soft_parameters w = r.calc_beta();

  double M1, M2, M1p;
  double betaM1, betaM2, betaM1p;

  M1 = r.get_MassB(); betaM1 = w.get_MassB();
  M2 = r.get_MassWB(); betaM2 = w.get_MassWB();
  M1p = r.get_MassBp(); betaM1p = w.get_MassBp();

  double g1, g2, gN;
  double  betag1, betag2, betagN, betalambda, betaTlambda;
  g1 = r.get_g1(); betag1 = w.get_g1();
  g2 = r.get_g2(); betag2 = w.get_g2();
  gN = r.get_gN(); betagN = w.get_gN();
  betalambda = w.get_Lambdax();
  betaTlambda = w.get_TLambdax();

  Eigen::Matrix<double,3,3> YuAdj, YdAdj, YeAdj, KappaAdj, BetaYuAdj, BetaYdAdj, BetaYeAdj, BetaKappaAdj;
  Eigen::Matrix<double,3,3> TYu, TYd, TYe, TKappa, BetaTYu, BetaTYd, BetaTYe, BetaTKappa;
  
  YuAdj = r.get_Yu().adjoint();
  YdAdj = r.get_Yd().adjoint();
  YeAdj = r.get_Ye().adjoint();
  KappaAdj = r.get_Kappa().adjoint();

  BetaYuAdj = w.get_Yu().adjoint();
  BetaYdAdj = w.get_Yd().adjoint();
  BetaYeAdj = w.get_Ye().adjoint();
  BetaKappaAdj = w.get_Kappa().adjoint();

  TYu = r.get_TYu(); BetaTYu = w.get_TYu();
  TYd = r.get_TYd(); BetaTYd = w.get_TYd();
  TYe = r.get_TYe(); BetaTYe = w.get_TYe();
  TKappa = r.get_TKappa(); BetaTKappa = w.get_TKappa();

  Eigen::Matrix<double,2,2> Lambda12Adj, BetaLambda12Adj;
  Eigen::Matrix<double,2,2> TLambda12, BetaTLambda12;

  Lambda12Adj = r.get_Lambda12().adjoint();
  BetaLambda12Adj = w.get_Lambda12().adjoint();
  TLambda12 = r.get_TLambda12();
  BetaTLambda12 = w.get_TLambda12();

  double beta_Alambda = (betaTlambda - Alambda * betalambda)/lambda; 

  double coeff = 2.4 * g1 * betag1 * M1 + 1.2 * Sqr(g1) * betaM1;
  coeff += 12.0 * g2 * betag2 * M2 + 6.0 * Sqr(g2) * betaM2;
  coeff += (8.0 * gN * betagN * M1p + 4.0 * Sqr(gN) * betaM1p)* (Sqr(QH1p) + Sqr(QH2p) + Sqr(QSp));
  coeff += 6.0 * (BetaYdAdj*TYd).trace() + 6.0 * (YdAdj*BetaTYd).trace();
  coeff += 6.0 * (BetaYuAdj*TYu).trace() + 6.0 * (YuAdj*BetaTYu).trace();
  coeff += 2.0 * (BetaYeAdj*TYe).trace() + 2.0 * (YeAdj*BetaTYe).trace();
  coeff += 6.0 * (BetaKappaAdj*TKappa).trace() + 6.0 * (KappaAdj*BetaTKappa).trace();
  coeff += 4.0 * (BetaLambda12Adj*TLambda12).trace() + 4.0 * (Lambda12Adj*BetaTLambda12).trace();
  coeff += 16.0 * lambda * betalambda * Alambda + 8.0 * Sqr(lambda) * beta_Alambda;

  coeff = 0.5*oneOver16PiSqr*coeff;

  return coeff;
}

double doCalcAtLogCoeff(genericE6SSM_soft_parameters r, int nLps)
{
  if (nLps == 1 || nLps == 2)
    {
      r.set_loops(nLps);    
    }
  else
    {
      cerr << "WARNING: can only do 1 or 2 loops: returning nLps = 1 loop result." << endl;
      r.set_loops(1);      
    }

  genericE6SSM_soft_parameters w = r.calc_beta();

  // Note that A_t must be calculated from the beta functions 
  // for y_t and T_{y_t}
  double yt, Tyt, At;
  
  yt = r.get_Yu(2,2);
  Tyt = r.get_TYu(2,2);
  
  if (Abs(Tyt) < EPSTOL) 
    {
      At = 0.0;
    }
  if (Abs(yt) < 1.0e-100)
    {
      ostringstream ii;
      ii << "WARNING: trying to calculate A_t where y_t coupling is " <<
	Abs(yt) << endl;
      throw ii.str();
    }
  else
    {
      At = Tyt/yt;
    }

  double beta_At = (w.get_TYu(2,2) - At * w.get_Yu(2,2))/yt; 

  return beta_At;
}

double doCalcAtLogSqCoeff(genericE6SSM_soft_parameters r, int nLogs)
{
  // Assumes first and second generation Yukawas vanish, 
  // consistent with the previous pMSSM study
  r.set_Yu(0,0,0.0);
  r.set_Yu(1,1,0.0);

  r.set_Yd(0,0,0.0);
  r.set_Ye(1,1,0.0);

  r.set_Ye(0,0,0.0);
  r.set_Ye(1,1,0.0);

  // For consistency the corresponding A terms should also vanish
  r.set_TYu(0,0,0.0);
  r.set_TYu(1,1,0.0);

  r.set_TYd(0,0,0.0);
  r.set_TYd(1,1,0.0);

  r.set_TYe(0,0,0.0);
  r.set_TYe(1,1,0.0);

  if (nLogs == 0)
    {
      return 0.0;
    }

  if (nLogs != 1)
    {
      cerr << "WARNING: cannot yet do more than LL calculation: using nLogs = 1." << endl;
    }

  r.set_loops(1);

  double QQp = r.get_input().QQp;
  double Qup = r.get_input().Qup;
  double Qdp = r.get_input().Qdp;
  double QLp = r.get_input().QLp;
  double Qep = r.get_input().Qep;
  double QH1p = r.get_input().QH1p;
  double QH2p = r.get_input().QH2p;
  double QSp = r.get_input().QSp;
  double QHpp = r.get_input().QHpp;
  double QHpbarp = r.get_input().QHpbarp;
  double QDxp = r.get_input().QDxp;
  double QDxbarp = r.get_input().QDxbarp;

  // Note that A_t must be calculated from the beta functions 
  // for y_t and T_{y_t}
  double yt, Tyt, At;
  
  yt = r.get_Yu(2,2);
  Tyt = r.get_TYu(2,2);
  
  if (Abs(Tyt) < EPSTOL) 
    {
      At = 0.0;
    }
  if (Abs(yt) < 1.0e-100)
    {
      ostringstream ii;
      ii << "WARNING: trying to calculate A_t where y_t coupling is " <<
	Abs(yt) << endl;
      throw ii.str();
    }
  else
    {
      At = Tyt/yt;
    }

  double yb, Tyb, Ab;
  
  yb = r.get_Yd(2,2);
  Tyb = r.get_TYd(2,2);
  
  if (Abs(Tyb) < EPSTOL) 
    {
      Ab = 0.0;
    }
  if (Abs(yb) < 1.0e-100)
    {
      ostringstream ii;
      ii << "WARNING: trying to calculate A_b where y_b coupling is " <<
	Abs(yb) << endl;
      throw ii.str();
    }
  else
    {
      Ab = Tyb/yb;
    }

  // Note that A_lambda must be calculated from the beta functions 
  // for lambda and T_lambda
  double lambda, Tlambda, Alambda;
  
  lambda = r.get_Lambdax();
  Tlambda = r.get_TLambdax();
  
  if (Abs(Tlambda) < EPSTOL) 
    {
      Alambda = 0.0;
    }
  if (Abs(lambda) < 1.0e-100)
    {
      ostringstream ii;
      ii << "WARNING: trying to calculate A_lambda(3) where lambda(3) coupling is " <<
	Abs(lambda) << endl;
      throw ii.str();
    }
  else
    {
      Alambda = Tlambda/lambda;
    }

  genericE6SSM_soft_parameters w = r.calc_beta();

  double M1, M2, M3, M1p;
  double betaM1, betaM2, betaM3, betaM1p;

  M1 = r.get_MassB(); betaM1 = w.get_MassB();
  M2 = r.get_MassWB(); betaM2 = w.get_MassWB();
  M3 = r.get_MassG(); betaM3 = w.get_MassG();
  M1p = r.get_MassBp(); betaM1p = w.get_MassBp();

  double g1, g2, g3, gN;
  double betag1, betag2, betag3, betagN;
  double betalambda, betayt, betayb;

  g1 = r.get_g1(); betag1 = w.get_g1();
  g2 = r.get_g2(); betag2 = w.get_g2();
  g3 = r.get_g3(); betag3 = w.get_g3();
  gN = r.get_gN(); betagN = w.get_gN();

  betalambda = w.get_Lambdax();
  betayt = w.get_Yu(2,2);
  betayb = w.get_Yd(2,2);

  double betaAt = (w.get_TYu(2,2) - At * w.get_Yu(2,2))/yt; 
  double betaAb = (w.get_TYd(2,2) - Ab * w.get_Yd(2,2))/yb; 
  double betaAlambda = (w.get_TLambdax() - Alambda * w.get_Lambdax())/lambda; 

  double coeff = 4.0 * lambda * betalambda * Alambda + 2.0 * Sqr(lambda) * betaAlambda;
  coeff += 24.0 * yt * betayt * At + 12.0 * Sqr(yt) * betaAt;
  coeff += 4.0 * yb * betayb * Ab + 2.0 * Sqr(yb) * betaAb;
  coeff += (52.0/15.0) * g1 * betag1 * M1 + (26.0/15.0) * Sqr(g1) * betaM1;
  coeff += 12.0 * g2 * betag2 * M2 + 6.0 * Sqr(g2) * betaM2;
  coeff += (64.0/3.0) * g3 * betag3 * M3 + (32.0/3.0) * Sqr(g3) * betaM3;
  coeff += (8.0 * gN * betagN * M1p + 4.0 * Sqr(gN) * betaM1p) * (Sqr(QQp) + Sqr(QH2p) + Sqr(Qup));

  coeff = 0.5*oneOver16PiSqr*coeff;

  return coeff;
}

// Functions for getting approximate derivatives of low scale parameters w.r.t high
// scale parameters. Returns the vector
// [ dlambda/dp dAlambda/dp dm_Hd^2/dp dm_Hu^2/dp dm_s^2/dp dm_Ql^2/dp dm_uR^2/dp dA_t/dp]^T. 
// Since there is a very large number of couplings in the E6SSM, it it not very practical
// to pass everything in as a vector of parameters, so we pass the object at MX instead.
// Because only a small subset of the high scale parameters have a significant contribution to
// the tuning, and to save time, for now we only calculate fine tuning for the third generation parameters
// listed below.
Eigen::Matrix<double,8,1> doCalcMh1SquaredDerivs(genericE6SSM_soft_parameters r, double ms, int gen,
				      bool & hasError)
{
  if (gen < 1 || gen > 3)
    {
      ostringstream ii;
      ii << "WARNING: can only do calculation of derivatives for m_{H1i}^2, i = 1, 2, 3, not m_{H1" << gen << "}^2" << endl;
      throw ii.str(); 
    }

  double mx = r.get_scale();

  // Taylor approximation to the RGE solution is assumed to be about the values at MX.
  double t = log(ms/mx);

  Eigen::Matrix<double,8,1> derivs;

  double dlambdadp, dAlambdadp, dmHdSqdp, dmHuSqdp, dmSSqdp, dmQlSqdp, dmUrSqdp, dAtdp;

  // These are zero for all soft squared masses
  dAlambdadp = 0.0;
  dAtdp = 0.0;
  dlambdadp = 0.0;

  // Calculate each of the remaining non-zero derivatives:
  // For derivatives of soft masses w.r.t. soft masses, take advantage
  // of the fact that beta functions are linear in the soft masses,
  // so can calculate derivatives by simple rise/run formula.
  int nLps = r.get_loops();
  int nLogs = 1;
  double shift = 1.0;

  double firstmHdSqLogCoeff = doCalcMh1SquaredLogCoeff(r, nLps);
  double firstmHdSqLogSqCoeff = doCalcMh1SquaredLogSqCoeff(r, nLogs);

  double firstmHuSqLogCoeff = doCalcMh2SquaredLogCoeff(r, nLps);
  double firstmHuSqLogSqCoeff = doCalcMh2SquaredLogSqCoeff(r, nLogs);

  double firstmSSqLogCoeff = doCalcMsSquaredLogCoeff(r, nLps);
  double firstmSSqLogSqCoeff = doCalcMsSquaredLogSqCoeff(r, nLogs);

  double firstmQlSqLogCoeff = doCalcmqL3SquaredLogCoeff(r, nLps);
  double firstmQlSqLogSqCoeff = doCalcmqL3SquaredLogSqCoeff(r, nLogs);

  double firstmUrSqLogCoeff = doCalcmtRSquaredLogCoeff(r, nLps);
  double firstmUrSqLogSqCoeff = doCalcmtRSquaredLogSqCoeff(r, nLogs);

  if (gen == 3)
    {
      r.set_mHd2(r.get_mHd2()+shift);
    }
  else
    {
      r.set_mH1I2(gen-1,gen-1, r.get_mH1I2(gen-1,gen-1)+shift);
    }

  double secondmHdSqLogCoeff = doCalcMh1SquaredLogCoeff(r, nLps);
  double secondmHdSqLogSqCoeff = doCalcMh1SquaredLogSqCoeff(r, nLogs);

  double secondmHuSqLogCoeff = doCalcMh2SquaredLogCoeff(r, nLps);
  double secondmHuSqLogSqCoeff = doCalcMh2SquaredLogSqCoeff(r, nLogs);

  double secondmSSqLogCoeff = doCalcMsSquaredLogCoeff(r, nLps);
  double secondmSSqLogSqCoeff = doCalcMsSquaredLogSqCoeff(r, nLogs);

  double secondmQlSqLogCoeff = doCalcmqL3SquaredLogCoeff(r, nLps);
  double secondmQlSqLogSqCoeff = doCalcmqL3SquaredLogSqCoeff(r, nLogs);

  double secondmUrSqLogCoeff = doCalcmtRSquaredLogCoeff(r, nLps);
  double secondmUrSqLogSqCoeff = doCalcmtRSquaredLogSqCoeff(r, nLogs);

  double dmHdSqOtCoeffdp = (secondmHdSqLogCoeff-firstmHdSqLogCoeff)/shift;
  double dmHdSqOt2Coeffdp = (secondmHdSqLogSqCoeff-firstmHdSqLogSqCoeff)/shift;

  double dmHuSqOtCoeffdp = (secondmHuSqLogCoeff-firstmHuSqLogCoeff)/shift;
  double dmHuSqOt2Coeffdp = (secondmHuSqLogSqCoeff-firstmHuSqLogSqCoeff)/shift;

  double dmSSqOtCoeffdp = (secondmSSqLogCoeff-firstmSSqLogCoeff)/shift;
  double dmSSqOt2Coeffdp = (secondmSSqLogSqCoeff-firstmSSqLogSqCoeff)/shift;

  double dmQlSqOtCoeffdp = (secondmQlSqLogCoeff-firstmQlSqLogCoeff)/shift;
  double dmQlSqOt2Coeffdp = (secondmQlSqLogSqCoeff-firstmQlSqLogSqCoeff)/shift;

  double dmUrSqOtCoeffdp = (secondmUrSqLogCoeff-firstmUrSqLogCoeff)/shift;
  double dmUrSqOt2Coeffdp = (secondmUrSqLogSqCoeff-firstmUrSqLogSqCoeff)/shift;

  derivs(0) = dlambdadp;
  derivs(1) = dAlambdadp;
  derivs(2) = KroneckerDelta(gen, 3) + t * dmHdSqOtCoeffdp + Sqr(t) * dmHdSqOt2Coeffdp;
  derivs(3) = t * dmHuSqOtCoeffdp + Sqr(t) * dmHuSqOt2Coeffdp;
  derivs(4) = t * dmSSqOtCoeffdp + Sqr(t) * dmSSqOt2Coeffdp;
  derivs(5) = t * dmQlSqOtCoeffdp + Sqr(t) * dmQlSqOt2Coeffdp;
  derivs(6) = t * dmUrSqOtCoeffdp + Sqr(t) * dmUrSqOt2Coeffdp;
  derivs(7) = dAtdp;

  return derivs;

}

Eigen::Matrix<double,8,1> doCalcMh2SquaredDerivs(genericE6SSM_soft_parameters r, double ms, int gen,
				      bool & hasError)
{
  if (gen < 1 || gen > 3)
    {
      ostringstream ii;
      ii << "WARNING: can only do calculation of derivatives for m_{H2i}^2, i = 1, 2, 3, not m_{H2" << gen << "}^2" << endl;
      throw ii.str(); 
    }

  double mx = r.get_scale();

  // Taylor approximation to the RGE solution is assumed to be about the values at MX.
  double t = log(ms/mx);

  Eigen::Matrix<double,8,1> derivs;

  double dlambdadp, dAlambdadp, dmHdSqdp, dmHuSqdp, dmSSqdp, dmQlSqdp, dmUrSqdp, dAtdp;

  // These are zero for all soft squared masses
  dAlambdadp = 0.0;
  dAtdp = 0.0;
  dlambdadp = 0.0;

  // Calculate each of the remaining non-zero derivatives:
  // For derivatives of soft masses w.r.t. soft masses, take advantage
  // of the fact that beta functions are linear in the soft masses,
  // so can calculate derivatives by simple rise/run formula.
  int nLps = r.get_loops();
  int nLogs = 1;
  double shift = 1.0;

  double firstmHdSqLogCoeff = doCalcMh1SquaredLogCoeff(r, nLps);
  double firstmHdSqLogSqCoeff = doCalcMh1SquaredLogSqCoeff(r, nLogs);

  double firstmHuSqLogCoeff = doCalcMh2SquaredLogCoeff(r, nLps);
  double firstmHuSqLogSqCoeff = doCalcMh2SquaredLogSqCoeff(r, nLogs);

  double firstmSSqLogCoeff = doCalcMsSquaredLogCoeff(r, nLps);
  double firstmSSqLogSqCoeff = doCalcMsSquaredLogSqCoeff(r, nLogs);

  double firstmQlSqLogCoeff = doCalcmqL3SquaredLogCoeff(r, nLps);
  double firstmQlSqLogSqCoeff = doCalcmqL3SquaredLogSqCoeff(r, nLogs);

  double firstmUrSqLogCoeff = doCalcmtRSquaredLogCoeff(r, nLps);
  double firstmUrSqLogSqCoeff = doCalcmtRSquaredLogSqCoeff(r, nLogs);

  if (gen == 3)
    {
      r.set_mHu2(r.get_mHu2()+shift);
    }
  else
    {
      r.set_mH2I2(gen-1,gen-1, r.get_mH2I2(gen-1,gen-1)+shift);
    }

  double secondmHdSqLogCoeff = doCalcMh1SquaredLogCoeff(r, nLps);
  double secondmHdSqLogSqCoeff = doCalcMh1SquaredLogSqCoeff(r, nLogs);

  double secondmHuSqLogCoeff = doCalcMh2SquaredLogCoeff(r, nLps);
  double secondmHuSqLogSqCoeff = doCalcMh2SquaredLogSqCoeff(r, nLogs);

  double secondmSSqLogCoeff = doCalcMsSquaredLogCoeff(r, nLps);
  double secondmSSqLogSqCoeff = doCalcMsSquaredLogSqCoeff(r, nLogs);

  double secondmQlSqLogCoeff = doCalcmqL3SquaredLogCoeff(r, nLps);
  double secondmQlSqLogSqCoeff = doCalcmqL3SquaredLogSqCoeff(r, nLogs);

  double secondmUrSqLogCoeff = doCalcmtRSquaredLogCoeff(r, nLps);
  double secondmUrSqLogSqCoeff = doCalcmtRSquaredLogSqCoeff(r, nLogs);

  double dmHdSqOtCoeffdp = (secondmHdSqLogCoeff-firstmHdSqLogCoeff)/shift;
  double dmHdSqOt2Coeffdp = (secondmHdSqLogSqCoeff-firstmHdSqLogSqCoeff)/shift;

  double dmHuSqOtCoeffdp = (secondmHuSqLogCoeff-firstmHuSqLogCoeff)/shift;
  double dmHuSqOt2Coeffdp = (secondmHuSqLogSqCoeff-firstmHuSqLogSqCoeff)/shift;

  double dmSSqOtCoeffdp = (secondmSSqLogCoeff-firstmSSqLogCoeff)/shift;
  double dmSSqOt2Coeffdp = (secondmSSqLogSqCoeff-firstmSSqLogSqCoeff)/shift;

  double dmQlSqOtCoeffdp = (secondmQlSqLogCoeff-firstmQlSqLogCoeff)/shift;
  double dmQlSqOt2Coeffdp = (secondmQlSqLogSqCoeff-firstmQlSqLogSqCoeff)/shift;

  double dmUrSqOtCoeffdp = (secondmUrSqLogCoeff-firstmUrSqLogCoeff)/shift;
  double dmUrSqOt2Coeffdp = (secondmUrSqLogSqCoeff-firstmUrSqLogSqCoeff)/shift;

  derivs(0) = dlambdadp;
  derivs(1) = dAlambdadp;
  derivs(2) = t * dmHdSqOtCoeffdp + Sqr(t) * dmHdSqOt2Coeffdp;
  derivs(3) = KroneckerDelta(gen, 3) + t * dmHuSqOtCoeffdp + Sqr(t) * dmHuSqOt2Coeffdp;
  derivs(4) = t * dmSSqOtCoeffdp + Sqr(t) * dmSSqOt2Coeffdp;
  derivs(5) = t * dmQlSqOtCoeffdp + Sqr(t) * dmQlSqOt2Coeffdp;
  derivs(6) = t * dmUrSqOtCoeffdp + Sqr(t) * dmUrSqOt2Coeffdp;
  derivs(7) = dAtdp;

  return derivs;


}

Eigen::Matrix<double,8,1> doCalcMsSquaredDerivs(genericE6SSM_soft_parameters r, double ms, int gen,
				      bool & hasError)
{
  if (gen < 1 || gen > 3)
    {
      ostringstream ii;
      ii << "WARNING: can only do calculation of derivatives for m_{si}^2, i = 1, 2, 3, not m_{s" << gen << "}^2" << endl;
      throw ii.str(); 
    }

  double mx = r.get_scale();

  // Taylor approximation to the RGE solution is assumed to be about the values at MX.
  double t = log(ms/mx);

  Eigen::Matrix<double,8,1> derivs;

  double dlambdadp, dAlambdadp, dmHdSqdp, dmHuSqdp, dmSSqdp, dmQlSqdp, dmUrSqdp, dAtdp;

  // These are zero for all soft squared masses
  dAlambdadp = 0.0;
  dAtdp = 0.0;
  dlambdadp = 0.0;

  // Calculate each of the remaining non-zero derivatives:
  // For derivatives of soft masses w.r.t. soft masses, take advantage
  // of the fact that beta functions are linear in the soft masses,
  // so can calculate derivatives by simple rise/run formula.
  int nLps = r.get_loops();
  int nLogs = 1;
  double shift = 1.0;

  double firstmHdSqLogCoeff = doCalcMh1SquaredLogCoeff(r, nLps);
  double firstmHdSqLogSqCoeff = doCalcMh1SquaredLogSqCoeff(r, nLogs);

  double firstmHuSqLogCoeff = doCalcMh2SquaredLogCoeff(r, nLps);
  double firstmHuSqLogSqCoeff = doCalcMh2SquaredLogSqCoeff(r, nLogs);

  double firstmSSqLogCoeff = doCalcMsSquaredLogCoeff(r, nLps);
  double firstmSSqLogSqCoeff = doCalcMsSquaredLogSqCoeff(r, nLogs);

  double firstmQlSqLogCoeff = doCalcmqL3SquaredLogCoeff(r, nLps);
  double firstmQlSqLogSqCoeff = doCalcmqL3SquaredLogSqCoeff(r, nLogs);

  double firstmUrSqLogCoeff = doCalcmtRSquaredLogCoeff(r, nLps);
  double firstmUrSqLogSqCoeff = doCalcmtRSquaredLogSqCoeff(r, nLogs);

  if (gen == 3)
    {
      r.set_ms2(r.get_ms2()+shift);
    }
  else
    {
      r.set_msI2(gen-1,gen-1, r.get_msI2(gen-1,gen-1)+shift);
    }

  double secondmHdSqLogCoeff = doCalcMh1SquaredLogCoeff(r, nLps);
  double secondmHdSqLogSqCoeff = doCalcMh1SquaredLogSqCoeff(r, nLogs);

  double secondmHuSqLogCoeff = doCalcMh2SquaredLogCoeff(r, nLps);
  double secondmHuSqLogSqCoeff = doCalcMh2SquaredLogSqCoeff(r, nLogs);

  double secondmSSqLogCoeff = doCalcMsSquaredLogCoeff(r, nLps);
  double secondmSSqLogSqCoeff = doCalcMsSquaredLogSqCoeff(r, nLogs);

  double secondmQlSqLogCoeff = doCalcmqL3SquaredLogCoeff(r, nLps);
  double secondmQlSqLogSqCoeff = doCalcmqL3SquaredLogSqCoeff(r, nLogs);

  double secondmUrSqLogCoeff = doCalcmtRSquaredLogCoeff(r, nLps);
  double secondmUrSqLogSqCoeff = doCalcmtRSquaredLogSqCoeff(r, nLogs);

  double dmHdSqOtCoeffdp = (secondmHdSqLogCoeff-firstmHdSqLogCoeff)/shift;
  double dmHdSqOt2Coeffdp = (secondmHdSqLogSqCoeff-firstmHdSqLogSqCoeff)/shift;

  double dmHuSqOtCoeffdp = (secondmHuSqLogCoeff-firstmHuSqLogCoeff)/shift;
  double dmHuSqOt2Coeffdp = (secondmHuSqLogSqCoeff-firstmHuSqLogSqCoeff)/shift;

  double dmSSqOtCoeffdp = (secondmSSqLogCoeff-firstmSSqLogCoeff)/shift;
  double dmSSqOt2Coeffdp = (secondmSSqLogSqCoeff-firstmSSqLogSqCoeff)/shift;

  double dmQlSqOtCoeffdp = (secondmQlSqLogCoeff-firstmQlSqLogCoeff)/shift;
  double dmQlSqOt2Coeffdp = (secondmQlSqLogSqCoeff-firstmQlSqLogSqCoeff)/shift;

  double dmUrSqOtCoeffdp = (secondmUrSqLogCoeff-firstmUrSqLogCoeff)/shift;
  double dmUrSqOt2Coeffdp = (secondmUrSqLogSqCoeff-firstmUrSqLogSqCoeff)/shift;

  derivs(0) = dlambdadp;
  derivs(1) = dAlambdadp;
  derivs(2) = t * dmHdSqOtCoeffdp + Sqr(t) * dmHdSqOt2Coeffdp;
  derivs(3) = t * dmHuSqOtCoeffdp + Sqr(t) * dmHuSqOt2Coeffdp;
  derivs(4) = KroneckerDelta(gen,3) + t * dmSSqOtCoeffdp + Sqr(t) * dmSSqOt2Coeffdp;
  derivs(5) = t * dmQlSqOtCoeffdp + Sqr(t) * dmQlSqOt2Coeffdp;
  derivs(6) = t * dmUrSqOtCoeffdp + Sqr(t) * dmUrSqOt2Coeffdp;
  derivs(7) = dAtdp;

  return derivs;

}

// Currently only 1-loop contributions are included in the derivatives for A_lambda
Eigen::Matrix<double,8,1> doCalcLambdaDerivs(genericE6SSM_soft_parameters r, double ms, int gen,
				bool & hasError)
{
  if (gen != 3)
    {
      ostringstream ii;
      ii << "WARNING: can only do calculation of derivatives for lambda(3), not lambda(" << gen << ")" << endl;
      throw ii.str(); 
    }

  double mx = r.get_scale();

  // Taylor approximation to the RGE solution is assumed to be about the values at MX.
  double t = log(ms/mx);

  Eigen::Matrix<double,8,1> derivs;

  genericE6SSM_input_parameters input = r.get_input();

  double g1 = r.get_g1();
  double g2 = r.get_g2();
  double g3 = r.get_g3();
  double gN = r.get_gN();

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
  
  double mHd2 = r.get_mHd2();
  double mHu2 = r.get_mHu2();
  double ms2 = r.get_ms2();

  Eigen::Matrix<double,3,3> Yd, Yu, Ye, kappa;  
  Eigen::Matrix<double,2,2> lambda12;

  Yd = r.get_Yd(); Yu = r.get_Yu(); Ye = r.get_Ye();

  kappa = r.get_Kappa();
  lambda12 = r.get_Lambda12();

  double M1 = r.get_MassB();
  double M2 = r.get_MassWB();
  double M3 = r.get_MassG();
  double M1p = r.get_MassBp();  

  double dlambdadp, dAlambdadp, dmHdSqdp, dmHuSqdp, dmSSqdp, dmQlSqdp, dmUrSqdp, dAtdp;


  // Calculate each of the remaining non-zero derivatives:

  // 1-loop betas
  double dbeta_1loop_lambdadp = -0.6 * Sqr(g1) - 3.0 * Sqr(g2) - 2.0 * Sqr(gN)*(Sqr(input.QH1p)+Sqr(input.QH2p) + Sqr(input.QSp))
    + 12.0 * Sqr(lambda) + 3.0 * (Yd*(Yd.adjoint())).trace() + 3.0 * (Yu*(Yu.adjoint())).trace()
    + (Ye*(Ye.adjoint())).trace() + 3.0 * (kappa*(kappa.adjoint())).trace() + 2.0 * (lambda12*(lambda12.adjoint())).trace();
  double dbeta_1loop_Alambdadp = 16.0*lambda*Alambda;
  double dbeta_1loop_mHdSqdp = 4.0 * lambda * (mHd2 + mHu2 + ms2 + Sqr(Alambda) );
  double dbeta_1loop_mHuSqdp = 4.0 * lambda * (mHd2 + mHu2 + ms2 + Sqr(Alambda) );
  double dbeta_1loop_mSSqdp = 8.0 * lambda * (mHd2 + mHu2 + ms2 + Sqr(Alambda) );
  double dbeta_1loop_mQlSqdp = 0.0;
  double dbeta_1loop_mUrSqdp = 0.0;
  // Assume first and second generation, and off-diagonal mixings, all vanish for simplicity
  double dbeta_1loop_Atdp = 4.0*lambda*Alambda;

  Eigen::Matrix<double,3,3> Td = r.get_TYd();
  Eigen::Matrix<double,3,3> Tu = r.get_TYu();
  Eigen::Matrix<double,3,3> Te = r.get_TYe();

  // Extra 1-loop derivatives needed below
  Eigen::Matrix<double,3,3> dbeta_1loop_Yddp = 2.0 * lambda * Yd;
  Eigen::Matrix<double,3,3> dbeta_1loop_kappadp = 4.0 * lambda * kappa;
  Eigen::Matrix<double,2,2> dbeta_1loop_lambda12dp = 4.0 * lambda * lambda12;
  Eigen::Matrix<double,3,3> dbeta_1loop_Yudp = 2.0 * lambda * Yu;
  Eigen::Matrix<double,3,3> dbeta_1loop_Yedp = 2.0 * lambda * Ye; 
  Eigen::Matrix<double,3,3> dbeta_1loop_Tddp = 2.0 * lambda * Td + 4.0 * lambda * Alambda * Yd;
  Eigen::Matrix<double,3,3> dbeta_1loop_Tudp = 2.0 * lambda * Tu + 4.0 * lambda * Alambda * Yu;

  int savedLoops = r.get_loops();
  r.set_loops(1);

  genericE6SSM_soft_parameters w = r.calc_beta();

  r.set_loops(savedLoops);

  double QQp = r.get_input().QQp;
  double Qup = r.get_input().Qup;
  double Qdp = r.get_input().Qdp;
  double QLp = r.get_input().QLp;
  double Qep = r.get_input().Qep;
  double QH1p = r.get_input().QH1p;
  double QH2p = r.get_input().QH2p;
  double QSp = r.get_input().QSp;
  double QHpp = r.get_input().QHpp;
  double QHpbarp = r.get_input().QHpbarp;
  double QDxp = r.get_input().QDxp;
  double QDxbarp = r.get_input().QDxbarp;

  // TODO:: (1-loop)^2 derivatives (which can be expressed in terms of the above)
  double dbeta_1loop2_lambdadp = 0.0;
  double dbeta_1loop2_Alambdadp = 0.0;
  double dbeta_1loop2_mHdSqdp = 0.0;
  double dbeta_1loop2_mHuSqdp = 0.0;
  double dbeta_1loop2_mSSqdp = 0.0;
  double dbeta_1loop2_mQlSqdp = 0.0;
  double dbeta_1loop2_mUrSqdp = 0.0;
  double dbeta_1loop2_Atdp = 0.0;  
  
  int nLogs = 1;

  if (nLogs > 0)
    {
      // Add in leading log contributions
      dbeta_1loop2_lambdadp += (-1.2 * g1 * w.get_g1() - 6.0 * g2 * w.get_g2() - 4.0 * gN * w.get_gN() * (Sqr(QH1p) + Sqr(QH2p) + Sqr(QSp))
				+ 8.0 * lambda * w.get_Lambdax()
				+ 3.0 * (w.get_Yd()*(Yd.adjoint())).trace() + 3.0 * (Yd * (w.get_Yd().adjoint())).trace()
				+ 3.0 * (w.get_Yu()*(Yu.adjoint())).trace() + 3.0 * (Yu * (w.get_Yu().adjoint())).trace()
				+ (w.get_Ye()*(Ye.adjoint())).trace() + (Ye * (w.get_Ye().adjoint())).trace()
				+ 3.0 * (w.get_Kappa()*(kappa.adjoint())).trace() + 3.0 * (kappa * (w.get_Kappa().adjoint())).trace()
				+ 2.0 * (w.get_Lambda12()*(lambda12.adjoint())).trace() + 2.0 * (lambda12 * (w.get_Lambda12().adjoint())).trace())/oneOver16PiSqr;
      dbeta_1loop2_lambdadp += lambda * (8.0 * w.get_Lambdax()/oneOver16PiSqr + 8.0 * lambda * dbeta_1loop_lambdadp
      					 + 3.0 * (dbeta_1loop_Yddp*(Yd.adjoint())).trace() + 3.0 * (Yd*(dbeta_1loop_Yddp.adjoint())).trace()
      					 + 3.0 * (dbeta_1loop_Yudp*(Yu.adjoint())).trace() + 3.0 * (Yu*(dbeta_1loop_Yudp.adjoint())).trace()
      					 + (dbeta_1loop_Yedp*(Ye.adjoint())).trace() + (Ye*(dbeta_1loop_Yedp.adjoint())).trace()
      					 + 3.0 * (dbeta_1loop_kappadp*(kappa.adjoint())).trace() + 3.0 * (kappa*(dbeta_1loop_kappadp.adjoint())).trace()
      					 + 2.0 * (dbeta_1loop_lambda12dp*(lambda12.adjoint())).trace() + 2.0 * (lambda12*(dbeta_1loop_lambda12dp.adjoint())).trace());
      dbeta_1loop2_lambdadp += dbeta_1loop_lambdadp * (-0.6*g1*g1-3.0*g2*g2-2.0*gN*gN*(Sqr(QH1p)+Sqr(QH2p)+Sqr(QSp))+4.0*lambda*lambda
      						       +3.0*(Yd*(Yd.adjoint())).trace()+3.0*(Yu*(Yu.adjoint())).trace()
      						       +(Ye*(Ye.adjoint())).trace()+3.0*(kappa*(kappa.adjoint())).trace()
      						       +2.0*(lambda12*(lambda12.adjoint())).trace());
      dbeta_1loop2_lambdadp += 8.0 * lambda * w.get_Lambdax()/oneOver16PiSqr;
      
      dbeta_1loop2_mQlSqdp += 2.0 * dbeta_1loop_mHdSqdp * (Yd.adjoint() * Yd)(2,2) + 2.0 * mHd2 * (dbeta_1loop_Yddp.adjoint() * Yd)(2,2)
	+ 2.0 * mHd2 * (Yd.adjoint() * dbeta_1loop_Yddp)(2,2) + 2.0 * dbeta_1loop_mHuSqdp * (Yu.adjoint() * Yu)(2,2)
	+ 2.0 * mHu2 * (dbeta_1loop_Yudp.adjoint() * Yu)(2,2) + 2.0 * mHu2 * (Yu.adjoint() * dbeta_1loop_Yudp)(2,2)
	+ 2.0 * (dbeta_1loop_Tddp.adjoint() * Td)(2,2) + 2.0 * (Td.adjoint() * dbeta_1loop_Tddp)(2,2)
	+ 2.0 * (dbeta_1loop_Tudp.adjoint() * Tu)(2,2) + 2.0 * (Tu.adjoint() * dbeta_1loop_Tudp)(2,2)
	+ (r.get_mq2()*(dbeta_1loop_Yudp.adjoint())*Yu)(2,2) + (r.get_mq2()*(Yu.adjoint())*dbeta_1loop_Yudp)(2,2)
	+ (r.get_mq2()*(dbeta_1loop_Yddp.adjoint())*Yd)(2,2) + (r.get_mq2()*(Yd.adjoint())*dbeta_1loop_Yddp)(2,2)
	+ 2.0 * (dbeta_1loop_Yddp.adjoint()*r.get_md2()*Yd)(2,2) + 2.0 * (Yd.adjoint()*r.get_md2()*dbeta_1loop_Yddp)(2,2)
	+ (dbeta_1loop_Yddp.adjoint()*Yd*r.get_mq2())(2,2) +(Yd.adjoint()*dbeta_1loop_Yddp*r.get_mq2())(2,2)
	+ 2.0 * (dbeta_1loop_Yudp.adjoint()*r.get_mu2()*Yu)(2,2) + 2.0 * (Yu.adjoint() * r.get_mu2() * dbeta_1loop_Yudp)(2,2)
	+ (dbeta_1loop_Yudp.adjoint()*Yu*r.get_mq2())(2,2) + (Yu.adjoint()*dbeta_1loop_Yudp*r.get_mq2())(2,2);


    }

  // TODO:: 2-loop derivatives
  double dbeta_2loop_lambdadp = 0.0;
  double dbeta_2loop_Alambdadp = 0.0;
  double dbeta_2loop_mHdSqdp = 0.0;
  double dbeta_2loop_mHuSqdp = 0.0;
  double dbeta_2loop_mSSqdp = 0.0;
  double dbeta_2loop_mQlSqdp = 0.0;
  double dbeta_2loop_mUrSqdp = 0.0;
  double dbeta_2loop_Atdp = 0.0;  
  
  if (r.get_loops() > 1)
    {
      // Add in 2-loop contributions
      dbeta_2loop_lambdadp += -0.02 * (-297.0*Sqr(g1*g1)-90.0*Sqr(g1*g2)-825.0*Sqr(g2*g2)+180.0*Sqr(g1*gN)*Qdp*QH1p
      			      	       +180.0*Sqr(g1*gN)*QDxbarp*QH1p-180.0*Sqr(g1*gN)*QDxp*QH1p+180.0*Sqr(g1*gN)*Qep*QH1p
				       -240.0*Sqr(g1*gN*QH1p)-300.0*Sqr(g2*gN*QH1p)-900.0*Sqr(gN*gN*Qdp*QH1p)
				       -900.0*Sqr(gN*gN*QDxbarp*QH1p)-900.0*Sqr(gN*gN*QDxp*QH1p)-300.0*Sqr(gN*gN*Qep*QH1p)
				       -800.0*Sqr(gN*gN*QH1p*QH1p)-180.0*Sqr(g1*gN)*Qdp*QH2p-180.0*Sqr(g1*gN)*QDxbarp*QH2p
				       +180.0*Sqr(g1*gN)*QDxp*QH2p-180.0*Sqr(g1*gN)*Qep*QH2p+360.0*Sqr(g1*gN)*QH1p*QH2p
				       -240.0*Sqr(g1*gN*QH2p)-300.0*Sqr(g2*gN*QH2p)-900.0*Sqr(gN*gN*Qdp*QH2p)
				       -900.0*Sqr(gN*gN*QDxbarp*QH2p)-900.0*Sqr(gN*gN*QDxp*QH2p)-300.0*Sqr(gN*gN*Qep*QH2p)
				       -1200.0*Sqr(gN*gN*QH1p*QH2p)-800.0*Sqr(gN*gN*QH2p*QH2p)+60.0*Sqr(g1*gN)*QH1p*QHpbarp
				       -60.0*Sqr(g1*gN)*QH2p*QHpbarp-200.0*Sqr(gN*gN*QH1p*QHpbarp)-200.0*Sqr(gN*gN*QH2p*QHpbarp)
				       -60.0*Sqr(g1*gN)*QH1p*QHpp+60.0*Sqr(g1*gN)*QH2p*QHpp-200.0*Sqr(gN*gN*QH1p*QHpp)
				       -200.0*Sqr(gN*gN*QH2p*QHpp)-180.0*Sqr(g1*gN)*QH1p*QLp+180.0*Sqr(g1*gN)*QH2p*QLp
				       -600.0*Sqr(gN*gN*QH1p*QLp)-600.0*Sqr(gN*gN*QH2p*QLp)+180.0*Sqr(g1*gN)*QH1p*QQp
				       -180.0*Sqr(g1*gN)*QH2p*QQp-1800.0*Sqr(gN*gN*QH1p*QQp)-1800.0*Sqr(gN*gN*QH2p*QQp)
				       -900.0*Sqr(gN*gN*Qdp*QSp)-900.0*Sqr(gN*gN*QDxbarp*QSp)-900.0*Sqr(gN*gN*QDxp*QSp)
				       -300.0*Sqr(gN*gN*Qep*QSp)-900.0*Sqr(gN*gN*QH1p*QSp)-900.0*Sqr(gN*gN*QH2p*QSp)
				       -200.0*Sqr(gN*gN*QHpbarp*QSp)-200.0*Sqr(gN*gN*QHpp*QSp)-600.0*Sqr(gN*gN*QLp*QSp)
				       -1800.0*Sqr(gN*gN*QQp*QSp)-500.0*Sqr(gN*gN*QSp*QSp)-360.0*Sqr(g1*gN)*QH1p*Qup
				       +360.0*Sqr(g1*gN)*QH2p*Qup-900.0*Sqr(gN*gN*QH1p*Qup)-900.0*Sqr(gN*gN*QH2p*Qup)
				       -900.0*Sqr(gN*gN*QSp*Qup)+500.0*Sqr(lambda*lambda)
				       +20.0*(-5.0*(3.0*Sqr(gN)*(-Sqr(QH1p)+Sqr(Qdp)+Sqr(QQp))+8.0*Sqr(g3))+Sqr(g1))*(Yd*(Yd.adjoint())).trace()
				       -20.0*(3.0*Sqr(g1)+5.0*Sqr(gN)*(Sqr(Qep)+Sqr(QLp)-Sqr(QH1p)))*(Ye*(Ye.adjoint())).trace()
				       -20.0*(2.0*Sqr(g1)+40.0*Sqr(g3)+15.0*Sqr(gN)*(Sqr(QQp)-Sqr(QH2p)+Sqr(Qup)))*(Yu*(Yu.adjoint())).trace()
				       -20.0*(2.0*Sqr(g1)+40.0*Sqr(g3)+15.0*Sqr(gN)*(Sqr(QDxbarp)+Sqr(QDxp)-Sqr(QSp)))*(kappa*(kappa.adjoint())).trace()
				       -10.0*Sqr(lambda)*(6.0*Sqr(g1)+30.0*Sqr(g2)+20.0*Sqr(gN)*(Sqr(QH1p)+Sqr(QH2p))
							  -45.0*(Yd*(Yd.adjoint())).trace()-15.0*(Ye*(Ye.adjoint())).trace()
							  -45.0*(Yu*(Yu.adjoint())).trace()-30.0*(kappa*(kappa.adjoint())).trace()
							  -20.0*(lambda12*(lambda12.adjoint())).trace())
				       -20.0*(3.0*Sqr(g1)+15.0*Sqr(g2)+10.0*Sqr(gN)*(Sqr(QH1p)+Sqr(QH2p)-Sqr(QSp)))*(lambda12*(lambda12.adjoint())).trace()
				       +450.0*(Yd*(Yd.adjoint())*Yd*(Yd.adjoint())).trace()
				       +300.0*(Yd*(Yu.adjoint())*Yu*(Yd.adjoint())).trace()
				       +150.0*(Ye*(Ye.adjoint())*Ye*(Ye.adjoint())).trace()
				       +450.0*(Yu*(Yu.adjoint())*Yu*(Yu.adjoint())).trace()
				       +300.0*(kappa*(kappa.adjoint())*kappa*(kappa.adjoint())).trace()
				       +200.0*(lambda12*(lambda12.adjoint())*lambda12*(lambda12.adjoint())).trace());
      dbeta_2loop_lambdadp += -0.02 * lambda * (2000.0*Sqr(lambda)*lambda
						-20.0*lambda*(6.0*Sqr(g1)+30.0*Sqr(g2)+20.0*Sqr(gN)*(Sqr(QH1p)+Sqr(QH2p))
							  -45.0*(Yd*(Yd.adjoint())).trace()-15.0*(Ye*(Ye.adjoint())).trace()
							  -45.0*(Yu*(Yu.adjoint())).trace()-30.0*(kappa*(kappa.adjoint())).trace()
							  -20.0*(lambda12*(lambda12.adjoint())).trace()));

      dbeta_2loop_Alambdadp = 0.0;
      dbeta_2loop_mHdSqdp = 0.0;
      dbeta_2loop_mHuSqdp = 0.0;
      dbeta_2loop_mSSqdp = -0.8 * (-6.0 * Sqr(g1) * lambda * Sqr(Alambda) - 30.0 * Sqr(g2) * lambda * Sqr(Alambda)
				   -20.0 * Sqr(gN) * lambda * Sqr(Alambda) * (Sqr(QH1p) + Sqr(QH2p) - Sqr(QSp))
				   + 80.0 * Sqr(lambda)*lambda * (mHd2 + mHu2 + ms2) + 6.0 * Sqr(g1) * M1 * lambda * Alambda
				   + 30.0 * Sqr(g2) * M2 * lambda * Alambda 
				   + 20.0 * Sqr(gN) * M1p * lambda * Alambda * (Sqr(QH1p) + Sqr(QH2p) - Sqr(QSp))
				   - Sqr(gN) * QSp * (-20.0 * lambda * (mHd2 * QH1p + mHu2 * QH2p + ms2 * QSp))
				   + 30.0 * lambda * Sqr(Alambda) * (Yd*(Yd.adjoint())).trace()
				   + 10.0 * lambda *Sqr(Alambda) * (Ye*(Ye.adjoint())).trace()
				   + 30.0 * lambda * Sqr(Alambda) * (Yu*(Yu.adjoint())).trace()
				   + 30.0 * lambda * Alambda * (Yd.adjoint()*Td).trace()
				   + 10.0 * lambda * Alambda * (Ye.adjoint()*Te).trace()
				   + 30.0 * lambda * Alambda *(Yu.adjoint()*Tu).trace()
				   - 5.0 * Sqr(gN) * M1p * (2.0 * (Sqr(QH1p) + Sqr(QH2p)-Sqr(QSp))*(2.0 * M1p * lambda - Tlambda)
							    + 2.0 * (Sqr(QH1p) + Sqr(QH2p) - Sqr(QSp))*lambda*(2.0*M1p-Alambda))
				   + 160.0 * Sqr(lambda*Alambda)*lambda
				   + 2.0 * lambda * (-3.0 * Sqr(g1) * mHd2 -15.0 * Sqr(g2) * mHd2 - 3.0 * Sqr(g1) * mHu2
						     -15.0 * Sqr(g2) * mHu2 - 3.0 * Sqr(g1) * ms2 - 15.0 * Sqr(g2) * ms2
						     -10.0*Sqr(gN)*(Sqr(QH1p)+Sqr(QH2p)-Sqr(QSp))*(mHd2+mHu2+ms2)
						     +3.0*Sqr(g1)*M1*(-2.0*M1+Alambda)+15.0*Sqr(g2)*M2*(-2.0*M2+Alambda)
						     +30.0*mHd2*(Yd*(Yd.adjoint())).trace()+15.0*mHu2*(Yd*(Yd.adjoint())).trace()
						     +15.0*ms2*(Yd*(Yd.adjoint())).trace()+10.0*mHd2*(Ye*(Ye.adjoint())).trace()
						     +5.0*mHu2*(Ye*(Ye.adjoint())).trace()+5.0*ms2*(Ye*(Ye.adjoint())).trace()
						     +15.0*mHd2*(Yu*(Yu.adjoint())).trace()+30.0*mHu2*(Yu*(Yu.adjoint())).trace()
						     +15.0*ms2*(Yu*(Yu.adjoint())).trace()+15.0*Alambda*(Yd*(Td.adjoint())).trace()
						     +15.0*(Td*(Td.adjoint())).trace()+5.0*(Te*(Te.adjoint())).trace()
						     +5.0*Alambda*(Ye*(Te.adjoint())).trace()+15.0*Alambda*(Yu*(Tu.adjoint())).trace()
						     +15.0*(Tu*(Tu.adjoint())).trace()+15.0*(r.get_md2()*Yd*(Yd.adjoint())).trace()
						     +5.0*(r.get_me2()*Ye*(Ye.adjoint())).trace()+5.0*(r.get_ml2()*(Ye.adjoint())*Ye).trace()
						     +15.0*(r.get_mq2()*(Yd.adjoint())*Yd).trace()+15.0*(r.get_mq2()*(Yu.adjoint())*Yu).trace()
						     +15.0*(r.get_mu2()*Yu*(Yu.adjoint())).trace()));
      dbeta_2loop_mQlSqdp += -8.0*lambda*mHd2*(Yd.adjoint()*Yd)(2,2)-4.0*mHu2*lambda*(Yd.adjoint()*Yd)(2,2)
	-4.0*lambda*ms2*(Yd.adjoint()*Yd)(2,2)-4.0*lambda*Sqr(Alambda)*(Yd.adjoint()*Yd)(2,2)
	-4.0*lambda*Alambda*(Yd.adjoint()*r.get_TYd())(2,2)-4.0*lambda*mHd2*(Yu.adjoint()*Yu)(2,2)
	-8.0*lambda*mHu2*(Yu.adjoint()*Yu)(2,2)-4.0*lambda*ms2*(Yu.adjoint()*Yu)(2,2)-4.0*lambda*Sqr(Alambda)*(Yu.adjoint()*Yu)(2,2)
	-4.0*lambda*Alambda*(Yu.adjoint()*r.get_TYu())(2,2)-4.0*lambda*(r.get_TYd().adjoint()*r.get_TYd())(2,2)
	-4.0*lambda*(r.get_TYu().adjoint()*r.get_TYu())(2,2)-2.0*lambda*(r.get_mq2()*Yd.adjoint()*Yd)(2,2)
	-2.0*lambda*(r.get_mq2()*Yu.adjoint()*Yu)(2,2)-4.0*lambda*(Yd.adjoint()*r.get_md2()*Yd)(2,2)
	-2.0*lambda*(Yd.adjoint()*Yd*r.get_mq2())(2,2)-4.0*lambda*(Yu.adjoint()*r.get_mu2()*Yu)(2,2)
	-2.0*lambda*(Yu.adjoint()*Yu*r.get_mq2())(2,2)-4.0*lambda*Alambda*(r.get_TYd().adjoint()*Yd)(2,2)
	-4.0*lambda*Alambda*(r.get_TYu().adjoint()*Yu)(2,2)+(1.0/75.0)*Sqr(g1)*(60.0*lambda*(mHd2-mHu2))
	+0.8*QQp*Sqr(gN)*(-20.0*lambda*(mHd2*QH1p+mHu2*QH2p+ms2*QSp));
      dbeta_2loop_mUrSqdp = 0.0;
      dbeta_2loop_Atdp = 0.0;  
				       
    }

  derivs(0) = 1.0 + t * (oneOver16PiSqr * dbeta_1loop_lambdadp + Sqr(oneOver16PiSqr) * dbeta_2loop_lambdadp) 
    + 0.5 * Sqr(t * oneOver16PiSqr) * dbeta_1loop2_lambdadp;
  derivs(1) = t * (oneOver16PiSqr * dbeta_1loop_Alambdadp + Sqr(oneOver16PiSqr) * dbeta_2loop_Alambdadp) 
    + 0.5 * Sqr(t * oneOver16PiSqr) * dbeta_1loop2_Alambdadp;   
  derivs(2) = t * (oneOver16PiSqr * dbeta_1loop_mHdSqdp + Sqr(oneOver16PiSqr) * dbeta_2loop_mHdSqdp) 
    + 0.5 * Sqr(t * oneOver16PiSqr) * dbeta_1loop2_mHdSqdp;
  derivs(3) = t * (oneOver16PiSqr * dbeta_1loop_mHuSqdp + Sqr(oneOver16PiSqr) * dbeta_2loop_mHuSqdp) 
    + 0.5 * Sqr(t * oneOver16PiSqr) * dbeta_1loop2_mHuSqdp;
  derivs(4) = t * (oneOver16PiSqr * dbeta_1loop_mSSqdp + Sqr(oneOver16PiSqr) * dbeta_2loop_mSSqdp) 
    + 0.5 * Sqr(t * oneOver16PiSqr) * dbeta_1loop2_mSSqdp;
  derivs(5) = t * (oneOver16PiSqr * dbeta_1loop_mQlSqdp + Sqr(oneOver16PiSqr) * dbeta_2loop_mQlSqdp) 
    + 0.5 * Sqr(t * oneOver16PiSqr) * dbeta_1loop2_mQlSqdp;
  derivs(6) = t * (oneOver16PiSqr * dbeta_1loop_mUrSqdp + Sqr(oneOver16PiSqr) * dbeta_2loop_mUrSqdp) 
    + 0.5 * Sqr(t * oneOver16PiSqr) * dbeta_1loop2_mUrSqdp;
  derivs(7) = t * (oneOver16PiSqr * dbeta_1loop_Atdp + Sqr(oneOver16PiSqr) * dbeta_2loop_Atdp) 
    + 0.5 * Sqr(t * oneOver16PiSqr) * dbeta_1loop2_Atdp;	
  
  return derivs;
}

Eigen::Matrix<double,8,1> doCalcGauginoDerivs(genericE6SSM_soft_parameters r, double ms, int whichGaugino,
				 bool & hasError)
{
  if (whichGaugino < 1 || whichGaugino > 4)
    {
      ostringstream ii;
      ii << "WARNING: can only do calculation of derivatives for M_1, M_2, M_3 and M_1'" << endl;
      throw ii.str(); 
    }

  double mx = r.get_scale();

  // Taylor approximation to the RGE solution is assumed to be about the values at MX.
  double t = log(ms/mx);

  Eigen::Matrix<double,8,1> derivs;

  genericE6SSM_input_parameters input = r.get_input();

  double g1 = r.get_g1();
  double g2 = r.get_g2();
  double g3 = r.get_g3();
  double gN = r.get_gN();

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
  
  double mHd2 = r.get_mHd2();
  double mHu2 = r.get_mHu2();
  double ms2 = r.get_ms2();

  Eigen::Matrix<double,3,3> Yd, Yu, Ye, TYd, TYu, TYe;  

  Yd = r.get_Yd(); Yu = r.get_Yu(); Ye = r.get_Ye();
  TYd = r.get_TYd(); TYu = r.get_TYu(); TYe = r.get_TYe();  

  double M1 = r.get_MassB();
  double M2 = r.get_MassWB();
  double M3 = r.get_MassG();
  double M1p = r.get_MassBp();

  double dlambdadp, dAlambdadp, dmHdSqdp, dmHuSqdp, dmSSqdp, dmQlSqdp, dmUrSqdp, dAtdp;

  switch(whichGaugino)
    {

    case 1:
      {

	// 1-loop betas
	double dbeta_1loop_lambdadp = 0.0;
	double dbeta_1loop_Alambdadp = 1.2 * Sqr(g1);
	double dbeta_1loop_mHdSqdp = -2.4 * Sqr(g1)*M1;
	double dbeta_1loop_mHuSqdp = -2.4 * Sqr(g1)*M1;
	double dbeta_1loop_mSSqdp = 0.0;
	double dbeta_1loop_mQlSqdp = -(4.0/15.0)*Sqr(g1)*M1;
	double dbeta_1loop_mUrSqdp = -(64.0/15.0)*Sqr(g1)*M1;
	// Assume first and second generation, and off-diagonal mixings, all vanish for simplicity
	double dbeta_1loop_Atdp = (26.0/15.0)*Sqr(g1);
	
	// TODO:: (1-loop)^2 derivatives (which can be expressed in terms of the above)
	double dbeta_1loop2_lambdadp = 0.0;
	double dbeta_1loop2_Alambdadp = 0.0;
	double dbeta_1loop2_mHdSqdp = 0.0;
	double dbeta_1loop2_mHuSqdp = 0.0;
	double dbeta_1loop2_mSSqdp = 0.0;
	double dbeta_1loop2_mQlSqdp = 0.0;
	double dbeta_1loop2_mUrSqdp = 0.0;
	double dbeta_1loop2_Atdp = 0.0;  

	int nLogs = 1;
	
	if (nLogs > 0)
	  {
	    // Add in leading log contributions
	  }
	
	// TODO:: 2-loop derivatives
	double dbeta_2loop_lambdadp = 0.0;
	double dbeta_2loop_Alambdadp = 0.0;
	double dbeta_2loop_mHdSqdp = 0.0;
	double dbeta_2loop_mHuSqdp = 0.0;
	double dbeta_2loop_mSSqdp = 0.0;
	double dbeta_2loop_mQlSqdp = 0.0;
	double dbeta_2loop_mUrSqdp = 0.0;
	double dbeta_2loop_Atdp = 0.0;  

	if (r.get_loops() > 1)
	  {
	    // Add in 2-loop contributions
	  }
	
	derivs(0) = t * (oneOver16PiSqr * dbeta_1loop_lambdadp + Sqr(oneOver16PiSqr) * dbeta_2loop_lambdadp) 
	  + 0.5 * Sqr(t * oneOver16PiSqr) * dbeta_1loop2_lambdadp;
	derivs(1) = t * (oneOver16PiSqr * dbeta_1loop_Alambdadp + Sqr(oneOver16PiSqr) * dbeta_2loop_Alambdadp)
	  + 0.5 * Sqr(t * oneOver16PiSqr) * dbeta_1loop2_Alambdadp;   
	derivs(2) = t * (oneOver16PiSqr * dbeta_1loop_mHdSqdp + Sqr(oneOver16PiSqr) * dbeta_2loop_mHdSqdp) 
	  + 0.5 * Sqr(t * oneOver16PiSqr) * dbeta_1loop2_mHdSqdp;
	derivs(3) = t * (oneOver16PiSqr * dbeta_1loop_mHuSqdp + Sqr(oneOver16PiSqr) * dbeta_2loop_mHuSqdp) 
	  + 0.5 * Sqr(t * oneOver16PiSqr) * dbeta_1loop2_mHuSqdp;
	derivs(4) = t * (oneOver16PiSqr * dbeta_1loop_mSSqdp + Sqr(oneOver16PiSqr) * dbeta_2loop_mSSqdp) 
	  + 0.5 * Sqr(t * oneOver16PiSqr) * dbeta_1loop2_mSSqdp;
	derivs(5) = t * (oneOver16PiSqr * dbeta_1loop_mQlSqdp + Sqr(oneOver16PiSqr) * dbeta_2loop_mQlSqdp) 
	  + 0.5 * Sqr(t * oneOver16PiSqr) * dbeta_1loop2_mQlSqdp;
	derivs(6) = t * (oneOver16PiSqr * dbeta_1loop_mUrSqdp + Sqr(oneOver16PiSqr) * dbeta_2loop_mUrSqdp) 
	  + 0.5 * Sqr(t * oneOver16PiSqr) * dbeta_1loop2_mUrSqdp;
	derivs(7) = t * (oneOver16PiSqr * dbeta_1loop_Atdp + Sqr(oneOver16PiSqr) * dbeta_2loop_Atdp)
	  + 0.5 * Sqr(t * oneOver16PiSqr) * dbeta_1loop2_Atdp;	

	break;		
      }
    case 2:
      {
	// 1-loop betas
	double dbeta_1loop_lambdadp = 0.0;
	double dbeta_1loop_Alambdadp = 6.0 * Sqr(g2);
	double dbeta_1loop_mHdSqdp = -12.0 * Sqr(g2) * M2;
	double dbeta_1loop_mHuSqdp = -12.0 * Sqr(g2) * M2;
	double dbeta_1loop_mSSqdp = 0.0;
	double dbeta_1loop_mQlSqdp = -12.0 * Sqr(g2) * M2;
	double dbeta_1loop_mUrSqdp = -12.0 * Sqr(g2) * M2;
	// Assume first and second generation, and off-diagonal mixings, all vanish for simplicity
	double dbeta_1loop_Atdp = 6.0 * Sqr(g2);

	// TODO:: (1-loop)^2 derivatives (which can be expressed in terms of the above)
	double dbeta_1loop2_lambdadp = 0.0;
	double dbeta_1loop2_Alambdadp = 0.0;
	double dbeta_1loop2_mHdSqdp = 0.0;
	double dbeta_1loop2_mHuSqdp = 0.0;
	double dbeta_1loop2_mSSqdp = 0.0;
	double dbeta_1loop2_mQlSqdp = 0.0;
	double dbeta_1loop2_mUrSqdp = 0.0;
	double dbeta_1loop2_Atdp = 0.0;  
	
	int nLogs = 1;
	
	if (nLogs > 0)
	  {
	    // Add in leading log contributions
	  }
	
	// TODO:: 2-loop derivatives
	double dbeta_2loop_lambdadp = 0.0;
	double dbeta_2loop_Alambdadp = 0.0;
	double dbeta_2loop_mHdSqdp = 0.0;
	double dbeta_2loop_mHuSqdp = 0.0;
	double dbeta_2loop_mSSqdp = 0.0;
	double dbeta_2loop_mQlSqdp = 0.0;
	double dbeta_2loop_mUrSqdp = 0.0;
	double dbeta_2loop_Atdp = 0.0;  
	
	if (r.get_loops() > 1)
	  {
	    // Add in 2-loop contributions
	  }

	derivs(0) = t * (oneOver16PiSqr * dbeta_1loop_lambdadp + Sqr(oneOver16PiSqr) * dbeta_2loop_lambdadp) 
	  + 0.5 * Sqr(t * oneOver16PiSqr) * dbeta_1loop2_lambdadp;
	derivs(1) = t * (oneOver16PiSqr * dbeta_1loop_Alambdadp + Sqr(oneOver16PiSqr) * dbeta_2loop_Alambdadp) 
	  + 0.5 * Sqr(t * oneOver16PiSqr) * dbeta_1loop2_Alambdadp;   
	derivs(2) = t * (oneOver16PiSqr * dbeta_1loop_mHdSqdp + Sqr(oneOver16PiSqr) * dbeta_2loop_mHdSqdp) 
	  + 0.5 * Sqr(t * oneOver16PiSqr) * dbeta_1loop2_mHdSqdp;
	derivs(3) = t * (oneOver16PiSqr * dbeta_1loop_mHuSqdp + Sqr(oneOver16PiSqr) * dbeta_2loop_mHuSqdp) 
	  + 0.5 * Sqr(t * oneOver16PiSqr) * dbeta_1loop2_mHuSqdp;
	derivs(4) = t * (oneOver16PiSqr * dbeta_1loop_mSSqdp + Sqr(oneOver16PiSqr) * dbeta_2loop_mSSqdp) 
	  + 0.5 * Sqr(t * oneOver16PiSqr) * dbeta_1loop2_mSSqdp;
	derivs(5) = t * (oneOver16PiSqr * dbeta_1loop_mQlSqdp + Sqr(oneOver16PiSqr) * dbeta_2loop_mQlSqdp) 
	  + 0.5 * Sqr(t * oneOver16PiSqr) * dbeta_1loop2_mQlSqdp;
	derivs(6) = t * (oneOver16PiSqr * dbeta_1loop_mUrSqdp + Sqr(oneOver16PiSqr) * dbeta_2loop_mUrSqdp) 
	  + 0.5 * Sqr(t * oneOver16PiSqr) * dbeta_1loop2_mUrSqdp;
	derivs(7) = t * (oneOver16PiSqr * dbeta_1loop_Atdp + Sqr(oneOver16PiSqr) * dbeta_2loop_Atdp) 
	  + 0.5 * Sqr(t * oneOver16PiSqr) * dbeta_1loop2_Atdp;

	break;
      }
    case 3:
      {
	// 1-loop betas
	double dbeta_1loop_lambdadp = 0.0;
	double dbeta_1loop_Alambdadp = 0.0;
	double dbeta_1loop_mHdSqdp = 0.0;
	double dbeta_1loop_mHuSqdp = 0.0;
	double dbeta_1loop_mSSqdp = 0.0;
	double dbeta_1loop_mQlSqdp = -(64.0/3.0)*Sqr(g3)*M3;
	double dbeta_1loop_mUrSqdp = -(64.0/3.0)*Sqr(g3)*M3;
	// Assume first and second generation, and off-diagonal mixings, all vanish for simplicity
	double dbeta_1loop_Atdp =  (32.0/3.0) * Sqr(g3);
	
	// TODO:: (1-loop)^2 derivatives (which can be expressed in terms of the above)
	double dbeta_1loop2_lambdadp = 0.0;
	double dbeta_1loop2_Alambdadp = 0.0;
	double dbeta_1loop2_mHdSqdp = 0.0;
	double dbeta_1loop2_mHuSqdp = 0.0;
	double dbeta_1loop2_mSSqdp = 0.0;
	double dbeta_1loop2_mQlSqdp = 0.0;
	double dbeta_1loop2_mUrSqdp = 0.0;
	double dbeta_1loop2_Atdp = 0.0;  

	int nLogs = 1;
	
	if (nLogs > 0)
	  {
	    // Add in leading log contributions
	  }
	
	// TODO:: 2-loop derivatives
	double dbeta_2loop_lambdadp = 0.0;
	double dbeta_2loop_Alambdadp = 0.0;
	double dbeta_2loop_mHdSqdp = 0.0;
	double dbeta_2loop_mHuSqdp = 0.0;
	double dbeta_2loop_mSSqdp = 0.0;
	double dbeta_2loop_mQlSqdp = 0.0;
	double dbeta_2loop_mUrSqdp = 0.0;
	double dbeta_2loop_Atdp = 0.0;  
	
	if (r.get_loops() > 1)
	  {
	    // Add in 2-loop contributions
	  }

	derivs(0) = t * (oneOver16PiSqr * dbeta_1loop_lambdadp + Sqr(oneOver16PiSqr) * dbeta_2loop_lambdadp) 
	  + 0.5 * Sqr(t * oneOver16PiSqr) * dbeta_1loop2_lambdadp;
	derivs(1) = t * (oneOver16PiSqr * dbeta_1loop_Alambdadp + Sqr(oneOver16PiSqr) * dbeta_2loop_Alambdadp) 
	  + 0.5 * Sqr(t * oneOver16PiSqr) * dbeta_1loop2_Alambdadp;   
	derivs(2) = t * (oneOver16PiSqr * dbeta_1loop_mHdSqdp + Sqr(oneOver16PiSqr) * dbeta_2loop_mHdSqdp)
	  + 0.5 * Sqr(t * oneOver16PiSqr) * dbeta_1loop2_mHdSqdp;
	derivs(3) = t * (oneOver16PiSqr * dbeta_1loop_mHuSqdp + Sqr(oneOver16PiSqr) * dbeta_2loop_mHuSqdp) 
	  + 0.5 * Sqr(t * oneOver16PiSqr) * dbeta_1loop2_mHuSqdp;
	derivs(4) = t * (oneOver16PiSqr * dbeta_1loop_mSSqdp + Sqr(oneOver16PiSqr) * dbeta_2loop_mSSqdp) 
	  + 0.5 * Sqr(t * oneOver16PiSqr) * dbeta_1loop2_mSSqdp;
	derivs(5) = t * (oneOver16PiSqr * dbeta_1loop_mQlSqdp + Sqr(oneOver16PiSqr) * dbeta_2loop_mQlSqdp)
	  + 0.5 * Sqr(t * oneOver16PiSqr) * dbeta_1loop2_mQlSqdp;
	derivs(6) = t * (oneOver16PiSqr * dbeta_1loop_mUrSqdp + Sqr(oneOver16PiSqr) * dbeta_2loop_mUrSqdp)
	  + 0.5 * Sqr(t * oneOver16PiSqr) * dbeta_1loop2_mUrSqdp;
	derivs(7) = t * (oneOver16PiSqr * dbeta_1loop_Atdp + Sqr(oneOver16PiSqr) * dbeta_2loop_Atdp) 
	  + 0.5 * Sqr(t * oneOver16PiSqr) * dbeta_1loop2_Atdp;

	break;
      }
    case 4:
      {
	// 1-loop betas
	double dbeta_1loop_lambdadp = 0.0;
	double dbeta_1loop_Alambdadp = 4.0 * Sqr(gN) * (Sqr(input.QH1p) + Sqr(input.QH2p) + Sqr(input.QSp));
	double dbeta_1loop_mHdSqdp = -16.0 * Sqr(gN*input.QH1p) * M1p;
	double dbeta_1loop_mHuSqdp = -16.0 * Sqr(gN*input.QH2p) * M1p;
	double dbeta_1loop_mSSqdp = -16.0 * Sqr(gN*input.QSp) * M1p;
	double dbeta_1loop_mQlSqdp = -16.0 * Sqr(gN*input.QQp) * M1p;
	double dbeta_1loop_mUrSqdp = -16.0 * Sqr(gN*input.Qup) * M1p;
	// Assume first and second generation, and off-diagonal mixings, all vanish for simplicity
	double dbeta_1loop_Atdp =  4.0 * Sqr(gN) * (Sqr(input.QQp) + Sqr(input.QH2p) + Sqr(input.Qup));
	
	// TODO:: (1-loop)^2 derivatives (which can be expressed in terms of the above)
	double dbeta_1loop2_lambdadp = 0.0;
	double dbeta_1loop2_Alambdadp = 0.0;
	double dbeta_1loop2_mHdSqdp = 0.0;
	double dbeta_1loop2_mHuSqdp = 0.0;
	double dbeta_1loop2_mSSqdp = 0.0;
	double dbeta_1loop2_mQlSqdp = 0.0;
	double dbeta_1loop2_mUrSqdp = 0.0;
	double dbeta_1loop2_Atdp = 0.0;  
	
	int nLogs = 1;
	
	if (nLogs > 0)
	  {
	    // Add in leading log contributions
	  }
	
	// TODO:: 2-loop derivatives
	double dbeta_2loop_lambdadp = 0.0;
	double dbeta_2loop_Alambdadp = 0.0;
	double dbeta_2loop_mHdSqdp = 0.0;
	double dbeta_2loop_mHuSqdp = 0.0;
	double dbeta_2loop_mSSqdp = 0.0;
	double dbeta_2loop_mQlSqdp = 0.0;
	double dbeta_2loop_mUrSqdp = 0.0;
	double dbeta_2loop_Atdp = 0.0;  
	
	if (r.get_loops() > 1)
	  {
	    // Add in 2-loop contributions
	  }

	derivs(0) = t * (oneOver16PiSqr * dbeta_1loop_lambdadp + Sqr(oneOver16PiSqr) * dbeta_2loop_lambdadp) 
	  + 0.5 * Sqr(t * oneOver16PiSqr) * dbeta_1loop2_lambdadp;
	derivs(1) = t * (oneOver16PiSqr * dbeta_1loop_Alambdadp + Sqr(oneOver16PiSqr) * dbeta_2loop_Alambdadp) 
	  + 0.5 * Sqr(t * oneOver16PiSqr) * dbeta_1loop2_Alambdadp;   
	derivs(2) = t * (oneOver16PiSqr * dbeta_1loop_mHdSqdp + Sqr(oneOver16PiSqr) * dbeta_2loop_mHdSqdp) 
	  + 0.5 * Sqr(t * oneOver16PiSqr) * dbeta_1loop2_mHdSqdp;
	derivs(3) = t * (oneOver16PiSqr * dbeta_1loop_mHuSqdp + Sqr(oneOver16PiSqr) * dbeta_2loop_mHuSqdp) 
	  + 0.5 * Sqr(t * oneOver16PiSqr) * dbeta_1loop2_mHuSqdp;
	derivs(4) = t * (oneOver16PiSqr * dbeta_1loop_mSSqdp + Sqr(oneOver16PiSqr) * dbeta_2loop_mSSqdp) 
	  + 0.5 * Sqr(t * oneOver16PiSqr) * dbeta_1loop2_mSSqdp;
	derivs(5) = t * (oneOver16PiSqr * dbeta_1loop_mQlSqdp + Sqr(oneOver16PiSqr) * dbeta_2loop_mQlSqdp) 
	  + 0.5 * Sqr(t * oneOver16PiSqr) * dbeta_1loop2_mQlSqdp;
	derivs(6) = t * (oneOver16PiSqr * dbeta_1loop_mUrSqdp + Sqr(oneOver16PiSqr) * dbeta_2loop_mUrSqdp)
	  + 0.5 * Sqr(t * oneOver16PiSqr) * dbeta_1loop2_mUrSqdp;
	derivs(7) = t * (oneOver16PiSqr * dbeta_1loop_Atdp + Sqr(oneOver16PiSqr) * dbeta_2loop_Atdp) 
	  + 0.5 * Sqr(t * oneOver16PiSqr) * dbeta_1loop2_Atdp;

	break;
      }
    default:
      {
	ostringstream ii;
	ii << "ERROR: unrecognised soft gaugino mass M_" << whichGaugino << " requested" << endl;
	throw ii.str();
      }
    }

  return derivs;
}

Eigen::Matrix<double,8,1> doCalcSoftAuDerivs(genericE6SSM_soft_parameters r, double ms, int m, int n,
				bool & hasError)
{
  if (m != 3 || n != 3)
    {
      ostringstream ii;
      ii << "WARNING: can only do calculation of derivatives for Au(2,2), not Au(" << m-1 << ", " << n-1 << ")" << endl;
      throw ii.str(); 
    }

  double mx = r.get_scale();

  // Taylor approximation to the RGE solution is assumed to be about the values at MX.
  double t = log(ms/mx);

  Eigen::Matrix<double,8,1> derivs;

  genericE6SSM_input_parameters input = r.get_input();

  double g1 = r.get_g1();
  double g2 = r.get_g2();
  double g3 = r.get_g3();
  double gN = r.get_gN();

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

  Eigen::Matrix<double,3,3> Yd, Yu, Ye, TYu;  

  Yd = r.get_Yd(); Yu = r.get_Yu(); Ye = r.get_Ye();

  TYu = r.get_TYu();
  
  double At;
  double TYt = r.get_TYu(2,2);
  
  if (Abs(TYt) < EPSTOL)
    {
      At = 0.0;
    }
  else if (Abs(Yu(2,2)) < 1.0e-100)
    {
      ostringstream ii;
      ii << "WARNING: trying to calculate A_t where y_t coupling is " <<
	Abs(Yu(2,2)) << endl;
      throw ii.str();
    }
  else
    {
      At = TYt/Yu(2,2);
    }  

  double mHd2 = r.get_mHd2();
  double mHu2 = r.get_mHu2();
  double ms2 = r.get_ms2();

  double M1 = r.get_MassB();
  double M2 = r.get_MassWB();
  double M3 = r.get_MassG();
  double M1p = r.get_MassBp();

  double dlambdadp, dAlambdadp, dmHdSqdp, dmHuSqdp, dmSSqdp, dmQlSqdp, dmUrSqdp, dAtdp;

  // 1-loop betas
  double dbeta_1loop_lambdadp = 0.0;
  double dbeta_1loop_Alambdadp = 6.0 * Sqr (Yu(2,2));
  double dbeta_1loop_mHdSqdp = 0.0;
  double dbeta_1loop_mHuSqdp = 12.0 * Sqr(Yu(2,2)) * At;
  double dbeta_1loop_mSSqdp = 0.0;
  double dbeta_1loop_mQlSqdp = 4.0 * Sqr(Yu(2,2)) * At;
  double dbeta_1loop_mUrSqdp = 8.0 * Sqr(Yu(2,2)) * At;
  // Assume first and second generation, and off-diagonal mixings, all vanish for simplicity
  double dbeta_1loop_Atdp =  12.0 * Sqr(Yu(2,2));
  
  // extra 1-loop betas needed below
  double dbeta_1loop_g1dp = 0.0;
  double dbeta_1loop_g2dp = 0.0;
  double dbeta_1loop_g3dp = 0.0;
  double dbeta_1loop_gNdp = 0.0;
  double dbeta_1loop_ytdp = 0.0;
  double dbeta_1loop_ybdp = 0.0;
  double dbeta_1loop_ytaudp = 0.0;
  double dbeta_1loop_kappadp = 0.0;
  
  // TODO:: (1-loop)^2 derivatives (which can be expressed in terms of the above)
  double dbeta_1loop2_lambdadp = 0.0;
  double dbeta_1loop2_Alambdadp = 0.0;
  double dbeta_1loop2_mHdSqdp = 0.0;
  double dbeta_1loop2_mHuSqdp = 0.0;
  double dbeta_1loop2_mSSqdp = 0.0;
  double dbeta_1loop2_mQlSqdp = 0.0;
  double dbeta_1loop2_mUrSqdp = 0.0;
  double dbeta_1loop2_Atdp = 0.0;  
  
  int nLogs = 1;
  
  if (nLogs > 0)
    {
	    // Add in leading log contributions
    }	

  // TODO:: 2-loop derivatives
  double dbeta_2loop_lambdadp = 0.0;
  double dbeta_2loop_Alambdadp = 0.0;
  double dbeta_2loop_mHdSqdp = 0.0;
  double dbeta_2loop_mHuSqdp = 0.0;
  double dbeta_2loop_mSSqdp = 0.0;
  double dbeta_2loop_mQlSqdp = 0.0;
  double dbeta_2loop_mUrSqdp = 0.0;
  double dbeta_2loop_Atdp = 0.0;  

  if (r.get_loops() > 1)
    {
      // Add in 2-loop contributions
    }

  derivs(0) = t * (oneOver16PiSqr * dbeta_1loop_lambdadp + Sqr(oneOver16PiSqr) * dbeta_2loop_lambdadp) 
    + 0.5 * Sqr(t * oneOver16PiSqr) * dbeta_1loop2_lambdadp;
  derivs(1) = t * (oneOver16PiSqr * dbeta_1loop_Alambdadp + Sqr(oneOver16PiSqr) * dbeta_2loop_Alambdadp) 
    + 0.5 * Sqr(t * oneOver16PiSqr) * dbeta_1loop2_Alambdadp;   
  derivs(2) = t * (oneOver16PiSqr * dbeta_1loop_mHdSqdp + Sqr(oneOver16PiSqr) * dbeta_2loop_mHdSqdp) 
    + 0.5 * Sqr(t * oneOver16PiSqr) * dbeta_1loop2_mHdSqdp;
  derivs(3) = t * (oneOver16PiSqr * dbeta_1loop_mHuSqdp + Sqr(oneOver16PiSqr) * dbeta_2loop_mHuSqdp) 
    + 0.5 * Sqr(t * oneOver16PiSqr) * dbeta_1loop2_mHuSqdp;
  derivs(4) = t * (oneOver16PiSqr* dbeta_1loop_mSSqdp + Sqr(oneOver16PiSqr) * dbeta_2loop_mSSqdp) 
    + 0.5 * Sqr(t * oneOver16PiSqr) * dbeta_1loop2_mSSqdp;
  derivs(5) = t * (oneOver16PiSqr * dbeta_1loop_mQlSqdp + Sqr(oneOver16PiSqr) * dbeta_2loop_mQlSqdp) 
    + 0.5 * Sqr(t * oneOver16PiSqr) * dbeta_1loop2_mQlSqdp;
  derivs(6) = t * (oneOver16PiSqr * dbeta_1loop_mUrSqdp + Sqr(oneOver16PiSqr) * dbeta_2loop_mUrSqdp) 
    + 0.5 * Sqr(t * oneOver16PiSqr) * dbeta_1loop2_mUrSqdp;
  derivs(7) = 1.0 + t * (oneOver16PiSqr * dbeta_1loop_Atdp + Sqr(oneOver16PiSqr) * dbeta_2loop_Atdp) 
    + 0.5 * Sqr(t * oneOver16PiSqr) * dbeta_1loop2_Atdp;

  return derivs;

}

Eigen::Matrix<double,8,1> doCalcSoftAlambdaDerivs(genericE6SSM_soft_parameters r, double ms, int gen,
				     bool & hasError)
{
  if (gen != 3)
    {
      ostringstream ii;
      ii << "WARNING: can only do calculation of derivatives for Alambda(3), not lambda(" << gen << ")" << endl;
      throw ii.str(); 
    }

  double mx = r.get_scale();

  // Taylor approximation to the RGE solution is assumed to be about the values at MX.
  double t = log(ms/mx);

  Eigen::Matrix<double,8,1> derivs;

  genericE6SSM_input_parameters input = r.get_input();

  double g1 = r.get_g1();
  double g2 = r.get_g2();
  double g3 = r.get_g3();
  double gN = r.get_gN();

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
  
  double mHd2 = r.get_mHd2();
  double mHu2 = r.get_mHu2();
  double ms2 = r.get_ms2();

  Eigen::Matrix<double,3,3> Yd, Yu, Ye, TYu;  

  Yd = r.get_Yd(); Yu = r.get_Yu(); Ye = r.get_Ye();

  TYu = r.get_TYu();  

  double M1 = r.get_MassB();
  double M2 = r.get_MassWB();
  double M3 = r.get_MassG();
  double M1p = r.get_MassBp();

  double dlambdadp, dAlambdadp, dmHdSqdp, dmHuSqdp, dmSSqdp, dmQlSqdp, dmUrSqdp, dAtdp;

  // 1-loop betas
  double dbeta_1loop_lambdadp = 0.0;
  double dbeta_1loop_Alambdadp = 8.0 * Sqr(lambda);
  double dbeta_1loop_mHdSqdp = 4.0 * lambda * Tlambda;
  double dbeta_1loop_mHuSqdp = 4.0 * lambda * Tlambda;
  double dbeta_1loop_mSSqdp = 8.0 * lambda * Tlambda;
  double dbeta_1loop_mQlSqdp = 0.0;
  double dbeta_1loop_mUrSqdp = 0.0;
  // Assume first and second generation, and off-diagonal mixings, all vanish for simplicity
  double dbeta_1loop_Atdp =  2.0 * Sqr(lambda);
  
  // extra 1-loop betas needed below
  double dbeta_1loop_g1dp = 0.0;
  double dbeta_1loop_g2dp = 0.0;
  double dbeta_1loop_g3dp = 0.0;
  double dbeta_1loop_gNdp = 0.0;
  double dbeta_1loop_ytdp = 0.0;
  double dbeta_1loop_ybdp = 0.0;
  double dbeta_1loop_ytaudp = 0.0;
  double dbeta_1loop_kappadp = 0.0;
  
  // TODO:: (1-loop)^2 derivatives (which can be expressed in terms of the above)
  double dbeta_1loop2_lambdadp = 0.0;
  double dbeta_1loop2_Alambdadp = 0.0;
  double dbeta_1loop2_mHdSqdp = 0.0;
  double dbeta_1loop2_mHuSqdp = 0.0;
  double dbeta_1loop2_mSSqdp = 0.0;
  double dbeta_1loop2_mQlSqdp = 0.0;
  double dbeta_1loop2_mUrSqdp = 0.0;
  double dbeta_1loop2_Atdp = 0.0;  
  
  int nLogs = 1;
  
  if (nLogs > 0)
    {
	    // Add in leading log contributions
    }

  // TODO:: 2-loop derivatives
  double dbeta_2loop_lambdadp = 0.0;
  double dbeta_2loop_Alambdadp = 0.0;
  double dbeta_2loop_mHdSqdp = 0.0;
  double dbeta_2loop_mHuSqdp = 0.0;
  double dbeta_2loop_mSSqdp = 0.0;
  double dbeta_2loop_mQlSqdp = 0.0;
  double dbeta_2loop_mUrSqdp = 0.0;
  double dbeta_2loop_Atdp = 0.0;  

  if (r.get_loops() > 1)
    {
      // Add in 2-loop contributions
    }

  derivs(0) = t * (oneOver16PiSqr * dbeta_1loop_lambdadp + Sqr(oneOver16PiSqr) * dbeta_2loop_lambdadp) 
    + 0.5 * Sqr(t * oneOver16PiSqr) * dbeta_1loop2_lambdadp;
  derivs(1) = 1.0 + t * (oneOver16PiSqr * dbeta_1loop_Alambdadp + Sqr(oneOver16PiSqr) * dbeta_2loop_Alambdadp) 
    + 0.5 * Sqr(t * oneOver16PiSqr) * dbeta_1loop2_Alambdadp;   
  derivs(2) = t * (oneOver16PiSqr * dbeta_1loop_mHdSqdp + Sqr(oneOver16PiSqr) * dbeta_2loop_mHdSqdp) 
    + 0.5 * Sqr(t * oneOver16PiSqr) * dbeta_1loop2_mHdSqdp;
  derivs(3) = t * (oneOver16PiSqr * dbeta_1loop_mHuSqdp + Sqr(oneOver16PiSqr) * dbeta_2loop_mHuSqdp) 
    + 0.5 * Sqr(t * oneOver16PiSqr) * dbeta_1loop2_mHuSqdp;
  derivs(4) = t * (oneOver16PiSqr * dbeta_1loop_mSSqdp + Sqr(oneOver16PiSqr) * dbeta_2loop_mSSqdp) 
    + 0.5 * Sqr(t * oneOver16PiSqr) * dbeta_1loop2_mSSqdp;
  derivs(5) = t * (oneOver16PiSqr * dbeta_1loop_mQlSqdp + Sqr(oneOver16PiSqr) * dbeta_2loop_mQlSqdp) 
    + 0.5 * Sqr(t * oneOver16PiSqr) * dbeta_1loop2_mQlSqdp;
  derivs(6) = t * (oneOver16PiSqr * dbeta_1loop_mUrSqdp + Sqr(oneOver16PiSqr) * dbeta_2loop_mUrSqdp) 
    + 0.5 * Sqr(t * oneOver16PiSqr) * dbeta_1loop2_mUrSqdp;
  derivs(7) = t * (oneOver16PiSqr * dbeta_1loop_Atdp + Sqr(oneOver16PiSqr) * dbeta_2loop_Atdp) 
    + 0.5 * Sqr(t * oneOver16PiSqr) * dbeta_1loop2_Atdp;

  return derivs;
}

Eigen::Matrix<double,8,1> doCalcMq2Derivs(genericE6SSM_soft_parameters r, double ms, int m, int n,
			     bool & hasError)
{
  if (m < 1 || m > 3 || n < 1 || n > 3)
    {
      ostringstream ii;
      ii << "WARNING: can only do calculation of derivatives for mq2(i,j), i, j = 1, 2, 3, not mq2(" << m << ", " << n << ")" << endl;
      throw ii.str(); 
    }

  double mx = r.get_scale();

  // Taylor approximation to the RGE solution is assumed to be about the values at MX.
  double t = log(ms/mx);

  Eigen::Matrix<double,8,1> derivs;

  double dlambdadp, dAlambdadp, dmHdSqdp, dmHuSqdp, dmSSqdp, dmQlSqdp, dmUrSqdp, dAtdp;

  // These are zero for all soft squared masses
  dAlambdadp = 0.0;
  dAtdp = 0.0;
  dlambdadp = 0.0;

  // Calculate each of the remaining non-zero derivatives:
  // For derivatives of soft masses w.r.t. soft masses, take advantage
  // of the fact that beta functions are linear in the soft masses,
  // so can calculate derivatives by simple rise/run formula.
  int nLps = r.get_loops();
  int nLogs = 1;
  double shift = 1.0;

  double firstmHdSqLogCoeff = doCalcMh1SquaredLogCoeff(r, nLps);
  double firstmHdSqLogSqCoeff = doCalcMh1SquaredLogSqCoeff(r, nLogs);

  double firstmHuSqLogCoeff = doCalcMh2SquaredLogCoeff(r, nLps);
  double firstmHuSqLogSqCoeff = doCalcMh2SquaredLogSqCoeff(r, nLogs);

  double firstmSSqLogCoeff = doCalcMsSquaredLogCoeff(r, nLps);
  double firstmSSqLogSqCoeff = doCalcMsSquaredLogSqCoeff(r, nLogs);

  double firstmQlSqLogCoeff = doCalcmqL3SquaredLogCoeff(r, nLps);
  double firstmQlSqLogSqCoeff = doCalcmqL3SquaredLogSqCoeff(r, nLogs);

  double firstmUrSqLogCoeff = doCalcmtRSquaredLogCoeff(r, nLps);
  double firstmUrSqLogSqCoeff = doCalcmtRSquaredLogSqCoeff(r, nLogs);

  r.set_mq2(m-1,n-1, r.get_mq2(m-1,n-1)+shift);

  double secondmHdSqLogCoeff = doCalcMh1SquaredLogCoeff(r, nLps);
  double secondmHdSqLogSqCoeff = doCalcMh1SquaredLogSqCoeff(r, nLogs);

  double secondmHuSqLogCoeff = doCalcMh2SquaredLogCoeff(r, nLps);
  double secondmHuSqLogSqCoeff = doCalcMh2SquaredLogSqCoeff(r, nLogs);

  double secondmSSqLogCoeff = doCalcMsSquaredLogCoeff(r, nLps);
  double secondmSSqLogSqCoeff = doCalcMsSquaredLogSqCoeff(r, nLogs);

  double secondmQlSqLogCoeff = doCalcmqL3SquaredLogCoeff(r, nLps);
  double secondmQlSqLogSqCoeff = doCalcmqL3SquaredLogSqCoeff(r, nLogs);

  double secondmUrSqLogCoeff = doCalcmtRSquaredLogCoeff(r, nLps);
  double secondmUrSqLogSqCoeff = doCalcmtRSquaredLogSqCoeff(r, nLogs);

  double dmHdSqOtCoeffdp = (secondmHdSqLogCoeff-firstmHdSqLogCoeff)/shift;
  double dmHdSqOt2Coeffdp = (secondmHdSqLogSqCoeff-firstmHdSqLogSqCoeff)/shift;

  double dmHuSqOtCoeffdp = (secondmHuSqLogCoeff-firstmHuSqLogCoeff)/shift;
  double dmHuSqOt2Coeffdp = (secondmHuSqLogSqCoeff-firstmHuSqLogSqCoeff)/shift;

  double dmSSqOtCoeffdp = (secondmSSqLogCoeff-firstmSSqLogCoeff)/shift;
  double dmSSqOt2Coeffdp = (secondmSSqLogSqCoeff-firstmSSqLogSqCoeff)/shift;

  double dmQlSqOtCoeffdp = (secondmQlSqLogCoeff-firstmQlSqLogCoeff)/shift;
  double dmQlSqOt2Coeffdp = (secondmQlSqLogSqCoeff-firstmQlSqLogSqCoeff)/shift;

  double dmUrSqOtCoeffdp = (secondmUrSqLogCoeff-firstmUrSqLogCoeff)/shift;
  double dmUrSqOt2Coeffdp = (secondmUrSqLogSqCoeff-firstmUrSqLogSqCoeff)/shift;

  derivs(0) = dlambdadp;
  derivs(1) = dAlambdadp;
  derivs(2) = t * dmHdSqOtCoeffdp + Sqr(t) * dmHdSqOt2Coeffdp;
  derivs(3) = t * dmHuSqOtCoeffdp + Sqr(t) * dmHuSqOt2Coeffdp;
  derivs(4) = t * dmSSqOtCoeffdp + Sqr(t) * dmSSqOt2Coeffdp;
  derivs(5) = KroneckerDelta(m,3)*KroneckerDelta(n,3) + t * dmQlSqOtCoeffdp + Sqr(t) * dmQlSqOt2Coeffdp;
  derivs(6) = t * dmUrSqOtCoeffdp + Sqr(t) * dmUrSqOt2Coeffdp;
  derivs(7) = dAtdp;

  return derivs;

}

Eigen::Matrix<double,8,1> doCalcMu2Derivs(genericE6SSM_soft_parameters r, double ms, int m, int n,
			     bool & hasError)
{
  if (m < 1 || m > 3 || n < 1 || n > 3)
    {
      ostringstream ii;
      ii << "WARNING: can only do calculation of derivatives for mu2(i,j), i, j = 1, 2, 3, not mu2(" << m << ", " << n << ")" << endl;
      throw ii.str(); 
    }

  double mx = r.get_scale();

  // Taylor approximation to the RGE solution is assumed to be about the values at MX.
  double t = log(ms/mx);

  Eigen::Matrix<double,8,1> derivs;

  double dlambdadp, dAlambdadp, dmHdSqdp, dmHuSqdp, dmSSqdp, dmQlSqdp, dmUrSqdp, dAtdp;

  // These are zero for all soft squared masses
  dAlambdadp = 0.0;
  dAtdp = 0.0;
  dlambdadp = 0.0;

  // Calculate each of the remaining non-zero derivatives:
  // For derivatives of soft masses w.r.t. soft masses, take advantage
  // of the fact that beta functions are linear in the soft masses,
  // so can calculate derivatives by simple rise/run formula.
  int nLps = r.get_loops();
  int nLogs = 1;
  double shift = 1.0;

  double firstmHdSqLogCoeff = doCalcMh1SquaredLogCoeff(r, nLps);
  double firstmHdSqLogSqCoeff = doCalcMh1SquaredLogSqCoeff(r, nLogs);

  double firstmHuSqLogCoeff = doCalcMh2SquaredLogCoeff(r, nLps);
  double firstmHuSqLogSqCoeff = doCalcMh2SquaredLogSqCoeff(r, nLogs);

  double firstmSSqLogCoeff = doCalcMsSquaredLogCoeff(r, nLps);
  double firstmSSqLogSqCoeff = doCalcMsSquaredLogSqCoeff(r, nLogs);

  double firstmQlSqLogCoeff = doCalcmqL3SquaredLogCoeff(r, nLps);
  double firstmQlSqLogSqCoeff = doCalcmqL3SquaredLogSqCoeff(r, nLogs);

  double firstmUrSqLogCoeff = doCalcmtRSquaredLogCoeff(r, nLps);
  double firstmUrSqLogSqCoeff = doCalcmtRSquaredLogSqCoeff(r, nLogs);

  r.set_mu2(m-1,n-1, r.get_mu2(m-1,n-1)+shift);

  double secondmHdSqLogCoeff = doCalcMh1SquaredLogCoeff(r, nLps);
  double secondmHdSqLogSqCoeff = doCalcMh1SquaredLogSqCoeff(r, nLogs);

  double secondmHuSqLogCoeff = doCalcMh2SquaredLogCoeff(r, nLps);
  double secondmHuSqLogSqCoeff = doCalcMh2SquaredLogSqCoeff(r, nLogs);

  double secondmSSqLogCoeff = doCalcMsSquaredLogCoeff(r, nLps);
  double secondmSSqLogSqCoeff = doCalcMsSquaredLogSqCoeff(r, nLogs);

  double secondmQlSqLogCoeff = doCalcmqL3SquaredLogCoeff(r, nLps);
  double secondmQlSqLogSqCoeff = doCalcmqL3SquaredLogSqCoeff(r, nLogs);

  double secondmUrSqLogCoeff = doCalcmtRSquaredLogCoeff(r, nLps);
  double secondmUrSqLogSqCoeff = doCalcmtRSquaredLogSqCoeff(r, nLogs);

  double dmHdSqOtCoeffdp = (secondmHdSqLogCoeff-firstmHdSqLogCoeff)/shift;
  double dmHdSqOt2Coeffdp = (secondmHdSqLogSqCoeff-firstmHdSqLogSqCoeff)/shift;

  double dmHuSqOtCoeffdp = (secondmHuSqLogCoeff-firstmHuSqLogCoeff)/shift;
  double dmHuSqOt2Coeffdp = (secondmHuSqLogSqCoeff-firstmHuSqLogSqCoeff)/shift;

  double dmSSqOtCoeffdp = (secondmSSqLogCoeff-firstmSSqLogCoeff)/shift;
  double dmSSqOt2Coeffdp = (secondmSSqLogSqCoeff-firstmSSqLogSqCoeff)/shift;

  double dmQlSqOtCoeffdp = (secondmQlSqLogCoeff-firstmQlSqLogCoeff)/shift;
  double dmQlSqOt2Coeffdp = (secondmQlSqLogSqCoeff-firstmQlSqLogSqCoeff)/shift;

  double dmUrSqOtCoeffdp = (secondmUrSqLogCoeff-firstmUrSqLogCoeff)/shift;
  double dmUrSqOt2Coeffdp = (secondmUrSqLogSqCoeff-firstmUrSqLogSqCoeff)/shift;

  derivs(0) = dlambdadp;
  derivs(1) = dAlambdadp;
  derivs(2) = t * dmHdSqOtCoeffdp + Sqr(t) * dmHdSqOt2Coeffdp;
  derivs(3) = t * dmHuSqOtCoeffdp + Sqr(t) * dmHuSqOt2Coeffdp;
  derivs(4) = t * dmSSqOtCoeffdp + Sqr(t) * dmSSqOt2Coeffdp;
  derivs(5) = t * dmQlSqOtCoeffdp + Sqr(t) * dmQlSqOt2Coeffdp;
  derivs(6) = KroneckerDelta(m,3)*KroneckerDelta(n,3) + t * dmUrSqOtCoeffdp + Sqr(t) * dmUrSqOt2Coeffdp;
  derivs(7) = dAtdp;

  return derivs;
}

// Parameters that we calculate the fine tuning for are stored in vectors
// in the order /{lambda, M1, M2, M3, M1', At, Alambda, mHd2, mHu2, ms2, mq3Sq, mtRSq\}
  Eigen::Matrix<double,3,1> doCalcRHSTuningVector_ESSM_Approx(flexiblesusy::genericE6SSM_soft_parameters modelAtMsusy, 
							  void (*ftBCatMX)(flexiblesusy::genericE6SSM_soft_parameters &, Eigen::ArrayXd &), 
							  Eigen::VectorXd const & vevs, tuning_parameters i, double mx, 
							  bool & hasError, Eigen::ArrayXd const & modelParsAtMx)
{
  Eigen::VectorXd rhsVec = vevs;

  if (vevs.size() != 3)
    {
      ostringstream kk;
      kk << "# WARNING: incorrect number of VEVs provided to doCalcRHSTuningVector_ESSM_Approx: exiting." << endl;
      throw kk.str();
    }
  else
    {
      double v1 = vevs(0);
      double v2 = vevs(1);
      double tb = v2/v1;
      double v = Sqrt(v1*v1+v2*v2);

      double s = vevs(2);
      
      genericE6SSM_input_parameters input = modelAtMsusy.get_input();
      
      double g1 = modelAtMsusy.get_g1();
      double g2 = modelAtMsusy.get_g2();
      double g3 = modelAtMsusy.get_g3();
      double gN = modelAtMsusy.get_gN();
      
      double Tlambda = modelAtMsusy.get_TLambdax();
      double lambda = modelAtMsusy.get_Lambdax();
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
      
      
      
      // First calculate the elements of the matrix that appears on the RHS in general, being derivatives of the EWSB
      // conditions wrt the EW scale parameters. For the EWSB conditions used in our study, this matrix has the form:
      // 
      //     [ df1/dlambda df1/dAlambda df1/dm_Hd^2 df1/dm_Hu^2 df1/dm_s^2 df1/dm_Ql^2 df1/dm_uR^2 df1/dA_t ]
      //     [ df2/dlambda df2/dAlambda df2/dm_Hd^2 df2/dm_Hu^2 df2/dm_s^2 df2/dm_Ql^2 df2/dm_uR^2 df2/dA_t ]
      //     [ df3/dlambda df3/dAlambda df3/dm_Hd^2 df3/dm_Hu^2 df3/dm_s^2 df3/dm_Ql^2 df3/dm_uR^2 df3/dA_t ]
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

	  EWderivs(0,0) += doCalcd2DeltaVdLambdadv1(modelAtMsusy, s, tb);
	  EWderivs(1,0) += doCalcd2DeltaVdLambdadv2(modelAtMsusy, s, tb);
	  EWderivs(2,0) += doCalcd2DeltaVdLambdadv3(modelAtMsusy, s, tb);
	  
	  EWderivs(0,5) += doCalcd2DeltaVdmQlsqdv1(modelAtMsusy, s, tb);
	  EWderivs(1,5) += doCalcd2DeltaVdmQlsqdv2(modelAtMsusy, s, tb);
	  EWderivs(2,5) += doCalcd2DeltaVdmQlsqdv3(modelAtMsusy, s, tb);
	  
	  EWderivs(0,6) += doCalcd2DeltaVdmUrsqdv1(modelAtMsusy, s, tb);
	  EWderivs(1,6) += doCalcd2DeltaVdmUrsqdv2(modelAtMsusy, s, tb);
	  EWderivs(2,6) += doCalcd2DeltaVdmUrsqdv3(modelAtMsusy, s, tb);
	  
	  EWderivs(0,7) += doCalcd2DeltaVdAtdv1(modelAtMsusy, s, tb);
	  EWderivs(1,7) += doCalcd2DeltaVdAtdv2(modelAtMsusy, s, tb);
	  EWderivs(2,7) += doCalcd2DeltaVdAtdv3(modelAtMsusy, s, tb);


	}    

      // Then calculate the vector that stores the derivatives of the EW scale parameters wrt the high scale parameter
      // of interest. In our study this has the form
      // [ dlambda/dp dAlambda/dp dm_Hd^2/dp dm_Hu^2/dp dm_s^2/dp dm_Ql^2/dp dm_uR^2/dp dA_t/dp]^T.

      // Get the scale M_SUSY - note that the model is assumed to be provided at this scale
      double ms = modelAtMsusy.get_scale();      

      Eigen::Matrix<double,8,1> paramDerivs;

      if (fabs(ms-mx) < TOLERANCE)
	{
	  switch(i)
	    {
	    case Au3: 
	      {
		paramDerivs(0) = 0.0;
		paramDerivs(1) = 0.0;
		paramDerivs(2) = 0.0;
		paramDerivs(3) = 0.0;
		paramDerivs(4) = 0.0;
		paramDerivs(5) = 0.0;
		paramDerivs(6) = 0.0;
		paramDerivs(7) = -1.0;
		break;
	      }
	    case mH13Sq:
	      {
		paramDerivs(0) = 0.0;
		paramDerivs(1) = 0.0;
		paramDerivs(2) = -1.0;
		paramDerivs(3) = 0.0;
		paramDerivs(4) = 0.0;
		paramDerivs(5) = 0.0;
		paramDerivs(6) = 0.0;
		paramDerivs(7) = 0.0;
		break;
	      }
	    case mH23Sq: 
	      {
		paramDerivs(0) = 0.0;
		paramDerivs(1) = 0.0;
		paramDerivs(2) = 0.0;
		paramDerivs(3) = -1.0;
		paramDerivs(4) = 0.0;
		paramDerivs(5) = 0.0;
		paramDerivs(6) = 0.0;
		paramDerivs(7) = 0.0;
		break;
	      }
	    case mS3Sq:
	      {
		paramDerivs(0) = 0.0;
		paramDerivs(1) = 0.0;
		paramDerivs(2) = 0.0;
		paramDerivs(3) = 0.0;
		paramDerivs(4) = -1.0;
		paramDerivs(5) = 0.0;
		paramDerivs(6) = 0.0;
		paramDerivs(7) = 0.0;
		break;
	      }
	    case lam3: 
	      {
		paramDerivs(0) = -1.0;
		paramDerivs(1) = 0.0;
		paramDerivs(2) = 0.0;
		paramDerivs(3) = 0.0;
		paramDerivs(4) = 0.0;
		paramDerivs(5) = 0.0;
		paramDerivs(6) = 0.0;
		paramDerivs(7) = 0.0;
		break;
	      }
	    case Alam3:
	      {
		paramDerivs(0) = 0.0;
		paramDerivs(1) = -1.0;
		paramDerivs(2) = 0.0;
		paramDerivs(3) = 0.0;
		paramDerivs(4) = 0.0;
		paramDerivs(5) = 0.0;
		paramDerivs(6) = 0.0;
		paramDerivs(7) = 0.0;
		break;
	      }
	    case mqL3Sq: 
	      {
		paramDerivs(0) = 0.0;
		paramDerivs(1) = 0.0;
		paramDerivs(2) = 0.0;
		paramDerivs(3) = 0.0;
		paramDerivs(4) = 0.0;
		paramDerivs(5) = -1.0;
		paramDerivs(6) = 0.0;
		paramDerivs(7) = 0.0;
		break;
	      }
	    case mtRSq:
	      {
		paramDerivs(0) = 0.0;
		paramDerivs(1) = 0.0;
		paramDerivs(2) = 0.0;
		paramDerivs(3) = 0.0;
		paramDerivs(4) = 0.0;
		paramDerivs(5) = 0.0;
		paramDerivs(6) = -1.0;
		paramDerivs(7) = 0.0;
		break;
	      }
	    default:
	      {
		paramDerivs(0) = 0.0;
		paramDerivs(1) = 0.0;
		paramDerivs(2) = 0.0;
		paramDerivs(3) = 0.0;
		paramDerivs(4) = 0.0;
		paramDerivs(5) = 0.0;
		paramDerivs(6) = 0.0;
		paramDerivs(7) = 0.0;
		break;
	      }
	    }
	}
      else
	{
	  // Here we use the approximate expression for each of the derivatives
	  // of the low scale parameters w.r.t the high scale parameter of interest.

	  genericE6SSM_soft_parameters r = genericE6SSM_soft_parameters(modelAtMsusy.get_input());
	  r.set(modelParsAtMx);

	  switch(i)
	    {
	    case lam3:
	      {
		paramDerivs = doCalcLambdaDerivs(r, ms, 3, hasError);
		break;
	      }
	    case M1:
	      {
		paramDerivs = doCalcGauginoDerivs(r, ms, 1, hasError);
		break;
	      }
	    case M2:
	      {
		paramDerivs = doCalcGauginoDerivs(r, ms, 2, hasError);
		break;
	      }
	    case M3:
	      {
		paramDerivs = doCalcGauginoDerivs(r, ms, 3, hasError);
		break;
	      }
	    case M1p:
	      {
		paramDerivs = doCalcGauginoDerivs(r, ms, 4, hasError);
		break;
	      }
	    case Au3:
	      {
		paramDerivs = doCalcSoftAuDerivs(r, ms, 3, 3, hasError);
		break;
	      }
	    case Alam3:
	      {
		paramDerivs = doCalcSoftAlambdaDerivs(r, ms, 3, hasError);
		break;
	      }
	    case mH13Sq:
	      {
		paramDerivs = doCalcMh1SquaredDerivs(r, ms, 3, hasError);
		break;
	      }
	    case mH23Sq:
	      {
		paramDerivs = doCalcMh2SquaredDerivs(r, ms, 3, hasError);
		break;
	      }
	    case mS3Sq:
	      {
		paramDerivs = doCalcMsSquaredDerivs(r, ms, 3, hasError);
		break;
	      }
	    case mqL3Sq:
	      {
		paramDerivs = doCalcMq2Derivs(r, ms, 3, 3, hasError);
		break;
	      }
	    case mtRSq:
	      {
		paramDerivs = doCalcMu2Derivs(r, ms, 3, 3, hasError);
		break;
	      }
	    default:
	      {
		cerr << "WARNING: unrecognised ESSM parameter " << i << " requested: ignoring it." << endl;
		hasError = true;
		break;
	      }
	    }
	  // Note additional minus sign needed
	  paramDerivs *= -1.0;
	}
      
   
      // Multiply the two together and return the result (note the additional minus sign!).
      rhsVec = EWderivs*paramDerivs; 
    }

  return rhsVec;
}

void pE6SSMftBCs(flexiblesusy::genericE6SSM_soft_parameters & model, Eigen::ArrayXd & tuningPars)
{

  // Set all first and second generation Yukawas and A terms to
  // zero, consistent with the approach taken in the pMSSM
  model.set_Yu(0,0,0.0);
  model.set_Yu(1,1,0.0);
  model.set_Yd(0,0,0.0);
  model.set_Yd(1,1,0.0);
  model.set_Ye(0,0,0.0);
  model.set_Ye(1,1,0.0);

  model.set_TYu(0,0,0.0);
  model.set_TYu(1,1,0.0);
  model.set_TYd(0,0,0.0);
  model.set_TYd(1,1,0.0);
  model.set_TYe(0,0,0.0);
  model.set_TYe(1,1,0.0);

  // Make sure off-diagonal elements vanish
  model.set_Yu(0,1,0.0); model.set_Yu(0,2,0.0);
  model.set_Yu(1,0,0.0); model.set_Yu(1,2,0.0);
  model.set_Yu(2,0,0.0); model.set_Yu(2,1,0.0);

  model.set_Yd(0,1,0.0); model.set_Yd(0,2,0.0);
  model.set_Yd(1,0,0.0); model.set_Yd(1,2,0.0);
  model.set_Yd(2,0,0.0); model.set_Yd(2,1,0.0);

  model.set_Ye(0,1,0.0); model.set_Ye(0,2,0.0);
  model.set_Ye(1,0,0.0); model.set_Ye(1,2,0.0);
  model.set_Ye(2,0,0.0); model.set_Ye(2,1,0.0);

  model.set_TYu(0,1,0.0); model.set_TYu(0,2,0.0);
  model.set_TYu(1,0,0.0); model.set_TYu(1,2,0.0);
  model.set_TYu(2,0,0.0); model.set_TYu(2,1,0.0);

  model.set_TYd(0,1,0.0); model.set_TYd(0,2,0.0);
  model.set_TYd(1,0,0.0); model.set_TYd(1,2,0.0);
  model.set_TYd(2,0,0.0); model.set_TYd(2,1,0.0);

  model.set_TYe(0,1,0.0); model.set_TYe(0,2,0.0);
  model.set_TYe(1,0,0.0); model.set_TYe(1,2,0.0);
  model.set_TYe(2,0,0.0); model.set_TYe(2,1,0.0);

  // Make sure off-diagonal soft masses also vanish
  model.set_mq2(0,1,0.0); model.set_mq2(0,2,0.0);
  model.set_mq2(1,0,0.0); model.set_mq2(1,2,0.0);
  model.set_mq2(2,0,0.0); model.set_mq2(2,1,0.0);

  model.set_mu2(0,1,0.0); model.set_mu2(0,2,0.0);
  model.set_mu2(1,0,0.0); model.set_mu2(1,2,0.0);
  model.set_mu2(2,0,0.0); model.set_mu2(2,1,0.0);

  model.set_md2(0,1,0.0); model.set_md2(0,2,0.0);
  model.set_md2(1,0,0.0); model.set_md2(1,2,0.0);
  model.set_md2(2,0,0.0); model.set_md2(2,1,0.0);

  model.set_ml2(0,1,0.0); model.set_ml2(0,2,0.0);
  model.set_ml2(1,0,0.0); model.set_ml2(1,2,0.0);
  model.set_ml2(2,0,0.0); model.set_ml2(2,1,0.0);

  model.set_me2(0,1,0.0); model.set_me2(0,2,0.0);
  model.set_me2(1,0,0.0); model.set_me2(1,2,0.0);
  model.set_me2(2,0,0.0); model.set_me2(2,1,0.0);

  // Off-diagonal exotic Yukawas should vanish
  model.set_Kappa(0,1,0.0); model.set_Kappa(0,2,0.0);
  model.set_Kappa(1,0,0.0); model.set_Kappa(1,2,0.0);
  model.set_Kappa(2,0,0.0); model.set_Kappa(2,1,0.0);

  model.set_TKappa(0,1,0.0); model.set_TKappa(0,2,0.0);
  model.set_TKappa(1,0,0.0); model.set_TKappa(1,2,0.0);
  model.set_TKappa(2,0,0.0); model.set_TKappa(2,1,0.0);

  model.set_Lambda12(0,1,0.0); model.set_Lambda12(1,0,0.0);

  model.set_mH1I2(0,1,0.0); model.set_mH1I2(1,0,0.0);
  model.set_mH2I2(0,1,0.0); model.set_mH2I2(1,0,0.0);
  model.set_msI2(0,1,0.0); model.set_msI2(1,0,0.0);

  // Set the remaining parameters
  model.set_Lambdax(tuningPars(tuning_parameters::lam3));
  model.set_MassB(tuningPars(tuning_parameters::M1));
  model.set_MassWB(tuningPars(tuning_parameters::M2));
  model.set_MassG(tuningPars(tuning_parameters::M3));
  model.set_MassBp(tuningPars(tuning_parameters::M1p));
  model.set_TYu(2,2,model.get_Yu(2,2)*tuningPars(tuning_parameters::Au3));
  model.set_TLambdax(tuningPars(tuning_parameters::Alam3)*tuningPars(tuning_parameters::lam3));
  model.set_mHd2(tuningPars(tuning_parameters::mH13Sq));
  model.set_mHu2(tuningPars(tuning_parameters::mH23Sq));
  model.set_ms2(tuningPars(tuning_parameters::mS3Sq));
  model.set_mq2(2,2,tuningPars(tuning_parameters::mqL3Sq));
  model.set_mu2(2,2,tuningPars(tuning_parameters::mtRSq));
}

  // Functions for numerically solving the EWSB conditions using GSL routines. 
  // They are basically identical to there equivalents in the FlexibleSUSY
  // generated spectrum generator, but using the simplified tadpole equations.
  // Note, however, that those routines solve for the soft masses, whereas
  // these solve for the VEVs.
int ewsb_conditions(const gsl_vector* x, void* params, gsl_vector* f)
{
  const int number_of_ewsb_equations = 3;

   if (contains_nan(x, number_of_ewsb_equations)) {
      for (std::size_t i = 0; i < number_of_ewsb_equations; ++i)
         gsl_vector_set(f, i, std::numeric_limits<double>::max());
      return GSL_EDOM;
   }

   genericE6SSM_soft_parameters* model
      = static_cast<genericE6SSM_soft_parameters*>(params);

   double tadpole[number_of_ewsb_equations];

   model->set_vd(gsl_vector_get(x, 0));
   model->set_vu(gsl_vector_get(x, 1));
   model->set_vs(gsl_vector_get(x, 2));

   tadpole[0] = ESSM_EWSBCondition1(*model);
   tadpole[1] = ESSM_EWSBCondition2(*model);
   tadpole[2] = ESSM_EWSBCondition3(*model);

   for (std::size_t i = 0; i < number_of_ewsb_equations; ++i)
      gsl_vector_set(f, i, tadpole[i]);

   return GSL_SUCCESS;
}

  int solve_for_vevs_iteratively(genericE6SSM_soft_parameters* model, int nints, double intprecis)
{
   const gsl_multiroot_fsolver_type* solvers[] = {
      gsl_multiroot_fsolver_hybrid, gsl_multiroot_fsolver_hybrids,
         gsl_multiroot_fsolver_broyden

   };

   double x_init[3];
   vevs_initial_guess(*model, x_init);

   int status;
   for (std::size_t i = 0; i < sizeof(solvers)/sizeof(*solvers); ++i) {
     status = solve_for_vevs_iteratively_with(solvers[i], x_init, model, nints, intprecis);
      if (status == GSL_SUCCESS) {
         break;
      }
   }

   if (status != GSL_SUCCESS) {
     cerr << "WARNING: problem solving EWSB conditions: tunings may be unreliable" << endl;
   } 
   return status;
}


  void vevs_initial_guess(genericE6SSM_soft_parameters model, double x_init[3])
{
  x_init[0] = model.get_vd();
  x_init[1] = model.get_vu();
  x_init[2] = model.get_vs();

}

  int solve_for_vevs_iteratively_with(const gsl_multiroot_fsolver_type* solver,
				      const double x_init[3],
				      genericE6SSM_soft_parameters * model, 
				      int number_of_ewsb_iterations, 
				      double ewsb_iteration_precision)
{

   Root_finder<3> root_finder(ewsb_conditions,
                              model,
                              number_of_ewsb_iterations,
                              ewsb_iteration_precision);
   root_finder.set_solver_type(solver);
   const int status = root_finder.find_root(x_init);

   return status;
}


// Variables used for getting information to the functions used in numerically calculating the
// fine tuning.
static genericE6SSM_soft_parameters *tempsoftTuning; // < a SoftParsMssm object given at the input scale MX
static double ftMsusy; // < the value of M_{SUSY} for the above object
static int ftParChoice; // < index labelling the current parameter we are varying
static Eigen::ArrayXd ftPars; // < vector containing the parameters varied as part of the fine tuning
static void (*currentftBC)(genericE6SSM_soft_parameters & , Eigen::ArrayXd &);


// A function for getting the predicted running value of M_Z^2, including 1-loop top and
// stop tadpole contributions. Used for numerically calculating the fine tuning in the
// E6SSM. Passed to SOFTSUSY's calcDerivative routine as part of the fine tuning calculation.
double predpE6SSMMzSqRun(double parVal)
{

  // Save copies of all initial values
  genericE6SSM_soft_parameters savedObject(*tempsoftTuning);
  Eigen::ArrayXd savedPars = ftPars;
  double saved_scale = tempsoftTuning->get_scale();

  // Update parameter values
  ftPars(ftParChoice) = parVal;
  currentftBC(*tempsoftTuning, ftPars);

  // Recalculate VEVs at M_{SUSY} using Newton's method. 
  // Solve iteratively using routines provided in the GSL
  // library instead of writing our own like in the MSSM,
  // as they are likely to be more reliable.

  int maxIters = 100;
  double tol = 1.0e-2;

  double ms = ftMsusy; // < guess for M_{SUSY}.
  Eigen::Matrix<double,3,1> vevs(3);

  int sing = 0;

  tempsoftTuning->run_to(ms, PRECISION);
  
  // If we don't vary M_{SUSY}, we just need to run down
  // and recalculate the VEVs.
  sing = solve_for_vevs_iteratively(tempsoftTuning, maxIters, tol);

  if (sing != GSL_SUCCESS)
    {
      cerr << "WARNING: error encountered in solving E6SSM EWSB conditions." << endl;
    }
    
  double v1 = tempsoftTuning->get_vd();
  double v2 = tempsoftTuning->get_vu();
  double s = tempsoftTuning->get_vs();

  double tb = v2/v1;

  // Calculate the new M_Z^2
  double t1Ov1 = 0.0;
  double t2Ov2 = 0.0;

  if (INCLUDE1LPTADPOLES)
    {
      t1Ov1 = -doCalcTadpoleESSMH1(*tempsoftTuning, s, tb);
      t2Ov2 = -doCalcTadpoleESSMH2(*tempsoftTuning, s, tb);
    }

  double mHdSq = tempsoftTuning->get_mHd2();
  double mHuSq = tempsoftTuning->get_mHu2();
  double mSSq = tempsoftTuning->get_ms2();
  double lambda = tempsoftTuning->get_Lambdax();
  
  double g1 = tempsoftTuning->get_g1();
  double g2 = tempsoftTuning->get_g2();
  double gbar = Sqrt(g2*g2+0.6*g1*g1);
  double gdash_1 = tempsoftTuning->get_gN();
  
  genericE6SSM_input_parameters input = tempsoftTuning->get_input();      
  
  // We neglect U(1) mixing for now.
  double Qtilde_1 = input.QH1p;
  double Qtilde_2 = input.QH2p;
  double Qtilde_s = input.QSp;
  
  double MzSq = -Sqr(lambda*s)+(2.0*(mHdSq+t1Ov1-tb*tb*(mHuSq+t2Ov2)))/(tb*tb-1.0)
    + Sqr(gdash_1)*(Qtilde_1*v1*v1+Qtilde_2*v2*v2+Qtilde_s*s*s)*(Qtilde_1-Qtilde_2*tb*tb)/(tb*tb-1.0);
  

  // Reset objects
  *tempsoftTuning = savedObject;
  tempsoftTuning->set_scale(saved_scale);
  ftPars = savedPars;

  return MzSq;

}

// Calculates the fine tuning in the E6SSM numerically.
// The object r is assumed to be provided at MX. The fine tuning is calculated
// as p / M_Z^2 \frac{\partial M_Z^2}{\partial p} where p is a set of parameters. M_Z^2
// is calculated using the tree level EWSB conditions + one-loop stop and top tadpole contributions
// at the scale M_SUSY. The function BCatMX is used to set the parameters at the high input scale MX.
 Eigen::VectorXd doCalcESSMTuningNumerically(genericE6SSM_soft_parameters r, double ms, double mx, 
					  Eigen::ArrayXd pars,
					  void (*BCatMX)(genericE6SSM_soft_parameters & , Eigen::ArrayXd &))
 {

   int nPars = pars.size();

   Eigen::VectorXd tunings;

   tunings.resize(nPars);

   // Check that we are at the right scale
   if (fabs(r.get_scale()-mx) > TOLERANCE)
    {
      cerr << "WARNING: object provided to doCalcESSMTuningNumerically should be at scale MX = " << mx
	   << ", not Q = " << r.get_scale() << ": skipping calculation." << endl;
      for (int i = 0; i < nPars; i++)
	{
	  tunings(i) = -numberOfTheBeast; // negative values indicate error, since tunings are positive by definition
	}
    }
  else
    {
      // Now for each parameter, estimate the numerical derivative by varying its
      // value at the input scale MX, running down to M_SUSY and recalculating M_Z.

      ftMsusy = ms;
      currentftBC = BCatMX;
      ftPars = pars;
      tempsoftTuning = &r;

      double epsilon = 1.0e-5;
      double temp;
      double h;
      double deriv;
      double err;

      // Initial value of M_Z^2
      genericE6SSM_soft_parameters w = r;
      w.run_to(ms, PRECISION);

      double v1 = w.get_vd();
      double v2 = w.get_vu();
      double s = w.get_vs();
      double tb = v2/v1;

      double t1Ov1 = 0.0;
      double t2Ov2 = 0.0;

      if (INCLUDE1LPTADPOLES)
	{
	  t1Ov1 = calculateDeltaTadpole1(w, s, tb);//-doCalcTadpoleESSMH1(w, s, tb);
	  t2Ov2 = calculateDeltaTadpole2(w,s,tb);//-doCalcTadpoleESSMH2(w, s, tb);
	}
      
      double mHdSq = w.get_mHd2();
      double mHuSq = w.get_mHu2();
      double mSSq = w.get_ms2();
      double lambda = w.get_Lambdax();
      
      double g1 = w.get_g1();
      double g2 = w.get_g2();
      double gbar = Sqrt(g2*g2+0.6*g1*g1);
      double gdash_1 = w.get_gN();

      genericE6SSM_input_parameters input = w.get_input();      

      // We neglect U(1) mixing for now.
      double Qtilde_1 = input.QH1p;
      double Qtilde_2 = input.QH2p;
      double Qtilde_s = input.QSp;

      double refMzSq = -Sqr(lambda*s)+(2.0*(mHdSq+t1Ov1-tb*tb*(mHuSq+t2Ov2)))/(tb*tb-1.0)
	+ Sqr(gdash_1)*(Qtilde_1*v1*v1+Qtilde_2*v2*v2+Qtilde_s*s*s)*(Qtilde_1-Qtilde_2*tb*tb)/(tb*tb-1.0);

      // Loop over the provided parameters and calculate the fine tuning for each
      for (int i = 0; i < nPars; i++)
	{
	  ftParChoice = i;

	  // Initial estimate for step size h.
	  h = epsilon*Abs(pars(i));
	  if (h == 0.0) h = epsilon;
	  temp = pars(i);
	  pars(i) = temp + h;
	  h = pars(i) - temp;

	  deriv = calcDerivative(predpE6SSMMzSqRun, temp, h, &err);

	  if (pars(i) > TOLERANCE && fabs(err / deriv) > 1.0) 
	    {
	      tunings(i) = -numberOfTheBeast; // tuning is inaccurate, so flag using negative value
	    }
	  else
	    {
	      tunings(i) = fabs((temp/refMzSq)*deriv); //DH::NOTE
	    }

	}

    }
  return tunings;
}

  Eigen::Matrix<double,tuning_parameters::NUMESSMTUNINGPARS,1> doCalcESSMTuningApprox(genericE6SSM_soft_parameters r, double ms, double mx, 
					  Eigen::ArrayXd pars, bool & hasTuningProblem, bool useApproxSolns)
  {
    Eigen::Matrix<double,tuning_parameters::NUMESSMTUNINGPARS,1> tunings;

    // Check that we are at the right scale
    if (fabs(r.get_scale()-mx) > TOLERANCE)
      {
	cerr << "WARNING: object provided to doCalcESSMTuningApprox should be at scale MX = " << mx
	     << ", not Q = " << r.get_scale() << ": skipping calculation." << endl;
	for (int i = 0; i < tuning_parameters::NUMESSMTUNINGPARS; i++)
	  {
	    tunings(i) = -numberOfTheBeast; // negative values indicate error, since tunings are positive by definition
	  }
    }
    else
    {
      // Save values of tuning parameters at MX
      double lambdaAtMx = r.get_Lambdax();
      double TlambdaAtMx = r.get_TLambdax();

      double AlambdaAtMx;
      
      if (Abs(TlambdaAtMx) < EPSTOL) 
	{
	  AlambdaAtMx = 0.0;
	}
      else if (Abs(lambdaAtMx) < 1.0e-100)
	{
	  ostringstream ii;
	  ii << "WARNING: trying to calculate A_lambda where lambda3 coupling is " <<
	    Abs(lambdaAtMx) << endl;
	  throw ii.str();
	}
      else
	{
	  AlambdaAtMx = TlambdaAtMx/lambdaAtMx;
	}

      double mHdSqAtMx = r.get_mHd2();
      double mHuSqAtMx = r.get_mHu2();
      double mSSqAtMx = r.get_ms2();
      double mqL3SqAtMx = r.get_mq2(2,2);
      double mtRSqAtMx = r.get_mu2(2,2);

      double ytAtMx = r.get_Yu(2,2);
      double TytAtMx = r.get_TYu(2,2);

      double AtAtMx;
      
      if (Abs(TytAtMx) < EPSTOL) 
	{
	  AtAtMx = 0.0;
	}
      else if (Abs(ytAtMx) < 1.0e-100)
	{
	  ostringstream ii;
	  ii << "WARNING: trying to calculate A_t where y_t coupling is " <<
	    Abs(ytAtMx) << endl;
	  throw ii.str();
	}
      else
	{
	  AtAtMx = TytAtMx/ytAtMx;
	}

      double M1AtMx = r.get_MassB();
      double M2AtMx = r.get_MassWB();
      double M3AtMx = r.get_MassG();
      double M1pAtMx = r.get_MassBp();

      // Calculate vectors of derivatives at MX, either using approximate analytic solutions 
      // or numerically (for later: change this so more easily extensible to general parameter set)
      Eigen::Matrix<double,8,1> rhs_vec_lambda, rhs_vec_Alambda, rhs_vec_mHdSq, rhs_vec_mHuSq, rhs_vec_mSSq;
      Eigen::Matrix<double,8,1> rhs_vec_mqL3Sq, rhs_vec_mtRSq, rhs_vec_At, rhs_vec_M1, rhs_vec_M2, rhs_vec_M3, rhs_vec_M1p;

      int nLogs = 1; //< include leading logarithms? Default is 1 (i.e. yes) for numerical calculation

      if (useApproxSolns)
	{

	  rhs_vec_lambda = doCalcLambdaDerivs(r, ms, 3, hasTuningProblem);
	  rhs_vec_Alambda = doCalcSoftAlambdaDerivs(r, ms, 3, hasTuningProblem);
	  rhs_vec_mHdSq = doCalcMh1SquaredDerivs(r, ms, 3, hasTuningProblem);
	  rhs_vec_mHuSq = doCalcMh2SquaredDerivs(r, ms, 3, hasTuningProblem);
	  rhs_vec_mSSq = doCalcMsSquaredDerivs(r, ms, 3, hasTuningProblem);
	  rhs_vec_mqL3Sq = doCalcMq2Derivs(r, ms, 3, 3, hasTuningProblem);
	  rhs_vec_mtRSq = doCalcMu2Derivs(r, ms, 3, 3, hasTuningProblem);
	  rhs_vec_At = doCalcSoftAuDerivs(r, ms, 3, 3, hasTuningProblem);
	  rhs_vec_M1 = doCalcGauginoDerivs(r, ms, 1, hasTuningProblem);
	  rhs_vec_M2 = doCalcGauginoDerivs(r, ms, 2, hasTuningProblem);
	  rhs_vec_M3 = doCalcGauginoDerivs(r, ms, 3, hasTuningProblem);
	  rhs_vec_M1p = doCalcGauginoDerivs(r, ms, 4, hasTuningProblem);
	}
      else
	{
	  rhs_vec_lambda = doCalcNumericDerivs(r, ms, tuning_parameters::lam3, r.get_loops(), nLogs);
	  rhs_vec_Alambda = doCalcNumericDerivs(r, ms, tuning_parameters::Alam3, r.get_loops(), nLogs);
	  // For the soft masses, the analytic and numerical routines are identical...
	  rhs_vec_mHdSq = doCalcMh1SquaredDerivs(r, ms, 3, hasTuningProblem);
	  rhs_vec_mHuSq = doCalcMh2SquaredDerivs(r, ms, 3, hasTuningProblem);
	  rhs_vec_mSSq = doCalcMsSquaredDerivs(r, ms, 3, hasTuningProblem);
	  rhs_vec_mqL3Sq = doCalcMq2Derivs(r, ms, 3, 3, hasTuningProblem);
	  rhs_vec_mtRSq = doCalcMu2Derivs(r, ms, 3, 3, hasTuningProblem);
	  rhs_vec_At = doCalcNumericDerivs(r, ms, tuning_parameters::Au3, r.get_loops(), nLogs);
	  rhs_vec_M1 = doCalcNumericDerivs(r, ms, tuning_parameters::M1, r.get_loops(), nLogs);
	  rhs_vec_M2 = doCalcNumericDerivs(r, ms, tuning_parameters::M2, r.get_loops(), nLogs);
	  rhs_vec_M3 = doCalcNumericDerivs(r, ms, tuning_parameters::M3, r.get_loops(), nLogs);
	  rhs_vec_M1p = doCalcNumericDerivs(r, ms, tuning_parameters::M1p, r.get_loops(), nLogs);
	}

      // Run to M_{SUSY} and calculate the matrices appearing on the left and
      // right hand sides (which are the same for all parameters)
      r.run_to(ms);

      double v1 = r.get_vd();
      double v2 = r.get_vu();
      double tb = v2/v1;
      double v = Sqrt(v1*v1+v2*v2);

      double s = r.get_vs();
      
      genericE6SSM_input_parameters input = r.get_input();
      
      double g1 = r.get_g1();
      double g2 = r.get_g2();
      double g3 = r.get_g3();
      double gN = r.get_gN();
      
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
      
      // First calculate the elements of the matrix that appears on the RHS in general, being derivatives of the EWSB
      // conditions wrt the EW scale parameters. For the EWSB conditions used in our study, this matrix has the form:
      // 
      //     [ df1/dlambda df1/dAlambda df1/dm_Hd^2 df1/dm_Hu^2 df1/dm_s^2 df1/dm_Ql^2 df1/dm_uR^2 df1/dA_t ]
      //     [ df2/dlambda df2/dAlambda df2/dm_Hd^2 df2/dm_Hu^2 df2/dm_s^2 df2/dm_Ql^2 df2/dm_uR^2 df2/dA_t ]
      //     [ df3/dlambda df3/dAlambda df3/dm_Hd^2 df3/dm_Hu^2 df3/dm_s^2 df3/dm_Ql^2 df3/dm_uR^2 df3/dA_t ]
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
	  EWderivs(0,0) += doCalcd2DeltaVdLambdadv1(r, s, tb);
	  EWderivs(1,0) += doCalcd2DeltaVdLambdadv2(r, s, tb);
	  EWderivs(2,0) += doCalcd2DeltaVdLambdadv3(r, s, tb);
	  
	  EWderivs(0,5) += doCalcd2DeltaVdmQlsqdv1(r, s, tb);
	  EWderivs(1,5) += doCalcd2DeltaVdmQlsqdv2(r, s, tb);
	  EWderivs(2,5) += doCalcd2DeltaVdmQlsqdv3(r, s, tb);
	  
	  EWderivs(0,6) += doCalcd2DeltaVdmUrsqdv1(r, s, tb);
	  EWderivs(1,6) += doCalcd2DeltaVdmUrsqdv2(r, s, tb);
	  EWderivs(2,6) += doCalcd2DeltaVdmUrsqdv3(r, s, tb);
	  
	  EWderivs(0,7) += doCalcd2DeltaVdAtdv1(r, s, tb);
	  EWderivs(1,7) += doCalcd2DeltaVdAtdv2(r, s, tb);
	  EWderivs(2,7) += doCalcd2DeltaVdAtdv3(r, s, tb);
	}    

      // Also save the VEVs and couplings at M_SUSY
      Eigen::Matrix<double,3,1> vevs;
      vevs(0) = v1; vevs(1) = v2; vevs(2) = s;

      Eigen::Matrix<double,3,3> lhs = doCalcLHSTuningMatrix(r, vevs);

      // Solve using QR routine for each parameter
      Eigen::Matrix<double,3,1> dvevs_dlambda, dvevs_dAlambda, dvevs_dmHdSq, dvevs_dmHuSq, dvevs_dmSSq;
      Eigen::Matrix<double,3,1> dvevs_dmqL3Sq, dvevs_dmtRSq, dvevs_dAt, dvevs_dM1, dvevs_dM2, dvevs_dM3, dvevs_dM1p;

      dvevs_dlambda = lhs.fullPivHouseholderQr().solve(-EWderivs*rhs_vec_lambda);
      dvevs_dAlambda = lhs.fullPivHouseholderQr().solve(-EWderivs*rhs_vec_Alambda);
      dvevs_dmHdSq = lhs.fullPivHouseholderQr().solve(-EWderivs*rhs_vec_mHdSq);
      dvevs_dmHuSq = lhs.fullPivHouseholderQr().solve(-EWderivs*rhs_vec_mHuSq);
      dvevs_dmSSq = lhs.fullPivHouseholderQr().solve(-EWderivs*rhs_vec_mSSq);
      dvevs_dmqL3Sq = lhs.fullPivHouseholderQr().solve(-EWderivs*rhs_vec_mqL3Sq);
      dvevs_dmtRSq = lhs.fullPivHouseholderQr().solve(-EWderivs*rhs_vec_mtRSq);
      dvevs_dAt = lhs.fullPivHouseholderQr().solve(-EWderivs*rhs_vec_At);
      dvevs_dM1 = lhs.fullPivHouseholderQr().solve(-EWderivs*rhs_vec_M1);
      dvevs_dM2 = lhs.fullPivHouseholderQr().solve(-EWderivs*rhs_vec_M2);
      dvevs_dM3 = lhs.fullPivHouseholderQr().solve(-EWderivs*rhs_vec_M3);
      dvevs_dM1p = lhs.fullPivHouseholderQr().solve(-EWderivs*rhs_vec_M1p);

      // Check solutions are valid
      bool has_lambda_solution = (lhs*dvevs_dlambda).isApprox(-EWderivs*rhs_vec_lambda);
      bool has_Alambda_solution = (lhs*dvevs_dAlambda).isApprox(-EWderivs*rhs_vec_Alambda);
      bool has_mHdSq_solution = (lhs*dvevs_dmHdSq).isApprox(-EWderivs*rhs_vec_mHdSq);
      bool has_mHuSq_solution = (lhs*dvevs_dmHuSq).isApprox(-EWderivs*rhs_vec_mHuSq);
      bool has_mSSq_solution = (lhs*dvevs_dmSSq).isApprox(-EWderivs*rhs_vec_mSSq);
      bool has_mqL3Sq_solution = (lhs*dvevs_dmqL3Sq).isApprox(-EWderivs*rhs_vec_mqL3Sq);
      bool has_mtRSq_solution = (lhs*dvevs_dmtRSq).isApprox(-EWderivs*rhs_vec_mtRSq);
      bool has_At_solution = (lhs*dvevs_dAt).isApprox(-EWderivs*rhs_vec_At);
      bool has_M1_solution = (lhs*dvevs_dM1).isApprox(-EWderivs*rhs_vec_M1);
      bool has_M2_solution = (lhs*dvevs_dM2).isApprox(-EWderivs*rhs_vec_M2);
      bool has_M3_solution = (lhs*dvevs_dM3).isApprox(-EWderivs*rhs_vec_M3);
      bool has_M1p_solution = (lhs*dvevs_dM1p).isApprox(-EWderivs*rhs_vec_M1p);


      if (!has_lambda_solution || !has_Alambda_solution || !has_mHdSq_solution
	  || !has_mHuSq_solution || !has_mSSq_solution || !has_mqL3Sq_solution
	  || !has_mtRSq_solution || !has_At_solution || !has_M1_solution
	  || !has_M2_solution || !has_M3_solution || !has_M1p_solution)
	{
	  cerr << "WARNING: problem calculating fine tunings: calculated fine tunings are inaccurate" << endl;
	  hasTuningProblem = true;
	}

      // Calculate fine tuning


      tunings(tuning_parameters::lam3) = Abs(doCalcdLogMzSqdLogParam(r, lambdaAtMx, vevs, dvevs_dlambda));
      tunings(tuning_parameters::Alam3) = Abs(doCalcdLogMzSqdLogParam(r, AlambdaAtMx, vevs, dvevs_dAlambda));
      tunings(tuning_parameters::mH13Sq) = Abs(doCalcdLogMzSqdLogParam(r, mHdSqAtMx, vevs, dvevs_dmHdSq));
      tunings(tuning_parameters::mH23Sq) = Abs(doCalcdLogMzSqdLogParam(r, mHuSqAtMx, vevs, dvevs_dmHuSq));
      tunings(tuning_parameters::mS3Sq) = Abs(doCalcdLogMzSqdLogParam(r, mSSqAtMx, vevs, dvevs_dmSSq));
      tunings(tuning_parameters::mqL3Sq) = Abs(doCalcdLogMzSqdLogParam(r, mqL3SqAtMx, vevs, dvevs_dmqL3Sq));
      tunings(tuning_parameters::mtRSq) = Abs(doCalcdLogMzSqdLogParam(r, mtRSqAtMx, vevs, dvevs_dmtRSq));
      tunings(tuning_parameters::Au3) = Abs(doCalcdLogMzSqdLogParam(r, AtAtMx, vevs, dvevs_dAt));
      tunings(tuning_parameters::M1) = Abs(doCalcdLogMzSqdLogParam(r, M1AtMx, vevs, dvevs_dM1));
      tunings(tuning_parameters::M2) = Abs(doCalcdLogMzSqdLogParam(r, M2AtMx, vevs, dvevs_dM2));
      tunings(tuning_parameters::M3) = Abs(doCalcdLogMzSqdLogParam(r, M3AtMx, vevs, dvevs_dM3));
      tunings(tuning_parameters::M1p) = Abs(doCalcdLogMzSqdLogParam(r, M1pAtMx, vevs, dvevs_dM1p));

    }

    return tunings;
  }

  // Useful functions for calculating derivatives of the 
  // low scale parameters numerically
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

  // Calculate the fine tuning using the matrix approach, but 
  // computing the derivatives of the high scale parameters
  // numerically (not using the approximate solutions at all).
  // If there are no issues with the calculation of the matrices, 
  // this should agree exactly with the results of doCalcESSMTuningNumerically.
  Eigen::Matrix<double,tuning_parameters::NUMESSMTUNINGPARS,1> doCalcESSMTuningSemianalytic(genericE6SSM_soft_parameters r, 
  											    double ms, double mx, 
  											    bool & hasTuningProblem)
  {
    // Object is provided at MX
    Eigen::Matrix<double,tuning_parameters::NUMESSMTUNINGPARS,1> tunings;

    // Check that we are at the right scale
    if (fabs(r.get_scale()-mx) > TOLERANCE)
      {
  	cerr << "WARNING: object provided to doCalcESSMTuningSemianalytic should be at scale MX = " << mx
  	     << ", not Q = " << r.get_scale() << ": skipping calculation." << endl;
  	for (int i = 0; i < tuning_parameters::NUMESSMTUNINGPARS; i++)
  	  {
  	    tunings(i) = -numberOfTheBeast; // negative values indicate error, since tunings are positive by definition
  	  }
  	return tunings;
      }

    // Save values of tuning parameters at MX
    double lambdaAtMx = r.get_Lambdax();
    double TlambdaAtMx = r.get_TLambdax();
    
    double AlambdaAtMx;
    
    if (Abs(TlambdaAtMx) < EPSTOL) 
      {
  	AlambdaAtMx = 0.0;
      }
    else if (Abs(lambdaAtMx) < 1.0e-100)
      {
  	ostringstream ii;
  	ii << "WARNING: trying to calculate A_lambda where lambda3 coupling is " <<
  	  Abs(lambdaAtMx) << endl;
  	throw ii.str();
      }
    else
      {
  	AlambdaAtMx = TlambdaAtMx/lambdaAtMx;
      }
    
    double mHdSqAtMx = r.get_mHd2();
    double mHuSqAtMx = r.get_mHu2();
    double mSSqAtMx = r.get_ms2();
    double mqL3SqAtMx = r.get_mq2(2,2);
    double mtRSqAtMx = r.get_mu2(2,2);
    
    double ytAtMx = r.get_Yu(2,2);
    double TytAtMx = r.get_TYu(2,2);
    
    double AtAtMx;
    
    if (Abs(TytAtMx) < EPSTOL) 
      {
  	AtAtMx = 0.0;
      }
    else if (Abs(ytAtMx) < 1.0e-100)
      {
  	ostringstream ii;
  	ii << "WARNING: trying to calculate A_t where y_t coupling is " <<
  	  Abs(ytAtMx) << endl;
  	throw ii.str();
      }
    else
      {
  	AtAtMx = TytAtMx/ytAtMx;
      }
    
    double M1AtMx = r.get_MassB();
    double M2AtMx = r.get_MassWB();
    double M3AtMx = r.get_MassG();
    double M1pAtMx = r.get_MassBp();
    
    // Calculate the derivatives of the low scale parameters w.r.t the
    // high scale input parameters
    Eigen::Matrix<double,8,1> rhs_vec_lambda, rhs_vec_Alambda, rhs_vec_mHdSq, rhs_vec_mHuSq, rhs_vec_mSSq;
    Eigen::Matrix<double,8,1> rhs_vec_mqL3Sq, rhs_vec_mtRSq, rhs_vec_At, rhs_vec_M1, rhs_vec_M2, rhs_vec_M3, rhs_vec_M1p;

    Eigen::Matrix<double,tuning_parameters::NUMESSMTUNINGPARS,1> params;

    params(tuning_parameters::lam3) = lambdaAtMx;
    params(tuning_parameters::Alam3) = AlambdaAtMx;
    params(tuning_parameters::mH13Sq) = mHdSqAtMx;
    params(tuning_parameters::mH23Sq) = mHuSqAtMx;
    params(tuning_parameters::mS3Sq) = mSSqAtMx;
    params(tuning_parameters::mqL3Sq) = mqL3SqAtMx;
    params(tuning_parameters::mtRSq) = mtRSqAtMx;
    params(tuning_parameters::Au3) = AtAtMx;
    params(tuning_parameters::M1) = M1AtMx;
    params(tuning_parameters::M2) = M2AtMx;
    params(tuning_parameters::M3) = M3AtMx;
    params(tuning_parameters::M1p) = M1pAtMx;

    rhs_vec_lambda = doCalcdLowScaledHighScale(r, ms, params, tuning_parameters::lam3, pE6SSMftBCs, hasTuningProblem);
    rhs_vec_Alambda = doCalcdLowScaledHighScale(r, ms, params, tuning_parameters::Alam3, pE6SSMftBCs, hasTuningProblem);
    rhs_vec_mHdSq = doCalcdLowScaledHighScale(r, ms, params, tuning_parameters::mH13Sq, pE6SSMftBCs, hasTuningProblem);
    rhs_vec_mHuSq = doCalcdLowScaledHighScale(r, ms, params, tuning_parameters::mH23Sq, pE6SSMftBCs, hasTuningProblem);
    rhs_vec_mSSq = doCalcdLowScaledHighScale(r, ms, params, tuning_parameters::mS3Sq, pE6SSMftBCs, hasTuningProblem);
    rhs_vec_mqL3Sq = doCalcdLowScaledHighScale(r, ms, params, tuning_parameters::mqL3Sq, pE6SSMftBCs, hasTuningProblem);
    rhs_vec_mtRSq = doCalcdLowScaledHighScale(r, ms, params, tuning_parameters::mtRSq, pE6SSMftBCs, hasTuningProblem);
    rhs_vec_At = doCalcdLowScaledHighScale(r, ms, params, tuning_parameters::Au3, pE6SSMftBCs, hasTuningProblem);
    rhs_vec_M1 = doCalcdLowScaledHighScale(r, ms, params, tuning_parameters::M1, pE6SSMftBCs, hasTuningProblem);
    rhs_vec_M2 = doCalcdLowScaledHighScale(r, ms, params, tuning_parameters::M2, pE6SSMftBCs, hasTuningProblem);
    rhs_vec_M3 = doCalcdLowScaledHighScale(r, ms, params, tuning_parameters::M3, pE6SSMftBCs, hasTuningProblem);
    rhs_vec_M1p = doCalcdLowScaledHighScale(r, ms, params, tuning_parameters::M1p, pE6SSMftBCs, hasTuningProblem);

    // Run to M_{SUSY} and calculate the matrices appearing on the left and
    // right hand sides (which are the same for all parameters)
    r.run_to(ms);
    
    double v1 = r.get_vd();
    double v2 = r.get_vu();
    double tb = v2/v1;
    double v = Sqrt(v1*v1+v2*v2);
    
    double s = r.get_vs();
    
    genericE6SSM_input_parameters input = r.get_input();
    
    double g1 = r.get_g1();
    double g2 = r.get_g2();
    double g3 = r.get_g3();
    double gN = r.get_gN();
    
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
    
    // First calculate the elements of the matrix that appears on the RHS in general, being derivatives of the EWSB
    // conditions wrt the EW scale parameters. For the EWSB conditions used in our study, this matrix has the form:
    // 
    //     [ df1/dlambda df1/dAlambda df1/dm_Hd^2 df1/dm_Hu^2 df1/dm_s^2 df1/dm_Ql^2 df1/dm_uR^2 df1/dA_t ]
    //     [ df2/dlambda df2/dAlambda df2/dm_Hd^2 df2/dm_Hu^2 df2/dm_s^2 df2/dm_Ql^2 df2/dm_uR^2 df2/dA_t ]
    //     [ df3/dlambda df3/dAlambda df3/dm_Hd^2 df3/dm_Hu^2 df3/dm_s^2 df3/dm_Ql^2 df3/dm_uR^2 df3/dA_t ]
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
  	EWderivs(0,0) += doCalcd2DeltaVdLambdadv1(r, s, tb);
  	EWderivs(1,0) += doCalcd2DeltaVdLambdadv2(r, s, tb);
  	EWderivs(2,0) += doCalcd2DeltaVdLambdadv3(r, s, tb);
	
  	EWderivs(0,5) += doCalcd2DeltaVdmQlsqdv1(r, s, tb);
  	EWderivs(1,5) += doCalcd2DeltaVdmQlsqdv2(r, s, tb);
  	EWderivs(2,5) += doCalcd2DeltaVdmQlsqdv3(r, s, tb);
	
  	EWderivs(0,6) += doCalcd2DeltaVdmUrsqdv1(r, s, tb);
  	EWderivs(1,6) += doCalcd2DeltaVdmUrsqdv2(r, s, tb);
  	EWderivs(2,6) += doCalcd2DeltaVdmUrsqdv3(r, s, tb);
	
  	EWderivs(0,7) += doCalcd2DeltaVdAtdv1(r, s, tb);
  	EWderivs(1,7) += doCalcd2DeltaVdAtdv2(r, s, tb);
  	EWderivs(2,7) += doCalcd2DeltaVdAtdv3(r, s, tb);
      }    
    
    // Also save the VEVs and couplings at M_SUSY
    Eigen::Matrix<double,3,1> vevs;
    vevs(0) = v1; vevs(1) = v2; vevs(2) = s;
    
    Eigen::Matrix<double,3,3> lhs = doCalcLHSTuningMatrix(r, vevs);
    
    // Solve using QR routine for each parameter
    Eigen::Matrix<double,3,1> dvevs_dlambda, dvevs_dAlambda, dvevs_dmHdSq, dvevs_dmHuSq, dvevs_dmSSq;
    Eigen::Matrix<double,3,1> dvevs_dmqL3Sq, dvevs_dmtRSq, dvevs_dAt, dvevs_dM1, dvevs_dM2, dvevs_dM3, dvevs_dM1p;
    
    dvevs_dlambda = lhs.fullPivHouseholderQr().solve(-EWderivs*rhs_vec_lambda);
    dvevs_dAlambda = lhs.fullPivHouseholderQr().solve(-EWderivs*rhs_vec_Alambda);
    dvevs_dmHdSq = lhs.fullPivHouseholderQr().solve(-EWderivs*rhs_vec_mHdSq);
    dvevs_dmHuSq = lhs.fullPivHouseholderQr().solve(-EWderivs*rhs_vec_mHuSq);
    dvevs_dmSSq = lhs.fullPivHouseholderQr().solve(-EWderivs*rhs_vec_mSSq);
    dvevs_dmqL3Sq = lhs.fullPivHouseholderQr().solve(-EWderivs*rhs_vec_mqL3Sq);
    dvevs_dmtRSq = lhs.fullPivHouseholderQr().solve(-EWderivs*rhs_vec_mtRSq);
    dvevs_dAt = lhs.fullPivHouseholderQr().solve(-EWderivs*rhs_vec_At);
    dvevs_dM1 = lhs.fullPivHouseholderQr().solve(-EWderivs*rhs_vec_M1);
    dvevs_dM2 = lhs.fullPivHouseholderQr().solve(-EWderivs*rhs_vec_M2);
    dvevs_dM3 = lhs.fullPivHouseholderQr().solve(-EWderivs*rhs_vec_M3);
    dvevs_dM1p = lhs.fullPivHouseholderQr().solve(-EWderivs*rhs_vec_M1p);
    
    // Check solutions are valid
    bool has_lambda_solution = (lhs*dvevs_dlambda).isApprox(-EWderivs*rhs_vec_lambda);
    bool has_Alambda_solution = (lhs*dvevs_dAlambda).isApprox(-EWderivs*rhs_vec_Alambda);
    bool has_mHdSq_solution = (lhs*dvevs_dmHdSq).isApprox(-EWderivs*rhs_vec_mHdSq);
    bool has_mHuSq_solution = (lhs*dvevs_dmHuSq).isApprox(-EWderivs*rhs_vec_mHuSq);
    bool has_mSSq_solution = (lhs*dvevs_dmSSq).isApprox(-EWderivs*rhs_vec_mSSq);
    bool has_mqL3Sq_solution = (lhs*dvevs_dmqL3Sq).isApprox(-EWderivs*rhs_vec_mqL3Sq);
    bool has_mtRSq_solution = (lhs*dvevs_dmtRSq).isApprox(-EWderivs*rhs_vec_mtRSq);
    bool has_At_solution = (lhs*dvevs_dAt).isApprox(-EWderivs*rhs_vec_At);
    bool has_M1_solution = (lhs*dvevs_dM1).isApprox(-EWderivs*rhs_vec_M1);
    bool has_M2_solution = (lhs*dvevs_dM2).isApprox(-EWderivs*rhs_vec_M2);
    bool has_M3_solution = (lhs*dvevs_dM3).isApprox(-EWderivs*rhs_vec_M3);
    bool has_M1p_solution = (lhs*dvevs_dM1p).isApprox(-EWderivs*rhs_vec_M1p);
    
    
    if (!has_lambda_solution || !has_Alambda_solution || !has_mHdSq_solution
  	|| !has_mHuSq_solution || !has_mSSq_solution || !has_mqL3Sq_solution
  	|| !has_mtRSq_solution || !has_At_solution || !has_M1_solution
  	|| !has_M2_solution || !has_M3_solution || !has_M1p_solution)
      {
  	cerr << "WARNING: problem calculating fine tunings: calculated fine tunings are inaccurate" << endl;
  	hasTuningProblem = true;
      }
    
    // Calculate fine tuning
 

    tunings(tuning_parameters::lam3) = Abs(doCalcdLogMzSqdLogParam(r, lambdaAtMx, vevs, dvevs_dlambda));
    tunings(tuning_parameters::Alam3) = Abs(doCalcdLogMzSqdLogParam(r, AlambdaAtMx, vevs, dvevs_dAlambda));
    tunings(tuning_parameters::mH13Sq) = Abs(doCalcdLogMzSqdLogParam(r, mHdSqAtMx, vevs, dvevs_dmHdSq));
    tunings(tuning_parameters::mH23Sq) = Abs(doCalcdLogMzSqdLogParam(r, mHuSqAtMx, vevs, dvevs_dmHuSq));
    tunings(tuning_parameters::mS3Sq) = Abs(doCalcdLogMzSqdLogParam(r, mSSqAtMx, vevs, dvevs_dmSSq));
    tunings(tuning_parameters::mqL3Sq) = Abs(doCalcdLogMzSqdLogParam(r, mqL3SqAtMx, vevs, dvevs_dmqL3Sq));
    tunings(tuning_parameters::mtRSq) = Abs(doCalcdLogMzSqdLogParam(r, mtRSqAtMx, vevs, dvevs_dmtRSq));
    tunings(tuning_parameters::Au3) = Abs(doCalcdLogMzSqdLogParam(r, AtAtMx, vevs, dvevs_dAt));
    tunings(tuning_parameters::M1) = Abs(doCalcdLogMzSqdLogParam(r, M1AtMx, vevs, dvevs_dM1));
    tunings(tuning_parameters::M2) = Abs(doCalcdLogMzSqdLogParam(r, M2AtMx, vevs, dvevs_dM2));
    tunings(tuning_parameters::M3) = Abs(doCalcdLogMzSqdLogParam(r, M3AtMx, vevs, dvevs_dM3));
    tunings(tuning_parameters::M1p) = Abs(doCalcdLogMzSqdLogParam(r, M1pAtMx, vevs, dvevs_dM1p));
    
    return tunings;
    
  }


static genericE6SSM_soft_parameters *tempessm; // < an E6SSM object given at the input scale MX
static double modelMsusy; // < the value of M_{SUSY} for the above object
static int essmParChoice; // < index labelling the current parameter we are varying
static Eigen::ArrayXd tuningPars; // < vector containing the parameters varied as part of the fine tuning
static void (*currentessmBC)(genericE6SSM_soft_parameters & , Eigen::ArrayXd &);

  double doCalcLowScaleLambda(double par_val)
  {
    // Get copy of current model
    genericE6SSM_soft_parameters model(*tempessm);

    // Update parameter values
    Eigen::ArrayXd params = tuningPars;

    params(essmParChoice) = par_val;

    currentessmBC(model, params);

    model.run_to(modelMsusy, PRECISION);

    return model.get_Lambdax();

  }

  double doCalcLowScaleAlambda(double par_val)
  {
    // Get copy of current model
    genericE6SSM_soft_parameters model(*tempessm);

    // Update parameter values
    Eigen::ArrayXd params = tuningPars;

    params(essmParChoice) = par_val;

    currentessmBC(model, params);

    model.run_to(modelMsusy, PRECISION);

    double lambda = model.get_Lambdax();
    double Tlambda = model.get_TLambdax();
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

    return Alambda;

  }

  double doCalcLowScaleMh1Squared(double par_val)
  {
    // Get copy of current model
    genericE6SSM_soft_parameters model(*tempessm);

    // Update parameter values
    Eigen::ArrayXd params = tuningPars;

    params(essmParChoice) = par_val;

    currentessmBC(model, params);

    model.run_to(modelMsusy, PRECISION);

    return model.get_mHd2();

  }

  double doCalcLowScaleMh2Squared(double par_val)
  {
    // Get copy of current model
    genericE6SSM_soft_parameters model(*tempessm);

    // Update parameter values
    Eigen::ArrayXd params = tuningPars;

    params(essmParChoice) = par_val;

    currentessmBC(model, params);

    model.run_to(modelMsusy, PRECISION);

    return model.get_mHu2();

  }

  double doCalcLowScaleMsSquared(double par_val)
  {
    // Get copy of current model
    genericE6SSM_soft_parameters model(*tempessm);

    // Update parameter values
    Eigen::ArrayXd params = tuningPars;

    params(essmParChoice) = par_val;

    currentessmBC(model, params);

    model.run_to(modelMsusy, PRECISION);

    return model.get_ms2();

  }

  double doCalcLowScaleMqL3Squared(double par_val)
  {
    // Get copy of current model
    genericE6SSM_soft_parameters model(*tempessm);

    // Update parameter values
    Eigen::ArrayXd params = tuningPars;

    params(essmParChoice) = par_val;

    currentessmBC(model, params);

    model.run_to(modelMsusy, PRECISION);

    return model.get_mq2(2,2);

  }

  double doCalcLowScaleMtRSquared(double par_val)
  {
    // Get copy of current model
    genericE6SSM_soft_parameters model(*tempessm);

    // Update parameter values
    Eigen::ArrayXd params = tuningPars;

    params(essmParChoice) = par_val;

    currentessmBC(model, params);

    model.run_to(modelMsusy, PRECISION);

    return model.get_mu2(2,2);

  }

  double doCalcLowScaleAt(double par_val)
  {
    // Get copy of current model
    genericE6SSM_soft_parameters model(*tempessm);

    // Update parameter values
    Eigen::ArrayXd params = tuningPars;

    params(essmParChoice) = par_val;

    currentessmBC(model, params);

    model.run_to(modelMsusy, PRECISION);

    double yt = model.get_Yu(2,2);
    double Tyt = model.get_TYu(2,2);
    
    double At;
    
    if (Abs(Tyt) < EPSTOL) 
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
  	At = Tyt/yt;
      }

    return At;

  }

  Eigen::Matrix<double,8,1> doCalcdLowScaledHighScale(genericE6SSM_soft_parameters r, double ms, 
						      Eigen::MatrixXd params, unsigned i,
						      void (*BCatMX)(genericE6SSM_soft_parameters & , Eigen::ArrayXd &), 
						      bool & hasTuningProblem)
  {

    essmParChoice = i;
    tuningPars = params;
    modelMsusy = ms;
    currentessmBC = BCatMX;
    tempessm = &r;

    double deriv, h, temp, err;
    
    double epsilon = 1.0e-5;
    
    // Initial estimate for step size h.
    h = epsilon*Abs(params(i));
    if (h == 0.0) h = epsilon;
    temp = params(i);
    params(i) = temp + h;
    h = params(i) - temp;
    
    Eigen::Matrix<double,8,1> derivs, errvec;
    
    derivs(0) = calcDerivative(doCalcLowScaleLambda, temp, h, &err); errvec(0) = err;
    derivs(1) = calcDerivative(doCalcLowScaleAlambda, temp, h, &err); errvec(1) = err;
    derivs(2) = calcDerivative(doCalcLowScaleMh1Squared, temp, h, &err); errvec(2) = err;
    derivs(3) = calcDerivative(doCalcLowScaleMh2Squared, temp, h, &err); errvec(3) = err;
    derivs(4) = calcDerivative(doCalcLowScaleMsSquared, temp, h, &err); errvec(4) = err;
    derivs(5) = calcDerivative(doCalcLowScaleMqL3Squared, temp, h, &err); errvec(5) = err;
    derivs(6) = calcDerivative(doCalcLowScaleMtRSquared, temp, h, &err); errvec(6) = err;
    derivs(7) = calcDerivative(doCalcLowScaleAt, temp, h, &err); errvec(7) = err;

    for (int j = 0; j < 8; j++)
      {
	if (Abs(errvec(j)) > 1.0e-8 && Abs(errvec(j)/derivs(j)) > 1.0)
	  {
	    cerr << "WARNING: problem calculating derivatives of low scale parameters" << endl;
	    derivs(j) = numberOfTheBeast;
	    hasTuningProblem = true;
	  }
      }
    
    return derivs;
  }

} // namespace essm_tuning_utils
