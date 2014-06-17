/*
  Implementations of the functions defined in essmtuningutils.h
 */

#include "essmtuningutils.h"

using namespace flexiblesusy;

// A function to calculate the 3x3 matrix appearing on the LHS of the calculation
// of the derivatives of the VEVs wrt the parameters in the E6SSM. Note that this
// assumes that the E6SSM object is already set to have its renormalisation
// scale set at M_SUSY.
DoubleMatrix doCalcLHSTuningMatrix(genericE6SSM_soft_parameters essmSusy, DoubleVector const & vevs)
{
  
  DoubleMatrix lhsMat(3,3);
  
  if (vevs.displayEnd()-vevs.displayStart()+1 < 3)
    {
      cerr << "WARNING: incorrect number of VEVs supplied to function: skipping calculating LHS matrix." << endl;
    }
  else
    {
      
      double v1 = vevs(vevs.displayStart());
      double v2 = vevs(vevs.displayStart()+1);
      double s = vevs(vevs.displayEnd());
      double tb = v2/v1;
      double v = Sqrt(v1*v1+v2*v2);

      double g1 = essmSusy.get_g1();
      double g2 = essmSusy.get_g2();
      double gbar = Sqrt(g2*g2+0.6*g1*g1);
      double gdash_1 = essmSusy.get_gN();

      double cb = 1.0/sqrt(1.0+tb*tb);
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
      
            
      lhsMat(1,1) = df1dv1;
      lhsMat(1,2) = df1dv2;
      lhsMat(1,3) = df1ds;
      lhsMat(2,1) = df2dv1;
      lhsMat(2,2) = df2dv2;
      lhsMat(2,3) = df2ds;
      lhsMat(3,1) = df3dv1;
      lhsMat(3,2) = df3dv2;
      lhsMat(3,3) = df3ds;

    }
  
  return lhsMat;

}

// A function to calculate d log(M_Z^2)/d log(p) using the value of the parameter p
// and the already calculated derivatives d v_i/ d p_j. For the 
double doCalcdLogMzSqdLogParam(genericE6SSM_soft_parameters r, double p, DoubleVector const & vevs, DoubleVector const & dVevsdp)
{

  double deriv = 0.0;

  if ((vevs.displayEnd()-vevs.displayStart()+1 < 3) || (dVevsdp.displayEnd()-dVevsdp.displayStart()+1 < 3))
    {
      cerr << "WARNING: incorrect number of VEVs supplied to function: skipping calculating derivative." << endl;
    }
  else
    {
      double g1 = r.get_g1();
      double g2 = r.get_g2();
      double gbar = sqrt(g2*g2+0.6*g1*g1);
      
      double v1 = vevs(vevs.displayStart());
      double v2 = vevs(vevs.displayStart()+1);

      double v = sqrt(v1*v1+v2*v2);

      // Neglect any neutral mixing for now, but later may want to include contribution from it.      
      deriv = (2.0*p/(v*v))*(v1*dVevsdp(dVevsdp.displayStart())+v2*dVevsdp(dVevsdp.displayEnd()));
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
      f1 = m1Sq + 0.5*lambda*lambda*(v2Sq+s*s)-lambda*Alambda*s*tb/sqrt(2.0)
	+gbar*gbar*v*v*c2b/8.0 + (0.5*Qtilde_1*gdash_1*gdash_1*
				  (Qtilde_1*v*v*cSqb + Qtilde_2*v*v*sSqb + Qtilde_s*s*s));
    }
  else
    {
      f1 = m1Sq + 0.5*lambda*lambda*(v2Sq+s*s)-lambda*Alambda*s*tb/sqrt(2.0)
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
      f2 = m2Sq + 0.5*lambda*lambda*(v1Sq+s*s)-lambda*Alambda*s/(sqrt(2.0)*tb)
	-gbar*gbar*v*v*c2b/8.0 + (0.5*Qtilde_2*gdash_1*gdash_1*
				  (Qtilde_1*v*v*cSqb + Qtilde_2*v*v*sSqb + Qtilde_s*s*s));
    }
  else
    {
      f2 = m2Sq + 0.5*lambda*lambda*(v1Sq+s*s)-lambda*Alambda*s/(sqrt(2.0)*tb)
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
      f3 = msSq + 0.5*lambda*lambda*v*v-lambda*Alambda*v*v*s2b/(2.0*sqrt(2.0)*s)
	+0.5*Qtilde_s*gdash_1*gdash_1*(Qtilde_1*v*v*cSqb + Qtilde_2*v*v*sSqb + Qtilde_s*s*s);
    }
  else
    {
      f3 = msSq + 0.5*lambda*lambda*v*v-lambda*Alambda*v*v*s2b/(2.0*sqrt(2.0)*s)
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
static double MsusyScaleFactor = 1.0e3;
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
	      r.run_to(ms);
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

	      mSSq = (1.0/s)*(lambda*Alambda*v*v*s2b/(2.0*Sqrt(2.0))-0.5*lambda*lambda*v*v*s
			      -0.5*gdash_1*gdash_1*Qtilde_s*s*(Qtilde_1*v*v*cSqb+Qtilde_2*v*v*sSqb+Qtilde_s*s*s));
	      mH1Sq = lambda*Alambda*s*tb/Sqrt(2.0)-0.5*lambda*lambda*(s*s+v*v*sSqb)-gbar*gbar*v*v*c2b/8.0
		-0.5*Qtilde_1*gdash_1*gdash_1*(Qtilde_1*v*v*cSqb+Qtilde_2*v*v*sSqb+Qtilde_s*s*s);
	      mH2Sq = lambda*Alambda*s/(Sqrt(2.0)*tb)-0.5*lambda*lambda*(s*s+v*v*cSqb)+gbar*gbar*v*v*c2b/8.0
		-0.5*Qtilde_2*gdash_1*gdash_1*(Qtilde_1*v*v*cSqb+Qtilde_2*v*v*sSqb+Qtilde_s*s*s);

	    }
	  else
	    {

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

	  double t_run = log(mx/ms);

	  double mHuSq_guess, mHdSq_guess, mSSq_guess;

	  // If t is too large, this is an unreliable estimate, so just use a default value.
	  const int TLIMIT = 10;
	  if (t_run >= TLIMIT)
	    {
	      mHdSqScaleFactor = r.get_mHd2();
	      mHuSqScaleFactor = r.get_mHu2();
	      mSSqScaleFactor = r.get_ms2();

	      mHuSq_guess = mHuSqScaleFactor;
	      mHdSq_guess = mHdSqScaleFactor;

	    }
	  else
	    {
	      // The Taylor series expansion is about the solution at M_{SUSY} initially.
	      genericE6SSM_soft_parameters w = r; // copy object to avoid numerical errors
	      w.run_to(ms);
	      
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
	      
	      // Construct required coefficients in Taylor series.
	      double mHdSqLogCoeff = doCalcMh1SquaredLogCoeff(w, nLps);
	      double mHdSqLogSqCoeff = doCalcMh1SquaredLogSqCoeff(w, 1);
	      double mHuSqLogCoeff = doCalcMh2SquaredLogCoeff(w, nLps);
	      double mHuSqLogSqCoeff = doCalcMh2SquaredLogSqCoeff(w, 1);
	      double mSSqLogCoeff = doCalcMsSquaredLogCoeff(w, nLps);
	      double mSSqLogSqCoeff = doCalcMsSquaredLogSqCoeff(w, 1);
	      
	      mHuSq_guess = mHuSqInit + t_run*mHuSqLogCoeff + Sqr(t_run)*mHuSqLogSqCoeff;
	      mHdSq_guess = mHdSqInit + t_run*mHdSqLogCoeff + Sqr(t_run)*mHdSqLogSqCoeff;
	      mSSq_guess = mSSqInit + t_run*mSSqLogCoeff + Sqr(t_run)*mSSqLogSqCoeff;
	      
	      mHuSqScaleFactor = mHuSq_guess;
	      mHdSqScaleFactor = mHdSq_guess;
	      mSSqScaleFactor = mSSq_guess;

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

  double f = (ms/sqrt(mstop(1)*mstop(2)))-1.0;

  return f;
}

static genericE6SSM_soft_parameters *tempsoftessm;
void ESSM_EWSB_Shooter_Functions(DoubleVector const & parVals, DoubleVector & f)
{
  genericE6SSM_soft_parameters s = *tempsoftessm;

  s.set_mHd2(parVals(1)*mHdSqScaleFactor);
  s.set_mHu2(parVals(2)*mHuSqScaleFactor);
  s.set_ms2(parVals(3)*mSSqScaleFactor);

  s.run_to(parVals(4)*MsusyScaleFactor);

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

  w.run_to(estimate(4)*MsusyScaleFactor);

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
  double f4 = estimate(4)*MsusyScaleFactor/sqrt(mstop(1)*mstop(2))-1.0;

  // cout << "f1/mHdSq Scale Factor = " << f1/mHdSqScaleFactor << endl;
  // cout << "f2/mHdSq Scale Factor = " << f1/mHuSqScaleFactor << endl;
  // cout << "f3 = " << f3 << endl;

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
    cout << " yt*v2/(sqrt(2.0)) = " <<  yt*v2/(sqrt(2.0)) << endl;
    cout << "mtop = " << mtop << endl;
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

  mstopsq(1) = 0.5*(mQlsq +  mUrsq + 2*mtop*mtop +0.125*sqr(g2)*(sqr(v1) - sqr(v2)) +0.125*3.0*sqr(g1)*(sqr(v1) - sqr(v2))/5.0 + oneO40*sqr(g1p)*(-3.0*sqr(v1) - 2.0*sqr(v2) + 5.0*sqr(s))  + sqrt(sqr(mQlsq- mUrsq+0.125*sqr(g2)*(sqr(v1) - sqr(v2)) - 0.125*sqr(g1)*(sqr(v1) - sqr(v2))) + 4.0*mtop*mtop*sqr(At - lambda(3)*s*oneOrt2*v1/v2)));

  if (mstopsq(1) >= 0.0)
    {
  mstop(1) = sqrt(0.5*(mQlsq + mUrsq + 2*mtop*mtop + 0.125*sqr(g2)*(sqr(v1) - sqr(v2)) + 0.125*3.0*sqr(g1)*(sqr(v1) - sqr(v2))/5.0 + oneO40*sqr(g1p)*(-3.0*sqr(v1) - 2.0*sqr(v2) + 5.0*sqr(s))  + sqrt(sqr(mQlsq- mUrsq+0.125*sqr(g2)*(sqr(v1) - sqr(v2)) - 0.125*sqr(g1)*(sqr(v1) - sqr(v2))) + 4.0*mtop*mtop*sqr(At - lambda(3)*s*oneOrt2*v1/v2))));
    }
  else
    {
      mstop(1) = -sqrt(Abs(mstopsq(1)));
    }

  mstopsq(2) = 0.5*(mQlsq + mUrsq + 2*mtop*mtop +0.125*sqr(g2)*(sqr(v1) - sqr(v2)) +0.125*3.0*sqr(g1)*(sqr(v1) - sqr(v2))/5.0 + oneO40*sqr(g1p)*(-3.0*sqr(v1) - 2.0*sqr(v2) + 5.0*sqr(s))  - sqrt(sqr(mQlsq- mUrsq+0.125*sqr(g2)*(sqr(v1) - sqr(v2)) - 0.125*sqr(g1)*(sqr(v1) - sqr(v2))) + 4.0*mtop*mtop*sqr(At - lambda(3)*s*oneOrt2*v1/v2)));

  if (mstopsq(2) >= 0.0)
    { 
  mstop(2) = sqrt(0.5*(mQlsq + mUrsq + 2*mtop*mtop +0.125*sqr(g2)*(sqr(v1) - sqr(v2)) +0.125*3.0*sqr(g1)*(sqr(v1) - sqr(v2))/5.0 + oneO40*sqr(g1p)*(-3.0*sqr(v1) - 2.0*sqr(v2) + 5.0*sqr(s))  - sqrt(sqr(mQlsq- mUrsq+0.125*sqr(g2)*(sqr(v1) - sqr(v2)) - 0.125*sqr(g1)*(sqr(v1) - sqr(v2))) + 4.0*mtop*mtop*sqr(At - lambda(3)*s*oneOrt2*v1/v2))));
    }
  else
    {
      mstop(2) = -sqrt(Abs(mstopsq(2)));
    }
 
  if(speak){								   
    cout << " mstop(1) = " << mstop(1) << endl;
    cout << " mstop(2) = " << mstop(2) << endl;
  }
 
  int gen=1;
  for(gen=1; gen< 4; gen++){
    mD1sq(gen) = 0.5*(mDsq(gen) + mDbarsq(gen) + kappa(gen)*kappa(gen)*s*s - 2.5*oneO40*sqr(g1p)*(-3.0*sqr(v1) - 2.0*sqr(v2) + 5.0*sqr(s))  
		      - sqrt(sqr(mDsq(gen)- mDbarsq(gen) + 0.5*oneO40*sqr(g1p)*(-3.0*sqr(v1) - 2.0*sqr(v2) + 5.0*sqr(s))   + 0.1*sqr(g1)*(sqr(v1) - sqr(v2))) + 4.0*sqr(Akappa(gen)*kappa(gen)*s*oneOrt2 - 0.5*kappa(gen)*lambda(3)*v1*v2)));

    mD2sq(gen) = 0.5*(mDsq(gen) + mDbarsq(gen) + kappa(gen)*kappa(gen)*s*s - 2.5*oneO40*sqr(g1p)*(-3.0*sqr(v1) - 2.0*sqr(v2) + 5.0*sqr(s))  
		      + sqrt(sqr(mDsq(gen)- mDbarsq(gen) + 0.5*oneO40*sqr(g1p)*(-3.0*sqr(v1) - 2.0*sqr(v2) + 5.0*sqr(s))   + 0.1*sqr(g1)*(sqr(v1) - sqr(v2))) + 4.0*sqr(Akappa(gen)*kappa(gen)*s*oneOrt2 - 0.5*kappa(gen)*lambda(3)*v1*v2)));
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
  double Delmtop,Delmstop1,Delmstop2;  
  Delmtop = 0;
  
  Delmstop1 = 0.5*(0.5*( 8.0*Sqr(mtop)*(At - lambda(3)*s*oneOrt2*v1/v2)*(-lambda(3)*v1*oneOrt2/(v2*s) ) )) /(   Sqrt(Sqr(mQlsq- mUrsq)  + 4.0*mtop*mtop*Sqr(At - lambda(3)*s*oneOrt2*v1/v2))  );
  
  Delmstop2 = 0.5*( -0.5*( 8.0*Sqr(mtop)*(At - lambda(3)*s*oneOrt2*v1/v2)*(-lambda(3)*v1*oneOrt2/(v2*s) ) )) /(   Sqrt(Sqr(mQlsq- mUrsq)+ 4.0*mtop*mtop*Sqr(At - lambda(3)*s*oneOrt2*v1/v2))  );
  
  
  double delta =   threeO32pisq*(-2.0*a0Peter(mstopsq(1),q)*Delmstop1-2.0*a0Peter(mstopsq(2),q)*Delmstop2 + 4.0*a0Peter(Sqr(mtop),q)*Delmtop);
  //cout <<"delta = " << delta << endl;
 
  return delta;
}

void HiggsMasses(genericE6SSM_soft_parameters & r, double s, double tb, DoubleVector & mstop, DoubleVector & mstopsq, int WhatCorrections, bool speak, bool Bugspeak, DoubleVector & bounds, int & ExpValid, DoubleVector & mhout, DoubleMatrix & mhmix, DoubleMatrix & msq, int & sing) {

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
   if(speak) cout << "tachyonics stop masses so i'm not going to bother doing higgs mass corrections.  Point already ruled out. " << endl;

     // Dylan:: I have modified this so that in the event of tachyonic stops, the returned Higgs masses are negative
     // as a flag.
     mhout.set(1, -1.0);
     mhout.set(2, -1.0);
     mhout.set(3, -1.0);
     ExpValid = 30;
     return;
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
	 cout << " Ut = " << Ut << endl;
	 cout << "Should be like 1lp HE11 = " <<  HE11*(1 - 3*Sqr(r.get_Yu(2, 2))*l/(8.0*Sqr(PI)))+ 3.0 * Sqr(Sqr(r.get_Yu( 2, 2)) )*Sqr(vev)*Sqr(Sqr(Sin(ArcTan(tb))))/(8.0*Sqr(PI))*(0.5*Ut + l) << endl;
       }
  
     if(Bugspeak) cout << "Should be likeDelMh11 = " <<  HE11*( - 3*Sqr(r.get_Yu( 2, 2))*l)/(8.0*Sqr(PI))+ 3.0 * Sqr(Sqr(r.get_Yu( 2, 2)) )*Sqr(vev)*Sqr(Sqr(Sin(ArcTan(tb))))/(8.0*Sqr(PI))*(0.5*Ut + l) << endl;
     
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
	 cout << "tachyonics stop masses so i'm not going to bother doing higgs mass corrections.  Point already ruled out. " << endl;
	 
	 // Dylan:: I have modified this so that in the event of tachyonic stops, the returned Higgs masses are negative
	 // as a flag.
	 mhout.set(1, -1.0);
	 mhout.set(2, -1.0);
	 mhout.set(3, -1.0);
	 ExpValid = 30;
       }  
   }
 else
   {

     if (MH_even.diagonaliseSym(mixMH, mhphysq) >  TOLERANCE * 1.0e-3) 
       { 
	 cout << "Warning:  accuracy bad in CP-even Higgs diagonalisation" << endl;
	 // Flag to indicate poor accuracy
	 sing = 1;
       }
     DoubleVector mhphy(3);
     
     if (mhphysq(1) < 0. || mhphysq(2) < 0. || mhphysq(3) < 0.) 
       {
	 
	 ExpValid = 30;
	 
	 if (mhphysq(1) < 0.0)
	   {
	     mhphy.set(1, -Sqrt(-mhphysq(1)));
	   }
	 else
	   {
	     mhphy.set(1, Sqrt(mhphysq(1)));
	   }
	 
	 if (mhphysq(2) < 0.0)
	   {
	     mhphy.set(2, -Sqrt(-mhphysq(2)));
	   }
	 else
	   {
	     mhphy.set(2, Sqrt(mhphysq(2)));
	   }
	 
	 if (mhphysq(3) < 0.0)
	   {
	     mhphy.set(3, -Sqrt(-mhphysq(3)));
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
     
     mhout = mhphy;    
     mhmix = mixMH; // NOTE ADDITIONAL OUTPUT TO SAVE MIXING AS WELL
     msq = MH_even;
     
     if(speak)cout <<  "2lp light_higgs = " <<mhphy(1) << endl;
     if(speak)cout <<  "higgs = " <<mhphy << endl;
     
   }
 
 
 // Have been using 120 GeV and 130 GeV as limits. Note additional check to make sure 
 // not already flagged as tachyonic
 if (mhout.display(1) >= bounds.display(1) && mhout.display(1) <= bounds.display(2) && ExpValid != 30)
   {
     ExpValid = 0;
   }
 
 
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
double mAsq_OneLoop(genericE6SSM_soft_parameters essmSusy, double s, double tb)
{

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
  double mtop = 165;//yt*v2/Sqrt(2.0);
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

  double sigma11 = r.get_mHu2()-r.get_mHd2()+r.get_mHpbar2()-r.get_mHp2()+r.get_md2().trace()
    -r.get_mDx2().trace()+r.get_mDxbar2().trace()+r.get_me2().trace()-r.get_mH1I2().trace()
    +r.get_mH2I2().trace()-r.get_ml2().trace()+r.get_mq2().trace()-2.0*r.get_mu2().trace();

  double sigma14 = 2.0*QH1p*r.get_mHd2()+2.0*QH2p*r.get_mHu2()+2.0*QHpbarp*r.get_mHpbar2()
    +2.0*QHpp*r.get_mHp2()+QSp*r.get_ms2()+3.0*Qdp*r.get_md2().trace()+3.0*QDxp*r.get_mDx2().trace()
    +3.0*QDxbarp*r.get_mDxbar2().trace()+Qep*r.get_me2().trace()+2.0*QH1p*r.get_mH1I2().trace()
    +2.0*QH2p*r.get_mH2I2().trace()+2.0*QLp*r.get_ml2().trace()+6.0*QQp*r.get_mq2().trace()+QSp*r.get_msI2().trace()
    +3.0*Qup*r.get_mu2().trace();

  genericE6SSM_soft_parameters w = r.calc_beta();

  double beta_sigma11 = w.get_mHu2()-w.get_mHd2()+w.get_mHpbar2()-w.get_mHp2()+w.get_md2().trace()
    -w.get_mDx2().trace()+w.get_mDxbar2().trace()+w.get_me2().trace()-w.get_mH1I2().trace()
    +w.get_mH2I2().trace()-w.get_ml2().trace()+w.get_mq2().trace()-2.0*w.get_mu2().trace();

  double beta_sigma14 = 2.0*QH1p*w.get_mHd2()+2.0*QH2p*w.get_mHu2()+2.0*QHpbarp*w.get_mHpbar2()
    +2.0*QHpp*w.get_mHp2()+QSp*w.get_ms2()+3.0*Qdp*w.get_md2().trace()+3.0*QDxp*w.get_mDx2().trace()
    +3.0*QDxbarp*w.get_mDxbar2().trace()+Qep*w.get_me2().trace()+2.0*QH1p*w.get_mH1I2().trace()
    +2.0*QH2p*w.get_mH2I2().trace()+2.0*QLp*w.get_ml2().trace()+6.0*QQp*w.get_mq2().trace()+QSp*w.get_msI2().trace()
    +3.0*Qup*w.get_mu2().trace();

  double coeff = -2.4 * (r.get_g1()*w.get_g1()*Sqr(r.get_MassB()) + Sqr(r.get_g1())*r.get_MassB()*w.get_MassB());
  coeff = coeff - 16.0 * Sqr(QH2p) * (r.get_gN()*w.get_gN()*Sqr(r.get_MassBp()) + Sqr(r.get_gN())*r.get_MassBp()*w.get_MassBp());
  coeff = coeff - 12.0 * (r.get_g2()*w.get_g2()*Sqr(r.get_MassWB()) + Sqr(r.get_g2())*r.get_MassWB()*w.get_MassWB());
  coeff = coeff + 4.0 * r.get_Lambdax()*w.get_Lambdax()*(r.get_mHd2()+r.get_mHu2()+r.get_ms2());
  coeff = coeff + 2.0 * Sqr(r.get_Lambdax())*(w.get_mHd2()+w.get_mHu2()+w.get_ms2()) + 4.0 * r.get_TLambdax()*w.get_TLambdax();
  coeff = coeff + 12.0 * r.get_Yu(2,2)*w.get_Yu(2,2)*(r.get_mHu2()+r.get_mq2(2,2)+r.get_mu2(2,2));
  coeff = coeff + 6.0 * Sqr(r.get_Yu(2,2))*(w.get_mHu2()+w.get_mq2(2,2)+w.get_mu2(2,2)) + 12.0*r.get_TYu(2,2)*w.get_TYu(2,2);
  coeff = coeff + 1.2 * r.get_g1()*w.get_g1()*sigma11 + 0.6 * Sqr(r.get_g1())*beta_sigma11;
  coeff = coeff + 4.0 * QH2p * r.get_gN()*w.get_gN()*sigma14 + 2.0 * QH2p * Sqr(r.get_gN()) * beta_sigma14;

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

  double sigma11 = r.get_mHu2()-r.get_mHd2()+r.get_mHpbar2()-r.get_mHp2()+r.get_md2().trace()
    -r.get_mDx2().trace()+r.get_mDxbar2().trace()+r.get_me2().trace()-r.get_mH1I2().trace()
    +r.get_mH2I2().trace()-r.get_ml2().trace()+r.get_mq2().trace()-2.0*r.get_mu2().trace();

  double sigma14 = 2.0*QH1p*r.get_mHd2()+2.0*QH2p*r.get_mHu2()+2.0*QHpbarp*r.get_mHpbar2()
    +2.0*QHpp*r.get_mHp2()+QSp*r.get_ms2()+3.0*Qdp*r.get_md2().trace()+3.0*QDxp*r.get_mDx2().trace()
    +3.0*QDxbarp*r.get_mDxbar2().trace()+Qep*r.get_me2().trace()+2.0*QH1p*r.get_mH1I2().trace()
    +2.0*QH2p*r.get_mH2I2().trace()+2.0*QLp*r.get_ml2().trace()+6.0*QQp*r.get_mq2().trace()+QSp*r.get_msI2().trace()
    +3.0*Qup*r.get_mu2().trace();

  genericE6SSM_soft_parameters w = r.calc_beta();

  double beta_sigma11 = w.get_mHu2()-w.get_mHd2()+w.get_mHpbar2()-w.get_mHp2()+w.get_md2().trace()
    -w.get_mDx2().trace()+w.get_mDxbar2().trace()+w.get_me2().trace()-w.get_mH1I2().trace()
    +w.get_mH2I2().trace()-w.get_ml2().trace()+w.get_mq2().trace()-2.0*w.get_mu2().trace();

  double beta_sigma14 = 2.0*QH1p*w.get_mHd2()+2.0*QH2p*w.get_mHu2()+2.0*QHpbarp*w.get_mHpbar2()
    +2.0*QHpp*w.get_mHp2()+QSp*w.get_ms2()+3.0*Qdp*w.get_md2().trace()+3.0*QDxp*w.get_mDx2().trace()
    +3.0*QDxbarp*w.get_mDxbar2().trace()+Qep*w.get_me2().trace()+2.0*QH1p*w.get_mH1I2().trace()
    +2.0*QH2p*w.get_mH2I2().trace()+2.0*QLp*w.get_ml2().trace()+6.0*QQp*w.get_mq2().trace()+QSp*w.get_msI2().trace()
    +3.0*Qup*w.get_mu2().trace();

  double coeff = -2.4 * (r.get_g1()*w.get_g1()*Sqr(r.get_MassB()) + Sqr(r.get_g1())*r.get_MassB()*w.get_MassB());
  coeff = coeff - 16.0 * Sqr(QH1p) * (r.get_gN()*w.get_gN()*Sqr(r.get_MassBp()) + Sqr(r.get_gN())*r.get_MassBp()*w.get_MassBp());
  coeff = coeff - 12.0 * (r.get_g2()*w.get_g2()*Sqr(r.get_MassWB()) + Sqr(r.get_g2())*r.get_MassWB()*w.get_MassWB());
  coeff = coeff + 4.0 * r.get_Lambdax()*w.get_Lambdax()*(r.get_mHd2()+r.get_mHu2()+r.get_ms2());
  coeff = coeff + 2.0 * Sqr(r.get_Lambdax())*(w.get_mHd2()+w.get_mHu2()+w.get_ms2()) + 4.0 * r.get_TLambdax()*w.get_TLambdax();
  coeff = coeff + 12.0 * r.get_Yd(2,2)*w.get_Yd(2,2)*(r.get_mHd2()+r.get_mq2(2,2)+r.get_md2(2,2));
  coeff = coeff + 6.0 * Sqr(r.get_Yd(2,2))*(w.get_mHd2()+w.get_mq2(2,2)+w.get_md2(2,2)) + 12.0*r.get_TYd(2,2)*w.get_TYd(2,2);
  coeff = coeff + 4.0 * r.get_Ye(2,2)*w.get_Ye(2,2)*(r.get_mHd2()+r.get_ml2(2,2)+r.get_me2(2,2));
  coeff = coeff + 2.0 * Sqr(r.get_Ye(2,2))*(w.get_mHd2()+w.get_ml2(2,2)+w.get_me2(2,2)) + 4.0*r.get_TYe(2,2)*w.get_TYe(2,2);
  coeff = coeff - 1.2 * r.get_g1()*w.get_g1()*sigma11 - 0.6 * Sqr(r.get_g1())*beta_sigma11;
  coeff = coeff + 4.0 * QH1p * r.get_gN()*w.get_gN()*sigma14 + 2.0 * QH1p * Sqr(r.get_gN()) * beta_sigma14;

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

  // For consistency the correspond A terms should also vanish
  r.set_TYu(0,0,0.0);
  r.set_TYu(1,1,0.0);

  r.set_TYd(0,0,0.0);
  r.set_TYd(1,1,0.0);

  r.set_TYe(0,0,0.0);
  r.set_TYe(1,1,0.0);

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

  double sigma14 = 2.0*QH1p*r.get_mHd2()+2.0*QH2p*r.get_mHu2()+2.0*QHpbarp*r.get_mHpbar2()
    +2.0*QHpp*r.get_mHp2()+QSp*r.get_ms2()+3.0*Qdp*r.get_md2().trace()+3.0*QDxp*r.get_mDx2().trace()
    +3.0*QDxbarp*r.get_mDxbar2().trace()+Qep*r.get_me2().trace()+2.0*QH1p*r.get_mH1I2().trace()
    +2.0*QH2p*r.get_mH2I2().trace()+2.0*QLp*r.get_ml2().trace()+6.0*QQp*r.get_mq2().trace()+QSp*r.get_msI2().trace()
    +3.0*Qup*r.get_mu2().trace();

  genericE6SSM_soft_parameters w = r.calc_beta();

  double beta_sigma14 = 2.0*QH1p*w.get_mHd2()+2.0*QH2p*w.get_mHu2()+2.0*QHpbarp*w.get_mHpbar2()
    +2.0*QHpp*w.get_mHp2()+QSp*w.get_ms2()+3.0*Qdp*w.get_md2().trace()+3.0*QDxp*w.get_mDx2().trace()
    +3.0*QDxbarp*w.get_mDxbar2().trace()+Qep*w.get_me2().trace()+2.0*QH1p*w.get_mH1I2().trace()
    +2.0*QH2p*w.get_mH2I2().trace()+2.0*QLp*w.get_ml2().trace()+6.0*QQp*w.get_mq2().trace()+QSp*w.get_msI2().trace()
    +3.0*Qup*w.get_mu2().trace();


  double coeff = -16.0 * Sqr(QSp) * (r.get_gN()*w.get_gN()*Sqr(r.get_MassBp()) + Sqr(r.get_gN())*r.get_MassBp()*w.get_MassBp());
  coeff = coeff + 8.0 * r.get_Lambdax()*w.get_Lambdax()*(r.get_mHd2()+r.get_mHu2()+r.get_ms2());
  coeff = coeff + 4.0 * Sqr(r.get_Lambdax())*(w.get_mHd2()+w.get_mHu2()+w.get_ms2()) + 8.0 * r.get_TLambdax()*w.get_TLambdax();
  coeff = coeff + 6.0*w.get_ms2()*(r.get_TKappa()*(r.get_TKappa().adjoint())).trace() + 6.0 * r.get_ms2()*(w.get_TKappa()*(r.get_TKappa().adjoint())).trace() + 6.0*r.get_ms2()*(r.get_TKappa()*(w.get_TKappa().adjoint())).trace();
  coeff = coeff + 4.0*w.get_ms2()*(r.get_Lambda12()*(r.get_Lambda12().adjoint())).trace()+4.0*r.get_ms2()*(w.get_Lambda12()*(r.get_Lambda12().adjoint())).trace()+4.0*r.get_ms2()*(r.get_Lambda12()*(w.get_Lambda12().adjoint())).trace();
  coeff = coeff + 6.0*(w.get_TKappa()*(r.get_TKappa().adjoint())).trace()+6.0*(r.get_TKappa()*(w.get_TKappa().adjoint())).trace();
  coeff = coeff + 4.0*(w.get_TLambda12()*(r.get_TLambda12().adjoint())).trace()+4.0*(r.get_TLambda12()*(w.get_TLambda12().adjoint())).trace();
  coeff = coeff + 4.0*(w.get_mH1I2()*(r.get_Lambda12().adjoint())*r.get_Lambda12() + r.get_mH1I2()*(w.get_Lambda12().adjoint())*r.get_Lambda12() + r.get_mH1I2()*(r.get_Lambda12().adjoint())*w.get_Lambda12()).trace();
  coeff = coeff + 6.0*(w.get_Kappa()*(r.get_Kappa().adjoint())*r.get_mDx2() + r.get_Kappa()*(w.get_Kappa().adjoint())*r.get_mDx2() + r.get_Kappa()*(r.get_Kappa().adjoint())*w.get_mDx2()).trace();
  coeff = coeff + 6.0*(w.get_Kappa()*r.get_mDxbar2()*(r.get_Kappa().adjoint()) + r.get_Kappa()*w.get_mDxbar2()*(r.get_Kappa().adjoint()) + r.get_Kappa()*r.get_mDxbar2()*(w.get_Kappa().adjoint())).trace();
  coeff = coeff + 4.0*(w.get_Lambda12()*(r.get_Lambda12().adjoint())*r.get_mH2I2() + r.get_Lambda12()*(w.get_Lambda12().adjoint())*r.get_mH2I2() + r.get_Lambda12()*(r.get_Lambda12().adjoint())*w.get_mH2I2()).trace();
  coeff = coeff + 4.0 * QSp * r.get_gN()*w.get_gN()*sigma14 + 2.0 * QSp * Sqr(r.get_gN()) * beta_sigma14;

  coeff = 0.5*oneOver16PiSqr*coeff;

  return coeff;
}
