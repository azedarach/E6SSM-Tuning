/*
  Implementations of the functions defined in essmtuningutils.h
 */

#include "essmtuningutils.h"

// ESSM_EWSBConditioni, i = 1, 2, calculates the value of the 
// EWSB condition for m_H_i^2 and m_S^2.
// Inputs:
//     SoftParsMssm const & r = pMSSM model object
//     DoubleVector const & vevs = the values of v1 and v2
double ESSM_EWSBCondition1(flexiblesusy::genericE6SSM_soft_parameters const & r)
{
  using namespace flexiblesusy;

  double m1Sq = r.get_mHd2();

  double Tlambda = r.get_TLambdax();
  double lambda = r.get_Lambdax();
  double Alambda;

  if (Abs(Tlambda) < EPSTOL) Alambda = 0.0;

  if (Abs(lambda) < 1.0e-100)
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
	-doCalcTadpoleESSMH1(r);
    }

  return f1;

}

double ESSM_EWSBCondition2(flexiblesusy::genericE6SSM_soft_parameters const & r)
{
  using namespace flexiblesusy;

  double m2Sq = r.get_mHu2();

  double Tlambda = r.get_TLambdax();
  double lambda = r.get_Lambdax();
  double Alambda;

  if (Abs(Tlambda) < EPSTOL) Alambda = 0.0;

  if (Abs(lambda) < 1.0e-100)
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
	-doCalcTadpoleESSMH2(r);
    }

  return f2;
}

double ESSM_EWSBCondition3(flexiblesusy::genericE6SSM_soft_parameters const & r)
{
  using namespace flexiblesusy;

  double msSq = r.get_ms2();

  double Tlambda = r.get_TLambdax();
  double lambda = r.get_Lambdax();
  double Alambda;

  if (Abs(Tlambda) < EPSTOL) Alambda = 0.0;

  if (Abs(lambda) < 1.0e-100)
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
	-doCalcTadpolesESSMS(r);
    }
  
  return f3;

}
