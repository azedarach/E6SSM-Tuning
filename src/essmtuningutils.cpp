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

      genericE6SSM_input_parameters input = r.get_input();
      
      // Because we are neglecting U(1) mixing, 
      // we will for now approximate the effective charges
      // by the U(1)' charges
      double Qtilde_1 = input.QH1p;
      double Qtilde_2 = input.QH2p;
      double Qtilde_s = input.QSp;
      
      
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


//Should be part of softsusy but you may not want to link to the file which
// uses this.  If you do just comment this out.  
double ccbSqrt(double f){ return Sqrt(Abs(f)); }

// Function a0 defined in paper.
// Note sign difference to paper - a0 term
// appears with opposite sign below for some reason.
double a0Peter(double mSq, double Q)
{
  return mSq*(log(mSq/(Q*Q))-1);
}

// Calculates the masses of the stops and exotic quarks for tadpoles.  the former are always included in our tadpole contributions, the latter could be, but provide small contributions in comparison to errors due to threshold treatment.  
void physical_ESSM(genericE6SSM_soft_parameters r,DoubleVector & mstop, DoubleVector & mstopsq, DoubleVector & mD1sq, DoubleVector & mD2sq, double s, double tb)
{
  bool speak = false;

  double yt = r.get_Yu(3, 3);

  double v1 = r.get_vd();
  double v2 = r.get_vu();

  double vev = Sqrt(v1*v1+v2*v2);

  double oneOrt2 = 1/(sqrt(2.0));
  double mtop = yt*v2/(sqrt(2.0));

  if(USEMTOFMT){
    mtop = 165;
    yt = mtop/(v2*oneOrt2);
    r.set_Yu(3, 3, yt);
    
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
      lambda(i) = r.get_Lambda12(i,i);
    }

  for (int i = 1; i <= 3; i++)
    {
      kappa(i) = r.get_Kappa(i,i);
      mDsq(i) = r.get_mDx2(i,i);
      mDbarsq(i) = r.get_mDxbar2(i,i);
    }
 
  DoubleVector Akappa(3);
  DoubleVector Alambda(3);

  DoubleVector Tlambda(3), Tkappa(3);

  Tlambda(3) = r.get_TLambdax();

  for (int i = 1; i <= 2; i++)
    {
      Tlambda(i) = r.get_TLambda12(i,i);

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
      ii << "WARNING: trying to calculate A_lambda where lambda3 coupling is " <<
	Abs(lambda) << endl;
      throw ii.str();
    }
  else
    {
      Alambda(3) = Tlambda(3)/lambda(3);
    }

  double At;
  double TYt = r.get_TYu(3,3);

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
      Tkappa(i) = r.get_Tkappa(i,i);
      
      if (AbsTkappa(i,i) < EPSTOL)
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
	  Akappa(i) = Tkappa(i,i)/kappa(i);
	}
    }

  double mQlsq = r.get_mq2(3,3);
  double mUrsq =  r.get_mu2(3,3);

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
      mstop(1) = -sqrt(fabs(mstopsq(1)));
    }

  mstopsq(2) = 0.5*(mQlsq + mUrsq + 2*mtop*mtop +0.125*sqr(g2)*(sqr(v1) - sqr(v2)) +0.125*3.0*sqr(g1)*(sqr(v1) - sqr(v2))/5.0 + oneO40*sqr(g1p)*(-3.0*sqr(v1) - 2.0*sqr(v2) + 5.0*sqr(s))  - sqrt(sqr(mQlsq- mUrsq+0.125*sqr(g2)*(sqr(v1) - sqr(v2)) - 0.125*sqr(g1)*(sqr(v1) - sqr(v2))) + 4.0*mtop*mtop*sqr(At - lambda(3)*s*oneOrt2*v1/v2)));

  if (mstopsq(2) >= 0.0)
    { 
  mstop(2) = sqrt(0.5*(mQlsq + mUrsq + 2*mtop*mtop +0.125*sqr(g2)*(sqr(v1) - sqr(v2)) +0.125*3.0*sqr(g1)*(sqr(v1) - sqr(v2))/5.0 + oneO40*sqr(g1p)*(-3.0*sqr(v1) - 2.0*sqr(v2) + 5.0*sqr(s))  - sqrt(sqr(mQlsq- mUrsq+0.125*sqr(g2)*(sqr(v1) - sqr(v2)) - 0.125*sqr(g1)*(sqr(v1) - sqr(v2))) + 4.0*mtop*mtop*sqr(At - lambda(3)*s*oneOrt2*v1/v2))));
    }
  else
    {
      mstop(2) = -sqrt(fabs(mstopsq(2)));
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
