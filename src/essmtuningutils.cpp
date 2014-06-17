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
double ESSM_EWSBCondition1(genericE6SSM_soft_parameters const & r)
{

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
	-doCalcTadpolesESSMS(r, s, tb);
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

  double oneOrt2 = 1/(Sqrt(2.0));
  double mtop = yt*v2/(Sqrt(2.0));

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
      ii << "WARNING: trying to calculate A_lambda(3) where lambda(3) coupling is " <<
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

//This version neglects U(1) D-terms and was written for comparison with Romans program
void physical_ESSM_Roman(genericE6SSM_soft_parameters r,DoubleVector & mstop, DoubleVector & mstopsq, DoubleVector & mD1sq, DoubleVector & mD2sq, double s, double tb){

  bool speak = false;
  double yt = r.get_Yu(3, 3);

  double v1 = r.get_vd();
  double v2 = r.get_vu();

  double vev = Sqrt(v1*v1+v2*v2);

  double oneOrt2 = 1/(Sqrt(2.0));
  double mtop = yt*v2/(Sqrt(2.0));
  mtop = 165;
  yt = mtop/(v2*oneOrt2);
  r.set_Yu(3, 3, yt);

  if(speak){
    cout << "mtop = " << mtop << endl;
    cout << "yt = " << yt << endl; 
  }
  
  double oneO40 = 1.0/(40.0);

  DoubleVector lambda(3), kappa(3), Akappa(3), Alambda(3), mDsq3(3), mDbarsq(3);
  
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
      ii << "WARNING: trying to calculate A_lambda(3) where lambda(3) coupling is " <<
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
  double yt = r.get_Yu(3, 3);

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
    r.set_Yu(3, 3, yt);
  }

  DoubleVector lambda(3), kappa(3), Akappa(3), Alambda(3), mDsq3(3), mDbarsq(3);

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
      ii << "WARNING: trying to calculate A_lambda(3) where lambda(3) coupling is " <<
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

  
  DoubleVector mstop(2), mstopsq(2);//stop1,stop2
  DoubleVector mD1sq(3), mD2sq(3);//3 gens

  double mQlsq = r.get_mq2(3,3);
  double mUrsq =  r.get_mu2(3,3);

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

  double yt = r.get_Yu(3, 3);

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
    r.set_Yu(3, 3, yt);    
  }

  DoubleVector lambda(3), kappa(3), Akappa(3), Alambda(3), mDsq3(3), mDbarsq(3);

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
      ii << "WARNING: trying to calculate A_lambda(3) where lambda(3) coupling is " <<
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

  DoubleVector mstop(2), mstopsq(2);//stop1,stop2
  DoubleVector mD1sq(3), mD2sq(3);//3 gens

  double mQlsq = r.get_mq2(3,3);
  double mUrsq =  r.get_mu2(3,3);

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
  double yt = r.get_Yu(3, 3);

  double v1 = r.get_vd();
  double v2 = r.get_vu();

  double vev = Sqrt(v1*v1+v2*v2);

  double oneOrt2 = 1/(Sqrt(2.0));
  double mtop = yt*v2*oneOrt2;
  mtop = 165; 
  yt = mtop/(v2*oneOrt2);
  r.set_Yu(3, 3, yt);

  double q = r.get_scale();
  q = 165; //GeV. Fudgeing to match Romans code

  DoubleVector lambda(3), kappa(3), Akappa(3), Alambda(3), mDsq3(3), mDbarsq(3);

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
      ii << "WARNING: trying to calculate A_lambda(3) where lambda(3) coupling is " <<
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

  DoubleVector mstop(2), mstopsq(2);//stop1,stop2
  DoubleVector mD1sq(3), mD2sq(3);//3 gens
  double mQlsq = r.get_mq2(3,3);
  double mUrsq =  r.get_mu2(3,3);
  
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
  double yt = r.get_Yu(3, 3);

  double v1 = r.get_vd();
  double v2 = r.get_vu();

  double vev = Sqrt(v1*v1+v2*v2);

  double oneOrt2 = 1/(Sqrt(2.0));
  double mtop = yt*v2*oneOrt2;
  //mtop = 165; 
  yt = mtop/(v2*oneOrt2);
  r.set_Yu(3, 3, yt);
  
  if(USEMTOFMT){
    mtop = 165;
    yt = mtop/(v2*oneOrt2);
    r.set_Yu(3, 3, yt);
  }
 
  double q = r.get_scale();
  //q = 165; //GeV. Fudgeing to match Romans code

  DoubleVector lambda(3), kappa(3), Akappa(3), Alambda(3), mDsq3(3), mDbarsq(3);

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
      ii << "WARNING: trying to calculate A_lambda(3) where lambda(3) coupling is " <<
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

  DoubleVector mstop(2), mstopsq(2);//stop1,stop2
  DoubleVector mD1sq(3), mD2sq(3);//3 gens
  double mQlsq = r.get_mq2(3,3);
  double mUrsq =  r.get_mu2(3,3);
  
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
   double yt = r.get_Yu(3, 3);
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
     r.set_Yu(3, 3, yt);
   }
 
   double threeO32pisq = 3.0/(32.0*Sqr(PI));
   //double oneOrt2 = 1/(Sqrt(2));			      
   double top = - 6.0 * Sqr(yt) * a0Peter(mtop, q)/ (16.0 * Sqr(PI));
   double topPeter = - 6.0 * Sqr(yt) * a0Peter(Sqr(mtop), q)/ (16.0 * Sqr(PI));
   double oneO40 = 1.0/(40.0);

   DoubleVector lambda(3), kappa(3), Akappa(3), Alambda(3), mDsq3(3), mDbarsq(3);
   
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
       ii << "WARNING: trying to calculate A_lambda(3) where lambda(3) coupling is " <<
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
   double Delmtop, Delmstop1,Delmstop2; 

   
   DoubleVector mstop(2), mstopsq(2);//stop1,stop2
   DoubleVector mD1sq(3), mD2sq(3);//3 gens
   
   physical_ESSM(r,mstop, mstopsq, mD1sq,mD2sq, s, tb); 
   Delmtop = Sqr(r.get_Yu(3,3));
 
   Delmstop1 = 0.5*(2.0*Sqr(r.get_Yu(3,3))- 0.6*0.25*Sqr(g1) -  0.25*Sqr(g2) - 0.1*Sqr(g1p)
		    +0.5*( 2.0*(mQlsq- mUrsq+0.125*Sqr(g2)*(Sqr(v1) - Sqr(v2)) - 0.125*Sqr(g1)*(Sqr(v1) - Sqr(v2)))*0.25*(Sqr(g1)- Sqr(g2)) + 4.0*Sqr( At - lambda(3)*s*oneOrt2*v1/v2)* Sqr(r.get_Yu(3,3)) + 8.0*Sqr(mtop)*(At - lambda(3)*s*oneOrt2*v1/v2)*(lambda(3)*s*v1*oneOrt2/(v2*v2*v2) )              ) /(   Sqrt(Sqr(mQlsq- mUrsq+0.125*Sqr(g2)*(Sqr(v1) - Sqr(v2)) - 0.125*Sqr(g1)*(Sqr(v1) - Sqr(v2))) + 4.0*mtop*mtop*Sqr(At - lambda(3)*s*oneOrt2*v1/v2))  ) );
 
 
   Delmstop2 = 0.5*(2.0*Sqr(r.get_Yu(3,3))- 0.6*0.25*Sqr(g1) -  0.25*Sqr(g2) - 0.1*Sqr(g1p)
		    -0.5*( 2.0*(mQlsq- mUrsq+0.125*Sqr(g2)*(Sqr(v1) - Sqr(v2)) - 0.125*Sqr(g1)*(Sqr(v1) - Sqr(v2)))*0.25*(Sqr(g1)- Sqr(g2)) + 4.0*Sqr( At - lambda(3)*s*oneOrt2*v1/v2)* Sqr(r.get_Yu(3,3)) + 8.0*Sqr(mtop)*(At - lambda(3)*s*oneOrt2*v1/v2)*(lambda(3)*s*v1*oneOrt2/(v2*v2*v2) )              ) /(   Sqrt(Sqr(mQlsq- mUrsq+0.125*Sqr(g2)*(Sqr(v1) - Sqr(v2)) - 0.125*Sqr(g1)*(Sqr(v1) - Sqr(v2))) + 4.0*mtop*mtop*Sqr(At - lambda(3)*s*oneOrt2*v1/v2))  ) );
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
  double yt = r.get_Yu(3,3);
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
    r.set_Yu(3, 3, yt);
  }
  
  double threeO32pisq = 3.0/(32.0*Sqr(PI));
  //double oneOrt2 = 1/(Sqrt(2));			      
  double oneO40 = 1.0/(40.0);
  
  DoubleVector lambda(3), kappa(3), Akappa(3), Alambda(3), mDsq3(3), mDbarsq(3);
  
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
      ii << "WARNING: trying to calculate A_lambda(3) where lambda(3) coupling is " <<
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

  double Delmtop, Delmstop1,Delmstop2; 

  double mQlsq = r.get_mq2(3,3);
  double mUrsq =  r.get_mu2(3,3);

  double g1 = r.get_g1();
  double g2 = r.get_g2();
  double g1p = r.get_gN();
  //getting stop masses
  DoubleVector mstop(2), mstopsq(2);//stop1,stop2
  DoubleVector mD1sq(3), mD2sq(3);//3 gens
  physical_ESSM(r,mstop, mstopsq, mD1sq,mD2sq, s, tb); 
  Delmtop = Sqr(r.get_Yu(3,3));
  // cout << "Delmtop = " << Delmtop << endl; 
  Delmstop1 = 0.5*(2.0*Sqr(r.get_Yu(3,3))- 0.6*0.25*Sqr(g1) -  0.25*Sqr(g2) - 0.1*Sqr(g1p)
		   +0.5*( 2.0*(mQlsq- mUrsq+0.125*Sqr(g2)*(Sqr(v1) - Sqr(v2)) - 0.125*Sqr(g1)*(Sqr(v1) - Sqr(v2)))*0.25*(Sqr(g1)- Sqr(g2)) + 4.0*Sqr( At - lambda(3)*s*oneOrt2*v1/v2)* Sqr(r.get_Yu(3,3)) + 8.0*Sqr(mtop)*(At - lambda(3)*s*oneOrt2*v1/v2)*(lambda(3)*s*v1*oneOrt2/(v2*v2*v2) )              ) /(   Sqrt(Sqr(mQlsq- mUrsq+0.125*Sqr(g2)*(Sqr(v1) - Sqr(v2)) - 0.125*Sqr(g1)*(Sqr(v1) - Sqr(v2))) + 4.0*mtop*mtop*Sqr(At - lambda(3)*s*oneOrt2*v1/v2))  ) );
      
  Delmstop2 = 0.5*(2.0*Sqr(r.get_Yu(3,3))- 0.6*0.25*Sqr(g1) -  0.25*Sqr(g2) - 0.1*Sqr(g1p)
		   -0.5*( 2.0*(mQlsq- mUrsq+0.125*Sqr(g2)*(Sqr(v1) - Sqr(v2)) - 0.125*Sqr(g1)*(Sqr(v1) - Sqr(v2)))*0.25*(Sqr(g1)- Sqr(g2)) + 4.0*Sqr( At - lambda(3)*s*oneOrt2*v1/v2)* Sqr(r.get_Yu(3,3)) + 8.0*Sqr(mtop)*(At - lambda(3)*s*oneOrt2*v1/v2)*(lambda(3)*s*v1*oneOrt2/(v2*v2*v2) )              ) /(   Sqrt(Sqr(mQlsq- mUrsq+0.125*Sqr(g2)*(Sqr(v1) - Sqr(v2)) - 0.125*Sqr(g1)*(Sqr(v1) - Sqr(v2))) + 4.0*mtop*mtop*Sqr(At - lambda(3)*s*oneOrt2*v1/v2))  ) );
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
  double yt = r.get_Yu(3,3);
  // cout << "yt = " << yt << endl;
  double v1 = r.get_vd()
    double v2 = r.get_vu();
 
  double vev = Sqrt(v1*v1+v2*v2);

  double oneOrt2 = 1/(Sqrt(2.0));
  double mtop = yt*v2/(Sqrt(2.0));
  mtop = 165;
  yt = mtop/(v2*oneOrt2);
  r.set_Yu(3,3,yt);
  double q = r.get_scale();
  q = 165; //GeV. Fi=udgeing to match Romans code
  double threeO32pisq = 3.0/(32.0*Sqr(PI));
  //double oneOrt2 = 1/(Sqrt(2));			      

  double oneO40 = 1.0/(40.0);

  DoubleVector lambda(3), kappa(3), Akappa(3), Alambda(3), mDsq3(3), mDbarsq(3);
  
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
      ii << "WARNING: trying to calculate A_lambda(3) where lambda(3) coupling is " <<
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
  double Delmtop, Delmstop1,Delmstop2; 
  //   cout << "in H2 tads" << endl; 
  
  //getting stop masses
  DoubleVector mstop(2), mstopsq(2);//stop1,stop2
  DoubleVector mD1sq(3), mD2sq(3);//3 gens
  physical_ESSM_Roman(r,mstop, mstopsq, mD1sq,mD2sq, s, tb); 
  Delmtop = Sqr(r.get_Yu(3,3));
  //cout << "Delmtop = " << Delmtop << endl; 
  
  Delmstop1 = 0.5*(2.0*Sqr(r.get_Yu(3,3))  +0.5*(  4.0*Sqr( At - lambda(3)*s*oneOrt2*v1/v2)* Sqr(r.get_Yu(3,3))+ 8.0*Sqr(mtop)*(At - lambda(3)*s*oneOrt2*v1/v2)*(lambda(3)*s*v1*oneOrt2/(v2*v2*v2) )             ) /(   Sqrt(Sqr(mQlsq- mUrsq) + 4.0*mtop*mtop*Sqr(At - lambda(3)*s*oneOrt2*v1/v2))  ) ) ;
  
  Delmstop2 = 0.5*(2.0*Sqr(r.get_Yu(3,3))
		   -0.5*(  4.0*Sqr( At - lambda(3)*s*oneOrt2*v1/v2)* Sqr(r.get_Yu(3,3)) + 8.0*Sqr(mtop)*(At - lambda(3)*s*oneOrt2*v1/v2)*(lambda(3)*s*v1*oneOrt2/(v2*v2*v2) )              ) /(   Sqrt(Sqr(mQlsq- mUrsq)+ 4.0*mtop*mtop*Sqr(At - lambda(3)*s*oneOrt2*v1/v2))       ));
  
  
  //checked and debuugged to this point 28/6/07 
  double delta =   threeO32pisq*(-2.0*a0Peter(mstopsq(1),q)*Delmstop1-2.0*a0Peter(mstopsq(2),q)*Delmstop2 + 4.0*a0Peter(Sqr(mtop),q)*Delmtop);
  //cout <<"delta = " << delta << endl;
  return delta;
}

double doCalcTadpoleESSMH2_Roman_atQ(genericE6SSM_soft_parameters r, double s , double tb  ){
  double yt = r.get_Yu(3,3);
  // cout << "yt = " << yt << endl; 
  double v1 = r.get_vd();
  double v2 = r.get_vu();

  double vev = Sqrt(v1*v1+v2*v2);//246; //NOTE CHANGE

  double oneOrt2 = 1/(Sqrt(2.0));
  double mtop = yt*v2/(Sqrt(2.0));
  //  mtop = 165;
  yt = mtop/(v2*oneOrt2);
  r.set_Yu(3,3,yt);
  double q = r.get_scale();
  
  if(USEMTOFMT){
    cout << "using mtop = " << mtop << endl; 
    mtop = 165;
    yt = mtop/(v2*oneOrt2);
    r.set_Yu(3,3,yt);
    
  }
  
  double threeO32pisq = 3.0/(32.0*Sqr(PI));
  //double oneOrt2 = 1/(Sqrt(2));			      

  double oneO40 = 1.0/(40.0);

  DoubleVector lambda(3), kappa(3), Akappa(3), Alambda(3), mDsq3(3), mDbarsq(3);
  
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
      ii << "WARNING: trying to calculate A_lambda(3) where lambda(3) coupling is " <<
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

  double Delmtop, Delmstop1,Delmstop2; 
  double mQlsq = r.get_mq2(3,3);
  double mUrsq =  r.get_mu2(3,3);
  
  double g1 = r.get_g1();
  double g2 = r.get_g2();
  double g1p = r.get_gN();

  
  DoubleVector mstop(2), mstopsq(2);//stop1,stop2
  DoubleVector mD1sq(3), mD2sq(3);//3 gens
  
  physical_ESSM_Roman(r,mstop, mstopsq, mD1sq,mD2sq, s, tb);
   
  Delmtop = Sqr(r.get_Yu(3,3));
  //cout << "Delmtop = " << Delmtop << endl; 
  
  Delmstop1 = 0.5*(2.0*Sqr(r.get_Yu(3,3))  +0.5*(  4.0*Sqr( At - lambda(3)*s*oneOrt2*v1/v2)* Sqr(r.get_Yu(3,3))+ 8.0*Sqr(mtop)*(At - lambda(3)*s*oneOrt2*v1/v2)*(lambda(3)*s*v1*oneOrt2/(v2*v2*v2) )             ) /(   Sqrt(Sqr(mQlsq- mUrsq) + 4.0*mtop*mtop*Sqr(At - lambda(3)*s*oneOrt2*v1/v2))  ) ) ;
  
  Delmstop2 = 0.5*(2.0*Sqr(r.get_Yu(3,3))
		   -0.5*(  4.0*Sqr( At - lambda(3)*s*oneOrt2*v1/v2)* Sqr(r.get_Yu(3,3)) + 8.0*Sqr(mtop)*(At - lambda(3)*s*oneOrt2*v1/v2)*(lambda(3)*s*v1*oneOrt2/(v2*v2*v2) )              ) /(   Sqrt(Sqr(mQlsq- mUrsq)+ 4.0*mtop*mtop*Sqr(At - lambda(3)*s*oneOrt2*v1/v2))       ));
  
  
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
  double yt = r.get_Yu(3,3);

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
    r.set_Yu(3,3,yt);
  }

  DoubleVector lambda(3), kappa(3), Akappa(3), Alambda(3), mDsq3(3), mDbarsq(3);
  
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
      ii << "WARNING: trying to calculate A_lambda(3) where lambda(3) coupling is " <<
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
  double yt = r.get_Yu(3,3);

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
    r.set_Yu(3,3,yt);
    
  }

  if(speak){
   cout << "q = " << q << endl;
  cout << "mtop = " << mtop << endl;
  }


  DoubleVector lambda(3), kappa(3), Akappa(3), Alambda(3), mDsq3(3), mDbarsq(3);
  
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
      ii << "WARNING: trying to calculate A_lambda(3) where lambda(3) coupling is " <<
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
  double yt = r.get_Yu(3,3);

  double v1 = r.get_vd();
  double v2 = r.get_vu();

  double vev = Sqrt(v1*v1+v2*v2);

  ////double oneOrt2 = 1/(Sqrt(2.0));
  double mtop = yt*v2*oneOrt2;
  mtop = 165;
  
  yt = mtop/(v2*oneOrt2);
  r.set_Yu(3,3,yt);
  //cout << "mtop = " << mtop << endl;
  double q = r.get_scale();
  q = 165; //GeV. Fudgeing to match Romans code
  //cout << "q = " << q << endl;
  
  DoubleVector lambda(3), kappa(3), Akappa(3), Alambda(3), mDsq3(3), mDbarsq(3);
  
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
      ii << "WARNING: trying to calculate A_lambda(3) where lambda(3) coupling is " <<
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
  double yt = r.get_Yu(3,3);

  double v1 = r.get_vd();
  double v2 = r.get_vu();

  double vev = Sqrt(v1*v1+v2*v2);

  ////double oneOrt2 = 1/(Sqrt(2.0));
  double mtop = yt*v2*oneOrt2;
  //mtop = 165;
  
  yt = mtop/(v2*oneOrt2);
  r.set_Yu(3,3,yt);
  //cout << "mtop = " << mtop << endl;
  double q = r.get_scale();
  //q = 165; //GeV. Fudgeing to match Romans code
  cout << "q = " << q << endl;
  
  if(USEMTOFMT){
    cout << "using mtop = " << mtop << endl; 
    mtop = 165;
    yt = mtop/(v2*oneOrt2);
    r.set_Yu(3,3,yt);
  }
  
  DoubleVector lambda(3), kappa(3), Akappa(3), Alambda(3), mDsq3(3), mDbarsq(3);
  
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
      ii << "WARNING: trying to calculate A_lambda(3) where lambda(3) coupling is " <<
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

  double mQLsq = r.get_mq2(3,3);
  double mURsq = r.get_mu2(3,3);
  double g2 = r.get_g2();
  double g1 = r.get_g1();

  double yt = r.get_Yu(3,3); //cout << "yt = " << yt << endl;
  double g1p = r.get_gN();
  double mu_eff = mylambda(3)*s/(Sqrt(2.0));
  double gbar = Sqrt(3.0*Sqr(g1)/5.0 + Sqr(g2)) ; 


  double v1 = r.get_vd();
  double v2 = r.get_vu();
 
  double vev = Sqrt(v1*v1+v2*v2);

  double oneOrt2 = 1.0/Sqrt(2.0);

  DoubleVector mylambda(3), mykappa(3), myAkappa(3), myAlambda(3), mDsq3(3), mDbarsq(3);
  
  mylambda(3) = r.get_Lambdax();
  
  for (int i = 1; i <= 2; i++)
    {
      mylambda(i) = r.get_Lambda12(i,i);
    }
  
  for (int i = 1; i <= 3; i++)
    {
      mykappa(i) = r.get_Kappa(i,i);
      mDsq(i) = r.get_mDx2(i,i);
      mDbarsq(i) = r.get_mDxbar2(i,i);
    }
  
  
  DoubleVector Tlambda(3), Tkappa(3);
  
  Tlambda(3) = r.get_TLambdax();
  
  for (int i = 1; i <= 2; i++)
    {
      Tlambda(i) = r.get_TLambda12(i,i);
      
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
	Abs(mylambda) << endl;
      throw ii.str();
    }
  else
    {
      myAlambda(3) = Tlambda(3)/mylambda(3);
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
	  myAkappa(i) = Tkappa(i,i)/mykappa(i);
	}
    }


  if( WhatCorrections ==2 || WhatCorrections == 4){
    r.set_scale(165);
  }

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
 r.set_Yu(3,3, yt);
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
 
 
 DelMh11prime = 3.0/(32.0*Sqr(PI))*(2.0*(2.0*log(mstop(1)/r.get_scale())*Sqr( mstop1derriv_1 )  + a0Peter(mstopsq(1),r.get_scale())* DoubleStop1derriv11)   +2.0*(2.0*log(mstop(2)/r.get_scale())*Sqr( mstop2derriv_1 )  + a0Peter(mstopsq(2),r.get_scale())* DoubleStop2derriv11));
 
 DelMh12prime = 3.0/(32.0*Sqr(PI))*(2.0*(2.0*log(mstop(1)/r.get_scale())*(mstop1derriv_1*mstop1derriv_2)+  a0Peter(mstopsq(1),r.get_scale())* DoubleStop1derriv12) + 2.0*( 2.0*log(mstop(2)/r.get_scale())*mstop2derriv_1*mstop2derriv_2 +  a0Peter(mstopsq(2),r.get_scale())* DoubleStop2derriv12));
 
 DelMh13prime = 3.0/(32.0*Sqr(PI))*(2.0*(2.0*log(mstop(1)/r.get_scale())*(mstop1derriv_1*mstop1derriv_s)+  a0Peter(mstopsq(1),r.get_scale())* DoubleStop1derriv13) + 2.0*( 2.0*log(mstop(2)/r.get_scale())*mstop2derriv_1*mstop2derriv_s +  a0Peter(mstopsq(2),r.get_scale())* DoubleStop2derriv13));
 
 DelMh22prime = 3.0/(32.0*Sqr(PI))*(2.0*(2.0*log(mstop(1)/r.get_scale())*Sqr( mstop1derriv_2 )  + a0Peter(mstopsq(1),r.get_scale())* DoubleStop1derriv22)   +2.0*(2.0*log(mstop(2)/r.get_scale())*Sqr( mstop2derriv_2 )  + a0Peter(mstopsq(2),r.get_scale())* DoubleStop2derriv22) - 4.0*(2.0*log(mtop/r.get_scale())*Sqr(Sqr(yt))*Sqr(v2) +a0Peter(Sqr(mtop),r.get_scale())*Sqr(yt)));
 
 DelMh23prime = 3.0/(32.0*Sqr(PI))*(2.0*(2.0*log(mstop(1)/r.get_scale())*(mstop1derriv_2*mstop1derriv_s) +  a0Peter(mstopsq(1), r.get_scale())* DoubleStop1derriv23) + 2.0*( 2.0*log(mstop(2)/r.get_scale())*mstop2derriv_2*mstop2derriv_s +  a0Peter(mstopsq(2),r.get_scale())* DoubleStop2derriv23));
 
 //cout << " DelMh23prime = " <<  DelMh23prime  << endl;
 DelMh33prime = 3.0/(32.0*Sqr(PI))*(2.0*(2.0*log(mstop(1)/r.get_scale())*Sqr( mstop1derriv_s )  + a0Peter(mstopsq(1),r.get_scale())* DoubleStop1derriv33)   +2.0*(2.0*log(mstop(2)/r.get_scale())*Sqr( mstop2derriv_s )  + a0Peter(mstopsq(2),r.get_scale())* DoubleStop2derriv33));
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
 DelMh11  = Sqr(cos(atan(tb))) * DelMh11prime + 2*tb/(1+tb*tb)*DelMh12prime + Sqr(sin(atan(tb)))* DelMh22prime;
 DelMh22  = Sqr(sin(atan(tb))) * DelMh11prime - 2*tb/(1+tb*tb)*DelMh12prime + Sqr(cos(atan(tb)))* DelMh22prime;
 DelMh33 = DelMh33prime;
 DelMh12  = (Sqr(cos(atan(tb))) -Sqr(sin(atan(tb)))) * DelMh12prime  + tb/(1+tb*tb)*( DelMh22prime - DelMh11prime);
 DelMh13 = cos(atan(tb))* DelMh13prime + sin(atan(tb))* DelMh23prime; 
 DelMh23 = cos(atan(tb))* DelMh23prime - sin(atan(tb))* DelMh13prime; 

 //entries
 double HE11 = 2.0*Sqr(mylambda(3))*v1*v2*tb/(1 + tb*tb) + 0.25*Sqr(gbar)*(Sqr(v1) - Sqr(v2))*(1-tb*tb)/(1 + tb*tb) + Sqr(g1p)*Sqr(v1*Q1eff*cos(atan(tb)) + v2*sin(atan(tb))*Q2eff) + DelMh11;
 
 const bool Includeleadtwoloop = true; // True = use two loop. Should use two loop in final results?
 if(Includeleadtwoloop)
   {
     //Remove one loop because we will replace it with equivelent terms obtained from rge evolution
     HE11 = HE11 - DelMh11;
     //if(speak) cout << "Tree (hopefully!) HE11 = " << HE11 << endl;
     double l = log(mstop(1)*mstop(2)/(Sqr(mtop)));
     //if(speak) cout << " l = " << l << endl;  
     double Xt = At - mu_eff/tb;
     //if(speak)cout << "Xt = " << Xt << endl;
     double Ut = 2.0*Sqr(Xt)/(mstop(1)*mstop(2))*(1 - Sqr(Xt)/(12.0*mstop(1)*mstop(2)));
     if(Bugspeak)
       { 
	 cout << " Ut = " << Ut << endl;
	 cout << "Should be like 1lp HE11 = " <<  HE11*(1 - 3*Sqr(r.displayYukawaElement(YU, 3, 3))*l/(8.0*Sqr(PI)))+ 3.0 * Sqr(Sqr(r.displayYukawaElement(YU, 3, 3)) )*Sqr(vev)*Sqr(Sqr(sin(atan(tb))))/(8.0*Sqr(PI))*(0.5*Ut + l) << endl;
       }
  
     if(Bugspeak) cout << "Should be likeDelMh11 = " <<  HE11*( - 3*Sqr(r.displayYukawaElement(YU, 3, 3))*l)/(8.0*Sqr(PI))+ 3.0 * Sqr(Sqr(r.displayYukawaElement(YU, 3, 3)) )*Sqr(vev)*Sqr(Sqr(sin(atan(tb))))/(8.0*Sqr(PI))*(0.5*Ut + l) << endl;
     
     HE11 = HE11*(1 - 3*Sqr(r.displayYukawaElement(YU, 3, 3))*l/(8.0*Sqr(PI)))+ 3.0 * Sqr(Sqr(r.displayYukawaElement(YU, 3, 3)) )*Sqr(vev)*Sqr(Sqr(sin(atan(tb))))/(8.0*Sqr(PI))*(0.5*Ut + l + 1/(16*Sqr(PI))*(1.5*Sqr(r.displayYukawaElement(YU, 3, 3)) - 8.0*Sqr(r.displayGaugeCoupling(3)))*(Ut + l)*l);
     
     //if(speak) cout << "After 2lp HE11 = " << HE11 << endl; 
     
   }
 double HE12 = (0.25*Sqr(mylambda(3)) - 0.125*Sqr(gbar))*Sqr(vev)* 4.0*tb*(1 - Sqr(tb))/(Sqr(1.0+tb*tb)) +   Sqr(g1p)*( Sqr(v1)*Q1eff + Sqr(v2)*Q2eff)*(Q2eff - Q1eff)*tb/(1+tb*tb)  + DelMh12;
 double HE13 = -mylambda(3)*myAlambda(3)*v2*Sqrt(2)*cos(atan(tb)) + Sqr( mylambda(3))*s*vev +  Sqr(g1p)*( Sqr(v1)*Q1eff + Sqr(v2)*Q2eff)*QSeff*s/(vev) + DelMh13;
 double HE22 = mylambda(3)*myAlambda(3)*Sqrt(2.0)*s*(1+tb*tb)/(2.0*tb) - (0.5*Sqr(mylambda(3)) - 0.25*Sqr(gbar))*Sqr(vev)* 4.0*tb*tb/(Sqr(1.0+tb*tb)) +Sqr(g1p)*Sqr(vev)*tb*tb/(Sqr(1.0+tb*tb))*Sqr(Q2eff - Q1eff) + DelMh22 ;
 double HE23 = - mylambda(3)*myAlambda(3)/(Sqrt(2))*vev*(1 - tb*tb)/(tb*tb +1.0) + Sqr(g1p)*(Q2eff - Q1eff)*QSeff*v2*s*cos(atan(tb)) + DelMh23;
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

  double lambda = r.get_Lambdax();
  double Tlambda;
  double Alambda;
  
  Tlambda = r.get_TLambdax();
  
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
  double yt = essmSusy.get_Yu(3, 3);
  
  double lambda = r.get_Lambdax();
  double Tlambda;
  double Alambda;
  
  Tlambda = r.get_TLambdax();
  
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

  double mtop = yt*v2/(Sqrt(2.0));
  double oneOrt2 = 1.0/Sqrt(2.0);

  if(USEMTOFMT){
    mtop = 165;
    yt = mtop/(v2*oneOrt2);
    essmSusy.set_Yu(3, 3, yt);    
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
  double fVal = (1.0/(m1Sq-m2Sq))*(m1Sq*log(m1Sq/(Q*Q))-m2Sq*log(m2Sq/(Q*Q))-m1Sq+m2Sq);

  return fVal;
}


