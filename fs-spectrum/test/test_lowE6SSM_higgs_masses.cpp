// ====================================================================
// Test suite for methods used to calculate the Higgs mass in 
// E6 models with general charges
// ====================================================================

#include <limits>
#include <random>
#include <chrono>
#include <Eigen/Core>

#include "ew_input.hpp"
#include "wrappers.hpp"
#include "lowE6SSM_two_scale_model.hpp"
#include "lowE6SSM_two_scale_ew_derivs.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_lowE6SSM_higgs_masses

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

using namespace std;
using namespace flexiblesusy;
using namespace softsusy;

int const UFBPROBLEM = 1; //< unbounded from below flag
int const CCBPROBLEM = 2; //< charge/colour breaking minimum flag

int const EWSBPROBLEM = 3; //< problem with iteration

int const WRONGVACUUM = 4; //< problem with unphysical vacuum

int const SQUARKTACHYON = 5; //< stop/sbottom tachyon
int const VECTORBOSONTACHYON = 6; //< vector boson tachyon
int const TADPOLESPROBLEM = 15; //< problem calculating 1-loop tadpoles

int const HIGGSPROBLEM = 10; //< problems with calculating the Higgs mass
int const NOTEXPVALID = 30; //< point not experimentally valid
int const HIGGSTACHYON = 31; //< ruled out because tachyonic
int const POLEHIGGSTACHYON = 32; //< ruled out because physical Higgs is tachyon

int const NUMERICALPROBLEM = 666; //< for serious numerical problems

int const TUNINGERROR = 35;

double const HIGGSCENT = 125.0; //< rough central value for Higgs mass, GeV
double const HIGGSERROR = 7.0; //< theory error in Higgs calculation, GeV

bool const USEMTOFMT = false;

lowE6SSM_input_parameters get_test_inputs(double theta)
{
   lowE6SSM_input_parameters input;

   input.TanBeta = 10.0;

   input.QQp = lowE6SSM_info::get_e6_charge(lowE6SSM_info::QL, theta);
   input.QLp = lowE6SSM_info::get_e6_charge(lowE6SSM_info::LL, theta);
   input.QH1p = lowE6SSM_info::get_e6_charge(lowE6SSM_info::Hd, theta);
   input.QH2p = lowE6SSM_info::get_e6_charge(lowE6SSM_info::Hu, theta);
   input.Qdp = lowE6SSM_info::get_e6_charge(lowE6SSM_info::dR, theta);
   input.Qup = lowE6SSM_info::get_e6_charge(lowE6SSM_info::uR, theta);
   input.Qep = lowE6SSM_info::get_e6_charge(lowE6SSM_info::eR, theta);
   input.QSp = lowE6SSM_info::get_e6_charge(lowE6SSM_info::S, theta);
   input.QDxp = lowE6SSM_info::get_e6_charge(lowE6SSM_info::Dx, theta);
   input.QDxbarp = lowE6SSM_info::get_e6_charge(lowE6SSM_info::Dxbar, theta);
   input.QHpp = lowE6SSM_info::get_e6_charge(lowE6SSM_info::Hp, theta);
   input.QHpbarp = lowE6SSM_info::get_e6_charge(lowE6SSM_info::Hpbar, theta);

   input.KappaInput(0,0) = 0.5;
   input.KappaInput(1,1) = 0.5;
   input.KappaInput(2,2) = 0.5;

   input.Lambda12Input(0,0) = 0.1;
   input.Lambda12Input(1,1) = 0.1;

   input.LambdaxInput = 0.5;

   input.MuPrInput = 10000.;

   input.gNInput = 0.45;

   input.vsInput = 8000.;

   input.TYdInput(0,0) = 0.;
   input.TYdInput(1,1) = 0.;
   input.TYdInput(2,2) = 5000.;

   input.TYeInput(0,0) = 0.;
   input.TYeInput(1,1) = 0.;
   input.TYeInput(2,2) = 5000.;

   input.TKappaInput(0,0) = 1000.;
   input.TKappaInput(1,1) = 1000.;
   input.TKappaInput(2,2) = 1000.;

   input.TLambda12Input(0,0) = 200.;
   input.TLambda12Input(1,1) = 200.;

   input.TLambdaxInput = 5000.;

   input.TYuInput(0,0) = 0.;
   input.TYuInput(1,1) = 0.;
   input.TYuInput(2,2) = 2500.;

   input.BMuPrInput = 10000.;

   input.mq2Input(0,0) = 2.5e7;
   input.mq2Input(1,1) = 2.5e7;
   input.mq2Input(2,2) = 1.5e6;

   input.ml2Input(0,0) = 2.5e7;
   input.ml2Input(1,1) = 2.5e7;
   input.ml2Input(2,2) = 2.5e7;

   input.md2Input(0,0) = 2.5e7;
   input.md2Input(1,1) = 2.5e7;
   input.md2Input(2,2) = 2.5e7;

   input.mu2Input(0,0) = 2.5e7;
   input.mu2Input(1,1) = 2.5e7;
   input.mu2Input(2,2) = 1.5e6;

   input.me2Input(0,0) = 2.5e7;
   input.me2Input(1,1) = 2.5e7;
   input.me2Input(2,2) = 2.5e7;

   input.mH1I2Input(0,0) = 2.5e7;
   input.mH1I2Input(1,1) = 2.5e7;

   input.mH2I2Input(0,0) = 2.5e7;
   input.mH2I2Input(1,1) = 2.5e7;

   input.msI2Input(0,0) = 2.5e7;
   input.msI2Input(1,1) = 2.5e7;

   input.mDx2Input(0,0) = 2.5e7;
   input.mDx2Input(1,1) = 2.5e7;
   input.mDx2Input(2,2) = 2.5e7;

   input.mDxbar2Input(0,0) = 2.5e7;
   input.mDxbar2Input(1,1) = 2.5e7;
   input.mDxbar2Input(2,2) = 2.5e7;

   input.mHp2Input = 7000.0;
   input.mHpbar2Input = 6000.0;
  
   input.MassBInput = 100.;
   input.MassWBInput = 1000.;
   input.MassGInput = 2000.;
   input.MassBpInput = 100.;

   return input;
}

void set_susy_parameters_from_input(const lowE6SSM_input_parameters& input, 
                                    lowE6SSM<Two_scale>& model)
{
   model.set_input_parameters(input);

   model.set_Yd(0, 0, 0.);
   model.set_Yd(1, 1, 0.);
   model.set_Yd(2, 2, 0.131960806);

   model.set_Ye(0, 0, 0.);
   model.set_Ye(1, 1, 0.);
   model.set_Ye(2, 2, 0.0985967845);

   model.set_Kappa(input.KappaInput);
   model.set_Lambda12(input.Lambda12Input);
   model.set_Lambdax(input.LambdaxInput);

   model.set_Yu(0, 0, 0.);
   model.set_Yu(1, 1, 0.);
   model.set_Yu(2, 2, 0.886476119);

   model.set_MuPr(input.MuPrInput);
   model.set_g1(0.4710277548);
   model.set_g2(0.634883796);
   model.set_g3(1.01228547);
   model.set_gN(input.gNInput);
   model.set_vd(25.3566930);
   model.set_vu(240.478577);
   model.set_vs(input.vsInput);

}

void set_soft_parameters_from_input(const lowE6SSM_input_parameters& input,
                                    lowE6SSM<Two_scale>& model)
{
   model.set_TYd(input.TYdInput);
   model.set_TYe(input.TYeInput);
   model.set_TKappa(input.TKappaInput);
   model.set_TLambda12(input.TLambda12Input);
   model.set_TLambdax(input.TLambdaxInput);
   model.set_TYu(input.TYuInput);
   model.set_BMuPr(input.BMuPrInput);
   model.set_mq2(input.mq2Input);
   model.set_ml2(input.ml2Input);
   model.set_mHd2(2.61451136e8);
   model.set_mHu2(-4.76353957e6);
   model.set_md2(input.md2Input);
   model.set_mu2(input.mu2Input);
   model.set_me2(input.me2Input);
   model.set_ms2(-7.20751463e6);
   model.set_mH1I2(input.mH1I2Input);
   model.set_mH2I2(input.mH2I2Input);
   model.set_msI2(input.msI2Input);
   model.set_mDx2(input.mDx2Input);
   model.set_mDxbar2(input.mDxbar2Input);
   model.set_mHp2(input.mHp2Input);
   model.set_mHpbar2(input.mHpbar2Input);
   model.set_MassB(input.MassBInput);
   model.set_MassWB(input.MassWBInput);
   model.set_MassG(input.MassGInput);
   model.set_MassBp(input.MassBpInput);

}

void initialize_model(const lowE6SSM_input_parameters& input, 
                      lowE6SSM<Two_scale>& model)
{
   model.set_scale(2.66356559e3);
   set_susy_parameters_from_input(input, model);
   set_soft_parameters_from_input(input, model);
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
void physical_ESSM(lowE6SSM<Two_scale> r,DoubleVector & mstop, DoubleVector & mstopsq, DoubleVector & mD1sq, DoubleVector & mD2sq, double s, double tb)
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

  if (mstopsq(1) >= 0.0)
    {
  mstop(1) = Sqrt(0.5*(mQlsq + mUrsq + 2*mtop*mtop + 0.125*Sqr(g2)*(Sqr(v1) - Sqr(v2)) + 0.125*3.0*Sqr(g1)*(Sqr(v1) - Sqr(v2))/5.0 + oneO40*Sqr(g1p)*(-3.0*Sqr(v1) - 2.0*Sqr(v2) + 5.0*Sqr(s))  + Sqrt(Sqr(mQlsq- mUrsq+0.125*Sqr(g2)*(Sqr(v1) - Sqr(v2)) - 0.125*Sqr(g1)*(Sqr(v1) - Sqr(v2))) + 4.0*mtop*mtop*Sqr(At - lambda(3)*s*oneOrt2*v1/v2))));
    }
  else
    {
      mstop(1) = -Sqrt(Abs(mstopsq(1)));
    }

  mstopsq(2) = 0.5*(mQlsq + mUrsq + 2*mtop*mtop +0.125*Sqr(g2)*(Sqr(v1) - Sqr(v2)) +0.125*3.0*Sqr(g1)*(Sqr(v1) - Sqr(v2))/5.0 + oneO40*Sqr(g1p)*(-3.0*Sqr(v1) - 2.0*Sqr(v2) + 5.0*Sqr(s))  - Sqrt(Sqr(mQlsq- mUrsq+0.125*Sqr(g2)*(Sqr(v1) - Sqr(v2)) - 0.125*Sqr(g1)*(Sqr(v1) - Sqr(v2))) + 4.0*mtop*mtop*Sqr(At - lambda(3)*s*oneOrt2*v1/v2)));

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
void physical_ESSM_Roman(lowE6SSM<Two_scale> r,DoubleVector & mstop, DoubleVector & mstopsq, DoubleVector & mD1sq, DoubleVector & mD2sq, double s, double tb){

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
double doCalcTadpoleESSMH1(lowE6SSM<Two_scale> r,  double s , double tb)
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
double doCalcTadpoleESSMH1_atMt(lowE6SSM<Two_scale> r,  double s , double tb){
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
double doCalcTadpoleESSMH1_Roman(lowE6SSM<Two_scale> r,  double s , double tb){
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
double doCalcTadpoleESSMH1_Roman_atQ(lowE6SSM<Two_scale> r,  double s , double tb){
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
double doCalcTadpoleESSMH2(lowE6SSM<Two_scale> r, double s , double tb  ){
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
 
double doCalcTadpoleESSMH2_atMt(lowE6SSM<Two_scale> r, double s , double tb  ){
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

double doCalcTadpoleESSMH2_Roman(lowE6SSM<Two_scale> r, double s , double tb  ){
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

double doCalcTadpoleESSMH2_Roman_atQ(lowE6SSM<Two_scale> r, double s , double tb  ){
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
double doCalcTadpolesESSMS( lowE6SSM<Two_scale> r, double s , double tb  ){
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

double doCalcTadpolesESSMS_atMt( lowE6SSM<Two_scale> r, double s , double tb  )  {
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
   

double doCalcTadpolesESSMS_Roman( lowE6SSM<Two_scale> r, double s , double tb  )  {
     
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


double doCalcTadpolesESSMS_Roman_atQ( lowE6SSM<Two_scale> r, double s , double tb  ){
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

// Peter's Higgs code
bool HiggsMasses(lowE6SSM<Two_scale> const & model, double s, double tb, DoubleVector & mstop, DoubleVector & mstopsq, int WhatCorrections, bool speak, bool Bugspeak, DoubleVector & bounds, int & ExpValid, DoubleVector & mhout, DoubleMatrix & mhmix, DoubleMatrix & msq, int & sing) {


  // Copy object for calculation
   lowE6SSM<Two_scale> r = model;

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

  lowE6SSM_input_parameters input = r.get_input();

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
 r.set_TYu(2,2,yt*At);
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
 // cout << "Tadpole1 = " << doCalcTadpoleESSMH1(r,s,tb) << endl;
 // cout << "Tadpole2 = " << doCalcTadpoleESSMH2(r,s,tb) << endl;
 // cout << "Tadpole3 = " << doCalcTadpolesESSMS(r,s,tb) << endl;
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
     if(speak) cout << "Tree (hopefully!) HE11 = " << HE11 << endl;
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

// Higgs mass matrix agrees at tree level so far, two loop turned off,
// testing one loop
BOOST_AUTO_TEST_CASE( test_ewsb_tree_level_soln_soft_masses )
{
   const double tol = 1.0e-10;
   const double theta = ArcTan(Sqrt(15.));
   
   lowE6SSM_input_parameters inputs = get_test_inputs(theta);

   lowE6SSM<Two_scale> model;

   initialize_model(inputs, model);

   model.set_ewsb_loop_order(2);
   model.set_pole_mass_loop_order(2);

   lowE6SSM_ew_derivs ew_derivs(model);

   bool fs_tachyon = ew_derivs.calculate_MHiggs();
   
   DoubleVector mstop(2), mstopsq(2), bounds(2), mhout(3);
   DoubleMatrix mhmix(3,3), msq(3,3);

   bounds(1) = 120.0;
   bounds(2) = 130.0;

   const int what_corrections = 1;
   const bool speak = false;
   const bool bugspeak = false;
   int exp_valid;
   int sing;

   bool softsusy_tachyon = HiggsMasses(model, model.get_vs(), model.get_vu() / model.get_vd(), mstop, 
                                      mstopsq, what_corrections, speak, bugspeak, bounds, exp_valid, 
                                      mhout, mhmix, msq, sing);

   BOOST_CHECK_EQUAL(fs_tachyon, softsusy_tachyon);
   BOOST_CHECK_LE(Abs(ew_derivs.get_MHiggs()(0) - mhout(1)), tol);
   BOOST_CHECK_LE(Abs(ew_derivs.get_MHiggs()(1) - mhout(2)), tol);
   BOOST_CHECK_LE(Abs(ew_derivs.get_MHiggs()(2) - mhout(3)), tol);

}
