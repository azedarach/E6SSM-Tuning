/*
  Scan over E6SSM parameter space and output the Higgs spectrum at
  each point. Note that *all* parameters are scanned linearly at
  the moment, to do log scans either need to modify this file
  or just do multiple linear scans to cover the ranges you want
  in more detail.
 */

#include <iostream>
#include <cstdlib>
#include <map>

#include "essmtuningutils.h"

/*
  --------------------------------------------------------------
  E6SSM U(1)_\psi and U(1)_\chi charges
  --------------------------------------------------------------
*/

const double QQPsi = 1.0/(2.0*flexiblesusy::Sqrt(6.0));
const double QuPsi = 1.0/(2.0*flexiblesusy::Sqrt(6.0));
const double QdPsi = 1.0/(2.0*flexiblesusy::Sqrt(6.0));
const double QLPsi = 1.0/(2.0*flexiblesusy::Sqrt(6.0));
const double QePsi = 1.0/(2.0*flexiblesusy::Sqrt(6.0));
const double QSPsi = 4.0/(2.0*flexiblesusy::Sqrt(6.0));
const double QH1Psi = -1.0/flexiblesusy::Sqrt(6.0);
const double QH2Psi = -1.0/flexiblesusy::Sqrt(6.0);
const double QXPsi = -1.0/flexiblesusy::Sqrt(6.0);
const double QXbarPsi = -1.0/flexiblesusy::Sqrt(6.0);
const double QHPrPsi = 1.0/(2.0*flexiblesusy::Sqrt(6.0));
const double QHbarPrPsi = -1.0/(2.0*flexiblesusy::Sqrt(6.0));

const double QQChi = -1.0/(2.0*flexiblesusy::Sqrt(10.0));
const double QuChi = -1.0/(2.0*flexiblesusy::Sqrt(10.0));
const double QdChi = 3.0/(2.0*flexiblesusy::Sqrt(10.0));
const double QLChi = 3.0/(2.0*flexiblesusy::Sqrt(10.0));
const double QeChi = -1.0/(2.0*flexiblesusy::Sqrt(10.0));
const double QSChi = 0.0;
const double QH1Chi = -1.0/flexiblesusy::Sqrt(10.0);
const double QH2Chi = 1.0/flexiblesusy::Sqrt(10.0);
const double QXChi = 1.0/flexiblesusy::Sqrt(10.0);
const double QXbarChi = -1.0/flexiblesusy::Sqrt(10.0);
const double QHPrChi = 3.0/(2.0*flexiblesusy::Sqrt(10.0));
const double QHbarPrChi = -3.0/(2.0*flexiblesusy::Sqrt(10.0));

// This code only applies to a specific E6 type model
const double thetaE6 = flexiblesusy::ArcTan(flexiblesusy::Sqrt(15.0));


/*
  --------------------------------------------------------------
  For switching tan(beta) values, it is convenient to store the
  benchmark SM couplings and VEVs for the different tan(beta)
  values we consider. It is probably easiest just to do this in 
  a struct.
  --------------------------------------------------------------
*/
struct genericE6SSM_tanb_benchmark {
  double yt;
  double yb;
  double ytau;
  double g1;
  double g2;
  double g3;
  double gN;
  double vd;
  double vu;
  double vs;
  double q;

  genericE6SSM_tanb_benchmark()
    : yt(0), yb(0), ytau(0), g1(0), 
      g2(0), g3(0), gN(0), vd(0), vu(0), vs(0), q(0)
  {}

  genericE6SSM_tanb_benchmark(double yt_, double yb_, double ytau_, double g1_, 
			      double g2_, double g3_, double gN_, double vd_, 
			      double vu_, double vs_, double q_)
    : yt(yt_), yb(yb_), ytau(ytau_), g1(g1_), g2(g2_), g3(g3_), 
      gN(gN_), vd(vd_), vu(vu_), vs(vs_), q(q_)
  {}

};

int main(int argc, const char* argv[])
{
  using namespace flexiblesusy;
  using namespace essm_tuning_utils;
  typedef Two_scale algorithm_type;
  
  /*
    --------------------------------------------------------------
    SM Parameters
    --------------------------------------------------------------
  
  double poleM_t = 173.2; // top quark pole mass in GeV
  double mbmb = 4.16; // bottom quark running mass m_b(m_b)^MSbar in GeV
  double poleM_tau = 1.777;// tau pole mass in GeV
  double poleM_Z = 91.1876; // Z boson pole mass in GeV
  double G_F = 1.16637900e-5; // G_F^MSbar in GeV^-2.
  double alphasmz = 0.1193; // alpha_s(mz)^MSbar (strong coupling constant = g_3^2/(4*PI))
  double alphaemmzinv = 127.9568; // 1 / alpha_em(mz)^MSbar (electromagnetic coupling = e^2/(4*PI))

  
    --------------------------------------------------------------
  */
  
  /*
    --------------------------------------------------------------
    Benchmark tan(beta) values
    --------------------------------------------------------------
  */


  // M_Z' ~ 2.5 TeV benchmarks
  std::map<double, genericE6SSM_tanb_benchmark> tanb_benchmarks_lightZp;

  tanb_benchmarks_lightZp[2] = genericE6SSM_tanb_benchmark(0.959124964,0.0400071651,0.0216468839,0.4863671017,0.655920190,1.03253774,0.490545769,114.499068,206.720817,6473.40002,20000.0);
  tanb_benchmarks_lightZp[4] = genericE6SSM_tanb_benchmark(0.859951732,0.0734441866,0.0399972936,0.48621451,0.655880073,1.03211622,0.490468747,62.0254231,227.861139,6488.98714,20000.0);
  tanb_benchmarks_lightZp[6] = genericE6SSM_tanb_benchmark(0.841945783,0.108333702,0.0591237075,0.4861821009,0.655834175,1.03207227,0.490446235,41.9982217,232.330999,6490.78019,20000.0);
  tanb_benchmarks_lightZp[8] = genericE6SSM_tanb_benchmark(0.835707235,0.143738028,0.0785402143,0.4861642477,0.655801403,1.03206204,0.490431588,31.6439425,233.943395,6491.35943,20000.0);
  tanb_benchmarks_lightZp[10] = genericE6SSM_tanb_benchmark(0.832883879,0.179506721,0.0981595687,0.4861513481,0.655780111,1.03206201,0.490419902,25.3426343,234.696923,6491.58828,20000.0);
  tanb_benchmarks_lightZp[15] = genericE6SSM_tanb_benchmark(0.830271928,0.270491511,0.148048266,0.4861273795,0.655766615,1.03208667,0.490397035,16.8431261,235.444017,6491.72748,20000.0);
  tanb_benchmarks_lightZp[20] = genericE6SSM_tanb_benchmark(0.829606947,0.364129251,0.199340931,0.4861097226,0.655803693,1.03215025,0.490380135,12.5410697,235.702334,6491.68545,20000.0);
  tanb_benchmarks_lightZp[30] = genericE6SSM_tanb_benchmark(0.830223730,0.561326513,0.307827785,0.4861051809,0.656077992,1.03256720,0.490378735,8.16984550,235.845120,6491.53142,20000.0);
  tanb_benchmarks_lightZp[40] = genericE6SSM_tanb_benchmark(0.837056784,0.797297515,0.423453064,0.4865421102,0.660715078,1.03234294,0.490892040,5.90104956,235.398004,6488.05664,20000.0);
  tanb_benchmarks_lightZp[50] = genericE6SSM_tanb_benchmark(0.809680220,1.03497910,0.555805753,0.4869302761,0.656906847,1.09553359,0.490973612,4.52472815,236.857502,6689.94285,20000.0);

  /*
    --------------------------------------------------------------
   */

  /*
    --------------------------------------------------------------
    Input scale
    --------------------------------------------------------------
  */

  const double MX = 20000.0; // GeV

  /*
    --------------------------------------------------------------
  */

    /*
    --------------------------------------------------------------
    Values for parameters that are not scanned over
    --------------------------------------------------------------
   */

  // Select tan(beta) value
  double tanb = 10.0;

  // Get values from saved benchmarks
  std::map<double,genericE6SSM_tanb_benchmark>::iterator benchmark = tanb_benchmarks_lightZp.find(tanb);

  if (benchmark == tanb_benchmarks_lightZp.end()) 
    {
      cerr << "ERROR: no benchmark for tan(beta) = " << tanb << endl;
      return 1;
    }

  // SM Yukawas, assumed diagonal and
  // 1st and 2nd gen vanishing
  Eigen::Matrix<double,3,3> Yu, Yd, Ye;

  Yu(0,0) = 0.0;
  Yu(0,1) = 0.0;
  Yu(0,2) = 0.0;
  Yu(1,0) = 0.0;
  Yu(1,1) = 0.0;
  Yu(1,2) = 0.0;
  Yu(2,0) = 0.0;
  Yu(2,1) = 0.0;
  Yu(2,2) = benchmark->second.yt;

  Yd(0,0) = 0.0;
  Yd(0,1) = 0.0;
  Yd(0,2) = 0.0;
  Yd(1,0) = 0.0;
  Yd(1,1) = 0.0;
  Yd(1,2) = 0.0;
  Yd(2,0) = 0.0;
  Yd(2,1) = 0.0;
  Yd(2,2) = benchmark->second.yb;

  Ye(0,0) = 0.0;
  Ye(0,1) = 0.0;
  Ye(0,2) = 0.0;
  Ye(1,0) = 0.0;
  Ye(1,1) = 0.0;
  Ye(1,2) = 0.0;
  Ye(2,0) = 0.0;
  Ye(2,1) = 0.0;
  Ye(2,2) = benchmark->second.ytau;

  // SM gauge couplings
  double g1 = benchmark->second.g1;
  double g2 = benchmark->second.g2;
  double g3 = benchmark->second.g3;
  double g1p = benchmark->second.gN;

  // Higgs and singlet VEVs
  double v = Sqrt(Sqr(benchmark->second.vd)+Sqr(benchmark->second.vu)); // GeV
  double v1 = v/Sqrt(1.0+tanb*tanb);
  double v2 = v1*tanb; 
  double vs = benchmark->second.vs; // GeV

  // Gaugino soft masses
  double M1 = 300.0; // GeV
  double M1p = 300.0; // GeV
  double M3 = 2000.0; // GeV
  
  // 3rd-gen down and lepton soft A-terms
  double Ab = 0.0; // GeV
  double Atau = 0.0; // GeV
  double At = 0.0; // GeV, we set this later in the scan

  // Soft SM A-terms
  Eigen::Matrix<double,3,3> TYu, TYd, TYe;

  TYu(0,0) = 0.0;
  TYu(0,1) = 0.0;
  TYu(0,2) = 0.0;
  TYu(1,0) = 0.0;
  TYu(1,1) = 0.0;
  TYu(1,2) = 0.0;
  TYu(2,0) = 0.0;
  TYu(2,1) = 0.0;
  TYu(2,2) = Yu(2,2)*At;

  TYd(0,0) = 0.0;
  TYd(0,1) = 0.0;
  TYd(0,2) = 0.0;
  TYd(1,0) = 0.0;
  TYd(1,1) = 0.0;
  TYd(1,2) = 0.0;
  TYd(2,0) = 0.0;
  TYd(2,1) = 0.0;
  TYd(2,2) = Yd(2,2)*Ab;

  TYe(0,0) = 0.0;
  TYe(0,1) = 0.0;
  TYe(0,2) = 0.0;
  TYe(1,0) = 0.0;
  TYe(1,1) = 0.0;
  TYe(1,2) = 0.0;
  TYe(2,0) = 0.0;
  TYe(2,1) = 0.0;
  TYe(2,2) = Ye(2,2)*Atau;

  // Soft left-handed slepton masses
  double meL = 5000.0; // GeV
  double mmuL = 5000.0; // GeV
  double mtauL = 5000.0; // GeV
  Eigen::Matrix<double,3,3> mLlSq;
  mLlSq(0,0) = Sqr(meL);
  mLlSq(0,1) = 0.0;
  mLlSq(0,2) = 0.0;
  mLlSq(1,0) = 0.0;
  mLlSq(1,1) = Sqr(mmuL);
  mLlSq(1,2) = 0.0;
  mLlSq(2,0) = 0.0;
  mLlSq(2,1) = 0.0;
  mLlSq(2,2) = Sqr(mtauL);  

  // Soft right-handed slepton masses
  double meR = 5000.0; // GeV
  double mmuR = 5000.0; // GeV
  double mtauR = 5000.0; // GeV
  Eigen::Matrix<double,3,3> mErSq;
  mErSq(0,0) = Sqr(meR);
  mErSq(0,1) = 0.0;
  mErSq(0,2) = 0.0;
  mErSq(1,0) = 0.0;
  mErSq(1,1) = Sqr(mmuR);
  mErSq(1,2) = 0.0;
  mErSq(2,0) = 0.0;
  mErSq(2,1) = 0.0;
  mErSq(2,2) = Sqr(mtauR);  

  // Soft left-handed squark masses
  double mqL1 = 5000.0; // GeV
  double mqL2 = 5000.0; // GeV
  Eigen::Matrix<double,3,3> mQlSq;
  mQlSq(0,0) = Sqr(mqL1);
  mQlSq(0,1) = 0.0;
  mQlSq(0,2) = 0.0;
  mQlSq(1,0) = 0.0;
  mQlSq(1,1) = Sqr(mqL2);
  mQlSq(2,0) = 0.0;
  mQlSq(2,1) = 0.0;
  mQlSq(2,2) = 0.0; //< we set this in the scan later

  // Soft right-handed up-squark masses
  double muR = 5000.0; // GeV
  double mcR = 5000.0; // GeV
  Eigen::Matrix<double,3,3> mUrSq;
  mUrSq(0,0) = Sqr(muR);
  mUrSq(0,1) = 0.0;
  mUrSq(0,2) = 0.0;
  mUrSq(1,0) = 0.0;
  mUrSq(1,1) = Sqr(mcR);
  mUrSq(1,2) = 0.0;
  mUrSq(2,0) = 0.0;
  mUrSq(2,1) = 0.0;
  mUrSq(2,2) = 0.0; //< we set this in the scan later

  // Soft right-handed down-squark masses
  double mdR = 5000.0; // GeV
  double msR = 5000.0; // GeV
  double mbR = 5000.0; // GeV
  Eigen::Matrix<double,3,3> mDrSq;
  mDrSq(0,0) = Sqr(mdR);
  mDrSq(0,1) = 0.0;
  mDrSq(0,2) = 0.0;
  mDrSq(1,0) = 0.0;
  mDrSq(1,1) = Sqr(msR);
  mDrSq(1,2) = 0.0;
  mDrSq(2,0) = 0.0;
  mDrSq(2,1) = 0.0;
  mDrSq(2,2) = Sqr(mbR);

  // SUSY singlet-D-Dbar trilinear couplings
  double universal_kappa = 0.1;
  double kappa1 = universal_kappa;
  double kappa2 = universal_kappa;
  double kappa3 = universal_kappa;
  Eigen::Matrix<double,3,3> Kappa;
  Kappa(0,0) = kappa1;
  Kappa(0,1) = 0.0;
  Kappa(0,2) = 0.0;
  Kappa(1,0) = 0.0;
  Kappa(1,1) = kappa2;
  Kappa(1,2) = 0.0;
  Kappa(2,0) = 0.0;
  Kappa(2,1) = 0.0;
  Kappa(2,2) = kappa3; 

  // Soft singlet-D-Dbar A-terms
  double Akappa1 = 0.0; // GeV 
  double Akappa2 = 0.0; // GeV 
  double Akappa3 = 0.0; // GeV 
  Eigen::Matrix<double,3,3> TKappa;
  TKappa(0,0) = Kappa(0,0) * Akappa1;
  TKappa(0,1) = 0.0;
  TKappa(0,2) = 0.0;
  TKappa(1,0) = 0.0;
  TKappa(1,1) = Kappa(1,1) * Akappa2;
  TKappa(1,2) = 0.0;
  TKappa(2,0) = 0.0;
  TKappa(2,1) = 0.0;
  TKappa(2,2) = Kappa(2,2) * Akappa3;

  // SUSY singlet-inert H1-inert H2 trilinear couplings
  double lambda1 = 0.2; 
  double lambda2 = 0.2; 
  Eigen::Matrix<double,2,2> Lambda12;
  Lambda12(0,0) = lambda1;
  Lambda12(0,1) = 0.0;
  Lambda12(1,0) = 0.0;
  Lambda12(1,1) = lambda2;

  // Soft singlet-inert H1-inert H2 A-terms
  double Alambda1 = 0.0; // GeV 
  double Alambda2 = 0.0; // GeV 
  Eigen::Matrix<double,2,2> TLambda12;
  TLambda12(0,0) = Lambda12(0,0) * Alambda1;
  TLambda12(0,1) = 0.0;
  TLambda12(1,0) = 0.0;
  TLambda12(1,1) = Lambda12(1,1) * Alambda2; 

  // SUSY H'-H'bar bilinear coupling
  double mupr = 5000.0; // GeV 
  // Soft H'-H'bar bilinear coupling
  double Bmupr = 5000.0; // GeV^2 

  // Soft inert H1 scalar squared masses
  double mH11Sq = Sqr(5000.0); // GeV^2 
  double mH12Sq = Sqr(5000.0); // GeV^2 
  Eigen::Matrix<double,2,2> mHdISq;
  mHdISq(0,0) = mH11Sq;
  mHdISq(0,1) = 0.0;
  mHdISq(1,0) = 0.0;
  mHdISq(1,1) = mH12Sq;

  // Soft inert H2 scalar squared masses
  double mH21Sq = Sqr(5000.0); // GeV^2 
  double mH22Sq = Sqr(5000.0); // GeV^2 
  Eigen::Matrix<double,2,2> mHuISq;
  mHuISq(0,0) = mH21Sq;
  mHuISq(0,1) = 0.0;
  mHuISq(1,0) = 0.0;
  mHuISq(1,1) = mH22Sq;

  // Soft inert singlet scalar squared masses
  double mS1Sq = Sqr(5000.0); // GeV^2 
  double mS2Sq = Sqr(5000.0); // GeV^2 
  Eigen::Matrix<double,2,2> mSISq;
  mSISq(0,0) = mS1Sq;
  mSISq(0,1) = 0.0;
  mSISq(1,0) = 0.0;
  mSISq(1,1) = mS2Sq;

  // Soft H' scalar squared mass
  double mHpSq = Sqr(5000.0); // GeV^2
  // Soft H'bar scalar squared mass 
  double mHpbarSq = Sqr(5000.0); // GeV^2 

  // Soft D scalar squared masses
  double mDSq = Sqr(5000.0); // GeV^2
  Eigen::Matrix<double,3,3> mDxSq;
  mDxSq(0,0) = mDSq;
  mDxSq(0,1) = 0.0;
  mDxSq(0,2) = 0.0;
  mDxSq(1,0) = 0.0;
  mDxSq(1,1) = mDSq;
  mDxSq(1,2) = 0.0;
  mDxSq(2,0) = 0.0;
  mDxSq(2,1) = 0.0;
  mDxSq(2,2) = mDSq;

  // Soft Dbar scalar squared masses
  double mDbarSq = Sqr(5000.0); // GeV^2
  Eigen::Matrix<double,3,3> mDxbarSq;
  mDxbarSq(0,0) = mDbarSq;
  mDxbarSq(0,1) = 0.0;
  mDxbarSq(0,2) = 0.0;
  mDxbarSq(1,0) = 0.0;
  mDxbarSq(1,1) = mDbarSq;
  mDxbarSq(1,2) = 0.0;
  mDxbarSq(2,0) = 0.0;
  mDxbarSq(2,1) = 0.0;
  mDxbarSq(2,2) = mDbarSq;

  /*
    --------------------------------------------------------------
   */

  /*
    ----------------------------------------------------
    Scan parameters - upper and lower bounds, 
    number of points in each parameter direction,
    and increment step size.
    ----------------------------------------------------
   */
  
  // Number of points in each direction
  const unsigned lambda3_npts = 1;
  const unsigned Alambda3_npts = 10;
  const unsigned mqL3sq_npts = 5;
  const unsigned mtRsq_npts = 5;
  const unsigned At_npts = 5;
  const unsigned M2_npts = 1;
  
  // Lower bounds
  const double lambda3_low = 1.1;
  const double Alambda3_low = 30000.0; // GeV
  const double At_low = -1000.0; // GeV
  const double mqL3_low = 2000.0; // GeV
  const double mtR_low = 2000.0; // GeV
  
  const double mqL3sq_low = Sqr(mqL3_low);
  const double mtRsq_low = Sqr(mtR_low);
  
  const double M2_low = 100.0; // GeV

  // Upper bounds
  const double lambda3_up = 3.0; // GeV
  const double Alambda3_up = 50000.0; // GeV
  const double At_up = 10000.0; // GeV
  const double mqL3_up = 2000.0; // GeV
  const double mtR_up = 2000.0; // GeV
  
  const double mqL3sq_up = Sqr(mqL3_up);
  const double mtRsq_up = Sqr(mtR_up);
  
  const double M2_up = 2000.0; // GeV

  /*
    ----------------------------------------------------
   */

  /*
    ----------------------------------------------------
    Do the scan and write to standard output
    ----------------------------------------------------
  */

  // Work out increments
  double par_incr = 0.0;
  
  if (lambda3_npts == 1)
    {
      par_incr = 0.0;
    }
  else
    {
      par_incr = (lambda3_up-lambda3_low)/(lambda3_npts-1.0);
    }

  const double lambda3_incr = par_incr;

  par_incr = 0.0;
  if (Alambda3_npts == 1)
    {
      par_incr = 0.0;
    }
  else
    {
      par_incr = (Alambda3_up-Alambda3_low)/(Alambda3_npts-1.0);
    }

  const double Alambda3_incr = par_incr;

  par_incr = 0.0;
  if (mqL3sq_npts == 1)
    {
      par_incr = 0.0;
    }
  else
    {
      par_incr = (mqL3sq_up-mqL3sq_low)/(mqL3sq_npts-1.0);
    }

  const double mqL3sq_incr = par_incr;

  par_incr = 0.0;
  if (mtRsq_npts == 1)
    {
      par_incr = 0.0;
    }
  else
    {
      par_incr = (mtRsq_up-mtRsq_low)/(mtRsq_npts-1.0);
    }

  const double mtRsq_incr = par_incr;

  par_incr = 0.0;
  if (At_npts == 1)
    {
      par_incr = 0.0;
    }
  else
    {
      par_incr = (At_up-At_low)/(At_npts-1.0);
    }

  const double At_incr = par_incr;

  par_incr = 0.0;
  if (M2_npts == 1)
    {
      par_incr = 0.0;
    }
  else
    {
      par_incr = (M2_up-M2_low)/(M2_npts-1.0);
    }

  const double M2_incr = par_incr;

  // Charges that are needed to define model
  double cE6 = cos(thetaE6);
  double sE6 = sin(thetaE6);
  
  // E6 input parameters
  genericE6SSM_input_parameters input = genericE6SSM_input_parameters();

  input.QQp = QQChi*cE6+QQPsi*sE6;
  input.Qup = QuChi*cE6+QuPsi*sE6;
  input.Qdp = QdChi*cE6+QdPsi*sE6;
  input.QLp = QLChi*cE6+QLPsi*sE6;
  input.Qep = QeChi*cE6+QePsi*sE6;
  input.QSp = QSChi*cE6+QSPsi*sE6;
  input.QH1p = QH1Chi*cE6+QH1Psi*sE6;
  input.QH2p = QH2Chi*cE6+QH2Psi*sE6;
  input.QDxp = QXChi*cE6+QXPsi*sE6;
  input.QDxbarp = QXbarChi*cE6+QXbarPsi*sE6;
  input.QHpp = QHPrChi*cE6+QHPrPsi*sE6;
  input.QHpbarp = QHbarPrChi*cE6+QHbarPrPsi*sE6;

  // Scan over all points and print out Higgs masses (calculated
  // using FlexibleSUSY routines and approximate routines)
  double lambda3 = lambda3_low;
  double Alambda3 = Alambda3_low;
  double TLambda3 = lambda3*Alambda3;
  double mqL3sq = mqL3sq_low;
  double mtRsq = mtRsq_low;
  double M2 = M2_low;

  double mHdSq = 1.0e6, mHuSq = 1.0e6, mSSq = 1.0e6;
  double tol = 10.0;

  for (unsigned i_lambda3 = 1; i_lambda3 <= lambda3_npts; ++i_lambda3)
    {
      lambda3 = lambda3_low+(i_lambda3-1.0)*lambda3_incr;
  for (unsigned i_Alambda3 = 1; i_Alambda3 <= Alambda3_npts; ++i_Alambda3)
    {
      Alambda3 = Alambda3_low+(i_Alambda3-1.0)*Alambda3_incr;
  for (unsigned i_mqL3sq = 1; i_mqL3sq <= mqL3sq_npts; ++i_mqL3sq)
    {
      mqL3sq = mqL3sq_low+(i_mqL3sq-1.0)*mqL3sq_incr;
  for (unsigned i_mtRsq = 1; i_mtRsq <= mtRsq_npts; ++i_mtRsq)
    {
      mtRsq = mtRsq_low+(i_mtRsq-1.0)*mtRsq_incr;
  for (unsigned i_At = 1; i_At <= At_npts; ++i_At)
    {
      At = At_low+(i_At-1.0)*At_incr;
  for (unsigned i_M2 = 1; i_M2 <= M2_npts; ++i_M2)
    {
      M2 = M2_low+(i_M2-1.0)*M2_incr;

      // Update all variables
      mQlSq(2,2) = mqL3sq;
      mUrSq(2,2) = mtRsq;
      TLambda3 = lambda3*Alambda3;
      TYu(2,2) = Yu(2,2)*At;

      // Create a new model
      genericE6SSM_susy_parameters susy_model 
	= genericE6SSM_susy_parameters(MX, LOOPS, THRESH, input, Yd, Ye, Kappa, Lambda12, 
				       lambda3, Yu, mupr, g1, g2, g3, g1p, v1, v2, vs);

      genericE6SSM_soft_parameters soft_model
	= genericE6SSM_soft_parameters(susy_model, TYd, TYe, TKappa, TLambda12, TLambda3, 
				       TYu, Bmupr, mQlSq, mLlSq, mHdSq, mHuSq, mDrSq, 
				       mUrSq, mErSq, mSSq, mHdISq, mHuISq, mSISq, mDxSq, 
				       mDxbarSq, mHpSq, mHpbarSq, M1, M2, M3, M1p);

      genericE6SSM<algorithm_type> model = genericE6SSM<algorithm_type>(soft_model);

      // Calculate the Higgs masses. This is done using the 
      // approximate code (actually used in the scan) and
      // the FlexibleSUSY routines, for comparison
      
      // Approximate solution, note have to get M_{SUSY} using Newton iteration
      DoubleVector ewsb_guess(4);
      ewsb_guess(1) = model.get_mHd2();
      ewsb_guess(2) = model.get_mHu2();
      ewsb_guess(3) = model.get_ms2();
      ewsb_guess(4) = 2.0e3;

      bool hasEWSBProblem = ESSM_ImplementEWSBConstraints_SoftMasses(model, MX, 2.0e3, 
								     false, ewsb_guess, tol);

      model.set_mHd2(ewsb_guess(1));
      model.set_mHu2(ewsb_guess(2));
      model.set_ms2(ewsb_guess(3));
      
      model.run_to(ewsb_guess(4), PRECISION);

      DoubleVector mstop(2), mstopsq(2), mD1sq(3), mD2sq(3);
      physical_ESSM(model, mstop, mstopsq, mD1sq, mD2sq, model.get_vs(), model.get_vu()/model.get_vd());      

      // In HiggsMasses, sing = 1 indicates inaccurate Higgs mass.
      // ExpValid = 15 indicates TADPOLEPROBLEM.
      // ExpValid = 30 indicates not experimentally valid, but otherwise no problems
      // WhatCorrections selects which version to use (use 1)
      int sing = 0;
      int ExpValid = 0;
      const int WhatCorrections = 1;
      DoubleVector expBounds(2);

      expBounds(1) = HIGGSCENT - HIGGSERROR;
      expBounds(2) = HIGGSCENT + HIGGSERROR;

      DoubleVector mh(3);
      DoubleMatrix mhmix(3,3), msq(3,3);

      bool poleHiggsTachyons = HiggsMasses(model, model.get_vs(), model.get_vu()/model.get_vd(), mstop, mstopsq, WhatCorrections, 
					   false, false, expBounds, ExpValid, mh, mhmix, msq, sing);

      // FlexibleSUSY routines. Assume M_{SUSY} does not change substantially 
      // when using the exact routines for simplicity.
      model.solve_ewsb_tree_level();
      model.calculate_DRbar_parameters();
      model.calculate_MStop();
      // Check to see if M_{SUSY} has changed much
      double scale_change = Abs(ewsb_guess(4)-Sqrt(model.get_MStop()(0)*model.get_MStop()(1)))
	/(0.5*(ewsb_guess(4)+Sqrt(model.get_MStop()(0)*model.get_MStop()(1))));
      model.solve_ewsb();
      model.calculate_Mhh_pole();
      //model.calculate_MAh_pole();

      Problems<genericE6SSM_info::NUMBER_OF_PARTICLES> model_problems = model.get_problems();

      if (!model_problems.have_problem() && !model_problems.have_serious_problem() && !poleHiggsTachyons && !hasEWSBProblem)
	{
	  cout << lambda3 << " ";
	  cout << Alambda3 << " ";
	  cout << model.get_MVZ() << " ";
	  cout << model.get_MVZp() << " ";
	  cout << model.get_Mhh().minCoeff() << " ";
	  cout << mh(1) << " ";
	  cout << model.get_physical().Mhh.minCoeff() << " ";
	  cout << scale_change << endl;
	}
      else
	{
	  if (model_problems.have_problem() || model_problems.have_serious_problem())
	    {
	      cout << "# " << model_problems;
	      cout << "lambda = " << lambda3 << ", Alambda = " << Alambda3 << endl;
	    }
	  else if (hasEWSBProblem)
	    {
	      cout << "# Problems: impose approximate EWSB conditions failed" << endl;
	    }
	  else if (poleHiggsTachyons)
	    {
	      cout << "# Problems: tachyon approximate Higgs masses" << endl;
	    }
	}
    } //< M2 scan
    } //< At scan
    } //< mtRsq scan
    } //< mqL3sq scan
    } //< Alambda3 scan
    } //< lambda3 scan

  /*
    ----------------------------------------------------
   */

  return 0;
}
