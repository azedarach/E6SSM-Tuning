/*
  Prints running of model parameters in the E6SSM
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
  double M2 = 200.0; // GeV
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
  double mqL3 = 2000.0; // GeV
  Eigen::Matrix<double,3,3> mQlSq;
  mQlSq(0,0) = Sqr(mqL1);
  mQlSq(0,1) = 0.0;
  mQlSq(0,2) = 0.0;
  mQlSq(1,0) = 0.0;
  mQlSq(1,1) = Sqr(mqL2);
  mQlSq(1,2) = 0.0;
  mQlSq(2,0) = 0.0;
  mQlSq(2,1) = 0.0;
  mQlSq(2,2) = Sqr(mqL3); 

  // Soft right-handed up-squark masses
  double muR = 5000.0; // GeV
  double mcR = 5000.0; // GeV
  double mtR = 2000.0; // GeV
  Eigen::Matrix<double,3,3> mUrSq;
  mUrSq(0,0) = Sqr(muR);
  mUrSq(0,1) = 0.0;
  mUrSq(0,2) = 0.0;
  mUrSq(1,0) = 0.0;
  mUrSq(1,1) = Sqr(mcR);
  mUrSq(1,2) = 0.0;
  mUrSq(2,0) = 0.0;
  mUrSq(2,1) = 0.0;
  mUrSq(2,2) = Sqr(mtR); 

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
  double universal_kappa = 0.6;
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

  // SUSY singlet-H1-H2 trilinear coupling
  double lambda3 = 1.1;

  // Soft singlet-H1-H2 A-terms
  double Alambda3 = 0.0; // GeV
  double TLambda3 = lambda3*Alambda3;

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

  double mHdSq = 2.35172e9;
  double mHuSq = 1.5164e8;
  double mSSq = 3.10096e8;

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

  /*
    --------------------------------------------------------------
   */

  /*
    --------------------------------------------------------------
    Scales to evaluate at
    --------------------------------------------------------------
   */

  const double q_start = MX;
  const double q_end = MZ;
  const unsigned q_npts = 5000;
  double par_incr = 0.0;
  if (q_npts != 1)
    {
      par_incr = (q_end-q_start)/(q_npts-1.0);
    }
  const double q_incr = par_incr;

  /*
    --------------------------------------------------------------
  */
  cout << "# q lambda Alambda At m_Hu^2 mHd^2 m_s^2 m_q3^2 m_u3^2 m_stop1 m_stop2 is_tachy(1=T) M1 M2 M3 M1'" << endl;
  double scale;
  outputCharacteristics(8);
  for (unsigned i = 1; i <= q_npts; ++i)
    {
      scale = q_start+(i-1.0)*q_incr;

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

      model.run_to(scale, PRECISION);

      model.calculate_MStop();

      Problems<genericE6SSM_info::NUMBER_OF_PARTICLES> model_problems = model.get_problems();

      cout << model.get_scale() << " ";
      cout << model.get_Lambdax() << " ";
      cout << model.get_TLambdax()/model.get_Lambdax() << " ";
      cout << model.get_TYu(2,2)/model.get_Yu(2,2) << " ";
      cout << model.get_mHu2() << " ";
      cout << model.get_mHd2() << " ";
      cout << model.get_ms2() << " ";
      cout << model.get_mq2(2,2) << " ";
      cout << model.get_mu2(2,2) << " ";
      cout << model.get_MStop()(0) << " ";
      cout << model.get_MStop()(1) << " ";
      if (model_problems.is_tachyon(genericE6SSM_info::Su))
	{
	  cout << "1" << " ";
	}
      else
	{
	  cout << "0" << " ";
	}
      cout << model.get_MassB() << " ";
      cout << model.get_MassWB() << " ";
      cout << model.get_MassG() << " ";
      cout << model.get_MassBp() << " ";
      cout << ESSM_EWSBCondition1(model) << " ";
      cout << ESSM_EWSBCondition2(model) << " ";
      cout << ESSM_EWSBCondition3(model) << " ";
      cout << model.get_scale()/(Sqrt(model.get_MStop()(0)*model.get_MStop()(1)))-1.0 << " ";
      cout << model.get_Kappa(2,2) << endl;
    }

  return 0;
}
