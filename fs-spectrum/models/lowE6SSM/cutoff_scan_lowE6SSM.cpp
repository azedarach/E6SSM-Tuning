// ====================================================================
// Scanning code for calculating fine tuning in the E6SSM, 
// keeping the parameters fixed and varying the cut-off (i.e.
// input) scale. Reads data from standard input and writes to standard
// output, with parameter values and scan ranges specified by flags. 
//
// Standard usage would be:
//
// >$ cat input_file.txt
// # Input file suitable for cutoff_scan_lowE6SSM.x
// MXLL=20000.0                  # lower limit of MX
// MXUL=1.0e16                   # upper limit of MX 
// MXNPTS=10                     # number of MX values to use
// TBVAL=10                      # value of \tan\beta(M_Z)
// LAMBDAXVAL=0.230769           # value of \lambda(M_{input})
// ALAMBDAXVAL=3792.69           # value of A_\lambda(M_{input})
// ATVAL=-1438.45                # value of A_t(M_{input})
// MQLSQVAL=449655.0             # value of m_Q^2(M_{input})
// MURSQVAL=586207.0             # value of m_u^2(M_{input})
// M2VAL=1050.0                  # value of M_2(M_{input})
// VSVAL=6700.0                  # value of singlet vev s(MSUSY) 
// THETAVAL=1.318116072          # U(1) mixing angle
// INPUTSCALE=MSUSY              # set M_{input} = MX or MSUSY
// INPUTSCALEBC=UNCONSTRAINED    # set remaining parameters using default input 
//                               # (UNCONSTRAINED) or to values in BM2 of 
//                               # arXiv:1302.5291 [hep-ph] (CONSTRAINED)
// # End of input file
// >$ ./cutoff_scan_lowE6SSM.x < input_file.txt
//
// ====================================================================

#include "lowE6SSM_input_parameters.hpp"
#include "lowE6SSM_spectrum_generator.hpp"
#include "lowE6SSM_two_scale_tuning_calculator.hpp"
#include "lowE6SSM_slha_io.hpp"

#include "error.hpp"
#include "grid_scanner.hpp"
#include "logger.hpp"
#include "lowe.h"
#include "scan.hpp"
#include "wrappers.hpp"

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include <chrono>
#include <iostream>
#include <map>
#include <random>

using namespace flexiblesusy;

std::string ToUpper(const std::string & s)
{
   std::string result;
   for (std::size_t index = 0; index < s.length(); ++index) {
      char a = s[index];
      a = std::toupper(a);
      result = result + a;
   }

   return result;
}

Eigen::VectorXd fill_linear_values(std::size_t num_points, double lower, double upper)
{
   Eigen::VectorXd vec;
   vec.resize(num_points);

   double incr = 0.;

   if (num_points > 1)
      incr = (upper - lower) / (num_points - 1.);

   for (std::size_t i = 0; i < num_points; ++i) {
      vec(i) = lower + incr * i;
   }

   return vec;
}

Eigen::VectorXd fill_log_values(std::size_t num_points, double lower, double upper)
{
   Eigen::VectorXd vec;
   vec.resize(num_points);

   double incr = 0.;

   if (lower <= 0. || upper <= 0.) {
      ERROR("Bounds must be positive.");
      return vec;
   }

   if (num_points > 1)
      incr = (Log(upper) - Log(lower)) / (num_points - 1.);

   for (std::size_t i = 0; i < num_points; ++i) {
      vec(i) = std::exp(Log(lower) + incr * i);
   }

   return vec;
}

void initialize_e6_charges(double theta, lowE6SSM_input_parameters& input)
{
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
}

void apply_constrained_bcs(lowE6SSM_input_parameters& input, double Lambda0, double Kappa0, double m12, double m0, double Azero)
{
      input.KappaInput(0,0) = Kappa0;
      input.KappaInput(1,1) = Kappa0;
      input.KappaInput(2,2) = Kappa0;
      
      input.Lambda12Input(0,0) = Lambda0;
      input.Lambda12Input(1,1) = Lambda0;
      
      input.AYdInput(0,0) = Azero;
      input.AYdInput(1,1) = Azero;
      input.AYdInput(2,2) = Azero;

      input.AYeInput(0,0) = Azero;
      input.AYeInput(1,1) = Azero;
      input.AYeInput(2,2) = Azero;

      input.AYuInput(0,0) = Azero;
      input.AYuInput(1,1) = Azero;
      input.AYuInput(2,2) = Azero;    

      input.mq2Input(0,0) = Sqr(m0); // GeV^2
      input.mq2Input(1,1) = Sqr(m0); // GeV^2
      
      input.ml2Input(0,0) = Sqr(m0); // GeV^2
      input.ml2Input(1,1) = Sqr(m0); // GeV^2
      input.ml2Input(2,2) = Sqr(m0); // GeV^2
      
      input.md2Input(0,0) = Sqr(m0); // GeV^2
      input.md2Input(1,1) = Sqr(m0); // GeV^2
      input.md2Input(2,2) = Sqr(m0); // GeV^2
      
      input.mu2Input(0,0) = Sqr(m0); // GeV^2
      input.mu2Input(1,1) = Sqr(m0); // GeV^2
      
      input.me2Input(0,0) = Sqr(m0); // GeV^2
      input.me2Input(1,1) = Sqr(m0); // GeV^2
      input.me2Input(2,2) = Sqr(m0); // GeV^2
      
      input.mH1I2Input(0,0) = Sqr(m0); // GeV^2
      input.mH1I2Input(1,1) = Sqr(m0); // GeV^2
      
      input.mH2I2Input(0,0) = Sqr(m0); // GeV^2
      input.mH2I2Input(1,1) = Sqr(m0); // GeV^2
      
      input.msI2Input(0,0) = Sqr(m0); // GeV^2
      input.msI2Input(1,1) = Sqr(m0); // GeV^2
      
      input.mDx2Input(0,0) = Sqr(m0); // GeV^2
      input.mDx2Input(1,1) = Sqr(m0); // GeV^2
      input.mDx2Input(2,2) = Sqr(m0); // GeV^2
      
      input.mDxbar2Input(0,0) = Sqr(m0); // GeV^2
      input.mDxbar2Input(1,1) = Sqr(m0); // GeV^2
      input.mDxbar2Input(2,2) = Sqr(m0); // GeV^2
      
      input.mHp2Input = Sqr(m0); // GeV^2
      input.mHpbar2Input = Sqr(m0); // GeV^2
      
      input.MassBInput = m12; // GeV
      input.MassGInput = m12; // GeV
      input.MassBpInput = m12; // GeV
}

enum Input_set : unsigned { cE6SSM_BM1, cE6SSM_BM2, cE6SSM_BM3, cE6SSM_BM4, 
      cE6SSM_BM5, cE6SSM_BM6, DEFAULT};

lowE6SSM_input_parameters get_default_inputs(double theta, Input_set point, bool is_constrained)
{
   lowE6SSM_input_parameters input;

   initialize_e6_charges(theta, input);

   if (!is_constrained) {

      switch (point) {
      case cE6SSM_BM1: {
         input.KappaInput = Eigen::Matrix<double,3,3>::Zero();
         input.KappaInput(0,0) = 4.97590875e-01;
         input.KappaInput(1,1) = 4.97590875e-01;
         input.KappaInput(2,2) = 4.97590875e-01;
         
         input.Lambda12Input = Eigen::Matrix<double,2,2>::Zero();
         input.Lambda12Input(0,0) = 1.38615712e-01;
         input.Lambda12Input(1,1) = 1.38615712e-01;
         
         input.MuPrInput = 8.99224583e+02; // GeV
         
         input.AYdInput = Eigen::Matrix<double,3,3>::Zero();
         input.AYdInput(0,0) = -2.490779818e+03;
         input.AYdInput(1,1) = -2.490770835e+03;
         input.AYdInput(2,2) = -2.302482528e+03;

         input.AYeInput = Eigen::Matrix<double,3,3>::Zero();
         input.AYeInput(0,0) = -6.673042246e+02;
         input.AYeInput(1,1) = -6.672906294e+02;
         input.AYeInput(2,2) = -6.632629918e+02;

         input.TKappaInput = Eigen::Matrix<double,3,3>::Zero();
         input.TKappaInput(0,0) = -7.31087487e+02;
         input.TKappaInput(1,1) = -7.31087487e+02;
         input.TKappaInput(2,2) = -7.31087487e+02;

         input.TLambda12Input = Eigen::Matrix<double,2,2>::Zero();
         input.TLambda12Input(0,0) = -9.99712778;
         input.TLambda12Input(1,1) = -9.99712778;

         input.AYuInput = Eigen::Matrix<double,3,3>::Zero();
         input.AYuInput(0,0) = -1.937449477e+03;
         input.AYuInput(1,1) = -1.937439764e+03;         

         input.BMuPrInput = -4.34941611e+05; // GeV^2
         
         input.mq2Input = Eigen::Matrix<double,3,3>::Zero();
         input.mq2Input(0,0) = 5.11244065e+06; // GeV^2
         input.mq2Input(1,1) = 5.11239866e+06; // GeV^2
         
         input.ml2Input = Eigen::Matrix<double,3,3>::Zero();
         input.ml2Input(0,0) = 4.24332080e+06; // GeV^2
         input.ml2Input(1,1) = 4.24321704e+06; // GeV^2
         input.ml2Input(2,2) = 4.21249586e+06; // GeV^2
         
         input.md2Input = Eigen::Matrix<double,3,3>::Zero();
         input.md2Input(0,0) = 5.07403459e+06; // GeV^2
         input.md2Input(1,1) = 5.07399534e+06; // GeV^2
         input.md2Input(2,2) = 4.99723927e+06; // GeV^2
         
         input.mu2Input = Eigen::Matrix<double,3,3>::Zero();
         input.mu2Input(0,0) = 5.09983374e+06; // GeV^2
         input.mu2Input(1,1) = 5.09978819e+06; // GeV^2
         
         input.me2Input = Eigen::Matrix<double,3,3>::Zero();
         input.me2Input(0,0) = 4.09753546e+06; // GeV^2
         input.me2Input(1,1) = 4.09732609e+06; // GeV^2
         input.me2Input(2,2) = 4.03532513e+06; // GeV^2
         
         input.mH1I2Input = Eigen::Matrix<double,2,2>::Zero();
         input.mH1I2Input(0,0) = 4.12843671e+06; // GeV^2
         input.mH1I2Input(1,1) = 4.12843671e+06; // GeV^2
         
         input.mH2I2Input = Eigen::Matrix<double,2,2>::Zero();
         input.mH2I2Input(0,0) = 4.09182359e+06; // GeV^2
         input.mH2I2Input(1,1) = 4.09182359e+06; // GeV^2
         
         input.msI2Input = Eigen::Matrix<double,2,2>::Zero();
         input.msI2Input(0,0) = 4.19038586e+06; // GeV^2
         input.msI2Input(1,1) = 4.19038586e+06; // GeV^2
         
         input.mDx2Input = Eigen::Matrix<double,3,3>::Zero();
         input.mDx2Input(0,0) = 4.55380360e+06; // GeV^2
         input.mDx2Input(1,1) = 4.55380360e+06; // GeV^2
         input.mDx2Input(2,2) = 4.55380360e+06; // GeV^2
         
         input.mDxbar2Input = Eigen::Matrix<double,3,3>::Zero();
         input.mDxbar2Input(0,0) = 4.51577226e+06; // GeV^2
         input.mDxbar2Input(1,1) = 4.51577226e+06; // GeV^2
         input.mDxbar2Input(2,2) = 4.51577226e+06; // GeV^2
         
         input.mHp2Input = 4.24332080e+06; // GeV^2
         input.mHpbar2Input = 4.15145499e+06; // GeV^2
         
         input.MassBInput = 1.77855744e+02; // GeV
         input.MassGInput = 7.29569906e+02; // GeV
         input.MassBpInput = 1.79328012e+02; // GeV

         break;
      }
      case cE6SSM_BM2: {
         input.KappaInput = Eigen::Matrix<double,3,3>::Zero();
         input.KappaInput(0,0) = 5.17726837e-01;
         input.KappaInput(1,1) = 5.17726837e-01;
         input.KappaInput(2,2) = 5.17726837e-01;
         
         input.Lambda12Input = Eigen::Matrix<double,2,2>::Zero();
         input.Lambda12Input(0,0) = 1.33872435e-01;
         input.Lambda12Input(1,1) = 1.33872435e-01;
         
         input.MuPrInput = 8.97916020e+02; // GeV
         
         input.AYdInput = Eigen::Matrix<double,3,3>::Zero();
         input.AYdInput(0,0) = -1.848250101e+03;
         input.AYdInput(1,1) = -1.848244014e+03;
         input.AYdInput(2,2) = -1.720674215e+03;

         input.AYeInput = Eigen::Matrix<double,3,3>::Zero();
         input.AYeInput(0,0) = -8.820525215e+01;
         input.AYeInput(1,1) = -8.820621015e+01;
         input.AYeInput(2,2) = -8.849023484e+01;

         input.TKappaInput = Eigen::Matrix<double,3,3>::Zero();
         input.TKappaInput(0,0) = -7.19207879e+02;
         input.TKappaInput(1,1) = -7.19207879e+02;
         input.TKappaInput(2,2) = 7.19207879e+02; // note sign change

         input.TLambda12Input = Eigen::Matrix<double,2,2>::Zero();
         input.TLambda12Input(0,0) = -7.08483859;
         input.TLambda12Input(1,1) = -7.08483859;

         input.AYuInput = Eigen::Matrix<double,3,3>::Zero();
         input.AYuInput(0,0) = -1.467418151e+03;
         input.AYuInput(1,1) = -1.467411566e+03;         

         input.BMuPrInput = -4.20990207e+05; // GeV^2
         
         input.mq2Input = Eigen::Matrix<double,3,3>::Zero();
         input.mq2Input(0,0) = 5.76210034E+06;//4.93021964e+06; // GeV^2
         input.mq2Input(1,1) = 5.76205345E+06;//4.93017789e+06; // GeV^2
         
         input.ml2Input = Eigen::Matrix<double,3,3>::Zero();
         input.ml2Input(0,0) = 4.93965102E+06;//3.96022681e+06; // GeV^2
         input.ml2Input(1,1) = 4.93949052E+06;//3.96008182e+06; // GeV^2
         input.ml2Input(2,2) = 4.89189756E+06;//3.91708183e+06; // GeV^2
         
         input.md2Input = Eigen::Matrix<double,3,3>::Zero();
         input.md2Input(0,0) = 5.88363593E+06;//5.04926079e+06; // GeV^2
         input.md2Input(1,1) = 5.88358594E+06;//5.04921454e+06; // GeV^2
         input.md2Input(2,2) = 5.78399108E+06;//4.95687546e+06; // GeV^2
         
         input.mu2Input = Eigen::Matrix<double,3,3>::Zero();
         input.mu2Input(0,0) = 5.53534602E+06;//4.63796804e+06; // GeV^2
         input.mu2Input(1,1) = 5.53530132E+06;//4.63792998e+06; // GeV^2
         
         input.me2Input = Eigen::Matrix<double,3,3>::Zero();
         input.me2Input(0,0) = 5.21097266E+06;//4.25677859e+06; // GeV^2
         input.me2Input(1,1) = 5.21064896E+06;//4.25648615e+06; // GeV^2
         input.me2Input(2,2) = 5.11465282E+06;//4.16975022e+06; // GeV^2
         
         input.mH1I2Input = Eigen::Matrix<double,2,2>::Zero();
         input.mH1I2Input(0,0) = 4.46478081E+06;//3.40512919e+06; // GeV^2
         input.mH1I2Input(1,1) = 4.46478081E+06;//3.40512919e+06; // GeV^2
         
         input.mH2I2Input = Eigen::Matrix<double,2,2>::Zero();
         input.mH2I2Input(0,0) = 4.80681012E+06;//3.81203970e+06; // GeV^2
         input.mH2I2Input(1,1) = 4.80681012E+06;//3.81203970e+06; // GeV^2
         
         input.msI2Input = Eigen::Matrix<double,2,2>::Zero();
         input.msI2Input(0,0) = 5.28331572E+06;//4.36110508e+06; // GeV^2
         input.msI2Input(1,1) = 5.28331572E+06;//4.36110508e+06; // GeV^2
         
         input.mDx2Input = Eigen::Matrix<double,3,3>::Zero();
         input.mDx2Input(0,0) = 4.81205009E+06;//3.98014509e+06; // GeV^2
         input.mDx2Input(1,1) = 4.81205009E+06;//3.98014509e+06; // GeV^2
         input.mDx2Input(2,2) = 4.81205009E+06;//3.98014509e+06; // GeV^2
         
         input.mDxbar2Input = Eigen::Matrix<double,3,3>::Zero();
         input.mDxbar2Input(0,0) = 4.90389439E+06;//4.08695764e+06; // GeV^2
         input.mDxbar2Input(1,1) = 4.90389439E+06;//4.08695764e+06; // GeV^2
         input.mDxbar2Input(2,2) = 4.90389439E+06;//4.08695764e+06; // GeV^2
         
         input.mHp2Input = 4.93965103E+06;//3.96022682e+06; // GeV^2
         input.mHpbar2Input = 4.86889351E+06;//3.85796087e+06; // GeV^2
         
         input.MassBInput = 1.73364059e+02; // GeV
         input.MassGInput = 1200.0;//7.11480225e+02; // GeV
         input.MassBpInput = 1.75249244e+02; // GeV

         break;
      }
      case cE6SSM_BM3: {
         input.KappaInput = Eigen::Matrix<double,3,3>::Zero();
         input.KappaInput(0,0) = 5.40456640e-01;
         input.KappaInput(1,1) = 5.40456640e-01;
         input.KappaInput(2,2) = 5.40456640e-01;
         
         input.Lambda12Input = Eigen::Matrix<double,2,2>::Zero();
         input.Lambda12Input(0,0) = 1.30058044e-01;
         input.Lambda12Input(1,1) = 1.30058044e-01;
         
         input.MuPrInput = 8.95186503e+02; // GeV
         
         input.AYdInput = Eigen::Matrix<double,3,3>::Zero();
         input.AYdInput(0,0) = -1.674591506e+03;
         input.AYdInput(1,1) = -1.674586225e+03;
         input.AYdInput(2,2) = -1.563941660e+03;

         input.AYeInput = Eigen::Matrix<double,3,3>::Zero();
         input.AYeInput(0,0) = 6.850428646e+01;
         input.AYeInput(1,1) = 6.849939791e+01;
         input.AYeInput(2,2) = 6.704902596e+01;

         input.TKappaInput = Eigen::Matrix<double,3,3>::Zero();
         input.TKappaInput(0,0) = -7.15392982e+02;
         input.TKappaInput(1,1) = -7.15392982e+02;
         input.TKappaInput(2,2) = -7.15392982e+02;

         input.TLambda12Input = Eigen::Matrix<double,2,2>::Zero();
         input.TLambda12Input(0,0) = -1.77822908;
         input.TLambda12Input(1,1) = -1.77822908;

         input.AYuInput = Eigen::Matrix<double,3,3>::Zero();
         input.AYuInput(0,0) = -1.342018824e+03;
         input.AYuInput(1,1) = -1.342013117e+03;         

         input.BMuPrInput = -4.19344670e+05; // GeV^2
         
         input.mq2Input = Eigen::Matrix<double,3,3>::Zero();
         input.mq2Input(0,0) = 5.80064516e+06; // GeV^2
         input.mq2Input(1,1) = 5.80059582e+06; // GeV^2
         
         input.ml2Input = Eigen::Matrix<double,3,3>::Zero();
         input.ml2Input(0,0) = 4.91363547e+06; // GeV^2
         input.ml2Input(1,1) = 4.91345181e+06; // GeV^2
         input.ml2Input(2,2) = 4.85897362e+06; // GeV^2
         
         input.md2Input = Eigen::Matrix<double,3,3>::Zero();
         input.md2Input(0,0) = 6.00457567e+06; // GeV^2
         input.md2Input(1,1) = 6.00452031e+06; // GeV^2
         input.md2Input(2,2) = 5.89352103e+06; // GeV^2
         
         input.mu2Input = Eigen::Matrix<double,3,3>::Zero();
         input.mu2Input(0,0) = 5.45231071e+06; // GeV^2
         input.mu2Input(1,1) = 5.45226646e+06; // GeV^2
         
         input.me2Input = Eigen::Matrix<double,3,3>::Zero();
         input.me2Input(0,0) = 5.35382589e+06; // GeV^2
         input.me2Input(1,1) = 5.35345553e+06; // GeV^2
         input.me2Input(2,2) = 5.24358231e+06; // GeV^2
         
         input.mH1I2Input = Eigen::Matrix<double,2,2>::Zero();
         input.mH1I2Input(0,0) = 4.16848415e+06; // GeV^2
         input.mH1I2Input(1,1) = 4.16848415e+06; // GeV^2
         
         input.mH2I2Input = Eigen::Matrix<double,2,2>::Zero();
         input.mH2I2Input(0,0) = 4.71691136e+06; // GeV^2
         input.mH2I2Input(1,1) = 4.71691136e+06; // GeV^2
         
         input.msI2Input = Eigen::Matrix<double,2,2>::Zero();
         input.msI2Input(0,0) = 5.50135338e+06; // GeV^2
         input.msI2Input(1,1) = 5.50135338e+06; // GeV^2
         
         input.mDx2Input = Eigen::Matrix<double,3,3>::Zero();
         input.mDx2Input(0,0) = 4.61767894e+06; // GeV^2
         input.mDx2Input(1,1) = 4.61767894e+06; // GeV^2
         input.mDx2Input(2,2) = 4.61767894e+06; // GeV^2
         
         input.mDxbar2Input = Eigen::Matrix<double,3,3>::Zero();
         input.mDxbar2Input(0,0) = 4.75445797e+06; // GeV^2
         input.mDxbar2Input(1,1) = 4.75445797e+06; // GeV^2
         input.mDxbar2Input(2,2) = 4.75445797e+06; // GeV^2
         
         input.mHp2Input = 4.91363547e+06; // GeV^2
         input.mHpbar2Input = 4.76893729e+06; // GeV^2
         
         input.MassBInput = 1.75426673e+02; // GeV
         input.MassGInput = 7.18915663e+02; // GeV
         input.MassBpInput = 1.77622859e+02; // GeV

         break;
      }
      case cE6SSM_BM4: {
         input.KappaInput = Eigen::Matrix<double,3,3>::Zero();
         input.KappaInput(0,0) = 5.59205400e-01;
         input.KappaInput(1,1) = 5.59205400e-01;
         input.KappaInput(2,2) = 5.59205400e-01;
         
         input.Lambda12Input = Eigen::Matrix<double,2,2>::Zero();
         input.Lambda12Input(0,0) = 1.26549378e-01;
         input.Lambda12Input(1,1) = 1.26549378e-01;
         
         input.MuPrInput = 8.92667263e+02; // GeV
         
         input.AYdInput = Eigen::Matrix<double,3,3>::Zero();
         input.AYdInput(0,0) = -1.534805370e+03;
         input.AYdInput(1,1) = -1.534800750e+03;
         input.AYdInput(2,2) = -1.437704187e+03;

         input.AYeInput = Eigen::Matrix<double,3,3>::Zero();
         input.AYeInput(0,0) = 1.879962498e+02;
         input.AYeInput(1,1) = 1.879884168e+02;
         input.AYeInput(2,2) = 1.856639354e+02;

         input.TKappaInput = Eigen::Matrix<double,3,3>::Zero();
         input.TKappaInput(0,0) = -7.05914987e+02;
         input.TKappaInput(1,1) = -7.05914987e+02;
         input.TKappaInput(2,2) = -7.05914987e+02;

         input.TLambda12Input = Eigen::Matrix<double,2,2>::Zero();
         input.TLambda12Input(0,0) = 2.66628129;
         input.TLambda12Input(1,1) = 2.66628129;

         input.AYuInput = Eigen::Matrix<double,3,3>::Zero();
         input.AYuInput(0,0) = -1.240931187e+03;
         input.AYuInput(1,1) = -1.240926181e+03;         

         input.BMuPrInput = -4.16556053e+05; // GeV^2
         
         input.mq2Input = Eigen::Matrix<double,3,3>::Zero();
         input.mq2Input(0,0) = 6.83149172e+06; // GeV^2
         input.mq2Input(1,1) = 6.83143358e+06; // GeV^2
         
         input.ml2Input = Eigen::Matrix<double,3,3>::Zero();
         input.ml2Input(0,0) = 6.07064679e+06; // GeV^2
         input.ml2Input(1,1) = 6.07042307e+06; // GeV^2
         input.ml2Input(2,2) = 6.00404622e+06; // GeV^2
         
         input.md2Input = Eigen::Matrix<double,3,3>::Zero();
         input.md2Input(0,0) = 7.11860259e+06; // GeV^2
         input.md2Input(1,1) = 7.11853756e+06; // GeV^2
         input.md2Input(2,2) = 6.98791267e+06; // GeV^2
         
         input.mu2Input = Eigen::Matrix<double,3,3>::Zero();
         input.mu2Input(0,0) = 6.45403247e+06; // GeV^2
         input.mu2Input(1,1) = 6.45398014e+06; // GeV^2
         
         input.me2Input = Eigen::Matrix<double,3,3>::Zero();
         input.me2Input(0,0) = 6.62614580e+06; // GeV^2
         input.me2Input(1,1) = 6.62569469e+06; // GeV^2
         input.me2Input(2,2) = 6.49183582e+06; // GeV^2
         
         input.mH1I2Input = Eigen::Matrix<double,2,2>::Zero();
         input.mH1I2Input(0,0) = 5.15387181e+06; // GeV^2
         input.mH1I2Input(1,1) = 5.15387181e+06; // GeV^2
         
         input.mH2I2Input = Eigen::Matrix<double,2,2>::Zero();
         input.mH2I2Input(0,0) = 5.81449092e+06; // GeV^2
         input.mH2I2Input(1,1) = 5.81449092e+06; // GeV^2
         
         input.msI2Input = Eigen::Matrix<double,2,2>::Zero();
         input.msI2Input(0,0) = 6.82607466e+06; // GeV^2
         input.msI2Input(1,1) = 6.82607466e+06; // GeV^2
         
         input.mDx2Input = Eigen::Matrix<double,3,3>::Zero();
         input.mDx2Input(0,0) = 5.42247824e+06; // GeV^2
         input.mDx2Input(1,1) = 5.42247824e+06; // GeV^2
         input.mDx2Input(2,2) = 5.42247824e+06; // GeV^2
         
         input.mDxbar2Input = Eigen::Matrix<double,3,3>::Zero();
         input.mDxbar2Input(0,0) = 5.57557806e+06; // GeV^2
         input.mDxbar2Input(1,1) = 5.57557806e+06; // GeV^2
         input.mDxbar2Input(2,2) = 5.57557806e+06; // GeV^2
         
         input.mHp2Input = 6.07064679e+06; // GeV^2
         input.mHpbar2Input = 5.87438461e+06; // GeV^2
         
         input.MassBInput = 1.76840930e+02; // GeV
         input.MassGInput = 7.23372153e+02; // GeV
         input.MassBpInput = 1.79306030e+02; // GeV

         break;
      }
      case cE6SSM_BM5: {
         input.KappaInput = Eigen::Matrix<double,3,3>::Zero();
         input.KappaInput(0,0) = 5.74572082e-01;
         input.KappaInput(1,1) = 5.74572082e-01;
         input.KappaInput(2,2) = 5.74572082e-01;
         
         input.Lambda12Input = Eigen::Matrix<double,2,2>::Zero();
         input.Lambda12Input(0,0) = 1.23417951e-01;
         input.Lambda12Input(1,1) = 1.23417951e-01;
         
         input.MuPrInput = 8.90331564e+02; // GeV
         
         input.AYdInput = Eigen::Matrix<double,3,3>::Zero();
         input.AYdInput(0,0) = -1.451760308e+03;
         input.AYdInput(1,1) = -1.451756075e+03;
         input.AYdInput(2,2) = -1.362740988e+03;

         input.AYeInput = Eigen::Matrix<double,3,3>::Zero();
         input.AYeInput(0,0) = 2.539645262e+02;
         input.AYeInput(1,1) = 2.539551063e+02;
         input.AYeInput(2,2) = 2.511591168e+02;

         input.TKappaInput = Eigen::Matrix<double,3,3>::Zero();
         input.TKappaInput(0,0) = -6.95099021e+02;
         input.TKappaInput(1,1) = -6.95099021e+02;
         input.TKappaInput(2,2) = -6.95099021e+02;

         input.TLambda12Input = Eigen::Matrix<double,2,2>::Zero();
         input.TLambda12Input(0,0) = 6.35816408;
         input.TLambda12Input(1,1) = 6.35816408;

         input.AYuInput = Eigen::Matrix<double,3,3>::Zero();
         input.AYuInput(0,0) = -1.181009192e+03;
         input.AYuInput(1,1) = -1.181004605e+03;         

         input.BMuPrInput = -4.14328113e+05; // GeV^2
         
         input.mq2Input = Eigen::Matrix<double,3,3>::Zero();
         input.mq2Input(0,0) = 8.02588791e+06; // GeV^2
         input.mq2Input(1,1) = 8.02581997e+06; // GeV^2
         
         input.ml2Input = Eigen::Matrix<double,3,3>::Zero();
         input.ml2Input(0,0) = 7.42587809e+06; // GeV^2
         input.ml2Input(1,1) = 7.42561524e+06; // GeV^2
         input.ml2Input(2,2) = 7.34760648e+06; // GeV^2
         
         input.md2Input = Eigen::Matrix<double,3,3>::Zero();
         input.md2Input(0,0) = 8.38666241e+06; // GeV^2
         input.md2Input(1,1) = 8.38658764e+06; // GeV^2
         input.md2Input(2,2) = 8.23637106e+06; // GeV^2
         
         input.mu2Input = Eigen::Matrix<double,3,3>::Zero();
         input.mu2Input(0,0) = 7.65468038e+06; // GeV^2
         input.mu2Input(1,1) = 7.65461801e+06; // GeV^2
         
         input.me2Input = Eigen::Matrix<double,3,3>::Zero();
         input.me2Input(0,0) = 8.05127304e+06; // GeV^2
         input.me2Input(1,1) = 8.05074304e+06; // GeV^2
         input.me2Input(2,2) = 7.89343510e+06; // GeV^2
         
         input.mH1I2Input = Eigen::Matrix<double,2,2>::Zero();
         input.mH1I2Input(0,0) = 6.37427737e+06; // GeV^2
         input.mH1I2Input(1,1) = 6.37427737e+06; // GeV^2
         
         input.mH2I2Input = Eigen::Matrix<double,2,2>::Zero();
         input.mH2I2Input(0,0) = 7.10132273e+06; // GeV^2
         input.mH2I2Input(1,1) = 7.10132273e+06; // GeV^2
         
         input.msI2Input = Eigen::Matrix<double,2,2>::Zero();
         input.msI2Input(0,0) = 8.31074684e+06; // GeV^2
         input.msI2Input(1,1) = 8.31074684e+06; // GeV^2
         
         input.mDx2Input = Eigen::Matrix<double,3,3>::Zero();
         input.mDx2Input(0,0) = 6.40670836e+06; // GeV^2
         input.mDx2Input(1,1) = 6.40670836e+06; // GeV^2
         input.mDx2Input(2,2) = 6.40670836e+06; // GeV^2
         
         input.mDxbar2Input = Eigen::Matrix<double,3,3>::Zero();
         input.mDxbar2Input(0,0) = 6.55787123e+06; // GeV^2
         input.mDxbar2Input(1,1) = 6.55787123e+06; // GeV^2
         input.mDxbar2Input(2,2) = 6.55787123e+06; // GeV^2
         
         input.mHp2Input = 7.42587810e+06; // GeV^2
         input.mHpbar2Input = 7.17094327e+06; // GeV^2
         
         input.MassBInput = 1.78329759e+02; // GeV
         input.MassGInput = 7.27946082e+02; // GeV
         input.MassBpInput = 1.81019879e+02; // GeV

         break;
      }
      case cE6SSM_BM6: {
         input.KappaInput = Eigen::Matrix<double,3,3>::Zero();
         input.KappaInput(0,0) = 5.86868800e-01;
         input.KappaInput(1,1) = 5.86868800e-01;
         input.KappaInput(2,2) = 5.86868800e-01;
         
         input.Lambda12Input = Eigen::Matrix<double,2,2>::Zero();
         input.Lambda12Input(0,0) = 1.20654518e-01;
         input.Lambda12Input(1,1) = 1.20654518e-01;
         
         input.MuPrInput = 8.88171989e+02; // GeV
         
         input.AYdInput = Eigen::Matrix<double,3,3>::Zero();
         input.AYdInput(0,0) = -1.403993187e+03;
         input.AYdInput(1,1) = -1.403989188e+03;
         input.AYdInput(2,2) = -1.319866262e+03;

         input.AYeInput = Eigen::Matrix<double,3,3>::Zero();
         input.AYeInput(0,0) = 2.945757408e+02;
         input.AYeInput(1,1) = 2.945653497e+02;
         input.AYeInput(2,2) = 2.914802792e+02;

         input.TKappaInput = Eigen::Matrix<double,3,3>::Zero();
         input.TKappaInput(0,0) = -6.87449476e+02;
         input.TKappaInput(1,1) = -6.87449476e+02;
         input.TKappaInput(2,2) = -6.87449476e+02;

         input.TLambda12Input = Eigen::Matrix<double,2,2>::Zero();
         input.TLambda12Input(0,0) = 9.25383127;
         input.TLambda12Input(1,1) = 9.25383127;

         input.AYuInput = Eigen::Matrix<double,3,3>::Zero();
         input.AYuInput(0,0) = -1.147193708e+03;
         input.AYuInput(1,1) = -1.147189376e+03;         

         input.BMuPrInput = -4.14328042e+05; // GeV^2
         
         input.mq2Input = Eigen::Matrix<double,3,3>::Zero();
         input.mq2Input(0,0) = 9.34919293e+06; // GeV^2
         input.mq2Input(1,1) = 9.34911426e+06; // GeV^2
         
         input.ml2Input = Eigen::Matrix<double,3,3>::Zero();
         input.ml2Input(0,0) = 8.91467870e+06; // GeV^2
         input.ml2Input(1,1) = 8.91437556e+06; // GeV^2
         input.ml2Input(2,2) = 8.82437868e+06; // GeV^2
         
         input.md2Input = Eigen::Matrix<double,3,3>::Zero();
         input.md2Input(0,0) = 9.78090191e+06; // GeV^2
         input.md2Input(1,1) = 9.78081684e+06; // GeV^2
         input.md2Input(2,2) = 9.61004722e+06; // GeV^2
         
         input.mu2Input = Eigen::Matrix<double,3,3>::Zero();
         input.mu2Input(0,0) = 8.99713861e+06; // GeV^2
         input.mu2Input(1,1) = 8.99706487e+06; // GeV^2
         
         input.me2Input = Eigen::Matrix<double,3,3>::Zero();
         input.me2Input(0,0) = 9.59172862e+06; // GeV^2
         input.me2Input(1,1) = 9.59111740e+06; // GeV^2
         input.me2Input(2,2) = 9.40964183e+06; // GeV^2
         
         input.mH1I2Input = Eigen::Matrix<double,2,2>::Zero();
         input.mH1I2Input(0,0) = 7.73710607e+06; // GeV^2
         input.mH1I2Input(1,1) = 7.73710607e+06; // GeV^2
         
         input.mH2I2Input = Eigen::Matrix<double,2,2>::Zero();
         input.mH2I2Input(0,0) = 8.51452594e+06; // GeV^2
         input.mH2I2Input(1,1) = 8.51452594e+06; // GeV^2
         
         input.msI2Input = Eigen::Matrix<double,2,2>::Zero();
         input.msI2Input(0,0) = 9.91724208e+06; // GeV^2
         input.msI2Input(1,1) = 9.91724208e+06; // GeV^2
         
         input.mDx2Input = Eigen::Matrix<double,3,3>::Zero();
         input.mDx2Input(0,0) = 7.51154798e+06; // GeV^2
         input.mDx2Input(1,1) = 7.51154798e+06; // GeV^2
         input.mDx2Input(2,2) = 7.51154798e+06; // GeV^2
         
         input.mDxbar2Input = Eigen::Matrix<double,3,3>::Zero();
         input.mDxbar2Input(0,0) = 7.65307614e+06; // GeV^2
         input.mDxbar2Input(1,1) = 7.65307614e+06; // GeV^2
         input.mDxbar2Input(2,2) = 7.65307614e+06; // GeV^2
         
         input.mHp2Input = 8.91467871e+06; // GeV^2
         input.mHpbar2Input = 8.59464762e+06; // GeV^2
         
         input.MassBInput = 1.80587402e+02; // GeV
         input.MassGInput = 7.35634256e+02; // GeV
         input.MassBpInput = 1.83479012e+02; // GeV
         
         break;
      }
      case DEFAULT: { 
         input.KappaInput = Eigen::Matrix<double,3,3>::Zero();
         input.KappaInput(0,0) = 0.6;
         input.KappaInput(1,1) = 0.6;
         input.KappaInput(2,2) = 0.6;
         
         input.Lambda12Input = Eigen::Matrix<double,2,2>::Zero();
         input.Lambda12Input(0,0) = 0.2;
         input.Lambda12Input(1,1) = 0.2;
         
         input.MuPrInput = 3000.; // GeV
         
         input.TYdInput = Eigen::Matrix<double,3,3>::Zero();
         input.TYeInput = Eigen::Matrix<double,3,3>::Zero();
         input.TKappaInput = Eigen::Matrix<double,3,3>::Zero();
         input.TLambda12Input = Eigen::Matrix<double,2,2>::Zero();
         input.TYuInput = Eigen::Matrix<double,3,3>::Zero();
         
         input.BMuPrInput = 5000.; // GeV^2
         
         input.mq2Input = Eigen::Matrix<double,3,3>::Zero();
         input.mq2Input(0,0) = Sqr(5000.); // GeV^2
         input.mq2Input(1,1) = Sqr(5000.); // GeV^2
         
         input.ml2Input = Eigen::Matrix<double,3,3>::Zero();
         input.ml2Input(0,0) = Sqr(5000.); // GeV^2
         input.ml2Input(1,1) = Sqr(5000.); // GeV^2
         input.ml2Input(2,2) = Sqr(5000.); // GeV^2
         
         input.md2Input = Eigen::Matrix<double,3,3>::Zero();
         input.md2Input(0,0) = Sqr(5000.); // GeV^2
         input.md2Input(1,1) = Sqr(5000.); // GeV^2
         input.md2Input(2,2) = Sqr(5000.); // GeV^2
         
         input.mu2Input = Eigen::Matrix<double,3,3>::Zero();
         input.mu2Input(0,0) = Sqr(5000.); // GeV^2
         input.mu2Input(1,1) = Sqr(5000.); // GeV^2
         
         input.me2Input = Eigen::Matrix<double,3,3>::Zero();
         input.me2Input(0,0) = Sqr(5000.); // GeV^2
         input.me2Input(1,1) = Sqr(5000.); // GeV^2
         input.me2Input(2,2) = Sqr(5000.); // GeV^2
         
         input.mH1I2Input = Eigen::Matrix<double,2,2>::Zero();
         input.mH1I2Input(0,0) = Sqr(5000.); // GeV^2
         input.mH1I2Input(1,1) = Sqr(5000.); // GeV^2
         
         input.mH2I2Input = Eigen::Matrix<double,2,2>::Zero();
         input.mH2I2Input(0,0) = Sqr(5000.); // GeV^2
         input.mH2I2Input(1,1) = Sqr(5000.); // GeV^2
         
         input.msI2Input = Eigen::Matrix<double,2,2>::Zero();
         input.msI2Input(0,0) = Sqr(5000.); // GeV^2
         input.msI2Input(1,1) = Sqr(5000.); // GeV^2
         
         input.mDx2Input = Eigen::Matrix<double,3,3>::Zero();
         input.mDx2Input(0,0) = Sqr(5000.); // GeV^2
         input.mDx2Input(1,1) = Sqr(5000.); // GeV^2
         input.mDx2Input(2,2) = Sqr(5000.); // GeV^2
         
         input.mDxbar2Input = Eigen::Matrix<double,3,3>::Zero();
         input.mDxbar2Input(0,0) = Sqr(5000.); // GeV^2
         input.mDxbar2Input(1,1) = Sqr(5000.); // GeV^2
         input.mDxbar2Input(2,2) = Sqr(5000.); // GeV^2
         
         input.mHp2Input = Sqr(5000.); // GeV^2
         input.mHpbar2Input = Sqr(5000.); // GeV^2
         
         input.MassBInput = 300.; // GeV
         input.MassGInput = 1200.0;//2000.; // GeV
         input.MassBpInput = 300.; // GeV
         break;
      }
      }
   } else {
      const double Lambda0 = 0.1;
      const double Kappa0 = 0.1923;
      const double m0 = 1951.;
      const double m12 = 1003.;
      const double Azero = 500.;

      input.MuPrInput = 569.; // GeV
      input.BMuPrInput = 5000.; // GeV^2

      apply_constrained_bcs(input, Lambda0, Kappa0, m12, m0, Azero);
   }

   return input;
}

bool calculate_3rd_generation_stop_masses(const lowE6SSM<Two_scale>& model, double & msf1, double & msf2)
{
   bool stop_tachyon = false;

   const double Lambdax = model.get_Lambdax();
   const double Yu22 = model.get_Yu(2,2);
   const double g1 = model.get_g1();
   const double g2 = model.get_g2();
   const double gN = model.get_gN();
   const double TYu22 = model.get_TYu(2,2);
   const double vd = model.get_vd();
   const double vu = model.get_vu();
   const double vs = model.get_vs();
   const double mq222 = model.get_mq2(2,2);
   const double mu222 = model.get_mu2(2,2);

   const double gbar = Sqrt(g2 * g2 + 0.6 * g1 * g1);

   const double mt = Yu22 * vu / Sqrt(2.);

   const double mix = vu * TYu22 - Yu22 * Lambdax * vs * vd / Sqrt(2.);

   const lowE6SSM_input_parameters input = model.get_input();
   const double QH1p = input.QH1p;
   const double QH2p = input.QH2p;
   const double QSp = input.QSp;
   const double QQp = input.QQp;
   const double Qup = input.Qup;

   const double deltaQ = 0.5*Sqr(gN)*(QH1p*Sqr(vd)+QH2p*Sqr(vu)+QSp*Sqr(vs))*QQp;
   const double deltaU = 0.5*Sqr(gN)*(QH1p*Sqr(vd)+QH2p*Sqr(vu)+QSp*Sqr(vs))*Qup;

   msf1 = 0.5 * (mq222 + mu222 + 0.125 * gbar * gbar * (vd * vd - vu * vu) + 2.0 * mt * mt + deltaQ + deltaU
                 - Sqrt(Sqr(mq222 - mu222 + 0.125 * (g2 * g2 - g1 * g1) * (vd * vd - vu * vu) + deltaQ - deltaU)
                        + 2.0 * Sqr(mix)));

   msf2 = 0.5 * (mq222 + mu222 + 0.125 * gbar * gbar * (vd * vd - vu * vu) + 2.0 * mt * mt + deltaQ + deltaU
                 + Sqrt(Sqr(mq222 - mu222 + 0.125 * (g2 * g2 - g1 * g1) * (vd * vd - vu * vu) + deltaQ - deltaU)
                        + 2.0 * Sqr(mix)));

   if (msf1 < 0.0 || msf2 < 0.0) {
      stop_tachyon = true;
   } 

   msf1 = AbsSqrt(msf1);
   msf2 = AbsSqrt(msf2);

   return stop_tachyon;
}

double maximum_tuning(const std::map<lowE6SSM_info::Parameters,double>& tunings)
{
   double max = -1.0e6;

   for (std::map<lowE6SSM_info::Parameters,double>::const_iterator it = tunings.begin(),
           end = tunings.end(); it != end; ++it) {
      if (it->second > max) max = it->second;
   }

   return max;
}

// partially randomize: only parameters in the main scan are set to random values
void partially_randomize_input_parameters(lowE6SSM_input_parameters& input, std::default_random_engine& gen)
{
   const double LambdaxInput = input.LambdaxInput;  
   const double TLambdaxInput = input.TLambdaxInput;
   const double AYu22Input = input.AYuInput(2,2);
   const double mq222Input = input.mq2Input(2,2);
   const double mu222Input = input.mu2Input(2,2);
   const double MassWBInput = input.MassWBInput;

   double Lambdax_width = 0.1;
   double TLambdax_width = 0.1;
   double AYu22_width = 0.1;
   double mq222_width = 0.1;
   double mu222_width = 0.1;
   double MassWB_width = 0.1;

   if (Abs(TLambdaxInput) > 1.) {
      TLambdax_width = 0.1 * Abs(TLambdaxInput);
   }

   if (Abs(AYu22Input) > 1.) {
      AYu22_width = 0.1 * Abs(AYu22Input);
   }

   if (Abs(mq222Input) > 1.) {
      mq222_width = 0.1 * Abs(mq222Input);
   }

   if (Abs(mu222Input) > 1.) {
      mu222_width = 0.1 * Abs(mu222Input);
   }

   if (Abs(MassWBInput) > 1.) {
      MassWB_width = 0.1 * Abs(MassWBInput);
   }

   std::uniform_real_distribution<double> Lambdax_distribution(LambdaxInput - Lambdax_width, LambdaxInput + Lambdax_width);
   std::uniform_real_distribution<double> TLambdax_distribution(TLambdaxInput - TLambdax_width, TLambdaxInput + TLambdax_width);
   std::uniform_real_distribution<double> AYu22_distribution(AYu22Input - AYu22_width, AYu22Input + AYu22_width);
   std::uniform_real_distribution<double> mq222_distribution(mq222Input - mq222_width, mq222Input + mq222_width);
   std::uniform_real_distribution<double> mu222_distribution(mu222Input - mu222_width, mu222Input + mu222_width);
   std::uniform_real_distribution<double> MassWB_distribution(MassWBInput - MassWB_width, MassWBInput + MassWB_width);

   input.LambdaxInput = Lambdax_distribution(gen);
   input.TLambdaxInput = TLambdax_distribution(gen);
   input.AYuInput(2,2) = AYu22_distribution(gen);
   input.mq2Input(2,2) = mq222_distribution(gen);
   input.mu2Input(2,2) = mu222_distribution(gen);
   input.MassWBInput = MassWB_distribution(gen);
}

// constrained randomize: contrained model inputs (m0, m12 etc) are set to random values, except MuPr
// and BMuPr
void constrained_randomize_input_parameters(lowE6SSM_input_parameters& input, double Lambdax, double Lambda12, double Kappa0, double m12, double m0, double Azero, std::default_random_engine& gen)
{
   double Lambdax_width = 0.1;
   double Lambda12_width = 0.1;
   double Kappa0_width = 0.1;
   double m12_width = 0.1;
   double m0_width = 0.1;
   double Azero_width = 0.1;

   if (Abs(m12) > 1.) {
      m12_width = 0.1 * Abs(m12);
   }

   if (Abs(m0) > 1.) {
      m0_width = 0.1 * Abs(m0);
   }

   if (Abs(Azero) > 1.) {
      Azero_width = 0.1 * Abs(Azero);
   }

   std::uniform_real_distribution<double> Lambdax_distribution(Lambdax - Lambdax_width, Lambdax + Lambdax_width);
   std::uniform_real_distribution<double> Lambda12_distribution(Lambda12 - Lambda12_width, Lambda12 + Lambda12_width);
   std::uniform_real_distribution<double> Kappa0_distribution(Kappa0 - Kappa0_width, Kappa0 + Kappa0_width);
   std::uniform_real_distribution<double> m0_distribution(m0 - m0_width, m0_width + m0_width);
   std::uniform_real_distribution<double> m12_distribution(m12 - m12_width, m12 + m12_width);
   std::uniform_real_distribution<double> Azero_distribution(Azero - Azero_width, Azero + Azero_width);

   double Lambdax_value = Lambdax_distribution(gen);
   double Lambda12_value = Lambda12_distribution(gen);
   double Kappa0_value = Kappa0_distribution(gen);
   double m0_value = m0_distribution(gen);
   double m12_value = m12_distribution(gen);
   double Azero_value = Azero_distribution(gen);

   apply_constrained_bcs(input, Lambda12_value, Kappa0_value, m12_value, m0_value, Azero_value);

   input.LambdaxInput = Lambdax_value;
   input.TLambdaxInput = Lambdax_value * Azero_value;
   input.AYuInput(2,2) = Azero_value;
   input.mq2Input(2,2) = Sqr(m0_value);
   input.mu2Input(2,2) = Sqr(m0_value);
   input.MassWBInput = m12_value;
}

int main()
{
   typedef Two_scale algorithm_type;

   // for "microscans"
   unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
   std::default_random_engine generator(seed);

   // default number of cut-off scales
   int init_mx_npts = 1; //< for dealing with user input
   std::size_t mx_npts = 1; //< actual value

   // default lower bound
   double mx_lower = 20000.0;

   // default upper bound
   double mx_upper = 1.0e16;

   // default values for input parameters
   double TanBeta_value = 10.0;
   double Lambdax_value = 0.1;
   double ALambdax_value = 5000.0; // GeV
   double AYu22_value = 5000.0; // GeV
   double mq222_value = 4.0e4; // GeV^2
   double mu222_value = 4.0e4; // GeV^2
   double MassWB_value = 1050.0; // GeV
   double vs_value = 6700.0; // GeV
   double theta_value = ArcTan(Sqrt(15.));

   bool has_mx_npts = false;
   bool has_mx_lower = false;
   bool has_mx_upper = false;
   bool has_TanBeta_value = false;
   bool has_Lambdax_value = false;
   bool has_ALambdax_value = false;
   bool has_AYu22_value = false;
   bool has_mq222_value = false;
   bool has_mu222_value = false;
   bool has_MassWB_value = false;
   bool has_vs_value = false;
   bool has_theta_value = false;

   bool fix_at_msusy = false;
   bool use_constrained_bcs = false;

   // additional options for the cE6SSM
   double Lambda12_value = 0.1;
   double Kappa0_value = 0.1923;
   double m12_value = 1003.;
   double m0_value = 1951.;
   double Azero_value = 500.;

   bool has_Lambda12_value = false;
   bool has_Kappa0_value = false;
   bool has_m12_value = false;
   bool has_m0_value = false;
   bool has_Azero_value = false;

   // options to do small scan about point
   bool microscanning = false;
   bool partial_randomize = true;
   bool constrained_randomize = false;

   Input_set initial_input_set = DEFAULT;

   // read in scan parameters
   try {
      std::string line;
      std::size_t line_num = 0;

      while (getline(cin, line)) {
         line_num++;

         std::istringstream input(line);
         std::string word1, word2;

         input >> word1;

         if (word1.find("#") != 0) {
            // remove anything after the first comment character
            if (word1.find("#") != string::npos) {
               word1 = word1.substr(0, word1.find("#"));
            }

            // all valid options have an = in them, so check for that
            if (word1.find("=") == string::npos) {
               WARNING("Invalid option " + word1 + ": ignoring it.");
            } else {
               // divide into two strings: option and value
               word2 = word1.substr(word1.find("=") + 1, string::npos); //< value
               word1 = word1.substr(0, word1.find("=")); //< option
               word1 = ToUpper(word1);

               std::istringstream kk(word2);

               if (word1 == "MXLL") {
                  kk >> mx_lower;
                  has_mx_lower = true;
               } else if (word1 == "MXUL") {
                  kk >> mx_upper;
                  has_mx_upper = true;
               } else if (word1 == "MXNPTS") {
                  kk >> init_mx_npts;
                  has_mx_npts = true;
               } else if (word1 == "TBVAL") {
                  kk >> TanBeta_value;
                  has_TanBeta_value = true;
               } else if (word1 == "LAMBDAXVAL") {
                  kk >> Lambdax_value;
                  has_Lambdax_value = true;
               } else if (word1 == "ALAMBDAXVAL") {
                  kk >> ALambdax_value;
                  has_ALambdax_value = true;
               } else if (word1 == "ATVAL") {
                  kk >> AYu22_value;
                  has_AYu22_value = true;
               } else if (word1 == "MQLSQVAL") {
                  kk >> mq222_value;
                  has_mq222_value = true;
               } else if (word1 == "MURSQVAL") {
                  kk >> mu222_value;
                  has_mu222_value = true;
               } else if (word1 == "M2VAL") {
                  kk >> MassWB_value;
                  has_MassWB_value = true;
               } else if (word1 == "VSVAL") {
                  kk >> vs_value;
                  has_vs_value = true;
               } else if (word1 == "THETAVAL") {
                  kk >> theta_value;
                  has_theta_value = true;
               } else if (word1 == "LAMBDA12VAL") {
                  kk >> Lambda12_value;
                  has_Lambda12_value = true;
               } else if (word1 == "KAPPA0VAL") {
                  kk >> Kappa0_value;
                  has_Kappa0_value = true;
               } else if (word1 == "M12VAL") {
                  kk >> m12_value;
                  has_m12_value = true;
               } else if (word1 == "M0VAL") {
                  kk >> m0_value;
                  has_m0_value = true;
               } else if (word1 == "A0VAL") {
                  kk >> Azero_value;
                  has_Azero_value = true;
               } else if (word1 == "INPUTSCALE") {
                  if (ToUpper(word2) == "MSUSY") {
                     fix_at_msusy = true;
                  } else if (ToUpper(word2) == "MX") {
                     fix_at_msusy = false;
                  } else {
                     WARNING("Unrecognised value '" + word2 + "' for option " + word1 + ": ignoring it");
                  }
               } else if (word1 == "INPUTSCALEBC") {
                  if (ToUpper(word2) == "CONSTRAINED") {
                     use_constrained_bcs = true;
                  } else if (ToUpper(word2) == "UNCONSTRAINED") {
                     use_constrained_bcs = false;
                  } else {
                     WARNING("Unrecognised value '" + word2 + "' for option " + word1 + ": ignoring it");
                  }
               } else if (word1 == "POINTSCAN") {
                  if (ToUpper(word2) == "TRUE") {
                     microscanning = true;
                  } else if (ToUpper(word2) == "FALSE") {
                     microscanning = false;
                  } else {
                     WARNING("Unrecognised value '" + word2 + "' for option " + word1 + ": ignoring it");
                  }
               } else if (word1 == "RANDOMIZE") {
                  if (ToUpper(word2) == "PARTIAL") {
                     partial_randomize = true;
                     constrained_randomize = false;
                  } else if (ToUpper(word2) == "CONSTRAINED") {
                     partial_randomize = false;
                     constrained_randomize = true;
                  } else {
                     WARNING("Unrecognised value '" + word2 + "' for option " + word1 + ": ignoring it");
                  }
               } else if (word1 == "INPUTSET") {
                  if (ToUpper(word2) == "CE6SSMBM1") {
                     initial_input_set = cE6SSM_BM1;
                  } else if (ToUpper(word2) == "CE6SSMBM2") {
                     initial_input_set = cE6SSM_BM2;
                  } else if (ToUpper(word2) == "CE6SSMBM3") {
                     initial_input_set = cE6SSM_BM3;
                  } else if (ToUpper(word2) == "CE6SSMBM4") {
                     initial_input_set = cE6SSM_BM4;
                  } else if (ToUpper(word2) == "CE6SSMBM5") {
                     initial_input_set = cE6SSM_BM5;
                  } else if (ToUpper(word2) == "CE6SSMBM6") {
                     initial_input_set = cE6SSM_BM6;
                  } else if (ToUpper(word2) == "DEFAULT") {
                     initial_input_set = DEFAULT;
                  } else {
                     initial_input_set = DEFAULT;
                  }
               } else {
                  WARNING("Unrecognised option '" + word1 + "': ignoring it");
               }

            }
         } // if (word1.find("#") != 0)
      } // while (getline(cin, line))

      // check inputs are valid
      if (!has_mx_lower) {
         WARNING("MX lower limit not found: using default value");
      }

      if (!has_mx_upper) {
         WARNING("MX upper limit not found: using default value");
      }

      if (!has_mx_npts) {
         WARNING("MX number of points not found: using default value");
      }

      if (!has_TanBeta_value) {
         WARNING("TanBeta value not found: using default value");
      }

      if (!has_Lambdax_value) {
         WARNING("Lambdax value not found: using default value");
      }

      if (!has_ALambdax_value) {
         WARNING("ALambdax value not found: using default value");
      }

      if (!has_AYu22_value) {
         WARNING("AYu(2,2) value not found: using default value");
      }

      if (!has_mq222_value) {
         WARNING("mq2(2,2) value not found: using default value");
      }

      if (!has_mu222_value) {
         WARNING("mu2(2,2) value not found: using default value");
      }

      if (!has_MassWB_value) {
         WARNING("MassWB value not found: using default value");
      }

      if (!has_vs_value) {
         WARNING("vs value not found: using default value");
      }

      if (!has_theta_value) {
         WARNING("Theta value not found: using default value");
      }

      if (use_constrained_bcs) {
         if (!has_Lambda12_value) {
            WARNING("Constrained BCs requested but Lambda_{1,2} value not found: using default value");
         }

         if (!has_Kappa0_value) {
            WARNING("Constrained BCs requested but Kappa_{1,2,3} value not found: using default value");
         }

         if (!has_m12_value) {
            WARNING("Constrained BCs requested but M_{1/2} value not found: using default value");
         }

         if (!has_m0_value) {
            WARNING("Constrained BCs requested but m_0 value not found: using default value");
         }

         if (!has_Azero_value) {
            WARNING("Constrained BCs requested but A_0 value not found: using default value");
         }
      }

      double temp;
      if (mx_lower > mx_upper) {
         temp = mx_lower;
         mx_lower = mx_upper;
         mx_upper = temp;
      }

      // check for valid values of TanBeta
      if (TanBeta_value < 1.0 || TanBeta_value > 500.0) {
         WARNING("Invalid TanBeta value of TanBeta: using default value");
         TanBeta_value = 10.0;
      }

      // check for valid number of points for cut-off scale
      if (init_mx_npts < 1) {
         WARNING("Invalid number of MX points: using default 1");
         mx_npts = 1;
      } else {
         mx_npts = init_mx_npts;
      }

      // get scanned cut-off scale values
      Eigen::VectorXd mx_vals = fill_log_values(mx_npts, mx_lower, mx_upper);

      std::vector<std::size_t> scan_dimensions
         = {mx_npts};

      Grid_scanner scan(scan_dimensions);

      // print header line
      std::cerr << "# "
                << std::setw(12) << std::left << "MX(GeV)" << ' '
                << std::setw(12) << std::left << "vd" << ' '
                << std::setw(12) << std::left << "vu" << ' '
                << std::setw(12) << std::left << "vs" << ' '
                << std::setw(12) << std::left << "alpha(0,0)" << ' '
                << std::setw(12) << std::left << "alpha(0,1)" << ' '
                << std::setw(12) << std::left << "alpha(0,2)" << ' '
                << std::setw(12) << std::left << "alpha(0,3)" << ' '
                << std::setw(12) << std::left << "alpha(0,4)" << ' '
                << std::setw(12) << std::left << "alpha(0,5)" << ' '
                << std::setw(12) << std::left << "alpha(0,6)" << ' '
                << std::setw(12) << std::left << "alpha(0,7)" << ' '
                << std::setw(12) << std::left << "alpha(0,8)" << ' '
                << std::setw(12) << std::left << "alpha(0,9)" << ' '
                << std::setw(12) << std::left << "alpha(0,10)" << ' '
                << std::setw(12) << std::left << "alpha(0,11)" << ' '
                << std::setw(12) << std::left << "alpha(1,0)" << ' '
                << std::setw(12) << std::left << "alpha(1,1)" << ' '
                << std::setw(12) << std::left << "alpha(1,2)" << ' '
                << std::setw(12) << std::left << "alpha(1,3)" << ' '
                << std::setw(12) << std::left << "alpha(1,4)" << ' '
                << std::setw(12) << std::left << "alpha(1,5)" << ' '
                << std::setw(12) << std::left << "alpha(1,6)" << ' '
                << std::setw(12) << std::left << "alpha(1,7)" << ' '
                << std::setw(12) << std::left << "alpha(1,8)" << ' '
                << std::setw(12) << std::left << "alpha(1,9)" << ' '
                << std::setw(12) << std::left << "alpha(1,10)" << ' '
                << std::setw(12) << std::left << "alpha(1,11)" << ' '
                << std::setw(12) << std::left << "alpha(2,0)" << ' '
                << std::setw(12) << std::left << "alpha(2,1)" << ' '
                << std::setw(12) << std::left << "alpha(2,2)" << ' '
                << std::setw(12) << std::left << "alpha(2,3)" << ' '
                << std::setw(12) << std::left << "alpha(2,4)" << ' '
                << std::setw(12) << std::left << "alpha(2,5)" << ' '
                << std::setw(12) << std::left << "alpha(2,6)" << ' '
                << std::setw(12) << std::left << "alpha(2,7)" << ' '
                << std::setw(12) << std::left << "alpha(2,8)" << ' '
                << std::setw(12) << std::left << "alpha(2,9)" << ' '
                << std::setw(12) << std::left << "alpha(2,10)" << ' '
                << std::setw(12) << std::left << "alpha(2,11)" << ' '
                << std::setw(12) << std::left << "dLambdaxdLambdax" << ' '
                << std::setw(12) << std::left << "dYu22dLambdax" << ' '
                << std::setw(12) << std::left << "dg1dLambdax" << ' '
                << std::setw(12) << std::left << "dg2dLambdax" << ' '
                << std::setw(12) << std::left << "dgNdLambdax" << ' '
                << std::setw(12) << std::left << "dTLambdaxdLambdax" << ' '
                << std::setw(12) << std::left << "dTYu22dLambdax" << ' '
                << std::setw(12) << std::left << "dmq222dLambdax" << ' '
                << std::setw(12) << std::left << "dmHd2dLambdax" << ' '
                << std::setw(12) << std::left << "dmHu2dLambdax" << ' '
                << std::setw(12) << std::left << "dmu222dLambdax" << ' '
                << std::setw(12) << std::left << "dms2dLambdax" << ' '
                << std::setw(12) << std::left << "dvddLambdax" << ' '
                << std::setw(12) << std::left << "dvudLamdbax" << ' '
                << std::setw(12) << std::left << "dvsdLamdbax" << ' '
                << std::setw(12) << std::left << "dLambdaxdTLambdax" << ' '
                << std::setw(12) << std::left << "dYu22dTLambdax" << ' '
                << std::setw(12) << std::left << "dg1dTLambdax" << ' '
                << std::setw(12) << std::left << "dg2dTLambdax" << ' '
                << std::setw(12) << std::left << "dgNdTLambdax" << ' '
                << std::setw(12) << std::left << "dTLambdaxdTLambdax" << ' '
                << std::setw(12) << std::left << "dTYu22dTLambdax" << ' '
                << std::setw(12) << std::left << "dmq222dTLambdax" << ' '
                << std::setw(12) << std::left << "dmHd2dTLambdax" << ' '
                << std::setw(12) << std::left << "dmHu2dTLambdax" << ' '
                << std::setw(12) << std::left << "dmu222dTLambdax" << ' '
                << std::setw(12) << std::left << "dms2dTLambdax" << ' '
                << std::setw(12) << std::left << "dvddTLambdax" << ' '
                << std::setw(12) << std::left << "dvudTLamdbax" << ' '
                << std::setw(12) << std::left << "dvsdTLamdbax" << ' '
                << std::setw(12) << std::left << "dLambdaxdTYu22" << ' '
                << std::setw(12) << std::left << "dYu22dTYu22" << ' '
                << std::setw(12) << std::left << "dg1dTYu22" << ' '
                << std::setw(12) << std::left << "dg2dTYu22" << ' '
                << std::setw(12) << std::left << "dgNdTYu22" << ' '
                << std::setw(12) << std::left << "dTLambdaxdTYu22" << ' '
                << std::setw(12) << std::left << "dTYu22dTYu22" << ' '
                << std::setw(12) << std::left << "dmq222dTYu22" << ' '
                << std::setw(12) << std::left << "dmHd2dTYu22" << ' '
                << std::setw(12) << std::left << "dmHu2dTYu22" << ' '
                << std::setw(12) << std::left << "dmu222dTYu22" << ' '
                << std::setw(12) << std::left << "dms2dTYu22" << ' '
                << std::setw(12) << std::left << "dvddTYu22" << ' '
                << std::setw(12) << std::left << "dvudTYu22" << ' '
                << std::setw(12) << std::left << "dvsdTYu22" << ' '
                << std::setw(12) << std::left << "dLambdaxdmq222" << ' '
                << std::setw(12) << std::left << "dYu22dmq222" << ' '
                << std::setw(12) << std::left << "dg1dmq222" << ' '
                << std::setw(12) << std::left << "dg2dmq222" << ' '
                << std::setw(12) << std::left << "dgNdmq222" << ' '
                << std::setw(12) << std::left << "dTLambdaxdmq222" << ' '
                << std::setw(12) << std::left << "dTYu22dmq222" << ' '
                << std::setw(12) << std::left << "dmq222dmq222" << ' '
                << std::setw(12) << std::left << "dmHd2dmq222" << ' '
                << std::setw(12) << std::left << "dmHu2dmq222" << ' '
                << std::setw(12) << std::left << "dmu222dmq222" << ' '
                << std::setw(12) << std::left << "dms2dmq222" << ' '
                << std::setw(12) << std::left << "dvddmq222" << ' '
                << std::setw(12) << std::left << "dvudmq222" << ' '
                << std::setw(12) << std::left << "dvsdmq222" << ' '
                << std::setw(12) << std::left << "dLambdaxdmHd2" << ' '
                << std::setw(12) << std::left << "dYu22dmHd2" << ' '
                << std::setw(12) << std::left << "dg1dmHd2" << ' '
                << std::setw(12) << std::left << "dg2dmHd2" << ' '
                << std::setw(12) << std::left << "dgNdmHd2" << ' '
                << std::setw(12) << std::left << "dTLambdaxdmHd2" << ' '
                << std::setw(12) << std::left << "dTYu22dmHd2" << ' '
                << std::setw(12) << std::left << "dmq222dmHd2" << ' '
                << std::setw(12) << std::left << "dmHd2dmHd2" << ' '
                << std::setw(12) << std::left << "dmHu2dmHd2" << ' '
                << std::setw(12) << std::left << "dmu222dmHd2" << ' '
                << std::setw(12) << std::left << "dms2dmHd2" << ' '
                << std::setw(12) << std::left << "dvddmHd2" << ' '
                << std::setw(12) << std::left << "dvudmHd2" << ' '
                << std::setw(12) << std::left << "dvsdmHd2" << ' '
                << std::setw(12) << std::left << "dLambdaxdmHu2" << ' '
                << std::setw(12) << std::left << "dYu22dmHu2" << ' '
                << std::setw(12) << std::left << "dg1dmHu2" << ' '
                << std::setw(12) << std::left << "dg2dmHu2" << ' '
                << std::setw(12) << std::left << "dgNdmHu2" << ' '
                << std::setw(12) << std::left << "dTLambdaxdmHu2" << ' '
                << std::setw(12) << std::left << "dTYu22dmHu2" << ' '
                << std::setw(12) << std::left << "dmq222dmHu2" << ' '
                << std::setw(12) << std::left << "dmHd2dmHu2" << ' '
                << std::setw(12) << std::left << "dmHu2dmHu2" << ' '
                << std::setw(12) << std::left << "dmu222dmHu2" << ' '
                << std::setw(12) << std::left << "dms2dmHu2" << ' '
                << std::setw(12) << std::left << "dvddmHu2" << ' '
                << std::setw(12) << std::left << "dvudmHu2" << ' '
                << std::setw(12) << std::left << "dvsdmHu2" << ' '
                << std::setw(12) << std::left << "dLambdaxdmu222" << ' '
                << std::setw(12) << std::left << "dYu22dmu222" << ' '
                << std::setw(12) << std::left << "dg1dmu222" << ' '
                << std::setw(12) << std::left << "dg2dmu222" << ' '
                << std::setw(12) << std::left << "dgNdmu222" << ' '
                << std::setw(12) << std::left << "dTLambdaxdmu222" << ' '
                << std::setw(12) << std::left << "dTYu22dmu222" << ' '
                << std::setw(12) << std::left << "dmq222dmu222" << ' '
                << std::setw(12) << std::left << "dmHd2dmu222" << ' '
                << std::setw(12) << std::left << "dmHu2dmu222" << ' '
                << std::setw(12) << std::left << "dmu222dmu222" << ' '
                << std::setw(12) << std::left << "dms2dmu222" << ' '
                << std::setw(12) << std::left << "dvddmu222" << ' '
                << std::setw(12) << std::left << "dvudmu222" << ' '
                << std::setw(12) << std::left << "dvsdmu222" << ' '
                << std::setw(12) << std::left << "dLambdaxdms2" << ' '
                << std::setw(12) << std::left << "dYu22dms2" << ' '
                << std::setw(12) << std::left << "dg1dms2" << ' '
                << std::setw(12) << std::left << "dg2dms2" << ' '
                << std::setw(12) << std::left << "dgNdms2" << ' '
                << std::setw(12) << std::left << "dTLambdaxdms2" << ' '
                << std::setw(12) << std::left << "dTYu22dms2" << ' '
                << std::setw(12) << std::left << "dmq222dms2" << ' '
                << std::setw(12) << std::left << "dmHd2dms2" << ' '
                << std::setw(12) << std::left << "dmHu2dms2" << ' '
                << std::setw(12) << std::left << "dmu222dms2" << ' '
                << std::setw(12) << std::left << "dms2dms2" << ' '
                << std::setw(12) << std::left << "dvddms2" << ' '
                << std::setw(12) << std::left << "dvudms2" << ' '
                << std::setw(12) << std::left << "dvsdms2" << ' '
                << std::setw(12) << std::left << "dLambdaxdMassB" << ' '
                << std::setw(12) << std::left << "dYu22dMassB" << ' '
                << std::setw(12) << std::left << "dg1dMassB" << ' '
                << std::setw(12) << std::left << "dg2dMassB" << ' '
                << std::setw(12) << std::left << "dgNdMassB" << ' '
                << std::setw(12) << std::left << "dTLambdaxdMassB" << ' '
                << std::setw(12) << std::left << "dTYu22dMassB" << ' '
                << std::setw(12) << std::left << "dmq222dMassB" << ' '
                << std::setw(12) << std::left << "dmHd2dMassB" << ' '
                << std::setw(12) << std::left << "dmHu2dMassB" << ' '
                << std::setw(12) << std::left << "dmu222dMassB" << ' '
                << std::setw(12) << std::left << "dms2dMassB" << ' '
                << std::setw(12) << std::left << "dvddMassB" << ' '
                << std::setw(12) << std::left << "dvudMassB" << ' '
                << std::setw(12) << std::left << "dvsdMassB" << ' '
                << std::setw(12) << std::left << "dLambdaxdMassWB" << ' '
                << std::setw(12) << std::left << "dYu22dMassWB" << ' '
                << std::setw(12) << std::left << "dg1dMassWB" << ' '
                << std::setw(12) << std::left << "dg2dMassWB" << ' '
                << std::setw(12) << std::left << "dgNdMassWB" << ' '
                << std::setw(12) << std::left << "dTLambdaxdMassWB" << ' '
                << std::setw(12) << std::left << "dTYu22dMassWB" << ' '
                << std::setw(12) << std::left << "dmq222dMassWB" << ' '
                << std::setw(12) << std::left << "dmHd2dMassWB" << ' '
                << std::setw(12) << std::left << "dmHu2dMassWB" << ' '
                << std::setw(12) << std::left << "dmu222dMassWB" << ' '
                << std::setw(12) << std::left << "dms2dMassWB" << ' '
                << std::setw(12) << std::left << "dvddMassWB" << ' '
                << std::setw(12) << std::left << "dvudMassWB" << ' '
                << std::setw(12) << std::left << "dvsdMassWB" << ' '
                << std::setw(12) << std::left << "dLambdaxdMassG" << ' '
                << std::setw(12) << std::left << "dYu22dMassG" << ' '
                << std::setw(12) << std::left << "dg1dMassG" << ' '
                << std::setw(12) << std::left << "dg2dMassG" << ' '
                << std::setw(12) << std::left << "dgNdMassG" << ' '
                << std::setw(12) << std::left << "dTLambdaxdMassG" << ' '
                << std::setw(12) << std::left << "dTYu22dMassG" << ' '
                << std::setw(12) << std::left << "dmq222dMassG" << ' '
                << std::setw(12) << std::left << "dmHd2dMassG" << ' '
                << std::setw(12) << std::left << "dmHu2dMassG" << ' '
                << std::setw(12) << std::left << "dmu222dMassG" << ' '
                << std::setw(12) << std::left << "dms2dMassG" << ' '
                << std::setw(12) << std::left << "dvddMassG" << ' '
                << std::setw(12) << std::left << "dvudMassG" << ' '
                << std::setw(12) << std::left << "dvsdMassG" << ' '
                << std::setw(12) << std::left << "dLambdaxdMassBp" << ' '
                << std::setw(12) << std::left << "dYu22dMassBp" << ' '
                << std::setw(12) << std::left << "dg1dMassBp" << ' '
                << std::setw(12) << std::left << "dg2dMassBp" << ' '
                << std::setw(12) << std::left << "dgNdMassBp" << ' '
                << std::setw(12) << std::left << "dTLambdaxdMassBp" << ' '
                << std::setw(12) << std::left << "dTYu22dMassBp" << ' '
                << std::setw(12) << std::left << "dmq222dMassBp" << ' '
                << std::setw(12) << std::left << "dmHd2dMassBp" << ' '
                << std::setw(12) << std::left << "dmHu2dMassBp" << ' '
                << std::setw(12) << std::left << "dmu222dMassBp" << ' '
                << std::setw(12) << std::left << "dms2dMassBp" << ' '
                << std::setw(12) << std::left << "dvddMassBp" << ' '
                << std::setw(12) << std::left << "dvudMassBp" << ' '
                << std::setw(12) << std::left << "dvsdMassBp" << ' '
                << std::setw(12) << std::left << "D(Lambdax)" << ' '
                << std::setw(12) << std::left << "D(TLambdax)" << ' '
                << std::setw(12) << std::left << "D(TYu22)" << ' '
                << std::setw(12) << std::left << "D(mq222)" << ' '
                << std::setw(12) << std::left << "D(mHd2)" << ' '
                << std::setw(12) << std::left << "D(mHu2)" << ' '
                << std::setw(12) << std::left << "D(mu222)" << ' '
                << std::setw(12) << std::left << "D(ms2)" << ' '
                << std::setw(12) << std::left << "D(MassB)" << ' '
                << std::setw(12) << std::left << "D(MassWB)" << ' '
                << std::setw(12) << std::left << "D(MassG)" << ' '
                << std::setw(12) << std::left << "D(MassBp)" << ' '
                << std::setw(12) << std::left << "MassG(MX)/GeV" << ' '
                << std::setw(12) << std::left << "mHd2(MX)/GeV^2" << ' '
                << std::setw(12) << std::left << "mHu2(MX)/GeV^2" << ' '
                << std::setw(12) << std::left << "ms2(MX)/GeV^2" << ' '
                << std::setw(12) << std::left << "mq2(2,2)(MX)/GeV^2" << ' '
                << std::setw(12) << std::left << "mu2(2,2)(MX)/GeV^2" << ' '
                << std::setw(12) << std::left << "AYu(2,2)(MX)/GeV" << ' '
                << '\n';

      std::cout << "# "
                << std::setw(12) << std::left << "TanBeta" << ' '
                << std::setw(12) << std::left << "Lambdax(MX)" << ' '
                << std::setw(12) << std::left << "ALambdax(MX)/GeV" << ' '
                << std::setw(12) << std::left << "AYu22(MX)/GeV" << ' '
                << std::setw(12) << std::left << "mq222(MX)/GeV^2" << ' '
                << std::setw(12) << std::left << "mu222(MX)/GeV^2" << ' '
                << std::setw(12) << std::left << "MassWB(MX)/GeV" << ' '
                << std::setw(12) << std::left << "Mhh(1)/GeV" << ' '
                << std::setw(12) << std::left << "Mhh(2)/GeV" << ' '
                << std::setw(12) << std::left << "Mhh(3)/GeV" << ' '
                << std::setw(12) << std::left << "mHd2/GeV^2" << ' '
                << std::setw(12) << std::left << "mHu2/GeV^2" << ' '
                << std::setw(12) << std::left << "ms2/GeV^2" << ' '
                << std::setw(12) << std::left << "g1(MS)" << ' '
                << std::setw(12) << std::left << "g2(MS)" << ' '
                << std::setw(12) << std::left << "MS/Gev" << ' '
                << std::setw(12) << std::left << "MX/GeV" << ' '
                << std::setw(12) << std::left << "s/GeV" << ' '
                << std::setw(12) << std::left << "MVZ/GeV" << ' '
                << std::setw(12) << std::left << "MVZp/GeV" << ' '
                << std::setw(12) << std::left << "MSu(1)/GeV" << ' '
                << std::setw(12) << std::left << "MSu(2)/GeV" << ' '
                << std::setw(12) << std::left << "MSb(1)/GeV" << ' '
                << std::setw(12) << std::left << "MSb(2)/GeV" << ' '
                << std::setw(12) << std::left << "MCha(1)/GeV" << ' '
                << std::setw(12) << std::left << "MCha(2)/GeV" << ' '
                << std::setw(12) << std::left << "MChi(1)/GeV" << ' '
                << std::setw(12) << std::left << "MChi(2)/GeV" << ' '
                << std::setw(12) << std::left << "MChi(3)/GeV" << ' '
                << std::setw(12) << std::left << "MChi(4)/GeV" << ' '
                << std::setw(12) << std::left << "MChi(5)/GeV" << ' '
                << std::setw(12) << std::left << "MChi(6)/GeV" << ' '
                << std::setw(12) << std::left << "Mhh(1)DRbar/GeV" << ' '
                << std::setw(12) << std::left << "Mhh(2)DRbar/GeV" << ' '
                << std::setw(12) << std::left << "Mhh(3)DRbar/GeV" << ' '
                << std::setw(12) << std::left << "MAh(1)DRbar/GeV" << ' '
                << std::setw(12) << std::left << "MAh(2)DRbar/GeV" << ' '
                << std::setw(12) << std::left << "MAh(3)DRbar/GeV" << ' '
                << std::setw(12) << std::left << "D(Max)" << ' '
                << std::setw(12) << std::left << "D(Lambdax)" << ' '
                << std::setw(12) << std::left << "D(ALambdax)" << ' '
                << std::setw(12) << std::left << "D(AYu22)" << ' '
                << std::setw(12) << std::left << "D(mq222)" << ' '
                << std::setw(12) << std::left << "D(mHd2)" << ' '
                << std::setw(12) << std::left << "D(mHu2)" << ' '
                << std::setw(12) << std::left << "D(mu222)" << ' '
                << std::setw(12) << std::left << "D(ms2)" << ' '
                << std::setw(12) << std::left << "D(MassB)" << ' '
                << std::setw(12) << std::left << "D(MassWB)" << ' '
                << std::setw(12) << std::left << "D(MassG)" << ' '
                << std::setw(12) << std::left << "D(MassBp)" << ' '
                << std::setw(12) << std::left << "tuning_err" << ' '
                << std::setw(12) << std::left << "error"
                << '\n';
      
      double input_scale = mx_lower;

      if (fix_at_msusy) {
         input_scale = -1.;
      }


      QedQcd oneset;
      oneset.toMz();
      // DH::note
      std::size_t point_num = 0;
      std::vector<std::size_t> position;
      while (!scan.has_finished()) {
         point_num++;
         // std::cout << "At output line number " << point_num << "\n";
         // std::cerr << "At output line number " << point_num << "\n";
         // initialise spectrum generator and do iteration
         position = scan.get_position();
         double mx_value = mx_vals(position.at(0));

         // initialise inputs
         lowE6SSM_input_parameters input = get_default_inputs(theta_value, initial_input_set, use_constrained_bcs);
         input.TanBeta = TanBeta_value;
         input.LambdaxInput = Lambdax_value;
         input.TLambdaxInput = Lambdax_value * ALambdax_value;
         input.AYuInput(2,2) = AYu22_value;
         input.mq2Input(2,2) = mq222_value;
         input.mu2Input(2,2) = mu222_value;
         input.MassWBInput = MassWB_value;
         input.vsInput = vs_value;
         
         if (use_constrained_bcs) {
            apply_constrained_bcs(input, Lambda12_value, Kappa0_value, m12_value, m0_value, Azero_value);
         }

         if (microscanning) {
            mx_value = mx_vals(mx_npts - 1);
            // we want to make sure the central point
            // is always included for reference
            if (point_num > 1) {
               if (partial_randomize) {
                  partially_randomize_input_parameters(input, generator);
               } else if (constrained_randomize) {
                  constrained_randomize_input_parameters(input, Lambdax_value, Lambda12_value, Kappa0_value, m12_value, m0_value, 
                                                         Azero_value, generator);
               }
            }
         }

         if (!fix_at_msusy) 
            input_scale = mx_value;

         lowE6SSM_spectrum_generator<algorithm_type> spectrum_generator;

         spectrum_generator.set_precision_goal(1.0e-4);
         spectrum_generator.set_input_scale(input_scale);
         spectrum_generator.set_pole_mass_loop_order(2);
         spectrum_generator.set_ewsb_loop_order(2);
         spectrum_generator.set_beta_loop_order(2);
         spectrum_generator.set_threshold_corrections_loop_order(1);
         //spectrum_generator.set_parameter_output_scale(1.0e16);
         spectrum_generator.run(oneset, input);

         const lowE6SSM<algorithm_type>& model
            = spectrum_generator.get_model();
         const lowE6SSM_physical& pole_masses = model.get_physical();
         const Problems<lowE6SSM_info::NUMBER_OF_PARTICLES>& problems
            = spectrum_generator.get_problems();

         const bool error = problems.have_serious_problem();         
         //std::cout << "MGlu = " << pole_masses.MGlu << "\n";
         // if no problems, calculate tuning (numerically)
         std::map<lowE6SSM_info::Parameters,double> fine_tunings;
         double max_tuning = 0.;
         bool tuning_problem = false;
         if (!error) {

            // solve ewsb
            lowE6SSM_ew_derivs ew_derivs(model);

            double mHd2 = ew_derivs.get_model().get_mHd2();
            double mHu2 = ew_derivs.get_model().get_mHu2();
            double ms2 = ew_derivs.get_model().get_ms2();         

            lowE6SSM<Two_scale> tuning_model(model);
            tuning_model.set_mHd2(mHd2);
            tuning_model.set_mHu2(mHu2);
            tuning_model.set_ms2(ms2);

            lowE6SSM_tuning_calculator tuning_calc(tuning_model);
            tuning_calc.set_tuning_scale(tuning_model.get_scale());
            tuning_calc.set_input_scale(mx_value);
            tuning_calc.set_tuning_ewsb_loop_order(1);
            tuning_calc.set_tuning_beta_loop_order(2);
            try {
               tuning_problem = tuning_calc.calculate_fine_tunings_using_rge_derivs();
            } catch (const std::string & a) {
               std::cerr << "WARNING: serious numerical problem encountered in fine tuning calculation.\n";
               tuning_problem = true;
            } catch (const char* a) {
               std::cerr << "WARNING: serious numerical problem encountered in fine tuning calculation.\n";
               tuning_problem = true;
            } catch (...) {
               std::cerr << "WARNING: serious numerical problem encountered in fine tuning calculation.\n";
               tuning_problem = true;
            }
            fine_tunings = tuning_calc.get_fine_tunings();
            max_tuning = maximum_tuning(fine_tunings);

            //tuning_model.run_to(spectrum_generator.get_high_scale());
            //tuning_model.print(std::cout);

            // tuning_model.run_to(mx_value);
            // std::cout << "In tuning, Q = " << tuning_model.get_scale() << "\n";
            // std::cout << "In tuning, Sigma1 = " << tuning_model.get_mq2().trace() - 2.0 * tuning_model.get_mu2().trace() + tuning_model.get_md2().trace() + tuning_model.get_me2().trace() - tuning_model.get_ml2().trace() + tuning_model.get_mHu2() - tuning_model.get_mHu2() + tuning_model.get_mH2I2().trace() - tuning_model.get_mH1I2().trace() + tuning_model.get_mDxbar2().trace() - tuning_model.get_mDx2().trace() - tuning_model.get_mHp2() + tuning_model.get_mHpbar2() << "\n";

            // std::cout << "In tuning, Sigma1' = " << 6.0 * tuning_model.get_mq2().trace() + 3.0 * tuning_model.get_mu2().trace() + 6.0 * tuning_model.get_md2().trace() + tuning_model.get_me2().trace() + 4.0 * tuning_model.get_ml2().trace() - 4.0 * tuning_model.get_mHu2() - 4.0 * tuning_model.get_mH2I2().trace() - 6.0 * tuning_model.get_mHd2() - 6.0 * tuning_model.get_mH1I2().trace() + 5.0 * tuning_model.get_ms2() + 5.0 * tuning_model.get_msI2().trace() - 9.0 * tuning_model.get_mDxbar2().trace() - 6.0 * tuning_model.get_mDx2().trace() + 4.0 * tuning_model.get_mHp2() - 4.0 * tuning_model.get_mHpbar2() << "\n";
         }

         // write out
         //tuning_model.run_to(spectrum_generator.get_high_scale());
         //tuning_model.print(std::cout);
         // lowE6SSM_slha_io slha_io;
         // slha_io.set_spinfo(problems);
         // slha_io.set_sminputs(oneset);
         // slha_io.set_minpar(input);
         // slha_io.set_extpar(input);
         // if (!problems.have_serious_problem())
         //    slha_io.set_spectrum(model);
         // slha_io.write_to_stream(std::cout);
         std::cout << " "
                   << std::setw(12) << std::left << input.TanBeta << ' '
                   << std::setw(12) << std::left << input.LambdaxInput << ' '
                   << std::setw(12) << std::left << input.TLambdaxInput / input.LambdaxInput << ' '
                   << std::setw(12) << std::left << input.AYuInput(2,2) << ' '
                   << std::setw(12) << std::left << input.mq2Input(2,2) << ' '
                   << std::setw(12) << std::left << input.mu2Input(2,2) << ' '
                   << std::setw(12) << std::left << input.MassWBInput << ' '
                   << std::setw(12) << std::left << pole_masses.Mhh(0) << ' '
                   << std::setw(12) << std::left << pole_masses.Mhh(1) << ' '
                   << std::setw(12) << std::left << pole_masses.Mhh(2) << ' '
                   << std::setw(12) << std::left << model.get_mHd2() << ' '
                   << std::setw(12) << std::left << model.get_mHu2() << ' '
                   << std::setw(12) << std::left << model.get_ms2() << ' '
                   << std::setw(12) << std::left << model.get_g1() << ' '
                   << std::setw(12) << std::left << model.get_gN() << ' '
                   << std::setw(12) << std::left << Sqrt(model.get_MSu()(0) * model.get_MSu()(1)) << ' '
                   << std::setw(12) << std::left << mx_value << ' '
                   << std::setw(12) << std::left << input.vsInput << ' '
                   << std::setw(12) << std::left << model.get_MVZ() << ' '
                   << std::setw(12) << std::left << model.get_MVZp() << ' '
                   << std::setw(12) << std::left << model.get_MSu()(0) << ' '
                   << std::setw(12) << std::left << model.get_MSu()(1) << ' '
                   << std::setw(12) << std::left << model.get_MSd()(0) << ' '
                   << std::setw(12) << std::left << model.get_MSd()(1) << ' '
                   << std::setw(12) << std::left << model.get_MCha()(0) << ' '
                   << std::setw(12) << std::left << model.get_MCha()(1) << ' '
                   << std::setw(12) << std::left << model.get_MChi()(0)<< ' '
                   << std::setw(12) << std::left << model.get_MChi()(1)<< ' '
                   << std::setw(12) << std::left << model.get_MChi()(2)<< ' '
                   << std::setw(12) << std::left << model.get_MChi()(3)<< ' '
                   << std::setw(12) << std::left << model.get_MChi()(4)<< ' '
                   << std::setw(12) << std::left << model.get_MChi()(5)<< ' '
                   << std::setw(12) << std::left << model.get_Mhh()(0) << ' '
                   << std::setw(12) << std::left << model.get_Mhh()(1) << ' '
                   << std::setw(12) << std::left << model.get_Mhh()(2) << ' '
                   << std::setw(12) << std::left << model.get_MAh()(0) << ' '
                   << std::setw(12) << std::left << model.get_MAh()(1) << ' '
                   << std::setw(12) << std::left << model.get_MAh()(2) << ' '
                   << std::setw(12) << std::left << max_tuning << ' '
                   << std::setw(12) << std::left << fine_tunings[lowE6SSM_info::Lambdax] << ' '
                   << std::setw(12) << std::left << fine_tunings[lowE6SSM_info::TLambdax] << ' '
                   << std::setw(12) << std::left << fine_tunings[lowE6SSM_info::TYu22] << ' '
                   << std::setw(12) << std::left << fine_tunings[lowE6SSM_info::mq222] << ' '
                   << std::setw(12) << std::left << fine_tunings[lowE6SSM_info::mHd2] << ' '
                   << std::setw(12) << std::left << fine_tunings[lowE6SSM_info::mHu2] << ' '
                   << std::setw(12) << std::left << fine_tunings[lowE6SSM_info::mu222] << ' '
                   << std::setw(12) << std::left << fine_tunings[lowE6SSM_info::ms2] << ' '
                   << std::setw(12) << std::left << fine_tunings[lowE6SSM_info::MassB] << ' '
                   << std::setw(12) << std::left << fine_tunings[lowE6SSM_info::MassWB] << ' '
                   << std::setw(12) << std::left << fine_tunings[lowE6SSM_info::MassG] << ' '
                   << std::setw(12) << std::left << fine_tunings[lowE6SSM_info::MassBp] << ' '
                   << std::setw(12) << std::left << tuning_problem << ' '
                   << std::setw(12) << std::left << error << ' ';
         if (error || tuning_problem) {
            std::cout << "# " << problems;
            if (tuning_problem) {
               std::cout << ", tuning error\n"; 
            } else {
               std::cout << '\n';
            }
         } else {
            std::cout << '\n';
         }
         
         scan.step_forward();
      } // while (!scan.has_finished())

   }
   catch(const std::string & a) { std::cout << a; return -1; }
   catch(const char * a) { cout << a; return -1; }
   catch (const Error& error) {
      ERROR(error.what());
      return EXIT_FAILURE;
   }
   catch(...) { std::cout << "Unknown type of exception caught.\n"; return -1; }

   return 0;
}
