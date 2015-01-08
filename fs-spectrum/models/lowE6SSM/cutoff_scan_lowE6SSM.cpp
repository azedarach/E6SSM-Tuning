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

#include "error.hpp"
#include "grid_scanner.hpp"
#include "logger.hpp"
#include "lowe.h"
#include "scan.hpp"
#include "wrappers.hpp"

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include <iostream>
#include <map>

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

lowE6SSM_input_parameters get_default_inputs(double theta, bool is_constrained)
{
   lowE6SSM_input_parameters input;

   initialize_e6_charges(theta, input);

   if (!is_constrained) {
      input.KappaInput = Eigen::Matrix<double,3,3>::Zero();
      input.KappaInput(0,0) = 0.6;
      input.KappaInput(1,1) = 0.6;
      input.KappaInput(2,2) = 0.6;
      
      input.Lambda12Input = Eigen::Matrix<double,2,2>::Zero();
      input.Lambda12Input(0,0) = 0.2;
      input.Lambda12Input(1,1) = 0.2;
      
      input.MuPrInput = 5000.; // GeV
      
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
      input.MassGInput = 2000.; // GeV
      input.MassBpInput = 300.; // GeV
   } else {
      const double Lambda0 = 0.1;
      const double Kappa0 = 0.1923;
      const double m0 = 1951.;
      const double m12 = 1003.;
      const double Azero = 500.;

      input.KappaInput = Eigen::Matrix<double,3,3>::Zero();
      input.KappaInput(0,0) = Kappa0;
      input.KappaInput(1,1) = Kappa0;
      input.KappaInput(2,2) = Kappa0;
      
      input.Lambda12Input = Eigen::Matrix<double,2,2>::Zero();
      input.Lambda12Input(0,0) = Lambda0;
      input.Lambda12Input(1,1) = Lambda0;
      
      input.MuPrInput = 569.; // GeV
      
      input.TYdInput = Eigen::Matrix<double,3,3>::Zero();
      input.TYeInput = Eigen::Matrix<double,3,3>::Zero();
      input.TKappaInput = Eigen::Matrix<double,3,3>::Zero();
      input.TLambda12Input = Eigen::Matrix<double,2,2>::Zero();
      input.TYuInput = Eigen::Matrix<double,3,3>::Zero();

      input.AYdInput = Eigen::Matrix<double,3,3>::Zero();
      input.AYdInput(0,0) = Azero;
      input.AYdInput(1,1) = Azero;
      input.AYdInput(2,2) = Azero;

      input.AYeInput = Eigen::Matrix<double,3,3>::Zero();
      input.AYeInput(0,0) = Azero;
      input.AYeInput(1,1) = Azero;
      input.AYeInput(2,2) = Azero;

      input.AYuInput = Eigen::Matrix<double,3,3>::Zero();
      input.AYuInput(0,0) = Azero;
      input.AYuInput(1,1) = Azero;
      input.AYuInput(2,2) = Azero;    

      input.BMuPrInput = 5000.; // GeV^2
      
      input.mq2Input = Eigen::Matrix<double,3,3>::Zero();
      input.mq2Input(0,0) = Sqr(m0); // GeV^2
      input.mq2Input(1,1) = Sqr(m0); // GeV^2
      
      input.ml2Input = Eigen::Matrix<double,3,3>::Zero();
      input.ml2Input(0,0) = Sqr(m0); // GeV^2
      input.ml2Input(1,1) = Sqr(m0); // GeV^2
      input.ml2Input(2,2) = Sqr(m0); // GeV^2
      
      input.md2Input = Eigen::Matrix<double,3,3>::Zero();
      input.md2Input(0,0) = Sqr(m0); // GeV^2
      input.md2Input(1,1) = Sqr(m0); // GeV^2
      input.md2Input(2,2) = Sqr(m0); // GeV^2
      
      input.mu2Input = Eigen::Matrix<double,3,3>::Zero();
      input.mu2Input(0,0) = Sqr(m0); // GeV^2
      input.mu2Input(1,1) = Sqr(m0); // GeV^2
      
      input.me2Input = Eigen::Matrix<double,3,3>::Zero();
      input.me2Input(0,0) = Sqr(m0); // GeV^2
      input.me2Input(1,1) = Sqr(m0); // GeV^2
      input.me2Input(2,2) = Sqr(m0); // GeV^2
      
      input.mH1I2Input = Eigen::Matrix<double,2,2>::Zero();
      input.mH1I2Input(0,0) = Sqr(m0); // GeV^2
      input.mH1I2Input(1,1) = Sqr(m0); // GeV^2
      
      input.mH2I2Input = Eigen::Matrix<double,2,2>::Zero();
      input.mH2I2Input(0,0) = Sqr(m0); // GeV^2
      input.mH2I2Input(1,1) = Sqr(m0); // GeV^2
      
      input.msI2Input = Eigen::Matrix<double,2,2>::Zero();
      input.msI2Input(0,0) = Sqr(m0); // GeV^2
      input.msI2Input(1,1) = Sqr(m0); // GeV^2
      
      input.mDx2Input = Eigen::Matrix<double,3,3>::Zero();
      input.mDx2Input(0,0) = Sqr(m0); // GeV^2
      input.mDx2Input(1,1) = Sqr(m0); // GeV^2
      input.mDx2Input(2,2) = Sqr(m0); // GeV^2
      
      input.mDxbar2Input = Eigen::Matrix<double,3,3>::Zero();
      input.mDxbar2Input(0,0) = Sqr(m0); // GeV^2
      input.mDxbar2Input(1,1) = Sqr(m0); // GeV^2
      input.mDxbar2Input(2,2) = Sqr(m0); // GeV^2
      
      input.mHp2Input = Sqr(m0); // GeV^2
      input.mHpbar2Input = Sqr(m0); // GeV^2
      
      input.MassBInput = m12; // GeV
      input.MassGInput = m12; // GeV
      input.MassBpInput = m12; // GeV
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

int main()
{
   typedef Two_scale algorithm_type;

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

      // initialise inputs
      lowE6SSM_input_parameters input = get_default_inputs(theta_value, use_constrained_bcs);
      input.TanBeta = TanBeta_value;
      input.LambdaxInput = Lambdax_value;
      input.TLambdaxInput = Lambdax_value * ALambdax_value;
      input.AYuInput(2,2) = AYu22_value;
      input.mq2Input(2,2) = mq222_value;
      input.mu2Input(2,2) = mu222_value;
      input.MassWBInput = MassWB_value;
      input.vsInput = vs_value;
      
      double input_scale = mx_lower;

      if (fix_at_msusy) {
         input_scale = -1.;
      }

      QedQcd oneset;
      oneset.toMz();
      // DH::note
      std::size_t point_num = 1;
      std::vector<std::size_t> position;
      while (!scan.has_finished()) {
         point_num++;
         // std::cout << "At output line number " << point_num << "\n";
         // std::cerr << "At output line number " << point_num << "\n";
         // initialise spectrum generator and do iteration
         position = scan.get_position();
         double mx_value = mx_vals(position.at(0));

         if (!fix_at_msusy) 
            input_scale = mx_value;

         lowE6SSM_spectrum_generator<algorithm_type> spectrum_generator;

         spectrum_generator.set_precision_goal(1.0e-4);
         spectrum_generator.set_input_scale(input_scale);
         spectrum_generator.set_pole_mass_loop_order(2);
         spectrum_generator.set_ewsb_loop_order(2);
         spectrum_generator.set_beta_loop_order(2);
         spectrum_generator.set_threshold_corrections_loop_order(1);

         spectrum_generator.run(oneset, input);

         const lowE6SSM<algorithm_type>& model
            = spectrum_generator.get_model();
         const lowE6SSM_physical& pole_masses = model.get_physical();
         const Problems<lowE6SSM_info::NUMBER_OF_PARTICLES>& problems
            = spectrum_generator.get_problems();

         const bool error = problems.have_serious_problem();         

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
         }

         // write out
         std::cout << " "
                   << std::setw(12) << std::left << input.TanBeta << ' '
                   << std::setw(12) << std::left << Lambdax_value << ' '
                   << std::setw(12) << std::left << ALambdax_value << ' '
                   << std::setw(12) << std::left << AYu22_value << ' '
                   << std::setw(12) << std::left << mq222_value << ' '
                   << std::setw(12) << std::left << mu222_value << ' '
                   << std::setw(12) << std::left << MassWB_value << ' '
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
