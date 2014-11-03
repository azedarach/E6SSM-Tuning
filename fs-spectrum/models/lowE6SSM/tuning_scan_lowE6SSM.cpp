// ====================================================================
// Scanning code for calculating fine tuning in the E6SSM
// ====================================================================

#include "lowE6SSM_two_scale_tuning_calculator.hpp"

#include "error.hpp"
#include "grid_scanner.hpp"
#include "logger.hpp"
#include "wrappers.hpp"

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <sys/time.h>

using namespace flexiblesusy;
using namespace softsusy;

std::string ToUpper(const std::string & s) 
{
   std::string result;
   unsigned int index;
   for (index = 0; index < s.length(); ++index) {
      char a = s[index];
      a = std::toupper(a);
      result = result + a;
   }
   
   return result;
}

double get_wall_time()
{
   struct timeval time;
   if (gettimeofday(&time,NULL)) {
      return 0;
   }
   return (double)time.tv_sec + (double)time.tv_usec*0.000001;
}

double get_cpu_time()
{
   return (double)clock() / CLOCKS_PER_SEC;
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

Eigen::VectorXd fill_symmetric_log_values(std::size_t num_points, double lower, double upper)
{
   Eigen::VectorXd vec;
   vec.resize(2 * num_points);

   double incr = 0.;

   if (num_points > 1)
      incr = (upper - lower) / (num_points - 1.);

   for (std::size_t i = 0; i < num_points; ++i) {
      vec(num_points + i) = std::exp(lower + incr * i);
      vec(num_points - i - 1) = -vec(num_points + i);
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

void initialize_model_from_input(lowE6SSM_input_parameters input, lowE6SSM<Two_scale>& model)
{
   model.set_input_parameters(input);

   model.set_Kappa(input.KappaInput);
   model.set_Lambda12(input.Lambda12Input);
   model.set_Lambdax(input.LambdaxInput);

   model.set_MuPr(input.MuPrInput);

   model.set_gN(input.gNInput);
   model.set_vs(input.vsInput);

   model.set_TYd(input.TYdInput);
   model.set_TYe(input.TYeInput);
   model.set_TKappa(input.TKappaInput);
   model.set_TLambda12(input.TLambda12Input);
   model.set_TLambdax(input.TLambdaxInput);
   model.set_TYu(input.TYuInput);
   model.set_BMuPr(input.BMuPrInput);

   model.set_mq2(input.mq2Input);
   model.set_ml2(input.ml2Input);
   model.set_md2(input.md2Input);
   model.set_mu2(input.mu2Input);
   model.set_me2(input.me2Input);

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

lowE6SSM_input_parameters get_default_inputs()
{
   lowE6SSM_input_parameters input;

   const double theta = ArcTan(Sqrt(15.));

   initialize_e6_charges(theta, input);

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

bool calculate_3rd_generation_sbottom_masses(const lowE6SSM<Two_scale>& model, double & msf1, double & msf2)
{
   bool sbottom_tachyon = false;

   const double Lambdax = model.get_Lambdax();
   const double Yd22 = model.get_Yd(2,2);
   const double g1 = model.get_g1();
   const double g2 = model.get_g2();
   const double gN = model.get_gN();
   const double TYd22 = model.get_TYd(2,2);
   const double vd = model.get_vd();
   const double vu = model.get_vu();
   const double vs = model.get_vs();
   const double mq222 = model.get_mq2(2,2);
   const double md222 = model.get_md2(2,2);

   const double gbar = Sqrt(g2 * g2 + 0.6 * g1 * g1);

   const double mb = Yd22 * vd / Sqrt(2.);

   const double mix = vd * TYd22 - Yd22 * Lambdax * vs * vu / Sqrt(2.);

   const lowE6SSM_input_parameters input = model.get_input();
   const double QH1p = input.QH1p;
   const double QH2p = input.QH2p;
   const double QSp = input.QSp;
   const double QQp = input.QQp;
   const double Qdp = input.Qdp;

   const double deltaQ = 0.5*Sqr(gN)*(QH1p*Sqr(vd)+QH2p*Sqr(vu)+QSp*Sqr(vs))*QQp;
   const double deltaD = 0.5*Sqr(gN)*(QH1p*Sqr(vd)+QH2p*Sqr(vu)+QSp*Sqr(vs))*Qdp;

   msf1 = 0.5 * (mq222 + md222 - 0.125 * gbar * gbar * (vd * vd - vu * vu) + 2.0 * mb * mb + deltaQ + deltaD
                 - Sqrt(Sqr(mq222 - md222 - 0.125 * (g2 * g2 - 0.2 * g1 * g1) * (vd * vd - vu * vu) + deltaQ - deltaD)
                        + 2.0 * Sqr(mix)));

   msf2 = 0.5 * (mq222 + md222 - 0.125 * gbar * gbar * (vd * vd - vu * vu) + 2.0 * mb * mb + deltaQ + deltaD
                 + Sqrt(Sqr(mq222 - md222 - 0.125 * (g2 * g2 - 0.2 * g1 * g1) * (vd * vd - vu * vu) + deltaQ - deltaD)
                        + 2.0 * Sqr(mix)));

   if (msf1 < 0.0 || msf2 < 0.0) {
      sbottom_tachyon = true;
   } 

   msf1 = AbsSqrt(msf1);
   msf2 = AbsSqrt(msf2);

   return sbottom_tachyon;
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
   // default number of points in each direction
   std::size_t TanBeta_npts = 1;
   std::size_t Lambdax_npts = 1;
   std::size_t ALambdax_npts = 1;
   std::size_t AYu22_npts = 1;
   std::size_t mq222_npts = 1;
   std::size_t mu222_npts = 1;
   std::size_t MassWB_npts = 1;

   // default lower bounds
   double TanBeta_lower = 2.;
   double Lambdax_lower = 0.;
   double ALambdax_lower = 10.; // GeV
   double AYu22_lower = 10.; // GeV
   double mq222_lower = Sqr(200.); // GeV^2
   double mu222_lower = Sqr(200.); // GeV^2
   double MassWB_lower = 500.; // GeV

   // default upper bounds
   double TanBeta_upper = 50.;
   double Lambdax_upper = 3.;
   double ALambdax_upper = 10000.; // GeV
   double AYu22_upper = 10000.; // GeV
   double mq222_upper = Sqr(2000.); // GeV^2
   double mu222_upper = Sqr(2000.); // GeV^2
   double MassWB_upper = 2000.; // GeV

   // flags for using log scans
   bool use_log_scan_AYu22 = false;
   bool use_log_scan_ALambdax = false;

   // default values of Yukawa and gauge couplings, and VEVs
   Eigen::Matrix<double,3,3> YuInput = Eigen::Matrix<double,3,3>::Zero();
   Eigen::Matrix<double,3,3> YdInput = Eigen::Matrix<double,3,3>::Zero();
   Eigen::Matrix<double,3,3> YeInput = Eigen::Matrix<double,3,3>::Zero();

   double g1Input = 0.46;
   double g2Input = 0.64;
   double g3Input = 1.02;
   double gNInput = 0.45;

   double vInput = 246.; // GeV
   double vsInput = 6700.; //GeV

   double mxInput = 20000.0; // GeV

   double thetaInput = ArcTan(Sqrt(15.));

   // flags for error checking
   bool has_log_AYu22_option = false;
   bool has_non_log_AYu22_option = false;
   bool has_log_ALambdax_option = false;
   bool has_non_log_ALambdax_option = false;

   bool has_TanBeta_lower = false;
   bool has_Lambdax_lower = false;
   bool has_ALambdax_lower = false;
   bool has_AYu22_lower = false;
   bool has_mq222_lower = false;
   bool has_mu222_lower = false;
   bool has_MassWB_lower = false;

   bool has_TanBeta_upper = false;
   bool has_Lambdax_upper = false;
   bool has_ALambdax_upper = false;
   bool has_AYu22_upper = false;
   bool has_mq222_upper = false;
   bool has_mu222_upper = false;
   bool has_MassWB_upper = false;

   bool has_TanBeta_npts = false;
   bool has_Lambdax_npts = false;
   bool has_ALambdax_npts = false;
   bool has_AYu22_npts = false;
   bool has_mq222_npts = false;
   bool has_mu222_npts = false;
   bool has_MassWB_npts = false;

   bool has_YuInput00 = false;
   bool has_YuInput11 = false;
   bool has_YuInput22 = false;
   bool has_YdInput00 = false;
   bool has_YdInput11 = false;
   bool has_YdInput22 = false;
   bool has_YeInput00 = false;
   bool has_YeInput11 = false;
   bool has_YeInput22 = false;

   bool has_g1Input = false;
   bool has_g2Input = false;
   bool has_g3Input = false;
   bool has_gNInput = false;

   bool has_vInput = false;
   bool has_vsInput = false;
   bool has_mxInput = false;
   bool has_thetaInput = false;

   // read in scan parameters
   try {

      std::string line; 
      int lineNum = 0;

      while (getline(cin,line)) {
         lineNum++;

         std::istringstream input(line); 
         std::string word1, word2;
         input >> word1; 

         if (word1.find("#") != 0) { 
            // Remove anything after the first comment character
            if (word1.find("#") != string::npos) {
               word1 = word1.substr(0, word1.find("#"));
            }

            // All valid options have an = in them, so check for that
            if (word1.find("=") == string::npos) {
               std::cerr << "WARNING: invalid option '" << word1 << "' at line " << lineNum;
               std::cerr << ": ignoring it." << std::endl;
            } else {
               // Divide into two strings: option and value
               word2 = word1.substr(word1.find("=")+1, string::npos); //< value
               word1 = word1.substr(0, word1.find("=")); //< option
               word1 = ToUpper(word1);

               std::istringstream kk(word2);

               // Check against list of valid options
               if (word1 == "TBLL") {
                  kk >> TanBeta_lower;
                  has_TanBeta_lower = true;
               } else if (word1 == "TBUL") {
                  kk >> TanBeta_upper;
                  has_TanBeta_upper = true;
               } else if (word1 == "TBNPTS") {
                  kk >> TanBeta_npts;
                  has_TanBeta_npts = true;
               } else if (word1 == "LAMBDALL") {
                  kk >> Lambdax_lower;
                  has_Lambdax_lower = true;
               } else if (word1 == "LAMBDAUL") {
                  kk >> Lambdax_upper;
                  has_Lambdax_upper = true;
               } else if (word1 == "LAMBDANPTS") {
                  kk >> Lambdax_npts;
                  has_Lambdax_npts = true;
               } else if (word1 == "MQLSQLL") {
                  kk >> mq222_lower;
                  has_mq222_lower = true;
               } else if (word1 == "MQLSQUL") {
                  kk >> mq222_upper;
                  has_mq222_upper = true;
               } else if (word1 == "MQLSQNPTS") {
                  kk >> mq222_npts;
                  has_mq222_npts = true;
               } else if (word1 == "MURSQLL") {
                  kk >> mu222_lower;
                  has_mu222_lower = true;
               } else if (word1 == "MURSQUL") {
                  kk >> mu222_upper;
                  has_mu222_upper = true;
               } else if (word1 == "MURSQNPTS") {
                  kk >> mu222_npts;
                  has_mu222_npts = true;
               } else if (word1 == "ATLL") {
                  if (has_log_AYu22_option) {
                     std::cerr << "WARNING: log scan of AYu(2,2) already requested:";
                     std::cerr << " ignoring option '" << word1 << "' at line " << lineNum << "." << std::endl;
                  } else {
                     has_non_log_AYu22_option = true;
                     kk >> AYu22_lower;
                     has_AYu22_lower = true;
                  }
               } else if (word1 == "ATUL") {
                  if (has_log_AYu22_option) {
                     std::cerr << "WARNING: log scan of AYu(2,2) already requested:";
                     std::cerr << " ignoring option '" << word1 << "' at line " << lineNum << "." << std::endl;
                  } else {
                     has_non_log_AYu22_option = true;
                     kk >> AYu22_upper;
                     has_AYu22_upper = true;
                  }
               } else if (word1 == "ATNPTS") {
                  if (has_log_AYu22_option) {
                     std::cerr << "WARNING: log scan of AYu(2,2) already requested:";
                     std::cerr << " ignoring option '" << word1 << "' at line " << lineNum << "." << std::endl;
                  } else {
                     has_non_log_AYu22_option = true;
                     kk >> AYu22_npts;
                     has_AYu22_npts = true;
                  }
               } else if (word1 == "LOGATLL") {
                  if (has_non_log_AYu22_option) {
                     std::cerr << "WARNING: linear scan of AYu(2,2) already requested:";
                     std::cerr << " ignoring option '" << word1 << "' at line " << lineNum << "." << std::endl;
                  } else {
                     has_log_AYu22_option = true;
                     kk >> AYu22_lower;
                     has_AYu22_lower = true;
                  }
               } else if (word1 == "LOGATUL") {
                  if (has_non_log_AYu22_option) {
                     std::cerr << "WARNING: linear scan of AYu(2,2) already requested:";
                     std::cerr << " ignoring option '" << word1 << "' at line " << lineNum << "." << std::endl;
                  } else {
                     has_log_AYu22_option = true;
                     kk >> AYu22_upper;
                     has_AYu22_upper = true;
                  }
               } else if (word1 == "LOGATNPTS") {
                  if (has_non_log_AYu22_option) {
                     std::cerr << "WARNING: linear scan of AYu(2,2) already requested:";
                     std::cerr << " ignoring option '" << word1 << "' at line " << lineNum << "." << std::endl;
                  } else {
                     has_log_AYu22_option = true;
                     kk >> AYu22_npts;
                     has_AYu22_npts = true;
                  }
               } else if (word1 == "ALAMBDALL") {
                  if (has_log_ALambdax_option) {
                     std::cerr << "WARNING: log scan of ALambdax already requested:";
                     std::cerr << " ignoring option '" << word1 << "' at line " << lineNum << "." << std::endl;
                  } else {
                     has_non_log_ALambdax_option = true;
                     kk >> ALambdax_lower;
                     has_ALambdax_lower = true;
                  }
               } else if (word1 == "ALAMBDAUL") {
                  if (has_log_ALambdax_option) {
                     std::cerr << "WARNING: log scan of ALambdax already requested:";
                     std::cerr << " ignoring option '" << word1 << "' at line " << lineNum << "." << std::endl;
                  } else {
                     has_non_log_ALambdax_option = true;
                     kk >> ALambdax_upper;
                     has_ALambdax_upper = true;
                  }
               } else if (word1 == "ALAMBDANPTS") {
                  if (has_log_ALambdax_option) {
                     std::cerr << "WARNING: log scan of ALambdax already requested:";
                     std::cerr << " ignoring option '" << word1 << "' at line " << lineNum << "." << std::endl;
                  } else {
                     has_non_log_ALambdax_option = true;
                     kk >> ALambdax_npts;
                     has_ALambdax_npts = true;
                  }
               } else if (word1 == "LOGALAMBDALL") {
                  if (has_non_log_ALambdax_option) {
                     std::cerr << "WARNING: linear scan of ALambdax already requested:";
                     std::cerr << " ignoring option '" << word1 << "' at line " << lineNum << "." << std::endl;
                  } else {
                     has_log_ALambdax_option = true;
                     kk >> ALambdax_lower;
                     has_ALambdax_lower = true;
                  }
               } else if (word1 == "LOGALAMBDAUL") {
                  if (has_non_log_ALambdax_option) {
                     std::cerr << "WARNING: linear scan of ALambdax already requested:";
                     std::cerr << " ignoring option '" << word1 << "' at line " << lineNum << "." << std::endl;
                  } else {
                     has_log_ALambdax_option = true;
                     kk >> ALambdax_upper;
                     has_ALambdax_upper = true;
                  }
               } else if (word1 == "LOGALAMBDANPTS") {
                  if (has_non_log_ALambdax_option) {
                     std::cerr << "WARNING: linear scan of ALambdax already requested:";
                     std::cerr << " ignoring option '" << word1 << "' at line " << lineNum << "." << std::endl;
                  } else {
                     has_log_ALambdax_option = true;
                     kk >> ALambdax_npts;
                     has_ALambdax_npts = true;
                  }
               } else if (word1 == "M2LL") {
                  kk >> MassWB_lower;
                  has_MassWB_lower = true;
               } else if (word1 == "M2UL") {
                  kk >> MassWB_upper;
                  has_MassWB_upper = true;
               } else if (word1 == "M2NPTS") {
                  kk >> MassWB_npts;
                  has_MassWB_npts = true;
               } else if (word1 == "YU") {
                  kk >> YuInput(0,0);
                  has_YuInput00 = true;
               } else if (word1 == "YC") {
                  kk >> YuInput(1,1);
                  has_YuInput11 = true;
               } else if (word1 == "YT") {
                  kk >> YuInput(2,2);
                  has_YuInput22 = true;
               } else if (word1 == "YD") {
                  kk >> YdInput(0,0);
                  has_YdInput00 = true;
               } else if (word1 == "YS") {
                  kk >> YdInput(1,1);
                  has_YdInput11 = true;
               } else if (word1 == "YB") {
                  kk >> YdInput(2,2);
                  has_YdInput22 = true;
               } else if (word1 == "YE") {
                  kk >> YeInput(0,0);
                  has_YeInput00 = true;
               } else if (word1 == "YMU") {
                  kk >> YeInput(1,1);
                  has_YeInput11 = true;
               } else if (word1 == "YTAU") {
                  kk >> YeInput(2,2);
                  has_YeInput22 = true;
               } else if (word1 == "G1") {
                  kk >> g1Input;
                  has_g1Input = true;
               } else if (word1 == "GN") {
                  kk >> gNInput;
                  has_gNInput = true;
               } else if (word1 == "G2"){
                  kk >> g2Input;
                  has_g2Input = true;
               } else if (word1 == "G3") {
                  kk >> g3Input;
                  has_g3Input = true;
               } else if (word1 == "HVEV") {
                  kk >> vInput;
                  has_vInput = true;
               } else if (word1 == "SVEV") {
                  kk >> vsInput;
                  has_vsInput = true;
               } else if (word1 == "MX") {
                  kk >> mxInput;
                  has_mxInput = true;
               } else if (word1 == "THETA") {
                  kk >> thetaInput;
                  has_thetaInput = true;
               } else { 
                  std::cerr << "WARNING: unrecognised option '" << word1 << "' requested:";
                  std::cerr << " ignoring it." << std::endl;
               }
            }
         } //< if (word1.find("#") != 0)
      } //< while(getline(cin,line)

     // Check inputs are valid
      if (!has_YuInput00) {
         std::cerr << "WARNING: Yu(0,0) value not found: default value is " << YuInput(0,0) << "." << std::endl;
      }

      if (!has_YuInput11) {
         std::cerr << "WARNING: Yu(1,1) value not found: default value is " << YuInput(1,1) << "." << std::endl;
      }

      if (!has_YuInput22) {
         std::cerr << "WARNING: Yu(2,2) value not found: default value is " << YuInput(2,2) << "." << std::endl;
      }

      if (!has_YdInput00) {
         std::cerr << "WARNING: Yd(0,0) value not found: default value is " << YdInput(0,0) << "." << std::endl;
      }

      if (!has_YdInput11) {
         std::cerr << "WARNING: Yd(1,1) value not found: default value is " << YdInput(1,1) << "." << std::endl;
      }

      if (!has_YdInput22) {
         std::cerr << "WARNING: Yd(2,2) value not found: default value is " << YdInput(2,2) << "." << std::endl;
      }

      if (!has_YeInput00) {
         std::cerr << "WARNING: Ye(0,0) value not found: default value is " << YeInput(0,0) << "." << std::endl;
      }

      if (!has_YeInput11) {
         std::cerr << "WARNING: Ye(1,1) value not found: default value is " << YeInput(1,1) << "." << std::endl;
      }

      if (!has_YeInput22) {
         std::cerr << "WARNING: Ye(2,2) value not found: default value is " << YeInput(2,2) << "." << std::endl;
      }

      if (!has_g1Input) {
         std::cerr << "WARNING: g_1 value not found: default value is " << g1Input << "." << std::endl;
      }

      if (!has_gNInput) {
         std::cerr << "WARNING: g_N value not found: default value is " << gNInput << "." << std::endl;
      }

      if (!has_g2Input) {
         std::cerr << "WARNING: g_2 value not found: default value is " << g2Input << "." << std::endl;
      }

      if (!has_g3Input) {
         std::cerr << "WARNING: g_3 value not found: default value is " << g3Input << "." << std::endl;
      }

      if (!has_vInput) {
         std::cerr << "WARNING: Higgs vev v value not found: default value is " << vInput << "." << std::endl;
      }

      if (!has_vsInput) {
         std::cerr << "WARNING: Singlet vev s value not found: default value is " << vsInput << "." << std::endl;
      }

      if (!has_mxInput) {
         std::cerr << "WARNING: Input scale MX value not found: default value is " << mxInput << "." << std::endl;
      }

      if (!has_thetaInput) {
         std::cerr << "WARNING: U(1) mixing theta value not found: default value is " << thetaInput << "." << std::endl;
      }

      if (!has_TanBeta_lower) {
         std::cerr << "WARNING: lower TanBeta limit not found: using default value " << TanBeta_lower << "." << std::endl;
      }

      if (!has_TanBeta_upper) {
         std::cerr << "WARNING: upper TanBeta limit not found: using default value " << TanBeta_upper << "." << std::endl;
      }

      if (!has_TanBeta_npts) {
         std::cerr << "WARNING: number of TanBeta points not found: using default value " << TanBeta_npts << "." << std::endl;
      }

      if (!has_Lambdax_lower) {
         std::cerr << "WARNING: lower Lambdax limit not found: using default value " << Lambdax_lower << "." << std::endl;
      }

      if (!has_Lambdax_upper) {
         std::cerr << "WARNING: upper Lambdax limit not found: using default value " << Lambdax_upper << "." << std::endl;
      }

      if (!has_Lambdax_npts) {
         std::cerr << "WARNING: number of Lambdax points not found: using default value " << Lambdax_npts << "." << std::endl;
      }

      if (!has_mq222_lower) {
         std::cerr << "WARNING: lower mq2(2,2) limit not found: using default value " << mq222_lower << "." << std::endl;
      }

      if (!has_mq222_upper) {
         std::cerr << "WARNING: upper mq2(2,2) limit not found: using default value " << mq222_upper << "." << std::endl;
      }

      if (!has_mq222_npts) {
         std::cerr << "WARNING: number of mq2(2,2) points not found: using default value " << mq222_npts << "." << std::endl;
      }

      if (!has_mu222_lower) {
         std::cerr << "WARNING: lower mu2(2,2) limit not found: using default value " << mu222_lower << "." << std::endl;
      }

      if (!has_mu222_upper) {
         std::cerr << "WARNING: upper mu2(2,2) limit not found: using default value " << mu222_upper << "." << std::endl;
      }

      if (!has_mu222_npts) {
         std::cerr << "WARNING: number of mu2(2,2) points not found: using default value " << mu222_npts << "." << std::endl;
      }

      if (!has_log_AYu22_option && !has_non_log_AYu22_option) {
         std::cerr << "WARNING: AYu(2,2) scan type not defined: using default ";
         if (!use_log_scan_AYu22) {
            std::cerr << "linear scan." << std::endl;
         } else {
            std::cerr << "log scan." << std::endl;
         }
      } else if (has_log_AYu22_option) {
         use_log_scan_AYu22 = true;
      } else {
         use_log_scan_AYu22 = false;
      }

      if (!has_AYu22_lower) {
         std::cerr << "WARNING: lower AYu(2,2) limit not found: using default value ";
         if (use_log_scan_AYu22) {
            std::cerr << -AYu22_upper << "." << std::endl;
         } else {
            std::cerr << AYu22_lower << "." << std::endl;
         }
      }

      if (!has_AYu22_upper) {
         std::cerr << "WARNING: upper AYu(2,2) limit not found: using default value " << AYu22_upper << "." << std::endl;
      }

      if (!has_AYu22_npts) {
         std::cerr << "WARNING: number of AYu(2,2) points not found: using default value " << AYu22_npts << "." << std::endl;
      }

      if (!has_log_ALambdax_option && !has_non_log_ALambdax_option) {
         std::cerr << "WARNING: ALambdax scan type not defined: using default ";
         if (!use_log_scan_ALambdax) {
            std::cerr << "linear scan." << std::endl;
         } else {
            std::cerr << "log scan." << std::endl;
         }
      } else if (has_log_ALambdax_option) {
         use_log_scan_ALambdax = true;
      } else {
         use_log_scan_ALambdax = false;
      }

      if (!has_ALambdax_lower) {
         std::cerr << "WARNING: lower ALambdax limit not found: using default value ";
         if (use_log_scan_ALambdax) {
            std::cerr << -ALambdax_upper << "." << std::endl;
         } else {
            std::cerr << ALambdax_lower << "." << std::endl;
         }
      }

      if (!has_ALambdax_upper) {
         std::cerr << "WARNING: upper ALambdax limit not found: using default value " << ALambdax_upper << "." << std::endl;
      }

      if (!has_ALambdax_npts) {
         std::cerr << "WARNING: number of ALambdax points not found: using default value " << ALambdax_npts << "." << std::endl;
      }

      if (!has_MassWB_lower) {
         std::cerr << "WARNING: lower M2 limit not found: using default value " << MassWB_lower << "." << std::endl;
      }

      if (!has_MassWB_upper) {
         std::cerr << "WARNING: upper M2 limit not found: using default value " << MassWB_upper << "." << std::endl;
      }

      if (!has_MassWB_npts) {
         std::cerr << "WARNING: number of M2 points not found: using default value " << MassWB_npts << "." << std::endl;
      }

      // Check ordering of upper and lower bounds
      double temp;
      if (TanBeta_lower > TanBeta_upper) {
         temp = TanBeta_lower;
         TanBeta_lower = TanBeta_upper;
         TanBeta_upper = temp;
      }

      if (Lambdax_lower > Lambdax_upper) {
         temp = Lambdax_lower;
         Lambdax_lower = Lambdax_upper;
         Lambdax_upper = temp;
      }

      if (ALambdax_lower > ALambdax_upper) {
         temp = ALambdax_lower;
         ALambdax_lower = ALambdax_upper;
         ALambdax_upper = temp;
      }

      if (mq222_lower > mq222_upper) {
         temp = mq222_lower;
         mq222_lower = mq222_upper;
         mq222_upper = temp;
      }

      if (mu222_lower > mu222_upper) {
         temp = mu222_lower;
         mu222_lower = mu222_upper;
         mu222_upper = temp;
      }

      if (AYu22_lower > AYu22_upper) {
         temp = AYu22_lower;
         AYu22_lower = AYu22_upper;
         AYu22_upper = temp;
      }

      if (MassWB_lower > MassWB_upper) {
         temp = MassWB_lower;
         MassWB_lower = MassWB_upper;
         MassWB_upper = temp;
      }

     // Check for valid values of tan(beta)
      if (TanBeta_lower < 1.0 || TanBeta_lower > 500.0) {
         std::cerr << "WARNING: invalid lower tan(beta) limit of " << TanBeta_lower << ":";
         std::cerr << " using default 2.0." << std::endl;
         TanBeta_lower = 2.0;
      }

      if (TanBeta_upper < 1.0 || TanBeta_upper > 500.0) {
         std::cerr << "WARNING: invalid upper tan(beta) limit of " << TanBeta_upper << ":";
         std::cerr << " using default 50.0." << std::endl;
         TanBeta_upper = 50.0;
      }

      // Check for valid numbers of points (i.e. at least 1 for each)
      if (TanBeta_npts < 1) {
         cerr << "WARNING: invalid number of tan(beta) points " << TanBeta_npts << " requested:";
         cerr << " using default 1." << endl;
         TanBeta_npts = 1;
      }

      if (Lambdax_npts < 1) {
         cerr << "WARNING: invalid number of mu points " << Lambdax_npts << " requested:";
         cerr << " using default 1." << endl;
         Lambdax_npts = 1;
      }

      if (ALambdax_npts < 1) {
         cerr << "WARNING: invalid number of A_lambda3 points " << ALambdax_npts << " requested:";
         cerr << " using default 1." << endl;
         ALambdax_npts = 1;
      }

      if (mq222_npts < 1) {
         cerr << "WARNING: invalid number of m_Q3^2 points " << mq222_npts << " requested:";
         cerr << " using default 1." << endl;
         mq222_npts = 1;
      }

      if (mu222_npts < 1) {
         cerr << "WARNING: invalid number of m_u3^2 points " << mu222_npts << " requested:";
         cerr << " using default 1." << endl;
         mu222_npts = 1;
      }

      if (AYu22_npts < 1) {
         cerr << "WARNING: invalid number of A_t points " << AYu22_npts << " requested:";
         cerr << " using default 1." << endl;
         AYu22_npts = 1;
      }

      if (MassWB_npts < 1) {
         cerr << "WARNING: invalid number of M_2 points " << MassWB_npts << " requested:";
         cerr << " using default 1." << endl;
         MassWB_npts = 1;
      }

      // Get values of scanned parameters
      Eigen::VectorXd TanBeta_vals = fill_linear_values(TanBeta_npts, TanBeta_lower, TanBeta_upper);
      Eigen::VectorXd Lambdax_vals = fill_linear_values(Lambdax_npts, Lambdax_lower, Lambdax_upper);
      Eigen::VectorXd ALambdax_vals;
      if (use_log_scan_ALambdax) {
         ALambdax_vals = fill_symmetric_log_values(ALambdax_npts, ALambdax_lower, ALambdax_upper);
         ALambdax_npts = 2 * ALambdax_npts;
      } else {
         ALambdax_vals = fill_linear_values(ALambdax_npts, ALambdax_lower, ALambdax_upper);
      }
      Eigen::VectorXd AYu22_vals;
      if (use_log_scan_AYu22) {
         AYu22_vals = fill_symmetric_log_values(AYu22_npts, AYu22_lower, AYu22_upper);
         AYu22_npts = 2 * AYu22_npts;
      } else {
         AYu22_vals = fill_linear_values(AYu22_npts, AYu22_lower, AYu22_upper);
      }
      Eigen::VectorXd mq222_vals = fill_linear_values(mq222_npts, mq222_lower, mq222_upper);
      Eigen::VectorXd mu222_vals = fill_linear_values(mu222_npts, mu222_lower, mu222_upper);
      Eigen::VectorXd MassWB_vals = fill_linear_values(MassWB_npts, MassWB_lower, MassWB_upper);

      std::vector<std::size_t> scan_dimensions 
         = {TanBeta_npts, Lambdax_npts, ALambdax_npts, AYu22_npts, 
            mq222_npts, mu222_npts, MassWB_npts};

      Grid_scanner scan(scan_dimensions);

      // print header line
               std::cout << "# "
                         << std::setw(12) << std::left << "TanBeta" << ' '
                         << std::setw(12) << std::left << "Lambdax(MS)" << ' '
                         << std::setw(12) << std::left << "ALambdax(MS)/GeV" << ' '
                         << std::setw(12) << std::left << "AYu22(MS)/GeV" << ' '
                         << std::setw(12) << std::left << "mq222(MS)/GeV^2" << ' '
                         << std::setw(12) << std::left << "mu222(MS)/GeV^2" << ' '
                         << std::setw(12) << std::left << "MassWB(MS)/GeV" << ' '
                         << std::setw(12) << std::left << "Lambdax(MX)" << ' '
                         << std::setw(12) << std::left << "Kappa22(MX)" << ' '
                         << std::setw(12) << std::left << "Mhh(1)/GeV" << ' '
                         << std::setw(12) << std::left << "Mhh(2)/GeV" << ' '
                         << std::setw(12) << std::left << "Mhh(3)/GeV" << ' '
                         << std::setw(12) << std::left << "mHd2/GeV^2" << ' '
                         << std::setw(12) << std::left << "mHu2/GeV^2" << ' '
                         << std::setw(12) << std::left << "ms2/GeV^2" << ' '
                         << std::setw(12) << std::left << "f1" << ' '
                         << std::setw(12) << std::left << "f2" << ' '
                         << std::setw(12) << std::left << "f3" << ' '
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
                         << std::setw(12) << std::left << "tachyon(Su)" << ' '
                         << std::setw(12) << std::left << "tachyon(Sd)" << ' '
                         << std::setw(12) << std::left << "tachyon(VZ)" << ' '
                         << std::setw(12) << std::left << "tachyon(VZp)" << ' '
                         << std::setw(12) << std::left << "tachyon(hh)" << ' '
                         << std::setw(12) << std::left << "tachyon(Ah)" << ' '
                         << std::setw(12) << std::left << "tuning_err" << ' '
                         << std::setw(12) << std::left << "error"
                         << '\n';

      // get default inputs
      lowE6SSM_input_parameters input = get_default_inputs();

      // set requested values of vs and gN
      input.gNInput = gNInput;
      input.vsInput = vsInput;

      std::vector<std::size_t> position;
      while (!scan.has_finished()) {

         bool has_serious_problem = false;
         position = scan.get_position();
         input.TanBeta = TanBeta_vals(position.at(0));
         input.LambdaxInput = Lambdax_vals(position.at(1));
         input.TLambdaxInput = Lambdax_vals(position.at(1)) * ALambdax_vals(position.at(2));
         input.TYuInput(2,2) = YuInput(2,2) * AYu22_vals(position.at(3));
         input.mq2Input(2,2) = mq222_vals(position.at(4));
         input.mu2Input(2,2) = mu222_vals(position.at(5));
         input.MassWBInput = MassWB_vals(position.at(6));

         lowE6SSM<Two_scale> model;
         model.set_loops(2);
         model.set_ewsb_loop_order(2);

         initialize_e6_charges(thetaInput, input);

         initialize_model_from_input(input, model);
         
         model.set_g1(g1Input);
         model.set_g2(g2Input);
         model.set_g3(g3Input);

         model.set_Yu(YuInput);
         model.set_Yd(YdInput);
         model.set_Ye(YeInput);

         model.set_vd(vInput / Sqrt(1.0 + Sqr(TanBeta_vals(position.at(0)))));
         model.set_vu(vInput * TanBeta_vals(position.at(0)) / Sqrt(1.0 + Sqr(TanBeta_vals(position.at(0)))));

         // calculate stop masses, if they are tachyonic skip the rest
         // of the calculation
         double mst1 = 0.;
         double mst2 = 0.;
         bool stop_tachyon = calculate_3rd_generation_stop_masses(model, mst1, mst2);
         if (stop_tachyon) {
            model.get_problems().flag_tachyon(lowE6SSM_info::Su);
         } else {
            model.get_problems().unflag_tachyon(lowE6SSM_info::Su);
            double msusy = Sqrt(mst1 * mst2);
            model.set_scale(msusy);

            double msb1 = 0.;
            double msb2 = 0.;
            bool sbottom_tachyon = calculate_3rd_generation_sbottom_masses(model, msb1, msb2);
            if (sbottom_tachyon) {
               model.get_problems().flag_tachyon(lowE6SSM_info::Sd);
            } else {
               model.get_problems().unflag_tachyon(lowE6SSM_info::Sd);
            }

            // solve ewsb
            lowE6SSM_ew_derivs ew_derivs(model);

            double mHd2 = ew_derivs.get_model().get_mHd2();
            double mHu2 = ew_derivs.get_model().get_mHu2();
            double ms2 = ew_derivs.get_model().get_ms2();

            model.set_mHd2(mHd2);
            model.set_mHu2(mHu2);
            model.set_ms2(ms2);            

            // calculate spectrum
            bool higgs_tachyon = ew_derivs.calculate_MHiggs();
            if (higgs_tachyon) {
               model.get_problems().flag_tachyon(lowE6SSM_info::hh);
            } else {
               model.get_problems().unflag_tachyon(lowE6SSM_info::hh);
            }
            model.solve_ewsb_tree_level();

            model.calculate_MVZ();
            model.calculate_MVZp();
            model.calculate_MChi();
            model.calculate_MCha();
            model.calculate_Mhh();
            model.calculate_MAh();

            model.set_mHd2(mHd2);
            model.set_mHu2(mHu2);
            model.set_ms2(ms2);  

            // calculate fine tuning
            lowE6SSM_tuning_calculator tuning_calc(model);
            tuning_calc.set_tuning_scale(msusy);
            tuning_calc.set_input_scale(mxInput);

            bool tuning_problem;
            try {
               tuning_problem = tuning_calc.calculate_fine_tunings_approximately();
            } catch (const std::string & a) {
               std::cerr << "WARNING: serious numerical problem encountered in fine tuning calculation.\n";
               has_serious_problem = true;
            } catch (const char* a) {
               std::cerr << "WARNING: serious numerical problem encountered in fine tuning calculation.\n";
               has_serious_problem = true;
            } catch (...) {
               std::cerr << "WARNING: serious numerical problem encountered in fine tuning calculation.\n";
               has_serious_problem = true;
            }
            std::map<lowE6SSM_info::Parameters,double> fine_tunings = tuning_calc.get_fine_tunings();
            double max_tuning = maximum_tuning(fine_tunings);

            double Lambdax_at_MX = tuning_calc.get_Lambdax_at_input_scale();
            double Kappa22_at_MX = tuning_calc.get_Kappa_at_input_scale()(2,2);

            // print results
            if (!has_serious_problem) {
               const Problems<lowE6SSM_info::NUMBER_OF_PARTICLES>& problems = model.get_problems();
               std::cout << " "
                         << std::setw(12) << std::left << TanBeta_vals(position.at(0)) << ' '
                         << std::setw(12) << std::left << Lambdax_vals(position.at(1)) << ' '
                         << std::setw(12) << std::left << ALambdax_vals(position.at(2)) << ' '
                         << std::setw(12) << std::left << AYu22_vals(position.at(3)) << ' '
                         << std::setw(12) << std::left << mq222_vals(position.at(4)) << ' '
                         << std::setw(12) << std::left << mu222_vals(position.at(5)) << ' '
                         << std::setw(12) << std::left << MassWB_vals(position.at(6)) << ' '
                         << std::setw(12) << std::left << Lambdax_at_MX << ' '
                         << std::setw(12) << std::left << Kappa22_at_MX << ' '
                         << std::setw(12) << std::left << ew_derivs.get_MHiggs()(0) << ' '
                         << std::setw(12) << std::left << ew_derivs.get_MHiggs()(1) << ' '
                         << std::setw(12) << std::left << ew_derivs.get_MHiggs()(2) << ' '
                         << std::setw(12) << std::left << mHd2 << ' '
                         << std::setw(12) << std::left << mHu2 << ' '
                         << std::setw(12) << std::left << ms2 << ' '
                         << std::setw(12) << std::left << ew_derivs.get_ewsb_condition_1() << ' '
                         << std::setw(12) << std::left << ew_derivs.get_ewsb_condition_2() << ' '
                         << std::setw(12) << std::left << ew_derivs.get_ewsb_condition_3() << ' '
                         << std::setw(12) << std::left << msusy << ' '
                         << std::setw(12) << std::left << mxInput << ' '
                         << std::setw(12) << std::left << vsInput << ' '
                         << std::setw(12) << std::left << model.get_MVZ() << ' '
                         << std::setw(12) << std::left << model.get_MVZp() << ' '
                         << std::setw(12) << std::left << mst1 << ' '
                         << std::setw(12) << std::left << mst2 << ' '
                         << std::setw(12) << std::left << msb1 << ' '
                         << std::setw(12) << std::left << msb2 << ' '
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
                         << std::setw(12) << std::left << problems.is_tachyon(lowE6SSM_info::Su) << ' '
                         << std::setw(12) << std::left << problems.is_tachyon(lowE6SSM_info::Sd) << ' '
                         << std::setw(12) << std::left << problems.is_tachyon(lowE6SSM_info::VZ) << ' '
                         << std::setw(12) << std::left << problems.is_tachyon(lowE6SSM_info::VZp) << ' '
                         << std::setw(12) << std::left << problems.is_tachyon(lowE6SSM_info::hh) << ' '
                         << std::setw(12) << std::left << problems.is_tachyon(lowE6SSM_info::Ah) << ' '
                         << std::setw(12) << std::left << tuning_problem << ' '
                         << std::setw(12) << std::left << (problems.have_serious_problem() || has_serious_problem)
                         << '\n';
                  }
         }
      }

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
