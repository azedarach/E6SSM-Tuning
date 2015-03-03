// ====================================================================
// Calculates the fine tuning at a set of points centred on 
// read in points. For reference:
// ====================================================================

#include "lowE6SSM_input_parameters.hpp"
#include "lowE6SSM_spectrum_generator.hpp"
#include "lowE6SSM_two_scale_tuning_calculator.hpp"

#include "error.hpp"
#include "scan.hpp"
#include "lowe.h"

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include <chrono>
#include <iostream>
#include <map>
#include <random>

namespace flexiblesusy {
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
      const auto Yd = model.get_Yd();
      const auto Ye = model.get_Ye();
      const auto Yu = model.get_Yu();

      const auto AYd = input.AYdInput;
      const auto AYe = input.AYeInput;
      const auto AYu = input.AYuInput;

      model.set_input_parameters(input);

      model.set_Kappa(input.KappaInput);
      model.set_Lambda12(input.Lambda12Input);
      model.set_Lambdax(input.LambdaxInput);
      
      model.set_MuPr(input.MuPrInput);
      
      model.set_gN(input.gNInput);
      model.set_vs(input.vsInput);
      
      model.set_TYd((Yd.array() * AYd.array()).matrix());
      model.set_TYe((Ye.array() * AYe.array()).matrix());
      model.set_TKappa(input.TKappaInput);
      model.set_TLambda12(input.TLambda12Input);
      model.set_TLambdax(input.TLambdaxInput);
      model.set_TYu((Yu.array() * AYu.array()).matrix());
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
   
   lowE6SSM_input_parameters get_default_inputs(double theta)
   {
      lowE6SSM_input_parameters input;
      
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
   
   inline void trim(std::string& str)
   {
      std::size_t startpos = str.find_first_not_of(" \t\n\v\f\r");
      if (startpos != std::string::npos) str.erase(0, startpos);
      std::size_t endpos = str.find_last_not_of(" \t\n\v\f\r");
      if (endpos != std::string::npos) str.erase(endpos+1);
   }

   double get_random_Lambdax(const lowE6SSM_input_parameters& input, double Lambdax_width, std::default_random_engine& generator)
   {
      const double LambdaxInput = input.LambdaxInput;

      std::uniform_real_distribution<double> Lambdax_distribution(LambdaxInput - Lambdax_width, LambdaxInput + Lambdax_width);

      return Lambdax_distribution(generator);
   }

   double get_random_TLambdax(const lowE6SSM_input_parameters& input, double ALambdax_width, std::default_random_engine& generator)
   {
      const double LambdaxInput = input.LambdaxInput;
      const double TLambdaxInput = input.TLambdaxInput;
      const double ALambdaxInput = TLambdaxInput / LambdaxInput;

      std::uniform_real_distribution<double> ALambdax_distribution(ALambdaxInput - ALambdax_width, ALambdaxInput + ALambdax_width);

      return LambdaxInput * ALambdax_distribution(generator);
   }

   double get_random_AYu22(const lowE6SSM_input_parameters& input, double AYu22_width, std::default_random_engine& generator)
   {
      const double AYu22Input = input.AYuInput(2,2);

      std::uniform_real_distribution<double> AYu22_distribution(AYu22Input - AYu22_width, AYu22Input + AYu22_width);

      return AYu22_distribution(generator);
   }

   double get_random_mq222(const lowE6SSM_input_parameters& input, double mq222_width, std::default_random_engine& generator)
   {
      const double mq222Input = input.mq2Input(2,2);

      std::uniform_real_distribution<double> mq222_distribution(mq222Input - mq222_width, mq222Input + mq222_width);

      return mq222_distribution(generator);
   }

   double get_random_mu222(const lowE6SSM_input_parameters& input, double mu222_width, std::default_random_engine& generator)
   {
      const double mu222Input = input.mu2Input(2,2);

      std::uniform_real_distribution<double> mu222_distribution(mu222Input - mu222_width, mu222Input + mu222_width);

      return mu222_distribution(generator);
   }

}

int main()
{
   using namespace flexiblesusy;
   using namespace softsusy;

   unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
   std::default_random_engine generator(seed);

   // Note:
   // U(1)_N: theta = ArcTan(Sqrt(15.0))
   // U(1)_I: theta = PI + ArcTan(Sqrt(0.6))
   // U(1)_psi: theta = 0.5 * PI
   // U(1)_eta: theta = -ArcTan(Sqrt(5.0 / 3.0))
   const double theta = -ArcTan(Sqrt(5.0 / 3.0));

   const std::size_t num_points = 10;

   Eigen::Matrix<double,3,3> YdInput; 
   Eigen::Matrix<double,3,3> YeInput;
   Eigen::Matrix<double,3,3> YuInput;

   YdInput(0,0) = 0.0;
   YdInput(0,1) = 0.0;
   YdInput(0,2) = 0.0;
   YdInput(1,0) = 0.0;
   YdInput(1,1) = 0.0;
   YdInput(1,2) = 0.0;
   YdInput(2,0) = 0.0;
   YdInput(2,1) = 0.0;
   YdInput(2,2) = 1.22357994e-01;

   YeInput(0,0) = 0.0;
   YeInput(0,1) = 0.0;
   YeInput(0,2) = 0.0;
   YeInput(1,0) = 0.0;
   YeInput(1,1) = 0.0;
   YeInput(1,2) = 0.0;
   YeInput(2,0) = 0.0;
   YeInput(2,1) = 0.0;
   YeInput(2,2) = 9.88254141e-02;

   YuInput(0,0) = 0.0;
   YuInput(0,1) = 0.0;
   YuInput(0,2) = 0.0;
   YuInput(1,0) = 0.0;
   YuInput(1,1) = 0.0;
   YuInput(1,2) = 0.0;
   YuInput(2,0) = 0.0;
   YuInput(2,1) = 0.0;
   YuInput(2,2) = 8.74991501e-01;

   const double g1Input = 4.71880484e-01;
   const double g2Input = 6.36109014e-01;
   const double g3Input = 9.99300252e-01;
   const double gNInput = 5.10342303e-01;
   const double vInput = 2.41460403e+02;

   lowE6SSM_input_parameters input;
   QedQcd oneset;
   oneset.toMz();

   std::size_t TanBeta_col = 0;
   std::size_t Lambdax_col = 0;
   std::size_t ALambdax_col = 0;
   std::size_t AYu22_col = 0;
   std::size_t mq222_col = 0;
   std::size_t mu222_col = 0;
   std::size_t MassWB_col = 0;
   std::size_t MX_col = 0;
   std::size_t vs_col = 0;

   // header line
   std::string comment_line;
   std::getline(std::cin, comment_line);

   // discard leading '#'
   comment_line = comment_line.substr(1);

   // get field names
   std::vector<std::string> field_names;
   boost::split(field_names, comment_line, boost::is_any_of(" \n\t"), boost::token_compress_on);
   field_names.erase(field_names.begin());
   field_names.erase(field_names.end());
   for (std::size_t i = 0; i < field_names.size(); ++i) {
      trim(field_names[i]);
      if (field_names[i] == "TanBeta") {
         TanBeta_col = i;
      } else if (field_names[i] == "Lambdax(MS)") {
         Lambdax_col = i;
      } else if (field_names[i] == "ALambdax(MS)/GeV") {
         ALambdax_col = i;
      } else if (field_names[i] == "AYu22(MS)/GeV") {
         AYu22_col = i;
      } else if (field_names[i] == "mq222(MS)/GeV^2") {
         mq222_col = i;
      } else if (field_names[i] == "mu222(MS)/GeV^2") {
         mu222_col = i;
      } else if (field_names[i] == "MassWB(MS)/GeV") {
         MassWB_col = i;
      } else if (field_names[i] == "MX/GeV") {
         MX_col = i;
      } else if (field_names[i] == "s/GeV") {
         vs_col = i;
      }
   }
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
             << std::setw(12) << std::left << "tachyon(Hpm)" << ' '
             << std::setw(12) << std::left << "tuning_err" << ' '
             << std::setw(12) << std::left << "error"
             << '\n';

      // read in points from file and calculate tuning
      std::string line;
      while (std::getline(std::cin, line)) {

         input = get_default_inputs(theta);

         // set requested value of gN
         input.gNInput = gNInput;

         // set scanned parameter values
         std::vector<std::string> values;
         boost::split(values, line, boost::is_any_of(" \n\t"), boost::token_compress_on);
         values.erase(values.begin());
         values.erase(values.end());

         try {
            input.TanBeta = boost::lexical_cast<double>(values[TanBeta_col]);
         } catch (const boost::bad_lexical_cast& error) {
            WARNING("Ignoring invalid input '" + values[TanBeta_col] + "'");
            continue;
         }

         try {
            input.LambdaxInput = boost::lexical_cast<double>(values[Lambdax_col]);
            try {
               input.TLambdaxInput = boost::lexical_cast<double>(values[ALambdax_col]);
               input.TLambdaxInput *= input.LambdaxInput;

            } catch (const boost::bad_lexical_cast& error) {
               WARNING("Ignoring invalid input '" + values[ALambdax_col] + "'");
               continue;
            }
         } catch (const boost::bad_lexical_cast& error) {
            WARNING("Ignoring invalid input '" + values[Lambdax_col] + "'");
            continue;
         }
         
         try {
            input.AYuInput(2,2) = boost::lexical_cast<double>(values[AYu22_col]);
         } catch (const boost::bad_lexical_cast& error) {
            WARNING("Ignoring invalid input '" + values[AYu22_col] + "'");
            continue;
         }
         
         try {
            input.mq2Input(2,2) = boost::lexical_cast<double>(values[mq222_col]);
         } catch (const boost::bad_lexical_cast& error) {
            WARNING("Ignoring invalid input '" + values[mq222_col] + "'");
            continue;
         }
         
         try {
            input.mu2Input(2,2) = boost::lexical_cast<double>(values[mu222_col]);
         } catch (const boost::bad_lexical_cast& error) {
            WARNING("Ignoring invalid input '" + values[mu222_col] + "'");
            continue;
         }
         
         try {
            input.MassWBInput = boost::lexical_cast<double>(values[MassWB_col]);
         } catch (const boost::bad_lexical_cast& error) {
            WARNING("Ignoring invalid input '" + values[MassWB_col] + "'");
            continue;
         }

         double mx = 0.;
         try {
            mx = boost::lexical_cast<double>(values[MX_col]);
         } catch (const boost::bad_lexical_cast& error) {
            WARNING("Ignoring invalid input '" + values[MX_col] + "'");
            continue;
         }
         
         try {
            input.vsInput = boost::lexical_cast<double>(values[vs_col]);
         } catch (const boost::bad_lexical_cast& error) {
            WARNING("Ignoring invalid input '" + values[vs_col] + "'");
            continue;
         }
         
         // generate additional points from read-in point
         const double ALambdax_width = 1000.0;
         const double AYu22_width = 1000.0;
         const double Lambdax_width = 0.05;
         const double mq222_width = 1000.0;
         const double mu222_width = 1000.0;

         for (std::size_t i = 0; i < num_points; ++i) {
            lowE6SSM_input_parameters random_input = input;
            random_input.TLambdaxInput = get_random_TLambdax(input, ALambdax_width, generator);
            random_input.AYuInput(2,2) = get_random_AYu22(input, AYu22_width, generator);

            // additional randomised values
            random_input.LambdaxInput = get_random_Lambdax(input, Lambdax_width, generator);
            random_input.mq2Input(2,2) = get_random_mq222(input, mq222_width, generator);
            random_input.mu2Input(2,2) = get_random_mu222(input, mu222_width, generator);

            bool has_serious_problem = false;

            lowE6SSM<Two_scale> model;
            model.set_loops(2);
            model.set_ewsb_loop_order(2);

            model.set_g1(g1Input);
            model.set_g2(g2Input);
            model.set_g3(g3Input);

            model.set_Yu(YuInput);
            model.set_Yd(YdInput);
            model.set_Ye(YeInput);

            initialize_model_from_input(random_input, model);

            model.set_vd(vInput / Sqrt(1.0 + Sqr(input.TanBeta)));
            model.set_vu(vInput * input.TanBeta / Sqrt(1.0 + Sqr(input.TanBeta)));

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
               model.calculate_MHpm();

               model.set_mHd2(mHd2);
               model.set_mHu2(mHu2);
               model.set_ms2(ms2);  
               
               // calculate fine tuning
               lowE6SSM_tuning_calculator tuning_calc(model);
               tuning_calc.set_tuning_scale(msusy);
               tuning_calc.set_input_scale(mx);

               bool tuning_problem = false;
               try {
               tuning_problem = tuning_calc.calculate_fine_tunings_approximately();
               //tuning_problem = tuning_calc.calculate_fine_tunings_numerically();
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
               if (!has_serious_problem && !tuning_problem) {
                  const Problems<lowE6SSM_info::NUMBER_OF_PARTICLES>& problems = model.get_problems();
                  std::cout << " "
                            << std::setw(12) << std::left << random_input.TanBeta << ' '
                            << std::setw(12) << std::left << random_input.LambdaxInput << ' '
                            << std::setw(12) << std::left << random_input.TLambdaxInput / random_input.LambdaxInput << ' '
                            << std::setw(12) << std::left << random_input.AYuInput(2,2) << ' '
                            << std::setw(12) << std::left << random_input.mq2Input(2,2) << ' '
                            << std::setw(12) << std::left << random_input.mu2Input(2,2) << ' '
                            << std::setw(12) << std::left << random_input.MassWBInput << ' '
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
                            << std::setw(12) << std::left << mx << ' '
                            << std::setw(12) << std::left << random_input.vsInput << ' '
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
                            << std::setw(12) << std::left << problems.is_tachyon(lowE6SSM_info::Hpm) << ' '
                            << std::setw(12) << std::left << tuning_problem << ' '
                            << std::setw(12) << std::left << (problems.have_serious_problem() || has_serious_problem)
                            << '\n';
               }
               
            } // if (stop_tachyon)

         } // for (std::size_t i = 0; i < num_points; ++i)

      } // while(std::getline(std::cin, line))

      return 0;
}
