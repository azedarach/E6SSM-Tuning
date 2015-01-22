// ====================================================================
// Calculates the full spectrum and the fine tuning in the E6SSM
// ====================================================================

// File generated at Sun 24 Aug 2014 16:15:32

#include "lowE6SSM_input_parameters.hpp"
#include "lowE6SSM_spectrum_generator.hpp"
#include "lowE6SSM_two_scale_tuning_calculator.hpp"

#include "error.hpp"
#include "scan.hpp"
#include "lowe.h"

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include <iostream>
#include <map>

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

   double find_closest_TanBeta(double tb_val, const std::vector<double>& vals) 
   {
      std::vector<double> differences(vals);
      double min_diff = 1.0e5;
      std::size_t min_index = 0;
      for (std::size_t i = 0; i < differences.size(); ++i) {
         if (Abs(tb_val - vals[i]) < min_diff) {
            min_diff = Abs(tb_val - vals[i]);
            min_index = i;
         }
      }

      return vals[min_index];
   }

}

int main()
{
   using namespace flexiblesusy;
   using namespace softsusy;
   typedef Two_scale algorithm_type;

   lowE6SSM_input_parameters input;
   QedQcd oneset;
   oneset.toMz();

   lowE6SSM_spectrum_generator<algorithm_type> spectrum_generator;
   spectrum_generator.set_precision_goal(1.0e-4);
   spectrum_generator.set_max_iterations(0);         // 0 == automatic
   spectrum_generator.set_calculate_sm_masses(0);    // 0 == no
   spectrum_generator.set_parameter_output_scale(0); // 0 == susy scale

   std::vector<double> TanBeta_vals = {2, 4, 6, 8, 10, 15, 20, 30, 40, 50};

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
                << std::setw(12) << std::left << "Mhh(1)/GeV" << ' '
                << std::setw(12) << std::left << "Mhh(2)/GeV" << ' '
                << std::setw(12) << std::left << "Mhh(3)/GeV" << ' '
                << std::setw(12) << std::left << "mHd2/GeV^2" << ' '
                << std::setw(12) << std::left << "mHu2/GeV^2" << ' '
                << std::setw(12) << std::left << "ms2/GeV^2" << ' '
                << std::setw(12) << std::left << "g1" << ' '
                << std::setw(12) << std::left << "gN" << ' '
                << std::setw(12) << std::left << "MS/GeV" << ' '
                << std::setw(12) << std::left << "MX/GeV" << ' '
                << std::setw(12) << std::left << "s/GeV" << ' '
                << std::setw(12) << std::left << "MVZ/GeV" << ' '
                << std::setw(12) << std::left << "MVZp/GeV" << ' '
                << std::setw(12) << std::left << "MSu(1)/GeV" << ' '
                << std::setw(12) << std::left << "MSu(2)/GeV" << ' '
                << std::setw(12) << std::left << "MSd(1)/GeV" << ' '
                << std::setw(12) << std::left << "MSd(2)/GeV" << ' '
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
                << std::setw(12) << std::left << "error" << '\n';

// read in points from file and calculate spectrum
   std::string line;
   while (std::getline(std::cin, line)) {
      input = get_default_inputs();
      
      // set scanned parameter values
      std::vector<std::string> values;
      boost::split(values, line, boost::is_any_of(" \n\t"), boost::token_compress_on);
      values.erase(values.begin());
      values.erase(values.end());

      try {
         input.TanBeta = find_closest_TanBeta(boost::lexical_cast<double>(values[TanBeta_col]), TanBeta_vals);
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

      // calculate spectrum
      spectrum_generator.run(oneset, input);

      const lowE6SSM<algorithm_type>& model = spectrum_generator.get_model();
      const lowE6SSM_physical& pole_masses = model.get_physical();
      const Problems<lowE6SSM_info::NUMBER_OF_PARTICLES>& problems
         = spectrum_generator.get_problems();

      const bool error = problems.have_serious_problem();

      std::map<lowE6SSM_info::Parameters,double> fine_tunings;
      double max_tuning = 0;
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

         // calculate fine tuning
         lowE6SSM_tuning_calculator tuning_calc(tuning_model);
         tuning_calc.set_tuning_scale(Sqrt(model.get_MSu()(0) * model.get_MSu()(1)));
         tuning_calc.set_input_scale(mx);
         tuning_calc.set_tuning_ewsb_loop_order(1);
         tuning_calc.set_tuning_beta_loop_order(2);
         try {
            tuning_problem = tuning_calc.calculate_fine_tunings_approximately();
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
      
      // print results
      std::cout << " "
                << std::setw(12) << std::left << input.TanBeta << ' '
                << std::setw(12) << std::left << model.get_Lambdax() << ' '
                << std::setw(12) << std::left << model.get_TLambdax() / model.get_Lambdax() << ' '
                << std::setw(12) << std::left << model.get_TYu(2,2) / model.get_Yu(2,2) << ' '
                << std::setw(12) << std::left << model.get_mq2(2,2) << ' '
                << std::setw(12) << std::left << model.get_mu2(2,2) << ' '
                << std::setw(12) << std::left << model.get_MassWB() << ' '
                << std::setw(12) << std::left << pole_masses.Mhh(0) << ' '
                << std::setw(12) << std::left << pole_masses.Mhh(1) << ' '
                << std::setw(12) << std::left << pole_masses.Mhh(2) << ' '
                << std::setw(12) << std::left << model.get_mHd2() << ' '
                << std::setw(12) << std::left << model.get_mHu2() << ' '
                << std::setw(12) << std::left << model.get_ms2() << ' '
                << std::setw(12) << std::left << model.get_g1() << ' '
                << std::setw(12) << std::left << model.get_gN() << ' '
                << std::setw(12) << std::left << Sqrt(model.get_MSu()(0) * model.get_MSu()(1)) << ' '
                << std::setw(12) << std::left << mx << ' '
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
   }

   return 0;
}
