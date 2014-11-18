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
      input.AYuInput = Eigen::Matrix<double,3,3>::Zero();
      
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
      }
   }
   
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
      } catch (const boost::bad_lexical_cast& error) {
         WARNING("Ignoring invalid input '" + values[Lambdax_col] + "'");
         continue;
      }

      try {
         inputALambdaxInput = boost::lexical_cast<double>(values[ALambdax_col]);
      } catch (const boost::bad_lexical_cast& error) {
         WARNING("Ignoring invalid input '" + values[ALambdax_col] + "'");
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

      // calculate spectrum
      std::cout << "TanBeta = " << input.TanBeta << "\n";
      std::cout << "Lambdax = " << input.LambdaxInput << "\n";
      std::cout << "TLambdax = " << input.TLambdaxInput << "\n";
      std::cout << "TYu22 = " << input.TYuInput(2,2) << "\n";
      std::cout << "mq222 = " << input.mq2Input(2,2) << "\n";
      std::cout << "mu222 = " << input.mu2Input(2,2) << "\n";
      std::cout << "MassWB = " << input.MassWBInput << "\n";
      spectrum_generator.run(oneset, input);

      const lowE6SSM<algorithm_type>& model = spectrum_generator.get_model();
      const lowE6SSM_physical& pole_masses = model.get_physical();
      const Problems<lowE6SSM_info::NUMBER_OF_PARTICLES>& problems
         = spectrum_generator.get_problems();

      const bool error = problems.have_serious_problem();

      std::map<lowE6SSM_info::Parameters,double> fine_tunings;
      double max_tuning;
      if (!error) {
         // calculate fine tuning
         lowE6SSM_tuning_calculator tuning_calc(model);
         tuning_calc.set_tuning_scale(model.get_scale());
         tuning_calc.set_input_scale(mx);
         bool tuning_problem;
         try {
            tuning_problem = tuning_calc.calculate_fine_tunings_numerically();
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
      if (error) {
         std::cout << problems << "\n";
      }
   }

   return 0;
}
