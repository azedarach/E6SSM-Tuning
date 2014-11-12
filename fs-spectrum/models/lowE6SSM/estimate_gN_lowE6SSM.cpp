// ====================================================================
// Estimate g_N at the SUSY scale based on 1-loop gauge unification,
// (no threshold corrections).
// For comparison, for \theta_{E_6}, should have g0 ~
// ====================================================================

#include "lowE6SSM_input_parameters.hpp"
#include "lowE6SSM_slha_io.hpp"
#include "lowE6SSM_spectrum_generator.hpp"

#include "command_line_options.hpp"
#include "error.hpp"
#include "lowe.h"
#include "spectrum_generator_settings.hpp"

#include <iostream>
#include <cstdlib>

int main(int argc, const char* argv[])
{
   using namespace flexiblesusy;
   using namespace softsusy;
   typedef Two_scale algorithm_type;

   Command_line_options options(argc, argv);
   if (options.must_print_model_info())
      lowE6SSM_info::print(std::cout);
   if (options.must_exit())
      return options.status();

   const std::string rgflow_file(options.get_rgflow_file());
   const std::string slha_input_file(options.get_slha_input_file());
   const std::string slha_output_file(options.get_slha_output_file());
   const std::string spectrum_file(options.get_spectrum_file());
   lowE6SSM_slha_io slha_io;
   Spectrum_generator_settings spectrum_generator_settings;
   QedQcd oneset;
   lowE6SSM_input_parameters input;

   if (slha_input_file.empty()) {
      ERROR("No SLHA input file given!\n"
            "   Please provide one via the option --slha-input-file=");
      return EXIT_FAILURE;
   }

   try {
      slha_io.read_from_file(slha_input_file);
      slha_io.fill(oneset);
      slha_io.fill(input);
      slha_io.fill(spectrum_generator_settings);
   } catch (const Error& error) {
      ERROR(error.what());
      return EXIT_FAILURE;
   }

   oneset.toMz(); // run SM fermion masses to MZ

   lowE6SSM_spectrum_generator<algorithm_type> spectrum_generator;
   spectrum_generator.set_precision_goal(
      spectrum_generator_settings.get(Spectrum_generator_settings::precision));
   spectrum_generator.set_max_iterations(
      spectrum_generator_settings.get(Spectrum_generator_settings::max_iterations));
   spectrum_generator.set_calculate_sm_masses(
      spectrum_generator_settings.get(Spectrum_generator_settings::calculate_sm_masses) >= 1.0);
   spectrum_generator.set_input_scale(
      slha_io.get_input_scale());
   spectrum_generator.set_parameter_output_scale(
      slha_io.get_parameter_output_scale());
   spectrum_generator.set_pole_mass_loop_order(
      spectrum_generator_settings.get(Spectrum_generator_settings::pole_mass_loop_order));
   spectrum_generator.set_ewsb_loop_order(
      spectrum_generator_settings.get(Spectrum_generator_settings::ewsb_loop_order));
   spectrum_generator.set_beta_loop_order(1);
   spectrum_generator.set_threshold_corrections_loop_order(1);

   spectrum_generator.run(oneset, input);

   lowE6SSM<algorithm_type> model = spectrum_generator.get_model();
   const Problems<lowE6SSM_info::NUMBER_OF_PARTICLES>& problems
      = spectrum_generator.get_problems();

   if (!problems.have_serious_problem()) {
      const double susy_scale = model.get_scale();
      model.set_loops(1);

      // estimate MX in the 1-loop approximation
      lowE6SSM_soft_parameters betas = model.calc_beta();
      const double beta1 = betas.get_g1() / (oneOver16PiSqr * Power(model.get_g1(), 3));
      const double beta2 = betas.get_g2() / (oneOver16PiSqr * Power(model.get_g2(), 3));
      const double oneOverg12 = 1.0 / Sqr(model.get_g1());
      const double oneOverg22 = 1.0 / Sqr(model.get_g2());

      double MX = susy_scale * std::exp(0.5 * (oneOverg12 - oneOverg22) / (oneOver16PiSqr * (beta1 - beta2)));
      
      std::cout << " gauge = " << oneOverg12 - oneOverg22 << "\n";
      std::cout << "beta1 = " << beta1 << "\n";
      std::cout << "beta2 = " << beta2 << "\n";
      std::cout << "beta = " << beta1 - beta2 << "\n";
      std::cout << "factor = " << 0.5 * (oneOverg12 - oneOverg22) / (oneOver16PiSqr * (beta1 - beta2)) << "\n";
      std::cout << "MS = " << susy_scale << "\n";
      std::cout << "MX = " << MX << "\n";

      // attempt to unify gauge couplings by varying gN
      model.run_to(MX);
      std::cout << "Initially, gN(MX) = " << model.get_gN() << "\n";
      model.set_gN(model.get_g1());
      std::cout << "g1(MX) = " << model.get_g1() << "\n";
      std::cout << "g2(MX) = " << model.get_g2() << "\n";
      std::cout << "g3(MX) = " << model.get_g3() << "\n";
      std::cout << "gN(MX) = " << model.get_gN() << "\n";
      model.run_to(susy_scale);std::cout.precision(12);
      std::cout << "g1(MS) = " << model.get_g1() << "\n";
      std::cout << "g2(MS) = " << model.get_g2() << "\n";
      std::cout << "g3(MS) = " << model.get_g3() << "\n";
      std::cout << "gN(MS) = " << model.get_gN() << "\n";
   }

   const int exit_code = spectrum_generator.get_exit_code();

   return exit_code;
}
