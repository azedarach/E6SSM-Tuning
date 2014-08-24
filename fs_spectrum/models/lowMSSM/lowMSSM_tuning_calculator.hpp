// ====================================================================
// Class for calculating the fine tuning in the MSSM. Provides
// generic numerical routines as well as pMSSM-specific 
// routines based on an analytical calculation.
// ====================================================================

#ifndef lowMSSM_TUNING_CALCULATOR_H
#define lowMSSM_TUNING_CALCULATOR_H

#include "lowMSSM_two_scale_model.hpp"
#include "lowMSSM_two_scale_susy_scale_constraint.hpp"
#include "lowMSSM_two_scale_low_scale_constraint.hpp"
#include "lowMSSM_two_scale_convergence_tester.hpp"
#include "lowMSSM_two_scale_initial_guesser.hpp"
#include "lowMSSM_utilities.hpp"

#include "coupling_monitor.hpp"
#include "error.hpp"
#include "numerics.hpp"
#include "two_scale_running_precision.hpp"
#include "two_scale_solver.hpp"

namespace flexiblesusy {

template <class T>
class lowMSSM_tuning_calculator {
public:
   lowMSSM_tuning_calculator()
      : model()
      , input_scale(0.)
      , tuning_scale(0.)
      , precision_goal(1.0e-4)
      , max_iterations(0)
      , beta_loop_order(2)
      , threshold_corrections_loop_order(1)
   ~lowMSSM_spectrum_generator() {}

   double get_tuning_scale() const { return tuning_scale; }
   double get_input_scale()  const { return input_scale;  }
   const lowMSSM<T>& get_model() const { return model; }
   const Problems<lowMSSM_info::NUMBER_OF_PARTICLES>& get_problems() const {
      return model.get_problems();
   }
   int get_exit_code() const { return get_problems().have_serious_problem(); };
   void set_input_scale(double s) { input_scale = s; }
   void set_tuning_scale(double s) { tuning_scale = s; }
   void set_precision_goal(double precision_goal_) { precision_goal = precision_goal_; }
   void set_ewsb_loop_order(unsigned l) { model.set_ewsb_loop_order(l); }
   void set_beta_loop_order(unsigned l) { beta_loop_order = l; }
   void set_max_iterations(unsigned n) { max_iterations = n; }
   void set_threshold_corrections_loop_order(unsigned t) { threshold_corrections_loop_order = t; }

private:
   lowMSSM<T> model;
   double input_scale, tuning_scale;
   double precision_goal; ///< precision goal
   unsigned max_iterations; ///< maximum number of iterations
   unsigned beta_loop_order; ///< beta-function loop order
   unsigned threshold_corrections_loop_order; ///< threshold corrections loop order

  /// Useful helper methods in analytic tuning expressions

};

} // namespace flexiblesusy

#endif
