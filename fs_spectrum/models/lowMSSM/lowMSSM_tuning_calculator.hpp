// ====================================================================
// Class for calculating the fine tuning in the MSSM. Provides
// generic numerical routines as well as pMSSM-specific 
// routines based on an analytical calculation.
// ====================================================================

#ifndef lowMSSM_TUNING_CALCULATOR_H
#define lowMSSM_TUNING_CALCULATOR_H

#include "lowMSSM_two_scale_model.hpp"
#include "lowMSSM_utilities.hpp"

#include "error.hpp"
#include "numerics.hpp"

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

   /// Helper methods in analytic tuning expressions.
   double gbar() const; ///< for convenience
   double MQQ2() const; ///< for convenience
   double RQQ() const;  ///< for convenience
   double rt() const; ///< for convenience
   /// DH:: Note a0 has OPPOSITE sign convention to that used in cE6SSM paper,
   /// and therefore to that used in our expressions. Also A0 takes as input
   /// the mass, NOT the mass squared, and is evaluated at the current scale
   double deriv_d2DeltaV_dvd_dvd() const;
   double deriv_d2DeltaV_dvu_dvu() const;
   double deriv_d2DeltaV_dvu_dvd() const;

   /// Unit tests, to be removed later.

};

} // namespace flexiblesusy

#endif
