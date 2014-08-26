// ====================================================================
// Class for calculating the fine tuning in the MSSM. Provides
// generic numerical routines as well as pMSSM-specific 
// routines based on an analytical calculation.
// Notes:
//   - as of 26/8/2014, we currently only calculate fine tuning
//     neglecting family mixing (so stop matrix is simple 2x2)
// ====================================================================

#ifndef lowMSSM_TUNING_CALCULATOR_H
#define lowMSSM_TUNING_CALCULATOR_H

#include "lowMSSM_two_scale_model.hpp"
#include "lowMSSM_utilities.hpp"

#include "error.hpp"
#include "numerics.hpp"
#include <Eigen/Core>

namespace flexiblesusy {

class lowMSSM_tuning_calculator {
public:
   lowMSSM_tuning_calculator()
      : model()
      , input_scale(0.)
      , tuning_scale(0.)
      , precision_goal(1.0e-4)
      , max_iterations(0)
      , beta_loop_order(2)
      , threshold_corrections_loop_order(1) {}
   ~lowMSSM_tuning_calculator() {}

   double get_tuning_scale() const { return tuning_scale; }
   double get_input_scale()  const { return input_scale;  }
   const lowMSSM<Two_scale>& get_model() const { return model; }
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
   lowMSSM<Two_scale> model;
   double input_scale, tuning_scale;
   double precision_goal; ///< precision goal
   unsigned max_iterations; ///< maximum number of iterations
   unsigned beta_loop_order; ///< beta-function loop order
   unsigned threshold_corrections_loop_order; ///< threshold corrections loop order
   /// For now, while just using stops, save the result for later use
   Eigen::Array<double,2,1> MStop; 

   /// Helper methods in analytic tuning expressions.
   double gbar() const; 
   double MFtop_DRbar() const;
   double stop_mass_matrix_LL_entry() const; ///< note LL = (0,0) entry
   double stop_mass_matrix_RR_entry() const; ///< note RR = (1,1) entry
   double stop_mass_matrix_LR_entry() const; ///< note LR = (0,1) entry
   double MQQ2() const; 
   double RQQ() const;  
   double stop_discriminant() const;   

   void calculate_MStop();

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
