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
#include <map>

namespace flexiblesusy {
   
   template <class T>
   class lowMSSM_tuning_calculator {
   public:
      lowMSSM_tuning_calculator()
         : model(0)
         , input_scale(0.)
         , tuning_scale(0.)
         , precision_goal(1.0e-4)
         , max_iterations(0)
         , tuning_loop_order(1)
         , num_tuning_pars(0) {}
      ~lowMSSM_tuning_calculator() {}

      double get_tuning_scale() const { return tuning_scale; }
      double get_input_scale()  const { return input_scale;  }
      const Problems<lowMSSM_info::NUMBER_OF_PARTICLES> &get_problems() const {
         return model->get_problems();
      }
      const std::map<lowMSSM_info::Parameters,double>& get_fine_tunings() const
         {
            return fine_tunings;
         }
      void set_model(lowMSSM<T>* m) { model = m; }
      void set_input_scale(double s) { input_scale = s; }
      void set_tuning_scale(double s) { tuning_scale = s; }
      void set_tolerance(double t) { tolerance = t; }
      void set_tuning_loop_order(unsigned l) { tuning_loop_order = l; }
      void set_max_iterations(unsigned n) { max_iterations = n; }
      void add_fine_tuning_parameter(lowMSSM_info::Parameters p);

      void calculate_fine_tunings_numerically();
      void calculate_fine_tunings_approximately();

   private:

      lowMSSM<T>* model;
      double input_scale; ///< parameter input scale
      double tuning_scale; ///< scale to calculate tuning at
      double tolerance; ///< tolerance in numerical derivatives
      std::size_t max_iterations; ///< maximum number of iterations
      unsigned tuning_loop_order; ///< order of CW loop corrections included in tuning calculation (<= 1)
      std::size_t num_tuning_pars;
      std::map<lowMSSM_info::Parameters,double> fine_tunings;
   };

   template <class T>
   void lowMSSM_tuning_calculator<T>::add_fine_tuning_parameter(lowMSSM_info::Parameters p)
   {
      std::pair<std::map<lowMSSM_info::Parameters,double>,bool> result;
      result = fine_tunings.insert(std::pair<lowMSSM_info::Parameters,double>(p, 0.0));

      if (result.second == true) {
         ++num_tuning_pars;
      }
   }

   template <class T>
   void lowMSSM_tuning_calculator<T>::calculate_fine_tunings_numerically()
   {

   }

   template <class T>
   void lowMSSM_tuning_calculator<T>::calculate_fine_tunings_approximately()
   {
      lowMSSM_ewsb_conditions ewsb_cond;
      // Get unrotated CP-even Higgs mass matrix at given
      // loop order
      if (model->get_scale() == tuning_scale) {
         ewsb_cond = lowMSSM_ewsb_conditions(cast_model<lowMSSM<Two_scale> >(*model));
         model->run_to(input_scale);

      } else if (model->get_scale() == input_scale) {

      } else {

      }

      Eigen::Matrix<double,2,2> mass_matrix_hh = ewsb_cond.calculate_unrotated_mass_matrix_hh();


   }

} // namespace flexiblesusy

#endif
