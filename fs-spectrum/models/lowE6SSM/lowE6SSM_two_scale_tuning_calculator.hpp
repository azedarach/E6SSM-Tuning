// ====================================================================
// Class for calculating the fine tuning in an E6 model.
// Notes:
//   - as of 26/8/2014, we currently only calculate fine tuning
//     neglecting family mixing (so stop matrix is simple 2x2)
//   - to fit in with our previous work, the fine tuning labelled
//     TLambdax is really that associated with the SUGRA parameter
//     ALambdax.
// ====================================================================

#ifndef lowE6SSM_TUNING_CALCULATOR_H
#define lowE6SSM_TUNING_CALCULATOR_H

#include "lowE6SSM_two_scale_ew_derivs.hpp"
#include "lowE6SSM_two_scale_model.hpp"
#include "lowE6SSM_utilities.hpp"

#include "error.hpp"
#include "numerics.hpp"
#include <Eigen/Core>
#include <Eigen/LU>
#include <map>

namespace flexiblesusy { 

   class lowE6SSM_tuning_calculator {
   public:
      explicit lowE6SSM_tuning_calculator(const lowE6SSM<Two_scale>&);
      ~lowE6SSM_tuning_calculator();

      double get_tuning_scale() const { return tuning_scale; }
      double get_input_scale()  const { return input_scale;  }
      const Problems<lowE6SSM_info::NUMBER_OF_PARTICLES> &get_problems() const {
         return model.get_problems();
      }
      const std::map<lowE6SSM_info::Parameters,double>& get_fine_tunings() const
         {
            return fine_tunings;
         }
      std::vector<lowE6SSM_info::Parameters> get_fine_tuning_params() const;
      void set_model(const lowE6SSM<Two_scale>& m) { model = m; }
      void set_input_scale(double s) { input_scale = s; }
      void set_tuning_scale(double s) { tuning_scale = s; }
      void set_tolerance(double t) { tolerance = t; }
      void set_tuning_ewsb_loop_order(unsigned l) { tuning_ewsb_loop_order = l; model.set_ewsb_loop_order(l); }
      void set_tuning_beta_loop_order(unsigned l) { tuning_beta_loop_order = l; model.set_loops(l); }
      void set_max_iterations(unsigned n) { max_iterations = n; }

      /// Calculate the fine tunings. Returns true if there is a problem.
      bool calculate_fine_tunings_numerically();
      bool calculate_fine_tunings_approximately();

   private:
      static const std::size_t num_ewsb_eqs = 3;
      static const std::size_t num_tree_level_ewsb_pars = 8;
      static const std::size_t num_one_loop_ewsb_pars = 4;
      static const std::size_t num_tuning_pars = 12;

      const double underflow = 1.0e-100;

      struct numerical_deriv_pars {
         lowE6SSM<Two_scale> model;
         double high_scale;
         double low_scale;
         unsigned ewsb_loop_order;
         unsigned beta_loop_order;
         lowE6SSM_info::Parameters p;
      };

      lowE6SSM<Two_scale> model;
      double input_scale; ///< parameter input scale
      double tuning_scale; ///< scale to calculate tuning at
      double tolerance; ///< tolerance in numerical derivatives
      std::size_t max_iterations; ///< maximum number of iterations
      unsigned tuning_ewsb_loop_order; ///< order of CW loop corrections included in tuning calculation (<= 1)
      unsigned tuning_beta_loop_order; ///< order of beta functions included in tuning calculation (<= 2)
      std::map<lowE6SSM_info::Parameters,double> fine_tunings;
      std::map<lowE6SSM_info::Parameters,double> input_scale_pars;

      static double calculate_MVZ2(double x, void * params);

      void get_input_scale_pars();
      double calculate_fine_tuning(double par, double g1, double g2, double vd, double vu, 
                                   double dg1dp, double dg2dp, double dvddp, double dvudp) const;
      std::size_t get_g1_row() const;
      std::size_t get_g2_row() const;
      
      double get_deriv(lowE6SSM_info::Parameters p, const Eigen::Array<double,1,num_tuning_pars>& derivs) const;

      Eigen::Matrix<double,num_ewsb_eqs,Eigen::Dynamic> calculate_ewsb_parameter_derivs(const lowE6SSM_ew_derivs&) const;
      Eigen::Matrix<double,Eigen::Dynamic,num_tuning_pars> calculate_beta_derivs() const;

      Eigen::Matrix<double,Eigen::Dynamic,1> calculate_deriv_dlowscale_dLambdax() const;
      Eigen::Matrix<double,Eigen::Dynamic,1> calculate_deriv_dlowscale_dTLambdax() const;
      Eigen::Matrix<double,Eigen::Dynamic,1> calculate_deriv_dlowscale_dTYu22() const;
      Eigen::Matrix<double,Eigen::Dynamic,1> calculate_deriv_dlowscale_dmq222() const;
      Eigen::Matrix<double,Eigen::Dynamic,1> calculate_deriv_dlowscale_dmHd2() const;
      Eigen::Matrix<double,Eigen::Dynamic,1> calculate_deriv_dlowscale_dmHu2() const;
      Eigen::Matrix<double,Eigen::Dynamic,1> calculate_deriv_dlowscale_dmu222() const;
      Eigen::Matrix<double,Eigen::Dynamic,1> calculate_deriv_dlowscale_dms2() const;
      Eigen::Matrix<double,Eigen::Dynamic,1> calculate_deriv_dlowscale_dMassB() const;
      Eigen::Matrix<double,Eigen::Dynamic,1> calculate_deriv_dlowscale_dMassWB() const;
      Eigen::Matrix<double,Eigen::Dynamic,1> calculate_deriv_dlowscale_dMassG() const;
      Eigen::Matrix<double,Eigen::Dynamic,1> calculate_deriv_dlowscale_dMassBp() const;

   };

} // namespace flexiblesusy

#endif
