#include "lowE6SSM_two_scale_tuning_calculator.hpp"
#include "numerics.hpp"

#include <gsl/gsl_deriv.h>

namespace flexiblesusy {

   lowE6SSM_tuning_calculator::lowE6SSM_tuning_calculator(const lowE6SSM<Two_scale>& m)
         : model(m)
         , input_scale(0.)
         , tuning_scale(0.)
         , tolerance(1.0e-4)
         , max_iterations(100)
   {
      tuning_ewsb_loop_order = m.get_ewsb_loop_order();
      tuning_beta_loop_order = m.get_loops();

      fine_tunings[lowE6SSM_info::Lambdax] = 0.0;
      fine_tunings[lowE6SSM_info::TLambdax] = 0.0;
      fine_tunings[lowE6SSM_info::TYu22] = 0.0;
      fine_tunings[lowE6SSM_info::mq222] = 0.0;
      fine_tunings[lowE6SSM_info::mHd2] = 0.0;
      fine_tunings[lowE6SSM_info::mHu2] = 0.0;
      fine_tunings[lowE6SSM_info::mu222] = 0.0;
      fine_tunings[lowE6SSM_info::ms2] = 0.0;
      fine_tunings[lowE6SSM_info::MassB] = 0.0;
      fine_tunings[lowE6SSM_info::MassWB] = 0.0;
      fine_tunings[lowE6SSM_info::MassG] = 0.0;
      fine_tunings[lowE6SSM_info::MassBp] = 0.0;

      input_scale_pars[lowE6SSM_info::Lambdax] = 0.0;
      input_scale_pars[lowE6SSM_info::TLambdax] = 0.0;
      input_scale_pars[lowE6SSM_info::TYu22] = 0.0;
      input_scale_pars[lowE6SSM_info::mq222] = 0.0;
      input_scale_pars[lowE6SSM_info::mHd2] = 0.0;
      input_scale_pars[lowE6SSM_info::mHu2] = 0.0;
      input_scale_pars[lowE6SSM_info::mu222] = 0.0;
      input_scale_pars[lowE6SSM_info::ms2] = 0.0;
      input_scale_pars[lowE6SSM_info::MassB] = 0.0;
      input_scale_pars[lowE6SSM_info::MassWB] = 0.0;
      input_scale_pars[lowE6SSM_info::MassG] = 0.0;
      input_scale_pars[lowE6SSM_info::MassBp] = 0.0;
   }

   lowE6SSM_tuning_calculator::~lowE6SSM_tuning_calculator()
   {

   }

   double lowE6SSM_tuning_calculator::calculate_MVZ2(double x, void * params)
   {
      const numerical_deriv_pars* pars
         = static_cast<numerical_deriv_pars*>(params);

      lowE6SSM<Two_scale> model = pars->model;
      double high_scale = pars->high_scale;
      double low_scale = pars->low_scale;
      unsigned ewsb_loop_order = pars->ewsb_loop_order;
      unsigned beta_loop_order = pars->beta_loop_order;
      lowE6SSM_info::Parameters p = pars->p;

      model.set_loops(beta_loop_order);
      model.set_ewsb_loop_order(ewsb_loop_order);

      if (beta_loop_order == 0) {
         model.set_scale(high_scale);
      } else if (!is_equal_rel(high_scale, model.get_scale())) {
         model.run_to(high_scale);
      }

      if (p == lowE6SSM_info::Lambdax) {
         const double underflow = 1.0e-100;
         double Alambda;
         if (is_zero(model.get_TLambdax())) {
            Alambda = 0.0;
         } else if (Abs(model.get_Lambdax()) < underflow) {
            throw DivideByZeroError("in lowE6SSM_tuning_calculator::get_input_scale_pars");
         } else {
            Alambda = model.get_TLambdax() / model.get_Lambdax();
         }
         model.set_parameter(lowE6SSM_info::TLambdax, x * Alambda);
         model.set_parameter(p, x);
      } else if (p == lowE6SSM_info::TLambdax) {
         model.set_parameter(lowE6SSM_info::TLambdax, model.get_Lambdax() * x);
      } else if (p == lowE6SSM_info::TYu22) {
         model.set_parameter(lowE6SSM_info::TYu22, model.get_Yu(2, 2) * x);
      } else {
         model.set_parameter(p, x);
      }
      
      if (beta_loop_order == 0) {
         model.set_scale(low_scale);
      } else {
         model.run_to(low_scale);
      }

      const double g1 = model.get_g1();
      const double g2 = model.get_g2();
      const double gbar2 = Sqr(g2) + 0.6 * Sqr(g1);

      const double mHd2 = model.get_mHd2();
      const double mHu2 = model.get_mHu2();
      const double ms2 = model.get_ms2();

      lowE6SSM_ew_derivs derivs(model);

      derivs.get_model().set_mHd2(mHd2);
      derivs.get_model().set_mHu2(mHu2);
      derivs.get_model().set_ms2(ms2);

      derivs.solve_ewsb_conditions_for_vevs();

      const double vd = derivs.get_model().get_vd();
      const double vu = derivs.get_model().get_vu();
      const double vev2 = Sqr(vu) + Sqr(vd);

      double MVZ2 = 0.25 * gbar2 * vev2;
      
      return MVZ2;
   }

   bool lowE6SSM_tuning_calculator::calculate_fine_tunings_numerically()
   {
      bool tuning_problem = false;

      lowE6SSM_ew_derivs ew_derivs(model);
      
      double scale = model.get_scale();

      if (is_equal_rel(scale, input_scale)) {
         get_input_scale_pars();
      } else {
         if (tuning_beta_loop_order == 0) {
            model.set_scale(input_scale);
         } else {
            model.run_to(input_scale);
         }
         get_input_scale_pars();
      }

      const double epsilon = 1.0e-5;

      // Reference value of M_Z^2
      const double g1 = model.get_g1();
      const double g2 = model.get_g2();
      const double gbar2 = Sqr(g2) + 0.6 * Sqr(g1);
      const double vd = model.get_vd();
      const double vu = model.get_vu();
      const double vev2 = Sqr(vu) + Sqr(vd);
      const double MVZ2 = 0.25 * gbar2 * vev2;

      for (std::map<lowE6SSM_info::Parameters,double>::const_iterator it = input_scale_pars.begin(),
              end = input_scale_pars.end(); it != end; ++it) {
         numerical_deriv_pars pars = {model, input_scale, tuning_scale, 
                                      tuning_ewsb_loop_order, tuning_beta_loop_order, it->first};
         gsl_function calc_MVZ2 = {calculate_MVZ2, &pars};

         double numeric_deriv;
         double abs_err;
         double h;
         
         h = epsilon * Abs(it->second);
         if (h == 0.0) h = epsilon;
         double tmp_1 = it->second;
         double tmp_2 = tmp_1 + h;
         h = tmp_2 - tmp_1;
         
         gsl_deriv_central(&calc_MVZ2, it->second, h, &numeric_deriv, &abs_err);

         if (Abs(abs_err / numeric_deriv) > 1.0 && Abs(numeric_deriv) > 1.0e-10)
            tuning_problem = true;

         fine_tunings[it->first] = Abs((it->second / MVZ2) * numeric_deriv);
      }

      return tuning_problem;
   }

   std::vector<lowE6SSM_info::Parameters> lowE6SSM_tuning_calculator::get_fine_tuning_params() const
   {
      std::vector<lowE6SSM_info::Parameters> params;
      for (std::map<lowE6SSM_info::Parameters,double>::const_iterator it = fine_tunings.begin(),
              end = fine_tunings.end(); it != end; ++it) {
         params.push_back(it->first);
      }
      return params;
   }
   
   bool lowE6SSM_tuning_calculator::calculate_fine_tunings_approximately()
   {
      bool tuning_problem = false;
      
      lowE6SSM_ew_derivs ew_derivs(model);
      Eigen::Matrix<double,Eigen::Dynamic,num_tuning_pars> beta_derivs;
      
      double scale = model.get_scale();

      if (tuning_beta_loop_order == 0) {
         model.set_scale(input_scale);
         beta_derivs = calculate_beta_derivs();
         get_input_scale_pars();
         model.set_scale(tuning_scale);
         ew_derivs.set_model(model);
      } else {
         if (is_equal_rel(scale, tuning_scale)) {
            ew_derivs.set_model(model);
            model.run_to(input_scale);
            beta_derivs = calculate_beta_derivs();
            get_input_scale_pars();
         } else if (is_equal_rel(scale, input_scale)) {
            beta_derivs = calculate_beta_derivs();
            get_input_scale_pars();
            model.run_to(tuning_scale);
            ew_derivs.set_model(model);
         } else {
            model.run_to(input_scale);
            beta_derivs = calculate_beta_derivs();
            get_input_scale_pars();
            model.run_to(tuning_scale);
            ew_derivs.set_model(model);
         }
      }

      ew_derivs.set_ewsb_loop_order(tuning_ewsb_loop_order);
      
      Eigen::Matrix<double,num_ewsb_eqs,num_ewsb_eqs> mass_matrix = ew_derivs.calculate_unrotated_mass_matrix_hh();
      
      // N.B. extra minus sign
      Eigen::Matrix<double, num_ewsb_eqs, Eigen::Dynamic> ewsb_derivs = -1.0 * calculate_ewsb_parameter_derivs(ew_derivs);
      
      // Solve system. Use inverse directly since mass matrix is small in the E6SSM (3 x 3),
      // but note this has to change for matrices larger than 4 x 4.
      bool invertible;
      Eigen::Matrix<double,num_ewsb_eqs,num_ewsb_eqs> mass_matrix_inverse;
      mass_matrix.computeInverseWithCheck(mass_matrix_inverse, invertible);

      if (!invertible) {
         tuning_problem = true;
         return tuning_problem;
      }

      Eigen::Matrix<double,num_ewsb_eqs,num_tuning_pars> vev_derivs
         = mass_matrix_inverse * ewsb_derivs * beta_derivs;

      // Additional solution check here?

      double g1_at_tuning_scale = ew_derivs.get_model().get_g1();
      double g2_at_tuning_scale = ew_derivs.get_model().get_g2();
      double vd_at_tuning_scale = ew_derivs.get_model().get_vd();
      double vu_at_tuning_scale = ew_derivs.get_model().get_vu();
      const Eigen::Array<double,1,num_tuning_pars> g1_derivs = beta_derivs.row(get_g1_row());
      const Eigen::Array<double,1,num_tuning_pars> g2_derivs = beta_derivs.row(get_g2_row());
      const Eigen::Array<double,1,num_tuning_pars> vd_derivs = vev_derivs.row(0);
      const Eigen::Array<double,1,num_tuning_pars> vu_derivs = vev_derivs.row(1);

      for (std::map<lowE6SSM_info::Parameters,double>::iterator it = fine_tunings.begin(),
              end = fine_tunings.end(); it != end; ++it) {
         it->second = calculate_fine_tuning(input_scale_pars[it->first], g1_at_tuning_scale,
                                            g2_at_tuning_scale, vd_at_tuning_scale, vu_at_tuning_scale,
                                            get_deriv(it->first, g1_derivs), get_deriv(it->first, g2_derivs),
                                            get_deriv(it->first, vd_derivs), get_deriv(it->first, vu_derivs));
      }

      return tuning_problem;
   }

   void lowE6SSM_tuning_calculator::get_input_scale_pars()
   {
      input_scale_pars[lowE6SSM_info::Lambdax] = model.get_Lambdax();

      double Alambda;

      if (is_zero(model.get_TLambdax())) {
         Alambda = 0.0;
      } else if (Abs(model.get_Lambdax()) < underflow) {
         throw DivideByZeroError("in lowE6SSM_tuning_calculator::get_input_scale_pars");
      } else {
         Alambda = model.get_TLambdax() / model.get_Lambdax();
      }

      input_scale_pars[lowE6SSM_info::TLambdax] = Alambda;

      double At;
      if (is_zero(model.get_TYu(2, 2))) {
         At = 0.0;
      } else if (Abs(model.get_Yu(2, 2)) < underflow) {
         throw DivideByZeroError("in lowE6SSM_tuning_calculator::get_input_scale_pars");
      } else {
         At = model.get_TYu(2, 2) / model.get_Yu(2, 2);
      }

      input_scale_pars[lowE6SSM_info::TYu22] = At;
      input_scale_pars[lowE6SSM_info::mq222] = model.get_mq2(2, 2);
      input_scale_pars[lowE6SSM_info::mHd2] = model.get_mHd2();
      input_scale_pars[lowE6SSM_info::mHu2] = model.get_mHu2();
      input_scale_pars[lowE6SSM_info::mu222] = model.get_mu2(2, 2);
      input_scale_pars[lowE6SSM_info::ms2] = model.get_ms2();
      input_scale_pars[lowE6SSM_info::MassB] = model.get_MassB();
      input_scale_pars[lowE6SSM_info::MassWB] = model.get_MassWB();
      input_scale_pars[lowE6SSM_info::MassG] = model.get_MassG();
      input_scale_pars[lowE6SSM_info::MassBp] = model.get_MassBp();
   }
   
   double lowE6SSM_tuning_calculator::calculate_fine_tuning(double par, double g1, double g2, double vd, double vu, 
                                                            double dg1dp, double dg2dp, double dvddp, double dvudp) const
   {
      double gbar2 = Sqr(g2) + 0.6 * Sqr(g1);
      double v2 = Sqr(vd) + Sqr(vu);

      // If either of gbar2 or v2 are zero then our problems are worse than just
      // dividing by zero...
      if (is_zero(gbar2) || is_zero(v2)) {
         throw DivideByZeroError("in lowE6SSM_tuning_calculator::calculate_fine_tuning");
      } else {
         double tuning = 2.0 * par * (vd * dvddp + vu * dvudp) / v2;
         tuning += (2.0 * par * (0.6 * g1 * dg1dp + g2 * dg2dp) / gbar2);
         return Abs(tuning); 
      }
   }

   std::size_t lowE6SSM_tuning_calculator::get_g1_row() const
   {
      if (tuning_beta_loop_order == 0) {
         return 1;
      } else {
         return 2;
      }
   }

   std::size_t lowE6SSM_tuning_calculator::get_g2_row() const
   {
      if (tuning_beta_loop_order == 0) {
         return 2;
      } else {
         return 3;
      }
   }

   double lowE6SSM_tuning_calculator::get_deriv(lowE6SSM_info::Parameters p, const Eigen::Array<double,1,num_tuning_pars>& derivs) const
   {
      switch (p) {
      case lowE6SSM_info::Lambdax: {
         return derivs(0);
      }
      case lowE6SSM_info::TLambdax: {
         return derivs(1);
      }
      case lowE6SSM_info::TYu22: {
         return derivs(2);
      }
      case lowE6SSM_info::mq222: {
         return derivs(3);
      }
      case lowE6SSM_info::mHd2: {
         return derivs(4);
      }
      case lowE6SSM_info::mHu2: {
         return derivs(5);
      }
      case lowE6SSM_info::mu222: {
         return derivs(6);
      } 
      case lowE6SSM_info::ms2: {
         return derivs(7);
      }
      case lowE6SSM_info::MassB: {
         return derivs(8);
      }
      case lowE6SSM_info::MassWB: {
         return derivs(9);
      }
      case lowE6SSM_info::MassG: {
         return derivs(10);
      }
      case lowE6SSM_info::MassBp: {
         return derivs(11);
      }
      default: {
         return 0.0;
      }
      }
   }

   Eigen::Matrix<double,lowE6SSM_tuning_calculator::num_ewsb_eqs,Eigen::Dynamic> lowE6SSM_tuning_calculator::calculate_ewsb_parameter_derivs(const lowE6SSM_ew_derivs & derivs) const
   {

      lowE6SSM_ew_derivs ew_derivs(derivs);

      if (tuning_ewsb_loop_order == 0) {

         Eigen::Matrix<double,num_ewsb_eqs,num_tree_level_ewsb_pars> derivs;

         derivs.col(0) = ew_derivs.calculate_ewsb_parameter_derivs(lowE6SSM_info::Lambdax);
         derivs.col(1) = ew_derivs.calculate_ewsb_parameter_derivs(lowE6SSM_info::g1);
         derivs.col(2) = ew_derivs.calculate_ewsb_parameter_derivs(lowE6SSM_info::g2);
         derivs.col(3) = ew_derivs.calculate_ewsb_parameter_derivs(lowE6SSM_info::gN);
         derivs.col(4) = ew_derivs.calculate_ewsb_parameter_derivs(lowE6SSM_info::TLambdax);
         derivs.col(5) = ew_derivs.calculate_ewsb_parameter_derivs(lowE6SSM_info::mHd2);
         derivs.col(6) = ew_derivs.calculate_ewsb_parameter_derivs(lowE6SSM_info::mHu2);
         derivs.col(7) = ew_derivs.calculate_ewsb_parameter_derivs(lowE6SSM_info::ms2);

         return derivs;

      } else {

         Eigen::Matrix<double,num_ewsb_eqs,num_tree_level_ewsb_pars + num_one_loop_ewsb_pars> derivs;
         
         derivs.col(0) = ew_derivs.calculate_ewsb_parameter_derivs(lowE6SSM_info::Lambdax);
         derivs.col(1) = ew_derivs.calculate_ewsb_parameter_derivs(lowE6SSM_info::Yu22);
         derivs.col(2) = ew_derivs.calculate_ewsb_parameter_derivs(lowE6SSM_info::g1);
         derivs.col(3) = ew_derivs.calculate_ewsb_parameter_derivs(lowE6SSM_info::g2);
         derivs.col(4) = ew_derivs.calculate_ewsb_parameter_derivs(lowE6SSM_info::gN);
         derivs.col(5) = ew_derivs.calculate_ewsb_parameter_derivs(lowE6SSM_info::TLambdax);
         derivs.col(6) = ew_derivs.calculate_ewsb_parameter_derivs(lowE6SSM_info::TYu22);
         derivs.col(7) = ew_derivs.calculate_ewsb_parameter_derivs(lowE6SSM_info::mq222);
         derivs.col(8) = ew_derivs.calculate_ewsb_parameter_derivs(lowE6SSM_info::mHd2);
         derivs.col(9) = ew_derivs.calculate_ewsb_parameter_derivs(lowE6SSM_info::mHu2);
         derivs.col(10) = ew_derivs.calculate_ewsb_parameter_derivs(lowE6SSM_info::mu222);
         derivs.col(11) = ew_derivs.calculate_ewsb_parameter_derivs(lowE6SSM_info::ms2);

         return derivs;
      }
   }

   Eigen::Matrix<double,Eigen::Dynamic,lowE6SSM_tuning_calculator::num_tuning_pars> lowE6SSM_tuning_calculator::calculate_beta_derivs() const
   {
      Eigen::Matrix<double,Eigen::Dynamic,num_tuning_pars> derivs;
      if (tuning_ewsb_loop_order == 0) {
         derivs.resize(num_tree_level_ewsb_pars, num_tuning_pars);
      } else {
         derivs.resize(num_tree_level_ewsb_pars + num_one_loop_ewsb_pars, num_tuning_pars);
      }
      
      derivs.col(0) = calculate_deriv_dlowscale_dLambdax();
      derivs.col(1) = calculate_deriv_dlowscale_dTLambdax();
      derivs.col(2) = calculate_deriv_dlowscale_dTYu22();
      derivs.col(3) = calculate_deriv_dlowscale_dmq222();
      derivs.col(4) = calculate_deriv_dlowscale_dmHd2();
      derivs.col(5) = calculate_deriv_dlowscale_dmHu2();
      derivs.col(6) = calculate_deriv_dlowscale_dmu222();
      derivs.col(7) = calculate_deriv_dlowscale_dms2();
      derivs.col(8) = calculate_deriv_dlowscale_dMassB();
      derivs.col(9) = calculate_deriv_dlowscale_dMassWB();
      derivs.col(10) = calculate_deriv_dlowscale_dMassG();
      derivs.col(11) = calculate_deriv_dlowscale_dMassBp();

      return derivs;
   }

   Eigen::Matrix<double,Eigen::Dynamic,1> lowE6SSM_tuning_calculator::calculate_deriv_dlowscale_dLambdax() const
   {
      if (tuning_ewsb_loop_order == 0) {
         Eigen::Matrix<double,num_tree_level_ewsb_pars,1> derivs;
         if (tuning_beta_loop_order == 0) {
            derivs(0) = 1.;
            derivs(1) = 0.;
            derivs(2) = 0.;
            derivs(3) = 0.;
            double Alambda;

            if (is_zero(model.get_TLambdax())) {
               Alambda = 0.0;
            } else if (Abs(model.get_Lambdax()) < underflow) {
               throw DivideByZeroError("in lowE6SSM_tuning_calculator::get_input_scale_pars");
            } else {
               Alambda = model.get_TLambdax() / model.get_Lambdax();
            }

            derivs(4) = Alambda;
            derivs(5) = 0.;
            derivs(6) = 0.;
            derivs(7) = 0.;
         } else {

         }
         return derivs;
      } else {
         Eigen::Matrix<double,num_tree_level_ewsb_pars + num_one_loop_ewsb_pars,1> derivs;
         if (tuning_beta_loop_order == 0) {
            derivs(0) = 1.;
            derivs(1) = 0.;
            derivs(2) = 0.;
            derivs(3) = 0.;
            derivs(4) = 0.;
            double Alambda;

            if (is_zero(model.get_TLambdax())) {
               Alambda = 0.0;
            } else if (Abs(model.get_Lambdax()) < underflow) {
               throw DivideByZeroError("in lowE6SSM_tuning_calculator::get_input_scale_pars");
            } else {
               Alambda = model.get_TLambdax() / model.get_Lambdax();
            }

            derivs(5) = Alambda;
            derivs(6) = 0.;
            derivs(7) = 0.;
            derivs(8) = 0.;
            derivs(9) = 0.;
            derivs(10) = 0.;
            derivs(11) = 0.;
         } else {

         }
         return derivs;
      }
   }

   Eigen::Matrix<double,Eigen::Dynamic,1> lowE6SSM_tuning_calculator::calculate_deriv_dlowscale_dTLambdax() const
   {
      if (tuning_ewsb_loop_order == 0) {
         Eigen::Matrix<double,num_tree_level_ewsb_pars,1> derivs;
         if (tuning_beta_loop_order == 0) {
            derivs(0) = 0.;
            derivs(1) = 0.;
            derivs(2) = 0.;
            derivs(3) = 0.;
            derivs(4) = model.get_Lambdax();
            derivs(5) = 0.;
            derivs(6) = 0.;
            derivs(7) = 0.;
         } else {

         }
         return derivs;
      } else {
         Eigen::Matrix<double,num_tree_level_ewsb_pars + num_one_loop_ewsb_pars,1> derivs;
         if (tuning_beta_loop_order == 0) {
            derivs(0) = 0.;
            derivs(1) = 0.;
            derivs(2) = 0.;
            derivs(3) = 0.;
            derivs(4) = 0.;
            derivs(5) = model.get_Lambdax();
            derivs(6) = 0.;
            derivs(7) = 0.;
            derivs(8) = 0.;
            derivs(9) = 0.;
            derivs(10) = 0.;
            derivs(11) = 0.;
         } else {

         }
         return derivs;
      }
   }

   Eigen::Matrix<double,Eigen::Dynamic,1> lowE6SSM_tuning_calculator::calculate_deriv_dlowscale_dTYu22() const
   {
      if (tuning_ewsb_loop_order == 0) {
         Eigen::Matrix<double,num_tree_level_ewsb_pars,1> derivs;
         if (tuning_beta_loop_order == 0) {
            derivs(0) = 0.;
            derivs(1) = 0.;
            derivs(2) = 0.;
            derivs(3) = 0.;
            derivs(4) = 0.;
            derivs(5) = 0.;
            derivs(6) = 0.;
            derivs(7) = 0.;
         } else {

         }
         return derivs;
      } else {
         Eigen::Matrix<double,num_tree_level_ewsb_pars + num_one_loop_ewsb_pars,1> derivs;
         if (tuning_beta_loop_order == 0) {
            derivs(0) = 0.;
            derivs(1) = 0.;
            derivs(2) = 0.;
            derivs(3) = 0.;
            derivs(4) = 0.;
            derivs(5) = 0.;
            derivs(6) = model.get_Yu(2, 2);
            derivs(7) = 0.;
            derivs(8) = 0.;
            derivs(9) = 0.;
            derivs(10) = 0.;
            derivs(11) = 0.;
         } else {

         }
         return derivs;
      }
   }

   Eigen::Matrix<double,Eigen::Dynamic,1> lowE6SSM_tuning_calculator::calculate_deriv_dlowscale_dmq222() const
   {
      if (tuning_ewsb_loop_order == 0) {
         Eigen::Matrix<double,num_tree_level_ewsb_pars,1> derivs;
         if (tuning_beta_loop_order == 0) {
            derivs(0) = 0.;
            derivs(1) = 0.;
            derivs(2) = 0.;
            derivs(3) = 0.;
            derivs(4) = 0.;
            derivs(5) = 0.;
            derivs(6) = 0.;
            derivs(7) = 0.;
         } else {

         }
         return derivs;
      } else {
         Eigen::Matrix<double,num_tree_level_ewsb_pars + num_one_loop_ewsb_pars,1> derivs;
         if (tuning_beta_loop_order == 0) {
            derivs(0) = 0.;
            derivs(1) = 0.;
            derivs(2) = 0.;
            derivs(3) = 0.;
            derivs(4) = 0.;
            derivs(5) = 0.;
            derivs(6) = 0.;
            derivs(7) = 1.;
            derivs(8) = 0.;
            derivs(9) = 0.;
            derivs(10) = 0.;
            derivs(11) = 0.;
         } else {

         }
         return derivs;
      }
   }

   Eigen::Matrix<double,Eigen::Dynamic,1> lowE6SSM_tuning_calculator::calculate_deriv_dlowscale_dmHd2() const
   {
      if (tuning_ewsb_loop_order == 0) {
         Eigen::Matrix<double,num_tree_level_ewsb_pars,1> derivs;
         if (tuning_beta_loop_order == 0) {
            derivs(0) = 0.;
            derivs(1) = 0.;
            derivs(2) = 0.;
            derivs(3) = 0.;
            derivs(4) = 0.;
            derivs(5) = 1.;
            derivs(6) = 0.;
            derivs(7) = 0.;
         } else {

         }
         return derivs;
      } else {
         Eigen::Matrix<double,num_tree_level_ewsb_pars + num_one_loop_ewsb_pars,1> derivs;
         if (tuning_beta_loop_order == 0) {
            derivs(0) = 0.;
            derivs(1) = 0.;
            derivs(2) = 0.;
            derivs(3) = 0.;
            derivs(4) = 0.;
            derivs(5) = 0.;
            derivs(6) = 0.;
            derivs(7) = 0.;
            derivs(8) = 1.;
            derivs(9) = 0.;
            derivs(10) = 0.;
            derivs(11) = 0.;
         } else {

         }
         return derivs;
      }
   }

   Eigen::Matrix<double,Eigen::Dynamic,1> lowE6SSM_tuning_calculator::calculate_deriv_dlowscale_dmHu2() const
   {
      if (tuning_ewsb_loop_order == 0) {
         Eigen::Matrix<double,num_tree_level_ewsb_pars,1> derivs;
         if (tuning_beta_loop_order == 0) {
            derivs(0) = 0.;
            derivs(1) = 0.;
            derivs(2) = 0.;
            derivs(3) = 0.;
            derivs(4) = 0.;
            derivs(5) = 0.;
            derivs(6) = 1.;
            derivs(7) = 0.;
         } else {

         }
         return derivs;
      } else {
         Eigen::Matrix<double,num_tree_level_ewsb_pars + num_one_loop_ewsb_pars,1> derivs;
         if (tuning_beta_loop_order == 0) {
            derivs(0) = 0.;
            derivs(1) = 0.;
            derivs(2) = 0.;
            derivs(3) = 0.;
            derivs(4) = 0.;
            derivs(5) = 0.;
            derivs(6) = 0.;
            derivs(7) = 0.;
            derivs(8) = 0.;
            derivs(9) = 1.;
            derivs(10) = 0.;
            derivs(11) = 0.;
         } else {

         }
         return derivs;
      }
   }

   Eigen::Matrix<double,Eigen::Dynamic,1> lowE6SSM_tuning_calculator::calculate_deriv_dlowscale_dmu222() const
   {
      if (tuning_ewsb_loop_order == 0) {
         Eigen::Matrix<double,num_tree_level_ewsb_pars,1> derivs;
         if (tuning_beta_loop_order == 0) {
            derivs(0) = 0.;
            derivs(1) = 0.;
            derivs(2) = 0.;
            derivs(3) = 0.;
            derivs(4) = 0.;
            derivs(5) = 0.;
            derivs(6) = 0.;
            derivs(7) = 0.;
         } else {

         }
         return derivs;
      } else {
         Eigen::Matrix<double,num_tree_level_ewsb_pars + num_one_loop_ewsb_pars,1> derivs;
         if (tuning_beta_loop_order == 0) {
            derivs(0) = 0.;
            derivs(1) = 0.;
            derivs(2) = 0.;
            derivs(3) = 0.;
            derivs(4) = 0.;
            derivs(5) = 0.;
            derivs(6) = 0.;
            derivs(7) = 0.;
            derivs(8) = 0.;
            derivs(9) = 0.;
            derivs(10) = 1.;
            derivs(11) = 0.;
         } else {

         }
         return derivs;
      }
   }

   Eigen::Matrix<double,Eigen::Dynamic,1> lowE6SSM_tuning_calculator::calculate_deriv_dlowscale_dms2() const
   {
      if (tuning_ewsb_loop_order == 0) {
         Eigen::Matrix<double,num_tree_level_ewsb_pars,1> derivs;
         if (tuning_beta_loop_order == 0) {
            derivs(0) = 0.;
            derivs(1) = 0.;
            derivs(2) = 0.;
            derivs(3) = 0.;
            derivs(4) = 0.;
            derivs(5) = 0.;
            derivs(6) = 0.;
            derivs(7) = 1.;
         } else {

         }
         return derivs;
      } else {
         Eigen::Matrix<double,num_tree_level_ewsb_pars + num_one_loop_ewsb_pars,1> derivs;
         if (tuning_beta_loop_order == 0) {
            derivs(0) = 0.;
            derivs(1) = 0.;
            derivs(2) = 0.;
            derivs(3) = 0.;
            derivs(4) = 0.;
            derivs(5) = 0.;
            derivs(6) = 0.;
            derivs(7) = 0.;
            derivs(8) = 0.;
            derivs(9) = 0.;
            derivs(10) = 0.;
            derivs(11) = 1.;
         } else {

         }
         return derivs;
      }
   }

   Eigen::Matrix<double,Eigen::Dynamic,1> lowE6SSM_tuning_calculator::calculate_deriv_dlowscale_dMassB() const
   {
      if (tuning_ewsb_loop_order == 0) {
         Eigen::Matrix<double,num_tree_level_ewsb_pars,1> derivs;
         if (tuning_beta_loop_order == 0) {
            derivs(0) = 0.;
            derivs(1) = 0.;
            derivs(2) = 0.;
            derivs(3) = 0.;
            derivs(4) = 0.;
            derivs(5) = 0.;
            derivs(6) = 0.;
            derivs(7) = 0.;
         } else {

         }
         return derivs;
      } else {
         Eigen::Matrix<double,num_tree_level_ewsb_pars + num_one_loop_ewsb_pars,1> derivs;
         if (tuning_beta_loop_order == 0) {
            derivs(0) = 0.;
            derivs(1) = 0.;
            derivs(2) = 0.;
            derivs(3) = 0.;
            derivs(4) = 0.;
            derivs(5) = 0.;
            derivs(6) = 0.;
            derivs(7) = 0.;
            derivs(8) = 0.;
            derivs(9) = 0.;
            derivs(10) = 0.;
            derivs(11) = 0.;
         } else {

         }
         return derivs;
      }
   }

   Eigen::Matrix<double,Eigen::Dynamic,1> lowE6SSM_tuning_calculator::calculate_deriv_dlowscale_dMassWB() const
   {
      if (tuning_ewsb_loop_order == 0) {
         Eigen::Matrix<double,num_tree_level_ewsb_pars,1> derivs;
         if (tuning_beta_loop_order == 0) {
            derivs(0) = 0.;
            derivs(1) = 0.;
            derivs(2) = 0.;
            derivs(3) = 0.;
            derivs(4) = 0.;
            derivs(5) = 0.;
            derivs(6) = 0.;
            derivs(7) = 0.;
         } else {

         }
         return derivs;
      } else {
         Eigen::Matrix<double,num_tree_level_ewsb_pars + num_one_loop_ewsb_pars,1> derivs;
         if (tuning_beta_loop_order == 0) {
            derivs(0) = 0.;
            derivs(1) = 0.;
            derivs(2) = 0.;
            derivs(3) = 0.;
            derivs(4) = 0.;
            derivs(5) = 0.;
            derivs(6) = 0.;
            derivs(7) = 0.;
            derivs(8) = 0.;
            derivs(9) = 0.;
            derivs(10) = 0.;
            derivs(11) = 0.;
         } else {

         }
         return derivs;
      }
   }

   Eigen::Matrix<double,Eigen::Dynamic,1> lowE6SSM_tuning_calculator::calculate_deriv_dlowscale_dMassG() const
   {
      if (tuning_ewsb_loop_order == 0) {
         Eigen::Matrix<double,num_tree_level_ewsb_pars,1> derivs;
         if (tuning_beta_loop_order == 0) {
            derivs(0) = 0.;
            derivs(1) = 0.;
            derivs(2) = 0.;
            derivs(3) = 0.;
            derivs(4) = 0.;
            derivs(5) = 0.;
            derivs(6) = 0.;
            derivs(7) = 0.;
         } else {

         }
         return derivs;
      } else {
         Eigen::Matrix<double,num_tree_level_ewsb_pars + num_one_loop_ewsb_pars,1> derivs;
         if (tuning_beta_loop_order == 0) {
            derivs(0) = 0.;
            derivs(1) = 0.;
            derivs(2) = 0.;
            derivs(3) = 0.;
            derivs(4) = 0.;
            derivs(5) = 0.;
            derivs(6) = 0.;
            derivs(7) = 0.;
            derivs(8) = 0.;
            derivs(9) = 0.;
            derivs(10) = 0.;
            derivs(11) = 0.;
         } else {

         }
         return derivs;
      }
   }

   Eigen::Matrix<double,Eigen::Dynamic,1> lowE6SSM_tuning_calculator::calculate_deriv_dlowscale_dMassBp() const
   {
      if (tuning_ewsb_loop_order == 0) {
         Eigen::Matrix<double,num_tree_level_ewsb_pars,1> derivs;
         if (tuning_beta_loop_order == 0) {
            derivs(0) = 0.;
            derivs(1) = 0.;
            derivs(2) = 0.;
            derivs(3) = 0.;
            derivs(4) = 0.;
            derivs(5) = 0.;
            derivs(6) = 0.;
            derivs(7) = 0.;
         } else {

         }
         return derivs;
      } else {
         Eigen::Matrix<double,num_tree_level_ewsb_pars + num_one_loop_ewsb_pars,1> derivs;
         if (tuning_beta_loop_order == 0) {
            derivs(0) = 0.;
            derivs(1) = 0.;
            derivs(2) = 0.;
            derivs(3) = 0.;
            derivs(4) = 0.;
            derivs(5) = 0.;
            derivs(6) = 0.;
            derivs(7) = 0.;
            derivs(8) = 0.;
            derivs(9) = 0.;
            derivs(10) = 0.;
            derivs(11) = 0.;
         } else {

         }
         return derivs;
      }
   }

} // namespace flexiblesusy
