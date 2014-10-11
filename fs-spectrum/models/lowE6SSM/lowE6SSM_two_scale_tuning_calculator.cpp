#include "lowE6SSM_two_scale_tuning_calculator.hpp"
#include "numerics.hpp"

namespace flexiblesusy {

   lowE6SSM_tuning_calculator::lowE6SSM_tuning_calculator(const lowE6SSM<Two_scale>& m)
         : model(m)
         , input_scale(0.)
         , tuning_scale(0.)
         , tolerance(1.0e-4)
         , max_iterations(100)
         , tuning_ewsb_loop_order(1)
         , tuning_beta_loop_order(2)
   {
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

   bool lowE6SSM_tuning_calculator::calculate_fine_tunings_numerically()
   {
      return true;
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
      
      lowE6SSM_ew_derivs ew_derivs;
      Eigen::Matrix<double,num_ewsb_eqs,num_tuning_pars> beta_derivs;
      
      double scale = model.get_scale();

      const std::size_t num_pars = (tuning_ewsb_loop_order > 0 ? num_tree_level_ewsb_pars + num_one_loop_ewsb_pars
                                    : num_tree_level_ewsb_pars );

      if (is_equal_rel(scale, tuning_scale)) {
         ew_derivs.set_model(model);
         model.run_to(input_scale);
         beta_derivs = calculate_beta_derivs(num_pars);
         get_input_scale_pars();
      } else if (is_equal_rel(scale, input_scale)) {
         beta_derivs = calculate_beta_derivs(num_pars);
         get_input_scale_pars();
         model.run_to(tuning_scale);
         ew_derivs.set_model(model);
      } else {
         model.run_to(input_scale);
         beta_derivs = calculate_beta_derivs(num_pars);
         get_input_scale_pars();
         model.run_to(tuning_scale);
         ew_derivs.set_model(model);
      }

      ew_derivs.set_ewsb_loop_order(tuning_ewsb_loop_order);

      Eigen::Matrix<double,num_ewsb_eqs,num_ewsb_eqs> mass_matrix = ew_derivs.calculate_unrotated_mass_matrix_hh();

      // N.B. extra minus sign
      Eigen::Matrix<double, num_ewsb_eqs, Eigen::Dynamic> ewsb_derivs = -1.0 * ew_derivs.calculate_ewsb_parameter_derivs();

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

      double g1_at_tuning_scale = ew_derivs.get_g1();
      double g2_at_tuning_scale = ew_derivs.get_g2();
      double vd_at_tuning_scale = ew_derivs.get_vd();
      double vu_at_tuning_scale = ew_derivs.get_vu();
      const Eigen::Array<double,num_tuning_pars,1> g1_derivs = beta_derivs.row(get_g1_row());
      const Eigen::Array<double,num_tuning_pars,1> g2_derivs = beta_derivs.row(get_g2_row());
      const Eigen::Array<double,num_tuning_pars,1> vd_derivs = vev_derivs.row(0);
      const Eigen::Array<double,num_tuning_pars,1> vu_derivs = vev_derivs.row(1);

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
      const double underflow = 1.0e-100;

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

   double lowE6SSM_tuning_calculator::get_deriv(lowE6SSM_info::Parameters p, const Eigen::Array<double,num_tuning_pars,1>& derivs) const
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

   Eigen::Matrix<double,Eigen::Dynamic,lowE6SSM_tuning_calculator::num_tuning_pars> lowE6SSM_tuning_calculator::calculate_beta_derivs(std::size_t num_pars) const
   {
      Eigen::Matrix<double,Eigen::Dynamic,num_tuning_pars> derivs;

      derivs.col(0) = calculate_deriv_dlowscale_dLambdax(num_pars);
      derivs.col(1) = calculate_deriv_dlowscale_dTLambdax(num_pars);
      derivs.col(2) = calculate_deriv_dlowscale_dTYu22(num_pars);
      derivs.col(3) = calculate_deriv_dlowscale_dmq222(num_pars);
      derivs.col(4) = calculate_deriv_dlowscale_dmHd2(num_pars);
      derivs.col(5) = calculate_deriv_dlowscale_dmHu2(num_pars);
      derivs.col(6) = calculate_deriv_dlowscale_dmu222(num_pars);
      derivs.col(7) = calculate_deriv_dlowscale_dms2(num_pars);
      derivs.col(8) = calculate_deriv_dlowscale_dMassB(num_pars);
      derivs.col(9) = calculate_deriv_dlowscale_dMassWB(num_pars);
      derivs.col(10) = calculate_deriv_dlowscale_dMassG(num_pars);
      derivs.col(11) = calculate_deriv_dlowscale_dMassBp(num_pars);

      return derivs;
   }

   Eigen::VectorXd lowE6SSM_tuning_calculator::calculate_deriv_dlowscale_dLambdax(std::size_t num_pars) const
   {

   }

   Eigen::VectorXd lowE6SSM_tuning_calculator::calculate_deriv_dlowscale_dTLambdax(std::size_t num_pars) const
   {

   }

   Eigen::VectorXd lowE6SSM_tuning_calculator::calculate_deriv_dlowscale_dTYu22(std::size_t num_pars) const
   {

   }

   Eigen::VectorXd lowE6SSM_tuning_calculator::calculate_deriv_dlowscale_dmq222(std::size_t num_pars) const
   {

   }

   Eigen::VectorXd lowE6SSM_tuning_calculator::calculate_deriv_dlowscale_dmHd2(std::size_t num_pars) const
   {

   }

   Eigen::VectorXd lowE6SSM_tuning_calculator::calculate_deriv_dlowscale_dmHu2(std::size_t num_pars) const
   {

   }

   Eigen::VectorXd lowE6SSM_tuning_calculator::calculate_deriv_dlowscale_dmu222(std::size_t num_pars) const
   {

   }

   Eigen::VectorXd lowE6SSM_tuning_calculator::calculate_deriv_dlowscale_dms2(std::size_t num_pars) const
   {

   }

   Eigen::VectorXd lowE6SSM_tuning_calculator::calculate_deriv_dlowscale_dMassB(std::size_t num_pars) const
   {

   }

   Eigen::VectorXd lowE6SSM_tuning_calculator::calculate_deriv_dlowscale_dMassWB(std::size_t num_pars) const
   {

   }

   Eigen::VectorXd lowE6SSM_tuning_calculator::calculate_deriv_dlowscale_dMassG(std::size_t num_pars) const
   {

   }

   Eigen::VectorXd lowE6SSM_tuning_calculator::calculate_deriv_dlowscale_dMassBp(std::size_t num_pars) const
   {

   }
} // namespace flexiblesusy
