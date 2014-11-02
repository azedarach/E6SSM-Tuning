// ====================================================================
// Test suite for implementation of one-loop leading log
// coefficients
// ====================================================================

#include <limits>
#include <random>
#include <chrono>
#include <Eigen/Core>
#include <gsl/gsl_deriv.h>

#include "wrappers.hpp"
#include "lowE6SSM_two_scale_model.hpp"
#include "lowE6SSM_two_scale_ew_derivs.hpp"
#include "lowE6SSM_two_scale_tuning_calculator.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_lowE6SSM_leading_log_coeffs

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

using namespace std;
using namespace flexiblesusy;

lowE6SSM_input_parameters get_test_inputs(double theta)
{
   lowE6SSM_input_parameters input;

   input.TanBeta = 10.0;

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

   input.KappaInput(0,0) = 0.5;
   input.KappaInput(1,1) = 0.5;
   input.KappaInput(2,2) = 0.5;

   input.Lambda12Input(0,0) = 0.1;
   input.Lambda12Input(1,1) = 0.1;

   input.LambdaxInput = 0.5;

   input.MuPrInput = 10132.2;

   input.gNInput = 0.45;

   input.vsInput = 16234.2;

   input.TYdInput(0,0) = 132.4;
   input.TYdInput(1,1) = 123.5;
   input.TYdInput(2,2) = 567.2;

   input.TYeInput(0,0) = 99.2;
   input.TYeInput(1,1) = 140.5;
   input.TYeInput(2,2) = 342.2;

   input.TKappaInput(0,0) = 542.2;
   input.TKappaInput(1,1) = 673.3;
   input.TKappaInput(2,2) = 782.3;

   input.TLambda12Input(0,0) = 231.4;
   input.TLambda12Input(1,1) = 212.5;

   input.TLambdaxInput = 645.2;

   input.TYuInput(0,0) = 321.4;
   input.TYuInput(1,1) = 354.1;
   input.TYuInput(2,2) = 745.2;

   input.BMuPrInput = 11323.2;

   input.mq2Input(0,0) = 5000.0;
   input.mq2Input(1,1) = 5000.0;
   input.mq2Input(2,2) = 6231.7;

   input.ml2Input(0,0) = 3000.0;
   input.ml2Input(1,1) = 1000.0;
   input.ml2Input(2,2) = 3000.0;

   input.md2Input(0,0) = 1500.0;
   input.md2Input(1,1) = 2000.0;
   input.md2Input(2,2) = 2734.2;

   input.mu2Input(0,0) = 5500.0;
   input.mu2Input(1,1) = 5000.0;
   input.mu2Input(2,2) = 2341.5;

   input.me2Input(0,0) = 5000.0;
   input.me2Input(1,1) = 4300.0;
   input.me2Input(2,2) = 6000.0;

   input.mH1I2Input(0,0) = 3000.0;
   input.mH1I2Input(1,1) = 6000.0;

   input.mH2I2Input(0,0) = 4000.0;
   input.mH2I2Input(1,1) = 5000.0;

   input.msI2Input(0,0) = 4300.0;
   input.msI2Input(1,1) = 3210.0;

   input.mDx2Input(0,0) = 6000.0;
   input.mDx2Input(1,1) = 1523.1;
   input.mDx2Input(2,2) = 3000.0;

   input.mDxbar2Input(0,0) = 5000.0;
   input.mDxbar2Input(1,1) = 8921.2;
   input.mDxbar2Input(2,2) = 5000.0;

   input.mHp2Input = 7000.0;
   input.mHpbar2Input = 6000.0;
  
   input.MassBInput = 253.2;
   input.MassWBInput = 1324.2;
   input.MassGInput = 5329.3;
   input.MassBpInput = 300.0;

   return input;
}

void set_susy_parameters_from_input(const lowE6SSM_input_parameters& input, 
                                    lowE6SSM<Two_scale>& model)
{
   model.set_input_parameters(input);

   model.set_Yd(0, 0, 0.);
   model.set_Yd(1, 1, 0.);
   model.set_Yd(2, 2,0.1);

   model.set_Ye(0, 0, 0.);
   model.set_Ye(1, 1, 0.);
   model.set_Ye(2, 2, 0.05);

   model.set_Kappa(input.KappaInput);
   model.set_Lambda12(input.Lambda12Input);
   model.set_Lambdax(input.LambdaxInput);

   model.set_Yu(0, 0, 0.);
   model.set_Yu(1, 1, 0.);
   model.set_Yu(2, 2, 0.8);

   model.set_MuPr(input.MuPrInput);
   model.set_g1(0.46);
   model.set_g2(0.63);
   model.set_g3(1.02);
   model.set_gN(input.gNInput);
   model.set_vd(246.0 / Sqrt(1. + Sqr(input.TanBeta)));
   model.set_vu(246.0 * input.TanBeta / Sqrt(1. + Sqr(input.TanBeta)));
   model.set_vs(input.vsInput);

}

void set_soft_parameters_from_input(const lowE6SSM_input_parameters& input,
                                    lowE6SSM<Two_scale>& model)
{
   model.set_TYd(input.TYdInput);
   model.set_TYe(input.TYeInput);
   model.set_TKappa(input.TKappaInput);
   model.set_TLambda12(input.TLambda12Input);
   model.set_TLambdax(input.TLambdaxInput);
   model.set_TYu(input.TYuInput);
   model.set_BMuPr(input.BMuPrInput);
   model.set_mq2(input.mq2Input);
   model.set_ml2(input.ml2Input);
   model.set_mHd2(154253.2);
   model.set_mHu2(-23032.6);
   model.set_md2(input.md2Input);
   model.set_mu2(input.mu2Input);
   model.set_me2(input.me2Input);
   model.set_ms2(-123534.);
   model.set_mH1I2(input.mH1I2Input);
   model.set_mH2I2(input.mH2I2Input);
   model.set_msI2(input.msI2Input);
   model.set_mDx2(input.mDx2Input);
   model.set_mDxbar2(input.mDxbar2Input);
   model.set_mHp2(input.mHp2Input);
   model.set_mHpbar2(input.mHpbar2Input);
   model.set_MassB(input.MassBInput);
   model.set_MassWB(input.MassWBInput);
   model.set_MassG(input.MassGInput);
   model.set_MassBp(input.MassBpInput);

}

void initialize_model(const lowE6SSM_input_parameters& input, 
                      lowE6SSM<Two_scale>& model)
{
   set_susy_parameters_from_input(input, model);
   set_soft_parameters_from_input(input, model);
   model.set_scale(20000.0);
}

struct beta_params {
   lowE6SSM<Two_scale>* model;
};

double beta_Lambdax(double t, void * params)
{
   const beta_params* pars
      = static_cast<beta_params*>(params);

   lowE6SSM<Two_scale>* model = pars->model;

   model->run_to(std::exp(t));

   lowE6SSM_soft_parameters beta_fns = model->calc_beta();

   double beta_fn = beta_fns.get_Lambdax();

   return beta_fn;
}

double beta_Yu22(double t, void * params)
{
   const beta_params* pars
      = static_cast<beta_params*>(params);

   lowE6SSM<Two_scale>* model = pars->model;

   model->run_to(std::exp(t));

   lowE6SSM_soft_parameters beta_fns = model->calc_beta();

   double beta_fn = beta_fns.get_Yu(2,2);

   return beta_fn;
}

double beta_g1(double t, void * params)
{
   const beta_params* pars
      = static_cast<beta_params*>(params);

   lowE6SSM<Two_scale>* model = pars->model;

   model->run_to(std::exp(t));

   lowE6SSM_soft_parameters beta_fns = model->calc_beta();

   double beta_fn = beta_fns.get_g1();

   return beta_fn;
}

double beta_g2(double t, void * params)
{
   const beta_params* pars
      = static_cast<beta_params*>(params);

   lowE6SSM<Two_scale>* model = pars->model;

   model->run_to(std::exp(t));

   lowE6SSM_soft_parameters beta_fns = model->calc_beta();

   double beta_fn = beta_fns.get_g2();

   return beta_fn;
}

double beta_gN(double t, void * params)
{
   const beta_params* pars
      = static_cast<beta_params*>(params);

   lowE6SSM<Two_scale>* model = pars->model;

   model->run_to(std::exp(t));

   lowE6SSM_soft_parameters beta_fns = model->calc_beta();

   double beta_fn = beta_fns.get_gN();

   return beta_fn;
}

double beta_TLambdax(double t, void * params)
{
   const beta_params* pars
      = static_cast<beta_params*>(params);

   lowE6SSM<Two_scale>* model = pars->model;

   model->run_to(std::exp(t));

   lowE6SSM_soft_parameters beta_fns = model->calc_beta();

   double beta_fn = beta_fns.get_TLambdax();

   return beta_fn;
}

double beta_TYu22(double t, void * params)
{
   const beta_params* pars
      = static_cast<beta_params*>(params);

   lowE6SSM<Two_scale>* model = pars->model;

   model->run_to(std::exp(t));

   lowE6SSM_soft_parameters beta_fns = model->calc_beta();

   double beta_fn = beta_fns.get_TYu(2,2);

   return beta_fn;
}

double beta_mq222(double t, void * params)
{
   const beta_params* pars
      = static_cast<beta_params*>(params);

   lowE6SSM<Two_scale>* model = pars->model;

   model->run_to(std::exp(t));

   lowE6SSM_soft_parameters beta_fns = model->calc_beta();

   double beta_fn = beta_fns.get_mq2(2,2);

   return beta_fn;
}

double beta_mHd2(double t, void * params)
{
   const beta_params* pars
      = static_cast<beta_params*>(params);

   lowE6SSM<Two_scale>* model = pars->model;

   model->run_to(std::exp(t));

   lowE6SSM_soft_parameters beta_fns = model->calc_beta();

   double beta_fn = beta_fns.get_mHd2();

   return beta_fn;
}

double beta_mHu2(double t, void * params)
{
   const beta_params* pars
      = static_cast<beta_params*>(params);

   lowE6SSM<Two_scale>* model = pars->model;

   model->run_to(std::exp(t));

   lowE6SSM_soft_parameters beta_fns = model->calc_beta();

   double beta_fn = beta_fns.get_mHu2();

   return beta_fn;
}

double beta_mu222(double t, void * params)
{
   const beta_params* pars
      = static_cast<beta_params*>(params);

   lowE6SSM<Two_scale>* model = pars->model;

   model->run_to(std::exp(t));

   lowE6SSM_soft_parameters beta_fns = model->calc_beta();

   double beta_fn = beta_fns.get_mu2(2,2);

   return beta_fn;
}

double beta_ms2(double t, void * params)
{
   const beta_params* pars
      = static_cast<beta_params*>(params);

   lowE6SSM<Two_scale>* model = pars->model;

   model->run_to(std::exp(t));

   lowE6SSM_soft_parameters beta_fns = model->calc_beta();

   double beta_fn = beta_fns.get_ms2();

   return beta_fn;
}

BOOST_AUTO_TEST_CASE( test_leading_log_coeff_Lambdax )
{
   const double max_err = 1.0e-4;
   const double theta = ArcTan(Sqrt(15.));

   lowE6SSM_input_parameters inputs = get_test_inputs(theta);

   lowE6SSM<Two_scale> model;

   initialize_model(inputs, model);

   model.set_loops(1);

   lowE6SSM_tuning_calculator tuning_calc(model);

   double analytic_deriv = tuning_calc.calculate_leading_log_coefficient(lowE6SSM_info::Lambdax);

   beta_params pars = {&model};

   gsl_function func = {beta_Lambdax, &pars};

   double numeric_deriv;
   double abs_err;

   gsl_deriv_central(&func, Log(model.get_scale()), 1.0e-3, &numeric_deriv, &abs_err);

   BOOST_REQUIRE(Abs(abs_err) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv - numeric_deriv), abs_err);

}

BOOST_AUTO_TEST_CASE( test_leading_log_coeff_Yu22 )
{
   const double max_err = 1.0e-4;
   const double theta = ArcTan(Sqrt(15.));

   lowE6SSM_input_parameters inputs = get_test_inputs(theta);

   lowE6SSM<Two_scale> model;

   initialize_model(inputs, model);

   model.set_loops(1);

   lowE6SSM_tuning_calculator tuning_calc(model);

   double analytic_deriv = tuning_calc.calculate_leading_log_coefficient(lowE6SSM_info::Yu22);

   beta_params pars = {&model};

   gsl_function func = {beta_Yu22, &pars};

   double numeric_deriv;
   double abs_err;

   gsl_deriv_central(&func, Log(model.get_scale()), 1.0e-3, &numeric_deriv, &abs_err);

   BOOST_REQUIRE(Abs(abs_err) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv - numeric_deriv), abs_err);

}

BOOST_AUTO_TEST_CASE( test_leading_log_coeff_g1 )
{
   const double max_err = 1.0e-4;
   const double theta = ArcTan(Sqrt(15.));

   lowE6SSM_input_parameters inputs = get_test_inputs(theta);

   lowE6SSM<Two_scale> model;

   initialize_model(inputs, model);

   model.set_loops(1);

   lowE6SSM_tuning_calculator tuning_calc(model);

   double analytic_deriv = tuning_calc.calculate_leading_log_coefficient(lowE6SSM_info::g1);

   beta_params pars = {&model};

   gsl_function func = {beta_g1, &pars};

   double numeric_deriv;
   double abs_err;

   gsl_deriv_central(&func, Log(model.get_scale()), 1.0e-3, &numeric_deriv, &abs_err);

   BOOST_REQUIRE(Abs(abs_err) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv - numeric_deriv), abs_err);

}

BOOST_AUTO_TEST_CASE( test_leading_log_coeff_g2 )
{
   const double max_err = 1.0e-4;
   const double theta = ArcTan(Sqrt(15.));

   lowE6SSM_input_parameters inputs = get_test_inputs(theta);

   lowE6SSM<Two_scale> model;

   initialize_model(inputs, model);

   model.set_loops(1);

   lowE6SSM_tuning_calculator tuning_calc(model);

   double analytic_deriv = tuning_calc.calculate_leading_log_coefficient(lowE6SSM_info::g2);

   beta_params pars = {&model};

   gsl_function func = {beta_g2, &pars};

   double numeric_deriv;
   double abs_err;

   gsl_deriv_central(&func, Log(model.get_scale()), 1.0e-3, &numeric_deriv, &abs_err);

   BOOST_REQUIRE(Abs(abs_err) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv - numeric_deriv), abs_err);

}

BOOST_AUTO_TEST_CASE( test_leading_log_coeff_gN )
{
   const double max_err = 1.0e-4;
   const double theta = ArcTan(Sqrt(15.));

   lowE6SSM_input_parameters inputs = get_test_inputs(theta);

   lowE6SSM<Two_scale> model;

   initialize_model(inputs, model);

   model.set_loops(1);

   lowE6SSM_tuning_calculator tuning_calc(model);

   double analytic_deriv = tuning_calc.calculate_leading_log_coefficient(lowE6SSM_info::gN);

   beta_params pars = {&model};

   gsl_function func = {beta_gN, &pars};

   double numeric_deriv;
   double abs_err;

   gsl_deriv_central(&func, Log(model.get_scale()), 1.0e-3, &numeric_deriv, &abs_err);

   BOOST_REQUIRE(Abs(abs_err) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv - numeric_deriv), abs_err);

}

BOOST_AUTO_TEST_CASE( test_leading_log_coeff_TLambdax )
{
   const double max_err = 1.0e-4;
   const double theta = ArcTan(Sqrt(15.));

   lowE6SSM_input_parameters inputs = get_test_inputs(theta);

   lowE6SSM<Two_scale> model;

   initialize_model(inputs, model);

   model.set_loops(1);

   lowE6SSM_tuning_calculator tuning_calc(model);

   double analytic_deriv = tuning_calc.calculate_leading_log_coefficient(lowE6SSM_info::TLambdax);

   beta_params pars = {&model};

   gsl_function func = {beta_TLambdax, &pars};

   double numeric_deriv;
   double abs_err;

   gsl_deriv_central(&func, Log(model.get_scale()), 1.0e-3, &numeric_deriv, &abs_err);

   BOOST_REQUIRE(Abs(abs_err) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv - numeric_deriv), abs_err);

}

BOOST_AUTO_TEST_CASE( test_leading_log_coeff_TYu22 )
{
   const double max_err = 1.0e-4;
   const double theta = ArcTan(Sqrt(15.));

   lowE6SSM_input_parameters inputs = get_test_inputs(theta);

   lowE6SSM<Two_scale> model;

   initialize_model(inputs, model);

   model.set_loops(1);

   lowE6SSM_tuning_calculator tuning_calc(model);

   double analytic_deriv = tuning_calc.calculate_leading_log_coefficient(lowE6SSM_info::TYu22);

   beta_params pars = {&model};

   gsl_function func = {beta_TYu22, &pars};

   double numeric_deriv;
   double abs_err;

   gsl_deriv_central(&func, Log(model.get_scale()), 1.0e-3, &numeric_deriv, &abs_err);

   BOOST_REQUIRE(Abs(abs_err) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv - numeric_deriv), abs_err);

}

BOOST_AUTO_TEST_CASE( test_leading_log_coeff_mq222 )
{
   const double max_err = 1.0e-4;
   const double theta = ArcTan(Sqrt(15.));

   lowE6SSM_input_parameters inputs = get_test_inputs(theta);

   lowE6SSM<Two_scale> model;

   initialize_model(inputs, model);

   model.set_loops(1);

   lowE6SSM_tuning_calculator tuning_calc(model);

   double analytic_deriv = tuning_calc.calculate_leading_log_coefficient(lowE6SSM_info::mq222);

   beta_params pars = {&model};

   gsl_function func = {beta_mq222, &pars};

   double numeric_deriv;
   double abs_err;

   gsl_deriv_central(&func, Log(model.get_scale()), 1.0e-3, &numeric_deriv, &abs_err);

   BOOST_REQUIRE(Abs(abs_err) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv - numeric_deriv), abs_err);

}

BOOST_AUTO_TEST_CASE( test_leading_log_coeff_mHd2 )
{
   const double max_err = 1.0e-4;
   const double theta = ArcTan(Sqrt(15.));

   lowE6SSM_input_parameters inputs = get_test_inputs(theta);

   lowE6SSM<Two_scale> model;

   initialize_model(inputs, model);

   model.set_loops(1);

   lowE6SSM_tuning_calculator tuning_calc(model);

   double analytic_deriv = tuning_calc.calculate_leading_log_coefficient(lowE6SSM_info::mHd2);

   beta_params pars = {&model};

   gsl_function func = {beta_mHd2, &pars};

   double numeric_deriv;
   double abs_err;

   gsl_deriv_central(&func, Log(model.get_scale()), 1.0e-3, &numeric_deriv, &abs_err);

   BOOST_REQUIRE(Abs(abs_err) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv - numeric_deriv), abs_err);

}

BOOST_AUTO_TEST_CASE( test_leading_log_coeff_mHu2 )
{
   const double max_err = 1.0e-4;
   const double theta = ArcTan(Sqrt(15.));

   lowE6SSM_input_parameters inputs = get_test_inputs(theta);

   lowE6SSM<Two_scale> model;

   initialize_model(inputs, model);

   model.set_loops(1);

   lowE6SSM_tuning_calculator tuning_calc(model);

   double analytic_deriv = tuning_calc.calculate_leading_log_coefficient(lowE6SSM_info::mHu2);

   beta_params pars = {&model};

   gsl_function func = {beta_mHu2, &pars};

   double numeric_deriv;
   double abs_err;

   gsl_deriv_central(&func, Log(model.get_scale()), 1.0e-3, &numeric_deriv, &abs_err);

   BOOST_REQUIRE(Abs(abs_err) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv - numeric_deriv), abs_err);

}

BOOST_AUTO_TEST_CASE( test_leading_log_coeff_mu222 )
{
   const double max_err = 1.0e-4;
   const double theta = ArcTan(Sqrt(15.));

   lowE6SSM_input_parameters inputs = get_test_inputs(theta);

   lowE6SSM<Two_scale> model;

   initialize_model(inputs, model);

   model.set_loops(1);

   lowE6SSM_tuning_calculator tuning_calc(model);

   double analytic_deriv = tuning_calc.calculate_leading_log_coefficient(lowE6SSM_info::mu222);

   beta_params pars = {&model};

   gsl_function func = {beta_mu222, &pars};

   double numeric_deriv;
   double abs_err;

   gsl_deriv_central(&func, Log(model.get_scale()), 1.0e-3, &numeric_deriv, &abs_err);

   BOOST_REQUIRE(Abs(abs_err) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv - numeric_deriv), abs_err);

}

BOOST_AUTO_TEST_CASE( test_leading_log_coeff_ms2 )
{
   const double max_err = 1.0e-4;
   const double theta = ArcTan(Sqrt(15.));

   lowE6SSM_input_parameters inputs = get_test_inputs(theta);

   lowE6SSM<Two_scale> model;

   initialize_model(inputs, model);

   model.set_loops(1);

   lowE6SSM_tuning_calculator tuning_calc(model);

   double analytic_deriv = tuning_calc.calculate_leading_log_coefficient(lowE6SSM_info::ms2);

   beta_params pars = {&model};

   gsl_function func = {beta_ms2, &pars};

   double numeric_deriv;
   double abs_err;

   gsl_deriv_central(&func, Log(model.get_scale()), 1.0e-3, &numeric_deriv, &abs_err);

   BOOST_REQUIRE(Abs(abs_err) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv - numeric_deriv), abs_err);

}
