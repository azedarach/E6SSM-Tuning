// ====================================================================
// Test suite for methods used in computing derivatives of the 
// EWSB conditions
// ====================================================================

#include <limits>
#include <random>
#include <chrono>
#include <Eigen/Core>
#include <gsl/gsl_deriv.h>

#include "ew_input.hpp"
#include "wrappers.hpp"
#include "lowE6SSM_two_scale_model.hpp"
#include "lowE6SSM_two_scale_ew_derivs.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_lowE6SSM_ew_derivs

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
   model.set_scale(Electroweak_constants::MZ);
   set_susy_parameters_from_input(input, model);
   set_soft_parameters_from_input(input, model);
}

// Test solution of EWSB via soft masses at tree level
BOOST_AUTO_TEST_CASE( test_ewsb_tree_level_soln_soft_masses )
{
   const double tol = 1.0e-4;
   const double theta = ArcTan(Sqrt(15.));

   lowE6SSM_input_parameters inputs = get_test_inputs(theta);

   lowE6SSM<Two_scale> model;

   initialize_model(inputs, model);

   lowE6SSM_ew_derivs ew_derivs(model);

   ew_derivs.set_ewsb_loop_order(0);

   BOOST_CHECK_LE(Abs(ew_derivs.get_ewsb_condition_1()), tol);
   BOOST_CHECK_LE(Abs(ew_derivs.get_ewsb_condition_2()), tol);
   BOOST_CHECK_LE(Abs(ew_derivs.get_ewsb_condition_3()), tol);
}

// Test solution of EWSB via soft masses at 1-loop
BOOST_AUTO_TEST_CASE( test_ewsb_one_loop_soln_soft_masses )
{
   const double tol = 1.0e-4;
   const double theta = ArcTan(Sqrt(15.));

   lowE6SSM_input_parameters inputs = get_test_inputs(theta);

   lowE6SSM<Two_scale> model;

   initialize_model(inputs, model);

   lowE6SSM_ew_derivs ew_derivs(model);

   ew_derivs.set_ewsb_loop_order(1);

   BOOST_CHECK_LE(Abs(ew_derivs.get_ewsb_condition_1()), tol);
   BOOST_CHECK_LE(Abs(ew_derivs.get_ewsb_condition_2()), tol);
   BOOST_CHECK_LE(Abs(ew_derivs.get_ewsb_condition_3()), tol);
}

struct deltaV_params {
   lowE6SSM_ew_derivs* derivs;
   lowE6SSM_info::Parameters p; 
};

double deltaV(double x, void * params)
{
   const double oneOverSqrt2 = 1.0 / Sqrt(2.);

   const deltaV_params* pars
      = static_cast<deltaV_params*>(params);

   lowE6SSM_ew_derivs* derivs = pars->derivs;
   lowE6SSM_info::Parameters p = pars->p;

   derivs->get_model().set_parameter(p, x);

   double scale = derivs->get_model().get_scale();

   double mt = oneOverSqrt2 * derivs->get_model().get_Yu(2,2)
      * derivs->get_model().get_vu();

   double msf1 = 0.;
   double msf2 = 0.;
   double theta = 0.;

   derivs->get_model().calculate_MSu_3rd_generation(msf1, msf2, theta);

   double tmp_1 = Power(msf1, 4) * (Log(Sqr(msf1) / Sqr(scale)) - 1.5);
   double tmp_2 = Power(msf2, 4) * (Log(Sqr(msf2) / Sqr(scale)) - 1.5);
   double tmp_3 = -2.0 * Power(mt, 4) * (Log(Sqr(mt) / Sqr(scale)) - 1.5);

   double result = oneOver16PiSqr * 1.5 * (tmp_1 + tmp_2 + tmp_3);

   return result;
}

BOOST_AUTO_TEST_CASE( test_dDeltaV_dvd )
{
   const double max_err = 0.1;
   const double theta = ArcTan(Sqrt(15.));
   const lowE6SSM_info::Parameters p = lowE6SSM_info::vd;

   lowE6SSM_input_parameters inputs = get_test_inputs(theta);

   lowE6SSM<Two_scale> model;

   initialize_model(inputs, model);

   lowE6SSM_ew_derivs ew_derivs(model);

   double analytic_deriv = ew_derivs.deriv_dDeltaV_dparam(p);

   deltaV_params pars = {&ew_derivs, p };

   gsl_function func = { deltaV, &pars };

   double numeric_deriv;
   double abs_err;

   gsl_deriv_central(&func, model.get_parameter(p), 1.0, &numeric_deriv, &abs_err);

   BOOST_REQUIRE(Abs(abs_err) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv - numeric_deriv), abs_err);

}

BOOST_AUTO_TEST_CASE( test_dDeltaV_dvu )
{
   const double max_err = 0.1;
   const double theta = ArcTan(Sqrt(15.));
   const lowE6SSM_info::Parameters p = lowE6SSM_info::vu;

   lowE6SSM_input_parameters inputs = get_test_inputs(theta);

   lowE6SSM<Two_scale> model;

   initialize_model(inputs, model);

   lowE6SSM_ew_derivs ew_derivs(model);

   double analytic_deriv = ew_derivs.deriv_dDeltaV_dparam(p);

   deltaV_params pars = {&ew_derivs, p };

   gsl_function func = { deltaV, &pars };

   double numeric_deriv;
   double abs_err;

   gsl_deriv_central(&func, model.get_parameter(p), 1.0, &numeric_deriv, &abs_err);

   BOOST_REQUIRE(Abs(abs_err) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv - numeric_deriv), abs_err);

}

BOOST_AUTO_TEST_CASE( test_dDeltaV_dvs )
{
   const double max_err = 0.1;
   const double theta = ArcTan(Sqrt(15.));
   const lowE6SSM_info::Parameters p = lowE6SSM_info::vs;

   lowE6SSM_input_parameters inputs = get_test_inputs(theta);

   lowE6SSM<Two_scale> model;

   initialize_model(inputs, model);

   lowE6SSM_ew_derivs ew_derivs(model);

   double analytic_deriv = ew_derivs.deriv_dDeltaV_dparam(p);

   deltaV_params pars = {&ew_derivs, p };

   gsl_function func = { deltaV, &pars };

   double numeric_deriv;
   double abs_err;

   gsl_deriv_central(&func, model.get_parameter(p), 1.0, &numeric_deriv, &abs_err);

   BOOST_REQUIRE(Abs(abs_err) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv - numeric_deriv), abs_err);

}

struct ewsb_condition_params {
   lowE6SSM_ew_derivs* derivs;
   lowE6SSM_info::Parameters p;
};

double ewsb_condition_1(double x, void * params)
{
   const ewsb_condition_params* pars 
      = static_cast<ewsb_condition_params*>(params);

   lowE6SSM_ew_derivs* derivs = pars->derivs;
   lowE6SSM_info::Parameters p = pars->p;

   derivs->get_model().set_parameter(p, x);

   return derivs->get_ewsb_condition_1();
}

double ewsb_condition_2(double x, void * params)
{
   const ewsb_condition_params* pars 
      = static_cast<ewsb_condition_params*>(params);

   lowE6SSM_ew_derivs* derivs = pars->derivs;
   lowE6SSM_info::Parameters p = pars->p;

   derivs->get_model().set_parameter(p, x);

   return derivs->get_ewsb_condition_2();
}

double ewsb_condition_3(double x, void * params)
{
   const ewsb_condition_params* pars 
      = static_cast<ewsb_condition_params*>(params);

   lowE6SSM_ew_derivs* derivs = pars->derivs;
   lowE6SSM_info::Parameters p = pars->p;

   derivs->get_model().set_parameter(p, x);

   return derivs->get_ewsb_condition_3();
}

BOOST_AUTO_TEST_CASE( test_tree_level_dewsb_conditions_dLambdax )
{
   const double max_err = 0.01;
   const double theta = ArcTan(Sqrt(15.));
   const lowE6SSM_info::Parameters p = lowE6SSM_info::Lambdax;

   lowE6SSM_input_parameters inputs = get_test_inputs(theta);

   lowE6SSM<Two_scale> model;

   initialize_model(inputs, model);

   model.set_ewsb_loop_order(0);

   lowE6SSM_ew_derivs ew_derivs(model);

   Eigen::Matrix<double,3,1> analytic_deriv 
      = ew_derivs.calculate_ewsb_parameter_derivs(p);

   ewsb_condition_params pars = {&ew_derivs, p };

   gsl_function func_1 = { ewsb_condition_1, &pars };
   gsl_function func_2 = { ewsb_condition_2, &pars };
   gsl_function func_3 = { ewsb_condition_3, &pars };

   double numeric_deriv_1;
   double abs_err_1;

   gsl_deriv_central(&func_1, model.get_parameter(p), 1.0e-1, &numeric_deriv_1, &abs_err_1);

   BOOST_REQUIRE(Abs(abs_err_1) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(0) - numeric_deriv_1), abs_err_1);

   double numeric_deriv_2;
   double abs_err_2;

   gsl_deriv_central(&func_2, model.get_parameter(p), 1.0e-1, &numeric_deriv_2, &abs_err_2);

   BOOST_REQUIRE(Abs(abs_err_2) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(1) - numeric_deriv_2), abs_err_2);

   double numeric_deriv_3;
   double abs_err_3;

   gsl_deriv_central(&func_3, model.get_parameter(p), 1.0e-1, &numeric_deriv_3, &abs_err_3);

   BOOST_REQUIRE(Abs(abs_err_3) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(2) - numeric_deriv_3), abs_err_3);

}

BOOST_AUTO_TEST_CASE( test_tree_level_dewsb_conditions_dg1 )
{
   const double max_err = 0.01;
   const double theta = ArcTan(Sqrt(15.));
   const lowE6SSM_info::Parameters p = lowE6SSM_info::g1;

   lowE6SSM_input_parameters inputs = get_test_inputs(theta);

   lowE6SSM<Two_scale> model;

   initialize_model(inputs, model);

   model.set_ewsb_loop_order(0);

   lowE6SSM_ew_derivs ew_derivs(model);

   Eigen::Matrix<double,3,1> analytic_deriv 
      = ew_derivs.calculate_ewsb_parameter_derivs(p);

   ewsb_condition_params pars = {&ew_derivs, p };

   gsl_function func_1 = { ewsb_condition_1, &pars };
   gsl_function func_2 = { ewsb_condition_2, &pars };
   gsl_function func_3 = { ewsb_condition_3, &pars };

   double numeric_deriv_1;
   double abs_err_1;

   gsl_deriv_central(&func_1, model.get_parameter(p), 1.0e-2, &numeric_deriv_1, &abs_err_1);

   BOOST_REQUIRE(Abs(abs_err_1) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(0) - numeric_deriv_1), abs_err_1);

   double numeric_deriv_2;
   double abs_err_2;

   gsl_deriv_central(&func_2, model.get_parameter(p), 1.0e-3, &numeric_deriv_2, &abs_err_2);

   BOOST_REQUIRE(Abs(abs_err_2) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(1) - numeric_deriv_2), abs_err_2);

   double numeric_deriv_3;
   double abs_err_3;

   gsl_deriv_central(&func_3, model.get_parameter(p), 1.0e-1, &numeric_deriv_3, &abs_err_3);

   BOOST_REQUIRE(Abs(abs_err_3) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(2) - numeric_deriv_3), abs_err_3);

}

BOOST_AUTO_TEST_CASE( test_tree_level_dewsb_conditions_dg2 )
{
   const double max_err = 0.01;
   const double theta = ArcTan(Sqrt(15.));
   const lowE6SSM_info::Parameters p = lowE6SSM_info::g2;

   lowE6SSM_input_parameters inputs = get_test_inputs(theta);

   lowE6SSM<Two_scale> model;

   initialize_model(inputs, model);

   model.set_ewsb_loop_order(0);

   lowE6SSM_ew_derivs ew_derivs(model);

   Eigen::Matrix<double,3,1> analytic_deriv 
      = ew_derivs.calculate_ewsb_parameter_derivs(p);

   ewsb_condition_params pars = {&ew_derivs, p };

   gsl_function func_1 = { ewsb_condition_1, &pars };
   gsl_function func_2 = { ewsb_condition_2, &pars };
   gsl_function func_3 = { ewsb_condition_3, &pars };

   double numeric_deriv_1;
   double abs_err_1;

   gsl_deriv_central(&func_1, model.get_parameter(p), 1.0e-3, &numeric_deriv_1, &abs_err_1);

   BOOST_REQUIRE(Abs(abs_err_1) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(0) - numeric_deriv_1), abs_err_1);

   double numeric_deriv_2;
   double abs_err_2;

   gsl_deriv_central(&func_2, model.get_parameter(p), 1.0e-4, &numeric_deriv_2, &abs_err_2);

   BOOST_REQUIRE(Abs(abs_err_2) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(1) - numeric_deriv_2), abs_err_2);

   double numeric_deriv_3;
   double abs_err_3;

   gsl_deriv_central(&func_3, model.get_parameter(p), 1.0e-1, &numeric_deriv_3, &abs_err_3);

   BOOST_REQUIRE(Abs(abs_err_3) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(2) - numeric_deriv_3), abs_err_3);

}

BOOST_AUTO_TEST_CASE( test_tree_level_dewsb_conditions_dgN )
{
   const double max_err = 0.01;
   const double theta = ArcTan(Sqrt(15.));
   const lowE6SSM_info::Parameters p = lowE6SSM_info::gN;

   lowE6SSM_input_parameters inputs = get_test_inputs(theta);

   lowE6SSM<Two_scale> model;

   initialize_model(inputs, model);

   model.set_ewsb_loop_order(0);

   lowE6SSM_ew_derivs ew_derivs(model);

   Eigen::Matrix<double,3,1> analytic_deriv 
      = ew_derivs.calculate_ewsb_parameter_derivs(p);

   ewsb_condition_params pars = {&ew_derivs, p };

   gsl_function func_1 = { ewsb_condition_1, &pars };
   gsl_function func_2 = { ewsb_condition_2, &pars };
   gsl_function func_3 = { ewsb_condition_3, &pars };

   double numeric_deriv_1;
   double abs_err_1;

   gsl_deriv_central(&func_1, model.get_parameter(p), 1.0e-1, &numeric_deriv_1, &abs_err_1);

   BOOST_REQUIRE(Abs(abs_err_1) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(0) - numeric_deriv_1), abs_err_1);

   double numeric_deriv_2;
   double abs_err_2;

   gsl_deriv_central(&func_2, model.get_parameter(p), 1.0e-1, &numeric_deriv_2, &abs_err_2);

   BOOST_REQUIRE(Abs(abs_err_2) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(1) - numeric_deriv_2), abs_err_2);

   double numeric_deriv_3;
   double abs_err_3;

   gsl_deriv_central(&func_3, model.get_parameter(p), 1.0e-1, &numeric_deriv_3, &abs_err_3);

   BOOST_REQUIRE(Abs(abs_err_3) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(2) - numeric_deriv_3), abs_err_3);

}

BOOST_AUTO_TEST_CASE( test_tree_level_dewsb_conditions_dvd )
{
   const double max_err = 0.01;
   const double theta = ArcTan(Sqrt(15.));
   const lowE6SSM_info::Parameters p = lowE6SSM_info::vd;

   lowE6SSM_input_parameters inputs = get_test_inputs(theta);

   lowE6SSM<Two_scale> model;

   initialize_model(inputs, model);

   model.set_ewsb_loop_order(0);

   lowE6SSM_ew_derivs ew_derivs(model);

   Eigen::Matrix<double,3,1> analytic_deriv 
      = ew_derivs.calculate_ewsb_parameter_derivs(p);

   ewsb_condition_params pars = {&ew_derivs, p };

   gsl_function func_1 = { ewsb_condition_1, &pars };
   gsl_function func_2 = { ewsb_condition_2, &pars };
   gsl_function func_3 = { ewsb_condition_3, &pars };

   double numeric_deriv_1;
   double abs_err_1;

   gsl_deriv_central(&func_1, model.get_parameter(p), 1.0e-2, &numeric_deriv_1, &abs_err_1);

   BOOST_REQUIRE(Abs(abs_err_1) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(0) - numeric_deriv_1), abs_err_1);

   double numeric_deriv_2;
   double abs_err_2;

   gsl_deriv_central(&func_2, model.get_parameter(p), 1.0e-2, &numeric_deriv_2, &abs_err_2);

   BOOST_REQUIRE(Abs(abs_err_2) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(1) - numeric_deriv_2), abs_err_2);

   double numeric_deriv_3;
   double abs_err_3;

   gsl_deriv_central(&func_3, model.get_parameter(p), 1.0e-2, &numeric_deriv_3, &abs_err_3);

   BOOST_REQUIRE(Abs(abs_err_3) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(2) - numeric_deriv_3), abs_err_3);

}

BOOST_AUTO_TEST_CASE( test_tree_level_dewsb_conditions_dvu )
{
   const double max_err = 0.01;
   const double theta = ArcTan(Sqrt(15.));
   const lowE6SSM_info::Parameters p = lowE6SSM_info::vu;

   lowE6SSM_input_parameters inputs = get_test_inputs(theta);

   lowE6SSM<Two_scale> model;

   initialize_model(inputs, model);

   model.set_ewsb_loop_order(0);

   lowE6SSM_ew_derivs ew_derivs(model);

   Eigen::Matrix<double,3,1> analytic_deriv 
      = ew_derivs.calculate_ewsb_parameter_derivs(p);

   ewsb_condition_params pars = {&ew_derivs, p };

   gsl_function func_1 = { ewsb_condition_1, &pars };
   gsl_function func_2 = { ewsb_condition_2, &pars };
   gsl_function func_3 = { ewsb_condition_3, &pars };

   double numeric_deriv_1;
   double abs_err_1;

   gsl_deriv_central(&func_1, model.get_parameter(p), 1.0e-3, &numeric_deriv_1, &abs_err_1);

   BOOST_REQUIRE(Abs(abs_err_1) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(0) - numeric_deriv_1), abs_err_1);

   double numeric_deriv_2;
   double abs_err_2;

   gsl_deriv_central(&func_2, model.get_parameter(p), 1.0e-3, &numeric_deriv_2, &abs_err_2);

   BOOST_REQUIRE(Abs(abs_err_2) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(1) - numeric_deriv_2), abs_err_2);

   double numeric_deriv_3;
   double abs_err_3;

   gsl_deriv_central(&func_3, model.get_parameter(p), 10.0, &numeric_deriv_3, &abs_err_3);

   BOOST_REQUIRE(Abs(abs_err_3) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(2) - numeric_deriv_3), abs_err_3);

}

BOOST_AUTO_TEST_CASE( test_tree_level_dewsb_conditions_dvs )
{
   const double max_err = 0.01;
   const double theta = ArcTan(Sqrt(15.));
   const lowE6SSM_info::Parameters p = lowE6SSM_info::vs;

   lowE6SSM_input_parameters inputs = get_test_inputs(theta);

   lowE6SSM<Two_scale> model;

   initialize_model(inputs, model);

   model.set_ewsb_loop_order(0);

   lowE6SSM_ew_derivs ew_derivs(model);

   Eigen::Matrix<double,3,1> analytic_deriv 
      = ew_derivs.calculate_ewsb_parameter_derivs(p);

   ewsb_condition_params pars = {&ew_derivs, p };

   gsl_function func_1 = { ewsb_condition_1, &pars };
   gsl_function func_2 = { ewsb_condition_2, &pars };
   gsl_function func_3 = { ewsb_condition_3, &pars };

   double numeric_deriv_1;
   double abs_err_1;

   gsl_deriv_central(&func_1, model.get_parameter(p), 1.0, &numeric_deriv_1, &abs_err_1);

   BOOST_REQUIRE(Abs(abs_err_1) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(0) - numeric_deriv_1), abs_err_1);

   double numeric_deriv_2;
   double abs_err_2;

   gsl_deriv_central(&func_2, model.get_parameter(p), 1.0e-1, &numeric_deriv_2, &abs_err_2);

   BOOST_REQUIRE(Abs(abs_err_2) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(1) - numeric_deriv_2), abs_err_2);

   double numeric_deriv_3;
   double abs_err_3;

   gsl_deriv_central(&func_3, model.get_parameter(p), 1.0e-1, &numeric_deriv_3, &abs_err_3);

   BOOST_REQUIRE(Abs(abs_err_3) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(2) - numeric_deriv_3), abs_err_3);

}

BOOST_AUTO_TEST_CASE( test_tree_level_dewsb_conditions_dTLambdax )
{
   const double max_err = 0.01;
   const double theta = ArcTan(Sqrt(15.));
   const lowE6SSM_info::Parameters p = lowE6SSM_info::TLambdax;

   lowE6SSM_input_parameters inputs = get_test_inputs(theta);

   lowE6SSM<Two_scale> model;

   initialize_model(inputs, model);

   model.set_ewsb_loop_order(0);

   lowE6SSM_ew_derivs ew_derivs(model);

   Eigen::Matrix<double,3,1> analytic_deriv 
      = ew_derivs.calculate_ewsb_parameter_derivs(p);

   ewsb_condition_params pars = {&ew_derivs, p };

   gsl_function func_1 = { ewsb_condition_1, &pars };
   gsl_function func_2 = { ewsb_condition_2, &pars };
   gsl_function func_3 = { ewsb_condition_3, &pars };

   double numeric_deriv_1;
   double abs_err_1;

   gsl_deriv_central(&func_1, model.get_parameter(p), 1.0e-1, &numeric_deriv_1, &abs_err_1);

   BOOST_REQUIRE(Abs(abs_err_1) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(0) - numeric_deriv_1), abs_err_1);

   double numeric_deriv_2;
   double abs_err_2;

   gsl_deriv_central(&func_2, model.get_parameter(p), 1.0e-1, &numeric_deriv_2, &abs_err_2);

   BOOST_REQUIRE(Abs(abs_err_2) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(1) - numeric_deriv_2), abs_err_2);

   double numeric_deriv_3;
   double abs_err_3;

   gsl_deriv_central(&func_3, model.get_parameter(p), 1.0e-1, &numeric_deriv_3, &abs_err_3);

   BOOST_REQUIRE(Abs(abs_err_3) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(2) - numeric_deriv_3), abs_err_3);

}

BOOST_AUTO_TEST_CASE( test_tree_level_dewsb_conditions_dmHd2 )
{
   const double max_err = 0.01;
   const double theta = ArcTan(Sqrt(15.));
   const lowE6SSM_info::Parameters p = lowE6SSM_info::mHd2;

   lowE6SSM_input_parameters inputs = get_test_inputs(theta);

   lowE6SSM<Two_scale> model;

   initialize_model(inputs, model);

   model.set_ewsb_loop_order(0);

   lowE6SSM_ew_derivs ew_derivs(model);

   Eigen::Matrix<double,3,1> analytic_deriv 
      = ew_derivs.calculate_ewsb_parameter_derivs(p);

   ewsb_condition_params pars = {&ew_derivs, p };

   gsl_function func_1 = { ewsb_condition_1, &pars };
   gsl_function func_2 = { ewsb_condition_2, &pars };
   gsl_function func_3 = { ewsb_condition_3, &pars };

   double numeric_deriv_1;
   double abs_err_1;

   gsl_deriv_central(&func_1, model.get_parameter(p), 1.0e-1, &numeric_deriv_1, &abs_err_1);

   BOOST_REQUIRE(Abs(abs_err_1) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(0) - numeric_deriv_1), abs_err_1);

   double numeric_deriv_2;
   double abs_err_2;

   gsl_deriv_central(&func_2, model.get_parameter(p), 1.0e-1, &numeric_deriv_2, &abs_err_2);

   BOOST_REQUIRE(Abs(abs_err_2) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(1) - numeric_deriv_2), abs_err_2);

   double numeric_deriv_3;
   double abs_err_3;

   gsl_deriv_central(&func_3, model.get_parameter(p), 1.0e-1, &numeric_deriv_3, &abs_err_3);

   BOOST_REQUIRE(Abs(abs_err_3) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(2) - numeric_deriv_3), abs_err_3);

}

BOOST_AUTO_TEST_CASE( test_tree_level_dewsb_conditions_dmHu2 )
{
   const double max_err = 0.01;
   const double theta = ArcTan(Sqrt(15.));
   const lowE6SSM_info::Parameters p = lowE6SSM_info::mHu2;

   lowE6SSM_input_parameters inputs = get_test_inputs(theta);

   lowE6SSM<Two_scale> model;

   initialize_model(inputs, model);

   model.set_ewsb_loop_order(0);

   lowE6SSM_ew_derivs ew_derivs(model);

   Eigen::Matrix<double,3,1> analytic_deriv 
      = ew_derivs.calculate_ewsb_parameter_derivs(p);

   ewsb_condition_params pars = {&ew_derivs, p };

   gsl_function func_1 = { ewsb_condition_1, &pars };
   gsl_function func_2 = { ewsb_condition_2, &pars };
   gsl_function func_3 = { ewsb_condition_3, &pars };

   double numeric_deriv_1;
   double abs_err_1;

   gsl_deriv_central(&func_1, model.get_parameter(p), 1.0e-1, &numeric_deriv_1, &abs_err_1);

   BOOST_REQUIRE(Abs(abs_err_1) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(0) - numeric_deriv_1), abs_err_1);

   double numeric_deriv_2;
   double abs_err_2;

   gsl_deriv_central(&func_2, model.get_parameter(p), 1.0e-1, &numeric_deriv_2, &abs_err_2);

   BOOST_REQUIRE(Abs(abs_err_2) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(1) - numeric_deriv_2), abs_err_2);

   double numeric_deriv_3;
   double abs_err_3;

   gsl_deriv_central(&func_3, model.get_parameter(p), 1.0e-1, &numeric_deriv_3, &abs_err_3);

   BOOST_REQUIRE(Abs(abs_err_3) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(2) - numeric_deriv_3), abs_err_3);

}

BOOST_AUTO_TEST_CASE( test_tree_level_dewsb_conditions_dms2 )
{
   const double max_err = 0.01;
   const double theta = ArcTan(Sqrt(15.));
   const lowE6SSM_info::Parameters p = lowE6SSM_info::ms2;

   lowE6SSM_input_parameters inputs = get_test_inputs(theta);

   lowE6SSM<Two_scale> model;

   initialize_model(inputs, model);

   model.set_ewsb_loop_order(0);

   lowE6SSM_ew_derivs ew_derivs(model);

   Eigen::Matrix<double,3,1> analytic_deriv 
      = ew_derivs.calculate_ewsb_parameter_derivs(p);

   ewsb_condition_params pars = {&ew_derivs, p };

   gsl_function func_1 = { ewsb_condition_1, &pars };
   gsl_function func_2 = { ewsb_condition_2, &pars };
   gsl_function func_3 = { ewsb_condition_3, &pars };

   double numeric_deriv_1;
   double abs_err_1;

   gsl_deriv_central(&func_1, model.get_parameter(p), 1.0e-1, &numeric_deriv_1, &abs_err_1);

   BOOST_REQUIRE(Abs(abs_err_1) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(0) - numeric_deriv_1), abs_err_1);

   double numeric_deriv_2;
   double abs_err_2;

   gsl_deriv_central(&func_2, model.get_parameter(p), 1.0e-1, &numeric_deriv_2, &abs_err_2);

   BOOST_REQUIRE(Abs(abs_err_2) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(1) - numeric_deriv_2), abs_err_2);

   double numeric_deriv_3;
   double abs_err_3;

   gsl_deriv_central(&func_3, model.get_parameter(p), 1.0e-1, &numeric_deriv_3, &abs_err_3);

   BOOST_REQUIRE(Abs(abs_err_3) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(2) - numeric_deriv_3), abs_err_3);

}

BOOST_AUTO_TEST_CASE( test_tree_level_dewsb_conditions_dmq222 )
{
   const double max_err = 0.01;
   const double theta = ArcTan(Sqrt(15.));
   const lowE6SSM_info::Parameters p = lowE6SSM_info::mq222;

   lowE6SSM_input_parameters inputs = get_test_inputs(theta);

   lowE6SSM<Two_scale> model;

   initialize_model(inputs, model);

   model.set_ewsb_loop_order(0);

   lowE6SSM_ew_derivs ew_derivs(model);

   Eigen::Matrix<double,3,1> analytic_deriv 
      = ew_derivs.calculate_ewsb_parameter_derivs(p);

   ewsb_condition_params pars = {&ew_derivs, p };

   gsl_function func_1 = { ewsb_condition_1, &pars };
   gsl_function func_2 = { ewsb_condition_2, &pars };
   gsl_function func_3 = { ewsb_condition_3, &pars };

   double numeric_deriv_1;
   double abs_err_1;

   gsl_deriv_central(&func_1, model.get_parameter(p), 1.0e-1, &numeric_deriv_1, &abs_err_1);

   BOOST_REQUIRE(Abs(abs_err_1) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(0) - numeric_deriv_1), abs_err_1);

   double numeric_deriv_2;
   double abs_err_2;

   gsl_deriv_central(&func_2, model.get_parameter(p), 1.0e-1, &numeric_deriv_2, &abs_err_2);

   BOOST_REQUIRE(Abs(abs_err_2) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(1) - numeric_deriv_2), abs_err_2);

   double numeric_deriv_3;
   double abs_err_3;

   gsl_deriv_central(&func_3, model.get_parameter(p), 1.0e-1, &numeric_deriv_3, &abs_err_3);

   BOOST_REQUIRE(Abs(abs_err_3) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(2) - numeric_deriv_3), abs_err_3);

}

BOOST_AUTO_TEST_CASE( test_tree_level_dewsb_conditions_dmu222 )
{
   const double max_err = 0.01;
   const double theta = ArcTan(Sqrt(15.));
   const lowE6SSM_info::Parameters p = lowE6SSM_info::mu222;

   lowE6SSM_input_parameters inputs = get_test_inputs(theta);

   lowE6SSM<Two_scale> model;

   initialize_model(inputs, model);

   model.set_ewsb_loop_order(0);

   lowE6SSM_ew_derivs ew_derivs(model);

   Eigen::Matrix<double,3,1> analytic_deriv 
      = ew_derivs.calculate_ewsb_parameter_derivs(p);

   ewsb_condition_params pars = {&ew_derivs, p };

   gsl_function func_1 = { ewsb_condition_1, &pars };
   gsl_function func_2 = { ewsb_condition_2, &pars };
   gsl_function func_3 = { ewsb_condition_3, &pars };

   double numeric_deriv_1;
   double abs_err_1;

   gsl_deriv_central(&func_1, model.get_parameter(p), 1.0e-1, &numeric_deriv_1, &abs_err_1);

   BOOST_REQUIRE(Abs(abs_err_1) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(0) - numeric_deriv_1), abs_err_1);

   double numeric_deriv_2;
   double abs_err_2;

   gsl_deriv_central(&func_2, model.get_parameter(p), 1.0e-1, &numeric_deriv_2, &abs_err_2);

   BOOST_REQUIRE(Abs(abs_err_2) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(1) - numeric_deriv_2), abs_err_2);

   double numeric_deriv_3;
   double abs_err_3;

   gsl_deriv_central(&func_3, model.get_parameter(p), 1.0e-1, &numeric_deriv_3, &abs_err_3);

   BOOST_REQUIRE(Abs(abs_err_3) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(2) - numeric_deriv_3), abs_err_3);

}

BOOST_AUTO_TEST_CASE( test_tree_level_dewsb_conditions_dTYu22 )
{
   const double max_err = 0.01;
   const double theta = ArcTan(Sqrt(15.));
   const lowE6SSM_info::Parameters p = lowE6SSM_info::TYu22;

   lowE6SSM_input_parameters inputs = get_test_inputs(theta);

   lowE6SSM<Two_scale> model;

   initialize_model(inputs, model);

   model.set_ewsb_loop_order(0);

   lowE6SSM_ew_derivs ew_derivs(model);

   Eigen::Matrix<double,3,1> analytic_deriv 
      = ew_derivs.calculate_ewsb_parameter_derivs(p);

   ewsb_condition_params pars = {&ew_derivs, p };

   gsl_function func_1 = { ewsb_condition_1, &pars };
   gsl_function func_2 = { ewsb_condition_2, &pars };
   gsl_function func_3 = { ewsb_condition_3, &pars };

   double numeric_deriv_1;
   double abs_err_1;

   gsl_deriv_central(&func_1, model.get_parameter(p), 1.0e-1, &numeric_deriv_1, &abs_err_1);

   BOOST_REQUIRE(Abs(abs_err_1) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(0) - numeric_deriv_1), abs_err_1);

   double numeric_deriv_2;
   double abs_err_2;

   gsl_deriv_central(&func_2, model.get_parameter(p), 1.0e-1, &numeric_deriv_2, &abs_err_2);

   BOOST_REQUIRE(Abs(abs_err_2) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(1) - numeric_deriv_2), abs_err_2);

   double numeric_deriv_3;
   double abs_err_3;

   gsl_deriv_central(&func_3, model.get_parameter(p), 1.0e-1, &numeric_deriv_3, &abs_err_3);

   BOOST_REQUIRE(Abs(abs_err_3) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(2) - numeric_deriv_3), abs_err_3);

}

BOOST_AUTO_TEST_CASE( test_tree_level_dewsb_conditions_dYu22 )
{
   const double max_err = 0.01;
   const double theta = ArcTan(Sqrt(15.));
   const lowE6SSM_info::Parameters p = lowE6SSM_info::Yu22;

   lowE6SSM_input_parameters inputs = get_test_inputs(theta);

   lowE6SSM<Two_scale> model;

   initialize_model(inputs, model);

   model.set_ewsb_loop_order(0);

   lowE6SSM_ew_derivs ew_derivs(model);

   Eigen::Matrix<double,3,1> analytic_deriv 
      = ew_derivs.calculate_ewsb_parameter_derivs(p);

   ewsb_condition_params pars = {&ew_derivs, p };

   gsl_function func_1 = { ewsb_condition_1, &pars };
   gsl_function func_2 = { ewsb_condition_2, &pars };
   gsl_function func_3 = { ewsb_condition_3, &pars };

   double numeric_deriv_1;
   double abs_err_1;

   gsl_deriv_central(&func_1, model.get_parameter(p), 1.0e-1, &numeric_deriv_1, &abs_err_1);

   BOOST_REQUIRE(Abs(abs_err_1) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(0) - numeric_deriv_1), abs_err_1);

   double numeric_deriv_2;
   double abs_err_2;

   gsl_deriv_central(&func_2, model.get_parameter(p), 1.0e-1, &numeric_deriv_2, &abs_err_2);

   BOOST_REQUIRE(Abs(abs_err_2) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(1) - numeric_deriv_2), abs_err_2);

   double numeric_deriv_3;
   double abs_err_3;

   gsl_deriv_central(&func_3, model.get_parameter(p), 1.0e-1, &numeric_deriv_3, &abs_err_3);

   BOOST_REQUIRE(Abs(abs_err_3) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(2) - numeric_deriv_3), abs_err_3);

}

// As above, but at 1-loop order
BOOST_AUTO_TEST_CASE( test_one_loop_dewsb_conditions_dLambdax )
{
   const double max_err = 0.01;
   const double theta = ArcTan(Sqrt(15.));
   const lowE6SSM_info::Parameters p = lowE6SSM_info::Lambdax;

   lowE6SSM_input_parameters inputs = get_test_inputs(theta);

   lowE6SSM<Two_scale> model;

   initialize_model(inputs, model);

   model.set_ewsb_loop_order(1);

   lowE6SSM_ew_derivs ew_derivs(model);

   Eigen::Matrix<double,3,1> analytic_deriv 
      = ew_derivs.calculate_ewsb_parameter_derivs(p);

   ewsb_condition_params pars = {&ew_derivs, p };

   gsl_function func_1 = { ewsb_condition_1, &pars };
   gsl_function func_2 = { ewsb_condition_2, &pars };
   gsl_function func_3 = { ewsb_condition_3, &pars };

   double numeric_deriv_1;
   double abs_err_1;

   gsl_deriv_central(&func_1, model.get_parameter(p), 1.0e-1, &numeric_deriv_1, &abs_err_1);

   BOOST_REQUIRE(Abs(abs_err_1) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(0) - numeric_deriv_1), abs_err_1);

   double numeric_deriv_2;
   double abs_err_2;

   gsl_deriv_central(&func_2, model.get_parameter(p), 1.0e-1, &numeric_deriv_2, &abs_err_2);

   BOOST_REQUIRE(Abs(abs_err_2) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(1) - numeric_deriv_2), abs_err_2);

   double numeric_deriv_3;
   double abs_err_3;

   gsl_deriv_central(&func_3, model.get_parameter(p), 1.0e-1, &numeric_deriv_3, &abs_err_3);

   BOOST_REQUIRE(Abs(abs_err_3) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(2) - numeric_deriv_3), abs_err_3);

}

BOOST_AUTO_TEST_CASE( test_one_loop_dewsb_conditions_dg1 )
{
   const double max_err = 0.05;
   const double theta = ArcTan(Sqrt(15.));
   const lowE6SSM_info::Parameters p = lowE6SSM_info::g1;

   lowE6SSM_input_parameters inputs = get_test_inputs(theta);

   lowE6SSM<Two_scale> model;

   initialize_model(inputs, model);

   model.set_ewsb_loop_order(1);

   lowE6SSM_ew_derivs ew_derivs(model);

   Eigen::Matrix<double,3,1> analytic_deriv 
      = ew_derivs.calculate_ewsb_parameter_derivs(p);

   ewsb_condition_params pars = {&ew_derivs, p };

   gsl_function func_1 = { ewsb_condition_1, &pars };
   gsl_function func_2 = { ewsb_condition_2, &pars };
   gsl_function func_3 = { ewsb_condition_3, &pars };

   double numeric_deriv_1;
   double abs_err_1;

   gsl_deriv_central(&func_1, model.get_parameter(p), 1.0e-3, &numeric_deriv_1, &abs_err_1);

   BOOST_REQUIRE(Abs(abs_err_1) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(0) - numeric_deriv_1), abs_err_1);

   double numeric_deriv_2;
   double abs_err_2;

   gsl_deriv_central(&func_2, model.get_parameter(p), 1.0e-3, &numeric_deriv_2, &abs_err_2);

   BOOST_REQUIRE(Abs(abs_err_2) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(1) - numeric_deriv_2), abs_err_2);

   double numeric_deriv_3;
   double abs_err_3;

   gsl_deriv_central(&func_3, model.get_parameter(p), 1.0e-3, &numeric_deriv_3, &abs_err_3);

   BOOST_REQUIRE(Abs(abs_err_3) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(2) - numeric_deriv_3), abs_err_3);

}

BOOST_AUTO_TEST_CASE( test_one_loop_dewsb_conditions_dg2 )
{
   const double max_err = 0.01;
   const double theta = ArcTan(Sqrt(15.));
   const lowE6SSM_info::Parameters p = lowE6SSM_info::g2;

   lowE6SSM_input_parameters inputs = get_test_inputs(theta);

   lowE6SSM<Two_scale> model;

   initialize_model(inputs, model);

   model.set_ewsb_loop_order(1);

   lowE6SSM_ew_derivs ew_derivs(model);

   Eigen::Matrix<double,3,1> analytic_deriv 
      = ew_derivs.calculate_ewsb_parameter_derivs(p);

   ewsb_condition_params pars = {&ew_derivs, p };

   gsl_function func_1 = { ewsb_condition_1, &pars };
   gsl_function func_2 = { ewsb_condition_2, &pars };
   gsl_function func_3 = { ewsb_condition_3, &pars };

   double numeric_deriv_1;
   double abs_err_1;

   gsl_deriv_central(&func_1, model.get_parameter(p), 1.0e-4, &numeric_deriv_1, &abs_err_1);

   BOOST_REQUIRE(Abs(abs_err_1) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(0) - numeric_deriv_1), abs_err_1);

   double numeric_deriv_2;
   double abs_err_2;

   gsl_deriv_central(&func_2, model.get_parameter(p), 1.0e-4, &numeric_deriv_2, &abs_err_2);

   BOOST_REQUIRE(Abs(abs_err_2) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(1) - numeric_deriv_2), abs_err_2);

   double numeric_deriv_3;
   double abs_err_3;

   gsl_deriv_central(&func_3, model.get_parameter(p), 1.0e-1, &numeric_deriv_3, &abs_err_3);

   BOOST_REQUIRE(Abs(abs_err_3) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(2) - numeric_deriv_3), abs_err_3);

}

BOOST_AUTO_TEST_CASE( test_one_loop_dewsb_conditions_dgN )
{
   const double max_err = 20.;
   const double theta = ArcTan(Sqrt(15.));
   const lowE6SSM_info::Parameters p = lowE6SSM_info::gN;

   lowE6SSM_input_parameters inputs = get_test_inputs(theta);

   lowE6SSM<Two_scale> model;

   initialize_model(inputs, model);

   model.set_ewsb_loop_order(1);

   lowE6SSM_ew_derivs ew_derivs(model);

   Eigen::Matrix<double,3,1> analytic_deriv 
      = ew_derivs.calculate_ewsb_parameter_derivs(p);

   ewsb_condition_params pars = {&ew_derivs, p };

   gsl_function func_1 = { ewsb_condition_1, &pars };
   gsl_function func_2 = { ewsb_condition_2, &pars };
   gsl_function func_3 = { ewsb_condition_3, &pars };

   double numeric_deriv_1;
   double abs_err_1;

   gsl_deriv_central(&func_1, model.get_parameter(p), 1.1e-4, &numeric_deriv_1, &abs_err_1);

   BOOST_REQUIRE(Abs(abs_err_1) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(0) - numeric_deriv_1), abs_err_1);

   double numeric_deriv_2;
   double abs_err_2;

   gsl_deriv_central(&func_2, model.get_parameter(p), 1.0e-1, &numeric_deriv_2, &abs_err_2);

   BOOST_REQUIRE(Abs(abs_err_2) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(1) - numeric_deriv_2), abs_err_2);

   double numeric_deriv_3;
   double abs_err_3;

   gsl_deriv_central(&func_3, model.get_parameter(p), 1.0e-3, &numeric_deriv_3, &abs_err_3);

   BOOST_REQUIRE(Abs(abs_err_3) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(2) - numeric_deriv_3), abs_err_3);

}

BOOST_AUTO_TEST_CASE( test_one_loop_dewsb_conditions_dvd )
{
   const double max_err = 0.01;
   const double theta = ArcTan(Sqrt(15.));
   const lowE6SSM_info::Parameters p = lowE6SSM_info::vd;

   lowE6SSM_input_parameters inputs = get_test_inputs(theta);

   lowE6SSM<Two_scale> model;

   initialize_model(inputs, model);

   model.set_ewsb_loop_order(1);

   lowE6SSM_ew_derivs ew_derivs(model);

   Eigen::Matrix<double,3,1> analytic_deriv 
      = ew_derivs.calculate_ewsb_parameter_derivs(p);

   ewsb_condition_params pars = {&ew_derivs, p };

   gsl_function func_1 = { ewsb_condition_1, &pars };
   gsl_function func_2 = { ewsb_condition_2, &pars };
   gsl_function func_3 = { ewsb_condition_3, &pars };

   double numeric_deriv_1;
   double abs_err_1;

   gsl_deriv_central(&func_1, model.get_parameter(p), 1.0e-2, &numeric_deriv_1, &abs_err_1);

   BOOST_REQUIRE(Abs(abs_err_1) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(0) - numeric_deriv_1), abs_err_1);

   double numeric_deriv_2;
   double abs_err_2;

   gsl_deriv_central(&func_2, model.get_parameter(p), 1.0e-2, &numeric_deriv_2, &abs_err_2);

   BOOST_REQUIRE(Abs(abs_err_2) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(1) - numeric_deriv_2), abs_err_2);

   double numeric_deriv_3;
   double abs_err_3;

   gsl_deriv_central(&func_3, model.get_parameter(p), 1.0e-2, &numeric_deriv_3, &abs_err_3);

   BOOST_REQUIRE(Abs(abs_err_3) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(2) - numeric_deriv_3), abs_err_3);

}

BOOST_AUTO_TEST_CASE( test_one_loop_dewsb_conditions_dvu )
{
   const double max_err = 0.01;
   const double theta = ArcTan(Sqrt(15.));
   const lowE6SSM_info::Parameters p = lowE6SSM_info::vu;

   lowE6SSM_input_parameters inputs = get_test_inputs(theta);

   lowE6SSM<Two_scale> model;

   initialize_model(inputs, model);

   model.set_ewsb_loop_order(1);

   lowE6SSM_ew_derivs ew_derivs(model);

   Eigen::Matrix<double,3,1> analytic_deriv 
      = ew_derivs.calculate_ewsb_parameter_derivs(p);

   ewsb_condition_params pars = {&ew_derivs, p };

   gsl_function func_1 = { ewsb_condition_1, &pars };
   gsl_function func_2 = { ewsb_condition_2, &pars };
   gsl_function func_3 = { ewsb_condition_3, &pars };

   double numeric_deriv_1;
   double abs_err_1;

   gsl_deriv_central(&func_1, model.get_parameter(p), 1.0e-3, &numeric_deriv_1, &abs_err_1);

   BOOST_REQUIRE(Abs(abs_err_1) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(0) - numeric_deriv_1), abs_err_1);

   double numeric_deriv_2;
   double abs_err_2;

   gsl_deriv_central(&func_2, model.get_parameter(p), 1.0, &numeric_deriv_2, &abs_err_2);

   BOOST_REQUIRE(Abs(abs_err_2) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(1) - numeric_deriv_2), abs_err_2);

   double numeric_deriv_3;
   double abs_err_3;

   gsl_deriv_central(&func_3, model.get_parameter(p), 10.0, &numeric_deriv_3, &abs_err_3);

   BOOST_REQUIRE(Abs(abs_err_3) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(2) - numeric_deriv_3), abs_err_3);

}

BOOST_AUTO_TEST_CASE( test_one_loop_dewsb_conditions_dvs )
{
   const double max_err = 0.01;
   const double theta = ArcTan(Sqrt(15.));
   const lowE6SSM_info::Parameters p = lowE6SSM_info::vs;

   lowE6SSM_input_parameters inputs = get_test_inputs(theta);

   lowE6SSM<Two_scale> model;

   initialize_model(inputs, model);

   model.set_ewsb_loop_order(1);

   lowE6SSM_ew_derivs ew_derivs(model);

   Eigen::Matrix<double,3,1> analytic_deriv 
      = ew_derivs.calculate_ewsb_parameter_derivs(p);

   ewsb_condition_params pars = {&ew_derivs, p };

   gsl_function func_1 = { ewsb_condition_1, &pars };
   gsl_function func_2 = { ewsb_condition_2, &pars };
   gsl_function func_3 = { ewsb_condition_3, &pars };

   double numeric_deriv_1;
   double abs_err_1;

   gsl_deriv_central(&func_1, model.get_parameter(p), 1.0, &numeric_deriv_1, &abs_err_1);

   BOOST_REQUIRE(Abs(abs_err_1) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(0) - numeric_deriv_1), abs_err_1);

   double numeric_deriv_2;
   double abs_err_2;

   gsl_deriv_central(&func_2, model.get_parameter(p), 1.0e-1, &numeric_deriv_2, &abs_err_2);

   BOOST_REQUIRE(Abs(abs_err_2) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(1) - numeric_deriv_2), abs_err_2);

   double numeric_deriv_3;
   double abs_err_3;

   gsl_deriv_central(&func_3, model.get_parameter(p), 1.0e-1, &numeric_deriv_3, &abs_err_3);

   BOOST_REQUIRE(Abs(abs_err_3) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(2) - numeric_deriv_3), abs_err_3);

}

BOOST_AUTO_TEST_CASE( test_one_loop_dewsb_conditions_dTLambdax )
{
   const double max_err = 0.01;
   const double theta = ArcTan(Sqrt(15.));
   const lowE6SSM_info::Parameters p = lowE6SSM_info::TLambdax;

   lowE6SSM_input_parameters inputs = get_test_inputs(theta);

   lowE6SSM<Two_scale> model;

   initialize_model(inputs, model);

   model.set_ewsb_loop_order(1);

   lowE6SSM_ew_derivs ew_derivs(model);

   Eigen::Matrix<double,3,1> analytic_deriv 
      = ew_derivs.calculate_ewsb_parameter_derivs(p);

   ewsb_condition_params pars = {&ew_derivs, p };

   gsl_function func_1 = { ewsb_condition_1, &pars };
   gsl_function func_2 = { ewsb_condition_2, &pars };
   gsl_function func_3 = { ewsb_condition_3, &pars };

   double numeric_deriv_1;
   double abs_err_1;

   gsl_deriv_central(&func_1, model.get_parameter(p), 1.0e-1, &numeric_deriv_1, &abs_err_1);

   BOOST_REQUIRE(Abs(abs_err_1) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(0) - numeric_deriv_1), abs_err_1);

   double numeric_deriv_2;
   double abs_err_2;

   gsl_deriv_central(&func_2, model.get_parameter(p), 1.0e-1, &numeric_deriv_2, &abs_err_2);

   BOOST_REQUIRE(Abs(abs_err_2) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(1) - numeric_deriv_2), abs_err_2);

   double numeric_deriv_3;
   double abs_err_3;

   gsl_deriv_central(&func_3, model.get_parameter(p), 1.0e-1, &numeric_deriv_3, &abs_err_3);

   BOOST_REQUIRE(Abs(abs_err_3) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(2) - numeric_deriv_3), abs_err_3);

}

BOOST_AUTO_TEST_CASE( test_one_loop_dewsb_conditions_dmHd2 )
{
   const double max_err = 0.01;
   const double theta = ArcTan(Sqrt(15.));
   const lowE6SSM_info::Parameters p = lowE6SSM_info::mHd2;

   lowE6SSM_input_parameters inputs = get_test_inputs(theta);

   lowE6SSM<Two_scale> model;

   initialize_model(inputs, model);

   model.set_ewsb_loop_order(1);

   lowE6SSM_ew_derivs ew_derivs(model);

   Eigen::Matrix<double,3,1> analytic_deriv 
      = ew_derivs.calculate_ewsb_parameter_derivs(p);

   ewsb_condition_params pars = {&ew_derivs, p };

   gsl_function func_1 = { ewsb_condition_1, &pars };
   gsl_function func_2 = { ewsb_condition_2, &pars };
   gsl_function func_3 = { ewsb_condition_3, &pars };

   double numeric_deriv_1;
   double abs_err_1;

   gsl_deriv_central(&func_1, model.get_parameter(p), 1.0e-1, &numeric_deriv_1, &abs_err_1);

   BOOST_REQUIRE(Abs(abs_err_1) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(0) - numeric_deriv_1), abs_err_1);

   double numeric_deriv_2;
   double abs_err_2;

   gsl_deriv_central(&func_2, model.get_parameter(p), 1.0e-1, &numeric_deriv_2, &abs_err_2);

   BOOST_REQUIRE(Abs(abs_err_2) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(1) - numeric_deriv_2), abs_err_2);

   double numeric_deriv_3;
   double abs_err_3;

   gsl_deriv_central(&func_3, model.get_parameter(p), 1.0e-1, &numeric_deriv_3, &abs_err_3);

   BOOST_REQUIRE(Abs(abs_err_3) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(2) - numeric_deriv_3), abs_err_3);

}

BOOST_AUTO_TEST_CASE( test_one_loop_dewsb_conditions_dmHu2 )
{
   const double max_err = 0.01;
   const double theta = ArcTan(Sqrt(15.));
   const lowE6SSM_info::Parameters p = lowE6SSM_info::mHu2;

   lowE6SSM_input_parameters inputs = get_test_inputs(theta);

   lowE6SSM<Two_scale> model;

   initialize_model(inputs, model);

   model.set_ewsb_loop_order(1);

   lowE6SSM_ew_derivs ew_derivs(model);

   Eigen::Matrix<double,3,1> analytic_deriv 
      = ew_derivs.calculate_ewsb_parameter_derivs(p);

   ewsb_condition_params pars = {&ew_derivs, p };

   gsl_function func_1 = { ewsb_condition_1, &pars };
   gsl_function func_2 = { ewsb_condition_2, &pars };
   gsl_function func_3 = { ewsb_condition_3, &pars };

   double numeric_deriv_1;
   double abs_err_1;

   gsl_deriv_central(&func_1, model.get_parameter(p), 1.0e-1, &numeric_deriv_1, &abs_err_1);

   BOOST_REQUIRE(Abs(abs_err_1) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(0) - numeric_deriv_1), abs_err_1);

   double numeric_deriv_2;
   double abs_err_2;

   gsl_deriv_central(&func_2, model.get_parameter(p), 1.0e-1, &numeric_deriv_2, &abs_err_2);

   BOOST_REQUIRE(Abs(abs_err_2) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(1) - numeric_deriv_2), abs_err_2);

   double numeric_deriv_3;
   double abs_err_3;

   gsl_deriv_central(&func_3, model.get_parameter(p), 1.0e-1, &numeric_deriv_3, &abs_err_3);

   BOOST_REQUIRE(Abs(abs_err_3) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(2) - numeric_deriv_3), abs_err_3);

}

BOOST_AUTO_TEST_CASE( test_one_loop_dewsb_conditions_dms2 )
{
   const double max_err = 0.01;
   const double theta = ArcTan(Sqrt(15.));
   const lowE6SSM_info::Parameters p = lowE6SSM_info::ms2;

   lowE6SSM_input_parameters inputs = get_test_inputs(theta);

   lowE6SSM<Two_scale> model;

   initialize_model(inputs, model);

   model.set_ewsb_loop_order(1);

   lowE6SSM_ew_derivs ew_derivs(model);

   Eigen::Matrix<double,3,1> analytic_deriv 
      = ew_derivs.calculate_ewsb_parameter_derivs(p);

   ewsb_condition_params pars = {&ew_derivs, p };

   gsl_function func_1 = { ewsb_condition_1, &pars };
   gsl_function func_2 = { ewsb_condition_2, &pars };
   gsl_function func_3 = { ewsb_condition_3, &pars };

   double numeric_deriv_1;
   double abs_err_1;

   gsl_deriv_central(&func_1, model.get_parameter(p), 1.0e-1, &numeric_deriv_1, &abs_err_1);

   BOOST_REQUIRE(Abs(abs_err_1) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(0) - numeric_deriv_1), abs_err_1);

   double numeric_deriv_2;
   double abs_err_2;

   gsl_deriv_central(&func_2, model.get_parameter(p), 1.0e-1, &numeric_deriv_2, &abs_err_2);

   BOOST_REQUIRE(Abs(abs_err_2) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(1) - numeric_deriv_2), abs_err_2);

   double numeric_deriv_3;
   double abs_err_3;

   gsl_deriv_central(&func_3, model.get_parameter(p), 1.0e-1, &numeric_deriv_3, &abs_err_3);

   BOOST_REQUIRE(Abs(abs_err_3) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(2) - numeric_deriv_3), abs_err_3);

}

BOOST_AUTO_TEST_CASE( test_one_loop_dewsb_conditions_dmq222 )
{
   const double max_err = 0.01;
   const double theta = ArcTan(Sqrt(15.));
   const lowE6SSM_info::Parameters p = lowE6SSM_info::mq222;

   lowE6SSM_input_parameters inputs = get_test_inputs(theta);

   lowE6SSM<Two_scale> model;

   initialize_model(inputs, model);

   model.set_ewsb_loop_order(1);

   lowE6SSM_ew_derivs ew_derivs(model);

   Eigen::Matrix<double,3,1> analytic_deriv 
      = ew_derivs.calculate_ewsb_parameter_derivs(p);

   ewsb_condition_params pars = {&ew_derivs, p };

   gsl_function func_1 = { ewsb_condition_1, &pars };
   gsl_function func_2 = { ewsb_condition_2, &pars };
   gsl_function func_3 = { ewsb_condition_3, &pars };

   double numeric_deriv_1;
   double abs_err_1;

   gsl_deriv_central(&func_1, model.get_parameter(p), 1.0e-1, &numeric_deriv_1, &abs_err_1);

   BOOST_REQUIRE(Abs(abs_err_1) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(0) - numeric_deriv_1), abs_err_1);

   double numeric_deriv_2;
   double abs_err_2;

   gsl_deriv_central(&func_2, model.get_parameter(p), 1.0e-1, &numeric_deriv_2, &abs_err_2);

   BOOST_REQUIRE(Abs(abs_err_2) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(1) - numeric_deriv_2), abs_err_2);

   double numeric_deriv_3;
   double abs_err_3;

   gsl_deriv_central(&func_3, model.get_parameter(p), 1.0e-1, &numeric_deriv_3, &abs_err_3);

   BOOST_REQUIRE(Abs(abs_err_3) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(2) - numeric_deriv_3), abs_err_3);

}

BOOST_AUTO_TEST_CASE( test_one_loop_dewsb_conditions_dmu222 )
{
   const double max_err = 0.01;
   const double theta = ArcTan(Sqrt(15.));
   const lowE6SSM_info::Parameters p = lowE6SSM_info::mu222;

   lowE6SSM_input_parameters inputs = get_test_inputs(theta);

   lowE6SSM<Two_scale> model;

   initialize_model(inputs, model);

   model.set_ewsb_loop_order(1);

   lowE6SSM_ew_derivs ew_derivs(model);

   Eigen::Matrix<double,3,1> analytic_deriv 
      = ew_derivs.calculate_ewsb_parameter_derivs(p);

   ewsb_condition_params pars = {&ew_derivs, p };

   gsl_function func_1 = { ewsb_condition_1, &pars };
   gsl_function func_2 = { ewsb_condition_2, &pars };
   gsl_function func_3 = { ewsb_condition_3, &pars };

   double numeric_deriv_1;
   double abs_err_1;

   gsl_deriv_central(&func_1, model.get_parameter(p), 1.0e-1, &numeric_deriv_1, &abs_err_1);

   BOOST_REQUIRE(Abs(abs_err_1) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(0) - numeric_deriv_1), abs_err_1);

   double numeric_deriv_2;
   double abs_err_2;

   gsl_deriv_central(&func_2, model.get_parameter(p), 1.0e-1, &numeric_deriv_2, &abs_err_2);

   BOOST_REQUIRE(Abs(abs_err_2) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(1) - numeric_deriv_2), abs_err_2);

   double numeric_deriv_3;
   double abs_err_3;

   gsl_deriv_central(&func_3, model.get_parameter(p), 1.0e-2, &numeric_deriv_3, &abs_err_3);

   BOOST_REQUIRE(Abs(abs_err_3) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(2) - numeric_deriv_3), abs_err_3);

}

BOOST_AUTO_TEST_CASE( test_one_loop_dewsb_conditions_dTYu22 )
{
   const double max_err = 0.01;
   const double theta = ArcTan(Sqrt(15.));
   const lowE6SSM_info::Parameters p = lowE6SSM_info::TYu22;

   lowE6SSM_input_parameters inputs = get_test_inputs(theta);

   lowE6SSM<Two_scale> model;

   initialize_model(inputs, model);

   model.set_ewsb_loop_order(1);

   lowE6SSM_ew_derivs ew_derivs(model);

   Eigen::Matrix<double,3,1> analytic_deriv 
      = ew_derivs.calculate_ewsb_parameter_derivs(p);

   ewsb_condition_params pars = {&ew_derivs, p };

   gsl_function func_1 = { ewsb_condition_1, &pars };
   gsl_function func_2 = { ewsb_condition_2, &pars };
   gsl_function func_3 = { ewsb_condition_3, &pars };

   double numeric_deriv_1;
   double abs_err_1;

   gsl_deriv_central(&func_1, model.get_parameter(p), 1.0e-3, &numeric_deriv_1, &abs_err_1);

   BOOST_REQUIRE(Abs(abs_err_1) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(0) - numeric_deriv_1), abs_err_1);

   double numeric_deriv_2;
   double abs_err_2;

   gsl_deriv_central(&func_2, model.get_parameter(p), 1.0e-1, &numeric_deriv_2, &abs_err_2);

   BOOST_REQUIRE(Abs(abs_err_2) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(1) - numeric_deriv_2), abs_err_2);

   double numeric_deriv_3;
   double abs_err_3;

   gsl_deriv_central(&func_3, model.get_parameter(p), 1.0e-2, &numeric_deriv_3, &abs_err_3);

   BOOST_REQUIRE(Abs(abs_err_3) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(2) - numeric_deriv_3), abs_err_3);

}

BOOST_AUTO_TEST_CASE( test_one_loop_dewsb_conditions_dYu22 )
{
   const double max_err = 0.01;
   const double theta = ArcTan(Sqrt(15.));
   const lowE6SSM_info::Parameters p = lowE6SSM_info::Yu22;

   lowE6SSM_input_parameters inputs = get_test_inputs(theta);

   lowE6SSM<Two_scale> model;

   initialize_model(inputs, model);

   model.set_ewsb_loop_order(1);

   lowE6SSM_ew_derivs ew_derivs(model);

   Eigen::Matrix<double,3,1> analytic_deriv 
      = ew_derivs.calculate_ewsb_parameter_derivs(p);

   ewsb_condition_params pars = {&ew_derivs, p };

   gsl_function func_1 = { ewsb_condition_1, &pars };
   gsl_function func_2 = { ewsb_condition_2, &pars };
   gsl_function func_3 = { ewsb_condition_3, &pars };

   double numeric_deriv_1;
   double abs_err_1;

   gsl_deriv_central(&func_1, model.get_parameter(p), 0.5, &numeric_deriv_1, &abs_err_1);

   BOOST_REQUIRE(Abs(abs_err_1) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(0) - numeric_deriv_1), abs_err_1);

   double numeric_deriv_2;
   double abs_err_2;

   gsl_deriv_central(&func_2, model.get_parameter(p), 1.0e-2, &numeric_deriv_2, &abs_err_2);

   BOOST_REQUIRE(Abs(abs_err_2) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(1) - numeric_deriv_2), abs_err_2);

   double numeric_deriv_3;
   double abs_err_3;

   gsl_deriv_central(&func_3, model.get_parameter(p), 1.0e-3, &numeric_deriv_3, &abs_err_3);

   BOOST_REQUIRE(Abs(abs_err_3) < max_err);
   BOOST_CHECK_LE(Abs(analytic_deriv(2) - numeric_deriv_3), abs_err_3);

}
