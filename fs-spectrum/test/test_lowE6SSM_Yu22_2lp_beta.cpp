// ====================================================================
// Test suite for checking coefficients of various terms in the
// 2-loop betas
// ====================================================================

#include <limits>
#include <random>
#include <chrono>
#include <Eigen/Core>

#include "wrappers.hpp"
#include "lowE6SSM_two_scale_model.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_lowE6SSM_Yu22_2lp_beta

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

   input.LambdaxInput = 0.3;

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
   input.MassGInput = 2000.0;
   input.MassBpInput = 300.0;

   return input;
}

void set_susy_parameters_from_input(const lowE6SSM_input_parameters& input, 
                                    lowE6SSM<Two_scale>& model)
{
   model.set_input_parameters(input);

   model.set_Yd(0, 0, 0.);
   model.set_Yd(1, 1, 0.);
   model.set_Yd(2, 2, 0.1);

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
}

BOOST_AUTO_TEST_CASE( test_2lp_beta_Yu22_g1_gN_coeff )
{
   const double max_err = 1.0e-5;
   const double theta = ArcTan(Sqrt(15.));

   lowE6SSM_input_parameters inputs = get_test_inputs(theta);

   lowE6SSM<Two_scale> model;

   initialize_model(inputs, model);

   // ensure only non-zero terms involve g1, gN and Yu22
   std::vector<lowE6SSM_info::Parameters> parameters_to_zero
      = {lowE6SSM_info::Yd00, lowE6SSM_info::Yd01, lowE6SSM_info::Yd02
         , lowE6SSM_info::Yd10, lowE6SSM_info::Yd11, lowE6SSM_info::Yd12
         , lowE6SSM_info::Yd20, lowE6SSM_info::Yd21, lowE6SSM_info::Yd22, 
         lowE6SSM_info::Ye00, lowE6SSM_info::Ye01, lowE6SSM_info::Ye02, 
         lowE6SSM_info::Ye10, lowE6SSM_info::Ye11, lowE6SSM_info::Ye12, 
         lowE6SSM_info::Ye20, lowE6SSM_info::Ye21, lowE6SSM_info::Ye22, 
         lowE6SSM_info::Yu00, lowE6SSM_info::Yu01, lowE6SSM_info::Yu02, 
         lowE6SSM_info::Yu10, lowE6SSM_info::Yu11, lowE6SSM_info::Yu12, 
         lowE6SSM_info::Yu20, lowE6SSM_info::Yu21, lowE6SSM_info::Kappa00
         , lowE6SSM_info::Kappa01, lowE6SSM_info::Kappa02, lowE6SSM_info::Kappa10
         , lowE6SSM_info::Kappa11, lowE6SSM_info::Kappa12, lowE6SSM_info::Kappa20
         , lowE6SSM_info::Kappa21, lowE6SSM_info::Kappa22, lowE6SSM_info::Lambda1200
         , lowE6SSM_info::Lambda1201, lowE6SSM_info::Lambda1210, lowE6SSM_info::Lambda1211
         , lowE6SSM_info::MuPr, lowE6SSM_info::Lambdax, lowE6SSM_info::g2
         , lowE6SSM_info::g3, lowE6SSM_info::vd, lowE6SSM_info::vu
         , lowE6SSM_info::vs};

   for (std::size_t i = 0; i < parameters_to_zero.size(); ++i) {
      model.set_parameter(parameters_to_zero.at(i), 0.);
   }

   // use 6 combinations of (Yu22, g1, gN) to obtain the 
   // coefficients of the various non-zero terms
   const std::size_t num_terms = 6;

   unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
   std::default_random_engine generator(seed);
   std::uniform_real_distribution<double> distribution(0., 1.);

   lowE6SSM_soft_parameters betas_1lp;
   lowE6SSM_soft_parameters betas_2lp;
   double yu22_2lp_beta;

   Eigen::Matrix<double, num_terms, num_terms> coefficients;
   Eigen::VectorXd rhs(num_terms);

   for (std::size_t i = 0; i < num_terms; ++i) {
      model.set_Yu(2, 2, distribution(generator));
      model.set_gN(distribution(generator));
      model.set_g1(distribution(generator));

      model.set_loops(1);

      betas_1lp = model.calc_beta();

      model.set_loops(2);

      betas_2lp = model.calc_beta();

      yu22_2lp_beta = (betas_2lp.get_Yu(2,2) - betas_1lp.get_Yu(2,2))
         / (twoLoop * model.get_Yu(2,2));

      // add row to matrix of values
      coefficients(i, 0) = Power(model.get_Yu(2,2), 4);
      coefficients(i, 1) = Power(model.get_gN(), 4);
      coefficients(i, 2) = Power(model.get_g1(), 4);
      coefficients(i, 3) = Sqr(model.get_Yu(2,2)) * Sqr(model.get_gN());
      coefficients(i, 4) = Sqr(model.get_Yu(2,2)) * Sqr(model.get_g1());
      coefficients(i, 5) = Sqr(model.get_gN()) * Sqr(model.get_g1());

      rhs(i) = yu22_2lp_beta;
   }

   // solve the system
   Eigen::VectorXd x = coefficients.fullPivHouseholderQr().solve(rhs);

   // expected values
   const double charges_sum = 6.0 * (Sqr(inputs.QH1p) + Sqr(inputs.QH2p) + Sqr(inputs.QLp))
      + 9.0 * (Sqr(inputs.Qdp) + Sqr(inputs.QDxp) + Sqr(inputs.QDxbarp) + Sqr(inputs.Qup)
               + 2.0 * Sqr(inputs.QQp))
      + 3.0 * (Sqr(inputs.Qep) + Sqr(inputs.QSp)) + 2.0 * (Sqr(inputs.QHpp) + Sqr(inputs.QHpbarp));

   const double mixed_charges_sum = 3.0 * (-inputs.QH1p + inputs.QH2p + inputs.Qdp - inputs.QDxp
                                           + inputs.Qep - inputs.QLp + inputs.QQp - 2.0 * inputs.Qup
                                           + inputs.QDxbarp) + inputs.QHpbarp - inputs.QHpp;

   const double yt4_coeff = -22.;
   const double g14_coeff = 3913.0 / 450.0;
   const double gN4_coeff = 2.0 * (2.0 * (Power(inputs.QH2p, 4) + Power(inputs.QQp, 4) + Power(inputs.Qup, 4)) + charges_sum * (Sqr(inputs.QH2p) + Sqr(inputs.QQp) + Sqr(inputs.Qup)));
   const double yt2gN2_coeff = 4.0 * (2.0 * Sqr(inputs.QQp) + Sqr(inputs.Qup));
   const double yt2g12_coeff = 1.2;
   const double gN2g12_coeff = 0.4 * (3.0 * Sqr(inputs.QH2p) + Sqr(inputs.QQp) / 3.0 + 16.0 * Sqr(inputs.Qup) / 3.0 + mixed_charges_sum * (3.0 * inputs.QH2p + inputs.QQp - 4.0 * inputs.Qup));

   BOOST_REQUIRE((coefficients * x).isApprox(rhs));
   BOOST_CHECK_LE(Abs(x(0) - yt4_coeff), max_err);
   BOOST_CHECK_LE(Abs(x(1) - gN4_coeff), max_err);
   BOOST_CHECK_LE(Abs(x(2) - g14_coeff), max_err);
   BOOST_CHECK_LE(Abs(x(3) - yt2gN2_coeff), max_err);
   BOOST_CHECK_LE(Abs(x(4) - yt2g12_coeff), max_err);
   BOOST_CHECK_LE(Abs(x(5) - gN2g12_coeff), max_err);

}
