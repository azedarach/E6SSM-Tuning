// ====================================================================
// Test suite for checking coefficients of various terms in the
// 2-loop beta for \lambda in the old E6SSM code
// ====================================================================

#include <limits>
#include <random>
#include <chrono>
#include <Eigen/Core>
#include <Eigen/Dense>

#include "wrappers.hpp"
#include "softsusy_essmsusy.h"
#include "softsusy_essmsoftpars.h"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_E6SpecGen_Lambdax_2lp_beta

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

using namespace std;
using namespace essmsoftsusy;

EssmSusy initialize_model()
{
   DoubleMatrix u(3,3);
   DoubleMatrix d(3,3);
   DoubleMatrix e(3,3);

   u(3,3) = 0.8;

   d(3,3) = 0.1;

   e(3,3) = 0.05;

   DoubleVector gauge(3);

   gauge(1) = 0.46;
   gauge(2) = 0.63;
   gauge(3) = 1.02;

   const double mu = 0.;
   const double tb = 10.;
   const double Q = 20000.;
   const int loops = 2;
   const int thresholds = 2;
   const double v = 246.;

   DoubleVector lambda(3);
   DoubleVector kappa(3);
   DoubleVector MN_SUSY(3);

   lambda(1) = 0.1;
   lambda(2) = 0.1;
   lambda(3) = 0.3;

   kappa(1) = 0.5;
   kappa(2) = 0.5;
   kappa(3) = 0.5;

   MN_SUSY(1) = 100000.;
   MN_SUSY(2) = 100000.;
   MN_SUSY(3) = 100000.;

   const double hN = 0.;
   const double mu_0 = 0.;
   const double gN = 0.46;
   const double g11 = 0.;

   return EssmSusy(u, d, e, gauge, mu, tb, Q, loops, thresholds,
                   v, lambda, kappa, MN_SUSY, hN, mu_0,
                   gN, g11);
}

BOOST_AUTO_TEST_CASE( test_2lp_beta_Lambdax_gN_Yd22_coeff )
{
   const double max_err = 1.0e-5;

   EssmSusy model(initialize_model());

   // ensure only non-zero terms involve gN, Yd(2,2) and
   // Lambdax
   DoubleMatrix zero_matrix(3,3);
   DoubleVector zero_vector(3);

   model.setYukawaMatrix(YU, zero_matrix);
   model.setYukawaMatrix(YD, zero_matrix);
   model.setYukawaMatrix(YE, zero_matrix);
   model.setAllGauge(zero_vector);
   model.setSusyMu(0.);
   model.seth_N(0.);
   model.setg_11(0.);
   model.setmu_0(0.);
   model.setkappa(zero_vector);
   model.setlambda(zero_vector);

   // use 6 combinations of (Lambdax, Yd(2,2), gN) to obtain the
   // coefficients of the various non-zero terms
   const std::size_t num_terms = 6;

   unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
   std::default_random_engine generator(seed);
   std::uniform_real_distribution<double> distribution(0., 1.);

   Eigen::Matrix<double, num_terms, num_terms> coefficients;
   Eigen::VectorXd rhs(num_terms);

   DoubleVector lambda(3);
   EssmSusy betas_1lp;
   EssmSusy betas_2lp;
   double lambdax_2lp_beta;

   for (std::size_t i = 0; i < num_terms; ++i) {
      lambda(3) = distribution(generator);
      model.setlambda(lambda);
      model.setYukawaElement(YD, 3, 3, distribution(generator));
      model.setgash_1(distribution(generator));

      model.setLoops(1);

      betas_1lp.set(model.beta());

      model.setLoops(2);

      betas_2lp.set(model.beta());

      lambdax_2lp_beta = (betas_2lp.displaylambda()(3) - betas_1lp.displaylambda()(3)) 
         / (flexiblesusy::twoLoop * model.displaylambda()(3));

      coefficients(i, 0) = flexiblesusy::Power(model.displaylambda()(3), 4);
      coefficients(i, 1) = flexiblesusy::Power(model.displaygdash_1(), 4);
      coefficients(i, 2) = flexiblesusy::Power(model.displayYukawaElement(YD, 3, 3), 4);
      coefficients(i, 3) = flexiblesusy::Sqr(model.displaylambda()(3)) * flexiblesusy::Sqr(model.displaygdash_1());
      coefficients(i, 4) = flexiblesusy::Sqr(model.displaylambda()(3)) * flexiblesusy::Sqr(model.displayYukawaElement(YD, 3, 3));
      coefficients(i, 5) = flexiblesusy::Sqr(model.displaygdash_1()) * flexiblesusy::Sqr(model.displayYukawaElement(YD, 3, 3));

      rhs(i) = lambdax_2lp_beta;
   }

   // solve the system
   Eigen::VectorXd x = coefficients.fullPivHouseholderQr().solve(rhs);

   // expected values (from arXiv:0904.2169v3 [hep-ph])
   const double lambdax4_coeff = -10.;
   const double gN4_coeff = 19.665;
   const double yb4_coeff = -9.;
   const double lambdax2gN2_coeff = 1.3;
   const double lambdax2yb2_coeff = -9.;
   const double gN2yb2_coeff = -0.2;

   BOOST_REQUIRE((coefficients * x).isApprox(rhs));
   BOOST_CHECK_LE(flexiblesusy::Abs(x(0) - lambdax4_coeff), max_err);
   BOOST_CHECK_LE(flexiblesusy::Abs(x(1) - gN4_coeff), max_err);
   BOOST_CHECK_LE(flexiblesusy::Abs(x(2) - yb4_coeff), max_err);
   BOOST_CHECK_LE(flexiblesusy::Abs(x(3) - lambdax2gN2_coeff), max_err);
   BOOST_CHECK_LE(flexiblesusy::Abs(x(4) - lambdax2yb2_coeff), max_err);
   BOOST_CHECK_LE(flexiblesusy::Abs(x(5) - gN2yb2_coeff), max_err);

}

BOOST_AUTO_TEST_CASE( test_2lp_beta_Lambdax_gN_g1_coeff )
{
   const double max_err = 1.0e-5;

   EssmSusy model(initialize_model());

   // ensure only non-zero terms involve gN, g1 and
   // Lambdax
   DoubleMatrix zero_matrix(3,3);
   DoubleVector zero_vector(3);

   model.setYukawaMatrix(YU, zero_matrix);
   model.setYukawaMatrix(YD, zero_matrix);
   model.setYukawaMatrix(YE, zero_matrix);
   model.setAllGauge(zero_vector);
   model.setSusyMu(0.);
   model.seth_N(0.);
   model.setg_11(0.);
   model.setmu_0(0.);
   model.setkappa(zero_vector);
   model.setlambda(zero_vector);

   // use 6 combinations of (Lambdax, g1, gN) to obtain the
   // coefficients of the various non-zero terms
   const std::size_t num_terms = 6;

   unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
   std::default_random_engine generator(seed);
   std::uniform_real_distribution<double> distribution(0., 1.);

   Eigen::Matrix<double, num_terms, num_terms> coefficients;
   Eigen::VectorXd rhs(num_terms);

   DoubleVector lambda(3);
   EssmSusy betas_1lp;
   EssmSusy betas_2lp;
   double lambdax_2lp_beta;

   for (std::size_t i = 0; i < num_terms; ++i) {
      lambda(3) = distribution(generator);
      model.setlambda(lambda);
      model.setGaugeCoupling(1, distribution(generator));
      model.setgash_1(distribution(generator));

      model.setLoops(1);

      betas_1lp.set(model.beta());

      model.setLoops(2);

      betas_2lp.set(model.beta());

      lambdax_2lp_beta = (betas_2lp.displaylambda()(3) - betas_1lp.displaylambda()(3)) 
         / (flexiblesusy::twoLoop * model.displaylambda()(3));

      coefficients(i, 0) = flexiblesusy::Power(model.displaylambda()(3), 4);
      coefficients(i, 1) = flexiblesusy::Power(model.displaygdash_1(), 4);
      coefficients(i, 2) = flexiblesusy::Power(model.displayGaugeCoupling(1), 4);
      coefficients(i, 3) = flexiblesusy::Sqr(model.displaylambda()(3)) * flexiblesusy::Sqr(model.displaygdash_1());
      coefficients(i, 4) = flexiblesusy::Sqr(model.displaylambda()(3)) * flexiblesusy::Sqr(model.displayGaugeCoupling(1));
      coefficients(i, 5) = flexiblesusy::Sqr(model.displaygdash_1()) * flexiblesusy::Sqr(model.displayGaugeCoupling(1));

      rhs(i) = lambdax_2lp_beta;
   }

   // solve the system
   Eigen::VectorXd x = coefficients.fullPivHouseholderQr().solve(rhs);

   // expected values (from arXiv:0904.2169v3 [hep-ph])
   const double lambdax4_coeff = -10.;
   const double gN4_coeff = 19.665;
   const double g14_coeff = 5.94;
   const double lambdax2gN2_coeff = 1.3;
   const double lambdax2g12_coeff = 1.2;
   const double gN2g12_coeff = 0.39;

   BOOST_REQUIRE((coefficients * x).isApprox(rhs));
   BOOST_CHECK_LE(flexiblesusy::Abs(x(0) - lambdax4_coeff), max_err);
   BOOST_CHECK_LE(flexiblesusy::Abs(x(1) - gN4_coeff), max_err);
   BOOST_CHECK_LE(flexiblesusy::Abs(x(2) - g14_coeff), max_err);
   BOOST_CHECK_LE(flexiblesusy::Abs(x(3) - lambdax2gN2_coeff), max_err);
   BOOST_CHECK_LE(flexiblesusy::Abs(x(4) - lambdax2g12_coeff), max_err);
   BOOST_CHECK_LE(flexiblesusy::Abs(x(5) - gN2g12_coeff), max_err);

}
