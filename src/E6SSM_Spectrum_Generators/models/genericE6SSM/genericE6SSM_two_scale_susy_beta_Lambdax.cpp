// ====================================================================
// This file is part of FlexibleSUSY.
//
// FlexibleSUSY is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// FlexibleSUSY is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with FlexibleSUSY.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

// File generated at Sun 15 Jun 2014 19:15:45

#include "genericE6SSM_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of Lambdax.
 *
 * @return one-loop beta function
 */
double genericE6SSM_susy_parameters::calc_beta_Lambdax_one_loop(const Susy_traces& susy_traces) const
{
   const auto QH1p = INPUT(QH1p);
   const auto QH2p = INPUT(QH2p);
   const auto QSp = INPUT(QSp);
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;


   double beta_Lambdax;

   beta_Lambdax = oneOver16PiSqr*(3*traceKappaAdjKappa*Lambdax + 2*
      traceLambda12AdjLambda12*Lambdax + 3*traceYdAdjYd*Lambdax + traceYeAdjYe*
      Lambdax + 3*traceYuAdjYu*Lambdax - 0.6*Lambdax*Sqr(g1) - 3*Lambdax*Sqr(g2
      ) - 2*Lambdax*Sqr(gN)*Sqr(QH1p) - 2*Lambdax*Sqr(gN)*Sqr(QH2p) - 2*Lambdax
      *Sqr(gN)*Sqr(QSp) + 4*Conj(Lambdax)*Sqr(Lambdax));


   return beta_Lambdax;
}

/**
 * Calculates the two-loop beta function of Lambdax.
 *
 * @return two-loop beta function
 */
double genericE6SSM_susy_parameters::calc_beta_Lambdax_two_loop(const Susy_traces& susy_traces) const
{
   const auto Qdp = INPUT(Qdp);
   const auto QH1p = INPUT(QH1p);
   const auto QDxbarp = INPUT(QDxbarp);
   const auto QDxp = INPUT(QDxp);
   const auto Qep = INPUT(Qep);
   const auto QH2p = INPUT(QH2p);
   const auto QHpbarp = INPUT(QHpbarp);
   const auto QHpp = INPUT(QHpp);
   const auto QLp = INPUT(QLp);
   const auto QQp = INPUT(QQp);
   const auto QSp = INPUT(QSp);
   const auto Qup = INPUT(Qup);
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;
   const double traceKappaAdjKappaKappaAdjKappa =
      TRACE_STRUCT.traceKappaAdjKappaKappaAdjKappa;
   const double traceLambda12AdjLambda12Lambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12Lambda12AdjLambda12;


   double beta_Lambdax;

   beta_Lambdax = -0.02*twoLoop*Lambdax*(-297*Power(g1,4) - 825*Power(g2,
      4) - 800*Power(gN,4)*Power(QH1p,4) - 800*Power(gN,4)*Power(QH2p,4) - 500*
      Power(gN,4)*Power(QSp,4) + 300*traceKappaAdjKappaKappaAdjKappa + 200*
      traceLambda12AdjLambda12Lambda12AdjLambda12 + 450*traceYdAdjYdYdAdjYd +
      300*traceYdAdjYuYuAdjYd + 150*traceYeAdjYeYeAdjYe + 450*
      traceYuAdjYuYuAdjYu - 40*traceKappaAdjKappa*Sqr(g1) - 60*
      traceLambda12AdjLambda12*Sqr(g1) - 60*traceYeAdjYe*Sqr(g1) - 40*
      traceYuAdjYu*Sqr(g1) - 300*traceLambda12AdjLambda12*Sqr(g2) - 90*Sqr(g1)*
      Sqr(g2) - 800*traceKappaAdjKappa*Sqr(g3) - 800*traceYuAdjYu*Sqr(g3) + 180
      *Qdp*QH1p*Sqr(g1)*Sqr(gN) + 180*QDxbarp*QH1p*Sqr(g1)*Sqr(gN) - 180*QDxp*
      QH1p*Sqr(g1)*Sqr(gN) + 180*Qep*QH1p*Sqr(g1)*Sqr(gN) - 180*Qdp*QH2p*Sqr(g1
      )*Sqr(gN) - 180*QDxbarp*QH2p*Sqr(g1)*Sqr(gN) + 180*QDxp*QH2p*Sqr(g1)*Sqr(
      gN) - 180*Qep*QH2p*Sqr(g1)*Sqr(gN) + 360*QH1p*QH2p*Sqr(g1)*Sqr(gN) + 60*
      QH1p*QHpbarp*Sqr(g1)*Sqr(gN) - 60*QH2p*QHpbarp*Sqr(g1)*Sqr(gN) - 60*QH1p*
      QHpp*Sqr(g1)*Sqr(gN) + 60*QH2p*QHpp*Sqr(g1)*Sqr(gN) - 180*QH1p*QLp*Sqr(g1
      )*Sqr(gN) + 180*QH2p*QLp*Sqr(g1)*Sqr(gN) + 180*QH1p*QQp*Sqr(g1)*Sqr(gN) -
      180*QH2p*QQp*Sqr(g1)*Sqr(gN) - 360*QH1p*Qup*Sqr(g1)*Sqr(gN) + 360*QH2p*
      Qup*Sqr(g1)*Sqr(gN) - 300*traceKappaAdjKappa*Sqr(gN)*Sqr(QDxbarp) - 300*
      traceKappaAdjKappa*Sqr(gN)*Sqr(QDxp) - 100*traceYeAdjYe*Sqr(gN)*Sqr(Qep)
      - 200*traceLambda12AdjLambda12*Sqr(gN)*Sqr(QH1p) + 100*traceYeAdjYe*Sqr(
      gN)*Sqr(QH1p) - 240*Sqr(g1)*Sqr(gN)*Sqr(QH1p) - 300*Sqr(g2)*Sqr(gN)*Sqr(
      QH1p) - 900*Power(gN,4)*Sqr(Qdp)*Sqr(QH1p) - 900*Power(gN,4)*Sqr(QDxbarp)
      *Sqr(QH1p) - 900*Power(gN,4)*Sqr(QDxp)*Sqr(QH1p) - 300*Power(gN,4)*Sqr(
      Qep)*Sqr(QH1p) - 200*traceLambda12AdjLambda12*Sqr(gN)*Sqr(QH2p) + 300*
      traceYuAdjYu*Sqr(gN)*Sqr(QH2p) - 240*Sqr(g1)*Sqr(gN)*Sqr(QH2p) - 300*Sqr(
      g2)*Sqr(gN)*Sqr(QH2p) - 900*Power(gN,4)*Sqr(Qdp)*Sqr(QH2p) - 900*Power(gN
      ,4)*Sqr(QDxbarp)*Sqr(QH2p) - 900*Power(gN,4)*Sqr(QDxp)*Sqr(QH2p) - 300*
      Power(gN,4)*Sqr(Qep)*Sqr(QH2p) - 1200*Power(gN,4)*Sqr(QH1p)*Sqr(QH2p) -
      10*AbsSqr(Lambdax)*(-30*traceKappaAdjKappa - 20*traceLambda12AdjLambda12
      - 45*traceYdAdjYd - 15*traceYeAdjYe - 45*traceYuAdjYu + 6*Sqr(g1) + 30*
      Sqr(g2) + 20*Sqr(gN)*Sqr(QH1p) + 20*Sqr(gN)*Sqr(QH2p)) - 200*Power(gN,4)*
      Sqr(QH1p)*Sqr(QHpbarp) - 200*Power(gN,4)*Sqr(QH2p)*Sqr(QHpbarp) - 200*
      Power(gN,4)*Sqr(QH1p)*Sqr(QHpp) - 200*Power(gN,4)*Sqr(QH2p)*Sqr(QHpp) -
      100*traceYeAdjYe*Sqr(gN)*Sqr(QLp) - 600*Power(gN,4)*Sqr(QH1p)*Sqr(QLp) -
      600*Power(gN,4)*Sqr(QH2p)*Sqr(QLp) - 300*traceYuAdjYu*Sqr(gN)*Sqr(QQp) -
      1800*Power(gN,4)*Sqr(QH1p)*Sqr(QQp) - 1800*Power(gN,4)*Sqr(QH2p)*Sqr(QQp)
      + 20*traceYdAdjYd*(Sqr(g1) - 5*(8*Sqr(g3) + 3*Sqr(gN)*(Sqr(Qdp) - Sqr(
      QH1p) + Sqr(QQp)))) + 300*traceKappaAdjKappa*Sqr(gN)*Sqr(QSp) + 200*
      traceLambda12AdjLambda12*Sqr(gN)*Sqr(QSp) - 900*Power(gN,4)*Sqr(Qdp)*Sqr(
      QSp) - 900*Power(gN,4)*Sqr(QDxbarp)*Sqr(QSp) - 900*Power(gN,4)*Sqr(QDxp)*
      Sqr(QSp) - 300*Power(gN,4)*Sqr(Qep)*Sqr(QSp) - 900*Power(gN,4)*Sqr(QH1p)*
      Sqr(QSp) - 900*Power(gN,4)*Sqr(QH2p)*Sqr(QSp) - 200*Power(gN,4)*Sqr(
      QHpbarp)*Sqr(QSp) - 200*Power(gN,4)*Sqr(QHpp)*Sqr(QSp) - 600*Power(gN,4)*
      Sqr(QLp)*Sqr(QSp) - 1800*Power(gN,4)*Sqr(QQp)*Sqr(QSp) - 300*traceYuAdjYu
      *Sqr(gN)*Sqr(Qup) - 900*Power(gN,4)*Sqr(QH1p)*Sqr(Qup) - 900*Power(gN,4)*
      Sqr(QH2p)*Sqr(Qup) - 900*Power(gN,4)*Sqr(QSp)*Sqr(Qup) + 500*Sqr(Conj(
      Lambdax))*Sqr(Lambdax));


   return beta_Lambdax;
}

} // namespace flexiblesusy
