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

// File generated at Sun 15 Jun 2014 19:15:46

#include "genericE6SSM_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of vd.
 *
 * @return one-loop beta function
 */
double genericE6SSM_susy_parameters::calc_beta_vd_one_loop(const Susy_traces& susy_traces) const
{
   const auto QH1p = INPUT(QH1p);
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   double beta_vd;

   beta_vd = 0.1*oneOver16PiSqr*vd*(-30*traceYdAdjYd - 10*traceYeAdjYe -
      10*AbsSqr(Lambdax) + 3*Sqr(g1) + 15*Sqr(g2) + 20*Sqr(gN)*Sqr(QH1p));


   return beta_vd;
}

/**
 * Calculates the two-loop beta function of vd.
 *
 * @return two-loop beta function
 */
double genericE6SSM_susy_parameters::calc_beta_vd_two_loop(const Susy_traces& susy_traces) const
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


   double beta_vd;

   beta_vd = -0.005*twoLoop*vd*(297*Power(g1,4) + 725*Power(g2,4) + 1600*
      Power(gN,4)*Power(QH1p,4) - 1800*traceYdAdjYdYdAdjYd - 600*
      traceYdAdjYuYuAdjYd - 600*traceYeAdjYeYeAdjYe + 300*traceYeAdjYe*Sqr(g1)
      + 300*traceYeAdjYe*Sqr(g2) + 90*Sqr(g1)*Sqr(g2) - 360*Qdp*QH1p*Sqr(g1)*
      Sqr(gN) - 360*QDxbarp*QH1p*Sqr(g1)*Sqr(gN) + 360*QDxp*QH1p*Sqr(g1)*Sqr(gN
      ) - 360*Qep*QH1p*Sqr(g1)*Sqr(gN) - 360*QH1p*QH2p*Sqr(g1)*Sqr(gN) - 120*
      QH1p*QHpbarp*Sqr(g1)*Sqr(gN) + 120*QH1p*QHpp*Sqr(g1)*Sqr(gN) + 360*QH1p*
      QLp*Sqr(g1)*Sqr(gN) - 360*QH1p*QQp*Sqr(g1)*Sqr(gN) + 720*QH1p*Qup*Sqr(g1)
      *Sqr(gN) + 400*traceYeAdjYe*Sqr(gN)*Sqr(Qep) + 480*Sqr(g1)*Sqr(gN)*Sqr(
      QH1p) + 600*Sqr(g2)*Sqr(gN)*Sqr(QH1p) + 1800*Power(gN,4)*Sqr(Qdp)*Sqr(
      QH1p) + 1800*Power(gN,4)*Sqr(QDxbarp)*Sqr(QH1p) + 1800*Power(gN,4)*Sqr(
      QDxp)*Sqr(QH1p) + 600*Power(gN,4)*Sqr(Qep)*Sqr(QH1p) + 1200*Power(gN,4)*
      Sqr(QH1p)*Sqr(QH2p) + 400*Power(gN,4)*Sqr(QH1p)*Sqr(QHpbarp) + 400*Power(
      gN,4)*Sqr(QH1p)*Sqr(QHpp) + 400*traceYeAdjYe*Sqr(gN)*Sqr(QLp) + 1200*
      Power(gN,4)*Sqr(QH1p)*Sqr(QLp) + 3600*Power(gN,4)*Sqr(QH1p)*Sqr(QQp) +
      100*traceYdAdjYd*(Sqr(g1) + 9*Sqr(g2) + 32*Sqr(g3) + 12*Sqr(gN)*Sqr(Qdp)
      + 12*Sqr(gN)*Sqr(QQp)) + 600*Power(gN,4)*Sqr(QH1p)*Sqr(QSp) + 20*AbsSqr(
      Lambdax)*(-30*traceKappaAdjKappa - 20*traceLambda12AdjLambda12 - 30*
      traceYuAdjYu + 3*Sqr(g1) + 15*Sqr(g2) + 20*Sqr(gN)*Sqr(QH2p) + 20*Sqr(gN)
      *Sqr(QSp)) + 1800*Power(gN,4)*Sqr(QH1p)*Sqr(Qup) - 600*Sqr(Conj(Lambdax))
      *Sqr(Lambdax));


   return beta_vd;
}

} // namespace flexiblesusy
