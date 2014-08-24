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

// File generated at Sun 24 Aug 2014 16:08:38

#include "lowE6SSM_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of vs.
 *
 * @return one-loop beta function
 */
double lowE6SSM_susy_parameters::calc_beta_vs_one_loop(const Susy_traces& susy_traces) const
{
   const auto QSp = INPUT(QSp);
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;


   double beta_vs;

   beta_vs = oneOver16PiSqr*vs*(-3*traceKappaAdjKappa - 2*
      traceLambda12AdjLambda12 - 2*AbsSqr(Lambdax) + 2*Sqr(gN)*Sqr(QSp));


   return beta_vs;
}

/**
 * Calculates the two-loop beta function of vs.
 *
 * @return two-loop beta function
 */
double lowE6SSM_susy_parameters::calc_beta_vs_two_loop(const Susy_traces& susy_traces) const
{
   const auto Qdp = INPUT(Qdp);
   const auto QSp = INPUT(QSp);
   const auto QDxbarp = INPUT(QDxbarp);
   const auto QDxp = INPUT(QDxp);
   const auto Qep = INPUT(Qep);
   const auto QH1p = INPUT(QH1p);
   const auto QH2p = INPUT(QH2p);
   const auto QHpbarp = INPUT(QHpbarp);
   const auto QHpp = INPUT(QHpp);
   const auto QLp = INPUT(QLp);
   const auto QQp = INPUT(QQp);
   const auto Qup = INPUT(Qup);
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;
   const double traceKappaAdjKappaKappaAdjKappa =
      TRACE_STRUCT.traceKappaAdjKappaKappaAdjKappa;
   const double traceLambda12AdjLambda12Lambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12Lambda12AdjLambda12;


   double beta_vs;

   beta_vs = -0.2*twoLoop*vs*(25*Power(gN,4)*Power(QSp,4) - 30*
      traceKappaAdjKappaKappaAdjKappa - 20*
      traceLambda12AdjLambda12Lambda12AdjLambda12 + 6*traceLambda12AdjLambda12*
      Sqr(g1) + 30*traceLambda12AdjLambda12*Sqr(g2) + traceKappaAdjKappa*(4*Sqr
      (g1) + 10*(8*Sqr(g3) + 3*Sqr(gN)*(Sqr(QDxbarp) + Sqr(QDxp)))) + 20*
      traceLambda12AdjLambda12*Sqr(gN)*Sqr(QH1p) + 20*traceLambda12AdjLambda12*
      Sqr(gN)*Sqr(QH2p) + 2*AbsSqr(Lambdax)*(-15*traceYdAdjYd - 5*traceYeAdjYe
      - 15*traceYuAdjYu + 3*Sqr(g1) + 15*Sqr(g2) + 10*Sqr(gN)*Sqr(QH1p) + 10*
      Sqr(gN)*Sqr(QH2p)) + 45*Power(gN,4)*Sqr(Qdp)*Sqr(QSp) + 45*Power(gN,4)*
      Sqr(QDxbarp)*Sqr(QSp) + 45*Power(gN,4)*Sqr(QDxp)*Sqr(QSp) + 15*Power(gN,4
      )*Sqr(Qep)*Sqr(QSp) + 30*Power(gN,4)*Sqr(QH1p)*Sqr(QSp) + 30*Power(gN,4)*
      Sqr(QH2p)*Sqr(QSp) + 10*Power(gN,4)*Sqr(QHpbarp)*Sqr(QSp) + 10*Power(gN,4
      )*Sqr(QHpp)*Sqr(QSp) + 30*Power(gN,4)*Sqr(QLp)*Sqr(QSp) + 90*Power(gN,4)*
      Sqr(QQp)*Sqr(QSp) + 45*Power(gN,4)*Sqr(QSp)*Sqr(Qup) - 20*Sqr(Conj(
      Lambdax))*Sqr(Lambdax));


   return beta_vs;
}

} // namespace flexiblesusy
