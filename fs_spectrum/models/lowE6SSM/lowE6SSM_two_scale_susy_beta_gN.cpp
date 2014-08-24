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

// File generated at Sun 24 Aug 2014 16:08:37

#include "lowE6SSM_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of gN.
 *
 * @return one-loop beta function
 */
double lowE6SSM_susy_parameters::calc_beta_gN_one_loop(const Susy_traces& susy_traces) const
{
   const auto Qdp = INPUT(Qdp);
   const auto QDxbarp = INPUT(QDxbarp);
   const auto QDxp = INPUT(QDxp);
   const auto Qep = INPUT(Qep);
   const auto QH1p = INPUT(QH1p);
   const auto QH2p = INPUT(QH2p);
   const auto QHpbarp = INPUT(QHpbarp);
   const auto QHpp = INPUT(QHpp);
   const auto QLp = INPUT(QLp);
   const auto QQp = INPUT(QQp);
   const auto QSp = INPUT(QSp);
   const auto Qup = INPUT(Qup);


   double beta_gN;

   beta_gN = Power(gN,3)*oneOver16PiSqr*(9*Sqr(Qdp) + 9*Sqr(QDxbarp) + 9*
      Sqr(QDxp) + 3*Sqr(Qep) + 6*Sqr(QH1p) + 6*Sqr(QH2p) + 2*Sqr(QHpbarp) + 2*
      Sqr(QHpp) + 6*Sqr(QLp) + 18*Sqr(QQp) + 3*Sqr(QSp) + 9*Sqr(Qup));


   return beta_gN;
}

/**
 * Calculates the two-loop beta function of gN.
 *
 * @return two-loop beta function
 */
double lowE6SSM_susy_parameters::calc_beta_gN_two_loop(const Susy_traces& susy_traces) const
{
   const auto Qdp = INPUT(Qdp);
   const auto QDxbarp = INPUT(QDxbarp);
   const auto QDxp = INPUT(QDxp);
   const auto Qep = INPUT(Qep);
   const auto QH1p = INPUT(QH1p);
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


   double beta_gN;

   beta_gN = 0.4*Power(gN,3)*twoLoop*(90*Power(Qdp,4)*Sqr(gN) + 90*Power(
      QDxbarp,4)*Sqr(gN) + 90*Power(QDxp,4)*Sqr(gN) + 30*Power(Qep,4)*Sqr(gN) +
      60*Power(QH1p,4)*Sqr(gN) + 60*Power(QH2p,4)*Sqr(gN) + 20*Power(QHpbarp,4
      )*Sqr(gN) + 20*Power(QHpp,4)*Sqr(gN) + 60*Power(QLp,4)*Sqr(gN) + 180*
      Power(QQp,4)*Sqr(gN) + 30*Power(QSp,4)*Sqr(gN) + 90*Power(Qup,4)*Sqr(gN)
      + 6*Sqr(g1)*Sqr(Qdp) + 120*Sqr(g3)*Sqr(Qdp) - 15*traceKappaAdjKappa*Sqr(
      QDxbarp) + 6*Sqr(g1)*Sqr(QDxbarp) + 120*Sqr(g3)*Sqr(QDxbarp) - 15*
      traceKappaAdjKappa*Sqr(QDxp) + 6*Sqr(g1)*Sqr(QDxp) + 120*Sqr(g3)*Sqr(QDxp
      ) - 10*traceYeAdjYe*Sqr(Qep) + 18*Sqr(g1)*Sqr(Qep) - 10*
      traceLambda12AdjLambda12*Sqr(QH1p) - 10*traceYeAdjYe*Sqr(QH1p) + 9*Sqr(g1
      )*Sqr(QH1p) + 45*Sqr(g2)*Sqr(QH1p) - 10*traceLambda12AdjLambda12*Sqr(QH2p
      ) - 30*traceYuAdjYu*Sqr(QH2p) + 9*Sqr(g1)*Sqr(QH2p) + 45*Sqr(g2)*Sqr(QH2p
      ) + 3*Sqr(g1)*Sqr(QHpbarp) + 15*Sqr(g2)*Sqr(QHpbarp) + 3*Sqr(g1)*Sqr(QHpp
      ) + 15*Sqr(g2)*Sqr(QHpp) - 10*traceYeAdjYe*Sqr(QLp) + 9*Sqr(g1)*Sqr(QLp)
      + 45*Sqr(g2)*Sqr(QLp) - 30*traceYuAdjYu*Sqr(QQp) + 3*Sqr(g1)*Sqr(QQp) +
      135*Sqr(g2)*Sqr(QQp) + 240*Sqr(g3)*Sqr(QQp) - 30*traceYdAdjYd*(Sqr(Qdp) +
      Sqr(QH1p) + Sqr(QQp)) - 15*traceKappaAdjKappa*Sqr(QSp) - 10*
      traceLambda12AdjLambda12*Sqr(QSp) - 10*AbsSqr(Lambdax)*(Sqr(QH1p) + Sqr(
      QH2p) + Sqr(QSp)) - 30*traceYuAdjYu*Sqr(Qup) + 24*Sqr(g1)*Sqr(Qup) + 120*
      Sqr(g3)*Sqr(Qup));


   return beta_gN;
}

} // namespace flexiblesusy
