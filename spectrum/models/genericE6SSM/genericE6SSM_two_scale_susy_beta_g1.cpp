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
 * Calculates the one-loop beta function of g1.
 *
 * @return one-loop beta function
 */
double genericE6SSM_susy_parameters::calc_beta_g1_one_loop(const Susy_traces& susy_traces) const
{


   double beta_g1;

   beta_g1 = 9.6*Power(g1,3)*oneOver16PiSqr;


   return beta_g1;
}

/**
 * Calculates the two-loop beta function of g1.
 *
 * @return two-loop beta function
 */
double genericE6SSM_susy_parameters::calc_beta_g1_two_loop(const Susy_traces& susy_traces) const
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
   const auto Qup = INPUT(Qup);
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;


   double beta_g1;

   beta_g1 = 0.08*Power(g1,3)*twoLoop*(-10*traceKappaAdjKappa - 15*
      traceLambda12AdjLambda12 - 35*traceYdAdjYd - 45*traceYeAdjYe - 65*
      traceYuAdjYu - 15*AbsSqr(Lambdax) + 117*Sqr(g1) + 135*Sqr(g2) + 300*Sqr(
      g3) + 30*Sqr(gN)*Sqr(Qdp) + 30*Sqr(gN)*Sqr(QDxbarp) + 30*Sqr(gN)*Sqr(QDxp
      ) + 90*Sqr(gN)*Sqr(Qep) + 45*Sqr(gN)*Sqr(QH1p) + 45*Sqr(gN)*Sqr(QH2p) +
      15*Sqr(gN)*Sqr(QHpbarp) + 15*Sqr(gN)*Sqr(QHpp) + 45*Sqr(gN)*Sqr(QLp) + 15
      *Sqr(gN)*Sqr(QQp) + 120*Sqr(gN)*Sqr(Qup));


   return beta_g1;
}

} // namespace flexiblesusy
