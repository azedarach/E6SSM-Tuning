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
 * Calculates the one-loop beta function of g2.
 *
 * @return one-loop beta function
 */
double lowE6SSM_susy_parameters::calc_beta_g2_one_loop(const Susy_traces& susy_traces) const
{


   double beta_g2;

   beta_g2 = 4*Power(g2,3)*oneOver16PiSqr;


   return beta_g2;
}

/**
 * Calculates the two-loop beta function of g2.
 *
 * @return two-loop beta function
 */
double lowE6SSM_susy_parameters::calc_beta_g2_two_loop(const Susy_traces& susy_traces) const
{
   const auto QH1p = INPUT(QH1p);
   const auto QH2p = INPUT(QH2p);
   const auto QHpbarp = INPUT(QHpbarp);
   const auto QHpp = INPUT(QHpp);
   const auto QLp = INPUT(QLp);
   const auto QQp = INPUT(QQp);
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;


   double beta_g2;

   beta_g2 = 0.4*Power(g2,3)*twoLoop*(-5*traceLambda12AdjLambda12 - 15*
      traceYdAdjYd - 5*traceYeAdjYe - 15*traceYuAdjYu - 5*AbsSqr(Lambdax) + 9*
      Sqr(g1) + 115*Sqr(g2) + 60*Sqr(g3) + 15*Sqr(gN)*Sqr(QH1p) + 15*Sqr(gN)*
      Sqr(QH2p) + 5*Sqr(gN)*Sqr(QHpbarp) + 5*Sqr(gN)*Sqr(QHpp) + 15*Sqr(gN)*Sqr
      (QLp) + 45*Sqr(gN)*Sqr(QQp));


   return beta_g2;
}

} // namespace flexiblesusy
