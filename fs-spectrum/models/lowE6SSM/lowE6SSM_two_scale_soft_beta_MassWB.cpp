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

// File generated at Sun 24 Aug 2014 16:10:22

#include "lowE6SSM_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of MassWB.
 *
 * @return one-loop beta function
 */
double lowE6SSM_soft_parameters::calc_beta_MassWB_one_loop(const Soft_traces& soft_traces) const
{


   double beta_MassWB;

   beta_MassWB = 8*MassWB*oneOver16PiSqr*Sqr(g2);


   return beta_MassWB;
}

/**
 * Calculates the two-loop beta function of MassWB.
 *
 * @return two-loop beta function
 */
double lowE6SSM_soft_parameters::calc_beta_MassWB_two_loop(const Soft_traces& soft_traces) const
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
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceAdjLambda12TLambda12 =
      TRACE_STRUCT.traceAdjLambda12TLambda12;


   double beta_MassWB;

   beta_MassWB = 0.8*twoLoop*Sqr(g2)*(5*traceAdjLambda12TLambda12 + 15*
      traceAdjYdTYd + 5*traceAdjYeTYe + 15*traceAdjYuTYu - 5*MassWB*
      traceLambda12AdjLambda12 - 15*MassWB*traceYdAdjYd - 5*MassWB*traceYeAdjYe
      - 15*MassWB*traceYuAdjYu + 9*MassB*Sqr(g1) + 9*MassWB*Sqr(g1) + 230*
      MassWB*Sqr(g2) + 60*MassG*Sqr(g3) + 60*MassWB*Sqr(g3) + 15*MassBp*Sqr(gN)
      *Sqr(QH1p) + 15*MassWB*Sqr(gN)*Sqr(QH1p) + 15*MassBp*Sqr(gN)*Sqr(QH2p) +
      15*MassWB*Sqr(gN)*Sqr(QH2p) + 5*MassBp*Sqr(gN)*Sqr(QHpbarp) + 5*MassWB*
      Sqr(gN)*Sqr(QHpbarp) + 5*MassBp*Sqr(gN)*Sqr(QHpp) + 5*MassWB*Sqr(gN)*Sqr(
      QHpp) + 15*MassBp*Sqr(gN)*Sqr(QLp) + 15*MassWB*Sqr(gN)*Sqr(QLp) + 45*
      MassBp*Sqr(gN)*Sqr(QQp) + 45*MassWB*Sqr(gN)*Sqr(QQp) + Conj(Lambdax)*(-5*
      MassWB*Lambdax + 5*TLambdax));


   return beta_MassWB;
}

} // namespace flexiblesusy
