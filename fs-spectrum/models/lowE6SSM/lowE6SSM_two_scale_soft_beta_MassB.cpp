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

// File generated at Sun 24 Aug 2014 16:10:21

#include "lowE6SSM_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of MassB.
 *
 * @return one-loop beta function
 */
double lowE6SSM_soft_parameters::calc_beta_MassB_one_loop(const Soft_traces& soft_traces) const
{


   double beta_MassB;

   beta_MassB = 19.2*MassB*oneOver16PiSqr*Sqr(g1);


   return beta_MassB;
}

/**
 * Calculates the two-loop beta function of MassB.
 *
 * @return two-loop beta function
 */
double lowE6SSM_soft_parameters::calc_beta_MassB_two_loop(const Soft_traces& soft_traces) const
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
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceAdjKappaTKappa = TRACE_STRUCT.traceAdjKappaTKappa;
   const double traceAdjLambda12TLambda12 =
      TRACE_STRUCT.traceAdjLambda12TLambda12;


   double beta_MassB;

   beta_MassB = 0.16*twoLoop*Sqr(g1)*(10*traceAdjKappaTKappa + 15*
      traceAdjLambda12TLambda12 + 35*traceAdjYdTYd + 45*traceAdjYeTYe + 65*
      traceAdjYuTYu - 10*MassB*traceKappaAdjKappa - 15*MassB*
      traceLambda12AdjLambda12 - 35*MassB*traceYdAdjYd - 45*MassB*traceYeAdjYe
      - 65*MassB*traceYuAdjYu + 234*MassB*Sqr(g1) + 135*MassB*Sqr(g2) + 135*
      MassWB*Sqr(g2) + 300*MassB*Sqr(g3) + 300*MassG*Sqr(g3) + 30*MassB*Sqr(gN)
      *Sqr(Qdp) + 30*MassBp*Sqr(gN)*Sqr(Qdp) + 30*MassB*Sqr(gN)*Sqr(QDxbarp) +
      30*MassBp*Sqr(gN)*Sqr(QDxbarp) + 30*MassB*Sqr(gN)*Sqr(QDxp) + 30*MassBp*
      Sqr(gN)*Sqr(QDxp) + 90*MassB*Sqr(gN)*Sqr(Qep) + 90*MassBp*Sqr(gN)*Sqr(Qep
      ) + 45*MassB*Sqr(gN)*Sqr(QH1p) + 45*MassBp*Sqr(gN)*Sqr(QH1p) + 45*MassB*
      Sqr(gN)*Sqr(QH2p) + 45*MassBp*Sqr(gN)*Sqr(QH2p) + 15*MassB*Sqr(gN)*Sqr(
      QHpbarp) + 15*MassBp*Sqr(gN)*Sqr(QHpbarp) + 15*MassB*Sqr(gN)*Sqr(QHpp) +
      15*MassBp*Sqr(gN)*Sqr(QHpp) + 45*MassB*Sqr(gN)*Sqr(QLp) + 45*MassBp*Sqr(
      gN)*Sqr(QLp) + 15*MassB*Sqr(gN)*Sqr(QQp) + 15*MassBp*Sqr(gN)*Sqr(QQp) +
      120*MassB*Sqr(gN)*Sqr(Qup) + 120*MassBp*Sqr(gN)*Sqr(Qup) - 15*Conj(
      Lambdax)*(MassB*Lambdax - TLambdax));


   return beta_MassB;
}

} // namespace flexiblesusy
