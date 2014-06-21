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

// File generated at Sun 15 Jun 2014 19:16:44

#include "genericE6SSM_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of MassG.
 *
 * @return one-loop beta function
 */
double genericE6SSM_soft_parameters::calc_beta_MassG_one_loop(const Soft_traces& soft_traces) const
{


   double beta_MassG;

   beta_MassG = 0;


   return beta_MassG;
}

/**
 * Calculates the two-loop beta function of MassG.
 *
 * @return two-loop beta function
 */
double genericE6SSM_soft_parameters::calc_beta_MassG_two_loop(const Soft_traces& soft_traces) const
{
   const auto Qdp = INPUT(Qdp);
   const auto QDxbarp = INPUT(QDxbarp);
   const auto QDxp = INPUT(QDxp);
   const auto QQp = INPUT(QQp);
   const auto Qup = INPUT(Qup);
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceAdjKappaTKappa = TRACE_STRUCT.traceAdjKappaTKappa;


   double beta_MassG;

   beta_MassG = 2*twoLoop*Sqr(g3)*(2*traceAdjKappaTKappa + 4*
      traceAdjYdTYd + 4*traceAdjYuTYu - 2*MassG*traceKappaAdjKappa - 4*MassG*
      traceYdAdjYd - 4*MassG*traceYuAdjYu + 3*MassB*Sqr(g1) + 3*MassG*Sqr(g1) +
      9*MassG*Sqr(g2) + 9*MassWB*Sqr(g2) + 96*MassG*Sqr(g3) + 6*MassBp*Sqr(gN)
      *Sqr(Qdp) + 6*MassG*Sqr(gN)*Sqr(Qdp) + 6*MassBp*Sqr(gN)*Sqr(QDxbarp) + 6*
      MassG*Sqr(gN)*Sqr(QDxbarp) + 6*MassBp*Sqr(gN)*Sqr(QDxp) + 6*MassG*Sqr(gN)
      *Sqr(QDxp) + 12*MassBp*Sqr(gN)*Sqr(QQp) + 12*MassG*Sqr(gN)*Sqr(QQp) + 6*
      MassBp*Sqr(gN)*Sqr(Qup) + 6*MassG*Sqr(gN)*Sqr(Qup));


   return beta_MassG;
}

} // namespace flexiblesusy
