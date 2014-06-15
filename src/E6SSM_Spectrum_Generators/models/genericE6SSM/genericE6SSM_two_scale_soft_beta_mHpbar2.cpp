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

// File generated at Sun 15 Jun 2014 19:16:43

#include "genericE6SSM_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of mHpbar2.
 *
 * @return one-loop beta function
 */
double genericE6SSM_soft_parameters::calc_beta_mHpbar2_one_loop(const Soft_traces& soft_traces) const
{
   const auto QHpbarp = INPUT(QHpbarp);
   const double Tr11 = TRACE_STRUCT.Tr11;
   const double Tr14 = TRACE_STRUCT.Tr14;


   double beta_mHpbar2;

   beta_mHpbar2 = oneOver16PiSqr*(0.7745966692414834*g1*Tr11 + 2*gN*
      QHpbarp*Tr14 - 1.2*AbsSqr(MassB)*Sqr(g1) - 6*AbsSqr(MassWB)*Sqr(g2) - 8*
      AbsSqr(MassBp)*Sqr(gN)*Sqr(QHpbarp));


   return beta_mHpbar2;
}

/**
 * Calculates the two-loop beta function of mHpbar2.
 *
 * @return two-loop beta function
 */
double genericE6SSM_soft_parameters::calc_beta_mHpbar2_two_loop(const Soft_traces& soft_traces) const
{
   const auto QHpbarp = INPUT(QHpbarp);
   const auto Qdp = INPUT(Qdp);
   const auto QDxbarp = INPUT(QDxbarp);
   const auto QDxp = INPUT(QDxp);
   const auto Qep = INPUT(Qep);
   const auto QH1p = INPUT(QH1p);
   const auto QH2p = INPUT(QH2p);
   const auto QHpp = INPUT(QHpp);
   const auto QLp = INPUT(QLp);
   const auto QQp = INPUT(QQp);
   const auto Qup = INPUT(Qup);
   const auto QSp = INPUT(QSp);
   const double Tr2U111 = TRACE_STRUCT.Tr2U111;
   const double Tr2U114 = TRACE_STRUCT.Tr2U114;
   const double Tr31 = TRACE_STRUCT.Tr31;
   const double Tr22 = TRACE_STRUCT.Tr22;
   const double Tr2U141 = TRACE_STRUCT.Tr2U141;
   const double Tr2U144 = TRACE_STRUCT.Tr2U144;
   const double Tr34 = TRACE_STRUCT.Tr34;


   double beta_mHpbar2;

   beta_mHpbar2 = 0.04*twoLoop*(3*Conj(MassB)*Sqr(g1)*(297*MassB*Sqr(g1)
      + 5*(3*(2*MassB + MassWB)*Sqr(g2) + 4*(2*MassB + MassBp)*QHpbarp*(3*Qdp +
      3*QDxbarp - 3*QDxp + 3*Qep - 3*QH1p + 3*QH2p + 2*QHpbarp - QHpp - 3*QLp
      + 3*QQp - 6*Qup)*Sqr(gN))) + 5*(30*Power(g2,4)*Tr22 + 15.491933384829668*
      g1*gN*QHpbarp*Tr2U114 + 15.491933384829668*g1*gN*QHpbarp*Tr2U141 +
      15.491933384829668*g1*Tr31 + 40*gN*QHpbarp*Tr34 + 435*Power(g2,4)*AbsSqr(
      MassWB) + 6*Tr2U111*Sqr(g1) + 18*AbsSqr(MassWB)*Sqr(g1)*Sqr(g2) + 9*MassB
      *Conj(MassWB)*Sqr(g1)*Sqr(g2) + 40*Tr2U144*Sqr(gN)*Sqr(QHpbarp) + 120*
      AbsSqr(MassWB)*Sqr(g2)*Sqr(gN)*Sqr(QHpbarp) + 60*MassBp*Conj(MassWB)*Sqr(
      g2)*Sqr(gN)*Sqr(QHpbarp) + 12*QHpbarp*Conj(MassBp)*Sqr(gN)*((MassB + 2*
      MassBp)*(3*Qdp + 3*QDxbarp - 3*QDxp + 3*Qep - 3*QH1p + 3*QH2p + 2*QHpbarp
      - QHpp - 3*QLp + 3*QQp - 6*Qup)*Sqr(g1) + 5*QHpbarp*((2*MassBp + MassWB)
      *Sqr(g2) + 2*MassBp*Sqr(gN)*(9*Sqr(Qdp) + 9*Sqr(QDxbarp) + 9*Sqr(QDxp) +
      3*Sqr(Qep) + 6*Sqr(QH1p) + 6*Sqr(QH2p) + 4*Sqr(QHpbarp) + 2*Sqr(QHpp) + 6
      *Sqr(QLp) + 18*Sqr(QQp) + 3*Sqr(QSp) + 9*Sqr(Qup))))));


   return beta_mHpbar2;
}

} // namespace flexiblesusy
