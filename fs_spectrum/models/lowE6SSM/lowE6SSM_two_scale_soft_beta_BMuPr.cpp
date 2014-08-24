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

// File generated at Sun 24 Aug 2014 16:10:07

#include "lowE6SSM_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of BMuPr.
 *
 * @return one-loop beta function
 */
double lowE6SSM_soft_parameters::calc_beta_BMuPr_one_loop(const Soft_traces& soft_traces) const
{
   const auto QHpbarp = INPUT(QHpbarp);
   const auto QHpp = INPUT(QHpp);


   double beta_BMuPr;

   beta_BMuPr = oneOver16PiSqr*(1.2*MassB*MuPr*Sqr(g1) + 6*MassWB*MuPr*
      Sqr(g2) + 4*MassBp*MuPr*Sqr(gN)*Sqr(QHpbarp) + 4*MassBp*MuPr*Sqr(gN)*Sqr(
      QHpp) - 0.2*BMuPr*(3*Sqr(g1) + 5*(3*Sqr(g2) + 2*Sqr(gN)*(Sqr(QHpbarp) +
      Sqr(QHpp)))));


   return beta_BMuPr;
}

/**
 * Calculates the two-loop beta function of BMuPr.
 *
 * @return two-loop beta function
 */
double lowE6SSM_soft_parameters::calc_beta_BMuPr_two_loop(const Soft_traces& soft_traces) const
{
   const auto QHpbarp = INPUT(QHpbarp);
   const auto QHpp = INPUT(QHpp);
   const auto QDxp = INPUT(QDxp);
   const auto Qep = INPUT(Qep);
   const auto QH1p = INPUT(QH1p);
   const auto QH2p = INPUT(QH2p);
   const auto Qdp = INPUT(Qdp);
   const auto QDxbarp = INPUT(QDxbarp);
   const auto QLp = INPUT(QLp);
   const auto QQp = INPUT(QQp);
   const auto QSp = INPUT(QSp);
   const auto Qup = INPUT(Qup);


   double beta_BMuPr;

   beta_BMuPr = 0.02*twoLoop*(BMuPr*(297*Power(g1,4) + 30*Sqr(g1)*(3*Sqr(
      g2) + 2*Sqr(gN)*(Sqr(QHpbarp) + Sqr(QHpp))) + 25*(33*Power(g2,4) + 12*Sqr
      (g2)*Sqr(gN)*(Sqr(QHpbarp) + Sqr(QHpp)) + 4*Power(gN,4)*(4*Power(QHpbarp,
      4) + 4*Power(QHpp,4) + 9*Sqr(QDxp)*Sqr(QHpbarp) + 3*Sqr(Qep)*Sqr(QHpbarp)
      + 6*Sqr(QH1p)*Sqr(QHpbarp) + 6*Sqr(QH2p)*Sqr(QHpbarp) + 9*Sqr(QDxp)*Sqr(
      QHpp) + 3*Sqr(Qep)*Sqr(QHpp) + 6*Sqr(QH1p)*Sqr(QHpp) + 6*Sqr(QH2p)*Sqr(
      QHpp) + 4*Sqr(QHpbarp)*Sqr(QHpp) + 9*Sqr(Qdp)*(Sqr(QHpbarp) + Sqr(QHpp))
      + 9*Sqr(QDxbarp)*(Sqr(QHpbarp) + Sqr(QHpp)) + 6*Sqr(QHpbarp)*Sqr(QLp) + 6
      *Sqr(QHpp)*Sqr(QLp) + 18*Sqr(QHpbarp)*Sqr(QQp) + 18*Sqr(QHpp)*Sqr(QQp) +
      3*Sqr(QHpbarp)*Sqr(QSp) + 3*Sqr(QHpp)*Sqr(QSp) + 9*Sqr(QHpbarp)*Sqr(Qup)
      + 9*Sqr(QHpp)*Sqr(Qup)))) - 4*MuPr*(297*Power(g1,4)*MassB + 15*Sqr(g1)*(3
      *(MassB + MassWB)*Sqr(g2) + 2*(MassB + MassBp)*Sqr(gN)*(Sqr(QHpbarp) +
      Sqr(QHpp))) + 25*(33*Power(g2,4)*MassWB + 6*(MassBp + MassWB)*Sqr(g2)*Sqr
      (gN)*(Sqr(QHpbarp) + Sqr(QHpp)) + 4*Power(gN,4)*MassBp*(4*Power(QHpbarp,4
      ) + 4*Power(QHpp,4) + 9*Sqr(QDxp)*Sqr(QHpbarp) + 3*Sqr(Qep)*Sqr(QHpbarp)
      + 6*Sqr(QH1p)*Sqr(QHpbarp) + 6*Sqr(QH2p)*Sqr(QHpbarp) + 9*Sqr(QDxp)*Sqr(
      QHpp) + 3*Sqr(Qep)*Sqr(QHpp) + 6*Sqr(QH1p)*Sqr(QHpp) + 6*Sqr(QH2p)*Sqr(
      QHpp) + 4*Sqr(QHpbarp)*Sqr(QHpp) + 9*Sqr(Qdp)*(Sqr(QHpbarp) + Sqr(QHpp))
      + 9*Sqr(QDxbarp)*(Sqr(QHpbarp) + Sqr(QHpp)) + 6*Sqr(QHpbarp)*Sqr(QLp) + 6
      *Sqr(QHpp)*Sqr(QLp) + 18*Sqr(QHpbarp)*Sqr(QQp) + 18*Sqr(QHpp)*Sqr(QQp) +
      3*Sqr(QHpbarp)*Sqr(QSp) + 3*Sqr(QHpp)*Sqr(QSp) + 9*Sqr(QHpbarp)*Sqr(Qup)
      + 9*Sqr(QHpp)*Sqr(Qup)))));


   return beta_BMuPr;
}

} // namespace flexiblesusy
