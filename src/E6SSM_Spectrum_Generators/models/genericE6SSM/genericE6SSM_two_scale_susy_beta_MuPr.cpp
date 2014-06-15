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

// File generated at Sun 15 Jun 2014 19:15:45

#include "genericE6SSM_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of MuPr.
 *
 * @return one-loop beta function
 */
double genericE6SSM_susy_parameters::calc_beta_MuPr_one_loop(const Susy_traces& susy_traces) const
{
   const auto QHpbarp = INPUT(QHpbarp);
   const auto QHpp = INPUT(QHpp);


   double beta_MuPr;

   beta_MuPr = -0.2*oneOver16PiSqr*MuPr*(3*Sqr(g1) + 5*(3*Sqr(g2) + 2*Sqr
      (gN)*(Sqr(QHpbarp) + Sqr(QHpp))));


   return beta_MuPr;
}

/**
 * Calculates the two-loop beta function of MuPr.
 *
 * @return two-loop beta function
 */
double genericE6SSM_susy_parameters::calc_beta_MuPr_two_loop(const Susy_traces& susy_traces) const
{
   const auto QDxbarp = INPUT(QDxbarp);
   const auto QHpbarp = INPUT(QHpbarp);
   const auto QDxp = INPUT(QDxp);
   const auto Qep = INPUT(Qep);
   const auto QH1p = INPUT(QH1p);
   const auto QH2p = INPUT(QH2p);
   const auto Qdp = INPUT(Qdp);
   const auto QHpp = INPUT(QHpp);
   const auto QLp = INPUT(QLp);
   const auto QQp = INPUT(QQp);
   const auto Qup = INPUT(Qup);
   const auto QSp = INPUT(QSp);


   double beta_MuPr;

   beta_MuPr = 0.02*twoLoop*MuPr*(297*Power(g1,4) + 30*Sqr(g1)*(3*Sqr(g2)
      + 2*Sqr(gN)*(3*QDxbarp*QHpbarp - 3*QDxp*QHpbarp + 3*Qep*QHpbarp - 3*QH1p
      *QHpbarp + 3*QH2p*QHpbarp + 3*Qdp*(QHpbarp - QHpp) - 3*QDxbarp*QHpp + 3*
      QDxp*QHpp - 3*Qep*QHpp + 3*QH1p*QHpp - 3*QH2p*QHpp - 2*QHpbarp*QHpp - 3*
      QHpbarp*QLp + 3*QHpp*QLp + 3*QHpbarp*QQp - 3*QHpp*QQp - 6*QHpbarp*Qup + 6
      *QHpp*Qup + 2*Sqr(QHpbarp) + 2*Sqr(QHpp))) + 25*(33*Power(g2,4) + 12*Sqr(
      g2)*Sqr(gN)*(Sqr(QHpbarp) + Sqr(QHpp)) + 4*Power(gN,4)*(4*Power(QHpbarp,4
      ) + 4*Power(QHpp,4) + 9*Sqr(QDxp)*Sqr(QHpbarp) + 3*Sqr(Qep)*Sqr(QHpbarp)
      + 6*Sqr(QH1p)*Sqr(QHpbarp) + 6*Sqr(QH2p)*Sqr(QHpbarp) + 9*Sqr(QDxp)*Sqr(
      QHpp) + 3*Sqr(Qep)*Sqr(QHpp) + 6*Sqr(QH1p)*Sqr(QHpp) + 6*Sqr(QH2p)*Sqr(
      QHpp) + 4*Sqr(QHpbarp)*Sqr(QHpp) + 9*Sqr(Qdp)*(Sqr(QHpbarp) + Sqr(QHpp))
      + 9*Sqr(QDxbarp)*(Sqr(QHpbarp) + Sqr(QHpp)) + 6*Sqr(QHpbarp)*Sqr(QLp) + 6
      *Sqr(QHpp)*Sqr(QLp) + 18*Sqr(QHpbarp)*Sqr(QQp) + 18*Sqr(QHpp)*Sqr(QQp) +
      3*Sqr(QHpbarp)*Sqr(QSp) + 3*Sqr(QHpp)*Sqr(QSp) + 9*Sqr(QHpbarp)*Sqr(Qup)
      + 9*Sqr(QHpp)*Sqr(Qup))));


   return beta_MuPr;
}

} // namespace flexiblesusy
