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
 * Calculates the one-loop beta function of MassBp.
 *
 * @return one-loop beta function
 */
double genericE6SSM_soft_parameters::calc_beta_MassBp_one_loop(const Soft_traces& soft_traces) const
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


   double beta_MassBp;

   beta_MassBp = 2*MassBp*oneOver16PiSqr*Sqr(gN)*(9*Sqr(Qdp) + 9*Sqr(
      QDxbarp) + 9*Sqr(QDxp) + 3*Sqr(Qep) + 6*Sqr(QH1p) + 6*Sqr(QH2p) + 2*Sqr(
      QHpbarp) + 2*Sqr(QHpp) + 6*Sqr(QLp) + 18*Sqr(QQp) + 3*Sqr(QSp) + 9*Sqr(
      Qup));


   return beta_MassBp;
}

/**
 * Calculates the two-loop beta function of MassBp.
 *
 * @return two-loop beta function
 */
double genericE6SSM_soft_parameters::calc_beta_MassBp_two_loop(const Soft_traces& soft_traces) const
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
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceAdjKappaTKappa = TRACE_STRUCT.traceAdjKappaTKappa;
   const double traceAdjLambda12TLambda12 =
      TRACE_STRUCT.traceAdjLambda12TLambda12;


   double beta_MassBp;

   beta_MassBp = 0.8*twoLoop*Sqr(gN)*(180*MassBp*Power(Qdp,4)*Sqr(gN) +
      180*MassBp*Power(QDxbarp,4)*Sqr(gN) + 180*MassBp*Power(QDxp,4)*Sqr(gN) +
      60*MassBp*Power(Qep,4)*Sqr(gN) + 120*MassBp*Power(QH1p,4)*Sqr(gN) + 120*
      MassBp*Power(QH2p,4)*Sqr(gN) + 40*MassBp*Power(QHpbarp,4)*Sqr(gN) + 40*
      MassBp*Power(QHpp,4)*Sqr(gN) + 120*MassBp*Power(QLp,4)*Sqr(gN) + 360*
      MassBp*Power(QQp,4)*Sqr(gN) + 60*MassBp*Power(QSp,4)*Sqr(gN) + 180*MassBp
      *Power(Qup,4)*Sqr(gN) + 30*traceAdjYdTYd*Sqr(Qdp) + 6*MassB*Sqr(g1)*Sqr(
      Qdp) + 6*MassBp*Sqr(g1)*Sqr(Qdp) + 120*MassBp*Sqr(g3)*Sqr(Qdp) + 120*
      MassG*Sqr(g3)*Sqr(Qdp) + 15*traceAdjKappaTKappa*Sqr(QDxbarp) - 15*MassBp*
      traceKappaAdjKappa*Sqr(QDxbarp) + 6*MassB*Sqr(g1)*Sqr(QDxbarp) + 6*MassBp
      *Sqr(g1)*Sqr(QDxbarp) + 120*MassBp*Sqr(g3)*Sqr(QDxbarp) + 120*MassG*Sqr(
      g3)*Sqr(QDxbarp) + 15*traceAdjKappaTKappa*Sqr(QDxp) - 15*MassBp*
      traceKappaAdjKappa*Sqr(QDxp) + 6*MassB*Sqr(g1)*Sqr(QDxp) + 6*MassBp*Sqr(
      g1)*Sqr(QDxp) + 120*MassBp*Sqr(g3)*Sqr(QDxp) + 120*MassG*Sqr(g3)*Sqr(QDxp
      ) + 10*traceAdjYeTYe*Sqr(Qep) - 10*MassBp*traceYeAdjYe*Sqr(Qep) + 18*
      MassB*Sqr(g1)*Sqr(Qep) + 18*MassBp*Sqr(g1)*Sqr(Qep) + 10*
      traceAdjLambda12TLambda12*Sqr(QH1p) + 30*traceAdjYdTYd*Sqr(QH1p) + 10*
      traceAdjYeTYe*Sqr(QH1p) - 10*MassBp*traceLambda12AdjLambda12*Sqr(QH1p) -
      10*MassBp*traceYeAdjYe*Sqr(QH1p) + 9*MassB*Sqr(g1)*Sqr(QH1p) + 9*MassBp*
      Sqr(g1)*Sqr(QH1p) + 45*MassBp*Sqr(g2)*Sqr(QH1p) + 45*MassWB*Sqr(g2)*Sqr(
      QH1p) + 10*traceAdjLambda12TLambda12*Sqr(QH2p) + 30*traceAdjYuTYu*Sqr(
      QH2p) - 10*MassBp*traceLambda12AdjLambda12*Sqr(QH2p) - 30*MassBp*
      traceYuAdjYu*Sqr(QH2p) + 9*MassB*Sqr(g1)*Sqr(QH2p) + 9*MassBp*Sqr(g1)*Sqr
      (QH2p) + 45*MassBp*Sqr(g2)*Sqr(QH2p) + 45*MassWB*Sqr(g2)*Sqr(QH2p) + 3*
      MassB*Sqr(g1)*Sqr(QHpbarp) + 3*MassBp*Sqr(g1)*Sqr(QHpbarp) + 15*MassBp*
      Sqr(g2)*Sqr(QHpbarp) + 15*MassWB*Sqr(g2)*Sqr(QHpbarp) + 3*MassB*Sqr(g1)*
      Sqr(QHpp) + 3*MassBp*Sqr(g1)*Sqr(QHpp) + 15*MassBp*Sqr(g2)*Sqr(QHpp) + 15
      *MassWB*Sqr(g2)*Sqr(QHpp) + 10*traceAdjYeTYe*Sqr(QLp) - 10*MassBp*
      traceYeAdjYe*Sqr(QLp) + 9*MassB*Sqr(g1)*Sqr(QLp) + 9*MassBp*Sqr(g1)*Sqr(
      QLp) + 45*MassBp*Sqr(g2)*Sqr(QLp) + 45*MassWB*Sqr(g2)*Sqr(QLp) + 30*
      traceAdjYdTYd*Sqr(QQp) + 30*traceAdjYuTYu*Sqr(QQp) - 30*MassBp*
      traceYuAdjYu*Sqr(QQp) + 3*MassB*Sqr(g1)*Sqr(QQp) + 3*MassBp*Sqr(g1)*Sqr(
      QQp) + 135*MassBp*Sqr(g2)*Sqr(QQp) + 135*MassWB*Sqr(g2)*Sqr(QQp) + 240*
      MassBp*Sqr(g3)*Sqr(QQp) + 240*MassG*Sqr(g3)*Sqr(QQp) - 30*MassBp*
      traceYdAdjYd*(Sqr(Qdp) + Sqr(QH1p) + Sqr(QQp)) + 15*traceAdjKappaTKappa*
      Sqr(QSp) + 10*traceAdjLambda12TLambda12*Sqr(QSp) - 15*MassBp*
      traceKappaAdjKappa*Sqr(QSp) - 10*MassBp*traceLambda12AdjLambda12*Sqr(QSp)
      + 30*traceAdjYuTYu*Sqr(Qup) - 30*MassBp*traceYuAdjYu*Sqr(Qup) + 24*MassB
      *Sqr(g1)*Sqr(Qup) + 24*MassBp*Sqr(g1)*Sqr(Qup) + 120*MassBp*Sqr(g3)*Sqr(
      Qup) + 120*MassG*Sqr(g3)*Sqr(Qup) - 10*Conj(Lambdax)*(Sqr(QH1p) + Sqr(
      QH2p) + Sqr(QSp))*(MassBp*Lambdax - TLambdax));


   return beta_MassBp;
}

} // namespace flexiblesusy
