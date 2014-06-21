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

// File generated at Sun 15 Jun 2014 19:16:41

#include "genericE6SSM_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of ms2.
 *
 * @return one-loop beta function
 */
double genericE6SSM_soft_parameters::calc_beta_ms2_one_loop(const Soft_traces& soft_traces) const
{
   const auto QSp = INPUT(QSp);
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;
   const double traceconjTKappaTpTKappa =
      TRACE_STRUCT.traceconjTKappaTpTKappa;
   const double traceconjTLambda12TpTLambda12 =
      TRACE_STRUCT.traceconjTLambda12TpTLambda12;
   const double tracemH1I2AdjLambda12Lambda12 =
      TRACE_STRUCT.tracemH1I2AdjLambda12Lambda12;
   const double traceKappaAdjKappaconjmDx2 =
      TRACE_STRUCT.traceKappaAdjKappaconjmDx2;
   const double traceKappaconjmDxbar2AdjKappa =
      TRACE_STRUCT.traceKappaconjmDxbar2AdjKappa;
   const double traceLambda12AdjLambda12conjmH2I2 =
      TRACE_STRUCT.traceLambda12AdjLambda12conjmH2I2;
   const double Tr14 = TRACE_STRUCT.Tr14;


   double beta_ms2;

   beta_ms2 = oneOver16PiSqr*(2*gN*QSp*Tr14 + 6*traceconjTKappaTpTKappa +
      4*traceconjTLambda12TpTLambda12 + 6*ms2*traceKappaAdjKappa + 6*
      traceKappaAdjKappaconjmDx2 + 6*traceKappaconjmDxbar2AdjKappa + 4*ms2*
      traceLambda12AdjLambda12 + 4*traceLambda12AdjLambda12conjmH2I2 + 4*
      tracemH1I2AdjLambda12Lambda12 + 4*(mHd2 + mHu2 + ms2)*AbsSqr(Lambdax) + 4
      *AbsSqr(TLambdax) - 8*AbsSqr(MassBp)*Sqr(gN)*Sqr(QSp));


   return beta_ms2;
}

/**
 * Calculates the two-loop beta function of ms2.
 *
 * @return two-loop beta function
 */
double genericE6SSM_soft_parameters::calc_beta_ms2_two_loop(const Soft_traces& soft_traces) const
{
   const auto QH1p = INPUT(QH1p);
   const auto QH2p = INPUT(QH2p);
   const auto QSp = INPUT(QSp);
   const auto QDxbarp = INPUT(QDxbarp);
   const auto QDxp = INPUT(QDxp);
   const auto Qdp = INPUT(Qdp);
   const auto Qep = INPUT(Qep);
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
   const double traceconjTKappaTpKappa =
      TRACE_STRUCT.traceconjTKappaTpKappa;
   const double traceconjTKappaTpTKappa =
      TRACE_STRUCT.traceconjTKappaTpTKappa;
   const double traceconjTLambda12TpLambda12 =
      TRACE_STRUCT.traceconjTLambda12TpLambda12;
   const double traceconjTLambda12TpTLambda12 =
      TRACE_STRUCT.traceconjTLambda12TpTLambda12;
   const double tracemH1I2AdjLambda12Lambda12 =
      TRACE_STRUCT.tracemH1I2AdjLambda12Lambda12;
   const double traceconjTYdTpYd = TRACE_STRUCT.traceconjTYdTpYd;
   const double traceconjTYdTpTYd = TRACE_STRUCT.traceconjTYdTpTYd;
   const double traceconjTYeTpYe = TRACE_STRUCT.traceconjTYeTpYe;
   const double traceconjTYeTpTYe = TRACE_STRUCT.traceconjTYeTpTYe;
   const double traceconjTYuTpYu = TRACE_STRUCT.traceconjTYuTpYu;
   const double traceconjTYuTpTYu = TRACE_STRUCT.traceconjTYuTpTYu;
   const double tracemd2YdAdjYd = TRACE_STRUCT.tracemd2YdAdjYd;
   const double traceme2YeAdjYe = TRACE_STRUCT.traceme2YeAdjYe;
   const double traceml2AdjYeYe = TRACE_STRUCT.traceml2AdjYeYe;
   const double tracemq2AdjYdYd = TRACE_STRUCT.tracemq2AdjYdYd;
   const double tracemq2AdjYuYu = TRACE_STRUCT.tracemq2AdjYuYu;
   const double tracemu2YuAdjYu = TRACE_STRUCT.tracemu2YuAdjYu;
   const double traceKappaAdjKappaconjmDx2 =
      TRACE_STRUCT.traceKappaAdjKappaconjmDx2;
   const double traceKappaconjmDxbar2AdjKappa =
      TRACE_STRUCT.traceKappaconjmDxbar2AdjKappa;
   const double traceLambda12AdjLambda12conjmH2I2 =
      TRACE_STRUCT.traceLambda12AdjLambda12conjmH2I2;
   const double traceKappaAdjKappaKappaAdjKappa =
      TRACE_STRUCT.traceKappaAdjKappaKappaAdjKappa;
   const double traceKappaAdjKappaTKappaAdjTKappa =
      TRACE_STRUCT.traceKappaAdjKappaTKappaAdjTKappa;
   const double traceKappaAdjTKappaTKappaAdjKappa =
      TRACE_STRUCT.traceKappaAdjTKappaTKappaAdjKappa;
   const double traceLambda12AdjLambda12Lambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12Lambda12AdjLambda12;
   const double traceLambda12AdjLambda12TLambda12AdjTLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12TLambda12AdjTLambda12;
   const double traceLambda12AdjTLambda12TLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjTLambda12TLambda12AdjLambda12;
   const double tracemH1I2AdjLambda12Lambda12AdjLambda12Lambda12 =
      TRACE_STRUCT.tracemH1I2AdjLambda12Lambda12AdjLambda12Lambda12;
   const double traceKappaAdjKappaKappaAdjKappaconjmDx2 =
      TRACE_STRUCT.traceKappaAdjKappaKappaAdjKappaconjmDx2;
   const double traceKappaAdjKappaKappaconjmDxbar2AdjKappa =
      TRACE_STRUCT.traceKappaAdjKappaKappaconjmDxbar2AdjKappa;
   const double traceKappaAdjKappaconjmDx2KappaAdjKappa =
      TRACE_STRUCT.traceKappaAdjKappaconjmDx2KappaAdjKappa;
   const double traceKappaconjmDxbar2AdjKappaKappaAdjKappa =
      TRACE_STRUCT.traceKappaconjmDxbar2AdjKappaKappaAdjKappa;
   const double traceLambda12AdjLambda12Lambda12AdjLambda12conjmH2I2 =
      TRACE_STRUCT.traceLambda12AdjLambda12Lambda12AdjLambda12conjmH2I2;
   const double traceLambda12AdjLambda12conjmH2I2Lambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12conjmH2I2Lambda12AdjLambda12;
   const double Tr2U144 = TRACE_STRUCT.Tr2U144;
   const double Tr34 = TRACE_STRUCT.Tr34;


   double beta_ms2;

   beta_ms2 = -0.8*twoLoop*(-10*gN*QSp*Tr34 + 15*
      traceKappaAdjKappaconjmDx2KappaAdjKappa + 30*ms2*
      traceKappaAdjKappaKappaAdjKappa + 15*
      traceKappaAdjKappaKappaAdjKappaconjmDx2 + 15*
      traceKappaAdjKappaKappaconjmDxbar2AdjKappa + 30*
      traceKappaAdjKappaTKappaAdjTKappa + 30*traceKappaAdjTKappaTKappaAdjKappa
      + 15*traceKappaconjmDxbar2AdjKappaKappaAdjKappa + 10*
      traceLambda12AdjLambda12conjmH2I2Lambda12AdjLambda12 + 20*ms2*
      traceLambda12AdjLambda12Lambda12AdjLambda12 + 10*
      traceLambda12AdjLambda12Lambda12AdjLambda12conjmH2I2 + 20*
      traceLambda12AdjLambda12TLambda12AdjTLambda12 + 20*
      traceLambda12AdjTLambda12TLambda12AdjLambda12 + 20*
      tracemH1I2AdjLambda12Lambda12AdjLambda12Lambda12 + 15*traceYdAdjYd*AbsSqr
      (TLambdax) + 5*traceYeAdjYe*AbsSqr(TLambdax) + 15*traceYuAdjYu*AbsSqr(
      TLambdax) + 15*traceAdjYdTYd*Conj(TLambdax)*Lambdax + 5*traceAdjYeTYe*
      Conj(TLambdax)*Lambdax + 15*traceAdjYuTYu*Conj(TLambdax)*Lambdax + 2*
      MassB*traceconjTKappaTpKappa*Sqr(g1) - 2*traceconjTKappaTpTKappa*Sqr(g1)
      + 3*MassB*traceconjTLambda12TpLambda12*Sqr(g1) - 3*
      traceconjTLambda12TpTLambda12*Sqr(g1) - 2*ms2*traceKappaAdjKappa*Sqr(g1)
      - 2*traceKappaAdjKappaconjmDx2*Sqr(g1) - 2*traceKappaconjmDxbar2AdjKappa*
      Sqr(g1) - 3*ms2*traceLambda12AdjLambda12*Sqr(g1) - 3*
      traceLambda12AdjLambda12conjmH2I2*Sqr(g1) - 3*
      tracemH1I2AdjLambda12Lambda12*Sqr(g1) - 4*traceKappaAdjKappa*AbsSqr(MassB
      )*Sqr(g1) - 6*traceLambda12AdjLambda12*AbsSqr(MassB)*Sqr(g1) - 3*AbsSqr(
      TLambdax)*Sqr(g1) + 2*traceAdjKappaTKappa*Conj(MassB)*Sqr(g1) + 3*
      traceAdjLambda12TLambda12*Conj(MassB)*Sqr(g1) + 3*MassB*Conj(TLambdax)*
      Lambdax*Sqr(g1) + 15*MassWB*traceconjTLambda12TpLambda12*Sqr(g2) - 15*
      traceconjTLambda12TpTLambda12*Sqr(g2) - 15*ms2*traceLambda12AdjLambda12*
      Sqr(g2) - 15*traceLambda12AdjLambda12conjmH2I2*Sqr(g2) - 15*
      tracemH1I2AdjLambda12Lambda12*Sqr(g2) - 30*traceLambda12AdjLambda12*
      AbsSqr(MassWB)*Sqr(g2) - 15*AbsSqr(TLambdax)*Sqr(g2) + 15*
      traceAdjLambda12TLambda12*Conj(MassWB)*Sqr(g2) + 15*MassWB*Conj(TLambdax)
      *Lambdax*Sqr(g2) + 40*MassG*traceconjTKappaTpKappa*Sqr(g3) - 40*
      traceconjTKappaTpTKappa*Sqr(g3) - 40*ms2*traceKappaAdjKappa*Sqr(g3) - 40*
      traceKappaAdjKappaconjmDx2*Sqr(g3) - 40*traceKappaconjmDxbar2AdjKappa*Sqr
      (g3) - 80*traceKappaAdjKappa*AbsSqr(MassG)*Sqr(g3) + 40*
      traceAdjKappaTKappa*Conj(MassG)*Sqr(g3) + 15*MassBp*
      traceconjTKappaTpKappa*Sqr(gN)*Sqr(QDxbarp) - 15*traceconjTKappaTpTKappa*
      Sqr(gN)*Sqr(QDxbarp) - 15*ms2*traceKappaAdjKappa*Sqr(gN)*Sqr(QDxbarp) -
      15*traceKappaAdjKappaconjmDx2*Sqr(gN)*Sqr(QDxbarp) - 15*
      traceKappaconjmDxbar2AdjKappa*Sqr(gN)*Sqr(QDxbarp) + 15*MassBp*
      traceconjTKappaTpKappa*Sqr(gN)*Sqr(QDxp) - 15*traceconjTKappaTpTKappa*Sqr
      (gN)*Sqr(QDxp) - 15*ms2*traceKappaAdjKappa*Sqr(gN)*Sqr(QDxp) - 15*
      traceKappaAdjKappaconjmDx2*Sqr(gN)*Sqr(QDxp) - 15*
      traceKappaconjmDxbar2AdjKappa*Sqr(gN)*Sqr(QDxp) + 10*MassBp*
      traceconjTLambda12TpLambda12*Sqr(gN)*Sqr(QH1p) - 10*
      traceconjTLambda12TpTLambda12*Sqr(gN)*Sqr(QH1p) - 10*ms2*
      traceLambda12AdjLambda12*Sqr(gN)*Sqr(QH1p) - 10*
      traceLambda12AdjLambda12conjmH2I2*Sqr(gN)*Sqr(QH1p) - 10*
      tracemH1I2AdjLambda12Lambda12*Sqr(gN)*Sqr(QH1p) - 10*AbsSqr(TLambdax)*Sqr
      (gN)*Sqr(QH1p) + 10*MassBp*Conj(TLambdax)*Lambdax*Sqr(gN)*Sqr(QH1p) + 10*
      MassBp*traceconjTLambda12TpLambda12*Sqr(gN)*Sqr(QH2p) - 10*
      traceconjTLambda12TpTLambda12*Sqr(gN)*Sqr(QH2p) - 10*ms2*
      traceLambda12AdjLambda12*Sqr(gN)*Sqr(QH2p) - 10*
      traceLambda12AdjLambda12conjmH2I2*Sqr(gN)*Sqr(QH2p) - 10*
      tracemH1I2AdjLambda12Lambda12*Sqr(gN)*Sqr(QH2p) - 10*AbsSqr(TLambdax)*Sqr
      (gN)*Sqr(QH2p) + 10*MassBp*Conj(TLambdax)*Lambdax*Sqr(gN)*Sqr(QH2p) - 10*
      Tr2U144*Sqr(gN)*Sqr(QSp) - 15*MassBp*traceconjTKappaTpKappa*Sqr(gN)*Sqr(
      QSp) + 15*traceconjTKappaTpTKappa*Sqr(gN)*Sqr(QSp) - 10*MassBp*
      traceconjTLambda12TpLambda12*Sqr(gN)*Sqr(QSp) + 10*
      traceconjTLambda12TpTLambda12*Sqr(gN)*Sqr(QSp) + 15*ms2*
      traceKappaAdjKappa*Sqr(gN)*Sqr(QSp) + 15*traceKappaAdjKappaconjmDx2*Sqr(
      gN)*Sqr(QSp) + 15*traceKappaconjmDxbar2AdjKappa*Sqr(gN)*Sqr(QSp) + 10*ms2
      *traceLambda12AdjLambda12*Sqr(gN)*Sqr(QSp) + 10*
      traceLambda12AdjLambda12conjmH2I2*Sqr(gN)*Sqr(QSp) + 10*
      tracemH1I2AdjLambda12Lambda12*Sqr(gN)*Sqr(QSp) + 10*AbsSqr(TLambdax)*Sqr(
      gN)*Sqr(QSp) - 10*MassBp*Conj(TLambdax)*Lambdax*Sqr(gN)*Sqr(QSp) + 20*(
      mHd2 + mHu2 + ms2)*Sqr(Conj(Lambdax))*Sqr(Lambdax) - 5*Conj(MassBp)*Sqr(
      gN)*(30*MassBp*Power(QSp,4)*Sqr(gN) - 3*traceAdjKappaTKappa*Sqr(QDxbarp)
      - 3*traceAdjKappaTKappa*Sqr(QDxp) - 2*traceAdjLambda12TLambda12*Sqr(QH1p)
      + 4*MassBp*traceLambda12AdjLambda12*Sqr(QH1p) - 2*
      traceAdjLambda12TLambda12*Sqr(QH2p) + 4*MassBp*traceLambda12AdjLambda12*
      Sqr(QH2p) + 6*MassBp*traceKappaAdjKappa*(Sqr(QDxbarp) + Sqr(QDxp) - Sqr(
      QSp)) + 3*traceAdjKappaTKappa*Sqr(QSp) + 2*traceAdjLambda12TLambda12*Sqr(
      QSp) - 4*MassBp*traceLambda12AdjLambda12*Sqr(QSp) + 54*MassBp*Sqr(gN)*Sqr
      (Qdp)*Sqr(QSp) + 54*MassBp*Sqr(gN)*Sqr(QDxbarp)*Sqr(QSp) + 54*MassBp*Sqr(
      gN)*Sqr(QDxp)*Sqr(QSp) + 18*MassBp*Sqr(gN)*Sqr(Qep)*Sqr(QSp) + 36*MassBp*
      Sqr(gN)*Sqr(QH1p)*Sqr(QSp) + 36*MassBp*Sqr(gN)*Sqr(QH2p)*Sqr(QSp) + 12*
      MassBp*Sqr(gN)*Sqr(QHpbarp)*Sqr(QSp) + 12*MassBp*Sqr(gN)*Sqr(QHpp)*Sqr(
      QSp) + 36*MassBp*Sqr(gN)*Sqr(QLp)*Sqr(QSp) + 108*MassBp*Sqr(gN)*Sqr(QQp)*
      Sqr(QSp) + 54*MassBp*Sqr(gN)*Sqr(QSp)*Sqr(Qup) + 2*Conj(Lambdax)*(Sqr(
      QH1p) + Sqr(QH2p) - Sqr(QSp))*(2*MassBp*Lambdax - TLambdax)) + Conj(
      Lambdax)*(15*traceconjTYdTpTYd*Lambdax + 5*traceconjTYeTpTYe*Lambdax + 15
      *traceconjTYuTpTYu*Lambdax + 15*tracemd2YdAdjYd*Lambdax + 5*
      traceme2YeAdjYe*Lambdax + 5*traceml2AdjYeYe*Lambdax + 15*tracemq2AdjYdYd*
      Lambdax + 15*tracemq2AdjYuYu*Lambdax + 15*tracemu2YuAdjYu*Lambdax + 30*
      mHd2*traceYdAdjYd*Lambdax + 15*mHu2*traceYdAdjYd*Lambdax + 15*ms2*
      traceYdAdjYd*Lambdax + 10*mHd2*traceYeAdjYe*Lambdax + 5*mHu2*traceYeAdjYe
      *Lambdax + 5*ms2*traceYeAdjYe*Lambdax + 15*mHd2*traceYuAdjYu*Lambdax + 30
      *mHu2*traceYuAdjYu*Lambdax + 15*ms2*traceYuAdjYu*Lambdax + 40*AbsSqr(
      TLambdax)*Lambdax - 3*mHd2*Lambdax*Sqr(g1) - 3*mHu2*Lambdax*Sqr(g1) - 3*
      ms2*Lambdax*Sqr(g1) - 15*mHd2*Lambdax*Sqr(g2) - 15*mHu2*Lambdax*Sqr(g2) -
      15*ms2*Lambdax*Sqr(g2) - 10*mHd2*Lambdax*Sqr(gN)*Sqr(QH1p) - 10*mHu2*
      Lambdax*Sqr(gN)*Sqr(QH1p) - 10*ms2*Lambdax*Sqr(gN)*Sqr(QH1p) - 10*mHd2*
      Lambdax*Sqr(gN)*Sqr(QH2p) - 10*mHu2*Lambdax*Sqr(gN)*Sqr(QH2p) - 10*ms2*
      Lambdax*Sqr(gN)*Sqr(QH2p) + 10*mHd2*Lambdax*Sqr(gN)*Sqr(QSp) + 10*mHu2*
      Lambdax*Sqr(gN)*Sqr(QSp) + 10*ms2*Lambdax*Sqr(gN)*Sqr(QSp) + 15*
      traceconjTYdTpYd*TLambdax + 5*traceconjTYeTpYe*TLambdax + 15*
      traceconjTYuTpYu*TLambdax + 3*Conj(MassB)*Sqr(g1)*(-2*MassB*Lambdax +
      TLambdax) + 15*Conj(MassWB)*Sqr(g2)*(-2*MassWB*Lambdax + TLambdax)));


   return beta_ms2;
}

} // namespace flexiblesusy
