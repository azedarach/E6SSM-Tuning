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

// File generated at Sun 24 Aug 2014 16:10:12

#include "lowE6SSM_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of mHd2.
 *
 * @return one-loop beta function
 */
double lowE6SSM_soft_parameters::calc_beta_mHd2_one_loop(const Soft_traces& soft_traces) const
{
   const auto QH1p = INPUT(QH1p);
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceconjTYdTpTYd = TRACE_STRUCT.traceconjTYdTpTYd;
   const double traceconjTYeTpTYe = TRACE_STRUCT.traceconjTYeTpTYe;
   const double tracemd2YdAdjYd = TRACE_STRUCT.tracemd2YdAdjYd;
   const double traceme2YeAdjYe = TRACE_STRUCT.traceme2YeAdjYe;
   const double traceml2AdjYeYe = TRACE_STRUCT.traceml2AdjYeYe;
   const double tracemq2AdjYdYd = TRACE_STRUCT.tracemq2AdjYdYd;
   const double Tr11 = TRACE_STRUCT.Tr11;
   const double Tr14 = TRACE_STRUCT.Tr14;


   double beta_mHd2;

   beta_mHd2 = oneOver16PiSqr*(-0.7745966692414834*g1*Tr11 + 2*gN*QH1p*
      Tr14 + 6*traceconjTYdTpTYd + 2*traceconjTYeTpTYe + 6*tracemd2YdAdjYd + 2*
      traceme2YeAdjYe + 2*traceml2AdjYeYe + 6*tracemq2AdjYdYd + 6*mHd2*
      traceYdAdjYd + 2*mHd2*traceYeAdjYe + 2*mHd2*AbsSqr(Lambdax) + 2*mHu2*
      AbsSqr(Lambdax) + 2*ms2*AbsSqr(Lambdax) + 2*AbsSqr(TLambdax) - 1.2*AbsSqr
      (MassB)*Sqr(g1) - 6*AbsSqr(MassWB)*Sqr(g2) - 8*AbsSqr(MassBp)*Sqr(gN)*Sqr
      (QH1p));


   return beta_mHd2;
}

/**
 * Calculates the two-loop beta function of mHd2.
 *
 * @return two-loop beta function
 */
double lowE6SSM_soft_parameters::calc_beta_mHd2_two_loop(const Soft_traces& soft_traces) const
{
   const auto QH1p = INPUT(QH1p);
   const auto QH2p = INPUT(QH2p);
   const auto QSp = INPUT(QSp);
   const auto Qdp = INPUT(Qdp);
   const auto QQp = INPUT(QQp);
   const auto Qep = INPUT(Qep);
   const auto QLp = INPUT(QLp);
   const auto QDxbarp = INPUT(QDxbarp);
   const auto QDxp = INPUT(QDxp);
   const auto QHpbarp = INPUT(QHpbarp);
   const auto QHpp = INPUT(QHpp);
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
   const double traceconjTYdTpYd = TRACE_STRUCT.traceconjTYdTpYd;
   const double traceconjTYdTpTYd = TRACE_STRUCT.traceconjTYdTpTYd;
   const double traceconjTYeTpYe = TRACE_STRUCT.traceconjTYeTpYe;
   const double traceconjTYeTpTYe = TRACE_STRUCT.traceconjTYeTpTYe;
   const double traceconjTYuTpYu = TRACE_STRUCT.traceconjTYuTpYu;
   const double traceconjTYuTpTYu = TRACE_STRUCT.traceconjTYuTpTYu;
   const double traceconjTKappaTpKappa =
      TRACE_STRUCT.traceconjTKappaTpKappa;
   const double traceconjTKappaTpTKappa =
      TRACE_STRUCT.traceconjTKappaTpTKappa;
   const double traceconjTLambda12TpLambda12 =
      TRACE_STRUCT.traceconjTLambda12TpLambda12;
   const double traceconjTLambda12TpTLambda12 =
      TRACE_STRUCT.traceconjTLambda12TpTLambda12;
   const double tracemd2YdAdjYd = TRACE_STRUCT.tracemd2YdAdjYd;
   const double traceme2YeAdjYe = TRACE_STRUCT.traceme2YeAdjYe;
   const double tracemH1I2AdjLambda12Lambda12 =
      TRACE_STRUCT.tracemH1I2AdjLambda12Lambda12;
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
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYdTYdAdjTYd =
      TRACE_STRUCT.traceYdAdjYdTYdAdjTYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYdAdjYuTYuAdjTYd =
      TRACE_STRUCT.traceYdAdjYuTYuAdjTYd;
   const double traceYdAdjTYdTYdAdjYd =
      TRACE_STRUCT.traceYdAdjTYdTYdAdjYd;
   const double traceYdAdjTYuTYuAdjYd =
      TRACE_STRUCT.traceYdAdjTYuTYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYeAdjYeTYeAdjTYe =
      TRACE_STRUCT.traceYeAdjYeTYeAdjTYe;
   const double traceYeAdjTYeTYeAdjYe =
      TRACE_STRUCT.traceYeAdjTYeTYeAdjYe;
   const double traceYuAdjYdTYdAdjTYu =
      TRACE_STRUCT.traceYuAdjYdTYdAdjTYu;
   const double traceYuAdjTYdTYdAdjYu =
      TRACE_STRUCT.traceYuAdjTYdTYdAdjYu;
   const double tracemd2YdAdjYdYdAdjYd =
      TRACE_STRUCT.tracemd2YdAdjYdYdAdjYd;
   const double tracemd2YdAdjYuYuAdjYd =
      TRACE_STRUCT.tracemd2YdAdjYuYuAdjYd;
   const double traceme2YeAdjYeYeAdjYe =
      TRACE_STRUCT.traceme2YeAdjYeYeAdjYe;
   const double traceml2AdjYeYeAdjYeYe =
      TRACE_STRUCT.traceml2AdjYeYeAdjYeYe;
   const double tracemq2AdjYdYdAdjYdYd =
      TRACE_STRUCT.tracemq2AdjYdYdAdjYdYd;
   const double tracemq2AdjYdYdAdjYuYu =
      TRACE_STRUCT.tracemq2AdjYdYdAdjYuYu;
   const double tracemq2AdjYuYuAdjYdYd =
      TRACE_STRUCT.tracemq2AdjYuYuAdjYdYd;
   const double tracemu2YuAdjYdYdAdjYu =
      TRACE_STRUCT.tracemu2YuAdjYdYdAdjYu;
   const double Tr2U111 = TRACE_STRUCT.Tr2U111;
   const double Tr2U114 = TRACE_STRUCT.Tr2U114;
   const double Tr31 = TRACE_STRUCT.Tr31;
   const double Tr22 = TRACE_STRUCT.Tr22;
   const double Tr2U141 = TRACE_STRUCT.Tr2U141;
   const double Tr2U144 = TRACE_STRUCT.Tr2U144;
   const double Tr34 = TRACE_STRUCT.Tr34;


   double beta_mHd2;

   beta_mHd2 = twoLoop*(6*Power(g2,4)*Tr22 - 3.0983866769659336*g1*gN*
      QH1p*Tr2U114 - 3.0983866769659336*g1*gN*QH1p*Tr2U141 - 3.0983866769659336
      *g1*Tr31 + 8*gN*QH1p*Tr34 - 36*tracemd2YdAdjYdYdAdjYd - 6*
      tracemd2YdAdjYuYuAdjYd - 12*traceme2YeAdjYeYeAdjYe - 12*
      traceml2AdjYeYeAdjYeYe - 36*tracemq2AdjYdYdAdjYdYd - 6*
      tracemq2AdjYdYdAdjYuYu - 6*tracemq2AdjYuYuAdjYdYd - 6*
      tracemu2YuAdjYdYdAdjYu - 36*traceYdAdjTYdTYdAdjYd - 6*
      traceYdAdjTYuTYuAdjYd - 36*traceYdAdjYdTYdAdjTYd - 36*mHd2*
      traceYdAdjYdYdAdjYd - 6*traceYdAdjYuTYuAdjTYd - 6*mHd2*
      traceYdAdjYuYuAdjYd - 6*mHu2*traceYdAdjYuYuAdjYd - 12*
      traceYeAdjTYeTYeAdjYe - 12*traceYeAdjYeTYeAdjTYe - 12*mHd2*
      traceYeAdjYeYeAdjYe - 6*traceYuAdjTYdTYdAdjYu - 6*traceYuAdjYdTYdAdjTYu +
      87*Power(g2,4)*AbsSqr(MassWB) - 6*traceconjTKappaTpTKappa*AbsSqr(Lambdax
      ) - 4*traceconjTLambda12TpTLambda12*AbsSqr(Lambdax) - 6*traceconjTYuTpTYu
      *AbsSqr(Lambdax) - 6*mHd2*traceKappaAdjKappa*AbsSqr(Lambdax) - 6*mHu2*
      traceKappaAdjKappa*AbsSqr(Lambdax) - 12*ms2*traceKappaAdjKappa*AbsSqr(
      Lambdax) - 6*traceKappaAdjKappaconjmDx2*AbsSqr(Lambdax) - 6*
      traceKappaconjmDxbar2AdjKappa*AbsSqr(Lambdax) - 4*mHd2*
      traceLambda12AdjLambda12*AbsSqr(Lambdax) - 4*mHu2*
      traceLambda12AdjLambda12*AbsSqr(Lambdax) - 8*ms2*traceLambda12AdjLambda12
      *AbsSqr(Lambdax) - 4*traceLambda12AdjLambda12conjmH2I2*AbsSqr(Lambdax) -
      4*tracemH1I2AdjLambda12Lambda12*AbsSqr(Lambdax) - 6*tracemq2AdjYuYu*
      AbsSqr(Lambdax) - 6*tracemu2YuAdjYu*AbsSqr(Lambdax) - 6*mHd2*traceYuAdjYu
      *AbsSqr(Lambdax) - 12*mHu2*traceYuAdjYu*AbsSqr(Lambdax) - 6*ms2*
      traceYuAdjYu*AbsSqr(Lambdax) - 6*traceKappaAdjKappa*AbsSqr(TLambdax) - 4*
      traceLambda12AdjLambda12*AbsSqr(TLambdax) - 6*traceYuAdjYu*AbsSqr(
      TLambdax) - 24*AbsSqr(Lambdax)*AbsSqr(TLambdax) - 6*traceAdjKappaTKappa*
      Conj(TLambdax)*Lambdax - 4*traceAdjLambda12TLambda12*Conj(TLambdax)*
      Lambdax - 6*traceAdjYuTYu*Conj(TLambdax)*Lambdax + 1.2*Tr2U111*Sqr(g1) -
      0.8*traceconjTYdTpTYd*Sqr(g1) + 0.8*MassB*traceconjTYdTpYd*Sqr(g1) + 2.4*
      traceconjTYeTpTYe*Sqr(g1) - 2.4*MassB*traceconjTYeTpYe*Sqr(g1) - 0.8*
      tracemd2YdAdjYd*Sqr(g1) + 2.4*traceme2YeAdjYe*Sqr(g1) + 2.4*
      traceml2AdjYeYe*Sqr(g1) - 0.8*tracemq2AdjYdYd*Sqr(g1) - 0.8*mHd2*
      traceYdAdjYd*Sqr(g1) + 2.4*mHd2*traceYeAdjYe*Sqr(g1) + 3.6*AbsSqr(MassWB)
      *Sqr(g1)*Sqr(g2) + 1.8*MassB*Conj(MassWB)*Sqr(g1)*Sqr(g2) + 32*
      traceconjTYdTpTYd*Sqr(g3) - 32*MassG*traceconjTYdTpYd*Sqr(g3) + 32*
      tracemd2YdAdjYd*Sqr(g3) + 32*tracemq2AdjYdYd*Sqr(g3) + 32*mHd2*
      traceYdAdjYd*Sqr(g3) + 64*traceYdAdjYd*AbsSqr(MassG)*Sqr(g3) - 32*
      traceAdjYdTYd*Conj(MassG)*Sqr(g3) + 12*traceconjTYdTpTYd*Sqr(gN)*Sqr(Qdp)
      - 12*MassBp*traceconjTYdTpYd*Sqr(gN)*Sqr(Qdp) + 12*tracemd2YdAdjYd*Sqr(
      gN)*Sqr(Qdp) + 12*tracemq2AdjYdYd*Sqr(gN)*Sqr(Qdp) + 12*mHd2*traceYdAdjYd
      *Sqr(gN)*Sqr(Qdp) + 4*traceconjTYeTpTYe*Sqr(gN)*Sqr(Qep) - 4*MassBp*
      traceconjTYeTpYe*Sqr(gN)*Sqr(Qep) + 4*traceme2YeAdjYe*Sqr(gN)*Sqr(Qep) +
      4*traceml2AdjYeYe*Sqr(gN)*Sqr(Qep) + 4*mHd2*traceYeAdjYe*Sqr(gN)*Sqr(Qep)
      + 8*Tr2U144*Sqr(gN)*Sqr(QH1p) - 12*traceconjTYdTpTYd*Sqr(gN)*Sqr(QH1p) +
      12*MassBp*traceconjTYdTpYd*Sqr(gN)*Sqr(QH1p) - 4*traceconjTYeTpTYe*Sqr(
      gN)*Sqr(QH1p) + 4*MassBp*traceconjTYeTpYe*Sqr(gN)*Sqr(QH1p) - 12*
      tracemd2YdAdjYd*Sqr(gN)*Sqr(QH1p) - 4*traceme2YeAdjYe*Sqr(gN)*Sqr(QH1p) -
      4*traceml2AdjYeYe*Sqr(gN)*Sqr(QH1p) - 12*tracemq2AdjYdYd*Sqr(gN)*Sqr(
      QH1p) - 12*mHd2*traceYdAdjYd*Sqr(gN)*Sqr(QH1p) - 4*mHd2*traceYeAdjYe*Sqr(
      gN)*Sqr(QH1p) - 4*mHd2*AbsSqr(Lambdax)*Sqr(gN)*Sqr(QH1p) - 4*mHu2*AbsSqr(
      Lambdax)*Sqr(gN)*Sqr(QH1p) - 4*ms2*AbsSqr(Lambdax)*Sqr(gN)*Sqr(QH1p) - 4*
      AbsSqr(TLambdax)*Sqr(gN)*Sqr(QH1p) + 4*MassBp*Conj(TLambdax)*Lambdax*Sqr(
      gN)*Sqr(QH1p) + 24*AbsSqr(MassWB)*Sqr(g2)*Sqr(gN)*Sqr(QH1p) + 12*MassBp*
      Conj(MassWB)*Sqr(g2)*Sqr(gN)*Sqr(QH1p) + 0.04*Conj(MassB)*Sqr(g1)*(20*
      traceAdjYdTYd - 60*traceAdjYeTYe - 40*MassB*traceYdAdjYd + 120*MassB*
      traceYeAdjYe + 891*MassB*Sqr(g1) + 90*MassB*Sqr(g2) + 45*MassWB*Sqr(g2) -
      360*MassB*Qdp*QH1p*Sqr(gN) - 180*MassBp*Qdp*QH1p*Sqr(gN) - 360*MassB*
      QDxbarp*QH1p*Sqr(gN) - 180*MassBp*QDxbarp*QH1p*Sqr(gN) + 360*MassB*QDxp*
      QH1p*Sqr(gN) + 180*MassBp*QDxp*QH1p*Sqr(gN) - 360*MassB*Qep*QH1p*Sqr(gN)
      - 180*MassBp*Qep*QH1p*Sqr(gN) - 360*MassB*QH1p*QH2p*Sqr(gN) - 180*MassBp*
      QH1p*QH2p*Sqr(gN) - 120*MassB*QH1p*QHpbarp*Sqr(gN) - 60*MassBp*QH1p*
      QHpbarp*Sqr(gN) + 120*MassB*QH1p*QHpp*Sqr(gN) + 60*MassBp*QH1p*QHpp*Sqr(
      gN) + 360*MassB*QH1p*QLp*Sqr(gN) + 180*MassBp*QH1p*QLp*Sqr(gN) - 360*
      MassB*QH1p*QQp*Sqr(gN) - 180*MassBp*QH1p*QQp*Sqr(gN) + 720*MassB*QH1p*Qup
      *Sqr(gN) + 360*MassBp*QH1p*Qup*Sqr(gN) + 480*MassB*Sqr(gN)*Sqr(QH1p) +
      240*MassBp*Sqr(gN)*Sqr(QH1p)) + 4*mHd2*AbsSqr(Lambdax)*Sqr(gN)*Sqr(QH2p)
      + 4*mHu2*AbsSqr(Lambdax)*Sqr(gN)*Sqr(QH2p) + 4*ms2*AbsSqr(Lambdax)*Sqr(gN
      )*Sqr(QH2p) + 4*AbsSqr(TLambdax)*Sqr(gN)*Sqr(QH2p) - 4*MassBp*Conj(
      TLambdax)*Lambdax*Sqr(gN)*Sqr(QH2p) + 4*traceconjTYeTpTYe*Sqr(gN)*Sqr(QLp
      ) - 4*MassBp*traceconjTYeTpYe*Sqr(gN)*Sqr(QLp) + 4*traceme2YeAdjYe*Sqr(gN
      )*Sqr(QLp) + 4*traceml2AdjYeYe*Sqr(gN)*Sqr(QLp) + 4*mHd2*traceYeAdjYe*Sqr
      (gN)*Sqr(QLp) + 12*traceconjTYdTpTYd*Sqr(gN)*Sqr(QQp) - 12*MassBp*
      traceconjTYdTpYd*Sqr(gN)*Sqr(QQp) + 12*tracemd2YdAdjYd*Sqr(gN)*Sqr(QQp) +
      12*tracemq2AdjYdYd*Sqr(gN)*Sqr(QQp) + 12*mHd2*traceYdAdjYd*Sqr(gN)*Sqr(
      QQp) + 4*mHd2*AbsSqr(Lambdax)*Sqr(gN)*Sqr(QSp) + 4*mHu2*AbsSqr(Lambdax)*
      Sqr(gN)*Sqr(QSp) + 4*ms2*AbsSqr(Lambdax)*Sqr(gN)*Sqr(QSp) + 4*AbsSqr(
      TLambdax)*Sqr(gN)*Sqr(QSp) - 4*MassBp*Conj(TLambdax)*Lambdax*Sqr(gN)*Sqr(
      QSp) - 12*mHd2*Sqr(Conj(Lambdax))*Sqr(Lambdax) - 12*mHu2*Sqr(Conj(Lambdax
      ))*Sqr(Lambdax) - 12*ms2*Sqr(Conj(Lambdax))*Sqr(Lambdax) + 0.8*Conj(
      MassBp)*Sqr(gN)*(-9*MassB*Qdp*QH1p*Sqr(g1) - 18*MassBp*Qdp*QH1p*Sqr(g1) -
      9*MassB*QDxbarp*QH1p*Sqr(g1) - 18*MassBp*QDxbarp*QH1p*Sqr(g1) + 9*MassB*
      QDxp*QH1p*Sqr(g1) + 18*MassBp*QDxp*QH1p*Sqr(g1) - 9*MassB*Qep*QH1p*Sqr(g1
      ) - 18*MassBp*Qep*QH1p*Sqr(g1) - 9*MassB*QH1p*QH2p*Sqr(g1) - 18*MassBp*
      QH1p*QH2p*Sqr(g1) - 3*MassB*QH1p*QHpbarp*Sqr(g1) - 6*MassBp*QH1p*QHpbarp*
      Sqr(g1) + 3*MassB*QH1p*QHpp*Sqr(g1) + 6*MassBp*QH1p*QHpp*Sqr(g1) + 9*
      MassB*QH1p*QLp*Sqr(g1) + 18*MassBp*QH1p*QLp*Sqr(g1) - 9*MassB*QH1p*QQp*
      Sqr(g1) - 18*MassBp*QH1p*QQp*Sqr(g1) + 18*MassB*QH1p*Qup*Sqr(g1) + 36*
      MassBp*QH1p*Qup*Sqr(g1) + 240*MassBp*Power(QH1p,4)*Sqr(gN) - 15*
      traceAdjYdTYd*Sqr(Qdp) - 5*traceAdjYeTYe*Sqr(Qep) + 10*MassBp*
      traceYeAdjYe*Sqr(Qep) + 15*traceAdjYdTYd*Sqr(QH1p) + 5*traceAdjYeTYe*Sqr(
      QH1p) - 10*MassBp*traceYeAdjYe*Sqr(QH1p) + 12*MassB*Sqr(g1)*Sqr(QH1p) +
      24*MassBp*Sqr(g1)*Sqr(QH1p) + 30*MassBp*Sqr(g2)*Sqr(QH1p) + 15*MassWB*Sqr
      (g2)*Sqr(QH1p) + 270*MassBp*Sqr(gN)*Sqr(Qdp)*Sqr(QH1p) + 270*MassBp*Sqr(
      gN)*Sqr(QDxbarp)*Sqr(QH1p) + 270*MassBp*Sqr(gN)*Sqr(QDxp)*Sqr(QH1p) + 90*
      MassBp*Sqr(gN)*Sqr(Qep)*Sqr(QH1p) + 180*MassBp*Sqr(gN)*Sqr(QH1p)*Sqr(QH2p
      ) + 60*MassBp*Sqr(gN)*Sqr(QH1p)*Sqr(QHpbarp) + 60*MassBp*Sqr(gN)*Sqr(QH1p
      )*Sqr(QHpp) - 5*traceAdjYeTYe*Sqr(QLp) + 10*MassBp*traceYeAdjYe*Sqr(QLp)
      + 180*MassBp*Sqr(gN)*Sqr(QH1p)*Sqr(QLp) - 15*traceAdjYdTYd*Sqr(QQp) + 540
      *MassBp*Sqr(gN)*Sqr(QH1p)*Sqr(QQp) + 30*MassBp*traceYdAdjYd*(Sqr(Qdp) -
      Sqr(QH1p) + Sqr(QQp)) + 90*MassBp*Sqr(gN)*Sqr(QH1p)*Sqr(QSp) + 270*MassBp
      *Sqr(gN)*Sqr(QH1p)*Sqr(Qup) - 5*Conj(Lambdax)*(Sqr(QH1p) - Sqr(QH2p) -
      Sqr(QSp))*(2*MassBp*Lambdax - TLambdax)) - 6*traceconjTKappaTpKappa*Conj(
      Lambdax)*TLambdax - 4*traceconjTLambda12TpLambda12*Conj(Lambdax)*TLambdax
      - 6*traceconjTYuTpYu*Conj(Lambdax)*TLambdax);


   return beta_mHd2;
}

} // namespace flexiblesusy
