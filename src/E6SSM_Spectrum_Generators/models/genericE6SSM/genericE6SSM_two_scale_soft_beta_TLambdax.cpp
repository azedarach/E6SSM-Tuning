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

// File generated at Sun 15 Jun 2014 19:16:33

#include "genericE6SSM_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of TLambdax.
 *
 * @return one-loop beta function
 */
double genericE6SSM_soft_parameters::calc_beta_TLambdax_one_loop(const Soft_traces& soft_traces) const
{
   const auto QH1p = INPUT(QH1p);
   const auto QH2p = INPUT(QH2p);
   const auto QSp = INPUT(QSp);
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


   double beta_TLambdax;

   beta_TLambdax = oneOver16PiSqr*(6*traceAdjKappaTKappa*Lambdax + 4*
      traceAdjLambda12TLambda12*Lambdax + 6*traceAdjYdTYd*Lambdax + 2*
      traceAdjYeTYe*Lambdax + 6*traceAdjYuTYu*Lambdax + 1.2*MassB*Lambdax*Sqr(
      g1) + 6*MassWB*Lambdax*Sqr(g2) + 4*MassBp*Lambdax*Sqr(gN)*Sqr(QH1p) + 4*
      MassBp*Lambdax*Sqr(gN)*Sqr(QH2p) + 4*MassBp*Lambdax*Sqr(gN)*Sqr(QSp) + (3
      *traceKappaAdjKappa + 2*traceLambda12AdjLambda12 + 3*traceYdAdjYd +
      traceYeAdjYe + 3*traceYuAdjYu + 12*AbsSqr(Lambdax) - 0.6*Sqr(g1) - 3*Sqr(
      g2) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QSp))*
      TLambdax);


   return beta_TLambdax;
}

/**
 * Calculates the two-loop beta function of TLambdax.
 *
 * @return two-loop beta function
 */
double genericE6SSM_soft_parameters::calc_beta_TLambdax_two_loop(const Soft_traces& soft_traces) const
{
   const auto QH1p = INPUT(QH1p);
   const auto QH2p = INPUT(QH2p);
   const auto Qdp = INPUT(Qdp);
   const auto QDxbarp = INPUT(QDxbarp);
   const auto QDxp = INPUT(QDxp);
   const auto Qep = INPUT(Qep);
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
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;
   const double traceKappaAdjKappaKappaAdjKappa =
      TRACE_STRUCT.traceKappaAdjKappaKappaAdjKappa;
   const double traceLambda12AdjLambda12Lambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12Lambda12AdjLambda12;
   const double traceYdAdjYdTYdAdjYd = TRACE_STRUCT.traceYdAdjYdTYdAdjYd;
   const double traceYdAdjYuTYuAdjYd = TRACE_STRUCT.traceYdAdjYuTYuAdjYd;
   const double traceYeAdjYeTYeAdjYe = TRACE_STRUCT.traceYeAdjYeTYeAdjYe;
   const double traceYuAdjYdTYdAdjYu = TRACE_STRUCT.traceYuAdjYdTYdAdjYu;
   const double traceYuAdjYuTYuAdjYu = TRACE_STRUCT.traceYuAdjYuTYuAdjYu;
   const double traceKappaAdjKappaTKappaAdjKappa =
      TRACE_STRUCT.traceKappaAdjKappaTKappaAdjKappa;
   const double traceLambda12AdjLambda12TLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12TLambda12AdjLambda12;


   double beta_TLambdax;

   beta_TLambdax = twoLoop*(-0.08*Lambdax*(297*Power(g1,4)*MassB + 825*
      Power(g2,4)*MassWB + 800*Power(gN,4)*MassBp*Power(QH1p,4) + 800*Power(gN,
      4)*MassBp*Power(QH2p,4) + 500*Power(gN,4)*MassBp*Power(QSp,4) + 300*
      traceKappaAdjKappaTKappaAdjKappa + 200*
      traceLambda12AdjLambda12TLambda12AdjLambda12 + 450*traceYdAdjYdTYdAdjYd +
      150*traceYdAdjYuTYuAdjYd + 150*traceYeAdjYeTYeAdjYe + 150*
      traceYuAdjYdTYdAdjYu + 450*traceYuAdjYuTYuAdjYu - 20*traceAdjKappaTKappa*
      Sqr(g1) - 30*traceAdjLambda12TLambda12*Sqr(g1) + 10*traceAdjYdTYd*Sqr(g1)
      - 30*traceAdjYeTYe*Sqr(g1) - 20*traceAdjYuTYu*Sqr(g1) + 20*MassB*
      traceKappaAdjKappa*Sqr(g1) + 30*MassB*traceLambda12AdjLambda12*Sqr(g1) +
      20*MassB*traceYuAdjYu*Sqr(g1) - 150*traceAdjLambda12TLambda12*Sqr(g2) +
      150*MassWB*traceLambda12AdjLambda12*Sqr(g2) + 45*MassB*Sqr(g1)*Sqr(g2) +
      45*MassWB*Sqr(g1)*Sqr(g2) - 400*traceAdjKappaTKappa*Sqr(g3) - 400*
      traceAdjYdTYd*Sqr(g3) - 400*traceAdjYuTYu*Sqr(g3) + 400*MassG*
      traceKappaAdjKappa*Sqr(g3) + 400*MassG*traceYuAdjYu*Sqr(g3) - 150*
      traceAdjYdTYd*Sqr(gN)*Sqr(Qdp) - 150*traceAdjKappaTKappa*Sqr(gN)*Sqr(
      QDxbarp) + 150*MassBp*traceKappaAdjKappa*Sqr(gN)*Sqr(QDxbarp) - 150*
      traceAdjKappaTKappa*Sqr(gN)*Sqr(QDxp) + 150*MassBp*traceKappaAdjKappa*Sqr
      (gN)*Sqr(QDxp) - 50*traceAdjYeTYe*Sqr(gN)*Sqr(Qep) - 100*
      traceAdjLambda12TLambda12*Sqr(gN)*Sqr(QH1p) + 150*traceAdjYdTYd*Sqr(gN)*
      Sqr(QH1p) + 50*traceAdjYeTYe*Sqr(gN)*Sqr(QH1p) + 100*MassBp*
      traceLambda12AdjLambda12*Sqr(gN)*Sqr(QH1p) + 30*MassB*Sqr(g1)*Sqr(gN)*Sqr
      (QH1p) + 30*MassBp*Sqr(g1)*Sqr(gN)*Sqr(QH1p) + 150*MassBp*Sqr(g2)*Sqr(gN)
      *Sqr(QH1p) + 150*MassWB*Sqr(g2)*Sqr(gN)*Sqr(QH1p) + 900*Power(gN,4)*
      MassBp*Sqr(Qdp)*Sqr(QH1p) + 900*Power(gN,4)*MassBp*Sqr(QDxbarp)*Sqr(QH1p)
      + 900*Power(gN,4)*MassBp*Sqr(QDxp)*Sqr(QH1p) + 300*Power(gN,4)*MassBp*
      Sqr(Qep)*Sqr(QH1p) - 100*traceAdjLambda12TLambda12*Sqr(gN)*Sqr(QH2p) +
      150*traceAdjYuTYu*Sqr(gN)*Sqr(QH2p) + 100*MassBp*traceLambda12AdjLambda12
      *Sqr(gN)*Sqr(QH2p) - 150*MassBp*traceYuAdjYu*Sqr(gN)*Sqr(QH2p) + 30*MassB
      *Sqr(g1)*Sqr(gN)*Sqr(QH2p) + 30*MassBp*Sqr(g1)*Sqr(gN)*Sqr(QH2p) + 150*
      MassBp*Sqr(g2)*Sqr(gN)*Sqr(QH2p) + 150*MassWB*Sqr(g2)*Sqr(gN)*Sqr(QH2p) +
      900*Power(gN,4)*MassBp*Sqr(Qdp)*Sqr(QH2p) + 900*Power(gN,4)*MassBp*Sqr(
      QDxbarp)*Sqr(QH2p) + 900*Power(gN,4)*MassBp*Sqr(QDxp)*Sqr(QH2p) + 300*
      Power(gN,4)*MassBp*Sqr(Qep)*Sqr(QH2p) + 1200*Power(gN,4)*MassBp*Sqr(QH1p)
      *Sqr(QH2p) + 200*Power(gN,4)*MassBp*Sqr(QH1p)*Sqr(QHpbarp) + 200*Power(gN
      ,4)*MassBp*Sqr(QH2p)*Sqr(QHpbarp) + 200*Power(gN,4)*MassBp*Sqr(QH1p)*Sqr(
      QHpp) + 200*Power(gN,4)*MassBp*Sqr(QH2p)*Sqr(QHpp) - 50*traceAdjYeTYe*Sqr
      (gN)*Sqr(QLp) + 600*Power(gN,4)*MassBp*Sqr(QH1p)*Sqr(QLp) + 600*Power(gN,
      4)*MassBp*Sqr(QH2p)*Sqr(QLp) + 10*traceYeAdjYe*(3*MassB*Sqr(g1) + 5*
      MassBp*Sqr(gN)*(Sqr(Qep) - Sqr(QH1p) + Sqr(QLp))) - 150*traceAdjYdTYd*Sqr
      (gN)*Sqr(QQp) - 150*traceAdjYuTYu*Sqr(gN)*Sqr(QQp) + 150*MassBp*
      traceYuAdjYu*Sqr(gN)*Sqr(QQp) + 1800*Power(gN,4)*MassBp*Sqr(QH1p)*Sqr(QQp
      ) + 1800*Power(gN,4)*MassBp*Sqr(QH2p)*Sqr(QQp) - 10*traceYdAdjYd*(MassB*
      Sqr(g1) - 5*(8*MassG*Sqr(g3) + 3*MassBp*Sqr(gN)*(Sqr(Qdp) - Sqr(QH1p) +
      Sqr(QQp)))) + 150*traceAdjKappaTKappa*Sqr(gN)*Sqr(QSp) + 100*
      traceAdjLambda12TLambda12*Sqr(gN)*Sqr(QSp) - 150*MassBp*
      traceKappaAdjKappa*Sqr(gN)*Sqr(QSp) - 100*MassBp*traceLambda12AdjLambda12
      *Sqr(gN)*Sqr(QSp) + 900*Power(gN,4)*MassBp*Sqr(Qdp)*Sqr(QSp) + 900*Power(
      gN,4)*MassBp*Sqr(QDxbarp)*Sqr(QSp) + 900*Power(gN,4)*MassBp*Sqr(QDxp)*Sqr
      (QSp) + 300*Power(gN,4)*MassBp*Sqr(Qep)*Sqr(QSp) + 900*Power(gN,4)*MassBp
      *Sqr(QH1p)*Sqr(QSp) + 900*Power(gN,4)*MassBp*Sqr(QH2p)*Sqr(QSp) + 200*
      Power(gN,4)*MassBp*Sqr(QHpbarp)*Sqr(QSp) + 200*Power(gN,4)*MassBp*Sqr(
      QHpp)*Sqr(QSp) + 600*Power(gN,4)*MassBp*Sqr(QLp)*Sqr(QSp) + 1800*Power(gN
      ,4)*MassBp*Sqr(QQp)*Sqr(QSp) - 150*traceAdjYuTYu*Sqr(gN)*Sqr(Qup) + 150*
      MassBp*traceYuAdjYu*Sqr(gN)*Sqr(Qup) + 900*Power(gN,4)*MassBp*Sqr(QH1p)*
      Sqr(Qup) + 900*Power(gN,4)*MassBp*Sqr(QH2p)*Sqr(Qup) + 900*Power(gN,4)*
      MassBp*Sqr(QSp)*Sqr(Qup)) + (5.94*Power(g1,4) + 16.5*Power(g2,4) + 16*
      Power(gN,4)*Power(QH1p,4) + 16*Power(gN,4)*Power(QH2p,4) + 10*Power(gN,4)
      *Power(QSp,4) - 6*traceKappaAdjKappaKappaAdjKappa - 4*
      traceLambda12AdjLambda12Lambda12AdjLambda12 - 9*traceYdAdjYdYdAdjYd - 6*
      traceYdAdjYuYuAdjYd - 3*traceYeAdjYeYeAdjYe - 9*traceYuAdjYuYuAdjYu + 0.8
      *traceKappaAdjKappa*Sqr(g1) + 1.2*traceLambda12AdjLambda12*Sqr(g1) + 0.8*
      traceYuAdjYu*Sqr(g1) + 6*traceLambda12AdjLambda12*Sqr(g2) + 1.8*Sqr(g1)*
      Sqr(g2) + 16*traceKappaAdjKappa*Sqr(g3) + 16*traceYuAdjYu*Sqr(g3) + 6*
      traceKappaAdjKappa*Sqr(gN)*Sqr(QDxbarp) + 6*traceKappaAdjKappa*Sqr(gN)*
      Sqr(QDxp) + 4*traceLambda12AdjLambda12*Sqr(gN)*Sqr(QH1p) + 1.2*Sqr(g1)*
      Sqr(gN)*Sqr(QH1p) + 6*Sqr(g2)*Sqr(gN)*Sqr(QH1p) + 18*Power(gN,4)*Sqr(Qdp)
      *Sqr(QH1p) + 18*Power(gN,4)*Sqr(QDxbarp)*Sqr(QH1p) + 18*Power(gN,4)*Sqr(
      QDxp)*Sqr(QH1p) + 6*Power(gN,4)*Sqr(Qep)*Sqr(QH1p) + 4*
      traceLambda12AdjLambda12*Sqr(gN)*Sqr(QH2p) - 6*traceYuAdjYu*Sqr(gN)*Sqr(
      QH2p) + 1.2*Sqr(g1)*Sqr(gN)*Sqr(QH2p) + 6*Sqr(g2)*Sqr(gN)*Sqr(QH2p) + 18*
      Power(gN,4)*Sqr(Qdp)*Sqr(QH2p) + 18*Power(gN,4)*Sqr(QDxbarp)*Sqr(QH2p) +
      18*Power(gN,4)*Sqr(QDxp)*Sqr(QH2p) + 6*Power(gN,4)*Sqr(Qep)*Sqr(QH2p) +
      24*Power(gN,4)*Sqr(QH1p)*Sqr(QH2p) + 4*Power(gN,4)*Sqr(QH1p)*Sqr(QHpbarp)
      + 4*Power(gN,4)*Sqr(QH2p)*Sqr(QHpbarp) + 4*Power(gN,4)*Sqr(QH1p)*Sqr(
      QHpp) + 4*Power(gN,4)*Sqr(QH2p)*Sqr(QHpp) + 12*Power(gN,4)*Sqr(QH1p)*Sqr(
      QLp) + 12*Power(gN,4)*Sqr(QH2p)*Sqr(QLp) + 0.4*traceYeAdjYe*(3*Sqr(g1) +
      5*Sqr(gN)*(Sqr(Qep) - Sqr(QH1p) + Sqr(QLp))) + 6*traceYuAdjYu*Sqr(gN)*Sqr
      (QQp) + 36*Power(gN,4)*Sqr(QH1p)*Sqr(QQp) + 36*Power(gN,4)*Sqr(QH2p)*Sqr(
      QQp) - 0.4*traceYdAdjYd*(Sqr(g1) - 5*(8*Sqr(g3) + 3*Sqr(gN)*(Sqr(Qdp) -
      Sqr(QH1p) + Sqr(QQp)))) - 6*traceKappaAdjKappa*Sqr(gN)*Sqr(QSp) - 4*
      traceLambda12AdjLambda12*Sqr(gN)*Sqr(QSp) + 18*Power(gN,4)*Sqr(Qdp)*Sqr(
      QSp) + 18*Power(gN,4)*Sqr(QDxbarp)*Sqr(QSp) + 18*Power(gN,4)*Sqr(QDxp)*
      Sqr(QSp) + 6*Power(gN,4)*Sqr(Qep)*Sqr(QSp) + 18*Power(gN,4)*Sqr(QH1p)*Sqr
      (QSp) + 18*Power(gN,4)*Sqr(QH2p)*Sqr(QSp) + 4*Power(gN,4)*Sqr(QHpbarp)*
      Sqr(QSp) + 4*Power(gN,4)*Sqr(QHpp)*Sqr(QSp) + 12*Power(gN,4)*Sqr(QLp)*Sqr
      (QSp) + 36*Power(gN,4)*Sqr(QQp)*Sqr(QSp) + 6*traceYuAdjYu*Sqr(gN)*Sqr(Qup
      ) + 18*Power(gN,4)*Sqr(QH1p)*Sqr(Qup) + 18*Power(gN,4)*Sqr(QH2p)*Sqr(Qup)
      + 18*Power(gN,4)*Sqr(QSp)*Sqr(Qup))*TLambdax - 50*Sqr(Conj(Lambdax))*Sqr
      (Lambdax)*TLambdax - 0.2*AbsSqr(Lambdax)*(2*Lambdax*(30*
      traceAdjKappaTKappa + 20*traceAdjLambda12TLambda12 + 45*traceAdjYdTYd +
      15*traceAdjYeTYe + 45*traceAdjYuTYu + 6*MassB*Sqr(g1) + 30*MassWB*Sqr(g2)
      + 20*MassBp*Sqr(gN)*Sqr(QH1p) + 20*MassBp*Sqr(gN)*Sqr(QH2p)) - 3*(-30*
      traceKappaAdjKappa - 20*traceLambda12AdjLambda12 - 45*traceYdAdjYd - 15*
      traceYeAdjYe - 45*traceYuAdjYu + 6*Sqr(g1) + 30*Sqr(g2) + 20*Sqr(gN)*Sqr(
      QH1p) + 20*Sqr(gN)*Sqr(QH2p))*TLambdax));


   return beta_TLambdax;
}

} // namespace flexiblesusy
