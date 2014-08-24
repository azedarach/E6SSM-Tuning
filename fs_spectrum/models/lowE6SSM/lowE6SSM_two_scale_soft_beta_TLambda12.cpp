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

// File generated at Sun 24 Aug 2014 16:10:05

#include "lowE6SSM_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of TLambda12.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,2,2> lowE6SSM_soft_parameters::calc_beta_TLambda12_one_loop(const Soft_traces& soft_traces) const
{
   const auto QH1p = INPUT(QH1p);
   const auto QH2p = INPUT(QH2p);
   const auto QSp = INPUT(QSp);
   const double traceAdjKappaTKappa = TRACE_STRUCT.traceAdjKappaTKappa;
   const double traceAdjLambda12TLambda12 =
      TRACE_STRUCT.traceAdjLambda12TLambda12;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;


   Eigen::Matrix<double,2,2> beta_TLambda12;

   beta_TLambda12 = oneOver16PiSqr*(Lambda12*(6*traceAdjKappaTKappa + 4*
      traceAdjLambda12TLambda12 + 1.2*MassB*Sqr(g1) + 6*MassWB*Sqr(g2) + 4*
      MassBp*Sqr(gN)*Sqr(QH1p) + 4*MassBp*Sqr(gN)*Sqr(QH2p) + 4*MassBp*Sqr(gN)*
      Sqr(QSp) + 4*Conj(Lambdax)*TLambdax) + 3*traceKappaAdjKappa*TLambda12 + 2
      *traceLambda12AdjLambda12*TLambda12 + 2*AbsSqr(Lambdax)*TLambda12 - 0.6*
      Sqr(g1)*TLambda12 - 3*Sqr(g2)*TLambda12 - 2*Sqr(gN)*Sqr(QH1p)*TLambda12 -
      2*Sqr(gN)*Sqr(QH2p)*TLambda12 - 2*Sqr(gN)*Sqr(QSp)*TLambda12 + 3*(
      Lambda12*(Lambda12).adjoint()*TLambda12) + 3*(TLambda12*(Lambda12)
      .adjoint()*Lambda12));


   return beta_TLambda12;
}

/**
 * Calculates the two-loop beta function of TLambda12.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,2,2> lowE6SSM_soft_parameters::calc_beta_TLambda12_two_loop(const Soft_traces& soft_traces) const
{
   const auto QH1p = INPUT(QH1p);
   const auto Qdp = INPUT(Qdp);
   const auto QDxbarp = INPUT(QDxbarp);
   const auto QDxp = INPUT(QDxp);
   const auto Qep = INPUT(Qep);
   const auto QH2p = INPUT(QH2p);
   const auto QHpbarp = INPUT(QHpbarp);
   const auto QHpp = INPUT(QHpp);
   const auto QLp = INPUT(QLp);
   const auto QQp = INPUT(QQp);
   const auto QSp = INPUT(QSp);
   const auto Qup = INPUT(Qup);
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceAdjKappaTKappa = TRACE_STRUCT.traceAdjKappaTKappa;
   const double traceAdjLambda12TLambda12 =
      TRACE_STRUCT.traceAdjLambda12TLambda12;
   const double traceKappaAdjKappaTKappaAdjKappa =
      TRACE_STRUCT.traceKappaAdjKappaTKappaAdjKappa;
   const double traceLambda12AdjLambda12TLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12TLambda12AdjLambda12;
   const double traceKappaAdjKappaKappaAdjKappa =
      TRACE_STRUCT.traceKappaAdjKappaKappaAdjKappa;
   const double traceLambda12AdjLambda12Lambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12Lambda12AdjLambda12;


   Eigen::Matrix<double,2,2> beta_TLambda12;

   beta_TLambda12 = twoLoop*(-0.08*Lambda12*(297*Power(g1,4)*MassB + 825*
      Power(g2,4)*MassWB + 800*Power(gN,4)*MassBp*Power(QH1p,4) + 800*Power(gN,
      4)*MassBp*Power(QH2p,4) + 500*Power(gN,4)*MassBp*Power(QSp,4) + 300*
      traceKappaAdjKappaTKappaAdjKappa + 200*
      traceLambda12AdjLambda12TLambda12AdjLambda12 - 20*traceAdjKappaTKappa*Sqr
      (g1) - 30*traceAdjLambda12TLambda12*Sqr(g1) + 30*MassB*
      traceLambda12AdjLambda12*Sqr(g1) - 150*traceAdjLambda12TLambda12*Sqr(g2)
      + 150*MassWB*traceLambda12AdjLambda12*Sqr(g2) + 45*MassB*Sqr(g1)*Sqr(g2)
      + 45*MassWB*Sqr(g1)*Sqr(g2) - 400*traceAdjKappaTKappa*Sqr(g3) - 150*
      traceAdjKappaTKappa*Sqr(gN)*Sqr(QDxbarp) - 150*traceAdjKappaTKappa*Sqr(gN
      )*Sqr(QDxp) - 100*traceAdjLambda12TLambda12*Sqr(gN)*Sqr(QH1p) + 100*
      MassBp*traceLambda12AdjLambda12*Sqr(gN)*Sqr(QH1p) + 30*MassB*Sqr(g1)*Sqr(
      gN)*Sqr(QH1p) + 30*MassBp*Sqr(g1)*Sqr(gN)*Sqr(QH1p) + 150*MassBp*Sqr(g2)*
      Sqr(gN)*Sqr(QH1p) + 150*MassWB*Sqr(g2)*Sqr(gN)*Sqr(QH1p) + 900*Power(gN,4
      )*MassBp*Sqr(Qdp)*Sqr(QH1p) + 900*Power(gN,4)*MassBp*Sqr(QDxbarp)*Sqr(
      QH1p) + 900*Power(gN,4)*MassBp*Sqr(QDxp)*Sqr(QH1p) + 300*Power(gN,4)*
      MassBp*Sqr(Qep)*Sqr(QH1p) - 100*traceAdjLambda12TLambda12*Sqr(gN)*Sqr(
      QH2p) + 100*MassBp*traceLambda12AdjLambda12*Sqr(gN)*Sqr(QH2p) + 30*MassB*
      Sqr(g1)*Sqr(gN)*Sqr(QH2p) + 30*MassBp*Sqr(g1)*Sqr(gN)*Sqr(QH2p) + 150*
      MassBp*Sqr(g2)*Sqr(gN)*Sqr(QH2p) + 150*MassWB*Sqr(g2)*Sqr(gN)*Sqr(QH2p) +
      900*Power(gN,4)*MassBp*Sqr(Qdp)*Sqr(QH2p) + 900*Power(gN,4)*MassBp*Sqr(
      QDxbarp)*Sqr(QH2p) + 900*Power(gN,4)*MassBp*Sqr(QDxp)*Sqr(QH2p) + 300*
      Power(gN,4)*MassBp*Sqr(Qep)*Sqr(QH2p) + 1200*Power(gN,4)*MassBp*Sqr(QH1p)
      *Sqr(QH2p) + 200*Power(gN,4)*MassBp*Sqr(QH1p)*Sqr(QHpbarp) + 200*Power(gN
      ,4)*MassBp*Sqr(QH2p)*Sqr(QHpbarp) + 200*Power(gN,4)*MassBp*Sqr(QH1p)*Sqr(
      QHpp) + 200*Power(gN,4)*MassBp*Sqr(QH2p)*Sqr(QHpp) + 600*Power(gN,4)*
      MassBp*Sqr(QH1p)*Sqr(QLp) + 600*Power(gN,4)*MassBp*Sqr(QH2p)*Sqr(QLp) +
      1800*Power(gN,4)*MassBp*Sqr(QH1p)*Sqr(QQp) + 1800*Power(gN,4)*MassBp*Sqr(
      QH2p)*Sqr(QQp) + 10*traceKappaAdjKappa*(2*MassB*Sqr(g1) + 5*(8*MassG*Sqr(
      g3) + 3*MassBp*Sqr(gN)*(Sqr(QDxbarp) + Sqr(QDxp) - Sqr(QSp)))) + 150*
      traceAdjKappaTKappa*Sqr(gN)*Sqr(QSp) + 100*traceAdjLambda12TLambda12*Sqr(
      gN)*Sqr(QSp) - 100*MassBp*traceLambda12AdjLambda12*Sqr(gN)*Sqr(QSp) + 900
      *Power(gN,4)*MassBp*Sqr(Qdp)*Sqr(QSp) + 900*Power(gN,4)*MassBp*Sqr(
      QDxbarp)*Sqr(QSp) + 900*Power(gN,4)*MassBp*Sqr(QDxp)*Sqr(QSp) + 300*Power
      (gN,4)*MassBp*Sqr(Qep)*Sqr(QSp) + 900*Power(gN,4)*MassBp*Sqr(QH1p)*Sqr(
      QSp) + 900*Power(gN,4)*MassBp*Sqr(QH2p)*Sqr(QSp) + 200*Power(gN,4)*MassBp
      *Sqr(QHpbarp)*Sqr(QSp) + 200*Power(gN,4)*MassBp*Sqr(QHpp)*Sqr(QSp) + 600*
      Power(gN,4)*MassBp*Sqr(QLp)*Sqr(QSp) + 1800*Power(gN,4)*MassBp*Sqr(QQp)*
      Sqr(QSp) + 900*Power(gN,4)*MassBp*Sqr(QH1p)*Sqr(Qup) + 900*Power(gN,4)*
      MassBp*Sqr(QH2p)*Sqr(Qup) + 900*Power(gN,4)*MassBp*Sqr(QSp)*Sqr(Qup) +
      200*Lambdax*Sqr(Conj(Lambdax))*TLambdax + 10*Conj(Lambdax)*(Lambdax*(15*
      traceAdjYdTYd + 5*traceAdjYeTYe + 15*traceAdjYuTYu + 3*MassB*Sqr(g1) + 15
      *MassWB*Sqr(g2) + 10*MassBp*Sqr(gN)*Sqr(QH1p) + 10*MassBp*Sqr(gN)*Sqr(
      QH2p) - 10*MassBp*Sqr(gN)*Sqr(QSp)) + (15*traceYdAdjYd + 5*traceYeAdjYe +
      15*traceYuAdjYu - 3*Sqr(g1) - 15*Sqr(g2) - 10*Sqr(gN)*Sqr(QH1p) - 10*Sqr
      (gN)*Sqr(QH2p) + 10*Sqr(gN)*Sqr(QSp))*TLambdax)) + 5.94*Power(g1,4)*
      TLambda12 + 16.5*Power(g2,4)*TLambda12 + 16*Power(gN,4)*Power(QH1p,4)*
      TLambda12 + 16*Power(gN,4)*Power(QH2p,4)*TLambda12 + 10*Power(gN,4)*Power
      (QSp,4)*TLambda12 - 6*traceKappaAdjKappaKappaAdjKappa*TLambda12 - 4*
      traceLambda12AdjLambda12Lambda12AdjLambda12*TLambda12 - 6*traceYdAdjYd*
      AbsSqr(Lambdax)*TLambda12 - 2*traceYeAdjYe*AbsSqr(Lambdax)*TLambda12 - 6*
      traceYuAdjYu*AbsSqr(Lambdax)*TLambda12 + 0.8*traceKappaAdjKappa*Sqr(g1)*
      TLambda12 + 1.2*traceLambda12AdjLambda12*Sqr(g1)*TLambda12 + 1.2*AbsSqr(
      Lambdax)*Sqr(g1)*TLambda12 + 6*traceLambda12AdjLambda12*Sqr(g2)*TLambda12
      + 6*AbsSqr(Lambdax)*Sqr(g2)*TLambda12 + 1.8*Sqr(g1)*Sqr(g2)*TLambda12 +
      16*traceKappaAdjKappa*Sqr(g3)*TLambda12 + 6*traceKappaAdjKappa*Sqr(gN)*
      Sqr(QDxbarp)*TLambda12 + 6*traceKappaAdjKappa*Sqr(gN)*Sqr(QDxp)*TLambda12
      + 4*traceLambda12AdjLambda12*Sqr(gN)*Sqr(QH1p)*TLambda12 + 4*AbsSqr(
      Lambdax)*Sqr(gN)*Sqr(QH1p)*TLambda12 + 1.2*Sqr(g1)*Sqr(gN)*Sqr(QH1p)*
      TLambda12 + 6*Sqr(g2)*Sqr(gN)*Sqr(QH1p)*TLambda12 + 18*Power(gN,4)*Sqr(
      Qdp)*Sqr(QH1p)*TLambda12 + 18*Power(gN,4)*Sqr(QDxbarp)*Sqr(QH1p)*
      TLambda12 + 18*Power(gN,4)*Sqr(QDxp)*Sqr(QH1p)*TLambda12 + 6*Power(gN,4)*
      Sqr(Qep)*Sqr(QH1p)*TLambda12 + 4*traceLambda12AdjLambda12*Sqr(gN)*Sqr(
      QH2p)*TLambda12 + 4*AbsSqr(Lambdax)*Sqr(gN)*Sqr(QH2p)*TLambda12 + 1.2*Sqr
      (g1)*Sqr(gN)*Sqr(QH2p)*TLambda12 + 6*Sqr(g2)*Sqr(gN)*Sqr(QH2p)*TLambda12
      + 18*Power(gN,4)*Sqr(Qdp)*Sqr(QH2p)*TLambda12 + 18*Power(gN,4)*Sqr(
      QDxbarp)*Sqr(QH2p)*TLambda12 + 18*Power(gN,4)*Sqr(QDxp)*Sqr(QH2p)*
      TLambda12 + 6*Power(gN,4)*Sqr(Qep)*Sqr(QH2p)*TLambda12 + 24*Power(gN,4)*
      Sqr(QH1p)*Sqr(QH2p)*TLambda12 + 4*Power(gN,4)*Sqr(QH1p)*Sqr(QHpbarp)*
      TLambda12 + 4*Power(gN,4)*Sqr(QH2p)*Sqr(QHpbarp)*TLambda12 + 4*Power(gN,4
      )*Sqr(QH1p)*Sqr(QHpp)*TLambda12 + 4*Power(gN,4)*Sqr(QH2p)*Sqr(QHpp)*
      TLambda12 + 12*Power(gN,4)*Sqr(QH1p)*Sqr(QLp)*TLambda12 + 12*Power(gN,4)*
      Sqr(QH2p)*Sqr(QLp)*TLambda12 + 36*Power(gN,4)*Sqr(QH1p)*Sqr(QQp)*
      TLambda12 + 36*Power(gN,4)*Sqr(QH2p)*Sqr(QQp)*TLambda12 - 6*
      traceKappaAdjKappa*Sqr(gN)*Sqr(QSp)*TLambda12 - 4*
      traceLambda12AdjLambda12*Sqr(gN)*Sqr(QSp)*TLambda12 - 4*AbsSqr(Lambdax)*
      Sqr(gN)*Sqr(QSp)*TLambda12 + 18*Power(gN,4)*Sqr(Qdp)*Sqr(QSp)*TLambda12 +
      18*Power(gN,4)*Sqr(QDxbarp)*Sqr(QSp)*TLambda12 + 18*Power(gN,4)*Sqr(QDxp
      )*Sqr(QSp)*TLambda12 + 6*Power(gN,4)*Sqr(Qep)*Sqr(QSp)*TLambda12 + 18*
      Power(gN,4)*Sqr(QH1p)*Sqr(QSp)*TLambda12 + 18*Power(gN,4)*Sqr(QH2p)*Sqr(
      QSp)*TLambda12 + 4*Power(gN,4)*Sqr(QHpbarp)*Sqr(QSp)*TLambda12 + 4*Power(
      gN,4)*Sqr(QHpp)*Sqr(QSp)*TLambda12 + 12*Power(gN,4)*Sqr(QLp)*Sqr(QSp)*
      TLambda12 + 36*Power(gN,4)*Sqr(QQp)*Sqr(QSp)*TLambda12 + 18*Power(gN,4)*
      Sqr(QH1p)*Sqr(Qup)*TLambda12 + 18*Power(gN,4)*Sqr(QH2p)*Sqr(Qup)*
      TLambda12 + 18*Power(gN,4)*Sqr(QSp)*Sqr(Qup)*TLambda12 - 4*Sqr(Conj(
      Lambdax))*Sqr(Lambdax)*TLambda12 - 4*(3*traceAdjKappaTKappa + 2*(
      traceAdjLambda12TLambda12 + MassBp*Sqr(gN)*Sqr(QSp)) + 2*Conj(Lambdax)*
      TLambdax)*(Lambda12*(Lambda12).adjoint()*Lambda12) - 9*traceKappaAdjKappa
      *(Lambda12*(Lambda12).adjoint()*TLambda12) - 6*traceLambda12AdjLambda12*(
      Lambda12*(Lambda12).adjoint()*TLambda12) - 6*AbsSqr(Lambdax)*(Lambda12*(
      Lambda12).adjoint()*TLambda12) - 2*Sqr(gN)*Sqr(QH1p)*(Lambda12*(Lambda12)
      .adjoint()*TLambda12) + 2*Sqr(gN)*Sqr(QH2p)*(Lambda12*(Lambda12).adjoint(
      )*TLambda12) + 6*Sqr(gN)*Sqr(QSp)*(Lambda12*(Lambda12).adjoint()*
      TLambda12) - 9*traceKappaAdjKappa*(TLambda12*(Lambda12).adjoint()*
      Lambda12) - 6*traceLambda12AdjLambda12*(TLambda12*(Lambda12).adjoint()*
      Lambda12) - 6*AbsSqr(Lambdax)*(TLambda12*(Lambda12).adjoint()*Lambda12) +
      2*Sqr(gN)*Sqr(QH1p)*(TLambda12*(Lambda12).adjoint()*Lambda12) - 2*Sqr(gN
      )*Sqr(QH2p)*(TLambda12*(Lambda12).adjoint()*Lambda12) + 6*Sqr(gN)*Sqr(QSp
      )*(TLambda12*(Lambda12).adjoint()*Lambda12) - 3*(Lambda12*(Lambda12)
      .adjoint()*Lambda12*(Lambda12).adjoint()*TLambda12) - 4*(Lambda12*(
      Lambda12).adjoint()*TLambda12*(Lambda12).adjoint()*Lambda12) - 3*(
      TLambda12*(Lambda12).adjoint()*Lambda12*(Lambda12).adjoint()*Lambda12));


   return beta_TLambda12;
}

} // namespace flexiblesusy
