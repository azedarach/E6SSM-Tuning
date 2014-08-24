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

// File generated at Sun 24 Aug 2014 16:10:02

#include "lowE6SSM_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of TYd.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> lowE6SSM_soft_parameters::calc_beta_TYd_one_loop(const Soft_traces& soft_traces) const
{
   const auto Qdp = INPUT(Qdp);
   const auto QH1p = INPUT(QH1p);
   const auto QQp = INPUT(QQp);
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   Eigen::Matrix<double,3,3> beta_TYd;

   beta_TYd = oneOver16PiSqr*(3*traceYdAdjYd*TYd + traceYeAdjYe*TYd +
      AbsSqr(Lambdax)*TYd - 0.4666666666666667*Sqr(g1)*TYd - 3*Sqr(g2)*TYd -
      5.333333333333333*Sqr(g3)*TYd - 2*Sqr(gN)*Sqr(Qdp)*TYd - 2*Sqr(gN)*Sqr(
      QH1p)*TYd - 2*Sqr(gN)*Sqr(QQp)*TYd + Yd*(6*traceAdjYdTYd + 2*
      traceAdjYeTYe + 0.9333333333333333*MassB*Sqr(g1) + 6*MassWB*Sqr(g2) +
      10.666666666666666*MassG*Sqr(g3) + 4*MassBp*Sqr(gN)*Sqr(Qdp) + 4*MassBp*
      Sqr(gN)*Sqr(QH1p) + 4*MassBp*Sqr(gN)*Sqr(QQp) + 2*Conj(Lambdax)*TLambdax)
      + 4*(Yd*Yd.adjoint()*TYd) + 2*(Yd*Yu.adjoint()*TYu) + 5*(TYd*Yd.adjoint(
      )*Yd) + TYd*Yu.adjoint()*Yu);


   return beta_TYd;
}

/**
 * Calculates the two-loop beta function of TYd.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> lowE6SSM_soft_parameters::calc_beta_TYd_two_loop(const Soft_traces& soft_traces) const
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
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceAdjKappaTKappa = TRACE_STRUCT.traceAdjKappaTKappa;
   const double traceAdjLambda12TLambda12 =
      TRACE_STRUCT.traceAdjLambda12TLambda12;
   const double traceYdAdjYdTYdAdjYd = TRACE_STRUCT.traceYdAdjYdTYdAdjYd;
   const double traceYdAdjYuTYuAdjYd = TRACE_STRUCT.traceYdAdjYuTYuAdjYd;
   const double traceYeAdjYeTYeAdjYe = TRACE_STRUCT.traceYeAdjYeTYeAdjYe;
   const double traceYuAdjYdTYdAdjYu = TRACE_STRUCT.traceYuAdjYdTYdAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;


   Eigen::Matrix<double,3,3> beta_TYd;

   beta_TYd = twoLoop*(4.588888888888889*Power(g1,4)*TYd + 16.5*Power(g2,
      4)*TYd + 14.222222222222221*Power(g3,4)*TYd + 22*Power(gN,4)*Power(Qdp,4)
      *TYd + 16*Power(gN,4)*Power(QH1p,4)*TYd + 40*Power(gN,4)*Power(QQp,4)*TYd
      - 9*traceYdAdjYdYdAdjYd*TYd - 3*traceYdAdjYuYuAdjYd*TYd - 3*
      traceYeAdjYeYeAdjYe*TYd - 3*traceKappaAdjKappa*AbsSqr(Lambdax)*TYd - 2*
      traceLambda12AdjLambda12*AbsSqr(Lambdax)*TYd - 3*traceYuAdjYu*AbsSqr(
      Lambdax)*TYd - 0.4*traceYdAdjYd*Sqr(g1)*TYd + 1.2*traceYeAdjYe*Sqr(g1)*
      TYd + Sqr(g1)*Sqr(g2)*TYd + 16*traceYdAdjYd*Sqr(g3)*TYd +
      0.8888888888888888*Sqr(g1)*Sqr(g3)*TYd + 8*Sqr(g2)*Sqr(g3)*TYd + 6*
      traceYdAdjYd*Sqr(gN)*Sqr(Qdp)*TYd + 0.5333333333333333*Sqr(g1)*Sqr(gN)*
      Sqr(Qdp)*TYd + 10.666666666666666*Sqr(g3)*Sqr(gN)*Sqr(Qdp)*TYd + 18*Power
      (gN,4)*Sqr(Qdp)*Sqr(QDxbarp)*TYd + 18*Power(gN,4)*Sqr(Qdp)*Sqr(QDxp)*TYd
      + 2*traceYeAdjYe*Sqr(gN)*Sqr(Qep)*TYd + 6*Power(gN,4)*Sqr(Qdp)*Sqr(Qep)*
      TYd - 6*traceYdAdjYd*Sqr(gN)*Sqr(QH1p)*TYd - 2*traceYeAdjYe*Sqr(gN)*Sqr(
      QH1p)*TYd - 2*AbsSqr(Lambdax)*Sqr(gN)*Sqr(QH1p)*TYd + 1.2*Sqr(g1)*Sqr(gN)
      *Sqr(QH1p)*TYd + 6*Sqr(g2)*Sqr(gN)*Sqr(QH1p)*TYd + 30*Power(gN,4)*Sqr(Qdp
      )*Sqr(QH1p)*TYd + 18*Power(gN,4)*Sqr(QDxbarp)*Sqr(QH1p)*TYd + 18*Power(gN
      ,4)*Sqr(QDxp)*Sqr(QH1p)*TYd + 6*Power(gN,4)*Sqr(Qep)*Sqr(QH1p)*TYd + 2*
      AbsSqr(Lambdax)*Sqr(gN)*Sqr(QH2p)*TYd + 12*Power(gN,4)*Sqr(Qdp)*Sqr(QH2p)
      *TYd + 12*Power(gN,4)*Sqr(QH1p)*Sqr(QH2p)*TYd + 4*Power(gN,4)*Sqr(Qdp)*
      Sqr(QHpbarp)*TYd + 4*Power(gN,4)*Sqr(QH1p)*Sqr(QHpbarp)*TYd + 4*Power(gN,
      4)*Sqr(Qdp)*Sqr(QHpp)*TYd + 4*Power(gN,4)*Sqr(QH1p)*Sqr(QHpp)*TYd + 2*
      traceYeAdjYe*Sqr(gN)*Sqr(QLp)*TYd + 12*Power(gN,4)*Sqr(Qdp)*Sqr(QLp)*TYd
      + 12*Power(gN,4)*Sqr(QH1p)*Sqr(QLp)*TYd + 6*traceYdAdjYd*Sqr(gN)*Sqr(QQp)
      *TYd + 0.13333333333333333*Sqr(g1)*Sqr(gN)*Sqr(QQp)*TYd + 6*Sqr(g2)*Sqr(
      gN)*Sqr(QQp)*TYd + 10.666666666666666*Sqr(g3)*Sqr(gN)*Sqr(QQp)*TYd + 54*
      Power(gN,4)*Sqr(Qdp)*Sqr(QQp)*TYd + 18*Power(gN,4)*Sqr(QDxbarp)*Sqr(QQp)*
      TYd + 18*Power(gN,4)*Sqr(QDxp)*Sqr(QQp)*TYd + 6*Power(gN,4)*Sqr(Qep)*Sqr(
      QQp)*TYd + 48*Power(gN,4)*Sqr(QH1p)*Sqr(QQp)*TYd + 12*Power(gN,4)*Sqr(
      QH2p)*Sqr(QQp)*TYd + 4*Power(gN,4)*Sqr(QHpbarp)*Sqr(QQp)*TYd + 4*Power(gN
      ,4)*Sqr(QHpp)*Sqr(QQp)*TYd + 12*Power(gN,4)*Sqr(QLp)*Sqr(QQp)*TYd + 2*
      AbsSqr(Lambdax)*Sqr(gN)*Sqr(QSp)*TYd + 6*Power(gN,4)*Sqr(Qdp)*Sqr(QSp)*
      TYd + 6*Power(gN,4)*Sqr(QH1p)*Sqr(QSp)*TYd + 6*Power(gN,4)*Sqr(QQp)*Sqr(
      QSp)*TYd + 18*Power(gN,4)*Sqr(Qdp)*Sqr(Qup)*TYd + 18*Power(gN,4)*Sqr(QH1p
      )*Sqr(Qup)*TYd + 18*Power(gN,4)*Sqr(QQp)*Sqr(Qup)*TYd - 3*Sqr(Conj(
      Lambdax))*Sqr(Lambdax)*TYd - 0.044444444444444446*Yd*(413*Power(g1,4)*
      MassB + 1280*Power(g3,4)*MassG + 1485*Power(g2,4)*MassWB + 1980*Power(gN,
      4)*MassBp*Power(Qdp,4) + 1440*Power(gN,4)*MassBp*Power(QH1p,4) + 3600*
      Power(gN,4)*MassBp*Power(QQp,4) + 810*traceYdAdjYdTYdAdjYd + 135*
      traceYdAdjYuTYuAdjYd + 270*traceYeAdjYeTYeAdjYe + 135*
      traceYuAdjYdTYdAdjYu + 18*traceAdjYdTYd*Sqr(g1) - 54*traceAdjYeTYe*Sqr(g1
      ) + 54*MassB*traceYeAdjYe*Sqr(g1) + 45*MassB*Sqr(g1)*Sqr(g2) + 45*MassWB*
      Sqr(g1)*Sqr(g2) - 720*traceAdjYdTYd*Sqr(g3) + 40*MassB*Sqr(g1)*Sqr(g3) +
      40*MassG*Sqr(g1)*Sqr(g3) + 360*MassG*Sqr(g2)*Sqr(g3) + 360*MassWB*Sqr(g2)
      *Sqr(g3) - 270*traceAdjYdTYd*Sqr(gN)*Sqr(Qdp) + 24*MassB*Sqr(g1)*Sqr(gN)*
      Sqr(Qdp) + 24*MassBp*Sqr(g1)*Sqr(gN)*Sqr(Qdp) + 480*MassBp*Sqr(g3)*Sqr(gN
      )*Sqr(Qdp) + 480*MassG*Sqr(g3)*Sqr(gN)*Sqr(Qdp) + 1620*Power(gN,4)*MassBp
      *Sqr(Qdp)*Sqr(QDxbarp) + 1620*Power(gN,4)*MassBp*Sqr(Qdp)*Sqr(QDxp) - 90*
      traceAdjYeTYe*Sqr(gN)*Sqr(Qep) + 90*MassBp*traceYeAdjYe*Sqr(gN)*Sqr(Qep)
      + 540*Power(gN,4)*MassBp*Sqr(Qdp)*Sqr(Qep) + 270*traceAdjYdTYd*Sqr(gN)*
      Sqr(QH1p) + 90*traceAdjYeTYe*Sqr(gN)*Sqr(QH1p) - 90*MassBp*traceYeAdjYe*
      Sqr(gN)*Sqr(QH1p) + 54*MassB*Sqr(g1)*Sqr(gN)*Sqr(QH1p) + 54*MassBp*Sqr(g1
      )*Sqr(gN)*Sqr(QH1p) + 270*MassBp*Sqr(g2)*Sqr(gN)*Sqr(QH1p) + 270*MassWB*
      Sqr(g2)*Sqr(gN)*Sqr(QH1p) + 2700*Power(gN,4)*MassBp*Sqr(Qdp)*Sqr(QH1p) +
      1620*Power(gN,4)*MassBp*Sqr(QDxbarp)*Sqr(QH1p) + 1620*Power(gN,4)*MassBp*
      Sqr(QDxp)*Sqr(QH1p) + 540*Power(gN,4)*MassBp*Sqr(Qep)*Sqr(QH1p) + 1080*
      Power(gN,4)*MassBp*Sqr(Qdp)*Sqr(QH2p) + 1080*Power(gN,4)*MassBp*Sqr(QH1p)
      *Sqr(QH2p) + 360*Power(gN,4)*MassBp*Sqr(Qdp)*Sqr(QHpbarp) + 360*Power(gN,
      4)*MassBp*Sqr(QH1p)*Sqr(QHpbarp) + 360*Power(gN,4)*MassBp*Sqr(Qdp)*Sqr(
      QHpp) + 360*Power(gN,4)*MassBp*Sqr(QH1p)*Sqr(QHpp) - 90*traceAdjYeTYe*Sqr
      (gN)*Sqr(QLp) + 90*MassBp*traceYeAdjYe*Sqr(gN)*Sqr(QLp) + 1080*Power(gN,4
      )*MassBp*Sqr(Qdp)*Sqr(QLp) + 1080*Power(gN,4)*MassBp*Sqr(QH1p)*Sqr(QLp) -
      270*traceAdjYdTYd*Sqr(gN)*Sqr(QQp) + 6*MassB*Sqr(g1)*Sqr(gN)*Sqr(QQp) +
      6*MassBp*Sqr(g1)*Sqr(gN)*Sqr(QQp) + 270*MassBp*Sqr(g2)*Sqr(gN)*Sqr(QQp) +
      270*MassWB*Sqr(g2)*Sqr(gN)*Sqr(QQp) + 480*MassBp*Sqr(g3)*Sqr(gN)*Sqr(QQp
      ) + 480*MassG*Sqr(g3)*Sqr(gN)*Sqr(QQp) + 4860*Power(gN,4)*MassBp*Sqr(Qdp)
      *Sqr(QQp) + 1620*Power(gN,4)*MassBp*Sqr(QDxbarp)*Sqr(QQp) + 1620*Power(gN
      ,4)*MassBp*Sqr(QDxp)*Sqr(QQp) + 540*Power(gN,4)*MassBp*Sqr(Qep)*Sqr(QQp)
      + 4320*Power(gN,4)*MassBp*Sqr(QH1p)*Sqr(QQp) + 1080*Power(gN,4)*MassBp*
      Sqr(QH2p)*Sqr(QQp) + 360*Power(gN,4)*MassBp*Sqr(QHpbarp)*Sqr(QQp) + 360*
      Power(gN,4)*MassBp*Sqr(QHpp)*Sqr(QQp) + 1080*Power(gN,4)*MassBp*Sqr(QLp)*
      Sqr(QQp) - 18*traceYdAdjYd*(MassB*Sqr(g1) - 5*(8*MassG*Sqr(g3) + 3*MassBp
      *Sqr(gN)*(Sqr(Qdp) - Sqr(QH1p) + Sqr(QQp)))) + 540*Power(gN,4)*MassBp*Sqr
      (Qdp)*Sqr(QSp) + 540*Power(gN,4)*MassBp*Sqr(QH1p)*Sqr(QSp) + 540*Power(gN
      ,4)*MassBp*Sqr(QQp)*Sqr(QSp) + 1620*Power(gN,4)*MassBp*Sqr(Qdp)*Sqr(Qup)
      + 1620*Power(gN,4)*MassBp*Sqr(QH1p)*Sqr(Qup) + 1620*Power(gN,4)*MassBp*
      Sqr(QQp)*Sqr(Qup) + 270*Lambdax*Sqr(Conj(Lambdax))*TLambdax - 45*Conj(
      Lambdax)*(Lambdax*(-3*traceAdjKappaTKappa - 2*traceAdjLambda12TLambda12 -
      3*traceAdjYuTYu + 2*MassBp*Sqr(gN)*Sqr(QH1p) - 2*MassBp*Sqr(gN)*Sqr(QH2p
      ) - 2*MassBp*Sqr(gN)*Sqr(QSp)) - (3*traceKappaAdjKappa + 2*
      traceLambda12AdjLambda12 + 3*traceYuAdjYu + 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(
      gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QSp))*TLambdax)) - 0.4*(45*traceAdjYdTYd +
      15*traceAdjYeTYe + 4*MassB*Sqr(g1) + 30*MassWB*Sqr(g2) - 10*MassBp*Sqr(gN
      )*Sqr(Qdp) + 30*MassBp*Sqr(gN)*Sqr(QH1p) + 10*MassBp*Sqr(gN)*Sqr(QQp) +
      15*Conj(Lambdax)*TLambdax)*(Yd*Yd.adjoint()*Yd) - 12*traceYdAdjYd*(Yd*
      Yd.adjoint()*TYd) - 4*traceYeAdjYe*(Yd*Yd.adjoint()*TYd) - 4*AbsSqr(
      Lambdax)*(Yd*Yd.adjoint()*TYd) + 1.2*Sqr(g1)*(Yd*Yd.adjoint()*TYd) + 6*
      Sqr(g2)*(Yd*Yd.adjoint()*TYd) + 8*Sqr(gN)*Sqr(QH1p)*(Yd*Yd.adjoint()*TYd)
      - 6*traceAdjYuTYu*(Yd*Yu.adjoint()*Yu) - 1.6*MassB*Sqr(g1)*(Yd*
      Yu.adjoint()*Yu) - 4*MassBp*Sqr(gN)*Sqr(QH2p)*(Yd*Yu.adjoint()*Yu) + 4*
      MassBp*Sqr(gN)*Sqr(QQp)*(Yd*Yu.adjoint()*Yu) - 4*MassBp*Sqr(gN)*Sqr(Qup)*
      (Yd*Yu.adjoint()*Yu) - 2*Conj(Lambdax)*TLambdax*(Yd*Yu.adjoint()*Yu) - 6*
      traceYuAdjYu*(Yd*Yu.adjoint()*TYu) - 2*AbsSqr(Lambdax)*(Yd*Yu.adjoint()*
      TYu) + 1.6*Sqr(g1)*(Yd*Yu.adjoint()*TYu) + 4*Sqr(gN)*Sqr(QH2p)*(Yd*
      Yu.adjoint()*TYu) - 4*Sqr(gN)*Sqr(QQp)*(Yd*Yu.adjoint()*TYu) + 4*Sqr(gN)*
      Sqr(Qup)*(Yd*Yu.adjoint()*TYu) - 15*traceYdAdjYd*(TYd*Yd.adjoint()*Yd) -
      5*traceYeAdjYe*(TYd*Yd.adjoint()*Yd) - 5*AbsSqr(Lambdax)*(TYd*Yd.adjoint(
      )*Yd) + 1.2*Sqr(g1)*(TYd*Yd.adjoint()*Yd) + 12*Sqr(g2)*(TYd*Yd.adjoint()*
      Yd) - 6*Sqr(gN)*Sqr(Qdp)*(TYd*Yd.adjoint()*Yd) + 10*Sqr(gN)*Sqr(QH1p)*(
      TYd*Yd.adjoint()*Yd) + 6*Sqr(gN)*Sqr(QQp)*(TYd*Yd.adjoint()*Yd) - 3*
      traceYuAdjYu*(TYd*Yu.adjoint()*Yu) - AbsSqr(Lambdax)*(TYd*Yu.adjoint()*Yu
      ) + 0.8*Sqr(g1)*(TYd*Yu.adjoint()*Yu) + 2*Sqr(gN)*Sqr(QH2p)*(TYd*
      Yu.adjoint()*Yu) - 2*Sqr(gN)*Sqr(QQp)*(TYd*Yu.adjoint()*Yu) + 2*Sqr(gN)*
      Sqr(Qup)*(TYd*Yu.adjoint()*Yu) - 6*(Yd*Yd.adjoint()*Yd*Yd.adjoint()*TYd)
      - 8*(Yd*Yd.adjoint()*TYd*Yd.adjoint()*Yd) - 2*(Yd*Yu.adjoint()*Yu*
      Yd.adjoint()*TYd) - 4*(Yd*Yu.adjoint()*Yu*Yu.adjoint()*TYu) - 4*(Yd*
      Yu.adjoint()*TYu*Yd.adjoint()*Yd) - 4*(Yd*Yu.adjoint()*TYu*Yu.adjoint()*
      Yu) - 6*(TYd*Yd.adjoint()*Yd*Yd.adjoint()*Yd) - 4*(TYd*Yu.adjoint()*Yu*
      Yd.adjoint()*Yd) - 2*(TYd*Yu.adjoint()*Yu*Yu.adjoint()*Yu));


   return beta_TYd;
}

} // namespace flexiblesusy
