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
 * Calculates the one-loop beta function of TYu.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> lowE6SSM_soft_parameters::calc_beta_TYu_one_loop(const Soft_traces& soft_traces) const
{
   const auto QH2p = INPUT(QH2p);
   const auto QQp = INPUT(QQp);
   const auto Qup = INPUT(Qup);
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   Eigen::Matrix<double,3,3> beta_TYu;

   beta_TYu = oneOver16PiSqr*(3*traceYuAdjYu*TYu + AbsSqr(Lambdax)*TYu -
      0.8666666666666667*Sqr(g1)*TYu - 3*Sqr(g2)*TYu - 5.333333333333333*Sqr(g3
      )*TYu - 2*Sqr(gN)*Sqr(QH2p)*TYu - 2*Sqr(gN)*Sqr(QQp)*TYu - 2*Sqr(gN)*Sqr(
      Qup)*TYu + Yu*(6*traceAdjYuTYu + 1.7333333333333334*MassB*Sqr(g1) + 6*
      MassWB*Sqr(g2) + 10.666666666666666*MassG*Sqr(g3) + 4*MassBp*Sqr(gN)*Sqr(
      QH2p) + 4*MassBp*Sqr(gN)*Sqr(QQp) + 4*MassBp*Sqr(gN)*Sqr(Qup) + 2*Conj(
      Lambdax)*TLambdax) + 2*(Yu*Yd.adjoint()*TYd) + 4*(Yu*Yu.adjoint()*TYu) +
      TYu*Yd.adjoint()*Yd + 5*(TYu*Yu.adjoint()*Yu));


   return beta_TYu;
}

/**
 * Calculates the two-loop beta function of TYu.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> lowE6SSM_soft_parameters::calc_beta_TYu_two_loop(const Soft_traces& soft_traces) const
{
   const auto QH2p = INPUT(QH2p);
   const auto Qdp = INPUT(Qdp);
   const auto QDxbarp = INPUT(QDxbarp);
   const auto QDxp = INPUT(QDxp);
   const auto Qep = INPUT(Qep);
   const auto QH1p = INPUT(QH1p);
   const auto QHpbarp = INPUT(QHpbarp);
   const auto QHpp = INPUT(QHpp);
   const auto QLp = INPUT(QLp);
   const auto QQp = INPUT(QQp);
   const auto QSp = INPUT(QSp);
   const auto Qup = INPUT(Qup);
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjKappaTKappa = TRACE_STRUCT.traceAdjKappaTKappa;
   const double traceAdjLambda12TLambda12 =
      TRACE_STRUCT.traceAdjLambda12TLambda12;
   const double traceYdAdjYuTYuAdjYd = TRACE_STRUCT.traceYdAdjYuTYuAdjYd;
   const double traceYuAdjYdTYdAdjYu = TRACE_STRUCT.traceYuAdjYdTYdAdjYu;
   const double traceYuAdjYuTYuAdjYu = TRACE_STRUCT.traceYuAdjYuTYuAdjYu;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   Eigen::Matrix<double,3,3> beta_TYu;

   beta_TYu = twoLoop*(8.695555555555556*Power(g1,4)*TYu + 16.5*Power(g2,
      4)*TYu + 14.222222222222221*Power(g3,4)*TYu + 16*Power(gN,4)*Power(QH2p,4
      )*TYu + 40*Power(gN,4)*Power(QQp,4)*TYu + 22*Power(gN,4)*Power(Qup,4)*TYu
      - 3*traceYdAdjYuYuAdjYd*TYu - 9*traceYuAdjYuYuAdjYu*TYu - 3*
      traceKappaAdjKappa*AbsSqr(Lambdax)*TYu - 2*traceLambda12AdjLambda12*
      AbsSqr(Lambdax)*TYu - 3*traceYdAdjYd*AbsSqr(Lambdax)*TYu - traceYeAdjYe*
      AbsSqr(Lambdax)*TYu + 0.8*traceYuAdjYu*Sqr(g1)*TYu + Sqr(g1)*Sqr(g2)*TYu
      + 16*traceYuAdjYu*Sqr(g3)*TYu + 3.022222222222222*Sqr(g1)*Sqr(g3)*TYu + 8
      *Sqr(g2)*Sqr(g3)*TYu + 2*AbsSqr(Lambdax)*Sqr(gN)*Sqr(QH1p)*TYu - 6*
      traceYuAdjYu*Sqr(gN)*Sqr(QH2p)*TYu - 2*AbsSqr(Lambdax)*Sqr(gN)*Sqr(QH2p)*
      TYu + 1.2*Sqr(g1)*Sqr(gN)*Sqr(QH2p)*TYu + 6*Sqr(g2)*Sqr(gN)*Sqr(QH2p)*TYu
      + 18*Power(gN,4)*Sqr(Qdp)*Sqr(QH2p)*TYu + 18*Power(gN,4)*Sqr(QDxbarp)*
      Sqr(QH2p)*TYu + 18*Power(gN,4)*Sqr(QDxp)*Sqr(QH2p)*TYu + 6*Power(gN,4)*
      Sqr(Qep)*Sqr(QH2p)*TYu + 12*Power(gN,4)*Sqr(QH1p)*Sqr(QH2p)*TYu + 4*Power
      (gN,4)*Sqr(QH2p)*Sqr(QHpbarp)*TYu + 4*Power(gN,4)*Sqr(QH2p)*Sqr(QHpp)*TYu
      + 12*Power(gN,4)*Sqr(QH2p)*Sqr(QLp)*TYu + 6*traceYuAdjYu*Sqr(gN)*Sqr(QQp
      )*TYu + 0.13333333333333333*Sqr(g1)*Sqr(gN)*Sqr(QQp)*TYu + 6*Sqr(g2)*Sqr(
      gN)*Sqr(QQp)*TYu + 10.666666666666666*Sqr(g3)*Sqr(gN)*Sqr(QQp)*TYu + 18*
      Power(gN,4)*Sqr(Qdp)*Sqr(QQp)*TYu + 18*Power(gN,4)*Sqr(QDxbarp)*Sqr(QQp)*
      TYu + 18*Power(gN,4)*Sqr(QDxp)*Sqr(QQp)*TYu + 6*Power(gN,4)*Sqr(Qep)*Sqr(
      QQp)*TYu + 12*Power(gN,4)*Sqr(QH1p)*Sqr(QQp)*TYu + 48*Power(gN,4)*Sqr(
      QH2p)*Sqr(QQp)*TYu + 4*Power(gN,4)*Sqr(QHpbarp)*Sqr(QQp)*TYu + 4*Power(gN
      ,4)*Sqr(QHpp)*Sqr(QQp)*TYu + 12*Power(gN,4)*Sqr(QLp)*Sqr(QQp)*TYu + 2*
      AbsSqr(Lambdax)*Sqr(gN)*Sqr(QSp)*TYu + 6*Power(gN,4)*Sqr(QH2p)*Sqr(QSp)*
      TYu + 6*Power(gN,4)*Sqr(QQp)*Sqr(QSp)*TYu + 6*traceYuAdjYu*Sqr(gN)*Sqr(
      Qup)*TYu + 2.1333333333333333*Sqr(g1)*Sqr(gN)*Sqr(Qup)*TYu +
      10.666666666666666*Sqr(g3)*Sqr(gN)*Sqr(Qup)*TYu + 18*Power(gN,4)*Sqr(Qdp)
      *Sqr(Qup)*TYu + 18*Power(gN,4)*Sqr(QDxbarp)*Sqr(Qup)*TYu + 18*Power(gN,4)
      *Sqr(QDxp)*Sqr(Qup)*TYu + 6*Power(gN,4)*Sqr(Qep)*Sqr(Qup)*TYu + 12*Power(
      gN,4)*Sqr(QH1p)*Sqr(Qup)*TYu + 30*Power(gN,4)*Sqr(QH2p)*Sqr(Qup)*TYu + 4*
      Power(gN,4)*Sqr(QHpbarp)*Sqr(Qup)*TYu + 4*Power(gN,4)*Sqr(QHpp)*Sqr(Qup)*
      TYu + 12*Power(gN,4)*Sqr(QLp)*Sqr(Qup)*TYu + 54*Power(gN,4)*Sqr(QQp)*Sqr(
      Qup)*TYu + 6*Power(gN,4)*Sqr(QSp)*Sqr(Qup)*TYu - 3*Sqr(Conj(Lambdax))*Sqr
      (Lambdax)*TYu - 0.008888888888888889*Yu*(3913*Power(g1,4)*MassB + 6400*
      Power(g3,4)*MassG + 7425*Power(g2,4)*MassWB + 7200*Power(gN,4)*MassBp*
      Power(QH2p,4) + 18000*Power(gN,4)*MassBp*Power(QQp,4) + 9900*Power(gN,4)*
      MassBp*Power(Qup,4) + 675*traceYdAdjYuTYuAdjYd + 675*traceYuAdjYdTYdAdjYu
      + 4050*traceYuAdjYuTYuAdjYu - 180*traceAdjYuTYu*Sqr(g1) + 225*MassB*Sqr(
      g1)*Sqr(g2) + 225*MassWB*Sqr(g1)*Sqr(g2) - 3600*traceAdjYuTYu*Sqr(g3) +
      680*MassB*Sqr(g1)*Sqr(g3) + 680*MassG*Sqr(g1)*Sqr(g3) + 1800*MassG*Sqr(g2
      )*Sqr(g3) + 1800*MassWB*Sqr(g2)*Sqr(g3) + 1350*traceAdjYuTYu*Sqr(gN)*Sqr(
      QH2p) + 270*MassB*Sqr(g1)*Sqr(gN)*Sqr(QH2p) + 270*MassBp*Sqr(g1)*Sqr(gN)*
      Sqr(QH2p) + 1350*MassBp*Sqr(g2)*Sqr(gN)*Sqr(QH2p) + 1350*MassWB*Sqr(g2)*
      Sqr(gN)*Sqr(QH2p) + 8100*Power(gN,4)*MassBp*Sqr(Qdp)*Sqr(QH2p) + 8100*
      Power(gN,4)*MassBp*Sqr(QDxbarp)*Sqr(QH2p) + 8100*Power(gN,4)*MassBp*Sqr(
      QDxp)*Sqr(QH2p) + 2700*Power(gN,4)*MassBp*Sqr(Qep)*Sqr(QH2p) + 5400*Power
      (gN,4)*MassBp*Sqr(QH1p)*Sqr(QH2p) + 1800*Power(gN,4)*MassBp*Sqr(QH2p)*Sqr
      (QHpbarp) + 1800*Power(gN,4)*MassBp*Sqr(QH2p)*Sqr(QHpp) + 5400*Power(gN,4
      )*MassBp*Sqr(QH2p)*Sqr(QLp) - 1350*traceAdjYuTYu*Sqr(gN)*Sqr(QQp) + 30*
      MassB*Sqr(g1)*Sqr(gN)*Sqr(QQp) + 30*MassBp*Sqr(g1)*Sqr(gN)*Sqr(QQp) +
      1350*MassBp*Sqr(g2)*Sqr(gN)*Sqr(QQp) + 1350*MassWB*Sqr(g2)*Sqr(gN)*Sqr(
      QQp) + 2400*MassBp*Sqr(g3)*Sqr(gN)*Sqr(QQp) + 2400*MassG*Sqr(g3)*Sqr(gN)*
      Sqr(QQp) + 8100*Power(gN,4)*MassBp*Sqr(Qdp)*Sqr(QQp) + 8100*Power(gN,4)*
      MassBp*Sqr(QDxbarp)*Sqr(QQp) + 8100*Power(gN,4)*MassBp*Sqr(QDxp)*Sqr(QQp)
      + 2700*Power(gN,4)*MassBp*Sqr(Qep)*Sqr(QQp) + 5400*Power(gN,4)*MassBp*
      Sqr(QH1p)*Sqr(QQp) + 21600*Power(gN,4)*MassBp*Sqr(QH2p)*Sqr(QQp) + 1800*
      Power(gN,4)*MassBp*Sqr(QHpbarp)*Sqr(QQp) + 1800*Power(gN,4)*MassBp*Sqr(
      QHpp)*Sqr(QQp) + 5400*Power(gN,4)*MassBp*Sqr(QLp)*Sqr(QQp) + 2700*Power(
      gN,4)*MassBp*Sqr(QH2p)*Sqr(QSp) + 2700*Power(gN,4)*MassBp*Sqr(QQp)*Sqr(
      QSp) - 1350*traceAdjYuTYu*Sqr(gN)*Sqr(Qup) + 480*MassB*Sqr(g1)*Sqr(gN)*
      Sqr(Qup) + 480*MassBp*Sqr(g1)*Sqr(gN)*Sqr(Qup) + 2400*MassBp*Sqr(g3)*Sqr(
      gN)*Sqr(Qup) + 2400*MassG*Sqr(g3)*Sqr(gN)*Sqr(Qup) + 8100*Power(gN,4)*
      MassBp*Sqr(Qdp)*Sqr(Qup) + 8100*Power(gN,4)*MassBp*Sqr(QDxbarp)*Sqr(Qup)
      + 8100*Power(gN,4)*MassBp*Sqr(QDxp)*Sqr(Qup) + 2700*Power(gN,4)*MassBp*
      Sqr(Qep)*Sqr(Qup) + 5400*Power(gN,4)*MassBp*Sqr(QH1p)*Sqr(Qup) + 13500*
      Power(gN,4)*MassBp*Sqr(QH2p)*Sqr(Qup) + 1800*Power(gN,4)*MassBp*Sqr(
      QHpbarp)*Sqr(Qup) + 1800*Power(gN,4)*MassBp*Sqr(QHpp)*Sqr(Qup) + 5400*
      Power(gN,4)*MassBp*Sqr(QLp)*Sqr(Qup) + 24300*Power(gN,4)*MassBp*Sqr(QQp)*
      Sqr(Qup) + 2700*Power(gN,4)*MassBp*Sqr(QSp)*Sqr(Qup) + 90*traceYuAdjYu*(2
      *MassB*Sqr(g1) + 5*(8*MassG*Sqr(g3) + 3*MassBp*Sqr(gN)*(-Sqr(QH2p) + Sqr(
      QQp) + Sqr(Qup)))) + 1350*Lambdax*Sqr(Conj(Lambdax))*TLambdax + 225*Conj(
      Lambdax)*(Lambdax*(3*traceAdjKappaTKappa + 2*traceAdjLambda12TLambda12 +
      3*traceAdjYdTYd + traceAdjYeTYe + 2*MassBp*Sqr(gN)*Sqr(QH1p) - 2*MassBp*
      Sqr(gN)*Sqr(QH2p) + 2*MassBp*Sqr(gN)*Sqr(QSp)) + (3*traceKappaAdjKappa +
      2*traceLambda12AdjLambda12 + 3*traceYdAdjYd + traceYeAdjYe - 2*Sqr(gN)*
      Sqr(QH1p) + 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QSp))*TLambdax)) + (-6*
      traceAdjYdTYd - 2*traceAdjYeTYe - 0.8*MassB*Sqr(g1) - 4*MassBp*Sqr(gN)*
      Sqr(Qdp) - 4*MassBp*Sqr(gN)*Sqr(QH1p) + 4*MassBp*Sqr(gN)*Sqr(QQp) - 2*
      Conj(Lambdax)*TLambdax)*(Yu*Yd.adjoint()*Yd) - 6*traceYdAdjYd*(Yu*
      Yd.adjoint()*TYd) - 2*traceYeAdjYe*(Yu*Yd.adjoint()*TYd) - 2*AbsSqr(
      Lambdax)*(Yu*Yd.adjoint()*TYd) + 0.8*Sqr(g1)*(Yu*Yd.adjoint()*TYd) + 4*
      Sqr(gN)*Sqr(Qdp)*(Yu*Yd.adjoint()*TYd) + 4*Sqr(gN)*Sqr(QH1p)*(Yu*
      Yd.adjoint()*TYd) - 4*Sqr(gN)*Sqr(QQp)*(Yu*Yd.adjoint()*TYd) - 18*
      traceAdjYuTYu*(Yu*Yu.adjoint()*Yu) - 0.8*MassB*Sqr(g1)*(Yu*Yu.adjoint()*
      Yu) - 12*MassWB*Sqr(g2)*(Yu*Yu.adjoint()*Yu) - 12*MassBp*Sqr(gN)*Sqr(QH2p
      )*(Yu*Yu.adjoint()*Yu) - 4*MassBp*Sqr(gN)*Sqr(QQp)*(Yu*Yu.adjoint()*Yu) +
      4*MassBp*Sqr(gN)*Sqr(Qup)*(Yu*Yu.adjoint()*Yu) - 6*Conj(Lambdax)*
      TLambdax*(Yu*Yu.adjoint()*Yu) - 12*traceYuAdjYu*(Yu*Yu.adjoint()*TYu) - 4
      *AbsSqr(Lambdax)*(Yu*Yu.adjoint()*TYu) + 1.2*Sqr(g1)*(Yu*Yu.adjoint()*TYu
      ) + 6*Sqr(g2)*(Yu*Yu.adjoint()*TYu) + 8*Sqr(gN)*Sqr(QH2p)*(Yu*Yu.adjoint(
      )*TYu) - 3*traceYdAdjYd*(TYu*Yd.adjoint()*Yd) - traceYeAdjYe*(TYu*
      Yd.adjoint()*Yd) - AbsSqr(Lambdax)*(TYu*Yd.adjoint()*Yd) + 0.4*Sqr(g1)*(
      TYu*Yd.adjoint()*Yd) + 2*Sqr(gN)*Sqr(Qdp)*(TYu*Yd.adjoint()*Yd) + 2*Sqr(
      gN)*Sqr(QH1p)*(TYu*Yd.adjoint()*Yd) - 2*Sqr(gN)*Sqr(QQp)*(TYu*Yd.adjoint(
      )*Yd) - 15*traceYuAdjYu*(TYu*Yu.adjoint()*Yu) - 5*AbsSqr(Lambdax)*(TYu*
      Yu.adjoint()*Yu) + 12*Sqr(g2)*(TYu*Yu.adjoint()*Yu) + 10*Sqr(gN)*Sqr(QH2p
      )*(TYu*Yu.adjoint()*Yu) + 6*Sqr(gN)*Sqr(QQp)*(TYu*Yu.adjoint()*Yu) - 6*
      Sqr(gN)*Sqr(Qup)*(TYu*Yu.adjoint()*Yu) - 4*(Yu*Yd.adjoint()*Yd*Yd.adjoint
      ()*TYd) - 2*(Yu*Yd.adjoint()*Yd*Yu.adjoint()*TYu) - 4*(Yu*Yd.adjoint()*
      TYd*Yd.adjoint()*Yd) - 4*(Yu*Yd.adjoint()*TYd*Yu.adjoint()*Yu) - 6*(Yu*
      Yu.adjoint()*Yu*Yu.adjoint()*TYu) - 8*(Yu*Yu.adjoint()*TYu*Yu.adjoint()*
      Yu) - 2*(TYu*Yd.adjoint()*Yd*Yd.adjoint()*Yd) - 4*(TYu*Yd.adjoint()*Yd*
      Yu.adjoint()*Yu) - 6*(TYu*Yu.adjoint()*Yu*Yu.adjoint()*Yu));


   return beta_TYu;
}

} // namespace flexiblesusy
