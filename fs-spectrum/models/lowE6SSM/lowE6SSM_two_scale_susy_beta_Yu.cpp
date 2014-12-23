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

// File generated at Sun 24 Aug 2014 16:08:37

#include "lowE6SSM_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of Yu.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> lowE6SSM_susy_parameters::calc_beta_Yu_one_loop(const Susy_traces& susy_traces) const
{
   const auto QH2p = INPUT(QH2p);
   const auto QQp = INPUT(QQp);
   const auto Qup = INPUT(Qup);
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = oneOver16PiSqr*(Yu*(3*traceYuAdjYu + AbsSqr(Lambdax) -
      0.8666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) - 2*
      Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup)) + Yu*
      Yd.adjoint()*Yd + 3*(Yu*Yu.adjoint()*Yu));


   return beta_Yu;
}

/**
 * Calculates the two-loop beta function of Yu.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> lowE6SSM_susy_parameters::calc_beta_Yu_two_loop(const Susy_traces& susy_traces) const
{
   const auto Qdp = INPUT(Qdp);
   const auto QH2p = INPUT(QH2p);
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
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = twoLoop*(Yu*(8.695555555555556*Power(g1,4) + 16.5*Power(g2,4
      ) + 14.222222222222221*Power(g3,4) + 16*Power(gN,4)*Power(QH2p,4) + 40*
      Power(gN,4)*Power(QQp,4) + 22*Power(gN,4)*Power(Qup,4) - 3*
      traceYdAdjYuYuAdjYd - 9*traceYuAdjYuYuAdjYu + Sqr(g1)*Sqr(g2) +
      3.022222222222222*Sqr(g1)*Sqr(g3) + 8*Sqr(g2)*Sqr(g3) + 3.6*Qdp*QH2p*Sqr(
      g1)*Sqr(gN) + 3.6*QDxbarp*QH2p*Sqr(g1)*Sqr(gN) - 3.6*QDxp*QH2p*Sqr(g1)*
      Sqr(gN) + 3.6*Qep*QH2p*Sqr(g1)*Sqr(gN) - 3.6*QH1p*QH2p*Sqr(g1)*Sqr(gN) +
      1.2*QH2p*QHpbarp*Sqr(g1)*Sqr(gN) - 1.2*QH2p*QHpp*Sqr(g1)*Sqr(gN) - 3.6*
      QH2p*QLp*Sqr(g1)*Sqr(gN) + 1.2*Qdp*QQp*Sqr(g1)*Sqr(gN) + 1.2*QDxbarp*QQp*
      Sqr(g1)*Sqr(gN) - 1.2*QDxp*QQp*Sqr(g1)*Sqr(gN) + 1.2*Qep*QQp*Sqr(g1)*Sqr(
      gN) - 1.2*QH1p*QQp*Sqr(g1)*Sqr(gN) + 4.8*QH2p*QQp*Sqr(g1)*Sqr(gN) + 0.4*
      QHpbarp*QQp*Sqr(g1)*Sqr(gN) - 0.4*QHpp*QQp*Sqr(g1)*Sqr(gN) - 1.2*QLp*QQp*
      Sqr(g1)*Sqr(gN) - 4.8*Qdp*Qup*Sqr(g1)*Sqr(gN) - 4.8*QDxbarp*Qup*Sqr(g1)*
      Sqr(gN) + 4.8*QDxp*Qup*Sqr(g1)*Sqr(gN) - 4.8*Qep*Qup*Sqr(g1)*Sqr(gN) +
      4.8*QH1p*Qup*Sqr(g1)*Sqr(gN) - 12*QH2p*Qup*Sqr(g1)*Sqr(gN) - 1.6*QHpbarp*
      Qup*Sqr(g1)*Sqr(gN) + 1.6*QHpp*Qup*Sqr(g1)*Sqr(gN) + 4.8*QLp*Qup*Sqr(g1)*
      Sqr(gN) - 7.2*QQp*Qup*Sqr(g1)*Sqr(gN) + 4.8*Sqr(g1)*Sqr(gN)*Sqr(QH2p) + 6
      *Sqr(g2)*Sqr(gN)*Sqr(QH2p) + 18*Power(gN,4)*Sqr(Qdp)*Sqr(QH2p) + 18*Power
      (gN,4)*Sqr(QDxbarp)*Sqr(QH2p) + 18*Power(gN,4)*Sqr(QDxp)*Sqr(QH2p) + 6*
      Power(gN,4)*Sqr(Qep)*Sqr(QH2p) + 12*Power(gN,4)*Sqr(QH1p)*Sqr(QH2p) + 4*
      Power(gN,4)*Sqr(QH2p)*Sqr(QHpbarp) + 4*Power(gN,4)*Sqr(QH2p)*Sqr(QHpp) +
      12*Power(gN,4)*Sqr(QH2p)*Sqr(QLp) + 1.3333333333333333*Sqr(g1)*Sqr(gN)*
      Sqr(QQp) + 6*Sqr(g2)*Sqr(gN)*Sqr(QQp) + 10.666666666666666*Sqr(g3)*Sqr(gN
      )*Sqr(QQp) + 18*Power(gN,4)*Sqr(Qdp)*Sqr(QQp) + 18*Power(gN,4)*Sqr(
      QDxbarp)*Sqr(QQp) + 18*Power(gN,4)*Sqr(QDxp)*Sqr(QQp) + 6*Power(gN,4)*Sqr
      (Qep)*Sqr(QQp) + 12*Power(gN,4)*Sqr(QH1p)*Sqr(QQp) + 48*Power(gN,4)*Sqr(
      QH2p)*Sqr(QQp) + 4*Power(gN,4)*Sqr(QHpbarp)*Sqr(QQp) + 4*Power(gN,4)*Sqr(
      QHpp)*Sqr(QQp) + 12*Power(gN,4)*Sqr(QLp)*Sqr(QQp) + 6*Power(gN,4)*Sqr(
      QH2p)*Sqr(QSp) + 6*Power(gN,4)*Sqr(QQp)*Sqr(QSp) + AbsSqr(Lambdax)*(-3*
      traceKappaAdjKappa - 2*traceLambda12AdjLambda12 - 3*traceYdAdjYd -
      traceYeAdjYe + 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QH2p) + 2*Sqr(gN)*Sqr(
      QSp)) + 11.733333333333333*Sqr(g1)*Sqr(gN)*Sqr(Qup) + 10.666666666666666*
      Sqr(g3)*Sqr(gN)*Sqr(Qup) + 18*Power(gN,4)*Sqr(Qdp)*Sqr(Qup) + 18*Power(gN
      ,4)*Sqr(QDxbarp)*Sqr(Qup) + 18*Power(gN,4)*Sqr(QDxp)*Sqr(Qup) + 6*Power(
      gN,4)*Sqr(Qep)*Sqr(Qup) + 12*Power(gN,4)*Sqr(QH1p)*Sqr(Qup) + 30*Power(gN
      ,4)*Sqr(QH2p)*Sqr(Qup) + 4*Power(gN,4)*Sqr(QHpbarp)*Sqr(Qup) + 4*Power(gN
      ,4)*Sqr(QHpp)*Sqr(Qup) + 12*Power(gN,4)*Sqr(QLp)*Sqr(Qup) + 54*Power(gN,4
      )*Sqr(QQp)*Sqr(Qup) + 6*Power(gN,4)*Sqr(QSp)*Sqr(Qup) + 0.4*traceYuAdjYu*
      (2*Sqr(g1) + 5*(8*Sqr(g3) + 3*Sqr(gN)*(-Sqr(QH2p) + Sqr(QQp) + Sqr(Qup)))
      ) - 3*Sqr(Conj(Lambdax))*Sqr(Lambdax)) + (-3*traceYdAdjYd - traceYeAdjYe
      - AbsSqr(Lambdax) + 0.4*Sqr(g1) + 2*Sqr(gN)*Sqr(Qdp) + 2*Sqr(gN)*Sqr(QH1p
      ) - 2*Sqr(gN)*Sqr(QQp))*(Yu*Yd.adjoint()*Yd) - 9*traceYuAdjYu*(Yu*
      Yu.adjoint()*Yu) - 3*AbsSqr(Lambdax)*(Yu*Yu.adjoint()*Yu) + 0.4*Sqr(g1)*(
      Yu*Yu.adjoint()*Yu) + 6*Sqr(g2)*(Yu*Yu.adjoint()*Yu) + 6*Sqr(gN)*Sqr(QH2p
      )*(Yu*Yu.adjoint()*Yu) + 2*Sqr(gN)*Sqr(QQp)*(Yu*Yu.adjoint()*Yu) - 2*Sqr(
      gN)*Sqr(Qup)*(Yu*Yu.adjoint()*Yu) - 2*(Yu*Yd.adjoint()*Yd*Yd.adjoint()*Yd
      ) - 2*(Yu*Yd.adjoint()*Yd*Yu.adjoint()*Yu) - 4*(Yu*Yu.adjoint()*Yu*
      Yu.adjoint()*Yu));


   return beta_Yu;
}

} // namespace flexiblesusy