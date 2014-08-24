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

// File generated at Sun 24 Aug 2014 16:08:35

#include "lowE6SSM_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of Yd.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> lowE6SSM_susy_parameters::calc_beta_Yd_one_loop(const Susy_traces& susy_traces) const
{
   const auto Qdp = INPUT(Qdp);
   const auto QH1p = INPUT(QH1p);
   const auto QQp = INPUT(QQp);
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   Eigen::Matrix<double,3,3> beta_Yd;

   beta_Yd = oneOver16PiSqr*(Yd*(3*traceYdAdjYd + traceYeAdjYe + AbsSqr(
      Lambdax) - 0.4666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr
      (g3) - 2*Sqr(gN)*Sqr(Qdp) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QQp)) + 3
      *(Yd*Yd.adjoint()*Yd) + Yd*Yu.adjoint()*Yu);


   return beta_Yd;
}

/**
 * Calculates the two-loop beta function of Yd.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> lowE6SSM_susy_parameters::calc_beta_Yd_two_loop(const Susy_traces& susy_traces) const
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
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;


   Eigen::Matrix<double,3,3> beta_Yd;

   beta_Yd = twoLoop*(Yd*(4.588888888888889*Power(g1,4) + 16.5*Power(g2,4
      ) + 14.222222222222221*Power(g3,4) + 22*Power(gN,4)*Power(Qdp,4) + 16*
      Power(gN,4)*Power(QH1p,4) + 40*Power(gN,4)*Power(QQp,4) - 9*
      traceYdAdjYdYdAdjYd - 3*traceYdAdjYuYuAdjYd - 3*traceYeAdjYeYeAdjYe + 1.2
      *traceYeAdjYe*Sqr(g1) + Sqr(g1)*Sqr(g2) + 0.8888888888888888*Sqr(g1)*Sqr(
      g3) + 8*Sqr(g2)*Sqr(g3) + 2.4*Qdp*QDxbarp*Sqr(g1)*Sqr(gN) - 2.4*Qdp*QDxp*
      Sqr(g1)*Sqr(gN) + 2.4*Qdp*Qep*Sqr(g1)*Sqr(gN) - 6*Qdp*QH1p*Sqr(g1)*Sqr(gN
      ) - 3.6*QDxbarp*QH1p*Sqr(g1)*Sqr(gN) + 3.6*QDxp*QH1p*Sqr(g1)*Sqr(gN) -
      3.6*Qep*QH1p*Sqr(g1)*Sqr(gN) + 2.4*Qdp*QH2p*Sqr(g1)*Sqr(gN) - 3.6*QH1p*
      QH2p*Sqr(g1)*Sqr(gN) + 0.8*Qdp*QHpbarp*Sqr(g1)*Sqr(gN) - 1.2*QH1p*QHpbarp
      *Sqr(g1)*Sqr(gN) - 0.8*Qdp*QHpp*Sqr(g1)*Sqr(gN) + 1.2*QH1p*QHpp*Sqr(g1)*
      Sqr(gN) - 2.4*Qdp*QLp*Sqr(g1)*Sqr(gN) + 3.6*QH1p*QLp*Sqr(g1)*Sqr(gN) +
      3.6*Qdp*QQp*Sqr(g1)*Sqr(gN) + 1.2*QDxbarp*QQp*Sqr(g1)*Sqr(gN) - 1.2*QDxp*
      QQp*Sqr(g1)*Sqr(gN) + 1.2*Qep*QQp*Sqr(g1)*Sqr(gN) - 4.8*QH1p*QQp*Sqr(g1)*
      Sqr(gN) + 1.2*QH2p*QQp*Sqr(g1)*Sqr(gN) + 0.4*QHpbarp*QQp*Sqr(g1)*Sqr(gN)
      - 0.4*QHpp*QQp*Sqr(g1)*Sqr(gN) - 1.2*QLp*QQp*Sqr(g1)*Sqr(gN) - 4.8*Qdp*
      Qup*Sqr(g1)*Sqr(gN) + 7.2*QH1p*Qup*Sqr(g1)*Sqr(gN) - 2.4*QQp*Qup*Sqr(g1)*
      Sqr(gN) + 2.933333333333333*Sqr(g1)*Sqr(gN)*Sqr(Qdp) + 10.666666666666666
      *Sqr(g3)*Sqr(gN)*Sqr(Qdp) + 18*Power(gN,4)*Sqr(Qdp)*Sqr(QDxbarp) + 18*
      Power(gN,4)*Sqr(Qdp)*Sqr(QDxp) + 2*traceYeAdjYe*Sqr(gN)*Sqr(Qep) + 6*
      Power(gN,4)*Sqr(Qdp)*Sqr(Qep) - 2*traceYeAdjYe*Sqr(gN)*Sqr(QH1p) + 4.8*
      Sqr(g1)*Sqr(gN)*Sqr(QH1p) + 6*Sqr(g2)*Sqr(gN)*Sqr(QH1p) + 30*Power(gN,4)*
      Sqr(Qdp)*Sqr(QH1p) + 18*Power(gN,4)*Sqr(QDxbarp)*Sqr(QH1p) + 18*Power(gN,
      4)*Sqr(QDxp)*Sqr(QH1p) + 6*Power(gN,4)*Sqr(Qep)*Sqr(QH1p) + 12*Power(gN,4
      )*Sqr(Qdp)*Sqr(QH2p) + 12*Power(gN,4)*Sqr(QH1p)*Sqr(QH2p) + 4*Power(gN,4)
      *Sqr(Qdp)*Sqr(QHpbarp) + 4*Power(gN,4)*Sqr(QH1p)*Sqr(QHpbarp) + 4*Power(
      gN,4)*Sqr(Qdp)*Sqr(QHpp) + 4*Power(gN,4)*Sqr(QH1p)*Sqr(QHpp) + 2*
      traceYeAdjYe*Sqr(gN)*Sqr(QLp) + 12*Power(gN,4)*Sqr(Qdp)*Sqr(QLp) + 12*
      Power(gN,4)*Sqr(QH1p)*Sqr(QLp) + 1.3333333333333333*Sqr(g1)*Sqr(gN)*Sqr(
      QQp) + 6*Sqr(g2)*Sqr(gN)*Sqr(QQp) + 10.666666666666666*Sqr(g3)*Sqr(gN)*
      Sqr(QQp) + 54*Power(gN,4)*Sqr(Qdp)*Sqr(QQp) + 18*Power(gN,4)*Sqr(QDxbarp)
      *Sqr(QQp) + 18*Power(gN,4)*Sqr(QDxp)*Sqr(QQp) + 6*Power(gN,4)*Sqr(Qep)*
      Sqr(QQp) + 48*Power(gN,4)*Sqr(QH1p)*Sqr(QQp) + 12*Power(gN,4)*Sqr(QH2p)*
      Sqr(QQp) + 4*Power(gN,4)*Sqr(QHpbarp)*Sqr(QQp) + 4*Power(gN,4)*Sqr(QHpp)*
      Sqr(QQp) + 12*Power(gN,4)*Sqr(QLp)*Sqr(QQp) - 0.4*traceYdAdjYd*(Sqr(g1) -
      5*(8*Sqr(g3) + 3*Sqr(gN)*(Sqr(Qdp) - Sqr(QH1p) + Sqr(QQp)))) + 6*Power(
      gN,4)*Sqr(Qdp)*Sqr(QSp) + 6*Power(gN,4)*Sqr(QH1p)*Sqr(QSp) + 6*Power(gN,4
      )*Sqr(QQp)*Sqr(QSp) - AbsSqr(Lambdax)*(3*traceKappaAdjKappa + 2*
      traceLambda12AdjLambda12 + 3*traceYuAdjYu + 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(
      gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QSp)) + 18*Power(gN,4)*Sqr(Qdp)*Sqr(Qup) +
      18*Power(gN,4)*Sqr(QH1p)*Sqr(Qup) + 18*Power(gN,4)*Sqr(QQp)*Sqr(Qup) - 3*
      Sqr(Conj(Lambdax))*Sqr(Lambdax)) + (-9*traceYdAdjYd - 3*traceYeAdjYe - 3*
      AbsSqr(Lambdax) + 0.8*Sqr(g1) + 6*Sqr(g2) - 2*Sqr(gN)*Sqr(Qdp) + 6*Sqr(gN
      )*Sqr(QH1p) + 2*Sqr(gN)*Sqr(QQp))*(Yd*Yd.adjoint()*Yd) - 3*traceYuAdjYu*(
      Yd*Yu.adjoint()*Yu) - AbsSqr(Lambdax)*(Yd*Yu.adjoint()*Yu) + 0.8*Sqr(g1)*
      (Yd*Yu.adjoint()*Yu) + 2*Sqr(gN)*Sqr(QH2p)*(Yd*Yu.adjoint()*Yu) - 2*Sqr(
      gN)*Sqr(QQp)*(Yd*Yu.adjoint()*Yu) + 2*Sqr(gN)*Sqr(Qup)*(Yd*Yu.adjoint()*
      Yu) - 4*(Yd*Yd.adjoint()*Yd*Yd.adjoint()*Yd) - 2*(Yd*Yu.adjoint()*Yu*
      Yd.adjoint()*Yd) - 2*(Yd*Yu.adjoint()*Yu*Yu.adjoint()*Yu));


   return beta_Yd;
}

} // namespace flexiblesusy
