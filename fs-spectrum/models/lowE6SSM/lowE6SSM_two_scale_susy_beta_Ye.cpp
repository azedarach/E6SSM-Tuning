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
 * Calculates the one-loop beta function of Ye.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> lowE6SSM_susy_parameters::calc_beta_Ye_one_loop(const Susy_traces& susy_traces) const
{
   const auto Qep = INPUT(Qep);
   const auto QH1p = INPUT(QH1p);
   const auto QLp = INPUT(QLp);
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   Eigen::Matrix<double,3,3> beta_Ye;

   beta_Ye = oneOver16PiSqr*(Ye*(3*traceYdAdjYd + traceYeAdjYe + AbsSqr(
      Lambdax) - 1.8*Sqr(g1) - 3*Sqr(g2) - 2*Sqr(gN)*Sqr(Qep) - 2*Sqr(gN)*Sqr(
      QH1p) - 2*Sqr(gN)*Sqr(QLp)) + 3*(Ye*Ye.adjoint()*Ye));


   return beta_Ye;
}

/**
 * Calculates the two-loop beta function of Ye.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> lowE6SSM_susy_parameters::calc_beta_Ye_two_loop(const Susy_traces& susy_traces) const
{
   const auto Qdp = INPUT(Qdp);
   const auto Qep = INPUT(Qep);
   const auto QDxbarp = INPUT(QDxbarp);
   const auto QDxp = INPUT(QDxp);
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


   Eigen::Matrix<double,3,3> beta_Ye;

   beta_Ye = twoLoop*(0.1*Ye*(189*Power(g1,4) + 165*Power(g2,4) + 100*
      Power(gN,4)*Power(Qep,4) + 160*Power(gN,4)*Power(QH1p,4) + 160*Power(gN,4
      )*Power(QLp,4) - 90*traceYdAdjYdYdAdjYd - 30*traceYdAdjYuYuAdjYd - 30*
      traceYeAdjYeYeAdjYe + 12*traceYeAdjYe*Sqr(g1) + 18*Sqr(g1)*Sqr(g2) + 72*
      Qdp*Qep*Sqr(g1)*Sqr(gN) + 72*QDxbarp*Qep*Sqr(g1)*Sqr(gN) - 72*QDxp*Qep*
      Sqr(g1)*Sqr(gN) - 36*Qdp*QH1p*Sqr(g1)*Sqr(gN) - 36*QDxbarp*QH1p*Sqr(g1)*
      Sqr(gN) + 36*QDxp*QH1p*Sqr(g1)*Sqr(gN) - 108*Qep*QH1p*Sqr(g1)*Sqr(gN) +
      72*Qep*QH2p*Sqr(g1)*Sqr(gN) - 36*QH1p*QH2p*Sqr(g1)*Sqr(gN) + 24*Qep*
      QHpbarp*Sqr(g1)*Sqr(gN) - 12*QH1p*QHpbarp*Sqr(g1)*Sqr(gN) - 24*Qep*QHpp*
      Sqr(g1)*Sqr(gN) + 12*QH1p*QHpp*Sqr(g1)*Sqr(gN) - 36*Qdp*QLp*Sqr(g1)*Sqr(
      gN) - 36*QDxbarp*QLp*Sqr(g1)*Sqr(gN) + 36*QDxp*QLp*Sqr(g1)*Sqr(gN) - 108*
      Qep*QLp*Sqr(g1)*Sqr(gN) + 72*QH1p*QLp*Sqr(g1)*Sqr(gN) - 36*QH2p*QLp*Sqr(
      g1)*Sqr(gN) - 12*QHpbarp*QLp*Sqr(g1)*Sqr(gN) + 12*QHpp*QLp*Sqr(g1)*Sqr(gN
      ) + 72*Qep*QQp*Sqr(g1)*Sqr(gN) - 36*QH1p*QQp*Sqr(g1)*Sqr(gN) - 36*QLp*QQp
      *Sqr(g1)*Sqr(gN) - 144*Qep*Qup*Sqr(g1)*Sqr(gN) + 72*QH1p*Qup*Sqr(g1)*Sqr(
      gN) + 72*QLp*Qup*Sqr(g1)*Sqr(gN) + 20*traceYeAdjYe*Sqr(gN)*Sqr(Qep) + 120
      *Sqr(g1)*Sqr(gN)*Sqr(Qep) + 180*Power(gN,4)*Sqr(Qdp)*Sqr(Qep) + 180*Power
      (gN,4)*Sqr(QDxbarp)*Sqr(Qep) + 180*Power(gN,4)*Sqr(QDxp)*Sqr(Qep) - 20*
      traceYeAdjYe*Sqr(gN)*Sqr(QH1p) + 48*Sqr(g1)*Sqr(gN)*Sqr(QH1p) + 60*Sqr(g2
      )*Sqr(gN)*Sqr(QH1p) + 180*Power(gN,4)*Sqr(Qdp)*Sqr(QH1p) + 180*Power(gN,4
      )*Sqr(QDxbarp)*Sqr(QH1p) + 180*Power(gN,4)*Sqr(QDxp)*Sqr(QH1p) + 180*
      Power(gN,4)*Sqr(Qep)*Sqr(QH1p) + 120*Power(gN,4)*Sqr(Qep)*Sqr(QH2p) + 120
      *Power(gN,4)*Sqr(QH1p)*Sqr(QH2p) + 40*Power(gN,4)*Sqr(Qep)*Sqr(QHpbarp) +
      40*Power(gN,4)*Sqr(QH1p)*Sqr(QHpbarp) + 40*Power(gN,4)*Sqr(Qep)*Sqr(QHpp
      ) + 40*Power(gN,4)*Sqr(QH1p)*Sqr(QHpp) + 20*traceYeAdjYe*Sqr(gN)*Sqr(QLp)
      + 48*Sqr(g1)*Sqr(gN)*Sqr(QLp) + 60*Sqr(g2)*Sqr(gN)*Sqr(QLp) + 180*Power(
      gN,4)*Sqr(Qdp)*Sqr(QLp) + 180*Power(gN,4)*Sqr(QDxbarp)*Sqr(QLp) + 180*
      Power(gN,4)*Sqr(QDxp)*Sqr(QLp) + 180*Power(gN,4)*Sqr(Qep)*Sqr(QLp) + 240*
      Power(gN,4)*Sqr(QH1p)*Sqr(QLp) + 120*Power(gN,4)*Sqr(QH2p)*Sqr(QLp) + 40*
      Power(gN,4)*Sqr(QHpbarp)*Sqr(QLp) + 40*Power(gN,4)*Sqr(QHpp)*Sqr(QLp) +
      360*Power(gN,4)*Sqr(Qep)*Sqr(QQp) + 360*Power(gN,4)*Sqr(QH1p)*Sqr(QQp) +
      360*Power(gN,4)*Sqr(QLp)*Sqr(QQp) - 4*traceYdAdjYd*(Sqr(g1) - 5*(8*Sqr(g3
      ) + 3*Sqr(gN)*(Sqr(Qdp) - Sqr(QH1p) + Sqr(QQp)))) + 60*Power(gN,4)*Sqr(
      Qep)*Sqr(QSp) + 60*Power(gN,4)*Sqr(QH1p)*Sqr(QSp) + 60*Power(gN,4)*Sqr(
      QLp)*Sqr(QSp) - 10*AbsSqr(Lambdax)*(3*traceKappaAdjKappa + 2*
      traceLambda12AdjLambda12 + 3*traceYuAdjYu + 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(
      gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QSp)) + 180*Power(gN,4)*Sqr(Qep)*Sqr(Qup) +
      180*Power(gN,4)*Sqr(QH1p)*Sqr(Qup) + 180*Power(gN,4)*Sqr(QLp)*Sqr(Qup) -
      30*Sqr(Conj(Lambdax))*Sqr(Lambdax)) + (-9*traceYdAdjYd - 3*traceYeAdjYe
      - 3*AbsSqr(Lambdax) + 6*Sqr(g2) - 2*Sqr(gN)*Sqr(Qep) + 6*Sqr(gN)*Sqr(QH1p
      ) + 2*Sqr(gN)*Sqr(QLp))*(Ye*Ye.adjoint()*Ye) - 4*(Ye*Ye.adjoint()*Ye*
      Ye.adjoint()*Ye));


   return beta_Ye;
}

} // namespace flexiblesusy
