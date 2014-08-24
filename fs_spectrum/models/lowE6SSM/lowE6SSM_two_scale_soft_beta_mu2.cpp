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

// File generated at Sun 24 Aug 2014 16:10:15

#include "lowE6SSM_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of mu2.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> lowE6SSM_soft_parameters::calc_beta_mu2_one_loop(const Soft_traces& soft_traces) const
{
   const auto Qup = INPUT(Qup);
   const double Tr11 = TRACE_STRUCT.Tr11;
   const double Tr14 = TRACE_STRUCT.Tr14;


   Eigen::Matrix<double,3,3> beta_mu2;

   beta_mu2 = oneOver16PiSqr*(4*mHu2*(Yu*Yu.adjoint()) + 4*(TYu*(TYu)
      .adjoint()) + 2*(mu2*Yu*Yu.adjoint()) + 4*(Yu*mq2*Yu.adjoint()) + 2*(Yu*
      Yu.adjoint()*mu2) - 1.0327955589886444*g1*Tr11*UNITMATRIX(3) + 2*gN*Qup*
      Tr14*UNITMATRIX(3) - 2.1333333333333333*AbsSqr(MassB)*Sqr(g1)*UNITMATRIX(
      3) - 10.666666666666666*AbsSqr(MassG)*Sqr(g3)*UNITMATRIX(3) - 8*AbsSqr(
      MassBp)*Sqr(gN)*Sqr(Qup)*UNITMATRIX(3));


   return beta_mu2;
}

/**
 * Calculates the two-loop beta function of mu2.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> lowE6SSM_soft_parameters::calc_beta_mu2_two_loop(const Soft_traces& soft_traces) const
{
   const auto Qup = INPUT(Qup);
   const auto QH2p = INPUT(QH2p);
   const auto QQp = INPUT(QQp);
   const auto Qdp = INPUT(Qdp);
   const auto QDxbarp = INPUT(QDxbarp);
   const auto QDxp = INPUT(QDxp);
   const auto Qep = INPUT(Qep);
   const auto QH1p = INPUT(QH1p);
   const auto QHpbarp = INPUT(QHpbarp);
   const auto QHpp = INPUT(QHpp);
   const auto QLp = INPUT(QLp);
   const auto QSp = INPUT(QSp);
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceconjTYuTpTYu = TRACE_STRUCT.traceconjTYuTpTYu;
   const double tracemq2AdjYuYu = TRACE_STRUCT.tracemq2AdjYuYu;
   const double tracemu2YuAdjYu = TRACE_STRUCT.tracemu2YuAdjYu;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceconjTYuTpYu = TRACE_STRUCT.traceconjTYuTpYu;
   const double Tr2U111 = TRACE_STRUCT.Tr2U111;
   const double Tr2U114 = TRACE_STRUCT.Tr2U114;
   const double Tr31 = TRACE_STRUCT.Tr31;
   const double Tr23 = TRACE_STRUCT.Tr23;
   const double Tr2U141 = TRACE_STRUCT.Tr2U141;
   const double Tr2U144 = TRACE_STRUCT.Tr2U144;
   const double Tr34 = TRACE_STRUCT.Tr34;


   Eigen::Matrix<double,3,3> beta_mu2;

   beta_mu2 = twoLoop*(-12*traceconjTYuTpTYu*(Yu*Yu.adjoint()) - 12*
      tracemq2AdjYuYu*(Yu*Yu.adjoint()) - 12*tracemu2YuAdjYu*(Yu*Yu.adjoint())
      - 24*mHu2*traceYuAdjYu*(Yu*Yu.adjoint()) - 4*mHd2*AbsSqr(Lambdax)*(Yu*
      Yu.adjoint()) - 8*mHu2*AbsSqr(Lambdax)*(Yu*Yu.adjoint()) - 4*ms2*AbsSqr(
      Lambdax)*(Yu*Yu.adjoint()) - 4*AbsSqr(TLambdax)*(Yu*Yu.adjoint()) - 0.8*
      mHu2*Sqr(g1)*(Yu*Yu.adjoint()) + 12*mHu2*Sqr(g2)*(Yu*Yu.adjoint()) + 24*
      AbsSqr(MassWB)*Sqr(g2)*(Yu*Yu.adjoint()) + 8*mHu2*Sqr(gN)*Sqr(QH2p)*(Yu*
      Yu.adjoint()) + 8*mHu2*Sqr(gN)*Sqr(QQp)*(Yu*Yu.adjoint()) - 8*mHu2*Sqr(gN
      )*Sqr(Qup)*(Yu*Yu.adjoint()) - 12*traceAdjYuTYu*(Yu*(TYu).adjoint()) +
      0.8*MassB*Sqr(g1)*(Yu*(TYu).adjoint()) - 12*MassWB*Sqr(g2)*(Yu*(TYu)
      .adjoint()) - 8*MassBp*Sqr(gN)*Sqr(QH2p)*(Yu*(TYu).adjoint()) - 8*MassBp*
      Sqr(gN)*Sqr(QQp)*(Yu*(TYu).adjoint()) + 8*MassBp*Sqr(gN)*Sqr(Qup)*(Yu*(
      TYu).adjoint()) - 4*Conj(Lambdax)*TLambdax*(Yu*(TYu).adjoint()) - 12*
      traceconjTYuTpYu*(TYu*Yu.adjoint()) - 4*Conj(TLambdax)*Lambdax*(TYu*
      Yu.adjoint()) - 12*Conj(MassWB)*Sqr(g2)*(TYu*Yu.adjoint()) - 12*
      traceYuAdjYu*(TYu*(TYu).adjoint()) - 4*AbsSqr(Lambdax)*(TYu*(TYu).adjoint
      ()) - 0.8*Sqr(g1)*(TYu*(TYu).adjoint()) + 12*Sqr(g2)*(TYu*(TYu).adjoint()
      ) + 8*Sqr(gN)*Sqr(QH2p)*(TYu*(TYu).adjoint()) + 8*Sqr(gN)*Sqr(QQp)*(TYu*(
      TYu).adjoint()) - 8*Sqr(gN)*Sqr(Qup)*(TYu*(TYu).adjoint()) - 6*
      traceYuAdjYu*(mu2*Yu*Yu.adjoint()) - 2*AbsSqr(Lambdax)*(mu2*Yu*Yu.adjoint
      ()) - 0.4*Sqr(g1)*(mu2*Yu*Yu.adjoint()) + 6*Sqr(g2)*(mu2*Yu*Yu.adjoint())
      + 4*Sqr(gN)*Sqr(QH2p)*(mu2*Yu*Yu.adjoint()) + 4*Sqr(gN)*Sqr(QQp)*(mu2*Yu
      *Yu.adjoint()) - 4*Sqr(gN)*Sqr(Qup)*(mu2*Yu*Yu.adjoint()) - 12*
      traceYuAdjYu*(Yu*mq2*Yu.adjoint()) - 4*AbsSqr(Lambdax)*(Yu*mq2*Yu.adjoint
      ()) - 0.8*Sqr(g1)*(Yu*mq2*Yu.adjoint()) + 12*Sqr(g2)*(Yu*mq2*Yu.adjoint()
      ) + 8*Sqr(gN)*Sqr(QH2p)*(Yu*mq2*Yu.adjoint()) + 8*Sqr(gN)*Sqr(QQp)*(Yu*
      mq2*Yu.adjoint()) - 8*Sqr(gN)*Sqr(Qup)*(Yu*mq2*Yu.adjoint()) - 6*
      traceYuAdjYu*(Yu*Yu.adjoint()*mu2) - 2*AbsSqr(Lambdax)*(Yu*Yu.adjoint()*
      mu2) - 0.4*Sqr(g1)*(Yu*Yu.adjoint()*mu2) + 6*Sqr(g2)*(Yu*Yu.adjoint()*mu2
      ) + 4*Sqr(gN)*Sqr(QH2p)*(Yu*Yu.adjoint()*mu2) + 4*Sqr(gN)*Sqr(QQp)*(Yu*
      Yu.adjoint()*mu2) - 4*Sqr(gN)*Sqr(Qup)*(Yu*Yu.adjoint()*mu2) - 4*mHd2*(Yu
      *Yd.adjoint()*Yd*Yu.adjoint()) - 4*mHu2*(Yu*Yd.adjoint()*Yd*Yu.adjoint())
      - 4*(Yu*Yd.adjoint()*TYd*(TYu).adjoint()) - 8*mHu2*(Yu*Yu.adjoint()*Yu*
      Yu.adjoint()) - 4*(Yu*Yu.adjoint()*TYu*(TYu).adjoint()) - 4*(Yu*(TYd)
      .adjoint()*TYd*Yu.adjoint()) - 4*(Yu*(TYu).adjoint()*TYu*Yu.adjoint()) -
      4*(TYu*Yd.adjoint()*Yd*(TYu).adjoint()) - 4*(TYu*Yu.adjoint()*Yu*(TYu)
      .adjoint()) - 4*(TYu*(TYd).adjoint()*Yd*Yu.adjoint()) - 4*(TYu*(TYu)
      .adjoint()*Yu*Yu.adjoint()) - 2*(mu2*Yu*Yd.adjoint()*Yd*Yu.adjoint()) - 2
      *(mu2*Yu*Yu.adjoint()*Yu*Yu.adjoint()) - 4*(Yu*mq2*Yd.adjoint()*Yd*
      Yu.adjoint()) - 4*(Yu*mq2*Yu.adjoint()*Yu*Yu.adjoint()) - 4*(Yu*
      Yd.adjoint()*md2*Yd*Yu.adjoint()) - 4*(Yu*Yd.adjoint()*Yd*mq2*Yu.adjoint(
      )) - 2*(Yu*Yd.adjoint()*Yd*Yu.adjoint()*mu2) - 4*(Yu*Yu.adjoint()*mu2*Yu*
      Yu.adjoint()) - 4*(Yu*Yu.adjoint()*Yu*mq2*Yu.adjoint()) - 2*(Yu*
      Yu.adjoint()*Yu*Yu.adjoint()*mu2) + 10.666666666666666*Power(g3,4)*Tr23*
      UNITMATRIX(3) - 4.131182235954578*g1*gN*Qup*Tr2U114*UNITMATRIX(3) -
      4.131182235954578*g1*gN*Qup*Tr2U141*UNITMATRIX(3) - 4.131182235954578*g1*
      Tr31*UNITMATRIX(3) + 8*gN*Qup*Tr34*UNITMATRIX(3) + 53.333333333333336*
      Power(g3,4)*AbsSqr(MassG)*UNITMATRIX(3) + 2.1333333333333333*Tr2U111*Sqr(
      g1)*UNITMATRIX(3) + 11.377777777777778*AbsSqr(MassG)*Sqr(g1)*Sqr(g3)*
      UNITMATRIX(3) + 5.688888888888889*MassB*Conj(MassG)*Sqr(g1)*Sqr(g3)*
      UNITMATRIX(3) + 8*Tr2U144*Sqr(gN)*Sqr(Qup)*UNITMATRIX(3) +
      42.666666666666664*AbsSqr(MassG)*Sqr(g3)*Sqr(gN)*Sqr(Qup)*UNITMATRIX(3) +
      21.333333333333332*MassBp*Conj(MassG)*Sqr(g3)*Sqr(gN)*Sqr(Qup)*
      UNITMATRIX(3) + 0.017777777777777778*Conj(MassB)*Sqr(g1)*(45*(-2*MassB*(
      Yu*Yu.adjoint()) + TYu*Yu.adjoint()) + 4*(912*MassB*Sqr(g1) + 5*(16*(2*
      MassB + MassG)*Sqr(g3) - 3*(2*MassB + MassBp)*(9*Qdp + 9*QDxbarp - 9*QDxp
      + 9*Qep - 9*QH1p + 9*QH2p + 3*QHpbarp - 3*QHpp - 9*QLp + 9*QQp - 22*Qup)
      *Qup*Sqr(gN)))*UNITMATRIX(3)) + 0.5333333333333333*Conj(MassBp)*Sqr(gN)*(
      15*(Sqr(QH2p) + Sqr(QQp) - Sqr(Qup))*(2*MassBp*(Yu*Yu.adjoint()) - TYu*
      Yu.adjoint()) + Qup*(-2*(MassB + 2*MassBp)*(9*Qdp + 9*QDxbarp - 9*QDxp +
      9*Qep - 9*QH1p + 9*QH2p + 3*QHpbarp - 3*QHpp - 9*QLp + 9*QQp - 22*Qup)*
      Sqr(g1) + 5*Qup*(8*(2*MassBp + MassG)*Sqr(g3) + 9*MassBp*Sqr(gN)*(9*Sqr(
      Qdp) + 9*Sqr(QDxbarp) + 9*Sqr(QDxp) + 3*Sqr(Qep) + 6*Sqr(QH1p) + 6*Sqr(
      QH2p) + 2*Sqr(QHpbarp) + 2*Sqr(QHpp) + 6*Sqr(QLp) + 18*Sqr(QQp) + 3*Sqr(
      QSp) + 11*Sqr(Qup))))*UNITMATRIX(3)));


   return beta_mu2;
}

} // namespace flexiblesusy
