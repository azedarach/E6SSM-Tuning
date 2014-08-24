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

// File generated at Sun 24 Aug 2014 16:10:14

#include "lowE6SSM_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of md2.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> lowE6SSM_soft_parameters::calc_beta_md2_one_loop(const Soft_traces& soft_traces) const
{
   const auto Qdp = INPUT(Qdp);
   const double Tr11 = TRACE_STRUCT.Tr11;
   const double Tr14 = TRACE_STRUCT.Tr14;


   Eigen::Matrix<double,3,3> beta_md2;

   beta_md2 = oneOver16PiSqr*(4*mHd2*(Yd*Yd.adjoint()) + 4*(TYd*(TYd)
      .adjoint()) + 2*(md2*Yd*Yd.adjoint()) + 4*(Yd*mq2*Yd.adjoint()) + 2*(Yd*
      Yd.adjoint()*md2) + 0.5163977794943222*g1*Tr11*UNITMATRIX(3) + 2*gN*Qdp*
      Tr14*UNITMATRIX(3) - 0.5333333333333333*AbsSqr(MassB)*Sqr(g1)*UNITMATRIX(
      3) - 10.666666666666666*AbsSqr(MassG)*Sqr(g3)*UNITMATRIX(3) - 8*AbsSqr(
      MassBp)*Sqr(gN)*Sqr(Qdp)*UNITMATRIX(3));


   return beta_md2;
}

/**
 * Calculates the two-loop beta function of md2.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> lowE6SSM_soft_parameters::calc_beta_md2_two_loop(const Soft_traces& soft_traces) const
{
   const auto Qdp = INPUT(Qdp);
   const auto QH1p = INPUT(QH1p);
   const auto QQp = INPUT(QQp);
   const auto QDxbarp = INPUT(QDxbarp);
   const auto QDxp = INPUT(QDxp);
   const auto Qep = INPUT(Qep);
   const auto QH2p = INPUT(QH2p);
   const auto QHpbarp = INPUT(QHpbarp);
   const auto QHpp = INPUT(QHpp);
   const auto QLp = INPUT(QLp);
   const auto Qup = INPUT(Qup);
   const auto QSp = INPUT(QSp);
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceconjTYdTpTYd = TRACE_STRUCT.traceconjTYdTpTYd;
   const double traceconjTYeTpTYe = TRACE_STRUCT.traceconjTYeTpTYe;
   const double tracemd2YdAdjYd = TRACE_STRUCT.tracemd2YdAdjYd;
   const double traceme2YeAdjYe = TRACE_STRUCT.traceme2YeAdjYe;
   const double traceml2AdjYeYe = TRACE_STRUCT.traceml2AdjYeYe;
   const double tracemq2AdjYdYd = TRACE_STRUCT.tracemq2AdjYdYd;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceconjTYdTpYd = TRACE_STRUCT.traceconjTYdTpYd;
   const double traceconjTYeTpYe = TRACE_STRUCT.traceconjTYeTpYe;
   const double Tr2U111 = TRACE_STRUCT.Tr2U111;
   const double Tr2U114 = TRACE_STRUCT.Tr2U114;
   const double Tr31 = TRACE_STRUCT.Tr31;
   const double Tr23 = TRACE_STRUCT.Tr23;
   const double Tr2U141 = TRACE_STRUCT.Tr2U141;
   const double Tr2U144 = TRACE_STRUCT.Tr2U144;
   const double Tr34 = TRACE_STRUCT.Tr34;


   Eigen::Matrix<double,3,3> beta_md2;

   beta_md2 = twoLoop*(-12*traceconjTYdTpTYd*(Yd*Yd.adjoint()) - 4*
      traceconjTYeTpTYe*(Yd*Yd.adjoint()) - 12*tracemd2YdAdjYd*(Yd*Yd.adjoint()
      ) - 4*traceme2YeAdjYe*(Yd*Yd.adjoint()) - 4*traceml2AdjYeYe*(Yd*
      Yd.adjoint()) - 12*tracemq2AdjYdYd*(Yd*Yd.adjoint()) - 24*mHd2*
      traceYdAdjYd*(Yd*Yd.adjoint()) - 8*mHd2*traceYeAdjYe*(Yd*Yd.adjoint()) -
      8*mHd2*AbsSqr(Lambdax)*(Yd*Yd.adjoint()) - 4*mHu2*AbsSqr(Lambdax)*(Yd*
      Yd.adjoint()) - 4*ms2*AbsSqr(Lambdax)*(Yd*Yd.adjoint()) - 4*AbsSqr(
      TLambdax)*(Yd*Yd.adjoint()) + 0.8*mHd2*Sqr(g1)*(Yd*Yd.adjoint()) + 12*
      mHd2*Sqr(g2)*(Yd*Yd.adjoint()) + 24*AbsSqr(MassWB)*Sqr(g2)*(Yd*Yd.adjoint
      ()) - 8*mHd2*Sqr(gN)*Sqr(Qdp)*(Yd*Yd.adjoint()) + 8*mHd2*Sqr(gN)*Sqr(QH1p
      )*(Yd*Yd.adjoint()) + 8*mHd2*Sqr(gN)*Sqr(QQp)*(Yd*Yd.adjoint()) - 12*
      traceAdjYdTYd*(Yd*(TYd).adjoint()) - 4*traceAdjYeTYe*(Yd*(TYd).adjoint())
      - 0.8*MassB*Sqr(g1)*(Yd*(TYd).adjoint()) - 12*MassWB*Sqr(g2)*(Yd*(TYd)
      .adjoint()) + 8*MassBp*Sqr(gN)*Sqr(Qdp)*(Yd*(TYd).adjoint()) - 8*MassBp*
      Sqr(gN)*Sqr(QH1p)*(Yd*(TYd).adjoint()) - 8*MassBp*Sqr(gN)*Sqr(QQp)*(Yd*(
      TYd).adjoint()) - 4*Conj(Lambdax)*TLambdax*(Yd*(TYd).adjoint()) - 12*
      traceconjTYdTpYd*(TYd*Yd.adjoint()) - 4*traceconjTYeTpYe*(TYd*Yd.adjoint(
      )) - 4*Conj(TLambdax)*Lambdax*(TYd*Yd.adjoint()) - 12*Conj(MassWB)*Sqr(g2
      )*(TYd*Yd.adjoint()) - 12*traceYdAdjYd*(TYd*(TYd).adjoint()) - 4*
      traceYeAdjYe*(TYd*(TYd).adjoint()) - 4*AbsSqr(Lambdax)*(TYd*(TYd).adjoint
      ()) + 0.8*Sqr(g1)*(TYd*(TYd).adjoint()) + 12*Sqr(g2)*(TYd*(TYd).adjoint()
      ) - 8*Sqr(gN)*Sqr(Qdp)*(TYd*(TYd).adjoint()) + 8*Sqr(gN)*Sqr(QH1p)*(TYd*(
      TYd).adjoint()) + 8*Sqr(gN)*Sqr(QQp)*(TYd*(TYd).adjoint()) - 6*
      traceYdAdjYd*(md2*Yd*Yd.adjoint()) - 2*traceYeAdjYe*(md2*Yd*Yd.adjoint())
      - 2*AbsSqr(Lambdax)*(md2*Yd*Yd.adjoint()) + 0.4*Sqr(g1)*(md2*Yd*
      Yd.adjoint()) + 6*Sqr(g2)*(md2*Yd*Yd.adjoint()) - 4*Sqr(gN)*Sqr(Qdp)*(md2
      *Yd*Yd.adjoint()) + 4*Sqr(gN)*Sqr(QH1p)*(md2*Yd*Yd.adjoint()) + 4*Sqr(gN)
      *Sqr(QQp)*(md2*Yd*Yd.adjoint()) - 12*traceYdAdjYd*(Yd*mq2*Yd.adjoint()) -
      4*traceYeAdjYe*(Yd*mq2*Yd.adjoint()) - 4*AbsSqr(Lambdax)*(Yd*mq2*
      Yd.adjoint()) + 0.8*Sqr(g1)*(Yd*mq2*Yd.adjoint()) + 12*Sqr(g2)*(Yd*mq2*
      Yd.adjoint()) - 8*Sqr(gN)*Sqr(Qdp)*(Yd*mq2*Yd.adjoint()) + 8*Sqr(gN)*Sqr(
      QH1p)*(Yd*mq2*Yd.adjoint()) + 8*Sqr(gN)*Sqr(QQp)*(Yd*mq2*Yd.adjoint()) -
      6*traceYdAdjYd*(Yd*Yd.adjoint()*md2) - 2*traceYeAdjYe*(Yd*Yd.adjoint()*
      md2) - 2*AbsSqr(Lambdax)*(Yd*Yd.adjoint()*md2) + 0.4*Sqr(g1)*(Yd*
      Yd.adjoint()*md2) + 6*Sqr(g2)*(Yd*Yd.adjoint()*md2) - 4*Sqr(gN)*Sqr(Qdp)*
      (Yd*Yd.adjoint()*md2) + 4*Sqr(gN)*Sqr(QH1p)*(Yd*Yd.adjoint()*md2) + 4*Sqr
      (gN)*Sqr(QQp)*(Yd*Yd.adjoint()*md2) - 8*mHd2*(Yd*Yd.adjoint()*Yd*
      Yd.adjoint()) - 4*(Yd*Yd.adjoint()*TYd*(TYd).adjoint()) - 4*mHd2*(Yd*
      Yu.adjoint()*Yu*Yd.adjoint()) - 4*mHu2*(Yd*Yu.adjoint()*Yu*Yd.adjoint())
      - 4*(Yd*Yu.adjoint()*TYu*(TYd).adjoint()) - 4*(Yd*(TYd).adjoint()*TYd*
      Yd.adjoint()) - 4*(Yd*(TYu).adjoint()*TYu*Yd.adjoint()) - 4*(TYd*
      Yd.adjoint()*Yd*(TYd).adjoint()) - 4*(TYd*Yu.adjoint()*Yu*(TYd).adjoint()
      ) - 4*(TYd*(TYd).adjoint()*Yd*Yd.adjoint()) - 4*(TYd*(TYu).adjoint()*Yu*
      Yd.adjoint()) - 2*(md2*Yd*Yd.adjoint()*Yd*Yd.adjoint()) - 2*(md2*Yd*
      Yu.adjoint()*Yu*Yd.adjoint()) - 4*(Yd*mq2*Yd.adjoint()*Yd*Yd.adjoint()) -
      4*(Yd*mq2*Yu.adjoint()*Yu*Yd.adjoint()) - 4*(Yd*Yd.adjoint()*md2*Yd*
      Yd.adjoint()) - 4*(Yd*Yd.adjoint()*Yd*mq2*Yd.adjoint()) - 2*(Yd*
      Yd.adjoint()*Yd*Yd.adjoint()*md2) - 4*(Yd*Yu.adjoint()*mu2*Yu*Yd.adjoint(
      )) - 4*(Yd*Yu.adjoint()*Yu*mq2*Yd.adjoint()) - 2*(Yd*Yu.adjoint()*Yu*
      Yd.adjoint()*md2) + 10.666666666666666*Power(g3,4)*Tr23*UNITMATRIX(3) +
      2.065591117977289*g1*gN*Qdp*Tr2U114*UNITMATRIX(3) + 2.065591117977289*g1*
      gN*Qdp*Tr2U141*UNITMATRIX(3) + 2.065591117977289*g1*Tr31*UNITMATRIX(3) +
      8*gN*Qdp*Tr34*UNITMATRIX(3) + 53.333333333333336*Power(g3,4)*AbsSqr(MassG
      )*UNITMATRIX(3) + 0.5333333333333333*Tr2U111*Sqr(g1)*UNITMATRIX(3) +
      2.8444444444444446*AbsSqr(MassG)*Sqr(g1)*Sqr(g3)*UNITMATRIX(3) +
      1.4222222222222223*MassB*Conj(MassG)*Sqr(g1)*Sqr(g3)*UNITMATRIX(3) + 8*
      Tr2U144*Sqr(gN)*Sqr(Qdp)*UNITMATRIX(3) + 42.666666666666664*AbsSqr(MassG)
      *Sqr(g3)*Sqr(gN)*Sqr(Qdp)*UNITMATRIX(3) + 21.333333333333332*MassBp*Conj(
      MassG)*Sqr(g3)*Sqr(gN)*Sqr(Qdp)*UNITMATRIX(3) + 0.017777777777777778*Conj
      (MassB)*Sqr(g1)*(90*MassB*(Yd*Yd.adjoint()) - 45*(TYd*Yd.adjoint()) + 2*(
      438*MassB*Sqr(g1) + 5*(8*(2*MassB + MassG)*Sqr(g3) + 3*(2*MassB + MassBp)
      *Qdp*(11*Qdp + 3*(3*QDxbarp - 3*QDxp + 3*Qep - 3*QH1p + 3*QH2p + QHpbarp
      - QHpp - 3*QLp + 3*QQp - 6*Qup))*Sqr(gN)))*UNITMATRIX(3)) +
      0.5333333333333333*Conj(MassBp)*Sqr(gN)*(-15*(Sqr(Qdp) - Sqr(QH1p) - Sqr(
      QQp))*(2*MassBp*(Yd*Yd.adjoint()) - TYd*Yd.adjoint()) + Qdp*((MassB + 2*
      MassBp)*(11*Qdp + 3*(3*QDxbarp - 3*QDxp + 3*Qep - 3*QH1p + 3*QH2p +
      QHpbarp - QHpp - 3*QLp + 3*QQp - 6*Qup))*Sqr(g1) + 5*Qdp*(8*(2*MassBp +
      MassG)*Sqr(g3) + 9*MassBp*Sqr(gN)*(11*Sqr(Qdp) + 9*Sqr(QDxbarp) + 9*Sqr(
      QDxp) + 3*Sqr(Qep) + 6*Sqr(QH1p) + 6*Sqr(QH2p) + 2*Sqr(QHpbarp) + 2*Sqr(
      QHpp) + 6*Sqr(QLp) + 18*Sqr(QQp) + 3*Sqr(QSp) + 9*Sqr(Qup))))*UNITMATRIX(
      3)));


   return beta_md2;
}

} // namespace flexiblesusy
