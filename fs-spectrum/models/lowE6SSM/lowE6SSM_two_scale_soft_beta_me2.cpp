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
 * Calculates the one-loop beta function of me2.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> lowE6SSM_soft_parameters::calc_beta_me2_one_loop(const Soft_traces& soft_traces) const
{
   const auto Qep = INPUT(Qep);
   const double Tr11 = TRACE_STRUCT.Tr11;
   const double Tr14 = TRACE_STRUCT.Tr14;


   Eigen::Matrix<double,3,3> beta_me2;

   beta_me2 = oneOver16PiSqr*(4*mHd2*(Ye*Ye.adjoint()) + 4*(TYe*(TYe)
      .adjoint()) + 2*(me2*Ye*Ye.adjoint()) + 4*(Ye*ml2*Ye.adjoint()) + 2*(Ye*
      Ye.adjoint()*me2) + 1.5491933384829668*g1*Tr11*UNITMATRIX(3) + 2*gN*Qep*
      Tr14*UNITMATRIX(3) - 4.8*AbsSqr(MassB)*Sqr(g1)*UNITMATRIX(3) - 8*AbsSqr(
      MassBp)*Sqr(gN)*Sqr(Qep)*UNITMATRIX(3));


   return beta_me2;
}

/**
 * Calculates the two-loop beta function of me2.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> lowE6SSM_soft_parameters::calc_beta_me2_two_loop(const Soft_traces& soft_traces) const
{
   const auto Qep = INPUT(Qep);
   const auto Qdp = INPUT(Qdp);
   const auto QDxbarp = INPUT(QDxbarp);
   const auto QDxp = INPUT(QDxp);
   const auto QH1p = INPUT(QH1p);
   const auto QH2p = INPUT(QH2p);
   const auto QHpbarp = INPUT(QHpbarp);
   const auto QHpp = INPUT(QHpp);
   const auto QLp = INPUT(QLp);
   const auto QQp = INPUT(QQp);
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
   const double Tr2U141 = TRACE_STRUCT.Tr2U141;
   const double Tr2U144 = TRACE_STRUCT.Tr2U144;
   const double Tr34 = TRACE_STRUCT.Tr34;


   Eigen::Matrix<double,3,3> beta_me2;

   beta_me2 = 0.08*twoLoop*(6*Conj(MassB)*Sqr(g1)*(5*(-2*MassB*(Ye*
      Ye.adjoint()) + TYe*Ye.adjoint()) + 2*(162*MassB*Sqr(g1) + 5*(2*MassB +
      MassBp)*Qep*(3*Qdp + 3*QDxbarp - 3*QDxp + 5*Qep - 3*QH1p + 3*QH2p +
      QHpbarp - QHpp - 3*QLp + 3*QQp - 6*Qup)*Sqr(gN))*UNITMATRIX(3)) + 5*(-30*
      traceconjTYdTpTYd*(Ye*Ye.adjoint()) - 10*traceconjTYeTpTYe*(Ye*Ye.adjoint
      ()) - 30*tracemd2YdAdjYd*(Ye*Ye.adjoint()) - 10*traceme2YeAdjYe*(Ye*
      Ye.adjoint()) - 10*traceml2AdjYeYe*(Ye*Ye.adjoint()) - 30*tracemq2AdjYdYd
      *(Ye*Ye.adjoint()) - 60*mHd2*traceYdAdjYd*(Ye*Ye.adjoint()) - 20*mHd2*
      traceYeAdjYe*(Ye*Ye.adjoint()) - 20*mHd2*AbsSqr(Lambdax)*(Ye*Ye.adjoint()
      ) - 10*mHu2*AbsSqr(Lambdax)*(Ye*Ye.adjoint()) - 10*ms2*AbsSqr(Lambdax)*(
      Ye*Ye.adjoint()) - 10*AbsSqr(TLambdax)*(Ye*Ye.adjoint()) - 6*mHd2*Sqr(g1)
      *(Ye*Ye.adjoint()) + 30*mHd2*Sqr(g2)*(Ye*Ye.adjoint()) + 60*AbsSqr(MassWB
      )*Sqr(g2)*(Ye*Ye.adjoint()) - 20*mHd2*Sqr(gN)*Sqr(Qep)*(Ye*Ye.adjoint())
      + 20*mHd2*Sqr(gN)*Sqr(QH1p)*(Ye*Ye.adjoint()) + 20*mHd2*Sqr(gN)*Sqr(QLp)*
      (Ye*Ye.adjoint()) - 30*traceAdjYdTYd*(Ye*(TYe).adjoint()) - 10*
      traceAdjYeTYe*(Ye*(TYe).adjoint()) + 6*MassB*Sqr(g1)*(Ye*(TYe).adjoint())
      - 30*MassWB*Sqr(g2)*(Ye*(TYe).adjoint()) + 20*MassBp*Sqr(gN)*Sqr(Qep)*(
      Ye*(TYe).adjoint()) - 20*MassBp*Sqr(gN)*Sqr(QH1p)*(Ye*(TYe).adjoint()) -
      20*MassBp*Sqr(gN)*Sqr(QLp)*(Ye*(TYe).adjoint()) - 10*Conj(Lambdax)*
      TLambdax*(Ye*(TYe).adjoint()) - 30*traceconjTYdTpYd*(TYe*Ye.adjoint()) -
      10*traceconjTYeTpYe*(TYe*Ye.adjoint()) - 10*Conj(TLambdax)*Lambdax*(TYe*
      Ye.adjoint()) - 30*Conj(MassWB)*Sqr(g2)*(TYe*Ye.adjoint()) - 30*
      traceYdAdjYd*(TYe*(TYe).adjoint()) - 10*traceYeAdjYe*(TYe*(TYe).adjoint()
      ) - 10*AbsSqr(Lambdax)*(TYe*(TYe).adjoint()) - 6*Sqr(g1)*(TYe*(TYe)
      .adjoint()) + 30*Sqr(g2)*(TYe*(TYe).adjoint()) - 20*Sqr(gN)*Sqr(Qep)*(TYe
      *(TYe).adjoint()) + 20*Sqr(gN)*Sqr(QH1p)*(TYe*(TYe).adjoint()) + 20*Sqr(
      gN)*Sqr(QLp)*(TYe*(TYe).adjoint()) - 15*traceYdAdjYd*(me2*Ye*Ye.adjoint()
      ) - 5*traceYeAdjYe*(me2*Ye*Ye.adjoint()) - 5*AbsSqr(Lambdax)*(me2*Ye*
      Ye.adjoint()) - 3*Sqr(g1)*(me2*Ye*Ye.adjoint()) + 15*Sqr(g2)*(me2*Ye*
      Ye.adjoint()) - 10*Sqr(gN)*Sqr(Qep)*(me2*Ye*Ye.adjoint()) + 10*Sqr(gN)*
      Sqr(QH1p)*(me2*Ye*Ye.adjoint()) + 10*Sqr(gN)*Sqr(QLp)*(me2*Ye*Ye.adjoint(
      )) - 30*traceYdAdjYd*(Ye*ml2*Ye.adjoint()) - 10*traceYeAdjYe*(Ye*ml2*
      Ye.adjoint()) - 10*AbsSqr(Lambdax)*(Ye*ml2*Ye.adjoint()) - 6*Sqr(g1)*(Ye*
      ml2*Ye.adjoint()) + 30*Sqr(g2)*(Ye*ml2*Ye.adjoint()) - 20*Sqr(gN)*Sqr(Qep
      )*(Ye*ml2*Ye.adjoint()) + 20*Sqr(gN)*Sqr(QH1p)*(Ye*ml2*Ye.adjoint()) + 20
      *Sqr(gN)*Sqr(QLp)*(Ye*ml2*Ye.adjoint()) - 15*traceYdAdjYd*(Ye*Ye.adjoint(
      )*me2) - 5*traceYeAdjYe*(Ye*Ye.adjoint()*me2) - 5*AbsSqr(Lambdax)*(Ye*
      Ye.adjoint()*me2) - 3*Sqr(g1)*(Ye*Ye.adjoint()*me2) + 15*Sqr(g2)*(Ye*
      Ye.adjoint()*me2) - 10*Sqr(gN)*Sqr(Qep)*(Ye*Ye.adjoint()*me2) + 10*Sqr(gN
      )*Sqr(QH1p)*(Ye*Ye.adjoint()*me2) + 10*Sqr(gN)*Sqr(QLp)*(Ye*Ye.adjoint()*
      me2) - 20*mHd2*(Ye*Ye.adjoint()*Ye*Ye.adjoint()) - 10*(Ye*Ye.adjoint()*
      TYe*(TYe).adjoint()) - 10*(Ye*(TYe).adjoint()*TYe*Ye.adjoint()) - 10*(TYe
      *Ye.adjoint()*Ye*(TYe).adjoint()) - 10*(TYe*(TYe).adjoint()*Ye*Ye.adjoint
      ()) - 5*(me2*Ye*Ye.adjoint()*Ye*Ye.adjoint()) - 10*(Ye*ml2*Ye.adjoint()*
      Ye*Ye.adjoint()) - 10*(Ye*Ye.adjoint()*me2*Ye*Ye.adjoint()) - 10*(Ye*
      Ye.adjoint()*Ye*ml2*Ye.adjoint()) - 5*(Ye*Ye.adjoint()*Ye*Ye.adjoint()*
      me2) + 4*(3.872983346207417*g1*(gN*Qep*(Tr2U114 + Tr2U141) + Tr31) + 5*gN
      *Qep*(gN*Qep*Tr2U144 + Tr34) + 3*Tr2U111*Sqr(g1))*UNITMATRIX(3) + 4*Conj(
      MassBp)*Sqr(gN)*(-5*(Sqr(Qep) - Sqr(QH1p) - Sqr(QLp))*(2*MassBp*(Ye*
      Ye.adjoint()) - TYe*Ye.adjoint()) + 3*Qep*((MassB + 2*MassBp)*(3*Qdp + 3*
      QDxbarp - 3*QDxp + 5*Qep - 3*QH1p + 3*QH2p + QHpbarp - QHpp - 3*QLp + 3*
      QQp - 6*Qup)*Sqr(g1) + 5*MassBp*Qep*Sqr(gN)*(9*Sqr(Qdp) + 9*Sqr(QDxbarp)
      + 9*Sqr(QDxp) + 5*Sqr(Qep) + 6*Sqr(QH1p) + 6*Sqr(QH2p) + 2*Sqr(QHpbarp) +
      2*Sqr(QHpp) + 6*Sqr(QLp) + 18*Sqr(QQp) + 3*Sqr(QSp) + 9*Sqr(Qup)))*
      UNITMATRIX(3))));


   return beta_me2;
}

} // namespace flexiblesusy
