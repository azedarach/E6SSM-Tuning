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

// File generated at Sun 15 Jun 2014 19:16:30

#include "genericE6SSM_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of TYe.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> genericE6SSM_soft_parameters::calc_beta_TYe_one_loop(const Soft_traces& soft_traces) const
{
   const auto Qep = INPUT(Qep);
   const auto QH1p = INPUT(QH1p);
   const auto QLp = INPUT(QLp);
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   Eigen::Matrix<double,3,3> beta_TYe;

   beta_TYe = oneOver16PiSqr*(3*traceYdAdjYd*TYe + traceYeAdjYe*TYe +
      AbsSqr(Lambdax)*TYe - 1.8*Sqr(g1)*TYe - 3*Sqr(g2)*TYe - 2*Sqr(gN)*Sqr(Qep
      )*TYe - 2*Sqr(gN)*Sqr(QH1p)*TYe - 2*Sqr(gN)*Sqr(QLp)*TYe + Ye*(6*
      traceAdjYdTYd + 2*traceAdjYeTYe + 3.6*MassB*Sqr(g1) + 6*MassWB*Sqr(g2) +
      4*MassBp*Sqr(gN)*Sqr(Qep) + 4*MassBp*Sqr(gN)*Sqr(QH1p) + 4*MassBp*Sqr(gN)
      *Sqr(QLp) + 2*Conj(Lambdax)*TLambdax) + 4*(Ye*Ye.adjoint()*TYe) + 5*(TYe*
      Ye.adjoint()*Ye));


   return beta_TYe;
}

/**
 * Calculates the two-loop beta function of TYe.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> genericE6SSM_soft_parameters::calc_beta_TYe_two_loop(const Soft_traces& soft_traces) const
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


   Eigen::Matrix<double,3,3> beta_TYe;

   beta_TYe = twoLoop*(18.9*Power(g1,4)*TYe + 16.5*Power(g2,4)*TYe + 10*
      Power(gN,4)*Power(Qep,4)*TYe + 16*Power(gN,4)*Power(QH1p,4)*TYe + 16*
      Power(gN,4)*Power(QLp,4)*TYe - 9*traceYdAdjYdYdAdjYd*TYe - 3*
      traceYdAdjYuYuAdjYd*TYe - 3*traceYeAdjYeYeAdjYe*TYe - 3*
      traceKappaAdjKappa*AbsSqr(Lambdax)*TYe - 2*traceLambda12AdjLambda12*
      AbsSqr(Lambdax)*TYe - 3*traceYuAdjYu*AbsSqr(Lambdax)*TYe - 0.4*
      traceYdAdjYd*Sqr(g1)*TYe + 1.2*traceYeAdjYe*Sqr(g1)*TYe + 1.8*Sqr(g1)*Sqr
      (g2)*TYe + 16*traceYdAdjYd*Sqr(g3)*TYe + 6*traceYdAdjYd*Sqr(gN)*Sqr(Qdp)*
      TYe + 2*traceYeAdjYe*Sqr(gN)*Sqr(Qep)*TYe + 4.8*Sqr(g1)*Sqr(gN)*Sqr(Qep)*
      TYe + 18*Power(gN,4)*Sqr(Qdp)*Sqr(Qep)*TYe + 18*Power(gN,4)*Sqr(QDxbarp)*
      Sqr(Qep)*TYe + 18*Power(gN,4)*Sqr(QDxp)*Sqr(Qep)*TYe - 6*traceYdAdjYd*Sqr
      (gN)*Sqr(QH1p)*TYe - 2*traceYeAdjYe*Sqr(gN)*Sqr(QH1p)*TYe - 2*AbsSqr(
      Lambdax)*Sqr(gN)*Sqr(QH1p)*TYe + 1.2*Sqr(g1)*Sqr(gN)*Sqr(QH1p)*TYe + 6*
      Sqr(g2)*Sqr(gN)*Sqr(QH1p)*TYe + 18*Power(gN,4)*Sqr(Qdp)*Sqr(QH1p)*TYe +
      18*Power(gN,4)*Sqr(QDxbarp)*Sqr(QH1p)*TYe + 18*Power(gN,4)*Sqr(QDxp)*Sqr(
      QH1p)*TYe + 18*Power(gN,4)*Sqr(Qep)*Sqr(QH1p)*TYe + 2*AbsSqr(Lambdax)*Sqr
      (gN)*Sqr(QH2p)*TYe + 12*Power(gN,4)*Sqr(Qep)*Sqr(QH2p)*TYe + 12*Power(gN,
      4)*Sqr(QH1p)*Sqr(QH2p)*TYe + 4*Power(gN,4)*Sqr(Qep)*Sqr(QHpbarp)*TYe + 4*
      Power(gN,4)*Sqr(QH1p)*Sqr(QHpbarp)*TYe + 4*Power(gN,4)*Sqr(Qep)*Sqr(QHpp)
      *TYe + 4*Power(gN,4)*Sqr(QH1p)*Sqr(QHpp)*TYe + 2*traceYeAdjYe*Sqr(gN)*Sqr
      (QLp)*TYe + 1.2*Sqr(g1)*Sqr(gN)*Sqr(QLp)*TYe + 6*Sqr(g2)*Sqr(gN)*Sqr(QLp)
      *TYe + 18*Power(gN,4)*Sqr(Qdp)*Sqr(QLp)*TYe + 18*Power(gN,4)*Sqr(QDxbarp)
      *Sqr(QLp)*TYe + 18*Power(gN,4)*Sqr(QDxp)*Sqr(QLp)*TYe + 18*Power(gN,4)*
      Sqr(Qep)*Sqr(QLp)*TYe + 24*Power(gN,4)*Sqr(QH1p)*Sqr(QLp)*TYe + 12*Power(
      gN,4)*Sqr(QH2p)*Sqr(QLp)*TYe + 4*Power(gN,4)*Sqr(QHpbarp)*Sqr(QLp)*TYe +
      4*Power(gN,4)*Sqr(QHpp)*Sqr(QLp)*TYe + 6*traceYdAdjYd*Sqr(gN)*Sqr(QQp)*
      TYe + 36*Power(gN,4)*Sqr(Qep)*Sqr(QQp)*TYe + 36*Power(gN,4)*Sqr(QH1p)*Sqr
      (QQp)*TYe + 36*Power(gN,4)*Sqr(QLp)*Sqr(QQp)*TYe + 2*AbsSqr(Lambdax)*Sqr(
      gN)*Sqr(QSp)*TYe + 6*Power(gN,4)*Sqr(Qep)*Sqr(QSp)*TYe + 6*Power(gN,4)*
      Sqr(QH1p)*Sqr(QSp)*TYe + 6*Power(gN,4)*Sqr(QLp)*Sqr(QSp)*TYe + 18*Power(
      gN,4)*Sqr(Qep)*Sqr(Qup)*TYe + 18*Power(gN,4)*Sqr(QH1p)*Sqr(Qup)*TYe + 18*
      Power(gN,4)*Sqr(QLp)*Sqr(Qup)*TYe - 3*Sqr(Conj(Lambdax))*Sqr(Lambdax)*TYe
      - 0.4*Ye*(189*Power(g1,4)*MassB + 165*Power(g2,4)*MassWB + 100*Power(gN,
      4)*MassBp*Power(Qep,4) + 160*Power(gN,4)*MassBp*Power(QH1p,4) + 160*Power
      (gN,4)*MassBp*Power(QLp,4) + 90*traceYdAdjYdTYdAdjYd + 15*
      traceYdAdjYuTYuAdjYd + 30*traceYeAdjYeTYeAdjYe + 15*traceYuAdjYdTYdAdjYu
      + 2*traceAdjYdTYd*Sqr(g1) - 6*traceAdjYeTYe*Sqr(g1) + 6*MassB*
      traceYeAdjYe*Sqr(g1) + 9*MassB*Sqr(g1)*Sqr(g2) + 9*MassWB*Sqr(g1)*Sqr(g2)
      - 80*traceAdjYdTYd*Sqr(g3) - 30*traceAdjYdTYd*Sqr(gN)*Sqr(Qdp) - 10*
      traceAdjYeTYe*Sqr(gN)*Sqr(Qep) + 10*MassBp*traceYeAdjYe*Sqr(gN)*Sqr(Qep)
      + 24*MassB*Sqr(g1)*Sqr(gN)*Sqr(Qep) + 24*MassBp*Sqr(g1)*Sqr(gN)*Sqr(Qep)
      + 180*Power(gN,4)*MassBp*Sqr(Qdp)*Sqr(Qep) + 180*Power(gN,4)*MassBp*Sqr(
      QDxbarp)*Sqr(Qep) + 180*Power(gN,4)*MassBp*Sqr(QDxp)*Sqr(Qep) + 30*
      traceAdjYdTYd*Sqr(gN)*Sqr(QH1p) + 10*traceAdjYeTYe*Sqr(gN)*Sqr(QH1p) - 10
      *MassBp*traceYeAdjYe*Sqr(gN)*Sqr(QH1p) + 6*MassB*Sqr(g1)*Sqr(gN)*Sqr(QH1p
      ) + 6*MassBp*Sqr(g1)*Sqr(gN)*Sqr(QH1p) + 30*MassBp*Sqr(g2)*Sqr(gN)*Sqr(
      QH1p) + 30*MassWB*Sqr(g2)*Sqr(gN)*Sqr(QH1p) + 180*Power(gN,4)*MassBp*Sqr(
      Qdp)*Sqr(QH1p) + 180*Power(gN,4)*MassBp*Sqr(QDxbarp)*Sqr(QH1p) + 180*
      Power(gN,4)*MassBp*Sqr(QDxp)*Sqr(QH1p) + 180*Power(gN,4)*MassBp*Sqr(Qep)*
      Sqr(QH1p) + 120*Power(gN,4)*MassBp*Sqr(Qep)*Sqr(QH2p) + 120*Power(gN,4)*
      MassBp*Sqr(QH1p)*Sqr(QH2p) + 40*Power(gN,4)*MassBp*Sqr(Qep)*Sqr(QHpbarp)
      + 40*Power(gN,4)*MassBp*Sqr(QH1p)*Sqr(QHpbarp) + 40*Power(gN,4)*MassBp*
      Sqr(Qep)*Sqr(QHpp) + 40*Power(gN,4)*MassBp*Sqr(QH1p)*Sqr(QHpp) - 10*
      traceAdjYeTYe*Sqr(gN)*Sqr(QLp) + 10*MassBp*traceYeAdjYe*Sqr(gN)*Sqr(QLp)
      + 6*MassB*Sqr(g1)*Sqr(gN)*Sqr(QLp) + 6*MassBp*Sqr(g1)*Sqr(gN)*Sqr(QLp) +
      30*MassBp*Sqr(g2)*Sqr(gN)*Sqr(QLp) + 30*MassWB*Sqr(g2)*Sqr(gN)*Sqr(QLp) +
      180*Power(gN,4)*MassBp*Sqr(Qdp)*Sqr(QLp) + 180*Power(gN,4)*MassBp*Sqr(
      QDxbarp)*Sqr(QLp) + 180*Power(gN,4)*MassBp*Sqr(QDxp)*Sqr(QLp) + 180*Power
      (gN,4)*MassBp*Sqr(Qep)*Sqr(QLp) + 240*Power(gN,4)*MassBp*Sqr(QH1p)*Sqr(
      QLp) + 120*Power(gN,4)*MassBp*Sqr(QH2p)*Sqr(QLp) + 40*Power(gN,4)*MassBp*
      Sqr(QHpbarp)*Sqr(QLp) + 40*Power(gN,4)*MassBp*Sqr(QHpp)*Sqr(QLp) - 30*
      traceAdjYdTYd*Sqr(gN)*Sqr(QQp) + 360*Power(gN,4)*MassBp*Sqr(Qep)*Sqr(QQp)
      + 360*Power(gN,4)*MassBp*Sqr(QH1p)*Sqr(QQp) + 360*Power(gN,4)*MassBp*Sqr
      (QLp)*Sqr(QQp) + traceYdAdjYd*(-2*MassB*Sqr(g1) + 10*(8*MassG*Sqr(g3) + 3
      *MassBp*Sqr(gN)*(Sqr(Qdp) - Sqr(QH1p) + Sqr(QQp)))) + 60*Power(gN,4)*
      MassBp*Sqr(Qep)*Sqr(QSp) + 60*Power(gN,4)*MassBp*Sqr(QH1p)*Sqr(QSp) + 60*
      Power(gN,4)*MassBp*Sqr(QLp)*Sqr(QSp) + 180*Power(gN,4)*MassBp*Sqr(Qep)*
      Sqr(Qup) + 180*Power(gN,4)*MassBp*Sqr(QH1p)*Sqr(Qup) + 180*Power(gN,4)*
      MassBp*Sqr(QLp)*Sqr(Qup) + 30*Lambdax*Sqr(Conj(Lambdax))*TLambdax - 5*
      Conj(Lambdax)*(Lambdax*(-3*traceAdjKappaTKappa - 2*
      traceAdjLambda12TLambda12 - 3*traceAdjYuTYu + 2*MassBp*Sqr(gN)*Sqr(QH1p)
      - 2*MassBp*Sqr(gN)*Sqr(QH2p) - 2*MassBp*Sqr(gN)*Sqr(QSp)) - (3*
      traceKappaAdjKappa + 2*traceLambda12AdjLambda12 + 3*traceYuAdjYu + 2*Sqr(
      gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QSp))*TLambdax)) - 2*
      (9*traceAdjYdTYd + 3*traceAdjYeTYe + 6*MassWB*Sqr(g2) - 2*MassBp*Sqr(gN)*
      Sqr(Qep) + 6*MassBp*Sqr(gN)*Sqr(QH1p) + 2*MassBp*Sqr(gN)*Sqr(QLp) + 3*
      Conj(Lambdax)*TLambdax)*(Ye*Ye.adjoint()*Ye) - 12*traceYdAdjYd*(Ye*
      Ye.adjoint()*TYe) - 4*traceYeAdjYe*(Ye*Ye.adjoint()*TYe) - 4*AbsSqr(
      Lambdax)*(Ye*Ye.adjoint()*TYe) + 1.2*Sqr(g1)*(Ye*Ye.adjoint()*TYe) + 6*
      Sqr(g2)*(Ye*Ye.adjoint()*TYe) + 8*Sqr(gN)*Sqr(QH1p)*(Ye*Ye.adjoint()*TYe)
      - 15*traceYdAdjYd*(TYe*Ye.adjoint()*Ye) - 5*traceYeAdjYe*(TYe*Ye.adjoint
      ()*Ye) - 5*AbsSqr(Lambdax)*(TYe*Ye.adjoint()*Ye) - 1.2*Sqr(g1)*(TYe*
      Ye.adjoint()*Ye) + 12*Sqr(g2)*(TYe*Ye.adjoint()*Ye) - 6*Sqr(gN)*Sqr(Qep)*
      (TYe*Ye.adjoint()*Ye) + 10*Sqr(gN)*Sqr(QH1p)*(TYe*Ye.adjoint()*Ye) + 6*
      Sqr(gN)*Sqr(QLp)*(TYe*Ye.adjoint()*Ye) - 6*(Ye*Ye.adjoint()*Ye*Ye.adjoint
      ()*TYe) - 8*(Ye*Ye.adjoint()*TYe*Ye.adjoint()*Ye) - 6*(TYe*Ye.adjoint()*
      Ye*Ye.adjoint()*Ye));


   return beta_TYe;
}

} // namespace flexiblesusy
