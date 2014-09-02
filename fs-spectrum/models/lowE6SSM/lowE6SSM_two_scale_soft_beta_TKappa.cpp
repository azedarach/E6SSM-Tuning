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

// File generated at Sun 24 Aug 2014 16:10:04

#include "lowE6SSM_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of TKappa.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> lowE6SSM_soft_parameters::calc_beta_TKappa_one_loop(const Soft_traces& soft_traces) const
{
   const auto QDxbarp = INPUT(QDxbarp);
   const auto QDxp = INPUT(QDxp);
   const auto QSp = INPUT(QSp);
   const double traceAdjKappaTKappa = TRACE_STRUCT.traceAdjKappaTKappa;
   const double traceAdjLambda12TLambda12 =
      TRACE_STRUCT.traceAdjLambda12TLambda12;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;


   Eigen::Matrix<double,3,3> beta_TKappa;

   beta_TKappa = oneOver16PiSqr*(3*traceKappaAdjKappa*TKappa + 2*
      traceLambda12AdjLambda12*TKappa + 2*AbsSqr(Lambdax)*TKappa -
      0.26666666666666666*Sqr(g1)*TKappa - 5.333333333333333*Sqr(g3)*TKappa - 2
      *Sqr(gN)*Sqr(QDxbarp)*TKappa - 2*Sqr(gN)*Sqr(QDxp)*TKappa - 2*Sqr(gN)*Sqr
      (QSp)*TKappa + Kappa*(6*traceAdjKappaTKappa + 4*traceAdjLambda12TLambda12
      + 0.5333333333333333*MassB*Sqr(g1) + 10.666666666666666*MassG*Sqr(g3) +
      4*MassBp*Sqr(gN)*Sqr(QDxbarp) + 4*MassBp*Sqr(gN)*Sqr(QDxp) + 4*MassBp*Sqr
      (gN)*Sqr(QSp) + 4*Conj(Lambdax)*TLambdax) + 3*(Kappa*(Kappa).adjoint()*
      TKappa) + 3*(TKappa*(Kappa).adjoint()*Kappa));


   return beta_TKappa;
}

/**
 * Calculates the two-loop beta function of TKappa.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> lowE6SSM_soft_parameters::calc_beta_TKappa_two_loop(const Soft_traces& soft_traces) const
{
   const auto QDxbarp = INPUT(QDxbarp);
   const auto Qdp = INPUT(Qdp);
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


   Eigen::Matrix<double,3,3> beta_TKappa;

   beta_TKappa = twoLoop*(2.5955555555555554*Power(g1,4)*TKappa +
      14.222222222222221*Power(g3,4)*TKappa + 22*Power(gN,4)*Power(QDxbarp,4)*
      TKappa + 22*Power(gN,4)*Power(QDxp,4)*TKappa + 10*Power(gN,4)*Power(QSp,4
      )*TKappa - 6*traceKappaAdjKappaKappaAdjKappa*TKappa - 4*
      traceLambda12AdjLambda12Lambda12AdjLambda12*TKappa - 6*traceYdAdjYd*
      AbsSqr(Lambdax)*TKappa - 2*traceYeAdjYe*AbsSqr(Lambdax)*TKappa - 6*
      traceYuAdjYu*AbsSqr(Lambdax)*TKappa + 0.8*traceKappaAdjKappa*Sqr(g1)*
      TKappa + 1.2*traceLambda12AdjLambda12*Sqr(g1)*TKappa + 1.2*AbsSqr(Lambdax
      )*Sqr(g1)*TKappa + 6*traceLambda12AdjLambda12*Sqr(g2)*TKappa + 6*AbsSqr(
      Lambdax)*Sqr(g2)*TKappa + 16*traceKappaAdjKappa*Sqr(g3)*TKappa +
      1.4222222222222223*Sqr(g1)*Sqr(g3)*TKappa + 6*traceKappaAdjKappa*Sqr(gN)*
      Sqr(QDxbarp)*TKappa + 0.5333333333333333*Sqr(g1)*Sqr(gN)*Sqr(QDxbarp)*
      TKappa + 10.666666666666666*Sqr(g3)*Sqr(gN)*Sqr(QDxbarp)*TKappa + 18*
      Power(gN,4)*Sqr(Qdp)*Sqr(QDxbarp)*TKappa + 6*traceKappaAdjKappa*Sqr(gN)*
      Sqr(QDxp)*TKappa + 0.5333333333333333*Sqr(g1)*Sqr(gN)*Sqr(QDxp)*TKappa +
      10.666666666666666*Sqr(g3)*Sqr(gN)*Sqr(QDxp)*TKappa + 18*Power(gN,4)*Sqr(
      Qdp)*Sqr(QDxp)*TKappa + 36*Power(gN,4)*Sqr(QDxbarp)*Sqr(QDxp)*TKappa + 6*
      Power(gN,4)*Sqr(QDxbarp)*Sqr(Qep)*TKappa + 6*Power(gN,4)*Sqr(QDxp)*Sqr(
      Qep)*TKappa + 4*traceLambda12AdjLambda12*Sqr(gN)*Sqr(QH1p)*TKappa + 4*
      AbsSqr(Lambdax)*Sqr(gN)*Sqr(QH1p)*TKappa + 12*Power(gN,4)*Sqr(QDxbarp)*
      Sqr(QH1p)*TKappa + 12*Power(gN,4)*Sqr(QDxp)*Sqr(QH1p)*TKappa + 4*
      traceLambda12AdjLambda12*Sqr(gN)*Sqr(QH2p)*TKappa + 4*AbsSqr(Lambdax)*Sqr
      (gN)*Sqr(QH2p)*TKappa + 12*Power(gN,4)*Sqr(QDxbarp)*Sqr(QH2p)*TKappa + 12
      *Power(gN,4)*Sqr(QDxp)*Sqr(QH2p)*TKappa + 4*Power(gN,4)*Sqr(QDxbarp)*Sqr(
      QHpbarp)*TKappa + 4*Power(gN,4)*Sqr(QDxp)*Sqr(QHpbarp)*TKappa + 4*Power(
      gN,4)*Sqr(QDxbarp)*Sqr(QHpp)*TKappa + 4*Power(gN,4)*Sqr(QDxp)*Sqr(QHpp)*
      TKappa + 12*Power(gN,4)*Sqr(QDxbarp)*Sqr(QLp)*TKappa + 12*Power(gN,4)*Sqr
      (QDxp)*Sqr(QLp)*TKappa + 36*Power(gN,4)*Sqr(QDxbarp)*Sqr(QQp)*TKappa + 36
      *Power(gN,4)*Sqr(QDxp)*Sqr(QQp)*TKappa - 6*traceKappaAdjKappa*Sqr(gN)*Sqr
      (QSp)*TKappa - 4*traceLambda12AdjLambda12*Sqr(gN)*Sqr(QSp)*TKappa - 4*
      AbsSqr(Lambdax)*Sqr(gN)*Sqr(QSp)*TKappa + 18*Power(gN,4)*Sqr(Qdp)*Sqr(QSp
      )*TKappa + 24*Power(gN,4)*Sqr(QDxbarp)*Sqr(QSp)*TKappa + 24*Power(gN,4)*
      Sqr(QDxp)*Sqr(QSp)*TKappa + 6*Power(gN,4)*Sqr(Qep)*Sqr(QSp)*TKappa + 12*
      Power(gN,4)*Sqr(QH1p)*Sqr(QSp)*TKappa + 12*Power(gN,4)*Sqr(QH2p)*Sqr(QSp)
      *TKappa + 4*Power(gN,4)*Sqr(QHpbarp)*Sqr(QSp)*TKappa + 4*Power(gN,4)*Sqr(
      QHpp)*Sqr(QSp)*TKappa + 12*Power(gN,4)*Sqr(QLp)*Sqr(QSp)*TKappa + 36*
      Power(gN,4)*Sqr(QQp)*Sqr(QSp)*TKappa + 18*Power(gN,4)*Sqr(QDxbarp)*Sqr(
      Qup)*TKappa + 18*Power(gN,4)*Sqr(QDxp)*Sqr(Qup)*TKappa + 18*Power(gN,4)*
      Sqr(QSp)*Sqr(Qup)*TKappa - 4*Sqr(Conj(Lambdax))*Sqr(Lambdax)*TKappa -
      0.017777777777777778*Kappa*(584*Power(g1,4)*MassB + 3200*Power(g3,4)*
      MassG + 4950*Power(gN,4)*MassBp*Power(QDxbarp,4) + 4950*Power(gN,4)*
      MassBp*Power(QDxp,4) + 2250*Power(gN,4)*MassBp*Power(QSp,4) + 1350*
      traceKappaAdjKappaTKappaAdjKappa + 900*
      traceLambda12AdjLambda12TLambda12AdjLambda12 - 90*traceAdjKappaTKappa*Sqr
      (g1) - 135*traceAdjLambda12TLambda12*Sqr(g1) + 135*MassB*
      traceLambda12AdjLambda12*Sqr(g1) - 675*traceAdjLambda12TLambda12*Sqr(g2)
      + 675*MassWB*traceLambda12AdjLambda12*Sqr(g2) - 1800*traceAdjKappaTKappa*
      Sqr(g3) + 160*MassB*Sqr(g1)*Sqr(g3) + 160*MassG*Sqr(g1)*Sqr(g3) - 675*
      traceAdjKappaTKappa*Sqr(gN)*Sqr(QDxbarp) + 60*MassB*Sqr(g1)*Sqr(gN)*Sqr(
      QDxbarp) + 60*MassBp*Sqr(g1)*Sqr(gN)*Sqr(QDxbarp) + 1200*MassBp*Sqr(g3)*
      Sqr(gN)*Sqr(QDxbarp) + 1200*MassG*Sqr(g3)*Sqr(gN)*Sqr(QDxbarp) + 4050*
      Power(gN,4)*MassBp*Sqr(Qdp)*Sqr(QDxbarp) - 675*traceAdjKappaTKappa*Sqr(gN
      )*Sqr(QDxp) + 60*MassB*Sqr(g1)*Sqr(gN)*Sqr(QDxp) + 60*MassBp*Sqr(g1)*Sqr(
      gN)*Sqr(QDxp) + 1200*MassBp*Sqr(g3)*Sqr(gN)*Sqr(QDxp) + 1200*MassG*Sqr(g3
      )*Sqr(gN)*Sqr(QDxp) + 4050*Power(gN,4)*MassBp*Sqr(Qdp)*Sqr(QDxp) + 8100*
      Power(gN,4)*MassBp*Sqr(QDxbarp)*Sqr(QDxp) + 1350*Power(gN,4)*MassBp*Sqr(
      QDxbarp)*Sqr(Qep) + 1350*Power(gN,4)*MassBp*Sqr(QDxp)*Sqr(Qep) - 450*
      traceAdjLambda12TLambda12*Sqr(gN)*Sqr(QH1p) + 450*MassBp*
      traceLambda12AdjLambda12*Sqr(gN)*Sqr(QH1p) + 2700*Power(gN,4)*MassBp*Sqr(
      QDxbarp)*Sqr(QH1p) + 2700*Power(gN,4)*MassBp*Sqr(QDxp)*Sqr(QH1p) - 450*
      traceAdjLambda12TLambda12*Sqr(gN)*Sqr(QH2p) + 450*MassBp*
      traceLambda12AdjLambda12*Sqr(gN)*Sqr(QH2p) + 2700*Power(gN,4)*MassBp*Sqr(
      QDxbarp)*Sqr(QH2p) + 2700*Power(gN,4)*MassBp*Sqr(QDxp)*Sqr(QH2p) + 900*
      Power(gN,4)*MassBp*Sqr(QDxbarp)*Sqr(QHpbarp) + 900*Power(gN,4)*MassBp*Sqr
      (QDxp)*Sqr(QHpbarp) + 900*Power(gN,4)*MassBp*Sqr(QDxbarp)*Sqr(QHpp) + 900
      *Power(gN,4)*MassBp*Sqr(QDxp)*Sqr(QHpp) + 2700*Power(gN,4)*MassBp*Sqr(
      QDxbarp)*Sqr(QLp) + 2700*Power(gN,4)*MassBp*Sqr(QDxp)*Sqr(QLp) + 8100*
      Power(gN,4)*MassBp*Sqr(QDxbarp)*Sqr(QQp) + 8100*Power(gN,4)*MassBp*Sqr(
      QDxp)*Sqr(QQp) + 45*traceKappaAdjKappa*(2*MassB*Sqr(g1) + 5*(8*MassG*Sqr(
      g3) + 3*MassBp*Sqr(gN)*(Sqr(QDxbarp) + Sqr(QDxp) - Sqr(QSp)))) + 675*
      traceAdjKappaTKappa*Sqr(gN)*Sqr(QSp) + 450*traceAdjLambda12TLambda12*Sqr(
      gN)*Sqr(QSp) - 450*MassBp*traceLambda12AdjLambda12*Sqr(gN)*Sqr(QSp) +
      4050*Power(gN,4)*MassBp*Sqr(Qdp)*Sqr(QSp) + 5400*Power(gN,4)*MassBp*Sqr(
      QDxbarp)*Sqr(QSp) + 5400*Power(gN,4)*MassBp*Sqr(QDxp)*Sqr(QSp) + 1350*
      Power(gN,4)*MassBp*Sqr(Qep)*Sqr(QSp) + 2700*Power(gN,4)*MassBp*Sqr(QH1p)*
      Sqr(QSp) + 2700*Power(gN,4)*MassBp*Sqr(QH2p)*Sqr(QSp) + 900*Power(gN,4)*
      MassBp*Sqr(QHpbarp)*Sqr(QSp) + 900*Power(gN,4)*MassBp*Sqr(QHpp)*Sqr(QSp)
      + 2700*Power(gN,4)*MassBp*Sqr(QLp)*Sqr(QSp) + 8100*Power(gN,4)*MassBp*Sqr
      (QQp)*Sqr(QSp) + 4050*Power(gN,4)*MassBp*Sqr(QDxbarp)*Sqr(Qup) + 4050*
      Power(gN,4)*MassBp*Sqr(QDxp)*Sqr(Qup) + 4050*Power(gN,4)*MassBp*Sqr(QSp)*
      Sqr(Qup) + 900*Lambdax*Sqr(Conj(Lambdax))*TLambdax + 45*Conj(Lambdax)*(
      Lambdax*(15*traceAdjYdTYd + 5*traceAdjYeTYe + 15*traceAdjYuTYu + 3*MassB*
      Sqr(g1) + 15*MassWB*Sqr(g2) + 10*MassBp*Sqr(gN)*Sqr(QH1p) + 10*MassBp*Sqr
      (gN)*Sqr(QH2p) - 10*MassBp*Sqr(gN)*Sqr(QSp)) + (15*traceYdAdjYd + 5*
      traceYeAdjYe + 15*traceYuAdjYu - 3*Sqr(g1) - 15*Sqr(g2) - 10*Sqr(gN)*Sqr(
      QH1p) - 10*Sqr(gN)*Sqr(QH2p) + 10*Sqr(gN)*Sqr(QSp))*TLambdax)) - 4*(3*
      traceAdjKappaTKappa + 2*(traceAdjLambda12TLambda12 + MassBp*Sqr(gN)*Sqr(
      QSp)) + 2*Conj(Lambdax)*TLambdax)*(Kappa*(Kappa).adjoint()*Kappa) - 9*
      traceKappaAdjKappa*(Kappa*(Kappa).adjoint()*TKappa) - 6*
      traceLambda12AdjLambda12*(Kappa*(Kappa).adjoint()*TKappa) - 6*AbsSqr(
      Lambdax)*(Kappa*(Kappa).adjoint()*TKappa) - 2*Sqr(gN)*Sqr(QDxbarp)*(Kappa
      *(Kappa).adjoint()*TKappa) + 2*Sqr(gN)*Sqr(QDxp)*(Kappa*(Kappa).adjoint()
      *TKappa) + 6*Sqr(gN)*Sqr(QSp)*(Kappa*(Kappa).adjoint()*TKappa) - 9*
      traceKappaAdjKappa*(TKappa*(Kappa).adjoint()*Kappa) - 6*
      traceLambda12AdjLambda12*(TKappa*(Kappa).adjoint()*Kappa) - 6*AbsSqr(
      Lambdax)*(TKappa*(Kappa).adjoint()*Kappa) + 2*Sqr(gN)*Sqr(QDxbarp)*(
      TKappa*(Kappa).adjoint()*Kappa) - 2*Sqr(gN)*Sqr(QDxp)*(TKappa*(Kappa)
      .adjoint()*Kappa) + 6*Sqr(gN)*Sqr(QSp)*(TKappa*(Kappa).adjoint()*Kappa) -
      3*(Kappa*(Kappa).adjoint()*Kappa*(Kappa).adjoint()*TKappa) - 4*(Kappa*(
      Kappa).adjoint()*TKappa*(Kappa).adjoint()*Kappa) - 3*(TKappa*(Kappa)
      .adjoint()*Kappa*(Kappa).adjoint()*Kappa));


   return beta_TKappa;
}

} // namespace flexiblesusy
