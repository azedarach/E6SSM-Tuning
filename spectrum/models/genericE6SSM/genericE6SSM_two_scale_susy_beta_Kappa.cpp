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

// File generated at Sun 15 Jun 2014 19:15:44

#include "genericE6SSM_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of Kappa.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> genericE6SSM_susy_parameters::calc_beta_Kappa_one_loop(const Susy_traces& susy_traces) const
{
   const auto QDxbarp = INPUT(QDxbarp);
   const auto QDxp = INPUT(QDxp);
   const auto QSp = INPUT(QSp);
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;


   Eigen::Matrix<double,3,3> beta_Kappa;

   beta_Kappa = oneOver16PiSqr*(Kappa*(3*traceKappaAdjKappa + 2*
      traceLambda12AdjLambda12 + 2*AbsSqr(Lambdax) - 0.26666666666666666*Sqr(g1
      ) - 5.333333333333333*Sqr(g3) - 2*Sqr(gN)*Sqr(QDxbarp) - 2*Sqr(gN)*Sqr(
      QDxp) - 2*Sqr(gN)*Sqr(QSp)) + 2*(Kappa*(Kappa).adjoint()*Kappa));


   return beta_Kappa;
}

/**
 * Calculates the two-loop beta function of Kappa.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> genericE6SSM_susy_parameters::calc_beta_Kappa_two_loop(const Susy_traces& susy_traces) const
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
   const double traceKappaAdjKappaKappaAdjKappa =
      TRACE_STRUCT.traceKappaAdjKappaKappaAdjKappa;
   const double traceLambda12AdjLambda12Lambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12Lambda12AdjLambda12;


   Eigen::Matrix<double,3,3> beta_Kappa;

   beta_Kappa = twoLoop*(Kappa*(2.5955555555555554*Power(g1,4) +
      14.222222222222221*Power(g3,4) + 22*Power(gN,4)*Power(QDxbarp,4) + 22*
      Power(gN,4)*Power(QDxp,4) + 10*Power(gN,4)*Power(QSp,4) - 6*
      traceKappaAdjKappaKappaAdjKappa - 4*
      traceLambda12AdjLambda12Lambda12AdjLambda12 + 1.2*
      traceLambda12AdjLambda12*Sqr(g1) + 6*traceLambda12AdjLambda12*Sqr(g2) +
      1.4222222222222223*Sqr(g1)*Sqr(g3) + 2.4*Qdp*QDxbarp*Sqr(g1)*Sqr(gN) -
      2.4*Qdp*QDxp*Sqr(g1)*Sqr(gN) - 4.8*QDxbarp*QDxp*Sqr(g1)*Sqr(gN) + 2.4*
      QDxbarp*Qep*Sqr(g1)*Sqr(gN) - 2.4*QDxp*Qep*Sqr(g1)*Sqr(gN) - 2.4*QDxbarp*
      QH1p*Sqr(g1)*Sqr(gN) + 2.4*QDxp*QH1p*Sqr(g1)*Sqr(gN) + 2.4*QDxbarp*QH2p*
      Sqr(g1)*Sqr(gN) - 2.4*QDxp*QH2p*Sqr(g1)*Sqr(gN) + 0.8*QDxbarp*QHpbarp*Sqr
      (g1)*Sqr(gN) - 0.8*QDxp*QHpbarp*Sqr(g1)*Sqr(gN) - 0.8*QDxbarp*QHpp*Sqr(g1
      )*Sqr(gN) + 0.8*QDxp*QHpp*Sqr(g1)*Sqr(gN) - 2.4*QDxbarp*QLp*Sqr(g1)*Sqr(
      gN) + 2.4*QDxp*QLp*Sqr(g1)*Sqr(gN) + 2.4*QDxbarp*QQp*Sqr(g1)*Sqr(gN) -
      2.4*QDxp*QQp*Sqr(g1)*Sqr(gN) - 4.8*QDxbarp*Qup*Sqr(g1)*Sqr(gN) + 4.8*QDxp
      *Qup*Sqr(g1)*Sqr(gN) + 2.933333333333333*Sqr(g1)*Sqr(gN)*Sqr(QDxbarp) +
      10.666666666666666*Sqr(g3)*Sqr(gN)*Sqr(QDxbarp) + 18*Power(gN,4)*Sqr(Qdp)
      *Sqr(QDxbarp) + 2.933333333333333*Sqr(g1)*Sqr(gN)*Sqr(QDxp) +
      10.666666666666666*Sqr(g3)*Sqr(gN)*Sqr(QDxp) + 18*Power(gN,4)*Sqr(Qdp)*
      Sqr(QDxp) + 36*Power(gN,4)*Sqr(QDxbarp)*Sqr(QDxp) + 6*Power(gN,4)*Sqr(
      QDxbarp)*Sqr(Qep) + 6*Power(gN,4)*Sqr(QDxp)*Sqr(Qep) + 4*
      traceLambda12AdjLambda12*Sqr(gN)*Sqr(QH1p) + 12*Power(gN,4)*Sqr(QDxbarp)*
      Sqr(QH1p) + 12*Power(gN,4)*Sqr(QDxp)*Sqr(QH1p) + 4*
      traceLambda12AdjLambda12*Sqr(gN)*Sqr(QH2p) + 12*Power(gN,4)*Sqr(QDxbarp)*
      Sqr(QH2p) + 12*Power(gN,4)*Sqr(QDxp)*Sqr(QH2p) + 4*Power(gN,4)*Sqr(
      QDxbarp)*Sqr(QHpbarp) + 4*Power(gN,4)*Sqr(QDxp)*Sqr(QHpbarp) + 4*Power(gN
      ,4)*Sqr(QDxbarp)*Sqr(QHpp) + 4*Power(gN,4)*Sqr(QDxp)*Sqr(QHpp) + 12*Power
      (gN,4)*Sqr(QDxbarp)*Sqr(QLp) + 12*Power(gN,4)*Sqr(QDxp)*Sqr(QLp) + 36*
      Power(gN,4)*Sqr(QDxbarp)*Sqr(QQp) + 36*Power(gN,4)*Sqr(QDxp)*Sqr(QQp) +
      0.4*traceKappaAdjKappa*(2*Sqr(g1) + 5*(8*Sqr(g3) + 3*Sqr(gN)*(Sqr(QDxbarp
      ) + Sqr(QDxp) - Sqr(QSp)))) - 4*traceLambda12AdjLambda12*Sqr(gN)*Sqr(QSp)
      + 18*Power(gN,4)*Sqr(Qdp)*Sqr(QSp) + 24*Power(gN,4)*Sqr(QDxbarp)*Sqr(QSp
      ) + 24*Power(gN,4)*Sqr(QDxp)*Sqr(QSp) + 6*Power(gN,4)*Sqr(Qep)*Sqr(QSp) +
      12*Power(gN,4)*Sqr(QH1p)*Sqr(QSp) + 12*Power(gN,4)*Sqr(QH2p)*Sqr(QSp) +
      4*Power(gN,4)*Sqr(QHpbarp)*Sqr(QSp) + 4*Power(gN,4)*Sqr(QHpp)*Sqr(QSp) +
      12*Power(gN,4)*Sqr(QLp)*Sqr(QSp) + 36*Power(gN,4)*Sqr(QQp)*Sqr(QSp) +
      Conj(Lambdax)*(-6*traceYdAdjYd*Lambdax - 2*traceYeAdjYe*Lambdax - 6*
      traceYuAdjYu*Lambdax + 1.2*Lambdax*Sqr(g1) + 6*Lambdax*Sqr(g2) + 4*
      Lambdax*Sqr(gN)*Sqr(QH1p) + 4*Lambdax*Sqr(gN)*Sqr(QH2p) - 4*Lambdax*Sqr(
      gN)*Sqr(QSp)) + 18*Power(gN,4)*Sqr(QDxbarp)*Sqr(Qup) + 18*Power(gN,4)*Sqr
      (QDxp)*Sqr(Qup) + 18*Power(gN,4)*Sqr(QSp)*Sqr(Qup) - 4*Sqr(Conj(Lambdax))
      *Sqr(Lambdax)) + (-6*traceKappaAdjKappa - 4*traceLambda12AdjLambda12 - 4*
      AbsSqr(Lambdax) + 4*Sqr(gN)*Sqr(QSp))*(Kappa*(Kappa).adjoint()*Kappa) - 2
      *(Kappa*(Kappa).adjoint()*Kappa*(Kappa).adjoint()*Kappa));


   return beta_Kappa;
}

} // namespace flexiblesusy
