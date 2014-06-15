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

// File generated at Sun 15 Jun 2014 19:15:45

#include "genericE6SSM_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of Lambda12.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,2,2> genericE6SSM_susy_parameters::calc_beta_Lambda12_one_loop(const Susy_traces& susy_traces) const
{
   const auto QH1p = INPUT(QH1p);
   const auto QH2p = INPUT(QH2p);
   const auto QSp = INPUT(QSp);
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;


   Eigen::Matrix<double,2,2> beta_Lambda12;

   beta_Lambda12 = oneOver16PiSqr*(Lambda12*(3*traceKappaAdjKappa + 2*
      traceLambda12AdjLambda12 + 2*AbsSqr(Lambdax) - 0.6*Sqr(g1) - 3*Sqr(g2) -
      2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QSp)) + 2*(
      Lambda12*(Lambda12).adjoint()*Lambda12));


   return beta_Lambda12;
}

/**
 * Calculates the two-loop beta function of Lambda12.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,2,2> genericE6SSM_susy_parameters::calc_beta_Lambda12_two_loop(const Susy_traces& susy_traces) const
{
   const auto Qdp = INPUT(Qdp);
   const auto QH1p = INPUT(QH1p);
   const auto QDxbarp = INPUT(QDxbarp);
   const auto QDxp = INPUT(QDxp);
   const auto Qep = INPUT(Qep);
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


   Eigen::Matrix<double,2,2> beta_Lambda12;

   beta_Lambda12 = twoLoop*(Lambda12*(5.94*Power(g1,4) + 16.5*Power(g2,4)
      + 16*Power(gN,4)*Power(QH1p,4) + 16*Power(gN,4)*Power(QH2p,4) + 10*Power
      (gN,4)*Power(QSp,4) - 6*traceKappaAdjKappaKappaAdjKappa - 4*
      traceLambda12AdjLambda12Lambda12AdjLambda12 + 1.2*
      traceLambda12AdjLambda12*Sqr(g1) + 6*traceLambda12AdjLambda12*Sqr(g2) +
      1.8*Sqr(g1)*Sqr(g2) - 3.6*Qdp*QH1p*Sqr(g1)*Sqr(gN) - 3.6*QDxbarp*QH1p*Sqr
      (g1)*Sqr(gN) + 3.6*QDxp*QH1p*Sqr(g1)*Sqr(gN) - 3.6*Qep*QH1p*Sqr(g1)*Sqr(
      gN) + 3.6*Qdp*QH2p*Sqr(g1)*Sqr(gN) + 3.6*QDxbarp*QH2p*Sqr(g1)*Sqr(gN) -
      3.6*QDxp*QH2p*Sqr(g1)*Sqr(gN) + 3.6*Qep*QH2p*Sqr(g1)*Sqr(gN) - 7.2*QH1p*
      QH2p*Sqr(g1)*Sqr(gN) - 1.2*QH1p*QHpbarp*Sqr(g1)*Sqr(gN) + 1.2*QH2p*
      QHpbarp*Sqr(g1)*Sqr(gN) + 1.2*QH1p*QHpp*Sqr(g1)*Sqr(gN) - 1.2*QH2p*QHpp*
      Sqr(g1)*Sqr(gN) + 3.6*QH1p*QLp*Sqr(g1)*Sqr(gN) - 3.6*QH2p*QLp*Sqr(g1)*Sqr
      (gN) - 3.6*QH1p*QQp*Sqr(g1)*Sqr(gN) + 3.6*QH2p*QQp*Sqr(g1)*Sqr(gN) + 7.2*
      QH1p*Qup*Sqr(g1)*Sqr(gN) - 7.2*QH2p*Qup*Sqr(g1)*Sqr(gN) + 4*
      traceLambda12AdjLambda12*Sqr(gN)*Sqr(QH1p) + 4.8*Sqr(g1)*Sqr(gN)*Sqr(QH1p
      ) + 6*Sqr(g2)*Sqr(gN)*Sqr(QH1p) + 18*Power(gN,4)*Sqr(Qdp)*Sqr(QH1p) + 18*
      Power(gN,4)*Sqr(QDxbarp)*Sqr(QH1p) + 18*Power(gN,4)*Sqr(QDxp)*Sqr(QH1p) +
      6*Power(gN,4)*Sqr(Qep)*Sqr(QH1p) + 4*traceLambda12AdjLambda12*Sqr(gN)*
      Sqr(QH2p) + 4.8*Sqr(g1)*Sqr(gN)*Sqr(QH2p) + 6*Sqr(g2)*Sqr(gN)*Sqr(QH2p) +
      18*Power(gN,4)*Sqr(Qdp)*Sqr(QH2p) + 18*Power(gN,4)*Sqr(QDxbarp)*Sqr(QH2p
      ) + 18*Power(gN,4)*Sqr(QDxp)*Sqr(QH2p) + 6*Power(gN,4)*Sqr(Qep)*Sqr(QH2p)
      + 24*Power(gN,4)*Sqr(QH1p)*Sqr(QH2p) + 4*Power(gN,4)*Sqr(QH1p)*Sqr(
      QHpbarp) + 4*Power(gN,4)*Sqr(QH2p)*Sqr(QHpbarp) + 4*Power(gN,4)*Sqr(QH1p)
      *Sqr(QHpp) + 4*Power(gN,4)*Sqr(QH2p)*Sqr(QHpp) + 12*Power(gN,4)*Sqr(QH1p)
      *Sqr(QLp) + 12*Power(gN,4)*Sqr(QH2p)*Sqr(QLp) + 36*Power(gN,4)*Sqr(QH1p)*
      Sqr(QQp) + 36*Power(gN,4)*Sqr(QH2p)*Sqr(QQp) + 0.4*traceKappaAdjKappa*(2*
      Sqr(g1) + 5*(8*Sqr(g3) + 3*Sqr(gN)*(Sqr(QDxbarp) + Sqr(QDxp) - Sqr(QSp)))
      ) - 4*traceLambda12AdjLambda12*Sqr(gN)*Sqr(QSp) + 18*Power(gN,4)*Sqr(Qdp)
      *Sqr(QSp) + 18*Power(gN,4)*Sqr(QDxbarp)*Sqr(QSp) + 18*Power(gN,4)*Sqr(
      QDxp)*Sqr(QSp) + 6*Power(gN,4)*Sqr(Qep)*Sqr(QSp) + 18*Power(gN,4)*Sqr(
      QH1p)*Sqr(QSp) + 18*Power(gN,4)*Sqr(QH2p)*Sqr(QSp) + 4*Power(gN,4)*Sqr(
      QHpbarp)*Sqr(QSp) + 4*Power(gN,4)*Sqr(QHpp)*Sqr(QSp) + 12*Power(gN,4)*Sqr
      (QLp)*Sqr(QSp) + 36*Power(gN,4)*Sqr(QQp)*Sqr(QSp) + Conj(Lambdax)*(-6*
      traceYdAdjYd*Lambdax - 2*traceYeAdjYe*Lambdax - 6*traceYuAdjYu*Lambdax +
      1.2*Lambdax*Sqr(g1) + 6*Lambdax*Sqr(g2) + 4*Lambdax*Sqr(gN)*Sqr(QH1p) + 4
      *Lambdax*Sqr(gN)*Sqr(QH2p) - 4*Lambdax*Sqr(gN)*Sqr(QSp)) + 18*Power(gN,4)
      *Sqr(QH1p)*Sqr(Qup) + 18*Power(gN,4)*Sqr(QH2p)*Sqr(Qup) + 18*Power(gN,4)*
      Sqr(QSp)*Sqr(Qup) - 4*Sqr(Conj(Lambdax))*Sqr(Lambdax)) + (-6*
      traceKappaAdjKappa - 4*traceLambda12AdjLambda12 - 4*AbsSqr(Lambdax) + 4*
      Sqr(gN)*Sqr(QSp))*(Lambda12*(Lambda12).adjoint()*Lambda12) - 2*(Lambda12*
      (Lambda12).adjoint()*Lambda12*(Lambda12).adjoint()*Lambda12));


   return beta_Lambda12;
}

} // namespace flexiblesusy
