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

// File generated at Sun 24 Aug 2014 16:10:20

#include "lowE6SSM_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of mDx2.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> lowE6SSM_soft_parameters::calc_beta_mDx2_one_loop(const Soft_traces& soft_traces) const
{
   const auto QDxp = INPUT(QDxp);
   const double Tr11 = TRACE_STRUCT.Tr11;
   const double Tr14 = TRACE_STRUCT.Tr14;


   Eigen::Matrix<double,3,3> beta_mDx2;

   beta_mDx2 = oneOver16PiSqr*(2*ms2*(Kappa.conjugate()*(Kappa).transpose
      ()) + 2*(TKappa.conjugate()*(TKappa).transpose()) + mDx2*Kappa.conjugate(
      )*(Kappa).transpose() + 2*(Kappa.conjugate()*mDxbar2*(Kappa).transpose())
      + Kappa.conjugate()*(Kappa).transpose()*mDx2 - 0.5163977794943222*g1*
      Tr11*UNITMATRIX(3) + 2*gN*QDxp*Tr14*UNITMATRIX(3) - 0.5333333333333333*
      AbsSqr(MassB)*Sqr(g1)*UNITMATRIX(3) - 10.666666666666666*AbsSqr(MassG)*
      Sqr(g3)*UNITMATRIX(3) - 8*AbsSqr(MassBp)*Sqr(gN)*Sqr(QDxp)*UNITMATRIX(3))
      ;


   return beta_mDx2;
}

/**
 * Calculates the two-loop beta function of mDx2.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> lowE6SSM_soft_parameters::calc_beta_mDx2_two_loop(const Soft_traces& soft_traces) const
{
   const auto QDxp = INPUT(QDxp);
   const auto Qdp = INPUT(Qdp);
   const auto QDxbarp = INPUT(QDxbarp);
   const auto Qep = INPUT(Qep);
   const auto QH1p = INPUT(QH1p);
   const auto QH2p = INPUT(QH2p);
   const auto QHpbarp = INPUT(QHpbarp);
   const auto QHpp = INPUT(QHpp);
   const auto QLp = INPUT(QLp);
   const auto QQp = INPUT(QQp);
   const auto Qup = INPUT(Qup);
   const auto QSp = INPUT(QSp);
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;
   const double traceconjTKappaTpTKappa =
      TRACE_STRUCT.traceconjTKappaTpTKappa;
   const double traceconjTLambda12TpTLambda12 =
      TRACE_STRUCT.traceconjTLambda12TpTLambda12;
   const double tracemH1I2AdjLambda12Lambda12 =
      TRACE_STRUCT.tracemH1I2AdjLambda12Lambda12;
   const double traceKappaAdjKappaconjmDx2 =
      TRACE_STRUCT.traceKappaAdjKappaconjmDx2;
   const double traceKappaconjmDxbar2AdjKappa =
      TRACE_STRUCT.traceKappaconjmDxbar2AdjKappa;
   const double traceLambda12AdjLambda12conjmH2I2 =
      TRACE_STRUCT.traceLambda12AdjLambda12conjmH2I2;
   const double traceconjTKappaTpKappa =
      TRACE_STRUCT.traceconjTKappaTpKappa;
   const double traceconjTLambda12TpLambda12 =
      TRACE_STRUCT.traceconjTLambda12TpLambda12;
   const double traceAdjKappaTKappa = TRACE_STRUCT.traceAdjKappaTKappa;
   const double traceAdjLambda12TLambda12 =
      TRACE_STRUCT.traceAdjLambda12TLambda12;
   const double Tr2U111 = TRACE_STRUCT.Tr2U111;
   const double Tr2U114 = TRACE_STRUCT.Tr2U114;
   const double Tr31 = TRACE_STRUCT.Tr31;
   const double Tr23 = TRACE_STRUCT.Tr23;
   const double Tr2U141 = TRACE_STRUCT.Tr2U141;
   const double Tr2U144 = TRACE_STRUCT.Tr2U144;
   const double Tr34 = TRACE_STRUCT.Tr34;


   Eigen::Matrix<double,3,3> beta_mDx2;

   beta_mDx2 = twoLoop*(-6*traceconjTKappaTpTKappa*(Kappa.conjugate()*(
      Kappa).transpose()) - 4*traceconjTLambda12TpTLambda12*(Kappa.conjugate()*
      (Kappa).transpose()) - 12*ms2*traceKappaAdjKappa*(Kappa.conjugate()*(
      Kappa).transpose()) - 6*traceKappaAdjKappaconjmDx2*(Kappa.conjugate()*(
      Kappa).transpose()) - 6*traceKappaconjmDxbar2AdjKappa*(Kappa.conjugate()*
      (Kappa).transpose()) - 8*ms2*traceLambda12AdjLambda12*(Kappa.conjugate()*
      (Kappa).transpose()) - 4*traceLambda12AdjLambda12conjmH2I2*(
      Kappa.conjugate()*(Kappa).transpose()) - 4*tracemH1I2AdjLambda12Lambda12*
      (Kappa.conjugate()*(Kappa).transpose()) - 4*mHd2*AbsSqr(Lambdax)*(
      Kappa.conjugate()*(Kappa).transpose()) - 4*mHu2*AbsSqr(Lambdax)*(
      Kappa.conjugate()*(Kappa).transpose()) - 8*ms2*AbsSqr(Lambdax)*(
      Kappa.conjugate()*(Kappa).transpose()) - 4*AbsSqr(TLambdax)*(
      Kappa.conjugate()*(Kappa).transpose()) + 4*ms2*Sqr(gN)*Sqr(QDxbarp)*(
      Kappa.conjugate()*(Kappa).transpose()) - 4*ms2*Sqr(gN)*Sqr(QDxp)*(
      Kappa.conjugate()*(Kappa).transpose()) + 4*ms2*Sqr(gN)*Sqr(QSp)*(
      Kappa.conjugate()*(Kappa).transpose()) - 6*traceconjTKappaTpKappa*(
      Kappa.conjugate()*(TKappa).transpose()) - 4*traceconjTLambda12TpLambda12*
      (Kappa.conjugate()*(TKappa).transpose()) - 4*Conj(TLambdax)*Lambdax*(
      Kappa.conjugate()*(TKappa).transpose()) - 6*traceAdjKappaTKappa*(
      TKappa.conjugate()*(Kappa).transpose()) - 4*traceAdjLambda12TLambda12*(
      TKappa.conjugate()*(Kappa).transpose()) - 4*MassBp*Sqr(gN)*Sqr(QDxbarp)*(
      TKappa.conjugate()*(Kappa).transpose()) + 4*MassBp*Sqr(gN)*Sqr(QDxp)*(
      TKappa.conjugate()*(Kappa).transpose()) - 4*MassBp*Sqr(gN)*Sqr(QSp)*(
      TKappa.conjugate()*(Kappa).transpose()) - 4*Conj(Lambdax)*TLambdax*(
      TKappa.conjugate()*(Kappa).transpose()) - 6*traceKappaAdjKappa*(
      TKappa.conjugate()*(TKappa).transpose()) - 4*traceLambda12AdjLambda12*(
      TKappa.conjugate()*(TKappa).transpose()) - 4*AbsSqr(Lambdax)*(
      TKappa.conjugate()*(TKappa).transpose()) + 4*Sqr(gN)*Sqr(QDxbarp)*(
      TKappa.conjugate()*(TKappa).transpose()) - 4*Sqr(gN)*Sqr(QDxp)*(
      TKappa.conjugate()*(TKappa).transpose()) + 4*Sqr(gN)*Sqr(QSp)*(
      TKappa.conjugate()*(TKappa).transpose()) - 3*traceKappaAdjKappa*(mDx2*
      Kappa.conjugate()*(Kappa).transpose()) - 2*traceLambda12AdjLambda12*(mDx2
      *Kappa.conjugate()*(Kappa).transpose()) - 2*AbsSqr(Lambdax)*(mDx2*
      Kappa.conjugate()*(Kappa).transpose()) + 2*Sqr(gN)*Sqr(QDxbarp)*(mDx2*
      Kappa.conjugate()*(Kappa).transpose()) - 2*Sqr(gN)*Sqr(QDxp)*(mDx2*
      Kappa.conjugate()*(Kappa).transpose()) + 2*Sqr(gN)*Sqr(QSp)*(mDx2*
      Kappa.conjugate()*(Kappa).transpose()) - 6*traceKappaAdjKappa*(
      Kappa.conjugate()*mDxbar2*(Kappa).transpose()) - 4*
      traceLambda12AdjLambda12*(Kappa.conjugate()*mDxbar2*(Kappa).transpose())
      - 4*AbsSqr(Lambdax)*(Kappa.conjugate()*mDxbar2*(Kappa).transpose()) + 4*
      Sqr(gN)*Sqr(QDxbarp)*(Kappa.conjugate()*mDxbar2*(Kappa).transpose()) - 4*
      Sqr(gN)*Sqr(QDxp)*(Kappa.conjugate()*mDxbar2*(Kappa).transpose()) + 4*Sqr
      (gN)*Sqr(QSp)*(Kappa.conjugate()*mDxbar2*(Kappa).transpose()) - 3*
      traceKappaAdjKappa*(Kappa.conjugate()*(Kappa).transpose()*mDx2) - 2*
      traceLambda12AdjLambda12*(Kappa.conjugate()*(Kappa).transpose()*mDx2) - 2
      *AbsSqr(Lambdax)*(Kappa.conjugate()*(Kappa).transpose()*mDx2) + 2*Sqr(gN)
      *Sqr(QDxbarp)*(Kappa.conjugate()*(Kappa).transpose()*mDx2) - 2*Sqr(gN)*
      Sqr(QDxp)*(Kappa.conjugate()*(Kappa).transpose()*mDx2) + 2*Sqr(gN)*Sqr(
      QSp)*(Kappa.conjugate()*(Kappa).transpose()*mDx2) - 4*ms2*(
      Kappa.conjugate()*(Kappa).transpose()*Kappa.conjugate()*(Kappa).transpose
      ()) - 2*(Kappa.conjugate()*(Kappa).transpose()*TKappa.conjugate()*(TKappa
      ).transpose()) - 2*(Kappa.conjugate()*(TKappa).transpose()*
      TKappa.conjugate()*(Kappa).transpose()) - 2*(TKappa.conjugate()*(Kappa)
      .transpose()*Kappa.conjugate()*(TKappa).transpose()) - 2*(
      TKappa.conjugate()*(TKappa).transpose()*Kappa.conjugate()*(Kappa)
      .transpose()) - mDx2*Kappa.conjugate()*(Kappa).transpose()*
      Kappa.conjugate()*(Kappa).transpose() - 2*(Kappa.conjugate()*mDxbar2*(
      Kappa).transpose()*Kappa.conjugate()*(Kappa).transpose()) - 2*(
      Kappa.conjugate()*(Kappa).transpose()*mDx2*Kappa.conjugate()*(Kappa)
      .transpose()) - 2*(Kappa.conjugate()*(Kappa).transpose()*Kappa.conjugate(
      )*mDxbar2*(Kappa).transpose()) - Kappa.conjugate()*(Kappa).transpose()*
      Kappa.conjugate()*(Kappa).transpose()*mDx2 + 10.666666666666666*Power(g3,
      4)*Tr23*UNITMATRIX(3) - 2.065591117977289*g1*gN*QDxp*Tr2U114*UNITMATRIX(3
      ) - 2.065591117977289*g1*gN*QDxp*Tr2U141*UNITMATRIX(3) -
      2.065591117977289*g1*Tr31*UNITMATRIX(3) + 8*gN*QDxp*Tr34*UNITMATRIX(3) +
      53.333333333333336*Power(g3,4)*AbsSqr(MassG)*UNITMATRIX(3) +
      0.5333333333333333*Tr2U111*Sqr(g1)*UNITMATRIX(3) + 2.8444444444444446*
      AbsSqr(MassG)*Sqr(g1)*Sqr(g3)*UNITMATRIX(3) + 1.4222222222222223*MassB*
      Conj(MassG)*Sqr(g1)*Sqr(g3)*UNITMATRIX(3) + 0.035555555555555556*Conj(
      MassB)*Sqr(g1)*(438*MassB*Sqr(g1) + 5*(8*(2*MassB + MassG)*Sqr(g3) + 3*(2
      *MassB + MassBp)*QDxp*(-9*Qdp - 9*QDxbarp + 11*QDxp - 9*Qep + 9*QH1p - 9*
      QH2p - 3*QHpbarp + 3*QHpp + 9*QLp - 9*QQp + 18*Qup)*Sqr(gN)))*UNITMATRIX(
      3) + 8*Tr2U144*Sqr(gN)*Sqr(QDxp)*UNITMATRIX(3) + 42.666666666666664*
      AbsSqr(MassG)*Sqr(g3)*Sqr(gN)*Sqr(QDxp)*UNITMATRIX(3) +
      21.333333333333332*MassBp*Conj(MassG)*Sqr(g3)*Sqr(gN)*Sqr(QDxp)*
      UNITMATRIX(3) + 0.26666666666666666*Conj(MassBp)*Sqr(gN)*(15*(Sqr(QDxbarp
      ) - Sqr(QDxp) + Sqr(QSp))*(2*MassBp*(Kappa.conjugate()*(Kappa).transpose(
      )) - Kappa.conjugate()*(TKappa).transpose()) + 2*QDxp*(-((MassB + 2*
      MassBp)*(9*Qdp + 9*QDxbarp - 11*QDxp + 9*Qep - 9*QH1p + 9*QH2p + 3*
      QHpbarp - 3*QHpp - 9*QLp + 9*QQp - 18*Qup)*Sqr(g1)) + 5*QDxp*(8*(2*MassBp
      + MassG)*Sqr(g3) + 9*MassBp*Sqr(gN)*(9*Sqr(Qdp) + 9*Sqr(QDxbarp) + 11*
      Sqr(QDxp) + 3*Sqr(Qep) + 6*Sqr(QH1p) + 6*Sqr(QH2p) + 2*Sqr(QHpbarp) + 2*
      Sqr(QHpp) + 6*Sqr(QLp) + 18*Sqr(QQp) + 3*Sqr(QSp) + 9*Sqr(Qup))))*
      UNITMATRIX(3)));


   return beta_mDx2;
}

} // namespace flexiblesusy
