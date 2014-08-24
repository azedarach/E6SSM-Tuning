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

// File generated at Sun 24 Aug 2014 16:10:21

#include "lowE6SSM_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of mDxbar2.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> lowE6SSM_soft_parameters::calc_beta_mDxbar2_one_loop(const Soft_traces& soft_traces) const
{
   const auto QDxbarp = INPUT(QDxbarp);
   const double Tr11 = TRACE_STRUCT.Tr11;
   const double Tr14 = TRACE_STRUCT.Tr14;


   Eigen::Matrix<double,3,3> beta_mDxbar2;

   beta_mDxbar2 = oneOver16PiSqr*(2*ms2*((Kappa).transpose()*
      Kappa.conjugate()) + 2*((TKappa).transpose()*TKappa.conjugate()) +
      mDxbar2*(Kappa).transpose()*Kappa.conjugate() + 2*((Kappa).transpose()*
      mDx2*Kappa.conjugate()) + (Kappa).transpose()*Kappa.conjugate()*mDxbar2 +
      0.5163977794943222*g1*Tr11*UNITMATRIX(3) + 2*gN*QDxbarp*Tr14*UNITMATRIX(
      3) - 0.5333333333333333*AbsSqr(MassB)*Sqr(g1)*UNITMATRIX(3) -
      10.666666666666666*AbsSqr(MassG)*Sqr(g3)*UNITMATRIX(3) - 8*AbsSqr(MassBp)
      *Sqr(gN)*Sqr(QDxbarp)*UNITMATRIX(3));


   return beta_mDxbar2;
}

/**
 * Calculates the two-loop beta function of mDxbar2.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> lowE6SSM_soft_parameters::calc_beta_mDxbar2_two_loop(const Soft_traces& soft_traces) const
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
   const double traceAdjKappaTKappa = TRACE_STRUCT.traceAdjKappaTKappa;
   const double traceAdjLambda12TLambda12 =
      TRACE_STRUCT.traceAdjLambda12TLambda12;
   const double traceconjTKappaTpKappa =
      TRACE_STRUCT.traceconjTKappaTpKappa;
   const double traceconjTLambda12TpLambda12 =
      TRACE_STRUCT.traceconjTLambda12TpLambda12;
   const double Tr2U111 = TRACE_STRUCT.Tr2U111;
   const double Tr2U114 = TRACE_STRUCT.Tr2U114;
   const double Tr31 = TRACE_STRUCT.Tr31;
   const double Tr23 = TRACE_STRUCT.Tr23;
   const double Tr2U141 = TRACE_STRUCT.Tr2U141;
   const double Tr2U144 = TRACE_STRUCT.Tr2U144;
   const double Tr34 = TRACE_STRUCT.Tr34;


   Eigen::Matrix<double,3,3> beta_mDxbar2;

   beta_mDxbar2 = twoLoop*(-6*traceconjTKappaTpTKappa*((Kappa).transpose(
      )*Kappa.conjugate()) - 4*traceconjTLambda12TpTLambda12*((Kappa).transpose
      ()*Kappa.conjugate()) - 12*ms2*traceKappaAdjKappa*((Kappa).transpose()*
      Kappa.conjugate()) - 6*traceKappaAdjKappaconjmDx2*((Kappa).transpose()*
      Kappa.conjugate()) - 6*traceKappaconjmDxbar2AdjKappa*((Kappa).transpose()
      *Kappa.conjugate()) - 8*ms2*traceLambda12AdjLambda12*((Kappa).transpose()
      *Kappa.conjugate()) - 4*traceLambda12AdjLambda12conjmH2I2*((Kappa)
      .transpose()*Kappa.conjugate()) - 4*tracemH1I2AdjLambda12Lambda12*((Kappa
      ).transpose()*Kappa.conjugate()) - 4*mHd2*AbsSqr(Lambdax)*((Kappa)
      .transpose()*Kappa.conjugate()) - 4*mHu2*AbsSqr(Lambdax)*((Kappa)
      .transpose()*Kappa.conjugate()) - 8*ms2*AbsSqr(Lambdax)*((Kappa)
      .transpose()*Kappa.conjugate()) - 4*AbsSqr(TLambdax)*((Kappa).transpose()
      *Kappa.conjugate()) - 4*ms2*Sqr(gN)*Sqr(QDxbarp)*((Kappa).transpose()*
      Kappa.conjugate()) + 4*ms2*Sqr(gN)*Sqr(QDxp)*((Kappa).transpose()*
      Kappa.conjugate()) + 4*ms2*Sqr(gN)*Sqr(QSp)*((Kappa).transpose()*
      Kappa.conjugate()) - 6*traceAdjKappaTKappa*((Kappa).transpose()*
      TKappa.conjugate()) - 4*traceAdjLambda12TLambda12*((Kappa).transpose()*
      TKappa.conjugate()) + 4*MassBp*Sqr(gN)*Sqr(QDxbarp)*((Kappa).transpose()*
      TKappa.conjugate()) - 4*MassBp*Sqr(gN)*Sqr(QDxp)*((Kappa).transpose()*
      TKappa.conjugate()) - 4*MassBp*Sqr(gN)*Sqr(QSp)*((Kappa).transpose()*
      TKappa.conjugate()) - 4*Conj(Lambdax)*TLambdax*((Kappa).transpose()*
      TKappa.conjugate()) - 6*traceconjTKappaTpKappa*((TKappa).transpose()*
      Kappa.conjugate()) - 4*traceconjTLambda12TpLambda12*((TKappa).transpose()
      *Kappa.conjugate()) - 4*Conj(TLambdax)*Lambdax*((TKappa).transpose()*
      Kappa.conjugate()) - 6*traceKappaAdjKappa*((TKappa).transpose()*
      TKappa.conjugate()) - 4*traceLambda12AdjLambda12*((TKappa).transpose()*
      TKappa.conjugate()) - 4*AbsSqr(Lambdax)*((TKappa).transpose()*
      TKappa.conjugate()) - 4*Sqr(gN)*Sqr(QDxbarp)*((TKappa).transpose()*
      TKappa.conjugate()) + 4*Sqr(gN)*Sqr(QDxp)*((TKappa).transpose()*
      TKappa.conjugate()) + 4*Sqr(gN)*Sqr(QSp)*((TKappa).transpose()*
      TKappa.conjugate()) - 3*traceKappaAdjKappa*(mDxbar2*(Kappa).transpose()*
      Kappa.conjugate()) - 2*traceLambda12AdjLambda12*(mDxbar2*(Kappa)
      .transpose()*Kappa.conjugate()) - 2*AbsSqr(Lambdax)*(mDxbar2*(Kappa)
      .transpose()*Kappa.conjugate()) - 2*Sqr(gN)*Sqr(QDxbarp)*(mDxbar2*(Kappa)
      .transpose()*Kappa.conjugate()) + 2*Sqr(gN)*Sqr(QDxp)*(mDxbar2*(Kappa)
      .transpose()*Kappa.conjugate()) + 2*Sqr(gN)*Sqr(QSp)*(mDxbar2*(Kappa)
      .transpose()*Kappa.conjugate()) - 6*traceKappaAdjKappa*((Kappa).transpose
      ()*mDx2*Kappa.conjugate()) - 4*traceLambda12AdjLambda12*((Kappa)
      .transpose()*mDx2*Kappa.conjugate()) - 4*AbsSqr(Lambdax)*((Kappa)
      .transpose()*mDx2*Kappa.conjugate()) - 4*Sqr(gN)*Sqr(QDxbarp)*((Kappa)
      .transpose()*mDx2*Kappa.conjugate()) + 4*Sqr(gN)*Sqr(QDxp)*((Kappa)
      .transpose()*mDx2*Kappa.conjugate()) + 4*Sqr(gN)*Sqr(QSp)*((Kappa)
      .transpose()*mDx2*Kappa.conjugate()) - 3*traceKappaAdjKappa*((Kappa)
      .transpose()*Kappa.conjugate()*mDxbar2) - 2*traceLambda12AdjLambda12*((
      Kappa).transpose()*Kappa.conjugate()*mDxbar2) - 2*AbsSqr(Lambdax)*((Kappa
      ).transpose()*Kappa.conjugate()*mDxbar2) - 2*Sqr(gN)*Sqr(QDxbarp)*((Kappa
      ).transpose()*Kappa.conjugate()*mDxbar2) + 2*Sqr(gN)*Sqr(QDxp)*((Kappa)
      .transpose()*Kappa.conjugate()*mDxbar2) + 2*Sqr(gN)*Sqr(QSp)*((Kappa)
      .transpose()*Kappa.conjugate()*mDxbar2) - 4*ms2*((Kappa).transpose()*
      Kappa.conjugate()*(Kappa).transpose()*Kappa.conjugate()) - 2*((Kappa)
      .transpose()*Kappa.conjugate()*(TKappa).transpose()*TKappa.conjugate()) -
      2*((Kappa).transpose()*TKappa.conjugate()*(TKappa).transpose()*
      Kappa.conjugate()) - 2*((TKappa).transpose()*Kappa.conjugate()*(Kappa)
      .transpose()*TKappa.conjugate()) - 2*((TKappa).transpose()*
      TKappa.conjugate()*(Kappa).transpose()*Kappa.conjugate()) - mDxbar2*(
      Kappa).transpose()*Kappa.conjugate()*(Kappa).transpose()*Kappa.conjugate(
      ) - 2*((Kappa).transpose()*mDx2*Kappa.conjugate()*(Kappa).transpose()*
      Kappa.conjugate()) - 2*((Kappa).transpose()*Kappa.conjugate()*mDxbar2*(
      Kappa).transpose()*Kappa.conjugate()) - 2*((Kappa).transpose()*
      Kappa.conjugate()*(Kappa).transpose()*mDx2*Kappa.conjugate()) - (Kappa)
      .transpose()*Kappa.conjugate()*(Kappa).transpose()*Kappa.conjugate()*
      mDxbar2 + 10.666666666666666*Power(g3,4)*Tr23*UNITMATRIX(3) +
      2.065591117977289*g1*gN*QDxbarp*Tr2U114*UNITMATRIX(3) + 2.065591117977289
      *g1*gN*QDxbarp*Tr2U141*UNITMATRIX(3) + 2.065591117977289*g1*Tr31*
      UNITMATRIX(3) + 8*gN*QDxbarp*Tr34*UNITMATRIX(3) + 53.333333333333336*
      Power(g3,4)*AbsSqr(MassG)*UNITMATRIX(3) + 0.5333333333333333*Tr2U111*Sqr(
      g1)*UNITMATRIX(3) + 2.8444444444444446*AbsSqr(MassG)*Sqr(g1)*Sqr(g3)*
      UNITMATRIX(3) + 1.4222222222222223*MassB*Conj(MassG)*Sqr(g1)*Sqr(g3)*
      UNITMATRIX(3) + 0.035555555555555556*Conj(MassB)*Sqr(g1)*(438*MassB*Sqr(
      g1) + 5*(8*(2*MassB + MassG)*Sqr(g3) + 3*(2*MassB + MassBp)*QDxbarp*(9*
      Qdp + 11*QDxbarp - 9*QDxp + 9*Qep - 9*QH1p + 9*QH2p + 3*QHpbarp - 3*QHpp
      - 9*QLp + 9*QQp - 18*Qup)*Sqr(gN)))*UNITMATRIX(3) + 8*Tr2U144*Sqr(gN)*Sqr
      (QDxbarp)*UNITMATRIX(3) + 42.666666666666664*AbsSqr(MassG)*Sqr(g3)*Sqr(gN
      )*Sqr(QDxbarp)*UNITMATRIX(3) + 21.333333333333332*MassBp*Conj(MassG)*Sqr(
      g3)*Sqr(gN)*Sqr(QDxbarp)*UNITMATRIX(3) + 0.26666666666666666*Conj(MassBp)
      *Sqr(gN)*(-15*(Sqr(QDxbarp) - Sqr(QDxp) - Sqr(QSp))*(2*MassBp*((Kappa)
      .transpose()*Kappa.conjugate()) - (TKappa).transpose()*Kappa.conjugate())
      + 2*QDxbarp*((MassB + 2*MassBp)*(9*Qdp + 11*QDxbarp - 9*QDxp + 9*Qep - 9
      *QH1p + 9*QH2p + 3*QHpbarp - 3*QHpp - 9*QLp + 9*QQp - 18*Qup)*Sqr(g1) + 5
      *QDxbarp*(8*(2*MassBp + MassG)*Sqr(g3) + 9*MassBp*Sqr(gN)*(9*Sqr(Qdp) +
      11*Sqr(QDxbarp) + 9*Sqr(QDxp) + 3*Sqr(Qep) + 6*Sqr(QH1p) + 6*Sqr(QH2p) +
      2*Sqr(QHpbarp) + 2*Sqr(QHpp) + 6*Sqr(QLp) + 18*Sqr(QQp) + 3*Sqr(QSp) + 9*
      Sqr(Qup))))*UNITMATRIX(3)));


   return beta_mDxbar2;
}

} // namespace flexiblesusy
