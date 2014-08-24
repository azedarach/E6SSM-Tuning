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

// File generated at Sun 24 Aug 2014 16:10:19

#include "lowE6SSM_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of msI2.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,2,2> lowE6SSM_soft_parameters::calc_beta_msI2_one_loop(const Soft_traces& soft_traces) const
{
   const auto QSp = INPUT(QSp);
   const double Tr14 = TRACE_STRUCT.Tr14;


   Eigen::Matrix<double,2,2> beta_msI2;

   beta_msI2 = 2*gN*oneOver16PiSqr*QSp*(Tr14 - 4*gN*QSp*AbsSqr(MassBp))*
      UNITMATRIX(2);


   return beta_msI2;
}

/**
 * Calculates the two-loop beta function of msI2.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,2,2> lowE6SSM_soft_parameters::calc_beta_msI2_two_loop(const Soft_traces& soft_traces) const
{
   const auto QSp = INPUT(QSp);
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
   const auto Qup = INPUT(Qup);
   const double Tr2U144 = TRACE_STRUCT.Tr2U144;
   const double Tr34 = TRACE_STRUCT.Tr34;


   Eigen::Matrix<double,2,2> beta_msI2;

   beta_msI2 = 8*gN*QSp*twoLoop*(gN*QSp*Tr2U144 + Tr34 + 3*Power(gN,3)*
      QSp*AbsSqr(MassBp)*(9*Sqr(Qdp) + 9*Sqr(QDxbarp) + 9*Sqr(QDxp) + 3*Sqr(Qep
      ) + 6*Sqr(QH1p) + 6*Sqr(QH2p) + 2*Sqr(QHpbarp) + 2*Sqr(QHpp) + 6*Sqr(QLp)
      + 18*Sqr(QQp) + 5*Sqr(QSp) + 9*Sqr(Qup)))*UNITMATRIX(2);


   return beta_msI2;
}

} // namespace flexiblesusy
