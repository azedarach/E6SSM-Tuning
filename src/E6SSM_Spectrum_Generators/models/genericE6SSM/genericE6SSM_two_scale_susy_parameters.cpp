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

#include <iostream>

namespace flexiblesusy {

#define CLASSNAME genericE6SSM_susy_parameters
#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

genericE6SSM_susy_parameters::genericE6SSM_susy_parameters(const genericE6SSM_input_parameters& input_)
   : Beta_function()
   , Yd(Eigen::Matrix<double,3,3>::Zero()), Ye(Eigen::Matrix<double,3,3>::Zero(
   )), Kappa(Eigen::Matrix<double,3,3>::Zero()), Lambda12(Eigen::Matrix<double,
   2,2>::Zero()), Lambdax(0), Yu(Eigen::Matrix<double,3,3>::Zero()), MuPr(0),
   g1(0), g2(0), g3(0), gN(0), vd(0), vu(0), vs(0)

   , input(input_)
{
   set_number_of_parameters(numberOfParameters);
}

genericE6SSM_susy_parameters::genericE6SSM_susy_parameters(
   double scale_, double loops_, double thresholds_,
   const genericE6SSM_input_parameters& input_
   , const Eigen::Matrix<double,3,3>& Yd_, const Eigen::Matrix<double,3,3>& Ye_
   , const Eigen::Matrix<double,3,3>& Kappa_, const Eigen::Matrix<double,2,2>&
   Lambda12_, double Lambdax_, const Eigen::Matrix<double,3,3>& Yu_, double
   MuPr_, double g1_, double g2_, double g3_, double gN_, double vd_, double
   vu_, double vs_

)
   : Beta_function()
   , Yd(Yd_), Ye(Ye_), Kappa(Kappa_), Lambda12(Lambda12_), Lambdax(Lambdax_),
   Yu(Yu_), MuPr(MuPr_), g1(g1_), g2(g2_), g3(g3_), gN(gN_), vd(vd_), vu(vu_),
   vs(vs_)

   , input(input_)
{
   set_number_of_parameters(numberOfParameters);
   set_scale(scale_);
   set_loops(loops_);
   set_thresholds(thresholds_);
}

Eigen::ArrayXd genericE6SSM_susy_parameters::beta() const
{
   return calc_beta().get();
}

genericE6SSM_susy_parameters genericE6SSM_susy_parameters::calc_beta() const
{
   Susy_traces susy_traces;
   calc_susy_traces(susy_traces);

   Eigen::Matrix<double,3,3> beta_Yd(calc_beta_Yd_one_loop(TRACE_STRUCT));
   Eigen::Matrix<double,3,3> beta_Ye(calc_beta_Ye_one_loop(TRACE_STRUCT));
   Eigen::Matrix<double,3,3> beta_Kappa(calc_beta_Kappa_one_loop(TRACE_STRUCT))
      ;
   Eigen::Matrix<double,2,2> beta_Lambda12(calc_beta_Lambda12_one_loop(
      TRACE_STRUCT));
   double beta_Lambdax(calc_beta_Lambdax_one_loop(TRACE_STRUCT));
   Eigen::Matrix<double,3,3> beta_Yu(calc_beta_Yu_one_loop(TRACE_STRUCT));
   double beta_MuPr(calc_beta_MuPr_one_loop(TRACE_STRUCT));
   double beta_g1(calc_beta_g1_one_loop(TRACE_STRUCT));
   double beta_g2(calc_beta_g2_one_loop(TRACE_STRUCT));
   double beta_g3(calc_beta_g3_one_loop(TRACE_STRUCT));
   double beta_gN(calc_beta_gN_one_loop(TRACE_STRUCT));
   double beta_vd(calc_beta_vd_one_loop(TRACE_STRUCT));
   double beta_vu(calc_beta_vu_one_loop(TRACE_STRUCT));
   double beta_vs(calc_beta_vs_one_loop(TRACE_STRUCT));

   if (get_loops() > 1) {
      beta_Yd += calc_beta_Yd_two_loop(TRACE_STRUCT);
      beta_Ye += calc_beta_Ye_two_loop(TRACE_STRUCT);
      beta_Kappa += calc_beta_Kappa_two_loop(TRACE_STRUCT);
      beta_Lambda12 += calc_beta_Lambda12_two_loop(TRACE_STRUCT);
      beta_Lambdax += calc_beta_Lambdax_two_loop(TRACE_STRUCT);
      beta_Yu += calc_beta_Yu_two_loop(TRACE_STRUCT);
      beta_MuPr += calc_beta_MuPr_two_loop(TRACE_STRUCT);
      beta_g1 += calc_beta_g1_two_loop(TRACE_STRUCT);
      beta_g2 += calc_beta_g2_two_loop(TRACE_STRUCT);
      beta_g3 += calc_beta_g3_two_loop(TRACE_STRUCT);
      beta_gN += calc_beta_gN_two_loop(TRACE_STRUCT);
      beta_vd += calc_beta_vd_two_loop(TRACE_STRUCT);
      beta_vu += calc_beta_vu_two_loop(TRACE_STRUCT);
      beta_vs += calc_beta_vs_two_loop(TRACE_STRUCT);

   }


   return genericE6SSM_susy_parameters(get_scale(), get_loops(), get_thresholds(), input,
                    beta_Yd, beta_Ye, beta_Kappa, beta_Lambda12, beta_Lambdax, beta_Yu, beta_MuPr, beta_g1, beta_g2, beta_g3, beta_gN, beta_vd, beta_vu, beta_vs);
}

void genericE6SSM_susy_parameters::clear()
{
   reset();
   Yd = Eigen::Matrix<double,3,3>::Zero();
   Ye = Eigen::Matrix<double,3,3>::Zero();
   Kappa = Eigen::Matrix<double,3,3>::Zero();
   Lambda12 = Eigen::Matrix<double,2,2>::Zero();
   Lambdax = 0.;
   Yu = Eigen::Matrix<double,3,3>::Zero();
   MuPr = 0.;
   g1 = 0.;
   g2 = 0.;
   g3 = 0.;
   gN = 0.;
   vd = 0.;
   vu = 0.;
   vs = 0.;

}

Eigen::Matrix<double,3,3> CLASSNAME::get_SqSq() const
{
   const double twoLoop = oneOver16PiSqr * oneOver16PiSqr;
   Eigen::Matrix<double,3,3> anomDim;
   const auto QQp = INPUT(QQp);
   const auto QH2p = INPUT(QH2p);
   const auto Qup = INPUT(Qup);
   const auto Qdp = INPUT(Qdp);
   const auto QH1p = INPUT(QH1p);
   const auto QDxbarp = INPUT(QDxbarp);
   const auto QDxp = INPUT(QDxp);
   const auto Qep = INPUT(Qep);
   const auto QHpbarp = INPUT(QHpbarp);
   const auto QHpp = INPUT(QHpp);
   const auto QLp = INPUT(QLp);
   const auto QSp = INPUT(QSp);

   anomDim = oneOver16PiSqr*(Yd.adjoint()*Yd + Yu.adjoint()*Yu -
      0.03333333333333333*(Sqr(g1) + 45*Sqr(g2) + 80*Sqr(g3) + 60*Sqr(gN)*Sqr(
      QQp))*UNITMATRIX(3));

   if (get_loops() > 1) {
      anomDim += twoLoop*(-(AbsSqr(Lambdax)*(Yu.adjoint()*Yu)) + 0.8*
         Sqr(g1)*(Yu.adjoint()*Yu) + 2*Sqr(gN)*Sqr(QH2p)*(Yu.adjoint()*Yu) - 2*
         Sqr(gN)*Sqr(QQp)*(Yu.adjoint()*Yu) + 2*Sqr(gN)*Sqr(Qup)*(Yu.adjoint()*
         Yu) - 2*(Yd.adjoint()*Yd*Yd.adjoint()*Yd) - 2*(Yu.adjoint()*Yu*
         Yu.adjoint()*Yu) + Yd.adjoint()*Yd*(-AbsSqr(Lambdax) + 0.4*Sqr(g1) + 2
         *Sqr(gN)*Sqr(Qdp) + 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QQp) - 3*(Yd*
         Yd.adjoint()).trace() - (Ye*Ye.adjoint()).trace()) - 3*(Yu.adjoint()*
         Yu)*(Yu*Yu.adjoint()).trace() + 0.0011111111111111111*(289*Power(g1,4)
         + 10*Sqr(g1)*(9*Sqr(g2) + 4*(4*Sqr(g3) + 3*QQp*(9*Qdp + 9*QDxbarp - 9
         *QDxp + 9*Qep - 9*QH1p + 9*QH2p + 3*QHpbarp - 3*QHpp - 9*QLp + 10*QQp
         - 18*Qup)*Sqr(gN))) + 25*(297*Power(g2,4) + 72*Sqr(g2)*(4*Sqr(g3) + 3*
         Sqr(gN)*Sqr(QQp)) + 8*(32*Power(g3,4) + 48*Sqr(g3)*Sqr(gN)*Sqr(QQp) +
         9*Power(gN,4)*Sqr(QQp)*(9*Sqr(Qdp) + 9*Sqr(QDxbarp) + 9*Sqr(QDxp) + 3*
         Sqr(Qep) + 6*Sqr(QH1p) + 6*Sqr(QH2p) + 2*Sqr(QHpbarp) + 2*Sqr(QHpp) +
         6*Sqr(QLp) + 20*Sqr(QQp) + 3*Sqr(QSp) + 9*Sqr(Qup)))))*UNITMATRIX(3));
   }

   return anomDim;
}

Eigen::Matrix<double,3,3> CLASSNAME::get_SlSl() const
{
   const double twoLoop = oneOver16PiSqr * oneOver16PiSqr;
   Eigen::Matrix<double,3,3> anomDim;
   const auto QLp = INPUT(QLp);
   const auto Qep = INPUT(Qep);
   const auto QH1p = INPUT(QH1p);
   const auto Qdp = INPUT(Qdp);
   const auto QDxbarp = INPUT(QDxbarp);
   const auto QDxp = INPUT(QDxp);
   const auto QH2p = INPUT(QH2p);
   const auto QHpbarp = INPUT(QHpbarp);
   const auto QHpp = INPUT(QHpp);
   const auto QQp = INPUT(QQp);
   const auto Qup = INPUT(Qup);
   const auto QSp = INPUT(QSp);

   anomDim = oneOver16PiSqr*(Ye.adjoint()*Ye - 0.1*(3*Sqr(g1) + 15*Sqr(g2
      ) + 20*Sqr(gN)*Sqr(QLp))*UNITMATRIX(3));

   if (get_loops() > 1) {
      anomDim += twoLoop*(-2*(Ye.adjoint()*Ye*Ye.adjoint()*Ye) +
         Ye.adjoint()*Ye*(-AbsSqr(Lambdax) + 1.2*Sqr(g1) + 2*Sqr(gN)*Sqr(Qep) +
         2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QLp) - 3*(Yd*Yd.adjoint()).trace(
         ) - (Ye*Ye.adjoint()).trace()) + 0.01*(297*Power(g1,4) + 30*Sqr(g1)*(3
         *Sqr(g2) + 4*QLp*(-3*Qdp - 3*QDxbarp + 3*QDxp - 3*Qep + 3*QH1p - 3*
         QH2p - QHpbarp + QHpp + 4*QLp - 3*QQp + 6*Qup)*Sqr(gN)) + 25*(33*Power
         (g2,4) + 24*Sqr(g2)*Sqr(gN)*Sqr(QLp) + 8*Power(gN,4)*Sqr(QLp)*(9*Sqr(
         Qdp) + 9*Sqr(QDxbarp) + 9*Sqr(QDxp) + 3*Sqr(Qep) + 6*Sqr(QH1p) + 6*Sqr
         (QH2p) + 2*Sqr(QHpbarp) + 2*Sqr(QHpp) + 8*Sqr(QLp) + 18*Sqr(QQp) + 3*
         Sqr(QSp) + 9*Sqr(Qup))))*UNITMATRIX(3));
   }

   return anomDim;
}

double CLASSNAME::get_SHdSHd() const
{
   const double twoLoop = oneOver16PiSqr * oneOver16PiSqr;
   double anomDim = 0;
   const auto QH1p = INPUT(QH1p);
   const auto Qdp = INPUT(Qdp);
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

   anomDim = oneOver16PiSqr*(AbsSqr(Lambdax) - 0.3*Sqr(g1) - 1.5*Sqr(g2)
      - 2*Sqr(gN)*Sqr(QH1p) + 3*(Yd*Yd.adjoint()).trace() + (Ye*Ye.adjoint())
      .trace());

   if (get_loops() > 1) {
      anomDim += twoLoop*(2.97*Power(g1,4) + 8.25*Power(g2,4) + 16*
         Power(gN,4)*Power(QH1p,4) + 0.9*Sqr(g1)*Sqr(g2) - 3.6*Qdp*QH1p*Sqr(g1)
         *Sqr(gN) - 3.6*QDxbarp*QH1p*Sqr(g1)*Sqr(gN) + 3.6*QDxp*QH1p*Sqr(g1)*
         Sqr(gN) - 3.6*Qep*QH1p*Sqr(g1)*Sqr(gN) - 3.6*QH1p*QH2p*Sqr(g1)*Sqr(gN)
         - 1.2*QH1p*QHpbarp*Sqr(g1)*Sqr(gN) + 1.2*QH1p*QHpp*Sqr(g1)*Sqr(gN) +
         3.6*QH1p*QLp*Sqr(g1)*Sqr(gN) - 3.6*QH1p*QQp*Sqr(g1)*Sqr(gN) + 7.2*QH1p
         *Qup*Sqr(g1)*Sqr(gN) + 4.8*Sqr(g1)*Sqr(gN)*Sqr(QH1p) + 6*Sqr(g2)*Sqr(
         gN)*Sqr(QH1p) + 18*Power(gN,4)*Sqr(Qdp)*Sqr(QH1p) + 18*Power(gN,4)*Sqr
         (QDxbarp)*Sqr(QH1p) + 18*Power(gN,4)*Sqr(QDxp)*Sqr(QH1p) + 6*Power(gN,
         4)*Sqr(Qep)*Sqr(QH1p) + 12*Power(gN,4)*Sqr(QH1p)*Sqr(QH2p) + 4*Power(
         gN,4)*Sqr(QH1p)*Sqr(QHpbarp) + 4*Power(gN,4)*Sqr(QH1p)*Sqr(QHpp) + 12*
         Power(gN,4)*Sqr(QH1p)*Sqr(QLp) + 36*Power(gN,4)*Sqr(QH1p)*Sqr(QQp) + 6
         *Power(gN,4)*Sqr(QH1p)*Sqr(QSp) + 18*Power(gN,4)*Sqr(QH1p)*Sqr(Qup) -
         3*Sqr(Conj(Lambdax))*Sqr(Lambdax) - 0.4*(Sqr(g1) - 5*(8*Sqr(g3) + 3*
         Sqr(gN)*(Sqr(Qdp) - Sqr(QH1p) + Sqr(QQp))))*(Yd*Yd.adjoint()).trace()
         + 1.2*Sqr(g1)*(Ye*Ye.adjoint()).trace() + 2*Sqr(gN)*Sqr(Qep)*(Ye*
         Ye.adjoint()).trace() - 2*Sqr(gN)*Sqr(QH1p)*(Ye*Ye.adjoint()).trace()
         + 2*Sqr(gN)*Sqr(QLp)*(Ye*Ye.adjoint()).trace() - AbsSqr(Lambdax)*(2*
         Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QSp) + 3*(Yu*
         Yu.adjoint()).trace() + 3*(Kappa*(Kappa).adjoint()).trace() + 2*(
         Lambda12*(Lambda12).adjoint()).trace()) - 9*(Yd*Yd.adjoint()*Yd*
         Yd.adjoint()).trace() - 3*(Yd*Yu.adjoint()*Yu*Yd.adjoint()).trace() -
         3*(Ye*Ye.adjoint()*Ye*Ye.adjoint()).trace());
   }

   return anomDim;
}

double CLASSNAME::get_SHuSHu() const
{
   const double twoLoop = oneOver16PiSqr * oneOver16PiSqr;
   double anomDim = 0;
   const auto QH2p = INPUT(QH2p);
   const auto Qdp = INPUT(Qdp);
   const auto QDxbarp = INPUT(QDxbarp);
   const auto QDxp = INPUT(QDxp);
   const auto Qep = INPUT(Qep);
   const auto QH1p = INPUT(QH1p);
   const auto QHpbarp = INPUT(QHpbarp);
   const auto QHpp = INPUT(QHpp);
   const auto QLp = INPUT(QLp);
   const auto QQp = INPUT(QQp);
   const auto QSp = INPUT(QSp);
   const auto Qup = INPUT(Qup);

   anomDim = oneOver16PiSqr*(AbsSqr(Lambdax) - 0.3*Sqr(g1) - 1.5*Sqr(g2)
      - 2*Sqr(gN)*Sqr(QH2p) + 3*(Yu*Yu.adjoint()).trace());

   if (get_loops() > 1) {
      anomDim += twoLoop*(2.97*Power(g1,4) + 8.25*Power(g2,4) + 16*
         Power(gN,4)*Power(QH2p,4) + 0.9*Sqr(g1)*Sqr(g2) + 3.6*Qdp*QH2p*Sqr(g1)
         *Sqr(gN) + 3.6*QDxbarp*QH2p*Sqr(g1)*Sqr(gN) - 3.6*QDxp*QH2p*Sqr(g1)*
         Sqr(gN) + 3.6*Qep*QH2p*Sqr(g1)*Sqr(gN) - 3.6*QH1p*QH2p*Sqr(g1)*Sqr(gN)
         + 1.2*QH2p*QHpbarp*Sqr(g1)*Sqr(gN) - 1.2*QH2p*QHpp*Sqr(g1)*Sqr(gN) -
         3.6*QH2p*QLp*Sqr(g1)*Sqr(gN) + 3.6*QH2p*QQp*Sqr(g1)*Sqr(gN) - 7.2*QH2p
         *Qup*Sqr(g1)*Sqr(gN) + 4.8*Sqr(g1)*Sqr(gN)*Sqr(QH2p) + 6*Sqr(g2)*Sqr(
         gN)*Sqr(QH2p) + 18*Power(gN,4)*Sqr(Qdp)*Sqr(QH2p) + 18*Power(gN,4)*Sqr
         (QDxbarp)*Sqr(QH2p) + 18*Power(gN,4)*Sqr(QDxp)*Sqr(QH2p) + 6*Power(gN,
         4)*Sqr(Qep)*Sqr(QH2p) + 12*Power(gN,4)*Sqr(QH1p)*Sqr(QH2p) + 4*Power(
         gN,4)*Sqr(QH2p)*Sqr(QHpbarp) + 4*Power(gN,4)*Sqr(QH2p)*Sqr(QHpp) + 12*
         Power(gN,4)*Sqr(QH2p)*Sqr(QLp) + 36*Power(gN,4)*Sqr(QH2p)*Sqr(QQp) + 6
         *Power(gN,4)*Sqr(QH2p)*Sqr(QSp) + 18*Power(gN,4)*Sqr(QH2p)*Sqr(Qup) -
         3*Sqr(Conj(Lambdax))*Sqr(Lambdax) + 0.4*(2*Sqr(g1) + 5*(8*Sqr(g3) + 3*
         Sqr(gN)*(-Sqr(QH2p) + Sqr(QQp) + Sqr(Qup))))*(Yu*Yu.adjoint()).trace()
         + AbsSqr(Lambdax)*(2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QH2p) + 2*Sqr(
         gN)*Sqr(QSp) - 3*(Yd*Yd.adjoint()).trace() - (Ye*Ye.adjoint()).trace()
         - 3*(Kappa*(Kappa).adjoint()).trace() - 2*(Lambda12*(Lambda12)
         .adjoint()).trace()) - 3*(Yd*Yu.adjoint()*Yu*Yd.adjoint()).trace() - 9
         *(Yu*Yu.adjoint()*Yu*Yu.adjoint()).trace());
   }

   return anomDim;
}

Eigen::Matrix<double,3,3> CLASSNAME::get_SdRSdR() const
{
   const double twoLoop = oneOver16PiSqr * oneOver16PiSqr;
   Eigen::Matrix<double,3,3> anomDim;
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

   anomDim = oneOver16PiSqr*(2*(Yd.conjugate()*Yd.transpose()) -
      0.13333333333333333*(Sqr(g1) + 20*Sqr(g3) + 15*Sqr(gN)*Sqr(Qdp))*
      UNITMATRIX(3));

   if (get_loops() > 1) {
      anomDim += twoLoop*(-2*(Yd.conjugate()*Yd.transpose()*
         Yd.conjugate()*Yd.transpose() + Yd.conjugate()*Yu.transpose()*
         Yu.conjugate()*Yd.transpose()) + Yd.conjugate()*Yd.transpose()*(-2*
         AbsSqr(Lambdax) + 0.4*Sqr(g1) + 6*Sqr(g2) - 4*Sqr(gN)*Sqr(Qdp) + 4*Sqr
         (gN)*Sqr(QH1p) + 4*Sqr(gN)*Sqr(QQp) - 6*(Yd*Yd.adjoint()).trace() - 2*
         (Ye*Ye.adjoint()).trace()) + 0.008888888888888889*(146*Power(g1,4) +
         10*Sqr(g1)*(8*Sqr(g3) + 3*Qdp*(11*Qdp + 3*(3*QDxbarp - 3*QDxp + 3*Qep
         - 3*QH1p + 3*QH2p + QHpbarp - QHpp - 3*QLp + 3*QQp - 6*Qup))*Sqr(gN))
         + 25*(32*Power(g3,4) + 48*Sqr(g3)*Sqr(gN)*Sqr(Qdp) + 9*Power(gN,4)*Sqr
         (Qdp)*(11*Sqr(Qdp) + 9*Sqr(QDxbarp) + 9*Sqr(QDxp) + 3*Sqr(Qep) + 6*Sqr
         (QH1p) + 6*Sqr(QH2p) + 2*Sqr(QHpbarp) + 2*Sqr(QHpp) + 6*Sqr(QLp) + 18*
         Sqr(QQp) + 3*Sqr(QSp) + 9*Sqr(Qup))))*UNITMATRIX(3));
   }

   return anomDim;
}

Eigen::Matrix<double,3,3> CLASSNAME::get_SuRSuR() const
{
   const double twoLoop = oneOver16PiSqr * oneOver16PiSqr;
   Eigen::Matrix<double,3,3> anomDim;
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

   anomDim = oneOver16PiSqr*(2*(Yu.conjugate()*Yu.transpose()) -
      0.13333333333333333*(4*Sqr(g1) + 20*Sqr(g3) + 15*Sqr(gN)*Sqr(Qup))*
      UNITMATRIX(3));

   if (get_loops() > 1) {
      anomDim += twoLoop*(-2*(Yu.conjugate()*Yd.transpose()*
         Yd.conjugate()*Yu.transpose() + Yu.conjugate()*Yu.transpose()*
         Yu.conjugate()*Yu.transpose()) + Yu.conjugate()*Yu.transpose()*(-2*
         AbsSqr(Lambdax) - 0.4*Sqr(g1) + 6*Sqr(g2) + 4*Sqr(gN)*Sqr(QH2p) + 4*
         Sqr(gN)*Sqr(QQp) - 4*Sqr(gN)*Sqr(Qup) - 6*(Yu*Yu.adjoint()).trace()) +
         0.008888888888888889*(608*Power(g1,4) + 20*Sqr(g1)*(16*Sqr(g3) - 3*(9
         *Qdp + 9*QDxbarp - 9*QDxp + 9*Qep - 9*QH1p + 9*QH2p + 3*QHpbarp - 3*
         QHpp - 9*QLp + 9*QQp - 22*Qup)*Qup*Sqr(gN)) + 25*(32*Power(g3,4) + 48*
         Sqr(g3)*Sqr(gN)*Sqr(Qup) + 9*Power(gN,4)*Sqr(Qup)*(9*Sqr(Qdp) + 9*Sqr(
         QDxbarp) + 9*Sqr(QDxp) + 3*Sqr(Qep) + 6*Sqr(QH1p) + 6*Sqr(QH2p) + 2*
         Sqr(QHpbarp) + 2*Sqr(QHpp) + 6*Sqr(QLp) + 18*Sqr(QQp) + 3*Sqr(QSp) +
         11*Sqr(Qup))))*UNITMATRIX(3));
   }

   return anomDim;
}

Eigen::Matrix<double,3,3> CLASSNAME::get_SeRSeR() const
{
   const double twoLoop = oneOver16PiSqr * oneOver16PiSqr;
   Eigen::Matrix<double,3,3> anomDim;
   const auto Qep = INPUT(Qep);
   const auto QH1p = INPUT(QH1p);
   const auto QLp = INPUT(QLp);
   const auto Qdp = INPUT(Qdp);
   const auto QDxbarp = INPUT(QDxbarp);
   const auto QDxp = INPUT(QDxp);
   const auto QH2p = INPUT(QH2p);
   const auto QHpbarp = INPUT(QHpbarp);
   const auto QHpp = INPUT(QHpp);
   const auto QQp = INPUT(QQp);
   const auto Qup = INPUT(Qup);
   const auto QSp = INPUT(QSp);

   anomDim = oneOver16PiSqr*(2*(Ye.conjugate()*Ye.transpose()) - 0.4*(3*
      Sqr(g1) + 5*Sqr(gN)*Sqr(Qep))*UNITMATRIX(3));

   if (get_loops() > 1) {
      anomDim += twoLoop*(-2*(Ye.conjugate()*Ye.transpose()*
         Ye.conjugate()*Ye.transpose()) + Ye.conjugate()*Ye.transpose()*(-2*
         AbsSqr(Lambdax) - 1.2*Sqr(g1) + 6*Sqr(g2) - 4*Sqr(gN)*Sqr(Qep) + 4*Sqr
         (gN)*Sqr(QH1p) + 4*Sqr(gN)*Sqr(QLp) - 6*(Yd*Yd.adjoint()).trace() - 2*
         (Ye*Ye.adjoint()).trace()) + 0.08*(162*Power(g1,4) + 30*Qep*(3*Qdp + 3
         *QDxbarp - 3*QDxp + 5*Qep - 3*QH1p + 3*QH2p + QHpbarp - QHpp - 3*QLp +
         3*QQp - 6*Qup)*Sqr(g1)*Sqr(gN) + 25*Power(gN,4)*Sqr(Qep)*(9*Sqr(Qdp)
         + 9*Sqr(QDxbarp) + 9*Sqr(QDxp) + 5*Sqr(Qep) + 6*Sqr(QH1p) + 6*Sqr(QH2p
         ) + 2*Sqr(QHpbarp) + 2*Sqr(QHpp) + 6*Sqr(QLp) + 18*Sqr(QQp) + 3*Sqr(
         QSp) + 9*Sqr(Qup)))*UNITMATRIX(3));
   }

   return anomDim;
}

double CLASSNAME::get_SsRSsR() const
{
   const double twoLoop = oneOver16PiSqr * oneOver16PiSqr;
   double anomDim = 0;
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

   anomDim = oneOver16PiSqr*(2*AbsSqr(Lambdax) - 2*Sqr(gN)*Sqr(QSp) + 3*(
      Kappa*(Kappa).adjoint()).trace() + 2*(Lambda12*(Lambda12).adjoint())
      .trace());

   if (get_loops() > 1) {
      anomDim += twoLoop*(10*Power(gN,4)*Power(QSp,4) + 18*Power(gN,4)
         *Sqr(Qdp)*Sqr(QSp) + 18*Power(gN,4)*Sqr(QDxbarp)*Sqr(QSp) + 18*Power(
         gN,4)*Sqr(QDxp)*Sqr(QSp) + 6*Power(gN,4)*Sqr(Qep)*Sqr(QSp) + 12*Power(
         gN,4)*Sqr(QH1p)*Sqr(QSp) + 12*Power(gN,4)*Sqr(QH2p)*Sqr(QSp) + 4*Power
         (gN,4)*Sqr(QHpbarp)*Sqr(QSp) + 4*Power(gN,4)*Sqr(QHpp)*Sqr(QSp) + 12*
         Power(gN,4)*Sqr(QLp)*Sqr(QSp) + 36*Power(gN,4)*Sqr(QQp)*Sqr(QSp) + 18*
         Power(gN,4)*Sqr(QSp)*Sqr(Qup) - 4*Sqr(Conj(Lambdax))*Sqr(Lambdax) +
         Conj(Lambdax)*(1.2*Lambdax*Sqr(g1) + 6*Lambdax*Sqr(g2) + 4*Lambdax*Sqr
         (gN)*Sqr(QH1p) + 4*Lambdax*Sqr(gN)*Sqr(QH2p) - 4*Lambdax*Sqr(gN)*Sqr(
         QSp) - 6*Lambdax*(Yd*Yd.adjoint()).trace() - 2*Lambdax*(Ye*Ye.adjoint(
         )).trace() - 6*Lambdax*(Yu*Yu.adjoint()).trace()) + 0.4*(2*Sqr(g1) + 5
         *(8*Sqr(g3) + 3*Sqr(gN)*(Sqr(QDxbarp) + Sqr(QDxp) - Sqr(QSp))))*(Kappa
         *(Kappa).adjoint()).trace() + 1.2*Sqr(g1)*(Lambda12*(Lambda12).adjoint
         ()).trace() + 6*Sqr(g2)*(Lambda12*(Lambda12).adjoint()).trace() + 4*
         Sqr(gN)*Sqr(QH1p)*(Lambda12*(Lambda12).adjoint()).trace() + 4*Sqr(gN)*
         Sqr(QH2p)*(Lambda12*(Lambda12).adjoint()).trace() - 4*Sqr(gN)*Sqr(QSp)
         *(Lambda12*(Lambda12).adjoint()).trace() - 6*(Kappa*(Kappa).adjoint()*
         Kappa*(Kappa).adjoint()).trace() - 4*(Lambda12*(Lambda12).adjoint()*
         Lambda12*(Lambda12).adjoint()).trace());
   }

   return anomDim;
}

Eigen::Matrix<double,2,2> CLASSNAME::get_SH1ISH1I() const
{
   const double twoLoop = oneOver16PiSqr * oneOver16PiSqr;
   Eigen::Matrix<double,2,2> anomDim;
   const auto QH1p = INPUT(QH1p);
   const auto QH2p = INPUT(QH2p);
   const auto QSp = INPUT(QSp);
   const auto Qdp = INPUT(Qdp);
   const auto QDxbarp = INPUT(QDxbarp);
   const auto QDxp = INPUT(QDxp);
   const auto Qep = INPUT(Qep);
   const auto QHpbarp = INPUT(QHpbarp);
   const auto QHpp = INPUT(QHpp);
   const auto QLp = INPUT(QLp);
   const auto QQp = INPUT(QQp);
   const auto Qup = INPUT(Qup);

   anomDim = oneOver16PiSqr*((Lambda12).adjoint()*Lambda12 - 0.1*(3*Sqr(
      g1) + 15*Sqr(g2) + 20*Sqr(gN)*Sqr(QH1p))*UNITMATRIX(2));

   if (get_loops() > 1) {
      anomDim += twoLoop*(-((Lambda12).adjoint()*Lambda12*(Lambda12)
         .adjoint()*Lambda12) - (Lambda12).adjoint()*Lambda12*(2*AbsSqr(Lambdax
         ) + 3*(Kappa*(Kappa).adjoint()).trace() + 2*(Sqr(gN)*Sqr(QH1p) - Sqr(
         gN)*Sqr(QH2p) - Sqr(gN)*Sqr(QSp) + (Lambda12*(Lambda12).adjoint())
         .trace())) + 0.01*(297*Power(g1,4) + 30*Sqr(g1)*(3*Sqr(g2) + 4*QH1p*(
         -3*Qdp - 3*QDxbarp + 3*QDxp - 3*Qep + 4*QH1p - 3*QH2p - QHpbarp + QHpp
         + 3*QLp - 3*QQp + 6*Qup)*Sqr(gN)) + 25*(33*Power(g2,4) + 24*Sqr(g2)*
         Sqr(gN)*Sqr(QH1p) + 8*Power(gN,4)*Sqr(QH1p)*(9*Sqr(Qdp) + 9*Sqr(
         QDxbarp) + 9*Sqr(QDxp) + 3*Sqr(Qep) + 8*Sqr(QH1p) + 6*Sqr(QH2p) + 2*
         Sqr(QHpbarp) + 2*Sqr(QHpp) + 6*Sqr(QLp) + 18*Sqr(QQp) + 3*Sqr(QSp) + 9
         *Sqr(Qup))))*UNITMATRIX(2));
   }

   return anomDim;
}

Eigen::Matrix<double,2,2> CLASSNAME::get_SH2ISH2I() const
{
   const double twoLoop = oneOver16PiSqr * oneOver16PiSqr;
   Eigen::Matrix<double,2,2> anomDim;
   const auto QH2p = INPUT(QH2p);
   const auto QH1p = INPUT(QH1p);
   const auto QSp = INPUT(QSp);
   const auto Qdp = INPUT(Qdp);
   const auto QDxbarp = INPUT(QDxbarp);
   const auto QDxp = INPUT(QDxp);
   const auto Qep = INPUT(Qep);
   const auto QHpbarp = INPUT(QHpbarp);
   const auto QHpp = INPUT(QHpp);
   const auto QLp = INPUT(QLp);
   const auto QQp = INPUT(QQp);
   const auto Qup = INPUT(Qup);

   anomDim = oneOver16PiSqr*(Lambda12.conjugate()*(Lambda12).transpose()
      - 0.1*(3*Sqr(g1) + 15*Sqr(g2) + 20*Sqr(gN)*Sqr(QH2p))*UNITMATRIX(2));

   if (get_loops() > 1) {
      anomDim += twoLoop*(-(Lambda12.conjugate()*(Lambda12).transpose(
         )*Lambda12.conjugate()*(Lambda12).transpose()) + Lambda12.conjugate()*
         (Lambda12).transpose()*(-2*AbsSqr(Lambdax) + 2*Sqr(gN)*Sqr(QH1p) - 2*
         Sqr(gN)*Sqr(QH2p) + 2*Sqr(gN)*Sqr(QSp) - 3*(Kappa*(Kappa).adjoint())
         .trace() - 2*(Lambda12*(Lambda12).adjoint()).trace()) + 0.01*(297*
         Power(g1,4) + 30*Sqr(g1)*(3*Sqr(g2) + 4*QH2p*(3*Qdp + 3*QDxbarp - 3*
         QDxp + 3*Qep - 3*QH1p + 4*QH2p + QHpbarp - QHpp - 3*QLp + 3*QQp - 6*
         Qup)*Sqr(gN)) + 25*(33*Power(g2,4) + 24*Sqr(g2)*Sqr(gN)*Sqr(QH2p) + 8*
         Power(gN,4)*Sqr(QH2p)*(9*Sqr(Qdp) + 9*Sqr(QDxbarp) + 9*Sqr(QDxp) + 3*
         Sqr(Qep) + 6*Sqr(QH1p) + 8*Sqr(QH2p) + 2*Sqr(QHpbarp) + 2*Sqr(QHpp) +
         6*Sqr(QLp) + 18*Sqr(QQp) + 3*Sqr(QSp) + 9*Sqr(Qup))))*UNITMATRIX(2));
   }

   return anomDim;
}

Eigen::Matrix<double,2,2> CLASSNAME::get_SsIRSsIR() const
{
   const double twoLoop = oneOver16PiSqr * oneOver16PiSqr;
   Eigen::Matrix<double,2,2> anomDim;
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

   anomDim = -2*oneOver16PiSqr*Sqr(gN)*Sqr(QSp)*UNITMATRIX(2);

   if (get_loops() > 1) {
      anomDim += 2*Power(gN,4)*twoLoop*Sqr(QSp)*(9*Sqr(Qdp) + 9*Sqr(
         QDxbarp) + 9*Sqr(QDxp) + 3*Sqr(Qep) + 6*Sqr(QH1p) + 6*Sqr(QH2p) + 2*
         Sqr(QHpbarp) + 2*Sqr(QHpp) + 6*Sqr(QLp) + 18*Sqr(QQp) + 5*Sqr(QSp) + 9
         *Sqr(Qup))*UNITMATRIX(2);
   }

   return anomDim;
}

Eigen::Matrix<double,3,3> CLASSNAME::get_SDxLSDxL() const
{
   const double twoLoop = oneOver16PiSqr * oneOver16PiSqr;
   Eigen::Matrix<double,3,3> anomDim;
   const auto QDxp = INPUT(QDxp);
   const auto QDxbarp = INPUT(QDxbarp);
   const auto QSp = INPUT(QSp);
   const auto Qdp = INPUT(Qdp);
   const auto Qep = INPUT(Qep);
   const auto QH1p = INPUT(QH1p);
   const auto QH2p = INPUT(QH2p);
   const auto QHpbarp = INPUT(QHpbarp);
   const auto QHpp = INPUT(QHpp);
   const auto QLp = INPUT(QLp);
   const auto QQp = INPUT(QQp);
   const auto Qup = INPUT(Qup);

   anomDim = oneOver16PiSqr*(Kappa.conjugate()*(Kappa).transpose() -
      0.13333333333333333*(Sqr(g1) + 20*Sqr(g3) + 15*Sqr(gN)*Sqr(QDxp))*
      UNITMATRIX(3));

   if (get_loops() > 1) {
      anomDim += twoLoop*(-(Kappa.conjugate()*(Kappa).transpose()*
         Kappa.conjugate()*(Kappa).transpose()) + Kappa.conjugate()*(Kappa)
         .transpose()*(-2*AbsSqr(Lambdax) + 2*Sqr(gN)*Sqr(QDxbarp) - 2*Sqr(gN)*
         Sqr(QDxp) + 2*Sqr(gN)*Sqr(QSp) - 3*(Kappa*(Kappa).adjoint()).trace() -
         2*(Lambda12*(Lambda12).adjoint()).trace()) + 0.008888888888888889*(
         146*Power(g1,4) + 10*Sqr(g1)*(8*Sqr(g3) + 3*QDxp*(-9*Qdp - 9*QDxbarp +
         11*QDxp - 9*Qep + 9*QH1p - 9*QH2p - 3*QHpbarp + 3*QHpp + 9*QLp - 9*
         QQp + 18*Qup)*Sqr(gN)) + 25*(32*Power(g3,4) + 48*Sqr(g3)*Sqr(gN)*Sqr(
         QDxp) + 9*Power(gN,4)*Sqr(QDxp)*(9*Sqr(Qdp) + 9*Sqr(QDxbarp) + 11*Sqr(
         QDxp) + 3*Sqr(Qep) + 6*Sqr(QH1p) + 6*Sqr(QH2p) + 2*Sqr(QHpbarp) + 2*
         Sqr(QHpp) + 6*Sqr(QLp) + 18*Sqr(QQp) + 3*Sqr(QSp) + 9*Sqr(Qup))))*
         UNITMATRIX(3));
   }

   return anomDim;
}

Eigen::Matrix<double,3,3> CLASSNAME::get_SDxbarRSDxbarR() const
{
   const double twoLoop = oneOver16PiSqr * oneOver16PiSqr;
   Eigen::Matrix<double,3,3> anomDim;
   const auto QDxbarp = INPUT(QDxbarp);
   const auto QDxp = INPUT(QDxp);
   const auto QSp = INPUT(QSp);
   const auto Qdp = INPUT(Qdp);
   const auto Qep = INPUT(Qep);
   const auto QH1p = INPUT(QH1p);
   const auto QH2p = INPUT(QH2p);
   const auto QHpbarp = INPUT(QHpbarp);
   const auto QHpp = INPUT(QHpp);
   const auto QLp = INPUT(QLp);
   const auto QQp = INPUT(QQp);
   const auto Qup = INPUT(Qup);

   anomDim = oneOver16PiSqr*((Kappa).adjoint()*Kappa -
      0.13333333333333333*(Sqr(g1) + 20*Sqr(g3) + 15*Sqr(gN)*Sqr(QDxbarp))*
      UNITMATRIX(3));

   if (get_loops() > 1) {
      anomDim += twoLoop*(-((Kappa).adjoint()*Kappa*(Kappa).adjoint()*
         Kappa) - (Kappa).adjoint()*Kappa*(2*AbsSqr(Lambdax) + 3*(Kappa*(Kappa)
         .adjoint()).trace() + 2*(Sqr(gN)*Sqr(QDxbarp) - Sqr(gN)*Sqr(QDxp) -
         Sqr(gN)*Sqr(QSp) + (Lambda12*(Lambda12).adjoint()).trace())) +
         0.008888888888888889*(146*Power(g1,4) + 10*Sqr(g1)*(8*Sqr(g3) + 3*
         QDxbarp*(9*Qdp + 11*QDxbarp - 9*QDxp + 9*Qep - 9*QH1p + 9*QH2p + 3*
         QHpbarp - 3*QHpp - 9*QLp + 9*QQp - 18*Qup)*Sqr(gN)) + 25*(32*Power(g3,
         4) + 48*Sqr(g3)*Sqr(gN)*Sqr(QDxbarp) + 9*Power(gN,4)*Sqr(QDxbarp)*(9*
         Sqr(Qdp) + 11*Sqr(QDxbarp) + 9*Sqr(QDxp) + 3*Sqr(Qep) + 6*Sqr(QH1p) +
         6*Sqr(QH2p) + 2*Sqr(QHpbarp) + 2*Sqr(QHpp) + 6*Sqr(QLp) + 18*Sqr(QQp)
         + 3*Sqr(QSp) + 9*Sqr(Qup))))*UNITMATRIX(3));
   }

   return anomDim;
}

double CLASSNAME::get_SHpSHp() const
{
   const double twoLoop = oneOver16PiSqr * oneOver16PiSqr;
   double anomDim = 0;
   const auto QHpp = INPUT(QHpp);
   const auto Qdp = INPUT(Qdp);
   const auto QDxbarp = INPUT(QDxbarp);
   const auto QDxp = INPUT(QDxp);
   const auto Qep = INPUT(Qep);
   const auto QH1p = INPUT(QH1p);
   const auto QH2p = INPUT(QH2p);
   const auto QHpbarp = INPUT(QHpbarp);
   const auto QLp = INPUT(QLp);
   const auto QQp = INPUT(QQp);
   const auto Qup = INPUT(Qup);
   const auto QSp = INPUT(QSp);

   anomDim = oneOver16PiSqr*(-0.3*Sqr(g1) - 1.5*Sqr(g2) - 2*Sqr(gN)*Sqr(
      QHpp));

   if (get_loops() > 1) {
      anomDim += twoLoop*(2.97*Power(g1,4) + 8.25*Power(g2,4) + 0.3*
         Sqr(g1)*(3*Sqr(g2) - 4*QHpp*(3*Qdp + 3*QDxbarp - 3*QDxp + 3*Qep - 3*
         QH1p + 3*QH2p + QHpbarp - 2*QHpp - 3*QLp + 3*QQp - 6*Qup)*Sqr(gN)) + 6
         *Sqr(g2)*Sqr(gN)*Sqr(QHpp) + 2*Power(gN,4)*Sqr(QHpp)*(9*Sqr(Qdp) + 9*
         Sqr(QDxbarp) + 9*Sqr(QDxp) + 3*Sqr(Qep) + 6*Sqr(QH1p) + 6*Sqr(QH2p) +
         2*Sqr(QHpbarp) + 4*Sqr(QHpp) + 6*Sqr(QLp) + 18*Sqr(QQp) + 3*Sqr(QSp) +
         9*Sqr(Qup)));
   }

   return anomDim;
}

double CLASSNAME::get_SHpbarSHpbar() const
{
   const double twoLoop = oneOver16PiSqr * oneOver16PiSqr;
   double anomDim = 0;
   const auto QHpbarp = INPUT(QHpbarp);
   const auto Qdp = INPUT(Qdp);
   const auto QDxbarp = INPUT(QDxbarp);
   const auto QDxp = INPUT(QDxp);
   const auto Qep = INPUT(Qep);
   const auto QH1p = INPUT(QH1p);
   const auto QH2p = INPUT(QH2p);
   const auto QHpp = INPUT(QHpp);
   const auto QLp = INPUT(QLp);
   const auto QQp = INPUT(QQp);
   const auto Qup = INPUT(Qup);
   const auto QSp = INPUT(QSp);

   anomDim = oneOver16PiSqr*(-0.3*Sqr(g1) - 1.5*Sqr(g2) - 2*Sqr(gN)*Sqr(
      QHpbarp));

   if (get_loops() > 1) {
      anomDim += twoLoop*(2.97*Power(g1,4) + 8.25*Power(g2,4) + 0.3*
         Sqr(g1)*(3*Sqr(g2) + 4*QHpbarp*(3*Qdp + 3*QDxbarp - 3*QDxp + 3*Qep - 3
         *QH1p + 3*QH2p + 2*QHpbarp - QHpp - 3*QLp + 3*QQp - 6*Qup)*Sqr(gN)) +
         6*Sqr(g2)*Sqr(gN)*Sqr(QHpbarp) + 2*Power(gN,4)*Sqr(QHpbarp)*(9*Sqr(Qdp
         ) + 9*Sqr(QDxbarp) + 9*Sqr(QDxp) + 3*Sqr(Qep) + 6*Sqr(QH1p) + 6*Sqr(
         QH2p) + 4*Sqr(QHpbarp) + 2*Sqr(QHpp) + 6*Sqr(QLp) + 18*Sqr(QQp) + 3*
         Sqr(QSp) + 9*Sqr(Qup)));
   }

   return anomDim;
}


const Eigen::ArrayXd genericE6SSM_susy_parameters::get() const
{
   Eigen::ArrayXd pars(numberOfParameters);

   pars(0) = Yd(0,0);
   pars(1) = Yd(0,1);
   pars(2) = Yd(0,2);
   pars(3) = Yd(1,0);
   pars(4) = Yd(1,1);
   pars(5) = Yd(1,2);
   pars(6) = Yd(2,0);
   pars(7) = Yd(2,1);
   pars(8) = Yd(2,2);
   pars(9) = Ye(0,0);
   pars(10) = Ye(0,1);
   pars(11) = Ye(0,2);
   pars(12) = Ye(1,0);
   pars(13) = Ye(1,1);
   pars(14) = Ye(1,2);
   pars(15) = Ye(2,0);
   pars(16) = Ye(2,1);
   pars(17) = Ye(2,2);
   pars(18) = Kappa(0,0);
   pars(19) = Kappa(0,1);
   pars(20) = Kappa(0,2);
   pars(21) = Kappa(1,0);
   pars(22) = Kappa(1,1);
   pars(23) = Kappa(1,2);
   pars(24) = Kappa(2,0);
   pars(25) = Kappa(2,1);
   pars(26) = Kappa(2,2);
   pars(27) = Lambda12(0,0);
   pars(28) = Lambda12(0,1);
   pars(29) = Lambda12(1,0);
   pars(30) = Lambda12(1,1);
   pars(31) = Lambdax;
   pars(32) = Yu(0,0);
   pars(33) = Yu(0,1);
   pars(34) = Yu(0,2);
   pars(35) = Yu(1,0);
   pars(36) = Yu(1,1);
   pars(37) = Yu(1,2);
   pars(38) = Yu(2,0);
   pars(39) = Yu(2,1);
   pars(40) = Yu(2,2);
   pars(41) = MuPr;
   pars(42) = g1;
   pars(43) = g2;
   pars(44) = g3;
   pars(45) = gN;
   pars(46) = vd;
   pars(47) = vu;
   pars(48) = vs;


   return pars;
}

void genericE6SSM_susy_parameters::print(std::ostream& ostr) const
{
   ostr << "susy parameters:\n";
   ostr << "Yd = " << Yd << '\n';
   ostr << "Ye = " << Ye << '\n';
   ostr << "Kappa = " << Kappa << '\n';
   ostr << "Lambda12 = " << Lambda12 << '\n';
   ostr << "Lambdax = " << Lambdax << '\n';
   ostr << "Yu = " << Yu << '\n';
   ostr << "MuPr = " << MuPr << '\n';
   ostr << "g1 = " << g1 << '\n';
   ostr << "g2 = " << g2 << '\n';
   ostr << "g3 = " << g3 << '\n';
   ostr << "gN = " << gN << '\n';
   ostr << "vd = " << vd << '\n';
   ostr << "vu = " << vu << '\n';
   ostr << "vs = " << vs << '\n';

}

void genericE6SSM_susy_parameters::set(const Eigen::ArrayXd& v)
{
   Yd(0,0) = v(0);
   Yd(0,1) = v(1);
   Yd(0,2) = v(2);
   Yd(1,0) = v(3);
   Yd(1,1) = v(4);
   Yd(1,2) = v(5);
   Yd(2,0) = v(6);
   Yd(2,1) = v(7);
   Yd(2,2) = v(8);
   Ye(0,0) = v(9);
   Ye(0,1) = v(10);
   Ye(0,2) = v(11);
   Ye(1,0) = v(12);
   Ye(1,1) = v(13);
   Ye(1,2) = v(14);
   Ye(2,0) = v(15);
   Ye(2,1) = v(16);
   Ye(2,2) = v(17);
   Kappa(0,0) = v(18);
   Kappa(0,1) = v(19);
   Kappa(0,2) = v(20);
   Kappa(1,0) = v(21);
   Kappa(1,1) = v(22);
   Kappa(1,2) = v(23);
   Kappa(2,0) = v(24);
   Kappa(2,1) = v(25);
   Kappa(2,2) = v(26);
   Lambda12(0,0) = v(27);
   Lambda12(0,1) = v(28);
   Lambda12(1,0) = v(29);
   Lambda12(1,1) = v(30);
   Lambdax = v(31);
   Yu(0,0) = v(32);
   Yu(0,1) = v(33);
   Yu(0,2) = v(34);
   Yu(1,0) = v(35);
   Yu(1,1) = v(36);
   Yu(1,2) = v(37);
   Yu(2,0) = v(38);
   Yu(2,1) = v(39);
   Yu(2,2) = v(40);
   MuPr = v(41);
   g1 = v(42);
   g2 = v(43);
   g3 = v(44);
   gN = v(45);
   vd = v(46);
   vu = v(47);
   vs = v(48);

}

const genericE6SSM_input_parameters& genericE6SSM_susy_parameters::get_input() const
{
   return input;
}

void genericE6SSM_susy_parameters::set_input_parameters(const genericE6SSM_input_parameters& input_)
{
   input = input_;
}

void genericE6SSM_susy_parameters::calc_susy_traces(Susy_traces& susy_traces) const
{
   TRACE_STRUCT.traceYdAdjYd = (Yd*Yd.adjoint()).trace();
   TRACE_STRUCT.traceYeAdjYe = (Ye*Ye.adjoint()).trace();
   TRACE_STRUCT.traceYuAdjYu = (Yu*Yu.adjoint()).trace();
   TRACE_STRUCT.traceKappaAdjKappa = (Kappa*(Kappa).adjoint()).trace();
   TRACE_STRUCT.traceLambda12AdjLambda12 = (Lambda12*(Lambda12).adjoint())
      .trace();
   TRACE_STRUCT.traceYdAdjYdYdAdjYd = (Yd*Yd.adjoint()*Yd*Yd.adjoint()).trace()
      ;
   TRACE_STRUCT.traceYdAdjYuYuAdjYd = (Yd*Yu.adjoint()*Yu*Yd.adjoint()).trace()
      ;
   TRACE_STRUCT.traceYeAdjYeYeAdjYe = (Ye*Ye.adjoint()*Ye*Ye.adjoint()).trace()
      ;
   TRACE_STRUCT.traceKappaAdjKappaKappaAdjKappa = (Kappa*(Kappa).adjoint()*
      Kappa*(Kappa).adjoint()).trace();
   TRACE_STRUCT.traceLambda12AdjLambda12Lambda12AdjLambda12 = (Lambda12*(
      Lambda12).adjoint()*Lambda12*(Lambda12).adjoint()).trace();
   TRACE_STRUCT.traceYuAdjYuYuAdjYu = (Yu*Yu.adjoint()*Yu*Yu.adjoint()).trace()
      ;

}

std::ostream& operator<<(std::ostream& ostr, const genericE6SSM_susy_parameters& susy_pars)
{
   susy_pars.print(std::cout);
   return ostr;
}

} // namespace flexiblesusy
