#include "lowE6SSM_two_scale_tuning_calculator.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

double lowE6SSM_tuning_calculator::leading_log_coefficient_Yu22() const
{
   const lowE6SSM_input_parameters inputs = model.get_input();

   const auto QQp = inputs.QQp;
   const auto QLp = inputs.QLp;
   const auto QH1p = inputs.QH1p;
   const auto QH2p = inputs.QH2p;
   const auto Qdp = inputs.Qdp;
   const auto Qup = inputs.Qup;
   const auto Qep = inputs.Qep;
   const auto QSp = inputs.QSp;
   const auto QDxp = inputs.QDxp;
   const auto QDxbarp = inputs.QDxbarp;
   const auto QHpp = inputs.QHpp;
   const auto QHpbarp = inputs.QHpbarp;
   
   const double Yd00 = model.get_Yd(0,0);
   const double Yd01 = model.get_Yd(0,1);
   const double Yd02 = model.get_Yd(0,2);
   const double Yd10 = model.get_Yd(1,0);
   const double Yd11 = model.get_Yd(1,1);
   const double Yd12 = model.get_Yd(1,2);
   const double Yd20 = model.get_Yd(2,0);
   const double Yd21 = model.get_Yd(2,1);
   const double Yd22 = model.get_Yd(2,2);
   const double Ye00 = model.get_Ye(0,0);
   const double Ye01 = model.get_Ye(0,1);
   const double Ye02 = model.get_Ye(0,2);
   const double Ye10 = model.get_Ye(1,0);
   const double Ye11 = model.get_Ye(1,1);
   const double Ye12 = model.get_Ye(1,2);
   const double Ye20 = model.get_Ye(2,0);
   const double Ye21 = model.get_Ye(2,1);
   const double Ye22 = model.get_Ye(2,2);
   const double Kappa00 = model.get_Kappa(0,0);
   const double Kappa01 = model.get_Kappa(0,1);
   const double Kappa02 = model.get_Kappa(0,2);
   const double Kappa10 = model.get_Kappa(1,0);
   const double Kappa11 = model.get_Kappa(1,1);
   const double Kappa12 = model.get_Kappa(1,2);
   const double Kappa20 = model.get_Kappa(2,0);
   const double Kappa21 = model.get_Kappa(2,1);
   const double Kappa22 = model.get_Kappa(2,2);
   const double Lambda1200 = model.get_Lambda12(0,0);
   const double Lambda1201 = model.get_Lambda12(0,1);
   const double Lambda1210 = model.get_Lambda12(1,0);
   const double Lambda1211 = model.get_Lambda12(1,1);
   const double Lambdax = model.get_Lambdax();
   const double Yu00 = model.get_Yu(0,0);
   const double Yu01 = model.get_Yu(0,1);
   const double Yu02 = model.get_Yu(0,2);
   const double Yu10 = model.get_Yu(1,0);
   const double Yu11 = model.get_Yu(1,1);
   const double Yu12 = model.get_Yu(1,2);
   const double Yu20 = model.get_Yu(2,0);
   const double Yu21 = model.get_Yu(2,1);
   const double Yu22 = model.get_Yu(2,2);
   const double g1 = model.get_g1();
   const double g2 = model.get_g2();
   const double g3 = model.get_g3();
   const double gN = model.get_gN();

   double coeff = -16.64*Power(g1,4)*Yu22 - 24*Power(g2,4)*Yu22 + Power(gN,3)*
      Yu22*(9*Sqr(Qdp) + 9*Sqr(QDxbarp) + 9*Sqr(QDxp) + 3*Sqr(Qep) + 6*Sqr(QH1p) +
      6*Sqr(QH2p) + 2*Sqr(QHpbarp) + 2*Sqr(QHpp) + 6*Sqr(QLp) + 18*Sqr(QQp) + 3*
      Sqr(QSp) + 9*Sqr(Qup))*(-4*gN*Sqr(QH2p) - 4*gN*Sqr(QQp) - 4*gN*Sqr(Qup)) +
      Yd02*Yu20*(Yu00*(Yd00*Yu00 + Yd01*Yu01 + Yd02*Yu02) + Yu10*(Yd00*Yu10 + Yd01
      *Yu11 + Yd02*Yu12) + Yu20*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22) + 3*(Yd10*(
      Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + Yd20*(Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22
      ) + Yd00*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02))) + Yd00*(-0.4666666666666667*
      Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*
      Sqr(Qdp) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QQp) + 3*(Sqr(Yd00) + Sqr(
      Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21
      ) + Sqr(Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) +
      Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22))) + Yd02*Yu21*(Yu01*(Yd00*
      Yu00 + Yd01*Yu01 + Yd02*Yu02) + Yu11*(Yd00*Yu10 + Yd01*Yu11 + Yd02*Yu12) +
      Yu21*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22) + 3*(Yd11*(Yd00*Yd10 + Yd01*Yd11 +
      Yd02*Yd12) + Yd21*(Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22) + Yd01*(Sqr(Yd00) +
      Sqr(Yd01) + Sqr(Yd02))) + Yd01*(-0.4666666666666667*Sqr(g1) - 3*Sqr(g2) -
      5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(Qdp) - 2*Sqr(gN)*
      Sqr(QH1p) - 2*Sqr(gN)*Sqr(QQp) + 3*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(
      Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(
      Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20
      ) + Sqr(Ye21) + Sqr(Ye22))) + (Yd00*Yu20 + Yd01*Yu21 + 2*Yd02*Yu22)*(Yu02*(
      Yd00*Yu00 + Yd01*Yu01 + Yd02*Yu02) + Yu12*(Yd00*Yu10 + Yd01*Yu11 + Yd02*Yu12
      ) + Yu22*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22) + 3*(Yd12*(Yd00*Yd10 + Yd01*
      Yd11 + Yd02*Yd12) + Yd22*(Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22) + Yd02*(Sqr(
      Yd00) + Sqr(Yd01) + Sqr(Yd02))) + Yd02*(-0.4666666666666667*Sqr(g1) - 3*Sqr(
      g2) - 5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(Qdp) - 2*Sqr(
      gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QQp) + 3*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) +
      Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr
      (Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(
      Ye20) + Sqr(Ye21) + Sqr(Ye22))) + Yd12*Yu20*(Yu00*(Yd10*Yu00 + Yd11*Yu01 +
      Yd12*Yu02) + Yu10*(Yd10*Yu10 + Yd11*Yu11 + Yd12*Yu12) + Yu20*(Yd10*Yu20 +
      Yd11*Yu21 + Yd12*Yu22) + 3*(Yd00*(Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + Yd20*
      (Yd10*Yd20 + Yd11*Yd21 + Yd12*Yd22) + Yd10*(Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12
      ))) + Yd10*(-0.4666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(
      g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(Qdp) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*
      Sqr(QQp) + 3*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) +
      Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr
      (Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(
      Ye22))) + Yd12*Yu21*(Yu01*(Yd10*Yu00 + Yd11*Yu01 + Yd12*Yu02) + Yu11*(Yd10*
      Yu10 + Yd11*Yu11 + Yd12*Yu12) + Yu21*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22) + 3
      *(Yd01*(Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + Yd21*(Yd10*Yd20 + Yd11*Yd21 +
      Yd12*Yd22) + Yd11*(Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12))) + Yd11*(
      -0.4666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(
      Lambdax) - 2*Sqr(gN)*Sqr(Qdp) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QQp) + 3
      *(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) +
      Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr
      (Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22))) + (Yd10
      *Yu20 + Yd11*Yu21 + 2*Yd12*Yu22)*(Yu02*(Yd10*Yu00 + Yd11*Yu01 + Yd12*Yu02) +
      Yu12*(Yd10*Yu10 + Yd11*Yu11 + Yd12*Yu12) + Yu22*(Yd10*Yu20 + Yd11*Yu21 +
      Yd12*Yu22) + 3*(Yd02*(Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + Yd22*(Yd10*Yd20 +
      Yd11*Yd21 + Yd12*Yd22) + Yd12*(Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12))) + Yd12*(
      -0.4666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(
      Lambdax) - 2*Sqr(gN)*Sqr(Qdp) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QQp) + 3
      *(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) +
      Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr
      (Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22))) + Yd22*
      Yu20*(Yu00*(Yd20*Yu00 + Yd21*Yu01 + Yd22*Yu02) + Yu10*(Yd20*Yu10 + Yd21*Yu11
      + Yd22*Yu12) + Yu20*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22) + 3*(Yd00*(Yd00*
      Yd20 + Yd01*Yd21 + Yd02*Yd22) + Yd10*(Yd10*Yd20 + Yd11*Yd21 + Yd12*Yd22) +
      Yd20*(Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22))) + Yd20*(-0.4666666666666667*Sqr(g1
      ) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(Qdp
      ) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QQp) + 3*(Sqr(Yd00) + Sqr(Yd01) +
      Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(
      Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(
      Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22))) + Yd22*Yu21*(Yu01*(Yd20*Yu00 +
      Yd21*Yu01 + Yd22*Yu02) + Yu11*(Yd20*Yu10 + Yd21*Yu11 + Yd22*Yu12) + Yu21*(
      Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22) + 3*(Yd01*(Yd00*Yd20 + Yd01*Yd21 + Yd02*
      Yd22) + Yd11*(Yd10*Yd20 + Yd11*Yd21 + Yd12*Yd22) + Yd21*(Sqr(Yd20) + Sqr(
      Yd21) + Sqr(Yd22))) + Yd21*(-0.4666666666666667*Sqr(g1) - 3*Sqr(g2) -
      5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(Qdp) - 2*Sqr(gN)*
      Sqr(QH1p) - 2*Sqr(gN)*Sqr(QQp) + 3*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(
      Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(
      Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20
      ) + Sqr(Ye21) + Sqr(Ye22))) + (Yd20*Yu20 + Yd21*Yu21 + 2*Yd22*Yu22)*(Yu02*(
      Yd20*Yu00 + Yd21*Yu01 + Yd22*Yu02) + Yu12*(Yd20*Yu10 + Yd21*Yu11 + Yd22*Yu12
      ) + Yu22*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22) + 3*(Yd02*(Yd00*Yd20 + Yd01*
      Yd21 + Yd02*Yd22) + Yd12*(Yd10*Yd20 + Yd11*Yd21 + Yd12*Yd22) + Yd22*(Sqr(
      Yd20) + Sqr(Yd21) + Sqr(Yd22))) + Yd22*(-0.4666666666666667*Sqr(g1) - 3*Sqr(
      g2) - 5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(Qdp) - 2*Sqr(
      gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QQp) + 3*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) +
      Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr
      (Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(
      Ye20) + Sqr(Ye21) + Sqr(Ye22))) + 2*Lambdax*Yu22*(4*Power(Lambdax,3) - 0.6*
      Lambdax*Sqr(g1) - 3*Lambdax*Sqr(g2) + 3*Lambdax*(Sqr(Kappa00) + Sqr(Kappa01)
      + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20)
      + Sqr(Kappa21) + Sqr(Kappa22)) + 2*Lambdax*(Sqr(Lambda1200) + Sqr(Lambda1201
      ) + Sqr(Lambda1210) + Sqr(Lambda1211)) - 2*Lambdax*Sqr(gN)*Sqr(QH1p) - 2*
      Lambdax*Sqr(gN)*Sqr(QH2p) - 2*Lambdax*Sqr(gN)*Sqr(QSp) + 3*Lambdax*(Sqr(Yd00
      ) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) +
      Sqr(Yd21) + Sqr(Yd22)) + Lambdax*(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(
      Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22)) + 3*
      Lambdax*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(
      Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + (3*Yu02*Yu20 + 6*Yu00*Yu22)*(
      Yd00*(Yd00*Yu00 + Yd01*Yu01 + Yd02*Yu02) + Yd10*(Yd10*Yu00 + Yd11*Yu01 +
      Yd12*Yu02) + Yd20*(Yd20*Yu00 + Yd21*Yu01 + Yd22*Yu02) + 3*(Yu10*(Yu00*Yu10 +
      Yu01*Yu11 + Yu02*Yu12) + Yu20*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu00*(
      Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02))) + Yu00*(-0.8666666666666667*Sqr(g1) - 3*
      Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2
      *Sqr(gN)*Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02
      ) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))))
      + (3*Yu02*Yu21 + 6*Yu01*Yu22)*(Yd01*(Yd00*Yu00 + Yd01*Yu01 + Yd02*Yu02) +
      Yd11*(Yd10*Yu00 + Yd11*Yu01 + Yd12*Yu02) + Yd21*(Yd20*Yu00 + Yd21*Yu01 +
      Yd22*Yu02) + 3*(Yu11*(Yu00*Yu10 + Yu01*Yu11 + Yu02*Yu12) + Yu21*(Yu00*Yu20 +
      Yu01*Yu21 + Yu02*Yu22) + Yu01*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02))) + Yu01*(
      -0.8666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(
      Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3
      *(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) +
      Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) + (6*Yu02*Yu22 + 3*(Yu00*Yu20 + Yu01*
      Yu21 + 2*Yu02*Yu22))*(Yd02*(Yd00*Yu00 + Yd01*Yu01 + Yd02*Yu02) + Yd12*(Yd10*
      Yu00 + Yd11*Yu01 + Yd12*Yu02) + Yd22*(Yd20*Yu00 + Yd21*Yu01 + Yd22*Yu02) + 3
      *(Yu12*(Yu00*Yu10 + Yu01*Yu11 + Yu02*Yu12) + Yu22*(Yu00*Yu20 + Yu01*Yu21 +
      Yu02*Yu22) + Yu02*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02))) + Yu02*(
      -0.8666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(
      Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3
      *(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) +
      Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) + (3*Yu12*Yu20 + 6*Yu10*Yu22)*(Yd00*(
      Yd00*Yu10 + Yd01*Yu11 + Yd02*Yu12) + Yd10*(Yd10*Yu10 + Yd11*Yu11 + Yd12*Yu12
      ) + Yd20*(Yd20*Yu10 + Yd21*Yu11 + Yd22*Yu12) + 3*(Yu00*(Yu00*Yu10 + Yu01*
      Yu11 + Yu02*Yu12) + Yu20*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + Yu10*(Sqr(
      Yu10) + Sqr(Yu11) + Sqr(Yu12))) + Yu10*(-0.8666666666666667*Sqr(g1) - 3*Sqr(
      g2) - 5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr
      (gN)*Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) +
      Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) + (
      3*Yu12*Yu21 + 6*Yu11*Yu22)*(Yd01*(Yd00*Yu10 + Yd01*Yu11 + Yd02*Yu12) + Yd11*
      (Yd10*Yu10 + Yd11*Yu11 + Yd12*Yu12) + Yd21*(Yd20*Yu10 + Yd21*Yu11 + Yd22*
      Yu12) + 3*(Yu01*(Yu00*Yu10 + Yu01*Yu11 + Yu02*Yu12) + Yu21*(Yu10*Yu20 + Yu11
      *Yu21 + Yu12*Yu22) + Yu11*(Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12))) + Yu11*(
      -0.8666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(
      Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3
      *(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) +
      Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) + (6*Yu12*Yu22 + 3*(Yu10*Yu20 + Yu11*
      Yu21 + 2*Yu12*Yu22))*(Yd02*(Yd00*Yu10 + Yd01*Yu11 + Yd02*Yu12) + Yd12*(Yd10*
      Yu10 + Yd11*Yu11 + Yd12*Yu12) + Yd22*(Yd20*Yu10 + Yd21*Yu11 + Yd22*Yu12) + 3
      *(Yu02*(Yu00*Yu10 + Yu01*Yu11 + Yu02*Yu12) + Yu22*(Yu10*Yu20 + Yu11*Yu21 +
      Yu12*Yu22) + Yu12*(Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12))) + Yu12*(
      -0.8666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(
      Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3
      *(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) +
      Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) + (Yd00*Yd02 + Yd10*Yd12 + Yd20*Yd22 +
      6*Yu20*Yu22 + 3*(Yu00*Yu02 + Yu10*Yu12 + 2*Yu20*Yu22))*(Yd00*(Yd00*Yu20 +
      Yd01*Yu21 + Yd02*Yu22) + Yd10*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22) + Yd20*(
      Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22) + 3*(Yu00*(Yu00*Yu20 + Yu01*Yu21 + Yu02*
      Yu22) + Yu10*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + Yu20*(Sqr(Yu20) + Sqr(
      Yu21) + Sqr(Yu22))) + Yu20*(-0.8666666666666667*Sqr(g1) - 3*Sqr(g2) -
      5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*
      Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(
      Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) + (Yd01
      *Yd02 + Yd11*Yd12 + Yd21*Yd22 + 6*Yu21*Yu22 + 3*(Yu01*Yu02 + Yu11*Yu12 + 2*
      Yu21*Yu22))*(Yd01*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22) + Yd11*(Yd10*Yu20 +
      Yd11*Yu21 + Yd12*Yu22) + Yd21*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22) + 3*(Yu01*
      (Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu11*(Yu10*Yu20 + Yu11*Yu21 + Yu12*
      Yu22) + Yu21*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + Yu21*(
      -0.8666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(
      Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3
      *(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) +
      Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) + (-0.8666666666666667*Sqr(g1) - 3*Sqr(
      g2) - 5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr
      (gN)*Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + Sqr(Yd02) + Sqr(Yd12) + Sqr(Yd22) + 6*
      Sqr(Yu22) + 3*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) +
      Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)) + 3*(Sqr(Yu02) + Sqr(Yu12) +
      Sqr(Yu20) + Sqr(Yu21) + 3*Sqr(Yu22)))*(Yd02*(Yd00*Yu20 + Yd01*Yu21 + Yd02*
      Yu22) + Yd12*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22) + Yd22*(Yd20*Yu20 + Yd21*
      Yu21 + Yd22*Yu22) + 3*(Yu02*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu12*(Yu10
      *Yu20 + Yu11*Yu21 + Yu12*Yu22) + Yu22*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) +
      Yu22*(-0.8666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) +
      Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup
      ) + 3*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12)
      + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))));


   return twoLoop*coeff;
}

} // namespace flexiblesusy
