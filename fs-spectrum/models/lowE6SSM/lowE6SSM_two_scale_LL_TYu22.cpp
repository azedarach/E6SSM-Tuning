#include "lowE6SSM_two_scale_tuning_calculator.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

double lowE6SSM_tuning_calculator::leading_log_coefficient_TYu22() const
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

   const double TYd00 = model.get_TYd(0,0);
   const double TYd01 = model.get_TYd(0,1);
   const double TYd02 = model.get_TYd(0,2);
   const double TYd10 = model.get_TYd(1,0);
   const double TYd11 = model.get_TYd(1,1);
   const double TYd12 = model.get_TYd(1,2);
   const double TYd20 = model.get_TYd(2,0);
   const double TYd21 = model.get_TYd(2,1);
   const double TYd22 = model.get_TYd(2,2);
   const double TYe00 = model.get_TYe(0,0);
   const double TYe01 = model.get_TYe(0,1);
   const double TYe02 = model.get_TYe(0,2);
   const double TYe10 = model.get_TYe(1,0);
   const double TYe11 = model.get_TYe(1,1);
   const double TYe12 = model.get_TYe(1,2);
   const double TYe20 = model.get_TYe(2,0);
   const double TYe21 = model.get_TYe(2,1);
   const double TYe22 = model.get_TYe(2,2);
   const double TKappa00 = model.get_TKappa(0,0);
   const double TKappa01 = model.get_TKappa(0,1);
   const double TKappa02 = model.get_TKappa(0,2);
   const double TKappa10 = model.get_TKappa(1,0);
   const double TKappa11 = model.get_TKappa(1,1);
   const double TKappa12 = model.get_TKappa(1,2);
   const double TKappa20 = model.get_TKappa(2,0);
   const double TKappa21 = model.get_TKappa(2,1);
   const double TKappa22 = model.get_TKappa(2,2);
   const double TLambda1200 = model.get_TLambda12(0,0);
   const double TLambda1201 = model.get_TLambda12(0,1);
   const double TLambda1210 = model.get_TLambda12(1,0);
   const double TLambda1211 = model.get_TLambda12(1,1);
   const double TLambdax = model.get_TLambdax();
   const double TYu00 = model.get_TYu(0,0);
   const double TYu01 = model.get_TYu(0,1);
   const double TYu02 = model.get_TYu(0,2);
   const double TYu10 = model.get_TYu(1,0);
   const double TYu11 = model.get_TYu(1,1);
   const double TYu12 = model.get_TYu(1,2);
   const double TYu20 = model.get_TYu(2,0);
   const double TYu21 = model.get_TYu(2,1);
   const double TYu22 = model.get_TYu(2,2);
   const double MassB = model.get_MassB();
   const double MassWB = model.get_MassWB();
   const double MassG = model.get_MassG();
   const double MassBp = model.get_MassBp();

   double coeff = 33.28*Power(g1,4)*MassB*Yu22 + 48*Power(g2,4)*MassWB*Yu22 +
      9.6*Power(g1,3)*(-1.7333333333333334*g1*TYu22 + 3.466666666666667*g1*MassB*
      Yu22) + 4*Power(g2,3)*(-6*g2*TYu22 + 12*g2*MassWB*Yu22) + 2*MassBp*Yu22*Sqr(
      gN)*(9*Sqr(Qdp) + 9*Sqr(QDxbarp) + 9*Sqr(QDxp) + 3*Sqr(Qep) + 6*Sqr(QH1p) +
      6*Sqr(QH2p) + 2*Sqr(QHpbarp) + 2*Sqr(QHpp) + 6*Sqr(QLp) + 18*Sqr(QQp) + 3*
      Sqr(QSp) + 9*Sqr(Qup))*(4*Sqr(gN)*Sqr(QH2p) + 4*Sqr(gN)*Sqr(QQp) + 4*Sqr(gN)
      *Sqr(Qup)) + Power(gN,3)*(9*Sqr(Qdp) + 9*Sqr(QDxbarp) + 9*Sqr(QDxp) + 3*Sqr(
      Qep) + 6*Sqr(QH1p) + 6*Sqr(QH2p) + 2*Sqr(QHpbarp) + 2*Sqr(QHpp) + 6*Sqr(QLp)
      + 18*Sqr(QQp) + 3*Sqr(QSp) + 9*Sqr(Qup))*(-4*gN*TYu22*Sqr(QH2p) - 4*gN*
      TYu22*Sqr(QQp) - 4*gN*TYu22*Sqr(Qup) + Yu22*(8*gN*MassBp*Sqr(QH2p) + 8*gN*
      MassBp*Sqr(QQp) + 8*gN*MassBp*Sqr(Qup))) + 2*(Yd00*Yu20 + Yd01*Yu21 + Yd02*
      Yu22)*(5*(Yd02*(TYd00*Yd00 + TYd01*Yd01 + TYd02*Yd02) + Yd12*(TYd00*Yd10 +
      TYd01*Yd11 + TYd02*Yd12) + Yd22*(TYd00*Yd20 + TYd01*Yd21 + TYd02*Yd22)) +
      Yu02*(TYd00*Yu00 + TYd01*Yu01 + TYd02*Yu02) + Yu12*(TYd00*Yu10 + TYd01*Yu11
      + TYd02*Yu12) + Yu22*(TYd00*Yu20 + TYd01*Yu21 + TYd02*Yu22) + 2*(TYu02*(Yd00
      *Yu00 + Yd01*Yu01 + Yd02*Yu02) + TYu12*(Yd00*Yu10 + Yd01*Yu11 + Yd02*Yu12) +
      TYu22*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22)) - 0.4666666666666667*TYd02*Sqr(
      g1) - 3*TYd02*Sqr(g2) - 5.333333333333333*TYd02*Sqr(g3) + TYd02*Sqr(Lambdax)
      - 2*TYd02*Sqr(gN)*Sqr(Qdp) - 2*TYd02*Sqr(gN)*Sqr(QH1p) - 2*TYd02*Sqr(gN)*
      Sqr(QQp) + Yd02*(2*Lambdax*TLambdax + 6*(TYd00*Yd00 + TYd01*Yd01 + TYd02*
      Yd02 + TYd10*Yd10 + TYd11*Yd11 + TYd12*Yd12 + TYd20*Yd20 + TYd21*Yd21 +
      TYd22*Yd22) + 2*(TYe00*Ye00 + TYe01*Ye01 + TYe02*Ye02 + TYe10*Ye10 + TYe11*
      Ye11 + TYe12*Ye12 + TYe20*Ye20 + TYe21*Ye21 + TYe22*Ye22) +
      0.9333333333333333*MassB*Sqr(g1) + 6*MassWB*Sqr(g2) + 10.666666666666666*
      MassG*Sqr(g3) + 4*MassBp*Sqr(gN)*Sqr(Qdp) + 4*MassBp*Sqr(gN)*Sqr(QH1p) + 4*
      MassBp*Sqr(gN)*Sqr(QQp)) + 4*(TYd12*(Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) +
      TYd22*(Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22) + TYd02*(Sqr(Yd00) + Sqr(Yd01) +
      Sqr(Yd02))) + 3*TYd02*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(
      Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + TYd02*(Sqr(Ye00) +
      Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(
      Ye21) + Sqr(Ye22))) + 2*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22)*(5*(Yd02*(TYd10*
      Yd00 + TYd11*Yd01 + TYd12*Yd02) + Yd12*(TYd10*Yd10 + TYd11*Yd11 + TYd12*Yd12
      ) + Yd22*(TYd10*Yd20 + TYd11*Yd21 + TYd12*Yd22)) + Yu02*(TYd10*Yu00 + TYd11*
      Yu01 + TYd12*Yu02) + Yu12*(TYd10*Yu10 + TYd11*Yu11 + TYd12*Yu12) + Yu22*(
      TYd10*Yu20 + TYd11*Yu21 + TYd12*Yu22) + 2*(TYu02*(Yd10*Yu00 + Yd11*Yu01 +
      Yd12*Yu02) + TYu12*(Yd10*Yu10 + Yd11*Yu11 + Yd12*Yu12) + TYu22*(Yd10*Yu20 +
      Yd11*Yu21 + Yd12*Yu22)) - 0.4666666666666667*TYd12*Sqr(g1) - 3*TYd12*Sqr(g2)
      - 5.333333333333333*TYd12*Sqr(g3) + TYd12*Sqr(Lambdax) - 2*TYd12*Sqr(gN)*
      Sqr(Qdp) - 2*TYd12*Sqr(gN)*Sqr(QH1p) - 2*TYd12*Sqr(gN)*Sqr(QQp) + Yd12*(2*
      Lambdax*TLambdax + 6*(TYd00*Yd00 + TYd01*Yd01 + TYd02*Yd02 + TYd10*Yd10 +
      TYd11*Yd11 + TYd12*Yd12 + TYd20*Yd20 + TYd21*Yd21 + TYd22*Yd22) + 2*(TYe00*
      Ye00 + TYe01*Ye01 + TYe02*Ye02 + TYe10*Ye10 + TYe11*Ye11 + TYe12*Ye12 +
      TYe20*Ye20 + TYe21*Ye21 + TYe22*Ye22) + 0.9333333333333333*MassB*Sqr(g1) + 6
      *MassWB*Sqr(g2) + 10.666666666666666*MassG*Sqr(g3) + 4*MassBp*Sqr(gN)*Sqr(
      Qdp) + 4*MassBp*Sqr(gN)*Sqr(QH1p) + 4*MassBp*Sqr(gN)*Sqr(QQp)) + 4*(TYd02*(
      Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + TYd22*(Yd10*Yd20 + Yd11*Yd21 + Yd12*
      Yd22) + TYd12*(Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12))) + 3*TYd12*(Sqr(Yd00) +
      Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(
      Yd21) + Sqr(Yd22)) + TYd12*(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) +
      Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22))) + 2*(Yd20*Yu20 +
      Yd21*Yu21 + Yd22*Yu22)*(5*(Yd02*(TYd20*Yd00 + TYd21*Yd01 + TYd22*Yd02) +
      Yd12*(TYd20*Yd10 + TYd21*Yd11 + TYd22*Yd12) + Yd22*(TYd20*Yd20 + TYd21*Yd21
      + TYd22*Yd22)) + Yu02*(TYd20*Yu00 + TYd21*Yu01 + TYd22*Yu02) + Yu12*(TYd20*
      Yu10 + TYd21*Yu11 + TYd22*Yu12) + Yu22*(TYd20*Yu20 + TYd21*Yu21 + TYd22*Yu22
      ) + 2*(TYu02*(Yd20*Yu00 + Yd21*Yu01 + Yd22*Yu02) + TYu12*(Yd20*Yu10 + Yd21*
      Yu11 + Yd22*Yu12) + TYu22*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22)) -
      0.4666666666666667*TYd22*Sqr(g1) - 3*TYd22*Sqr(g2) - 5.333333333333333*TYd22
      *Sqr(g3) + TYd22*Sqr(Lambdax) - 2*TYd22*Sqr(gN)*Sqr(Qdp) - 2*TYd22*Sqr(gN)*
      Sqr(QH1p) - 2*TYd22*Sqr(gN)*Sqr(QQp) + Yd22*(2*Lambdax*TLambdax + 6*(TYd00*
      Yd00 + TYd01*Yd01 + TYd02*Yd02 + TYd10*Yd10 + TYd11*Yd11 + TYd12*Yd12 +
      TYd20*Yd20 + TYd21*Yd21 + TYd22*Yd22) + 2*(TYe00*Ye00 + TYe01*Ye01 + TYe02*
      Ye02 + TYe10*Ye10 + TYe11*Ye11 + TYe12*Ye12 + TYe20*Ye20 + TYe21*Ye21 +
      TYe22*Ye22) + 0.9333333333333333*MassB*Sqr(g1) + 6*MassWB*Sqr(g2) +
      10.666666666666666*MassG*Sqr(g3) + 4*MassBp*Sqr(gN)*Sqr(Qdp) + 4*MassBp*Sqr(
      gN)*Sqr(QH1p) + 4*MassBp*Sqr(gN)*Sqr(QQp)) + 3*TYd22*(Sqr(Yd00) + Sqr(Yd01)
      + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) +
      Sqr(Yd22)) + 4*(TYd02*(Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22) + TYd12*(Yd10*Yd20
      + Yd11*Yd21 + Yd12*Yd22) + TYd22*(Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22))) +
      TYd22*(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12)
      + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22))) + (TYu20*Yd02 + 2*TYd02*Yu20)*(Yu00*(
      Yd00*Yu00 + Yd01*Yu01 + Yd02*Yu02) + Yu10*(Yd00*Yu10 + Yd01*Yu11 + Yd02*Yu12
      ) + Yu20*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22) + 3*(Yd10*(Yd00*Yd10 + Yd01*
      Yd11 + Yd02*Yd12) + Yd20*(Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22) + Yd00*(Sqr(
      Yd00) + Sqr(Yd01) + Sqr(Yd02))) + Yd00*(-0.4666666666666667*Sqr(g1) - 3*Sqr(
      g2) - 5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(Qdp) - 2*Sqr(
      gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QQp) + 3*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) +
      Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr
      (Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(
      Ye20) + Sqr(Ye21) + Sqr(Ye22))) + (TYu21*Yd02 + 2*TYd02*Yu21)*(Yu01*(Yd00*
      Yu00 + Yd01*Yu01 + Yd02*Yu02) + Yu11*(Yd00*Yu10 + Yd01*Yu11 + Yd02*Yu12) +
      Yu21*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22) + 3*(Yd11*(Yd00*Yd10 + Yd01*Yd11 +
      Yd02*Yd12) + Yd21*(Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22) + Yd01*(Sqr(Yd00) +
      Sqr(Yd01) + Sqr(Yd02))) + Yd01*(-0.4666666666666667*Sqr(g1) - 3*Sqr(g2) -
      5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(Qdp) - 2*Sqr(gN)*
      Sqr(QH1p) - 2*Sqr(gN)*Sqr(QQp) + 3*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(
      Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(
      Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20
      ) + Sqr(Ye21) + Sqr(Ye22))) + (TYu20*Yd00 + TYu21*Yd01 + 2*TYu22*Yd02 + 2*
      TYd02*Yu22)*(Yu02*(Yd00*Yu00 + Yd01*Yu01 + Yd02*Yu02) + Yu12*(Yd00*Yu10 +
      Yd01*Yu11 + Yd02*Yu12) + Yu22*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22) + 3*(Yd12*
      (Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + Yd22*(Yd00*Yd20 + Yd01*Yd21 + Yd02*
      Yd22) + Yd02*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02))) + Yd02*(
      -0.4666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(
      Lambdax) - 2*Sqr(gN)*Sqr(Qdp) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QQp) + 3
      *(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) +
      Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr
      (Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22))) + (
      TYu20*Yd12 + 2*TYd12*Yu20)*(Yu00*(Yd10*Yu00 + Yd11*Yu01 + Yd12*Yu02) + Yu10*
      (Yd10*Yu10 + Yd11*Yu11 + Yd12*Yu12) + Yu20*(Yd10*Yu20 + Yd11*Yu21 + Yd12*
      Yu22) + 3*(Yd00*(Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + Yd20*(Yd10*Yd20 + Yd11
      *Yd21 + Yd12*Yd22) + Yd10*(Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12))) + Yd10*(
      -0.4666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(
      Lambdax) - 2*Sqr(gN)*Sqr(Qdp) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QQp) + 3
      *(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) +
      Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr
      (Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22))) + (
      TYu21*Yd12 + 2*TYd12*Yu21)*(Yu01*(Yd10*Yu00 + Yd11*Yu01 + Yd12*Yu02) + Yu11*
      (Yd10*Yu10 + Yd11*Yu11 + Yd12*Yu12) + Yu21*(Yd10*Yu20 + Yd11*Yu21 + Yd12*
      Yu22) + 3*(Yd01*(Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + Yd21*(Yd10*Yd20 + Yd11
      *Yd21 + Yd12*Yd22) + Yd11*(Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12))) + Yd11*(
      -0.4666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(
      Lambdax) - 2*Sqr(gN)*Sqr(Qdp) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QQp) + 3
      *(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) +
      Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr
      (Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22))) + (
      TYu20*Yd10 + TYu21*Yd11 + 2*TYu22*Yd12 + 2*TYd12*Yu22)*(Yu02*(Yd10*Yu00 +
      Yd11*Yu01 + Yd12*Yu02) + Yu12*(Yd10*Yu10 + Yd11*Yu11 + Yd12*Yu12) + Yu22*(
      Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22) + 3*(Yd02*(Yd00*Yd10 + Yd01*Yd11 + Yd02*
      Yd12) + Yd22*(Yd10*Yd20 + Yd11*Yd21 + Yd12*Yd22) + Yd12*(Sqr(Yd10) + Sqr(
      Yd11) + Sqr(Yd12))) + Yd12*(-0.4666666666666667*Sqr(g1) - 3*Sqr(g2) -
      5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(Qdp) - 2*Sqr(gN)*
      Sqr(QH1p) - 2*Sqr(gN)*Sqr(QQp) + 3*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(
      Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(
      Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20
      ) + Sqr(Ye21) + Sqr(Ye22))) + (TYu20*Yd22 + 2*TYd22*Yu20)*(Yu00*(Yd20*Yu00 +
      Yd21*Yu01 + Yd22*Yu02) + Yu10*(Yd20*Yu10 + Yd21*Yu11 + Yd22*Yu12) + Yu20*(
      Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22) + 3*(Yd00*(Yd00*Yd20 + Yd01*Yd21 + Yd02*
      Yd22) + Yd10*(Yd10*Yd20 + Yd11*Yd21 + Yd12*Yd22) + Yd20*(Sqr(Yd20) + Sqr(
      Yd21) + Sqr(Yd22))) + Yd20*(-0.4666666666666667*Sqr(g1) - 3*Sqr(g2) -
      5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(Qdp) - 2*Sqr(gN)*
      Sqr(QH1p) - 2*Sqr(gN)*Sqr(QQp) + 3*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(
      Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(
      Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20
      ) + Sqr(Ye21) + Sqr(Ye22))) + (TYu21*Yd22 + 2*TYd22*Yu21)*(Yu01*(Yd20*Yu00 +
      Yd21*Yu01 + Yd22*Yu02) + Yu11*(Yd20*Yu10 + Yd21*Yu11 + Yd22*Yu12) + Yu21*(
      Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22) + 3*(Yd01*(Yd00*Yd20 + Yd01*Yd21 + Yd02*
      Yd22) + Yd11*(Yd10*Yd20 + Yd11*Yd21 + Yd12*Yd22) + Yd21*(Sqr(Yd20) + Sqr(
      Yd21) + Sqr(Yd22))) + Yd21*(-0.4666666666666667*Sqr(g1) - 3*Sqr(g2) -
      5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(Qdp) - 2*Sqr(gN)*
      Sqr(QH1p) - 2*Sqr(gN)*Sqr(QQp) + 3*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(
      Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(
      Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20
      ) + Sqr(Ye21) + Sqr(Ye22))) + (TYu20*Yd20 + TYu21*Yd21 + 2*TYu22*Yd22 + 2*
      TYd22*Yu22)*(Yu02*(Yd20*Yu00 + Yd21*Yu01 + Yd22*Yu02) + Yu12*(Yd20*Yu10 +
      Yd21*Yu11 + Yd22*Yu12) + Yu22*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22) + 3*(Yd02*
      (Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22) + Yd12*(Yd10*Yd20 + Yd11*Yd21 + Yd12*
      Yd22) + Yd22*(Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22))) + Yd22*(
      -0.4666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(
      Lambdax) - 2*Sqr(gN)*Sqr(Qdp) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QQp) + 3
      *(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) +
      Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr
      (Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22))) + (2*
      Lambdax*TYu22 + 2*TLambdax*Yu22)*(4*Power(Lambdax,3) - 0.6*Lambdax*Sqr(g1) -
      3*Lambdax*Sqr(g2) + 3*Lambdax*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) +
      Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) +
      Sqr(Kappa22)) + 2*Lambdax*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(
      Lambda1210) + Sqr(Lambda1211)) - 2*Lambdax*Sqr(gN)*Sqr(QH1p) - 2*Lambdax*Sqr
      (gN)*Sqr(QH2p) - 2*Lambdax*Sqr(gN)*Sqr(QSp) + 3*Lambdax*(Sqr(Yd00) + Sqr(
      Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21
      ) + Sqr(Yd22)) + Lambdax*(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) +
      Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22)) + 3*Lambdax*(Sqr(
      Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20
      ) + Sqr(Yu21) + Sqr(Yu22))) + 6*Yu00*Yu22*(Yd00*(TYu00*Yd00 + TYu01*Yd01 +
      TYu02*Yd02) + Yd10*(TYu00*Yd10 + TYu01*Yd11 + TYu02*Yd12) + Yd20*(TYu00*Yd20
      + TYu01*Yd21 + TYu02*Yd22) + 2*(TYd00*(Yd00*Yu00 + Yd01*Yu01 + Yd02*Yu02) +
      TYd10*(Yd10*Yu00 + Yd11*Yu01 + Yd12*Yu02) + TYd20*(Yd20*Yu00 + Yd21*Yu01 +
      Yd22*Yu02)) + 5*(Yu00*(TYu00*Yu00 + TYu01*Yu01 + TYu02*Yu02) + Yu10*(TYu00*
      Yu10 + TYu01*Yu11 + TYu02*Yu12) + Yu20*(TYu00*Yu20 + TYu01*Yu21 + TYu02*Yu22
      )) - 0.8666666666666667*TYu00*Sqr(g1) - 3*TYu00*Sqr(g2) - 5.333333333333333*
      TYu00*Sqr(g3) + TYu00*Sqr(Lambdax) - 2*TYu00*Sqr(gN)*Sqr(QH2p) - 2*TYu00*Sqr
      (gN)*Sqr(QQp) - 2*TYu00*Sqr(gN)*Sqr(Qup) + Yu00*(2*Lambdax*TLambdax + 6*(
      TYu00*Yu00 + TYu01*Yu01 + TYu02*Yu02 + TYu10*Yu10 + TYu11*Yu11 + TYu12*Yu12
      + TYu20*Yu20 + TYu21*Yu21 + TYu22*Yu22) + 1.7333333333333334*MassB*Sqr(g1) +
      6*MassWB*Sqr(g2) + 10.666666666666666*MassG*Sqr(g3) + 4*MassBp*Sqr(gN)*Sqr(
      QH2p) + 4*MassBp*Sqr(gN)*Sqr(QQp) + 4*MassBp*Sqr(gN)*Sqr(Qup)) + 4*(TYu10*(
      Yu00*Yu10 + Yu01*Yu11 + Yu02*Yu12) + TYu20*(Yu00*Yu20 + Yu01*Yu21 + Yu02*
      Yu22) + TYu00*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02))) + 3*TYu00*(Sqr(Yu00) +
      Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(
      Yu21) + Sqr(Yu22))) + 6*Yu01*Yu22*(Yd01*(TYu00*Yd00 + TYu01*Yd01 + TYu02*
      Yd02) + Yd11*(TYu00*Yd10 + TYu01*Yd11 + TYu02*Yd12) + Yd21*(TYu00*Yd20 +
      TYu01*Yd21 + TYu02*Yd22) + 2*(TYd01*(Yd00*Yu00 + Yd01*Yu01 + Yd02*Yu02) +
      TYd11*(Yd10*Yu00 + Yd11*Yu01 + Yd12*Yu02) + TYd21*(Yd20*Yu00 + Yd21*Yu01 +
      Yd22*Yu02)) + 5*(Yu01*(TYu00*Yu00 + TYu01*Yu01 + TYu02*Yu02) + Yu11*(TYu00*
      Yu10 + TYu01*Yu11 + TYu02*Yu12) + Yu21*(TYu00*Yu20 + TYu01*Yu21 + TYu02*Yu22
      )) - 0.8666666666666667*TYu01*Sqr(g1) - 3*TYu01*Sqr(g2) - 5.333333333333333*
      TYu01*Sqr(g3) + TYu01*Sqr(Lambdax) - 2*TYu01*Sqr(gN)*Sqr(QH2p) - 2*TYu01*Sqr
      (gN)*Sqr(QQp) - 2*TYu01*Sqr(gN)*Sqr(Qup) + Yu01*(2*Lambdax*TLambdax + 6*(
      TYu00*Yu00 + TYu01*Yu01 + TYu02*Yu02 + TYu10*Yu10 + TYu11*Yu11 + TYu12*Yu12
      + TYu20*Yu20 + TYu21*Yu21 + TYu22*Yu22) + 1.7333333333333334*MassB*Sqr(g1) +
      6*MassWB*Sqr(g2) + 10.666666666666666*MassG*Sqr(g3) + 4*MassBp*Sqr(gN)*Sqr(
      QH2p) + 4*MassBp*Sqr(gN)*Sqr(QQp) + 4*MassBp*Sqr(gN)*Sqr(Qup)) + 4*(TYu11*(
      Yu00*Yu10 + Yu01*Yu11 + Yu02*Yu12) + TYu21*(Yu00*Yu20 + Yu01*Yu21 + Yu02*
      Yu22) + TYu01*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02))) + 3*TYu01*(Sqr(Yu00) +
      Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(
      Yu21) + Sqr(Yu22))) + (6*Yu02*Yu22 + 4*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22))*
      (Yd02*(TYu00*Yd00 + TYu01*Yd01 + TYu02*Yd02) + Yd12*(TYu00*Yd10 + TYu01*Yd11
      + TYu02*Yd12) + Yd22*(TYu00*Yd20 + TYu01*Yd21 + TYu02*Yd22) + 2*(TYd02*(
      Yd00*Yu00 + Yd01*Yu01 + Yd02*Yu02) + TYd12*(Yd10*Yu00 + Yd11*Yu01 + Yd12*
      Yu02) + TYd22*(Yd20*Yu00 + Yd21*Yu01 + Yd22*Yu02)) + 5*(Yu02*(TYu00*Yu00 +
      TYu01*Yu01 + TYu02*Yu02) + Yu12*(TYu00*Yu10 + TYu01*Yu11 + TYu02*Yu12) +
      Yu22*(TYu00*Yu20 + TYu01*Yu21 + TYu02*Yu22)) - 0.8666666666666667*TYu02*Sqr(
      g1) - 3*TYu02*Sqr(g2) - 5.333333333333333*TYu02*Sqr(g3) + TYu02*Sqr(Lambdax)
      - 2*TYu02*Sqr(gN)*Sqr(QH2p) - 2*TYu02*Sqr(gN)*Sqr(QQp) - 2*TYu02*Sqr(gN)*
      Sqr(Qup) + Yu02*(2*Lambdax*TLambdax + 6*(TYu00*Yu00 + TYu01*Yu01 + TYu02*
      Yu02 + TYu10*Yu10 + TYu11*Yu11 + TYu12*Yu12 + TYu20*Yu20 + TYu21*Yu21 +
      TYu22*Yu22) + 1.7333333333333334*MassB*Sqr(g1) + 6*MassWB*Sqr(g2) +
      10.666666666666666*MassG*Sqr(g3) + 4*MassBp*Sqr(gN)*Sqr(QH2p) + 4*MassBp*Sqr
      (gN)*Sqr(QQp) + 4*MassBp*Sqr(gN)*Sqr(Qup)) + 4*(TYu12*(Yu00*Yu10 + Yu01*Yu11
      + Yu02*Yu12) + TYu22*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + TYu02*(Sqr(Yu00)
      + Sqr(Yu01) + Sqr(Yu02))) + 3*TYu02*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) +
      Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + 6*
      Yu10*Yu22*(Yd00*(TYu10*Yd00 + TYu11*Yd01 + TYu12*Yd02) + Yd10*(TYu10*Yd10 +
      TYu11*Yd11 + TYu12*Yd12) + Yd20*(TYu10*Yd20 + TYu11*Yd21 + TYu12*Yd22) + 2*(
      TYd00*(Yd00*Yu10 + Yd01*Yu11 + Yd02*Yu12) + TYd10*(Yd10*Yu10 + Yd11*Yu11 +
      Yd12*Yu12) + TYd20*(Yd20*Yu10 + Yd21*Yu11 + Yd22*Yu12)) + 5*(Yu00*(TYu10*
      Yu00 + TYu11*Yu01 + TYu12*Yu02) + Yu10*(TYu10*Yu10 + TYu11*Yu11 + TYu12*Yu12
      ) + Yu20*(TYu10*Yu20 + TYu11*Yu21 + TYu12*Yu22)) - 0.8666666666666667*TYu10*
      Sqr(g1) - 3*TYu10*Sqr(g2) - 5.333333333333333*TYu10*Sqr(g3) + TYu10*Sqr(
      Lambdax) - 2*TYu10*Sqr(gN)*Sqr(QH2p) - 2*TYu10*Sqr(gN)*Sqr(QQp) - 2*TYu10*
      Sqr(gN)*Sqr(Qup) + Yu10*(2*Lambdax*TLambdax + 6*(TYu00*Yu00 + TYu01*Yu01 +
      TYu02*Yu02 + TYu10*Yu10 + TYu11*Yu11 + TYu12*Yu12 + TYu20*Yu20 + TYu21*Yu21
      + TYu22*Yu22) + 1.7333333333333334*MassB*Sqr(g1) + 6*MassWB*Sqr(g2) +
      10.666666666666666*MassG*Sqr(g3) + 4*MassBp*Sqr(gN)*Sqr(QH2p) + 4*MassBp*Sqr
      (gN)*Sqr(QQp) + 4*MassBp*Sqr(gN)*Sqr(Qup)) + 4*(TYu00*(Yu00*Yu10 + Yu01*Yu11
      + Yu02*Yu12) + TYu20*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + TYu10*(Sqr(Yu10)
      + Sqr(Yu11) + Sqr(Yu12))) + 3*TYu10*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) +
      Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + 6*
      Yu11*Yu22*(Yd01*(TYu10*Yd00 + TYu11*Yd01 + TYu12*Yd02) + Yd11*(TYu10*Yd10 +
      TYu11*Yd11 + TYu12*Yd12) + Yd21*(TYu10*Yd20 + TYu11*Yd21 + TYu12*Yd22) + 2*(
      TYd01*(Yd00*Yu10 + Yd01*Yu11 + Yd02*Yu12) + TYd11*(Yd10*Yu10 + Yd11*Yu11 +
      Yd12*Yu12) + TYd21*(Yd20*Yu10 + Yd21*Yu11 + Yd22*Yu12)) + 5*(Yu01*(TYu10*
      Yu00 + TYu11*Yu01 + TYu12*Yu02) + Yu11*(TYu10*Yu10 + TYu11*Yu11 + TYu12*Yu12
      ) + Yu21*(TYu10*Yu20 + TYu11*Yu21 + TYu12*Yu22)) - 0.8666666666666667*TYu11*
      Sqr(g1) - 3*TYu11*Sqr(g2) - 5.333333333333333*TYu11*Sqr(g3) + TYu11*Sqr(
      Lambdax) - 2*TYu11*Sqr(gN)*Sqr(QH2p) - 2*TYu11*Sqr(gN)*Sqr(QQp) - 2*TYu11*
      Sqr(gN)*Sqr(Qup) + Yu11*(2*Lambdax*TLambdax + 6*(TYu00*Yu00 + TYu01*Yu01 +
      TYu02*Yu02 + TYu10*Yu10 + TYu11*Yu11 + TYu12*Yu12 + TYu20*Yu20 + TYu21*Yu21
      + TYu22*Yu22) + 1.7333333333333334*MassB*Sqr(g1) + 6*MassWB*Sqr(g2) +
      10.666666666666666*MassG*Sqr(g3) + 4*MassBp*Sqr(gN)*Sqr(QH2p) + 4*MassBp*Sqr
      (gN)*Sqr(QQp) + 4*MassBp*Sqr(gN)*Sqr(Qup)) + 4*(TYu01*(Yu00*Yu10 + Yu01*Yu11
      + Yu02*Yu12) + TYu21*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + TYu11*(Sqr(Yu10)
      + Sqr(Yu11) + Sqr(Yu12))) + 3*TYu11*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) +
      Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + (6
      *Yu12*Yu22 + 4*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22))*(Yd02*(TYu10*Yd00 +
      TYu11*Yd01 + TYu12*Yd02) + Yd12*(TYu10*Yd10 + TYu11*Yd11 + TYu12*Yd12) +
      Yd22*(TYu10*Yd20 + TYu11*Yd21 + TYu12*Yd22) + 2*(TYd02*(Yd00*Yu10 + Yd01*
      Yu11 + Yd02*Yu12) + TYd12*(Yd10*Yu10 + Yd11*Yu11 + Yd12*Yu12) + TYd22*(Yd20*
      Yu10 + Yd21*Yu11 + Yd22*Yu12)) + 5*(Yu02*(TYu10*Yu00 + TYu11*Yu01 + TYu12*
      Yu02) + Yu12*(TYu10*Yu10 + TYu11*Yu11 + TYu12*Yu12) + Yu22*(TYu10*Yu20 +
      TYu11*Yu21 + TYu12*Yu22)) - 0.8666666666666667*TYu12*Sqr(g1) - 3*TYu12*Sqr(
      g2) - 5.333333333333333*TYu12*Sqr(g3) + TYu12*Sqr(Lambdax) - 2*TYu12*Sqr(gN)
      *Sqr(QH2p) - 2*TYu12*Sqr(gN)*Sqr(QQp) - 2*TYu12*Sqr(gN)*Sqr(Qup) + Yu12*(2*
      Lambdax*TLambdax + 6*(TYu00*Yu00 + TYu01*Yu01 + TYu02*Yu02 + TYu10*Yu10 +
      TYu11*Yu11 + TYu12*Yu12 + TYu20*Yu20 + TYu21*Yu21 + TYu22*Yu22) +
      1.7333333333333334*MassB*Sqr(g1) + 6*MassWB*Sqr(g2) + 10.666666666666666*
      MassG*Sqr(g3) + 4*MassBp*Sqr(gN)*Sqr(QH2p) + 4*MassBp*Sqr(gN)*Sqr(QQp) + 4*
      MassBp*Sqr(gN)*Sqr(Qup)) + 4*(TYu02*(Yu00*Yu10 + Yu01*Yu11 + Yu02*Yu12) +
      TYu22*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + TYu12*(Sqr(Yu10) + Sqr(Yu11) +
      Sqr(Yu12))) + 3*TYu12*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(
      Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + (Yd00*Yd02 + Yd10*
      Yd12 + Yd20*Yd22 + 6*Yu20*Yu22 + 5*(Yu00*Yu02 + Yu10*Yu12 + Yu20*Yu22))*(
      Yd00*(TYu20*Yd00 + TYu21*Yd01 + TYu22*Yd02) + Yd10*(TYu20*Yd10 + TYu21*Yd11
      + TYu22*Yd12) + Yd20*(TYu20*Yd20 + TYu21*Yd21 + TYu22*Yd22) + 5*(Yu00*(TYu20
      *Yu00 + TYu21*Yu01 + TYu22*Yu02) + Yu10*(TYu20*Yu10 + TYu21*Yu11 + TYu22*
      Yu12) + Yu20*(TYu20*Yu20 + TYu21*Yu21 + TYu22*Yu22)) + 2*(TYd00*(Yd00*Yu20 +
      Yd01*Yu21 + Yd02*Yu22) + TYd10*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22) + TYd20*
      (Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22)) - 0.8666666666666667*TYu20*Sqr(g1) - 3*
      TYu20*Sqr(g2) - 5.333333333333333*TYu20*Sqr(g3) + TYu20*Sqr(Lambdax) - 2*
      TYu20*Sqr(gN)*Sqr(QH2p) - 2*TYu20*Sqr(gN)*Sqr(QQp) - 2*TYu20*Sqr(gN)*Sqr(Qup
      ) + Yu20*(2*Lambdax*TLambdax + 6*(TYu00*Yu00 + TYu01*Yu01 + TYu02*Yu02 +
      TYu10*Yu10 + TYu11*Yu11 + TYu12*Yu12 + TYu20*Yu20 + TYu21*Yu21 + TYu22*Yu22)
      + 1.7333333333333334*MassB*Sqr(g1) + 6*MassWB*Sqr(g2) + 10.666666666666666*
      MassG*Sqr(g3) + 4*MassBp*Sqr(gN)*Sqr(QH2p) + 4*MassBp*Sqr(gN)*Sqr(QQp) + 4*
      MassBp*Sqr(gN)*Sqr(Qup)) + 3*TYu20*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(
      Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)) + 4*(
      TYu00*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + TYu10*(Yu10*Yu20 + Yu11*Yu21 +
      Yu12*Yu22) + TYu20*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) + (Yd01*Yd02 + Yd11
      *Yd12 + Yd21*Yd22 + 6*Yu21*Yu22 + 5*(Yu01*Yu02 + Yu11*Yu12 + Yu21*Yu22))*(
      Yd01*(TYu20*Yd00 + TYu21*Yd01 + TYu22*Yd02) + Yd11*(TYu20*Yd10 + TYu21*Yd11
      + TYu22*Yd12) + Yd21*(TYu20*Yd20 + TYu21*Yd21 + TYu22*Yd22) + 5*(Yu01*(TYu20
      *Yu00 + TYu21*Yu01 + TYu22*Yu02) + Yu11*(TYu20*Yu10 + TYu21*Yu11 + TYu22*
      Yu12) + Yu21*(TYu20*Yu20 + TYu21*Yu21 + TYu22*Yu22)) + 2*(TYd01*(Yd00*Yu20 +
      Yd01*Yu21 + Yd02*Yu22) + TYd11*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22) + TYd21*
      (Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22)) - 0.8666666666666667*TYu21*Sqr(g1) - 3*
      TYu21*Sqr(g2) - 5.333333333333333*TYu21*Sqr(g3) + TYu21*Sqr(Lambdax) - 2*
      TYu21*Sqr(gN)*Sqr(QH2p) - 2*TYu21*Sqr(gN)*Sqr(QQp) - 2*TYu21*Sqr(gN)*Sqr(Qup
      ) + Yu21*(2*Lambdax*TLambdax + 6*(TYu00*Yu00 + TYu01*Yu01 + TYu02*Yu02 +
      TYu10*Yu10 + TYu11*Yu11 + TYu12*Yu12 + TYu20*Yu20 + TYu21*Yu21 + TYu22*Yu22)
      + 1.7333333333333334*MassB*Sqr(g1) + 6*MassWB*Sqr(g2) + 10.666666666666666*
      MassG*Sqr(g3) + 4*MassBp*Sqr(gN)*Sqr(QH2p) + 4*MassBp*Sqr(gN)*Sqr(QQp) + 4*
      MassBp*Sqr(gN)*Sqr(Qup)) + 3*TYu21*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(
      Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)) + 4*(
      TYu01*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + TYu11*(Yu10*Yu20 + Yu11*Yu21 +
      Yu12*Yu22) + TYu21*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) + (
      -0.8666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(
      Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) +
      Sqr(Yd02) + Sqr(Yd12) + Sqr(Yd22) + 6*Sqr(Yu22) + 5*(Sqr(Yu02) + Sqr(Yu12) +
      Sqr(Yu22)) + 4*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)) + 3*(Sqr(Yu00) + Sqr(
      Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21
      ) + Sqr(Yu22)))*(Yd02*(TYu20*Yd00 + TYu21*Yd01 + TYu22*Yd02) + Yd12*(TYu20*
      Yd10 + TYu21*Yd11 + TYu22*Yd12) + Yd22*(TYu20*Yd20 + TYu21*Yd21 + TYu22*Yd22
      ) + 5*(Yu02*(TYu20*Yu00 + TYu21*Yu01 + TYu22*Yu02) + Yu12*(TYu20*Yu10 +
      TYu21*Yu11 + TYu22*Yu12) + Yu22*(TYu20*Yu20 + TYu21*Yu21 + TYu22*Yu22)) + 2*
      (TYd02*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22) + TYd12*(Yd10*Yu20 + Yd11*Yu21 +
      Yd12*Yu22) + TYd22*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22)) - 0.8666666666666667
      *TYu22*Sqr(g1) - 3*TYu22*Sqr(g2) - 5.333333333333333*TYu22*Sqr(g3) + TYu22*
      Sqr(Lambdax) - 2*TYu22*Sqr(gN)*Sqr(QH2p) - 2*TYu22*Sqr(gN)*Sqr(QQp) - 2*
      TYu22*Sqr(gN)*Sqr(Qup) + Yu22*(2*Lambdax*TLambdax + 6*(TYu00*Yu00 + TYu01*
      Yu01 + TYu02*Yu02 + TYu10*Yu10 + TYu11*Yu11 + TYu12*Yu12 + TYu20*Yu20 +
      TYu21*Yu21 + TYu22*Yu22) + 1.7333333333333334*MassB*Sqr(g1) + 6*MassWB*Sqr(
      g2) + 10.666666666666666*MassG*Sqr(g3) + 4*MassBp*Sqr(gN)*Sqr(QH2p) + 4*
      MassBp*Sqr(gN)*Sqr(QQp) + 4*MassBp*Sqr(gN)*Sqr(Qup)) + 3*TYu22*(Sqr(Yu00) +
      Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(
      Yu21) + Sqr(Yu22)) + 4*(TYu02*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + TYu12*(
      Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + TYu22*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22
      )))) + (6*TYu22*Yu00 + 5*TYu20*Yu02 + 4*TYu02*Yu20 + 6*TYu00*Yu22)*(Yd00*(
      Yd00*Yu00 + Yd01*Yu01 + Yd02*Yu02) + Yd10*(Yd10*Yu00 + Yd11*Yu01 + Yd12*Yu02
      ) + Yd20*(Yd20*Yu00 + Yd21*Yu01 + Yd22*Yu02) + 3*(Yu10*(Yu00*Yu10 + Yu01*
      Yu11 + Yu02*Yu12) + Yu20*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu00*(Sqr(
      Yu00) + Sqr(Yu01) + Sqr(Yu02))) + Yu00*(-0.8666666666666667*Sqr(g1) - 3*Sqr(
      g2) - 5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr
      (gN)*Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) +
      Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) + (
      6*TYu22*Yu01 + 5*TYu21*Yu02 + 4*TYu02*Yu21 + 6*TYu01*Yu22)*(Yd01*(Yd00*Yu00
      + Yd01*Yu01 + Yd02*Yu02) + Yd11*(Yd10*Yu00 + Yd11*Yu01 + Yd12*Yu02) + Yd21*(
      Yd20*Yu00 + Yd21*Yu01 + Yd22*Yu02) + 3*(Yu11*(Yu00*Yu10 + Yu01*Yu11 + Yu02*
      Yu12) + Yu21*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu01*(Sqr(Yu00) + Sqr(
      Yu01) + Sqr(Yu02))) + Yu01*(-0.8666666666666667*Sqr(g1) - 3*Sqr(g2) -
      5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*
      Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(
      Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) + (6*
      TYu22*Yu02 + 5*(TYu20*Yu00 + TYu21*Yu01 + 2*TYu22*Yu02) + 10*TYu02*Yu22)*(
      Yd02*(Yd00*Yu00 + Yd01*Yu01 + Yd02*Yu02) + Yd12*(Yd10*Yu00 + Yd11*Yu01 +
      Yd12*Yu02) + Yd22*(Yd20*Yu00 + Yd21*Yu01 + Yd22*Yu02) + 3*(Yu12*(Yu00*Yu10 +
      Yu01*Yu11 + Yu02*Yu12) + Yu22*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu02*(
      Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02))) + Yu02*(-0.8666666666666667*Sqr(g1) - 3*
      Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2
      *Sqr(gN)*Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02
      ) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))))
      + (6*TYu22*Yu10 + 5*TYu20*Yu12 + 4*TYu12*Yu20 + 6*TYu10*Yu22)*(Yd00*(Yd00*
      Yu10 + Yd01*Yu11 + Yd02*Yu12) + Yd10*(Yd10*Yu10 + Yd11*Yu11 + Yd12*Yu12) +
      Yd20*(Yd20*Yu10 + Yd21*Yu11 + Yd22*Yu12) + 3*(Yu00*(Yu00*Yu10 + Yu01*Yu11 +
      Yu02*Yu12) + Yu20*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + Yu10*(Sqr(Yu10) +
      Sqr(Yu11) + Sqr(Yu12))) + Yu10*(-0.8666666666666667*Sqr(g1) - 3*Sqr(g2) -
      5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*
      Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(
      Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) + (6*
      TYu22*Yu11 + 5*TYu21*Yu12 + 4*TYu12*Yu21 + 6*TYu11*Yu22)*(Yd01*(Yd00*Yu10 +
      Yd01*Yu11 + Yd02*Yu12) + Yd11*(Yd10*Yu10 + Yd11*Yu11 + Yd12*Yu12) + Yd21*(
      Yd20*Yu10 + Yd21*Yu11 + Yd22*Yu12) + 3*(Yu01*(Yu00*Yu10 + Yu01*Yu11 + Yu02*
      Yu12) + Yu21*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + Yu11*(Sqr(Yu10) + Sqr(
      Yu11) + Sqr(Yu12))) + Yu11*(-0.8666666666666667*Sqr(g1) - 3*Sqr(g2) -
      5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*
      Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(
      Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) + (6*
      TYu22*Yu12 + 5*(TYu20*Yu10 + TYu21*Yu11 + 2*TYu22*Yu12) + 10*TYu12*Yu22)*(
      Yd02*(Yd00*Yu10 + Yd01*Yu11 + Yd02*Yu12) + Yd12*(Yd10*Yu10 + Yd11*Yu11 +
      Yd12*Yu12) + Yd22*(Yd20*Yu10 + Yd21*Yu11 + Yd22*Yu12) + 3*(Yu02*(Yu00*Yu10 +
      Yu01*Yu11 + Yu02*Yu12) + Yu22*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + Yu12*(
      Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12))) + Yu12*(-0.8666666666666667*Sqr(g1) - 3*
      Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2
      *Sqr(gN)*Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02
      ) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))))
      + (2*(TYd02*Yd00 + TYd12*Yd10 + TYd22*Yd20) + 6*TYu22*Yu20 + 4*(TYu02*Yu00
      + TYu12*Yu10 + 2*TYu22*Yu20) + 11*TYu20*Yu22)*(Yd00*(Yd00*Yu20 + Yd01*Yu21 +
      Yd02*Yu22) + Yd10*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22) + Yd20*(Yd20*Yu20 +
      Yd21*Yu21 + Yd22*Yu22) + 3*(Yu00*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu10*
      (Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + Yu20*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22
      ))) + Yu20*(-0.8666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(
      g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QQp) - 2*Sqr(gN)*
      Sqr(Qup) + 3*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) +
      Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) + (2*(TYd02*Yd01 + TYd12*
      Yd11 + TYd22*Yd21) + 6*TYu22*Yu21 + 4*(TYu02*Yu01 + TYu12*Yu11 + 2*TYu22*
      Yu21) + 11*TYu21*Yu22)*(Yd01*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22) + Yd11*(
      Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22) + Yd21*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22
      ) + 3*(Yu01*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu11*(Yu10*Yu20 + Yu11*
      Yu21 + Yu12*Yu22) + Yu21*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + Yu21*(
      -0.8666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(
      Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3
      *(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) +
      Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) + (2*Lambdax*TLambdax + 2*(TYd02*Yd02 +
      TYd12*Yd12 + TYd22*Yd22) + 12*TYu22*Yu22 + 6*(TYu00*Yu00 + TYu01*Yu01 +
      TYu02*Yu02 + TYu10*Yu10 + TYu11*Yu11 + TYu12*Yu12 + TYu20*Yu20 + TYu21*Yu21
      + TYu22*Yu22) + 4*(TYu02*Yu02 + TYu12*Yu12 + 2*TYu22*Yu22) + 5*(TYu20*Yu20 +
      TYu21*Yu21 + 2*TYu22*Yu22) + 1.7333333333333334*MassB*Sqr(g1) + 6*MassWB*
      Sqr(g2) + 10.666666666666666*MassG*Sqr(g3) + 4*MassBp*Sqr(gN)*Sqr(QH2p) + 4*
      MassBp*Sqr(gN)*Sqr(QQp) + 4*MassBp*Sqr(gN)*Sqr(Qup))*(Yd02*(Yd00*Yu20 + Yd01
      *Yu21 + Yd02*Yu22) + Yd12*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22) + Yd22*(Yd20*
      Yu20 + Yd21*Yu21 + Yd22*Yu22) + 3*(Yu02*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22)
      + Yu12*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + Yu22*(Sqr(Yu20) + Sqr(Yu21) +
      Sqr(Yu22))) + Yu22*(-0.8666666666666667*Sqr(g1) - 3*Sqr(g2) -
      5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*
      Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(
      Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) + 2*
      Lambdax*Yu22*(6*Lambdax*(Kappa00*TKappa00 + Kappa01*TKappa01 + Kappa02*
      TKappa02 + Kappa10*TKappa10 + Kappa11*TKappa11 + Kappa12*TKappa12 + Kappa20*
      TKappa20 + Kappa21*TKappa21 + Kappa22*TKappa22) + 4*Lambdax*(Lambda1200*
      TLambda1200 + Lambda1201*TLambda1201 + Lambda1210*TLambda1210 + Lambda1211*
      TLambda1211) + 6*Lambdax*(TYd00*Yd00 + TYd01*Yd01 + TYd02*Yd02 + TYd10*Yd10
      + TYd11*Yd11 + TYd12*Yd12 + TYd20*Yd20 + TYd21*Yd21 + TYd22*Yd22) + 2*
      Lambdax*(TYe00*Ye00 + TYe01*Ye01 + TYe02*Ye02 + TYe10*Ye10 + TYe11*Ye11 +
      TYe12*Ye12 + TYe20*Ye20 + TYe21*Ye21 + TYe22*Ye22) + 6*Lambdax*(TYu00*Yu00 +
      TYu01*Yu01 + TYu02*Yu02 + TYu10*Yu10 + TYu11*Yu11 + TYu12*Yu12 + TYu20*Yu20
      + TYu21*Yu21 + TYu22*Yu22) + 1.2*Lambdax*MassB*Sqr(g1) + 6*Lambdax*MassWB*
      Sqr(g2) + 4*Lambdax*MassBp*Sqr(gN)*Sqr(QH1p) + 4*Lambdax*MassBp*Sqr(gN)*Sqr(
      QH2p) + 4*Lambdax*MassBp*Sqr(gN)*Sqr(QSp) + TLambdax*(-0.6*Sqr(g1) - 3*Sqr(
      g2) + 3*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(
      Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*(
      Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 12*
      Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QSp
      ) + 3*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12)
      + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) +
      Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22) + 3*(
      Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(
      Yu20) + Sqr(Yu21) + Sqr(Yu22))));


   return twoLoop*coeff;
}

} // namespace flexiblesusy
