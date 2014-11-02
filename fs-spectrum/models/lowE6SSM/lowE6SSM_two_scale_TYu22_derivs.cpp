#include "lowE6SSM_two_scale_tuning_calculator.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

// one loop
double lowE6SSM_tuning_calculator::deriv_dbeta_one_loop_TYu22_dLambdax() const 
{

   const double Lambdax = model.get_Lambdax();
   const double Yu22 = model.get_Yu(2,2);

   const double ALambdax = get_ALambdax();
   const double TYu22 = model.get_TYu(2,2);

   double deriv = 2*Lambdax*TYu22 + 4*ALambdax*Lambdax*Yu22;


   return oneOver16PiSqr*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_one_loop_TYu22_dALambdax() const 
{

   const double Lambdax = model.get_Lambdax();
   const double Yu22 = model.get_Yu(2,2);

   double deriv = 2*Yu22*Sqr(Lambdax);


   return oneOver16PiSqr*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_one_loop_TYu22_dAYu22() const 
{
   const lowE6SSM_input_parameters inputs = model.get_input();

   const auto QQp = inputs.QQp;
   const auto QH2p = inputs.QH2p;
   const auto Qup = inputs.Qup;
   
   const double Yd02 = model.get_Yd(0,2);
   const double Yd12 = model.get_Yd(1,2);
   const double Yd22 = model.get_Yd(2,2);
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

   double deriv = 6*Power(Yu22,3) - 0.8666666666666667*Yu22*Sqr(g1) - 3*Yu22*
      Sqr(g2) - 5.333333333333333*Yu22*Sqr(g3) + Yu22*Sqr(Lambdax) - 2*Yu22*Sqr(gN
      )*Sqr(QH2p) - 2*Yu22*Sqr(gN)*Sqr(QQp) - 2*Yu22*Sqr(gN)*Sqr(Qup) + Yu22*Sqr(
      Yd02) + Yu22*Sqr(Yd12) + Yu22*Sqr(Yd22) + 5*(Power(Yu22,3) + Yu22*Sqr(Yu02)
      + Yu22*Sqr(Yu12)) + 4*Yu22*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)) + 3*Yu22*(Sqr
      (Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(
      Yu20) + Sqr(Yu21) + Sqr(Yu22));


   return oneOver16PiSqr*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_one_loop_TYu22_dMassB() const 
{

   const double Yu22 = model.get_Yu(2,2);
   const double g1 = model.get_g1();

   double deriv = 1.7333333333333334*Yu22*Sqr(g1);


   return oneOver16PiSqr*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_one_loop_TYu22_dMassWB() const 
{

   const double Yu22 = model.get_Yu(2,2);
   const double g2 = model.get_g2();

   double deriv = 6*Yu22*Sqr(g2);


   return oneOver16PiSqr*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_one_loop_TYu22_dMassG() const 
{

   const double Yu22 = model.get_Yu(2,2);
   const double g3 = model.get_g3();

   double deriv = 10.666666666666666*Yu22*Sqr(g3);


   return oneOver16PiSqr*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_one_loop_TYu22_dMassBp() const 
{
   const lowE6SSM_input_parameters inputs = model.get_input();

   const auto QQp = inputs.QQp;
   const auto QH2p = inputs.QH2p;
   const auto Qup = inputs.Qup;

   const double Yu22 = model.get_Yu(2,2);
   const double gN = model.get_gN();

   double deriv = Yu22*(4*Sqr(gN)*Sqr(QH2p) + 4*Sqr(gN)*Sqr(QQp) + 4*Sqr(gN)*
      Sqr(Qup));


   return oneOver16PiSqr*deriv;
}

// two loop
double lowE6SSM_tuning_calculator::deriv_dbeta_two_loop_TYu22_dLambdax() const 
{
   const lowE6SSM_input_parameters inputs = model.get_input();

   const auto QH1p = inputs.QH1p;
   const auto QH2p = inputs.QH2p;
   const auto QSp = inputs.QSp;
   
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
   const double ALambdax = get_ALambdax();
   const double TYu02 = model.get_TYu(0,2);
   const double TYu12 = model.get_TYu(1,2);
   const double TYu20 = model.get_TYu(2,0);
   const double TYu21 = model.get_TYu(2,1);
   const double TYu22 = model.get_TYu(2,2);
   const double MassBp = model.get_MassBp();

   double deriv = -12*Power(Lambdax,3)*TYu22 - 2*Lambdax*(Yd02*(TYu20*Yd00 +
      TYu21*Yd01 + TYu22*Yd02) + Yd12*(TYu20*Yd10 + TYu21*Yd11 + TYu22*Yd12) +
      Yd22*(TYu20*Yd20 + TYu21*Yd21 + TYu22*Yd22)) - 10*Lambdax*(Yu02*(TYu20*Yu00
      + TYu21*Yu01 + TYu22*Yu02) + Yu12*(TYu20*Yu10 + TYu21*Yu11 + TYu22*Yu12) +
      Yu22*(TYu20*Yu20 + TYu21*Yu21 + TYu22*Yu22)) - 4*Lambdax*(TYd02*(Yd00*Yu20 +
      Yd01*Yu21 + Yd02*Yu22) + TYd12*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22) + TYd22*
      (Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22)) - 4*ALambdax*Lambdax*(Yd02*(Yd00*Yu20 +
      Yd01*Yu21 + Yd02*Yu22) + Yd12*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22) + Yd22*(
      Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22)) - 6*Lambdax*TYu22*(Sqr(Kappa00) + Sqr(
      Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(
      Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) - 4*Lambdax*TYu22*(Sqr(Lambda1200) +
      Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 4*Lambdax*TYu22*Sqr(
      gN)*Sqr(QH1p) - 4*Lambdax*TYu22*Sqr(gN)*Sqr(QH2p) + 4*Lambdax*TYu22*Sqr(gN)*
      Sqr(QSp) - 6*Lambdax*TYu22*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) +
      Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) - 2*Lambdax*TYu22
      *(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) +
      Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22)) - 0.008888888888888889*Yu22*(5400*
      ALambdax*Power(Lambdax,3) + 225*Lambdax*(3*(Kappa00*TKappa00 + Kappa01*
      TKappa01 + Kappa02*TKappa02 + Kappa10*TKappa10 + Kappa11*TKappa11 + Kappa12*
      TKappa12 + Kappa20*TKappa20 + Kappa21*TKappa21 + Kappa22*TKappa22) + 2*(
      Lambda1200*TLambda1200 + Lambda1201*TLambda1201 + Lambda1210*TLambda1210 +
      Lambda1211*TLambda1211) + 3*(TYd00*Yd00 + TYd01*Yd01 + TYd02*Yd02 + TYd10*
      Yd10 + TYd11*Yd11 + TYd12*Yd12 + TYd20*Yd20 + TYd21*Yd21 + TYd22*Yd22) +
      TYe00*Ye00 + TYe01*Ye01 + TYe02*Ye02 + TYe10*Ye10 + TYe11*Ye11 + TYe12*Ye12
      + TYe20*Ye20 + TYe21*Ye21 + TYe22*Ye22 + 2*MassBp*Sqr(gN)*Sqr(QH1p) - 2*
      MassBp*Sqr(gN)*Sqr(QH2p) + 2*MassBp*Sqr(gN)*Sqr(QSp) + ALambdax*(3*(Sqr(
      Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(
      Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*(Sqr(Lambda1200)
      + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) - 2*Sqr(gN)*Sqr(QH1p)
      + 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QSp) + 3*(Sqr(Yd00) + Sqr(Yd01) + Sqr
      (Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(
      Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(
      Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22))) + 225*(Lambdax*(3*(Kappa00*
      TKappa00 + Kappa01*TKappa01 + Kappa02*TKappa02 + Kappa10*TKappa10 + Kappa11*
      TKappa11 + Kappa12*TKappa12 + Kappa20*TKappa20 + Kappa21*TKappa21 + Kappa22*
      TKappa22) + 2*(Lambda1200*TLambda1200 + Lambda1201*TLambda1201 + Lambda1210*
      TLambda1210 + Lambda1211*TLambda1211) + 3*(TYd00*Yd00 + TYd01*Yd01 + TYd02*
      Yd02 + TYd10*Yd10 + TYd11*Yd11 + TYd12*Yd12 + TYd20*Yd20 + TYd21*Yd21 +
      TYd22*Yd22) + TYe00*Ye00 + TYe01*Ye01 + TYe02*Ye02 + TYe10*Ye10 + TYe11*Ye11
      + TYe12*Ye12 + TYe20*Ye20 + TYe21*Ye21 + TYe22*Ye22 + 2*MassBp*Sqr(gN)*Sqr(
      QH1p) - 2*MassBp*Sqr(gN)*Sqr(QH2p) + 2*MassBp*Sqr(gN)*Sqr(QSp)) + ALambdax*
      Lambdax*(3*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(
      Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*(
      Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) - 2*
      Sqr(gN)*Sqr(QH1p) + 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QSp) + 3*(Sqr(Yd00)
      + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) +
      Sqr(Yd21) + Sqr(Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr
      (Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22)))) - 8*Lambdax*(TYu02
      *(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + TYu12*(Yu10*Yu20 + Yu11*Yu21 + Yu12*
      Yu22) + TYu22*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) - 12*ALambdax*Lambdax*(
      Yu02*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu12*(Yu10*Yu20 + Yu11*Yu21 +
      Yu12*Yu22) + Yu22*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_two_loop_TYu22_dALambdax() const 
{
   const lowE6SSM_input_parameters inputs = model.get_input();

   const auto QH1p = inputs.QH1p;
   const auto QH2p = inputs.QH2p;
   const auto QSp = inputs.QSp;
   
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
   const double gN = model.get_gN();

   double deriv = -2*(Yd02*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22) + Yd12*(Yd10*
      Yu20 + Yd11*Yu21 + Yd12*Yu22) + Yd22*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22))*
      Sqr(Lambdax) - 0.008888888888888889*Yu22*(1350*Power(Lambdax,4) + 225*Sqr(
      Lambdax)*(3*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr
      (Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*(
      Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) - 2*
      Sqr(gN)*Sqr(QH1p) + 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QSp) + 3*(Sqr(Yd00)
      + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) +
      Sqr(Yd21) + Sqr(Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr
      (Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22))) - 6*Sqr(Lambdax)*(
      Yu02*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu12*(Yu10*Yu20 + Yu11*Yu21 +
      Yu12*Yu22) + Yu22*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_two_loop_TYu22_dAYu22() const 
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

   double deriv = 8.695555555555556*Power(g1,4)*Yu22 + 16.5*Power(g2,4)*Yu22 +
      14.222222222222221*Power(g3,4)*Yu22 - 3*Power(Lambdax,4)*Yu22 + 16*Power(gN,
      4)*Power(QH2p,4)*Yu22 + 40*Power(gN,4)*Power(QQp,4)*Yu22 + 22*Power(gN,4)*
      Power(Qup,4)*Yu22 - 2*Yu22*(Yu20*(Yd00*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22) +
      Yd10*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22) + Yd20*(Yd20*Yu20 + Yd21*Yu21 +
      Yd22*Yu22)) + Yu21*(Yd01*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22) + Yd11*(Yd10*
      Yu20 + Yd11*Yu21 + Yd12*Yu22) + Yd21*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22)) +
      Yu22*(Yd02*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22) + Yd12*(Yd10*Yu20 + Yd11*Yu21
      + Yd12*Yu22) + Yd22*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22))) - 3*Yu22*(Yd00*(
      Yu00*(Yd00*Yu00 + Yd01*Yu01 + Yd02*Yu02) + Yu10*(Yd00*Yu10 + Yd01*Yu11 +
      Yd02*Yu12) + Yu20*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22)) + Yd01*(Yu01*(Yd00*
      Yu00 + Yd01*Yu01 + Yd02*Yu02) + Yu11*(Yd00*Yu10 + Yd01*Yu11 + Yd02*Yu12) +
      Yu21*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22)) + Yd02*(Yu02*(Yd00*Yu00 + Yd01*
      Yu01 + Yd02*Yu02) + Yu12*(Yd00*Yu10 + Yd01*Yu11 + Yd02*Yu12) + Yu22*(Yd00*
      Yu20 + Yd01*Yu21 + Yd02*Yu22)) + Yd10*(Yu00*(Yd10*Yu00 + Yd11*Yu01 + Yd12*
      Yu02) + Yu10*(Yd10*Yu10 + Yd11*Yu11 + Yd12*Yu12) + Yu20*(Yd10*Yu20 + Yd11*
      Yu21 + Yd12*Yu22)) + Yd11*(Yu01*(Yd10*Yu00 + Yd11*Yu01 + Yd12*Yu02) + Yu11*(
      Yd10*Yu10 + Yd11*Yu11 + Yd12*Yu12) + Yu21*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22
      )) + Yd12*(Yu02*(Yd10*Yu00 + Yd11*Yu01 + Yd12*Yu02) + Yu12*(Yd10*Yu10 + Yd11
      *Yu11 + Yd12*Yu12) + Yu22*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22)) + Yd20*(Yu00*
      (Yd20*Yu00 + Yd21*Yu01 + Yd22*Yu02) + Yu10*(Yd20*Yu10 + Yd21*Yu11 + Yd22*
      Yu12) + Yu20*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22)) + Yd21*(Yu01*(Yd20*Yu00 +
      Yd21*Yu01 + Yd22*Yu02) + Yu11*(Yd20*Yu10 + Yd21*Yu11 + Yd22*Yu12) + Yu21*(
      Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22)) + Yd22*(Yu02*(Yd20*Yu00 + Yd21*Yu01 +
      Yd22*Yu02) + Yu12*(Yd20*Yu10 + Yd21*Yu11 + Yd22*Yu12) + Yu22*(Yd20*Yu20 +
      Yd21*Yu21 + Yd22*Yu22))) + Yu22*Sqr(g1)*Sqr(g2) + 3.022222222222222*Yu22*Sqr
      (g1)*Sqr(g3) + 8*Yu22*Sqr(g2)*Sqr(g3) - 3*Yu22*(Sqr(Kappa00) + Sqr(Kappa01)
      + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) +
      Sqr(Kappa21) + Sqr(Kappa22))*Sqr(Lambdax) - 2*Yu22*(Sqr(Lambda1200) + Sqr(
      Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211))*Sqr(Lambdax) + 2*Yu22*Sqr(
      gN)*Sqr(Lambdax)*Sqr(QH1p) + 1.2*Yu22*Sqr(g1)*Sqr(gN)*Sqr(QH2p) + 6*Yu22*Sqr
      (g2)*Sqr(gN)*Sqr(QH2p) - 2*Yu22*Sqr(gN)*Sqr(Lambdax)*Sqr(QH2p) + 18*Power(gN
      ,4)*Yu22*Sqr(Qdp)*Sqr(QH2p) + 18*Power(gN,4)*Yu22*Sqr(QDxbarp)*Sqr(QH2p) +
      18*Power(gN,4)*Yu22*Sqr(QDxp)*Sqr(QH2p) + 6*Power(gN,4)*Yu22*Sqr(Qep)*Sqr(
      QH2p) + 12*Power(gN,4)*Yu22*Sqr(QH1p)*Sqr(QH2p) + 4*Power(gN,4)*Yu22*Sqr(
      QH2p)*Sqr(QHpbarp) + 4*Power(gN,4)*Yu22*Sqr(QH2p)*Sqr(QHpp) + 12*Power(gN,4)
      *Yu22*Sqr(QH2p)*Sqr(QLp) + 0.13333333333333333*Yu22*Sqr(g1)*Sqr(gN)*Sqr(QQp)
      + 6*Yu22*Sqr(g2)*Sqr(gN)*Sqr(QQp) + 10.666666666666666*Yu22*Sqr(g3)*Sqr(gN)
      *Sqr(QQp) + 18*Power(gN,4)*Yu22*Sqr(Qdp)*Sqr(QQp) + 18*Power(gN,4)*Yu22*Sqr(
      QDxbarp)*Sqr(QQp) + 18*Power(gN,4)*Yu22*Sqr(QDxp)*Sqr(QQp) + 6*Power(gN,4)*
      Yu22*Sqr(Qep)*Sqr(QQp) + 12*Power(gN,4)*Yu22*Sqr(QH1p)*Sqr(QQp) + 48*Power(
      gN,4)*Yu22*Sqr(QH2p)*Sqr(QQp) + 4*Power(gN,4)*Yu22*Sqr(QHpbarp)*Sqr(QQp) + 4
      *Power(gN,4)*Yu22*Sqr(QHpp)*Sqr(QQp) + 12*Power(gN,4)*Yu22*Sqr(QLp)*Sqr(QQp)
      + 2*Yu22*Sqr(gN)*Sqr(Lambdax)*Sqr(QSp) + 6*Power(gN,4)*Yu22*Sqr(QH2p)*Sqr(
      QSp) + 6*Power(gN,4)*Yu22*Sqr(QQp)*Sqr(QSp) + 2.1333333333333333*Yu22*Sqr(g1
      )*Sqr(gN)*Sqr(Qup) + 10.666666666666666*Yu22*Sqr(g3)*Sqr(gN)*Sqr(Qup) + 18*
      Power(gN,4)*Yu22*Sqr(Qdp)*Sqr(Qup) + 18*Power(gN,4)*Yu22*Sqr(QDxbarp)*Sqr(
      Qup) + 18*Power(gN,4)*Yu22*Sqr(QDxp)*Sqr(Qup) + 6*Power(gN,4)*Yu22*Sqr(Qep)*
      Sqr(Qup) + 12*Power(gN,4)*Yu22*Sqr(QH1p)*Sqr(Qup) + 30*Power(gN,4)*Yu22*Sqr(
      QH2p)*Sqr(Qup) + 4*Power(gN,4)*Yu22*Sqr(QHpbarp)*Sqr(Qup) + 4*Power(gN,4)*
      Yu22*Sqr(QHpp)*Sqr(Qup) + 12*Power(gN,4)*Yu22*Sqr(QLp)*Sqr(Qup) + 54*Power(
      gN,4)*Yu22*Sqr(QQp)*Sqr(Qup) + 6*Power(gN,4)*Yu22*Sqr(QSp)*Sqr(Qup) - 3*Yu22
      *Sqr(Lambdax)*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) +
      Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + 0.4*Sqr(g1)*(Yu22*Sqr(Yd02)
      + Yu22*Sqr(Yd12) + Yu22*Sqr(Yd22)) - Sqr(Lambdax)*(Yu22*Sqr(Yd02) + Yu22*
      Sqr(Yd12) + Yu22*Sqr(Yd22)) + 2*Sqr(gN)*Sqr(Qdp)*(Yu22*Sqr(Yd02) + Yu22*Sqr(
      Yd12) + Yu22*Sqr(Yd22)) + 2*Sqr(gN)*Sqr(QH1p)*(Yu22*Sqr(Yd02) + Yu22*Sqr(
      Yd12) + Yu22*Sqr(Yd22)) - 2*Sqr(gN)*Sqr(QQp)*(Yu22*Sqr(Yd02) + Yu22*Sqr(Yd12
      ) + Yu22*Sqr(Yd22)) - 3*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr
      (Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22))*(Yu22*Sqr(Yd02) +
      Yu22*Sqr(Yd12) + Yu22*Sqr(Yd22)) - 2*(Yd02*(Yd00*(Yd00*Yd02*Yu22 + Yd10*Yd12
      *Yu22 + Yd20*Yd22*Yu22) + Yd01*(Yd01*Yd02*Yu22 + Yd11*Yd12*Yu22 + Yd21*Yd22*
      Yu22) + Yd02*(Yu22*Sqr(Yd02) + Yu22*Sqr(Yd12) + Yu22*Sqr(Yd22))) + Yd12*(
      Yd10*(Yd00*Yd02*Yu22 + Yd10*Yd12*Yu22 + Yd20*Yd22*Yu22) + Yd11*(Yd01*Yd02*
      Yu22 + Yd11*Yd12*Yu22 + Yd21*Yd22*Yu22) + Yd12*(Yu22*Sqr(Yd02) + Yu22*Sqr(
      Yd12) + Yu22*Sqr(Yd22))) + Yd22*(Yd20*(Yd00*Yd02*Yu22 + Yd10*Yd12*Yu22 +
      Yd20*Yd22*Yu22) + Yd21*(Yd01*Yd02*Yu22 + Yd11*Yd12*Yu22 + Yd21*Yd22*Yu22) +
      Yd22*(Yu22*Sqr(Yd02) + Yu22*Sqr(Yd12) + Yu22*Sqr(Yd22)))) - 4*(Yu02*(Yu00*(
      Yd00*Yd02*Yu22 + Yd10*Yd12*Yu22 + Yd20*Yd22*Yu22) + Yu01*(Yd01*Yd02*Yu22 +
      Yd11*Yd12*Yu22 + Yd21*Yd22*Yu22) + Yu02*(Yu22*Sqr(Yd02) + Yu22*Sqr(Yd12) +
      Yu22*Sqr(Yd22))) + Yu12*(Yu10*(Yd00*Yd02*Yu22 + Yd10*Yd12*Yu22 + Yd20*Yd22*
      Yu22) + Yu11*(Yd01*Yd02*Yu22 + Yd11*Yd12*Yu22 + Yd21*Yd22*Yu22) + Yu12*(Yu22
      *Sqr(Yd02) + Yu22*Sqr(Yd12) + Yu22*Sqr(Yd22))) + Yu22*(Yu20*(Yd00*Yd02*Yu22
      + Yd10*Yd12*Yu22 + Yd20*Yd22*Yu22) + Yu21*(Yd01*Yd02*Yu22 + Yd11*Yd12*Yu22 +
      Yd21*Yd22*Yu22) + Yu22*(Yu22*Sqr(Yd02) + Yu22*Sqr(Yd12) + Yu22*Sqr(Yd22))))
      - Yu22*Sqr(Lambdax)*(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(
      Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22)) - (Yu22*Sqr(Yd02) +
      Yu22*Sqr(Yd12) + Yu22*Sqr(Yd22))*(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(
      Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22)) + 12*Sqr(
      g2)*(Power(Yu22,3) + Yu22*Sqr(Yu02) + Yu22*Sqr(Yu12)) - 5*Sqr(Lambdax)*(
      Power(Yu22,3) + Yu22*Sqr(Yu02) + Yu22*Sqr(Yu12)) + 10*Sqr(gN)*Sqr(QH2p)*(
      Power(Yu22,3) + Yu22*Sqr(Yu02) + Yu22*Sqr(Yu12)) + 6*Sqr(gN)*Sqr(QQp)*(Power
      (Yu22,3) + Yu22*Sqr(Yu02) + Yu22*Sqr(Yu12)) - 6*Sqr(gN)*Sqr(Qup)*(Power(Yu22
      ,3) + Yu22*Sqr(Yu02) + Yu22*Sqr(Yu12)) + 1.2*Yu22*Sqr(g1)*(Sqr(Yu20) + Sqr(
      Yu21) + Sqr(Yu22)) + 6*Yu22*Sqr(g2)*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)) - 4*
      Yu22*Sqr(Lambdax)*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)) + 8*Yu22*Sqr(gN)*Sqr(
      QH2p)*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)) + 0.8*Yu22*Sqr(g1)*(Sqr(Yu00) +
      Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(
      Yu21) + Sqr(Yu22)) + 16*Yu22*Sqr(g3)*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) +
      Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)) - 6*
      Yu22*Sqr(gN)*Sqr(QH2p)*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(
      Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)) + 6*Yu22*Sqr(gN)*Sqr(
      QQp)*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12)
      + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)) + 6*Yu22*Sqr(gN)*Sqr(Qup)*(Sqr(Yu00) +
      Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(
      Yu21) + Sqr(Yu22)) - 15*(Power(Yu22,3) + Yu22*Sqr(Yu02) + Yu22*Sqr(Yu12))*(
      Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(
      Yu20) + Sqr(Yu21) + Sqr(Yu22)) - 12*Yu22*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))
      *(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) +
      Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)) - 18*Sqr(Yu22)*(Yu02*(Yu00*Yu20 + Yu01*
      Yu21 + Yu02*Yu22) + Yu12*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + Yu22*(Sqr(
      Yu20) + Sqr(Yu21) + Sqr(Yu22))) - 8*(Power(Yu22,3)*(Sqr(Yu20) + Sqr(Yu21) +
      Sqr(Yu22)) + Yu22*Sqr(Yu02)*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)) + Yu22*Sqr(
      Yu12)*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) - 6*Yu22*(Yu20*(Yu00*(Yu00*Yu20 +
      Yu01*Yu21 + Yu02*Yu22) + Yu10*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + Yu20*(
      Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + Yu21*(Yu01*(Yu00*Yu20 + Yu01*Yu21 +
      Yu02*Yu22) + Yu11*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + Yu21*(Sqr(Yu20) +
      Sqr(Yu21) + Sqr(Yu22))) + Yu22*(Yu02*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) +
      Yu12*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + Yu22*(Sqr(Yu20) + Sqr(Yu21) + Sqr
      (Yu22)))) - 9*Yu22*(Yu00*(Yu10*(Yu00*Yu10 + Yu01*Yu11 + Yu02*Yu12) + Yu20*(
      Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu00*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02)
      )) + Yu01*(Yu11*(Yu00*Yu10 + Yu01*Yu11 + Yu02*Yu12) + Yu21*(Yu00*Yu20 + Yu01
      *Yu21 + Yu02*Yu22) + Yu01*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02))) + Yu02*(Yu12*
      (Yu00*Yu10 + Yu01*Yu11 + Yu02*Yu12) + Yu22*(Yu00*Yu20 + Yu01*Yu21 + Yu02*
      Yu22) + Yu02*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02))) + Yu10*(Yu00*(Yu00*Yu10 +
      Yu01*Yu11 + Yu02*Yu12) + Yu20*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + Yu10*(
      Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12))) + Yu11*(Yu01*(Yu00*Yu10 + Yu01*Yu11 +
      Yu02*Yu12) + Yu21*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + Yu11*(Sqr(Yu10) +
      Sqr(Yu11) + Sqr(Yu12))) + Yu12*(Yu02*(Yu00*Yu10 + Yu01*Yu11 + Yu02*Yu12) +
      Yu22*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + Yu12*(Sqr(Yu10) + Sqr(Yu11) + Sqr
      (Yu12))) + Yu20*(Yu00*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu10*(Yu10*Yu20
      + Yu11*Yu21 + Yu12*Yu22) + Yu20*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + Yu21*
      (Yu01*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu11*(Yu10*Yu20 + Yu11*Yu21 +
      Yu12*Yu22) + Yu21*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + Yu22*(Yu02*(Yu00*
      Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu12*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) +
      Yu22*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) - 0.008888888888888889*Yu22*(675*
      (Yd02*Yu22*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22) + Yd12*Yu22*(Yd10*Yu20 + Yd11
      *Yu21 + Yd12*Yu22) + Yd22*Yu22*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22)) - 180*
      Sqr(g1)*Sqr(Yu22) - 3600*Sqr(g3)*Sqr(Yu22) + 1350*Sqr(gN)*Sqr(QH2p)*Sqr(Yu22
      ) - 1350*Sqr(gN)*Sqr(QQp)*Sqr(Yu22) - 1350*Sqr(gN)*Sqr(Qup)*Sqr(Yu22) + 4050
      *(Yu02*Yu22*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu12*Yu22*(Yu10*Yu20 +
      Yu11*Yu21 + Yu12*Yu22) + Sqr(Yu22)*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) - 6
      *(Yu02*(Yu02*(Power(Yu22,3) + Yu22*Sqr(Yu02) + Yu22*Sqr(Yu12)) + Yu00*(Yu00*
      Yu02*Yu22 + Yu10*Yu12*Yu22 + Yu20*Sqr(Yu22)) + Yu01*(Yu01*Yu02*Yu22 + Yu11*
      Yu12*Yu22 + Yu21*Sqr(Yu22))) + Yu12*(Yu12*(Power(Yu22,3) + Yu22*Sqr(Yu02) +
      Yu22*Sqr(Yu12)) + Yu10*(Yu00*Yu02*Yu22 + Yu10*Yu12*Yu22 + Yu20*Sqr(Yu22)) +
      Yu11*(Yu01*Yu02*Yu22 + Yu11*Yu12*Yu22 + Yu21*Sqr(Yu22))) + Yu22*(Yu22*(Power
      (Yu22,3) + Yu22*Sqr(Yu02) + Yu22*Sqr(Yu12)) + Yu20*(Yu00*Yu02*Yu22 + Yu10*
      Yu12*Yu22 + Yu20*Sqr(Yu22)) + Yu21*(Yu01*Yu02*Yu22 + Yu11*Yu12*Yu22 + Yu21*
      Sqr(Yu22))));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_two_loop_TYu22_dMassB() const 
{
   const lowE6SSM_input_parameters inputs = model.get_input();

   const auto QQp = inputs.QQp;
   const auto QH2p = inputs.QH2p;
   const auto Qup = inputs.Qup;
   
   const double Yd00 = model.get_Yd(0,0);
   const double Yd01 = model.get_Yd(0,1);
   const double Yd02 = model.get_Yd(0,2);
   const double Yd10 = model.get_Yd(1,0);
   const double Yd11 = model.get_Yd(1,1);
   const double Yd12 = model.get_Yd(1,2);
   const double Yd20 = model.get_Yd(2,0);
   const double Yd21 = model.get_Yd(2,1);
   const double Yd22 = model.get_Yd(2,2);
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

   double deriv = -0.8*(Yd02*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22) + Yd12*(Yd10*
      Yu20 + Yd11*Yu21 + Yd12*Yu22) + Yd22*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22))*
      Sqr(g1) - 0.8*Sqr(g1)*(Yu02*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu12*(Yu10
      *Yu20 + Yu11*Yu21 + Yu12*Yu22) + Yu22*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) -
      0.008888888888888889*Yu22*(3913*Power(g1,4) + 225*Sqr(g1)*Sqr(g2) + 680*Sqr
      (g1)*Sqr(g3) + 270*Sqr(g1)*Sqr(gN)*Sqr(QH2p) + 30*Sqr(g1)*Sqr(gN)*Sqr(QQp) +
      480*Sqr(g1)*Sqr(gN)*Sqr(Qup) + 180*Sqr(g1)*(Sqr(Yu00) + Sqr(Yu01) + Sqr(
      Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22
      )));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_two_loop_TYu22_dMassWB() const 
{
   const lowE6SSM_input_parameters inputs = model.get_input();

   const auto QQp = inputs.QQp;
   const auto QH2p = inputs.QH2p;

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

   double deriv = -0.008888888888888889*Yu22*(7425*Power(g2,4) + 225*Sqr(g1)*
      Sqr(g2) + 1800*Sqr(g2)*Sqr(g3) + 1350*Sqr(g2)*Sqr(gN)*Sqr(QH2p) + 1350*Sqr(
      g2)*Sqr(gN)*Sqr(QQp)) - 12*Sqr(g2)*(Yu02*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22)
      + Yu12*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + Yu22*(Sqr(Yu20) + Sqr(Yu21) +
      Sqr(Yu22)));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_two_loop_TYu22_dMassG() const 
{
   const lowE6SSM_input_parameters inputs = model.get_input();

   const auto QQp = inputs.QQp;
   const auto Qup = inputs.Qup;
   
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

   double deriv = -0.008888888888888889*Yu22*(6400*Power(g3,4) + 680*Sqr(g1)*
      Sqr(g3) + 1800*Sqr(g2)*Sqr(g3) + 2400*Sqr(g3)*Sqr(gN)*Sqr(QQp) + 2400*Sqr(g3
      )*Sqr(gN)*Sqr(Qup) + 3600*Sqr(g3)*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(
      Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_two_loop_TYu22_dMassBp() const 
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

   double deriv = (Yd02*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22) + Yd12*(Yd10*Yu20 +
      Yd11*Yu21 + Yd12*Yu22) + Yd22*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22))*(-4*Sqr(
      gN)*Sqr(Qdp) - 4*Sqr(gN)*Sqr(QH1p) + 4*Sqr(gN)*Sqr(QQp)) - 12*Sqr(gN)*Sqr(
      QH2p)*(Yu02*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu12*(Yu10*Yu20 + Yu11*
      Yu21 + Yu12*Yu22) + Yu22*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) - 4*Sqr(gN)*
      Sqr(QQp)*(Yu02*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu12*(Yu10*Yu20 + Yu11*
      Yu21 + Yu12*Yu22) + Yu22*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + 4*Sqr(gN)*
      Sqr(Qup)*(Yu02*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu12*(Yu10*Yu20 + Yu11*
      Yu21 + Yu12*Yu22) + Yu22*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) -
      0.008888888888888889*Yu22*(7200*Power(gN,4)*Power(QH2p,4) + 18000*Power(gN,4
      )*Power(QQp,4) + 9900*Power(gN,4)*Power(Qup,4) + 270*Sqr(g1)*Sqr(gN)*Sqr(
      QH2p) + 1350*Sqr(g2)*Sqr(gN)*Sqr(QH2p) + 8100*Power(gN,4)*Sqr(Qdp)*Sqr(QH2p)
      + 8100*Power(gN,4)*Sqr(QDxbarp)*Sqr(QH2p) + 8100*Power(gN,4)*Sqr(QDxp)*Sqr(
      QH2p) + 2700*Power(gN,4)*Sqr(Qep)*Sqr(QH2p) + 5400*Power(gN,4)*Sqr(QH1p)*Sqr
      (QH2p) + 1800*Power(gN,4)*Sqr(QH2p)*Sqr(QHpbarp) + 1800*Power(gN,4)*Sqr(QH2p
      )*Sqr(QHpp) + 5400*Power(gN,4)*Sqr(QH2p)*Sqr(QLp) + 30*Sqr(g1)*Sqr(gN)*Sqr(
      QQp) + 1350*Sqr(g2)*Sqr(gN)*Sqr(QQp) + 2400*Sqr(g3)*Sqr(gN)*Sqr(QQp) + 8100*
      Power(gN,4)*Sqr(Qdp)*Sqr(QQp) + 8100*Power(gN,4)*Sqr(QDxbarp)*Sqr(QQp) +
      8100*Power(gN,4)*Sqr(QDxp)*Sqr(QQp) + 2700*Power(gN,4)*Sqr(Qep)*Sqr(QQp) +
      5400*Power(gN,4)*Sqr(QH1p)*Sqr(QQp) + 21600*Power(gN,4)*Sqr(QH2p)*Sqr(QQp) +
      1800*Power(gN,4)*Sqr(QHpbarp)*Sqr(QQp) + 1800*Power(gN,4)*Sqr(QHpp)*Sqr(QQp
      ) + 5400*Power(gN,4)*Sqr(QLp)*Sqr(QQp) + 2700*Power(gN,4)*Sqr(QH2p)*Sqr(QSp)
      + 2700*Power(gN,4)*Sqr(QQp)*Sqr(QSp) + 225*Sqr(Lambdax)*(2*Sqr(gN)*Sqr(QH1p
      ) - 2*Sqr(gN)*Sqr(QH2p) + 2*Sqr(gN)*Sqr(QSp)) + 480*Sqr(g1)*Sqr(gN)*Sqr(Qup)
      + 2400*Sqr(g3)*Sqr(gN)*Sqr(Qup) + 8100*Power(gN,4)*Sqr(Qdp)*Sqr(Qup) + 8100
      *Power(gN,4)*Sqr(QDxbarp)*Sqr(Qup) + 8100*Power(gN,4)*Sqr(QDxp)*Sqr(Qup) +
      2700*Power(gN,4)*Sqr(Qep)*Sqr(Qup) + 5400*Power(gN,4)*Sqr(QH1p)*Sqr(Qup) +
      13500*Power(gN,4)*Sqr(QH2p)*Sqr(Qup) + 1800*Power(gN,4)*Sqr(QHpbarp)*Sqr(Qup
      ) + 1800*Power(gN,4)*Sqr(QHpp)*Sqr(Qup) + 5400*Power(gN,4)*Sqr(QLp)*Sqr(Qup)
      + 24300*Power(gN,4)*Sqr(QQp)*Sqr(Qup) + 2700*Power(gN,4)*Sqr(QSp)*Sqr(Qup)
      + 1350*Sqr(gN)*(-Sqr(QH2p) + Sqr(QQp) + Sqr(Qup))*(Sqr(Yu00) + Sqr(Yu01) +
      Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(
      Yu22)));


   return twoLoop*deriv;
}

// leading log
double lowE6SSM_tuning_calculator::deriv_dlead_log_one_loop_TYu22_dLambdax() const 
{
   const lowE6SSM_input_parameters inputs = model.get_input();

   const auto QQp = inputs.QQp;
   const auto QH1p = inputs.QH1p;
   const auto QH2p = inputs.QH2p;
   const auto Qup = inputs.Qup;
   const auto QSp = inputs.QSp;
   
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
   const double ALambdax = get_ALambdax();
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

   double deriv = 2*Lambdax*Yd00*(TYu20*Yd02 + 2*TYd02*Yu20) + 2*Lambdax*Yd10*(
      TYu20*Yd12 + 2*TYd12*Yu20) + 2*Lambdax*Yd20*(TYu20*Yd22 + 2*TYd22*Yu20) + 2*
      Lambdax*Yd01*(TYu21*Yd02 + 2*TYd02*Yu21) + 2*Lambdax*Yd11*(TYu21*Yd12 + 2*
      TYd12*Yu21) + 2*Lambdax*Yd21*(TYu21*Yd22 + 2*TYd22*Yu21) + 6*Yu00*(2*Lambdax
      *TYu00 + 4*ALambdax*Lambdax*Yu00)*Yu22 + 6*Yu01*(2*Lambdax*TYu01 + 4*
      ALambdax*Lambdax*Yu01)*Yu22 + 6*Yu10*(2*Lambdax*TYu10 + 4*ALambdax*Lambdax*
      Yu10)*Yu22 + 6*Yu11*(2*Lambdax*TYu11 + 4*ALambdax*Lambdax*Yu11)*Yu22 + 2*
      Lambdax*Yd02*(TYu20*Yd00 + TYu21*Yd01 + 2*TYu22*Yd02 + 2*TYd02*Yu22) + 2*
      Lambdax*Yd12*(TYu20*Yd10 + TYu21*Yd11 + 2*TYu22*Yd12 + 2*TYd12*Yu22) + 2*
      Lambdax*Yd22*(TYu20*Yd20 + TYu21*Yd21 + 2*TYu22*Yd22 + 2*TYd22*Yu22) + 2*
      Lambdax*Yu00*(6*TYu22*Yu00 + 5*TYu20*Yu02 + 4*TYu02*Yu20 + 6*TYu00*Yu22) + 2
      *Lambdax*Yu01*(6*TYu22*Yu01 + 5*TYu21*Yu02 + 4*TYu02*Yu21 + 6*TYu01*Yu22) +
      2*Lambdax*Yu02*(6*TYu22*Yu02 + 5*(TYu20*Yu00 + TYu21*Yu01 + 2*TYu22*Yu02) +
      10*TYu02*Yu22) + 2*Lambdax*Yu10*(6*TYu22*Yu10 + 5*TYu20*Yu12 + 4*TYu12*Yu20
      + 6*TYu10*Yu22) + 2*Lambdax*Yu11*(6*TYu22*Yu11 + 5*TYu21*Yu12 + 4*TYu12*Yu21
      + 6*TYu11*Yu22) + 2*Lambdax*Yu12*(6*TYu22*Yu12 + 5*(TYu20*Yu10 + TYu21*Yu11
      + 2*TYu22*Yu12) + 10*TYu12*Yu22) + 2*Lambdax*Yu20*(2*(TYd02*Yd00 + TYd12*
      Yd10 + TYd22*Yd20) + 6*TYu22*Yu20 + 4*(TYu02*Yu00 + TYu12*Yu10 + 2*TYu22*
      Yu20) + 11*TYu20*Yu22) + 2*Lambdax*Yu21*(2*(TYd02*Yd01 + TYd12*Yd11 + TYd22*
      Yd21) + 6*TYu22*Yu21 + 4*(TYu02*Yu01 + TYu12*Yu11 + 2*TYu22*Yu21) + 11*TYu21
      *Yu22) + 2*(2*Lambdax*TYd02 + 4*ALambdax*Lambdax*Yd02)*(Yd00*Yu20 + Yd01*
      Yu21 + Yd02*Yu22) + 2*(2*Lambdax*TYd12 + 4*ALambdax*Lambdax*Yd12)*(Yd10*Yu20
      + Yd11*Yu21 + Yd12*Yu22) + 2*(2*Lambdax*TYd22 + 4*ALambdax*Lambdax*Yd22)*(
      Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22) + (2*Lambdax*TYu02 + 4*ALambdax*Lambdax*
      Yu02)*(6*Yu02*Yu22 + 4*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22)) + (2*Lambdax*
      TYu12 + 4*ALambdax*Lambdax*Yu12)*(6*Yu12*Yu22 + 4*(Yu10*Yu20 + Yu11*Yu21 +
      Yu12*Yu22)) + (2*Lambdax*TYu20 + 4*ALambdax*Lambdax*Yu20)*(Yd00*Yd02 + Yd10*
      Yd12 + Yd20*Yd22 + 6*Yu20*Yu22 + 5*(Yu00*Yu02 + Yu10*Yu12 + Yu20*Yu22)) + (2
      *Lambdax*TYu21 + 4*ALambdax*Lambdax*Yu21)*(Yd01*Yd02 + Yd11*Yd12 + Yd21*Yd22
      + 6*Yu21*Yu22 + 5*(Yu01*Yu02 + Yu11*Yu12 + Yu21*Yu22)) + 2*Lambdax*Yu22*(2*
      (TYd02*Yd02 + TYd12*Yd12 + TYd22*Yd22) + 12*TYu22*Yu22 + 6*(TYu00*Yu00 +
      TYu01*Yu01 + TYu02*Yu02 + TYu10*Yu10 + TYu11*Yu11 + TYu12*Yu12 + TYu20*Yu20
      + TYu21*Yu21 + TYu22*Yu22) + 4*(TYu02*Yu02 + TYu12*Yu12 + 2*TYu22*Yu22) + 5*
      (TYu20*Yu20 + TYu21*Yu21 + 2*TYu22*Yu22) + 1.7333333333333334*MassB*Sqr(g1)
      + 6*MassWB*Sqr(g2) + 10.666666666666666*MassG*Sqr(g3) + 2*ALambdax*Sqr(
      Lambdax) + 4*MassBp*Sqr(gN)*Sqr(QH2p) + 4*MassBp*Sqr(gN)*Sqr(QQp) + 4*MassBp
      *Sqr(gN)*Sqr(Qup)) + (2*Lambdax*TYu22 + 2*ALambdax*Lambdax*Yu22)*(-0.6*Sqr(
      g1) - 3*Sqr(g2) + 3*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(
      Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(
      Kappa22)) + 2*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(
      Lambda1211)) + 12*Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QH2p) -
      2*Sqr(gN)*Sqr(QSp) + 3*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr
      (Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(Ye00) + Sqr(
      Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21
      ) + Sqr(Ye22) + 3*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11)
      + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + (2*Lambdax*TYu22 + 4*
      ALambdax*Lambdax*Yu22)*(-0.8666666666666667*Sqr(g1) - 3*Sqr(g2) -
      5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*
      Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + Sqr(Yd02) + Sqr(Yd12) + Sqr(Yd22) + 6*Sqr(
      Yu22) + 5*(Sqr(Yu02) + Sqr(Yu12) + Sqr(Yu22)) + 4*(Sqr(Yu20) + Sqr(Yu21) +
      Sqr(Yu22)) + 3*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) +
      Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + (2*TYu22 + 2*ALambdax*Yu22
      )*(4*Power(Lambdax,3) - 0.6*Lambdax*Sqr(g1) - 3*Lambdax*Sqr(g2) + 3*Lambdax*
      (Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) +
      Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*Lambdax*(Sqr(
      Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) - 2*
      Lambdax*Sqr(gN)*Sqr(QH1p) - 2*Lambdax*Sqr(gN)*Sqr(QH2p) - 2*Lambdax*Sqr(gN)*
      Sqr(QSp) + 3*Lambdax*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(
      Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Lambdax*(Sqr(Ye00)
      + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) +
      Sqr(Ye21) + Sqr(Ye22)) + 3*Lambdax*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(
      Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + 2*
      Lambdax*(Yd02*(TYu20*Yd00 + TYu21*Yd01 + TYu22*Yd02) + Yd12*(TYu20*Yd10 +
      TYu21*Yd11 + TYu22*Yd12) + Yd22*(TYu20*Yd20 + TYu21*Yd21 + TYu22*Yd22) + 5*(
      Yu02*(TYu20*Yu00 + TYu21*Yu01 + TYu22*Yu02) + Yu12*(TYu20*Yu10 + TYu21*Yu11
      + TYu22*Yu12) + Yu22*(TYu20*Yu20 + TYu21*Yu21 + TYu22*Yu22)) + 2*(TYd02*(
      Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22) + TYd12*(Yd10*Yu20 + Yd11*Yu21 + Yd12*
      Yu22) + TYd22*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22)) - 0.8666666666666667*
      TYu22*Sqr(g1) - 3*TYu22*Sqr(g2) - 5.333333333333333*TYu22*Sqr(g3) + TYu22*
      Sqr(Lambdax) - 2*TYu22*Sqr(gN)*Sqr(QH2p) - 2*TYu22*Sqr(gN)*Sqr(QQp) - 2*
      TYu22*Sqr(gN)*Sqr(Qup) + Yu22*(6*(TYu00*Yu00 + TYu01*Yu01 + TYu02*Yu02 +
      TYu10*Yu10 + TYu11*Yu11 + TYu12*Yu12 + TYu20*Yu20 + TYu21*Yu21 + TYu22*Yu22)
      + 1.7333333333333334*MassB*Sqr(g1) + 6*MassWB*Sqr(g2) + 10.666666666666666*
      MassG*Sqr(g3) + 2*ALambdax*Sqr(Lambdax) + 4*MassBp*Sqr(gN)*Sqr(QH2p) + 4*
      MassBp*Sqr(gN)*Sqr(QQp) + 4*MassBp*Sqr(gN)*Sqr(Qup)) + 3*TYu22*(Sqr(Yu00) +
      Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(
      Yu21) + Sqr(Yu22)) + 4*(TYu02*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + TYu12*(
      Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + TYu22*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22
      )))) + 4*ALambdax*Lambdax*(Yd02*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22) + Yd12*(
      Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22) + Yd22*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22
      ) + 3*(Yu02*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu12*(Yu10*Yu20 + Yu11*
      Yu21 + Yu12*Yu22) + Yu22*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + Yu22*(
      -0.8666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(
      Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3
      *(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) +
      Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) + 2*Lambdax*Yu22*(6*(Kappa00*TKappa00 +
      Kappa01*TKappa01 + Kappa02*TKappa02 + Kappa10*TKappa10 + Kappa11*TKappa11 +
      Kappa12*TKappa12 + Kappa20*TKappa20 + Kappa21*TKappa21 + Kappa22*TKappa22)
      + 4*(Lambda1200*TLambda1200 + Lambda1201*TLambda1201 + Lambda1210*
      TLambda1210 + Lambda1211*TLambda1211) + 6*(TYd00*Yd00 + TYd01*Yd01 + TYd02*
      Yd02 + TYd10*Yd10 + TYd11*Yd11 + TYd12*Yd12 + TYd20*Yd20 + TYd21*Yd21 +
      TYd22*Yd22) + 2*(TYe00*Ye00 + TYe01*Ye01 + TYe02*Ye02 + TYe10*Ye10 + TYe11*
      Ye11 + TYe12*Ye12 + TYe20*Ye20 + TYe21*Ye21 + TYe22*Ye22) + 6*(TYu00*Yu00 +
      TYu01*Yu01 + TYu02*Yu02 + TYu10*Yu10 + TYu11*Yu11 + TYu12*Yu12 + TYu20*Yu20
      + TYu21*Yu21 + TYu22*Yu22) + 1.2*MassB*Sqr(g1) + 6*MassWB*Sqr(g2) + 24*
      ALambdax*Sqr(Lambdax) + 4*MassBp*Sqr(gN)*Sqr(QH1p) + 4*MassBp*Sqr(gN)*Sqr(
      QH2p) + 4*MassBp*Sqr(gN)*Sqr(QSp) + ALambdax*(-0.6*Sqr(g1) - 3*Sqr(g2) + 3*(
      Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) +
      Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*(Sqr(
      Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 12*Sqr(
      Lambdax) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QSp) +
      3*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) +
      Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr
      (Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22) + 3*(Sqr(
      Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20
      ) + Sqr(Yu21) + Sqr(Yu22)))) + 2*Yu22*(6*Lambdax*(Kappa00*TKappa00 + Kappa01
      *TKappa01 + Kappa02*TKappa02 + Kappa10*TKappa10 + Kappa11*TKappa11 + Kappa12
      *TKappa12 + Kappa20*TKappa20 + Kappa21*TKappa21 + Kappa22*TKappa22) + 4*
      Lambdax*(Lambda1200*TLambda1200 + Lambda1201*TLambda1201 + Lambda1210*
      TLambda1210 + Lambda1211*TLambda1211) + 6*Lambdax*(TYd00*Yd00 + TYd01*Yd01 +
      TYd02*Yd02 + TYd10*Yd10 + TYd11*Yd11 + TYd12*Yd12 + TYd20*Yd20 + TYd21*Yd21
      + TYd22*Yd22) + 2*Lambdax*(TYe00*Ye00 + TYe01*Ye01 + TYe02*Ye02 + TYe10*
      Ye10 + TYe11*Ye11 + TYe12*Ye12 + TYe20*Ye20 + TYe21*Ye21 + TYe22*Ye22) + 6*
      Lambdax*(TYu00*Yu00 + TYu01*Yu01 + TYu02*Yu02 + TYu10*Yu10 + TYu11*Yu11 +
      TYu12*Yu12 + TYu20*Yu20 + TYu21*Yu21 + TYu22*Yu22) + 1.2*Lambdax*MassB*Sqr(
      g1) + 6*Lambdax*MassWB*Sqr(g2) + 4*Lambdax*MassBp*Sqr(gN)*Sqr(QH1p) + 4*
      Lambdax*MassBp*Sqr(gN)*Sqr(QH2p) + 4*Lambdax*MassBp*Sqr(gN)*Sqr(QSp) +
      ALambdax*Lambdax*(-0.6*Sqr(g1) - 3*Sqr(g2) + 3*(Sqr(Kappa00) + Sqr(Kappa01)
      + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) +
      Sqr(Kappa21) + Sqr(Kappa22)) + 2*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(
      Lambda1210) + Sqr(Lambda1211)) + 12*Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH1p) - 2*
      Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QSp) + 3*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02
      ) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) +
      Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr
      (Ye20) + Sqr(Ye21) + Sqr(Ye22) + 3*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(
      Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dlead_log_one_loop_TYu22_dALambdax() const 
{
   const lowE6SSM_input_parameters inputs = model.get_input();

   const auto QQp = inputs.QQp;
   const auto QH1p = inputs.QH1p;
   const auto QH2p = inputs.QH2p;
   const auto Qup = inputs.Qup;
   const auto QSp = inputs.QSp;
   
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

   double deriv = 4*Yd02*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22)*Sqr(Lambdax) + 4*
      Yd12*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22)*Sqr(Lambdax) + 4*Yd22*(Yd20*Yu20 +
      Yd21*Yu21 + Yd22*Yu22)*Sqr(Lambdax) + 2*Yu02*(6*Yu02*Yu22 + 4*(Yu00*Yu20 +
      Yu01*Yu21 + Yu02*Yu22))*Sqr(Lambdax) + 2*Yu12*(6*Yu12*Yu22 + 4*(Yu10*Yu20 +
      Yu11*Yu21 + Yu12*Yu22))*Sqr(Lambdax) + 2*Yu20*(Yd00*Yd02 + Yd10*Yd12 + Yd20*
      Yd22 + 6*Yu20*Yu22 + 5*(Yu00*Yu02 + Yu10*Yu12 + Yu20*Yu22))*Sqr(Lambdax) + 2
      *Yu21*(Yd01*Yd02 + Yd11*Yd12 + Yd21*Yd22 + 6*Yu21*Yu22 + 5*(Yu01*Yu02 + Yu11
      *Yu12 + Yu21*Yu22))*Sqr(Lambdax) + 12*Yu22*Sqr(Lambdax)*Sqr(Yu00) + 12*Yu22*
      Sqr(Lambdax)*Sqr(Yu01) + 12*Yu22*Sqr(Lambdax)*Sqr(Yu10) + 12*Yu22*Sqr(
      Lambdax)*Sqr(Yu11) + 2*Yu22*Sqr(Lambdax)*(-0.6*Sqr(g1) - 3*Sqr(g2) + 3*(Sqr(
      Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(
      Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*(Sqr(Lambda1200)
      + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 12*Sqr(Lambdax) - 2
      *Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QSp) + 3*(Sqr(Yd00)
      + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) +
      Sqr(Yd21) + Sqr(Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr
      (Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22) + 3*(Sqr(Yu00) + Sqr(
      Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21
      ) + Sqr(Yu22))) + 2*Yu22*Sqr(Lambdax)*(-0.8666666666666667*Sqr(g1) - 3*Sqr(
      g2) - 5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr
      (gN)*Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + Sqr(Yd02) + Sqr(Yd12) + Sqr(Yd22) + 6*
      Sqr(Yu22) + 5*(Sqr(Yu02) + Sqr(Yu12) + Sqr(Yu22)) + 4*(Sqr(Yu20) + Sqr(Yu21)
      + Sqr(Yu22)) + 3*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11)
      + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + 2*Lambdax*Yu22*(4*Power
      (Lambdax,3) - 0.6*Lambdax*Sqr(g1) - 3*Lambdax*Sqr(g2) + 3*Lambdax*(Sqr(
      Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(
      Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*Lambdax*(Sqr(
      Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) - 2*
      Lambdax*Sqr(gN)*Sqr(QH1p) - 2*Lambdax*Sqr(gN)*Sqr(QH2p) - 2*Lambdax*Sqr(gN)*
      Sqr(QSp) + 3*Lambdax*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(
      Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Lambdax*(Sqr(Ye00)
      + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) +
      Sqr(Ye21) + Sqr(Ye22)) + 3*Lambdax*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(
      Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + 2*Sqr(
      Lambdax)*(Yd02*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22) + Yd12*(Yd10*Yu20 + Yd11*
      Yu21 + Yd12*Yu22) + Yd22*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22) + 3*(Yu02*(Yu00
      *Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu12*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) +
      Yu22*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + Yu22*(-0.8666666666666667*Sqr(g1
      ) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(
      QH2p) - 2*Sqr(gN)*Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3*(Sqr(Yu00) + Sqr(Yu01) +
      Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr
      (Yu22))));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dlead_log_one_loop_TYu22_dAYu22() const 
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

   double deriv = -16.64*Power(g1,4)*Yu22 - 24*Power(g2,4)*Yu22 + 12*Power(Yu22
      ,3)*Sqr(Lambdax) + Power(gN,3)*(9*Sqr(Qdp) + 9*Sqr(QDxbarp) + 9*Sqr(QDxp) +
      3*Sqr(Qep) + 6*Sqr(QH1p) + 6*Sqr(QH2p) + 2*Sqr(QHpbarp) + 2*Sqr(QHpp) + 6*
      Sqr(QLp) + 18*Sqr(QQp) + 3*Sqr(QSp) + 9*Sqr(Qup))*(-4*gN*Yu22*Sqr(QH2p) - 4*
      gN*Yu22*Sqr(QQp) - 4*gN*Yu22*Sqr(Qup)) + 2*Yd02*Yu22*(Yu02*(Yd00*Yu00 + Yd01
      *Yu01 + Yd02*Yu02) + Yu12*(Yd00*Yu10 + Yd01*Yu11 + Yd02*Yu12) + Yu22*(Yd00*
      Yu20 + Yd01*Yu21 + Yd02*Yu22) + 3*(Yd12*(Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12)
      + Yd22*(Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22) + Yd02*(Sqr(Yd00) + Sqr(Yd01) +
      Sqr(Yd02))) + Yd02*(-0.4666666666666667*Sqr(g1) - 3*Sqr(g2) -
      5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(Qdp) - 2*Sqr(gN)*
      Sqr(QH1p) - 2*Sqr(gN)*Sqr(QQp) + 3*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(
      Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(
      Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20
      ) + Sqr(Ye21) + Sqr(Ye22))) + 2*Yd12*Yu22*(Yu02*(Yd10*Yu00 + Yd11*Yu01 +
      Yd12*Yu02) + Yu12*(Yd10*Yu10 + Yd11*Yu11 + Yd12*Yu12) + Yu22*(Yd10*Yu20 +
      Yd11*Yu21 + Yd12*Yu22) + 3*(Yd02*(Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + Yd22*
      (Yd10*Yd20 + Yd11*Yd21 + Yd12*Yd22) + Yd12*(Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12
      ))) + Yd12*(-0.4666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(
      g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(Qdp) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*
      Sqr(QQp) + 3*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) +
      Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr
      (Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(
      Ye22))) + 2*Yd22*Yu22*(Yu02*(Yd20*Yu00 + Yd21*Yu01 + Yd22*Yu02) + Yu12*(Yd20
      *Yu10 + Yd21*Yu11 + Yd22*Yu12) + Yu22*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22) +
      3*(Yd02*(Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22) + Yd12*(Yd10*Yd20 + Yd11*Yd21 +
      Yd12*Yd22) + Yd22*(Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22))) + Yd22*(
      -0.4666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(
      Lambdax) - 2*Sqr(gN)*Sqr(Qdp) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QQp) + 3
      *(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) +
      Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr
      (Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22))) + 36*
      Power(Yu22,3)*Sqr(Yu00) + 36*Power(Yu22,3)*Sqr(Yu01) + 36*Power(Yu22,3)*Sqr(
      Yu10) + 36*Power(Yu22,3)*Sqr(Yu11) + (6*Yu02*Yu22 + 4*(Yu00*Yu20 + Yu01*Yu21
      + Yu02*Yu22))*(4*Yu22*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + 6*Yu02*Sqr(Yu22
      )) + (6*Yu12*Yu22 + 4*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22))*(4*Yu22*(Yu10*
      Yu20 + Yu11*Yu21 + Yu12*Yu22) + 6*Yu12*Sqr(Yu22)) + 2*Lambdax*Yu22*(4*Power(
      Lambdax,3) - 0.6*Lambdax*Sqr(g1) - 3*Lambdax*Sqr(g2) + 3*Lambdax*(Sqr(
      Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(
      Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*Lambdax*(Sqr(
      Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) - 2*
      Lambdax*Sqr(gN)*Sqr(QH1p) - 2*Lambdax*Sqr(gN)*Sqr(QH2p) - 2*Lambdax*Sqr(gN)*
      Sqr(QSp) + 3*Lambdax*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(
      Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Lambdax*(Sqr(Ye00)
      + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) +
      Sqr(Ye21) + Sqr(Ye22)) + 3*Lambdax*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(
      Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + (
      -0.8666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(
      Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) +
      Sqr(Yd02) + Sqr(Yd12) + Sqr(Yd22) + 6*Sqr(Yu22) + 5*(Sqr(Yu02) + Sqr(Yu12) +
      Sqr(Yu22)) + 4*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)) + 3*(Sqr(Yu00) + Sqr(
      Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21
      ) + Sqr(Yu22)))*(6*Power(Yu22,3) - 0.8666666666666667*Yu22*Sqr(g1) - 3*Yu22*
      Sqr(g2) - 5.333333333333333*Yu22*Sqr(g3) + Yu22*Sqr(Lambdax) - 2*Yu22*Sqr(gN
      )*Sqr(QH2p) - 2*Yu22*Sqr(gN)*Sqr(QQp) - 2*Yu22*Sqr(gN)*Sqr(Qup) + Yu22*Sqr(
      Yd02) + Yu22*Sqr(Yd12) + Yu22*Sqr(Yd22) + 5*(Power(Yu22,3) + Yu22*Sqr(Yu02)
      + Yu22*Sqr(Yu12)) + 4*Yu22*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)) + 3*Yu22*(Sqr
      (Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(
      Yu20) + Sqr(Yu21) + Sqr(Yu22))) + (Yd00*Yd02 + Yd10*Yd12 + Yd20*Yd22 + 6*
      Yu20*Yu22 + 5*(Yu00*Yu02 + Yu10*Yu12 + Yu20*Yu22))*(Yd00*Yd02*Yu22 + Yd10*
      Yd12*Yu22 + Yd20*Yd22*Yu22 + 6*Yu20*Sqr(Yu22) + 5*(Yu00*Yu02*Yu22 + Yu10*
      Yu12*Yu22 + Yu20*Sqr(Yu22))) + (Yd01*Yd02 + Yd11*Yd12 + Yd21*Yd22 + 6*Yu21*
      Yu22 + 5*(Yu01*Yu02 + Yu11*Yu12 + Yu21*Yu22))*(Yd01*Yd02*Yu22 + Yd11*Yd12*
      Yu22 + Yd21*Yd22*Yu22 + 6*Yu21*Sqr(Yu22) + 5*(Yu01*Yu02*Yu22 + Yu11*Yu12*
      Yu22 + Yu21*Sqr(Yu22))) + 6*Yu00*Yu22*(Yd00*(Yd00*Yu00 + Yd01*Yu01 + Yd02*
      Yu02) + Yd10*(Yd10*Yu00 + Yd11*Yu01 + Yd12*Yu02) + Yd20*(Yd20*Yu00 + Yd21*
      Yu01 + Yd22*Yu02) + 3*(Yu10*(Yu00*Yu10 + Yu01*Yu11 + Yu02*Yu12) + Yu20*(Yu00
      *Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu00*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02))) +
      Yu00*(-0.8666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) +
      Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup
      ) + 3*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12)
      + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) + 6*Yu01*Yu22*(Yd01*(Yd00*Yu00 +
      Yd01*Yu01 + Yd02*Yu02) + Yd11*(Yd10*Yu00 + Yd11*Yu01 + Yd12*Yu02) + Yd21*(
      Yd20*Yu00 + Yd21*Yu01 + Yd22*Yu02) + 3*(Yu11*(Yu00*Yu10 + Yu01*Yu11 + Yu02*
      Yu12) + Yu21*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu01*(Sqr(Yu00) + Sqr(
      Yu01) + Sqr(Yu02))) + Yu01*(-0.8666666666666667*Sqr(g1) - 3*Sqr(g2) -
      5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*
      Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(
      Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) + 16*
      Yu02*Yu22*(Yd02*(Yd00*Yu00 + Yd01*Yu01 + Yd02*Yu02) + Yd12*(Yd10*Yu00 + Yd11
      *Yu01 + Yd12*Yu02) + Yd22*(Yd20*Yu00 + Yd21*Yu01 + Yd22*Yu02) + 3*(Yu12*(
      Yu00*Yu10 + Yu01*Yu11 + Yu02*Yu12) + Yu22*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22
      ) + Yu02*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02))) + Yu02*(-0.8666666666666667*
      Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*
      Sqr(QH2p) - 2*Sqr(gN)*Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3*(Sqr(Yu00) + Sqr(
      Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21
      ) + Sqr(Yu22)))) + 6*Yu10*Yu22*(Yd00*(Yd00*Yu10 + Yd01*Yu11 + Yd02*Yu12) +
      Yd10*(Yd10*Yu10 + Yd11*Yu11 + Yd12*Yu12) + Yd20*(Yd20*Yu10 + Yd21*Yu11 +
      Yd22*Yu12) + 3*(Yu00*(Yu00*Yu10 + Yu01*Yu11 + Yu02*Yu12) + Yu20*(Yu10*Yu20 +
      Yu11*Yu21 + Yu12*Yu22) + Yu10*(Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12))) + Yu10*(
      -0.8666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(
      Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3
      *(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) +
      Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) + 6*Yu11*Yu22*(Yd01*(Yd00*Yu10 + Yd01*
      Yu11 + Yd02*Yu12) + Yd11*(Yd10*Yu10 + Yd11*Yu11 + Yd12*Yu12) + Yd21*(Yd20*
      Yu10 + Yd21*Yu11 + Yd22*Yu12) + 3*(Yu01*(Yu00*Yu10 + Yu01*Yu11 + Yu02*Yu12)
      + Yu21*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + Yu11*(Sqr(Yu10) + Sqr(Yu11) +
      Sqr(Yu12))) + Yu11*(-0.8666666666666667*Sqr(g1) - 3*Sqr(g2) -
      5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*
      Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(
      Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) + 16*
      Yu12*Yu22*(Yd02*(Yd00*Yu10 + Yd01*Yu11 + Yd02*Yu12) + Yd12*(Yd10*Yu10 + Yd11
      *Yu11 + Yd12*Yu12) + Yd22*(Yd20*Yu10 + Yd21*Yu11 + Yd22*Yu12) + 3*(Yu02*(
      Yu00*Yu10 + Yu01*Yu11 + Yu02*Yu12) + Yu22*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22
      ) + Yu12*(Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12))) + Yu12*(-0.8666666666666667*
      Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*
      Sqr(QH2p) - 2*Sqr(gN)*Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3*(Sqr(Yu00) + Sqr(
      Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21
      ) + Sqr(Yu22)))) + 14*Yu20*Yu22*(Yd00*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22) +
      Yd10*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22) + Yd20*(Yd20*Yu20 + Yd21*Yu21 +
      Yd22*Yu22) + 3*(Yu00*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu10*(Yu10*Yu20 +
      Yu11*Yu21 + Yu12*Yu22) + Yu20*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + Yu20*(
      -0.8666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(
      Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3
      *(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) +
      Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) + 14*Yu21*Yu22*(Yd01*(Yd00*Yu20 + Yd01*
      Yu21 + Yd02*Yu22) + Yd11*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22) + Yd21*(Yd20*
      Yu20 + Yd21*Yu21 + Yd22*Yu22) + 3*(Yu01*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22)
      + Yu11*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + Yu21*(Sqr(Yu20) + Sqr(Yu21) +
      Sqr(Yu22))) + Yu21*(-0.8666666666666667*Sqr(g1) - 3*Sqr(g2) -
      5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*
      Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(
      Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) + 36*
      Sqr(Yu22)*(Yd02*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22) + Yd12*(Yd10*Yu20 + Yd11
      *Yu21 + Yd12*Yu22) + Yd22*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22) + 3*(Yu02*(
      Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu12*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22
      ) + Yu22*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + Yu22*(-0.8666666666666667*
      Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*
      Sqr(QH2p) - 2*Sqr(gN)*Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3*(Sqr(Yu00) + Sqr(
      Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21
      ) + Sqr(Yu22)))) + 4*Yu22*Sqr(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22) + 4*Yu22*
      Sqr(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22) + 4*Yu22*Sqr(Yd20*Yu20 + Yd21*Yu21 +
      Yd22*Yu22);


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dlead_log_one_loop_TYu22_dMassB() const 
{
   const lowE6SSM_input_parameters inputs = model.get_input();

   const auto QQp = inputs.QQp;
   const auto QH2p = inputs.QH2p;
   const auto Qup = inputs.Qup;
   
   const double Yd00 = model.get_Yd(0,0);
   const double Yd01 = model.get_Yd(0,1);
   const double Yd02 = model.get_Yd(0,2);
   const double Yd10 = model.get_Yd(1,0);
   const double Yd11 = model.get_Yd(1,1);
   const double Yd12 = model.get_Yd(1,2);
   const double Yd20 = model.get_Yd(2,0);
   const double Yd21 = model.get_Yd(2,1);
   const double Yd22 = model.get_Yd(2,2);
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

   double deriv = 66.56*Power(g1,4)*Yu22 + 1.8666666666666667*Yd02*(Yd00*Yu20 +
      Yd01*Yu21 + Yd02*Yu22)*Sqr(g1) + 1.8666666666666667*Yd12*(Yd10*Yu20 + Yd11*
      Yu21 + Yd12*Yu22)*Sqr(g1) + 1.8666666666666667*Yd22*(Yd20*Yu20 + Yd21*Yu21 +
      Yd22*Yu22)*Sqr(g1) + 1.7333333333333334*Yu02*(6*Yu02*Yu22 + 4*(Yu00*Yu20 +
      Yu01*Yu21 + Yu02*Yu22))*Sqr(g1) + 1.7333333333333334*Yu12*(6*Yu12*Yu22 + 4*(
      Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22))*Sqr(g1) + 1.7333333333333334*Yu20*(Yd00*
      Yd02 + Yd10*Yd12 + Yd20*Yd22 + 6*Yu20*Yu22 + 5*(Yu00*Yu02 + Yu10*Yu12 + Yu20
      *Yu22))*Sqr(g1) + 1.7333333333333334*Yu21*(Yd01*Yd02 + Yd11*Yd12 + Yd21*Yd22
      + 6*Yu21*Yu22 + 5*(Yu01*Yu02 + Yu11*Yu12 + Yu21*Yu22))*Sqr(g1) + 2.4*Yu22*
      Sqr(g1)*Sqr(Lambdax) + 10.4*Yu22*Sqr(g1)*Sqr(Yu00) + 10.4*Yu22*Sqr(g1)*Sqr(
      Yu01) + 10.4*Yu22*Sqr(g1)*Sqr(Yu10) + 10.4*Yu22*Sqr(g1)*Sqr(Yu11) +
      1.7333333333333334*Yu22*Sqr(g1)*(-0.8666666666666667*Sqr(g1) - 3*Sqr(g2) -
      5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*
      Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + Sqr(Yd02) + Sqr(Yd12) + Sqr(Yd22) + 6*Sqr(
      Yu22) + 5*(Sqr(Yu02) + Sqr(Yu12) + Sqr(Yu22)) + 4*(Sqr(Yu20) + Sqr(Yu21) +
      Sqr(Yu22)) + 3*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) +
      Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + 1.7333333333333334*Sqr(g1)
      *(Yd02*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22) + Yd12*(Yd10*Yu20 + Yd11*Yu21 +
      Yd12*Yu22) + Yd22*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22) + 3*(Yu02*(Yu00*Yu20 +
      Yu01*Yu21 + Yu02*Yu22) + Yu12*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + Yu22*(
      Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + Yu22*(-0.8666666666666667*Sqr(g1) - 3*
      Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2
      *Sqr(gN)*Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02
      ) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))))
      ;


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dlead_log_one_loop_TYu22_dMassWB() const 
{
   const lowE6SSM_input_parameters inputs = model.get_input();

   const auto QQp = inputs.QQp;
   const auto QH2p = inputs.QH2p;
   const auto Qup = inputs.Qup;
   
   const double Yd00 = model.get_Yd(0,0);
   const double Yd01 = model.get_Yd(0,1);
   const double Yd02 = model.get_Yd(0,2);
   const double Yd10 = model.get_Yd(1,0);
   const double Yd11 = model.get_Yd(1,1);
   const double Yd12 = model.get_Yd(1,2);
   const double Yd20 = model.get_Yd(2,0);
   const double Yd21 = model.get_Yd(2,1);
   const double Yd22 = model.get_Yd(2,2);
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

   double deriv = 96*Power(g2,4)*Yu22 + 12*Yd02*(Yd00*Yu20 + Yd01*Yu21 + Yd02*
      Yu22)*Sqr(g2) + 12*Yd12*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22)*Sqr(g2) + 12*
      Yd22*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22)*Sqr(g2) + 6*Yu02*(6*Yu02*Yu22 + 4*(
      Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22))*Sqr(g2) + 6*Yu12*(6*Yu12*Yu22 + 4*(Yu10*
      Yu20 + Yu11*Yu21 + Yu12*Yu22))*Sqr(g2) + 6*Yu20*(Yd00*Yd02 + Yd10*Yd12 +
      Yd20*Yd22 + 6*Yu20*Yu22 + 5*(Yu00*Yu02 + Yu10*Yu12 + Yu20*Yu22))*Sqr(g2) + 6
      *Yu21*(Yd01*Yd02 + Yd11*Yd12 + Yd21*Yd22 + 6*Yu21*Yu22 + 5*(Yu01*Yu02 + Yu11
      *Yu12 + Yu21*Yu22))*Sqr(g2) + 12*Yu22*Sqr(g2)*Sqr(Lambdax) + 36*Yu22*Sqr(g2)
      *Sqr(Yu00) + 36*Yu22*Sqr(g2)*Sqr(Yu01) + 36*Yu22*Sqr(g2)*Sqr(Yu10) + 36*Yu22
      *Sqr(g2)*Sqr(Yu11) + 6*Yu22*Sqr(g2)*(-0.8666666666666667*Sqr(g1) - 3*Sqr(g2)
      - 5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN
      )*Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + Sqr(Yd02) + Sqr(Yd12) + Sqr(Yd22) + 6*Sqr(
      Yu22) + 5*(Sqr(Yu02) + Sqr(Yu12) + Sqr(Yu22)) + 4*(Sqr(Yu20) + Sqr(Yu21) +
      Sqr(Yu22)) + 3*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) +
      Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + 6*Sqr(g2)*(Yd02*(Yd00*Yu20
      + Yd01*Yu21 + Yd02*Yu22) + Yd12*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22) + Yd22*
      (Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22) + 3*(Yu02*(Yu00*Yu20 + Yu01*Yu21 + Yu02*
      Yu22) + Yu12*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + Yu22*(Sqr(Yu20) + Sqr(
      Yu21) + Sqr(Yu22))) + Yu22*(-0.8666666666666667*Sqr(g1) - 3*Sqr(g2) -
      5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*
      Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(
      Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dlead_log_one_loop_TYu22_dMassG() const 
{
   const lowE6SSM_input_parameters inputs = model.get_input();

   const auto QQp = inputs.QQp;
   const auto QH2p = inputs.QH2p;
   const auto Qup = inputs.Qup;
   
   const double Yd00 = model.get_Yd(0,0);
   const double Yd01 = model.get_Yd(0,1);
   const double Yd02 = model.get_Yd(0,2);
   const double Yd10 = model.get_Yd(1,0);
   const double Yd11 = model.get_Yd(1,1);
   const double Yd12 = model.get_Yd(1,2);
   const double Yd20 = model.get_Yd(2,0);
   const double Yd21 = model.get_Yd(2,1);
   const double Yd22 = model.get_Yd(2,2);
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

   double deriv = 21.333333333333332*Yd02*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22)*
      Sqr(g3) + 21.333333333333332*Yd12*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22)*Sqr(g3
      ) + 21.333333333333332*Yd22*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22)*Sqr(g3) +
      10.666666666666666*Yu02*(6*Yu02*Yu22 + 4*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22)
      )*Sqr(g3) + 10.666666666666666*Yu12*(6*Yu12*Yu22 + 4*(Yu10*Yu20 + Yu11*Yu21
      + Yu12*Yu22))*Sqr(g3) + 10.666666666666666*Yu20*(Yd00*Yd02 + Yd10*Yd12 +
      Yd20*Yd22 + 6*Yu20*Yu22 + 5*(Yu00*Yu02 + Yu10*Yu12 + Yu20*Yu22))*Sqr(g3) +
      10.666666666666666*Yu21*(Yd01*Yd02 + Yd11*Yd12 + Yd21*Yd22 + 6*Yu21*Yu22 + 5
      *(Yu01*Yu02 + Yu11*Yu12 + Yu21*Yu22))*Sqr(g3) + 64*Yu22*Sqr(g3)*Sqr(Yu00) +
      64*Yu22*Sqr(g3)*Sqr(Yu01) + 64*Yu22*Sqr(g3)*Sqr(Yu10) + 64*Yu22*Sqr(g3)*Sqr(
      Yu11) + 10.666666666666666*Yu22*Sqr(g3)*(-0.8666666666666667*Sqr(g1) - 3*Sqr
      (g2) - 5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*
      Sqr(gN)*Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + Sqr(Yd02) + Sqr(Yd12) + Sqr(Yd22) +
      6*Sqr(Yu22) + 5*(Sqr(Yu02) + Sqr(Yu12) + Sqr(Yu22)) + 4*(Sqr(Yu20) + Sqr(
      Yu21) + Sqr(Yu22)) + 3*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(
      Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + 10.666666666666666
      *Sqr(g3)*(Yd02*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22) + Yd12*(Yd10*Yu20 + Yd11*
      Yu21 + Yd12*Yu22) + Yd22*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22) + 3*(Yu02*(Yu00
      *Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu12*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) +
      Yu22*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + Yu22*(-0.8666666666666667*Sqr(g1
      ) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(
      QH2p) - 2*Sqr(gN)*Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3*(Sqr(Yu00) + Sqr(Yu01) +
      Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr
      (Yu22))));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dlead_log_one_loop_TYu22_dMassBp() const 
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

   double deriv = 2*Yd02*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22)*(4*Sqr(gN)*Sqr(Qdp
      ) + 4*Sqr(gN)*Sqr(QH1p) + 4*Sqr(gN)*Sqr(QQp)) + 2*Yd12*(Yd10*Yu20 + Yd11*
      Yu21 + Yd12*Yu22)*(4*Sqr(gN)*Sqr(Qdp) + 4*Sqr(gN)*Sqr(QH1p) + 4*Sqr(gN)*Sqr(
      QQp)) + 2*Yd22*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22)*(4*Sqr(gN)*Sqr(Qdp) + 4*
      Sqr(gN)*Sqr(QH1p) + 4*Sqr(gN)*Sqr(QQp)) + 2*Lambdax*Yu22*(4*Lambdax*Sqr(gN)*
      Sqr(QH1p) + 4*Lambdax*Sqr(gN)*Sqr(QH2p) + 4*Lambdax*Sqr(gN)*Sqr(QSp)) +
      Power(gN,3)*Yu22*(9*Sqr(Qdp) + 9*Sqr(QDxbarp) + 9*Sqr(QDxp) + 3*Sqr(Qep) + 6
      *Sqr(QH1p) + 6*Sqr(QH2p) + 2*Sqr(QHpbarp) + 2*Sqr(QHpp) + 6*Sqr(QLp) + 18*
      Sqr(QQp) + 3*Sqr(QSp) + 9*Sqr(Qup))*(8*gN*Sqr(QH2p) + 8*gN*Sqr(QQp) + 8*gN*
      Sqr(Qup)) + Yu02*(6*Yu02*Yu22 + 4*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22))*(4*
      Sqr(gN)*Sqr(QH2p) + 4*Sqr(gN)*Sqr(QQp) + 4*Sqr(gN)*Sqr(Qup)) + Yu12*(6*Yu12*
      Yu22 + 4*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22))*(4*Sqr(gN)*Sqr(QH2p) + 4*Sqr(
      gN)*Sqr(QQp) + 4*Sqr(gN)*Sqr(Qup)) + Yu20*(Yd00*Yd02 + Yd10*Yd12 + Yd20*Yd22
      + 6*Yu20*Yu22 + 5*(Yu00*Yu02 + Yu10*Yu12 + Yu20*Yu22))*(4*Sqr(gN)*Sqr(QH2p)
      + 4*Sqr(gN)*Sqr(QQp) + 4*Sqr(gN)*Sqr(Qup)) + Yu21*(Yd01*Yd02 + Yd11*Yd12 +
      Yd21*Yd22 + 6*Yu21*Yu22 + 5*(Yu01*Yu02 + Yu11*Yu12 + Yu21*Yu22))*(4*Sqr(gN)*
      Sqr(QH2p) + 4*Sqr(gN)*Sqr(QQp) + 4*Sqr(gN)*Sqr(Qup)) + 2*Yu22*Sqr(gN)*(9*Sqr
      (Qdp) + 9*Sqr(QDxbarp) + 9*Sqr(QDxp) + 3*Sqr(Qep) + 6*Sqr(QH1p) + 6*Sqr(QH2p
      ) + 2*Sqr(QHpbarp) + 2*Sqr(QHpp) + 6*Sqr(QLp) + 18*Sqr(QQp) + 3*Sqr(QSp) + 9
      *Sqr(Qup))*(4*Sqr(gN)*Sqr(QH2p) + 4*Sqr(gN)*Sqr(QQp) + 4*Sqr(gN)*Sqr(Qup)) +
      6*Yu22*(4*Sqr(gN)*Sqr(QH2p) + 4*Sqr(gN)*Sqr(QQp) + 4*Sqr(gN)*Sqr(Qup))*Sqr(
      Yu00) + 6*Yu22*(4*Sqr(gN)*Sqr(QH2p) + 4*Sqr(gN)*Sqr(QQp) + 4*Sqr(gN)*Sqr(Qup
      ))*Sqr(Yu01) + 6*Yu22*(4*Sqr(gN)*Sqr(QH2p) + 4*Sqr(gN)*Sqr(QQp) + 4*Sqr(gN)*
      Sqr(Qup))*Sqr(Yu10) + 6*Yu22*(4*Sqr(gN)*Sqr(QH2p) + 4*Sqr(gN)*Sqr(QQp) + 4*
      Sqr(gN)*Sqr(Qup))*Sqr(Yu11) + Yu22*(4*Sqr(gN)*Sqr(QH2p) + 4*Sqr(gN)*Sqr(QQp)
      + 4*Sqr(gN)*Sqr(Qup))*(-0.8666666666666667*Sqr(g1) - 3*Sqr(g2) -
      5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*
      Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + Sqr(Yd02) + Sqr(Yd12) + Sqr(Yd22) + 6*Sqr(
      Yu22) + 5*(Sqr(Yu02) + Sqr(Yu12) + Sqr(Yu22)) + 4*(Sqr(Yu20) + Sqr(Yu21) +
      Sqr(Yu22)) + 3*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) +
      Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + (4*Sqr(gN)*Sqr(QH2p) + 4*
      Sqr(gN)*Sqr(QQp) + 4*Sqr(gN)*Sqr(Qup))*(Yd02*(Yd00*Yu20 + Yd01*Yu21 + Yd02*
      Yu22) + Yd12*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22) + Yd22*(Yd20*Yu20 + Yd21*
      Yu21 + Yd22*Yu22) + 3*(Yu02*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu12*(Yu10
      *Yu20 + Yu11*Yu21 + Yu12*Yu22) + Yu22*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) +
      Yu22*(-0.8666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) +
      Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup
      ) + 3*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12)
      + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))));


   return twoLoop*deriv;
}

} // namespace flexiblesusy
