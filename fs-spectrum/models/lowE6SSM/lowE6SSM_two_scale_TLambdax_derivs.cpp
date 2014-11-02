#include "lowE6SSM_two_scale_tuning_calculator.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

// one loop
double lowE6SSM_tuning_calculator::deriv_dbeta_one_loop_TLambdax_dLambdax() const 
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
   const double g1 = model.get_g1();
   const double g2 = model.get_g2();
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
   const double MassBp = model.get_MassBp();

   double deriv = 6*(Kappa00*TKappa00 + Kappa01*TKappa01 + Kappa02*TKappa02 +
      Kappa10*TKappa10 + Kappa11*TKappa11 + Kappa12*TKappa12 + Kappa20*TKappa20 +
      Kappa21*TKappa21 + Kappa22*TKappa22) + 4*(Lambda1200*TLambda1200 +
      Lambda1201*TLambda1201 + Lambda1210*TLambda1210 + Lambda1211*TLambda1211) +
      6*(TYd00*Yd00 + TYd01*Yd01 + TYd02*Yd02 + TYd10*Yd10 + TYd11*Yd11 + TYd12*
      Yd12 + TYd20*Yd20 + TYd21*Yd21 + TYd22*Yd22) + 2*(TYe00*Ye00 + TYe01*Ye01 +
      TYe02*Ye02 + TYe10*Ye10 + TYe11*Ye11 + TYe12*Ye12 + TYe20*Ye20 + TYe21*Ye21
      + TYe22*Ye22) + 6*(TYu00*Yu00 + TYu01*Yu01 + TYu02*Yu02 + TYu10*Yu10 + TYu11
      *Yu11 + TYu12*Yu12 + TYu20*Yu20 + TYu21*Yu21 + TYu22*Yu22) + 1.2*MassB*Sqr(
      g1) + 6*MassWB*Sqr(g2) + 24*ALambdax*Sqr(Lambdax) + 4*MassBp*Sqr(gN)*Sqr(
      QH1p) + 4*MassBp*Sqr(gN)*Sqr(QH2p) + 4*MassBp*Sqr(gN)*Sqr(QSp) + ALambdax*(
      -0.6*Sqr(g1) - 3*Sqr(g2) + 3*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) +
      Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) +
      Sqr(Kappa22)) + 2*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr
      (Lambda1211)) + 12*Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QH2p)
      - 2*Sqr(gN)*Sqr(QSp) + 3*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) +
      Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(Ye00) + Sqr
      (Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(
      Ye21) + Sqr(Ye22) + 3*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(
      Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)));


   return oneOver16PiSqr*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_one_loop_TLambdax_dALambdax() const 
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
   const double g1 = model.get_g1();
   const double g2 = model.get_g2();
   const double gN = model.get_gN();

   double deriv = Lambdax*(-0.6*Sqr(g1) - 3*Sqr(g2) + 3*(Sqr(Kappa00) + Sqr(
      Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(
      Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*(Sqr(Lambda1200) + Sqr(
      Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 12*Sqr(Lambdax) - 2*Sqr(
      gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QSp) + 3*(Sqr(Yd00) +
      Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(
      Yd21) + Sqr(Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(
      Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22) + 3*(Sqr(Yu00) + Sqr(
      Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21
      ) + Sqr(Yu22)));


   return oneOver16PiSqr*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_one_loop_TLambdax_dAYu22() const 
{

   const double Lambdax = model.get_Lambdax();
   const double Yu22 = model.get_Yu(2,2);

   double deriv = 6*Lambdax*Sqr(Yu22);


   return oneOver16PiSqr*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_one_loop_TLambdax_dMassB() const 
{

   const double Lambdax = model.get_Lambdax();
   const double g1 = model.get_g1();

   double deriv = 1.2*Lambdax*Sqr(g1);


   return oneOver16PiSqr*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_one_loop_TLambdax_dMassWB() const 
{

   const double Lambdax = model.get_Lambdax();
   const double g2 = model.get_g2();

   double deriv = 6*Lambdax*Sqr(g2);


   return oneOver16PiSqr*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_one_loop_TLambdax_dMassG() const 
{

   double deriv = 0;


   return oneOver16PiSqr*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_one_loop_TLambdax_dMassBp() const 
{
   const lowE6SSM_input_parameters inputs = model.get_input();

   const auto QH1p = inputs.QH1p;
   const auto QH2p = inputs.QH2p;
   const auto QSp = inputs.QSp;

   const double Lambdax = model.get_Lambdax();
   const double gN = model.get_gN();

   double deriv = 4*Lambdax*Sqr(gN)*Sqr(QH1p) + 4*Lambdax*Sqr(gN)*Sqr(QH2p) + 4
      *Lambdax*Sqr(gN)*Sqr(QSp);


   return oneOver16PiSqr*deriv;
}

// two loop
double lowE6SSM_tuning_calculator::deriv_dbeta_two_loop_TLambdax_dLambdax() const 
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

   double deriv = -250*ALambdax*Power(Lambdax,4) - 0.2*Sqr(Lambdax)*(2*(30*(
      Kappa00*TKappa00 + Kappa01*TKappa01 + Kappa02*TKappa02 + Kappa10*TKappa10 +
      Kappa11*TKappa11 + Kappa12*TKappa12 + Kappa20*TKappa20 + Kappa21*TKappa21 +
      Kappa22*TKappa22) + 20*(Lambda1200*TLambda1200 + Lambda1201*TLambda1201 +
      Lambda1210*TLambda1210 + Lambda1211*TLambda1211) + 45*(TYd00*Yd00 + TYd01*
      Yd01 + TYd02*Yd02 + TYd10*Yd10 + TYd11*Yd11 + TYd12*Yd12 + TYd20*Yd20 +
      TYd21*Yd21 + TYd22*Yd22) + 15*(TYe00*Ye00 + TYe01*Ye01 + TYe02*Ye02 + TYe10*
      Ye10 + TYe11*Ye11 + TYe12*Ye12 + TYe20*Ye20 + TYe21*Ye21 + TYe22*Ye22) + 45*
      (TYu00*Yu00 + TYu01*Yu01 + TYu02*Yu02 + TYu10*Yu10 + TYu11*Yu11 + TYu12*Yu12
      + TYu20*Yu20 + TYu21*Yu21 + TYu22*Yu22) + 6*MassB*Sqr(g1) + 30*MassWB*Sqr(
      g2) + 20*MassBp*Sqr(gN)*Sqr(QH1p) + 20*MassBp*Sqr(gN)*Sqr(QH2p)) - 3*
      ALambdax*(6*Sqr(g1) + 30*Sqr(g2) - 30*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(
      Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(
      Kappa21) + Sqr(Kappa22)) - 20*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(
      Lambda1210) + Sqr(Lambda1211)) + 20*Sqr(gN)*Sqr(QH1p) + 20*Sqr(gN)*Sqr(QH2p)
      - 45*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12)
      + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) - 15*(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02
      ) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22)) -
      45*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) +
      Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) - 0.4*Lambdax*(2*Lambdax*(30*(Kappa00*
      TKappa00 + Kappa01*TKappa01 + Kappa02*TKappa02 + Kappa10*TKappa10 + Kappa11*
      TKappa11 + Kappa12*TKappa12 + Kappa20*TKappa20 + Kappa21*TKappa21 + Kappa22*
      TKappa22) + 20*(Lambda1200*TLambda1200 + Lambda1201*TLambda1201 + Lambda1210
      *TLambda1210 + Lambda1211*TLambda1211) + 45*(TYd00*Yd00 + TYd01*Yd01 + TYd02
      *Yd02 + TYd10*Yd10 + TYd11*Yd11 + TYd12*Yd12 + TYd20*Yd20 + TYd21*Yd21 +
      TYd22*Yd22) + 15*(TYe00*Ye00 + TYe01*Ye01 + TYe02*Ye02 + TYe10*Ye10 + TYe11*
      Ye11 + TYe12*Ye12 + TYe20*Ye20 + TYe21*Ye21 + TYe22*Ye22) + 45*(TYu00*Yu00 +
      TYu01*Yu01 + TYu02*Yu02 + TYu10*Yu10 + TYu11*Yu11 + TYu12*Yu12 + TYu20*Yu20
      + TYu21*Yu21 + TYu22*Yu22) + 6*MassB*Sqr(g1) + 30*MassWB*Sqr(g2) + 20*
      MassBp*Sqr(gN)*Sqr(QH1p) + 20*MassBp*Sqr(gN)*Sqr(QH2p)) - 3*ALambdax*Lambdax
      *(6*Sqr(g1) + 30*Sqr(g2) - 30*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) +
      Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) +
      Sqr(Kappa22)) - 20*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) +
      Sqr(Lambda1211)) + 20*Sqr(gN)*Sqr(QH1p) + 20*Sqr(gN)*Sqr(QH2p) - 45*(Sqr(
      Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20
      ) + Sqr(Yd21) + Sqr(Yd22)) - 15*(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(
      Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22)) - 45*(Sqr
      (Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(
      Yu20) + Sqr(Yu21) + Sqr(Yu22)))) - 0.08*(297*Power(g1,4)*MassB + 825*Power(
      g2,4)*MassWB + 800*Power(gN,4)*MassBp*Power(QH1p,4) + 800*Power(gN,4)*MassBp
      *Power(QH2p,4) + 500*Power(gN,4)*MassBp*Power(QSp,4) + 150*(Yu00*(TYd00*(
      Yd00*Yu00 + Yd01*Yu01 + Yd02*Yu02) + TYd10*(Yd10*Yu00 + Yd11*Yu01 + Yd12*
      Yu02) + TYd20*(Yd20*Yu00 + Yd21*Yu01 + Yd22*Yu02)) + Yu01*(TYd01*(Yd00*Yu00
      + Yd01*Yu01 + Yd02*Yu02) + TYd11*(Yd10*Yu00 + Yd11*Yu01 + Yd12*Yu02) + TYd21
      *(Yd20*Yu00 + Yd21*Yu01 + Yd22*Yu02)) + Yu02*(TYd02*(Yd00*Yu00 + Yd01*Yu01 +
      Yd02*Yu02) + TYd12*(Yd10*Yu00 + Yd11*Yu01 + Yd12*Yu02) + TYd22*(Yd20*Yu00 +
      Yd21*Yu01 + Yd22*Yu02)) + Yu10*(TYd00*(Yd00*Yu10 + Yd01*Yu11 + Yd02*Yu12) +
      TYd10*(Yd10*Yu10 + Yd11*Yu11 + Yd12*Yu12) + TYd20*(Yd20*Yu10 + Yd21*Yu11 +
      Yd22*Yu12)) + Yu11*(TYd01*(Yd00*Yu10 + Yd01*Yu11 + Yd02*Yu12) + TYd11*(Yd10*
      Yu10 + Yd11*Yu11 + Yd12*Yu12) + TYd21*(Yd20*Yu10 + Yd21*Yu11 + Yd22*Yu12)) +
      Yu12*(TYd02*(Yd00*Yu10 + Yd01*Yu11 + Yd02*Yu12) + TYd12*(Yd10*Yu10 + Yd11*
      Yu11 + Yd12*Yu12) + TYd22*(Yd20*Yu10 + Yd21*Yu11 + Yd22*Yu12)) + Yu20*(TYd00
      *(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22) + TYd10*(Yd10*Yu20 + Yd11*Yu21 + Yd12*
      Yu22) + TYd20*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22)) + Yu21*(TYd01*(Yd00*Yu20
      + Yd01*Yu21 + Yd02*Yu22) + TYd11*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22) + TYd21
      *(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22)) + Yu22*(TYd02*(Yd00*Yu20 + Yd01*Yu21 +
      Yd02*Yu22) + TYd12*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22) + TYd22*(Yd20*Yu20 +
      Yd21*Yu21 + Yd22*Yu22))) + 150*(Yd00*(TYu00*(Yd00*Yu00 + Yd01*Yu01 + Yd02*
      Yu02) + TYu10*(Yd00*Yu10 + Yd01*Yu11 + Yd02*Yu12) + TYu20*(Yd00*Yu20 + Yd01*
      Yu21 + Yd02*Yu22)) + Yd01*(TYu01*(Yd00*Yu00 + Yd01*Yu01 + Yd02*Yu02) + TYu11
      *(Yd00*Yu10 + Yd01*Yu11 + Yd02*Yu12) + TYu21*(Yd00*Yu20 + Yd01*Yu21 + Yd02*
      Yu22)) + Yd02*(TYu02*(Yd00*Yu00 + Yd01*Yu01 + Yd02*Yu02) + TYu12*(Yd00*Yu10
      + Yd01*Yu11 + Yd02*Yu12) + TYu22*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22)) + Yd10
      *(TYu00*(Yd10*Yu00 + Yd11*Yu01 + Yd12*Yu02) + TYu10*(Yd10*Yu10 + Yd11*Yu11 +
      Yd12*Yu12) + TYu20*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22)) + Yd11*(TYu01*(Yd10
      *Yu00 + Yd11*Yu01 + Yd12*Yu02) + TYu11*(Yd10*Yu10 + Yd11*Yu11 + Yd12*Yu12) +
      TYu21*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22)) + Yd12*(TYu02*(Yd10*Yu00 + Yd11*
      Yu01 + Yd12*Yu02) + TYu12*(Yd10*Yu10 + Yd11*Yu11 + Yd12*Yu12) + TYu22*(Yd10*
      Yu20 + Yd11*Yu21 + Yd12*Yu22)) + Yd20*(TYu00*(Yd20*Yu00 + Yd21*Yu01 + Yd22*
      Yu02) + TYu10*(Yd20*Yu10 + Yd21*Yu11 + Yd22*Yu12) + TYu20*(Yd20*Yu20 + Yd21*
      Yu21 + Yd22*Yu22)) + Yd21*(TYu01*(Yd20*Yu00 + Yd21*Yu01 + Yd22*Yu02) + TYu11
      *(Yd20*Yu10 + Yd21*Yu11 + Yd22*Yu12) + TYu21*(Yd20*Yu20 + Yd21*Yu21 + Yd22*
      Yu22)) + Yd22*(TYu02*(Yd20*Yu00 + Yd21*Yu01 + Yd22*Yu02) + TYu12*(Yd20*Yu10
      + Yd21*Yu11 + Yd22*Yu12) + TYu22*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22))) - 20*
      (Kappa00*TKappa00 + Kappa01*TKappa01 + Kappa02*TKappa02 + Kappa10*TKappa10 +
      Kappa11*TKappa11 + Kappa12*TKappa12 + Kappa20*TKappa20 + Kappa21*TKappa21 +
      Kappa22*TKappa22)*Sqr(g1) - 30*(Lambda1200*TLambda1200 + Lambda1201*
      TLambda1201 + Lambda1210*TLambda1210 + Lambda1211*TLambda1211)*Sqr(g1) + 10*
      (TYd00*Yd00 + TYd01*Yd01 + TYd02*Yd02 + TYd10*Yd10 + TYd11*Yd11 + TYd12*Yd12
      + TYd20*Yd20 + TYd21*Yd21 + TYd22*Yd22)*Sqr(g1) - 30*(TYe00*Ye00 + TYe01*
      Ye01 + TYe02*Ye02 + TYe10*Ye10 + TYe11*Ye11 + TYe12*Ye12 + TYe20*Ye20 +
      TYe21*Ye21 + TYe22*Ye22)*Sqr(g1) - 20*(TYu00*Yu00 + TYu01*Yu01 + TYu02*Yu02
      + TYu10*Yu10 + TYu11*Yu11 + TYu12*Yu12 + TYu20*Yu20 + TYu21*Yu21 + TYu22*
      Yu22)*Sqr(g1) - 150*(Lambda1200*TLambda1200 + Lambda1201*TLambda1201 +
      Lambda1210*TLambda1210 + Lambda1211*TLambda1211)*Sqr(g2) + 45*MassB*Sqr(g1)*
      Sqr(g2) + 45*MassWB*Sqr(g1)*Sqr(g2) - 400*(Kappa00*TKappa00 + Kappa01*
      TKappa01 + Kappa02*TKappa02 + Kappa10*TKappa10 + Kappa11*TKappa11 + Kappa12*
      TKappa12 + Kappa20*TKappa20 + Kappa21*TKappa21 + Kappa22*TKappa22)*Sqr(g3) -
      400*(TYd00*Yd00 + TYd01*Yd01 + TYd02*Yd02 + TYd10*Yd10 + TYd11*Yd11 + TYd12
      *Yd12 + TYd20*Yd20 + TYd21*Yd21 + TYd22*Yd22)*Sqr(g3) - 400*(TYu00*Yu00 +
      TYu01*Yu01 + TYu02*Yu02 + TYu10*Yu10 + TYu11*Yu11 + TYu12*Yu12 + TYu20*Yu20
      + TYu21*Yu21 + TYu22*Yu22)*Sqr(g3) + 20*MassB*Sqr(g1)*(Sqr(Kappa00) + Sqr(
      Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(
      Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 400*MassG*Sqr(g3)*(Sqr(Kappa00) +
      Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) +
      Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 300*(Kappa00*((Kappa00*Kappa10
      + Kappa01*Kappa11 + Kappa02*Kappa12)*TKappa10 + (Kappa00*Kappa20 + Kappa01*
      Kappa21 + Kappa02*Kappa22)*TKappa20 + TKappa00*(Sqr(Kappa00) + Sqr(Kappa01)
      + Sqr(Kappa02))) + Kappa01*((Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*
      Kappa12)*TKappa11 + (Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)*
      TKappa21 + TKappa01*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02))) + Kappa02*
      ((Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12)*TKappa12 + (Kappa00*
      Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)*TKappa22 + TKappa02*(Sqr(
      Kappa00) + Sqr(Kappa01) + Sqr(Kappa02))) + Kappa10*((Kappa00*Kappa10 +
      Kappa01*Kappa11 + Kappa02*Kappa12)*TKappa00 + (Kappa10*Kappa20 + Kappa11*
      Kappa21 + Kappa12*Kappa22)*TKappa20 + TKappa10*(Sqr(Kappa10) + Sqr(Kappa11)
      + Sqr(Kappa12))) + Kappa11*((Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*
      Kappa12)*TKappa01 + (Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)*
      TKappa21 + TKappa11*(Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12))) + Kappa12*
      ((Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12)*TKappa02 + (Kappa10*
      Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)*TKappa22 + TKappa12*(Sqr(
      Kappa10) + Sqr(Kappa11) + Sqr(Kappa12))) + Kappa20*((Kappa00*Kappa20 +
      Kappa01*Kappa21 + Kappa02*Kappa22)*TKappa00 + (Kappa10*Kappa20 + Kappa11*
      Kappa21 + Kappa12*Kappa22)*TKappa10 + TKappa20*(Sqr(Kappa20) + Sqr(Kappa21)
      + Sqr(Kappa22))) + Kappa21*((Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*
      Kappa22)*TKappa01 + (Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)*
      TKappa11 + TKappa21*(Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22))) + Kappa22*
      ((Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)*TKappa02 + (Kappa10*
      Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)*TKappa12 + TKappa22*(Sqr(
      Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)))) + 30*MassB*Sqr(g1)*(Sqr(Lambda1200
      ) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 150*MassWB*Sqr(g2
      )*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) +
      200*(Lambda1200*((Lambda1200*Lambda1210 + Lambda1201*Lambda1211)*TLambda1210
      + TLambda1200*(Sqr(Lambda1200) + Sqr(Lambda1201))) + Lambda1201*((
      Lambda1200*Lambda1210 + Lambda1201*Lambda1211)*TLambda1211 + TLambda1201*(
      Sqr(Lambda1200) + Sqr(Lambda1201))) + Lambda1210*((Lambda1200*Lambda1210 +
      Lambda1201*Lambda1211)*TLambda1200 + TLambda1210*(Sqr(Lambda1210) + Sqr(
      Lambda1211))) + Lambda1211*((Lambda1200*Lambda1210 + Lambda1201*Lambda1211)*
      TLambda1201 + TLambda1211*(Sqr(Lambda1210) + Sqr(Lambda1211)))) - 150*(TYd00
      *Yd00 + TYd01*Yd01 + TYd02*Yd02 + TYd10*Yd10 + TYd11*Yd11 + TYd12*Yd12 +
      TYd20*Yd20 + TYd21*Yd21 + TYd22*Yd22)*Sqr(gN)*Sqr(Qdp) - 150*(Kappa00*
      TKappa00 + Kappa01*TKappa01 + Kappa02*TKappa02 + Kappa10*TKappa10 + Kappa11*
      TKappa11 + Kappa12*TKappa12 + Kappa20*TKappa20 + Kappa21*TKappa21 + Kappa22*
      TKappa22)*Sqr(gN)*Sqr(QDxbarp) + 150*MassBp*Sqr(gN)*(Sqr(Kappa00) + Sqr(
      Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(
      Kappa20) + Sqr(Kappa21) + Sqr(Kappa22))*Sqr(QDxbarp) - 150*(Kappa00*TKappa00
      + Kappa01*TKappa01 + Kappa02*TKappa02 + Kappa10*TKappa10 + Kappa11*TKappa11
      + Kappa12*TKappa12 + Kappa20*TKappa20 + Kappa21*TKappa21 + Kappa22*TKappa22
      )*Sqr(gN)*Sqr(QDxp) + 150*MassBp*Sqr(gN)*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(
      Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(
      Kappa21) + Sqr(Kappa22))*Sqr(QDxp) - 50*(TYe00*Ye00 + TYe01*Ye01 + TYe02*
      Ye02 + TYe10*Ye10 + TYe11*Ye11 + TYe12*Ye12 + TYe20*Ye20 + TYe21*Ye21 +
      TYe22*Ye22)*Sqr(gN)*Sqr(Qep) - 100*(Lambda1200*TLambda1200 + Lambda1201*
      TLambda1201 + Lambda1210*TLambda1210 + Lambda1211*TLambda1211)*Sqr(gN)*Sqr(
      QH1p) + 150*(TYd00*Yd00 + TYd01*Yd01 + TYd02*Yd02 + TYd10*Yd10 + TYd11*Yd11
      + TYd12*Yd12 + TYd20*Yd20 + TYd21*Yd21 + TYd22*Yd22)*Sqr(gN)*Sqr(QH1p) + 50*
      (TYe00*Ye00 + TYe01*Ye01 + TYe02*Ye02 + TYe10*Ye10 + TYe11*Ye11 + TYe12*Ye12
      + TYe20*Ye20 + TYe21*Ye21 + TYe22*Ye22)*Sqr(gN)*Sqr(QH1p) + 30*MassB*Sqr(g1
      )*Sqr(gN)*Sqr(QH1p) + 30*MassBp*Sqr(g1)*Sqr(gN)*Sqr(QH1p) + 150*MassBp*Sqr(
      g2)*Sqr(gN)*Sqr(QH1p) + 150*MassWB*Sqr(g2)*Sqr(gN)*Sqr(QH1p) + 100*MassBp*
      Sqr(gN)*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(
      Lambda1211))*Sqr(QH1p) + 900*Power(gN,4)*MassBp*Sqr(Qdp)*Sqr(QH1p) + 900*
      Power(gN,4)*MassBp*Sqr(QDxbarp)*Sqr(QH1p) + 900*Power(gN,4)*MassBp*Sqr(QDxp)
      *Sqr(QH1p) + 300*Power(gN,4)*MassBp*Sqr(Qep)*Sqr(QH1p) - 100*(Lambda1200*
      TLambda1200 + Lambda1201*TLambda1201 + Lambda1210*TLambda1210 + Lambda1211*
      TLambda1211)*Sqr(gN)*Sqr(QH2p) + 150*(TYu00*Yu00 + TYu01*Yu01 + TYu02*Yu02 +
      TYu10*Yu10 + TYu11*Yu11 + TYu12*Yu12 + TYu20*Yu20 + TYu21*Yu21 + TYu22*Yu22
      )*Sqr(gN)*Sqr(QH2p) + 30*MassB*Sqr(g1)*Sqr(gN)*Sqr(QH2p) + 30*MassBp*Sqr(g1)
      *Sqr(gN)*Sqr(QH2p) + 150*MassBp*Sqr(g2)*Sqr(gN)*Sqr(QH2p) + 150*MassWB*Sqr(
      g2)*Sqr(gN)*Sqr(QH2p) + 100*MassBp*Sqr(gN)*(Sqr(Lambda1200) + Sqr(Lambda1201
      ) + Sqr(Lambda1210) + Sqr(Lambda1211))*Sqr(QH2p) + 900*Power(gN,4)*MassBp*
      Sqr(Qdp)*Sqr(QH2p) + 900*Power(gN,4)*MassBp*Sqr(QDxbarp)*Sqr(QH2p) + 900*
      Power(gN,4)*MassBp*Sqr(QDxp)*Sqr(QH2p) + 300*Power(gN,4)*MassBp*Sqr(Qep)*Sqr
      (QH2p) + 1200*Power(gN,4)*MassBp*Sqr(QH1p)*Sqr(QH2p) + 200*Power(gN,4)*
      MassBp*Sqr(QH1p)*Sqr(QHpbarp) + 200*Power(gN,4)*MassBp*Sqr(QH2p)*Sqr(QHpbarp
      ) + 200*Power(gN,4)*MassBp*Sqr(QH1p)*Sqr(QHpp) + 200*Power(gN,4)*MassBp*Sqr(
      QH2p)*Sqr(QHpp) - 50*(TYe00*Ye00 + TYe01*Ye01 + TYe02*Ye02 + TYe10*Ye10 +
      TYe11*Ye11 + TYe12*Ye12 + TYe20*Ye20 + TYe21*Ye21 + TYe22*Ye22)*Sqr(gN)*Sqr(
      QLp) + 600*Power(gN,4)*MassBp*Sqr(QH1p)*Sqr(QLp) + 600*Power(gN,4)*MassBp*
      Sqr(QH2p)*Sqr(QLp) - 150*(TYd00*Yd00 + TYd01*Yd01 + TYd02*Yd02 + TYd10*Yd10
      + TYd11*Yd11 + TYd12*Yd12 + TYd20*Yd20 + TYd21*Yd21 + TYd22*Yd22)*Sqr(gN)*
      Sqr(QQp) - 150*(TYu00*Yu00 + TYu01*Yu01 + TYu02*Yu02 + TYu10*Yu10 + TYu11*
      Yu11 + TYu12*Yu12 + TYu20*Yu20 + TYu21*Yu21 + TYu22*Yu22)*Sqr(gN)*Sqr(QQp) +
      1800*Power(gN,4)*MassBp*Sqr(QH1p)*Sqr(QQp) + 1800*Power(gN,4)*MassBp*Sqr(
      QH2p)*Sqr(QQp) + 150*(Kappa00*TKappa00 + Kappa01*TKappa01 + Kappa02*TKappa02
      + Kappa10*TKappa10 + Kappa11*TKappa11 + Kappa12*TKappa12 + Kappa20*TKappa20
      + Kappa21*TKappa21 + Kappa22*TKappa22)*Sqr(gN)*Sqr(QSp) + 100*(Lambda1200*
      TLambda1200 + Lambda1201*TLambda1201 + Lambda1210*TLambda1210 + Lambda1211*
      TLambda1211)*Sqr(gN)*Sqr(QSp) - 150*MassBp*Sqr(gN)*(Sqr(Kappa00) + Sqr(
      Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(
      Kappa20) + Sqr(Kappa21) + Sqr(Kappa22))*Sqr(QSp) - 100*MassBp*Sqr(gN)*(Sqr(
      Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211))*Sqr(QSp)
      + 900*Power(gN,4)*MassBp*Sqr(Qdp)*Sqr(QSp) + 900*Power(gN,4)*MassBp*Sqr(
      QDxbarp)*Sqr(QSp) + 900*Power(gN,4)*MassBp*Sqr(QDxp)*Sqr(QSp) + 300*Power(gN
      ,4)*MassBp*Sqr(Qep)*Sqr(QSp) + 900*Power(gN,4)*MassBp*Sqr(QH1p)*Sqr(QSp) +
      900*Power(gN,4)*MassBp*Sqr(QH2p)*Sqr(QSp) + 200*Power(gN,4)*MassBp*Sqr(
      QHpbarp)*Sqr(QSp) + 200*Power(gN,4)*MassBp*Sqr(QHpp)*Sqr(QSp) + 600*Power(gN
      ,4)*MassBp*Sqr(QLp)*Sqr(QSp) + 1800*Power(gN,4)*MassBp*Sqr(QQp)*Sqr(QSp) -
      150*(TYu00*Yu00 + TYu01*Yu01 + TYu02*Yu02 + TYu10*Yu10 + TYu11*Yu11 + TYu12*
      Yu12 + TYu20*Yu20 + TYu21*Yu21 + TYu22*Yu22)*Sqr(gN)*Sqr(Qup) + 900*Power(gN
      ,4)*MassBp*Sqr(QH1p)*Sqr(Qup) + 900*Power(gN,4)*MassBp*Sqr(QH2p)*Sqr(Qup) +
      900*Power(gN,4)*MassBp*Sqr(QSp)*Sqr(Qup) - 10*(MassB*Sqr(g1) - 5*(8*MassG*
      Sqr(g3) + 3*MassBp*Sqr(gN)*(Sqr(Qdp) - Sqr(QH1p) + Sqr(QQp))))*(Sqr(Yd00) +
      Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(
      Yd21) + Sqr(Yd22)) + 450*(Yd00*(TYd10*(Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) +
      TYd20*(Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22) + TYd00*(Sqr(Yd00) + Sqr(Yd01) +
      Sqr(Yd02))) + Yd01*(TYd11*(Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + TYd21*(Yd00*
      Yd20 + Yd01*Yd21 + Yd02*Yd22) + TYd01*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02))) +
      Yd02*(TYd12*(Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + TYd22*(Yd00*Yd20 + Yd01*
      Yd21 + Yd02*Yd22) + TYd02*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02))) + Yd10*(TYd00
      *(Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + TYd20*(Yd10*Yd20 + Yd11*Yd21 + Yd12*
      Yd22) + TYd10*(Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12))) + Yd11*(TYd01*(Yd00*Yd10
      + Yd01*Yd11 + Yd02*Yd12) + TYd21*(Yd10*Yd20 + Yd11*Yd21 + Yd12*Yd22) + TYd11
      *(Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12))) + Yd12*(TYd02*(Yd00*Yd10 + Yd01*Yd11 +
      Yd02*Yd12) + TYd22*(Yd10*Yd20 + Yd11*Yd21 + Yd12*Yd22) + TYd12*(Sqr(Yd10) +
      Sqr(Yd11) + Sqr(Yd12))) + Yd20*(TYd00*(Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22) +
      TYd10*(Yd10*Yd20 + Yd11*Yd21 + Yd12*Yd22) + TYd20*(Sqr(Yd20) + Sqr(Yd21) +
      Sqr(Yd22))) + Yd21*(TYd01*(Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22) + TYd11*(Yd10*
      Yd20 + Yd11*Yd21 + Yd12*Yd22) + TYd21*(Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22))) +
      Yd22*(TYd02*(Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22) + TYd12*(Yd10*Yd20 + Yd11*
      Yd21 + Yd12*Yd22) + TYd22*(Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)))) + 10*(3*
      MassB*Sqr(g1) + 5*MassBp*Sqr(gN)*(Sqr(Qep) - Sqr(QH1p) + Sqr(QLp)))*(Sqr(
      Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20
      ) + Sqr(Ye21) + Sqr(Ye22)) + 150*(Ye00*(TYe10*(Ye00*Ye10 + Ye01*Ye11 + Ye02*
      Ye12) + TYe20*(Ye00*Ye20 + Ye01*Ye21 + Ye02*Ye22) + TYe00*(Sqr(Ye00) + Sqr(
      Ye01) + Sqr(Ye02))) + Ye01*(TYe11*(Ye00*Ye10 + Ye01*Ye11 + Ye02*Ye12) +
      TYe21*(Ye00*Ye20 + Ye01*Ye21 + Ye02*Ye22) + TYe01*(Sqr(Ye00) + Sqr(Ye01) +
      Sqr(Ye02))) + Ye02*(TYe12*(Ye00*Ye10 + Ye01*Ye11 + Ye02*Ye12) + TYe22*(Ye00*
      Ye20 + Ye01*Ye21 + Ye02*Ye22) + TYe02*(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02))) +
      Ye10*(TYe00*(Ye00*Ye10 + Ye01*Ye11 + Ye02*Ye12) + TYe20*(Ye10*Ye20 + Ye11*
      Ye21 + Ye12*Ye22) + TYe10*(Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12))) + Ye11*(TYe01
      *(Ye00*Ye10 + Ye01*Ye11 + Ye02*Ye12) + TYe21*(Ye10*Ye20 + Ye11*Ye21 + Ye12*
      Ye22) + TYe11*(Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12))) + Ye12*(TYe02*(Ye00*Ye10
      + Ye01*Ye11 + Ye02*Ye12) + TYe22*(Ye10*Ye20 + Ye11*Ye21 + Ye12*Ye22) + TYe12
      *(Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12))) + Ye20*(TYe00*(Ye00*Ye20 + Ye01*Ye21 +
      Ye02*Ye22) + TYe10*(Ye10*Ye20 + Ye11*Ye21 + Ye12*Ye22) + TYe20*(Sqr(Ye20) +
      Sqr(Ye21) + Sqr(Ye22))) + Ye21*(TYe01*(Ye00*Ye20 + Ye01*Ye21 + Ye02*Ye22) +
      TYe11*(Ye10*Ye20 + Ye11*Ye21 + Ye12*Ye22) + TYe21*(Sqr(Ye20) + Sqr(Ye21) +
      Sqr(Ye22))) + Ye22*(TYe02*(Ye00*Ye20 + Ye01*Ye21 + Ye02*Ye22) + TYe12*(Ye10*
      Ye20 + Ye11*Ye21 + Ye12*Ye22) + TYe22*(Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22))))
      + 20*MassB*Sqr(g1)*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11
      ) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)) + 400*MassG*Sqr(g3)*(Sqr(
      Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20
      ) + Sqr(Yu21) + Sqr(Yu22)) - 150*MassBp*Sqr(gN)*Sqr(QH2p)*(Sqr(Yu00) + Sqr(
      Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21
      ) + Sqr(Yu22)) + 150*MassBp*Sqr(gN)*Sqr(QQp)*(Sqr(Yu00) + Sqr(Yu01) + Sqr(
      Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22
      )) + 150*MassBp*Sqr(gN)*Sqr(Qup)*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(
      Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)) + 450*(
      Yu00*(TYu10*(Yu00*Yu10 + Yu01*Yu11 + Yu02*Yu12) + TYu20*(Yu00*Yu20 + Yu01*
      Yu21 + Yu02*Yu22) + TYu00*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02))) + Yu01*(TYu11
      *(Yu00*Yu10 + Yu01*Yu11 + Yu02*Yu12) + TYu21*(Yu00*Yu20 + Yu01*Yu21 + Yu02*
      Yu22) + TYu01*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02))) + Yu02*(TYu12*(Yu00*Yu10
      + Yu01*Yu11 + Yu02*Yu12) + TYu22*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + TYu02
      *(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02))) + Yu10*(TYu00*(Yu00*Yu10 + Yu01*Yu11 +
      Yu02*Yu12) + TYu20*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + TYu10*(Sqr(Yu10) +
      Sqr(Yu11) + Sqr(Yu12))) + Yu11*(TYu01*(Yu00*Yu10 + Yu01*Yu11 + Yu02*Yu12) +
      TYu21*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + TYu11*(Sqr(Yu10) + Sqr(Yu11) +
      Sqr(Yu12))) + Yu12*(TYu02*(Yu00*Yu10 + Yu01*Yu11 + Yu02*Yu12) + TYu22*(Yu10*
      Yu20 + Yu11*Yu21 + Yu12*Yu22) + TYu12*(Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12))) +
      Yu20*(TYu00*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + TYu10*(Yu10*Yu20 + Yu11*
      Yu21 + Yu12*Yu22) + TYu20*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + Yu21*(TYu01
      *(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + TYu11*(Yu10*Yu20 + Yu11*Yu21 + Yu12*
      Yu22) + TYu21*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + Yu22*(TYu02*(Yu00*Yu20
      + Yu01*Yu21 + Yu02*Yu22) + TYu12*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + TYu22
      *(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))))) + ALambdax*(5.94*Power(g1,4) + 16.5*
      Power(g2,4) + 16*Power(gN,4)*Power(QH1p,4) + 16*Power(gN,4)*Power(QH2p,4) +
      10*Power(gN,4)*Power(QSp,4) - 6*(Yd00*(Yu00*(Yd00*Yu00 + Yd01*Yu01 + Yd02*
      Yu02) + Yu10*(Yd00*Yu10 + Yd01*Yu11 + Yd02*Yu12) + Yu20*(Yd00*Yu20 + Yd01*
      Yu21 + Yd02*Yu22)) + Yd01*(Yu01*(Yd00*Yu00 + Yd01*Yu01 + Yd02*Yu02) + Yu11*(
      Yd00*Yu10 + Yd01*Yu11 + Yd02*Yu12) + Yu21*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22
      )) + Yd02*(Yu02*(Yd00*Yu00 + Yd01*Yu01 + Yd02*Yu02) + Yu12*(Yd00*Yu10 + Yd01
      *Yu11 + Yd02*Yu12) + Yu22*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22)) + Yd10*(Yu00*
      (Yd10*Yu00 + Yd11*Yu01 + Yd12*Yu02) + Yu10*(Yd10*Yu10 + Yd11*Yu11 + Yd12*
      Yu12) + Yu20*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22)) + Yd11*(Yu01*(Yd10*Yu00 +
      Yd11*Yu01 + Yd12*Yu02) + Yu11*(Yd10*Yu10 + Yd11*Yu11 + Yd12*Yu12) + Yu21*(
      Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22)) + Yd12*(Yu02*(Yd10*Yu00 + Yd11*Yu01 +
      Yd12*Yu02) + Yu12*(Yd10*Yu10 + Yd11*Yu11 + Yd12*Yu12) + Yu22*(Yd10*Yu20 +
      Yd11*Yu21 + Yd12*Yu22)) + Yd20*(Yu00*(Yd20*Yu00 + Yd21*Yu01 + Yd22*Yu02) +
      Yu10*(Yd20*Yu10 + Yd21*Yu11 + Yd22*Yu12) + Yu20*(Yd20*Yu20 + Yd21*Yu21 +
      Yd22*Yu22)) + Yd21*(Yu01*(Yd20*Yu00 + Yd21*Yu01 + Yd22*Yu02) + Yu11*(Yd20*
      Yu10 + Yd21*Yu11 + Yd22*Yu12) + Yu21*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22)) +
      Yd22*(Yu02*(Yd20*Yu00 + Yd21*Yu01 + Yd22*Yu02) + Yu12*(Yd20*Yu10 + Yd21*Yu11
      + Yd22*Yu12) + Yu22*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22))) + 1.8*Sqr(g1)*Sqr
      (g2) + 0.8*Sqr(g1)*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10
      ) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)
      ) + 16*Sqr(g3)*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) +
      Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) -
      6*(Kappa00*(Kappa10*(Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12) +
      Kappa20*(Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22) + Kappa00*(Sqr
      (Kappa00) + Sqr(Kappa01) + Sqr(Kappa02))) + Kappa01*(Kappa11*(Kappa00*
      Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12) + Kappa21*(Kappa00*Kappa20 +
      Kappa01*Kappa21 + Kappa02*Kappa22) + Kappa01*(Sqr(Kappa00) + Sqr(Kappa01) +
      Sqr(Kappa02))) + Kappa02*(Kappa12*(Kappa00*Kappa10 + Kappa01*Kappa11 +
      Kappa02*Kappa12) + Kappa22*(Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*
      Kappa22) + Kappa02*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02))) + Kappa10*(
      Kappa00*(Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12) + Kappa20*(
      Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22) + Kappa10*(Sqr(Kappa10)
      + Sqr(Kappa11) + Sqr(Kappa12))) + Kappa11*(Kappa01*(Kappa00*Kappa10 +
      Kappa01*Kappa11 + Kappa02*Kappa12) + Kappa21*(Kappa10*Kappa20 + Kappa11*
      Kappa21 + Kappa12*Kappa22) + Kappa11*(Sqr(Kappa10) + Sqr(Kappa11) + Sqr(
      Kappa12))) + Kappa12*(Kappa02*(Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*
      Kappa12) + Kappa22*(Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22) +
      Kappa12*(Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12))) + Kappa20*(Kappa00*(
      Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22) + Kappa10*(Kappa10*
      Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22) + Kappa20*(Sqr(Kappa20) + Sqr(
      Kappa21) + Sqr(Kappa22))) + Kappa21*(Kappa01*(Kappa00*Kappa20 + Kappa01*
      Kappa21 + Kappa02*Kappa22) + Kappa11*(Kappa10*Kappa20 + Kappa11*Kappa21 +
      Kappa12*Kappa22) + Kappa21*(Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22))) +
      Kappa22*(Kappa02*(Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22) +
      Kappa12*(Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22) + Kappa22*(Sqr
      (Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)))) + 1.2*Sqr(g1)*(Sqr(Lambda1200) +
      Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 6*Sqr(g2)*(Sqr(
      Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) - 4*(
      Lambda1200*(Lambda1210*(Lambda1200*Lambda1210 + Lambda1201*Lambda1211) +
      Lambda1200*(Sqr(Lambda1200) + Sqr(Lambda1201))) + Lambda1201*(Lambda1211*(
      Lambda1200*Lambda1210 + Lambda1201*Lambda1211) + Lambda1201*(Sqr(Lambda1200)
      + Sqr(Lambda1201))) + Lambda1210*(Lambda1200*(Lambda1200*Lambda1210 +
      Lambda1201*Lambda1211) + Lambda1210*(Sqr(Lambda1210) + Sqr(Lambda1211))) +
      Lambda1211*(Lambda1201*(Lambda1200*Lambda1210 + Lambda1201*Lambda1211) +
      Lambda1211*(Sqr(Lambda1210) + Sqr(Lambda1211)))) + 6*Sqr(gN)*(Sqr(Kappa00) +
      Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) +
      Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22))*Sqr(QDxbarp) + 6*Sqr(gN)*(Sqr(
      Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(
      Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22))*Sqr(QDxp) + 1.2*Sqr(
      g1)*Sqr(gN)*Sqr(QH1p) + 6*Sqr(g2)*Sqr(gN)*Sqr(QH1p) + 4*Sqr(gN)*(Sqr(
      Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211))*Sqr(QH1p)
      + 18*Power(gN,4)*Sqr(Qdp)*Sqr(QH1p) + 18*Power(gN,4)*Sqr(QDxbarp)*Sqr(QH1p)
      + 18*Power(gN,4)*Sqr(QDxp)*Sqr(QH1p) + 6*Power(gN,4)*Sqr(Qep)*Sqr(QH1p) +
      1.2*Sqr(g1)*Sqr(gN)*Sqr(QH2p) + 6*Sqr(g2)*Sqr(gN)*Sqr(QH2p) + 4*Sqr(gN)*(Sqr
      (Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211))*Sqr(QH2p
      ) + 18*Power(gN,4)*Sqr(Qdp)*Sqr(QH2p) + 18*Power(gN,4)*Sqr(QDxbarp)*Sqr(QH2p
      ) + 18*Power(gN,4)*Sqr(QDxp)*Sqr(QH2p) + 6*Power(gN,4)*Sqr(Qep)*Sqr(QH2p) +
      24*Power(gN,4)*Sqr(QH1p)*Sqr(QH2p) + 4*Power(gN,4)*Sqr(QH1p)*Sqr(QHpbarp) +
      4*Power(gN,4)*Sqr(QH2p)*Sqr(QHpbarp) + 4*Power(gN,4)*Sqr(QH1p)*Sqr(QHpp) + 4
      *Power(gN,4)*Sqr(QH2p)*Sqr(QHpp) + 12*Power(gN,4)*Sqr(QH1p)*Sqr(QLp) + 12*
      Power(gN,4)*Sqr(QH2p)*Sqr(QLp) + 36*Power(gN,4)*Sqr(QH1p)*Sqr(QQp) + 36*
      Power(gN,4)*Sqr(QH2p)*Sqr(QQp) - 6*Sqr(gN)*(Sqr(Kappa00) + Sqr(Kappa01) +
      Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) +
      Sqr(Kappa21) + Sqr(Kappa22))*Sqr(QSp) - 4*Sqr(gN)*(Sqr(Lambda1200) + Sqr(
      Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211))*Sqr(QSp) + 18*Power(gN,4)*
      Sqr(Qdp)*Sqr(QSp) + 18*Power(gN,4)*Sqr(QDxbarp)*Sqr(QSp) + 18*Power(gN,4)*
      Sqr(QDxp)*Sqr(QSp) + 6*Power(gN,4)*Sqr(Qep)*Sqr(QSp) + 18*Power(gN,4)*Sqr(
      QH1p)*Sqr(QSp) + 18*Power(gN,4)*Sqr(QH2p)*Sqr(QSp) + 4*Power(gN,4)*Sqr(
      QHpbarp)*Sqr(QSp) + 4*Power(gN,4)*Sqr(QHpp)*Sqr(QSp) + 12*Power(gN,4)*Sqr(
      QLp)*Sqr(QSp) + 36*Power(gN,4)*Sqr(QQp)*Sqr(QSp) + 18*Power(gN,4)*Sqr(QH1p)*
      Sqr(Qup) + 18*Power(gN,4)*Sqr(QH2p)*Sqr(Qup) + 18*Power(gN,4)*Sqr(QSp)*Sqr(
      Qup) - 0.4*(Sqr(g1) - 5*(8*Sqr(g3) + 3*Sqr(gN)*(Sqr(Qdp) - Sqr(QH1p) + Sqr(
      QQp))))*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(
      Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) - 9*(Yd00*(Yd10*(Yd00*Yd10 + Yd01
      *Yd11 + Yd02*Yd12) + Yd20*(Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22) + Yd00*(Sqr(
      Yd00) + Sqr(Yd01) + Sqr(Yd02))) + Yd01*(Yd11*(Yd00*Yd10 + Yd01*Yd11 + Yd02*
      Yd12) + Yd21*(Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22) + Yd01*(Sqr(Yd00) + Sqr(
      Yd01) + Sqr(Yd02))) + Yd02*(Yd12*(Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + Yd22*
      (Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22) + Yd02*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02
      ))) + Yd10*(Yd00*(Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + Yd20*(Yd10*Yd20 +
      Yd11*Yd21 + Yd12*Yd22) + Yd10*(Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12))) + Yd11*(
      Yd01*(Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + Yd21*(Yd10*Yd20 + Yd11*Yd21 +
      Yd12*Yd22) + Yd11*(Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12))) + Yd12*(Yd02*(Yd00*
      Yd10 + Yd01*Yd11 + Yd02*Yd12) + Yd22*(Yd10*Yd20 + Yd11*Yd21 + Yd12*Yd22) +
      Yd12*(Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12))) + Yd20*(Yd00*(Yd00*Yd20 + Yd01*
      Yd21 + Yd02*Yd22) + Yd10*(Yd10*Yd20 + Yd11*Yd21 + Yd12*Yd22) + Yd20*(Sqr(
      Yd20) + Sqr(Yd21) + Sqr(Yd22))) + Yd21*(Yd01*(Yd00*Yd20 + Yd01*Yd21 + Yd02*
      Yd22) + Yd11*(Yd10*Yd20 + Yd11*Yd21 + Yd12*Yd22) + Yd21*(Sqr(Yd20) + Sqr(
      Yd21) + Sqr(Yd22))) + Yd22*(Yd02*(Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22) + Yd12*
      (Yd10*Yd20 + Yd11*Yd21 + Yd12*Yd22) + Yd22*(Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22
      )))) + 0.4*(3*Sqr(g1) + 5*Sqr(gN)*(Sqr(Qep) - Sqr(QH1p) + Sqr(QLp)))*(Sqr(
      Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20
      ) + Sqr(Ye21) + Sqr(Ye22)) - 3*(Ye00*(Ye10*(Ye00*Ye10 + Ye01*Ye11 + Ye02*
      Ye12) + Ye20*(Ye00*Ye20 + Ye01*Ye21 + Ye02*Ye22) + Ye00*(Sqr(Ye00) + Sqr(
      Ye01) + Sqr(Ye02))) + Ye01*(Ye11*(Ye00*Ye10 + Ye01*Ye11 + Ye02*Ye12) + Ye21*
      (Ye00*Ye20 + Ye01*Ye21 + Ye02*Ye22) + Ye01*(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02
      ))) + Ye02*(Ye12*(Ye00*Ye10 + Ye01*Ye11 + Ye02*Ye12) + Ye22*(Ye00*Ye20 +
      Ye01*Ye21 + Ye02*Ye22) + Ye02*(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02))) + Ye10*(
      Ye00*(Ye00*Ye10 + Ye01*Ye11 + Ye02*Ye12) + Ye20*(Ye10*Ye20 + Ye11*Ye21 +
      Ye12*Ye22) + Ye10*(Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12))) + Ye11*(Ye01*(Ye00*
      Ye10 + Ye01*Ye11 + Ye02*Ye12) + Ye21*(Ye10*Ye20 + Ye11*Ye21 + Ye12*Ye22) +
      Ye11*(Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12))) + Ye12*(Ye02*(Ye00*Ye10 + Ye01*
      Ye11 + Ye02*Ye12) + Ye22*(Ye10*Ye20 + Ye11*Ye21 + Ye12*Ye22) + Ye12*(Sqr(
      Ye10) + Sqr(Ye11) + Sqr(Ye12))) + Ye20*(Ye00*(Ye00*Ye20 + Ye01*Ye21 + Ye02*
      Ye22) + Ye10*(Ye10*Ye20 + Ye11*Ye21 + Ye12*Ye22) + Ye20*(Sqr(Ye20) + Sqr(
      Ye21) + Sqr(Ye22))) + Ye21*(Ye01*(Ye00*Ye20 + Ye01*Ye21 + Ye02*Ye22) + Ye11*
      (Ye10*Ye20 + Ye11*Ye21 + Ye12*Ye22) + Ye21*(Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22
      ))) + Ye22*(Ye02*(Ye00*Ye20 + Ye01*Ye21 + Ye02*Ye22) + Ye12*(Ye10*Ye20 +
      Ye11*Ye21 + Ye12*Ye22) + Ye22*(Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22)))) + 0.8*
      Sqr(g1)*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(
      Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)) + 16*Sqr(g3)*(Sqr(Yu00) + Sqr(
      Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21
      ) + Sqr(Yu22)) - 6*Sqr(gN)*Sqr(QH2p)*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) +
      Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)) + 6*
      Sqr(gN)*Sqr(QQp)*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11)
      + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)) + 6*Sqr(gN)*Sqr(Qup)*(Sqr(
      Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20
      ) + Sqr(Yu21) + Sqr(Yu22)) - 9*(Yu00*(Yu10*(Yu00*Yu10 + Yu01*Yu11 + Yu02*
      Yu12) + Yu20*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu00*(Sqr(Yu00) + Sqr(
      Yu01) + Sqr(Yu02))) + Yu01*(Yu11*(Yu00*Yu10 + Yu01*Yu11 + Yu02*Yu12) + Yu21*
      (Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu01*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02
      ))) + Yu02*(Yu12*(Yu00*Yu10 + Yu01*Yu11 + Yu02*Yu12) + Yu22*(Yu00*Yu20 +
      Yu01*Yu21 + Yu02*Yu22) + Yu02*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02))) + Yu10*(
      Yu00*(Yu00*Yu10 + Yu01*Yu11 + Yu02*Yu12) + Yu20*(Yu10*Yu20 + Yu11*Yu21 +
      Yu12*Yu22) + Yu10*(Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12))) + Yu11*(Yu01*(Yu00*
      Yu10 + Yu01*Yu11 + Yu02*Yu12) + Yu21*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) +
      Yu11*(Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12))) + Yu12*(Yu02*(Yu00*Yu10 + Yu01*
      Yu11 + Yu02*Yu12) + Yu22*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + Yu12*(Sqr(
      Yu10) + Sqr(Yu11) + Sqr(Yu12))) + Yu20*(Yu00*(Yu00*Yu20 + Yu01*Yu21 + Yu02*
      Yu22) + Yu10*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + Yu20*(Sqr(Yu20) + Sqr(
      Yu21) + Sqr(Yu22))) + Yu21*(Yu01*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu11*
      (Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + Yu21*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22
      ))) + Yu22*(Yu02*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu12*(Yu10*Yu20 +
      Yu11*Yu21 + Yu12*Yu22) + Yu22*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_two_loop_TLambdax_dALambdax() const 
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

   double deriv = -50*Power(Lambdax,5) + 0.6*Power(Lambdax,3)*(6*Sqr(g1) + 30*
      Sqr(g2) - 30*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) +
      Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) -
      20*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) +
      20*Sqr(gN)*Sqr(QH1p) + 20*Sqr(gN)*Sqr(QH2p) - 45*(Sqr(Yd00) + Sqr(Yd01) +
      Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(
      Yd22)) - 15*(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr
      (Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22)) - 45*(Sqr(Yu00) + Sqr(Yu01) +
      Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(
      Yu22))) + Lambdax*(5.94*Power(g1,4) + 16.5*Power(g2,4) + 16*Power(gN,4)*
      Power(QH1p,4) + 16*Power(gN,4)*Power(QH2p,4) + 10*Power(gN,4)*Power(QSp,4) -
      6*(Yd00*(Yu00*(Yd00*Yu00 + Yd01*Yu01 + Yd02*Yu02) + Yu10*(Yd00*Yu10 + Yd01*
      Yu11 + Yd02*Yu12) + Yu20*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22)) + Yd01*(Yu01*(
      Yd00*Yu00 + Yd01*Yu01 + Yd02*Yu02) + Yu11*(Yd00*Yu10 + Yd01*Yu11 + Yd02*Yu12
      ) + Yu21*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22)) + Yd02*(Yu02*(Yd00*Yu00 + Yd01
      *Yu01 + Yd02*Yu02) + Yu12*(Yd00*Yu10 + Yd01*Yu11 + Yd02*Yu12) + Yu22*(Yd00*
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
      Yd21*Yu21 + Yd22*Yu22))) + 1.8*Sqr(g1)*Sqr(g2) + 0.8*Sqr(g1)*(Sqr(Kappa00) +
      Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) +
      Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 16*Sqr(g3)*(Sqr(Kappa00) + Sqr
      (Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(
      Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) - 6*(Kappa00*(Kappa10*(Kappa00*
      Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12) + Kappa20*(Kappa00*Kappa20 +
      Kappa01*Kappa21 + Kappa02*Kappa22) + Kappa00*(Sqr(Kappa00) + Sqr(Kappa01) +
      Sqr(Kappa02))) + Kappa01*(Kappa11*(Kappa00*Kappa10 + Kappa01*Kappa11 +
      Kappa02*Kappa12) + Kappa21*(Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*
      Kappa22) + Kappa01*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02))) + Kappa02*(
      Kappa12*(Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12) + Kappa22*(
      Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22) + Kappa02*(Sqr(Kappa00)
      + Sqr(Kappa01) + Sqr(Kappa02))) + Kappa10*(Kappa00*(Kappa00*Kappa10 +
      Kappa01*Kappa11 + Kappa02*Kappa12) + Kappa20*(Kappa10*Kappa20 + Kappa11*
      Kappa21 + Kappa12*Kappa22) + Kappa10*(Sqr(Kappa10) + Sqr(Kappa11) + Sqr(
      Kappa12))) + Kappa11*(Kappa01*(Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*
      Kappa12) + Kappa21*(Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22) +
      Kappa11*(Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12))) + Kappa12*(Kappa02*(
      Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12) + Kappa22*(Kappa10*
      Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22) + Kappa12*(Sqr(Kappa10) + Sqr(
      Kappa11) + Sqr(Kappa12))) + Kappa20*(Kappa00*(Kappa00*Kappa20 + Kappa01*
      Kappa21 + Kappa02*Kappa22) + Kappa10*(Kappa10*Kappa20 + Kappa11*Kappa21 +
      Kappa12*Kappa22) + Kappa20*(Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22))) +
      Kappa21*(Kappa01*(Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22) +
      Kappa11*(Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22) + Kappa21*(Sqr
      (Kappa20) + Sqr(Kappa21) + Sqr(Kappa22))) + Kappa22*(Kappa02*(Kappa00*
      Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22) + Kappa12*(Kappa10*Kappa20 +
      Kappa11*Kappa21 + Kappa12*Kappa22) + Kappa22*(Sqr(Kappa20) + Sqr(Kappa21) +
      Sqr(Kappa22)))) + 1.2*Sqr(g1)*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(
      Lambda1210) + Sqr(Lambda1211)) + 6*Sqr(g2)*(Sqr(Lambda1200) + Sqr(Lambda1201
      ) + Sqr(Lambda1210) + Sqr(Lambda1211)) - 4*(Lambda1200*(Lambda1210*(
      Lambda1200*Lambda1210 + Lambda1201*Lambda1211) + Lambda1200*(Sqr(Lambda1200)
      + Sqr(Lambda1201))) + Lambda1201*(Lambda1211*(Lambda1200*Lambda1210 +
      Lambda1201*Lambda1211) + Lambda1201*(Sqr(Lambda1200) + Sqr(Lambda1201))) +
      Lambda1210*(Lambda1200*(Lambda1200*Lambda1210 + Lambda1201*Lambda1211) +
      Lambda1210*(Sqr(Lambda1210) + Sqr(Lambda1211))) + Lambda1211*(Lambda1201*(
      Lambda1200*Lambda1210 + Lambda1201*Lambda1211) + Lambda1211*(Sqr(Lambda1210)
      + Sqr(Lambda1211)))) + 6*Sqr(gN)*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02
      ) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21)
      + Sqr(Kappa22))*Sqr(QDxbarp) + 6*Sqr(gN)*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr
      (Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(
      Kappa21) + Sqr(Kappa22))*Sqr(QDxp) + 1.2*Sqr(g1)*Sqr(gN)*Sqr(QH1p) + 6*Sqr(
      g2)*Sqr(gN)*Sqr(QH1p) + 4*Sqr(gN)*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(
      Lambda1210) + Sqr(Lambda1211))*Sqr(QH1p) + 18*Power(gN,4)*Sqr(Qdp)*Sqr(QH1p)
      + 18*Power(gN,4)*Sqr(QDxbarp)*Sqr(QH1p) + 18*Power(gN,4)*Sqr(QDxp)*Sqr(QH1p
      ) + 6*Power(gN,4)*Sqr(Qep)*Sqr(QH1p) + 1.2*Sqr(g1)*Sqr(gN)*Sqr(QH2p) + 6*Sqr
      (g2)*Sqr(gN)*Sqr(QH2p) + 4*Sqr(gN)*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(
      Lambda1210) + Sqr(Lambda1211))*Sqr(QH2p) + 18*Power(gN,4)*Sqr(Qdp)*Sqr(QH2p)
      + 18*Power(gN,4)*Sqr(QDxbarp)*Sqr(QH2p) + 18*Power(gN,4)*Sqr(QDxp)*Sqr(QH2p
      ) + 6*Power(gN,4)*Sqr(Qep)*Sqr(QH2p) + 24*Power(gN,4)*Sqr(QH1p)*Sqr(QH2p) +
      4*Power(gN,4)*Sqr(QH1p)*Sqr(QHpbarp) + 4*Power(gN,4)*Sqr(QH2p)*Sqr(QHpbarp)
      + 4*Power(gN,4)*Sqr(QH1p)*Sqr(QHpp) + 4*Power(gN,4)*Sqr(QH2p)*Sqr(QHpp) + 12
      *Power(gN,4)*Sqr(QH1p)*Sqr(QLp) + 12*Power(gN,4)*Sqr(QH2p)*Sqr(QLp) + 36*
      Power(gN,4)*Sqr(QH1p)*Sqr(QQp) + 36*Power(gN,4)*Sqr(QH2p)*Sqr(QQp) - 6*Sqr(
      gN)*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11
      ) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22))*Sqr(QSp) - 4*
      Sqr(gN)*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(
      Lambda1211))*Sqr(QSp) + 18*Power(gN,4)*Sqr(Qdp)*Sqr(QSp) + 18*Power(gN,4)*
      Sqr(QDxbarp)*Sqr(QSp) + 18*Power(gN,4)*Sqr(QDxp)*Sqr(QSp) + 6*Power(gN,4)*
      Sqr(Qep)*Sqr(QSp) + 18*Power(gN,4)*Sqr(QH1p)*Sqr(QSp) + 18*Power(gN,4)*Sqr(
      QH2p)*Sqr(QSp) + 4*Power(gN,4)*Sqr(QHpbarp)*Sqr(QSp) + 4*Power(gN,4)*Sqr(
      QHpp)*Sqr(QSp) + 12*Power(gN,4)*Sqr(QLp)*Sqr(QSp) + 36*Power(gN,4)*Sqr(QQp)*
      Sqr(QSp) + 18*Power(gN,4)*Sqr(QH1p)*Sqr(Qup) + 18*Power(gN,4)*Sqr(QH2p)*Sqr(
      Qup) + 18*Power(gN,4)*Sqr(QSp)*Sqr(Qup) - 0.4*(Sqr(g1) - 5*(8*Sqr(g3) + 3*
      Sqr(gN)*(Sqr(Qdp) - Sqr(QH1p) + Sqr(QQp))))*(Sqr(Yd00) + Sqr(Yd01) + Sqr(
      Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22
      )) - 9*(Yd00*(Yd10*(Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + Yd20*(Yd00*Yd20 +
      Yd01*Yd21 + Yd02*Yd22) + Yd00*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02))) + Yd01*(
      Yd11*(Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + Yd21*(Yd00*Yd20 + Yd01*Yd21 +
      Yd02*Yd22) + Yd01*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02))) + Yd02*(Yd12*(Yd00*
      Yd10 + Yd01*Yd11 + Yd02*Yd12) + Yd22*(Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22) +
      Yd02*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02))) + Yd10*(Yd00*(Yd00*Yd10 + Yd01*
      Yd11 + Yd02*Yd12) + Yd20*(Yd10*Yd20 + Yd11*Yd21 + Yd12*Yd22) + Yd10*(Sqr(
      Yd10) + Sqr(Yd11) + Sqr(Yd12))) + Yd11*(Yd01*(Yd00*Yd10 + Yd01*Yd11 + Yd02*
      Yd12) + Yd21*(Yd10*Yd20 + Yd11*Yd21 + Yd12*Yd22) + Yd11*(Sqr(Yd10) + Sqr(
      Yd11) + Sqr(Yd12))) + Yd12*(Yd02*(Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + Yd22*
      (Yd10*Yd20 + Yd11*Yd21 + Yd12*Yd22) + Yd12*(Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12
      ))) + Yd20*(Yd00*(Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22) + Yd10*(Yd10*Yd20 +
      Yd11*Yd21 + Yd12*Yd22) + Yd20*(Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22))) + Yd21*(
      Yd01*(Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22) + Yd11*(Yd10*Yd20 + Yd11*Yd21 +
      Yd12*Yd22) + Yd21*(Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22))) + Yd22*(Yd02*(Yd00*
      Yd20 + Yd01*Yd21 + Yd02*Yd22) + Yd12*(Yd10*Yd20 + Yd11*Yd21 + Yd12*Yd22) +
      Yd22*(Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)))) + 0.4*(3*Sqr(g1) + 5*Sqr(gN)*(Sqr
      (Qep) - Sqr(QH1p) + Sqr(QLp)))*(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10
      ) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22)) - 3*(Ye00*(
      Ye10*(Ye00*Ye10 + Ye01*Ye11 + Ye02*Ye12) + Ye20*(Ye00*Ye20 + Ye01*Ye21 +
      Ye02*Ye22) + Ye00*(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02))) + Ye01*(Ye11*(Ye00*
      Ye10 + Ye01*Ye11 + Ye02*Ye12) + Ye21*(Ye00*Ye20 + Ye01*Ye21 + Ye02*Ye22) +
      Ye01*(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02))) + Ye02*(Ye12*(Ye00*Ye10 + Ye01*
      Ye11 + Ye02*Ye12) + Ye22*(Ye00*Ye20 + Ye01*Ye21 + Ye02*Ye22) + Ye02*(Sqr(
      Ye00) + Sqr(Ye01) + Sqr(Ye02))) + Ye10*(Ye00*(Ye00*Ye10 + Ye01*Ye11 + Ye02*
      Ye12) + Ye20*(Ye10*Ye20 + Ye11*Ye21 + Ye12*Ye22) + Ye10*(Sqr(Ye10) + Sqr(
      Ye11) + Sqr(Ye12))) + Ye11*(Ye01*(Ye00*Ye10 + Ye01*Ye11 + Ye02*Ye12) + Ye21*
      (Ye10*Ye20 + Ye11*Ye21 + Ye12*Ye22) + Ye11*(Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12
      ))) + Ye12*(Ye02*(Ye00*Ye10 + Ye01*Ye11 + Ye02*Ye12) + Ye22*(Ye10*Ye20 +
      Ye11*Ye21 + Ye12*Ye22) + Ye12*(Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12))) + Ye20*(
      Ye00*(Ye00*Ye20 + Ye01*Ye21 + Ye02*Ye22) + Ye10*(Ye10*Ye20 + Ye11*Ye21 +
      Ye12*Ye22) + Ye20*(Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22))) + Ye21*(Ye01*(Ye00*
      Ye20 + Ye01*Ye21 + Ye02*Ye22) + Ye11*(Ye10*Ye20 + Ye11*Ye21 + Ye12*Ye22) +
      Ye21*(Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22))) + Ye22*(Ye02*(Ye00*Ye20 + Ye01*
      Ye21 + Ye02*Ye22) + Ye12*(Ye10*Ye20 + Ye11*Ye21 + Ye12*Ye22) + Ye22*(Sqr(
      Ye20) + Sqr(Ye21) + Sqr(Ye22)))) + 0.8*Sqr(g1)*(Sqr(Yu00) + Sqr(Yu01) + Sqr(
      Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22
      )) + 16*Sqr(g3)*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) +
      Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)) - 6*Sqr(gN)*Sqr(QH2p)*(Sqr(
      Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20
      ) + Sqr(Yu21) + Sqr(Yu22)) + 6*Sqr(gN)*Sqr(QQp)*(Sqr(Yu00) + Sqr(Yu01) + Sqr
      (Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(
      Yu22)) + 6*Sqr(gN)*Sqr(Qup)*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) +
      Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)) - 9*(Yu00*(Yu10*
      (Yu00*Yu10 + Yu01*Yu11 + Yu02*Yu12) + Yu20*(Yu00*Yu20 + Yu01*Yu21 + Yu02*
      Yu22) + Yu00*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02))) + Yu01*(Yu11*(Yu00*Yu10 +
      Yu01*Yu11 + Yu02*Yu12) + Yu21*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu01*(
      Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02))) + Yu02*(Yu12*(Yu00*Yu10 + Yu01*Yu11 +
      Yu02*Yu12) + Yu22*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu02*(Sqr(Yu00) +
      Sqr(Yu01) + Sqr(Yu02))) + Yu10*(Yu00*(Yu00*Yu10 + Yu01*Yu11 + Yu02*Yu12) +
      Yu20*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + Yu10*(Sqr(Yu10) + Sqr(Yu11) + Sqr
      (Yu12))) + Yu11*(Yu01*(Yu00*Yu10 + Yu01*Yu11 + Yu02*Yu12) + Yu21*(Yu10*Yu20
      + Yu11*Yu21 + Yu12*Yu22) + Yu11*(Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12))) + Yu12*
      (Yu02*(Yu00*Yu10 + Yu01*Yu11 + Yu02*Yu12) + Yu22*(Yu10*Yu20 + Yu11*Yu21 +
      Yu12*Yu22) + Yu12*(Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12))) + Yu20*(Yu00*(Yu00*
      Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu10*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) +
      Yu20*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + Yu21*(Yu01*(Yu00*Yu20 + Yu01*
      Yu21 + Yu02*Yu22) + Yu11*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + Yu21*(Sqr(
      Yu20) + Sqr(Yu21) + Sqr(Yu22))) + Yu22*(Yu02*(Yu00*Yu20 + Yu01*Yu21 + Yu02*
      Yu22) + Yu12*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + Yu22*(Sqr(Yu20) + Sqr(
      Yu21) + Sqr(Yu22)))));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_two_loop_TLambdax_dAYu22() const 
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
   const double g3 = model.get_g3();
   const double gN = model.get_gN();

   double deriv = -18*Power(Lambdax,3)*Sqr(Yu22) - 0.08*Lambdax*(150*(Yd02*Yu22
      *(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22) + Yd12*Yu22*(Yd10*Yu20 + Yd11*Yu21 +
      Yd12*Yu22) + Yd22*Yu22*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22)) - 20*Sqr(g1)*Sqr
      (Yu22) - 400*Sqr(g3)*Sqr(Yu22) + 150*Sqr(gN)*Sqr(QH2p)*Sqr(Yu22) - 150*Sqr(
      gN)*Sqr(QQp)*Sqr(Yu22) - 150*Sqr(gN)*Sqr(Qup)*Sqr(Yu22) + 450*(Yu02*Yu22*(
      Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu12*Yu22*(Yu10*Yu20 + Yu11*Yu21 + Yu12
      *Yu22) + Sqr(Yu22)*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_two_loop_TLambdax_dMassB() const 
{
   const lowE6SSM_input_parameters inputs = model.get_input();

   const auto QH1p = inputs.QH1p;
   const auto QH2p = inputs.QH2p;
   
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
   const double gN = model.get_gN();

   double deriv = -2.4*Power(Lambdax,3)*Sqr(g1) - 0.08*Lambdax*(297*Power(g1,4)
      + 45*Sqr(g1)*Sqr(g2) + 20*Sqr(g1)*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(
      Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(
      Kappa21) + Sqr(Kappa22)) + 30*Sqr(g1)*(Sqr(Lambda1200) + Sqr(Lambda1201) +
      Sqr(Lambda1210) + Sqr(Lambda1211)) + 30*Sqr(g1)*Sqr(gN)*Sqr(QH1p) + 30*Sqr(
      g1)*Sqr(gN)*Sqr(QH2p) - 10*Sqr(g1)*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(
      Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + 30*Sqr(
      g1)*(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) +
      Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22)) + 20*Sqr(g1)*(Sqr(Yu00) + Sqr(Yu01) +
      Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(
      Yu22)));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_two_loop_TLambdax_dMassWB() const 
{
   const lowE6SSM_input_parameters inputs = model.get_input();

   const auto QH1p = inputs.QH1p;
   const auto QH2p = inputs.QH2p;

   const double Lambda1200 = model.get_Lambda12(0,0);
   const double Lambda1201 = model.get_Lambda12(0,1);
   const double Lambda1210 = model.get_Lambda12(1,0);
   const double Lambda1211 = model.get_Lambda12(1,1);
   const double Lambdax = model.get_Lambdax();
   const double g1 = model.get_g1();
   const double g2 = model.get_g2();
   const double gN = model.get_gN();

   double deriv = -12*Power(Lambdax,3)*Sqr(g2) - 0.08*Lambdax*(825*Power(g2,4)
      + 45*Sqr(g1)*Sqr(g2) + 150*Sqr(g2)*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(
      Lambda1210) + Sqr(Lambda1211)) + 150*Sqr(g2)*Sqr(gN)*Sqr(QH1p) + 150*Sqr(g2)
      *Sqr(gN)*Sqr(QH2p));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_two_loop_TLambdax_dMassG() const 
{
   
   const double Yd00 = model.get_Yd(0,0);
   const double Yd01 = model.get_Yd(0,1);
   const double Yd02 = model.get_Yd(0,2);
   const double Yd10 = model.get_Yd(1,0);
   const double Yd11 = model.get_Yd(1,1);
   const double Yd12 = model.get_Yd(1,2);
   const double Yd20 = model.get_Yd(2,0);
   const double Yd21 = model.get_Yd(2,1);
   const double Yd22 = model.get_Yd(2,2);
   const double Kappa00 = model.get_Kappa(0,0);
   const double Kappa01 = model.get_Kappa(0,1);
   const double Kappa02 = model.get_Kappa(0,2);
   const double Kappa10 = model.get_Kappa(1,0);
   const double Kappa11 = model.get_Kappa(1,1);
   const double Kappa12 = model.get_Kappa(1,2);
   const double Kappa20 = model.get_Kappa(2,0);
   const double Kappa21 = model.get_Kappa(2,1);
   const double Kappa22 = model.get_Kappa(2,2);
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
   const double g3 = model.get_g3();

   double deriv = -0.08*Lambdax*(400*Sqr(g3)*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr
      (Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(
      Kappa21) + Sqr(Kappa22)) + 400*Sqr(g3)*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) +
      Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + 400
      *Sqr(g3)*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(
      Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_two_loop_TLambdax_dMassBp() const 
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
   const double gN = model.get_gN();

   double deriv = -0.4*Power(Lambdax,3)*(20*Sqr(gN)*Sqr(QH1p) + 20*Sqr(gN)*Sqr(
      QH2p)) - 0.08*Lambdax*(800*Power(gN,4)*Power(QH1p,4) + 800*Power(gN,4)*Power
      (QH2p,4) + 500*Power(gN,4)*Power(QSp,4) + 150*Sqr(gN)*(Sqr(Kappa00) + Sqr(
      Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(
      Kappa20) + Sqr(Kappa21) + Sqr(Kappa22))*Sqr(QDxbarp) + 150*Sqr(gN)*(Sqr(
      Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(
      Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22))*Sqr(QDxp) + 30*Sqr(g1
      )*Sqr(gN)*Sqr(QH1p) + 150*Sqr(g2)*Sqr(gN)*Sqr(QH1p) + 100*Sqr(gN)*(Sqr(
      Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211))*Sqr(QH1p)
      + 900*Power(gN,4)*Sqr(Qdp)*Sqr(QH1p) + 900*Power(gN,4)*Sqr(QDxbarp)*Sqr(
      QH1p) + 900*Power(gN,4)*Sqr(QDxp)*Sqr(QH1p) + 300*Power(gN,4)*Sqr(Qep)*Sqr(
      QH1p) + 30*Sqr(g1)*Sqr(gN)*Sqr(QH2p) + 150*Sqr(g2)*Sqr(gN)*Sqr(QH2p) + 100*
      Sqr(gN)*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(
      Lambda1211))*Sqr(QH2p) + 900*Power(gN,4)*Sqr(Qdp)*Sqr(QH2p) + 900*Power(gN,4
      )*Sqr(QDxbarp)*Sqr(QH2p) + 900*Power(gN,4)*Sqr(QDxp)*Sqr(QH2p) + 300*Power(
      gN,4)*Sqr(Qep)*Sqr(QH2p) + 1200*Power(gN,4)*Sqr(QH1p)*Sqr(QH2p) + 200*Power(
      gN,4)*Sqr(QH1p)*Sqr(QHpbarp) + 200*Power(gN,4)*Sqr(QH2p)*Sqr(QHpbarp) + 200*
      Power(gN,4)*Sqr(QH1p)*Sqr(QHpp) + 200*Power(gN,4)*Sqr(QH2p)*Sqr(QHpp) + 600*
      Power(gN,4)*Sqr(QH1p)*Sqr(QLp) + 600*Power(gN,4)*Sqr(QH2p)*Sqr(QLp) + 1800*
      Power(gN,4)*Sqr(QH1p)*Sqr(QQp) + 1800*Power(gN,4)*Sqr(QH2p)*Sqr(QQp) - 150*
      Sqr(gN)*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(
      Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22))*Sqr(
      QSp) - 100*Sqr(gN)*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) +
      Sqr(Lambda1211))*Sqr(QSp) + 900*Power(gN,4)*Sqr(Qdp)*Sqr(QSp) + 900*Power(gN
      ,4)*Sqr(QDxbarp)*Sqr(QSp) + 900*Power(gN,4)*Sqr(QDxp)*Sqr(QSp) + 300*Power(
      gN,4)*Sqr(Qep)*Sqr(QSp) + 900*Power(gN,4)*Sqr(QH1p)*Sqr(QSp) + 900*Power(gN,
      4)*Sqr(QH2p)*Sqr(QSp) + 200*Power(gN,4)*Sqr(QHpbarp)*Sqr(QSp) + 200*Power(gN
      ,4)*Sqr(QHpp)*Sqr(QSp) + 600*Power(gN,4)*Sqr(QLp)*Sqr(QSp) + 1800*Power(gN,4
      )*Sqr(QQp)*Sqr(QSp) + 900*Power(gN,4)*Sqr(QH1p)*Sqr(Qup) + 900*Power(gN,4)*
      Sqr(QH2p)*Sqr(Qup) + 900*Power(gN,4)*Sqr(QSp)*Sqr(Qup) + 150*Sqr(gN)*(Sqr(
      Qdp) - Sqr(QH1p) + Sqr(QQp))*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10)
      + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + 50*Sqr(gN)*(
      Sqr(Qep) - Sqr(QH1p) + Sqr(QLp))*(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(
      Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22)) - 150*Sqr
      (gN)*Sqr(QH2p)*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) +
      Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)) + 150*Sqr(gN)*Sqr(QQp)*(Sqr(
      Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20
      ) + Sqr(Yu21) + Sqr(Yu22)) + 150*Sqr(gN)*Sqr(Qup)*(Sqr(Yu00) + Sqr(Yu01) +
      Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(
      Yu22)));


   return twoLoop*deriv;
}

// leading log
double lowE6SSM_tuning_calculator::deriv_dlead_log_one_loop_TLambdax_dLambdax() const 
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

   double deriv = 23.04*Power(g1,4)*MassB + 9.6*Power(g1,3)*(-1.2*ALambdax*g1 +
      2.4*g1*MassB) + 48*Power(g2,4)*MassWB + 4*Power(g2,3)*(-6*ALambdax*g2 + 12*
      g2*MassWB) + 6*Kappa00*Lambdax*(8*ALambdax*Kappa00*Lambdax + 4*Lambdax*
      TKappa00) + 4*Kappa00*Lambdax*(6*ALambdax*Kappa00*Lambdax + 6*Lambdax*
      TKappa00) + 6*Kappa01*Lambdax*(8*ALambdax*Kappa01*Lambdax + 4*Lambdax*
      TKappa01) + 4*Kappa01*Lambdax*(6*ALambdax*Kappa01*Lambdax + 6*Lambdax*
      TKappa01) + 6*Kappa02*Lambdax*(8*ALambdax*Kappa02*Lambdax + 4*Lambdax*
      TKappa02) + 4*Kappa02*Lambdax*(6*ALambdax*Kappa02*Lambdax + 6*Lambdax*
      TKappa02) + 6*Kappa10*Lambdax*(8*ALambdax*Kappa10*Lambdax + 4*Lambdax*
      TKappa10) + 4*Kappa10*Lambdax*(6*ALambdax*Kappa10*Lambdax + 6*Lambdax*
      TKappa10) + 6*Kappa11*Lambdax*(8*ALambdax*Kappa11*Lambdax + 4*Lambdax*
      TKappa11) + 4*Kappa11*Lambdax*(6*ALambdax*Kappa11*Lambdax + 6*Lambdax*
      TKappa11) + 6*Kappa12*Lambdax*(8*ALambdax*Kappa12*Lambdax + 4*Lambdax*
      TKappa12) + 4*Kappa12*Lambdax*(6*ALambdax*Kappa12*Lambdax + 6*Lambdax*
      TKappa12) + 6*Kappa20*Lambdax*(8*ALambdax*Kappa20*Lambdax + 4*Lambdax*
      TKappa20) + 4*Kappa20*Lambdax*(6*ALambdax*Kappa20*Lambdax + 6*Lambdax*
      TKappa20) + 6*Kappa21*Lambdax*(8*ALambdax*Kappa21*Lambdax + 4*Lambdax*
      TKappa21) + 4*Kappa21*Lambdax*(6*ALambdax*Kappa21*Lambdax + 6*Lambdax*
      TKappa21) + 6*Kappa22*Lambdax*(8*ALambdax*Kappa22*Lambdax + 4*Lambdax*
      TKappa22) + 4*Kappa22*Lambdax*(6*ALambdax*Kappa22*Lambdax + 6*Lambdax*
      TKappa22) + 4*Lambda1200*Lambdax*(4*ALambdax*Lambda1200*Lambdax + 4*Lambdax*
      TLambda1200) + 4*Lambda1200*Lambdax*(8*ALambdax*Lambda1200*Lambdax + 4*
      Lambdax*TLambda1200) + 4*Lambda1201*Lambdax*(4*ALambdax*Lambda1201*Lambdax +
      4*Lambdax*TLambda1201) + 4*Lambda1201*Lambdax*(8*ALambdax*Lambda1201*
      Lambdax + 4*Lambdax*TLambda1201) + 4*Lambda1210*Lambdax*(4*ALambdax*
      Lambda1210*Lambdax + 4*Lambdax*TLambda1210) + 4*Lambda1210*Lambdax*(8*
      ALambdax*Lambda1210*Lambdax + 4*Lambdax*TLambda1210) + 4*Lambda1211*Lambdax*
      (4*ALambdax*Lambda1211*Lambdax + 4*Lambdax*TLambda1211) + 4*Lambda1211*
      Lambdax*(8*ALambdax*Lambda1211*Lambdax + 4*Lambdax*TLambda1211) + 6*Lambdax*
      Yd00*(2*Lambdax*TYd00 + 4*ALambdax*Lambdax*Yd00) + 2*Lambdax*Yd00*(6*Lambdax
      *TYd00 + 6*ALambdax*Lambdax*Yd00) + 6*Lambdax*Yd01*(2*Lambdax*TYd01 + 4*
      ALambdax*Lambdax*Yd01) + 2*Lambdax*Yd01*(6*Lambdax*TYd01 + 6*ALambdax*
      Lambdax*Yd01) + 6*Lambdax*Yd02*(2*Lambdax*TYd02 + 4*ALambdax*Lambdax*Yd02) +
      2*Lambdax*Yd02*(6*Lambdax*TYd02 + 6*ALambdax*Lambdax*Yd02) + 6*Lambdax*Yd10
      *(2*Lambdax*TYd10 + 4*ALambdax*Lambdax*Yd10) + 2*Lambdax*Yd10*(6*Lambdax*
      TYd10 + 6*ALambdax*Lambdax*Yd10) + 6*Lambdax*Yd11*(2*Lambdax*TYd11 + 4*
      ALambdax*Lambdax*Yd11) + 2*Lambdax*Yd11*(6*Lambdax*TYd11 + 6*ALambdax*
      Lambdax*Yd11) + 6*Lambdax*Yd12*(2*Lambdax*TYd12 + 4*ALambdax*Lambdax*Yd12) +
      2*Lambdax*Yd12*(6*Lambdax*TYd12 + 6*ALambdax*Lambdax*Yd12) + 6*Lambdax*Yd20
      *(2*Lambdax*TYd20 + 4*ALambdax*Lambdax*Yd20) + 2*Lambdax*Yd20*(6*Lambdax*
      TYd20 + 6*ALambdax*Lambdax*Yd20) + 6*Lambdax*Yd21*(2*Lambdax*TYd21 + 4*
      ALambdax*Lambdax*Yd21) + 2*Lambdax*Yd21*(6*Lambdax*TYd21 + 6*ALambdax*
      Lambdax*Yd21) + 6*Lambdax*Yd22*(2*Lambdax*TYd22 + 4*ALambdax*Lambdax*Yd22) +
      2*Lambdax*Yd22*(6*Lambdax*TYd22 + 6*ALambdax*Lambdax*Yd22) + 2*Lambdax*Ye00
      *(2*Lambdax*TYe00 + 2*ALambdax*Lambdax*Ye00) + 2*Lambdax*Ye00*(2*Lambdax*
      TYe00 + 4*ALambdax*Lambdax*Ye00) + 2*Lambdax*Ye01*(2*Lambdax*TYe01 + 2*
      ALambdax*Lambdax*Ye01) + 2*Lambdax*Ye01*(2*Lambdax*TYe01 + 4*ALambdax*
      Lambdax*Ye01) + 2*Lambdax*Ye02*(2*Lambdax*TYe02 + 2*ALambdax*Lambdax*Ye02) +
      2*Lambdax*Ye02*(2*Lambdax*TYe02 + 4*ALambdax*Lambdax*Ye02) + 2*Lambdax*Ye10
      *(2*Lambdax*TYe10 + 2*ALambdax*Lambdax*Ye10) + 2*Lambdax*Ye10*(2*Lambdax*
      TYe10 + 4*ALambdax*Lambdax*Ye10) + 2*Lambdax*Ye11*(2*Lambdax*TYe11 + 2*
      ALambdax*Lambdax*Ye11) + 2*Lambdax*Ye11*(2*Lambdax*TYe11 + 4*ALambdax*
      Lambdax*Ye11) + 2*Lambdax*Ye12*(2*Lambdax*TYe12 + 2*ALambdax*Lambdax*Ye12) +
      2*Lambdax*Ye12*(2*Lambdax*TYe12 + 4*ALambdax*Lambdax*Ye12) + 2*Lambdax*Ye20
      *(2*Lambdax*TYe20 + 2*ALambdax*Lambdax*Ye20) + 2*Lambdax*Ye20*(2*Lambdax*
      TYe20 + 4*ALambdax*Lambdax*Ye20) + 2*Lambdax*Ye21*(2*Lambdax*TYe21 + 2*
      ALambdax*Lambdax*Ye21) + 2*Lambdax*Ye21*(2*Lambdax*TYe21 + 4*ALambdax*
      Lambdax*Ye21) + 2*Lambdax*Ye22*(2*Lambdax*TYe22 + 2*ALambdax*Lambdax*Ye22) +
      2*Lambdax*Ye22*(2*Lambdax*TYe22 + 4*ALambdax*Lambdax*Ye22) + 6*Lambdax*Yu00
      *(2*Lambdax*TYu00 + 4*ALambdax*Lambdax*Yu00) + 2*Lambdax*Yu00*(6*Lambdax*
      TYu00 + 6*ALambdax*Lambdax*Yu00) + 6*Lambdax*Yu01*(2*Lambdax*TYu01 + 4*
      ALambdax*Lambdax*Yu01) + 2*Lambdax*Yu01*(6*Lambdax*TYu01 + 6*ALambdax*
      Lambdax*Yu01) + 6*Lambdax*Yu02*(2*Lambdax*TYu02 + 4*ALambdax*Lambdax*Yu02) +
      2*Lambdax*Yu02*(6*Lambdax*TYu02 + 6*ALambdax*Lambdax*Yu02) + 6*Lambdax*Yu10
      *(2*Lambdax*TYu10 + 4*ALambdax*Lambdax*Yu10) + 2*Lambdax*Yu10*(6*Lambdax*
      TYu10 + 6*ALambdax*Lambdax*Yu10) + 6*Lambdax*Yu11*(2*Lambdax*TYu11 + 4*
      ALambdax*Lambdax*Yu11) + 2*Lambdax*Yu11*(6*Lambdax*TYu11 + 6*ALambdax*
      Lambdax*Yu11) + 6*Lambdax*Yu12*(2*Lambdax*TYu12 + 4*ALambdax*Lambdax*Yu12) +
      2*Lambdax*Yu12*(6*Lambdax*TYu12 + 6*ALambdax*Lambdax*Yu12) + 6*Lambdax*Yu20
      *(2*Lambdax*TYu20 + 4*ALambdax*Lambdax*Yu20) + 2*Lambdax*Yu20*(6*Lambdax*
      TYu20 + 6*ALambdax*Lambdax*Yu20) + 6*Lambdax*Yu21*(2*Lambdax*TYu21 + 4*
      ALambdax*Lambdax*Yu21) + 2*Lambdax*Yu21*(6*Lambdax*TYu21 + 6*ALambdax*
      Lambdax*Yu21) + 6*Lambdax*Yu22*(2*Lambdax*TYu22 + 4*ALambdax*Lambdax*Yu22) +
      2*Lambdax*Yu22*(6*Lambdax*TYu22 + 6*ALambdax*Lambdax*Yu22) + (6*ALambdax*
      Kappa00 + 6*TKappa00)*(2*(Kappa10*(Kappa00*Kappa10 + Kappa01*Kappa11 +
      Kappa02*Kappa12) + Kappa20*(Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*
      Kappa22) + Kappa00*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02))) + Kappa00*(
      -0.26666666666666666*Sqr(g1) - 5.333333333333333*Sqr(g3) + 3*(Sqr(Kappa00) +
      Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) +
      Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*(Sqr(Lambda1200) + Sqr(
      Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 2*Sqr(Lambdax) - 2*Sqr(gN
      )*Sqr(QDxbarp) - 2*Sqr(gN)*Sqr(QDxp) - 2*Sqr(gN)*Sqr(QSp))) + (6*ALambdax*
      Kappa01 + 6*TKappa01)*(2*(Kappa11*(Kappa00*Kappa10 + Kappa01*Kappa11 +
      Kappa02*Kappa12) + Kappa21*(Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*
      Kappa22) + Kappa01*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02))) + Kappa01*(
      -0.26666666666666666*Sqr(g1) - 5.333333333333333*Sqr(g3) + 3*(Sqr(Kappa00) +
      Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) +
      Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*(Sqr(Lambda1200) + Sqr(
      Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 2*Sqr(Lambdax) - 2*Sqr(gN
      )*Sqr(QDxbarp) - 2*Sqr(gN)*Sqr(QDxp) - 2*Sqr(gN)*Sqr(QSp))) + (6*ALambdax*
      Kappa02 + 6*TKappa02)*(2*(Kappa12*(Kappa00*Kappa10 + Kappa01*Kappa11 +
      Kappa02*Kappa12) + Kappa22*(Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*
      Kappa22) + Kappa02*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02))) + Kappa02*(
      -0.26666666666666666*Sqr(g1) - 5.333333333333333*Sqr(g3) + 3*(Sqr(Kappa00) +
      Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) +
      Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*(Sqr(Lambda1200) + Sqr(
      Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 2*Sqr(Lambdax) - 2*Sqr(gN
      )*Sqr(QDxbarp) - 2*Sqr(gN)*Sqr(QDxp) - 2*Sqr(gN)*Sqr(QSp))) + (6*ALambdax*
      Kappa10 + 6*TKappa10)*(2*(Kappa00*(Kappa00*Kappa10 + Kappa01*Kappa11 +
      Kappa02*Kappa12) + Kappa20*(Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*
      Kappa22) + Kappa10*(Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12))) + Kappa10*(
      -0.26666666666666666*Sqr(g1) - 5.333333333333333*Sqr(g3) + 3*(Sqr(Kappa00) +
      Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) +
      Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*(Sqr(Lambda1200) + Sqr(
      Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 2*Sqr(Lambdax) - 2*Sqr(gN
      )*Sqr(QDxbarp) - 2*Sqr(gN)*Sqr(QDxp) - 2*Sqr(gN)*Sqr(QSp))) + (6*ALambdax*
      Kappa11 + 6*TKappa11)*(2*(Kappa01*(Kappa00*Kappa10 + Kappa01*Kappa11 +
      Kappa02*Kappa12) + Kappa21*(Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*
      Kappa22) + Kappa11*(Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12))) + Kappa11*(
      -0.26666666666666666*Sqr(g1) - 5.333333333333333*Sqr(g3) + 3*(Sqr(Kappa00) +
      Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) +
      Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*(Sqr(Lambda1200) + Sqr(
      Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 2*Sqr(Lambdax) - 2*Sqr(gN
      )*Sqr(QDxbarp) - 2*Sqr(gN)*Sqr(QDxp) - 2*Sqr(gN)*Sqr(QSp))) + (6*ALambdax*
      Kappa12 + 6*TKappa12)*(2*(Kappa02*(Kappa00*Kappa10 + Kappa01*Kappa11 +
      Kappa02*Kappa12) + Kappa22*(Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*
      Kappa22) + Kappa12*(Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12))) + Kappa12*(
      -0.26666666666666666*Sqr(g1) - 5.333333333333333*Sqr(g3) + 3*(Sqr(Kappa00) +
      Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) +
      Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*(Sqr(Lambda1200) + Sqr(
      Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 2*Sqr(Lambdax) - 2*Sqr(gN
      )*Sqr(QDxbarp) - 2*Sqr(gN)*Sqr(QDxp) - 2*Sqr(gN)*Sqr(QSp))) + (6*ALambdax*
      Kappa20 + 6*TKappa20)*(2*(Kappa00*(Kappa00*Kappa20 + Kappa01*Kappa21 +
      Kappa02*Kappa22) + Kappa10*(Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*
      Kappa22) + Kappa20*(Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22))) + Kappa20*(
      -0.26666666666666666*Sqr(g1) - 5.333333333333333*Sqr(g3) + 3*(Sqr(Kappa00) +
      Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) +
      Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*(Sqr(Lambda1200) + Sqr(
      Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 2*Sqr(Lambdax) - 2*Sqr(gN
      )*Sqr(QDxbarp) - 2*Sqr(gN)*Sqr(QDxp) - 2*Sqr(gN)*Sqr(QSp))) + (6*ALambdax*
      Kappa21 + 6*TKappa21)*(2*(Kappa01*(Kappa00*Kappa20 + Kappa01*Kappa21 +
      Kappa02*Kappa22) + Kappa11*(Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*
      Kappa22) + Kappa21*(Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22))) + Kappa21*(
      -0.26666666666666666*Sqr(g1) - 5.333333333333333*Sqr(g3) + 3*(Sqr(Kappa00) +
      Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) +
      Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*(Sqr(Lambda1200) + Sqr(
      Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 2*Sqr(Lambdax) - 2*Sqr(gN
      )*Sqr(QDxbarp) - 2*Sqr(gN)*Sqr(QDxp) - 2*Sqr(gN)*Sqr(QSp))) + (6*ALambdax*
      Kappa22 + 6*TKappa22)*(2*(Kappa02*(Kappa00*Kappa20 + Kappa01*Kappa21 +
      Kappa02*Kappa22) + Kappa12*(Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*
      Kappa22) + Kappa22*(Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22))) + Kappa22*(
      -0.26666666666666666*Sqr(g1) - 5.333333333333333*Sqr(g3) + 3*(Sqr(Kappa00) +
      Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) +
      Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*(Sqr(Lambda1200) + Sqr(
      Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 2*Sqr(Lambdax) - 2*Sqr(gN
      )*Sqr(QDxbarp) - 2*Sqr(gN)*Sqr(QDxp) - 2*Sqr(gN)*Sqr(QSp))) + (4*ALambdax*
      Lambda1200 + 4*TLambda1200)*(2*(Lambda1210*(Lambda1200*Lambda1210 +
      Lambda1201*Lambda1211) + Lambda1200*(Sqr(Lambda1200) + Sqr(Lambda1201))) +
      Lambda1200*(-0.6*Sqr(g1) - 3*Sqr(g2) + 3*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(
      Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(
      Kappa21) + Sqr(Kappa22)) + 2*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(
      Lambda1210) + Sqr(Lambda1211)) + 2*Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH1p) - 2*
      Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QSp))) + (4*ALambdax*Lambda1201 + 4*
      TLambda1201)*(2*(Lambda1211*(Lambda1200*Lambda1210 + Lambda1201*Lambda1211)
      + Lambda1201*(Sqr(Lambda1200) + Sqr(Lambda1201))) + Lambda1201*(-0.6*Sqr(g1)
      - 3*Sqr(g2) + 3*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10)
      + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22))
      + 2*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211))
      + 2*Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr
      (QSp))) + (4*ALambdax*Lambda1210 + 4*TLambda1210)*(2*(Lambda1200*(Lambda1200
      *Lambda1210 + Lambda1201*Lambda1211) + Lambda1210*(Sqr(Lambda1210) + Sqr(
      Lambda1211))) + Lambda1210*(-0.6*Sqr(g1) - 3*Sqr(g2) + 3*(Sqr(Kappa00) + Sqr
      (Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(
      Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*(Sqr(Lambda1200) + Sqr(
      Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 2*Sqr(Lambdax) - 2*Sqr(gN
      )*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QSp))) + (4*ALambdax*
      Lambda1211 + 4*TLambda1211)*(2*(Lambda1201*(Lambda1200*Lambda1210 +
      Lambda1201*Lambda1211) + Lambda1211*(Sqr(Lambda1210) + Sqr(Lambda1211))) +
      Lambda1211*(-0.6*Sqr(g1) - 3*Sqr(g2) + 3*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(
      Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(
      Kappa21) + Sqr(Kappa22)) + 2*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(
      Lambda1210) + Sqr(Lambda1211)) + 2*Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH1p) - 2*
      Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QSp))) + 6*Kappa00*(3*(Kappa00*(Kappa00*
      TKappa00 + Kappa01*TKappa01 + Kappa02*TKappa02) + Kappa10*(Kappa10*TKappa00
      + Kappa11*TKappa01 + Kappa12*TKappa02) + Kappa20*(Kappa20*TKappa00 + Kappa21
      *TKappa01 + Kappa22*TKappa02)) - 0.26666666666666666*TKappa00*Sqr(g1) -
      5.333333333333333*TKappa00*Sqr(g3) + 3*((Kappa00*Kappa10 + Kappa01*Kappa11 +
      Kappa02*Kappa12)*TKappa10 + (Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*
      Kappa22)*TKappa20 + TKappa00*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02))) +
      3*TKappa00*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr
      (Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*
      TKappa00*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(
      Lambda1211)) + 2*TKappa00*Sqr(Lambdax) - 2*TKappa00*Sqr(gN)*Sqr(QDxbarp) - 2
      *TKappa00*Sqr(gN)*Sqr(QDxp) - 2*TKappa00*Sqr(gN)*Sqr(QSp) + Kappa00*(6*(
      Kappa00*TKappa00 + Kappa01*TKappa01 + Kappa02*TKappa02 + Kappa10*TKappa10 +
      Kappa11*TKappa11 + Kappa12*TKappa12 + Kappa20*TKappa20 + Kappa21*TKappa21 +
      Kappa22*TKappa22) + 4*(Lambda1200*TLambda1200 + Lambda1201*TLambda1201 +
      Lambda1210*TLambda1210 + Lambda1211*TLambda1211) + 0.5333333333333333*MassB*
      Sqr(g1) + 10.666666666666666*MassG*Sqr(g3) + 4*ALambdax*Sqr(Lambdax) + 4*
      MassBp*Sqr(gN)*Sqr(QDxbarp) + 4*MassBp*Sqr(gN)*Sqr(QDxp) + 4*MassBp*Sqr(gN)*
      Sqr(QSp))) + 6*Kappa01*(3*(Kappa01*(Kappa00*TKappa00 + Kappa01*TKappa01 +
      Kappa02*TKappa02) + Kappa11*(Kappa10*TKappa00 + Kappa11*TKappa01 + Kappa12*
      TKappa02) + Kappa21*(Kappa20*TKappa00 + Kappa21*TKappa01 + Kappa22*TKappa02)
      ) - 0.26666666666666666*TKappa01*Sqr(g1) - 5.333333333333333*TKappa01*Sqr(g3
      ) + 3*((Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12)*TKappa11 + (
      Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)*TKappa21 + TKappa01*(
      Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02))) + 3*TKappa01*(Sqr(Kappa00) +
      Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) +
      Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*TKappa01*(Sqr(Lambda1200) +
      Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 2*TKappa01*Sqr(
      Lambdax) - 2*TKappa01*Sqr(gN)*Sqr(QDxbarp) - 2*TKappa01*Sqr(gN)*Sqr(QDxp) -
      2*TKappa01*Sqr(gN)*Sqr(QSp) + Kappa01*(6*(Kappa00*TKappa00 + Kappa01*
      TKappa01 + Kappa02*TKappa02 + Kappa10*TKappa10 + Kappa11*TKappa11 + Kappa12*
      TKappa12 + Kappa20*TKappa20 + Kappa21*TKappa21 + Kappa22*TKappa22) + 4*(
      Lambda1200*TLambda1200 + Lambda1201*TLambda1201 + Lambda1210*TLambda1210 +
      Lambda1211*TLambda1211) + 0.5333333333333333*MassB*Sqr(g1) +
      10.666666666666666*MassG*Sqr(g3) + 4*ALambdax*Sqr(Lambdax) + 4*MassBp*Sqr(gN
      )*Sqr(QDxbarp) + 4*MassBp*Sqr(gN)*Sqr(QDxp) + 4*MassBp*Sqr(gN)*Sqr(QSp))) +
      6*Kappa02*(3*(Kappa02*(Kappa00*TKappa00 + Kappa01*TKappa01 + Kappa02*
      TKappa02) + Kappa12*(Kappa10*TKappa00 + Kappa11*TKappa01 + Kappa12*TKappa02)
      + Kappa22*(Kappa20*TKappa00 + Kappa21*TKappa01 + Kappa22*TKappa02)) -
      0.26666666666666666*TKappa02*Sqr(g1) - 5.333333333333333*TKappa02*Sqr(g3) +
      3*((Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12)*TKappa12 + (Kappa00
      *Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)*TKappa22 + TKappa02*(Sqr(
      Kappa00) + Sqr(Kappa01) + Sqr(Kappa02))) + 3*TKappa02*(Sqr(Kappa00) + Sqr(
      Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(
      Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*TKappa02*(Sqr(Lambda1200) + Sqr(
      Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 2*TKappa02*Sqr(Lambdax) -
      2*TKappa02*Sqr(gN)*Sqr(QDxbarp) - 2*TKappa02*Sqr(gN)*Sqr(QDxp) - 2*TKappa02
      *Sqr(gN)*Sqr(QSp) + Kappa02*(6*(Kappa00*TKappa00 + Kappa01*TKappa01 +
      Kappa02*TKappa02 + Kappa10*TKappa10 + Kappa11*TKappa11 + Kappa12*TKappa12 +
      Kappa20*TKappa20 + Kappa21*TKappa21 + Kappa22*TKappa22) + 4*(Lambda1200*
      TLambda1200 + Lambda1201*TLambda1201 + Lambda1210*TLambda1210 + Lambda1211*
      TLambda1211) + 0.5333333333333333*MassB*Sqr(g1) + 10.666666666666666*MassG*
      Sqr(g3) + 4*ALambdax*Sqr(Lambdax) + 4*MassBp*Sqr(gN)*Sqr(QDxbarp) + 4*MassBp
      *Sqr(gN)*Sqr(QDxp) + 4*MassBp*Sqr(gN)*Sqr(QSp))) + 6*Kappa10*(3*(Kappa00*(
      Kappa00*TKappa10 + Kappa01*TKappa11 + Kappa02*TKappa12) + Kappa10*(Kappa10*
      TKappa10 + Kappa11*TKappa11 + Kappa12*TKappa12) + Kappa20*(Kappa20*TKappa10
      + Kappa21*TKappa11 + Kappa22*TKappa12)) - 0.26666666666666666*TKappa10*Sqr(
      g1) - 5.333333333333333*TKappa10*Sqr(g3) + 3*((Kappa00*Kappa10 + Kappa01*
      Kappa11 + Kappa02*Kappa12)*TKappa00 + (Kappa10*Kappa20 + Kappa11*Kappa21 +
      Kappa12*Kappa22)*TKappa20 + TKappa10*(Sqr(Kappa10) + Sqr(Kappa11) + Sqr(
      Kappa12))) + 3*TKappa10*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(
      Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(
      Kappa22)) + 2*TKappa10*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210)
      + Sqr(Lambda1211)) + 2*TKappa10*Sqr(Lambdax) - 2*TKappa10*Sqr(gN)*Sqr(
      QDxbarp) - 2*TKappa10*Sqr(gN)*Sqr(QDxp) - 2*TKappa10*Sqr(gN)*Sqr(QSp) +
      Kappa10*(6*(Kappa00*TKappa00 + Kappa01*TKappa01 + Kappa02*TKappa02 + Kappa10
      *TKappa10 + Kappa11*TKappa11 + Kappa12*TKappa12 + Kappa20*TKappa20 + Kappa21
      *TKappa21 + Kappa22*TKappa22) + 4*(Lambda1200*TLambda1200 + Lambda1201*
      TLambda1201 + Lambda1210*TLambda1210 + Lambda1211*TLambda1211) +
      0.5333333333333333*MassB*Sqr(g1) + 10.666666666666666*MassG*Sqr(g3) + 4*
      ALambdax*Sqr(Lambdax) + 4*MassBp*Sqr(gN)*Sqr(QDxbarp) + 4*MassBp*Sqr(gN)*Sqr
      (QDxp) + 4*MassBp*Sqr(gN)*Sqr(QSp))) + 6*Kappa11*(3*(Kappa01*(Kappa00*
      TKappa10 + Kappa01*TKappa11 + Kappa02*TKappa12) + Kappa11*(Kappa10*TKappa10
      + Kappa11*TKappa11 + Kappa12*TKappa12) + Kappa21*(Kappa20*TKappa10 + Kappa21
      *TKappa11 + Kappa22*TKappa12)) - 0.26666666666666666*TKappa11*Sqr(g1) -
      5.333333333333333*TKappa11*Sqr(g3) + 3*((Kappa00*Kappa10 + Kappa01*Kappa11 +
      Kappa02*Kappa12)*TKappa01 + (Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*
      Kappa22)*TKappa21 + TKappa11*(Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12))) +
      3*TKappa11*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr
      (Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*
      TKappa11*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(
      Lambda1211)) + 2*TKappa11*Sqr(Lambdax) - 2*TKappa11*Sqr(gN)*Sqr(QDxbarp) - 2
      *TKappa11*Sqr(gN)*Sqr(QDxp) - 2*TKappa11*Sqr(gN)*Sqr(QSp) + Kappa11*(6*(
      Kappa00*TKappa00 + Kappa01*TKappa01 + Kappa02*TKappa02 + Kappa10*TKappa10 +
      Kappa11*TKappa11 + Kappa12*TKappa12 + Kappa20*TKappa20 + Kappa21*TKappa21 +
      Kappa22*TKappa22) + 4*(Lambda1200*TLambda1200 + Lambda1201*TLambda1201 +
      Lambda1210*TLambda1210 + Lambda1211*TLambda1211) + 0.5333333333333333*MassB*
      Sqr(g1) + 10.666666666666666*MassG*Sqr(g3) + 4*ALambdax*Sqr(Lambdax) + 4*
      MassBp*Sqr(gN)*Sqr(QDxbarp) + 4*MassBp*Sqr(gN)*Sqr(QDxp) + 4*MassBp*Sqr(gN)*
      Sqr(QSp))) + 6*Kappa12*(3*(Kappa02*(Kappa00*TKappa10 + Kappa01*TKappa11 +
      Kappa02*TKappa12) + Kappa12*(Kappa10*TKappa10 + Kappa11*TKappa11 + Kappa12*
      TKappa12) + Kappa22*(Kappa20*TKappa10 + Kappa21*TKappa11 + Kappa22*TKappa12)
      ) - 0.26666666666666666*TKappa12*Sqr(g1) - 5.333333333333333*TKappa12*Sqr(g3
      ) + 3*((Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12)*TKappa02 + (
      Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)*TKappa22 + TKappa12*(
      Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12))) + 3*TKappa12*(Sqr(Kappa00) +
      Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) +
      Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*TKappa12*(Sqr(Lambda1200) +
      Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 2*TKappa12*Sqr(
      Lambdax) - 2*TKappa12*Sqr(gN)*Sqr(QDxbarp) - 2*TKappa12*Sqr(gN)*Sqr(QDxp) -
      2*TKappa12*Sqr(gN)*Sqr(QSp) + Kappa12*(6*(Kappa00*TKappa00 + Kappa01*
      TKappa01 + Kappa02*TKappa02 + Kappa10*TKappa10 + Kappa11*TKappa11 + Kappa12*
      TKappa12 + Kappa20*TKappa20 + Kappa21*TKappa21 + Kappa22*TKappa22) + 4*(
      Lambda1200*TLambda1200 + Lambda1201*TLambda1201 + Lambda1210*TLambda1210 +
      Lambda1211*TLambda1211) + 0.5333333333333333*MassB*Sqr(g1) +
      10.666666666666666*MassG*Sqr(g3) + 4*ALambdax*Sqr(Lambdax) + 4*MassBp*Sqr(gN
      )*Sqr(QDxbarp) + 4*MassBp*Sqr(gN)*Sqr(QDxp) + 4*MassBp*Sqr(gN)*Sqr(QSp))) +
      6*Kappa20*(3*(Kappa00*(Kappa00*TKappa20 + Kappa01*TKappa21 + Kappa02*
      TKappa22) + Kappa10*(Kappa10*TKappa20 + Kappa11*TKappa21 + Kappa12*TKappa22)
      + Kappa20*(Kappa20*TKappa20 + Kappa21*TKappa21 + Kappa22*TKappa22)) -
      0.26666666666666666*TKappa20*Sqr(g1) - 5.333333333333333*TKappa20*Sqr(g3) +
      3*TKappa20*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(
      Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 3*((
      Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)*TKappa00 + (Kappa10*
      Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)*TKappa10 + TKappa20*(Sqr(
      Kappa20) + Sqr(Kappa21) + Sqr(Kappa22))) + 2*TKappa20*(Sqr(Lambda1200) + Sqr
      (Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 2*TKappa20*Sqr(Lambdax)
      - 2*TKappa20*Sqr(gN)*Sqr(QDxbarp) - 2*TKappa20*Sqr(gN)*Sqr(QDxp) - 2*
      TKappa20*Sqr(gN)*Sqr(QSp) + Kappa20*(6*(Kappa00*TKappa00 + Kappa01*TKappa01
      + Kappa02*TKappa02 + Kappa10*TKappa10 + Kappa11*TKappa11 + Kappa12*TKappa12
      + Kappa20*TKappa20 + Kappa21*TKappa21 + Kappa22*TKappa22) + 4*(Lambda1200*
      TLambda1200 + Lambda1201*TLambda1201 + Lambda1210*TLambda1210 + Lambda1211*
      TLambda1211) + 0.5333333333333333*MassB*Sqr(g1) + 10.666666666666666*MassG*
      Sqr(g3) + 4*ALambdax*Sqr(Lambdax) + 4*MassBp*Sqr(gN)*Sqr(QDxbarp) + 4*MassBp
      *Sqr(gN)*Sqr(QDxp) + 4*MassBp*Sqr(gN)*Sqr(QSp))) + 6*Kappa21*(3*(Kappa01*(
      Kappa00*TKappa20 + Kappa01*TKappa21 + Kappa02*TKappa22) + Kappa11*(Kappa10*
      TKappa20 + Kappa11*TKappa21 + Kappa12*TKappa22) + Kappa21*(Kappa20*TKappa20
      + Kappa21*TKappa21 + Kappa22*TKappa22)) - 0.26666666666666666*TKappa21*Sqr(
      g1) - 5.333333333333333*TKappa21*Sqr(g3) + 3*TKappa21*(Sqr(Kappa00) + Sqr(
      Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(
      Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 3*((Kappa00*Kappa20 + Kappa01*
      Kappa21 + Kappa02*Kappa22)*TKappa01 + (Kappa10*Kappa20 + Kappa11*Kappa21 +
      Kappa12*Kappa22)*TKappa11 + TKappa21*(Sqr(Kappa20) + Sqr(Kappa21) + Sqr(
      Kappa22))) + 2*TKappa21*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210)
      + Sqr(Lambda1211)) + 2*TKappa21*Sqr(Lambdax) - 2*TKappa21*Sqr(gN)*Sqr(
      QDxbarp) - 2*TKappa21*Sqr(gN)*Sqr(QDxp) - 2*TKappa21*Sqr(gN)*Sqr(QSp) +
      Kappa21*(6*(Kappa00*TKappa00 + Kappa01*TKappa01 + Kappa02*TKappa02 + Kappa10
      *TKappa10 + Kappa11*TKappa11 + Kappa12*TKappa12 + Kappa20*TKappa20 + Kappa21
      *TKappa21 + Kappa22*TKappa22) + 4*(Lambda1200*TLambda1200 + Lambda1201*
      TLambda1201 + Lambda1210*TLambda1210 + Lambda1211*TLambda1211) +
      0.5333333333333333*MassB*Sqr(g1) + 10.666666666666666*MassG*Sqr(g3) + 4*
      ALambdax*Sqr(Lambdax) + 4*MassBp*Sqr(gN)*Sqr(QDxbarp) + 4*MassBp*Sqr(gN)*Sqr
      (QDxp) + 4*MassBp*Sqr(gN)*Sqr(QSp))) + 6*Kappa22*(3*(Kappa02*(Kappa00*
      TKappa20 + Kappa01*TKappa21 + Kappa02*TKappa22) + Kappa12*(Kappa10*TKappa20
      + Kappa11*TKappa21 + Kappa12*TKappa22) + Kappa22*(Kappa20*TKappa20 + Kappa21
      *TKappa21 + Kappa22*TKappa22)) - 0.26666666666666666*TKappa22*Sqr(g1) -
      5.333333333333333*TKappa22*Sqr(g3) + 3*TKappa22*(Sqr(Kappa00) + Sqr(Kappa01)
      + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20)
      + Sqr(Kappa21) + Sqr(Kappa22)) + 3*((Kappa00*Kappa20 + Kappa01*Kappa21 +
      Kappa02*Kappa22)*TKappa02 + (Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*
      Kappa22)*TKappa12 + TKappa22*(Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22))) +
      2*TKappa22*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(
      Lambda1211)) + 2*TKappa22*Sqr(Lambdax) - 2*TKappa22*Sqr(gN)*Sqr(QDxbarp) - 2
      *TKappa22*Sqr(gN)*Sqr(QDxp) - 2*TKappa22*Sqr(gN)*Sqr(QSp) + Kappa22*(6*(
      Kappa00*TKappa00 + Kappa01*TKappa01 + Kappa02*TKappa02 + Kappa10*TKappa10 +
      Kappa11*TKappa11 + Kappa12*TKappa12 + Kappa20*TKappa20 + Kappa21*TKappa21 +
      Kappa22*TKappa22) + 4*(Lambda1200*TLambda1200 + Lambda1201*TLambda1201 +
      Lambda1210*TLambda1210 + Lambda1211*TLambda1211) + 0.5333333333333333*MassB*
      Sqr(g1) + 10.666666666666666*MassG*Sqr(g3) + 4*ALambdax*Sqr(Lambdax) + 4*
      MassBp*Sqr(gN)*Sqr(QDxbarp) + 4*MassBp*Sqr(gN)*Sqr(QDxp) + 4*MassBp*Sqr(gN)*
      Sqr(QSp))) + 4*Lambda1200*(3*(Lambda1200*(Lambda1200*TLambda1200 +
      Lambda1201*TLambda1201) + Lambda1210*(Lambda1210*TLambda1200 + Lambda1211*
      TLambda1201)) - 0.6*TLambda1200*Sqr(g1) - 3*TLambda1200*Sqr(g2) + 3*
      TLambda1200*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr
      (Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 3*(
      (Lambda1200*Lambda1210 + Lambda1201*Lambda1211)*TLambda1210 + TLambda1200*(
      Sqr(Lambda1200) + Sqr(Lambda1201))) + 2*TLambda1200*(Sqr(Lambda1200) + Sqr(
      Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 2*TLambda1200*Sqr(Lambdax
      ) - 2*TLambda1200*Sqr(gN)*Sqr(QH1p) - 2*TLambda1200*Sqr(gN)*Sqr(QH2p) - 2*
      TLambda1200*Sqr(gN)*Sqr(QSp) + Lambda1200*(6*(Kappa00*TKappa00 + Kappa01*
      TKappa01 + Kappa02*TKappa02 + Kappa10*TKappa10 + Kappa11*TKappa11 + Kappa12*
      TKappa12 + Kappa20*TKappa20 + Kappa21*TKappa21 + Kappa22*TKappa22) + 4*(
      Lambda1200*TLambda1200 + Lambda1201*TLambda1201 + Lambda1210*TLambda1210 +
      Lambda1211*TLambda1211) + 1.2*MassB*Sqr(g1) + 6*MassWB*Sqr(g2) + 4*ALambdax*
      Sqr(Lambdax) + 4*MassBp*Sqr(gN)*Sqr(QH1p) + 4*MassBp*Sqr(gN)*Sqr(QH2p) + 4*
      MassBp*Sqr(gN)*Sqr(QSp))) + 4*Lambda1201*(3*(Lambda1201*(Lambda1200*
      TLambda1200 + Lambda1201*TLambda1201) + Lambda1211*(Lambda1210*TLambda1200 +
      Lambda1211*TLambda1201)) - 0.6*TLambda1201*Sqr(g1) - 3*TLambda1201*Sqr(g2)
      + 3*TLambda1201*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) +
      Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) +
      3*((Lambda1200*Lambda1210 + Lambda1201*Lambda1211)*TLambda1211 +
      TLambda1201*(Sqr(Lambda1200) + Sqr(Lambda1201))) + 2*TLambda1201*(Sqr(
      Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 2*
      TLambda1201*Sqr(Lambdax) - 2*TLambda1201*Sqr(gN)*Sqr(QH1p) - 2*TLambda1201*
      Sqr(gN)*Sqr(QH2p) - 2*TLambda1201*Sqr(gN)*Sqr(QSp) + Lambda1201*(6*(Kappa00*
      TKappa00 + Kappa01*TKappa01 + Kappa02*TKappa02 + Kappa10*TKappa10 + Kappa11*
      TKappa11 + Kappa12*TKappa12 + Kappa20*TKappa20 + Kappa21*TKappa21 + Kappa22*
      TKappa22) + 4*(Lambda1200*TLambda1200 + Lambda1201*TLambda1201 + Lambda1210*
      TLambda1210 + Lambda1211*TLambda1211) + 1.2*MassB*Sqr(g1) + 6*MassWB*Sqr(g2)
      + 4*ALambdax*Sqr(Lambdax) + 4*MassBp*Sqr(gN)*Sqr(QH1p) + 4*MassBp*Sqr(gN)*
      Sqr(QH2p) + 4*MassBp*Sqr(gN)*Sqr(QSp))) + 4*Lambda1210*(3*(Lambda1200*(
      Lambda1200*TLambda1210 + Lambda1201*TLambda1211) + Lambda1210*(Lambda1210*
      TLambda1210 + Lambda1211*TLambda1211)) - 0.6*TLambda1210*Sqr(g1) - 3*
      TLambda1210*Sqr(g2) + 3*TLambda1210*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(
      Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(
      Kappa21) + Sqr(Kappa22)) + 2*TLambda1210*(Sqr(Lambda1200) + Sqr(Lambda1201)
      + Sqr(Lambda1210) + Sqr(Lambda1211)) + 3*((Lambda1200*Lambda1210 +
      Lambda1201*Lambda1211)*TLambda1200 + TLambda1210*(Sqr(Lambda1210) + Sqr(
      Lambda1211))) + 2*TLambda1210*Sqr(Lambdax) - 2*TLambda1210*Sqr(gN)*Sqr(QH1p)
      - 2*TLambda1210*Sqr(gN)*Sqr(QH2p) - 2*TLambda1210*Sqr(gN)*Sqr(QSp) +
      Lambda1210*(6*(Kappa00*TKappa00 + Kappa01*TKappa01 + Kappa02*TKappa02 +
      Kappa10*TKappa10 + Kappa11*TKappa11 + Kappa12*TKappa12 + Kappa20*TKappa20 +
      Kappa21*TKappa21 + Kappa22*TKappa22) + 4*(Lambda1200*TLambda1200 +
      Lambda1201*TLambda1201 + Lambda1210*TLambda1210 + Lambda1211*TLambda1211) +
      1.2*MassB*Sqr(g1) + 6*MassWB*Sqr(g2) + 4*ALambdax*Sqr(Lambdax) + 4*MassBp*
      Sqr(gN)*Sqr(QH1p) + 4*MassBp*Sqr(gN)*Sqr(QH2p) + 4*MassBp*Sqr(gN)*Sqr(QSp)))
      + 4*Lambda1211*(3*(Lambda1201*(Lambda1200*TLambda1210 + Lambda1201*
      TLambda1211) + Lambda1211*(Lambda1210*TLambda1210 + Lambda1211*TLambda1211))
      - 0.6*TLambda1211*Sqr(g1) - 3*TLambda1211*Sqr(g2) + 3*TLambda1211*(Sqr(
      Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(
      Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*TLambda1211*(Sqr(
      Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 3*((
      Lambda1200*Lambda1210 + Lambda1201*Lambda1211)*TLambda1201 + TLambda1211*(
      Sqr(Lambda1210) + Sqr(Lambda1211))) + 2*TLambda1211*Sqr(Lambdax) - 2*
      TLambda1211*Sqr(gN)*Sqr(QH1p) - 2*TLambda1211*Sqr(gN)*Sqr(QH2p) - 2*
      TLambda1211*Sqr(gN)*Sqr(QSp) + Lambda1211*(6*(Kappa00*TKappa00 + Kappa01*
      TKappa01 + Kappa02*TKappa02 + Kappa10*TKappa10 + Kappa11*TKappa11 + Kappa12*
      TKappa12 + Kappa20*TKappa20 + Kappa21*TKappa21 + Kappa22*TKappa22) + 4*(
      Lambda1200*TLambda1200 + Lambda1201*TLambda1201 + Lambda1210*TLambda1210 +
      Lambda1211*TLambda1211) + 1.2*MassB*Sqr(g1) + 6*MassWB*Sqr(g2) + 4*ALambdax*
      Sqr(Lambdax) + 4*MassBp*Sqr(gN)*Sqr(QH1p) + 4*MassBp*Sqr(gN)*Sqr(QH2p) + 4*
      MassBp*Sqr(gN)*Sqr(QSp))) + 2*MassBp*Sqr(gN)*(4*Sqr(gN)*Sqr(QH1p) + 4*Sqr(gN
      )*Sqr(QH2p) + 4*Sqr(gN)*Sqr(QSp))*(9*Sqr(Qdp) + 9*Sqr(QDxbarp) + 9*Sqr(QDxp)
      + 3*Sqr(Qep) + 6*Sqr(QH1p) + 6*Sqr(QH2p) + 2*Sqr(QHpbarp) + 2*Sqr(QHpp) + 6
      *Sqr(QLp) + 18*Sqr(QQp) + 3*Sqr(QSp) + 9*Sqr(Qup)) + Power(gN,3)*(8*gN*
      MassBp*Sqr(QH1p) + 8*gN*MassBp*Sqr(QH2p) + 8*gN*MassBp*Sqr(QSp) + ALambdax*(
      -4*gN*Sqr(QH1p) - 4*gN*Sqr(QH2p) - 4*gN*Sqr(QSp)))*(9*Sqr(Qdp) + 9*Sqr(
      QDxbarp) + 9*Sqr(QDxp) + 3*Sqr(Qep) + 6*Sqr(QH1p) + 6*Sqr(QH2p) + 2*Sqr(
      QHpbarp) + 2*Sqr(QHpp) + 6*Sqr(QLp) + 18*Sqr(QQp) + 3*Sqr(QSp) + 9*Sqr(Qup))
      + 6*Yd00*(5*(Yd00*(TYd00*Yd00 + TYd01*Yd01 + TYd02*Yd02) + Yd10*(TYd00*Yd10
      + TYd01*Yd11 + TYd02*Yd12) + Yd20*(TYd00*Yd20 + TYd01*Yd21 + TYd02*Yd22)) +
      Yu00*(TYd00*Yu00 + TYd01*Yu01 + TYd02*Yu02) + Yu10*(TYd00*Yu10 + TYd01*Yu11
      + TYd02*Yu12) + Yu20*(TYd00*Yu20 + TYd01*Yu21 + TYd02*Yu22) + 2*(TYu00*(
      Yd00*Yu00 + Yd01*Yu01 + Yd02*Yu02) + TYu10*(Yd00*Yu10 + Yd01*Yu11 + Yd02*
      Yu12) + TYu20*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22)) - 0.4666666666666667*
      TYd00*Sqr(g1) - 3*TYd00*Sqr(g2) - 5.333333333333333*TYd00*Sqr(g3) + TYd00*
      Sqr(Lambdax) - 2*TYd00*Sqr(gN)*Sqr(Qdp) - 2*TYd00*Sqr(gN)*Sqr(QH1p) - 2*
      TYd00*Sqr(gN)*Sqr(QQp) + Yd00*(6*(TYd00*Yd00 + TYd01*Yd01 + TYd02*Yd02 +
      TYd10*Yd10 + TYd11*Yd11 + TYd12*Yd12 + TYd20*Yd20 + TYd21*Yd21 + TYd22*Yd22)
      + 2*(TYe00*Ye00 + TYe01*Ye01 + TYe02*Ye02 + TYe10*Ye10 + TYe11*Ye11 + TYe12
      *Ye12 + TYe20*Ye20 + TYe21*Ye21 + TYe22*Ye22) + 0.9333333333333333*MassB*Sqr
      (g1) + 6*MassWB*Sqr(g2) + 10.666666666666666*MassG*Sqr(g3) + 2*ALambdax*Sqr(
      Lambdax) + 4*MassBp*Sqr(gN)*Sqr(Qdp) + 4*MassBp*Sqr(gN)*Sqr(QH1p) + 4*MassBp
      *Sqr(gN)*Sqr(QQp)) + 4*(TYd10*(Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + TYd20*(
      Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22) + TYd00*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02
      ))) + 3*TYd00*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) +
      Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + TYd00*(Sqr(Ye00) + Sqr(Ye01
      ) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) +
      Sqr(Ye22))) + 6*Yd01*(5*(Yd01*(TYd00*Yd00 + TYd01*Yd01 + TYd02*Yd02) + Yd11*
      (TYd00*Yd10 + TYd01*Yd11 + TYd02*Yd12) + Yd21*(TYd00*Yd20 + TYd01*Yd21 +
      TYd02*Yd22)) + Yu01*(TYd00*Yu00 + TYd01*Yu01 + TYd02*Yu02) + Yu11*(TYd00*
      Yu10 + TYd01*Yu11 + TYd02*Yu12) + Yu21*(TYd00*Yu20 + TYd01*Yu21 + TYd02*Yu22
      ) + 2*(TYu01*(Yd00*Yu00 + Yd01*Yu01 + Yd02*Yu02) + TYu11*(Yd00*Yu10 + Yd01*
      Yu11 + Yd02*Yu12) + TYu21*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22)) -
      0.4666666666666667*TYd01*Sqr(g1) - 3*TYd01*Sqr(g2) - 5.333333333333333*TYd01
      *Sqr(g3) + TYd01*Sqr(Lambdax) - 2*TYd01*Sqr(gN)*Sqr(Qdp) - 2*TYd01*Sqr(gN)*
      Sqr(QH1p) - 2*TYd01*Sqr(gN)*Sqr(QQp) + Yd01*(6*(TYd00*Yd00 + TYd01*Yd01 +
      TYd02*Yd02 + TYd10*Yd10 + TYd11*Yd11 + TYd12*Yd12 + TYd20*Yd20 + TYd21*Yd21
      + TYd22*Yd22) + 2*(TYe00*Ye00 + TYe01*Ye01 + TYe02*Ye02 + TYe10*Ye10 + TYe11
      *Ye11 + TYe12*Ye12 + TYe20*Ye20 + TYe21*Ye21 + TYe22*Ye22) +
      0.9333333333333333*MassB*Sqr(g1) + 6*MassWB*Sqr(g2) + 10.666666666666666*
      MassG*Sqr(g3) + 2*ALambdax*Sqr(Lambdax) + 4*MassBp*Sqr(gN)*Sqr(Qdp) + 4*
      MassBp*Sqr(gN)*Sqr(QH1p) + 4*MassBp*Sqr(gN)*Sqr(QQp)) + 4*(TYd11*(Yd00*Yd10
      + Yd01*Yd11 + Yd02*Yd12) + TYd21*(Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22) + TYd01
      *(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02))) + 3*TYd01*(Sqr(Yd00) + Sqr(Yd01) + Sqr
      (Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(
      Yd22)) + TYd01*(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) +
      Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22))) + 6*Yd02*(5*(Yd02*(TYd00*
      Yd00 + TYd01*Yd01 + TYd02*Yd02) + Yd12*(TYd00*Yd10 + TYd01*Yd11 + TYd02*Yd12
      ) + Yd22*(TYd00*Yd20 + TYd01*Yd21 + TYd02*Yd22)) + Yu02*(TYd00*Yu00 + TYd01*
      Yu01 + TYd02*Yu02) + Yu12*(TYd00*Yu10 + TYd01*Yu11 + TYd02*Yu12) + Yu22*(
      TYd00*Yu20 + TYd01*Yu21 + TYd02*Yu22) + 2*(TYu02*(Yd00*Yu00 + Yd01*Yu01 +
      Yd02*Yu02) + TYu12*(Yd00*Yu10 + Yd01*Yu11 + Yd02*Yu12) + TYu22*(Yd00*Yu20 +
      Yd01*Yu21 + Yd02*Yu22)) - 0.4666666666666667*TYd02*Sqr(g1) - 3*TYd02*Sqr(g2)
      - 5.333333333333333*TYd02*Sqr(g3) + TYd02*Sqr(Lambdax) - 2*TYd02*Sqr(gN)*
      Sqr(Qdp) - 2*TYd02*Sqr(gN)*Sqr(QH1p) - 2*TYd02*Sqr(gN)*Sqr(QQp) + Yd02*(6*(
      TYd00*Yd00 + TYd01*Yd01 + TYd02*Yd02 + TYd10*Yd10 + TYd11*Yd11 + TYd12*Yd12
      + TYd20*Yd20 + TYd21*Yd21 + TYd22*Yd22) + 2*(TYe00*Ye00 + TYe01*Ye01 + TYe02
      *Ye02 + TYe10*Ye10 + TYe11*Ye11 + TYe12*Ye12 + TYe20*Ye20 + TYe21*Ye21 +
      TYe22*Ye22) + 0.9333333333333333*MassB*Sqr(g1) + 6*MassWB*Sqr(g2) +
      10.666666666666666*MassG*Sqr(g3) + 2*ALambdax*Sqr(Lambdax) + 4*MassBp*Sqr(gN
      )*Sqr(Qdp) + 4*MassBp*Sqr(gN)*Sqr(QH1p) + 4*MassBp*Sqr(gN)*Sqr(QQp)) + 4*(
      TYd12*(Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + TYd22*(Yd00*Yd20 + Yd01*Yd21 +
      Yd02*Yd22) + TYd02*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02))) + 3*TYd02*(Sqr(Yd00)
      + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) +
      Sqr(Yd21) + Sqr(Yd22)) + TYd02*(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10
      ) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22))) + 6*Yd10*(5*
      (Yd00*(TYd10*Yd00 + TYd11*Yd01 + TYd12*Yd02) + Yd10*(TYd10*Yd10 + TYd11*Yd11
      + TYd12*Yd12) + Yd20*(TYd10*Yd20 + TYd11*Yd21 + TYd12*Yd22)) + Yu00*(TYd10*
      Yu00 + TYd11*Yu01 + TYd12*Yu02) + Yu10*(TYd10*Yu10 + TYd11*Yu11 + TYd12*Yu12
      ) + Yu20*(TYd10*Yu20 + TYd11*Yu21 + TYd12*Yu22) + 2*(TYu00*(Yd10*Yu00 + Yd11
      *Yu01 + Yd12*Yu02) + TYu10*(Yd10*Yu10 + Yd11*Yu11 + Yd12*Yu12) + TYu20*(Yd10
      *Yu20 + Yd11*Yu21 + Yd12*Yu22)) - 0.4666666666666667*TYd10*Sqr(g1) - 3*TYd10
      *Sqr(g2) - 5.333333333333333*TYd10*Sqr(g3) + TYd10*Sqr(Lambdax) - 2*TYd10*
      Sqr(gN)*Sqr(Qdp) - 2*TYd10*Sqr(gN)*Sqr(QH1p) - 2*TYd10*Sqr(gN)*Sqr(QQp) +
      Yd10*(6*(TYd00*Yd00 + TYd01*Yd01 + TYd02*Yd02 + TYd10*Yd10 + TYd11*Yd11 +
      TYd12*Yd12 + TYd20*Yd20 + TYd21*Yd21 + TYd22*Yd22) + 2*(TYe00*Ye00 + TYe01*
      Ye01 + TYe02*Ye02 + TYe10*Ye10 + TYe11*Ye11 + TYe12*Ye12 + TYe20*Ye20 +
      TYe21*Ye21 + TYe22*Ye22) + 0.9333333333333333*MassB*Sqr(g1) + 6*MassWB*Sqr(
      g2) + 10.666666666666666*MassG*Sqr(g3) + 2*ALambdax*Sqr(Lambdax) + 4*MassBp*
      Sqr(gN)*Sqr(Qdp) + 4*MassBp*Sqr(gN)*Sqr(QH1p) + 4*MassBp*Sqr(gN)*Sqr(QQp)) +
      4*(TYd00*(Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + TYd20*(Yd10*Yd20 + Yd11*Yd21
      + Yd12*Yd22) + TYd10*(Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12))) + 3*TYd10*(Sqr(
      Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20
      ) + Sqr(Yd21) + Sqr(Yd22)) + TYd10*(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(
      Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22))) + 6*Yd11
      *(5*(Yd01*(TYd10*Yd00 + TYd11*Yd01 + TYd12*Yd02) + Yd11*(TYd10*Yd10 + TYd11*
      Yd11 + TYd12*Yd12) + Yd21*(TYd10*Yd20 + TYd11*Yd21 + TYd12*Yd22)) + Yu01*(
      TYd10*Yu00 + TYd11*Yu01 + TYd12*Yu02) + Yu11*(TYd10*Yu10 + TYd11*Yu11 +
      TYd12*Yu12) + Yu21*(TYd10*Yu20 + TYd11*Yu21 + TYd12*Yu22) + 2*(TYu01*(Yd10*
      Yu00 + Yd11*Yu01 + Yd12*Yu02) + TYu11*(Yd10*Yu10 + Yd11*Yu11 + Yd12*Yu12) +
      TYu21*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22)) - 0.4666666666666667*TYd11*Sqr(g1
      ) - 3*TYd11*Sqr(g2) - 5.333333333333333*TYd11*Sqr(g3) + TYd11*Sqr(Lambdax) -
      2*TYd11*Sqr(gN)*Sqr(Qdp) - 2*TYd11*Sqr(gN)*Sqr(QH1p) - 2*TYd11*Sqr(gN)*Sqr(
      QQp) + Yd11*(6*(TYd00*Yd00 + TYd01*Yd01 + TYd02*Yd02 + TYd10*Yd10 + TYd11*
      Yd11 + TYd12*Yd12 + TYd20*Yd20 + TYd21*Yd21 + TYd22*Yd22) + 2*(TYe00*Ye00 +
      TYe01*Ye01 + TYe02*Ye02 + TYe10*Ye10 + TYe11*Ye11 + TYe12*Ye12 + TYe20*Ye20
      + TYe21*Ye21 + TYe22*Ye22) + 0.9333333333333333*MassB*Sqr(g1) + 6*MassWB*Sqr
      (g2) + 10.666666666666666*MassG*Sqr(g3) + 2*ALambdax*Sqr(Lambdax) + 4*MassBp
      *Sqr(gN)*Sqr(Qdp) + 4*MassBp*Sqr(gN)*Sqr(QH1p) + 4*MassBp*Sqr(gN)*Sqr(QQp))
      + 4*(TYd01*(Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + TYd21*(Yd10*Yd20 + Yd11*
      Yd21 + Yd12*Yd22) + TYd11*(Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12))) + 3*TYd11*(
      Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(
      Yd20) + Sqr(Yd21) + Sqr(Yd22)) + TYd11*(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) +
      Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22))) + 6*
      Yd12*(5*(Yd02*(TYd10*Yd00 + TYd11*Yd01 + TYd12*Yd02) + Yd12*(TYd10*Yd10 +
      TYd11*Yd11 + TYd12*Yd12) + Yd22*(TYd10*Yd20 + TYd11*Yd21 + TYd12*Yd22)) +
      Yu02*(TYd10*Yu00 + TYd11*Yu01 + TYd12*Yu02) + Yu12*(TYd10*Yu10 + TYd11*Yu11
      + TYd12*Yu12) + Yu22*(TYd10*Yu20 + TYd11*Yu21 + TYd12*Yu22) + 2*(TYu02*(Yd10
      *Yu00 + Yd11*Yu01 + Yd12*Yu02) + TYu12*(Yd10*Yu10 + Yd11*Yu11 + Yd12*Yu12) +
      TYu22*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22)) - 0.4666666666666667*TYd12*Sqr(
      g1) - 3*TYd12*Sqr(g2) - 5.333333333333333*TYd12*Sqr(g3) + TYd12*Sqr(Lambdax)
      - 2*TYd12*Sqr(gN)*Sqr(Qdp) - 2*TYd12*Sqr(gN)*Sqr(QH1p) - 2*TYd12*Sqr(gN)*
      Sqr(QQp) + Yd12*(6*(TYd00*Yd00 + TYd01*Yd01 + TYd02*Yd02 + TYd10*Yd10 +
      TYd11*Yd11 + TYd12*Yd12 + TYd20*Yd20 + TYd21*Yd21 + TYd22*Yd22) + 2*(TYe00*
      Ye00 + TYe01*Ye01 + TYe02*Ye02 + TYe10*Ye10 + TYe11*Ye11 + TYe12*Ye12 +
      TYe20*Ye20 + TYe21*Ye21 + TYe22*Ye22) + 0.9333333333333333*MassB*Sqr(g1) + 6
      *MassWB*Sqr(g2) + 10.666666666666666*MassG*Sqr(g3) + 2*ALambdax*Sqr(Lambdax)
      + 4*MassBp*Sqr(gN)*Sqr(Qdp) + 4*MassBp*Sqr(gN)*Sqr(QH1p) + 4*MassBp*Sqr(gN)
      *Sqr(QQp)) + 4*(TYd02*(Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + TYd22*(Yd10*Yd20
      + Yd11*Yd21 + Yd12*Yd22) + TYd12*(Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12))) + 3*
      TYd12*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12)
      + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + TYd12*(Sqr(Ye00) + Sqr(Ye01) + Sqr(
      Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22
      ))) + 6*Yd20*(5*(Yd00*(TYd20*Yd00 + TYd21*Yd01 + TYd22*Yd02) + Yd10*(TYd20*
      Yd10 + TYd21*Yd11 + TYd22*Yd12) + Yd20*(TYd20*Yd20 + TYd21*Yd21 + TYd22*Yd22
      )) + Yu00*(TYd20*Yu00 + TYd21*Yu01 + TYd22*Yu02) + Yu10*(TYd20*Yu10 + TYd21*
      Yu11 + TYd22*Yu12) + Yu20*(TYd20*Yu20 + TYd21*Yu21 + TYd22*Yu22) + 2*(TYu00*
      (Yd20*Yu00 + Yd21*Yu01 + Yd22*Yu02) + TYu10*(Yd20*Yu10 + Yd21*Yu11 + Yd22*
      Yu12) + TYu20*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22)) - 0.4666666666666667*
      TYd20*Sqr(g1) - 3*TYd20*Sqr(g2) - 5.333333333333333*TYd20*Sqr(g3) + TYd20*
      Sqr(Lambdax) - 2*TYd20*Sqr(gN)*Sqr(Qdp) - 2*TYd20*Sqr(gN)*Sqr(QH1p) - 2*
      TYd20*Sqr(gN)*Sqr(QQp) + Yd20*(6*(TYd00*Yd00 + TYd01*Yd01 + TYd02*Yd02 +
      TYd10*Yd10 + TYd11*Yd11 + TYd12*Yd12 + TYd20*Yd20 + TYd21*Yd21 + TYd22*Yd22)
      + 2*(TYe00*Ye00 + TYe01*Ye01 + TYe02*Ye02 + TYe10*Ye10 + TYe11*Ye11 + TYe12
      *Ye12 + TYe20*Ye20 + TYe21*Ye21 + TYe22*Ye22) + 0.9333333333333333*MassB*Sqr
      (g1) + 6*MassWB*Sqr(g2) + 10.666666666666666*MassG*Sqr(g3) + 2*ALambdax*Sqr(
      Lambdax) + 4*MassBp*Sqr(gN)*Sqr(Qdp) + 4*MassBp*Sqr(gN)*Sqr(QH1p) + 4*MassBp
      *Sqr(gN)*Sqr(QQp)) + 3*TYd20*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10)
      + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + 4*(TYd00*(
      Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22) + TYd10*(Yd10*Yd20 + Yd11*Yd21 + Yd12*
      Yd22) + TYd20*(Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22))) + TYd20*(Sqr(Ye00) + Sqr(
      Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21
      ) + Sqr(Ye22))) + 6*Yd21*(5*(Yd01*(TYd20*Yd00 + TYd21*Yd01 + TYd22*Yd02) +
      Yd11*(TYd20*Yd10 + TYd21*Yd11 + TYd22*Yd12) + Yd21*(TYd20*Yd20 + TYd21*Yd21
      + TYd22*Yd22)) + Yu01*(TYd20*Yu00 + TYd21*Yu01 + TYd22*Yu02) + Yu11*(TYd20*
      Yu10 + TYd21*Yu11 + TYd22*Yu12) + Yu21*(TYd20*Yu20 + TYd21*Yu21 + TYd22*Yu22
      ) + 2*(TYu01*(Yd20*Yu00 + Yd21*Yu01 + Yd22*Yu02) + TYu11*(Yd20*Yu10 + Yd21*
      Yu11 + Yd22*Yu12) + TYu21*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22)) -
      0.4666666666666667*TYd21*Sqr(g1) - 3*TYd21*Sqr(g2) - 5.333333333333333*TYd21
      *Sqr(g3) + TYd21*Sqr(Lambdax) - 2*TYd21*Sqr(gN)*Sqr(Qdp) - 2*TYd21*Sqr(gN)*
      Sqr(QH1p) - 2*TYd21*Sqr(gN)*Sqr(QQp) + Yd21*(6*(TYd00*Yd00 + TYd01*Yd01 +
      TYd02*Yd02 + TYd10*Yd10 + TYd11*Yd11 + TYd12*Yd12 + TYd20*Yd20 + TYd21*Yd21
      + TYd22*Yd22) + 2*(TYe00*Ye00 + TYe01*Ye01 + TYe02*Ye02 + TYe10*Ye10 + TYe11
      *Ye11 + TYe12*Ye12 + TYe20*Ye20 + TYe21*Ye21 + TYe22*Ye22) +
      0.9333333333333333*MassB*Sqr(g1) + 6*MassWB*Sqr(g2) + 10.666666666666666*
      MassG*Sqr(g3) + 2*ALambdax*Sqr(Lambdax) + 4*MassBp*Sqr(gN)*Sqr(Qdp) + 4*
      MassBp*Sqr(gN)*Sqr(QH1p) + 4*MassBp*Sqr(gN)*Sqr(QQp)) + 3*TYd21*(Sqr(Yd00) +
      Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr
      (Yd21) + Sqr(Yd22)) + 4*(TYd01*(Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22) + TYd11*(
      Yd10*Yd20 + Yd11*Yd21 + Yd12*Yd22) + TYd21*(Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22
      ))) + TYd21*(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr
      (Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22))) + 6*Yd22*(5*(Yd02*(TYd20*Yd00 +
      TYd21*Yd01 + TYd22*Yd02) + Yd12*(TYd20*Yd10 + TYd21*Yd11 + TYd22*Yd12) +
      Yd22*(TYd20*Yd20 + TYd21*Yd21 + TYd22*Yd22)) + Yu02*(TYd20*Yu00 + TYd21*Yu01
      + TYd22*Yu02) + Yu12*(TYd20*Yu10 + TYd21*Yu11 + TYd22*Yu12) + Yu22*(TYd20*
      Yu20 + TYd21*Yu21 + TYd22*Yu22) + 2*(TYu02*(Yd20*Yu00 + Yd21*Yu01 + Yd22*
      Yu02) + TYu12*(Yd20*Yu10 + Yd21*Yu11 + Yd22*Yu12) + TYu22*(Yd20*Yu20 + Yd21*
      Yu21 + Yd22*Yu22)) - 0.4666666666666667*TYd22*Sqr(g1) - 3*TYd22*Sqr(g2) -
      5.333333333333333*TYd22*Sqr(g3) + TYd22*Sqr(Lambdax) - 2*TYd22*Sqr(gN)*Sqr(
      Qdp) - 2*TYd22*Sqr(gN)*Sqr(QH1p) - 2*TYd22*Sqr(gN)*Sqr(QQp) + Yd22*(6*(TYd00
      *Yd00 + TYd01*Yd01 + TYd02*Yd02 + TYd10*Yd10 + TYd11*Yd11 + TYd12*Yd12 +
      TYd20*Yd20 + TYd21*Yd21 + TYd22*Yd22) + 2*(TYe00*Ye00 + TYe01*Ye01 + TYe02*
      Ye02 + TYe10*Ye10 + TYe11*Ye11 + TYe12*Ye12 + TYe20*Ye20 + TYe21*Ye21 +
      TYe22*Ye22) + 0.9333333333333333*MassB*Sqr(g1) + 6*MassWB*Sqr(g2) +
      10.666666666666666*MassG*Sqr(g3) + 2*ALambdax*Sqr(Lambdax) + 4*MassBp*Sqr(gN
      )*Sqr(Qdp) + 4*MassBp*Sqr(gN)*Sqr(QH1p) + 4*MassBp*Sqr(gN)*Sqr(QQp)) + 3*
      TYd22*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12)
      + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + 4*(TYd02*(Yd00*Yd20 + Yd01*Yd21 +
      Yd02*Yd22) + TYd12*(Yd10*Yd20 + Yd11*Yd21 + Yd12*Yd22) + TYd22*(Sqr(Yd20) +
      Sqr(Yd21) + Sqr(Yd22))) + TYd22*(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(
      Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22))) + 2*Ye00
      *(5*(Ye00*(TYe00*Ye00 + TYe01*Ye01 + TYe02*Ye02) + Ye10*(TYe00*Ye10 + TYe01*
      Ye11 + TYe02*Ye12) + Ye20*(TYe00*Ye20 + TYe01*Ye21 + TYe02*Ye22)) - 1.8*
      TYe00*Sqr(g1) - 3*TYe00*Sqr(g2) + TYe00*Sqr(Lambdax) - 2*TYe00*Sqr(gN)*Sqr(
      Qep) - 2*TYe00*Sqr(gN)*Sqr(QH1p) - 2*TYe00*Sqr(gN)*Sqr(QLp) + Ye00*(6*(TYd00
      *Yd00 + TYd01*Yd01 + TYd02*Yd02 + TYd10*Yd10 + TYd11*Yd11 + TYd12*Yd12 +
      TYd20*Yd20 + TYd21*Yd21 + TYd22*Yd22) + 2*(TYe00*Ye00 + TYe01*Ye01 + TYe02*
      Ye02 + TYe10*Ye10 + TYe11*Ye11 + TYe12*Ye12 + TYe20*Ye20 + TYe21*Ye21 +
      TYe22*Ye22) + 3.6*MassB*Sqr(g1) + 6*MassWB*Sqr(g2) + 2*ALambdax*Sqr(Lambdax)
      + 4*MassBp*Sqr(gN)*Sqr(Qep) + 4*MassBp*Sqr(gN)*Sqr(QH1p) + 4*MassBp*Sqr(gN)
      *Sqr(QLp)) + 3*TYe00*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(
      Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + 4*(TYe10*(Ye00*Ye10
      + Ye01*Ye11 + Ye02*Ye12) + TYe20*(Ye00*Ye20 + Ye01*Ye21 + Ye02*Ye22) +
      TYe00*(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02))) + TYe00*(Sqr(Ye00) + Sqr(Ye01) +
      Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(
      Ye22))) + 2*Ye01*(5*(Ye01*(TYe00*Ye00 + TYe01*Ye01 + TYe02*Ye02) + Ye11*(
      TYe00*Ye10 + TYe01*Ye11 + TYe02*Ye12) + Ye21*(TYe00*Ye20 + TYe01*Ye21 +
      TYe02*Ye22)) - 1.8*TYe01*Sqr(g1) - 3*TYe01*Sqr(g2) + TYe01*Sqr(Lambdax) - 2*
      TYe01*Sqr(gN)*Sqr(Qep) - 2*TYe01*Sqr(gN)*Sqr(QH1p) - 2*TYe01*Sqr(gN)*Sqr(QLp
      ) + Ye01*(6*(TYd00*Yd00 + TYd01*Yd01 + TYd02*Yd02 + TYd10*Yd10 + TYd11*Yd11
      + TYd12*Yd12 + TYd20*Yd20 + TYd21*Yd21 + TYd22*Yd22) + 2*(TYe00*Ye00 + TYe01
      *Ye01 + TYe02*Ye02 + TYe10*Ye10 + TYe11*Ye11 + TYe12*Ye12 + TYe20*Ye20 +
      TYe21*Ye21 + TYe22*Ye22) + 3.6*MassB*Sqr(g1) + 6*MassWB*Sqr(g2) + 2*ALambdax
      *Sqr(Lambdax) + 4*MassBp*Sqr(gN)*Sqr(Qep) + 4*MassBp*Sqr(gN)*Sqr(QH1p) + 4*
      MassBp*Sqr(gN)*Sqr(QLp)) + 3*TYe01*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(
      Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + 4*(
      TYe11*(Ye00*Ye10 + Ye01*Ye11 + Ye02*Ye12) + TYe21*(Ye00*Ye20 + Ye01*Ye21 +
      Ye02*Ye22) + TYe01*(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02))) + TYe01*(Sqr(Ye00) +
      Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr
      (Ye21) + Sqr(Ye22))) + 2*Ye02*(5*(Ye02*(TYe00*Ye00 + TYe01*Ye01 + TYe02*Ye02
      ) + Ye12*(TYe00*Ye10 + TYe01*Ye11 + TYe02*Ye12) + Ye22*(TYe00*Ye20 + TYe01*
      Ye21 + TYe02*Ye22)) - 1.8*TYe02*Sqr(g1) - 3*TYe02*Sqr(g2) + TYe02*Sqr(
      Lambdax) - 2*TYe02*Sqr(gN)*Sqr(Qep) - 2*TYe02*Sqr(gN)*Sqr(QH1p) - 2*TYe02*
      Sqr(gN)*Sqr(QLp) + Ye02*(6*(TYd00*Yd00 + TYd01*Yd01 + TYd02*Yd02 + TYd10*
      Yd10 + TYd11*Yd11 + TYd12*Yd12 + TYd20*Yd20 + TYd21*Yd21 + TYd22*Yd22) + 2*(
      TYe00*Ye00 + TYe01*Ye01 + TYe02*Ye02 + TYe10*Ye10 + TYe11*Ye11 + TYe12*Ye12
      + TYe20*Ye20 + TYe21*Ye21 + TYe22*Ye22) + 3.6*MassB*Sqr(g1) + 6*MassWB*Sqr(
      g2) + 2*ALambdax*Sqr(Lambdax) + 4*MassBp*Sqr(gN)*Sqr(Qep) + 4*MassBp*Sqr(gN)
      *Sqr(QH1p) + 4*MassBp*Sqr(gN)*Sqr(QLp)) + 3*TYe02*(Sqr(Yd00) + Sqr(Yd01) +
      Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(
      Yd22)) + 4*(TYe12*(Ye00*Ye10 + Ye01*Ye11 + Ye02*Ye12) + TYe22*(Ye00*Ye20 +
      Ye01*Ye21 + Ye02*Ye22) + TYe02*(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02))) + TYe02*
      (Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr
      (Ye20) + Sqr(Ye21) + Sqr(Ye22))) + 2*Ye10*(5*(Ye00*(TYe10*Ye00 + TYe11*Ye01
      + TYe12*Ye02) + Ye10*(TYe10*Ye10 + TYe11*Ye11 + TYe12*Ye12) + Ye20*(TYe10*
      Ye20 + TYe11*Ye21 + TYe12*Ye22)) - 1.8*TYe10*Sqr(g1) - 3*TYe10*Sqr(g2) +
      TYe10*Sqr(Lambdax) - 2*TYe10*Sqr(gN)*Sqr(Qep) - 2*TYe10*Sqr(gN)*Sqr(QH1p) -
      2*TYe10*Sqr(gN)*Sqr(QLp) + Ye10*(6*(TYd00*Yd00 + TYd01*Yd01 + TYd02*Yd02 +
      TYd10*Yd10 + TYd11*Yd11 + TYd12*Yd12 + TYd20*Yd20 + TYd21*Yd21 + TYd22*Yd22)
      + 2*(TYe00*Ye00 + TYe01*Ye01 + TYe02*Ye02 + TYe10*Ye10 + TYe11*Ye11 + TYe12
      *Ye12 + TYe20*Ye20 + TYe21*Ye21 + TYe22*Ye22) + 3.6*MassB*Sqr(g1) + 6*MassWB
      *Sqr(g2) + 2*ALambdax*Sqr(Lambdax) + 4*MassBp*Sqr(gN)*Sqr(Qep) + 4*MassBp*
      Sqr(gN)*Sqr(QH1p) + 4*MassBp*Sqr(gN)*Sqr(QLp)) + 3*TYe10*(Sqr(Yd00) + Sqr(
      Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21
      ) + Sqr(Yd22)) + 4*(TYe00*(Ye00*Ye10 + Ye01*Ye11 + Ye02*Ye12) + TYe20*(Ye10*
      Ye20 + Ye11*Ye21 + Ye12*Ye22) + TYe10*(Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12))) +
      TYe10*(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12
      ) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22))) + 2*Ye11*(5*(Ye01*(TYe10*Ye00 +
      TYe11*Ye01 + TYe12*Ye02) + Ye11*(TYe10*Ye10 + TYe11*Ye11 + TYe12*Ye12) +
      Ye21*(TYe10*Ye20 + TYe11*Ye21 + TYe12*Ye22)) - 1.8*TYe11*Sqr(g1) - 3*TYe11*
      Sqr(g2) + TYe11*Sqr(Lambdax) - 2*TYe11*Sqr(gN)*Sqr(Qep) - 2*TYe11*Sqr(gN)*
      Sqr(QH1p) - 2*TYe11*Sqr(gN)*Sqr(QLp) + Ye11*(6*(TYd00*Yd00 + TYd01*Yd01 +
      TYd02*Yd02 + TYd10*Yd10 + TYd11*Yd11 + TYd12*Yd12 + TYd20*Yd20 + TYd21*Yd21
      + TYd22*Yd22) + 2*(TYe00*Ye00 + TYe01*Ye01 + TYe02*Ye02 + TYe10*Ye10 + TYe11
      *Ye11 + TYe12*Ye12 + TYe20*Ye20 + TYe21*Ye21 + TYe22*Ye22) + 3.6*MassB*Sqr(
      g1) + 6*MassWB*Sqr(g2) + 2*ALambdax*Sqr(Lambdax) + 4*MassBp*Sqr(gN)*Sqr(Qep)
      + 4*MassBp*Sqr(gN)*Sqr(QH1p) + 4*MassBp*Sqr(gN)*Sqr(QLp)) + 3*TYe11*(Sqr(
      Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20
      ) + Sqr(Yd21) + Sqr(Yd22)) + 4*(TYe01*(Ye00*Ye10 + Ye01*Ye11 + Ye02*Ye12) +
      TYe21*(Ye10*Ye20 + Ye11*Ye21 + Ye12*Ye22) + TYe11*(Sqr(Ye10) + Sqr(Ye11) +
      Sqr(Ye12))) + TYe11*(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(
      Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22))) + 2*Ye12*(5*(Ye02*(
      TYe10*Ye00 + TYe11*Ye01 + TYe12*Ye02) + Ye12*(TYe10*Ye10 + TYe11*Ye11 +
      TYe12*Ye12) + Ye22*(TYe10*Ye20 + TYe11*Ye21 + TYe12*Ye22)) - 1.8*TYe12*Sqr(
      g1) - 3*TYe12*Sqr(g2) + TYe12*Sqr(Lambdax) - 2*TYe12*Sqr(gN)*Sqr(Qep) - 2*
      TYe12*Sqr(gN)*Sqr(QH1p) - 2*TYe12*Sqr(gN)*Sqr(QLp) + Ye12*(6*(TYd00*Yd00 +
      TYd01*Yd01 + TYd02*Yd02 + TYd10*Yd10 + TYd11*Yd11 + TYd12*Yd12 + TYd20*Yd20
      + TYd21*Yd21 + TYd22*Yd22) + 2*(TYe00*Ye00 + TYe01*Ye01 + TYe02*Ye02 + TYe10
      *Ye10 + TYe11*Ye11 + TYe12*Ye12 + TYe20*Ye20 + TYe21*Ye21 + TYe22*Ye22) +
      3.6*MassB*Sqr(g1) + 6*MassWB*Sqr(g2) + 2*ALambdax*Sqr(Lambdax) + 4*MassBp*
      Sqr(gN)*Sqr(Qep) + 4*MassBp*Sqr(gN)*Sqr(QH1p) + 4*MassBp*Sqr(gN)*Sqr(QLp)) +
      3*TYe12*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(
      Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + 4*(TYe02*(Ye00*Ye10 + Ye01*Ye11
      + Ye02*Ye12) + TYe22*(Ye10*Ye20 + Ye11*Ye21 + Ye12*Ye22) + TYe12*(Sqr(Ye10)
      + Sqr(Ye11) + Sqr(Ye12))) + TYe12*(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(
      Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22))) + (2*
      TYe00 + 2*ALambdax*Ye00)*(3*(Ye10*(Ye00*Ye10 + Ye01*Ye11 + Ye02*Ye12) + Ye20
      *(Ye00*Ye20 + Ye01*Ye21 + Ye02*Ye22) + Ye00*(Sqr(Ye00) + Sqr(Ye01) + Sqr(
      Ye02))) + Ye00*(-1.8*Sqr(g1) - 3*Sqr(g2) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(Qep)
      - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QLp) + 3*(Sqr(Yd00) + Sqr(Yd01) + Sqr
      (Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(
      Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(
      Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22))) + (2*TYe01 + 2*ALambdax*Ye01)*(3
      *(Ye11*(Ye00*Ye10 + Ye01*Ye11 + Ye02*Ye12) + Ye21*(Ye00*Ye20 + Ye01*Ye21 +
      Ye02*Ye22) + Ye01*(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02))) + Ye01*(-1.8*Sqr(g1)
      - 3*Sqr(g2) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(Qep) - 2*Sqr(gN)*Sqr(QH1p) - 2*
      Sqr(gN)*Sqr(QLp) + 3*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(
      Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(Ye00) + Sqr(
      Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21
      ) + Sqr(Ye22))) + (2*TYe02 + 2*ALambdax*Ye02)*(3*(Ye12*(Ye00*Ye10 + Ye01*
      Ye11 + Ye02*Ye12) + Ye22*(Ye00*Ye20 + Ye01*Ye21 + Ye02*Ye22) + Ye02*(Sqr(
      Ye00) + Sqr(Ye01) + Sqr(Ye02))) + Ye02*(-1.8*Sqr(g1) - 3*Sqr(g2) + Sqr(
      Lambdax) - 2*Sqr(gN)*Sqr(Qep) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QLp) + 3
      *(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) +
      Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr
      (Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22))) + (2*
      TYe10 + 2*ALambdax*Ye10)*(3*(Ye00*(Ye00*Ye10 + Ye01*Ye11 + Ye02*Ye12) + Ye20
      *(Ye10*Ye20 + Ye11*Ye21 + Ye12*Ye22) + Ye10*(Sqr(Ye10) + Sqr(Ye11) + Sqr(
      Ye12))) + Ye10*(-1.8*Sqr(g1) - 3*Sqr(g2) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(Qep)
      - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QLp) + 3*(Sqr(Yd00) + Sqr(Yd01) + Sqr
      (Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(
      Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(
      Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22))) + (2*TYe11 + 2*ALambdax*Ye11)*(3
      *(Ye01*(Ye00*Ye10 + Ye01*Ye11 + Ye02*Ye12) + Ye21*(Ye10*Ye20 + Ye11*Ye21 +
      Ye12*Ye22) + Ye11*(Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12))) + Ye11*(-1.8*Sqr(g1)
      - 3*Sqr(g2) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(Qep) - 2*Sqr(gN)*Sqr(QH1p) - 2*
      Sqr(gN)*Sqr(QLp) + 3*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(
      Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(Ye00) + Sqr(
      Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21
      ) + Sqr(Ye22))) + (2*TYe12 + 2*ALambdax*Ye12)*(3*(Ye02*(Ye00*Ye10 + Ye01*
      Ye11 + Ye02*Ye12) + Ye22*(Ye10*Ye20 + Ye11*Ye21 + Ye12*Ye22) + Ye12*(Sqr(
      Ye10) + Sqr(Ye11) + Sqr(Ye12))) + Ye12*(-1.8*Sqr(g1) - 3*Sqr(g2) + Sqr(
      Lambdax) - 2*Sqr(gN)*Sqr(Qep) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QLp) + 3
      *(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) +
      Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr
      (Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22))) + (6*
      TYd00 + 6*ALambdax*Yd00)*(Yu00*(Yd00*Yu00 + Yd01*Yu01 + Yd02*Yu02) + Yu10*(
      Yd00*Yu10 + Yd01*Yu11 + Yd02*Yu12) + Yu20*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22
      ) + 3*(Yd10*(Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + Yd20*(Yd00*Yd20 + Yd01*
      Yd21 + Yd02*Yd22) + Yd00*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02))) + Yd00*(
      -0.4666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(
      Lambdax) - 2*Sqr(gN)*Sqr(Qdp) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QQp) + 3
      *(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) +
      Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr
      (Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22))) + (6*
      TYd01 + 6*ALambdax*Yd01)*(Yu01*(Yd00*Yu00 + Yd01*Yu01 + Yd02*Yu02) + Yu11*(
      Yd00*Yu10 + Yd01*Yu11 + Yd02*Yu12) + Yu21*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22
      ) + 3*(Yd11*(Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + Yd21*(Yd00*Yd20 + Yd01*
      Yd21 + Yd02*Yd22) + Yd01*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02))) + Yd01*(
      -0.4666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(
      Lambdax) - 2*Sqr(gN)*Sqr(Qdp) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QQp) + 3
      *(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) +
      Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr
      (Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22))) + (6*
      TYd02 + 6*ALambdax*Yd02)*(Yu02*(Yd00*Yu00 + Yd01*Yu01 + Yd02*Yu02) + Yu12*(
      Yd00*Yu10 + Yd01*Yu11 + Yd02*Yu12) + Yu22*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22
      ) + 3*(Yd12*(Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + Yd22*(Yd00*Yd20 + Yd01*
      Yd21 + Yd02*Yd22) + Yd02*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02))) + Yd02*(
      -0.4666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(
      Lambdax) - 2*Sqr(gN)*Sqr(Qdp) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QQp) + 3
      *(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) +
      Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr
      (Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22))) + (6*
      TYd10 + 6*ALambdax*Yd10)*(Yu00*(Yd10*Yu00 + Yd11*Yu01 + Yd12*Yu02) + Yu10*(
      Yd10*Yu10 + Yd11*Yu11 + Yd12*Yu12) + Yu20*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22
      ) + 3*(Yd00*(Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + Yd20*(Yd10*Yd20 + Yd11*
      Yd21 + Yd12*Yd22) + Yd10*(Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12))) + Yd10*(
      -0.4666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(
      Lambdax) - 2*Sqr(gN)*Sqr(Qdp) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QQp) + 3
      *(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) +
      Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr
      (Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22))) + (6*
      TYd11 + 6*ALambdax*Yd11)*(Yu01*(Yd10*Yu00 + Yd11*Yu01 + Yd12*Yu02) + Yu11*(
      Yd10*Yu10 + Yd11*Yu11 + Yd12*Yu12) + Yu21*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22
      ) + 3*(Yd01*(Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + Yd21*(Yd10*Yd20 + Yd11*
      Yd21 + Yd12*Yd22) + Yd11*(Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12))) + Yd11*(
      -0.4666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(
      Lambdax) - 2*Sqr(gN)*Sqr(Qdp) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QQp) + 3
      *(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) +
      Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr
      (Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22))) + (6*
      TYd12 + 6*ALambdax*Yd12)*(Yu02*(Yd10*Yu00 + Yd11*Yu01 + Yd12*Yu02) + Yu12*(
      Yd10*Yu10 + Yd11*Yu11 + Yd12*Yu12) + Yu22*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22
      ) + 3*(Yd02*(Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + Yd22*(Yd10*Yd20 + Yd11*
      Yd21 + Yd12*Yd22) + Yd12*(Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12))) + Yd12*(
      -0.4666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(
      Lambdax) - 2*Sqr(gN)*Sqr(Qdp) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QQp) + 3
      *(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) +
      Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr
      (Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22))) + (6*
      TYd20 + 6*ALambdax*Yd20)*(Yu00*(Yd20*Yu00 + Yd21*Yu01 + Yd22*Yu02) + Yu10*(
      Yd20*Yu10 + Yd21*Yu11 + Yd22*Yu12) + Yu20*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22
      ) + 3*(Yd00*(Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22) + Yd10*(Yd10*Yd20 + Yd11*
      Yd21 + Yd12*Yd22) + Yd20*(Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22))) + Yd20*(
      -0.4666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(
      Lambdax) - 2*Sqr(gN)*Sqr(Qdp) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QQp) + 3
      *(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) +
      Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr
      (Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22))) + (6*
      TYd21 + 6*ALambdax*Yd21)*(Yu01*(Yd20*Yu00 + Yd21*Yu01 + Yd22*Yu02) + Yu11*(
      Yd20*Yu10 + Yd21*Yu11 + Yd22*Yu12) + Yu21*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22
      ) + 3*(Yd01*(Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22) + Yd11*(Yd10*Yd20 + Yd11*
      Yd21 + Yd12*Yd22) + Yd21*(Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22))) + Yd21*(
      -0.4666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(
      Lambdax) - 2*Sqr(gN)*Sqr(Qdp) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QQp) + 3
      *(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) +
      Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr
      (Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22))) + (6*
      TYd22 + 6*ALambdax*Yd22)*(Yu02*(Yd20*Yu00 + Yd21*Yu01 + Yd22*Yu02) + Yu12*(
      Yd20*Yu10 + Yd21*Yu11 + Yd22*Yu12) + Yu22*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22
      ) + 3*(Yd02*(Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22) + Yd12*(Yd10*Yd20 + Yd11*
      Yd21 + Yd12*Yd22) + Yd22*(Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22))) + Yd22*(
      -0.4666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(
      Lambdax) - 2*Sqr(gN)*Sqr(Qdp) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QQp) + 3
      *(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) +
      Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr
      (Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22))) + 2*
      Ye20*(5*(Ye00*(TYe20*Ye00 + TYe21*Ye01 + TYe22*Ye02) + Ye10*(TYe20*Ye10 +
      TYe21*Ye11 + TYe22*Ye12) + Ye20*(TYe20*Ye20 + TYe21*Ye21 + TYe22*Ye22)) -
      1.8*TYe20*Sqr(g1) - 3*TYe20*Sqr(g2) + TYe20*Sqr(Lambdax) - 2*TYe20*Sqr(gN)*
      Sqr(Qep) - 2*TYe20*Sqr(gN)*Sqr(QH1p) - 2*TYe20*Sqr(gN)*Sqr(QLp) + Ye20*(6*(
      TYd00*Yd00 + TYd01*Yd01 + TYd02*Yd02 + TYd10*Yd10 + TYd11*Yd11 + TYd12*Yd12
      + TYd20*Yd20 + TYd21*Yd21 + TYd22*Yd22) + 2*(TYe00*Ye00 + TYe01*Ye01 + TYe02
      *Ye02 + TYe10*Ye10 + TYe11*Ye11 + TYe12*Ye12 + TYe20*Ye20 + TYe21*Ye21 +
      TYe22*Ye22) + 3.6*MassB*Sqr(g1) + 6*MassWB*Sqr(g2) + 2*ALambdax*Sqr(Lambdax)
      + 4*MassBp*Sqr(gN)*Sqr(Qep) + 4*MassBp*Sqr(gN)*Sqr(QH1p) + 4*MassBp*Sqr(gN)
      *Sqr(QLp)) + 3*TYe20*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(
      Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + TYe20*(Sqr(Ye00) +
      Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(
      Ye21) + Sqr(Ye22)) + 4*(TYe00*(Ye00*Ye20 + Ye01*Ye21 + Ye02*Ye22) + TYe10*(
      Ye10*Ye20 + Ye11*Ye21 + Ye12*Ye22) + TYe20*(Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22
      )))) + 2*Ye21*(5*(Ye01*(TYe20*Ye00 + TYe21*Ye01 + TYe22*Ye02) + Ye11*(TYe20*
      Ye10 + TYe21*Ye11 + TYe22*Ye12) + Ye21*(TYe20*Ye20 + TYe21*Ye21 + TYe22*Ye22
      )) - 1.8*TYe21*Sqr(g1) - 3*TYe21*Sqr(g2) + TYe21*Sqr(Lambdax) - 2*TYe21*Sqr(
      gN)*Sqr(Qep) - 2*TYe21*Sqr(gN)*Sqr(QH1p) - 2*TYe21*Sqr(gN)*Sqr(QLp) + Ye21*(
      6*(TYd00*Yd00 + TYd01*Yd01 + TYd02*Yd02 + TYd10*Yd10 + TYd11*Yd11 + TYd12*
      Yd12 + TYd20*Yd20 + TYd21*Yd21 + TYd22*Yd22) + 2*(TYe00*Ye00 + TYe01*Ye01 +
      TYe02*Ye02 + TYe10*Ye10 + TYe11*Ye11 + TYe12*Ye12 + TYe20*Ye20 + TYe21*Ye21
      + TYe22*Ye22) + 3.6*MassB*Sqr(g1) + 6*MassWB*Sqr(g2) + 2*ALambdax*Sqr(
      Lambdax) + 4*MassBp*Sqr(gN)*Sqr(Qep) + 4*MassBp*Sqr(gN)*Sqr(QH1p) + 4*MassBp
      *Sqr(gN)*Sqr(QLp)) + 3*TYe21*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10)
      + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + TYe21*(Sqr(
      Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20
      ) + Sqr(Ye21) + Sqr(Ye22)) + 4*(TYe01*(Ye00*Ye20 + Ye01*Ye21 + Ye02*Ye22) +
      TYe11*(Ye10*Ye20 + Ye11*Ye21 + Ye12*Ye22) + TYe21*(Sqr(Ye20) + Sqr(Ye21) +
      Sqr(Ye22)))) + 2*Ye22*(5*(Ye02*(TYe20*Ye00 + TYe21*Ye01 + TYe22*Ye02) + Ye12
      *(TYe20*Ye10 + TYe21*Ye11 + TYe22*Ye12) + Ye22*(TYe20*Ye20 + TYe21*Ye21 +
      TYe22*Ye22)) - 1.8*TYe22*Sqr(g1) - 3*TYe22*Sqr(g2) + TYe22*Sqr(Lambdax) - 2*
      TYe22*Sqr(gN)*Sqr(Qep) - 2*TYe22*Sqr(gN)*Sqr(QH1p) - 2*TYe22*Sqr(gN)*Sqr(QLp
      ) + Ye22*(6*(TYd00*Yd00 + TYd01*Yd01 + TYd02*Yd02 + TYd10*Yd10 + TYd11*Yd11
      + TYd12*Yd12 + TYd20*Yd20 + TYd21*Yd21 + TYd22*Yd22) + 2*(TYe00*Ye00 + TYe01
      *Ye01 + TYe02*Ye02 + TYe10*Ye10 + TYe11*Ye11 + TYe12*Ye12 + TYe20*Ye20 +
      TYe21*Ye21 + TYe22*Ye22) + 3.6*MassB*Sqr(g1) + 6*MassWB*Sqr(g2) + 2*ALambdax
      *Sqr(Lambdax) + 4*MassBp*Sqr(gN)*Sqr(Qep) + 4*MassBp*Sqr(gN)*Sqr(QH1p) + 4*
      MassBp*Sqr(gN)*Sqr(QLp)) + 3*TYe22*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(
      Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + TYe22*(
      Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(
      Ye20) + Sqr(Ye21) + Sqr(Ye22)) + 4*(TYe02*(Ye00*Ye20 + Ye01*Ye21 + Ye02*Ye22
      ) + TYe12*(Ye10*Ye20 + Ye11*Ye21 + Ye12*Ye22) + TYe22*(Sqr(Ye20) + Sqr(Ye21)
      + Sqr(Ye22)))) + (2*TYe20 + 2*ALambdax*Ye20)*(Ye20*(-1.8*Sqr(g1) - 3*Sqr(g2
      ) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(Qep) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(
      QLp) + 3*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(
      Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr(
      Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22
      )) + 3*(Ye00*(Ye00*Ye20 + Ye01*Ye21 + Ye02*Ye22) + Ye10*(Ye10*Ye20 + Ye11*
      Ye21 + Ye12*Ye22) + Ye20*(Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22)))) + (2*TYe21 +
      2*ALambdax*Ye21)*(Ye21*(-1.8*Sqr(g1) - 3*Sqr(g2) + Sqr(Lambdax) - 2*Sqr(gN)*
      Sqr(Qep) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QLp) + 3*(Sqr(Yd00) + Sqr(
      Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21
      ) + Sqr(Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) +
      Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22)) + 3*(Ye01*(Ye00*Ye20 + Ye01*
      Ye21 + Ye02*Ye22) + Ye11*(Ye10*Ye20 + Ye11*Ye21 + Ye12*Ye22) + Ye21*(Sqr(
      Ye20) + Sqr(Ye21) + Sqr(Ye22)))) + (2*TYe22 + 2*ALambdax*Ye22)*(Ye22*(-1.8*
      Sqr(g1) - 3*Sqr(g2) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(Qep) - 2*Sqr(gN)*Sqr(QH1p
      ) - 2*Sqr(gN)*Sqr(QLp) + 3*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) +
      Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(Ye00) + Sqr
      (Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(
      Ye21) + Sqr(Ye22)) + 3*(Ye02*(Ye00*Ye20 + Ye01*Ye21 + Ye02*Ye22) + Ye12*(
      Ye10*Ye20 + Ye11*Ye21 + Ye12*Ye22) + Ye22*(Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22)
      ))) + (6*(Kappa00*TKappa00 + Kappa01*TKappa01 + Kappa02*TKappa02 + Kappa10*
      TKappa10 + Kappa11*TKappa11 + Kappa12*TKappa12 + Kappa20*TKappa20 + Kappa21*
      TKappa21 + Kappa22*TKappa22) + 4*(Lambda1200*TLambda1200 + Lambda1201*
      TLambda1201 + Lambda1210*TLambda1210 + Lambda1211*TLambda1211) + 6*(TYd00*
      Yd00 + TYd01*Yd01 + TYd02*Yd02 + TYd10*Yd10 + TYd11*Yd11 + TYd12*Yd12 +
      TYd20*Yd20 + TYd21*Yd21 + TYd22*Yd22) + 2*(TYe00*Ye00 + TYe01*Ye01 + TYe02*
      Ye02 + TYe10*Ye10 + TYe11*Ye11 + TYe12*Ye12 + TYe20*Ye20 + TYe21*Ye21 +
      TYe22*Ye22) + 6*(TYu00*Yu00 + TYu01*Yu01 + TYu02*Yu02 + TYu10*Yu10 + TYu11*
      Yu11 + TYu12*Yu12 + TYu20*Yu20 + TYu21*Yu21 + TYu22*Yu22) + 1.2*MassB*Sqr(g1
      ) + 6*MassWB*Sqr(g2) + 24*ALambdax*Sqr(Lambdax) + 4*MassBp*Sqr(gN)*Sqr(QH1p)
      + 4*MassBp*Sqr(gN)*Sqr(QH2p) + 4*MassBp*Sqr(gN)*Sqr(QSp))*(-0.6*Sqr(g1) - 3
      *Sqr(g2) + 3*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) +
      Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) +
      2*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) +
      12*Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(
      QSp) + 3*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(
      Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr(
      Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22
      ) + 3*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12)
      + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + 48*ALambdax*Lambdax*(4*Power(
      Lambdax,3) - 0.6*Lambdax*Sqr(g1) - 3*Lambdax*Sqr(g2) + 3*Lambdax*(Sqr(
      Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(
      Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*Lambdax*(Sqr(
      Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) - 2*
      Lambdax*Sqr(gN)*Sqr(QH1p) - 2*Lambdax*Sqr(gN)*Sqr(QH2p) - 2*Lambdax*Sqr(gN)*
      Sqr(QSp) + 3*Lambdax*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(
      Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Lambdax*(Sqr(Ye00)
      + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) +
      Sqr(Ye21) + Sqr(Ye22)) + 3*Lambdax*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(
      Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + 6*Yu00
      *(Yd00*(TYu00*Yd00 + TYu01*Yd01 + TYu02*Yd02) + Yd10*(TYu00*Yd10 + TYu01*
      Yd11 + TYu02*Yd12) + Yd20*(TYu00*Yd20 + TYu01*Yd21 + TYu02*Yd22) + 2*(TYd00*
      (Yd00*Yu00 + Yd01*Yu01 + Yd02*Yu02) + TYd10*(Yd10*Yu00 + Yd11*Yu01 + Yd12*
      Yu02) + TYd20*(Yd20*Yu00 + Yd21*Yu01 + Yd22*Yu02)) + 5*(Yu00*(TYu00*Yu00 +
      TYu01*Yu01 + TYu02*Yu02) + Yu10*(TYu00*Yu10 + TYu01*Yu11 + TYu02*Yu12) +
      Yu20*(TYu00*Yu20 + TYu01*Yu21 + TYu02*Yu22)) - 0.8666666666666667*TYu00*Sqr(
      g1) - 3*TYu00*Sqr(g2) - 5.333333333333333*TYu00*Sqr(g3) + TYu00*Sqr(Lambdax)
      - 2*TYu00*Sqr(gN)*Sqr(QH2p) - 2*TYu00*Sqr(gN)*Sqr(QQp) - 2*TYu00*Sqr(gN)*
      Sqr(Qup) + Yu00*(6*(TYu00*Yu00 + TYu01*Yu01 + TYu02*Yu02 + TYu10*Yu10 +
      TYu11*Yu11 + TYu12*Yu12 + TYu20*Yu20 + TYu21*Yu21 + TYu22*Yu22) +
      1.7333333333333334*MassB*Sqr(g1) + 6*MassWB*Sqr(g2) + 10.666666666666666*
      MassG*Sqr(g3) + 2*ALambdax*Sqr(Lambdax) + 4*MassBp*Sqr(gN)*Sqr(QH2p) + 4*
      MassBp*Sqr(gN)*Sqr(QQp) + 4*MassBp*Sqr(gN)*Sqr(Qup)) + 4*(TYu10*(Yu00*Yu10 +
      Yu01*Yu11 + Yu02*Yu12) + TYu20*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + TYu00*
      (Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02))) + 3*TYu00*(Sqr(Yu00) + Sqr(Yu01) + Sqr(
      Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22
      ))) + 6*Yu01*(Yd01*(TYu00*Yd00 + TYu01*Yd01 + TYu02*Yd02) + Yd11*(TYu00*Yd10
      + TYu01*Yd11 + TYu02*Yd12) + Yd21*(TYu00*Yd20 + TYu01*Yd21 + TYu02*Yd22) +
      2*(TYd01*(Yd00*Yu00 + Yd01*Yu01 + Yd02*Yu02) + TYd11*(Yd10*Yu00 + Yd11*Yu01
      + Yd12*Yu02) + TYd21*(Yd20*Yu00 + Yd21*Yu01 + Yd22*Yu02)) + 5*(Yu01*(TYu00*
      Yu00 + TYu01*Yu01 + TYu02*Yu02) + Yu11*(TYu00*Yu10 + TYu01*Yu11 + TYu02*Yu12
      ) + Yu21*(TYu00*Yu20 + TYu01*Yu21 + TYu02*Yu22)) - 0.8666666666666667*TYu01*
      Sqr(g1) - 3*TYu01*Sqr(g2) - 5.333333333333333*TYu01*Sqr(g3) + TYu01*Sqr(
      Lambdax) - 2*TYu01*Sqr(gN)*Sqr(QH2p) - 2*TYu01*Sqr(gN)*Sqr(QQp) - 2*TYu01*
      Sqr(gN)*Sqr(Qup) + Yu01*(6*(TYu00*Yu00 + TYu01*Yu01 + TYu02*Yu02 + TYu10*
      Yu10 + TYu11*Yu11 + TYu12*Yu12 + TYu20*Yu20 + TYu21*Yu21 + TYu22*Yu22) +
      1.7333333333333334*MassB*Sqr(g1) + 6*MassWB*Sqr(g2) + 10.666666666666666*
      MassG*Sqr(g3) + 2*ALambdax*Sqr(Lambdax) + 4*MassBp*Sqr(gN)*Sqr(QH2p) + 4*
      MassBp*Sqr(gN)*Sqr(QQp) + 4*MassBp*Sqr(gN)*Sqr(Qup)) + 4*(TYu11*(Yu00*Yu10 +
      Yu01*Yu11 + Yu02*Yu12) + TYu21*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + TYu01*
      (Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02))) + 3*TYu01*(Sqr(Yu00) + Sqr(Yu01) + Sqr(
      Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22
      ))) + 6*Yu02*(Yd02*(TYu00*Yd00 + TYu01*Yd01 + TYu02*Yd02) + Yd12*(TYu00*Yd10
      + TYu01*Yd11 + TYu02*Yd12) + Yd22*(TYu00*Yd20 + TYu01*Yd21 + TYu02*Yd22) +
      2*(TYd02*(Yd00*Yu00 + Yd01*Yu01 + Yd02*Yu02) + TYd12*(Yd10*Yu00 + Yd11*Yu01
      + Yd12*Yu02) + TYd22*(Yd20*Yu00 + Yd21*Yu01 + Yd22*Yu02)) + 5*(Yu02*(TYu00*
      Yu00 + TYu01*Yu01 + TYu02*Yu02) + Yu12*(TYu00*Yu10 + TYu01*Yu11 + TYu02*Yu12
      ) + Yu22*(TYu00*Yu20 + TYu01*Yu21 + TYu02*Yu22)) - 0.8666666666666667*TYu02*
      Sqr(g1) - 3*TYu02*Sqr(g2) - 5.333333333333333*TYu02*Sqr(g3) + TYu02*Sqr(
      Lambdax) - 2*TYu02*Sqr(gN)*Sqr(QH2p) - 2*TYu02*Sqr(gN)*Sqr(QQp) - 2*TYu02*
      Sqr(gN)*Sqr(Qup) + Yu02*(6*(TYu00*Yu00 + TYu01*Yu01 + TYu02*Yu02 + TYu10*
      Yu10 + TYu11*Yu11 + TYu12*Yu12 + TYu20*Yu20 + TYu21*Yu21 + TYu22*Yu22) +
      1.7333333333333334*MassB*Sqr(g1) + 6*MassWB*Sqr(g2) + 10.666666666666666*
      MassG*Sqr(g3) + 2*ALambdax*Sqr(Lambdax) + 4*MassBp*Sqr(gN)*Sqr(QH2p) + 4*
      MassBp*Sqr(gN)*Sqr(QQp) + 4*MassBp*Sqr(gN)*Sqr(Qup)) + 4*(TYu12*(Yu00*Yu10 +
      Yu01*Yu11 + Yu02*Yu12) + TYu22*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + TYu02*
      (Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02))) + 3*TYu02*(Sqr(Yu00) + Sqr(Yu01) + Sqr(
      Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22
      ))) + 6*Yu10*(Yd00*(TYu10*Yd00 + TYu11*Yd01 + TYu12*Yd02) + Yd10*(TYu10*Yd10
      + TYu11*Yd11 + TYu12*Yd12) + Yd20*(TYu10*Yd20 + TYu11*Yd21 + TYu12*Yd22) +
      2*(TYd00*(Yd00*Yu10 + Yd01*Yu11 + Yd02*Yu12) + TYd10*(Yd10*Yu10 + Yd11*Yu11
      + Yd12*Yu12) + TYd20*(Yd20*Yu10 + Yd21*Yu11 + Yd22*Yu12)) + 5*(Yu00*(TYu10*
      Yu00 + TYu11*Yu01 + TYu12*Yu02) + Yu10*(TYu10*Yu10 + TYu11*Yu11 + TYu12*Yu12
      ) + Yu20*(TYu10*Yu20 + TYu11*Yu21 + TYu12*Yu22)) - 0.8666666666666667*TYu10*
      Sqr(g1) - 3*TYu10*Sqr(g2) - 5.333333333333333*TYu10*Sqr(g3) + TYu10*Sqr(
      Lambdax) - 2*TYu10*Sqr(gN)*Sqr(QH2p) - 2*TYu10*Sqr(gN)*Sqr(QQp) - 2*TYu10*
      Sqr(gN)*Sqr(Qup) + Yu10*(6*(TYu00*Yu00 + TYu01*Yu01 + TYu02*Yu02 + TYu10*
      Yu10 + TYu11*Yu11 + TYu12*Yu12 + TYu20*Yu20 + TYu21*Yu21 + TYu22*Yu22) +
      1.7333333333333334*MassB*Sqr(g1) + 6*MassWB*Sqr(g2) + 10.666666666666666*
      MassG*Sqr(g3) + 2*ALambdax*Sqr(Lambdax) + 4*MassBp*Sqr(gN)*Sqr(QH2p) + 4*
      MassBp*Sqr(gN)*Sqr(QQp) + 4*MassBp*Sqr(gN)*Sqr(Qup)) + 4*(TYu00*(Yu00*Yu10 +
      Yu01*Yu11 + Yu02*Yu12) + TYu20*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + TYu10*
      (Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12))) + 3*TYu10*(Sqr(Yu00) + Sqr(Yu01) + Sqr(
      Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22
      ))) + 6*Yu11*(Yd01*(TYu10*Yd00 + TYu11*Yd01 + TYu12*Yd02) + Yd11*(TYu10*Yd10
      + TYu11*Yd11 + TYu12*Yd12) + Yd21*(TYu10*Yd20 + TYu11*Yd21 + TYu12*Yd22) +
      2*(TYd01*(Yd00*Yu10 + Yd01*Yu11 + Yd02*Yu12) + TYd11*(Yd10*Yu10 + Yd11*Yu11
      + Yd12*Yu12) + TYd21*(Yd20*Yu10 + Yd21*Yu11 + Yd22*Yu12)) + 5*(Yu01*(TYu10*
      Yu00 + TYu11*Yu01 + TYu12*Yu02) + Yu11*(TYu10*Yu10 + TYu11*Yu11 + TYu12*Yu12
      ) + Yu21*(TYu10*Yu20 + TYu11*Yu21 + TYu12*Yu22)) - 0.8666666666666667*TYu11*
      Sqr(g1) - 3*TYu11*Sqr(g2) - 5.333333333333333*TYu11*Sqr(g3) + TYu11*Sqr(
      Lambdax) - 2*TYu11*Sqr(gN)*Sqr(QH2p) - 2*TYu11*Sqr(gN)*Sqr(QQp) - 2*TYu11*
      Sqr(gN)*Sqr(Qup) + Yu11*(6*(TYu00*Yu00 + TYu01*Yu01 + TYu02*Yu02 + TYu10*
      Yu10 + TYu11*Yu11 + TYu12*Yu12 + TYu20*Yu20 + TYu21*Yu21 + TYu22*Yu22) +
      1.7333333333333334*MassB*Sqr(g1) + 6*MassWB*Sqr(g2) + 10.666666666666666*
      MassG*Sqr(g3) + 2*ALambdax*Sqr(Lambdax) + 4*MassBp*Sqr(gN)*Sqr(QH2p) + 4*
      MassBp*Sqr(gN)*Sqr(QQp) + 4*MassBp*Sqr(gN)*Sqr(Qup)) + 4*(TYu01*(Yu00*Yu10 +
      Yu01*Yu11 + Yu02*Yu12) + TYu21*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + TYu11*
      (Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12))) + 3*TYu11*(Sqr(Yu00) + Sqr(Yu01) + Sqr(
      Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22
      ))) + 6*Yu12*(Yd02*(TYu10*Yd00 + TYu11*Yd01 + TYu12*Yd02) + Yd12*(TYu10*Yd10
      + TYu11*Yd11 + TYu12*Yd12) + Yd22*(TYu10*Yd20 + TYu11*Yd21 + TYu12*Yd22) +
      2*(TYd02*(Yd00*Yu10 + Yd01*Yu11 + Yd02*Yu12) + TYd12*(Yd10*Yu10 + Yd11*Yu11
      + Yd12*Yu12) + TYd22*(Yd20*Yu10 + Yd21*Yu11 + Yd22*Yu12)) + 5*(Yu02*(TYu10*
      Yu00 + TYu11*Yu01 + TYu12*Yu02) + Yu12*(TYu10*Yu10 + TYu11*Yu11 + TYu12*Yu12
      ) + Yu22*(TYu10*Yu20 + TYu11*Yu21 + TYu12*Yu22)) - 0.8666666666666667*TYu12*
      Sqr(g1) - 3*TYu12*Sqr(g2) - 5.333333333333333*TYu12*Sqr(g3) + TYu12*Sqr(
      Lambdax) - 2*TYu12*Sqr(gN)*Sqr(QH2p) - 2*TYu12*Sqr(gN)*Sqr(QQp) - 2*TYu12*
      Sqr(gN)*Sqr(Qup) + Yu12*(6*(TYu00*Yu00 + TYu01*Yu01 + TYu02*Yu02 + TYu10*
      Yu10 + TYu11*Yu11 + TYu12*Yu12 + TYu20*Yu20 + TYu21*Yu21 + TYu22*Yu22) +
      1.7333333333333334*MassB*Sqr(g1) + 6*MassWB*Sqr(g2) + 10.666666666666666*
      MassG*Sqr(g3) + 2*ALambdax*Sqr(Lambdax) + 4*MassBp*Sqr(gN)*Sqr(QH2p) + 4*
      MassBp*Sqr(gN)*Sqr(QQp) + 4*MassBp*Sqr(gN)*Sqr(Qup)) + 4*(TYu02*(Yu00*Yu10 +
      Yu01*Yu11 + Yu02*Yu12) + TYu22*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + TYu12*
      (Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12))) + 3*TYu12*(Sqr(Yu00) + Sqr(Yu01) + Sqr(
      Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22
      ))) + 6*Yu20*(Yd00*(TYu20*Yd00 + TYu21*Yd01 + TYu22*Yd02) + Yd10*(TYu20*Yd10
      + TYu21*Yd11 + TYu22*Yd12) + Yd20*(TYu20*Yd20 + TYu21*Yd21 + TYu22*Yd22) +
      5*(Yu00*(TYu20*Yu00 + TYu21*Yu01 + TYu22*Yu02) + Yu10*(TYu20*Yu10 + TYu21*
      Yu11 + TYu22*Yu12) + Yu20*(TYu20*Yu20 + TYu21*Yu21 + TYu22*Yu22)) + 2*(TYd00
      *(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22) + TYd10*(Yd10*Yu20 + Yd11*Yu21 + Yd12*
      Yu22) + TYd20*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22)) - 0.8666666666666667*
      TYu20*Sqr(g1) - 3*TYu20*Sqr(g2) - 5.333333333333333*TYu20*Sqr(g3) + TYu20*
      Sqr(Lambdax) - 2*TYu20*Sqr(gN)*Sqr(QH2p) - 2*TYu20*Sqr(gN)*Sqr(QQp) - 2*
      TYu20*Sqr(gN)*Sqr(Qup) + Yu20*(6*(TYu00*Yu00 + TYu01*Yu01 + TYu02*Yu02 +
      TYu10*Yu10 + TYu11*Yu11 + TYu12*Yu12 + TYu20*Yu20 + TYu21*Yu21 + TYu22*Yu22)
      + 1.7333333333333334*MassB*Sqr(g1) + 6*MassWB*Sqr(g2) + 10.666666666666666*
      MassG*Sqr(g3) + 2*ALambdax*Sqr(Lambdax) + 4*MassBp*Sqr(gN)*Sqr(QH2p) + 4*
      MassBp*Sqr(gN)*Sqr(QQp) + 4*MassBp*Sqr(gN)*Sqr(Qup)) + 3*TYu20*(Sqr(Yu00) +
      Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(
      Yu21) + Sqr(Yu22)) + 4*(TYu00*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + TYu10*(
      Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + TYu20*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22
      )))) + 6*Yu21*(Yd01*(TYu20*Yd00 + TYu21*Yd01 + TYu22*Yd02) + Yd11*(TYu20*
      Yd10 + TYu21*Yd11 + TYu22*Yd12) + Yd21*(TYu20*Yd20 + TYu21*Yd21 + TYu22*Yd22
      ) + 5*(Yu01*(TYu20*Yu00 + TYu21*Yu01 + TYu22*Yu02) + Yu11*(TYu20*Yu10 +
      TYu21*Yu11 + TYu22*Yu12) + Yu21*(TYu20*Yu20 + TYu21*Yu21 + TYu22*Yu22)) + 2*
      (TYd01*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22) + TYd11*(Yd10*Yu20 + Yd11*Yu21 +
      Yd12*Yu22) + TYd21*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22)) - 0.8666666666666667
      *TYu21*Sqr(g1) - 3*TYu21*Sqr(g2) - 5.333333333333333*TYu21*Sqr(g3) + TYu21*
      Sqr(Lambdax) - 2*TYu21*Sqr(gN)*Sqr(QH2p) - 2*TYu21*Sqr(gN)*Sqr(QQp) - 2*
      TYu21*Sqr(gN)*Sqr(Qup) + Yu21*(6*(TYu00*Yu00 + TYu01*Yu01 + TYu02*Yu02 +
      TYu10*Yu10 + TYu11*Yu11 + TYu12*Yu12 + TYu20*Yu20 + TYu21*Yu21 + TYu22*Yu22)
      + 1.7333333333333334*MassB*Sqr(g1) + 6*MassWB*Sqr(g2) + 10.666666666666666*
      MassG*Sqr(g3) + 2*ALambdax*Sqr(Lambdax) + 4*MassBp*Sqr(gN)*Sqr(QH2p) + 4*
      MassBp*Sqr(gN)*Sqr(QQp) + 4*MassBp*Sqr(gN)*Sqr(Qup)) + 3*TYu21*(Sqr(Yu00) +
      Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(
      Yu21) + Sqr(Yu22)) + 4*(TYu01*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + TYu11*(
      Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + TYu21*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22
      )))) + 6*Yu22*(Yd02*(TYu20*Yd00 + TYu21*Yd01 + TYu22*Yd02) + Yd12*(TYu20*
      Yd10 + TYu21*Yd11 + TYu22*Yd12) + Yd22*(TYu20*Yd20 + TYu21*Yd21 + TYu22*Yd22
      ) + 5*(Yu02*(TYu20*Yu00 + TYu21*Yu01 + TYu22*Yu02) + Yu12*(TYu20*Yu10 +
      TYu21*Yu11 + TYu22*Yu12) + Yu22*(TYu20*Yu20 + TYu21*Yu21 + TYu22*Yu22)) + 2*
      (TYd02*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22) + TYd12*(Yd10*Yu20 + Yd11*Yu21 +
      Yd12*Yu22) + TYd22*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22)) - 0.8666666666666667
      *TYu22*Sqr(g1) - 3*TYu22*Sqr(g2) - 5.333333333333333*TYu22*Sqr(g3) + TYu22*
      Sqr(Lambdax) - 2*TYu22*Sqr(gN)*Sqr(QH2p) - 2*TYu22*Sqr(gN)*Sqr(QQp) - 2*
      TYu22*Sqr(gN)*Sqr(Qup) + Yu22*(6*(TYu00*Yu00 + TYu01*Yu01 + TYu02*Yu02 +
      TYu10*Yu10 + TYu11*Yu11 + TYu12*Yu12 + TYu20*Yu20 + TYu21*Yu21 + TYu22*Yu22)
      + 1.7333333333333334*MassB*Sqr(g1) + 6*MassWB*Sqr(g2) + 10.666666666666666*
      MassG*Sqr(g3) + 2*ALambdax*Sqr(Lambdax) + 4*MassBp*Sqr(gN)*Sqr(QH2p) + 4*
      MassBp*Sqr(gN)*Sqr(QQp) + 4*MassBp*Sqr(gN)*Sqr(Qup)) + 3*TYu22*(Sqr(Yu00) +
      Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(
      Yu21) + Sqr(Yu22)) + 4*(TYu02*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + TYu12*(
      Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + TYu22*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22
      )))) + (6*TYu00 + 6*ALambdax*Yu00)*(Yd00*(Yd00*Yu00 + Yd01*Yu01 + Yd02*Yu02)
      + Yd10*(Yd10*Yu00 + Yd11*Yu01 + Yd12*Yu02) + Yd20*(Yd20*Yu00 + Yd21*Yu01 +
      Yd22*Yu02) + 3*(Yu10*(Yu00*Yu10 + Yu01*Yu11 + Yu02*Yu12) + Yu20*(Yu00*Yu20 +
      Yu01*Yu21 + Yu02*Yu22) + Yu00*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02))) + Yu00*(
      -0.8666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(
      Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3
      *(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) +
      Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) + (6*TYu01 + 6*ALambdax*Yu01)*(Yd01*(
      Yd00*Yu00 + Yd01*Yu01 + Yd02*Yu02) + Yd11*(Yd10*Yu00 + Yd11*Yu01 + Yd12*Yu02
      ) + Yd21*(Yd20*Yu00 + Yd21*Yu01 + Yd22*Yu02) + 3*(Yu11*(Yu00*Yu10 + Yu01*
      Yu11 + Yu02*Yu12) + Yu21*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu01*(Sqr(
      Yu00) + Sqr(Yu01) + Sqr(Yu02))) + Yu01*(-0.8666666666666667*Sqr(g1) - 3*Sqr(
      g2) - 5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr
      (gN)*Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) +
      Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) + (
      6*TYu02 + 6*ALambdax*Yu02)*(Yd02*(Yd00*Yu00 + Yd01*Yu01 + Yd02*Yu02) + Yd12*
      (Yd10*Yu00 + Yd11*Yu01 + Yd12*Yu02) + Yd22*(Yd20*Yu00 + Yd21*Yu01 + Yd22*
      Yu02) + 3*(Yu12*(Yu00*Yu10 + Yu01*Yu11 + Yu02*Yu12) + Yu22*(Yu00*Yu20 + Yu01
      *Yu21 + Yu02*Yu22) + Yu02*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02))) + Yu02*(
      -0.8666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(
      Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3
      *(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) +
      Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) + (6*TYu10 + 6*ALambdax*Yu10)*(Yd00*(
      Yd00*Yu10 + Yd01*Yu11 + Yd02*Yu12) + Yd10*(Yd10*Yu10 + Yd11*Yu11 + Yd12*Yu12
      ) + Yd20*(Yd20*Yu10 + Yd21*Yu11 + Yd22*Yu12) + 3*(Yu00*(Yu00*Yu10 + Yu01*
      Yu11 + Yu02*Yu12) + Yu20*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + Yu10*(Sqr(
      Yu10) + Sqr(Yu11) + Sqr(Yu12))) + Yu10*(-0.8666666666666667*Sqr(g1) - 3*Sqr(
      g2) - 5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr
      (gN)*Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) +
      Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) + (
      6*TYu11 + 6*ALambdax*Yu11)*(Yd01*(Yd00*Yu10 + Yd01*Yu11 + Yd02*Yu12) + Yd11*
      (Yd10*Yu10 + Yd11*Yu11 + Yd12*Yu12) + Yd21*(Yd20*Yu10 + Yd21*Yu11 + Yd22*
      Yu12) + 3*(Yu01*(Yu00*Yu10 + Yu01*Yu11 + Yu02*Yu12) + Yu21*(Yu10*Yu20 + Yu11
      *Yu21 + Yu12*Yu22) + Yu11*(Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12))) + Yu11*(
      -0.8666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(
      Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3
      *(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) +
      Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) + (6*TYu12 + 6*ALambdax*Yu12)*(Yd02*(
      Yd00*Yu10 + Yd01*Yu11 + Yd02*Yu12) + Yd12*(Yd10*Yu10 + Yd11*Yu11 + Yd12*Yu12
      ) + Yd22*(Yd20*Yu10 + Yd21*Yu11 + Yd22*Yu12) + 3*(Yu02*(Yu00*Yu10 + Yu01*
      Yu11 + Yu02*Yu12) + Yu22*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + Yu12*(Sqr(
      Yu10) + Sqr(Yu11) + Sqr(Yu12))) + Yu12*(-0.8666666666666667*Sqr(g1) - 3*Sqr(
      g2) - 5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr
      (gN)*Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) +
      Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) + (
      6*TYu20 + 6*ALambdax*Yu20)*(Yd00*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22) + Yd10*
      (Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22) + Yd20*(Yd20*Yu20 + Yd21*Yu21 + Yd22*
      Yu22) + 3*(Yu00*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu10*(Yu10*Yu20 + Yu11
      *Yu21 + Yu12*Yu22) + Yu20*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + Yu20*(
      -0.8666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(
      Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3
      *(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) +
      Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) + (6*TYu21 + 6*ALambdax*Yu21)*(Yd01*(
      Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22) + Yd11*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22
      ) + Yd21*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22) + 3*(Yu01*(Yu00*Yu20 + Yu01*
      Yu21 + Yu02*Yu22) + Yu11*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + Yu21*(Sqr(
      Yu20) + Sqr(Yu21) + Sqr(Yu22))) + Yu21*(-0.8666666666666667*Sqr(g1) - 3*Sqr(
      g2) - 5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr
      (gN)*Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) +
      Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) + (
      6*TYu22 + 6*ALambdax*Yu22)*(Yd02*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22) + Yd12*
      (Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22) + Yd22*(Yd20*Yu20 + Yd21*Yu21 + Yd22*
      Yu22) + 3*(Yu02*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu12*(Yu10*Yu20 + Yu11
      *Yu21 + Yu12*Yu22) + Yu22*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + Yu22*(
      -0.8666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(
      Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3
      *(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) +
      Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) + (-0.6*Sqr(g1) - 3*Sqr(g2) + 3*(Sqr(
      Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(
      Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*(Sqr(Lambda1200)
      + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 12*Sqr(Lambdax) - 2
      *Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QSp) + 3*(Sqr(Yd00)
      + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) +
      Sqr(Yd21) + Sqr(Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr
      (Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22) + 3*(Sqr(Yu00) + Sqr(
      Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21
      ) + Sqr(Yu22)))*(6*(Kappa00*TKappa00 + Kappa01*TKappa01 + Kappa02*TKappa02 +
      Kappa10*TKappa10 + Kappa11*TKappa11 + Kappa12*TKappa12 + Kappa20*TKappa20 +
      Kappa21*TKappa21 + Kappa22*TKappa22) + 4*(Lambda1200*TLambda1200 +
      Lambda1201*TLambda1201 + Lambda1210*TLambda1210 + Lambda1211*TLambda1211) +
      6*(TYd00*Yd00 + TYd01*Yd01 + TYd02*Yd02 + TYd10*Yd10 + TYd11*Yd11 + TYd12*
      Yd12 + TYd20*Yd20 + TYd21*Yd21 + TYd22*Yd22) + 2*(TYe00*Ye00 + TYe01*Ye01 +
      TYe02*Ye02 + TYe10*Ye10 + TYe11*Ye11 + TYe12*Ye12 + TYe20*Ye20 + TYe21*Ye21
      + TYe22*Ye22) + 6*(TYu00*Yu00 + TYu01*Yu01 + TYu02*Yu02 + TYu10*Yu10 + TYu11
      *Yu11 + TYu12*Yu12 + TYu20*Yu20 + TYu21*Yu21 + TYu22*Yu22) + 1.2*MassB*Sqr(
      g1) + 6*MassWB*Sqr(g2) + 24*ALambdax*Sqr(Lambdax) + 4*MassBp*Sqr(gN)*Sqr(
      QH1p) + 4*MassBp*Sqr(gN)*Sqr(QH2p) + 4*MassBp*Sqr(gN)*Sqr(QSp) + ALambdax*(
      -0.6*Sqr(g1) - 3*Sqr(g2) + 3*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) +
      Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) +
      Sqr(Kappa22)) + 2*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr
      (Lambda1211)) + 12*Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QH2p)
      - 2*Sqr(gN)*Sqr(QSp) + 3*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) +
      Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(Ye00) + Sqr
      (Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(
      Ye21) + Sqr(Ye22) + 3*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(
      Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) + 24*Lambdax*(6*
      Lambdax*(Kappa00*TKappa00 + Kappa01*TKappa01 + Kappa02*TKappa02 + Kappa10*
      TKappa10 + Kappa11*TKappa11 + Kappa12*TKappa12 + Kappa20*TKappa20 + Kappa21*
      TKappa21 + Kappa22*TKappa22) + 4*Lambdax*(Lambda1200*TLambda1200 +
      Lambda1201*TLambda1201 + Lambda1210*TLambda1210 + Lambda1211*TLambda1211) +
      6*Lambdax*(TYd00*Yd00 + TYd01*Yd01 + TYd02*Yd02 + TYd10*Yd10 + TYd11*Yd11 +
      TYd12*Yd12 + TYd20*Yd20 + TYd21*Yd21 + TYd22*Yd22) + 2*Lambdax*(TYe00*Ye00 +
      TYe01*Ye01 + TYe02*Ye02 + TYe10*Ye10 + TYe11*Ye11 + TYe12*Ye12 + TYe20*Ye20
      + TYe21*Ye21 + TYe22*Ye22) + 6*Lambdax*(TYu00*Yu00 + TYu01*Yu01 + TYu02*
      Yu02 + TYu10*Yu10 + TYu11*Yu11 + TYu12*Yu12 + TYu20*Yu20 + TYu21*Yu21 +
      TYu22*Yu22) + 1.2*Lambdax*MassB*Sqr(g1) + 6*Lambdax*MassWB*Sqr(g2) + 4*
      Lambdax*MassBp*Sqr(gN)*Sqr(QH1p) + 4*Lambdax*MassBp*Sqr(gN)*Sqr(QH2p) + 4*
      Lambdax*MassBp*Sqr(gN)*Sqr(QSp) + ALambdax*Lambdax*(-0.6*Sqr(g1) - 3*Sqr(g2)
      + 3*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(
      Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*(
      Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 12*
      Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QSp
      ) + 3*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12)
      + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) +
      Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22) + 3*(
      Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(
      Yu20) + Sqr(Yu21) + Sqr(Yu22))));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dlead_log_one_loop_TLambdax_dALambdax() const 
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

   double deriv = -11.52*Power(g1,4)*Lambdax - 24*Power(g2,4)*Lambdax + 24*
      Power(Lambdax,3)*Sqr(Kappa00) + 24*Power(Lambdax,3)*Sqr(Kappa01) + 24*Power(
      Lambdax,3)*Sqr(Kappa02) + 24*Power(Lambdax,3)*Sqr(Kappa10) + 24*Power(
      Lambdax,3)*Sqr(Kappa11) + 24*Power(Lambdax,3)*Sqr(Kappa12) + 24*Power(
      Lambdax,3)*Sqr(Kappa20) + 24*Power(Lambdax,3)*Sqr(Kappa21) + 24*Power(
      Lambdax,3)*Sqr(Kappa22) + 16*Power(Lambdax,3)*Sqr(Lambda1200) + 16*Power(
      Lambdax,3)*Sqr(Lambda1201) + 16*Power(Lambdax,3)*Sqr(Lambda1210) + 16*Power(
      Lambdax,3)*Sqr(Lambda1211) + 6*Kappa00*Lambdax*(2*(Kappa10*(Kappa00*Kappa10
      + Kappa01*Kappa11 + Kappa02*Kappa12) + Kappa20*(Kappa00*Kappa20 + Kappa01*
      Kappa21 + Kappa02*Kappa22) + Kappa00*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(
      Kappa02))) + Kappa00*(-0.26666666666666666*Sqr(g1) - 5.333333333333333*Sqr(
      g3) + 3*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(
      Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*(
      Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 2*
      Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QDxbarp) - 2*Sqr(gN)*Sqr(QDxp) - 2*Sqr(gN)*Sqr(
      QSp))) + 6*Kappa01*Lambdax*(2*(Kappa11*(Kappa00*Kappa10 + Kappa01*Kappa11 +
      Kappa02*Kappa12) + Kappa21*(Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*
      Kappa22) + Kappa01*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02))) + Kappa01*(
      -0.26666666666666666*Sqr(g1) - 5.333333333333333*Sqr(g3) + 3*(Sqr(Kappa00) +
      Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) +
      Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*(Sqr(Lambda1200) + Sqr(
      Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 2*Sqr(Lambdax) - 2*Sqr(gN
      )*Sqr(QDxbarp) - 2*Sqr(gN)*Sqr(QDxp) - 2*Sqr(gN)*Sqr(QSp))) + 6*Kappa02*
      Lambdax*(2*(Kappa12*(Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12) +
      Kappa22*(Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22) + Kappa02*(Sqr
      (Kappa00) + Sqr(Kappa01) + Sqr(Kappa02))) + Kappa02*(-0.26666666666666666*
      Sqr(g1) - 5.333333333333333*Sqr(g3) + 3*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(
      Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(
      Kappa21) + Sqr(Kappa22)) + 2*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(
      Lambda1210) + Sqr(Lambda1211)) + 2*Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QDxbarp) - 2
      *Sqr(gN)*Sqr(QDxp) - 2*Sqr(gN)*Sqr(QSp))) + 6*Kappa10*Lambdax*(2*(Kappa00*(
      Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12) + Kappa20*(Kappa10*
      Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22) + Kappa10*(Sqr(Kappa10) + Sqr(
      Kappa11) + Sqr(Kappa12))) + Kappa10*(-0.26666666666666666*Sqr(g1) -
      5.333333333333333*Sqr(g3) + 3*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) +
      Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) +
      Sqr(Kappa22)) + 2*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr
      (Lambda1211)) + 2*Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QDxbarp) - 2*Sqr(gN)*Sqr(QDxp
      ) - 2*Sqr(gN)*Sqr(QSp))) + 6*Kappa11*Lambdax*(2*(Kappa01*(Kappa00*Kappa10 +
      Kappa01*Kappa11 + Kappa02*Kappa12) + Kappa21*(Kappa10*Kappa20 + Kappa11*
      Kappa21 + Kappa12*Kappa22) + Kappa11*(Sqr(Kappa10) + Sqr(Kappa11) + Sqr(
      Kappa12))) + Kappa11*(-0.26666666666666666*Sqr(g1) - 5.333333333333333*Sqr(
      g3) + 3*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(
      Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*(
      Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 2*
      Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QDxbarp) - 2*Sqr(gN)*Sqr(QDxp) - 2*Sqr(gN)*Sqr(
      QSp))) + 6*Kappa12*Lambdax*(2*(Kappa02*(Kappa00*Kappa10 + Kappa01*Kappa11 +
      Kappa02*Kappa12) + Kappa22*(Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*
      Kappa22) + Kappa12*(Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12))) + Kappa12*(
      -0.26666666666666666*Sqr(g1) - 5.333333333333333*Sqr(g3) + 3*(Sqr(Kappa00) +
      Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) +
      Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*(Sqr(Lambda1200) + Sqr(
      Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 2*Sqr(Lambdax) - 2*Sqr(gN
      )*Sqr(QDxbarp) - 2*Sqr(gN)*Sqr(QDxp) - 2*Sqr(gN)*Sqr(QSp))) + 6*Kappa20*
      Lambdax*(2*(Kappa00*(Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22) +
      Kappa10*(Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22) + Kappa20*(Sqr
      (Kappa20) + Sqr(Kappa21) + Sqr(Kappa22))) + Kappa20*(-0.26666666666666666*
      Sqr(g1) - 5.333333333333333*Sqr(g3) + 3*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(
      Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(
      Kappa21) + Sqr(Kappa22)) + 2*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(
      Lambda1210) + Sqr(Lambda1211)) + 2*Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QDxbarp) - 2
      *Sqr(gN)*Sqr(QDxp) - 2*Sqr(gN)*Sqr(QSp))) + 6*Kappa21*Lambdax*(2*(Kappa01*(
      Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22) + Kappa11*(Kappa10*
      Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22) + Kappa21*(Sqr(Kappa20) + Sqr(
      Kappa21) + Sqr(Kappa22))) + Kappa21*(-0.26666666666666666*Sqr(g1) -
      5.333333333333333*Sqr(g3) + 3*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) +
      Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) +
      Sqr(Kappa22)) + 2*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr
      (Lambda1211)) + 2*Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QDxbarp) - 2*Sqr(gN)*Sqr(QDxp
      ) - 2*Sqr(gN)*Sqr(QSp))) + 6*Kappa22*Lambdax*(2*(Kappa02*(Kappa00*Kappa20 +
      Kappa01*Kappa21 + Kappa02*Kappa22) + Kappa12*(Kappa10*Kappa20 + Kappa11*
      Kappa21 + Kappa12*Kappa22) + Kappa22*(Sqr(Kappa20) + Sqr(Kappa21) + Sqr(
      Kappa22))) + Kappa22*(-0.26666666666666666*Sqr(g1) - 5.333333333333333*Sqr(
      g3) + 3*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(
      Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*(
      Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 2*
      Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QDxbarp) - 2*Sqr(gN)*Sqr(QDxp) - 2*Sqr(gN)*Sqr(
      QSp))) + 4*Lambda1200*Lambdax*(2*(Lambda1210*(Lambda1200*Lambda1210 +
      Lambda1201*Lambda1211) + Lambda1200*(Sqr(Lambda1200) + Sqr(Lambda1201))) +
      Lambda1200*(-0.6*Sqr(g1) - 3*Sqr(g2) + 3*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(
      Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(
      Kappa21) + Sqr(Kappa22)) + 2*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(
      Lambda1210) + Sqr(Lambda1211)) + 2*Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH1p) - 2*
      Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QSp))) + 4*Lambda1201*Lambdax*(2*(
      Lambda1211*(Lambda1200*Lambda1210 + Lambda1201*Lambda1211) + Lambda1201*(Sqr
      (Lambda1200) + Sqr(Lambda1201))) + Lambda1201*(-0.6*Sqr(g1) - 3*Sqr(g2) + 3*
      (Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) +
      Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*(Sqr(
      Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 2*Sqr(
      Lambdax) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QSp)))
      + 4*Lambda1210*Lambdax*(2*(Lambda1200*(Lambda1200*Lambda1210 + Lambda1201*
      Lambda1211) + Lambda1210*(Sqr(Lambda1210) + Sqr(Lambda1211))) + Lambda1210*(
      -0.6*Sqr(g1) - 3*Sqr(g2) + 3*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) +
      Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) +
      Sqr(Kappa22)) + 2*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr
      (Lambda1211)) + 2*Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QH2p) -
      2*Sqr(gN)*Sqr(QSp))) + 4*Lambda1211*Lambdax*(2*(Lambda1201*(Lambda1200*
      Lambda1210 + Lambda1201*Lambda1211) + Lambda1211*(Sqr(Lambda1210) + Sqr(
      Lambda1211))) + Lambda1211*(-0.6*Sqr(g1) - 3*Sqr(g2) + 3*(Sqr(Kappa00) + Sqr
      (Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(
      Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*(Sqr(Lambda1200) + Sqr(
      Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 2*Sqr(Lambdax) - 2*Sqr(gN
      )*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QSp))) + Power(gN,3)*
      Lambdax*(-4*gN*Sqr(QH1p) - 4*gN*Sqr(QH2p) - 4*gN*Sqr(QSp))*(9*Sqr(Qdp) + 9*
      Sqr(QDxbarp) + 9*Sqr(QDxp) + 3*Sqr(Qep) + 6*Sqr(QH1p) + 6*Sqr(QH2p) + 2*Sqr(
      QHpbarp) + 2*Sqr(QHpp) + 6*Sqr(QLp) + 18*Sqr(QQp) + 3*Sqr(QSp) + 9*Sqr(Qup))
      + 12*Power(Lambdax,3)*Sqr(Yd00) + 12*Power(Lambdax,3)*Sqr(Yd01) + 12*Power(
      Lambdax,3)*Sqr(Yd02) + 12*Power(Lambdax,3)*Sqr(Yd10) + 12*Power(Lambdax,3)*
      Sqr(Yd11) + 12*Power(Lambdax,3)*Sqr(Yd12) + 12*Power(Lambdax,3)*Sqr(Yd20) +
      12*Power(Lambdax,3)*Sqr(Yd21) + 12*Power(Lambdax,3)*Sqr(Yd22) + 4*Power(
      Lambdax,3)*Sqr(Ye00) + 4*Power(Lambdax,3)*Sqr(Ye01) + 4*Power(Lambdax,3)*Sqr
      (Ye02) + 4*Power(Lambdax,3)*Sqr(Ye10) + 4*Power(Lambdax,3)*Sqr(Ye11) + 4*
      Power(Lambdax,3)*Sqr(Ye12) + 4*Power(Lambdax,3)*Sqr(Ye20) + 4*Power(Lambdax,
      3)*Sqr(Ye21) + 4*Power(Lambdax,3)*Sqr(Ye22) + 2*Lambdax*Ye00*(3*(Ye10*(Ye00*
      Ye10 + Ye01*Ye11 + Ye02*Ye12) + Ye20*(Ye00*Ye20 + Ye01*Ye21 + Ye02*Ye22) +
      Ye00*(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02))) + Ye00*(-1.8*Sqr(g1) - 3*Sqr(g2) +
      Sqr(Lambdax) - 2*Sqr(gN)*Sqr(Qep) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QLp
      ) + 3*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12)
      + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) +
      Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22))) + 2*
      Lambdax*Ye01*(3*(Ye11*(Ye00*Ye10 + Ye01*Ye11 + Ye02*Ye12) + Ye21*(Ye00*Ye20
      + Ye01*Ye21 + Ye02*Ye22) + Ye01*(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02))) + Ye01*
      (-1.8*Sqr(g1) - 3*Sqr(g2) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(Qep) - 2*Sqr(gN)*
      Sqr(QH1p) - 2*Sqr(gN)*Sqr(QLp) + 3*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(
      Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(
      Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20
      ) + Sqr(Ye21) + Sqr(Ye22))) + 2*Lambdax*Ye02*(3*(Ye12*(Ye00*Ye10 + Ye01*Ye11
      + Ye02*Ye12) + Ye22*(Ye00*Ye20 + Ye01*Ye21 + Ye02*Ye22) + Ye02*(Sqr(Ye00) +
      Sqr(Ye01) + Sqr(Ye02))) + Ye02*(-1.8*Sqr(g1) - 3*Sqr(g2) + Sqr(Lambdax) - 2
      *Sqr(gN)*Sqr(Qep) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QLp) + 3*(Sqr(Yd00)
      + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) +
      Sqr(Yd21) + Sqr(Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr
      (Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22))) + 2*Lambdax*Ye10*(3
      *(Ye00*(Ye00*Ye10 + Ye01*Ye11 + Ye02*Ye12) + Ye20*(Ye10*Ye20 + Ye11*Ye21 +
      Ye12*Ye22) + Ye10*(Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12))) + Ye10*(-1.8*Sqr(g1)
      - 3*Sqr(g2) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(Qep) - 2*Sqr(gN)*Sqr(QH1p) - 2*
      Sqr(gN)*Sqr(QLp) + 3*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(
      Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(Ye00) + Sqr(
      Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21
      ) + Sqr(Ye22))) + 2*Lambdax*Ye11*(3*(Ye01*(Ye00*Ye10 + Ye01*Ye11 + Ye02*Ye12
      ) + Ye21*(Ye10*Ye20 + Ye11*Ye21 + Ye12*Ye22) + Ye11*(Sqr(Ye10) + Sqr(Ye11) +
      Sqr(Ye12))) + Ye11*(-1.8*Sqr(g1) - 3*Sqr(g2) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr
      (Qep) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QLp) + 3*(Sqr(Yd00) + Sqr(Yd01)
      + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) +
      Sqr(Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr
      (Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22))) + 2*Lambdax*Ye12*(3*(Ye02*(Ye00
      *Ye10 + Ye01*Ye11 + Ye02*Ye12) + Ye22*(Ye10*Ye20 + Ye11*Ye21 + Ye12*Ye22) +
      Ye12*(Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12))) + Ye12*(-1.8*Sqr(g1) - 3*Sqr(g2) +
      Sqr(Lambdax) - 2*Sqr(gN)*Sqr(Qep) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QLp
      ) + 3*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12)
      + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) +
      Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22))) + 6*
      Lambdax*Yd00*(Yu00*(Yd00*Yu00 + Yd01*Yu01 + Yd02*Yu02) + Yu10*(Yd00*Yu10 +
      Yd01*Yu11 + Yd02*Yu12) + Yu20*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22) + 3*(Yd10*
      (Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + Yd20*(Yd00*Yd20 + Yd01*Yd21 + Yd02*
      Yd22) + Yd00*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02))) + Yd00*(
      -0.4666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(
      Lambdax) - 2*Sqr(gN)*Sqr(Qdp) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QQp) + 3
      *(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) +
      Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr
      (Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22))) + 6*
      Lambdax*Yd01*(Yu01*(Yd00*Yu00 + Yd01*Yu01 + Yd02*Yu02) + Yu11*(Yd00*Yu10 +
      Yd01*Yu11 + Yd02*Yu12) + Yu21*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22) + 3*(Yd11*
      (Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + Yd21*(Yd00*Yd20 + Yd01*Yd21 + Yd02*
      Yd22) + Yd01*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02))) + Yd01*(
      -0.4666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(
      Lambdax) - 2*Sqr(gN)*Sqr(Qdp) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QQp) + 3
      *(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) +
      Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr
      (Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22))) + 6*
      Lambdax*Yd02*(Yu02*(Yd00*Yu00 + Yd01*Yu01 + Yd02*Yu02) + Yu12*(Yd00*Yu10 +
      Yd01*Yu11 + Yd02*Yu12) + Yu22*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22) + 3*(Yd12*
      (Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + Yd22*(Yd00*Yd20 + Yd01*Yd21 + Yd02*
      Yd22) + Yd02*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02))) + Yd02*(
      -0.4666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(
      Lambdax) - 2*Sqr(gN)*Sqr(Qdp) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QQp) + 3
      *(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) +
      Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr
      (Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22))) + 6*
      Lambdax*Yd10*(Yu00*(Yd10*Yu00 + Yd11*Yu01 + Yd12*Yu02) + Yu10*(Yd10*Yu10 +
      Yd11*Yu11 + Yd12*Yu12) + Yu20*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22) + 3*(Yd00*
      (Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + Yd20*(Yd10*Yd20 + Yd11*Yd21 + Yd12*
      Yd22) + Yd10*(Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12))) + Yd10*(
      -0.4666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(
      Lambdax) - 2*Sqr(gN)*Sqr(Qdp) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QQp) + 3
      *(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) +
      Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr
      (Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22))) + 6*
      Lambdax*Yd11*(Yu01*(Yd10*Yu00 + Yd11*Yu01 + Yd12*Yu02) + Yu11*(Yd10*Yu10 +
      Yd11*Yu11 + Yd12*Yu12) + Yu21*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22) + 3*(Yd01*
      (Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + Yd21*(Yd10*Yd20 + Yd11*Yd21 + Yd12*
      Yd22) + Yd11*(Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12))) + Yd11*(
      -0.4666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(
      Lambdax) - 2*Sqr(gN)*Sqr(Qdp) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QQp) + 3
      *(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) +
      Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr
      (Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22))) + 6*
      Lambdax*Yd12*(Yu02*(Yd10*Yu00 + Yd11*Yu01 + Yd12*Yu02) + Yu12*(Yd10*Yu10 +
      Yd11*Yu11 + Yd12*Yu12) + Yu22*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22) + 3*(Yd02*
      (Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + Yd22*(Yd10*Yd20 + Yd11*Yd21 + Yd12*
      Yd22) + Yd12*(Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12))) + Yd12*(
      -0.4666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(
      Lambdax) - 2*Sqr(gN)*Sqr(Qdp) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QQp) + 3
      *(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) +
      Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr
      (Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22))) + 6*
      Lambdax*Yd20*(Yu00*(Yd20*Yu00 + Yd21*Yu01 + Yd22*Yu02) + Yu10*(Yd20*Yu10 +
      Yd21*Yu11 + Yd22*Yu12) + Yu20*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22) + 3*(Yd00*
      (Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22) + Yd10*(Yd10*Yd20 + Yd11*Yd21 + Yd12*
      Yd22) + Yd20*(Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22))) + Yd20*(
      -0.4666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(
      Lambdax) - 2*Sqr(gN)*Sqr(Qdp) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QQp) + 3
      *(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) +
      Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr
      (Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22))) + 6*
      Lambdax*Yd21*(Yu01*(Yd20*Yu00 + Yd21*Yu01 + Yd22*Yu02) + Yu11*(Yd20*Yu10 +
      Yd21*Yu11 + Yd22*Yu12) + Yu21*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22) + 3*(Yd01*
      (Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22) + Yd11*(Yd10*Yd20 + Yd11*Yd21 + Yd12*
      Yd22) + Yd21*(Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22))) + Yd21*(
      -0.4666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(
      Lambdax) - 2*Sqr(gN)*Sqr(Qdp) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QQp) + 3
      *(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) +
      Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr
      (Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22))) + 6*
      Lambdax*Yd22*(Yu02*(Yd20*Yu00 + Yd21*Yu01 + Yd22*Yu02) + Yu12*(Yd20*Yu10 +
      Yd21*Yu11 + Yd22*Yu12) + Yu22*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22) + 3*(Yd02*
      (Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22) + Yd12*(Yd10*Yd20 + Yd11*Yd21 + Yd12*
      Yd22) + Yd22*(Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22))) + Yd22*(
      -0.4666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(
      Lambdax) - 2*Sqr(gN)*Sqr(Qdp) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QQp) + 3
      *(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) +
      Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr
      (Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22))) + 2*
      Lambdax*Ye20*(Ye20*(-1.8*Sqr(g1) - 3*Sqr(g2) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(
      Qep) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QLp) + 3*(Sqr(Yd00) + Sqr(Yd01) +
      Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr
      (Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(
      Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22)) + 3*(Ye00*(Ye00*Ye20 + Ye01*Ye21
      + Ye02*Ye22) + Ye10*(Ye10*Ye20 + Ye11*Ye21 + Ye12*Ye22) + Ye20*(Sqr(Ye20) +
      Sqr(Ye21) + Sqr(Ye22)))) + 2*Lambdax*Ye21*(Ye21*(-1.8*Sqr(g1) - 3*Sqr(g2) +
      Sqr(Lambdax) - 2*Sqr(gN)*Sqr(Qep) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QLp)
      + 3*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12)
      + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) +
      Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22)) + 3*(
      Ye01*(Ye00*Ye20 + Ye01*Ye21 + Ye02*Ye22) + Ye11*(Ye10*Ye20 + Ye11*Ye21 +
      Ye12*Ye22) + Ye21*(Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22)))) + 2*Lambdax*Ye22*(
      Ye22*(-1.8*Sqr(g1) - 3*Sqr(g2) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(Qep) - 2*Sqr(
      gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QLp) + 3*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) +
      Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr
      (Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(
      Ye20) + Sqr(Ye21) + Sqr(Ye22)) + 3*(Ye02*(Ye00*Ye20 + Ye01*Ye21 + Ye02*Ye22)
      + Ye12*(Ye10*Ye20 + Ye11*Ye21 + Ye12*Ye22) + Ye22*(Sqr(Ye20) + Sqr(Ye21) +
      Sqr(Ye22)))) + 12*Power(Lambdax,3)*Sqr(Yu00) + 12*Power(Lambdax,3)*Sqr(Yu01)
      + 12*Power(Lambdax,3)*Sqr(Yu02) + 12*Power(Lambdax,3)*Sqr(Yu10) + 12*Power(
      Lambdax,3)*Sqr(Yu11) + 12*Power(Lambdax,3)*Sqr(Yu12) + 12*Power(Lambdax,3)*
      Sqr(Yu20) + 12*Power(Lambdax,3)*Sqr(Yu21) + 12*Power(Lambdax,3)*Sqr(Yu22) +
      24*Sqr(Lambdax)*(4*Power(Lambdax,3) - 0.6*Lambdax*Sqr(g1) - 3*Lambdax*Sqr(g2
      ) + 3*Lambdax*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) +
      Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) +
      2*Lambdax*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(
      Lambda1211)) - 2*Lambdax*Sqr(gN)*Sqr(QH1p) - 2*Lambdax*Sqr(gN)*Sqr(QH2p) - 2
      *Lambdax*Sqr(gN)*Sqr(QSp) + 3*Lambdax*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) +
      Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) +
      Lambdax*(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(
      Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22)) + 3*Lambdax*(Sqr(Yu00) + Sqr(Yu01
      ) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) +
      Sqr(Yu22))) + 6*Lambdax*Yu00*(Yd00*(Yd00*Yu00 + Yd01*Yu01 + Yd02*Yu02) +
      Yd10*(Yd10*Yu00 + Yd11*Yu01 + Yd12*Yu02) + Yd20*(Yd20*Yu00 + Yd21*Yu01 +
      Yd22*Yu02) + 3*(Yu10*(Yu00*Yu10 + Yu01*Yu11 + Yu02*Yu12) + Yu20*(Yu00*Yu20 +
      Yu01*Yu21 + Yu02*Yu22) + Yu00*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02))) + Yu00*(
      -0.8666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(
      Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3
      *(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) +
      Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) + 6*Lambdax*Yu01*(Yd01*(Yd00*Yu00 +
      Yd01*Yu01 + Yd02*Yu02) + Yd11*(Yd10*Yu00 + Yd11*Yu01 + Yd12*Yu02) + Yd21*(
      Yd20*Yu00 + Yd21*Yu01 + Yd22*Yu02) + 3*(Yu11*(Yu00*Yu10 + Yu01*Yu11 + Yu02*
      Yu12) + Yu21*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu01*(Sqr(Yu00) + Sqr(
      Yu01) + Sqr(Yu02))) + Yu01*(-0.8666666666666667*Sqr(g1) - 3*Sqr(g2) -
      5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*
      Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(
      Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) + 6*
      Lambdax*Yu02*(Yd02*(Yd00*Yu00 + Yd01*Yu01 + Yd02*Yu02) + Yd12*(Yd10*Yu00 +
      Yd11*Yu01 + Yd12*Yu02) + Yd22*(Yd20*Yu00 + Yd21*Yu01 + Yd22*Yu02) + 3*(Yu12*
      (Yu00*Yu10 + Yu01*Yu11 + Yu02*Yu12) + Yu22*(Yu00*Yu20 + Yu01*Yu21 + Yu02*
      Yu22) + Yu02*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02))) + Yu02*(
      -0.8666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(
      Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3
      *(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) +
      Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) + 6*Lambdax*Yu10*(Yd00*(Yd00*Yu10 +
      Yd01*Yu11 + Yd02*Yu12) + Yd10*(Yd10*Yu10 + Yd11*Yu11 + Yd12*Yu12) + Yd20*(
      Yd20*Yu10 + Yd21*Yu11 + Yd22*Yu12) + 3*(Yu00*(Yu00*Yu10 + Yu01*Yu11 + Yu02*
      Yu12) + Yu20*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + Yu10*(Sqr(Yu10) + Sqr(
      Yu11) + Sqr(Yu12))) + Yu10*(-0.8666666666666667*Sqr(g1) - 3*Sqr(g2) -
      5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*
      Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(
      Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) + 6*
      Lambdax*Yu11*(Yd01*(Yd00*Yu10 + Yd01*Yu11 + Yd02*Yu12) + Yd11*(Yd10*Yu10 +
      Yd11*Yu11 + Yd12*Yu12) + Yd21*(Yd20*Yu10 + Yd21*Yu11 + Yd22*Yu12) + 3*(Yu01*
      (Yu00*Yu10 + Yu01*Yu11 + Yu02*Yu12) + Yu21*(Yu10*Yu20 + Yu11*Yu21 + Yu12*
      Yu22) + Yu11*(Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12))) + Yu11*(
      -0.8666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(
      Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3
      *(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) +
      Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) + 6*Lambdax*Yu12*(Yd02*(Yd00*Yu10 +
      Yd01*Yu11 + Yd02*Yu12) + Yd12*(Yd10*Yu10 + Yd11*Yu11 + Yd12*Yu12) + Yd22*(
      Yd20*Yu10 + Yd21*Yu11 + Yd22*Yu12) + 3*(Yu02*(Yu00*Yu10 + Yu01*Yu11 + Yu02*
      Yu12) + Yu22*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + Yu12*(Sqr(Yu10) + Sqr(
      Yu11) + Sqr(Yu12))) + Yu12*(-0.8666666666666667*Sqr(g1) - 3*Sqr(g2) -
      5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*
      Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(
      Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) + 6*
      Lambdax*Yu20*(Yd00*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22) + Yd10*(Yd10*Yu20 +
      Yd11*Yu21 + Yd12*Yu22) + Yd20*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22) + 3*(Yu00*
      (Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu10*(Yu10*Yu20 + Yu11*Yu21 + Yu12*
      Yu22) + Yu20*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + Yu20*(
      -0.8666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(
      Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3
      *(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) +
      Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) + 6*Lambdax*Yu21*(Yd01*(Yd00*Yu20 +
      Yd01*Yu21 + Yd02*Yu22) + Yd11*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22) + Yd21*(
      Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22) + 3*(Yu01*(Yu00*Yu20 + Yu01*Yu21 + Yu02*
      Yu22) + Yu11*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + Yu21*(Sqr(Yu20) + Sqr(
      Yu21) + Sqr(Yu22))) + Yu21*(-0.8666666666666667*Sqr(g1) - 3*Sqr(g2) -
      5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*
      Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(
      Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) + 6*
      Lambdax*Yu22*(Yd02*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22) + Yd12*(Yd10*Yu20 +
      Yd11*Yu21 + Yd12*Yu22) + Yd22*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22) + 3*(Yu02*
      (Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu12*(Yu10*Yu20 + Yu11*Yu21 + Yu12*
      Yu22) + Yu22*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + Yu22*(
      -0.8666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(
      Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3
      *(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) +
      Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) + Lambdax*Sqr(-0.6*Power(g1,2) - 3*
      Power(g2,2) + 3*(Power(Kappa00,2) + Power(Kappa01,2) + Power(Kappa02,2) +
      Power(Kappa10,2) + Power(Kappa11,2) + Power(Kappa12,2) + Power(Kappa20,2) +
      Power(Kappa21,2) + Power(Kappa22,2)) + 2*(Power(Lambda1200,2) + Power(
      Lambda1201,2) + Power(Lambda1210,2) + Power(Lambda1211,2)) + 12*Power(
      Lambdax,2) - 2*Power(gN,2)*Power(QH1p,2) - 2*Power(gN,2)*Power(QH2p,2) - 2*
      Power(gN,2)*Power(QSp,2) + 3*(Power(Yd00,2) + Power(Yd01,2) + Power(Yd02,2)
      + Power(Yd10,2) + Power(Yd11,2) + Power(Yd12,2) + Power(Yd20,2) + Power(Yd21
      ,2) + Power(Yd22,2)) + Power(Ye00,2) + Power(Ye01,2) + Power(Ye02,2) + Power
      (Ye10,2) + Power(Ye11,2) + Power(Ye12,2) + Power(Ye20,2) + Power(Ye21,2) +
      Power(Ye22,2) + 3*(Power(Yu00,2) + Power(Yu01,2) + Power(Yu02,2) + Power(
      Yu10,2) + Power(Yu11,2) + Power(Yu12,2) + Power(Yu20,2) + Power(Yu21,2) +
      Power(Yu22,2)));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dlead_log_one_loop_TLambdax_dAYu22() const 
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

   double deriv = 12*Lambdax*Yd02*Yu22*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22) + 12
      *Lambdax*Yd12*Yu22*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22) + 12*Lambdax*Yd22*
      Yu22*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22) + 36*Lambdax*Sqr(Yu00)*Sqr(Yu22) +
      36*Lambdax*Sqr(Yu01)*Sqr(Yu22) + 36*Lambdax*Sqr(Yu10)*Sqr(Yu22) + 36*Lambdax
      *Sqr(Yu11)*Sqr(Yu22) + 6*Lambdax*Yu02*(4*Yu22*(Yu00*Yu20 + Yu01*Yu21 + Yu02*
      Yu22) + 6*Yu02*Sqr(Yu22)) + 6*Lambdax*Yu12*(4*Yu22*(Yu10*Yu20 + Yu11*Yu21 +
      Yu12*Yu22) + 6*Yu12*Sqr(Yu22)) + 6*Lambdax*Sqr(Yu22)*(-0.6*Sqr(g1) - 3*Sqr(
      g2) + 3*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(
      Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*(
      Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 12*
      Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QSp
      ) + 3*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12)
      + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) +
      Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22) + 3*(
      Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(
      Yu20) + Sqr(Yu21) + Sqr(Yu22))) + 6*Sqr(Yu22)*(4*Power(Lambdax,3) - 0.6*
      Lambdax*Sqr(g1) - 3*Lambdax*Sqr(g2) + 3*Lambdax*(Sqr(Kappa00) + Sqr(Kappa01)
      + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20)
      + Sqr(Kappa21) + Sqr(Kappa22)) + 2*Lambdax*(Sqr(Lambda1200) + Sqr(Lambda1201
      ) + Sqr(Lambda1210) + Sqr(Lambda1211)) - 2*Lambdax*Sqr(gN)*Sqr(QH1p) - 2*
      Lambdax*Sqr(gN)*Sqr(QH2p) - 2*Lambdax*Sqr(gN)*Sqr(QSp) + 3*Lambdax*(Sqr(Yd00
      ) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) +
      Sqr(Yd21) + Sqr(Yd22)) + Lambdax*(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(
      Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22)) + 3*
      Lambdax*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(
      Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + 6*Lambdax*Yu22*(6*Power(Yu22,3
      ) - 0.8666666666666667*Yu22*Sqr(g1) - 3*Yu22*Sqr(g2) - 5.333333333333333*
      Yu22*Sqr(g3) + Yu22*Sqr(Lambdax) - 2*Yu22*Sqr(gN)*Sqr(QH2p) - 2*Yu22*Sqr(gN)
      *Sqr(QQp) - 2*Yu22*Sqr(gN)*Sqr(Qup) + Yu22*Sqr(Yd02) + Yu22*Sqr(Yd12) + Yu22
      *Sqr(Yd22) + 5*(Power(Yu22,3) + Yu22*Sqr(Yu02) + Yu22*Sqr(Yu12)) + 4*Yu22*(
      Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)) + 3*Yu22*(Sqr(Yu00) + Sqr(Yu01) + Sqr(
      Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22
      ))) + 6*Lambdax*Yu20*(Yd00*Yd02*Yu22 + Yd10*Yd12*Yu22 + Yd20*Yd22*Yu22 + 6*
      Yu20*Sqr(Yu22) + 5*(Yu00*Yu02*Yu22 + Yu10*Yu12*Yu22 + Yu20*Sqr(Yu22))) + 6*
      Lambdax*Yu21*(Yd01*Yd02*Yu22 + Yd11*Yd12*Yu22 + Yd21*Yd22*Yu22 + 6*Yu21*Sqr(
      Yu22) + 5*(Yu01*Yu02*Yu22 + Yu11*Yu12*Yu22 + Yu21*Sqr(Yu22))) + 6*Lambdax*
      Yu22*(Yd02*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22) + Yd12*(Yd10*Yu20 + Yd11*Yu21
      + Yd12*Yu22) + Yd22*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22) + 3*(Yu02*(Yu00*
      Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu12*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) +
      Yu22*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + Yu22*(-0.8666666666666667*Sqr(g1
      ) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(
      QH2p) - 2*Sqr(gN)*Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3*(Sqr(Yu00) + Sqr(Yu01) +
      Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr
      (Yu22))));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dlead_log_one_loop_TLambdax_dMassB() const 
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
   const double g1 = model.get_g1();
   const double g2 = model.get_g2();
   const double gN = model.get_gN();

   double deriv = 46.08*Power(g1,4)*Lambdax + 3.2*Lambdax*Sqr(g1)*Sqr(Kappa00)
      + 3.2*Lambdax*Sqr(g1)*Sqr(Kappa01) + 3.2*Lambdax*Sqr(g1)*Sqr(Kappa02) + 3.2*
      Lambdax*Sqr(g1)*Sqr(Kappa10) + 3.2*Lambdax*Sqr(g1)*Sqr(Kappa11) + 3.2*
      Lambdax*Sqr(g1)*Sqr(Kappa12) + 3.2*Lambdax*Sqr(g1)*Sqr(Kappa20) + 3.2*
      Lambdax*Sqr(g1)*Sqr(Kappa21) + 3.2*Lambdax*Sqr(g1)*Sqr(Kappa22) + 4.8*
      Lambdax*Sqr(g1)*Sqr(Lambda1200) + 4.8*Lambdax*Sqr(g1)*Sqr(Lambda1201) + 4.8*
      Lambdax*Sqr(g1)*Sqr(Lambda1210) + 4.8*Lambdax*Sqr(g1)*Sqr(Lambda1211) + 5.6*
      Lambdax*Sqr(g1)*Sqr(Yd00) + 5.6*Lambdax*Sqr(g1)*Sqr(Yd01) + 5.6*Lambdax*Sqr(
      g1)*Sqr(Yd02) + 5.6*Lambdax*Sqr(g1)*Sqr(Yd10) + 5.6*Lambdax*Sqr(g1)*Sqr(Yd11
      ) + 5.6*Lambdax*Sqr(g1)*Sqr(Yd12) + 5.6*Lambdax*Sqr(g1)*Sqr(Yd20) + 5.6*
      Lambdax*Sqr(g1)*Sqr(Yd21) + 5.6*Lambdax*Sqr(g1)*Sqr(Yd22) + 7.2*Lambdax*Sqr(
      g1)*Sqr(Ye00) + 7.2*Lambdax*Sqr(g1)*Sqr(Ye01) + 7.2*Lambdax*Sqr(g1)*Sqr(Ye02
      ) + 7.2*Lambdax*Sqr(g1)*Sqr(Ye10) + 7.2*Lambdax*Sqr(g1)*Sqr(Ye11) + 7.2*
      Lambdax*Sqr(g1)*Sqr(Ye12) + 7.2*Lambdax*Sqr(g1)*Sqr(Ye20) + 7.2*Lambdax*Sqr(
      g1)*Sqr(Ye21) + 7.2*Lambdax*Sqr(g1)*Sqr(Ye22) + 10.4*Lambdax*Sqr(g1)*Sqr(
      Yu00) + 10.4*Lambdax*Sqr(g1)*Sqr(Yu01) + 10.4*Lambdax*Sqr(g1)*Sqr(Yu02) +
      10.4*Lambdax*Sqr(g1)*Sqr(Yu10) + 10.4*Lambdax*Sqr(g1)*Sqr(Yu11) + 10.4*
      Lambdax*Sqr(g1)*Sqr(Yu12) + 10.4*Lambdax*Sqr(g1)*Sqr(Yu20) + 10.4*Lambdax*
      Sqr(g1)*Sqr(Yu21) + 10.4*Lambdax*Sqr(g1)*Sqr(Yu22) + 1.2*Lambdax*Sqr(g1)*(
      -0.6*Sqr(g1) - 3*Sqr(g2) + 3*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) +
      Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) +
      Sqr(Kappa22)) + 2*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr
      (Lambda1211)) + 12*Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QH2p)
      - 2*Sqr(gN)*Sqr(QSp) + 3*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) +
      Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(Ye00) + Sqr
      (Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(
      Ye21) + Sqr(Ye22) + 3*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(
      Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + 1.2*Sqr(g1)*(4*
      Power(Lambdax,3) - 0.6*Lambdax*Sqr(g1) - 3*Lambdax*Sqr(g2) + 3*Lambdax*(Sqr(
      Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(
      Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*Lambdax*(Sqr(
      Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) - 2*
      Lambdax*Sqr(gN)*Sqr(QH1p) - 2*Lambdax*Sqr(gN)*Sqr(QH2p) - 2*Lambdax*Sqr(gN)*
      Sqr(QSp) + 3*Lambdax*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(
      Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Lambdax*(Sqr(Ye00)
      + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) +
      Sqr(Ye21) + Sqr(Ye22)) + 3*Lambdax*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(
      Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dlead_log_one_loop_TLambdax_dMassWB() const 
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
   const double g1 = model.get_g1();
   const double g2 = model.get_g2();
   const double gN = model.get_gN();

   double deriv = 96*Power(g2,4)*Lambdax + 24*Lambdax*Sqr(g2)*Sqr(Lambda1200) +
      24*Lambdax*Sqr(g2)*Sqr(Lambda1201) + 24*Lambdax*Sqr(g2)*Sqr(Lambda1210) +
      24*Lambdax*Sqr(g2)*Sqr(Lambda1211) + 36*Lambdax*Sqr(g2)*Sqr(Yd00) + 36*
      Lambdax*Sqr(g2)*Sqr(Yd01) + 36*Lambdax*Sqr(g2)*Sqr(Yd02) + 36*Lambdax*Sqr(g2
      )*Sqr(Yd10) + 36*Lambdax*Sqr(g2)*Sqr(Yd11) + 36*Lambdax*Sqr(g2)*Sqr(Yd12) +
      36*Lambdax*Sqr(g2)*Sqr(Yd20) + 36*Lambdax*Sqr(g2)*Sqr(Yd21) + 36*Lambdax*Sqr
      (g2)*Sqr(Yd22) + 12*Lambdax*Sqr(g2)*Sqr(Ye00) + 12*Lambdax*Sqr(g2)*Sqr(Ye01)
      + 12*Lambdax*Sqr(g2)*Sqr(Ye02) + 12*Lambdax*Sqr(g2)*Sqr(Ye10) + 12*Lambdax*
      Sqr(g2)*Sqr(Ye11) + 12*Lambdax*Sqr(g2)*Sqr(Ye12) + 12*Lambdax*Sqr(g2)*Sqr(
      Ye20) + 12*Lambdax*Sqr(g2)*Sqr(Ye21) + 12*Lambdax*Sqr(g2)*Sqr(Ye22) + 36*
      Lambdax*Sqr(g2)*Sqr(Yu00) + 36*Lambdax*Sqr(g2)*Sqr(Yu01) + 36*Lambdax*Sqr(g2
      )*Sqr(Yu02) + 36*Lambdax*Sqr(g2)*Sqr(Yu10) + 36*Lambdax*Sqr(g2)*Sqr(Yu11) +
      36*Lambdax*Sqr(g2)*Sqr(Yu12) + 36*Lambdax*Sqr(g2)*Sqr(Yu20) + 36*Lambdax*Sqr
      (g2)*Sqr(Yu21) + 36*Lambdax*Sqr(g2)*Sqr(Yu22) + 6*Lambdax*Sqr(g2)*(-0.6*Sqr(
      g1) - 3*Sqr(g2) + 3*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(
      Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(
      Kappa22)) + 2*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(
      Lambda1211)) + 12*Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QH2p) -
      2*Sqr(gN)*Sqr(QSp) + 3*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr
      (Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(Ye00) + Sqr(
      Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21
      ) + Sqr(Ye22) + 3*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11)
      + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + 6*Sqr(g2)*(4*Power(
      Lambdax,3) - 0.6*Lambdax*Sqr(g1) - 3*Lambdax*Sqr(g2) + 3*Lambdax*(Sqr(
      Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(
      Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*Lambdax*(Sqr(
      Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) - 2*
      Lambdax*Sqr(gN)*Sqr(QH1p) - 2*Lambdax*Sqr(gN)*Sqr(QH2p) - 2*Lambdax*Sqr(gN)*
      Sqr(QSp) + 3*Lambdax*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(
      Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Lambdax*(Sqr(Ye00)
      + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) +
      Sqr(Ye21) + Sqr(Ye22)) + 3*Lambdax*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(
      Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dlead_log_one_loop_TLambdax_dMassG() const 
{
   
   const double Yd00 = model.get_Yd(0,0);
   const double Yd01 = model.get_Yd(0,1);
   const double Yd02 = model.get_Yd(0,2);
   const double Yd10 = model.get_Yd(1,0);
   const double Yd11 = model.get_Yd(1,1);
   const double Yd12 = model.get_Yd(1,2);
   const double Yd20 = model.get_Yd(2,0);
   const double Yd21 = model.get_Yd(2,1);
   const double Yd22 = model.get_Yd(2,2);
   const double Kappa00 = model.get_Kappa(0,0);
   const double Kappa01 = model.get_Kappa(0,1);
   const double Kappa02 = model.get_Kappa(0,2);
   const double Kappa10 = model.get_Kappa(1,0);
   const double Kappa11 = model.get_Kappa(1,1);
   const double Kappa12 = model.get_Kappa(1,2);
   const double Kappa20 = model.get_Kappa(2,0);
   const double Kappa21 = model.get_Kappa(2,1);
   const double Kappa22 = model.get_Kappa(2,2);
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
   const double g3 = model.get_g3();

   double deriv = 64*Lambdax*Sqr(g3)*Sqr(Kappa00) + 64*Lambdax*Sqr(g3)*Sqr(
      Kappa01) + 64*Lambdax*Sqr(g3)*Sqr(Kappa02) + 64*Lambdax*Sqr(g3)*Sqr(Kappa10)
      + 64*Lambdax*Sqr(g3)*Sqr(Kappa11) + 64*Lambdax*Sqr(g3)*Sqr(Kappa12) + 64*
      Lambdax*Sqr(g3)*Sqr(Kappa20) + 64*Lambdax*Sqr(g3)*Sqr(Kappa21) + 64*Lambdax*
      Sqr(g3)*Sqr(Kappa22) + 64*Lambdax*Sqr(g3)*Sqr(Yd00) + 64*Lambdax*Sqr(g3)*Sqr
      (Yd01) + 64*Lambdax*Sqr(g3)*Sqr(Yd02) + 64*Lambdax*Sqr(g3)*Sqr(Yd10) + 64*
      Lambdax*Sqr(g3)*Sqr(Yd11) + 64*Lambdax*Sqr(g3)*Sqr(Yd12) + 64*Lambdax*Sqr(g3
      )*Sqr(Yd20) + 64*Lambdax*Sqr(g3)*Sqr(Yd21) + 64*Lambdax*Sqr(g3)*Sqr(Yd22) +
      64*Lambdax*Sqr(g3)*Sqr(Yu00) + 64*Lambdax*Sqr(g3)*Sqr(Yu01) + 64*Lambdax*Sqr
      (g3)*Sqr(Yu02) + 64*Lambdax*Sqr(g3)*Sqr(Yu10) + 64*Lambdax*Sqr(g3)*Sqr(Yu11)
      + 64*Lambdax*Sqr(g3)*Sqr(Yu12) + 64*Lambdax*Sqr(g3)*Sqr(Yu20) + 64*Lambdax*
      Sqr(g3)*Sqr(Yu21) + 64*Lambdax*Sqr(g3)*Sqr(Yu22);


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dlead_log_one_loop_TLambdax_dMassBp() const 
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
   const double gN = model.get_gN();

   double deriv = 6*Lambdax*Sqr(Kappa00)*(4*Sqr(gN)*Sqr(QDxbarp) + 4*Sqr(gN)*
      Sqr(QDxp) + 4*Sqr(gN)*Sqr(QSp)) + 6*Lambdax*Sqr(Kappa01)*(4*Sqr(gN)*Sqr(
      QDxbarp) + 4*Sqr(gN)*Sqr(QDxp) + 4*Sqr(gN)*Sqr(QSp)) + 6*Lambdax*Sqr(Kappa02
      )*(4*Sqr(gN)*Sqr(QDxbarp) + 4*Sqr(gN)*Sqr(QDxp) + 4*Sqr(gN)*Sqr(QSp)) + 6*
      Lambdax*Sqr(Kappa10)*(4*Sqr(gN)*Sqr(QDxbarp) + 4*Sqr(gN)*Sqr(QDxp) + 4*Sqr(
      gN)*Sqr(QSp)) + 6*Lambdax*Sqr(Kappa11)*(4*Sqr(gN)*Sqr(QDxbarp) + 4*Sqr(gN)*
      Sqr(QDxp) + 4*Sqr(gN)*Sqr(QSp)) + 6*Lambdax*Sqr(Kappa12)*(4*Sqr(gN)*Sqr(
      QDxbarp) + 4*Sqr(gN)*Sqr(QDxp) + 4*Sqr(gN)*Sqr(QSp)) + 6*Lambdax*Sqr(Kappa20
      )*(4*Sqr(gN)*Sqr(QDxbarp) + 4*Sqr(gN)*Sqr(QDxp) + 4*Sqr(gN)*Sqr(QSp)) + 6*
      Lambdax*Sqr(Kappa21)*(4*Sqr(gN)*Sqr(QDxbarp) + 4*Sqr(gN)*Sqr(QDxp) + 4*Sqr(
      gN)*Sqr(QSp)) + 6*Lambdax*Sqr(Kappa22)*(4*Sqr(gN)*Sqr(QDxbarp) + 4*Sqr(gN)*
      Sqr(QDxp) + 4*Sqr(gN)*Sqr(QSp)) + 4*Lambdax*Sqr(Lambda1200)*(4*Sqr(gN)*Sqr(
      QH1p) + 4*Sqr(gN)*Sqr(QH2p) + 4*Sqr(gN)*Sqr(QSp)) + 4*Lambdax*Sqr(Lambda1201
      )*(4*Sqr(gN)*Sqr(QH1p) + 4*Sqr(gN)*Sqr(QH2p) + 4*Sqr(gN)*Sqr(QSp)) + 4*
      Lambdax*Sqr(Lambda1210)*(4*Sqr(gN)*Sqr(QH1p) + 4*Sqr(gN)*Sqr(QH2p) + 4*Sqr(
      gN)*Sqr(QSp)) + 4*Lambdax*Sqr(Lambda1211)*(4*Sqr(gN)*Sqr(QH1p) + 4*Sqr(gN)*
      Sqr(QH2p) + 4*Sqr(gN)*Sqr(QSp)) + Power(gN,3)*(8*gN*Lambdax*Sqr(QH1p) + 8*gN
      *Lambdax*Sqr(QH2p) + 8*gN*Lambdax*Sqr(QSp))*(9*Sqr(Qdp) + 9*Sqr(QDxbarp) + 9
      *Sqr(QDxp) + 3*Sqr(Qep) + 6*Sqr(QH1p) + 6*Sqr(QH2p) + 2*Sqr(QHpbarp) + 2*Sqr
      (QHpp) + 6*Sqr(QLp) + 18*Sqr(QQp) + 3*Sqr(QSp) + 9*Sqr(Qup)) + 2*Sqr(gN)*(4*
      Lambdax*Sqr(gN)*Sqr(QH1p) + 4*Lambdax*Sqr(gN)*Sqr(QH2p) + 4*Lambdax*Sqr(gN)*
      Sqr(QSp))*(9*Sqr(Qdp) + 9*Sqr(QDxbarp) + 9*Sqr(QDxp) + 3*Sqr(Qep) + 6*Sqr(
      QH1p) + 6*Sqr(QH2p) + 2*Sqr(QHpbarp) + 2*Sqr(QHpp) + 6*Sqr(QLp) + 18*Sqr(QQp
      ) + 3*Sqr(QSp) + 9*Sqr(Qup)) + 6*Lambdax*(4*Sqr(gN)*Sqr(Qdp) + 4*Sqr(gN)*Sqr
      (QH1p) + 4*Sqr(gN)*Sqr(QQp))*Sqr(Yd00) + 6*Lambdax*(4*Sqr(gN)*Sqr(Qdp) + 4*
      Sqr(gN)*Sqr(QH1p) + 4*Sqr(gN)*Sqr(QQp))*Sqr(Yd01) + 6*Lambdax*(4*Sqr(gN)*Sqr
      (Qdp) + 4*Sqr(gN)*Sqr(QH1p) + 4*Sqr(gN)*Sqr(QQp))*Sqr(Yd02) + 6*Lambdax*(4*
      Sqr(gN)*Sqr(Qdp) + 4*Sqr(gN)*Sqr(QH1p) + 4*Sqr(gN)*Sqr(QQp))*Sqr(Yd10) + 6*
      Lambdax*(4*Sqr(gN)*Sqr(Qdp) + 4*Sqr(gN)*Sqr(QH1p) + 4*Sqr(gN)*Sqr(QQp))*Sqr(
      Yd11) + 6*Lambdax*(4*Sqr(gN)*Sqr(Qdp) + 4*Sqr(gN)*Sqr(QH1p) + 4*Sqr(gN)*Sqr(
      QQp))*Sqr(Yd12) + 6*Lambdax*(4*Sqr(gN)*Sqr(Qdp) + 4*Sqr(gN)*Sqr(QH1p) + 4*
      Sqr(gN)*Sqr(QQp))*Sqr(Yd20) + 6*Lambdax*(4*Sqr(gN)*Sqr(Qdp) + 4*Sqr(gN)*Sqr(
      QH1p) + 4*Sqr(gN)*Sqr(QQp))*Sqr(Yd21) + 6*Lambdax*(4*Sqr(gN)*Sqr(Qdp) + 4*
      Sqr(gN)*Sqr(QH1p) + 4*Sqr(gN)*Sqr(QQp))*Sqr(Yd22) + 2*Lambdax*(4*Sqr(gN)*Sqr
      (Qep) + 4*Sqr(gN)*Sqr(QH1p) + 4*Sqr(gN)*Sqr(QLp))*Sqr(Ye00) + 2*Lambdax*(4*
      Sqr(gN)*Sqr(Qep) + 4*Sqr(gN)*Sqr(QH1p) + 4*Sqr(gN)*Sqr(QLp))*Sqr(Ye01) + 2*
      Lambdax*(4*Sqr(gN)*Sqr(Qep) + 4*Sqr(gN)*Sqr(QH1p) + 4*Sqr(gN)*Sqr(QLp))*Sqr(
      Ye02) + 2*Lambdax*(4*Sqr(gN)*Sqr(Qep) + 4*Sqr(gN)*Sqr(QH1p) + 4*Sqr(gN)*Sqr(
      QLp))*Sqr(Ye10) + 2*Lambdax*(4*Sqr(gN)*Sqr(Qep) + 4*Sqr(gN)*Sqr(QH1p) + 4*
      Sqr(gN)*Sqr(QLp))*Sqr(Ye11) + 2*Lambdax*(4*Sqr(gN)*Sqr(Qep) + 4*Sqr(gN)*Sqr(
      QH1p) + 4*Sqr(gN)*Sqr(QLp))*Sqr(Ye12) + 2*Lambdax*(4*Sqr(gN)*Sqr(Qep) + 4*
      Sqr(gN)*Sqr(QH1p) + 4*Sqr(gN)*Sqr(QLp))*Sqr(Ye20) + 2*Lambdax*(4*Sqr(gN)*Sqr
      (Qep) + 4*Sqr(gN)*Sqr(QH1p) + 4*Sqr(gN)*Sqr(QLp))*Sqr(Ye21) + 2*Lambdax*(4*
      Sqr(gN)*Sqr(Qep) + 4*Sqr(gN)*Sqr(QH1p) + 4*Sqr(gN)*Sqr(QLp))*Sqr(Ye22) + 6*
      Lambdax*(4*Sqr(gN)*Sqr(QH2p) + 4*Sqr(gN)*Sqr(QQp) + 4*Sqr(gN)*Sqr(Qup))*Sqr(
      Yu00) + 6*Lambdax*(4*Sqr(gN)*Sqr(QH2p) + 4*Sqr(gN)*Sqr(QQp) + 4*Sqr(gN)*Sqr(
      Qup))*Sqr(Yu01) + 6*Lambdax*(4*Sqr(gN)*Sqr(QH2p) + 4*Sqr(gN)*Sqr(QQp) + 4*
      Sqr(gN)*Sqr(Qup))*Sqr(Yu02) + 6*Lambdax*(4*Sqr(gN)*Sqr(QH2p) + 4*Sqr(gN)*Sqr
      (QQp) + 4*Sqr(gN)*Sqr(Qup))*Sqr(Yu10) + 6*Lambdax*(4*Sqr(gN)*Sqr(QH2p) + 4*
      Sqr(gN)*Sqr(QQp) + 4*Sqr(gN)*Sqr(Qup))*Sqr(Yu11) + 6*Lambdax*(4*Sqr(gN)*Sqr(
      QH2p) + 4*Sqr(gN)*Sqr(QQp) + 4*Sqr(gN)*Sqr(Qup))*Sqr(Yu12) + 6*Lambdax*(4*
      Sqr(gN)*Sqr(QH2p) + 4*Sqr(gN)*Sqr(QQp) + 4*Sqr(gN)*Sqr(Qup))*Sqr(Yu20) + 6*
      Lambdax*(4*Sqr(gN)*Sqr(QH2p) + 4*Sqr(gN)*Sqr(QQp) + 4*Sqr(gN)*Sqr(Qup))*Sqr(
      Yu21) + 6*Lambdax*(4*Sqr(gN)*Sqr(QH2p) + 4*Sqr(gN)*Sqr(QQp) + 4*Sqr(gN)*Sqr(
      Qup))*Sqr(Yu22) + (4*Lambdax*Sqr(gN)*Sqr(QH1p) + 4*Lambdax*Sqr(gN)*Sqr(QH2p)
      + 4*Lambdax*Sqr(gN)*Sqr(QSp))*(-0.6*Sqr(g1) - 3*Sqr(g2) + 3*(Sqr(Kappa00) +
      Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) +
      Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*(Sqr(Lambda1200) + Sqr(
      Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 12*Sqr(Lambdax) - 2*Sqr(
      gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QSp) + 3*(Sqr(Yd00) +
      Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(
      Yd21) + Sqr(Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(
      Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22) + 3*(Sqr(Yu00) + Sqr(
      Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21
      ) + Sqr(Yu22))) + (4*Sqr(gN)*Sqr(QH1p) + 4*Sqr(gN)*Sqr(QH2p) + 4*Sqr(gN)*Sqr
      (QSp))*(4*Power(Lambdax,3) - 0.6*Lambdax*Sqr(g1) - 3*Lambdax*Sqr(g2) + 3*
      Lambdax*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(
      Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*
      Lambdax*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(
      Lambda1211)) - 2*Lambdax*Sqr(gN)*Sqr(QH1p) - 2*Lambdax*Sqr(gN)*Sqr(QH2p) - 2
      *Lambdax*Sqr(gN)*Sqr(QSp) + 3*Lambdax*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) +
      Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) +
      Lambdax*(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(
      Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22)) + 3*Lambdax*(Sqr(Yu00) + Sqr(Yu01
      ) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) +
      Sqr(Yu22)));


   return twoLoop*deriv;
}

} // namespace flexiblesusy
