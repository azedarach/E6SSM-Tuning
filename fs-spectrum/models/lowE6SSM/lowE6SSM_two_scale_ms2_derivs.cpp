#include "lowE6SSM_two_scale_tuning_calculator.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

// one loop
double lowE6SSM_tuning_calculator::deriv_dbeta_one_loop_ms2_dLambdax() const
{

   const double Lambdax = model.get_Lambdax();

   const double ALambdax = get_ALambdax();

   const double mHd2 = model.get_mHd2();
   const double mHu2 = model.get_mHu2();
   const double ms2 = model.get_ms2();

   double deriv = 8*Lambdax*(mHd2 + mHu2 + ms2) + 8*Lambdax*Sqr(ALambdax);


   return oneOver16PiSqr*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_one_loop_ms2_dALambdax() const
{

   const double Lambdax = model.get_Lambdax();

   const double ALambdax = get_ALambdax();

   double deriv = 8*ALambdax*Sqr(Lambdax);


   return oneOver16PiSqr*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_one_loop_ms2_dAYu22() const
{

   double deriv = 0;


   return oneOver16PiSqr*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_one_loop_ms2_dmq222() const
{
   const lowE6SSM_input_parameters inputs = model.get_input();

   const auto QQp = inputs.QQp;
   const auto QSp = inputs.QSp;
   
   const double gN = model.get_gN();

   double deriv = 12*QQp*QSp*Sqr(gN);


   return oneOver16PiSqr*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_one_loop_ms2_dmHd2() const
{
   const lowE6SSM_input_parameters inputs = model.get_input();

   const auto QH1p = inputs.QH1p;
   const auto QSp = inputs.QSp;
   
   const double Lambdax = model.get_Lambdax();
   const double gN = model.get_gN();

   double deriv = 4*QH1p*QSp*Sqr(gN) + 4*Sqr(Lambdax);


   return oneOver16PiSqr*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_one_loop_ms2_dmHu2() const
{
   const lowE6SSM_input_parameters inputs = model.get_input();

   const auto QH2p = inputs.QH2p;
   const auto QSp = inputs.QSp;
   
   const double Lambdax = model.get_Lambdax();
   const double gN = model.get_gN();

   double deriv = 4*QH2p*QSp*Sqr(gN) + 4*Sqr(Lambdax);


   return oneOver16PiSqr*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_one_loop_ms2_dmu222() const
{
   const lowE6SSM_input_parameters inputs = model.get_input();

   const auto Qup = inputs.Qup;
   const auto QSp = inputs.QSp;
   
   const double gN = model.get_gN();

   double deriv = 6*QSp*Qup*Sqr(gN);


   return oneOver16PiSqr*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_one_loop_ms2_dms2() const
{
   const lowE6SSM_input_parameters inputs = model.get_input();

   const auto QSp = inputs.QSp;
   
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
   const double gN = model.get_gN();

   double deriv = 6*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10)
      + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22))
      + 4*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211))
      + 4*Sqr(Lambdax) + 2*Sqr(gN)*Sqr(QSp);


   return oneOver16PiSqr*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_one_loop_ms2_dMassB() const
{

   double deriv = 0;


   return oneOver16PiSqr*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_one_loop_ms2_dMassWB() const
{

   double deriv = 0;


   return oneOver16PiSqr*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_one_loop_ms2_dMassG() const
{

   double deriv = 0;


   return oneOver16PiSqr*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_one_loop_ms2_dMassBp() const
{
   const lowE6SSM_input_parameters inputs = model.get_input();

   const auto QSp = inputs.QSp;
   
   const double gN = model.get_gN();

   const double MassBp = model.get_MassBp();

   double deriv = -16*MassBp*Sqr(gN)*Sqr(QSp);


   return oneOver16PiSqr*deriv;
}

// two loop
double lowE6SSM_tuning_calculator::deriv_dbeta_two_loop_ms2_dLambdax() const
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
   const double AYu22 = get_AYu22();
   const double MassB = model.get_MassB();
   const double MassWB = model.get_MassWB();
   const double MassG = model.get_MassG();
   const double MassBp = model.get_MassBp();

   const double mq200 = model.get_mq2(0,0);
   const double mq201 = model.get_mq2(0,1);
   const double mq202 = model.get_mq2(0,2);
   const double mq210 = model.get_mq2(1,0);
   const double mq211 = model.get_mq2(1,1);
   const double mq212 = model.get_mq2(1,2);
   const double mq220 = model.get_mq2(2,0);
   const double mq221 = model.get_mq2(2,1);
   const double mq222 = model.get_mq2(2,2);
   const double ml200 = model.get_ml2(0,0);
   const double ml201 = model.get_ml2(0,1);
   const double ml202 = model.get_ml2(0,2);
   const double ml210 = model.get_ml2(1,0);
   const double ml211 = model.get_ml2(1,1);
   const double ml212 = model.get_ml2(1,2);
   const double ml220 = model.get_ml2(2,0);
   const double ml221 = model.get_ml2(2,1);
   const double ml222 = model.get_ml2(2,2);
   const double mHd2 = model.get_mHd2();
   const double mHu2 = model.get_mHu2();
   const double md200 = model.get_md2(0,0);
   const double md201 = model.get_md2(0,1);
   const double md202 = model.get_md2(0,2);
   const double md210 = model.get_md2(1,0);
   const double md211 = model.get_md2(1,1);
   const double md212 = model.get_md2(1,2);
   const double md220 = model.get_md2(2,0);
   const double md221 = model.get_md2(2,1);
   const double md222 = model.get_md2(2,2);
   const double mu200 = model.get_mu2(0,0);
   const double mu201 = model.get_mu2(0,1);
   const double mu202 = model.get_mu2(0,2);
   const double mu210 = model.get_mu2(1,0);
   const double mu211 = model.get_mu2(1,1);
   const double mu212 = model.get_mu2(1,2);
   const double mu220 = model.get_mu2(2,0);
   const double mu221 = model.get_mu2(2,1);
   const double mu222 = model.get_mu2(2,2);
   const double me200 = model.get_me2(0,0);
   const double me201 = model.get_me2(0,1);
   const double me202 = model.get_me2(0,2);
   const double me210 = model.get_me2(1,0);
   const double me211 = model.get_me2(1,1);
   const double me212 = model.get_me2(1,2);
   const double me220 = model.get_me2(2,0);
   const double me221 = model.get_me2(2,1);
   const double me222 = model.get_me2(2,2);
   const double ms2 = model.get_ms2();
   const double mH1I200 = model.get_mH1I2(0,0);
   const double mH1I201 = model.get_mH1I2(0,1);
   const double mH1I210 = model.get_mH1I2(1,0);
   const double mH1I211 = model.get_mH1I2(1,1);
   const double mH2I200 = model.get_mH2I2(0,0);
   const double mH2I201 = model.get_mH2I2(0,1);
   const double mH2I210 = model.get_mH2I2(1,0);
   const double mH2I211 = model.get_mH2I2(1,1);
   const double msI200 = model.get_msI2(0,0);
   const double msI201 = model.get_msI2(0,1);
   const double msI210 = model.get_msI2(1,0);
   const double msI211 = model.get_msI2(1,1);
   const double mDx200 = model.get_mDx2(0,0);
   const double mDx201 = model.get_mDx2(0,1);
   const double mDx202 = model.get_mDx2(0,2);
   const double mDx210 = model.get_mDx2(1,0);
   const double mDx211 = model.get_mDx2(1,1);
   const double mDx212 = model.get_mDx2(1,2);
   const double mDx220 = model.get_mDx2(2,0);
   const double mDx221 = model.get_mDx2(2,1);
   const double mDx222 = model.get_mDx2(2,2);
   const double mDxbar200 = model.get_mDxbar2(0,0);
   const double mDxbar201 = model.get_mDxbar2(0,1);
   const double mDxbar202 = model.get_mDxbar2(0,2);
   const double mDxbar210 = model.get_mDxbar2(1,0);
   const double mDxbar211 = model.get_mDxbar2(1,1);
   const double mDxbar212 = model.get_mDxbar2(1,2);
   const double mDxbar220 = model.get_mDxbar2(2,0);
   const double mDxbar221 = model.get_mDxbar2(2,1);
   const double mDxbar222 = model.get_mDxbar2(2,2);
   const double mHp2 = model.get_mHp2();
   const double mHpbar2 = model.get_mHpbar2();

   double deriv = -0.8*(80*Power(Lambdax,3)*(mHd2 + mHu2 + ms2) + 45*ALambdax*
      Lambdax*(TYd00*Yd00 + TYd01*Yd01 + TYd02*Yd02 + TYd10*Yd10 + TYd11*Yd11 +
      TYd12*Yd12 + TYd20*Yd20 + TYd21*Yd21 + TYd22*Yd22) + 15*Lambdax*(Yd00*(md200
      *Yd00 + md201*Yd10 + md202*Yd20) + Yd10*(md210*Yd00 + md211*Yd10 + md212*
      Yd20) + Yd20*(md220*Yd00 + md221*Yd10 + md222*Yd20) + Yd01*(md200*Yd01 +
      md201*Yd11 + md202*Yd21) + Yd11*(md210*Yd01 + md211*Yd11 + md212*Yd21) +
      Yd21*(md220*Yd01 + md221*Yd11 + md222*Yd21) + Yd02*(md200*Yd02 + md201*Yd12
      + md202*Yd22) + Yd12*(md210*Yd02 + md211*Yd12 + md212*Yd22) + Yd22*(md220*
      Yd02 + md221*Yd12 + md222*Yd22)) + 15*Lambdax*(Yd00*(mq200*Yd00 + mq201*Yd01
      + mq202*Yd02) + Yd01*(mq210*Yd00 + mq211*Yd01 + mq212*Yd02) + Yd02*(mq220*
      Yd00 + mq221*Yd01 + mq222*Yd02) + Yd10*(mq200*Yd10 + mq201*Yd11 + mq202*Yd12
      ) + Yd11*(mq210*Yd10 + mq211*Yd11 + mq212*Yd12) + Yd12*(mq220*Yd10 + mq221*
      Yd11 + mq222*Yd12) + Yd20*(mq200*Yd20 + mq201*Yd21 + mq202*Yd22) + Yd21*(
      mq210*Yd20 + mq211*Yd21 + mq212*Yd22) + Yd22*(mq220*Yd20 + mq221*Yd21 +
      mq222*Yd22)) + 15*ALambdax*Lambdax*(TYe00*Ye00 + TYe01*Ye01 + TYe02*Ye02 +
      TYe10*Ye10 + TYe11*Ye11 + TYe12*Ye12 + TYe20*Ye20 + TYe21*Ye21 + TYe22*Ye22)
      + 5*Lambdax*(Ye00*(me200*Ye00 + me201*Ye10 + me202*Ye20) + Ye10*(me210*Ye00
      + me211*Ye10 + me212*Ye20) + Ye20*(me220*Ye00 + me221*Ye10 + me222*Ye20) +
      Ye01*(me200*Ye01 + me201*Ye11 + me202*Ye21) + Ye11*(me210*Ye01 + me211*Ye11
      + me212*Ye21) + Ye21*(me220*Ye01 + me221*Ye11 + me222*Ye21) + Ye02*(me200*
      Ye02 + me201*Ye12 + me202*Ye22) + Ye12*(me210*Ye02 + me211*Ye12 + me212*Ye22
      ) + Ye22*(me220*Ye02 + me221*Ye12 + me222*Ye22)) + 5*Lambdax*(Ye00*(ml200*
      Ye00 + ml201*Ye01 + ml202*Ye02) + Ye01*(ml210*Ye00 + ml211*Ye01 + ml212*Ye02
      ) + Ye02*(ml220*Ye00 + ml221*Ye01 + ml222*Ye02) + Ye10*(ml200*Ye10 + ml201*
      Ye11 + ml202*Ye12) + Ye11*(ml210*Ye10 + ml211*Ye11 + ml212*Ye12) + Ye12*(
      ml220*Ye10 + ml221*Ye11 + ml222*Ye12) + Ye20*(ml200*Ye20 + ml201*Ye21 +
      ml202*Ye22) + Ye21*(ml210*Ye20 + ml211*Ye21 + ml212*Ye22) + Ye22*(ml220*Ye20
      + ml221*Ye21 + ml222*Ye22)) + 45*ALambdax*Lambdax*(TYu00*Yu00 + TYu01*Yu01
      + TYu02*Yu02 + TYu10*Yu10 + TYu11*Yu11 + TYu12*Yu12 + TYu20*Yu20 + TYu21*
      Yu21 + TYu22*Yu22) + 15*Lambdax*(Yu00*(mq200*Yu00 + mq201*Yu01 + mq202*Yu02)
      + Yu01*(mq210*Yu00 + mq211*Yu01 + mq212*Yu02) + Yu02*(mq220*Yu00 + mq221*
      Yu01 + mq222*Yu02) + Yu10*(mq200*Yu10 + mq201*Yu11 + mq202*Yu12) + Yu11*(
      mq210*Yu10 + mq211*Yu11 + mq212*Yu12) + Yu12*(mq220*Yu10 + mq221*Yu11 +
      mq222*Yu12) + Yu20*(mq200*Yu20 + mq201*Yu21 + mq202*Yu22) + Yu21*(mq210*Yu20
      + mq211*Yu21 + mq212*Yu22) + Yu22*(mq220*Yu20 + mq221*Yu21 + mq222*Yu22)) +
      15*Lambdax*(Yu00*(mu200*Yu00 + mu201*Yu10 + mu202*Yu20) + Yu10*(mu210*Yu00
      + mu211*Yu10 + mu212*Yu20) + Yu20*(mu220*Yu00 + mu221*Yu10 + mu222*Yu20) +
      Yu01*(mu200*Yu01 + mu201*Yu11 + mu202*Yu21) + Yu11*(mu210*Yu01 + mu211*Yu11
      + mu212*Yu21) + Yu21*(mu220*Yu01 + mu221*Yu11 + mu222*Yu21) + Yu02*(mu200*
      Yu02 + mu201*Yu12 + mu202*Yu22) + Yu12*(mu210*Yu02 + mu211*Yu12 + mu212*Yu22
      ) + Yu22*(mu220*Yu02 + mu221*Yu12 + mu222*Yu22)) + 40*Power(Lambdax,3)*Sqr(
      ALambdax) + 6*ALambdax*Lambdax*MassB*Sqr(g1) + 3*MassB*(ALambdax*Lambdax - 2
      *Lambdax*MassB)*Sqr(g1) - 3*Lambdax*mHd2*Sqr(g1) - 3*Lambdax*mHu2*Sqr(g1) -
      3*Lambdax*ms2*Sqr(g1) - 6*Lambdax*Sqr(ALambdax)*Sqr(g1) + 30*ALambdax*
      Lambdax*MassWB*Sqr(g2) + 15*MassWB*(ALambdax*Lambdax - 2*Lambdax*MassWB)*Sqr
      (g2) - 15*Lambdax*mHd2*Sqr(g2) - 15*Lambdax*mHu2*Sqr(g2) - 15*Lambdax*ms2*
      Sqr(g2) - 30*Lambdax*Sqr(ALambdax)*Sqr(g2) + 20*Lambdax*QSp*(mHd2*QH1p +
      mHu2*QH2p + ms2*QSp)*Sqr(gN) + 20*ALambdax*Lambdax*MassBp*Sqr(gN)*Sqr(QH1p)
      - 10*Lambdax*mHd2*Sqr(gN)*Sqr(QH1p) - 10*Lambdax*mHu2*Sqr(gN)*Sqr(QH1p) - 10
      *Lambdax*ms2*Sqr(gN)*Sqr(QH1p) - 20*Lambdax*Sqr(ALambdax)*Sqr(gN)*Sqr(QH1p)
      + 20*ALambdax*Lambdax*MassBp*Sqr(gN)*Sqr(QH2p) - 10*Lambdax*mHd2*Sqr(gN)*Sqr
      (QH2p) - 10*Lambdax*mHu2*Sqr(gN)*Sqr(QH2p) - 10*Lambdax*ms2*Sqr(gN)*Sqr(QH2p
      ) - 20*Lambdax*Sqr(ALambdax)*Sqr(gN)*Sqr(QH2p) - 5*MassBp*Sqr(gN)*(2*Lambdax
      *(-ALambdax + 2*MassBp)*(Sqr(QH1p) + Sqr(QH2p) - Sqr(QSp)) + 2*(-(ALambdax*
      Lambdax) + 2*Lambdax*MassBp)*(Sqr(QH1p) + Sqr(QH2p) - Sqr(QSp))) - 20*
      ALambdax*Lambdax*MassBp*Sqr(gN)*Sqr(QSp) + 10*Lambdax*mHd2*Sqr(gN)*Sqr(QSp)
      + 10*Lambdax*mHu2*Sqr(gN)*Sqr(QSp) + 10*Lambdax*ms2*Sqr(gN)*Sqr(QSp) + 20*
      Lambdax*Sqr(ALambdax)*Sqr(gN)*Sqr(QSp) + 15*Lambdax*(Sqr(TYd00) + Sqr(TYd01)
      + Sqr(TYd02) + Sqr(TYd10) + Sqr(TYd11) + Sqr(TYd12) + Sqr(TYd20) + Sqr(
      TYd21) + Sqr(TYd22)) + 5*Lambdax*(Sqr(TYe00) + Sqr(TYe01) + Sqr(TYe02) + Sqr
      (TYe10) + Sqr(TYe11) + Sqr(TYe12) + Sqr(TYe20) + Sqr(TYe21) + Sqr(TYe22)) +
      15*Lambdax*(Sqr(TYu00) + Sqr(TYu01) + Sqr(TYu02) + Sqr(TYu10) + Sqr(TYu11) +
      Sqr(TYu12) + Sqr(TYu20) + Sqr(TYu21) + Sqr(TYu22)) + 30*Lambdax*mHd2*(Sqr(
      Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20
      ) + Sqr(Yd21) + Sqr(Yd22)) + 15*Lambdax*mHu2*(Sqr(Yd00) + Sqr(Yd01) + Sqr(
      Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22
      )) + 15*Lambdax*ms2*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(
      Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + 30*Lambdax*Sqr(
      ALambdax)*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(
      Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + 10*Lambdax*mHd2*(Sqr(Ye00) +
      Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(
      Ye21) + Sqr(Ye22)) + 5*Lambdax*mHu2*(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr
      (Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22)) + 5*
      Lambdax*ms2*(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr
      (Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22)) + 10*Lambdax*Sqr(ALambdax)*(Sqr(
      Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20
      ) + Sqr(Ye21) + Sqr(Ye22)) + 15*Lambdax*mHd2*(Sqr(Yu00) + Sqr(Yu01) + Sqr(
      Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22
      )) + 30*Lambdax*mHu2*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(
      Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)) + 15*Lambdax*ms2*(Sqr
      (Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(
      Yu20) + Sqr(Yu21) + Sqr(Yu22)) + 30*Lambdax*Sqr(ALambdax)*(Sqr(Yu00) + Sqr(
      Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21
      ) + Sqr(Yu22)) + Lambdax*(15*ALambdax*(TYd00*Yd00 + TYd01*Yd01 + TYd02*Yd02
      + TYd10*Yd10 + TYd11*Yd11 + TYd12*Yd12 + TYd20*Yd20 + TYd21*Yd21 + TYd22*
      Yd22) + 15*(Yd00*(md200*Yd00 + md201*Yd10 + md202*Yd20) + Yd10*(md210*Yd00 +
      md211*Yd10 + md212*Yd20) + Yd20*(md220*Yd00 + md221*Yd10 + md222*Yd20) +
      Yd01*(md200*Yd01 + md201*Yd11 + md202*Yd21) + Yd11*(md210*Yd01 + md211*Yd11
      + md212*Yd21) + Yd21*(md220*Yd01 + md221*Yd11 + md222*Yd21) + Yd02*(md200*
      Yd02 + md201*Yd12 + md202*Yd22) + Yd12*(md210*Yd02 + md211*Yd12 + md212*Yd22
      ) + Yd22*(md220*Yd02 + md221*Yd12 + md222*Yd22)) + 15*(Yd00*(mq200*Yd00 +
      mq201*Yd01 + mq202*Yd02) + Yd01*(mq210*Yd00 + mq211*Yd01 + mq212*Yd02) +
      Yd02*(mq220*Yd00 + mq221*Yd01 + mq222*Yd02) + Yd10*(mq200*Yd10 + mq201*Yd11
      + mq202*Yd12) + Yd11*(mq210*Yd10 + mq211*Yd11 + mq212*Yd12) + Yd12*(mq220*
      Yd10 + mq221*Yd11 + mq222*Yd12) + Yd20*(mq200*Yd20 + mq201*Yd21 + mq202*Yd22
      ) + Yd21*(mq210*Yd20 + mq211*Yd21 + mq212*Yd22) + Yd22*(mq220*Yd20 + mq221*
      Yd21 + mq222*Yd22)) + 5*ALambdax*(TYe00*Ye00 + TYe01*Ye01 + TYe02*Ye02 +
      TYe10*Ye10 + TYe11*Ye11 + TYe12*Ye12 + TYe20*Ye20 + TYe21*Ye21 + TYe22*Ye22)
      + 5*(Ye00*(me200*Ye00 + me201*Ye10 + me202*Ye20) + Ye10*(me210*Ye00 + me211
      *Ye10 + me212*Ye20) + Ye20*(me220*Ye00 + me221*Ye10 + me222*Ye20) + Ye01*(
      me200*Ye01 + me201*Ye11 + me202*Ye21) + Ye11*(me210*Ye01 + me211*Ye11 +
      me212*Ye21) + Ye21*(me220*Ye01 + me221*Ye11 + me222*Ye21) + Ye02*(me200*Ye02
      + me201*Ye12 + me202*Ye22) + Ye12*(me210*Ye02 + me211*Ye12 + me212*Ye22) +
      Ye22*(me220*Ye02 + me221*Ye12 + me222*Ye22)) + 5*(Ye00*(ml200*Ye00 + ml201*
      Ye01 + ml202*Ye02) + Ye01*(ml210*Ye00 + ml211*Ye01 + ml212*Ye02) + Ye02*(
      ml220*Ye00 + ml221*Ye01 + ml222*Ye02) + Ye10*(ml200*Ye10 + ml201*Ye11 +
      ml202*Ye12) + Ye11*(ml210*Ye10 + ml211*Ye11 + ml212*Ye12) + Ye12*(ml220*Ye10
      + ml221*Ye11 + ml222*Ye12) + Ye20*(ml200*Ye20 + ml201*Ye21 + ml202*Ye22) +
      Ye21*(ml210*Ye20 + ml211*Ye21 + ml212*Ye22) + Ye22*(ml220*Ye20 + ml221*Ye21
      + ml222*Ye22)) + 15*ALambdax*(TYu00*Yu00 + TYu01*Yu01 + TYu02*Yu02 + TYu10*
      Yu10 + TYu11*Yu11 + TYu12*Yu12 + TYu20*Yu20 + TYu21*Yu21 + TYu22*Yu22) + 15*
      (Yu00*(mq200*Yu00 + mq201*Yu01 + mq202*Yu02) + Yu01*(mq210*Yu00 + mq211*Yu01
      + mq212*Yu02) + Yu02*(mq220*Yu00 + mq221*Yu01 + mq222*Yu02) + Yu10*(mq200*
      Yu10 + mq201*Yu11 + mq202*Yu12) + Yu11*(mq210*Yu10 + mq211*Yu11 + mq212*Yu12
      ) + Yu12*(mq220*Yu10 + mq221*Yu11 + mq222*Yu12) + Yu20*(mq200*Yu20 + mq201*
      Yu21 + mq202*Yu22) + Yu21*(mq210*Yu20 + mq211*Yu21 + mq212*Yu22) + Yu22*(
      mq220*Yu20 + mq221*Yu21 + mq222*Yu22)) + 15*(Yu00*(mu200*Yu00 + mu201*Yu10 +
      mu202*Yu20) + Yu10*(mu210*Yu00 + mu211*Yu10 + mu212*Yu20) + Yu20*(mu220*
      Yu00 + mu221*Yu10 + mu222*Yu20) + Yu01*(mu200*Yu01 + mu201*Yu11 + mu202*Yu21
      ) + Yu11*(mu210*Yu01 + mu211*Yu11 + mu212*Yu21) + Yu21*(mu220*Yu01 + mu221*
      Yu11 + mu222*Yu21) + Yu02*(mu200*Yu02 + mu201*Yu12 + mu202*Yu22) + Yu12*(
      mu210*Yu02 + mu211*Yu12 + mu212*Yu22) + Yu22*(mu220*Yu02 + mu221*Yu12 +
      mu222*Yu22)) + 3*(ALambdax - 2*MassB)*MassB*Sqr(g1) - 3*mHd2*Sqr(g1) - 3*
      mHu2*Sqr(g1) - 3*ms2*Sqr(g1) + 15*(ALambdax - 2*MassWB)*MassWB*Sqr(g2) - 15*
      mHd2*Sqr(g2) - 15*mHu2*Sqr(g2) - 15*ms2*Sqr(g2) + 120*Sqr(ALambdax)*Sqr(
      Lambdax) - 10*mHd2*Sqr(gN)*Sqr(QH1p) - 10*mHu2*Sqr(gN)*Sqr(QH1p) - 10*ms2*
      Sqr(gN)*Sqr(QH1p) - 10*mHd2*Sqr(gN)*Sqr(QH2p) - 10*mHu2*Sqr(gN)*Sqr(QH2p) -
      10*ms2*Sqr(gN)*Sqr(QH2p) + 10*mHd2*Sqr(gN)*Sqr(QSp) + 10*mHu2*Sqr(gN)*Sqr(
      QSp) + 10*ms2*Sqr(gN)*Sqr(QSp) + 15*(Sqr(TYd00) + Sqr(TYd01) + Sqr(TYd02) +
      Sqr(TYd10) + Sqr(TYd11) + Sqr(TYd12) + Sqr(TYd20) + Sqr(TYd21) + Sqr(TYd22))
      + 5*(Sqr(TYe00) + Sqr(TYe01) + Sqr(TYe02) + Sqr(TYe10) + Sqr(TYe11) + Sqr(
      TYe12) + Sqr(TYe20) + Sqr(TYe21) + Sqr(TYe22)) + 15*(Sqr(TYu00) + Sqr(TYu01)
      + Sqr(TYu02) + Sqr(TYu10) + Sqr(TYu11) + Sqr(TYu12) + Sqr(TYu20) + Sqr(
      TYu21) + Sqr(TYu22)) + 30*mHd2*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10
      ) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + 15*mHu2*(
      Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(
      Yd20) + Sqr(Yd21) + Sqr(Yd22)) + 15*ms2*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) +
      Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + 10
      *mHd2*(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12)
      + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22)) + 5*mHu2*(Sqr(Ye00) + Sqr(Ye01) + Sqr(
      Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22
      )) + 5*ms2*(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(
      Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22)) + 15*mHd2*(Sqr(Yu00) + Sqr(Yu01)
      + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) +
      Sqr(Yu22)) + 30*mHu2*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(
      Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)) + 15*ms2*(Sqr(Yu00) +
      Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr
      (Yu21) + Sqr(Yu22))));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_two_loop_ms2_dALambdax() const
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
   const double AYu22 = get_AYu22();
   const double MassB = model.get_MassB();
   const double MassWB = model.get_MassWB();
   const double MassG = model.get_MassG();
   const double MassBp = model.get_MassBp();

   double deriv = -0.8*(Lambdax*(80*ALambdax*Power(Lambdax,3) + 15*Lambdax*(
      TYd00*Yd00 + TYd01*Yd01 + TYd02*Yd02 + TYd10*Yd10 + TYd11*Yd11 + TYd12*Yd12
      + TYd20*Yd20 + TYd21*Yd21 + TYd22*Yd22) + 5*Lambdax*(TYe00*Ye00 + TYe01*Ye01
      + TYe02*Ye02 + TYe10*Ye10 + TYe11*Ye11 + TYe12*Ye12 + TYe20*Ye20 + TYe21*
      Ye21 + TYe22*Ye22) + 15*Lambdax*(TYu00*Yu00 + TYu01*Yu01 + TYu02*Yu02 +
      TYu10*Yu10 + TYu11*Yu11 + TYu12*Yu12 + TYu20*Yu20 + TYu21*Yu21 + TYu22*Yu22)
      + 3*Lambdax*MassB*Sqr(g1) + 15*Lambdax*MassWB*Sqr(g2)) + 15*(TYd00*Yd00 +
      TYd01*Yd01 + TYd02*Yd02 + TYd10*Yd10 + TYd11*Yd11 + TYd12*Yd12 + TYd20*Yd20
      + TYd21*Yd21 + TYd22*Yd22)*Sqr(Lambdax) + 5*(TYe00*Ye00 + TYe01*Ye01 + TYe02
      *Ye02 + TYe10*Ye10 + TYe11*Ye11 + TYe12*Ye12 + TYe20*Ye20 + TYe21*Ye21 +
      TYe22*Ye22)*Sqr(Lambdax) + 15*(TYu00*Yu00 + TYu01*Yu01 + TYu02*Yu02 + TYu10*
      Yu10 + TYu11*Yu11 + TYu12*Yu12 + TYu20*Yu20 + TYu21*Yu21 + TYu22*Yu22)*Sqr(
      Lambdax) - 6*ALambdax*Sqr(g1)*Sqr(Lambdax) + 3*MassB*Sqr(g1)*Sqr(Lambdax) -
      30*ALambdax*Sqr(g2)*Sqr(Lambdax) + 15*MassWB*Sqr(g2)*Sqr(Lambdax) - 20*
      ALambdax*Sqr(gN)*Sqr(Lambdax)*Sqr(QH1p) + 10*MassBp*Sqr(gN)*Sqr(Lambdax)*Sqr
      (QH1p) - 20*ALambdax*Sqr(gN)*Sqr(Lambdax)*Sqr(QH2p) + 10*MassBp*Sqr(gN)*Sqr(
      Lambdax)*Sqr(QH2p) + 10*MassBp*Sqr(gN)*Sqr(Lambdax)*(Sqr(QH1p) + Sqr(QH2p) -
      Sqr(QSp)) + 20*ALambdax*Sqr(gN)*Sqr(Lambdax)*Sqr(QSp) - 10*MassBp*Sqr(gN)*
      Sqr(Lambdax)*Sqr(QSp) + 30*ALambdax*Sqr(Lambdax)*(Sqr(Yd00) + Sqr(Yd01) +
      Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(
      Yd22)) + 10*ALambdax*Sqr(Lambdax)*(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(
      Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22)) + 30*
      ALambdax*Sqr(Lambdax)*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(
      Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_two_loop_ms2_dAYu22() const
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
   const double AYu22 = get_AYu22();
   const double MassB = model.get_MassB();
   const double MassWB = model.get_MassWB();
   const double MassG = model.get_MassG();
   const double MassBp = model.get_MassBp();

   double deriv = -0.8*(15*Lambdax*TLambdax*Sqr(Yu22) + Lambdax*(30*AYu22*
      Lambdax*Sqr(Yu22) + 15*TLambdax*Sqr(Yu22)));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_two_loop_ms2_dmq222() const
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

   double deriv = -0.8*(-60*Power(gN,4)*Sqr(QQp)*Sqr(QSp) + Lambdax*(15*Lambdax
      *(Sqr(Yd02) + Sqr(Yd12) + Sqr(Yd22)) + 15*Lambdax*(Sqr(Yu02) + Sqr(Yu12) +
      Sqr(Yu22))) - QSp*Sqr(gN)*(QQp*Sqr(g1) + 45*QQp*Sqr(g2) + 80*QQp*Sqr(g3) +
      60*Power(QQp,3)*Sqr(gN) - 30*QQp*(Sqr(Yd02) + Sqr(Yd12) + Sqr(Yd22)) - 30*
      QQp*(Sqr(Yu02) + Sqr(Yu12) + Sqr(Yu22))));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_two_loop_ms2_dmHd2() const
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

   double deriv = -0.8*(20*Power(Lambdax,4) - 20*Power(gN,4)*Sqr(QH1p)*Sqr(QSp)
      - QSp*Sqr(gN)*(3*QH1p*Sqr(g1) + 15*QH1p*Sqr(g2) + 20*Power(QH1p,3)*Sqr(gN)
      - 10*QH1p*Sqr(Lambdax) - 30*QH1p*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(
      Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) - 10*QH1p
      *(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) +
      Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22))) + Lambdax*(-3*Lambdax*Sqr(g1) - 15*
      Lambdax*Sqr(g2) - 10*Lambdax*Sqr(gN)*Sqr(QH1p) - 10*Lambdax*Sqr(gN)*Sqr(QH2p
      ) + 10*Lambdax*Sqr(gN)*Sqr(QSp) + 30*Lambdax*(Sqr(Yd00) + Sqr(Yd01) + Sqr(
      Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22
      )) + 10*Lambdax*(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) +
      Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22)) + 15*Lambdax*(Sqr(Yu00) +
      Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(
      Yu21) + Sqr(Yu22))));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_two_loop_ms2_dmHu2() const
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

   double deriv = -0.8*(20*Power(Lambdax,4) - 20*Power(gN,4)*Sqr(QH2p)*Sqr(QSp)
      + Lambdax*(-3*Lambdax*Sqr(g1) - 15*Lambdax*Sqr(g2) - 10*Lambdax*Sqr(gN)*Sqr
      (QH1p) - 10*Lambdax*Sqr(gN)*Sqr(QH2p) + 10*Lambdax*Sqr(gN)*Sqr(QSp) + 15*
      Lambdax*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(
      Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + 5*Lambdax*(Sqr(Ye00) + Sqr(Ye01
      ) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) +
      Sqr(Ye22)) + 30*Lambdax*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr
      (Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) - QSp*Sqr(gN)*(3*
      QH2p*Sqr(g1) + 15*QH2p*Sqr(g2) + 20*Power(QH2p,3)*Sqr(gN) - 10*QH2p*Sqr(
      Lambdax) - 30*QH2p*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11
      ) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_two_loop_ms2_dmu222() const
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

   double deriv = -0.8*(-30*Power(gN,4)*Sqr(QSp)*Sqr(Qup) + 15*Sqr(Lambdax)*(
      Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)) - QSp*Sqr(gN)*(8*Qup*Sqr(g1) + 40*Qup*Sqr
      (g3) + 30*Power(Qup,3)*Sqr(gN) - 30*Qup*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))
      );


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_two_loop_ms2_dms2() const
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

   double deriv = -0.8*(20*Power(Lambdax,4) - 10*Power(gN,4)*Power(QSp,4) - 2*
      Sqr(g1)*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(
      Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) - 40*
      Sqr(g3)*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(
      Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 30*(
      Kappa00*(Kappa10*(Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12) +
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
      (Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)))) - 3*Sqr(g1)*(Sqr(Lambda1200) +
      Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) - 15*Sqr(g2)*(Sqr(
      Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 20*(
      Lambda1200*(Lambda1210*(Lambda1200*Lambda1210 + Lambda1201*Lambda1211) +
      Lambda1200*(Sqr(Lambda1200) + Sqr(Lambda1201))) + Lambda1201*(Lambda1211*(
      Lambda1200*Lambda1210 + Lambda1201*Lambda1211) + Lambda1201*(Sqr(Lambda1200)
      + Sqr(Lambda1201))) + Lambda1210*(Lambda1200*(Lambda1200*Lambda1210 +
      Lambda1201*Lambda1211) + Lambda1210*(Sqr(Lambda1210) + Sqr(Lambda1211))) +
      Lambda1211*(Lambda1201*(Lambda1200*Lambda1210 + Lambda1201*Lambda1211) +
      Lambda1211*(Sqr(Lambda1210) + Sqr(Lambda1211)))) - QSp*Sqr(gN)*(10*Power(QSp
      ,3)*Sqr(gN) - 15*QSp*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(
      Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(
      Kappa22)) - 10*QSp*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) +
      Sqr(Lambda1211)) - 10*QSp*Sqr(Lambdax)) - 15*Sqr(gN)*(Sqr(Kappa00) + Sqr(
      Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(
      Kappa20) + Sqr(Kappa21) + Sqr(Kappa22))*Sqr(QDxbarp) - 15*Sqr(gN)*(Sqr(
      Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(
      Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22))*Sqr(QDxp) - 10*Sqr(gN
      )*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211))*
      Sqr(QH1p) - 10*Sqr(gN)*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210)
      + Sqr(Lambda1211))*Sqr(QH2p) + 15*Sqr(gN)*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr
      (Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(
      Kappa21) + Sqr(Kappa22))*Sqr(QSp) + 10*Sqr(gN)*(Sqr(Lambda1200) + Sqr(
      Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211))*Sqr(QSp) + Lambdax*(-3*
      Lambdax*Sqr(g1) - 15*Lambdax*Sqr(g2) - 10*Lambdax*Sqr(gN)*Sqr(QH1p) - 10*
      Lambdax*Sqr(gN)*Sqr(QH2p) + 10*Lambdax*Sqr(gN)*Sqr(QSp) + 15*Lambdax*(Sqr(
      Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20
      ) + Sqr(Yd21) + Sqr(Yd22)) + 5*Lambdax*(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) +
      Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22)) + 15*
      Lambdax*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(
      Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_two_loop_ms2_dMassB() const
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

   double deriv = -0.8*(4*(Kappa00*TKappa00 + Kappa01*TKappa01 + Kappa02*
      TKappa02 + Kappa10*TKappa10 + Kappa11*TKappa11 + Kappa12*TKappa12 + Kappa20*
      TKappa20 + Kappa21*TKappa21 + Kappa22*TKappa22)*Sqr(g1) + 6*(Lambda1200*
      TLambda1200 + Lambda1201*TLambda1201 + Lambda1210*TLambda1210 + Lambda1211*
      TLambda1211)*Sqr(g1) + 3*Lambdax*TLambdax*Sqr(g1) + Lambdax*(-6*Lambdax*
      MassB*Sqr(g1) + 3*(-2*Lambdax*MassB + TLambdax)*Sqr(g1)) - 8*MassB*Sqr(g1)*(
      Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) +
      Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) - 12*MassB*Sqr(g1
      )*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_two_loop_ms2_dMassWB() const
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

   double deriv = -0.8*(30*(Lambda1200*TLambda1200 + Lambda1201*TLambda1201 +
      Lambda1210*TLambda1210 + Lambda1211*TLambda1211)*Sqr(g2) + 15*Lambdax*
      TLambdax*Sqr(g2) + Lambdax*(-30*Lambdax*MassWB*Sqr(g2) + 15*(-2*Lambdax*
      MassWB + TLambdax)*Sqr(g2)) - 60*MassWB*Sqr(g2)*(Sqr(Lambda1200) + Sqr(
      Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_two_loop_ms2_dMassG() const
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

   double deriv = -0.8*(80*(Kappa00*TKappa00 + Kappa01*TKappa01 + Kappa02*
      TKappa02 + Kappa10*TKappa10 + Kappa11*TKappa11 + Kappa12*TKappa12 + Kappa20*
      TKappa20 + Kappa21*TKappa21 + Kappa22*TKappa22)*Sqr(g3) - 160*MassG*Sqr(g3)*
      (Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) +
      Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_two_loop_ms2_dMassBp() const
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

   double deriv = -0.8*(15*(Kappa00*TKappa00 + Kappa01*TKappa01 + Kappa02*
      TKappa02 + Kappa10*TKappa10 + Kappa11*TKappa11 + Kappa12*TKappa12 + Kappa20*
      TKappa20 + Kappa21*TKappa21 + Kappa22*TKappa22)*Sqr(gN)*Sqr(QDxbarp) + 15*(
      Kappa00*TKappa00 + Kappa01*TKappa01 + Kappa02*TKappa02 + Kappa10*TKappa10 +
      Kappa11*TKappa11 + Kappa12*TKappa12 + Kappa20*TKappa20 + Kappa21*TKappa21 +
      Kappa22*TKappa22)*Sqr(gN)*Sqr(QDxp) + 10*(Lambda1200*TLambda1200 +
      Lambda1201*TLambda1201 + Lambda1210*TLambda1210 + Lambda1211*TLambda1211)*
      Sqr(gN)*Sqr(QH1p) + 10*Lambdax*TLambdax*Sqr(gN)*Sqr(QH1p) + 10*(Lambda1200*
      TLambda1200 + Lambda1201*TLambda1201 + Lambda1210*TLambda1210 + Lambda1211*
      TLambda1211)*Sqr(gN)*Sqr(QH2p) + 10*Lambdax*TLambdax*Sqr(gN)*Sqr(QH2p) - 15*
      (Kappa00*TKappa00 + Kappa01*TKappa01 + Kappa02*TKappa02 + Kappa10*TKappa10 +
      Kappa11*TKappa11 + Kappa12*TKappa12 + Kappa20*TKappa20 + Kappa21*TKappa21 +
      Kappa22*TKappa22)*Sqr(gN)*Sqr(QSp) - 10*(Lambda1200*TLambda1200 +
      Lambda1201*TLambda1201 + Lambda1210*TLambda1210 + Lambda1211*TLambda1211)*
      Sqr(gN)*Sqr(QSp) - 10*Lambdax*TLambdax*Sqr(gN)*Sqr(QSp) - 5*MassBp*Sqr(gN)*(
      30*Power(QSp,4)*Sqr(gN) + 4*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(
      Lambda1210) + Sqr(Lambda1211))*Sqr(QH1p) + 4*(Sqr(Lambda1200) + Sqr(
      Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211))*Sqr(QH2p) + 6*(Sqr(Kappa00)
      + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12)
      + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22))*(Sqr(QDxbarp) + Sqr(QDxp) -
      Sqr(QSp)) + 4*Sqr(Lambdax)*(Sqr(QH1p) + Sqr(QH2p) - Sqr(QSp)) - 4*(Sqr(
      Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211))*Sqr(QSp)
      + 54*Sqr(gN)*Sqr(Qdp)*Sqr(QSp) + 54*Sqr(gN)*Sqr(QDxbarp)*Sqr(QSp) + 54*Sqr(
      gN)*Sqr(QDxp)*Sqr(QSp) + 18*Sqr(gN)*Sqr(Qep)*Sqr(QSp) + 36*Sqr(gN)*Sqr(QH1p)
      *Sqr(QSp) + 36*Sqr(gN)*Sqr(QH2p)*Sqr(QSp) + 12*Sqr(gN)*Sqr(QHpbarp)*Sqr(QSp)
      + 12*Sqr(gN)*Sqr(QHpp)*Sqr(QSp) + 36*Sqr(gN)*Sqr(QLp)*Sqr(QSp) + 108*Sqr(gN
      )*Sqr(QQp)*Sqr(QSp) + 54*Sqr(gN)*Sqr(QSp)*Sqr(Qup)) - 5*Sqr(gN)*(30*MassBp*
      Power(QSp,4)*Sqr(gN) - 3*(Kappa00*TKappa00 + Kappa01*TKappa01 + Kappa02*
      TKappa02 + Kappa10*TKappa10 + Kappa11*TKappa11 + Kappa12*TKappa12 + Kappa20*
      TKappa20 + Kappa21*TKappa21 + Kappa22*TKappa22)*Sqr(QDxbarp) - 3*(Kappa00*
      TKappa00 + Kappa01*TKappa01 + Kappa02*TKappa02 + Kappa10*TKappa10 + Kappa11*
      TKappa11 + Kappa12*TKappa12 + Kappa20*TKappa20 + Kappa21*TKappa21 + Kappa22*
      TKappa22)*Sqr(QDxp) - 2*(Lambda1200*TLambda1200 + Lambda1201*TLambda1201 +
      Lambda1210*TLambda1210 + Lambda1211*TLambda1211)*Sqr(QH1p) + 4*MassBp*(Sqr(
      Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211))*Sqr(QH1p)
      - 2*(Lambda1200*TLambda1200 + Lambda1201*TLambda1201 + Lambda1210*
      TLambda1210 + Lambda1211*TLambda1211)*Sqr(QH2p) + 4*MassBp*(Sqr(Lambda1200)
      + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211))*Sqr(QH2p) + 6*MassBp*
      (Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) +
      Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22))*(Sqr(QDxbarp) +
      Sqr(QDxp) - Sqr(QSp)) + 2*Lambdax*(2*Lambdax*MassBp - TLambdax)*(Sqr(QH1p) +
      Sqr(QH2p) - Sqr(QSp)) + 3*(Kappa00*TKappa00 + Kappa01*TKappa01 + Kappa02*
      TKappa02 + Kappa10*TKappa10 + Kappa11*TKappa11 + Kappa12*TKappa12 + Kappa20*
      TKappa20 + Kappa21*TKappa21 + Kappa22*TKappa22)*Sqr(QSp) + 2*(Lambda1200*
      TLambda1200 + Lambda1201*TLambda1201 + Lambda1210*TLambda1210 + Lambda1211*
      TLambda1211)*Sqr(QSp) - 4*MassBp*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(
      Lambda1210) + Sqr(Lambda1211))*Sqr(QSp) + 54*MassBp*Sqr(gN)*Sqr(Qdp)*Sqr(QSp
      ) + 54*MassBp*Sqr(gN)*Sqr(QDxbarp)*Sqr(QSp) + 54*MassBp*Sqr(gN)*Sqr(QDxp)*
      Sqr(QSp) + 18*MassBp*Sqr(gN)*Sqr(Qep)*Sqr(QSp) + 36*MassBp*Sqr(gN)*Sqr(QH1p)
      *Sqr(QSp) + 36*MassBp*Sqr(gN)*Sqr(QH2p)*Sqr(QSp) + 12*MassBp*Sqr(gN)*Sqr(
      QHpbarp)*Sqr(QSp) + 12*MassBp*Sqr(gN)*Sqr(QHpp)*Sqr(QSp) + 36*MassBp*Sqr(gN)
      *Sqr(QLp)*Sqr(QSp) + 108*MassBp*Sqr(gN)*Sqr(QQp)*Sqr(QSp) + 54*MassBp*Sqr(gN
      )*Sqr(QSp)*Sqr(Qup)));


   return twoLoop*deriv;
}

// leading log
double lowE6SSM_tuning_calculator::deriv_dlead_log_one_loop_ms2_dLambdax() const
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
   const double AYu22 = get_AYu22();
   const double MassB = model.get_MassB();
   const double MassWB = model.get_MassWB();
   const double MassG = model.get_MassG();
   const double MassBp = model.get_MassBp();

   const double mq200 = model.get_mq2(0,0);
   const double mq201 = model.get_mq2(0,1);
   const double mq202 = model.get_mq2(0,2);
   const double mq210 = model.get_mq2(1,0);
   const double mq211 = model.get_mq2(1,1);
   const double mq212 = model.get_mq2(1,2);
   const double mq220 = model.get_mq2(2,0);
   const double mq221 = model.get_mq2(2,1);
   const double mq222 = model.get_mq2(2,2);
   const double ml200 = model.get_ml2(0,0);
   const double ml201 = model.get_ml2(0,1);
   const double ml202 = model.get_ml2(0,2);
   const double ml210 = model.get_ml2(1,0);
   const double ml211 = model.get_ml2(1,1);
   const double ml212 = model.get_ml2(1,2);
   const double ml220 = model.get_ml2(2,0);
   const double ml221 = model.get_ml2(2,1);
   const double ml222 = model.get_ml2(2,2);
   const double mHd2 = model.get_mHd2();
   const double mHu2 = model.get_mHu2();
   const double md200 = model.get_md2(0,0);
   const double md201 = model.get_md2(0,1);
   const double md202 = model.get_md2(0,2);
   const double md210 = model.get_md2(1,0);
   const double md211 = model.get_md2(1,1);
   const double md212 = model.get_md2(1,2);
   const double md220 = model.get_md2(2,0);
   const double md221 = model.get_md2(2,1);
   const double md222 = model.get_md2(2,2);
   const double mu200 = model.get_mu2(0,0);
   const double mu201 = model.get_mu2(0,1);
   const double mu202 = model.get_mu2(0,2);
   const double mu210 = model.get_mu2(1,0);
   const double mu211 = model.get_mu2(1,1);
   const double mu212 = model.get_mu2(1,2);
   const double mu220 = model.get_mu2(2,0);
   const double mu221 = model.get_mu2(2,1);
   const double mu222 = model.get_mu2(2,2);
   const double me200 = model.get_me2(0,0);
   const double me201 = model.get_me2(0,1);
   const double me202 = model.get_me2(0,2);
   const double me210 = model.get_me2(1,0);
   const double me211 = model.get_me2(1,1);
   const double me212 = model.get_me2(1,2);
   const double me220 = model.get_me2(2,0);
   const double me221 = model.get_me2(2,1);
   const double me222 = model.get_me2(2,2);
   const double ms2 = model.get_ms2();
   const double mH1I200 = model.get_mH1I2(0,0);
   const double mH1I201 = model.get_mH1I2(0,1);
   const double mH1I210 = model.get_mH1I2(1,0);
   const double mH1I211 = model.get_mH1I2(1,1);
   const double mH2I200 = model.get_mH2I2(0,0);
   const double mH2I201 = model.get_mH2I2(0,1);
   const double mH2I210 = model.get_mH2I2(1,0);
   const double mH2I211 = model.get_mH2I2(1,1);
   const double msI200 = model.get_msI2(0,0);
   const double msI201 = model.get_msI2(0,1);
   const double msI210 = model.get_msI2(1,0);
   const double msI211 = model.get_msI2(1,1);
   const double mDx200 = model.get_mDx2(0,0);
   const double mDx201 = model.get_mDx2(0,1);
   const double mDx202 = model.get_mDx2(0,2);
   const double mDx210 = model.get_mDx2(1,0);
   const double mDx211 = model.get_mDx2(1,1);
   const double mDx212 = model.get_mDx2(1,2);
   const double mDx220 = model.get_mDx2(2,0);
   const double mDx221 = model.get_mDx2(2,1);
   const double mDx222 = model.get_mDx2(2,2);
   const double mDxbar200 = model.get_mDxbar2(0,0);
   const double mDxbar201 = model.get_mDxbar2(0,1);
   const double mDxbar202 = model.get_mDxbar2(0,2);
   const double mDxbar210 = model.get_mDxbar2(1,0);
   const double mDxbar211 = model.get_mDxbar2(1,1);
   const double mDxbar212 = model.get_mDxbar2(1,2);
   const double mDxbar220 = model.get_mDxbar2(2,0);
   const double mDxbar221 = model.get_mDxbar2(2,1);
   const double mDxbar222 = model.get_mDxbar2(2,2);
   const double mHp2 = model.get_mHp2();
   const double mHpbar2 = model.get_mHpbar2();

   double deriv = 4*Kappa00*Lambdax*(6*(2*Kappa00*mDx200 + Kappa10*mDx201 +
      Kappa20*mDx202 + Kappa10*mDx210 + Kappa20*mDx220) + 6*(2*Kappa00*mDxbar200 +
      Kappa01*mDxbar201 + Kappa02*mDxbar202 + Kappa01*mDxbar210 + Kappa02*
      mDxbar220) + 12*Kappa00*ms2) + 4*Kappa01*Lambdax*(6*(2*Kappa01*mDx200 +
      Kappa11*mDx201 + Kappa21*mDx202 + Kappa11*mDx210 + Kappa21*mDx220) + 6*(
      Kappa00*mDxbar201 + Kappa00*mDxbar210 + 2*Kappa01*mDxbar211 + Kappa02*
      mDxbar212 + Kappa02*mDxbar221) + 12*Kappa01*ms2) + 4*Kappa02*Lambdax*(6*(2*
      Kappa02*mDx200 + Kappa12*mDx201 + Kappa22*mDx202 + Kappa12*mDx210 + Kappa22*
      mDx220) + 6*(Kappa00*mDxbar202 + Kappa01*mDxbar212 + Kappa00*mDxbar220 +
      Kappa01*mDxbar221 + 2*Kappa02*mDxbar222) + 12*Kappa02*ms2) + 4*Kappa10*
      Lambdax*(6*(Kappa00*mDx201 + Kappa00*mDx210 + 2*Kappa10*mDx211 + Kappa20*
      mDx212 + Kappa20*mDx221) + 6*(2*Kappa10*mDxbar200 + Kappa11*mDxbar201 +
      Kappa12*mDxbar202 + Kappa11*mDxbar210 + Kappa12*mDxbar220) + 12*Kappa10*ms2)
      + 4*Kappa11*Lambdax*(6*(Kappa01*mDx201 + Kappa01*mDx210 + 2*Kappa11*mDx211
      + Kappa21*mDx212 + Kappa21*mDx221) + 6*(Kappa10*mDxbar201 + Kappa10*
      mDxbar210 + 2*Kappa11*mDxbar211 + Kappa12*mDxbar212 + Kappa12*mDxbar221) +
      12*Kappa11*ms2) + 4*Kappa12*Lambdax*(6*(Kappa02*mDx201 + Kappa02*mDx210 + 2*
      Kappa12*mDx211 + Kappa22*mDx212 + Kappa22*mDx221) + 6*(Kappa10*mDxbar202 +
      Kappa11*mDxbar212 + Kappa10*mDxbar220 + Kappa11*mDxbar221 + 2*Kappa12*
      mDxbar222) + 12*Kappa12*ms2) + 4*Kappa20*Lambdax*(6*(Kappa00*mDx202 +
      Kappa10*mDx212 + Kappa00*mDx220 + Kappa10*mDx221 + 2*Kappa20*mDx222) + 6*(2*
      Kappa20*mDxbar200 + Kappa21*mDxbar201 + Kappa22*mDxbar202 + Kappa21*
      mDxbar210 + Kappa22*mDxbar220) + 12*Kappa20*ms2) + 4*Kappa21*Lambdax*(6*(
      Kappa01*mDx202 + Kappa11*mDx212 + Kappa01*mDx220 + Kappa11*mDx221 + 2*
      Kappa21*mDx222) + 6*(Kappa20*mDxbar201 + Kappa20*mDxbar210 + 2*Kappa21*
      mDxbar211 + Kappa22*mDxbar212 + Kappa22*mDxbar221) + 12*Kappa21*ms2) + 4*
      Kappa22*Lambdax*(6*(Kappa02*mDx202 + Kappa12*mDx212 + Kappa02*mDx220 +
      Kappa12*mDx221 + 2*Kappa22*mDx222) + 6*(Kappa20*mDxbar202 + Kappa21*
      mDxbar212 + Kappa20*mDxbar220 + Kappa21*mDxbar221 + 2*Kappa22*mDxbar222) +
      12*Kappa22*ms2) + 4*Lambda1200*Lambdax*(4*(2*Lambda1200*mH1I200 + Lambda1201
      *mH1I201 + Lambda1201*mH1I210) + 4*(2*Lambda1200*mH2I200 + Lambda1210*
      mH2I201 + Lambda1210*mH2I210) + 8*Lambda1200*ms2) + 4*Lambda1201*Lambdax*(4*
      (Lambda1200*mH1I201 + Lambda1200*mH1I210 + 2*Lambda1201*mH1I211) + 4*(2*
      Lambda1201*mH2I200 + Lambda1211*mH2I201 + Lambda1211*mH2I210) + 8*Lambda1201
      *ms2) + 4*Lambda1210*Lambdax*(4*(2*Lambda1210*mH1I200 + Lambda1211*mH1I201 +
      Lambda1211*mH1I210) + 4*(Lambda1200*mH2I201 + Lambda1200*mH2I210 + 2*
      Lambda1210*mH2I211) + 8*Lambda1210*ms2) + 4*Lambda1211*Lambdax*(4*(
      Lambda1210*mH1I201 + Lambda1210*mH1I210 + 2*Lambda1211*mH1I211) + 4*(
      Lambda1201*mH2I201 + Lambda1201*mH2I210 + 2*Lambda1211*mH2I211) + 8*
      Lambda1211*ms2) + 12*TKappa00*(8*ALambdax*Kappa00*Lambdax + 4*Lambdax*
      TKappa00) + 12*TKappa01*(8*ALambdax*Kappa01*Lambdax + 4*Lambdax*TKappa01) +
      12*TKappa02*(8*ALambdax*Kappa02*Lambdax + 4*Lambdax*TKappa02) + 12*TKappa10*
      (8*ALambdax*Kappa10*Lambdax + 4*Lambdax*TKappa10) + 12*TKappa11*(8*ALambdax*
      Kappa11*Lambdax + 4*Lambdax*TKappa11) + 12*TKappa12*(8*ALambdax*Kappa12*
      Lambdax + 4*Lambdax*TKappa12) + 12*TKappa20*(8*ALambdax*Kappa20*Lambdax + 4*
      Lambdax*TKappa20) + 12*TKappa21*(8*ALambdax*Kappa21*Lambdax + 4*Lambdax*
      TKappa21) + 12*TKappa22*(8*ALambdax*Kappa22*Lambdax + 4*Lambdax*TKappa22) +
      8*TLambda1200*(8*ALambdax*Lambda1200*Lambdax + 4*Lambdax*TLambda1200) + 8*
      TLambda1201*(8*ALambdax*Lambda1201*Lambdax + 4*Lambdax*TLambda1201) + 8*
      TLambda1210*(8*ALambdax*Lambda1210*Lambdax + 4*Lambdax*TLambda1210) + 8*
      TLambda1211*(8*ALambdax*Lambda1211*Lambdax + 4*Lambdax*TLambda1211) + (4*
      Lambdax*mHd2 + 4*Lambdax*mHu2 + 4*Lambdax*ms2 + 4*Lambdax*Sqr(ALambdax))*(4*
      QH1p*QSp*Sqr(gN) + 4*Sqr(Lambdax)) + (4*Lambdax*mHd2 + 4*Lambdax*mHu2 + 4*
      Lambdax*ms2 + 4*Lambdax*Sqr(ALambdax))*(4*QH2p*QSp*Sqr(gN) + 4*Sqr(Lambdax))
      + (8*Lambdax*(mHd2 + mHu2 + ms2) + 8*Lambdax*Sqr(ALambdax))*(6*(Sqr(Kappa00
      ) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12)
      + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 4*(Sqr(Lambda1200) + Sqr(
      Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 4*Sqr(Lambdax) + 2*Sqr(gN
      )*Sqr(QSp)) + 8*Lambdax*(6*(Kappa00*(Kappa00*mDxbar200 + Kappa01*mDxbar210 +
      Kappa02*mDxbar220) + Kappa10*(Kappa10*mDxbar200 + Kappa11*mDxbar210 +
      Kappa12*mDxbar220) + Kappa20*(Kappa20*mDxbar200 + Kappa21*mDxbar210 +
      Kappa22*mDxbar220) + Kappa01*(Kappa00*mDxbar201 + Kappa01*mDxbar211 +
      Kappa02*mDxbar221) + Kappa11*(Kappa10*mDxbar201 + Kappa11*mDxbar211 +
      Kappa12*mDxbar221) + Kappa21*(Kappa20*mDxbar201 + Kappa21*mDxbar211 +
      Kappa22*mDxbar221) + Kappa02*(Kappa00*mDxbar202 + Kappa01*mDxbar212 +
      Kappa02*mDxbar222) + Kappa12*(Kappa10*mDxbar202 + Kappa11*mDxbar212 +
      Kappa12*mDxbar222) + Kappa22*(Kappa20*mDxbar202 + Kappa21*mDxbar212 +
      Kappa22*mDxbar222)) + 4*(Lambda1200*(Lambda1200*mH1I200 + Lambda1201*mH1I201
      ) + Lambda1210*(Lambda1210*mH1I200 + Lambda1211*mH1I201) + Lambda1201*(
      Lambda1200*mH1I210 + Lambda1201*mH1I211) + Lambda1211*(Lambda1210*mH1I210 +
      Lambda1211*mH1I211)) + 2*QSp*(3*(md200 + md211 + md222)*Qdp + 3*(mDxbar200 +
      mDxbar211 + mDxbar222)*QDxbarp + 3*(mDx200 + mDx211 + mDx222)*QDxp + (me200
      + me211 + me222)*Qep + 2*(mH1I200 + mH1I211)*QH1p + 2*mHd2*QH1p + 2*(
      mH2I200 + mH2I211)*QH2p + 2*mHu2*QH2p + 2*mHpbar2*QHpbarp + 2*mHp2*QHpp + 2*
      (ml200 + ml211 + ml222)*QLp + 6*(mq200 + mq211 + mq222)*QQp + ms2*QSp + (
      msI200 + msI211)*QSp + 3*(mu200 + mu211 + mu222)*Qup)*Sqr(gN) + 6*ms2*(Sqr(
      Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(
      Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 6*((Kappa00*Kappa10
      + Kappa01*Kappa11 + Kappa02*Kappa12)*mDx201 + (Kappa00*Kappa20 + Kappa01*
      Kappa21 + Kappa02*Kappa22)*mDx202 + (Kappa00*Kappa10 + Kappa01*Kappa11 +
      Kappa02*Kappa12)*mDx210 + (Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*
      Kappa22)*mDx212 + (Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)*
      mDx220 + (Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)*mDx221 +
      mDx200*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02)) + mDx211*(Sqr(Kappa10) +
      Sqr(Kappa11) + Sqr(Kappa12)) + mDx222*(Sqr(Kappa20) + Sqr(Kappa21) + Sqr(
      Kappa22))) + 4*ms2*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) +
      Sqr(Lambda1211)) + 4*((Lambda1200*Lambda1210 + Lambda1201*Lambda1211)*
      mH2I201 + (Lambda1200*Lambda1210 + Lambda1201*Lambda1211)*mH2I210 + mH2I200*
      (Sqr(Lambda1200) + Sqr(Lambda1201)) + mH2I211*(Sqr(Lambda1210) + Sqr(
      Lambda1211))) + 4*(mHd2 + mHu2 + ms2)*Sqr(Lambdax) + 4*Sqr(ALambdax)*Sqr(
      Lambdax) - 8*Sqr(gN)*Sqr(MassBp)*Sqr(QSp) + 6*(Sqr(TKappa00) + Sqr(TKappa01)
      + Sqr(TKappa02) + Sqr(TKappa10) + Sqr(TKappa11) + Sqr(TKappa12) + Sqr(
      TKappa20) + Sqr(TKappa21) + Sqr(TKappa22)) + 4*(Sqr(TLambda1200) + Sqr(
      TLambda1201) + Sqr(TLambda1210) + Sqr(TLambda1211))) + 8*Lambdax*(6*(Yd00*(
      md200*Yd00 + md201*Yd10 + md202*Yd20) + Yd10*(md210*Yd00 + md211*Yd10 +
      md212*Yd20) + Yd20*(md220*Yd00 + md221*Yd10 + md222*Yd20) + Yd01*(md200*Yd01
      + md201*Yd11 + md202*Yd21) + Yd11*(md210*Yd01 + md211*Yd11 + md212*Yd21) +
      Yd21*(md220*Yd01 + md221*Yd11 + md222*Yd21) + Yd02*(md200*Yd02 + md201*Yd12
      + md202*Yd22) + Yd12*(md210*Yd02 + md211*Yd12 + md212*Yd22) + Yd22*(md220*
      Yd02 + md221*Yd12 + md222*Yd22)) + 6*(Yd00*(mq200*Yd00 + mq201*Yd01 + mq202*
      Yd02) + Yd01*(mq210*Yd00 + mq211*Yd01 + mq212*Yd02) + Yd02*(mq220*Yd00 +
      mq221*Yd01 + mq222*Yd02) + Yd10*(mq200*Yd10 + mq201*Yd11 + mq202*Yd12) +
      Yd11*(mq210*Yd10 + mq211*Yd11 + mq212*Yd12) + Yd12*(mq220*Yd10 + mq221*Yd11
      + mq222*Yd12) + Yd20*(mq200*Yd20 + mq201*Yd21 + mq202*Yd22) + Yd21*(mq210*
      Yd20 + mq211*Yd21 + mq212*Yd22) + Yd22*(mq220*Yd20 + mq221*Yd21 + mq222*Yd22
      )) + 2*(Ye00*(me200*Ye00 + me201*Ye10 + me202*Ye20) + Ye10*(me210*Ye00 +
      me211*Ye10 + me212*Ye20) + Ye20*(me220*Ye00 + me221*Ye10 + me222*Ye20) +
      Ye01*(me200*Ye01 + me201*Ye11 + me202*Ye21) + Ye11*(me210*Ye01 + me211*Ye11
      + me212*Ye21) + Ye21*(me220*Ye01 + me221*Ye11 + me222*Ye21) + Ye02*(me200*
      Ye02 + me201*Ye12 + me202*Ye22) + Ye12*(me210*Ye02 + me211*Ye12 + me212*Ye22
      ) + Ye22*(me220*Ye02 + me221*Ye12 + me222*Ye22)) + 2*(Ye00*(ml200*Ye00 +
      ml201*Ye01 + ml202*Ye02) + Ye01*(ml210*Ye00 + ml211*Ye01 + ml212*Ye02) +
      Ye02*(ml220*Ye00 + ml221*Ye01 + ml222*Ye02) + Ye10*(ml200*Ye10 + ml201*Ye11
      + ml202*Ye12) + Ye11*(ml210*Ye10 + ml211*Ye11 + ml212*Ye12) + Ye12*(ml220*
      Ye10 + ml221*Ye11 + ml222*Ye12) + Ye20*(ml200*Ye20 + ml201*Ye21 + ml202*Ye22
      ) + Ye21*(ml210*Ye20 + ml211*Ye21 + ml212*Ye22) + Ye22*(ml220*Ye20 + ml221*
      Ye21 + ml222*Ye22)) - 0.6*(md200 + md211 + md222 - mDx200 - mDx211 - mDx222
      + mDxbar200 + mDxbar211 + mDxbar222 + me200 + me211 + me222 - mH1I200 -
      mH1I211 + mH2I200 + mH2I211 - mHd2 - mHp2 + mHpbar2 + mHu2 - ml200 - ml211 -
      ml222 + mq200 + mq211 + mq222 - 2*(mu200 + mu211 + mu222))*Sqr(g1) + 2*QH1p
      *(3*(md200 + md211 + md222)*Qdp + 3*(mDxbar200 + mDxbar211 + mDxbar222)*
      QDxbarp + 3*(mDx200 + mDx211 + mDx222)*QDxp + (me200 + me211 + me222)*Qep +
      2*(mH1I200 + mH1I211)*QH1p + 2*mHd2*QH1p + 2*(mH2I200 + mH2I211)*QH2p + 2*
      mHu2*QH2p + 2*mHpbar2*QHpbarp + 2*mHp2*QHpp + 2*(ml200 + ml211 + ml222)*QLp
      + 6*(mq200 + mq211 + mq222)*QQp + ms2*QSp + (msI200 + msI211)*QSp + 3*(mu200
      + mu211 + mu222)*Qup)*Sqr(gN) + 2*mHd2*Sqr(Lambdax) + 2*mHu2*Sqr(Lambdax) +
      2*ms2*Sqr(Lambdax) + 2*Sqr(ALambdax)*Sqr(Lambdax) - 1.2*Sqr(g1)*Sqr(MassB)
      - 6*Sqr(g2)*Sqr(MassWB) - 8*Sqr(gN)*Sqr(MassBp)*Sqr(QH1p) + 6*(Sqr(TYd00) +
      Sqr(TYd01) + Sqr(TYd02) + Sqr(TYd10) + Sqr(TYd11) + Sqr(TYd12) + Sqr(TYd20)
      + Sqr(TYd21) + Sqr(TYd22)) + 2*(Sqr(TYe00) + Sqr(TYe01) + Sqr(TYe02) + Sqr(
      TYe10) + Sqr(TYe11) + Sqr(TYe12) + Sqr(TYe20) + Sqr(TYe21) + Sqr(TYe22)) + 6
      *mHd2*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12)
      + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + 2*mHd2*(Sqr(Ye00) + Sqr(Ye01) + Sqr(
      Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22
      ))) + 8*Lambdax*(mHd2 + mHu2 + ms2)*(-0.6*Sqr(g1) - 3*Sqr(g2) + 3*(Sqr(
      Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(
      Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*(Sqr(Lambda1200)
      + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 12*Sqr(Lambdax) - 2
      *Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QSp) + 3*(Sqr(Yd00)
      + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) +
      Sqr(Yd21) + Sqr(Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr
      (Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22) + 3*(Sqr(Yu00) + Sqr(
      Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21
      ) + Sqr(Yu22))) + 8*(mHd2 + mHu2 + ms2)*(4*Power(Lambdax,3) - 0.6*Lambdax*
      Sqr(g1) - 3*Lambdax*Sqr(g2) + 3*Lambdax*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(
      Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(
      Kappa21) + Sqr(Kappa22)) + 2*Lambdax*(Sqr(Lambda1200) + Sqr(Lambda1201) +
      Sqr(Lambda1210) + Sqr(Lambda1211)) - 2*Lambdax*Sqr(gN)*Sqr(QH1p) - 2*Lambdax
      *Sqr(gN)*Sqr(QH2p) - 2*Lambdax*Sqr(gN)*Sqr(QSp) + 3*Lambdax*(Sqr(Yd00) + Sqr
      (Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(
      Yd21) + Sqr(Yd22)) + Lambdax*(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10)
      + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22)) + 3*Lambdax*(
      Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(
      Yu20) + Sqr(Yu21) + Sqr(Yu22))) + 8*Lambdax*(6*(Yu00*(mq200*Yu00 + mq201*
      Yu01 + mq202*Yu02) + Yu01*(mq210*Yu00 + mq211*Yu01 + mq212*Yu02) + Yu02*(
      mq220*Yu00 + mq221*Yu01 + mq222*Yu02) + Yu10*(mq200*Yu10 + mq201*Yu11 +
      mq202*Yu12) + Yu11*(mq210*Yu10 + mq211*Yu11 + mq212*Yu12) + Yu12*(mq220*Yu10
      + mq221*Yu11 + mq222*Yu12) + Yu20*(mq200*Yu20 + mq201*Yu21 + mq202*Yu22) +
      Yu21*(mq210*Yu20 + mq211*Yu21 + mq212*Yu22) + Yu22*(mq220*Yu20 + mq221*Yu21
      + mq222*Yu22)) + 6*(Yu00*(mu200*Yu00 + mu201*Yu10 + mu202*Yu20) + Yu10*(
      mu210*Yu00 + mu211*Yu10 + mu212*Yu20) + Yu20*(mu220*Yu00 + mu221*Yu10 +
      mu222*Yu20) + Yu01*(mu200*Yu01 + mu201*Yu11 + mu202*Yu21) + Yu11*(mu210*Yu01
      + mu211*Yu11 + mu212*Yu21) + Yu21*(mu220*Yu01 + mu221*Yu11 + mu222*Yu21) +
      Yu02*(mu200*Yu02 + mu201*Yu12 + mu202*Yu22) + Yu12*(mu210*Yu02 + mu211*Yu12
      + mu212*Yu22) + Yu22*(mu220*Yu02 + mu221*Yu12 + mu222*Yu22)) + 0.6*(md200 +
      md211 + md222 - mDx200 - mDx211 - mDx222 + mDxbar200 + mDxbar211 + mDxbar222
      + me200 + me211 + me222 - mH1I200 - mH1I211 + mH2I200 + mH2I211 - mHd2 -
      mHp2 + mHpbar2 + mHu2 - ml200 - ml211 - ml222 + mq200 + mq211 + mq222 - 2*(
      mu200 + mu211 + mu222))*Sqr(g1) + 2*QH2p*(3*(md200 + md211 + md222)*Qdp + 3*
      (mDxbar200 + mDxbar211 + mDxbar222)*QDxbarp + 3*(mDx200 + mDx211 + mDx222)*
      QDxp + (me200 + me211 + me222)*Qep + 2*(mH1I200 + mH1I211)*QH1p + 2*mHd2*
      QH1p + 2*(mH2I200 + mH2I211)*QH2p + 2*mHu2*QH2p + 2*mHpbar2*QHpbarp + 2*mHp2
      *QHpp + 2*(ml200 + ml211 + ml222)*QLp + 6*(mq200 + mq211 + mq222)*QQp + ms2*
      QSp + (msI200 + msI211)*QSp + 3*(mu200 + mu211 + mu222)*Qup)*Sqr(gN) + 2*
      mHd2*Sqr(Lambdax) + 2*mHu2*Sqr(Lambdax) + 2*ms2*Sqr(Lambdax) + 2*Sqr(
      ALambdax)*Sqr(Lambdax) - 1.2*Sqr(g1)*Sqr(MassB) - 6*Sqr(g2)*Sqr(MassWB) - 8*
      Sqr(gN)*Sqr(MassBp)*Sqr(QH2p) + 6*(Sqr(TYu00) + Sqr(TYu01) + Sqr(TYu02) +
      Sqr(TYu10) + Sqr(TYu11) + Sqr(TYu12) + Sqr(TYu20) + Sqr(TYu21) + Sqr(TYu22))
      + 6*mHu2*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(
      Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + 8*ALambdax*Lambdax*(6*(Kappa00
      *TKappa00 + Kappa01*TKappa01 + Kappa02*TKappa02 + Kappa10*TKappa10 + Kappa11
      *TKappa11 + Kappa12*TKappa12 + Kappa20*TKappa20 + Kappa21*TKappa21 + Kappa22
      *TKappa22) + 4*(Lambda1200*TLambda1200 + Lambda1201*TLambda1201 + Lambda1210
      *TLambda1210 + Lambda1211*TLambda1211) + 6*(TYd00*Yd00 + TYd01*Yd01 + TYd02*
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
      ) + Sqr(Yu21) + Sqr(Yu22)))) + 8*ALambdax*(6*Lambdax*(Kappa00*TKappa00 +
      Kappa01*TKappa01 + Kappa02*TKappa02 + Kappa10*TKappa10 + Kappa11*TKappa11 +
      Kappa12*TKappa12 + Kappa20*TKappa20 + Kappa21*TKappa21 + Kappa22*TKappa22) +
      4*Lambdax*(Lambda1200*TLambda1200 + Lambda1201*TLambda1201 + Lambda1210*
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

double lowE6SSM_tuning_calculator::deriv_dlead_log_one_loop_ms2_dALambdax() const
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
   const double AYu22 = get_AYu22();
   const double MassB = model.get_MassB();
   const double MassWB = model.get_MassWB();
   const double MassG = model.get_MassG();
   const double MassBp = model.get_MassBp();

   double deriv = 48*Kappa00*TKappa00*Sqr(Lambdax) + 48*Kappa01*TKappa01*Sqr(
      Lambdax) + 48*Kappa02*TKappa02*Sqr(Lambdax) + 48*Kappa10*TKappa10*Sqr(
      Lambdax) + 48*Kappa11*TKappa11*Sqr(Lambdax) + 48*Kappa12*TKappa12*Sqr(
      Lambdax) + 48*Kappa20*TKappa20*Sqr(Lambdax) + 48*Kappa21*TKappa21*Sqr(
      Lambdax) + 48*Kappa22*TKappa22*Sqr(Lambdax) + 32*Lambda1200*TLambda1200*Sqr(
      Lambdax) + 32*Lambda1201*TLambda1201*Sqr(Lambdax) + 32*Lambda1210*
      TLambda1210*Sqr(Lambdax) + 32*Lambda1211*TLambda1211*Sqr(Lambdax) + 4*
      ALambdax*Sqr(Lambdax)*(4*QH1p*QSp*Sqr(gN) + 4*Sqr(Lambdax)) + 4*ALambdax*Sqr
      (Lambdax)*(4*QH2p*QSp*Sqr(gN) + 4*Sqr(Lambdax)) + 8*ALambdax*Sqr(Lambdax)*(6
      *(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) +
      Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 4*(Sqr(
      Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 4*Sqr(
      Lambdax) + 2*Sqr(gN)*Sqr(QSp)) + 8*ALambdax*Sqr(Lambdax)*(-0.6*Sqr(g1) - 3*
      Sqr(g2) + 3*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr
      (Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*(
      Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 12*
      Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QSp
      ) + 3*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12)
      + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) +
      Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22) + 3*(
      Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(
      Yu20) + Sqr(Yu21) + Sqr(Yu22))) + 8*Lambdax*(6*Lambdax*(Kappa00*TKappa00 +
      Kappa01*TKappa01 + Kappa02*TKappa02 + Kappa10*TKappa10 + Kappa11*TKappa11 +
      Kappa12*TKappa12 + Kappa20*TKappa20 + Kappa21*TKappa21 + Kappa22*TKappa22) +
      4*Lambdax*(Lambda1200*TLambda1200 + Lambda1201*TLambda1201 + Lambda1210*
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

double lowE6SSM_tuning_calculator::deriv_dlead_log_one_loop_ms2_dAYu22() const
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
   const double ALamdbax = get_ALambdax();
   const double TYu00 = model.get_TYu(0,0);
   const double TYu01 = model.get_TYu(0,1);
   const double TYu02 = model.get_TYu(0,2);
   const double TYu10 = model.get_TYu(1,0);
   const double TYu11 = model.get_TYu(1,1);
   const double TYu12 = model.get_TYu(1,2);
   const double TYu20 = model.get_TYu(2,0);
   const double TYu21 = model.get_TYu(2,1);
   const double TYu22 = model.get_TYu(2,2);
   const double AYu22 = get_AYu22();
   const double MassB = model.get_MassB();
   const double MassWB = model.get_MassWB();
   const double MassG = model.get_MassG();
   const double MassBp = model.get_MassBp();

   double deriv = 48*Lambdax*TLambdax*Sqr(Yu22) + 48*AYu22*QQp*QSp*Sqr(gN)*Sqr(
      Yu22) + 48*AYu22*QSp*Qup*Sqr(gN)*Sqr(Yu22) + 12*AYu22*(4*QH2p*QSp*Sqr(gN) +
      4*Sqr(Lambdax))*Sqr(Yu22);


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dlead_log_one_loop_ms2_dmq222() const
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

   double deriv = 48*Power(gN,4)*QQp*Power(QSp,3) + 6*Qep*QSp*Sqr(gN)*(1.2*Sqr(
      g1) + 12*Qep*QQp*Sqr(gN)) + 4*QHpbarp*QSp*Sqr(gN)*(0.6*Sqr(g1) + 12*QHpbarp*
      QQp*Sqr(gN)) + 4*QHpp*QSp*Sqr(gN)*(-0.6*Sqr(g1) + 12*QHpp*QQp*Sqr(gN)) + 12*
      QLp*QSp*Sqr(gN)*(-0.6*Sqr(g1) + 12*QLp*QQp*Sqr(gN)) + (-0.4*Sqr(g1) + 12*
      QDxp*QQp*Sqr(gN))*(6*QDxp*QSp*Sqr(gN) + 6*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr
      (Kappa02))) + (-0.4*Sqr(g1) + 12*QDxp*QQp*Sqr(gN))*(6*QDxp*QSp*Sqr(gN) + 6*(
      Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12))) + (0.4*Sqr(g1) + 12*QDxbarp*QQp
      *Sqr(gN))*(6*QDxbarp*QSp*Sqr(gN) + 6*(Sqr(Kappa00) + Sqr(Kappa10) + Sqr(
      Kappa20))) + (0.4*Sqr(g1) + 12*QDxbarp*QQp*Sqr(gN))*(6*QDxbarp*QSp*Sqr(gN) +
      6*(Sqr(Kappa01) + Sqr(Kappa11) + Sqr(Kappa21))) + (0.4*Sqr(g1) + 12*QDxbarp
      *QQp*Sqr(gN))*(6*QDxbarp*QSp*Sqr(gN) + 6*(Sqr(Kappa02) + Sqr(Kappa12) + Sqr(
      Kappa22))) + (-0.4*Sqr(g1) + 12*QDxp*QQp*Sqr(gN))*(6*QDxp*QSp*Sqr(gN) + 6*(
      Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22))) + (0.6*Sqr(g1) + 12*QH2p*QQp*
      Sqr(gN))*(4*QH2p*QSp*Sqr(gN) + 4*(Sqr(Lambda1200) + Sqr(Lambda1201))) + (
      -0.6*Sqr(g1) + 12*QH1p*QQp*Sqr(gN))*(4*QH1p*QSp*Sqr(gN) + 4*(Sqr(Lambda1200)
      + Sqr(Lambda1210))) + (-0.6*Sqr(g1) + 12*QH1p*QQp*Sqr(gN))*(4*QH1p*QSp*Sqr(
      gN) + 4*(Sqr(Lambda1201) + Sqr(Lambda1211))) + (0.6*Sqr(g1) + 12*QH2p*QQp*
      Sqr(gN))*(4*QH2p*QSp*Sqr(gN) + 4*(Sqr(Lambda1210) + Sqr(Lambda1211))) + 24*
      QQp*QSp*Sqr(gN)*(0.2*Sqr(g1) + 12*Sqr(gN)*Sqr(QQp)) + 12*QQp*QSp*Sqr(gN)*(6*
      (Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) +
      Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 4*(Sqr(
      Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 4*Sqr(
      Lambdax) + 2*Sqr(gN)*Sqr(QSp)) + 24*Power(gN,4)*QQp*QSp*(9*Sqr(Qdp) + 9*Sqr(
      QDxbarp) + 9*Sqr(QDxp) + 3*Sqr(Qep) + 6*Sqr(QH1p) + 6*Sqr(QH2p) + 2*Sqr(
      QHpbarp) + 2*Sqr(QHpp) + 6*Sqr(QLp) + 18*Sqr(QQp) + 3*Sqr(QSp) + 9*Sqr(Qup))
      + 6*Qdp*QSp*Sqr(gN)*(0.4*Sqr(g1) + 12*Qdp*QQp*Sqr(gN) + 4*Sqr(Yd02)) + 6*
      Qdp*QSp*Sqr(gN)*(0.4*Sqr(g1) + 12*Qdp*QQp*Sqr(gN) + 4*Sqr(Yd12)) + 6*Qdp*QSp
      *Sqr(gN)*(0.4*Sqr(g1) + 12*Qdp*QQp*Sqr(gN) + 4*Sqr(Yd22)) + (4*QH1p*QSp*Sqr(
      gN) + 4*Sqr(Lambdax))*(-0.6*Sqr(g1) + 12*QH1p*QQp*Sqr(gN) + 6*(Sqr(Yd02) +
      Sqr(Yd12) + Sqr(Yd22))) + 6*QSp*Qup*Sqr(gN)*(-0.8*Sqr(g1) + 12*QQp*Qup*Sqr(
      gN) + 4*Sqr(Yu02)) + 6*QSp*Qup*Sqr(gN)*(-0.8*Sqr(g1) + 12*QQp*Qup*Sqr(gN) +
      4*Sqr(Yu12)) + 12*QQp*QSp*Sqr(gN)*(0.2*Sqr(g1) + 12*Sqr(gN)*Sqr(QQp) + 2*Sqr
      (Yd02) + 2*Sqr(Yd12) + 2*Sqr(Yd22) + 2*Sqr(Yu02) + 2*Sqr(Yu12) + 2*Sqr(Yu22)
      ) + 6*QSp*Qup*Sqr(gN)*(-0.8*Sqr(g1) + 12*QQp*Qup*Sqr(gN) + 4*Sqr(Yu22)) + (4
      *QH2p*QSp*Sqr(gN) + 4*Sqr(Lambdax))*(0.6*Sqr(g1) + 12*QH2p*QQp*Sqr(gN) + 6*(
      Sqr(Yu02) + Sqr(Yu12) + Sqr(Yu22)));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dlead_log_one_loop_ms2_dmHd2() const
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

   double deriv = 16*Power(gN,4)*QH1p*Power(QSp,3) + 4*QHpbarp*QSp*Sqr(gN)*(
      -0.6*Sqr(g1) + 4*QH1p*QHpbarp*Sqr(gN)) + 4*QHpp*QSp*Sqr(gN)*(0.6*Sqr(g1) + 4
      *QH1p*QHpp*Sqr(gN)) + 18*QSp*Qup*Sqr(gN)*(0.8*Sqr(g1) + 4*QH1p*Qup*Sqr(gN))
      + (0.4*Sqr(g1) + 4*QDxp*QH1p*Sqr(gN))*(6*QDxp*QSp*Sqr(gN) + 6*(Sqr(Kappa00)
      + Sqr(Kappa01) + Sqr(Kappa02))) + (0.4*Sqr(g1) + 4*QDxp*QH1p*Sqr(gN))*(6*
      QDxp*QSp*Sqr(gN) + 6*(Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12))) + (-0.4*
      Sqr(g1) + 4*QDxbarp*QH1p*Sqr(gN))*(6*QDxbarp*QSp*Sqr(gN) + 6*(Sqr(Kappa00) +
      Sqr(Kappa10) + Sqr(Kappa20))) + (-0.4*Sqr(g1) + 4*QDxbarp*QH1p*Sqr(gN))*(6*
      QDxbarp*QSp*Sqr(gN) + 6*(Sqr(Kappa01) + Sqr(Kappa11) + Sqr(Kappa21))) + (
      -0.4*Sqr(g1) + 4*QDxbarp*QH1p*Sqr(gN))*(6*QDxbarp*QSp*Sqr(gN) + 6*(Sqr(
      Kappa02) + Sqr(Kappa12) + Sqr(Kappa22))) + (0.4*Sqr(g1) + 4*QDxp*QH1p*Sqr(gN
      ))*(6*QDxp*QSp*Sqr(gN) + 6*(Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22))) + (
      -0.6*Sqr(g1) + 4*QH1p*QH2p*Sqr(gN))*(4*QH2p*QSp*Sqr(gN) + 4*(Sqr(Lambda1200)
      + Sqr(Lambda1201))) + (-0.6*Sqr(g1) + 4*QH1p*QH2p*Sqr(gN))*(4*QH2p*QSp*Sqr(
      gN) + 4*(Sqr(Lambda1210) + Sqr(Lambda1211))) + (-0.6*Sqr(g1) + 4*QH1p*QH2p*
      Sqr(gN) + 2*Sqr(Lambdax))*(4*QH2p*QSp*Sqr(gN) + 4*Sqr(Lambdax)) + (4*QH1p*
      QSp*Sqr(gN) + 4*(Sqr(Lambda1200) + Sqr(Lambda1210)))*(0.6*Sqr(g1) + 4*Sqr(gN
      )*Sqr(QH1p)) + (4*QH1p*QSp*Sqr(gN) + 4*(Sqr(Lambda1201) + Sqr(Lambda1211)))*
      (0.6*Sqr(g1) + 4*Sqr(gN)*Sqr(QH1p)) + (4*QH1p*QSp*Sqr(gN) + 4*Sqr(Lambdax))*
      (6*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11)
      + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 4*(Sqr(
      Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 4*Sqr(
      Lambdax) + 2*Sqr(gN)*Sqr(QSp)) + 8*Power(gN,4)*QH1p*QSp*(9*Sqr(Qdp) + 9*Sqr(
      QDxbarp) + 9*Sqr(QDxp) + 3*Sqr(Qep) + 6*Sqr(QH1p) + 6*Sqr(QH2p) + 2*Sqr(
      QHpbarp) + 2*Sqr(QHpp) + 6*Sqr(QLp) + 18*Sqr(QQp) + 3*Sqr(QSp) + 9*Sqr(Qup))
      + 6*Qdp*QSp*Sqr(gN)*(-0.4*Sqr(g1) + 4*Qdp*QH1p*Sqr(gN) + 4*(Sqr(Yd00) + Sqr
      (Yd01) + Sqr(Yd02))) + 6*Qdp*QSp*Sqr(gN)*(-0.4*Sqr(g1) + 4*Qdp*QH1p*Sqr(gN)
      + 4*(Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12))) + 12*QQp*QSp*Sqr(gN)*(-0.2*Sqr(g1)
      + 4*QH1p*QQp*Sqr(gN) + 2*(Sqr(Yd00) + Sqr(Yd10) + Sqr(Yd20))) + 12*QQp*QSp*
      Sqr(gN)*(-0.2*Sqr(g1) + 4*QH1p*QQp*Sqr(gN) + 2*(Sqr(Yd01) + Sqr(Yd11) + Sqr(
      Yd21))) + 12*QQp*QSp*Sqr(gN)*(-0.2*Sqr(g1) + 4*QH1p*QQp*Sqr(gN) + 2*(Sqr(
      Yd02) + Sqr(Yd12) + Sqr(Yd22))) + 6*Qdp*QSp*Sqr(gN)*(-0.4*Sqr(g1) + 4*Qdp*
      QH1p*Sqr(gN) + 4*(Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22))) + 2*Qep*QSp*Sqr(gN)*(
      -1.2*Sqr(g1) + 4*Qep*QH1p*Sqr(gN) + 4*(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02))) +
      2*Qep*QSp*Sqr(gN)*(-1.2*Sqr(g1) + 4*Qep*QH1p*Sqr(gN) + 4*(Sqr(Ye10) + Sqr(
      Ye11) + Sqr(Ye12))) + 4*QLp*QSp*Sqr(gN)*(0.6*Sqr(g1) + 4*QH1p*QLp*Sqr(gN) +
      2*(Sqr(Ye00) + Sqr(Ye10) + Sqr(Ye20))) + 4*QLp*QSp*Sqr(gN)*(0.6*Sqr(g1) + 4*
      QH1p*QLp*Sqr(gN) + 2*(Sqr(Ye01) + Sqr(Ye11) + Sqr(Ye21))) + 4*QLp*QSp*Sqr(gN
      )*(0.6*Sqr(g1) + 4*QH1p*QLp*Sqr(gN) + 2*(Sqr(Ye02) + Sqr(Ye12) + Sqr(Ye22)))
      + 2*Qep*QSp*Sqr(gN)*(-1.2*Sqr(g1) + 4*Qep*QH1p*Sqr(gN) + 4*(Sqr(Ye20) + Sqr
      (Ye21) + Sqr(Ye22))) + (4*QH1p*QSp*Sqr(gN) + 4*Sqr(Lambdax))*(0.6*Sqr(g1) +
      2*Sqr(Lambdax) + 4*Sqr(gN)*Sqr(QH1p) + 6*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02)
      + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + 2
      *(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) +
      Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22))) + 8*Lambdax*(4*Power(Lambdax,3) - 0.6*
      Lambdax*Sqr(g1) - 3*Lambdax*Sqr(g2) + 3*Lambdax*(Sqr(Kappa00) + Sqr(Kappa01)
      + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20)
      + Sqr(Kappa21) + Sqr(Kappa22)) + 2*Lambdax*(Sqr(Lambda1200) + Sqr(Lambda1201
      ) + Sqr(Lambda1210) + Sqr(Lambda1211)) - 2*Lambdax*Sqr(gN)*Sqr(QH1p) - 2*
      Lambdax*Sqr(gN)*Sqr(QH2p) - 2*Lambdax*Sqr(gN)*Sqr(QSp) + 3*Lambdax*(Sqr(Yd00
      ) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) +
      Sqr(Yd21) + Sqr(Yd22)) + Lambdax*(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(
      Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22)) + 3*
      Lambdax*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(
      Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dlead_log_one_loop_ms2_dmHu2() const
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

   double deriv = 16*Power(gN,4)*QH2p*Power(QSp,3) + 18*Qdp*QSp*Sqr(gN)*(0.4*
      Sqr(g1) + 4*Qdp*QH2p*Sqr(gN)) + 6*Qep*QSp*Sqr(gN)*(1.2*Sqr(g1) + 4*Qep*QH2p*
      Sqr(gN)) + 4*QHpbarp*QSp*Sqr(gN)*(0.6*Sqr(g1) + 4*QH2p*QHpbarp*Sqr(gN)) + 4*
      QHpp*QSp*Sqr(gN)*(-0.6*Sqr(g1) + 4*QH2p*QHpp*Sqr(gN)) + 12*QLp*QSp*Sqr(gN)*(
      -0.6*Sqr(g1) + 4*QH2p*QLp*Sqr(gN)) + (-0.4*Sqr(g1) + 4*QDxp*QH2p*Sqr(gN))*(6
      *QDxp*QSp*Sqr(gN) + 6*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02))) + (-0.4*
      Sqr(g1) + 4*QDxp*QH2p*Sqr(gN))*(6*QDxp*QSp*Sqr(gN) + 6*(Sqr(Kappa10) + Sqr(
      Kappa11) + Sqr(Kappa12))) + (0.4*Sqr(g1) + 4*QDxbarp*QH2p*Sqr(gN))*(6*
      QDxbarp*QSp*Sqr(gN) + 6*(Sqr(Kappa00) + Sqr(Kappa10) + Sqr(Kappa20))) + (0.4
      *Sqr(g1) + 4*QDxbarp*QH2p*Sqr(gN))*(6*QDxbarp*QSp*Sqr(gN) + 6*(Sqr(Kappa01)
      + Sqr(Kappa11) + Sqr(Kappa21))) + (0.4*Sqr(g1) + 4*QDxbarp*QH2p*Sqr(gN))*(6*
      QDxbarp*QSp*Sqr(gN) + 6*(Sqr(Kappa02) + Sqr(Kappa12) + Sqr(Kappa22))) + (
      -0.4*Sqr(g1) + 4*QDxp*QH2p*Sqr(gN))*(6*QDxp*QSp*Sqr(gN) + 6*(Sqr(Kappa20) +
      Sqr(Kappa21) + Sqr(Kappa22))) + (-0.6*Sqr(g1) + 4*QH1p*QH2p*Sqr(gN))*(4*QH1p
      *QSp*Sqr(gN) + 4*(Sqr(Lambda1200) + Sqr(Lambda1210))) + (-0.6*Sqr(g1) + 4*
      QH1p*QH2p*Sqr(gN))*(4*QH1p*QSp*Sqr(gN) + 4*(Sqr(Lambda1201) + Sqr(Lambda1211
      ))) + (-0.6*Sqr(g1) + 4*QH1p*QH2p*Sqr(gN) + 2*Sqr(Lambdax))*(4*QH1p*QSp*Sqr(
      gN) + 4*Sqr(Lambdax)) + (4*QH2p*QSp*Sqr(gN) + 4*(Sqr(Lambda1200) + Sqr(
      Lambda1201)))*(0.6*Sqr(g1) + 4*Sqr(gN)*Sqr(QH2p)) + (4*QH2p*QSp*Sqr(gN) + 4*
      (Sqr(Lambda1210) + Sqr(Lambda1211)))*(0.6*Sqr(g1) + 4*Sqr(gN)*Sqr(QH2p)) + (
      4*QH2p*QSp*Sqr(gN) + 4*Sqr(Lambdax))*(6*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(
      Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(
      Kappa21) + Sqr(Kappa22)) + 4*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(
      Lambda1210) + Sqr(Lambda1211)) + 4*Sqr(Lambdax) + 2*Sqr(gN)*Sqr(QSp)) + 8*
      Power(gN,4)*QH2p*QSp*(9*Sqr(Qdp) + 9*Sqr(QDxbarp) + 9*Sqr(QDxp) + 3*Sqr(Qep)
      + 6*Sqr(QH1p) + 6*Sqr(QH2p) + 2*Sqr(QHpbarp) + 2*Sqr(QHpp) + 6*Sqr(QLp) +
      18*Sqr(QQp) + 3*Sqr(QSp) + 9*Sqr(Qup)) + 6*QSp*Qup*Sqr(gN)*(-0.8*Sqr(g1) + 4
      *QH2p*Qup*Sqr(gN) + 4*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02))) + 6*QSp*Qup*Sqr(
      gN)*(-0.8*Sqr(g1) + 4*QH2p*Qup*Sqr(gN) + 4*(Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12
      ))) + 12*QQp*QSp*Sqr(gN)*(0.2*Sqr(g1) + 4*QH2p*QQp*Sqr(gN) + 2*(Sqr(Yu00) +
      Sqr(Yu10) + Sqr(Yu20))) + 12*QQp*QSp*Sqr(gN)*(0.2*Sqr(g1) + 4*QH2p*QQp*Sqr(
      gN) + 2*(Sqr(Yu01) + Sqr(Yu11) + Sqr(Yu21))) + 12*QQp*QSp*Sqr(gN)*(0.2*Sqr(
      g1) + 4*QH2p*QQp*Sqr(gN) + 2*(Sqr(Yu02) + Sqr(Yu12) + Sqr(Yu22))) + 6*QSp*
      Qup*Sqr(gN)*(-0.8*Sqr(g1) + 4*QH2p*Qup*Sqr(gN) + 4*(Sqr(Yu20) + Sqr(Yu21) +
      Sqr(Yu22))) + (4*QH2p*QSp*Sqr(gN) + 4*Sqr(Lambdax))*(0.6*Sqr(g1) + 2*Sqr(
      Lambdax) + 4*Sqr(gN)*Sqr(QH2p) + 6*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(
      Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + 8*
      Lambdax*(4*Power(Lambdax,3) - 0.6*Lambdax*Sqr(g1) - 3*Lambdax*Sqr(g2) + 3*
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

double lowE6SSM_tuning_calculator::deriv_dlead_log_one_loop_ms2_dmu222() const
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

   double deriv = 24*Power(gN,4)*Power(QSp,3)*Qup + 18*Qdp*QSp*Sqr(gN)*(-0.8*
      Sqr(g1) + 6*Qdp*Qup*Sqr(gN)) + 6*Qep*QSp*Sqr(gN)*(-2.4*Sqr(g1) + 6*Qep*Qup*
      Sqr(gN)) + 4*QHpbarp*QSp*Sqr(gN)*(-1.2*Sqr(g1) + 6*QHpbarp*Qup*Sqr(gN)) + 4*
      QHpp*QSp*Sqr(gN)*(1.2*Sqr(g1) + 6*QHpp*Qup*Sqr(gN)) + 12*QLp*QSp*Sqr(gN)*(
      1.2*Sqr(g1) + 6*QLp*Qup*Sqr(gN)) + (0.8*Sqr(g1) + 6*QDxp*Qup*Sqr(gN))*(6*
      QDxp*QSp*Sqr(gN) + 6*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02))) + (0.8*
      Sqr(g1) + 6*QDxp*Qup*Sqr(gN))*(6*QDxp*QSp*Sqr(gN) + 6*(Sqr(Kappa10) + Sqr(
      Kappa11) + Sqr(Kappa12))) + (-0.8*Sqr(g1) + 6*QDxbarp*Qup*Sqr(gN))*(6*
      QDxbarp*QSp*Sqr(gN) + 6*(Sqr(Kappa00) + Sqr(Kappa10) + Sqr(Kappa20))) + (
      -0.8*Sqr(g1) + 6*QDxbarp*Qup*Sqr(gN))*(6*QDxbarp*QSp*Sqr(gN) + 6*(Sqr(
      Kappa01) + Sqr(Kappa11) + Sqr(Kappa21))) + (-0.8*Sqr(g1) + 6*QDxbarp*Qup*Sqr
      (gN))*(6*QDxbarp*QSp*Sqr(gN) + 6*(Sqr(Kappa02) + Sqr(Kappa12) + Sqr(Kappa22)
      )) + (0.8*Sqr(g1) + 6*QDxp*Qup*Sqr(gN))*(6*QDxp*QSp*Sqr(gN) + 6*(Sqr(Kappa20
      ) + Sqr(Kappa21) + Sqr(Kappa22))) + (-1.2*Sqr(g1) + 6*QH2p*Qup*Sqr(gN))*(4*
      QH2p*QSp*Sqr(gN) + 4*(Sqr(Lambda1200) + Sqr(Lambda1201))) + (1.2*Sqr(g1) + 6
      *QH1p*Qup*Sqr(gN))*(4*QH1p*QSp*Sqr(gN) + 4*(Sqr(Lambda1200) + Sqr(Lambda1210
      ))) + (1.2*Sqr(g1) + 6*QH1p*Qup*Sqr(gN))*(4*QH1p*QSp*Sqr(gN) + 4*(Sqr(
      Lambda1201) + Sqr(Lambda1211))) + (-1.2*Sqr(g1) + 6*QH2p*Qup*Sqr(gN))*(4*
      QH2p*QSp*Sqr(gN) + 4*(Sqr(Lambda1210) + Sqr(Lambda1211))) + (1.2*Sqr(g1) + 6
      *QH1p*Qup*Sqr(gN))*(4*QH1p*QSp*Sqr(gN) + 4*Sqr(Lambdax)) + 6*QSp*Qup*Sqr(gN)
      *(6*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11
      ) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 4*(Sqr(
      Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 4*Sqr(
      Lambdax) + 2*Sqr(gN)*Sqr(QSp)) + 12*Power(gN,4)*QSp*Qup*(9*Sqr(Qdp) + 9*Sqr(
      QDxbarp) + 9*Sqr(QDxp) + 3*Sqr(Qep) + 6*Sqr(QH1p) + 6*Sqr(QH2p) + 2*Sqr(
      QHpbarp) + 2*Sqr(QHpp) + 6*Sqr(QLp) + 18*Sqr(QQp) + 3*Sqr(QSp) + 9*Sqr(Qup))
      + 12*QSp*Qup*Sqr(gN)*(1.6*Sqr(g1) + 6*Sqr(gN)*Sqr(Qup)) + 12*QQp*QSp*Sqr(gN
      )*(-0.4*Sqr(g1) + 6*QQp*Qup*Sqr(gN) + 2*Sqr(Yu20)) + 12*QQp*QSp*Sqr(gN)*(
      -0.4*Sqr(g1) + 6*QQp*Qup*Sqr(gN) + 2*Sqr(Yu21)) + 12*QQp*QSp*Sqr(gN)*(-0.4*
      Sqr(g1) + 6*QQp*Qup*Sqr(gN) + 2*Sqr(Yu22)) + 6*QSp*Qup*Sqr(gN)*(1.6*Sqr(g1)
      + 6*Sqr(gN)*Sqr(Qup) + 4*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + (4*QH2p*QSp*
      Sqr(gN) + 4*Sqr(Lambdax))*(-1.2*Sqr(g1) + 6*QH2p*Qup*Sqr(gN) + 6*(Sqr(Yu20)
      + Sqr(Yu21) + Sqr(Yu22)));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dlead_log_one_loop_ms2_dms2() const
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

   double deriv = 8*Power(gN,4)*Power(QSp,4) + (2*QDxp*QSp*Sqr(gN) + 2*(Sqr(
      Kappa00) + Sqr(Kappa01) + Sqr(Kappa02)))*(6*QDxp*QSp*Sqr(gN) + 6*(Sqr(
      Kappa00) + Sqr(Kappa01) + Sqr(Kappa02))) + (2*QDxp*QSp*Sqr(gN) + 2*(Sqr(
      Kappa10) + Sqr(Kappa11) + Sqr(Kappa12)))*(6*QDxp*QSp*Sqr(gN) + 6*(Sqr(
      Kappa10) + Sqr(Kappa11) + Sqr(Kappa12))) + 24*Sqr(Kappa00*Kappa10 + Kappa01*
      Kappa11 + Kappa02*Kappa12) + (2*QDxbarp*QSp*Sqr(gN) + 2*(Sqr(Kappa00) + Sqr(
      Kappa10) + Sqr(Kappa20)))*(6*QDxbarp*QSp*Sqr(gN) + 6*(Sqr(Kappa00) + Sqr(
      Kappa10) + Sqr(Kappa20))) + (2*QDxbarp*QSp*Sqr(gN) + 2*(Sqr(Kappa01) + Sqr(
      Kappa11) + Sqr(Kappa21)))*(6*QDxbarp*QSp*Sqr(gN) + 6*(Sqr(Kappa01) + Sqr(
      Kappa11) + Sqr(Kappa21))) + 24*Sqr(Kappa00*Kappa01 + Kappa10*Kappa11 +
      Kappa20*Kappa21) + (2*QDxbarp*QSp*Sqr(gN) + 2*(Sqr(Kappa02) + Sqr(Kappa12) +
      Sqr(Kappa22)))*(6*QDxbarp*QSp*Sqr(gN) + 6*(Sqr(Kappa02) + Sqr(Kappa12) +
      Sqr(Kappa22))) + (2*QDxp*QSp*Sqr(gN) + 2*(Sqr(Kappa20) + Sqr(Kappa21) + Sqr(
      Kappa22)))*(6*QDxp*QSp*Sqr(gN) + 6*(Sqr(Kappa20) + Sqr(Kappa21) + Sqr(
      Kappa22))) + 24*Sqr(Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22) +
      24*Sqr(Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22) + 24*Sqr(Kappa00
      *Kappa02 + Kappa10*Kappa12 + Kappa20*Kappa22) + 24*Sqr(Kappa01*Kappa02 +
      Kappa11*Kappa12 + Kappa21*Kappa22) + (2*QH2p*QSp*Sqr(gN) + 2*(Sqr(Lambda1200
      ) + Sqr(Lambda1201)))*(4*QH2p*QSp*Sqr(gN) + 4*(Sqr(Lambda1200) + Sqr(
      Lambda1201))) + (2*QH1p*QSp*Sqr(gN) + 2*(Sqr(Lambda1200) + Sqr(Lambda1210)))
      *(4*QH1p*QSp*Sqr(gN) + 4*(Sqr(Lambda1200) + Sqr(Lambda1210))) + (2*QH1p*QSp*
      Sqr(gN) + 2*(Sqr(Lambda1201) + Sqr(Lambda1211)))*(4*QH1p*QSp*Sqr(gN) + 4*(
      Sqr(Lambda1201) + Sqr(Lambda1211))) + (2*QH2p*QSp*Sqr(gN) + 2*(Sqr(
      Lambda1210) + Sqr(Lambda1211)))*(4*QH2p*QSp*Sqr(gN) + 4*(Sqr(Lambda1210) +
      Sqr(Lambda1211))) + 16*Sqr(Lambda1200*Lambda1210 + Lambda1201*Lambda1211) +
      16*Sqr(Lambda1200*Lambda1201 + Lambda1210*Lambda1211) + (2*QH1p*QSp*Sqr(gN)
      + 2*Sqr(Lambdax))*(4*QH1p*QSp*Sqr(gN) + 4*Sqr(Lambdax)) + (2*QH2p*QSp*Sqr(gN
      ) + 2*Sqr(Lambdax))*(4*QH2p*QSp*Sqr(gN) + 4*Sqr(Lambdax)) + 36*Power(gN,4)*
      Sqr(Qdp)*Sqr(QSp) + 12*Power(gN,4)*Sqr(Qep)*Sqr(QSp) + 8*Power(gN,4)*Sqr(
      QHpbarp)*Sqr(QSp) + 8*Power(gN,4)*Sqr(QHpp)*Sqr(QSp) + 24*Power(gN,4)*Sqr(
      QLp)*Sqr(QSp) + 72*Power(gN,4)*Sqr(QQp)*Sqr(QSp) + 12*Kappa00*(2*(Kappa10*(
      Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12) + Kappa20*(Kappa00*
      Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22) + Kappa00*(Sqr(Kappa00) + Sqr(
      Kappa01) + Sqr(Kappa02))) + Kappa00*(-0.26666666666666666*Sqr(g1) -
      5.333333333333333*Sqr(g3) + 3*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) +
      Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) +
      Sqr(Kappa22)) + 2*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr
      (Lambda1211)) + 2*Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QDxbarp) - 2*Sqr(gN)*Sqr(QDxp
      ) - 2*Sqr(gN)*Sqr(QSp))) + 12*Kappa01*(2*(Kappa11*(Kappa00*Kappa10 + Kappa01
      *Kappa11 + Kappa02*Kappa12) + Kappa21*(Kappa00*Kappa20 + Kappa01*Kappa21 +
      Kappa02*Kappa22) + Kappa01*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02))) +
      Kappa01*(-0.26666666666666666*Sqr(g1) - 5.333333333333333*Sqr(g3) + 3*(Sqr(
      Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(
      Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*(Sqr(Lambda1200)
      + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 2*Sqr(Lambdax) - 2*
      Sqr(gN)*Sqr(QDxbarp) - 2*Sqr(gN)*Sqr(QDxp) - 2*Sqr(gN)*Sqr(QSp))) + 12*
      Kappa02*(2*(Kappa12*(Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12) +
      Kappa22*(Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22) + Kappa02*(Sqr
      (Kappa00) + Sqr(Kappa01) + Sqr(Kappa02))) + Kappa02*(-0.26666666666666666*
      Sqr(g1) - 5.333333333333333*Sqr(g3) + 3*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(
      Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(
      Kappa21) + Sqr(Kappa22)) + 2*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(
      Lambda1210) + Sqr(Lambda1211)) + 2*Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QDxbarp) - 2
      *Sqr(gN)*Sqr(QDxp) - 2*Sqr(gN)*Sqr(QSp))) + 12*Kappa10*(2*(Kappa00*(Kappa00*
      Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12) + Kappa20*(Kappa10*Kappa20 +
      Kappa11*Kappa21 + Kappa12*Kappa22) + Kappa10*(Sqr(Kappa10) + Sqr(Kappa11) +
      Sqr(Kappa12))) + Kappa10*(-0.26666666666666666*Sqr(g1) - 5.333333333333333*
      Sqr(g3) + 3*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr
      (Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*(
      Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 2*
      Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QDxbarp) - 2*Sqr(gN)*Sqr(QDxp) - 2*Sqr(gN)*Sqr(
      QSp))) + 12*Kappa11*(2*(Kappa01*(Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02
      *Kappa12) + Kappa21*(Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22) +
      Kappa11*(Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12))) + Kappa11*(
      -0.26666666666666666*Sqr(g1) - 5.333333333333333*Sqr(g3) + 3*(Sqr(Kappa00) +
      Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) +
      Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*(Sqr(Lambda1200) + Sqr(
      Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 2*Sqr(Lambdax) - 2*Sqr(gN
      )*Sqr(QDxbarp) - 2*Sqr(gN)*Sqr(QDxp) - 2*Sqr(gN)*Sqr(QSp))) + 12*Kappa12*(2*
      (Kappa02*(Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12) + Kappa22*(
      Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22) + Kappa12*(Sqr(Kappa10)
      + Sqr(Kappa11) + Sqr(Kappa12))) + Kappa12*(-0.26666666666666666*Sqr(g1) -
      5.333333333333333*Sqr(g3) + 3*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) +
      Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) +
      Sqr(Kappa22)) + 2*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr
      (Lambda1211)) + 2*Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QDxbarp) - 2*Sqr(gN)*Sqr(QDxp
      ) - 2*Sqr(gN)*Sqr(QSp))) + 12*Kappa20*(2*(Kappa00*(Kappa00*Kappa20 + Kappa01
      *Kappa21 + Kappa02*Kappa22) + Kappa10*(Kappa10*Kappa20 + Kappa11*Kappa21 +
      Kappa12*Kappa22) + Kappa20*(Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22))) +
      Kappa20*(-0.26666666666666666*Sqr(g1) - 5.333333333333333*Sqr(g3) + 3*(Sqr(
      Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(
      Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*(Sqr(Lambda1200)
      + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 2*Sqr(Lambdax) - 2*
      Sqr(gN)*Sqr(QDxbarp) - 2*Sqr(gN)*Sqr(QDxp) - 2*Sqr(gN)*Sqr(QSp))) + 12*
      Kappa21*(2*(Kappa01*(Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22) +
      Kappa11*(Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22) + Kappa21*(Sqr
      (Kappa20) + Sqr(Kappa21) + Sqr(Kappa22))) + Kappa21*(-0.26666666666666666*
      Sqr(g1) - 5.333333333333333*Sqr(g3) + 3*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(
      Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(
      Kappa21) + Sqr(Kappa22)) + 2*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(
      Lambda1210) + Sqr(Lambda1211)) + 2*Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QDxbarp) - 2
      *Sqr(gN)*Sqr(QDxp) - 2*Sqr(gN)*Sqr(QSp))) + 12*Kappa22*(2*(Kappa02*(Kappa00*
      Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22) + Kappa12*(Kappa10*Kappa20 +
      Kappa11*Kappa21 + Kappa12*Kappa22) + Kappa22*(Sqr(Kappa20) + Sqr(Kappa21) +
      Sqr(Kappa22))) + Kappa22*(-0.26666666666666666*Sqr(g1) - 5.333333333333333*
      Sqr(g3) + 3*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr
      (Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*(
      Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 2*
      Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QDxbarp) - 2*Sqr(gN)*Sqr(QDxp) - 2*Sqr(gN)*Sqr(
      QSp))) + 8*Lambda1200*(2*(Lambda1210*(Lambda1200*Lambda1210 + Lambda1201*
      Lambda1211) + Lambda1200*(Sqr(Lambda1200) + Sqr(Lambda1201))) + Lambda1200*(
      -0.6*Sqr(g1) - 3*Sqr(g2) + 3*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) +
      Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) +
      Sqr(Kappa22)) + 2*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr
      (Lambda1211)) + 2*Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QH2p) -
      2*Sqr(gN)*Sqr(QSp))) + 8*Lambda1201*(2*(Lambda1211*(Lambda1200*Lambda1210 +
      Lambda1201*Lambda1211) + Lambda1201*(Sqr(Lambda1200) + Sqr(Lambda1201))) +
      Lambda1201*(-0.6*Sqr(g1) - 3*Sqr(g2) + 3*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(
      Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(
      Kappa21) + Sqr(Kappa22)) + 2*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(
      Lambda1210) + Sqr(Lambda1211)) + 2*Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH1p) - 2*
      Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QSp))) + 8*Lambda1210*(2*(Lambda1200*(
      Lambda1200*Lambda1210 + Lambda1201*Lambda1211) + Lambda1210*(Sqr(Lambda1210)
      + Sqr(Lambda1211))) + Lambda1210*(-0.6*Sqr(g1) - 3*Sqr(g2) + 3*(Sqr(Kappa00
      ) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12)
      + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*(Sqr(Lambda1200) + Sqr(
      Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 2*Sqr(Lambdax) - 2*Sqr(gN
      )*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QSp))) + 8*Lambda1211*(2*(
      Lambda1201*(Lambda1200*Lambda1210 + Lambda1201*Lambda1211) + Lambda1211*(Sqr
      (Lambda1210) + Sqr(Lambda1211))) + Lambda1211*(-0.6*Sqr(g1) - 3*Sqr(g2) + 3*
      (Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) +
      Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*(Sqr(
      Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 2*Sqr(
      Lambdax) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QSp)))
      + Sqr(6*(Power(Kappa00,2) + Power(Kappa01,2) + Power(Kappa02,2) + Power(
      Kappa10,2) + Power(Kappa11,2) + Power(Kappa12,2) + Power(Kappa20,2) + Power(
      Kappa21,2) + Power(Kappa22,2)) + 4*(Power(Lambda1200,2) + Power(Lambda1201,2
      ) + Power(Lambda1210,2) + Power(Lambda1211,2)) + 4*Power(Lambdax,2) + 2*
      Power(gN,2)*Power(QSp,2)) + 36*Power(gN,4)*Sqr(QSp)*Sqr(Qup) + 4*Power(gN,4)
      *Sqr(QSp)*(9*Sqr(Qdp) + 9*Sqr(QDxbarp) + 9*Sqr(QDxp) + 3*Sqr(Qep) + 6*Sqr(
      QH1p) + 6*Sqr(QH2p) + 2*Sqr(QHpbarp) + 2*Sqr(QHpp) + 6*Sqr(QLp) + 18*Sqr(QQp
      ) + 3*Sqr(QSp) + 9*Sqr(Qup)) + 8*Lambdax*(4*Power(Lambdax,3) - 0.6*Lambdax*
      Sqr(g1) - 3*Lambdax*Sqr(g2) + 3*Lambdax*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(
      Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(
      Kappa21) + Sqr(Kappa22)) + 2*Lambdax*(Sqr(Lambda1200) + Sqr(Lambda1201) +
      Sqr(Lambda1210) + Sqr(Lambda1211)) - 2*Lambdax*Sqr(gN)*Sqr(QH1p) - 2*Lambdax
      *Sqr(gN)*Sqr(QH2p) - 2*Lambdax*Sqr(gN)*Sqr(QSp) + 3*Lambdax*(Sqr(Yd00) + Sqr
      (Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(
      Yd21) + Sqr(Yd22)) + Lambdax*(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10)
      + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22)) + 3*Lambdax*(
      Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(
      Yu20) + Sqr(Yu21) + Sqr(Yu22)));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dlead_log_one_loop_ms2_dMassB() const
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

   double deriv = 6.4*Kappa00*TKappa00*Sqr(g1) + 6.4*Kappa01*TKappa01*Sqr(g1) +
      6.4*Kappa02*TKappa02*Sqr(g1) + 6.4*Kappa10*TKappa10*Sqr(g1) + 6.4*Kappa11*
      TKappa11*Sqr(g1) + 6.4*Kappa12*TKappa12*Sqr(g1) + 6.4*Kappa20*TKappa20*Sqr(
      g1) + 6.4*Kappa21*TKappa21*Sqr(g1) + 6.4*Kappa22*TKappa22*Sqr(g1) + 9.6*
      Lambda1200*TLambda1200*Sqr(g1) + 9.6*Lambda1201*TLambda1201*Sqr(g1) + 9.6*
      Lambda1210*TLambda1210*Sqr(g1) + 9.6*Lambda1211*TLambda1211*Sqr(g1) + 9.6*
      Lambdax*TLambdax*Sqr(g1) - 19.2*MassB*Qdp*QSp*Sqr(g1)*Sqr(gN) - 57.6*MassB*
      Qep*QSp*Sqr(g1)*Sqr(gN) - 9.6*MassB*QHpbarp*QSp*Sqr(g1)*Sqr(gN) - 9.6*MassB*
      QHpp*QSp*Sqr(g1)*Sqr(gN) - 28.8*MassB*QLp*QSp*Sqr(g1)*Sqr(gN) - 9.6*MassB*
      QQp*QSp*Sqr(g1)*Sqr(gN) - 76.8*MassB*QSp*Qup*Sqr(g1)*Sqr(gN) -
      1.0666666666666667*MassB*Sqr(g1)*(6*QDxp*QSp*Sqr(gN) + 6*(Sqr(Kappa00) + Sqr
      (Kappa01) + Sqr(Kappa02))) - 1.0666666666666667*MassB*Sqr(g1)*(6*QDxp*QSp*
      Sqr(gN) + 6*(Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12))) -
      1.0666666666666667*MassB*Sqr(g1)*(6*QDxbarp*QSp*Sqr(gN) + 6*(Sqr(Kappa00) +
      Sqr(Kappa10) + Sqr(Kappa20))) - 1.0666666666666667*MassB*Sqr(g1)*(6*QDxbarp*
      QSp*Sqr(gN) + 6*(Sqr(Kappa01) + Sqr(Kappa11) + Sqr(Kappa21))) -
      1.0666666666666667*MassB*Sqr(g1)*(6*QDxbarp*QSp*Sqr(gN) + 6*(Sqr(Kappa02) +
      Sqr(Kappa12) + Sqr(Kappa22))) - 1.0666666666666667*MassB*Sqr(g1)*(6*QDxp*QSp
      *Sqr(gN) + 6*(Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22))) - 2.4*MassB*Sqr(
      g1)*(4*QH2p*QSp*Sqr(gN) + 4*(Sqr(Lambda1200) + Sqr(Lambda1201))) - 2.4*MassB
      *Sqr(g1)*(4*QH1p*QSp*Sqr(gN) + 4*(Sqr(Lambda1200) + Sqr(Lambda1210))) - 2.4*
      MassB*Sqr(g1)*(4*QH1p*QSp*Sqr(gN) + 4*(Sqr(Lambda1201) + Sqr(Lambda1211))) -
      2.4*MassB*Sqr(g1)*(4*QH2p*QSp*Sqr(gN) + 4*(Sqr(Lambda1210) + Sqr(Lambda1211
      ))) - 2.4*MassB*Sqr(g1)*(4*QH1p*QSp*Sqr(gN) + 4*Sqr(Lambdax)) - 2.4*MassB*
      Sqr(g1)*(4*QH2p*QSp*Sqr(gN) + 4*Sqr(Lambdax));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dlead_log_one_loop_ms2_dMassWB() const
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

   double deriv = 48*Lambda1200*TLambda1200*Sqr(g2) + 48*Lambda1201*TLambda1201
      *Sqr(g2) + 48*Lambda1210*TLambda1210*Sqr(g2) + 48*Lambda1211*TLambda1211*Sqr
      (g2) + 48*Lambdax*TLambdax*Sqr(g2) - 48*MassWB*QHpbarp*QSp*Sqr(g2)*Sqr(gN) -
      48*MassWB*QHpp*QSp*Sqr(g2)*Sqr(gN) - 144*MassWB*QLp*QSp*Sqr(g2)*Sqr(gN) -
      432*MassWB*QQp*QSp*Sqr(g2)*Sqr(gN) - 12*MassWB*Sqr(g2)*(4*QH2p*QSp*Sqr(gN) +
      4*(Sqr(Lambda1200) + Sqr(Lambda1201))) - 12*MassWB*Sqr(g2)*(4*QH1p*QSp*Sqr(
      gN) + 4*(Sqr(Lambda1200) + Sqr(Lambda1210))) - 12*MassWB*Sqr(g2)*(4*QH1p*QSp
      *Sqr(gN) + 4*(Sqr(Lambda1201) + Sqr(Lambda1211))) - 12*MassWB*Sqr(g2)*(4*
      QH2p*QSp*Sqr(gN) + 4*(Sqr(Lambda1210) + Sqr(Lambda1211))) - 12*MassWB*Sqr(g2
      )*(4*QH1p*QSp*Sqr(gN) + 4*Sqr(Lambdax)) - 12*MassWB*Sqr(g2)*(4*QH2p*QSp*Sqr(
      gN) + 4*Sqr(Lambdax));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dlead_log_one_loop_ms2_dMassG() const
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

   double deriv = 128*Kappa00*TKappa00*Sqr(g3) + 128*Kappa01*TKappa01*Sqr(g3) +
      128*Kappa02*TKappa02*Sqr(g3) + 128*Kappa10*TKappa10*Sqr(g3) + 128*Kappa11*
      TKappa11*Sqr(g3) + 128*Kappa12*TKappa12*Sqr(g3) + 128*Kappa20*TKappa20*Sqr(
      g3) + 128*Kappa21*TKappa21*Sqr(g3) + 128*Kappa22*TKappa22*Sqr(g3) - 384*
      MassG*Qdp*QSp*Sqr(g3)*Sqr(gN) - 768*MassG*QQp*QSp*Sqr(g3)*Sqr(gN) - 384*
      MassG*QSp*Qup*Sqr(g3)*Sqr(gN) - 21.333333333333332*MassG*Sqr(g3)*(6*QDxp*QSp
      *Sqr(gN) + 6*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02))) -
      21.333333333333332*MassG*Sqr(g3)*(6*QDxp*QSp*Sqr(gN) + 6*(Sqr(Kappa10) + Sqr
      (Kappa11) + Sqr(Kappa12))) - 21.333333333333332*MassG*Sqr(g3)*(6*QDxbarp*QSp
      *Sqr(gN) + 6*(Sqr(Kappa00) + Sqr(Kappa10) + Sqr(Kappa20))) -
      21.333333333333332*MassG*Sqr(g3)*(6*QDxbarp*QSp*Sqr(gN) + 6*(Sqr(Kappa01) +
      Sqr(Kappa11) + Sqr(Kappa21))) - 21.333333333333332*MassG*Sqr(g3)*(6*QDxbarp*
      QSp*Sqr(gN) + 6*(Sqr(Kappa02) + Sqr(Kappa12) + Sqr(Kappa22))) -
      21.333333333333332*MassG*Sqr(g3)*(6*QDxp*QSp*Sqr(gN) + 6*(Sqr(Kappa20) + Sqr
      (Kappa21) + Sqr(Kappa22)));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dlead_log_one_loop_ms2_dMassBp() const
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

   double deriv = -288*Power(gN,4)*MassBp*Power(Qdp,3)*QSp - 96*Power(gN,4)*
      MassBp*Power(Qep,3)*QSp - 64*Power(gN,4)*MassBp*Power(QHpbarp,3)*QSp - 64*
      Power(gN,4)*MassBp*Power(QHpp,3)*QSp - 192*Power(gN,4)*MassBp*Power(QLp,3)*
      QSp - 576*Power(gN,4)*MassBp*Power(QQp,3)*QSp - 64*Power(gN,4)*MassBp*Power(
      QSp,4) - 288*Power(gN,4)*MassBp*QSp*Power(Qup,3) - 16*MassBp*Sqr(gN)*(6*
      QDxbarp*QSp*Sqr(gN) + 6*(Sqr(Kappa00) + Sqr(Kappa10) + Sqr(Kappa20)))*Sqr(
      QDxbarp) - 16*MassBp*Sqr(gN)*(6*QDxbarp*QSp*Sqr(gN) + 6*(Sqr(Kappa01) + Sqr(
      Kappa11) + Sqr(Kappa21)))*Sqr(QDxbarp) - 16*MassBp*Sqr(gN)*(6*QDxbarp*QSp*
      Sqr(gN) + 6*(Sqr(Kappa02) + Sqr(Kappa12) + Sqr(Kappa22)))*Sqr(QDxbarp) - 16*
      MassBp*Sqr(gN)*(6*QDxp*QSp*Sqr(gN) + 6*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(
      Kappa02)))*Sqr(QDxp) - 16*MassBp*Sqr(gN)*(6*QDxp*QSp*Sqr(gN) + 6*(Sqr(
      Kappa10) + Sqr(Kappa11) + Sqr(Kappa12)))*Sqr(QDxp) - 16*MassBp*Sqr(gN)*(6*
      QDxp*QSp*Sqr(gN) + 6*(Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)))*Sqr(QDxp)
      - 16*MassBp*Sqr(gN)*(4*QH1p*QSp*Sqr(gN) + 4*(Sqr(Lambda1200) + Sqr(
      Lambda1210)))*Sqr(QH1p) - 16*MassBp*Sqr(gN)*(4*QH1p*QSp*Sqr(gN) + 4*(Sqr(
      Lambda1201) + Sqr(Lambda1211)))*Sqr(QH1p) - 16*MassBp*Sqr(gN)*(4*QH1p*QSp*
      Sqr(gN) + 4*Sqr(Lambdax))*Sqr(QH1p) - 16*MassBp*Sqr(gN)*(4*QH2p*QSp*Sqr(gN)
      + 4*(Sqr(Lambda1200) + Sqr(Lambda1201)))*Sqr(QH2p) - 16*MassBp*Sqr(gN)*(4*
      QH2p*QSp*Sqr(gN) + 4*(Sqr(Lambda1210) + Sqr(Lambda1211)))*Sqr(QH2p) - 16*
      MassBp*Sqr(gN)*(4*QH2p*QSp*Sqr(gN) + 4*Sqr(Lambdax))*Sqr(QH2p) - 16*MassBp*
      Sqr(gN)*Sqr(QSp)*(6*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(
      Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(
      Kappa22)) + 4*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(
      Lambda1211)) + 4*Sqr(Lambdax) + 2*Sqr(gN)*Sqr(QSp)) + 12*Kappa00*TKappa00*(4
      *Sqr(gN)*Sqr(QDxbarp) + 4*Sqr(gN)*Sqr(QDxp) + 4*Sqr(gN)*Sqr(QSp)) + 12*
      Kappa01*TKappa01*(4*Sqr(gN)*Sqr(QDxbarp) + 4*Sqr(gN)*Sqr(QDxp) + 4*Sqr(gN)*
      Sqr(QSp)) + 12*Kappa02*TKappa02*(4*Sqr(gN)*Sqr(QDxbarp) + 4*Sqr(gN)*Sqr(QDxp
      ) + 4*Sqr(gN)*Sqr(QSp)) + 12*Kappa10*TKappa10*(4*Sqr(gN)*Sqr(QDxbarp) + 4*
      Sqr(gN)*Sqr(QDxp) + 4*Sqr(gN)*Sqr(QSp)) + 12*Kappa11*TKappa11*(4*Sqr(gN)*Sqr
      (QDxbarp) + 4*Sqr(gN)*Sqr(QDxp) + 4*Sqr(gN)*Sqr(QSp)) + 12*Kappa12*TKappa12*
      (4*Sqr(gN)*Sqr(QDxbarp) + 4*Sqr(gN)*Sqr(QDxp) + 4*Sqr(gN)*Sqr(QSp)) + 12*
      Kappa20*TKappa20*(4*Sqr(gN)*Sqr(QDxbarp) + 4*Sqr(gN)*Sqr(QDxp) + 4*Sqr(gN)*
      Sqr(QSp)) + 12*Kappa21*TKappa21*(4*Sqr(gN)*Sqr(QDxbarp) + 4*Sqr(gN)*Sqr(QDxp
      ) + 4*Sqr(gN)*Sqr(QSp)) + 12*Kappa22*TKappa22*(4*Sqr(gN)*Sqr(QDxbarp) + 4*
      Sqr(gN)*Sqr(QDxp) + 4*Sqr(gN)*Sqr(QSp)) + 8*Lambda1200*TLambda1200*(4*Sqr(gN
      )*Sqr(QH1p) + 4*Sqr(gN)*Sqr(QH2p) + 4*Sqr(gN)*Sqr(QSp)) + 8*Lambda1201*
      TLambda1201*(4*Sqr(gN)*Sqr(QH1p) + 4*Sqr(gN)*Sqr(QH2p) + 4*Sqr(gN)*Sqr(QSp))
      + 8*Lambda1210*TLambda1210*(4*Sqr(gN)*Sqr(QH1p) + 4*Sqr(gN)*Sqr(QH2p) + 4*
      Sqr(gN)*Sqr(QSp)) + 8*Lambda1211*TLambda1211*(4*Sqr(gN)*Sqr(QH1p) + 4*Sqr(gN
      )*Sqr(QH2p) + 4*Sqr(gN)*Sqr(QSp)) + 8*TLambdax*(4*Lambdax*Sqr(gN)*Sqr(QH1p)
      + 4*Lambdax*Sqr(gN)*Sqr(QH2p) + 4*Lambdax*Sqr(gN)*Sqr(QSp)) - 96*Power(gN,4)
      *MassBp*Sqr(QSp)*(9*Sqr(Qdp) + 9*Sqr(QDxbarp) + 9*Sqr(QDxp) + 3*Sqr(Qep) + 6
      *Sqr(QH1p) + 6*Sqr(QH2p) + 2*Sqr(QHpbarp) + 2*Sqr(QHpp) + 6*Sqr(QLp) + 18*
      Sqr(QQp) + 3*Sqr(QSp) + 9*Sqr(Qup));


   return twoLoop*deriv;
}

} // namespace flexiblesusy
