#include "lowE6SSM_two_scale_tuning_calculator.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

// one loop
double lowE6SSM_tuning_calculator::deriv_dbeta_one_loop_mq222_dLambdax() const
{

   double deriv = 0;


   return oneOver16PiSqr*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_one_loop_mq222_dALambdax() const
{

   double deriv = 0;


   return oneOver16PiSqr*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_one_loop_mq222_dAYu22() const
{

   const double Yu22 = model.get_Yu(2,2);

   const double AYu22 = get_AYu22();

   double deriv = 4*AYu22*Sqr(Yu22);


   return oneOver16PiSqr*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_one_loop_mq222_dmq222() const
{
   const lowE6SSM_input_parameters inputs = model.get_input();

   const auto QQp = inputs.QQp;
   
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
   const double gN = model.get_gN();

   double deriv = 0.2*Sqr(g1) + 12*Sqr(gN)*Sqr(QQp) + 2*Sqr(Yd02) + 2*Sqr(Yd12)
      + 2*Sqr(Yd22) + 2*Sqr(Yu02) + 2*Sqr(Yu12) + 2*Sqr(Yu22);


   return oneOver16PiSqr*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_one_loop_mq222_dmHd2() const
{
   const lowE6SSM_input_parameters inputs = model.get_input();

   const auto QQp = inputs.QQp;
   const auto QH1p = inputs.QH1p;
   
   const double Yd02 = model.get_Yd(0,2);
   const double Yd12 = model.get_Yd(1,2);
   const double Yd22 = model.get_Yd(2,2);
   const double g1 = model.get_g1();
   const double gN = model.get_gN();

   double deriv = -0.2*Sqr(g1) + 4*QH1p*QQp*Sqr(gN) + 2*(Sqr(Yd02) + Sqr(Yd12)
      + Sqr(Yd22));


   return oneOver16PiSqr*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_one_loop_mq222_dmHu2() const
{
   const lowE6SSM_input_parameters inputs = model.get_input();

   const auto QQp = inputs.QQp;
   const auto QH2p = inputs.QH2p;
   
   const double Yu02 = model.get_Yu(0,2);
   const double Yu12 = model.get_Yu(1,2);
   const double Yu22 = model.get_Yu(2,2);
   const double g1 = model.get_g1();
   const double gN = model.get_gN();

   double deriv = 0.2*Sqr(g1) + 4*QH2p*QQp*Sqr(gN) + 2*(Sqr(Yu02) + Sqr(Yu12) +
      Sqr(Yu22));


   return oneOver16PiSqr*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_one_loop_mq222_dmu222() const
{
   const lowE6SSM_input_parameters inputs = model.get_input();

   const auto QQp = inputs.QQp;
   const auto Qup = inputs.Qup;
   
   const double Yu22 = model.get_Yu(2,2);
   const double g1 = model.get_g1();
   const double gN = model.get_gN();

   double deriv = -0.4*Sqr(g1) + 6*QQp*Qup*Sqr(gN) + 2*Sqr(Yu22);


   return oneOver16PiSqr*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_one_loop_mq222_dms2() const
{
   const lowE6SSM_input_parameters inputs = model.get_input();

   const auto QQp = inputs.QQp;
   const auto QSp = inputs.QSp;
   
   const double gN = model.get_gN();

   double deriv = 2*QQp*QSp*Sqr(gN);


   return oneOver16PiSqr*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_one_loop_mq222_dMassB() const
{

   const double g1 = model.get_g1();

   const double MassB = model.get_MassB();

   double deriv = -0.26666666666666666*MassB*Sqr(g1);


   return oneOver16PiSqr*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_one_loop_mq222_dMassWB() const
{

   const double g2 = model.get_g2();

   const double MassWB = model.get_MassWB();

   double deriv = -12*MassWB*Sqr(g2);


   return oneOver16PiSqr*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_one_loop_mq222_dMassG() const
{

   const double g3 = model.get_g3();

   const double MassG = model.get_MassG();

   double deriv = -21.333333333333332*MassG*Sqr(g3);


   return oneOver16PiSqr*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_one_loop_mq222_dMassBp() const
{
   const lowE6SSM_input_parameters inputs = model.get_input();

   const auto QQp = inputs.QQp;
   
   const double gN = model.get_gN();

   const double MassBp = model.get_MassBp();

   double deriv = -16*MassBp*Sqr(gN)*Sqr(QQp);


   return oneOver16PiSqr*deriv;
}

// two loop
double lowE6SSM_tuning_calculator::deriv_dbeta_two_loop_mq222_dLambdax() const
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

   double deriv = -8*ALambdax*Lambdax*(TYd02*Yd02 + TYd12*Yd12 + TYd22*Yd22) -
      4*Lambdax*(Yd02*(md200*Yd02 + md210*Yd12 + md220*Yd22) + Yd12*(md201*Yd02 +
      md211*Yd12 + md221*Yd22) + Yd22*(md202*Yd02 + md212*Yd12 + md222*Yd22)) - 2*
      Lambdax*(Yd02*(mq220*Yd00 + mq221*Yd01 + mq222*Yd02) + Yd12*(mq220*Yd10 +
      mq221*Yd11 + mq222*Yd12) + Yd22*(mq220*Yd20 + mq221*Yd21 + mq222*Yd22)) - 8*
      ALambdax*Lambdax*(TYu02*Yu02 + TYu12*Yu12 + TYu22*Yu22) - 2*Lambdax*(Yu02*(
      mq220*Yu00 + mq221*Yu01 + mq222*Yu02) + Yu12*(mq220*Yu10 + mq221*Yu11 +
      mq222*Yu12) + Yu22*(mq220*Yu20 + mq221*Yu21 + mq222*Yu22)) - 4*Lambdax*(Yu02
      *(mu200*Yu02 + mu210*Yu12 + mu220*Yu22) + Yu12*(mu201*Yu02 + mu211*Yu12 +
      mu221*Yu22) + Yu22*(mu202*Yu02 + mu212*Yu12 + mu222*Yu22)) + 0.8*Lambdax*(
      mHd2 - mHu2)*Sqr(g1) - 16*Lambdax*QQp*(mHd2*QH1p + mHu2*QH2p + ms2*QSp)*Sqr(
      gN) - 4*Lambdax*(Sqr(TYd02) + Sqr(TYd12) + Sqr(TYd22)) - 4*Lambdax*(Sqr(
      TYu02) + Sqr(TYu12) + Sqr(TYu22)) - 8*Lambdax*mHd2*(Sqr(Yd02) + Sqr(Yd12) +
      Sqr(Yd22)) - 4*Lambdax*mHu2*(Sqr(Yd02) + Sqr(Yd12) + Sqr(Yd22)) - 4*Lambdax*
      ms2*(Sqr(Yd02) + Sqr(Yd12) + Sqr(Yd22)) - 4*Lambdax*Sqr(ALambdax)*(Sqr(Yd02)
      + Sqr(Yd12) + Sqr(Yd22)) - 2*Lambdax*(mq202*(Yd00*Yd02 + Yd10*Yd12 + Yd20*
      Yd22) + mq212*(Yd01*Yd02 + Yd11*Yd12 + Yd21*Yd22) + mq222*(Sqr(Yd02) + Sqr(
      Yd12) + Sqr(Yd22))) - 4*Lambdax*mHd2*(Sqr(Yu02) + Sqr(Yu12) + Sqr(Yu22)) - 8
      *Lambdax*mHu2*(Sqr(Yu02) + Sqr(Yu12) + Sqr(Yu22)) - 4*Lambdax*ms2*(Sqr(Yu02)
      + Sqr(Yu12) + Sqr(Yu22)) - 4*Lambdax*Sqr(ALambdax)*(Sqr(Yu02) + Sqr(Yu12) +
      Sqr(Yu22)) - 2*Lambdax*(mq202*(Yu00*Yu02 + Yu10*Yu12 + Yu20*Yu22) + mq212*(
      Yu01*Yu02 + Yu11*Yu12 + Yu21*Yu22) + mq222*(Sqr(Yu02) + Sqr(Yu12) + Sqr(Yu22
      )));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_two_loop_mq222_dALambdax() const
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

   double deriv = -4*(TYd02*Yd02 + TYd12*Yd12 + TYd22*Yd22)*Sqr(Lambdax) - 4*(
      TYu02*Yu02 + TYu12*Yu12 + TYu22*Yu22)*Sqr(Lambdax) - 4*ALambdax*Sqr(Lambdax)
      *(Sqr(Yd02) + Sqr(Yd12) + Sqr(Yd22)) - 4*ALambdax*Sqr(Lambdax)*(Sqr(Yu02) +
      Sqr(Yu12) + Sqr(Yu22));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_two_loop_mq222_dAYu22() const
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

   double deriv = -4*Lambdax*TLambdax*Sqr(Yu22) + 3.2*AYu22*Sqr(g1)*Sqr(Yu22) -
      3.2*MassB*Sqr(g1)*Sqr(Yu22) - 4*AYu22*Sqr(Lambdax)*Sqr(Yu22) + 8*AYu22*Sqr(
      gN)*Sqr(QH2p)*Sqr(Yu22) - 4*MassBp*Sqr(gN)*Sqr(QH2p)*Sqr(Yu22) - 8*AYu22*Sqr
      (gN)*Sqr(QQp)*Sqr(Yu22) + 4*MassBp*Sqr(gN)*Sqr(QQp)*Sqr(Yu22) + 8*AYu22*Sqr(
      gN)*Sqr(Qup)*Sqr(Yu22) - 4*MassBp*Sqr(gN)*Sqr(Qup)*Sqr(Yu22) - 4*MassBp*Sqr(
      gN)*(Sqr(QH2p) - Sqr(QQp) + Sqr(Qup))*Sqr(Yu22) - 12*AYu22*Sqr(Yu22)*(Sqr(
      Yu02) + Sqr(Yu12) + Sqr(Yu22)) - 12*AYu22*Sqr(Yu22)*(Sqr(Yu00) + Sqr(Yu01) +
      Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr
      (Yu22)) - 12*Sqr(Yu22)*(TYu02*Yu02 + TYu12*Yu12 + AYu22*Sqr(Yu22)) - 12*Sqr(
      Yu22)*(TYu00*Yu00 + TYu01*Yu01 + TYu02*Yu02 + TYu10*Yu10 + TYu11*Yu11 +
      TYu12*Yu12 + TYu20*Yu20 + TYu21*Yu21 + AYu22*Sqr(Yu22)) - 4*(Yu22*(TYu20*
      Yu20*Yu22 + TYu21*Yu21*Yu22 + 2*AYu22*Power(Yu22,3)) + Yu02*(TYu20*Yu00*Yu22
      + TYu21*Yu01*Yu22 + 2*AYu22*Yu02*Sqr(Yu22)) + Yu12*(TYu20*Yu10*Yu22 + TYu21
      *Yu11*Yu22 + 2*AYu22*Yu12*Sqr(Yu22))) - 4*(AYu22*Sqr(Yu22)*(Sqr(Yu02) + Sqr(
      Yu12) + Sqr(Yu22)) + Yu22*(TYu20*(Yu00*Yu02 + Yu10*Yu12 + Yu20*Yu22) + TYu21
      *(Yu01*Yu02 + Yu11*Yu12 + Yu21*Yu22) + AYu22*Yu22*(Sqr(Yu02) + Sqr(Yu12) +
      Sqr(Yu22)))) - 4*(TYu02*Yu02*Sqr(Yu22) + TYu12*Yu12*Sqr(Yu22) + Yu22*(AYu22*
      Power(Yu22,3) + Yu22*(TYu02*Yu02 + TYu12*Yu12 + AYu22*Sqr(Yu22)))) - 4*(
      AYu22*Yu22*(Power(Yu22,3) + Yu22*Sqr(Yu20) + Yu22*Sqr(Yu21)) + TYu02*(Yu00*
      Yu20*Yu22 + Yu01*Yu21*Yu22 + Yu02*Sqr(Yu22)) + TYu12*(Yu10*Yu20*Yu22 + Yu11*
      Yu21*Yu22 + Yu12*Sqr(Yu22)) + Yu22*(Yu20*(TYu02*Yu00 + TYu12*Yu10 + AYu22*
      Yu20*Yu22) + Yu21*(TYu02*Yu01 + TYu12*Yu11 + AYu22*Yu21*Yu22) + Yu22*(TYu02*
      Yu02 + TYu12*Yu12 + AYu22*Sqr(Yu22))));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_two_loop_mq222_dmq222() const
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

   double deriv = 0.013333333333333334*Power(g1,4) + 9*Power(g2,4) +
      10.666666666666666*Power(g3,4) + 48*Power(gN,4)*Power(QQp,4) + 1.6*Sqr(g1)*
      Sqr(gN)*Sqr(QQp) + 0.8*Sqr(g1)*(Sqr(Yd02) + Sqr(Yd12) + Sqr(Yd22)) - 2*Sqr(
      Lambdax)*(Sqr(Yd02) + Sqr(Yd12) + Sqr(Yd22)) + 4*Sqr(gN)*Sqr(Qdp)*(Sqr(Yd02)
      + Sqr(Yd12) + Sqr(Yd22)) + 4*Sqr(gN)*Sqr(QH1p)*(Sqr(Yd02) + Sqr(Yd12) + Sqr
      (Yd22)) - 4*Sqr(gN)*Sqr(QQp)*(Sqr(Yd02) + Sqr(Yd12) + Sqr(Yd22)) - 6*(Sqr(
      Yd02) + Sqr(Yd12) + Sqr(Yd22))*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10
      ) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) - 4*(Sqr(Yd02
      )*(Sqr(Yd02) + Sqr(Yd12) + Sqr(Yd22)) + Sqr(Yd12)*(Sqr(Yd02) + Sqr(Yd12) +
      Sqr(Yd22)) + Sqr(Yd22)*(Sqr(Yd02) + Sqr(Yd12) + Sqr(Yd22))) - 4*(Yd02*(Yd00*
      (Yd00*Yd02 + Yd10*Yd12 + Yd20*Yd22) + Yd01*(Yd01*Yd02 + Yd11*Yd12 + Yd21*
      Yd22) + Yd02*(Sqr(Yd02) + Sqr(Yd12) + Sqr(Yd22))) + Yd12*(Yd10*(Yd00*Yd02 +
      Yd10*Yd12 + Yd20*Yd22) + Yd11*(Yd01*Yd02 + Yd11*Yd12 + Yd21*Yd22) + Yd12*(
      Sqr(Yd02) + Sqr(Yd12) + Sqr(Yd22))) + Yd22*(Yd20*(Yd00*Yd02 + Yd10*Yd12 +
      Yd20*Yd22) + Yd21*(Yd01*Yd02 + Yd11*Yd12 + Yd21*Yd22) + Yd22*(Sqr(Yd02) +
      Sqr(Yd12) + Sqr(Yd22)))) - 6*Sqr(Power(Yd02,2) + Power(Yd12,2) + Power(Yd22,
      2)) - 2*(Sqr(Yd02) + Sqr(Yd12) + Sqr(Yd22))*(Sqr(Ye00) + Sqr(Ye01) + Sqr(
      Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22
      )) + 1.6*Sqr(g1)*(Sqr(Yu02) + Sqr(Yu12) + Sqr(Yu22)) - 2*Sqr(Lambdax)*(Sqr(
      Yu02) + Sqr(Yu12) + Sqr(Yu22)) + 4*Sqr(gN)*Sqr(QH2p)*(Sqr(Yu02) + Sqr(Yu12)
      + Sqr(Yu22)) - 4*Sqr(gN)*Sqr(QQp)*(Sqr(Yu02) + Sqr(Yu12) + Sqr(Yu22)) + 4*
      Sqr(gN)*Sqr(Qup)*(Sqr(Yu02) + Sqr(Yu12) + Sqr(Yu22)) - 6*(Sqr(Yu02) + Sqr(
      Yu12) + Sqr(Yu22))*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11
      ) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)) + 0.013333333333333334*
      Sqr(g1)*(Sqr(g1) + 45*Sqr(g2) + 80*Sqr(g3) + 60*Sqr(gN)*Sqr(QQp) - 30*(Sqr(
      Yd02) + Sqr(Yd12) + Sqr(Yd22)) - 30*(Sqr(Yu02) + Sqr(Yu12) + Sqr(Yu22))) +
      0.8*QQp*Sqr(gN)*(QQp*Sqr(g1) + 45*QQp*Sqr(g2) + 80*QQp*Sqr(g3) + 60*Power(
      QQp,3)*Sqr(gN) - 30*QQp*(Sqr(Yd02) + Sqr(Yd12) + Sqr(Yd22)) - 30*QQp*(Sqr(
      Yu02) + Sqr(Yu12) + Sqr(Yu22))) - 4*(Sqr(Yu02)*(Sqr(Yu02) + Sqr(Yu12) + Sqr(
      Yu22)) + Sqr(Yu12)*(Sqr(Yu02) + Sqr(Yu12) + Sqr(Yu22)) + Sqr(Yu22)*(Sqr(Yu02
      ) + Sqr(Yu12) + Sqr(Yu22))) - 4*(Yu02*(Yu00*(Yu00*Yu02 + Yu10*Yu12 + Yu20*
      Yu22) + Yu01*(Yu01*Yu02 + Yu11*Yu12 + Yu21*Yu22) + Yu02*(Sqr(Yu02) + Sqr(
      Yu12) + Sqr(Yu22))) + Yu12*(Yu10*(Yu00*Yu02 + Yu10*Yu12 + Yu20*Yu22) + Yu11*
      (Yu01*Yu02 + Yu11*Yu12 + Yu21*Yu22) + Yu12*(Sqr(Yu02) + Sqr(Yu12) + Sqr(Yu22
      ))) + Yu22*(Yu20*(Yu00*Yu02 + Yu10*Yu12 + Yu20*Yu22) + Yu21*(Yu01*Yu02 +
      Yu11*Yu12 + Yu21*Yu22) + Yu22*(Sqr(Yu02) + Sqr(Yu12) + Sqr(Yu22)))) - 6*Sqr(
      Power(Yu02,2) + Power(Yu12,2) + Power(Yu22,2));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_two_loop_mq222_dmHd2() const
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

   double deriv = 0.04*Power(g1,4) + 3*Power(g2,4) - 1.6*QH1p*QQp*Sqr(g1)*Sqr(
      gN) + 16*Power(gN,4)*Sqr(QH1p)*Sqr(QQp) + 0.8*Sqr(g1)*(Sqr(Yd02) + Sqr(Yd12)
      + Sqr(Yd22)) - 4*Sqr(Lambdax)*(Sqr(Yd02) + Sqr(Yd12) + Sqr(Yd22)) + 4*Sqr(
      gN)*Sqr(Qdp)*(Sqr(Yd02) + Sqr(Yd12) + Sqr(Yd22)) + 4*Sqr(gN)*Sqr(QH1p)*(Sqr(
      Yd02) + Sqr(Yd12) + Sqr(Yd22)) - 4*Sqr(gN)*Sqr(QQp)*(Sqr(Yd02) + Sqr(Yd12) +
      Sqr(Yd22)) - 12*(Sqr(Yd02) + Sqr(Yd12) + Sqr(Yd22))*(Sqr(Yd00) + Sqr(Yd01)
      + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) +
      Sqr(Yd22)) - 8*(Yd02*(Yd00*(Yd00*Yd02 + Yd10*Yd12 + Yd20*Yd22) + Yd01*(Yd01*
      Yd02 + Yd11*Yd12 + Yd21*Yd22) + Yd02*(Sqr(Yd02) + Sqr(Yd12) + Sqr(Yd22))) +
      Yd12*(Yd10*(Yd00*Yd02 + Yd10*Yd12 + Yd20*Yd22) + Yd11*(Yd01*Yd02 + Yd11*Yd12
      + Yd21*Yd22) + Yd12*(Sqr(Yd02) + Sqr(Yd12) + Sqr(Yd22))) + Yd22*(Yd20*(Yd00
      *Yd02 + Yd10*Yd12 + Yd20*Yd22) + Yd21*(Yd01*Yd02 + Yd11*Yd12 + Yd21*Yd22) +
      Yd22*(Sqr(Yd02) + Sqr(Yd12) + Sqr(Yd22)))) - 4*(Sqr(Yd02) + Sqr(Yd12) + Sqr(
      Yd22))*(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12
      ) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22)) + 0.013333333333333334*Sqr(g1)*(-9*
      Sqr(g1) - 45*Sqr(g2) + 30*Sqr(Lambdax) - 60*Sqr(gN)*Sqr(QH1p) + 90*(Sqr(Yd00
      ) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) +
      Sqr(Yd21) + Sqr(Yd22)) + 30*(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) +
      Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22))) + 0.8*QQp*Sqr(
      gN)*(3*QH1p*Sqr(g1) + 15*QH1p*Sqr(g2) + 20*Power(QH1p,3)*Sqr(gN) - 10*QH1p*
      Sqr(Lambdax) - 30*QH1p*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(
      Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) - 10*QH1p*(Sqr(Ye00)
      + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) +
      Sqr(Ye21) + Sqr(Ye22))) - 2*Sqr(Lambdax)*(Sqr(Yu02) + Sqr(Yu12) + Sqr(Yu22))
      ;


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_two_loop_mq222_dmHu2() const
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

   double deriv = 0.04*Power(g1,4) + 3*Power(g2,4) + 1.6*QH2p*QQp*Sqr(g1)*Sqr(
      gN) + 16*Power(gN,4)*Sqr(QH2p)*Sqr(QQp) - 2*Sqr(Lambdax)*(Sqr(Yd02) + Sqr(
      Yd12) + Sqr(Yd22)) + 1.6*Sqr(g1)*(Sqr(Yu02) + Sqr(Yu12) + Sqr(Yu22)) - 4*Sqr
      (Lambdax)*(Sqr(Yu02) + Sqr(Yu12) + Sqr(Yu22)) + 4*Sqr(gN)*Sqr(QH2p)*(Sqr(
      Yu02) + Sqr(Yu12) + Sqr(Yu22)) - 4*Sqr(gN)*Sqr(QQp)*(Sqr(Yu02) + Sqr(Yu12) +
      Sqr(Yu22)) + 4*Sqr(gN)*Sqr(Qup)*(Sqr(Yu02) + Sqr(Yu12) + Sqr(Yu22)) - 12*(
      Sqr(Yu02) + Sqr(Yu12) + Sqr(Yu22))*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(
      Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)) +
      0.013333333333333334*Sqr(g1)*(9*Sqr(g1) + 45*Sqr(g2) - 30*Sqr(Lambdax) + 60*
      Sqr(gN)*Sqr(QH2p) - 90*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(
      Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + 0.8*QQp*Sqr(gN)*(3
      *QH2p*Sqr(g1) + 15*QH2p*Sqr(g2) + 20*Power(QH2p,3)*Sqr(gN) - 10*QH2p*Sqr(
      Lambdax) - 30*QH2p*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11
      ) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) - 8*(Yu02*(Yu00*(Yu00*
      Yu02 + Yu10*Yu12 + Yu20*Yu22) + Yu01*(Yu01*Yu02 + Yu11*Yu12 + Yu21*Yu22) +
      Yu02*(Sqr(Yu02) + Sqr(Yu12) + Sqr(Yu22))) + Yu12*(Yu10*(Yu00*Yu02 + Yu10*
      Yu12 + Yu20*Yu22) + Yu11*(Yu01*Yu02 + Yu11*Yu12 + Yu21*Yu22) + Yu12*(Sqr(
      Yu02) + Sqr(Yu12) + Sqr(Yu22))) + Yu22*(Yu20*(Yu00*Yu02 + Yu10*Yu12 + Yu20*
      Yu22) + Yu21*(Yu01*Yu02 + Yu11*Yu12 + Yu21*Yu22) + Yu22*(Sqr(Yu02) + Sqr(
      Yu12) + Sqr(Yu22))));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_two_loop_mq222_dmu222() const
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

   double deriv = 0.10666666666666667*Power(g1,4) + 5.333333333333333*Power(g3,
      4) - 3.2*QQp*Qup*Sqr(g1)*Sqr(gN) + 24*Power(gN,4)*Sqr(QQp)*Sqr(Qup) + 1.6*
      Sqr(g1)*Sqr(Yu22) - 2*Sqr(Lambdax)*Sqr(Yu22) + 4*Sqr(gN)*Sqr(QH2p)*Sqr(Yu22)
      - 4*Sqr(gN)*Sqr(QQp)*Sqr(Yu22) + 4*Sqr(gN)*Sqr(Qup)*Sqr(Yu22) - 6*(Sqr(Yu02
      ) + Sqr(Yu12) + Sqr(Yu22))*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)) - 6*Sqr(Yu22)
      *(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) +
      Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)) - 4*Yu22*(Yu20*(Yu00*Yu02 + Yu10*Yu12 +
      Yu20*Yu22) + Yu21*(Yu01*Yu02 + Yu11*Yu12 + Yu21*Yu22) + Yu22*(Sqr(Yu02) +
      Sqr(Yu12) + Sqr(Yu22))) + 0.013333333333333334*Sqr(g1)*(-32*Sqr(g1) - 160*
      Sqr(g3) - 120*Sqr(gN)*Sqr(Qup) + 120*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) +
      0.8*QQp*Sqr(gN)*(8*Qup*Sqr(g1) + 40*Qup*Sqr(g3) + 30*Power(Qup,3)*Sqr(gN) -
      30*Qup*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) - 4*(Yu22*(Power(Yu22,3) + Yu22*
      Sqr(Yu20) + Yu22*Sqr(Yu21)) + Yu02*(Yu00*Yu20*Yu22 + Yu01*Yu21*Yu22 + Yu02*
      Sqr(Yu22)) + Yu12*(Yu10*Yu20*Yu22 + Yu11*Yu21*Yu22 + Yu12*Sqr(Yu22)));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_two_loop_mq222_dms2() const
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

   double deriv = 0.8*QQp*Sqr(gN)*(10*Power(QSp,3)*Sqr(gN) - 15*QSp*(Sqr(
      Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(
      Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) - 10*QSp*(Sqr(
      Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) - 10*QSp*
      Sqr(Lambdax)) + 8*Power(gN,4)*Sqr(QQp)*Sqr(QSp) - 2*Sqr(Lambdax)*(Sqr(Yd02)
      + Sqr(Yd12) + Sqr(Yd22)) - 2*Sqr(Lambdax)*(Sqr(Yu02) + Sqr(Yu12) + Sqr(Yu22)
      );


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_two_loop_mq222_dMassB() const
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

   double deriv = -0.8*(TYd02*Yd02 + TYd12*Yd12 + TYd22*Yd22)*Sqr(g1) - 1.6*(
      TYu02*Yu02 + TYu12*Yu12 + TYu22*Yu22)*Sqr(g1) + 0.2*MassWB*Sqr(g1)*Sqr(g2) +
      0.35555555555555557*MassG*Sqr(g1)*Sqr(g3) + 0.26666666666666666*MassBp*QQp*
      (9*Qdp + 9*QDxbarp - 9*QDxp + 9*Qep - 9*QH1p + 9*QH2p + 3*QHpbarp - 3*QHpp -
      9*QLp + 10*QQp - 18*Qup)*Sqr(g1)*Sqr(gN) + 0.0044444444444444444*MassB*Sqr(
      g1)*(867*Sqr(g1) + 5*(18*Sqr(g2) + 4*(8*Sqr(g3) + 6*QQp*(9*Qdp + 9*QDxbarp -
      9*QDxp + 9*Qep - 9*QH1p + 9*QH2p + 3*QHpbarp - 3*QHpp - 9*QLp + 10*QQp - 18
      *Qup)*Sqr(gN))) + 180*(2*(Sqr(Yd02) + Sqr(Yd12) + Sqr(Yd22)) + 4*(Sqr(Yu02)
      + Sqr(Yu12) + Sqr(Yu22)))) + 0.0044444444444444444*Sqr(g1)*(867*MassB*Sqr(g1
      ) + 5*(9*(2*MassB + MassWB)*Sqr(g2) + 4*(4*(2*MassB + MassG)*Sqr(g3) + 3*(2*
      MassB + MassBp)*QQp*(9*Qdp + 9*QDxbarp - 9*QDxp + 9*Qep - 9*QH1p + 9*QH2p +
      3*QHpbarp - 3*QHpp - 9*QLp + 10*QQp - 18*Qup)*Sqr(gN))) + 180*(-(TYd02*Yd02)
      - TYd12*Yd12 - TYd22*Yd22 - 2*(TYu02*Yu02 + TYu12*Yu12 + TYu22*Yu22) + 2*
      MassB*(Sqr(Yd02) + Sqr(Yd12) + Sqr(Yd22)) + 4*MassB*(Sqr(Yu02) + Sqr(Yu12) +
      Sqr(Yu22))));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_two_loop_mq222_dMassWB() const
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

   double deriv = 174*Power(g2,4)*MassWB + 0.4*MassB*Sqr(g1)*Sqr(g2) + 0.8*
      MassWB*Sqr(g1)*Sqr(g2) + 32*MassG*Sqr(g2)*Sqr(g3) + 64*MassWB*Sqr(g2)*Sqr(g3
      ) + 24*MassBp*Sqr(g2)*Sqr(gN)*Sqr(QQp) + 48*MassWB*Sqr(g2)*Sqr(gN)*Sqr(QQp);


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_two_loop_mq222_dMassG() const
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

   double deriv = 106.66666666666667*Power(g3,4)*MassG + 0.7111111111111111*
      MassB*Sqr(g1)*Sqr(g3) + 1.4222222222222223*MassG*Sqr(g1)*Sqr(g3) + 64*MassG*
      Sqr(g2)*Sqr(g3) + 32*MassWB*Sqr(g2)*Sqr(g3) + 42.666666666666664*MassBp*Sqr(
      g3)*Sqr(gN)*Sqr(QQp) + 85.33333333333333*MassG*Sqr(g3)*Sqr(gN)*Sqr(QQp);


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_two_loop_mq222_dMassBp() const
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

   double deriv = 0.26666666666666666*MassB*QQp*(9*Qdp + 9*QDxbarp - 9*QDxp + 9
      *Qep - 9*QH1p + 9*QH2p + 3*QHpbarp - 3*QHpp - 9*QLp + 10*QQp - 18*Qup)*Sqr(
      g1)*Sqr(gN) - 4*(TYd02*Yd02 + TYd12*Yd12 + TYd22*Yd22)*Sqr(gN)*Sqr(Qdp) - 4*
      (TYd02*Yd02 + TYd12*Yd12 + TYd22*Yd22)*Sqr(gN)*Sqr(QH1p) - 4*(TYu02*Yu02 +
      TYu12*Yu12 + TYu22*Yu22)*Sqr(gN)*Sqr(QH2p) + 4*(TYd02*Yd02 + TYd12*Yd12 +
      TYd22*Yd22)*Sqr(gN)*Sqr(QQp) + 4*(TYu02*Yu02 + TYu12*Yu12 + TYu22*Yu22)*Sqr(
      gN)*Sqr(QQp) + 12*MassWB*Sqr(g2)*Sqr(gN)*Sqr(QQp) + 21.333333333333332*MassG
      *Sqr(g3)*Sqr(gN)*Sqr(QQp) - 4*(TYu02*Yu02 + TYu12*Yu12 + TYu22*Yu22)*Sqr(gN)
      *Sqr(Qup) + 0.26666666666666666*MassBp*Sqr(gN)*(QQp*(2*(9*Qdp + 9*QDxbarp -
      9*QDxp + 9*Qep - 9*QH1p + 9*QH2p + 3*QHpbarp - 3*QHpp - 9*QLp + 10*QQp - 18*
      Qup)*Sqr(g1) + 5*QQp*(18*Sqr(g2) + 2*(16*Sqr(g3) + 9*Sqr(gN)*(9*Sqr(Qdp) + 9
      *Sqr(QDxbarp) + 9*Sqr(QDxp) + 3*Sqr(Qep) + 6*Sqr(QH1p) + 6*Sqr(QH2p) + 2*Sqr
      (QHpbarp) + 2*Sqr(QHpp) + 6*Sqr(QLp) + 20*Sqr(QQp) + 3*Sqr(QSp) + 9*Sqr(Qup)
      )))) + 15*(2*(Sqr(Qdp) + Sqr(QH1p) - Sqr(QQp))*(Sqr(Yd02) + Sqr(Yd12) + Sqr(
      Yd22)) + 2*(Sqr(QH2p) - Sqr(QQp) + Sqr(Qup))*(Sqr(Yu02) + Sqr(Yu12) + Sqr(
      Yu22)))) + 0.26666666666666666*Sqr(gN)*(QQp*((MassB + 2*MassBp)*(9*Qdp + 9*
      QDxbarp - 9*QDxp + 9*Qep - 9*QH1p + 9*QH2p + 3*QHpbarp - 3*QHpp - 9*QLp + 10
      *QQp - 18*Qup)*Sqr(g1) + 5*QQp*(9*(2*MassBp + MassWB)*Sqr(g2) + 2*(8*(2*
      MassBp + MassG)*Sqr(g3) + 9*MassBp*Sqr(gN)*(9*Sqr(Qdp) + 9*Sqr(QDxbarp) + 9*
      Sqr(QDxp) + 3*Sqr(Qep) + 6*Sqr(QH1p) + 6*Sqr(QH2p) + 2*Sqr(QHpbarp) + 2*Sqr(
      QHpp) + 6*Sqr(QLp) + 20*Sqr(QQp) + 3*Sqr(QSp) + 9*Sqr(Qup))))) + 15*(-((
      TYd02*Yd02 + TYd12*Yd12 + TYd22*Yd22)*(Sqr(Qdp) + Sqr(QH1p) - Sqr(QQp))) + 2
      *MassBp*(Sqr(Qdp) + Sqr(QH1p) - Sqr(QQp))*(Sqr(Yd02) + Sqr(Yd12) + Sqr(Yd22)
      ) + (Sqr(QH2p) - Sqr(QQp) + Sqr(Qup))*(-(TYu02*Yu02) - TYu12*Yu12 - TYu22*
      Yu22 + 2*MassBp*(Sqr(Yu02) + Sqr(Yu12) + Sqr(Yu22)))));


   return twoLoop*deriv;
}

// leading log
double lowE6SSM_tuning_calculator::deriv_dlead_log_one_loop_mq222_dLambdax() const
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

   double deriv = 4*TYd02*(2*Lambdax*TYd02 + 4*ALambdax*Lambdax*Yd02) + 2*
      Lambdax*Yd00*(mq202*Yd02 + mq220*Yd02) + 2*Lambdax*Yd01*(mq212*Yd02 + mq221*
      Yd02) + 4*TYd12*(2*Lambdax*TYd12 + 4*ALambdax*Lambdax*Yd12) + 2*Lambdax*Yd10
      *(mq202*Yd12 + mq220*Yd12) + 2*Lambdax*Yd11*(mq212*Yd12 + mq221*Yd12) + 4*
      TYd22*(2*Lambdax*TYd22 + 4*ALambdax*Lambdax*Yd22) + 2*Lambdax*Yd20*(mq202*
      Yd22 + mq220*Yd22) + 2*Lambdax*Yd21*(mq212*Yd22 + mq221*Yd22) + 2*Lambdax*
      Yd02*(mq202*Yd00 + mq220*Yd00 + mq212*Yd01 + mq221*Yd01 + 4*mHd2*Yd02 + 4*
      mq222*Yd02 + 2*(2*md200*Yd02 + md201*Yd12 + md210*Yd12 + md202*Yd22 + md220*
      Yd22)) + 2*Lambdax*Yd12*(mq202*Yd10 + mq220*Yd10 + mq212*Yd11 + mq221*Yd11 +
      4*mHd2*Yd12 + 4*mq222*Yd12 + 2*(md201*Yd02 + md210*Yd02 + 2*md211*Yd12 +
      md212*Yd22 + md221*Yd22)) + 2*Lambdax*Yd22*(mq202*Yd20 + mq220*Yd20 + mq212*
      Yd21 + mq221*Yd21 + 4*mHd2*Yd22 + 4*mq222*Yd22 + 2*(md202*Yd02 + md220*Yd02
      + md212*Yd12 + md221*Yd12 + 2*md222*Yd22)) + 4*TYu02*(2*Lambdax*TYu02 + 4*
      ALambdax*Lambdax*Yu02) + 2*Lambdax*Yu00*(mq202*Yu02 + mq220*Yu02) + 2*
      Lambdax*Yu01*(mq212*Yu02 + mq221*Yu02) + 4*TYu12*(2*Lambdax*TYu12 + 4*
      ALambdax*Lambdax*Yu12) + 2*Lambdax*Yu10*(mq202*Yu12 + mq220*Yu12) + 2*
      Lambdax*Yu11*(mq212*Yu12 + mq221*Yu12) + 4*TYu22*(2*Lambdax*TYu22 + 4*
      ALambdax*Lambdax*Yu22) + 2*Lambdax*Yu20*(mq202*Yu22 + mq220*Yu22) + 2*
      Lambdax*Yu21*(mq212*Yu22 + mq221*Yu22) + 2*Lambdax*Yu02*(mq202*Yu00 + mq220*
      Yu00 + mq212*Yu01 + mq221*Yu01 + 4*mHu2*Yu02 + 4*mq222*Yu02 + 2*(2*mu200*
      Yu02 + mu201*Yu12 + mu210*Yu12 + mu202*Yu22 + mu220*Yu22)) + 2*Lambdax*Yu12*
      (mq202*Yu10 + mq220*Yu10 + mq212*Yu11 + mq221*Yu11 + 4*mHu2*Yu12 + 4*mq222*
      Yu12 + 2*(mu201*Yu02 + mu210*Yu02 + 2*mu211*Yu12 + mu212*Yu22 + mu221*Yu22))
      + 2*Lambdax*Yu22*(mq202*Yu20 + mq220*Yu20 + mq212*Yu21 + mq221*Yu21 + 4*
      mHu2*Yu22 + 4*mq222*Yu22 + 2*(mu202*Yu02 + mu220*Yu02 + mu212*Yu12 + mu221*
      Yu12 + 2*mu222*Yu22)) + 2*QQp*QSp*(8*Lambdax*(mHd2 + mHu2 + ms2) + 8*Lambdax
      *Sqr(ALambdax))*Sqr(gN) + (4*Lambdax*mHd2 + 4*Lambdax*mHu2 + 4*Lambdax*ms2 +
      4*Lambdax*Sqr(ALambdax))*(-0.2*Sqr(g1) + 4*QH1p*QQp*Sqr(gN) + 2*(Sqr(Yd02)
      + Sqr(Yd12) + Sqr(Yd22))) + (4*Lambdax*mHd2 + 4*Lambdax*mHu2 + 4*Lambdax*ms2
      + 4*Lambdax*Sqr(ALambdax))*(0.2*Sqr(g1) + 4*QH2p*QQp*Sqr(gN) + 2*(Sqr(Yu02)
      + Sqr(Yu12) + Sqr(Yu22)));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dlead_log_one_loop_mq222_dALambdax() const
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

   double deriv = 8*TYd02*Yd02*Sqr(Lambdax) + 8*TYd12*Yd12*Sqr(Lambdax) + 8*
      TYd22*Yd22*Sqr(Lambdax) + 8*TYu02*Yu02*Sqr(Lambdax) + 8*TYu12*Yu12*Sqr(
      Lambdax) + 8*TYu22*Yu22*Sqr(Lambdax) + 16*ALambdax*QQp*QSp*Sqr(gN)*Sqr(
      Lambdax) + 4*ALambdax*Sqr(Lambdax)*(-0.2*Sqr(g1) + 4*QH1p*QQp*Sqr(gN) + 2*(
      Sqr(Yd02) + Sqr(Yd12) + Sqr(Yd22))) + 4*ALambdax*Sqr(Lambdax)*(0.2*Sqr(g1) +
      4*QH2p*QQp*Sqr(gN) + 2*(Sqr(Yu02) + Sqr(Yu12) + Sqr(Yu22)));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dlead_log_one_loop_mq222_dAYu22() const
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

   double deriv = 8*TYd02*Yu22*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22) + 8*TYd12*
      Yu22*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22) + 8*TYd22*Yu22*(Yd20*Yu20 + Yd21*
      Yu21 + Yd22*Yu22) + 4*TYu20*Yu22*(Yd00*Yd02 + Yd10*Yd12 + Yd20*Yd22 + Yu00*
      Yu02 + Yu10*Yu12 + Yu20*Yu22) + 4*TYu21*Yu22*(Yd01*Yd02 + Yd11*Yd12 + Yd21*
      Yd22 + Yu01*Yu02 + Yu11*Yu12 + Yu21*Yu22) + 16*TYu02*Yu02*Sqr(Yu22) + 16*
      TYu12*Yu12*Sqr(Yu22) + 8*AYu22*Sqr(Yu22)*(-0.4*Sqr(g1) + 6*QQp*Qup*Sqr(gN) +
      2*Sqr(Yu22)) + 4*AYu22*Sqr(Yu22)*(0.2*Sqr(g1) + 12*Sqr(gN)*Sqr(QQp) + 2*Sqr
      (Yd02) + 2*Sqr(Yd12) + 2*Sqr(Yd22) + 2*Sqr(Yu02) + 2*Sqr(Yu12) + 2*Sqr(Yu22)
      ) + 4*TYu02*(4*Yu22*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + 6*Yu02*Sqr(Yu22))
      + 4*TYu12*(4*Yu22*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + 6*Yu12*Sqr(Yu22)) +
      12*AYu22*Sqr(Yu22)*(0.2*Sqr(g1) + 4*QH2p*QQp*Sqr(gN) + 2*(Sqr(Yu02) + Sqr(
      Yu12) + Sqr(Yu22))) + 4*AYu22*Yu22*(6*Power(Yu22,3) - 0.8666666666666667*
      Yu22*Sqr(g1) - 3*Yu22*Sqr(g2) - 5.333333333333333*Yu22*Sqr(g3) + Yu22*Sqr(
      Lambdax) - 2*Yu22*Sqr(gN)*Sqr(QH2p) - 2*Yu22*Sqr(gN)*Sqr(QQp) - 2*Yu22*Sqr(
      gN)*Sqr(Qup) + Yu22*Sqr(Yd02) + Yu22*Sqr(Yd12) + Yu22*Sqr(Yd22) + 5*(Power(
      Yu22,3) + Yu22*Sqr(Yu02) + Yu22*Sqr(Yu12)) + 4*Yu22*(Sqr(Yu20) + Sqr(Yu21) +
      Sqr(Yu22)) + 3*Yu22*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(
      Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + 4*Yu22*(Yd02*(
      TYu20*Yd00 + TYu21*Yd01 + AYu22*Yd02*Yu22) + Yd12*(TYu20*Yd10 + TYu21*Yd11 +
      AYu22*Yd12*Yu22) + Yd22*(TYu20*Yd20 + TYu21*Yd21 + AYu22*Yd22*Yu22) + 2*(
      TYd02*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22) + TYd12*(Yd10*Yu20 + Yd11*Yu21 +
      Yd12*Yu22) + TYd22*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22)) - 0.8666666666666667
      *AYu22*Yu22*Sqr(g1) - 3*AYu22*Yu22*Sqr(g2) - 5.333333333333333*AYu22*Yu22*
      Sqr(g3) + AYu22*Yu22*Sqr(Lambdax) - 2*AYu22*Yu22*Sqr(gN)*Sqr(QH2p) - 2*AYu22
      *Yu22*Sqr(gN)*Sqr(QQp) - 2*AYu22*Yu22*Sqr(gN)*Sqr(Qup) + 3*AYu22*Yu22*(Sqr(
      Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20
      ) + Sqr(Yu21) + Sqr(Yu22)) + 4*(TYu02*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) +
      TYu12*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + AYu22*Yu22*(Sqr(Yu20) + Sqr(Yu21
      ) + Sqr(Yu22))) + 5*(Yu02*(TYu20*Yu00 + TYu21*Yu01 + AYu22*Yu02*Yu22) + Yu12
      *(TYu20*Yu10 + TYu21*Yu11 + AYu22*Yu12*Yu22) + Yu22*(TYu20*Yu20 + TYu21*Yu21
      + AYu22*Sqr(Yu22))) + Yu22*(2*Lambdax*TLambdax + 1.7333333333333334*MassB*
      Sqr(g1) + 6*MassWB*Sqr(g2) + 10.666666666666666*MassG*Sqr(g3) + 4*MassBp*Sqr
      (gN)*Sqr(QH2p) + 4*MassBp*Sqr(gN)*Sqr(QQp) + 4*MassBp*Sqr(gN)*Sqr(Qup) + 6*(
      TYu00*Yu00 + TYu01*Yu01 + TYu02*Yu02 + TYu10*Yu10 + TYu11*Yu11 + TYu12*Yu12
      + TYu20*Yu20 + TYu21*Yu21 + AYu22*Sqr(Yu22))));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dlead_log_one_loop_mq222_dmq222() const
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

   double deriv = 3.84*Power(g1,4) + 3*(0.2*Sqr(g1) + 6*QDxbarp*QQp*Sqr(gN))*(
      0.4*Sqr(g1) + 12*QDxbarp*QQp*Sqr(gN)) + 3*(-0.2*Sqr(g1) + 6*QDxp*QQp*Sqr(gN)
      )*(-0.4*Sqr(g1) + 12*QDxp*QQp*Sqr(gN)) + 3*(0.2*Sqr(g1) + 2*Qep*QQp*Sqr(gN))
      *(1.2*Sqr(g1) + 12*Qep*QQp*Sqr(gN)) + 2*(-0.2*Sqr(g1) + 4*QH1p*QQp*Sqr(gN))*
      (-0.6*Sqr(g1) + 12*QH1p*QQp*Sqr(gN)) + 2*(0.2*Sqr(g1) + 4*QH2p*QQp*Sqr(gN))*
      (0.6*Sqr(g1) + 12*QH2p*QQp*Sqr(gN)) + (0.2*Sqr(g1) + 4*QHpbarp*QQp*Sqr(gN))*
      (0.6*Sqr(g1) + 12*QHpbarp*QQp*Sqr(gN)) + (-0.2*Sqr(g1) + 4*QHpp*QQp*Sqr(gN))
      *(-0.6*Sqr(g1) + 12*QHpp*QQp*Sqr(gN)) + 3*(-0.2*Sqr(g1) + 4*QLp*QQp*Sqr(gN))
      *(-0.6*Sqr(g1) + 12*QLp*QQp*Sqr(gN)) + 2*Sqr(0.2*Power(g1,2) + 12*Power(gN,2
      )*Power(QQp,2)) + 72*Power(gN,4)*Sqr(QQp)*Sqr(QSp) + 24*Power(gN,4)*Sqr(QQp)
      *(9*Sqr(Qdp) + 9*Sqr(QDxbarp) + 9*Sqr(QDxp) + 3*Sqr(Qep) + 6*Sqr(QH1p) + 6*
      Sqr(QH2p) + 2*Sqr(QHpbarp) + 2*Sqr(QHpp) + 6*Sqr(QLp) + 18*Sqr(QQp) + 3*Sqr(
      QSp) + 9*Sqr(Qup)) + (0.2*Sqr(g1) + 6*Qdp*QQp*Sqr(gN) + 2*Sqr(Yd02))*(0.4*
      Sqr(g1) + 12*Qdp*QQp*Sqr(gN) + 4*Sqr(Yd02)) + 16*Sqr(Yd02)*Sqr(Yd12) + (0.2*
      Sqr(g1) + 6*Qdp*QQp*Sqr(gN) + 2*Sqr(Yd12))*(0.4*Sqr(g1) + 12*Qdp*QQp*Sqr(gN)
      + 4*Sqr(Yd12)) + 16*Sqr(Yd02)*Sqr(Yd22) + 16*Sqr(Yd12)*Sqr(Yd22) + (0.2*Sqr
      (g1) + 6*Qdp*QQp*Sqr(gN) + 2*Sqr(Yd22))*(0.4*Sqr(g1) + 12*Qdp*QQp*Sqr(gN) +
      4*Sqr(Yd22)) + (-0.2*Sqr(g1) + 4*QH1p*QQp*Sqr(gN) + 2*(Sqr(Yd02) + Sqr(Yd12)
      + Sqr(Yd22)))*(-0.6*Sqr(g1) + 12*QH1p*QQp*Sqr(gN) + 6*(Sqr(Yd02) + Sqr(Yd12
      ) + Sqr(Yd22))) + 4*Yd02*(Yu02*(Yd00*Yu00 + Yd01*Yu01 + Yd02*Yu02) + Yu12*(
      Yd00*Yu10 + Yd01*Yu11 + Yd02*Yu12) + Yu22*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22
      ) + 3*(Yd12*(Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + Yd22*(Yd00*Yd20 + Yd01*
      Yd21 + Yd02*Yd22) + Yd02*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02))) + Yd02*(
      -0.4666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(
      Lambdax) - 2*Sqr(gN)*Sqr(Qdp) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QQp) + 3
      *(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) +
      Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr
      (Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22))) + 4*
      Yd12*(Yu02*(Yd10*Yu00 + Yd11*Yu01 + Yd12*Yu02) + Yu12*(Yd10*Yu10 + Yd11*Yu11
      + Yd12*Yu12) + Yu22*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22) + 3*(Yd02*(Yd00*
      Yd10 + Yd01*Yd11 + Yd02*Yd12) + Yd22*(Yd10*Yd20 + Yd11*Yd21 + Yd12*Yd22) +
      Yd12*(Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12))) + Yd12*(-0.4666666666666667*Sqr(g1
      ) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(Qdp
      ) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QQp) + 3*(Sqr(Yd00) + Sqr(Yd01) +
      Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(
      Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(
      Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22))) + 4*Yd22*(Yu02*(Yd20*Yu00 + Yd21
      *Yu01 + Yd22*Yu02) + Yu12*(Yd20*Yu10 + Yd21*Yu11 + Yd22*Yu12) + Yu22*(Yd20*
      Yu20 + Yd21*Yu21 + Yd22*Yu22) + 3*(Yd02*(Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22)
      + Yd12*(Yd10*Yd20 + Yd11*Yd21 + Yd12*Yd22) + Yd22*(Sqr(Yd20) + Sqr(Yd21) +
      Sqr(Yd22))) + Yd22*(-0.4666666666666667*Sqr(g1) - 3*Sqr(g2) -
      5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(Qdp) - 2*Sqr(gN)*
      Sqr(QH1p) - 2*Sqr(gN)*Sqr(QQp) + 3*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(
      Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(
      Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20
      ) + Sqr(Ye21) + Sqr(Ye22))) + (-0.4*Sqr(g1) + 6*QQp*Qup*Sqr(gN) + 2*Sqr(Yu02
      ))*(-0.8*Sqr(g1) + 12*QQp*Qup*Sqr(gN) + 4*Sqr(Yu02)) + 16*Sqr(Yu02)*Sqr(Yu12
      ) + (-0.4*Sqr(g1) + 6*QQp*Qup*Sqr(gN) + 2*Sqr(Yu12))*(-0.8*Sqr(g1) + 12*QQp*
      Qup*Sqr(gN) + 4*Sqr(Yu12)) + 16*Sqr(Yu02)*Sqr(Yu22) + 16*Sqr(Yu12)*Sqr(Yu22)
      + (-0.4*Sqr(g1) + 6*QQp*Qup*Sqr(gN) + 2*Sqr(Yu22))*(-0.8*Sqr(g1) + 12*QQp*
      Qup*Sqr(gN) + 4*Sqr(Yu22)) + (0.2*Sqr(g1) + 4*QH2p*QQp*Sqr(gN) + 2*(Sqr(Yu02
      ) + Sqr(Yu12) + Sqr(Yu22)))*(0.6*Sqr(g1) + 12*QH2p*QQp*Sqr(gN) + 6*(Sqr(Yu02
      ) + Sqr(Yu12) + Sqr(Yu22))) + 4*Yu02*(Yd02*(Yd00*Yu00 + Yd01*Yu01 + Yd02*
      Yu02) + Yd12*(Yd10*Yu00 + Yd11*Yu01 + Yd12*Yu02) + Yd22*(Yd20*Yu00 + Yd21*
      Yu01 + Yd22*Yu02) + 3*(Yu12*(Yu00*Yu10 + Yu01*Yu11 + Yu02*Yu12) + Yu22*(Yu00
      *Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu02*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02))) +
      Yu02*(-0.8666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) +
      Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup
      ) + 3*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12)
      + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) + 4*Yu12*(Yd02*(Yd00*Yu10 + Yd01*
      Yu11 + Yd02*Yu12) + Yd12*(Yd10*Yu10 + Yd11*Yu11 + Yd12*Yu12) + Yd22*(Yd20*
      Yu10 + Yd21*Yu11 + Yd22*Yu12) + 3*(Yu02*(Yu00*Yu10 + Yu01*Yu11 + Yu02*Yu12)
      + Yu22*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + Yu12*(Sqr(Yu10) + Sqr(Yu11) +
      Sqr(Yu12))) + Yu12*(-0.8666666666666667*Sqr(g1) - 3*Sqr(g2) -
      5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*
      Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(
      Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) + 4*
      Yu22*(Yd02*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22) + Yd12*(Yd10*Yu20 + Yd11*Yu21
      + Yd12*Yu22) + Yd22*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22) + 3*(Yu02*(Yu00*
      Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu12*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) +
      Yu22*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + Yu22*(-0.8666666666666667*Sqr(g1
      ) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(
      QH2p) - 2*Sqr(gN)*Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3*(Sqr(Yu00) + Sqr(Yu01) +
      Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr
      (Yu22)))) + 2*Sqr(Yd00*Yd02 + Yd10*Yd12 + Yd20*Yd22 + Yu00*Yu02 + Yu10*Yu12
      + Yu20*Yu22) + 2*Sqr(Yd01*Yd02 + Yd11*Yd12 + Yd21*Yd22 + Yu01*Yu02 + Yu11*
      Yu12 + Yu21*Yu22) + Sqr(0.2*Power(g1,2) + 12*Power(gN,2)*Power(QQp,2) + 2*
      Power(Yd02,2) + 2*Power(Yd12,2) + 2*Power(Yd22,2) + 2*Power(Yu02,2) + 2*
      Power(Yu12,2) + 2*Power(Yu22,2));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dlead_log_one_loop_mq222_dmHd2() const
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

   double deriv = -3.84*Power(g1,4) + 16*Yd02*Yd12*(Yd00*Yd10 + Yd01*Yd11 +
      Yd02*Yd12) + 16*Yd02*Yd22*(Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22) + 16*Yd12*Yd22
      *(Yd10*Yd20 + Yd11*Yd21 + Yd12*Yd22) + 4*(Yd00*Yd02 + Yd10*Yd12 + Yd20*Yd22)
      *(Yd00*Yd02 + Yd10*Yd12 + Yd20*Yd22 + Yu00*Yu02 + Yu10*Yu12 + Yu20*Yu22) + 4
      *(Yd01*Yd02 + Yd11*Yd12 + Yd21*Yd22)*(Yd01*Yd02 + Yd11*Yd12 + Yd21*Yd22 +
      Yu01*Yu02 + Yu11*Yu12 + Yu21*Yu22) + 3*(-0.4*Sqr(g1) + 4*QDxbarp*QH1p*Sqr(gN
      ))*(0.2*Sqr(g1) + 6*QDxbarp*QQp*Sqr(gN)) + 3*(0.4*Sqr(g1) + 4*QDxp*QH1p*Sqr(
      gN))*(-0.2*Sqr(g1) + 6*QDxp*QQp*Sqr(gN)) + 2*(-0.6*Sqr(g1) + 4*QH1p*QH2p*Sqr
      (gN))*(0.2*Sqr(g1) + 4*QH2p*QQp*Sqr(gN)) + (-0.6*Sqr(g1) + 4*QH1p*QHpbarp*
      Sqr(gN))*(0.2*Sqr(g1) + 4*QHpbarp*QQp*Sqr(gN)) + (0.6*Sqr(g1) + 4*QH1p*QHpp*
      Sqr(gN))*(-0.2*Sqr(g1) + 4*QHpp*QQp*Sqr(gN)) + 2*QQp*QSp*Sqr(gN)*(4*QH1p*QSp
      *Sqr(gN) + 4*Sqr(Lambdax)) + 2*(-0.2*Sqr(g1) + 4*QH1p*QQp*Sqr(gN))*(0.6*Sqr(
      g1) + 4*Sqr(gN)*Sqr(QH1p)) + 16*Power(gN,4)*QH1p*QQp*Sqr(QSp) + 8*Power(gN,4
      )*QH1p*QQp*(9*Sqr(Qdp) + 9*Sqr(QDxbarp) + 9*Sqr(QDxp) + 3*Sqr(Qep) + 6*Sqr(
      QH1p) + 6*Sqr(QH2p) + 2*Sqr(QHpbarp) + 2*Sqr(QHpp) + 6*Sqr(QLp) + 18*Sqr(QQp
      ) + 3*Sqr(QSp) + 9*Sqr(Qup)) + (0.2*Sqr(g1) + 6*Qdp*QQp*Sqr(gN) + 2*Sqr(Yd02
      ))*(-0.4*Sqr(g1) + 4*Qdp*QH1p*Sqr(gN) + 4*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02)
      )) + (0.2*Sqr(g1) + 6*Qdp*QQp*Sqr(gN) + 2*Sqr(Yd12))*(-0.4*Sqr(g1) + 4*Qdp*
      QH1p*Sqr(gN) + 4*(Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12))) + (0.2*Sqr(g1) + 12*
      Sqr(gN)*Sqr(QQp))*(-0.2*Sqr(g1) + 4*QH1p*QQp*Sqr(gN) + 2*(Sqr(Yd00) + Sqr(
      Yd10) + Sqr(Yd20))) + (0.2*Sqr(g1) + 12*Sqr(gN)*Sqr(QQp))*(-0.2*Sqr(g1) + 4*
      QH1p*QQp*Sqr(gN) + 2*(Sqr(Yd01) + Sqr(Yd11) + Sqr(Yd21))) + (0.2*Sqr(g1) + 6
      *Qdp*QQp*Sqr(gN) + 2*Sqr(Yd22))*(-0.4*Sqr(g1) + 4*Qdp*QH1p*Sqr(gN) + 4*(Sqr(
      Yd20) + Sqr(Yd21) + Sqr(Yd22))) + (0.2*Sqr(g1) + 2*Qep*QQp*Sqr(gN))*(-1.2*
      Sqr(g1) + 4*Qep*QH1p*Sqr(gN) + 4*(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02))) + (0.2
      *Sqr(g1) + 2*Qep*QQp*Sqr(gN))*(-1.2*Sqr(g1) + 4*Qep*QH1p*Sqr(gN) + 4*(Sqr(
      Ye10) + Sqr(Ye11) + Sqr(Ye12))) + (-0.2*Sqr(g1) + 4*QLp*QQp*Sqr(gN))*(0.6*
      Sqr(g1) + 4*QH1p*QLp*Sqr(gN) + 2*(Sqr(Ye00) + Sqr(Ye10) + Sqr(Ye20))) + (
      -0.2*Sqr(g1) + 4*QLp*QQp*Sqr(gN))*(0.6*Sqr(g1) + 4*QH1p*QLp*Sqr(gN) + 2*(Sqr
      (Ye01) + Sqr(Ye11) + Sqr(Ye21))) + (-0.2*Sqr(g1) + 4*QLp*QQp*Sqr(gN))*(0.6*
      Sqr(g1) + 4*QH1p*QLp*Sqr(gN) + 2*(Sqr(Ye02) + Sqr(Ye12) + Sqr(Ye22))) + (0.2
      *Sqr(g1) + 2*Qep*QQp*Sqr(gN))*(-1.2*Sqr(g1) + 4*Qep*QH1p*Sqr(gN) + 4*(Sqr(
      Ye20) + Sqr(Ye21) + Sqr(Ye22))) + (-0.2*Sqr(g1) + 4*QH1p*QQp*Sqr(gN) + 2*(
      Sqr(Yd02) + Sqr(Yd12) + Sqr(Yd22)))*(0.6*Sqr(g1) + 2*Sqr(Lambdax) + 4*Sqr(gN
      )*Sqr(QH1p) + 6*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) +
      Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + 2*(Sqr(Ye00) + Sqr(Ye01) +
      Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr
      (Ye22))) + 4*Yd02*(Yu02*(Yd00*Yu00 + Yd01*Yu01 + Yd02*Yu02) + Yu12*(Yd00*
      Yu10 + Yd01*Yu11 + Yd02*Yu12) + Yu22*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22) + 3
      *(Yd12*(Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + Yd22*(Yd00*Yd20 + Yd01*Yd21 +
      Yd02*Yd22) + Yd02*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02))) + Yd02*(
      -0.4666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(
      Lambdax) - 2*Sqr(gN)*Sqr(Qdp) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QQp) + 3
      *(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) +
      Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr
      (Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22))) + 4*
      Yd12*(Yu02*(Yd10*Yu00 + Yd11*Yu01 + Yd12*Yu02) + Yu12*(Yd10*Yu10 + Yd11*Yu11
      + Yd12*Yu12) + Yu22*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22) + 3*(Yd02*(Yd00*
      Yd10 + Yd01*Yd11 + Yd02*Yd12) + Yd22*(Yd10*Yd20 + Yd11*Yd21 + Yd12*Yd22) +
      Yd12*(Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12))) + Yd12*(-0.4666666666666667*Sqr(g1
      ) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(Qdp
      ) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QQp) + 3*(Sqr(Yd00) + Sqr(Yd01) +
      Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(
      Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(
      Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22))) + 4*Yd22*(Yu02*(Yd20*Yu00 + Yd21
      *Yu01 + Yd22*Yu02) + Yu12*(Yd20*Yu10 + Yd21*Yu11 + Yd22*Yu12) + Yu22*(Yd20*
      Yu20 + Yd21*Yu21 + Yd22*Yu22) + 3*(Yd02*(Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22)
      + Yd12*(Yd10*Yd20 + Yd11*Yd21 + Yd12*Yd22) + Yd22*(Sqr(Yd20) + Sqr(Yd21) +
      Sqr(Yd22))) + Yd22*(-0.4666666666666667*Sqr(g1) - 3*Sqr(g2) -
      5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(Qdp) - 2*Sqr(gN)*
      Sqr(QH1p) - 2*Sqr(gN)*Sqr(QQp) + 3*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(
      Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(
      Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20
      ) + Sqr(Ye21) + Sqr(Ye22))) + (0.8*Sqr(g1) + 4*QH1p*Qup*Sqr(gN))*(-0.4*Sqr(
      g1) + 6*QQp*Qup*Sqr(gN) + 2*Sqr(Yu02)) + (0.8*Sqr(g1) + 4*QH1p*Qup*Sqr(gN))*
      (-0.4*Sqr(g1) + 6*QQp*Qup*Sqr(gN) + 2*Sqr(Yu12)) + (0.8*Sqr(g1) + 4*QH1p*Qup
      *Sqr(gN))*(-0.4*Sqr(g1) + 6*QQp*Qup*Sqr(gN) + 2*Sqr(Yu22)) + (-0.2*Sqr(g1) +
      4*QH1p*QQp*Sqr(gN) + 2*(Sqr(Yd02) + Sqr(Yd12) + Sqr(Yd22)))*(0.2*Sqr(g1) +
      12*Sqr(gN)*Sqr(QQp) + 2*Sqr(Yd02) + 2*Sqr(Yd12) + 2*Sqr(Yd22) + 2*Sqr(Yu02)
      + 2*Sqr(Yu12) + 2*Sqr(Yu22)) + (-0.6*Sqr(g1) + 4*QH1p*QH2p*Sqr(gN) + 2*Sqr(
      Lambdax))*(0.2*Sqr(g1) + 4*QH2p*QQp*Sqr(gN) + 2*(Sqr(Yu02) + Sqr(Yu12) + Sqr
      (Yu22)));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dlead_log_one_loop_mq222_dmHu2() const
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

   double deriv = 3.84*Power(g1,4) + 16*Yu02*Yu12*(Yu00*Yu10 + Yu01*Yu11 + Yu02
      *Yu12) + 16*Yu02*Yu22*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + 16*Yu12*Yu22*(
      Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + 4*(Yu00*Yu02 + Yu10*Yu12 + Yu20*Yu22)*(
      Yd00*Yd02 + Yd10*Yd12 + Yd20*Yd22 + Yu00*Yu02 + Yu10*Yu12 + Yu20*Yu22) + 4*(
      Yu01*Yu02 + Yu11*Yu12 + Yu21*Yu22)*(Yd01*Yd02 + Yd11*Yd12 + Yd21*Yd22 + Yu01
      *Yu02 + Yu11*Yu12 + Yu21*Yu22) + 3*(0.4*Sqr(g1) + 4*QDxbarp*QH2p*Sqr(gN))*(
      0.2*Sqr(g1) + 6*QDxbarp*QQp*Sqr(gN)) + 3*(-0.4*Sqr(g1) + 4*QDxp*QH2p*Sqr(gN)
      )*(-0.2*Sqr(g1) + 6*QDxp*QQp*Sqr(gN)) + 3*(1.2*Sqr(g1) + 4*Qep*QH2p*Sqr(gN))
      *(0.2*Sqr(g1) + 2*Qep*QQp*Sqr(gN)) + 2*(-0.6*Sqr(g1) + 4*QH1p*QH2p*Sqr(gN))*
      (-0.2*Sqr(g1) + 4*QH1p*QQp*Sqr(gN)) + (0.6*Sqr(g1) + 4*QH2p*QHpbarp*Sqr(gN))
      *(0.2*Sqr(g1) + 4*QHpbarp*QQp*Sqr(gN)) + (-0.6*Sqr(g1) + 4*QH2p*QHpp*Sqr(gN)
      )*(-0.2*Sqr(g1) + 4*QHpp*QQp*Sqr(gN)) + 3*(-0.6*Sqr(g1) + 4*QH2p*QLp*Sqr(gN)
      )*(-0.2*Sqr(g1) + 4*QLp*QQp*Sqr(gN)) + 2*QQp*QSp*Sqr(gN)*(4*QH2p*QSp*Sqr(gN)
      + 4*Sqr(Lambdax)) + 2*(0.2*Sqr(g1) + 4*QH2p*QQp*Sqr(gN))*(0.6*Sqr(g1) + 4*
      Sqr(gN)*Sqr(QH2p)) + 16*Power(gN,4)*QH2p*QQp*Sqr(QSp) + 8*Power(gN,4)*QH2p*
      QQp*(9*Sqr(Qdp) + 9*Sqr(QDxbarp) + 9*Sqr(QDxp) + 3*Sqr(Qep) + 6*Sqr(QH1p) +
      6*Sqr(QH2p) + 2*Sqr(QHpbarp) + 2*Sqr(QHpp) + 6*Sqr(QLp) + 18*Sqr(QQp) + 3*
      Sqr(QSp) + 9*Sqr(Qup)) + (0.4*Sqr(g1) + 4*Qdp*QH2p*Sqr(gN))*(0.2*Sqr(g1) + 6
      *Qdp*QQp*Sqr(gN) + 2*Sqr(Yd02)) + (0.4*Sqr(g1) + 4*Qdp*QH2p*Sqr(gN))*(0.2*
      Sqr(g1) + 6*Qdp*QQp*Sqr(gN) + 2*Sqr(Yd12)) + (0.4*Sqr(g1) + 4*Qdp*QH2p*Sqr(
      gN))*(0.2*Sqr(g1) + 6*Qdp*QQp*Sqr(gN) + 2*Sqr(Yd22)) + (-0.6*Sqr(g1) + 4*
      QH1p*QH2p*Sqr(gN) + 2*Sqr(Lambdax))*(-0.2*Sqr(g1) + 4*QH1p*QQp*Sqr(gN) + 2*(
      Sqr(Yd02) + Sqr(Yd12) + Sqr(Yd22))) + (-0.4*Sqr(g1) + 6*QQp*Qup*Sqr(gN) + 2*
      Sqr(Yu02))*(-0.8*Sqr(g1) + 4*QH2p*Qup*Sqr(gN) + 4*(Sqr(Yu00) + Sqr(Yu01) +
      Sqr(Yu02))) + (-0.4*Sqr(g1) + 6*QQp*Qup*Sqr(gN) + 2*Sqr(Yu12))*(-0.8*Sqr(g1)
      + 4*QH2p*Qup*Sqr(gN) + 4*(Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12))) + (0.2*Sqr(g1
      ) + 12*Sqr(gN)*Sqr(QQp))*(0.2*Sqr(g1) + 4*QH2p*QQp*Sqr(gN) + 2*(Sqr(Yu00) +
      Sqr(Yu10) + Sqr(Yu20))) + (0.2*Sqr(g1) + 12*Sqr(gN)*Sqr(QQp))*(0.2*Sqr(g1) +
      4*QH2p*QQp*Sqr(gN) + 2*(Sqr(Yu01) + Sqr(Yu11) + Sqr(Yu21))) + (0.2*Sqr(g1)
      + 12*Sqr(gN)*Sqr(QQp) + 2*Sqr(Yd02) + 2*Sqr(Yd12) + 2*Sqr(Yd22) + 2*Sqr(Yu02
      ) + 2*Sqr(Yu12) + 2*Sqr(Yu22))*(0.2*Sqr(g1) + 4*QH2p*QQp*Sqr(gN) + 2*(Sqr(
      Yu02) + Sqr(Yu12) + Sqr(Yu22))) + (-0.4*Sqr(g1) + 6*QQp*Qup*Sqr(gN) + 2*Sqr(
      Yu22))*(-0.8*Sqr(g1) + 4*QH2p*Qup*Sqr(gN) + 4*(Sqr(Yu20) + Sqr(Yu21) + Sqr(
      Yu22))) + (0.2*Sqr(g1) + 4*QH2p*QQp*Sqr(gN) + 2*(Sqr(Yu02) + Sqr(Yu12) + Sqr
      (Yu22)))*(0.6*Sqr(g1) + 2*Sqr(Lambdax) + 4*Sqr(gN)*Sqr(QH2p) + 6*(Sqr(Yu00)
      + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) +
      Sqr(Yu21) + Sqr(Yu22))) + 4*Yu02*(Yd02*(Yd00*Yu00 + Yd01*Yu01 + Yd02*Yu02) +
      Yd12*(Yd10*Yu00 + Yd11*Yu01 + Yd12*Yu02) + Yd22*(Yd20*Yu00 + Yd21*Yu01 +
      Yd22*Yu02) + 3*(Yu12*(Yu00*Yu10 + Yu01*Yu11 + Yu02*Yu12) + Yu22*(Yu00*Yu20 +
      Yu01*Yu21 + Yu02*Yu22) + Yu02*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02))) + Yu02*(
      -0.8666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(
      Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3
      *(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) +
      Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) + 4*Yu12*(Yd02*(Yd00*Yu10 + Yd01*Yu11 +
      Yd02*Yu12) + Yd12*(Yd10*Yu10 + Yd11*Yu11 + Yd12*Yu12) + Yd22*(Yd20*Yu10 +
      Yd21*Yu11 + Yd22*Yu12) + 3*(Yu02*(Yu00*Yu10 + Yu01*Yu11 + Yu02*Yu12) + Yu22*
      (Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + Yu12*(Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12
      ))) + Yu12*(-0.8666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(
      g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QQp) - 2*Sqr(gN)*
      Sqr(Qup) + 3*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) +
      Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) + 4*Yu22*(Yd02*(Yd00*Yu20 +
      Yd01*Yu21 + Yd02*Yu22) + Yd12*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22) + Yd22*(
      Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22) + 3*(Yu02*(Yu00*Yu20 + Yu01*Yu21 + Yu02*
      Yu22) + Yu12*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + Yu22*(Sqr(Yu20) + Sqr(
      Yu21) + Sqr(Yu22))) + Yu22*(-0.8666666666666667*Sqr(g1) - 3*Sqr(g2) -
      5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*
      Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(
      Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dlead_log_one_loop_mq222_dmu222() const
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

   double deriv = -7.68*Power(g1,4) + 8*Yu02*Yu22*(Yu00*Yu20 + Yu01*Yu21 + Yu02
      *Yu22) + 8*Yu12*Yu22*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + 4*Yu20*Yu22*(Yd00
      *Yd02 + Yd10*Yd12 + Yd20*Yd22 + Yu00*Yu02 + Yu10*Yu12 + Yu20*Yu22) + 4*Yu21*
      Yu22*(Yd01*Yd02 + Yd11*Yd12 + Yd21*Yd22 + Yu01*Yu02 + Yu11*Yu12 + Yu21*Yu22)
      + 3*(0.2*Sqr(g1) + 6*QDxbarp*QQp*Sqr(gN))*(-0.8*Sqr(g1) + 6*QDxbarp*Qup*Sqr
      (gN)) + 3*(-0.2*Sqr(g1) + 6*QDxp*QQp*Sqr(gN))*(0.8*Sqr(g1) + 6*QDxp*Qup*Sqr(
      gN)) + 3*(0.2*Sqr(g1) + 2*Qep*QQp*Sqr(gN))*(-2.4*Sqr(g1) + 6*Qep*Qup*Sqr(gN)
      ) + 2*(-0.2*Sqr(g1) + 4*QH1p*QQp*Sqr(gN))*(1.2*Sqr(g1) + 6*QH1p*Qup*Sqr(gN))
      + 2*(0.2*Sqr(g1) + 4*QH2p*QQp*Sqr(gN))*(-1.2*Sqr(g1) + 6*QH2p*Qup*Sqr(gN))
      + (0.2*Sqr(g1) + 4*QHpbarp*QQp*Sqr(gN))*(-1.2*Sqr(g1) + 6*QHpbarp*Qup*Sqr(gN
      )) + (-0.2*Sqr(g1) + 4*QHpp*QQp*Sqr(gN))*(1.2*Sqr(g1) + 6*QHpp*Qup*Sqr(gN))
      + 3*(-0.2*Sqr(g1) + 4*QLp*QQp*Sqr(gN))*(1.2*Sqr(g1) + 6*QLp*Qup*Sqr(gN)) +
      36*Power(gN,4)*QQp*Qup*Sqr(QSp) + 12*Power(gN,4)*QQp*Qup*(9*Sqr(Qdp) + 9*Sqr
      (QDxbarp) + 9*Sqr(QDxp) + 3*Sqr(Qep) + 6*Sqr(QH1p) + 6*Sqr(QH2p) + 2*Sqr(
      QHpbarp) + 2*Sqr(QHpp) + 6*Sqr(QLp) + 18*Sqr(QQp) + 3*Sqr(QSp) + 9*Sqr(Qup))
      + (-0.8*Sqr(g1) + 6*Qdp*Qup*Sqr(gN))*(0.2*Sqr(g1) + 6*Qdp*QQp*Sqr(gN) + 2*
      Sqr(Yd02)) + (-0.8*Sqr(g1) + 6*Qdp*Qup*Sqr(gN))*(0.2*Sqr(g1) + 6*Qdp*QQp*Sqr
      (gN) + 2*Sqr(Yd12)) + (-0.8*Sqr(g1) + 6*Qdp*Qup*Sqr(gN))*(0.2*Sqr(g1) + 6*
      Qdp*QQp*Sqr(gN) + 2*Sqr(Yd22)) + (1.2*Sqr(g1) + 6*QH1p*Qup*Sqr(gN))*(-0.2*
      Sqr(g1) + 4*QH1p*QQp*Sqr(gN) + 2*(Sqr(Yd02) + Sqr(Yd12) + Sqr(Yd22))) + (1.6
      *Sqr(g1) + 6*Sqr(gN)*Sqr(Qup))*(-0.4*Sqr(g1) + 6*QQp*Qup*Sqr(gN) + 2*Sqr(
      Yu02)) + (1.6*Sqr(g1) + 6*Sqr(gN)*Sqr(Qup))*(-0.4*Sqr(g1) + 6*QQp*Qup*Sqr(gN
      ) + 2*Sqr(Yu12)) + (0.2*Sqr(g1) + 12*Sqr(gN)*Sqr(QQp))*(-0.4*Sqr(g1) + 6*QQp
      *Qup*Sqr(gN) + 2*Sqr(Yu20)) + (0.2*Sqr(g1) + 12*Sqr(gN)*Sqr(QQp))*(-0.4*Sqr(
      g1) + 6*QQp*Qup*Sqr(gN) + 2*Sqr(Yu21)) + (-0.4*Sqr(g1) + 6*QQp*Qup*Sqr(gN) +
      2*Sqr(Yu22))*(0.2*Sqr(g1) + 12*Sqr(gN)*Sqr(QQp) + 2*Sqr(Yd02) + 2*Sqr(Yd12)
      + 2*Sqr(Yd22) + 2*Sqr(Yu02) + 2*Sqr(Yu12) + 2*Sqr(Yu22)) + (-0.4*Sqr(g1) +
      6*QQp*Qup*Sqr(gN) + 2*Sqr(Yu22))*(1.6*Sqr(g1) + 6*Sqr(gN)*Sqr(Qup) + 4*(Sqr(
      Yu20) + Sqr(Yu21) + Sqr(Yu22))) + (0.2*Sqr(g1) + 4*QH2p*QQp*Sqr(gN) + 2*(Sqr
      (Yu02) + Sqr(Yu12) + Sqr(Yu22)))*(-1.2*Sqr(g1) + 6*QH2p*Qup*Sqr(gN) + 6*(Sqr
      (Yu20) + Sqr(Yu21) + Sqr(Yu22))) + 4*Yu22*(Yd02*(Yd00*Yu20 + Yd01*Yu21 +
      Yd02*Yu22) + Yd12*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22) + Yd22*(Yd20*Yu20 +
      Yd21*Yu21 + Yd22*Yu22) + 3*(Yu02*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu12*
      (Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + Yu22*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22
      ))) + Yu22*(-0.8666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(
      g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QQp) - 2*Sqr(gN)*
      Sqr(Qup) + 3*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) +
      Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dlead_log_one_loop_mq222_dms2() const
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

   double deriv = 8*Power(gN,4)*QQp*Power(QSp,3) + 6*Qep*QSp*Sqr(gN)*(0.2*Sqr(
      g1) + 2*Qep*QQp*Sqr(gN)) + 2*QHpbarp*QSp*Sqr(gN)*(0.2*Sqr(g1) + 4*QHpbarp*
      QQp*Sqr(gN)) + 2*QHpp*QSp*Sqr(gN)*(-0.2*Sqr(g1) + 4*QHpp*QQp*Sqr(gN)) + 6*
      QLp*QSp*Sqr(gN)*(-0.2*Sqr(g1) + 4*QLp*QQp*Sqr(gN)) + (-0.2*Sqr(g1) + 6*QDxp*
      QQp*Sqr(gN))*(2*QDxp*QSp*Sqr(gN) + 2*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(
      Kappa02))) + (-0.2*Sqr(g1) + 6*QDxp*QQp*Sqr(gN))*(2*QDxp*QSp*Sqr(gN) + 2*(
      Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12))) + (0.2*Sqr(g1) + 6*QDxbarp*QQp*
      Sqr(gN))*(2*QDxbarp*QSp*Sqr(gN) + 2*(Sqr(Kappa00) + Sqr(Kappa10) + Sqr(
      Kappa20))) + (0.2*Sqr(g1) + 6*QDxbarp*QQp*Sqr(gN))*(2*QDxbarp*QSp*Sqr(gN) +
      2*(Sqr(Kappa01) + Sqr(Kappa11) + Sqr(Kappa21))) + (0.2*Sqr(g1) + 6*QDxbarp*
      QQp*Sqr(gN))*(2*QDxbarp*QSp*Sqr(gN) + 2*(Sqr(Kappa02) + Sqr(Kappa12) + Sqr(
      Kappa22))) + (-0.2*Sqr(g1) + 6*QDxp*QQp*Sqr(gN))*(2*QDxp*QSp*Sqr(gN) + 2*(
      Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22))) + (0.2*Sqr(g1) + 4*QH2p*QQp*Sqr
      (gN))*(2*QH2p*QSp*Sqr(gN) + 2*(Sqr(Lambda1200) + Sqr(Lambda1201))) + (-0.2*
      Sqr(g1) + 4*QH1p*QQp*Sqr(gN))*(2*QH1p*QSp*Sqr(gN) + 2*(Sqr(Lambda1200) + Sqr
      (Lambda1210))) + (-0.2*Sqr(g1) + 4*QH1p*QQp*Sqr(gN))*(2*QH1p*QSp*Sqr(gN) + 2
      *(Sqr(Lambda1201) + Sqr(Lambda1211))) + (0.2*Sqr(g1) + 4*QH2p*QQp*Sqr(gN))*(
      2*QH2p*QSp*Sqr(gN) + 2*(Sqr(Lambda1210) + Sqr(Lambda1211))) + 4*QQp*QSp*Sqr(
      gN)*(0.2*Sqr(g1) + 12*Sqr(gN)*Sqr(QQp)) + 2*QQp*QSp*Sqr(gN)*(6*(Sqr(Kappa00)
      + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12)
      + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 4*(Sqr(Lambda1200) + Sqr(
      Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 4*Sqr(Lambdax) + 2*Sqr(gN
      )*Sqr(QSp)) + 4*Power(gN,4)*QQp*QSp*(9*Sqr(Qdp) + 9*Sqr(QDxbarp) + 9*Sqr(
      QDxp) + 3*Sqr(Qep) + 6*Sqr(QH1p) + 6*Sqr(QH2p) + 2*Sqr(QHpbarp) + 2*Sqr(QHpp
      ) + 6*Sqr(QLp) + 18*Sqr(QQp) + 3*Sqr(QSp) + 9*Sqr(Qup)) + 2*Qdp*QSp*Sqr(gN)*
      (0.2*Sqr(g1) + 6*Qdp*QQp*Sqr(gN) + 2*Sqr(Yd02)) + 2*Qdp*QSp*Sqr(gN)*(0.2*Sqr
      (g1) + 6*Qdp*QQp*Sqr(gN) + 2*Sqr(Yd12)) + 2*Qdp*QSp*Sqr(gN)*(0.2*Sqr(g1) + 6
      *Qdp*QQp*Sqr(gN) + 2*Sqr(Yd22)) + (2*QH1p*QSp*Sqr(gN) + 2*Sqr(Lambdax))*(
      -0.2*Sqr(g1) + 4*QH1p*QQp*Sqr(gN) + 2*(Sqr(Yd02) + Sqr(Yd12) + Sqr(Yd22))) +
      2*QSp*Qup*Sqr(gN)*(-0.4*Sqr(g1) + 6*QQp*Qup*Sqr(gN) + 2*Sqr(Yu02)) + 2*QSp*
      Qup*Sqr(gN)*(-0.4*Sqr(g1) + 6*QQp*Qup*Sqr(gN) + 2*Sqr(Yu12)) + 2*QSp*Qup*Sqr
      (gN)*(-0.4*Sqr(g1) + 6*QQp*Qup*Sqr(gN) + 2*Sqr(Yu22)) + 2*QQp*QSp*Sqr(gN)*(
      0.2*Sqr(g1) + 12*Sqr(gN)*Sqr(QQp) + 2*Sqr(Yd02) + 2*Sqr(Yd12) + 2*Sqr(Yd22)
      + 2*Sqr(Yu02) + 2*Sqr(Yu12) + 2*Sqr(Yu22)) + (2*QH2p*QSp*Sqr(gN) + 2*Sqr(
      Lambdax))*(0.2*Sqr(g1) + 4*QH2p*QQp*Sqr(gN) + 2*(Sqr(Yu02) + Sqr(Yu12) + Sqr
      (Yu22)));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dlead_log_one_loop_mq222_dMassB() const
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

   double deriv = -15.36*Power(g1,4)*MassB + 3.7333333333333334*TYd02*Yd02*Sqr(
      g1) + 3.7333333333333334*TYd12*Yd12*Sqr(g1) + 3.7333333333333334*TYd22*Yd22*
      Sqr(g1) + 6.933333333333334*TYu02*Yu02*Sqr(g1) + 6.933333333333334*TYu12*
      Yu12*Sqr(g1) + 6.933333333333334*TYu22*Yu22*Sqr(g1) - 3.2*MassB*Sqr(g1)*(0.2
      *Sqr(g1) + 6*QDxbarp*QQp*Sqr(gN)) - 3.2*MassB*Sqr(g1)*(-0.2*Sqr(g1) + 6*QDxp
      *QQp*Sqr(gN)) - 28.8*MassB*Sqr(g1)*(0.2*Sqr(g1) + 2*Qep*QQp*Sqr(gN)) - 4.8*
      MassB*Sqr(g1)*(-0.2*Sqr(g1) + 4*QH1p*QQp*Sqr(gN)) - 4.8*MassB*Sqr(g1)*(0.2*
      Sqr(g1) + 4*QH2p*QQp*Sqr(gN)) - 2.4*MassB*Sqr(g1)*(0.2*Sqr(g1) + 4*QHpbarp*
      QQp*Sqr(gN)) - 2.4*MassB*Sqr(g1)*(-0.2*Sqr(g1) + 4*QHpp*QQp*Sqr(gN)) - 7.2*
      MassB*Sqr(g1)*(-0.2*Sqr(g1) + 4*QLp*QQp*Sqr(gN)) - 0.5333333333333333*MassB*
      Sqr(g1)*(0.2*Sqr(g1) + 12*Sqr(gN)*Sqr(QQp)) - 1.0666666666666667*MassB*Sqr(
      g1)*(0.2*Sqr(g1) + 6*Qdp*QQp*Sqr(gN) + 2*Sqr(Yd02)) - 1.0666666666666667*
      MassB*Sqr(g1)*(0.2*Sqr(g1) + 6*Qdp*QQp*Sqr(gN) + 2*Sqr(Yd12)) -
      1.0666666666666667*MassB*Sqr(g1)*(0.2*Sqr(g1) + 6*Qdp*QQp*Sqr(gN) + 2*Sqr(
      Yd22)) - 2.4*MassB*Sqr(g1)*(-0.2*Sqr(g1) + 4*QH1p*QQp*Sqr(gN) + 2*(Sqr(Yd02)
      + Sqr(Yd12) + Sqr(Yd22))) - 4.266666666666667*MassB*Sqr(g1)*(-0.4*Sqr(g1) +
      6*QQp*Qup*Sqr(gN) + 2*Sqr(Yu02)) - 4.266666666666667*MassB*Sqr(g1)*(-0.4*
      Sqr(g1) + 6*QQp*Qup*Sqr(gN) + 2*Sqr(Yu12)) - 4.266666666666667*MassB*Sqr(g1)
      *(-0.4*Sqr(g1) + 6*QQp*Qup*Sqr(gN) + 2*Sqr(Yu22)) - 0.26666666666666666*
      MassB*Sqr(g1)*(0.2*Sqr(g1) + 12*Sqr(gN)*Sqr(QQp) + 2*Sqr(Yd02) + 2*Sqr(Yd12)
      + 2*Sqr(Yd22) + 2*Sqr(Yu02) + 2*Sqr(Yu12) + 2*Sqr(Yu22)) - 2.4*MassB*Sqr(g1
      )*(0.2*Sqr(g1) + 4*QH2p*QQp*Sqr(gN) + 2*(Sqr(Yu02) + Sqr(Yu12) + Sqr(Yu22)))
      ;


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dlead_log_one_loop_mq222_dMassWB() const
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

   double deriv = -288*Power(g2,4)*MassWB + 24*TYd02*Yd02*Sqr(g2) + 24*TYd12*
      Yd12*Sqr(g2) + 24*TYd22*Yd22*Sqr(g2) + 24*TYu02*Yu02*Sqr(g2) + 24*TYu12*Yu12
      *Sqr(g2) + 24*TYu22*Yu22*Sqr(g2) - 24*MassWB*Sqr(g2)*(-0.2*Sqr(g1) + 4*QH1p*
      QQp*Sqr(gN)) - 24*MassWB*Sqr(g2)*(0.2*Sqr(g1) + 4*QH2p*QQp*Sqr(gN)) - 12*
      MassWB*Sqr(g2)*(0.2*Sqr(g1) + 4*QHpbarp*QQp*Sqr(gN)) - 12*MassWB*Sqr(g2)*(
      -0.2*Sqr(g1) + 4*QHpp*QQp*Sqr(gN)) - 36*MassWB*Sqr(g2)*(-0.2*Sqr(g1) + 4*QLp
      *QQp*Sqr(gN)) - 24*MassWB*Sqr(g2)*(0.2*Sqr(g1) + 12*Sqr(gN)*Sqr(QQp)) - 12*
      MassWB*Sqr(g2)*(-0.2*Sqr(g1) + 4*QH1p*QQp*Sqr(gN) + 2*(Sqr(Yd02) + Sqr(Yd12)
      + Sqr(Yd22))) - 12*MassWB*Sqr(g2)*(0.2*Sqr(g1) + 12*Sqr(gN)*Sqr(QQp) + 2*
      Sqr(Yd02) + 2*Sqr(Yd12) + 2*Sqr(Yd22) + 2*Sqr(Yu02) + 2*Sqr(Yu12) + 2*Sqr(
      Yu22)) - 12*MassWB*Sqr(g2)*(0.2*Sqr(g1) + 4*QH2p*QQp*Sqr(gN) + 2*(Sqr(Yu02)
      + Sqr(Yu12) + Sqr(Yu22)));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dlead_log_one_loop_mq222_dMassG() const
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

   double deriv = 42.666666666666664*TYd02*Yd02*Sqr(g3) + 42.666666666666664*
      TYd12*Yd12*Sqr(g3) + 42.666666666666664*TYd22*Yd22*Sqr(g3) +
      42.666666666666664*TYu02*Yu02*Sqr(g3) + 42.666666666666664*TYu12*Yu12*Sqr(g3
      ) + 42.666666666666664*TYu22*Yu22*Sqr(g3) - 64*MassG*Sqr(g3)*(0.2*Sqr(g1) +
      6*QDxbarp*QQp*Sqr(gN)) - 64*MassG*Sqr(g3)*(-0.2*Sqr(g1) + 6*QDxp*QQp*Sqr(gN)
      ) - 42.666666666666664*MassG*Sqr(g3)*(0.2*Sqr(g1) + 12*Sqr(gN)*Sqr(QQp)) -
      21.333333333333332*MassG*Sqr(g3)*(0.2*Sqr(g1) + 6*Qdp*QQp*Sqr(gN) + 2*Sqr(
      Yd02)) - 21.333333333333332*MassG*Sqr(g3)*(0.2*Sqr(g1) + 6*Qdp*QQp*Sqr(gN) +
      2*Sqr(Yd12)) - 21.333333333333332*MassG*Sqr(g3)*(0.2*Sqr(g1) + 6*Qdp*QQp*
      Sqr(gN) + 2*Sqr(Yd22)) - 21.333333333333332*MassG*Sqr(g3)*(-0.4*Sqr(g1) + 6*
      QQp*Qup*Sqr(gN) + 2*Sqr(Yu02)) - 21.333333333333332*MassG*Sqr(g3)*(-0.4*Sqr(
      g1) + 6*QQp*Qup*Sqr(gN) + 2*Sqr(Yu12)) - 21.333333333333332*MassG*Sqr(g3)*(
      -0.4*Sqr(g1) + 6*QQp*Qup*Sqr(gN) + 2*Sqr(Yu22)) - 21.333333333333332*MassG*
      Sqr(g3)*(0.2*Sqr(g1) + 12*Sqr(gN)*Sqr(QQp) + 2*Sqr(Yd02) + 2*Sqr(Yd12) + 2*
      Sqr(Yd22) + 2*Sqr(Yu02) + 2*Sqr(Yu12) + 2*Sqr(Yu22));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dlead_log_one_loop_mq222_dMassBp() const
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

   double deriv = -96*Power(gN,4)*MassBp*QQp*Power(QSp,3) - 48*MassBp*Sqr(gN)*(
      0.2*Sqr(g1) + 6*QDxbarp*QQp*Sqr(gN))*Sqr(QDxbarp) - 48*MassBp*Sqr(gN)*(-0.2*
      Sqr(g1) + 6*QDxp*QQp*Sqr(gN))*Sqr(QDxp) - 48*MassBp*Sqr(gN)*(0.2*Sqr(g1) + 2
      *Qep*QQp*Sqr(gN))*Sqr(Qep) - 32*MassBp*Sqr(gN)*(-0.2*Sqr(g1) + 4*QH1p*QQp*
      Sqr(gN))*Sqr(QH1p) - 32*MassBp*Sqr(gN)*(0.2*Sqr(g1) + 4*QH2p*QQp*Sqr(gN))*
      Sqr(QH2p) - 16*MassBp*Sqr(gN)*(0.2*Sqr(g1) + 4*QHpbarp*QQp*Sqr(gN))*Sqr(
      QHpbarp) - 16*MassBp*Sqr(gN)*(-0.2*Sqr(g1) + 4*QHpp*QQp*Sqr(gN))*Sqr(QHpp) -
      48*MassBp*Sqr(gN)*(-0.2*Sqr(g1) + 4*QLp*QQp*Sqr(gN))*Sqr(QLp) + 4*TYd02*
      Yd02*(4*Sqr(gN)*Sqr(Qdp) + 4*Sqr(gN)*Sqr(QH1p) + 4*Sqr(gN)*Sqr(QQp)) + 4*
      TYd12*Yd12*(4*Sqr(gN)*Sqr(Qdp) + 4*Sqr(gN)*Sqr(QH1p) + 4*Sqr(gN)*Sqr(QQp)) +
      4*TYd22*Yd22*(4*Sqr(gN)*Sqr(Qdp) + 4*Sqr(gN)*Sqr(QH1p) + 4*Sqr(gN)*Sqr(QQp)
      ) - 32*MassBp*Sqr(gN)*Sqr(QQp)*(0.2*Sqr(g1) + 12*Sqr(gN)*Sqr(QQp)) - 96*
      Power(gN,4)*MassBp*Sqr(QQp)*(9*Sqr(Qdp) + 9*Sqr(QDxbarp) + 9*Sqr(QDxp) + 3*
      Sqr(Qep) + 6*Sqr(QH1p) + 6*Sqr(QH2p) + 2*Sqr(QHpbarp) + 2*Sqr(QHpp) + 6*Sqr(
      QLp) + 18*Sqr(QQp) + 3*Sqr(QSp) + 9*Sqr(Qup)) + 4*TYu02*Yu02*(4*Sqr(gN)*Sqr(
      QH2p) + 4*Sqr(gN)*Sqr(QQp) + 4*Sqr(gN)*Sqr(Qup)) + 4*TYu12*Yu12*(4*Sqr(gN)*
      Sqr(QH2p) + 4*Sqr(gN)*Sqr(QQp) + 4*Sqr(gN)*Sqr(Qup)) + 4*TYu22*Yu22*(4*Sqr(
      gN)*Sqr(QH2p) + 4*Sqr(gN)*Sqr(QQp) + 4*Sqr(gN)*Sqr(Qup)) - 16*MassBp*Sqr(gN)
      *Sqr(Qdp)*(0.2*Sqr(g1) + 6*Qdp*QQp*Sqr(gN) + 2*Sqr(Yd02)) - 16*MassBp*Sqr(gN
      )*Sqr(Qdp)*(0.2*Sqr(g1) + 6*Qdp*QQp*Sqr(gN) + 2*Sqr(Yd12)) - 16*MassBp*Sqr(
      gN)*Sqr(Qdp)*(0.2*Sqr(g1) + 6*Qdp*QQp*Sqr(gN) + 2*Sqr(Yd22)) - 16*MassBp*Sqr
      (gN)*Sqr(QH1p)*(-0.2*Sqr(g1) + 4*QH1p*QQp*Sqr(gN) + 2*(Sqr(Yd02) + Sqr(Yd12)
      + Sqr(Yd22))) - 16*MassBp*Sqr(gN)*Sqr(Qup)*(-0.4*Sqr(g1) + 6*QQp*Qup*Sqr(gN
      ) + 2*Sqr(Yu02)) - 16*MassBp*Sqr(gN)*Sqr(Qup)*(-0.4*Sqr(g1) + 6*QQp*Qup*Sqr(
      gN) + 2*Sqr(Yu12)) - 16*MassBp*Sqr(gN)*Sqr(Qup)*(-0.4*Sqr(g1) + 6*QQp*Qup*
      Sqr(gN) + 2*Sqr(Yu22)) - 16*MassBp*Sqr(gN)*Sqr(QQp)*(0.2*Sqr(g1) + 12*Sqr(gN
      )*Sqr(QQp) + 2*Sqr(Yd02) + 2*Sqr(Yd12) + 2*Sqr(Yd22) + 2*Sqr(Yu02) + 2*Sqr(
      Yu12) + 2*Sqr(Yu22)) - 16*MassBp*Sqr(gN)*Sqr(QH2p)*(0.2*Sqr(g1) + 4*QH2p*QQp
      *Sqr(gN) + 2*(Sqr(Yu02) + Sqr(Yu12) + Sqr(Yu22)));


   return twoLoop*deriv;
}

} // namespace flexiblesusy
