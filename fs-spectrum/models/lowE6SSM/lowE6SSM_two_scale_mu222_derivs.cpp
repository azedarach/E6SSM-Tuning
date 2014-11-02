#include "lowE6SSM_two_scale_tuning_calculator.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

// one loop
double lowE6SSM_tuning_calculator::deriv_dbeta_one_loop_mu222_dLambdax() const
{

   double deriv = 0;


   return oneOver16PiSqr*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_one_loop_mu222_dALambdax() const
{

   double deriv = 0;


   return oneOver16PiSqr*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_one_loop_mu222_dAYu22() const
{

   const double Yu22 = model.get_Yu(2,2);

   const double AYu22 = get_AYu22();

   double deriv = 8*AYu22*Sqr(Yu22);


   return oneOver16PiSqr*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_one_loop_mu222_dmq222() const
{
   const lowE6SSM_input_parameters inputs = model.get_input();

   const auto QQp = inputs.QQp;
   const auto Qup = inputs.Qup;
   
   const double Yu22 = model.get_Yu(2,2);
   const double g1 = model.get_g1();
   const double gN = model.get_gN();

   double deriv = -0.8*Sqr(g1) + 12*QQp*Qup*Sqr(gN) + 4*Sqr(Yu22);


   return oneOver16PiSqr*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_one_loop_mu222_dmHd2() const
{
   const lowE6SSM_input_parameters inputs = model.get_input();

   const auto QH1p = inputs.QH1p;
   const auto Qup = inputs.Qup;
   
   const double g1 = model.get_g1();
   const double gN = model.get_gN();

   double deriv = 0.8*Sqr(g1) + 4*QH1p*Qup*Sqr(gN);


   return oneOver16PiSqr*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_one_loop_mu222_dmHu2() const
{
   const lowE6SSM_input_parameters inputs = model.get_input();

   const auto QH2p = inputs.QH2p;
   const auto Qup = inputs.Qup;
   
   const double Yu20 = model.get_Yu(2,0);
   const double Yu21 = model.get_Yu(2,1);
   const double Yu22 = model.get_Yu(2,2);
   const double g1 = model.get_g1();
   const double gN = model.get_gN();

   double deriv = -0.8*Sqr(g1) + 4*QH2p*Qup*Sqr(gN) + 4*(Sqr(Yu20) + Sqr(Yu21)
      + Sqr(Yu22));


   return oneOver16PiSqr*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_one_loop_mu222_dmu222() const
{
   const lowE6SSM_input_parameters inputs = model.get_input();

   const auto Qup = inputs.Qup;
   
   const double Yu20 = model.get_Yu(2,0);
   const double Yu21 = model.get_Yu(2,1);
   const double Yu22 = model.get_Yu(2,2);
   const double g1 = model.get_g1();
   const double gN = model.get_gN();

   double deriv = 1.6*Sqr(g1) + 6*Sqr(gN)*Sqr(Qup) + 4*(Sqr(Yu20) + Sqr(Yu21) +
      Sqr(Yu22));


   return oneOver16PiSqr*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_one_loop_mu222_dms2() const
{
   const lowE6SSM_input_parameters inputs = model.get_input();

   const auto Qup = inputs.Qup;
   const auto QSp = inputs.QSp;
   
   const double gN = model.get_gN();

   double deriv = 2*QSp*Qup*Sqr(gN);


   return oneOver16PiSqr*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_one_loop_mu222_dMassB() const
{

   const double g1 = model.get_g1();

   const double MassB = model.get_MassB();

   double deriv = -4.266666666666667*MassB*Sqr(g1);


   return oneOver16PiSqr*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_one_loop_mu222_dMassWB() const
{

   double deriv = 0;


   return oneOver16PiSqr*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_one_loop_mu222_dMassG() const
{

   const double g3 = model.get_g3();

   const double MassG = model.get_MassG();

   double deriv = -21.333333333333332*MassG*Sqr(g3);


   return oneOver16PiSqr*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_one_loop_mu222_dMassBp() const
{
   const lowE6SSM_input_parameters inputs = model.get_input();

   const auto Qup = inputs.Qup;
   
   const double gN = model.get_gN();

   const double MassBp = model.get_MassBp();

   double deriv = -16*MassBp*Sqr(gN)*Sqr(Qup);


   return oneOver16PiSqr*deriv;
}

// two loop
double lowE6SSM_tuning_calculator::deriv_dbeta_two_loop_mu222_dLambdax() const
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

   double deriv = -16*ALambdax*Lambdax*(TYu20*Yu20 + TYu21*Yu21 + TYu22*Yu22) -
      8*Lambdax*(Yu20*(mq200*Yu20 + mq210*Yu21 + mq220*Yu22) + Yu21*(mq201*Yu20 +
      mq211*Yu21 + mq221*Yu22) + Yu22*(mq202*Yu20 + mq212*Yu21 + mq222*Yu22)) - 4
      *Lambdax*(Yu20*(mu220*Yu00 + mu221*Yu10 + mu222*Yu20) + Yu21*(mu220*Yu01 +
      mu221*Yu11 + mu222*Yu21) + Yu22*(mu220*Yu02 + mu221*Yu12 + mu222*Yu22)) -
      3.2*Lambdax*(mHd2 - mHu2)*Sqr(g1) - 16*Lambdax*(mHd2*QH1p + mHu2*QH2p + ms2*
      QSp)*Qup*Sqr(gN) - 8*Lambdax*(Sqr(TYu20) + Sqr(TYu21) + Sqr(TYu22)) - 8*
      Lambdax*mHd2*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)) - 16*Lambdax*mHu2*(Sqr(Yu20
      ) + Sqr(Yu21) + Sqr(Yu22)) - 8*Lambdax*ms2*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22
      )) - 8*Lambdax*Sqr(ALambdax)*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)) - 4*Lambdax
      *(mu202*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + mu212*(Yu10*Yu20 + Yu11*Yu21 +
      Yu12*Yu22) + mu222*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_two_loop_mu222_dALambdax() const
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

   double deriv = -8*(TYu20*Yu20 + TYu21*Yu21 + TYu22*Yu22)*Sqr(Lambdax) - 8*
      ALambdax*Sqr(Lambdax)*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_two_loop_mu222_dAYu22() const
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

   double deriv = -4*Yu22*(TYd02*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22) + TYd12*(
      Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22) + TYd22*(Yd20*Yu20 + Yd21*Yu21 + Yd22*
      Yu22)) - 4*(Yu20*(TYd02*Yd00*Yu22 + TYd12*Yd10*Yu22 + TYd22*Yd20*Yu22) +
      Yu21*(TYd02*Yd01*Yu22 + TYd12*Yd11*Yu22 + TYd22*Yd21*Yu22) + Yu22*(TYd02*
      Yd02*Yu22 + TYd12*Yd12*Yu22 + TYd22*Yd22*Yu22)) - 4*(TYu20*(Yd00*Yd02*Yu22 +
      Yd10*Yd12*Yu22 + Yd20*Yd22*Yu22) + TYu21*(Yd01*Yd02*Yu22 + Yd11*Yd12*Yu22 +
      Yd21*Yd22*Yu22) + Yu22*(Yd02*(TYu20*Yd00 + TYu21*Yd01 + AYu22*Yd02*Yu22) +
      Yd12*(TYu20*Yd10 + TYu21*Yd11 + AYu22*Yd12*Yu22) + Yd22*(TYu20*Yd20 + TYu21*
      Yd21 + AYu22*Yd22*Yu22)) + AYu22*Yu22*(Yu22*Sqr(Yd02) + Yu22*Sqr(Yd12) +
      Yu22*Sqr(Yd22))) - 8*Lambdax*TLambdax*Sqr(Yu22) - 1.6*AYu22*Sqr(g1)*Sqr(Yu22
      ) + 1.6*MassB*Sqr(g1)*Sqr(Yu22) + 24*AYu22*Sqr(g2)*Sqr(Yu22) - 24*MassWB*Sqr
      (g2)*Sqr(Yu22) - 8*AYu22*Sqr(Lambdax)*Sqr(Yu22) + 16*AYu22*Sqr(gN)*Sqr(QH2p)
      *Sqr(Yu22) - 8*MassBp*Sqr(gN)*Sqr(QH2p)*Sqr(Yu22) + 16*AYu22*Sqr(gN)*Sqr(QQp
      )*Sqr(Yu22) - 8*MassBp*Sqr(gN)*Sqr(QQp)*Sqr(Yu22) - 8*MassBp*Sqr(gN)*(Sqr(
      QH2p) + Sqr(QQp) - Sqr(Qup))*Sqr(Yu22) - 16*AYu22*Sqr(gN)*Sqr(Qup)*Sqr(Yu22)
      + 8*MassBp*Sqr(gN)*Sqr(Qup)*Sqr(Yu22) - 24*AYu22*Sqr(Yu22)*(Sqr(Yu20) + Sqr
      (Yu21) + Sqr(Yu22)) - 24*AYu22*Sqr(Yu22)*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02)
      + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)) -
      24*Sqr(Yu22)*(TYu20*Yu20 + TYu21*Yu21 + AYu22*Sqr(Yu22)) - 24*Sqr(Yu22)*(
      TYu00*Yu00 + TYu01*Yu01 + TYu02*Yu02 + TYu10*Yu10 + TYu11*Yu11 + TYu12*Yu12
      + TYu20*Yu20 + TYu21*Yu21 + AYu22*Sqr(Yu22)) - 4*(Yu22*(TYu02*Yu02*Yu22 +
      TYu12*Yu12*Yu22 + 2*AYu22*Power(Yu22,3)) + Yu20*(TYu02*Yu00*Yu22 + TYu12*
      Yu10*Yu22 + 2*AYu22*Yu20*Sqr(Yu22)) + Yu21*(TYu02*Yu01*Yu22 + TYu12*Yu11*
      Yu22 + 2*AYu22*Yu21*Sqr(Yu22))) - 4*(AYu22*Sqr(Yu22)*(Sqr(Yu20) + Sqr(Yu21)
      + Sqr(Yu22)) + Yu22*(TYu02*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + TYu12*(Yu10
      *Yu20 + Yu11*Yu21 + Yu12*Yu22) + AYu22*Yu22*(Sqr(Yu20) + Sqr(Yu21) + Sqr(
      Yu22)))) - 4*(TYu20*Yu20*Sqr(Yu22) + TYu21*Yu21*Sqr(Yu22) + Yu22*(AYu22*
      Power(Yu22,3) + Yu22*(TYu20*Yu20 + TYu21*Yu21 + AYu22*Sqr(Yu22)))) - 4*(
      AYu22*Yu22*(Power(Yu22,3) + Yu22*Sqr(Yu02) + Yu22*Sqr(Yu12)) + TYu20*(Yu00*
      Yu02*Yu22 + Yu10*Yu12*Yu22 + Yu20*Sqr(Yu22)) + TYu21*(Yu01*Yu02*Yu22 + Yu11*
      Yu12*Yu22 + Yu21*Sqr(Yu22)) + Yu22*(Yu02*(TYu20*Yu00 + TYu21*Yu01 + AYu22*
      Yu02*Yu22) + Yu12*(TYu20*Yu10 + TYu21*Yu11 + AYu22*Yu12*Yu22) + Yu22*(TYu20*
      Yu20 + TYu21*Yu21 + AYu22*Sqr(Yu22))));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_two_loop_mu222_dmq222() const
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

   double deriv = 0.21333333333333335*Power(g1,4) + 10.666666666666666*Power(g3
      ,4) - 4*Yu22*(Yd02*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22) + Yd12*(Yd10*Yu20 +
      Yd11*Yu21 + Yd12*Yu22) + Yd22*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22)) - 6.4*QQp
      *Qup*Sqr(g1)*Sqr(gN) + 48*Power(gN,4)*Sqr(QQp)*Sqr(Qup) - 4*(Yu20*(Yd00*Yd02
      *Yu22 + Yd10*Yd12*Yu22 + Yd20*Yd22*Yu22) + Yu21*(Yd01*Yd02*Yu22 + Yd11*Yd12*
      Yu22 + Yd21*Yd22*Yu22) + Yu22*(Yu22*Sqr(Yd02) + Yu22*Sqr(Yd12) + Yu22*Sqr(
      Yd22))) - 0.8*Sqr(g1)*Sqr(Yu22) + 12*Sqr(g2)*Sqr(Yu22) - 4*Sqr(Lambdax)*Sqr(
      Yu22) + 8*Sqr(gN)*Sqr(QH2p)*Sqr(Yu22) + 8*Sqr(gN)*Sqr(QQp)*Sqr(Yu22) - 8*Sqr
      (gN)*Sqr(Qup)*Sqr(Yu22) - 12*(Sqr(Yu02) + Sqr(Yu12) + Sqr(Yu22))*(Sqr(Yu20)
      + Sqr(Yu21) + Sqr(Yu22)) - 12*Sqr(Yu22)*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) +
      Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)) -
      0.05333333333333334*Sqr(g1)*(Sqr(g1) + 45*Sqr(g2) + 80*Sqr(g3) + 60*Sqr(gN)*
      Sqr(QQp) - 30*(Sqr(Yd02) + Sqr(Yd12) + Sqr(Yd22)) - 30*(Sqr(Yu02) + Sqr(Yu12
      ) + Sqr(Yu22))) + 0.8*Qup*Sqr(gN)*(QQp*Sqr(g1) + 45*QQp*Sqr(g2) + 80*QQp*Sqr
      (g3) + 60*Power(QQp,3)*Sqr(gN) - 30*QQp*(Sqr(Yd02) + Sqr(Yd12) + Sqr(Yd22))
      - 30*QQp*(Sqr(Yu02) + Sqr(Yu12) + Sqr(Yu22))) - 4*Yu22*(Yu02*(Yu00*Yu20 +
      Yu01*Yu21 + Yu02*Yu22) + Yu12*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + Yu22*(
      Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) - 4*(Yu22*(Power(Yu22,3) + Yu22*Sqr(Yu02
      ) + Yu22*Sqr(Yu12)) + Yu20*(Yu00*Yu02*Yu22 + Yu10*Yu12*Yu22 + Yu20*Sqr(Yu22)
      ) + Yu21*(Yu01*Yu02*Yu22 + Yu11*Yu12*Yu22 + Yu21*Sqr(Yu22)));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_two_loop_mu222_dmHd2() const
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

   double deriv = 0.64*Power(g1,4) - 4*(Yu20*(Yd00*(Yd00*Yu20 + Yd01*Yu21 +
      Yd02*Yu22) + Yd10*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22) + Yd20*(Yd20*Yu20 +
      Yd21*Yu21 + Yd22*Yu22)) + Yu21*(Yd01*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22) +
      Yd11*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22) + Yd21*(Yd20*Yu20 + Yd21*Yu21 +
      Yd22*Yu22)) + Yu22*(Yd02*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22) + Yd12*(Yd10*
      Yu20 + Yd11*Yu21 + Yd12*Yu22) + Yd22*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22))) +
      6.4*QH1p*Qup*Sqr(g1)*Sqr(gN) + 16*Power(gN,4)*Sqr(QH1p)*Sqr(Qup) -
      0.05333333333333334*Sqr(g1)*(-9*Sqr(g1) - 45*Sqr(g2) + 30*Sqr(Lambdax) - 60*
      Sqr(gN)*Sqr(QH1p) + 90*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(
      Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + 30*(Sqr(Ye00) + Sqr
      (Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(
      Ye21) + Sqr(Ye22))) + 0.8*Qup*Sqr(gN)*(3*QH1p*Sqr(g1) + 15*QH1p*Sqr(g2) + 20
      *Power(QH1p,3)*Sqr(gN) - 10*QH1p*Sqr(Lambdax) - 30*QH1p*(Sqr(Yd00) + Sqr(
      Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21
      ) + Sqr(Yd22)) - 10*QH1p*(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) +
      Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22))) - 4*Sqr(Lambdax)
      *(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_two_loop_mu222_dmHu2() const
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

   double deriv = 0.64*Power(g1,4) - 4*(Yu20*(Yd00*(Yd00*Yu20 + Yd01*Yu21 +
      Yd02*Yu22) + Yd10*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22) + Yd20*(Yd20*Yu20 +
      Yd21*Yu21 + Yd22*Yu22)) + Yu21*(Yd01*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22) +
      Yd11*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22) + Yd21*(Yd20*Yu20 + Yd21*Yu21 +
      Yd22*Yu22)) + Yu22*(Yd02*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22) + Yd12*(Yd10*
      Yu20 + Yd11*Yu21 + Yd12*Yu22) + Yd22*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22))) -
      6.4*QH2p*Qup*Sqr(g1)*Sqr(gN) + 16*Power(gN,4)*Sqr(QH2p)*Sqr(Qup) - 0.8*Sqr(
      g1)*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)) + 12*Sqr(g2)*(Sqr(Yu20) + Sqr(Yu21)
      + Sqr(Yu22)) - 8*Sqr(Lambdax)*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)) + 8*Sqr(gN
      )*Sqr(QH2p)*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)) + 8*Sqr(gN)*Sqr(QQp)*(Sqr(
      Yu20) + Sqr(Yu21) + Sqr(Yu22)) - 8*Sqr(gN)*Sqr(Qup)*(Sqr(Yu20) + Sqr(Yu21) +
      Sqr(Yu22)) - 24*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))*(Sqr(Yu00) + Sqr(Yu01)
      + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) +
      Sqr(Yu22)) - 0.05333333333333334*Sqr(g1)*(9*Sqr(g1) + 45*Sqr(g2) - 30*Sqr(
      Lambdax) + 60*Sqr(gN)*Sqr(QH2p) - 90*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) +
      Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) +
      0.8*Qup*Sqr(gN)*(3*QH2p*Sqr(g1) + 15*QH2p*Sqr(g2) + 20*Power(QH2p,3)*Sqr(gN)
      - 10*QH2p*Sqr(Lambdax) - 30*QH2p*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(
      Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) - 8*(
      Yu20*(Yu00*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu10*(Yu10*Yu20 + Yu11*Yu21
      + Yu12*Yu22) + Yu20*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + Yu21*(Yu01*(Yu00
      *Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu11*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) +
      Yu21*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + Yu22*(Yu02*(Yu00*Yu20 + Yu01*
      Yu21 + Yu02*Yu22) + Yu12*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + Yu22*(Sqr(
      Yu20) + Sqr(Yu21) + Sqr(Yu22))));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_two_loop_mu222_dmu222() const
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

   double deriv = 1.7066666666666668*Power(g1,4) + 5.333333333333333*Power(g3,4
      ) + 24*Power(gN,4)*Power(Qup,4) - 4*(Yu20*(Yd00*(Yd00*Yu20 + Yd01*Yu21 +
      Yd02*Yu22) + Yd10*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22) + Yd20*(Yd20*Yu20 +
      Yd21*Yu21 + Yd22*Yu22)) + Yu21*(Yd01*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22) +
      Yd11*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22) + Yd21*(Yd20*Yu20 + Yd21*Yu21 +
      Yd22*Yu22)) + Yu22*(Yd02*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22) + Yd12*(Yd10*
      Yu20 + Yd11*Yu21 + Yd12*Yu22) + Yd22*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22))) +
      12.8*Sqr(g1)*Sqr(gN)*Sqr(Qup) - 0.8*Sqr(g1)*(Sqr(Yu20) + Sqr(Yu21) + Sqr(
      Yu22)) + 12*Sqr(g2)*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)) - 4*Sqr(Lambdax)*(
      Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)) + 8*Sqr(gN)*Sqr(QH2p)*(Sqr(Yu20) + Sqr(
      Yu21) + Sqr(Yu22)) + 8*Sqr(gN)*Sqr(QQp)*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))
      - 8*Sqr(gN)*Sqr(Qup)*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)) - 12*(Sqr(Yu20) +
      Sqr(Yu21) + Sqr(Yu22))*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(
      Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)) - 0.05333333333333334
      *Sqr(g1)*(-32*Sqr(g1) - 160*Sqr(g3) - 120*Sqr(gN)*Sqr(Qup) + 120*(Sqr(Yu20)
      + Sqr(Yu21) + Sqr(Yu22))) + 0.8*Qup*Sqr(gN)*(8*Qup*Sqr(g1) + 40*Qup*Sqr(g3)
      + 30*Power(Qup,3)*Sqr(gN) - 30*Qup*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) - 4*
      (Sqr(Yu20)*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)) + Sqr(Yu21)*(Sqr(Yu20) + Sqr(
      Yu21) + Sqr(Yu22)) + Sqr(Yu22)*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) - 4*(
      Yu20*(Yu00*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu10*(Yu10*Yu20 + Yu11*Yu21
      + Yu12*Yu22) + Yu20*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + Yu21*(Yu01*(Yu00
      *Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu11*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) +
      Yu21*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + Yu22*(Yu02*(Yu00*Yu20 + Yu01*
      Yu21 + Yu02*Yu22) + Yu12*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + Yu22*(Sqr(
      Yu20) + Sqr(Yu21) + Sqr(Yu22)))) - 12*Sqr(Power(Yu20,2) + Power(Yu21,2) +
      Power(Yu22,2));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_two_loop_mu222_dms2() const
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

   double deriv = 0.8*Qup*Sqr(gN)*(10*Power(QSp,3)*Sqr(gN) - 15*QSp*(Sqr(
      Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(
      Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) - 10*QSp*(Sqr(
      Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) - 10*QSp*
      Sqr(Lambdax)) + 8*Power(gN,4)*Sqr(QSp)*Sqr(Qup) - 4*Sqr(Lambdax)*(Sqr(Yu20)
      + Sqr(Yu21) + Sqr(Yu22));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_two_loop_mu222_dMassB() const
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

   double deriv = 0.8*(TYu20*Yu20 + TYu21*Yu21 + TYu22*Yu22)*Sqr(g1) +
      5.688888888888889*MassG*Sqr(g1)*Sqr(g3) - 1.0666666666666667*MassBp*(9*Qdp +
      9*QDxbarp - 9*QDxp + 9*Qep - 9*QH1p + 9*QH2p + 3*QHpbarp - 3*QHpp - 9*QLp +
      9*QQp - 22*Qup)*Qup*Sqr(g1)*Sqr(gN) + 0.017777777777777778*MassB*Sqr(g1)*(4
      *(912*Sqr(g1) + 5*(32*Sqr(g3) - 6*(9*Qdp + 9*QDxbarp - 9*QDxp + 9*Qep - 9*
      QH1p + 9*QH2p + 3*QHpbarp - 3*QHpp - 9*QLp + 9*QQp - 22*Qup)*Qup*Sqr(gN))) -
      90*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + 0.017777777777777778*Sqr(g1)*(4*(
      912*MassB*Sqr(g1) + 5*(16*(2*MassB + MassG)*Sqr(g3) - 3*(2*MassB + MassBp)*(
      9*Qdp + 9*QDxbarp - 9*QDxp + 9*Qep - 9*QH1p + 9*QH2p + 3*QHpbarp - 3*QHpp -
      9*QLp + 9*QQp - 22*Qup)*Qup*Sqr(gN))) + 45*(TYu20*Yu20 + TYu21*Yu21 + TYu22*
      Yu22 - 2*MassB*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_two_loop_mu222_dMassWB() const
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

   double deriv = -24*(TYu20*Yu20 + TYu21*Yu21 + TYu22*Yu22)*Sqr(g2) + 48*
      MassWB*Sqr(g2)*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_two_loop_mu222_dMassG() const
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

   double deriv = 106.66666666666667*Power(g3,4)*MassG + 11.377777777777778*
      MassB*Sqr(g1)*Sqr(g3) + 22.755555555555556*MassG*Sqr(g1)*Sqr(g3) +
      42.666666666666664*MassBp*Sqr(g3)*Sqr(gN)*Sqr(Qup) + 85.33333333333333*MassG
      *Sqr(g3)*Sqr(gN)*Sqr(Qup);


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_two_loop_mu222_dMassBp() const
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

   double deriv = -1.0666666666666667*MassB*(9*Qdp + 9*QDxbarp - 9*QDxp + 9*Qep
      - 9*QH1p + 9*QH2p + 3*QHpbarp - 3*QHpp - 9*QLp + 9*QQp - 22*Qup)*Qup*Sqr(g1
      )*Sqr(gN) - 8*(TYu20*Yu20 + TYu21*Yu21 + TYu22*Yu22)*Sqr(gN)*Sqr(QH2p) - 8*(
      TYu20*Yu20 + TYu21*Yu21 + TYu22*Yu22)*Sqr(gN)*Sqr(QQp) + 8*(TYu20*Yu20 +
      TYu21*Yu21 + TYu22*Yu22)*Sqr(gN)*Sqr(Qup) + 21.333333333333332*MassG*Sqr(g3)
      *Sqr(gN)*Sqr(Qup) + 0.5333333333333333*MassBp*Sqr(gN)*(Qup*(-4*(9*Qdp + 9*
      QDxbarp - 9*QDxp + 9*Qep - 9*QH1p + 9*QH2p + 3*QHpbarp - 3*QHpp - 9*QLp + 9*
      QQp - 22*Qup)*Sqr(g1) + 5*Qup*(16*Sqr(g3) + 9*Sqr(gN)*(9*Sqr(Qdp) + 9*Sqr(
      QDxbarp) + 9*Sqr(QDxp) + 3*Sqr(Qep) + 6*Sqr(QH1p) + 6*Sqr(QH2p) + 2*Sqr(
      QHpbarp) + 2*Sqr(QHpp) + 6*Sqr(QLp) + 18*Sqr(QQp) + 3*Sqr(QSp) + 11*Sqr(Qup)
      ))) + 30*(Sqr(QH2p) + Sqr(QQp) - Sqr(Qup))*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22
      ))) + 0.5333333333333333*Sqr(gN)*(Qup*(-2*(MassB + 2*MassBp)*(9*Qdp + 9*
      QDxbarp - 9*QDxp + 9*Qep - 9*QH1p + 9*QH2p + 3*QHpbarp - 3*QHpp - 9*QLp + 9*
      QQp - 22*Qup)*Sqr(g1) + 5*Qup*(8*(2*MassBp + MassG)*Sqr(g3) + 9*MassBp*Sqr(
      gN)*(9*Sqr(Qdp) + 9*Sqr(QDxbarp) + 9*Sqr(QDxp) + 3*Sqr(Qep) + 6*Sqr(QH1p) +
      6*Sqr(QH2p) + 2*Sqr(QHpbarp) + 2*Sqr(QHpp) + 6*Sqr(QLp) + 18*Sqr(QQp) + 3*
      Sqr(QSp) + 11*Sqr(Qup)))) + 15*(Sqr(QH2p) + Sqr(QQp) - Sqr(Qup))*(-(TYu20*
      Yu20) - TYu21*Yu21 - TYu22*Yu22 + 2*MassBp*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22
      ))));


   return twoLoop*deriv;
}

// leading log
double lowE6SSM_tuning_calculator::deriv_dlead_log_one_loop_mu222_dLambdax() const
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

   double deriv = 8*TYu20*(2*Lambdax*TYu20 + 4*ALambdax*Lambdax*Yu20) + 2*
      Lambdax*Yu00*(2*mu202*Yu20 + 2*mu220*Yu20) + 2*Lambdax*Yu10*(2*mu212*Yu20 +
      2*mu221*Yu20) + 8*TYu21*(2*Lambdax*TYu21 + 4*ALambdax*Lambdax*Yu21) + 2*
      Lambdax*Yu01*(2*mu202*Yu21 + 2*mu220*Yu21) + 2*Lambdax*Yu11*(2*mu212*Yu21 +
      2*mu221*Yu21) + 8*TYu22*(2*Lambdax*TYu22 + 4*ALambdax*Lambdax*Yu22) + 2*
      Lambdax*Yu02*(2*mu202*Yu22 + 2*mu220*Yu22) + 2*Lambdax*Yu12*(2*mu212*Yu22 +
      2*mu221*Yu22) + 2*Lambdax*Yu20*(8*mHu2*Yu20 + 2*(mu202*Yu00 + mu212*Yu10 + 2
      *mu222*Yu20) + 2*(mu220*Yu00 + mu221*Yu10 + 2*mu222*Yu20) + 4*(2*mq200*Yu20
      + mq201*Yu21 + mq210*Yu21 + mq202*Yu22 + mq220*Yu22)) + 2*Lambdax*Yu21*(8*
      mHu2*Yu21 + 2*(mu202*Yu01 + mu212*Yu11 + 2*mu222*Yu21) + 2*(mu220*Yu01 +
      mu221*Yu11 + 2*mu222*Yu21) + 4*(mq201*Yu20 + mq210*Yu20 + 2*mq211*Yu21 +
      mq212*Yu22 + mq221*Yu22)) + 2*Lambdax*Yu22*(8*mHu2*Yu22 + 4*(mq202*Yu20 +
      mq220*Yu20 + mq212*Yu21 + mq221*Yu21 + 2*mq222*Yu22) + 2*(mu202*Yu02 + mu212
      *Yu12 + 2*mu222*Yu22) + 2*(mu220*Yu02 + mu221*Yu12 + 2*mu222*Yu22)) + 2*QSp*
      Qup*(8*Lambdax*(mHd2 + mHu2 + ms2) + 8*Lambdax*Sqr(ALambdax))*Sqr(gN) + (4*
      Lambdax*mHd2 + 4*Lambdax*mHu2 + 4*Lambdax*ms2 + 4*Lambdax*Sqr(ALambdax))*(
      0.8*Sqr(g1) + 4*QH1p*Qup*Sqr(gN)) + (4*Lambdax*mHd2 + 4*Lambdax*mHu2 + 4*
      Lambdax*ms2 + 4*Lambdax*Sqr(ALambdax))*(-0.8*Sqr(g1) + 4*QH2p*Qup*Sqr(gN) +
      4*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dlead_log_one_loop_mu222_dALambdax() const
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

   double deriv = 16*TYu20*Yu20*Sqr(Lambdax) + 16*TYu21*Yu21*Sqr(Lambdax) + 16*
      TYu22*Yu22*Sqr(Lambdax) + 16*ALambdax*QSp*Qup*Sqr(gN)*Sqr(Lambdax) + 4*
      ALambdax*(0.8*Sqr(g1) + 4*QH1p*Qup*Sqr(gN))*Sqr(Lambdax) + 4*ALambdax*Sqr(
      Lambdax)*(-0.8*Sqr(g1) + 4*QH2p*Qup*Sqr(gN) + 4*(Sqr(Yu20) + Sqr(Yu21) + Sqr
      (Yu22)));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dlead_log_one_loop_mu222_dAYu22() const
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

   double deriv = 16*TYu02*Yu22*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + 16*TYu12*
      Yu22*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + 16*TYu20*Yu20*Sqr(Yu22) + 16*
      TYu21*Yu21*Sqr(Yu22) + 4*AYu22*Sqr(Yu22)*(-0.8*Sqr(g1) + 12*QQp*Qup*Sqr(gN)
      + 4*Sqr(Yu22)) + 12*AYu22*Sqr(Yu22)*(-0.8*Sqr(g1) + 4*QH2p*Qup*Sqr(gN) + 4*(
      Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + 8*AYu22*Sqr(Yu22)*(1.6*Sqr(g1) + 6*Sqr
      (gN)*Sqr(Qup) + 4*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + 8*AYu22*Yu22*(6*
      Power(Yu22,3) - 0.8666666666666667*Yu22*Sqr(g1) - 3*Yu22*Sqr(g2) -
      5.333333333333333*Yu22*Sqr(g3) + Yu22*Sqr(Lambdax) - 2*Yu22*Sqr(gN)*Sqr(QH2p
      ) - 2*Yu22*Sqr(gN)*Sqr(QQp) - 2*Yu22*Sqr(gN)*Sqr(Qup) + Yu22*Sqr(Yd02) +
      Yu22*Sqr(Yd12) + Yu22*Sqr(Yd22) + 5*(Power(Yu22,3) + Yu22*Sqr(Yu02) + Yu22*
      Sqr(Yu12)) + 4*Yu22*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)) + 3*Yu22*(Sqr(Yu00)
      + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) +
      Sqr(Yu21) + Sqr(Yu22))) + 8*TYu20*(Yd00*Yd02*Yu22 + Yd10*Yd12*Yu22 + Yd20*
      Yd22*Yu22 + 6*Yu20*Sqr(Yu22) + 5*(Yu00*Yu02*Yu22 + Yu10*Yu12*Yu22 + Yu20*Sqr
      (Yu22))) + 8*TYu21*(Yd01*Yd02*Yu22 + Yd11*Yd12*Yu22 + Yd21*Yd22*Yu22 + 6*
      Yu21*Sqr(Yu22) + 5*(Yu01*Yu02*Yu22 + Yu11*Yu12*Yu22 + Yu21*Sqr(Yu22))) + 8*
      Yu22*(Yd02*(TYu20*Yd00 + TYu21*Yd01 + AYu22*Yd02*Yu22) + Yd12*(TYu20*Yd10 +
      TYu21*Yd11 + AYu22*Yd12*Yu22) + Yd22*(TYu20*Yd20 + TYu21*Yd21 + AYu22*Yd22*
      Yu22) + 2*(TYd02*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22) + TYd12*(Yd10*Yu20 +
      Yd11*Yu21 + Yd12*Yu22) + TYd22*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22)) -
      0.8666666666666667*AYu22*Yu22*Sqr(g1) - 3*AYu22*Yu22*Sqr(g2) -
      5.333333333333333*AYu22*Yu22*Sqr(g3) + AYu22*Yu22*Sqr(Lambdax) - 2*AYu22*
      Yu22*Sqr(gN)*Sqr(QH2p) - 2*AYu22*Yu22*Sqr(gN)*Sqr(QQp) - 2*AYu22*Yu22*Sqr(gN
      )*Sqr(Qup) + 3*AYu22*Yu22*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) +
      Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)) + 4*(TYu02*(Yu00*
      Yu20 + Yu01*Yu21 + Yu02*Yu22) + TYu12*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) +
      AYu22*Yu22*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + 5*(Yu02*(TYu20*Yu00 +
      TYu21*Yu01 + AYu22*Yu02*Yu22) + Yu12*(TYu20*Yu10 + TYu21*Yu11 + AYu22*Yu12*
      Yu22) + Yu22*(TYu20*Yu20 + TYu21*Yu21 + AYu22*Sqr(Yu22))) + Yu22*(2*Lambdax*
      TLambdax + 1.7333333333333334*MassB*Sqr(g1) + 6*MassWB*Sqr(g2) +
      10.666666666666666*MassG*Sqr(g3) + 4*MassBp*Sqr(gN)*Sqr(QH2p) + 4*MassBp*Sqr
      (gN)*Sqr(QQp) + 4*MassBp*Sqr(gN)*Sqr(Qup) + 6*(TYu00*Yu00 + TYu01*Yu01 +
      TYu02*Yu02 + TYu10*Yu10 + TYu11*Yu11 + TYu12*Yu12 + TYu20*Yu20 + TYu21*Yu21
      + AYu22*Sqr(Yu22))));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dlead_log_one_loop_mu222_dmq222() const
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

   double deriv = -15.36*Power(g1,4) + 16*Yu02*Yu22*(Yu00*Yu20 + Yu01*Yu21 +
      Yu02*Yu22) + 16*Yu12*Yu22*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + 8*Yu20*Yu22*
      (Yd00*Yd02 + Yd10*Yd12 + Yd20*Yd22 + Yu00*Yu02 + Yu10*Yu12 + Yu20*Yu22) + 8*
      Yu21*Yu22*(Yd01*Yd02 + Yd11*Yd12 + Yd21*Yd22 + Yu01*Yu02 + Yu11*Yu12 + Yu21*
      Yu22) + 3*(0.4*Sqr(g1) + 12*QDxbarp*QQp*Sqr(gN))*(-0.8*Sqr(g1) + 6*QDxbarp*
      Qup*Sqr(gN)) + 3*(-0.4*Sqr(g1) + 12*QDxp*QQp*Sqr(gN))*(0.8*Sqr(g1) + 6*QDxp*
      Qup*Sqr(gN)) + 3*(1.2*Sqr(g1) + 12*Qep*QQp*Sqr(gN))*(-0.8*Sqr(g1) + 2*Qep*
      Qup*Sqr(gN)) + 2*(-0.6*Sqr(g1) + 12*QH1p*QQp*Sqr(gN))*(0.8*Sqr(g1) + 4*QH1p*
      Qup*Sqr(gN)) + 2*(0.6*Sqr(g1) + 12*QH2p*QQp*Sqr(gN))*(-0.8*Sqr(g1) + 4*QH2p*
      Qup*Sqr(gN)) + (0.6*Sqr(g1) + 12*QHpbarp*QQp*Sqr(gN))*(-0.8*Sqr(g1) + 4*
      QHpbarp*Qup*Sqr(gN)) + (-0.6*Sqr(g1) + 12*QHpp*QQp*Sqr(gN))*(0.8*Sqr(g1) + 4
      *QHpp*Qup*Sqr(gN)) + 3*(-0.6*Sqr(g1) + 12*QLp*QQp*Sqr(gN))*(0.8*Sqr(g1) + 4*
      QLp*Qup*Sqr(gN)) + 72*Power(gN,4)*QQp*Qup*Sqr(QSp) + 24*Power(gN,4)*QQp*Qup*
      (9*Sqr(Qdp) + 9*Sqr(QDxbarp) + 9*Sqr(QDxp) + 3*Sqr(Qep) + 6*Sqr(QH1p) + 6*
      Sqr(QH2p) + 2*Sqr(QHpbarp) + 2*Sqr(QHpp) + 6*Sqr(QLp) + 18*Sqr(QQp) + 3*Sqr(
      QSp) + 9*Sqr(Qup)) + (-0.8*Sqr(g1) + 6*Qdp*Qup*Sqr(gN))*(0.4*Sqr(g1) + 12*
      Qdp*QQp*Sqr(gN) + 4*Sqr(Yd02)) + (-0.8*Sqr(g1) + 6*Qdp*Qup*Sqr(gN))*(0.4*Sqr
      (g1) + 12*Qdp*QQp*Sqr(gN) + 4*Sqr(Yd12)) + (-0.8*Sqr(g1) + 6*Qdp*Qup*Sqr(gN)
      )*(0.4*Sqr(g1) + 12*Qdp*QQp*Sqr(gN) + 4*Sqr(Yd22)) + (0.8*Sqr(g1) + 4*QH1p*
      Qup*Sqr(gN))*(-0.6*Sqr(g1) + 12*QH1p*QQp*Sqr(gN) + 6*(Sqr(Yd02) + Sqr(Yd12)
      + Sqr(Yd22))) + (1.6*Sqr(g1) + 6*Sqr(gN)*Sqr(Qup))*(-0.8*Sqr(g1) + 12*QQp*
      Qup*Sqr(gN) + 4*Sqr(Yu02)) + (1.6*Sqr(g1) + 6*Sqr(gN)*Sqr(Qup))*(-0.8*Sqr(g1
      ) + 12*QQp*Qup*Sqr(gN) + 4*Sqr(Yu12)) + (0.2*Sqr(g1) + 12*Sqr(gN)*Sqr(QQp))*
      (-0.8*Sqr(g1) + 12*QQp*Qup*Sqr(gN) + 4*Sqr(Yu20)) + (0.2*Sqr(g1) + 12*Sqr(gN
      )*Sqr(QQp))*(-0.8*Sqr(g1) + 12*QQp*Qup*Sqr(gN) + 4*Sqr(Yu21)) + (0.2*Sqr(g1)
      + 12*Sqr(gN)*Sqr(QQp) + 2*Sqr(Yd02) + 2*Sqr(Yd12) + 2*Sqr(Yd22) + 2*Sqr(
      Yu02) + 2*Sqr(Yu12) + 2*Sqr(Yu22))*(-0.8*Sqr(g1) + 12*QQp*Qup*Sqr(gN) + 4*
      Sqr(Yu22)) + (0.6*Sqr(g1) + 12*QH2p*QQp*Sqr(gN) + 6*(Sqr(Yu02) + Sqr(Yu12) +
      Sqr(Yu22)))*(-0.8*Sqr(g1) + 4*QH2p*Qup*Sqr(gN) + 4*(Sqr(Yu20) + Sqr(Yu21) +
      Sqr(Yu22))) + (-0.8*Sqr(g1) + 12*QQp*Qup*Sqr(gN) + 4*Sqr(Yu22))*(1.6*Sqr(g1
      ) + 6*Sqr(gN)*Sqr(Qup) + 4*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + 8*Yu22*(
      Yd02*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22) + Yd12*(Yd10*Yu20 + Yd11*Yu21 +
      Yd12*Yu22) + Yd22*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22) + 3*(Yu02*(Yu00*Yu20 +
      Yu01*Yu21 + Yu02*Yu22) + Yu12*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + Yu22*(
      Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + Yu22*(-0.8666666666666667*Sqr(g1) - 3*
      Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2
      *Sqr(gN)*Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02
      ) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))))
      ;


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dlead_log_one_loop_mu222_dmHd2() const
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

   double deriv = 15.36*Power(g1,4) + 16*(Yd00*Yd01 + Yd10*Yd11 + Yd20*Yd21)*
      Yu20*Yu21 + 16*(Yd00*Yd02 + Yd10*Yd12 + Yd20*Yd22)*Yu20*Yu22 + 16*(Yd01*Yd02
      + Yd11*Yd12 + Yd21*Yd22)*Yu21*Yu22 + 3*(-0.4*Sqr(g1) + 4*QDxbarp*QH1p*Sqr(
      gN))*(-0.8*Sqr(g1) + 6*QDxbarp*Qup*Sqr(gN)) + 3*(0.4*Sqr(g1) + 4*QDxp*QH1p*
      Sqr(gN))*(0.8*Sqr(g1) + 6*QDxp*Qup*Sqr(gN)) + 2*(-0.6*Sqr(g1) + 4*QH1p*QH2p*
      Sqr(gN))*(-0.8*Sqr(g1) + 4*QH2p*Qup*Sqr(gN)) + (-0.6*Sqr(g1) + 4*QH1p*
      QHpbarp*Sqr(gN))*(-0.8*Sqr(g1) + 4*QHpbarp*Qup*Sqr(gN)) + (0.6*Sqr(g1) + 4*
      QH1p*QHpp*Sqr(gN))*(0.8*Sqr(g1) + 4*QHpp*Qup*Sqr(gN)) + 2*QSp*Qup*Sqr(gN)*(4
      *QH1p*QSp*Sqr(gN) + 4*Sqr(Lambdax)) + 2*(0.8*Sqr(g1) + 4*QH1p*Qup*Sqr(gN))*(
      0.6*Sqr(g1) + 4*Sqr(gN)*Sqr(QH1p)) + 16*Power(gN,4)*QH1p*Qup*Sqr(QSp) + 8*
      Power(gN,4)*QH1p*Qup*(9*Sqr(Qdp) + 9*Sqr(QDxbarp) + 9*Sqr(QDxp) + 3*Sqr(Qep)
      + 6*Sqr(QH1p) + 6*Sqr(QH2p) + 2*Sqr(QHpbarp) + 2*Sqr(QHpp) + 6*Sqr(QLp) +
      18*Sqr(QQp) + 3*Sqr(QSp) + 9*Sqr(Qup)) + 2*(0.8*Sqr(g1) + 4*QH1p*Qup*Sqr(gN)
      )*(1.6*Sqr(g1) + 6*Sqr(gN)*Sqr(Qup)) + (-0.8*Sqr(g1) + 6*Qdp*Qup*Sqr(gN))*(
      -0.4*Sqr(g1) + 4*Qdp*QH1p*Sqr(gN) + 4*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02))) +
      (-0.8*Sqr(g1) + 6*Qdp*Qup*Sqr(gN))*(-0.4*Sqr(g1) + 4*Qdp*QH1p*Sqr(gN) + 4*(
      Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12))) + (-0.8*Sqr(g1) + 6*Qdp*Qup*Sqr(gN))*(
      -0.4*Sqr(g1) + 4*Qdp*QH1p*Sqr(gN) + 4*(Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22))) +
      (-0.8*Sqr(g1) + 2*Qep*Qup*Sqr(gN))*(-1.2*Sqr(g1) + 4*Qep*QH1p*Sqr(gN) + 4*(
      Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02))) + (-0.8*Sqr(g1) + 2*Qep*Qup*Sqr(gN))*(
      -1.2*Sqr(g1) + 4*Qep*QH1p*Sqr(gN) + 4*(Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12))) +
      (0.8*Sqr(g1) + 4*QLp*Qup*Sqr(gN))*(0.6*Sqr(g1) + 4*QH1p*QLp*Sqr(gN) + 2*(
      Sqr(Ye00) + Sqr(Ye10) + Sqr(Ye20))) + (0.8*Sqr(g1) + 4*QLp*Qup*Sqr(gN))*(0.6
      *Sqr(g1) + 4*QH1p*QLp*Sqr(gN) + 2*(Sqr(Ye01) + Sqr(Ye11) + Sqr(Ye21))) + (
      0.8*Sqr(g1) + 4*QLp*Qup*Sqr(gN))*(0.6*Sqr(g1) + 4*QH1p*QLp*Sqr(gN) + 2*(Sqr(
      Ye02) + Sqr(Ye12) + Sqr(Ye22))) + (-0.8*Sqr(g1) + 2*Qep*Qup*Sqr(gN))*(-1.2*
      Sqr(g1) + 4*Qep*QH1p*Sqr(gN) + 4*(Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22))) + (0.8
      *Sqr(g1) + 4*QH1p*Qup*Sqr(gN))*(0.6*Sqr(g1) + 2*Sqr(Lambdax) + 4*Sqr(gN)*Sqr
      (QH1p) + 6*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(
      Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + 2*(Sqr(Ye00) + Sqr(Ye01) + Sqr(
      Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22
      ))) + (-0.2*Sqr(g1) + 4*QH1p*QQp*Sqr(gN) + 2*(Sqr(Yd00) + Sqr(Yd10) + Sqr(
      Yd20)))*(-0.8*Sqr(g1) + 12*QQp*Qup*Sqr(gN) + 4*Sqr(Yu20)) + (-0.2*Sqr(g1) +
      4*QH1p*QQp*Sqr(gN) + 2*(Sqr(Yd01) + Sqr(Yd11) + Sqr(Yd21)))*(-0.8*Sqr(g1) +
      12*QQp*Qup*Sqr(gN) + 4*Sqr(Yu21)) + (-0.2*Sqr(g1) + 4*QH1p*QQp*Sqr(gN) + 2*(
      Sqr(Yd02) + Sqr(Yd12) + Sqr(Yd22)))*(-0.8*Sqr(g1) + 12*QQp*Qup*Sqr(gN) + 4*
      Sqr(Yu22)) + (-0.6*Sqr(g1) + 4*QH1p*QH2p*Sqr(gN) + 2*Sqr(Lambdax))*(-0.8*Sqr
      (g1) + 4*QH2p*Qup*Sqr(gN) + 4*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + (0.8*
      Sqr(g1) + 4*QH1p*Qup*Sqr(gN))*(1.6*Sqr(g1) + 6*Sqr(gN)*Sqr(Qup) + 4*(Sqr(
      Yu20) + Sqr(Yu21) + Sqr(Yu22)));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dlead_log_one_loop_mu222_dmHu2() const
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

   double deriv = -15.36*Power(g1,4) + 16*Yu20*Yu21*(Yu00*Yu01 + Yu10*Yu11 +
      Yu20*Yu21) + 16*Yu20*Yu22*(Yu00*Yu02 + Yu10*Yu12 + Yu20*Yu22) + 16*Yu21*Yu22
      *(Yu01*Yu02 + Yu11*Yu12 + Yu21*Yu22) + 3*(0.4*Sqr(g1) + 4*Qdp*QH2p*Sqr(gN))*
      (-0.8*Sqr(g1) + 6*Qdp*Qup*Sqr(gN)) + 3*(0.4*Sqr(g1) + 4*QDxbarp*QH2p*Sqr(gN)
      )*(-0.8*Sqr(g1) + 6*QDxbarp*Qup*Sqr(gN)) + 3*(-0.4*Sqr(g1) + 4*QDxp*QH2p*Sqr
      (gN))*(0.8*Sqr(g1) + 6*QDxp*Qup*Sqr(gN)) + 3*(1.2*Sqr(g1) + 4*Qep*QH2p*Sqr(
      gN))*(-0.8*Sqr(g1) + 2*Qep*Qup*Sqr(gN)) + 2*(-0.6*Sqr(g1) + 4*QH1p*QH2p*Sqr(
      gN))*(0.8*Sqr(g1) + 4*QH1p*Qup*Sqr(gN)) + (0.6*Sqr(g1) + 4*QH2p*QHpbarp*Sqr(
      gN))*(-0.8*Sqr(g1) + 4*QHpbarp*Qup*Sqr(gN)) + (-0.6*Sqr(g1) + 4*QH2p*QHpp*
      Sqr(gN))*(0.8*Sqr(g1) + 4*QHpp*Qup*Sqr(gN)) + 3*(-0.6*Sqr(g1) + 4*QH2p*QLp*
      Sqr(gN))*(0.8*Sqr(g1) + 4*QLp*Qup*Sqr(gN)) + (0.8*Sqr(g1) + 4*QH1p*Qup*Sqr(
      gN))*(-0.6*Sqr(g1) + 4*QH1p*QH2p*Sqr(gN) + 2*Sqr(Lambdax)) + 2*QSp*Qup*Sqr(
      gN)*(4*QH2p*QSp*Sqr(gN) + 4*Sqr(Lambdax)) + 2*(-0.8*Sqr(g1) + 4*QH2p*Qup*Sqr
      (gN))*(0.6*Sqr(g1) + 4*Sqr(gN)*Sqr(QH2p)) + 16*Power(gN,4)*QH2p*Qup*Sqr(QSp)
      + 8*Power(gN,4)*QH2p*Qup*(9*Sqr(Qdp) + 9*Sqr(QDxbarp) + 9*Sqr(QDxp) + 3*Sqr
      (Qep) + 6*Sqr(QH1p) + 6*Sqr(QH2p) + 2*Sqr(QHpbarp) + 2*Sqr(QHpp) + 6*Sqr(QLp
      ) + 18*Sqr(QQp) + 3*Sqr(QSp) + 9*Sqr(Qup)) + (1.6*Sqr(g1) + 6*Sqr(gN)*Sqr(
      Qup))*(-0.8*Sqr(g1) + 4*QH2p*Qup*Sqr(gN) + 4*(Sqr(Yu00) + Sqr(Yu01) + Sqr(
      Yu02))) + (1.6*Sqr(g1) + 6*Sqr(gN)*Sqr(Qup))*(-0.8*Sqr(g1) + 4*QH2p*Qup*Sqr(
      gN) + 4*(Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12))) + (-0.8*Sqr(g1) + 12*QQp*Qup*
      Sqr(gN) + 4*Sqr(Yu20))*(0.2*Sqr(g1) + 4*QH2p*QQp*Sqr(gN) + 2*(Sqr(Yu00) +
      Sqr(Yu10) + Sqr(Yu20))) + (-0.8*Sqr(g1) + 12*QQp*Qup*Sqr(gN) + 4*Sqr(Yu21))*
      (0.2*Sqr(g1) + 4*QH2p*QQp*Sqr(gN) + 2*(Sqr(Yu01) + Sqr(Yu11) + Sqr(Yu21))) +
      (-0.8*Sqr(g1) + 12*QQp*Qup*Sqr(gN) + 4*Sqr(Yu22))*(0.2*Sqr(g1) + 4*QH2p*QQp
      *Sqr(gN) + 2*(Sqr(Yu02) + Sqr(Yu12) + Sqr(Yu22))) + (-0.8*Sqr(g1) + 4*QH2p*
      Qup*Sqr(gN) + 4*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))*(1.6*Sqr(g1) + 6*Sqr(gN
      )*Sqr(Qup) + 4*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + (-0.8*Sqr(g1) + 4*QH2p
      *Qup*Sqr(gN) + 4*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))*(0.6*Sqr(g1) + 2*Sqr(
      Lambdax) + 4*Sqr(gN)*Sqr(QH2p) + 6*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(
      Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + 8*Yu20
      *(Yd00*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22) + Yd10*(Yd10*Yu20 + Yd11*Yu21 +
      Yd12*Yu22) + Yd20*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22) + 3*(Yu00*(Yu00*Yu20 +
      Yu01*Yu21 + Yu02*Yu22) + Yu10*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + Yu20*(
      Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + Yu20*(-0.8666666666666667*Sqr(g1) - 3*
      Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2
      *Sqr(gN)*Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02
      ) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))))
      + 8*Yu21*(Yd01*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22) + Yd11*(Yd10*Yu20 + Yd11
      *Yu21 + Yd12*Yu22) + Yd21*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22) + 3*(Yu01*(
      Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu11*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22
      ) + Yu21*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + Yu21*(-0.8666666666666667*
      Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*
      Sqr(QH2p) - 2*Sqr(gN)*Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3*(Sqr(Yu00) + Sqr(
      Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21
      ) + Sqr(Yu22)))) + 8*Yu22*(Yd02*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22) + Yd12*(
      Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22) + Yd22*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22
      ) + 3*(Yu02*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu12*(Yu10*Yu20 + Yu11*
      Yu21 + Yu12*Yu22) + Yu22*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + Yu22*(
      -0.8666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(
      Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3
      *(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) +
      Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) + 16*Sqr(Yu00*Yu20 + Yu01*Yu21 + Yu02*
      Yu22) + 16*Sqr(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22);


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dlead_log_one_loop_mu222_dmu222() const
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

   double deriv = 30.72*Power(g1,4) + 3*(-0.8*Sqr(g1) + 2*Qep*Qup*Sqr(gN))*(
      -2.4*Sqr(g1) + 6*Qep*Qup*Sqr(gN)) + 3*(0.8*Sqr(g1) + 4*QH1p*Qup*Sqr(gN))*(
      1.2*Sqr(g1) + 6*QH1p*Qup*Sqr(gN)) + 2*(-0.8*Sqr(g1) + 4*QH2p*Qup*Sqr(gN))*(
      -1.2*Sqr(g1) + 6*QH2p*Qup*Sqr(gN)) + (-0.8*Sqr(g1) + 4*QHpbarp*Qup*Sqr(gN))*
      (-1.2*Sqr(g1) + 6*QHpbarp*Qup*Sqr(gN)) + (0.8*Sqr(g1) + 4*QHpp*Qup*Sqr(gN))*
      (1.2*Sqr(g1) + 6*QHpp*Qup*Sqr(gN)) + 3*(0.8*Sqr(g1) + 4*QLp*Qup*Sqr(gN))*(
      1.2*Sqr(g1) + 6*QLp*Qup*Sqr(gN)) + 36*Power(gN,4)*Sqr(QSp)*Sqr(Qup) + 12*
      Power(gN,4)*Sqr(Qup)*(9*Sqr(Qdp) + 9*Sqr(QDxbarp) + 9*Sqr(QDxp) + 3*Sqr(Qep)
      + 6*Sqr(QH1p) + 6*Sqr(QH2p) + 2*Sqr(QHpbarp) + 2*Sqr(QHpp) + 6*Sqr(QLp) +
      18*Sqr(QQp) + 3*Sqr(QSp) + 9*Sqr(Qup)) + 3*Sqr(-0.8*Power(g1,2) + 6*Power(gN
      ,2)*Qdp*Qup) + 3*Sqr(-0.8*Power(g1,2) + 6*Power(gN,2)*QDxbarp*Qup) + 3*Sqr(
      0.8*Power(g1,2) + 6*Power(gN,2)*QDxp*Qup) + 2*Sqr(1.6*Power(g1,2) + 6*Power(
      gN,2)*Power(Qup,2)) + (-0.4*Sqr(g1) + 6*QQp*Qup*Sqr(gN) + 2*Sqr(Yu20))*(-0.8
      *Sqr(g1) + 12*QQp*Qup*Sqr(gN) + 4*Sqr(Yu20)) + 16*Sqr(Yu20)*Sqr(Yu21) + (
      -0.4*Sqr(g1) + 6*QQp*Qup*Sqr(gN) + 2*Sqr(Yu21))*(-0.8*Sqr(g1) + 12*QQp*Qup*
      Sqr(gN) + 4*Sqr(Yu21)) + 16*Sqr(Yu20)*Sqr(Yu22) + 16*Sqr(Yu21)*Sqr(Yu22) + (
      -0.4*Sqr(g1) + 6*QQp*Qup*Sqr(gN) + 2*Sqr(Yu22))*(-0.8*Sqr(g1) + 12*QQp*Qup*
      Sqr(gN) + 4*Sqr(Yu22)) + (-0.8*Sqr(g1) + 4*QH2p*Qup*Sqr(gN) + 4*(Sqr(Yu20) +
      Sqr(Yu21) + Sqr(Yu22)))*(-1.2*Sqr(g1) + 6*QH2p*Qup*Sqr(gN) + 6*(Sqr(Yu20) +
      Sqr(Yu21) + Sqr(Yu22))) + 8*Yu20*(Yd00*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22)
      + Yd10*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22) + Yd20*(Yd20*Yu20 + Yd21*Yu21 +
      Yd22*Yu22) + 3*(Yu00*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu10*(Yu10*Yu20 +
      Yu11*Yu21 + Yu12*Yu22) + Yu20*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + Yu20*(
      -0.8666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(
      Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3
      *(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) +
      Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) + 8*Yu21*(Yd01*(Yd00*Yu20 + Yd01*Yu21 +
      Yd02*Yu22) + Yd11*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22) + Yd21*(Yd20*Yu20 +
      Yd21*Yu21 + Yd22*Yu22) + 3*(Yu01*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu11*
      (Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + Yu21*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22
      ))) + Yu21*(-0.8666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(
      g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QQp) - 2*Sqr(gN)*
      Sqr(Qup) + 3*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) +
      Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) + 8*Yu22*(Yd02*(Yd00*Yu20 +
      Yd01*Yu21 + Yd02*Yu22) + Yd12*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22) + Yd22*(
      Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22) + 3*(Yu02*(Yu00*Yu20 + Yu01*Yu21 + Yu02*
      Yu22) + Yu12*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + Yu22*(Sqr(Yu20) + Sqr(
      Yu21) + Sqr(Yu22))) + Yu22*(-0.8666666666666667*Sqr(g1) - 3*Sqr(g2) -
      5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*
      Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(
      Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) + 8*Sqr
      (Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + 8*Sqr(Yu10*Yu20 + Yu11*Yu21 + Yu12*
      Yu22) + Sqr(1.6*Power(g1,2) + 6*Power(gN,2)*Power(Qup,2) + 4*(Power(Yu20,2)
      + Power(Yu21,2) + Power(Yu22,2)));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dlead_log_one_loop_mu222_dms2() const
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

   double deriv = 8*Power(gN,4)*Power(QSp,3)*Qup + 6*Qdp*QSp*Sqr(gN)*(-0.8*Sqr(
      g1) + 6*Qdp*Qup*Sqr(gN)) + 6*Qep*QSp*Sqr(gN)*(-0.8*Sqr(g1) + 2*Qep*Qup*Sqr(
      gN)) + 2*QHpbarp*QSp*Sqr(gN)*(-0.8*Sqr(g1) + 4*QHpbarp*Qup*Sqr(gN)) + 2*QHpp
      *QSp*Sqr(gN)*(0.8*Sqr(g1) + 4*QHpp*Qup*Sqr(gN)) + 6*QLp*QSp*Sqr(gN)*(0.8*Sqr
      (g1) + 4*QLp*Qup*Sqr(gN)) + (0.8*Sqr(g1) + 6*QDxp*Qup*Sqr(gN))*(2*QDxp*QSp*
      Sqr(gN) + 2*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02))) + (0.8*Sqr(g1) + 6
      *QDxp*Qup*Sqr(gN))*(2*QDxp*QSp*Sqr(gN) + 2*(Sqr(Kappa10) + Sqr(Kappa11) +
      Sqr(Kappa12))) + (-0.8*Sqr(g1) + 6*QDxbarp*Qup*Sqr(gN))*(2*QDxbarp*QSp*Sqr(
      gN) + 2*(Sqr(Kappa00) + Sqr(Kappa10) + Sqr(Kappa20))) + (-0.8*Sqr(g1) + 6*
      QDxbarp*Qup*Sqr(gN))*(2*QDxbarp*QSp*Sqr(gN) + 2*(Sqr(Kappa01) + Sqr(Kappa11)
      + Sqr(Kappa21))) + (-0.8*Sqr(g1) + 6*QDxbarp*Qup*Sqr(gN))*(2*QDxbarp*QSp*
      Sqr(gN) + 2*(Sqr(Kappa02) + Sqr(Kappa12) + Sqr(Kappa22))) + (0.8*Sqr(g1) + 6
      *QDxp*Qup*Sqr(gN))*(2*QDxp*QSp*Sqr(gN) + 2*(Sqr(Kappa20) + Sqr(Kappa21) +
      Sqr(Kappa22))) + (-0.8*Sqr(g1) + 4*QH2p*Qup*Sqr(gN))*(2*QH2p*QSp*Sqr(gN) + 2
      *(Sqr(Lambda1200) + Sqr(Lambda1201))) + (0.8*Sqr(g1) + 4*QH1p*Qup*Sqr(gN))*(
      2*QH1p*QSp*Sqr(gN) + 2*(Sqr(Lambda1200) + Sqr(Lambda1210))) + (0.8*Sqr(g1) +
      4*QH1p*Qup*Sqr(gN))*(2*QH1p*QSp*Sqr(gN) + 2*(Sqr(Lambda1201) + Sqr(
      Lambda1211))) + (-0.8*Sqr(g1) + 4*QH2p*Qup*Sqr(gN))*(2*QH2p*QSp*Sqr(gN) + 2*
      (Sqr(Lambda1210) + Sqr(Lambda1211))) + (0.8*Sqr(g1) + 4*QH1p*Qup*Sqr(gN))*(2
      *QH1p*QSp*Sqr(gN) + 2*Sqr(Lambdax)) + 2*QSp*Qup*Sqr(gN)*(6*(Sqr(Kappa00) +
      Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) +
      Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 4*(Sqr(Lambda1200) + Sqr(
      Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 4*Sqr(Lambdax) + 2*Sqr(gN
      )*Sqr(QSp)) + 4*Power(gN,4)*QSp*Qup*(9*Sqr(Qdp) + 9*Sqr(QDxbarp) + 9*Sqr(
      QDxp) + 3*Sqr(Qep) + 6*Sqr(QH1p) + 6*Sqr(QH2p) + 2*Sqr(QHpbarp) + 2*Sqr(QHpp
      ) + 6*Sqr(QLp) + 18*Sqr(QQp) + 3*Sqr(QSp) + 9*Sqr(Qup)) + 4*QSp*Qup*Sqr(gN)*
      (1.6*Sqr(g1) + 6*Sqr(gN)*Sqr(Qup)) + 2*QQp*QSp*Sqr(gN)*(-0.8*Sqr(g1) + 12*
      QQp*Qup*Sqr(gN) + 4*Sqr(Yu20)) + 2*QQp*QSp*Sqr(gN)*(-0.8*Sqr(g1) + 12*QQp*
      Qup*Sqr(gN) + 4*Sqr(Yu21)) + 2*QQp*QSp*Sqr(gN)*(-0.8*Sqr(g1) + 12*QQp*Qup*
      Sqr(gN) + 4*Sqr(Yu22)) + (2*QH2p*QSp*Sqr(gN) + 2*Sqr(Lambdax))*(-0.8*Sqr(g1)
      + 4*QH2p*Qup*Sqr(gN) + 4*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + 2*QSp*Qup*
      Sqr(gN)*(1.6*Sqr(g1) + 6*Sqr(gN)*Sqr(Qup) + 4*(Sqr(Yu20) + Sqr(Yu21) + Sqr(
      Yu22)));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dlead_log_one_loop_mu222_dMassB() const
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

   double deriv = -245.76*Power(g1,4)*MassB + 13.866666666666667*TYu20*Yu20*Sqr
      (g1) + 13.866666666666667*TYu21*Yu21*Sqr(g1) + 13.866666666666667*TYu22*Yu22
      *Sqr(g1) - 3.2*MassB*Sqr(g1)*(-0.8*Sqr(g1) + 6*Qdp*Qup*Sqr(gN)) - 3.2*MassB*
      Sqr(g1)*(-0.8*Sqr(g1) + 6*QDxbarp*Qup*Sqr(gN)) - 3.2*MassB*Sqr(g1)*(0.8*Sqr(
      g1) + 6*QDxp*Qup*Sqr(gN)) - 28.8*MassB*Sqr(g1)*(-0.8*Sqr(g1) + 2*Qep*Qup*Sqr
      (gN)) - 7.2*MassB*Sqr(g1)*(0.8*Sqr(g1) + 4*QH1p*Qup*Sqr(gN)) - 4.8*MassB*Sqr
      (g1)*(-0.8*Sqr(g1) + 4*QH2p*Qup*Sqr(gN)) - 2.4*MassB*Sqr(g1)*(-0.8*Sqr(g1) +
      4*QHpbarp*Qup*Sqr(gN)) - 2.4*MassB*Sqr(g1)*(0.8*Sqr(g1) + 4*QHpp*Qup*Sqr(gN
      )) - 7.2*MassB*Sqr(g1)*(0.8*Sqr(g1) + 4*QLp*Qup*Sqr(gN)) - 8.533333333333333
      *MassB*Sqr(g1)*(1.6*Sqr(g1) + 6*Sqr(gN)*Sqr(Qup)) - 0.26666666666666666*
      MassB*Sqr(g1)*(-0.8*Sqr(g1) + 12*QQp*Qup*Sqr(gN) + 4*Sqr(Yu20)) -
      0.26666666666666666*MassB*Sqr(g1)*(-0.8*Sqr(g1) + 12*QQp*Qup*Sqr(gN) + 4*Sqr
      (Yu21)) - 0.26666666666666666*MassB*Sqr(g1)*(-0.8*Sqr(g1) + 12*QQp*Qup*Sqr(
      gN) + 4*Sqr(Yu22)) - 2.4*MassB*Sqr(g1)*(-0.8*Sqr(g1) + 4*QH2p*Qup*Sqr(gN) +
      4*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) - 4.266666666666667*MassB*Sqr(g1)*(
      1.6*Sqr(g1) + 6*Sqr(gN)*Sqr(Qup) + 4*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dlead_log_one_loop_mu222_dMassWB() const
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

   double deriv = 48*TYu20*Yu20*Sqr(g2) + 48*TYu21*Yu21*Sqr(g2) + 48*TYu22*Yu22
      *Sqr(g2) - 36*MassWB*Sqr(g2)*(0.8*Sqr(g1) + 4*QH1p*Qup*Sqr(gN)) - 24*MassWB*
      Sqr(g2)*(-0.8*Sqr(g1) + 4*QH2p*Qup*Sqr(gN)) - 12*MassWB*Sqr(g2)*(-0.8*Sqr(g1
      ) + 4*QHpbarp*Qup*Sqr(gN)) - 12*MassWB*Sqr(g2)*(0.8*Sqr(g1) + 4*QHpp*Qup*Sqr
      (gN)) - 36*MassWB*Sqr(g2)*(0.8*Sqr(g1) + 4*QLp*Qup*Sqr(gN)) - 12*MassWB*Sqr(
      g2)*(-0.8*Sqr(g1) + 12*QQp*Qup*Sqr(gN) + 4*Sqr(Yu20)) - 12*MassWB*Sqr(g2)*(
      -0.8*Sqr(g1) + 12*QQp*Qup*Sqr(gN) + 4*Sqr(Yu21)) - 12*MassWB*Sqr(g2)*(-0.8*
      Sqr(g1) + 12*QQp*Qup*Sqr(gN) + 4*Sqr(Yu22)) - 12*MassWB*Sqr(g2)*(-0.8*Sqr(g1
      ) + 4*QH2p*Qup*Sqr(gN) + 4*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dlead_log_one_loop_mu222_dMassG() const
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

   double deriv = 85.33333333333333*TYu20*Yu20*Sqr(g3) + 85.33333333333333*
      TYu21*Yu21*Sqr(g3) + 85.33333333333333*TYu22*Yu22*Sqr(g3) - 64*MassG*Sqr(g3)
      *(-0.8*Sqr(g1) + 6*Qdp*Qup*Sqr(gN)) - 64*MassG*Sqr(g3)*(-0.8*Sqr(g1) + 6*
      QDxbarp*Qup*Sqr(gN)) - 64*MassG*Sqr(g3)*(0.8*Sqr(g1) + 6*QDxp*Qup*Sqr(gN)) -
      42.666666666666664*MassG*Sqr(g3)*(1.6*Sqr(g1) + 6*Sqr(gN)*Sqr(Qup)) -
      21.333333333333332*MassG*Sqr(g3)*(-0.8*Sqr(g1) + 12*QQp*Qup*Sqr(gN) + 4*Sqr(
      Yu20)) - 21.333333333333332*MassG*Sqr(g3)*(-0.8*Sqr(g1) + 12*QQp*Qup*Sqr(gN)
      + 4*Sqr(Yu21)) - 21.333333333333332*MassG*Sqr(g3)*(-0.8*Sqr(g1) + 12*QQp*
      Qup*Sqr(gN) + 4*Sqr(Yu22)) - 21.333333333333332*MassG*Sqr(g3)*(1.6*Sqr(g1) +
      6*Sqr(gN)*Sqr(Qup) + 4*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dlead_log_one_loop_mu222_dMassBp() const
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

   double deriv = -96*Power(gN,4)*MassBp*Power(QSp,3)*Qup - 48*MassBp*Sqr(gN)*(
      -0.8*Sqr(g1) + 6*Qdp*Qup*Sqr(gN))*Sqr(Qdp) - 48*MassBp*Sqr(gN)*(-0.8*Sqr(g1)
      + 6*QDxbarp*Qup*Sqr(gN))*Sqr(QDxbarp) - 48*MassBp*Sqr(gN)*(0.8*Sqr(g1) + 6*
      QDxp*Qup*Sqr(gN))*Sqr(QDxp) - 48*MassBp*Sqr(gN)*(-0.8*Sqr(g1) + 2*Qep*Qup*
      Sqr(gN))*Sqr(Qep) - 48*MassBp*Sqr(gN)*(0.8*Sqr(g1) + 4*QH1p*Qup*Sqr(gN))*Sqr
      (QH1p) - 32*MassBp*Sqr(gN)*(-0.8*Sqr(g1) + 4*QH2p*Qup*Sqr(gN))*Sqr(QH2p) -
      16*MassBp*Sqr(gN)*(-0.8*Sqr(g1) + 4*QHpbarp*Qup*Sqr(gN))*Sqr(QHpbarp) - 16*
      MassBp*Sqr(gN)*(0.8*Sqr(g1) + 4*QHpp*Qup*Sqr(gN))*Sqr(QHpp) - 48*MassBp*Sqr(
      gN)*(0.8*Sqr(g1) + 4*QLp*Qup*Sqr(gN))*Sqr(QLp) - 96*Power(gN,4)*MassBp*Sqr(
      Qup)*(9*Sqr(Qdp) + 9*Sqr(QDxbarp) + 9*Sqr(QDxp) + 3*Sqr(Qep) + 6*Sqr(QH1p) +
      6*Sqr(QH2p) + 2*Sqr(QHpbarp) + 2*Sqr(QHpp) + 6*Sqr(QLp) + 18*Sqr(QQp) + 3*
      Sqr(QSp) + 9*Sqr(Qup)) + 8*TYu20*Yu20*(4*Sqr(gN)*Sqr(QH2p) + 4*Sqr(gN)*Sqr(
      QQp) + 4*Sqr(gN)*Sqr(Qup)) + 8*TYu21*Yu21*(4*Sqr(gN)*Sqr(QH2p) + 4*Sqr(gN)*
      Sqr(QQp) + 4*Sqr(gN)*Sqr(Qup)) + 8*TYu22*Yu22*(4*Sqr(gN)*Sqr(QH2p) + 4*Sqr(
      gN)*Sqr(QQp) + 4*Sqr(gN)*Sqr(Qup)) - 32*MassBp*Sqr(gN)*Sqr(Qup)*(1.6*Sqr(g1)
      + 6*Sqr(gN)*Sqr(Qup)) - 16*MassBp*Sqr(gN)*Sqr(QQp)*(-0.8*Sqr(g1) + 12*QQp*
      Qup*Sqr(gN) + 4*Sqr(Yu20)) - 16*MassBp*Sqr(gN)*Sqr(QQp)*(-0.8*Sqr(g1) + 12*
      QQp*Qup*Sqr(gN) + 4*Sqr(Yu21)) - 16*MassBp*Sqr(gN)*Sqr(QQp)*(-0.8*Sqr(g1) +
      12*QQp*Qup*Sqr(gN) + 4*Sqr(Yu22)) - 16*MassBp*Sqr(gN)*Sqr(QH2p)*(-0.8*Sqr(g1
      ) + 4*QH2p*Qup*Sqr(gN) + 4*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) - 16*MassBp*
      Sqr(gN)*Sqr(Qup)*(1.6*Sqr(g1) + 6*Sqr(gN)*Sqr(Qup) + 4*(Sqr(Yu20) + Sqr(Yu21
      ) + Sqr(Yu22)));


   return twoLoop*deriv;
}

} // namespace flexiblesusy
