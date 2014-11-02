#include "lowE6SSM_two_scale_tuning_calculator.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

// one loop
double lowE6SSM_tuning_calculator::deriv_dbeta_one_loop_mHu2_dLambdax() const
{

   const double Lambdax = model.get_Lambdax();

   const double ALambdax = get_ALambdax();

   const double mHd2 = model.get_mHd2();
   const double mHu2 = model.get_mHu2();
   const double ms2 = model.get_ms2();

   double deriv = 4*Lambdax*mHd2 + 4*Lambdax*mHu2 + 4*Lambdax*ms2 + 4*Lambdax*
      Sqr(ALambdax);


   return oneOver16PiSqr*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_one_loop_mHu2_dALambdax() const
{

   const double Lambdax = model.get_Lambdax();

   const double ALambdax = get_ALambdax();

   double deriv = 4*ALambdax*Sqr(Lambdax);


   return oneOver16PiSqr*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_one_loop_mHu2_dAYu22() const
{

   const double Yu22 = model.get_Yu(2,2);

   const double AYu22 = get_AYu22();

   double deriv = 12*AYu22*Sqr(Yu22);


   return oneOver16PiSqr*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_one_loop_mHu2_dmq222() const
{
   const lowE6SSM_input_parameters inputs = model.get_input();

   const auto QQp = inputs.QQp;
   const auto QH2p = inputs.QH2p;
   
   const double Yu02 = model.get_Yu(0,2);
   const double Yu12 = model.get_Yu(1,2);
   const double Yu22 = model.get_Yu(2,2);
   const double g1 = model.get_g1();
   const double gN = model.get_gN();

   double deriv = 0.6*Sqr(g1) + 12*QH2p*QQp*Sqr(gN) + 6*(Sqr(Yu02) + Sqr(Yu12)
      + Sqr(Yu22));


   return oneOver16PiSqr*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_one_loop_mHu2_dmHd2() const
{
   const lowE6SSM_input_parameters inputs = model.get_input();

   const auto QH1p = inputs.QH1p;
   const auto QH2p = inputs.QH2p;
   
   const double Lambdax = model.get_Lambdax();
   const double g1 = model.get_g1();
   const double gN = model.get_gN();

   double deriv = -0.6*Sqr(g1) + 4*QH1p*QH2p*Sqr(gN) + 2*Sqr(Lambdax);


   return oneOver16PiSqr*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_one_loop_mHu2_dmHu2() const
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

   double deriv = 0.6*Sqr(g1) + 2*Sqr(Lambdax) + 4*Sqr(gN)*Sqr(QH2p) + 6*(Sqr(
      Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20
      ) + Sqr(Yu21) + Sqr(Yu22));


   return oneOver16PiSqr*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_one_loop_mHu2_dmu222() const
{
   const lowE6SSM_input_parameters inputs = model.get_input();

   const auto QH2p = inputs.QH2p;
   const auto Qup = inputs.Qup;

   const double Yu20 = model.get_Yu(2,0);
   const double Yu21 = model.get_Yu(2,1);
   const double Yu22 = model.get_Yu(2,2);
   const double g1 = model.get_g1();
   const double gN = model.get_gN();

   double deriv = -1.2*Sqr(g1) + 6*QH2p*Qup*Sqr(gN) + 6*(Sqr(Yu20) + Sqr(Yu21)
      + Sqr(Yu22));


   return oneOver16PiSqr*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_one_loop_mHu2_dms2() const
{
   const lowE6SSM_input_parameters inputs = model.get_input();

   const auto QH2p = inputs.QH2p;
   const auto QSp = inputs.QSp;
   
   const double Lambdax = model.get_Lambdax();
   const double gN = model.get_gN();

   double deriv = 2*QH2p*QSp*Sqr(gN) + 2*Sqr(Lambdax);


   return oneOver16PiSqr*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_one_loop_mHu2_dMassB() const
{

   const double g1 = model.get_g1();

   const double MassB = model.get_MassB();

   double deriv = -2.4*MassB*Sqr(g1);


   return oneOver16PiSqr*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_one_loop_mHu2_dMassWB() const
{

   const double g2 = model.get_g2();

   const double MassWB = model.get_MassWB();

   double deriv = -12*MassWB*Sqr(g2);


   return oneOver16PiSqr*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_one_loop_mHu2_dMassG() const
{

   double deriv = 0;


   return oneOver16PiSqr*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_one_loop_mHu2_dMassBp() const
{
   const lowE6SSM_input_parameters inputs = model.get_input();

   const auto QH2p = inputs.QH2p;
   
   const double gN = model.get_gN();

   const double MassBp = model.get_MassBp();

   double deriv = -16*MassBp*Sqr(gN)*Sqr(QH2p);


   return oneOver16PiSqr*deriv;
}

// two loop
double lowE6SSM_tuning_calculator::deriv_dbeta_two_loop_mHu2_dLambdax() const
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

   double deriv = -12*Lambdax*(Kappa00*(Kappa00*mDxbar200 + Kappa01*mDxbar210 +
      Kappa02*mDxbar220) + Kappa10*(Kappa10*mDxbar200 + Kappa11*mDxbar210 +
      Kappa12*mDxbar220) + Kappa20*(Kappa20*mDxbar200 + Kappa21*mDxbar210 +
      Kappa22*mDxbar220) + Kappa01*(Kappa00*mDxbar201 + Kappa01*mDxbar211 +
      Kappa02*mDxbar221) + Kappa11*(Kappa10*mDxbar201 + Kappa11*mDxbar211 +
      Kappa12*mDxbar221) + Kappa21*(Kappa20*mDxbar201 + Kappa21*mDxbar211 +
      Kappa22*mDxbar221) + Kappa02*(Kappa00*mDxbar202 + Kappa01*mDxbar212 +
      Kappa02*mDxbar222) + Kappa12*(Kappa10*mDxbar202 + Kappa11*mDxbar212 +
      Kappa12*mDxbar222) + Kappa22*(Kappa20*mDxbar202 + Kappa21*mDxbar212 +
      Kappa22*mDxbar222)) - 8*Lambdax*(Lambda1200*(Lambda1200*mH1I200 + Lambda1201
      *mH1I201) + Lambda1210*(Lambda1210*mH1I200 + Lambda1211*mH1I201) +
      Lambda1201*(Lambda1200*mH1I210 + Lambda1201*mH1I211) + Lambda1211*(
      Lambda1210*mH1I210 + Lambda1211*mH1I211)) - 48*Power(Lambdax,3)*mHd2 - 48*
      Power(Lambdax,3)*mHu2 - 48*Power(Lambdax,3)*ms2 - 24*ALambdax*Lambdax*(
      Kappa00*TKappa00 + Kappa01*TKappa01 + Kappa02*TKappa02 + Kappa10*TKappa10 +
      Kappa11*TKappa11 + Kappa12*TKappa12 + Kappa20*TKappa20 + Kappa21*TKappa21 +
      Kappa22*TKappa22) - 16*ALambdax*Lambdax*(Lambda1200*TLambda1200 + Lambda1201
      *TLambda1201 + Lambda1210*TLambda1210 + Lambda1211*TLambda1211) - 24*
      ALambdax*Lambdax*(TYd00*Yd00 + TYd01*Yd01 + TYd02*Yd02 + TYd10*Yd10 + TYd11*
      Yd11 + TYd12*Yd12 + TYd20*Yd20 + TYd21*Yd21 + TYd22*Yd22) - 12*Lambdax*(Yd00
      *(md200*Yd00 + md201*Yd10 + md202*Yd20) + Yd10*(md210*Yd00 + md211*Yd10 +
      md212*Yd20) + Yd20*(md220*Yd00 + md221*Yd10 + md222*Yd20) + Yd01*(md200*Yd01
      + md201*Yd11 + md202*Yd21) + Yd11*(md210*Yd01 + md211*Yd11 + md212*Yd21) +
      Yd21*(md220*Yd01 + md221*Yd11 + md222*Yd21) + Yd02*(md200*Yd02 + md201*Yd12
      + md202*Yd22) + Yd12*(md210*Yd02 + md211*Yd12 + md212*Yd22) + Yd22*(md220*
      Yd02 + md221*Yd12 + md222*Yd22)) - 12*Lambdax*(Yd00*(mq200*Yd00 + mq201*Yd01
      + mq202*Yd02) + Yd01*(mq210*Yd00 + mq211*Yd01 + mq212*Yd02) + Yd02*(mq220*
      Yd00 + mq221*Yd01 + mq222*Yd02) + Yd10*(mq200*Yd10 + mq201*Yd11 + mq202*Yd12
      ) + Yd11*(mq210*Yd10 + mq211*Yd11 + mq212*Yd12) + Yd12*(mq220*Yd10 + mq221*
      Yd11 + mq222*Yd12) + Yd20*(mq200*Yd20 + mq201*Yd21 + mq202*Yd22) + Yd21*(
      mq210*Yd20 + mq211*Yd21 + mq212*Yd22) + Yd22*(mq220*Yd20 + mq221*Yd21 +
      mq222*Yd22)) - 8*ALambdax*Lambdax*(TYe00*Ye00 + TYe01*Ye01 + TYe02*Ye02 +
      TYe10*Ye10 + TYe11*Ye11 + TYe12*Ye12 + TYe20*Ye20 + TYe21*Ye21 + TYe22*Ye22)
      - 4*Lambdax*(Ye00*(me200*Ye00 + me201*Ye10 + me202*Ye20) + Ye10*(me210*Ye00
      + me211*Ye10 + me212*Ye20) + Ye20*(me220*Ye00 + me221*Ye10 + me222*Ye20) +
      Ye01*(me200*Ye01 + me201*Ye11 + me202*Ye21) + Ye11*(me210*Ye01 + me211*Ye11
      + me212*Ye21) + Ye21*(me220*Ye01 + me221*Ye11 + me222*Ye21) + Ye02*(me200*
      Ye02 + me201*Ye12 + me202*Ye22) + Ye12*(me210*Ye02 + me211*Ye12 + me212*Ye22
      ) + Ye22*(me220*Ye02 + me221*Ye12 + me222*Ye22)) - 4*Lambdax*(Ye00*(ml200*
      Ye00 + ml201*Ye01 + ml202*Ye02) + Ye01*(ml210*Ye00 + ml211*Ye01 + ml212*Ye02
      ) + Ye02*(ml220*Ye00 + ml221*Ye01 + ml222*Ye02) + Ye10*(ml200*Ye10 + ml201*
      Ye11 + ml202*Ye12) + Ye11*(ml210*Ye10 + ml211*Ye11 + ml212*Ye12) + Ye12*(
      ml220*Ye10 + ml221*Ye11 + ml222*Ye12) + Ye20*(ml200*Ye20 + ml201*Ye21 +
      ml202*Ye22) + Ye21*(ml210*Ye20 + ml211*Ye21 + ml212*Ye22) + Ye22*(ml220*Ye20
      + ml221*Ye21 + ml222*Ye22)) - 96*Power(Lambdax,3)*Sqr(ALambdax) + 2.4*
      Lambdax*(mHd2 - mHu2)*Sqr(g1) - 16*Lambdax*QH2p*(mHd2*QH1p + mHu2*QH2p + ms2
      *QSp)*Sqr(gN) - 12*Lambdax*mHd2*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02)
      + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) +
      Sqr(Kappa22)) - 12*Lambdax*mHu2*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02)
      + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21)
      + Sqr(Kappa22)) - 24*Lambdax*ms2*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02)
      + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21)
      + Sqr(Kappa22)) - 12*Lambdax*Sqr(ALambdax)*(Sqr(Kappa00) + Sqr(Kappa01) +
      Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) +
      Sqr(Kappa21) + Sqr(Kappa22)) - 12*Lambdax*((Kappa00*Kappa10 + Kappa01*
      Kappa11 + Kappa02*Kappa12)*mDx201 + (Kappa00*Kappa20 + Kappa01*Kappa21 +
      Kappa02*Kappa22)*mDx202 + (Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*
      Kappa12)*mDx210 + (Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)*
      mDx212 + (Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)*mDx220 + (
      Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)*mDx221 + mDx200*(Sqr(
      Kappa00) + Sqr(Kappa01) + Sqr(Kappa02)) + mDx211*(Sqr(Kappa10) + Sqr(Kappa11
      ) + Sqr(Kappa12)) + mDx222*(Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22))) - 8
      *Lambdax*mHd2*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(
      Lambda1211)) - 8*Lambdax*mHu2*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(
      Lambda1210) + Sqr(Lambda1211)) - 16*Lambdax*ms2*(Sqr(Lambda1200) + Sqr(
      Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) - 8*Lambdax*Sqr(ALambdax)*(
      Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) - 8*
      Lambdax*((Lambda1200*Lambda1210 + Lambda1201*Lambda1211)*mH2I201 + (
      Lambda1200*Lambda1210 + Lambda1201*Lambda1211)*mH2I210 + mH2I200*(Sqr(
      Lambda1200) + Sqr(Lambda1201)) + mH2I211*(Sqr(Lambda1210) + Sqr(Lambda1211))
      ) - 8*ALambdax*Lambdax*MassBp*Sqr(gN)*Sqr(QH1p) + 8*Lambdax*mHd2*Sqr(gN)*Sqr
      (QH1p) + 8*Lambdax*mHu2*Sqr(gN)*Sqr(QH1p) + 8*Lambdax*ms2*Sqr(gN)*Sqr(QH1p)
      + 8*Lambdax*Sqr(ALambdax)*Sqr(gN)*Sqr(QH1p) + 8*ALambdax*Lambdax*MassBp*Sqr(
      gN)*Sqr(QH2p) - 8*Lambdax*mHd2*Sqr(gN)*Sqr(QH2p) - 8*Lambdax*mHu2*Sqr(gN)*
      Sqr(QH2p) - 8*Lambdax*ms2*Sqr(gN)*Sqr(QH2p) - 8*Lambdax*Sqr(ALambdax)*Sqr(gN
      )*Sqr(QH2p) - 8*ALambdax*Lambdax*MassBp*Sqr(gN)*Sqr(QSp) + 8*Lambdax*mHd2*
      Sqr(gN)*Sqr(QSp) + 8*Lambdax*mHu2*Sqr(gN)*Sqr(QSp) + 8*Lambdax*ms2*Sqr(gN)*
      Sqr(QSp) + 8*Lambdax*Sqr(ALambdax)*Sqr(gN)*Sqr(QSp) + 0.8*MassBp*Sqr(gN)*(5*
      Lambdax*(-ALambdax + 2*MassBp)*(Sqr(QH1p) - Sqr(QH2p) + Sqr(QSp)) + 5*(-(
      ALambdax*Lambdax) + 2*Lambdax*MassBp)*(Sqr(QH1p) - Sqr(QH2p) + Sqr(QSp))) -
      12*Lambdax*(Sqr(TKappa00) + Sqr(TKappa01) + Sqr(TKappa02) + Sqr(TKappa10) +
      Sqr(TKappa11) + Sqr(TKappa12) + Sqr(TKappa20) + Sqr(TKappa21) + Sqr(TKappa22
      )) - 8*Lambdax*(Sqr(TLambda1200) + Sqr(TLambda1201) + Sqr(TLambda1210) + Sqr
      (TLambda1211)) - 12*Lambdax*(Sqr(TYd00) + Sqr(TYd01) + Sqr(TYd02) + Sqr(
      TYd10) + Sqr(TYd11) + Sqr(TYd12) + Sqr(TYd20) + Sqr(TYd21) + Sqr(TYd22)) - 4
      *Lambdax*(Sqr(TYe00) + Sqr(TYe01) + Sqr(TYe02) + Sqr(TYe10) + Sqr(TYe11) +
      Sqr(TYe12) + Sqr(TYe20) + Sqr(TYe21) + Sqr(TYe22)) - 24*Lambdax*mHd2*(Sqr(
      Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20
      ) + Sqr(Yd21) + Sqr(Yd22)) - 12*Lambdax*mHu2*(Sqr(Yd00) + Sqr(Yd01) + Sqr(
      Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22
      )) - 12*Lambdax*ms2*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(
      Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) - 12*Lambdax*Sqr(
      ALambdax)*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(
      Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) - 8*Lambdax*mHd2*(Sqr(Ye00) + Sqr
      (Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(
      Ye21) + Sqr(Ye22)) - 4*Lambdax*mHu2*(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr
      (Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22)) - 4*
      Lambdax*ms2*(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr
      (Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22)) - 4*Lambdax*Sqr(ALambdax)*(Sqr(
      Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20
      ) + Sqr(Ye21) + Sqr(Ye22));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_two_loop_mHu2_dALambdax() const
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

   double deriv = -48*ALambdax*Power(Lambdax,4) - 12*(Kappa00*TKappa00 +
      Kappa01*TKappa01 + Kappa02*TKappa02 + Kappa10*TKappa10 + Kappa11*TKappa11 +
      Kappa12*TKappa12 + Kappa20*TKappa20 + Kappa21*TKappa21 + Kappa22*TKappa22)*
      Sqr(Lambdax) - 8*(Lambda1200*TLambda1200 + Lambda1201*TLambda1201 +
      Lambda1210*TLambda1210 + Lambda1211*TLambda1211)*Sqr(Lambdax) - 12*(TYd00*
      Yd00 + TYd01*Yd01 + TYd02*Yd02 + TYd10*Yd10 + TYd11*Yd11 + TYd12*Yd12 +
      TYd20*Yd20 + TYd21*Yd21 + TYd22*Yd22)*Sqr(Lambdax) - 4*(TYe00*Ye00 + TYe01*
      Ye01 + TYe02*Ye02 + TYe10*Ye10 + TYe11*Ye11 + TYe12*Ye12 + TYe20*Ye20 +
      TYe21*Ye21 + TYe22*Ye22)*Sqr(Lambdax) - 12*ALambdax*(Sqr(Kappa00) + Sqr(
      Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(
      Kappa20) + Sqr(Kappa21) + Sqr(Kappa22))*Sqr(Lambdax) - 8*ALambdax*(Sqr(
      Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211))*Sqr(
      Lambdax) + 8*ALambdax*Sqr(gN)*Sqr(Lambdax)*Sqr(QH1p) - 4*MassBp*Sqr(gN)*Sqr(
      Lambdax)*Sqr(QH1p) - 8*ALambdax*Sqr(gN)*Sqr(Lambdax)*Sqr(QH2p) + 4*MassBp*
      Sqr(gN)*Sqr(Lambdax)*Sqr(QH2p) + 8*ALambdax*Sqr(gN)*Sqr(Lambdax)*Sqr(QSp) -
      4*MassBp*Sqr(gN)*Sqr(Lambdax)*Sqr(QSp) - 4*MassBp*Sqr(gN)*Sqr(Lambdax)*(Sqr(
      QH1p) - Sqr(QH2p) + Sqr(QSp)) - 12*ALambdax*Sqr(Lambdax)*(Sqr(Yd00) + Sqr(
      Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21
      ) + Sqr(Yd22)) - 4*ALambdax*Sqr(Lambdax)*(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02)
      + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_two_loop_mHu2_dAYu22() const
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

   double deriv = -6*Yu22*(TYd02*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22) + TYd12*(
      Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22) + TYd22*(Yd20*Yu20 + Yd21*Yu21 + Yd22*
      Yu22)) - 6*(TYd02*Yu22*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22) + TYd12*Yu22*(
      Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22) + TYd22*Yu22*(Yd20*Yu20 + Yd21*Yu21 +
      Yd22*Yu22)) + 3.2*AYu22*Sqr(g1)*Sqr(Yu22) - 3.2*MassB*Sqr(g1)*Sqr(Yu22) + 64
      *AYu22*Sqr(g3)*Sqr(Yu22) - 64*MassG*Sqr(g3)*Sqr(Yu22) - 24*AYu22*Sqr(gN)*Sqr
      (QH2p)*Sqr(Yu22) + 12*MassBp*Sqr(gN)*Sqr(QH2p)*Sqr(Yu22) + 24*AYu22*Sqr(gN)*
      Sqr(QQp)*Sqr(Yu22) - 12*MassBp*Sqr(gN)*Sqr(QQp)*Sqr(Yu22) + 12*MassBp*Sqr(gN
      )*(Sqr(QH2p) - Sqr(QQp) - Sqr(Qup))*Sqr(Yu22) + 24*AYu22*Sqr(gN)*Sqr(Qup)*
      Sqr(Yu22) - 12*MassBp*Sqr(gN)*Sqr(Qup)*Sqr(Yu22) - 6*(TYu20*Yd00*Yd02*Yu22 +
      TYu21*Yd01*Yd02*Yu22 + TYu20*Yd10*Yd12*Yu22 + TYu21*Yd11*Yd12*Yu22 + TYu20*
      Yd20*Yd22*Yu22 + TYu21*Yd21*Yd22*Yu22 + Yd02*(Yu22*(TYu20*Yd00 + TYu21*Yd01
      + AYu22*Yd02*Yu22) + AYu22*Yd02*Sqr(Yu22)) + Yd12*(Yu22*(TYu20*Yd10 + TYu21*
      Yd11 + AYu22*Yd12*Yu22) + AYu22*Yd12*Sqr(Yu22)) + Yd22*(Yu22*(TYu20*Yd20 +
      TYu21*Yd21 + AYu22*Yd22*Yu22) + AYu22*Yd22*Sqr(Yu22))) - 36*(TYu02*Yu22*(
      Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + TYu12*Yu22*(Yu10*Yu20 + Yu11*Yu21 +
      Yu12*Yu22) + AYu22*Sqr(Yu22)*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)) + Yu22*(
      TYu02*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + TYu12*(Yu10*Yu20 + Yu11*Yu21 +
      Yu12*Yu22) + AYu22*Yu22*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) - 36*(TYu20*
      Yu00*Yu02*Yu22 + TYu21*Yu01*Yu02*Yu22 + TYu20*Yu10*Yu12*Yu22 + TYu21*Yu11*
      Yu12*Yu22 + TYu20*Yu20*Sqr(Yu22) + TYu21*Yu21*Sqr(Yu22) + Yu02*(Yu22*(TYu20*
      Yu00 + TYu21*Yu01 + AYu22*Yu02*Yu22) + AYu22*Yu02*Sqr(Yu22)) + Yu12*(Yu22*(
      TYu20*Yu10 + TYu21*Yu11 + AYu22*Yu12*Yu22) + AYu22*Yu12*Sqr(Yu22)) + Yu22*(
      AYu22*Power(Yu22,3) + Yu22*(TYu20*Yu20 + TYu21*Yu21 + AYu22*Sqr(Yu22))));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_two_loop_mHu2_dmq222() const
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

   double deriv = 0.12*Power(g1,4) + 9*Power(g2,4) + 4.8*QH2p*QQp*Sqr(g1)*Sqr(
      gN) + 48*Power(gN,4)*Sqr(QH2p)*Sqr(QQp) - 6*Sqr(Lambdax)*(Sqr(Yd02) + Sqr(
      Yd12) + Sqr(Yd22)) - 6*(Yu02*((Yd00*Yd02 + Yd10*Yd12 + Yd20*Yd22)*Yu00 + (
      Yd01*Yd02 + Yd11*Yd12 + Yd21*Yd22)*Yu01 + Yu02*(Sqr(Yd02) + Sqr(Yd12) + Sqr(
      Yd22))) + Yu12*((Yd00*Yd02 + Yd10*Yd12 + Yd20*Yd22)*Yu10 + (Yd01*Yd02 + Yd11
      *Yd12 + Yd21*Yd22)*Yu11 + Yu12*(Sqr(Yd02) + Sqr(Yd12) + Sqr(Yd22))) + Yu22*(
      (Yd00*Yd02 + Yd10*Yd12 + Yd20*Yd22)*Yu20 + (Yd01*Yd02 + Yd11*Yd12 + Yd21*
      Yd22)*Yu21 + Yu22*(Sqr(Yd02) + Sqr(Yd12) + Sqr(Yd22)))) + 1.6*Sqr(g1)*(Sqr(
      Yu02) + Sqr(Yu12) + Sqr(Yu22)) + 32*Sqr(g3)*(Sqr(Yu02) + Sqr(Yu12) + Sqr(
      Yu22)) - 12*Sqr(gN)*Sqr(QH2p)*(Sqr(Yu02) + Sqr(Yu12) + Sqr(Yu22)) + 12*Sqr(
      gN)*Sqr(QQp)*(Sqr(Yu02) + Sqr(Yu12) + Sqr(Yu22)) + 12*Sqr(gN)*Sqr(Qup)*(Sqr(
      Yu02) + Sqr(Yu12) + Sqr(Yu22)) + 0.04*Sqr(g1)*(Sqr(g1) + 45*Sqr(g2) + 80*Sqr
      (g3) + 60*Sqr(gN)*Sqr(QQp) - 30*(Sqr(Yd02) + Sqr(Yd12) + Sqr(Yd22)) - 30*(
      Sqr(Yu02) + Sqr(Yu12) + Sqr(Yu22))) + 0.8*QH2p*Sqr(gN)*(QQp*Sqr(g1) + 45*QQp
      *Sqr(g2) + 80*QQp*Sqr(g3) + 60*Power(QQp,3)*Sqr(gN) - 30*QQp*(Sqr(Yd02) +
      Sqr(Yd12) + Sqr(Yd22)) - 30*QQp*(Sqr(Yu02) + Sqr(Yu12) + Sqr(Yu22))) - 6*(
      Yd02*(Yd00*(Yu00*Yu02 + Yu10*Yu12 + Yu20*Yu22) + Yd01*(Yu01*Yu02 + Yu11*Yu12
      + Yu21*Yu22) + Yd02*(Sqr(Yu02) + Sqr(Yu12) + Sqr(Yu22))) + Yd12*(Yd10*(Yu00
      *Yu02 + Yu10*Yu12 + Yu20*Yu22) + Yd11*(Yu01*Yu02 + Yu11*Yu12 + Yu21*Yu22) +
      Yd12*(Sqr(Yu02) + Sqr(Yu12) + Sqr(Yu22))) + Yd22*(Yd20*(Yu00*Yu02 + Yu10*
      Yu12 + Yu20*Yu22) + Yd21*(Yu01*Yu02 + Yu11*Yu12 + Yu21*Yu22) + Yd22*(Sqr(
      Yu02) + Sqr(Yu12) + Sqr(Yu22)))) - 36*(Yu02*(Yu00*(Yu00*Yu02 + Yu10*Yu12 +
      Yu20*Yu22) + Yu01*(Yu01*Yu02 + Yu11*Yu12 + Yu21*Yu22) + Yu02*(Sqr(Yu02) +
      Sqr(Yu12) + Sqr(Yu22))) + Yu12*(Yu10*(Yu00*Yu02 + Yu10*Yu12 + Yu20*Yu22) +
      Yu11*(Yu01*Yu02 + Yu11*Yu12 + Yu21*Yu22) + Yu12*(Sqr(Yu02) + Sqr(Yu12) + Sqr
      (Yu22))) + Yu22*(Yu20*(Yu00*Yu02 + Yu10*Yu12 + Yu20*Yu22) + Yu21*(Yu01*Yu02
      + Yu11*Yu12 + Yu21*Yu22) + Yu22*(Sqr(Yu02) + Sqr(Yu12) + Sqr(Yu22))));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_two_loop_mHu2_dmHd2() const
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

   double deriv = 0.36*Power(g1,4) + 3*Power(g2,4) - 12*Power(Lambdax,4) - 6*(
      Yd00*(Yu00*(Yd00*Yu00 + Yd01*Yu01 + Yd02*Yu02) + Yu10*(Yd00*Yu10 + Yd01*Yu11
      + Yd02*Yu12) + Yu20*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22)) + Yd01*(Yu01*(Yd00
      *Yu00 + Yd01*Yu01 + Yd02*Yu02) + Yu11*(Yd00*Yu10 + Yd01*Yu11 + Yd02*Yu12) +
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
      Yd21*Yu21 + Yd22*Yu22))) - 4.8*QH1p*QH2p*Sqr(g1)*Sqr(gN) - 6*(Sqr(Kappa00) +
      Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) +
      Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22))*Sqr(Lambdax) - 4*(Sqr(Lambda1200
      ) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211))*Sqr(Lambdax) + 4*
      Sqr(gN)*Sqr(Lambdax)*Sqr(QH1p) - 4*Sqr(gN)*Sqr(Lambdax)*Sqr(QH2p) + 16*Power
      (gN,4)*Sqr(QH1p)*Sqr(QH2p) + 4*Sqr(gN)*Sqr(Lambdax)*Sqr(QSp) - 12*Sqr(
      Lambdax)*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(
      Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) - 4*Sqr(Lambdax)*(Sqr(Ye00) + Sqr
      (Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(
      Ye21) + Sqr(Ye22)) + 0.04*Sqr(g1)*(-9*Sqr(g1) - 45*Sqr(g2) + 30*Sqr(Lambdax)
      - 60*Sqr(gN)*Sqr(QH1p) + 90*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10)
      + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + 30*(Sqr(Ye00)
      + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) +
      Sqr(Ye21) + Sqr(Ye22))) + 0.8*QH2p*Sqr(gN)*(3*QH1p*Sqr(g1) + 15*QH1p*Sqr(g2)
      + 20*Power(QH1p,3)*Sqr(gN) - 10*QH1p*Sqr(Lambdax) - 30*QH1p*(Sqr(Yd00) +
      Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(
      Yd21) + Sqr(Yd22)) - 10*QH1p*(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10)
      + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22)));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_two_loop_mHu2_dmHu2() const
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

   double deriv = 0.36*Power(g1,4) + 3*Power(g2,4) - 12*Power(Lambdax,4) + 16*
      Power(gN,4)*Power(QH2p,4) - 6*(Yd00*(Yu00*(Yd00*Yu00 + Yd01*Yu01 + Yd02*Yu02
      ) + Yu10*(Yd00*Yu10 + Yd01*Yu11 + Yd02*Yu12) + Yu20*(Yd00*Yu20 + Yd01*Yu21 +
      Yd02*Yu22)) + Yd01*(Yu01*(Yd00*Yu00 + Yd01*Yu01 + Yd02*Yu02) + Yu11*(Yd00*
      Yu10 + Yd01*Yu11 + Yd02*Yu12) + Yu21*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22)) +
      Yd02*(Yu02*(Yd00*Yu00 + Yd01*Yu01 + Yd02*Yu02) + Yu12*(Yd00*Yu10 + Yd01*Yu11
      + Yd02*Yu12) + Yu22*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22)) + Yd10*(Yu00*(Yd10
      *Yu00 + Yd11*Yu01 + Yd12*Yu02) + Yu10*(Yd10*Yu10 + Yd11*Yu11 + Yd12*Yu12) +
      Yu20*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22)) + Yd11*(Yu01*(Yd10*Yu00 + Yd11*
      Yu01 + Yd12*Yu02) + Yu11*(Yd10*Yu10 + Yd11*Yu11 + Yd12*Yu12) + Yu21*(Yd10*
      Yu20 + Yd11*Yu21 + Yd12*Yu22)) + Yd12*(Yu02*(Yd10*Yu00 + Yd11*Yu01 + Yd12*
      Yu02) + Yu12*(Yd10*Yu10 + Yd11*Yu11 + Yd12*Yu12) + Yu22*(Yd10*Yu20 + Yd11*
      Yu21 + Yd12*Yu22)) + Yd20*(Yu00*(Yd20*Yu00 + Yd21*Yu01 + Yd22*Yu02) + Yu10*(
      Yd20*Yu10 + Yd21*Yu11 + Yd22*Yu12) + Yu20*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22
      )) + Yd21*(Yu01*(Yd20*Yu00 + Yd21*Yu01 + Yd22*Yu02) + Yu11*(Yd20*Yu10 + Yd21
      *Yu11 + Yd22*Yu12) + Yu21*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22)) + Yd22*(Yu02*
      (Yd20*Yu00 + Yd21*Yu01 + Yd22*Yu02) + Yu12*(Yd20*Yu10 + Yd21*Yu11 + Yd22*
      Yu12) + Yu22*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22))) - 6*(Sqr(Kappa00) + Sqr(
      Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(
      Kappa20) + Sqr(Kappa21) + Sqr(Kappa22))*Sqr(Lambdax) - 4*(Sqr(Lambda1200) +
      Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211))*Sqr(Lambdax) + 4*Sqr(gN
      )*Sqr(Lambdax)*Sqr(QH1p) + 4.8*Sqr(g1)*Sqr(gN)*Sqr(QH2p) - 4*Sqr(gN)*Sqr(
      Lambdax)*Sqr(QH2p) + 4*Sqr(gN)*Sqr(Lambdax)*Sqr(QSp) - 6*Sqr(Lambdax)*(Sqr(
      Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20
      ) + Sqr(Yd21) + Sqr(Yd22)) - 2*Sqr(Lambdax)*(Sqr(Ye00) + Sqr(Ye01) + Sqr(
      Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22
      )) + 1.6*Sqr(g1)*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11)
      + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)) + 32*Sqr(g3)*(Sqr(Yu00) +
      Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(
      Yu21) + Sqr(Yu22)) - 12*Sqr(gN)*Sqr(QH2p)*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02)
      + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)) +
      12*Sqr(gN)*Sqr(QQp)*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(
      Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)) + 12*Sqr(gN)*Sqr(Qup)
      *(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) +
      Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)) + 0.04*Sqr(g1)*(9*Sqr(g1) + 45*Sqr(g2) -
      30*Sqr(Lambdax) + 60*Sqr(gN)*Sqr(QH2p) - 90*(Sqr(Yu00) + Sqr(Yu01) + Sqr(
      Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22
      ))) + 0.8*QH2p*Sqr(gN)*(3*QH2p*Sqr(g1) + 15*QH2p*Sqr(g2) + 20*Power(QH2p,3)*
      Sqr(gN) - 10*QH2p*Sqr(Lambdax) - 30*QH2p*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02)
      + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) -
      36*(Yu00*(Yu10*(Yu00*Yu10 + Yu01*Yu11 + Yu02*Yu12) + Yu20*(Yu00*Yu20 + Yu01*
      Yu21 + Yu02*Yu22) + Yu00*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02))) + Yu01*(Yu11*(
      Yu00*Yu10 + Yu01*Yu11 + Yu02*Yu12) + Yu21*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22
      ) + Yu01*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02))) + Yu02*(Yu12*(Yu00*Yu10 + Yu01
      *Yu11 + Yu02*Yu12) + Yu22*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu02*(Sqr(
      Yu00) + Sqr(Yu01) + Sqr(Yu02))) + Yu10*(Yu00*(Yu00*Yu10 + Yu01*Yu11 + Yu02*
      Yu12) + Yu20*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + Yu10*(Sqr(Yu10) + Sqr(
      Yu11) + Sqr(Yu12))) + Yu11*(Yu01*(Yu00*Yu10 + Yu01*Yu11 + Yu02*Yu12) + Yu21*
      (Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + Yu11*(Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12
      ))) + Yu12*(Yu02*(Yu00*Yu10 + Yu01*Yu11 + Yu02*Yu12) + Yu22*(Yu10*Yu20 +
      Yu11*Yu21 + Yu12*Yu22) + Yu12*(Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12))) + Yu20*(
      Yu00*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu10*(Yu10*Yu20 + Yu11*Yu21 +
      Yu12*Yu22) + Yu20*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + Yu21*(Yu01*(Yu00*
      Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu11*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) +
      Yu21*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + Yu22*(Yu02*(Yu00*Yu20 + Yu01*
      Yu21 + Yu02*Yu22) + Yu12*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + Yu22*(Sqr(
      Yu20) + Sqr(Yu21) + Sqr(Yu22))));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_two_loop_mHu2_dmu222() const
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

   double deriv = 0.96*Power(g1,4) - 6*(Yu20*(Yd00*(Yd00*Yu20 + Yd01*Yu21 +
      Yd02*Yu22) + Yd10*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22) + Yd20*(Yd20*Yu20 +
      Yd21*Yu21 + Yd22*Yu22)) + Yu21*(Yd01*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22) +
      Yd11*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22) + Yd21*(Yd20*Yu20 + Yd21*Yu21 +
      Yd22*Yu22)) + Yu22*(Yd02*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22) + Yd12*(Yd10*
      Yu20 + Yd11*Yu21 + Yd12*Yu22) + Yd22*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22))) -
      9.6*QH2p*Qup*Sqr(g1)*Sqr(gN) + 24*Power(gN,4)*Sqr(QH2p)*Sqr(Qup) + 1.6*Sqr(
      g1)*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)) + 32*Sqr(g3)*(Sqr(Yu20) + Sqr(Yu21)
      + Sqr(Yu22)) - 12*Sqr(gN)*Sqr(QH2p)*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)) + 12
      *Sqr(gN)*Sqr(QQp)*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)) + 12*Sqr(gN)*Sqr(Qup)*
      (Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)) + 0.04*Sqr(g1)*(-32*Sqr(g1) - 160*Sqr(g3
      ) - 120*Sqr(gN)*Sqr(Qup) + 120*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + 0.8*
      QH2p*Sqr(gN)*(8*Qup*Sqr(g1) + 40*Qup*Sqr(g3) + 30*Power(Qup,3)*Sqr(gN) - 30*
      Qup*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) - 36*(Yu20*(Yu00*(Yu00*Yu20 + Yu01*
      Yu21 + Yu02*Yu22) + Yu10*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + Yu20*(Sqr(
      Yu20) + Sqr(Yu21) + Sqr(Yu22))) + Yu21*(Yu01*(Yu00*Yu20 + Yu01*Yu21 + Yu02*
      Yu22) + Yu11*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + Yu21*(Sqr(Yu20) + Sqr(
      Yu21) + Sqr(Yu22))) + Yu22*(Yu02*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu12*
      (Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + Yu22*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22
      ))));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_two_loop_mHu2_dms2() const
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

   double deriv = -12*Power(Lambdax,4) - 12*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(
      Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(
      Kappa21) + Sqr(Kappa22))*Sqr(Lambdax) - 8*(Sqr(Lambda1200) + Sqr(Lambda1201)
      + Sqr(Lambda1210) + Sqr(Lambda1211))*Sqr(Lambdax) + 0.8*QH2p*Sqr(gN)*(10*
      Power(QSp,3)*Sqr(gN) - 15*QSp*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) +
      Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) +
      Sqr(Kappa22)) - 10*QSp*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210)
      + Sqr(Lambda1211)) - 10*QSp*Sqr(Lambdax)) + 4*Sqr(gN)*Sqr(Lambdax)*Sqr(QH1p)
      - 4*Sqr(gN)*Sqr(Lambdax)*Sqr(QH2p) + 4*Sqr(gN)*Sqr(Lambdax)*Sqr(QSp) + 8*
      Power(gN,4)*Sqr(QH2p)*Sqr(QSp) - 6*Sqr(Lambdax)*(Sqr(Yd00) + Sqr(Yd01) + Sqr
      (Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(
      Yd22)) - 2*Sqr(Lambdax)*(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr
      (Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_two_loop_mHu2_dMassB() const
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

   double deriv = -1.6*(TYu00*Yu00 + TYu01*Yu01 + TYu02*Yu02 + TYu10*Yu10 +
      TYu11*Yu11 + TYu12*Yu12 + TYu20*Yu20 + TYu21*Yu21 + TYu22*Yu22)*Sqr(g1) +
      1.8*MassWB*Sqr(g1)*Sqr(g2) + 2.4*MassBp*QH2p*(3*Qdp + 3*QDxbarp - 3*QDxp + 3
      *Qep - 3*QH1p + 4*QH2p + QHpbarp - QHpp - 3*QLp + 3*QQp - 6*Qup)*Sqr(g1)*Sqr
      (gN) + 0.04*MassB*Sqr(g1)*(3*(297*Sqr(g1) + 5*(6*Sqr(g2) + 8*QH2p*(3*Qdp + 3
      *QDxbarp - 3*QDxp + 3*Qep - 3*QH1p + 4*QH2p + QHpbarp - QHpp - 3*QLp + 3*QQp
      - 6*Qup)*Sqr(gN))) + 80*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) +
      Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + 0.04*Sqr(g1)*(
      -40*(TYu00*Yu00 + TYu01*Yu01 + TYu02*Yu02 + TYu10*Yu10 + TYu11*Yu11 + TYu12*
      Yu12 + TYu20*Yu20 + TYu21*Yu21 + TYu22*Yu22) + 3*(297*MassB*Sqr(g1) + 5*(3*(
      2*MassB + MassWB)*Sqr(g2) + 4*(2*MassB + MassBp)*QH2p*(3*Qdp + 3*QDxbarp - 3
      *QDxp + 3*Qep - 3*QH1p + 4*QH2p + QHpbarp - QHpp - 3*QLp + 3*QQp - 6*Qup)*
      Sqr(gN))) + 80*MassB*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(
      Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_two_loop_mHu2_dMassWB() const
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

   double deriv = 174*Power(g2,4)*MassWB + 3.6*MassB*Sqr(g1)*Sqr(g2) + 7.2*
      MassWB*Sqr(g1)*Sqr(g2) + 24*MassBp*Sqr(g2)*Sqr(gN)*Sqr(QH2p) + 48*MassWB*Sqr
      (g2)*Sqr(gN)*Sqr(QH2p);


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_two_loop_mHu2_dMassG() const
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

   double deriv = -64*(TYu00*Yu00 + TYu01*Yu01 + TYu02*Yu02 + TYu10*Yu10 +
      TYu11*Yu11 + TYu12*Yu12 + TYu20*Yu20 + TYu21*Yu21 + TYu22*Yu22)*Sqr(g3) +
      128*MassG*Sqr(g3)*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11)
      + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_two_loop_mHu2_dMassBp() const
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

   double deriv = 2.4*MassB*QH2p*(3*Qdp + 3*QDxbarp - 3*QDxp + 3*Qep - 3*QH1p +
      4*QH2p + QHpbarp - QHpp - 3*QLp + 3*QQp - 6*Qup)*Sqr(g1)*Sqr(gN) - 4*
      Lambdax*TLambdax*Sqr(gN)*Sqr(QH1p) + 4*Lambdax*TLambdax*Sqr(gN)*Sqr(QH2p) +
      12*(TYu00*Yu00 + TYu01*Yu01 + TYu02*Yu02 + TYu10*Yu10 + TYu11*Yu11 + TYu12*
      Yu12 + TYu20*Yu20 + TYu21*Yu21 + TYu22*Yu22)*Sqr(gN)*Sqr(QH2p) + 12*MassWB*
      Sqr(g2)*Sqr(gN)*Sqr(QH2p) - 12*(TYu00*Yu00 + TYu01*Yu01 + TYu02*Yu02 + TYu10
      *Yu10 + TYu11*Yu11 + TYu12*Yu12 + TYu20*Yu20 + TYu21*Yu21 + TYu22*Yu22)*Sqr(
      gN)*Sqr(QQp) - 4*Lambdax*TLambdax*Sqr(gN)*Sqr(QSp) - 12*(TYu00*Yu00 + TYu01*
      Yu01 + TYu02*Yu02 + TYu10*Yu10 + TYu11*Yu11 + TYu12*Yu12 + TYu20*Yu20 +
      TYu21*Yu21 + TYu22*Yu22)*Sqr(gN)*Sqr(Qup) + 0.8*MassBp*Sqr(gN)*(10*Sqr(
      Lambdax)*(Sqr(QH1p) - Sqr(QH2p) + Sqr(QSp)) + 3*(QH2p*(2*(3*Qdp + 3*QDxbarp
      - 3*QDxp + 3*Qep - 3*QH1p + 4*QH2p + QHpbarp - QHpp - 3*QLp + 3*QQp - 6*Qup)
      *Sqr(g1) + 5*QH2p*(2*Sqr(g2) + 2*Sqr(gN)*(9*Sqr(Qdp) + 9*Sqr(QDxbarp) + 9*
      Sqr(QDxp) + 3*Sqr(Qep) + 6*Sqr(QH1p) + 8*Sqr(QH2p) + 2*Sqr(QHpbarp) + 2*Sqr(
      QHpp) + 6*Sqr(QLp) + 18*Sqr(QQp) + 3*Sqr(QSp) + 9*Sqr(Qup)))) - 10*(Sqr(QH2p
      ) - Sqr(QQp) - Sqr(Qup))*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) +
      Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) + 0.8*Sqr(gN)*(
      5*Lambdax*(2*Lambdax*MassBp - TLambdax)*(Sqr(QH1p) - Sqr(QH2p) + Sqr(QSp)) +
      3*(5*(TYu00*Yu00 + TYu01*Yu01 + TYu02*Yu02 + TYu10*Yu10 + TYu11*Yu11 +
      TYu12*Yu12 + TYu20*Yu20 + TYu21*Yu21 + TYu22*Yu22)*(Sqr(QH2p) - Sqr(QQp) -
      Sqr(Qup)) + QH2p*((MassB + 2*MassBp)*(3*Qdp + 3*QDxbarp - 3*QDxp + 3*Qep - 3
      *QH1p + 4*QH2p + QHpbarp - QHpp - 3*QLp + 3*QQp - 6*Qup)*Sqr(g1) + 5*QH2p*((
      2*MassBp + MassWB)*Sqr(g2) + 2*MassBp*Sqr(gN)*(9*Sqr(Qdp) + 9*Sqr(QDxbarp) +
      9*Sqr(QDxp) + 3*Sqr(Qep) + 6*Sqr(QH1p) + 8*Sqr(QH2p) + 2*Sqr(QHpbarp) + 2*
      Sqr(QHpp) + 6*Sqr(QLp) + 18*Sqr(QQp) + 3*Sqr(QSp) + 9*Sqr(Qup)))) - 10*
      MassBp*(Sqr(QH2p) - Sqr(QQp) - Sqr(Qup))*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02)
      + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))));


   return twoLoop*deriv;
}

// leading log
double lowE6SSM_tuning_calculator::deriv_dlead_log_one_loop_mHu2_dLambdax() const
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

   double deriv = 12*TYu00*(2*Lambdax*TYu00 + 4*ALambdax*Lambdax*Yu00) + 12*
      TYu01*(2*Lambdax*TYu01 + 4*ALambdax*Lambdax*Yu01) + 12*TYu02*(2*Lambdax*
      TYu02 + 4*ALambdax*Lambdax*Yu02) + 12*TYu10*(2*Lambdax*TYu10 + 4*ALambdax*
      Lambdax*Yu10) + 12*TYu11*(2*Lambdax*TYu11 + 4*ALambdax*Lambdax*Yu11) + 12*
      TYu12*(2*Lambdax*TYu12 + 4*ALambdax*Lambdax*Yu12) + 12*TYu20*(2*Lambdax*
      TYu20 + 4*ALambdax*Lambdax*Yu20) + 2*Lambdax*Yu00*(12*mHu2*Yu00 + 6*(2*mq200
      *Yu00 + mq201*Yu01 + mq210*Yu01 + mq202*Yu02 + mq220*Yu02) + 6*(2*mu200*Yu00
      + mu201*Yu10 + mu210*Yu10 + mu202*Yu20 + mu220*Yu20)) + 2*Lambdax*Yu10*(12*
      mHu2*Yu10 + 6*(2*mq200*Yu10 + mq201*Yu11 + mq210*Yu11 + mq202*Yu12 + mq220*
      Yu12) + 6*(mu201*Yu00 + mu210*Yu00 + 2*mu211*Yu10 + mu212*Yu20 + mu221*Yu20)
      ) + 12*TYu21*(2*Lambdax*TYu21 + 4*ALambdax*Lambdax*Yu21) + 2*Lambdax*Yu01*(
      12*mHu2*Yu01 + 6*(mq201*Yu00 + mq210*Yu00 + 2*mq211*Yu01 + mq212*Yu02 +
      mq221*Yu02) + 6*(2*mu200*Yu01 + mu201*Yu11 + mu210*Yu11 + mu202*Yu21 + mu220
      *Yu21)) + 2*Lambdax*Yu11*(12*mHu2*Yu11 + 6*(mq201*Yu10 + mq210*Yu10 + 2*
      mq211*Yu11 + mq212*Yu12 + mq221*Yu12) + 6*(mu201*Yu01 + mu210*Yu01 + 2*mu211
      *Yu11 + mu212*Yu21 + mu221*Yu21)) + 12*TYu22*(2*Lambdax*TYu22 + 4*ALambdax*
      Lambdax*Yu22) + 2*Lambdax*Yu20*(12*mHu2*Yu20 + 6*(mu202*Yu00 + mu220*Yu00 +
      mu212*Yu10 + mu221*Yu10 + 2*mu222*Yu20) + 6*(2*mq200*Yu20 + mq201*Yu21 +
      mq210*Yu21 + mq202*Yu22 + mq220*Yu22)) + 2*Lambdax*Yu21*(12*mHu2*Yu21 + 6*(
      mu202*Yu01 + mu220*Yu01 + mu212*Yu11 + mu221*Yu11 + 2*mu222*Yu21) + 6*(mq201
      *Yu20 + mq210*Yu20 + 2*mq211*Yu21 + mq212*Yu22 + mq221*Yu22)) + 2*Lambdax*
      Yu02*(12*mHu2*Yu02 + 6*(mq202*Yu00 + mq220*Yu00 + mq212*Yu01 + mq221*Yu01 +
      2*mq222*Yu02) + 6*(2*mu200*Yu02 + mu201*Yu12 + mu210*Yu12 + mu202*Yu22 +
      mu220*Yu22)) + 2*Lambdax*Yu12*(12*mHu2*Yu12 + 6*(mq202*Yu10 + mq220*Yu10 +
      mq212*Yu11 + mq221*Yu11 + 2*mq222*Yu12) + 6*(mu201*Yu02 + mu210*Yu02 + 2*
      mu211*Yu12 + mu212*Yu22 + mu221*Yu22)) + 2*Lambdax*Yu22*(12*mHu2*Yu22 + 6*(
      mq202*Yu20 + mq220*Yu20 + mq212*Yu21 + mq221*Yu21 + 2*mq222*Yu22) + 6*(mu202
      *Yu02 + mu220*Yu02 + mu212*Yu12 + mu221*Yu12 + 2*mu222*Yu22)) + (4*Lambdax*
      mHd2 + 4*Lambdax*mHu2 + 4*Lambdax*ms2 + 4*Lambdax*Sqr(ALambdax))*(-0.6*Sqr(
      g1) + 4*QH1p*QH2p*Sqr(gN) + 2*Sqr(Lambdax)) + (8*Lambdax*(mHd2 + mHu2 + ms2)
      + 8*Lambdax*Sqr(ALambdax))*(2*QH2p*QSp*Sqr(gN) + 2*Sqr(Lambdax)) + 4*
      Lambdax*(6*(Kappa00*(Kappa00*mDxbar200 + Kappa01*mDxbar210 + Kappa02*
      mDxbar220) + Kappa10*(Kappa10*mDxbar200 + Kappa11*mDxbar210 + Kappa12*
      mDxbar220) + Kappa20*(Kappa20*mDxbar200 + Kappa21*mDxbar210 + Kappa22*
      mDxbar220) + Kappa01*(Kappa00*mDxbar201 + Kappa01*mDxbar211 + Kappa02*
      mDxbar221) + Kappa11*(Kappa10*mDxbar201 + Kappa11*mDxbar211 + Kappa12*
      mDxbar221) + Kappa21*(Kappa20*mDxbar201 + Kappa21*mDxbar211 + Kappa22*
      mDxbar221) + Kappa02*(Kappa00*mDxbar202 + Kappa01*mDxbar212 + Kappa02*
      mDxbar222) + Kappa12*(Kappa10*mDxbar202 + Kappa11*mDxbar212 + Kappa12*
      mDxbar222) + Kappa22*(Kappa20*mDxbar202 + Kappa21*mDxbar212 + Kappa22*
      mDxbar222)) + 4*(Lambda1200*(Lambda1200*mH1I200 + Lambda1201*mH1I201) +
      Lambda1210*(Lambda1210*mH1I200 + Lambda1211*mH1I201) + Lambda1201*(
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
      TLambda1201) + Sqr(TLambda1210) + Sqr(TLambda1211))) + 4*Lambdax*(6*(Yd00*(
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
      ))) + (4*Lambdax*mHd2 + 4*Lambdax*mHu2 + 4*Lambdax*ms2)*(-0.6*Sqr(g1) - 3*
      Sqr(g2) + 3*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr
      (Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*(
      Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 12*
      Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QSp
      ) + 3*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12)
      + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) +
      Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22) + 3*(
      Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(
      Yu20) + Sqr(Yu21) + Sqr(Yu22))) + (4*Lambdax*mHd2 + 4*Lambdax*mHu2 + 4*
      Lambdax*ms2 + 4*Lambdax*Sqr(ALambdax))*(0.6*Sqr(g1) + 2*Sqr(Lambdax) + 4*Sqr
      (gN)*Sqr(QH2p) + 6*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11
      ) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + (4*mHd2 + 4*mHu2 + 4*
      ms2)*(4*Power(Lambdax,3) - 0.6*Lambdax*Sqr(g1) - 3*Lambdax*Sqr(g2) + 3*
      Lambdax*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(
      Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*
      Lambdax*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(
      Lambda1211)) - 2*Lambdax*Sqr(gN)*Sqr(QH1p) - 2*Lambdax*Sqr(gN)*Sqr(QH2p) - 2
      *Lambdax*Sqr(gN)*Sqr(QSp) + 3*Lambdax*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) +
      Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) +
      Lambdax*(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(
      Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22)) + 3*Lambdax*(Sqr(Yu00) + Sqr(Yu01
      ) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) +
      Sqr(Yu22))) + 4*Lambdax*(6*(Yu00*(mq200*Yu00 + mq201*Yu01 + mq202*Yu02) +
      Yu01*(mq210*Yu00 + mq211*Yu01 + mq212*Yu02) + Yu02*(mq220*Yu00 + mq221*Yu01
      + mq222*Yu02) + Yu10*(mq200*Yu10 + mq201*Yu11 + mq202*Yu12) + Yu11*(mq210*
      Yu10 + mq211*Yu11 + mq212*Yu12) + Yu12*(mq220*Yu10 + mq221*Yu11 + mq222*Yu12
      ) + Yu20*(mq200*Yu20 + mq201*Yu21 + mq202*Yu22) + Yu21*(mq210*Yu20 + mq211*
      Yu21 + mq212*Yu22) + Yu22*(mq220*Yu20 + mq221*Yu21 + mq222*Yu22)) + 6*(Yu00*
      (mu200*Yu00 + mu201*Yu10 + mu202*Yu20) + Yu10*(mu210*Yu00 + mu211*Yu10 +
      mu212*Yu20) + Yu20*(mu220*Yu00 + mu221*Yu10 + mu222*Yu20) + Yu01*(mu200*Yu01
      + mu201*Yu11 + mu202*Yu21) + Yu11*(mu210*Yu01 + mu211*Yu11 + mu212*Yu21) +
      Yu21*(mu220*Yu01 + mu221*Yu11 + mu222*Yu21) + Yu02*(mu200*Yu02 + mu201*Yu12
      + mu202*Yu22) + Yu12*(mu210*Yu02 + mu211*Yu12 + mu212*Yu22) + Yu22*(mu220*
      Yu02 + mu221*Yu12 + mu222*Yu22)) + 0.6*(md200 + md211 + md222 - mDx200 -
      mDx211 - mDx222 + mDxbar200 + mDxbar211 + mDxbar222 + me200 + me211 + me222
      - mH1I200 - mH1I211 + mH2I200 + mH2I211 - mHd2 - mHp2 + mHpbar2 + mHu2 -
      ml200 - ml211 - ml222 + mq200 + mq211 + mq222 - 2*(mu200 + mu211 + mu222))*
      Sqr(g1) + 2*QH2p*(3*(md200 + md211 + md222)*Qdp + 3*(mDxbar200 + mDxbar211 +
      mDxbar222)*QDxbarp + 3*(mDx200 + mDx211 + mDx222)*QDxp + (me200 + me211 +
      me222)*Qep + 2*(mH1I200 + mH1I211)*QH1p + 2*mHd2*QH1p + 2*(mH2I200 + mH2I211
      )*QH2p + 2*mHu2*QH2p + 2*mHpbar2*QHpbarp + 2*mHp2*QHpp + 2*(ml200 + ml211 +
      ml222)*QLp + 6*(mq200 + mq211 + mq222)*QQp + ms2*QSp + (msI200 + msI211)*QSp
      + 3*(mu200 + mu211 + mu222)*Qup)*Sqr(gN) + 2*mHd2*Sqr(Lambdax) + 2*mHu2*Sqr
      (Lambdax) + 2*ms2*Sqr(Lambdax) + 2*Sqr(ALambdax)*Sqr(Lambdax) - 1.2*Sqr(g1)*
      Sqr(MassB) - 6*Sqr(g2)*Sqr(MassWB) - 8*Sqr(gN)*Sqr(MassBp)*Sqr(QH2p) + 6*(
      Sqr(TYu00) + Sqr(TYu01) + Sqr(TYu02) + Sqr(TYu10) + Sqr(TYu11) + Sqr(TYu12)
      + Sqr(TYu20) + Sqr(TYu21) + Sqr(TYu22)) + 6*mHu2*(Sqr(Yu00) + Sqr(Yu01) +
      Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(
      Yu22))) + 4*ALambdax*Lambdax*(6*(Kappa00*TKappa00 + Kappa01*TKappa01 +
      Kappa02*TKappa02 + Kappa10*TKappa10 + Kappa11*TKappa11 + Kappa12*TKappa12 +
      Kappa20*TKappa20 + Kappa21*TKappa21 + Kappa22*TKappa22) + 4*(Lambda1200*
      TLambda1200 + Lambda1201*TLambda1201 + Lambda1210*TLambda1210 + Lambda1211*
      TLambda1211) + 6*(TYd00*Yd00 + TYd01*Yd01 + TYd02*Yd02 + TYd10*Yd10 + TYd11*
      Yd11 + TYd12*Yd12 + TYd20*Yd20 + TYd21*Yd21 + TYd22*Yd22) + 2*(TYe00*Ye00 +
      TYe01*Ye01 + TYe02*Ye02 + TYe10*Ye10 + TYe11*Ye11 + TYe12*Ye12 + TYe20*Ye20
      + TYe21*Ye21 + TYe22*Ye22) + 6*(TYu00*Yu00 + TYu01*Yu01 + TYu02*Yu02 + TYu10
      *Yu10 + TYu11*Yu11 + TYu12*Yu12 + TYu20*Yu20 + TYu21*Yu21 + TYu22*Yu22) +
      1.2*MassB*Sqr(g1) + 6*MassWB*Sqr(g2) + 24*ALambdax*Sqr(Lambdax) + 4*MassBp*
      Sqr(gN)*Sqr(QH1p) + 4*MassBp*Sqr(gN)*Sqr(QH2p) + 4*MassBp*Sqr(gN)*Sqr(QSp) +
      ALambdax*(-0.6*Sqr(g1) - 3*Sqr(g2) + 3*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(
      Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(
      Kappa21) + Sqr(Kappa22)) + 2*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(
      Lambda1210) + Sqr(Lambda1211)) + 12*Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH1p) - 2*
      Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QSp) + 3*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02
      ) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) +
      Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr
      (Ye20) + Sqr(Ye21) + Sqr(Ye22) + 3*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(
      Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) + 4*
      ALambdax*(6*Lambdax*(Kappa00*TKappa00 + Kappa01*TKappa01 + Kappa02*TKappa02
      + Kappa10*TKappa10 + Kappa11*TKappa11 + Kappa12*TKappa12 + Kappa20*TKappa20
      + Kappa21*TKappa21 + Kappa22*TKappa22) + 4*Lambdax*(Lambda1200*TLambda1200 +
      Lambda1201*TLambda1201 + Lambda1210*TLambda1210 + Lambda1211*TLambda1211) +
      6*Lambdax*(TYd00*Yd00 + TYd01*Yd01 + TYd02*Yd02 + TYd10*Yd10 + TYd11*Yd11 +
      TYd12*Yd12 + TYd20*Yd20 + TYd21*Yd21 + TYd22*Yd22) + 2*Lambdax*(TYe00*Ye00
      + TYe01*Ye01 + TYe02*Ye02 + TYe10*Ye10 + TYe11*Ye11 + TYe12*Ye12 + TYe20*
      Ye20 + TYe21*Ye21 + TYe22*Ye22) + 6*Lambdax*(TYu00*Yu00 + TYu01*Yu01 + TYu02
      *Yu02 + TYu10*Yu10 + TYu11*Yu11 + TYu12*Yu12 + TYu20*Yu20 + TYu21*Yu21 +
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

double lowE6SSM_tuning_calculator::deriv_dlead_log_one_loop_mHu2_dALambdax() const
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

   double deriv = 24*TYu00*Yu00*Sqr(Lambdax) + 24*TYu01*Yu01*Sqr(Lambdax) + 24*
      TYu02*Yu02*Sqr(Lambdax) + 24*TYu10*Yu10*Sqr(Lambdax) + 24*TYu11*Yu11*Sqr(
      Lambdax) + 24*TYu12*Yu12*Sqr(Lambdax) + 24*TYu20*Yu20*Sqr(Lambdax) + 24*
      TYu21*Yu21*Sqr(Lambdax) + 24*TYu22*Yu22*Sqr(Lambdax) + 4*ALambdax*Sqr(
      Lambdax)*(-0.6*Sqr(g1) + 4*QH1p*QH2p*Sqr(gN) + 2*Sqr(Lambdax)) + 8*ALambdax*
      Sqr(Lambdax)*(2*QH2p*QSp*Sqr(gN) + 2*Sqr(Lambdax)) + 4*ALambdax*Sqr(Lambdax)
      *(-0.6*Sqr(g1) - 3*Sqr(g2) + 3*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) +
      Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) +
      Sqr(Kappa22)) + 2*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr
      (Lambda1211)) + 12*Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QH2p)
      - 2*Sqr(gN)*Sqr(QSp) + 3*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) +
      Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(Ye00) + Sqr
      (Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(
      Ye21) + Sqr(Ye22) + 3*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(
      Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + 4*ALambdax*Sqr(
      Lambdax)*(0.6*Sqr(g1) + 2*Sqr(Lambdax) + 4*Sqr(gN)*Sqr(QH2p) + 6*(Sqr(Yu00)
      + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) +
      Sqr(Yu21) + Sqr(Yu22))) + 4*Lambdax*(6*Lambdax*(Kappa00*TKappa00 + Kappa01*
      TKappa01 + Kappa02*TKappa02 + Kappa10*TKappa10 + Kappa11*TKappa11 + Kappa12*
      TKappa12 + Kappa20*TKappa20 + Kappa21*TKappa21 + Kappa22*TKappa22) + 4*
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

double lowE6SSM_tuning_calculator::deriv_dlead_log_one_loop_mHu2_dAYu22() const
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

   double deriv = 48*TYu02*Yu22*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + 48*TYu12*
      Yu22*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + 24*TYu20*Yu22*(Yu00*Yu02 + Yu10*
      Yu12 + Yu20*Yu22) + 24*TYu21*Yu22*(Yu01*Yu02 + Yu11*Yu12 + Yu21*Yu22) + 24*
      Lambdax*TLambdax*Sqr(Yu22) + 72*TYu00*Yu00*Sqr(Yu22) + 72*TYu01*Yu01*Sqr(
      Yu22) + 72*TYu10*Yu10*Sqr(Yu22) + 72*TYu11*Yu11*Sqr(Yu22) + 12*TYu02*(4*Yu22
      *(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + 6*Yu02*Sqr(Yu22)) + 12*TYu12*(4*Yu22*
      (Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + 6*Yu12*Sqr(Yu22)) + 4*AYu22*Sqr(Yu22)*
      (0.6*Sqr(g1) + 12*QH2p*QQp*Sqr(gN) + 6*(Sqr(Yu02) + Sqr(Yu12) + Sqr(Yu22)))
      + 8*AYu22*Sqr(Yu22)*(-1.2*Sqr(g1) + 6*QH2p*Qup*Sqr(gN) + 6*(Sqr(Yu20) + Sqr(
      Yu21) + Sqr(Yu22))) + 12*AYu22*Sqr(Yu22)*(0.6*Sqr(g1) + 2*Sqr(Lambdax) + 4*
      Sqr(gN)*Sqr(QH2p) + 6*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(
      Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + 12*AYu22*Yu22*(6*
      Power(Yu22,3) - 0.8666666666666667*Yu22*Sqr(g1) - 3*Yu22*Sqr(g2) -
      5.333333333333333*Yu22*Sqr(g3) + Yu22*Sqr(Lambdax) - 2*Yu22*Sqr(gN)*Sqr(QH2p
      ) - 2*Yu22*Sqr(gN)*Sqr(QQp) - 2*Yu22*Sqr(gN)*Sqr(Qup) + Yu22*Sqr(Yd02) +
      Yu22*Sqr(Yd12) + Yu22*Sqr(Yd22) + 5*(Power(Yu22,3) + Yu22*Sqr(Yu02) + Yu22*
      Sqr(Yu12)) + 4*Yu22*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)) + 3*Yu22*(Sqr(Yu00)
      + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) +
      Sqr(Yu21) + Sqr(Yu22))) + 12*TYu20*(Yd00*Yd02*Yu22 + Yd10*Yd12*Yu22 + Yd20*
      Yd22*Yu22 + 6*Yu20*Sqr(Yu22) + 5*(Yu00*Yu02*Yu22 + Yu10*Yu12*Yu22 + Yu20*Sqr
      (Yu22))) + 12*TYu21*(Yd01*Yd02*Yu22 + Yd11*Yd12*Yu22 + Yd21*Yd22*Yu22 + 6*
      Yu21*Sqr(Yu22) + 5*(Yu01*Yu02*Yu22 + Yu11*Yu12*Yu22 + Yu21*Sqr(Yu22))) + 12*
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

double lowE6SSM_tuning_calculator::deriv_dlead_log_one_loop_mHu2_dmq222() const
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

   double deriv = 11.52*Power(g1,4) + 48*Yu02*Yu12*(Yu00*Yu10 + Yu01*Yu11 +
      Yu02*Yu12) + 48*Yu02*Yu22*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + 48*Yu12*Yu22
      *(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + 12*(Yu00*Yu02 + Yu10*Yu12 + Yu20*Yu22
      )*(Yd00*Yd02 + Yd10*Yd12 + Yd20*Yd22 + Yu00*Yu02 + Yu10*Yu12 + Yu20*Yu22) +
      12*(Yu01*Yu02 + Yu11*Yu12 + Yu21*Yu22)*(Yd01*Yd02 + Yd11*Yd12 + Yd21*Yd22 +
      Yu01*Yu02 + Yu11*Yu12 + Yu21*Yu22) + 3*(0.6*Sqr(g1) + 6*QDxbarp*QH2p*Sqr(gN)
      )*(0.4*Sqr(g1) + 12*QDxbarp*QQp*Sqr(gN)) + 3*(-0.6*Sqr(g1) + 6*QDxp*QH2p*Sqr
      (gN))*(-0.4*Sqr(g1) + 12*QDxp*QQp*Sqr(gN)) + 3*(0.6*Sqr(g1) + 2*Qep*QH2p*Sqr
      (gN))*(1.2*Sqr(g1) + 12*Qep*QQp*Sqr(gN)) + 2*(-0.6*Sqr(g1) + 4*QH1p*QH2p*Sqr
      (gN))*(-0.6*Sqr(g1) + 12*QH1p*QQp*Sqr(gN)) + (0.6*Sqr(g1) + 4*QH2p*QHpbarp*
      Sqr(gN))*(0.6*Sqr(g1) + 12*QHpbarp*QQp*Sqr(gN)) + (-0.6*Sqr(g1) + 4*QH2p*
      QHpp*Sqr(gN))*(-0.6*Sqr(g1) + 12*QHpp*QQp*Sqr(gN)) + 3*(-0.6*Sqr(g1) + 4*
      QH2p*QLp*Sqr(gN))*(-0.6*Sqr(g1) + 12*QLp*QQp*Sqr(gN)) + 12*QQp*QSp*Sqr(gN)*(
      2*QH2p*QSp*Sqr(gN) + 2*Sqr(Lambdax)) + 2*(0.6*Sqr(g1) + 12*QH2p*QQp*Sqr(gN))
      *(0.6*Sqr(g1) + 4*Sqr(gN)*Sqr(QH2p)) + 48*Power(gN,4)*QH2p*QQp*Sqr(QSp) + 24
      *Power(gN,4)*QH2p*QQp*(9*Sqr(Qdp) + 9*Sqr(QDxbarp) + 9*Sqr(QDxp) + 3*Sqr(Qep
      ) + 6*Sqr(QH1p) + 6*Sqr(QH2p) + 2*Sqr(QHpbarp) + 2*Sqr(QHpp) + 6*Sqr(QLp) +
      18*Sqr(QQp) + 3*Sqr(QSp) + 9*Sqr(Qup)) + (0.6*Sqr(g1) + 6*Qdp*QH2p*Sqr(gN))*
      (0.4*Sqr(g1) + 12*Qdp*QQp*Sqr(gN) + 4*Sqr(Yd02)) + (0.6*Sqr(g1) + 6*Qdp*QH2p
      *Sqr(gN))*(0.4*Sqr(g1) + 12*Qdp*QQp*Sqr(gN) + 4*Sqr(Yd12)) + (0.6*Sqr(g1) +
      6*Qdp*QH2p*Sqr(gN))*(0.4*Sqr(g1) + 12*Qdp*QQp*Sqr(gN) + 4*Sqr(Yd22)) + (-0.6
      *Sqr(g1) + 4*QH1p*QH2p*Sqr(gN) + 2*Sqr(Lambdax))*(-0.6*Sqr(g1) + 12*QH1p*QQp
      *Sqr(gN) + 6*(Sqr(Yd02) + Sqr(Yd12) + Sqr(Yd22))) + (-0.8*Sqr(g1) + 12*QQp*
      Qup*Sqr(gN) + 4*Sqr(Yu02))*(-1.2*Sqr(g1) + 6*QH2p*Qup*Sqr(gN) + 6*(Sqr(Yu00)
      + Sqr(Yu01) + Sqr(Yu02))) + (-0.8*Sqr(g1) + 12*QQp*Qup*Sqr(gN) + 4*Sqr(Yu12
      ))*(-1.2*Sqr(g1) + 6*QH2p*Qup*Sqr(gN) + 6*(Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12)
      )) + (0.2*Sqr(g1) + 12*Sqr(gN)*Sqr(QQp))*(0.6*Sqr(g1) + 12*QH2p*QQp*Sqr(gN)
      + 6*(Sqr(Yu00) + Sqr(Yu10) + Sqr(Yu20))) + (0.2*Sqr(g1) + 12*Sqr(gN)*Sqr(QQp
      ))*(0.6*Sqr(g1) + 12*QH2p*QQp*Sqr(gN) + 6*(Sqr(Yu01) + Sqr(Yu11) + Sqr(Yu21)
      )) + (0.2*Sqr(g1) + 12*Sqr(gN)*Sqr(QQp) + 2*Sqr(Yd02) + 2*Sqr(Yd12) + 2*Sqr(
      Yd22) + 2*Sqr(Yu02) + 2*Sqr(Yu12) + 2*Sqr(Yu22))*(0.6*Sqr(g1) + 12*QH2p*QQp*
      Sqr(gN) + 6*(Sqr(Yu02) + Sqr(Yu12) + Sqr(Yu22))) + (-0.8*Sqr(g1) + 12*QQp*
      Qup*Sqr(gN) + 4*Sqr(Yu22))*(-1.2*Sqr(g1) + 6*QH2p*Qup*Sqr(gN) + 6*(Sqr(Yu20)
      + Sqr(Yu21) + Sqr(Yu22))) + (0.6*Sqr(g1) + 12*QH2p*QQp*Sqr(gN) + 6*(Sqr(
      Yu02) + Sqr(Yu12) + Sqr(Yu22)))*(0.6*Sqr(g1) + 2*Sqr(Lambdax) + 4*Sqr(gN)*
      Sqr(QH2p) + 6*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) +
      Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + 12*Yu02*(Yd02*(Yd00*Yu00 +
      Yd01*Yu01 + Yd02*Yu02) + Yd12*(Yd10*Yu00 + Yd11*Yu01 + Yd12*Yu02) + Yd22*(
      Yd20*Yu00 + Yd21*Yu01 + Yd22*Yu02) + 3*(Yu12*(Yu00*Yu10 + Yu01*Yu11 + Yu02*
      Yu12) + Yu22*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu02*(Sqr(Yu00) + Sqr(
      Yu01) + Sqr(Yu02))) + Yu02*(-0.8666666666666667*Sqr(g1) - 3*Sqr(g2) -
      5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*
      Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(
      Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) + 12*
      Yu12*(Yd02*(Yd00*Yu10 + Yd01*Yu11 + Yd02*Yu12) + Yd12*(Yd10*Yu10 + Yd11*Yu11
      + Yd12*Yu12) + Yd22*(Yd20*Yu10 + Yd21*Yu11 + Yd22*Yu12) + 3*(Yu02*(Yu00*
      Yu10 + Yu01*Yu11 + Yu02*Yu12) + Yu22*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) +
      Yu12*(Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12))) + Yu12*(-0.8666666666666667*Sqr(g1
      ) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(
      QH2p) - 2*Sqr(gN)*Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3*(Sqr(Yu00) + Sqr(Yu01) +
      Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr
      (Yu22)))) + 12*Yu22*(Yd02*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22) + Yd12*(Yd10*
      Yu20 + Yd11*Yu21 + Yd12*Yu22) + Yd22*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22) + 3
      *(Yu02*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu12*(Yu10*Yu20 + Yu11*Yu21 +
      Yu12*Yu22) + Yu22*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + Yu22*(
      -0.8666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(
      Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3
      *(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) +
      Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dlead_log_one_loop_mHu2_dmHd2() const
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

   double deriv = -11.52*Power(g1,4) + 24*(Yd00*Yd01 + Yd10*Yd11 + Yd20*Yd21)*(
      Yu00*Yu01 + Yu10*Yu11 + Yu20*Yu21) + 24*(Yd00*Yd02 + Yd10*Yd12 + Yd20*Yd22)*
      (Yu00*Yu02 + Yu10*Yu12 + Yu20*Yu22) + 24*(Yd01*Yd02 + Yd11*Yd12 + Yd21*Yd22)
      *(Yu01*Yu02 + Yu11*Yu12 + Yu21*Yu22) + 3*(-0.4*Sqr(g1) + 4*QDxbarp*QH1p*Sqr(
      gN))*(0.6*Sqr(g1) + 6*QDxbarp*QH2p*Sqr(gN)) + 3*(0.4*Sqr(g1) + 4*QDxp*QH1p*
      Sqr(gN))*(-0.6*Sqr(g1) + 6*QDxp*QH2p*Sqr(gN)) + (-0.6*Sqr(g1) + 4*QH1p*
      QHpbarp*Sqr(gN))*(0.6*Sqr(g1) + 4*QH2p*QHpbarp*Sqr(gN)) + (0.6*Sqr(g1) + 4*
      QH1p*QHpp*Sqr(gN))*(-0.6*Sqr(g1) + 4*QH2p*QHpp*Sqr(gN)) + (2*QH2p*QSp*Sqr(gN
      ) + 2*Sqr(Lambdax))*(4*QH1p*QSp*Sqr(gN) + 4*Sqr(Lambdax)) + 2*(-0.6*Sqr(g1)
      + 4*QH1p*QH2p*Sqr(gN))*(0.6*Sqr(g1) + 4*Sqr(gN)*Sqr(QH1p)) + 2*(-0.6*Sqr(g1)
      + 4*QH1p*QH2p*Sqr(gN))*(0.6*Sqr(g1) + 4*Sqr(gN)*Sqr(QH2p)) + 16*Power(gN,4)
      *QH1p*QH2p*Sqr(QSp) + 8*Power(gN,4)*QH1p*QH2p*(9*Sqr(Qdp) + 9*Sqr(QDxbarp) +
      9*Sqr(QDxp) + 3*Sqr(Qep) + 6*Sqr(QH1p) + 6*Sqr(QH2p) + 2*Sqr(QHpbarp) + 2*
      Sqr(QHpp) + 6*Sqr(QLp) + 18*Sqr(QQp) + 3*Sqr(QSp) + 9*Sqr(Qup)) + (0.6*Sqr(
      g1) + 6*Qdp*QH2p*Sqr(gN))*(-0.4*Sqr(g1) + 4*Qdp*QH1p*Sqr(gN) + 4*(Sqr(Yd00)
      + Sqr(Yd01) + Sqr(Yd02))) + (0.6*Sqr(g1) + 6*Qdp*QH2p*Sqr(gN))*(-0.4*Sqr(g1)
      + 4*Qdp*QH1p*Sqr(gN) + 4*(Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12))) + (0.6*Sqr(g1
      ) + 6*Qdp*QH2p*Sqr(gN))*(-0.4*Sqr(g1) + 4*Qdp*QH1p*Sqr(gN) + 4*(Sqr(Yd20) +
      Sqr(Yd21) + Sqr(Yd22))) + (0.6*Sqr(g1) + 2*Qep*QH2p*Sqr(gN))*(-1.2*Sqr(g1) +
      4*Qep*QH1p*Sqr(gN) + 4*(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02))) + (0.6*Sqr(g1)
      + 2*Qep*QH2p*Sqr(gN))*(-1.2*Sqr(g1) + 4*Qep*QH1p*Sqr(gN) + 4*(Sqr(Ye10) +
      Sqr(Ye11) + Sqr(Ye12))) + (-0.6*Sqr(g1) + 4*QH2p*QLp*Sqr(gN))*(0.6*Sqr(g1) +
      4*QH1p*QLp*Sqr(gN) + 2*(Sqr(Ye00) + Sqr(Ye10) + Sqr(Ye20))) + (-0.6*Sqr(g1)
      + 4*QH2p*QLp*Sqr(gN))*(0.6*Sqr(g1) + 4*QH1p*QLp*Sqr(gN) + 2*(Sqr(Ye01) +
      Sqr(Ye11) + Sqr(Ye21))) + (-0.6*Sqr(g1) + 4*QH2p*QLp*Sqr(gN))*(0.6*Sqr(g1) +
      4*QH1p*QLp*Sqr(gN) + 2*(Sqr(Ye02) + Sqr(Ye12) + Sqr(Ye22))) + (0.6*Sqr(g1)
      + 2*Qep*QH2p*Sqr(gN))*(-1.2*Sqr(g1) + 4*Qep*QH1p*Sqr(gN) + 4*(Sqr(Ye20) +
      Sqr(Ye21) + Sqr(Ye22))) + (-0.6*Sqr(g1) + 4*QH1p*QH2p*Sqr(gN) + 2*Sqr(
      Lambdax))*(0.6*Sqr(g1) + 2*Sqr(Lambdax) + 4*Sqr(gN)*Sqr(QH1p) + 6*(Sqr(Yd00)
      + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) +
      Sqr(Yd21) + Sqr(Yd22)) + 2*(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) +
      Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22))) + (0.8*Sqr(g1) +
      4*QH1p*Qup*Sqr(gN))*(-1.2*Sqr(g1) + 6*QH2p*Qup*Sqr(gN) + 6*(Sqr(Yu00) + Sqr
      (Yu01) + Sqr(Yu02))) + (0.8*Sqr(g1) + 4*QH1p*Qup*Sqr(gN))*(-1.2*Sqr(g1) + 6*
      QH2p*Qup*Sqr(gN) + 6*(Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12))) + (-0.2*Sqr(g1) +
      4*QH1p*QQp*Sqr(gN) + 2*(Sqr(Yd00) + Sqr(Yd10) + Sqr(Yd20)))*(0.6*Sqr(g1) +
      12*QH2p*QQp*Sqr(gN) + 6*(Sqr(Yu00) + Sqr(Yu10) + Sqr(Yu20))) + (-0.2*Sqr(g1)
      + 4*QH1p*QQp*Sqr(gN) + 2*(Sqr(Yd01) + Sqr(Yd11) + Sqr(Yd21)))*(0.6*Sqr(g1)
      + 12*QH2p*QQp*Sqr(gN) + 6*(Sqr(Yu01) + Sqr(Yu11) + Sqr(Yu21))) + (-0.2*Sqr(
      g1) + 4*QH1p*QQp*Sqr(gN) + 2*(Sqr(Yd02) + Sqr(Yd12) + Sqr(Yd22)))*(0.6*Sqr(
      g1) + 12*QH2p*QQp*Sqr(gN) + 6*(Sqr(Yu02) + Sqr(Yu12) + Sqr(Yu22))) + (0.8*
      Sqr(g1) + 4*QH1p*Qup*Sqr(gN))*(-1.2*Sqr(g1) + 6*QH2p*Qup*Sqr(gN) + 6*(Sqr(
      Yu20) + Sqr(Yu21) + Sqr(Yu22))) + (-0.6*Sqr(g1) + 4*QH1p*QH2p*Sqr(gN) + 2*
      Sqr(Lambdax))*(0.6*Sqr(g1) + 2*Sqr(Lambdax) + 4*Sqr(gN)*Sqr(QH2p) + 6*(Sqr(
      Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20
      ) + Sqr(Yu21) + Sqr(Yu22))) + 4*Lambdax*(4*Power(Lambdax,3) - 0.6*Lambdax*
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

double lowE6SSM_tuning_calculator::deriv_dlead_log_one_loop_mHu2_dmHu2() const
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

   double deriv = 11.52*Power(g1,4) + 3*(0.4*Sqr(g1) + 4*Qdp*QH2p*Sqr(gN))*(0.6
      *Sqr(g1) + 6*Qdp*QH2p*Sqr(gN)) + 3*(0.4*Sqr(g1) + 4*QDxbarp*QH2p*Sqr(gN))*(
      0.6*Sqr(g1) + 6*QDxbarp*QH2p*Sqr(gN)) + 3*(-0.4*Sqr(g1) + 4*QDxp*QH2p*Sqr(gN
      ))*(-0.6*Sqr(g1) + 6*QDxp*QH2p*Sqr(gN)) + 3*(0.6*Sqr(g1) + 2*Qep*QH2p*Sqr(gN
      ))*(1.2*Sqr(g1) + 4*Qep*QH2p*Sqr(gN)) + (2*QH2p*QSp*Sqr(gN) + 2*Sqr(Lambdax)
      )*(4*QH2p*QSp*Sqr(gN) + 4*Sqr(Lambdax)) + 2*Sqr(-0.6*Power(g1,2) + 4*Power(
      gN,2)*QH1p*QH2p) + Sqr(-0.6*Power(g1,2) + 2*Power(Lambdax,2) + 4*Power(gN,2)
      *QH1p*QH2p) + 2*Sqr(0.6*Power(g1,2) + 4*Power(gN,2)*Power(QH2p,2)) + Sqr(0.6
      *Power(g1,2) + 4*Power(gN,2)*QH2p*QHpbarp) + Sqr(-0.6*Power(g1,2) + 4*Power(
      gN,2)*QH2p*QHpp) + 3*Sqr(-0.6*Power(g1,2) + 4*Power(gN,2)*QH2p*QLp) + 16*
      Power(gN,4)*Sqr(QH2p)*Sqr(QSp) + 8*Power(gN,4)*Sqr(QH2p)*(9*Sqr(Qdp) + 9*Sqr
      (QDxbarp) + 9*Sqr(QDxp) + 3*Sqr(Qep) + 6*Sqr(QH1p) + 6*Sqr(QH2p) + 2*Sqr(
      QHpbarp) + 2*Sqr(QHpp) + 6*Sqr(QLp) + 18*Sqr(QQp) + 3*Sqr(QSp) + 9*Sqr(Qup))
      + (-0.8*Sqr(g1) + 4*QH2p*Qup*Sqr(gN) + 4*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02)
      ))*(-1.2*Sqr(g1) + 6*QH2p*Qup*Sqr(gN) + 6*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02)
      )) + (-0.8*Sqr(g1) + 4*QH2p*Qup*Sqr(gN) + 4*(Sqr(Yu10) + Sqr(Yu11) + Sqr(
      Yu12)))*(-1.2*Sqr(g1) + 6*QH2p*Qup*Sqr(gN) + 6*(Sqr(Yu10) + Sqr(Yu11) + Sqr(
      Yu12))) + 48*Sqr(Yu00*Yu10 + Yu01*Yu11 + Yu02*Yu12) + (0.2*Sqr(g1) + 4*QH2p*
      QQp*Sqr(gN) + 2*(Sqr(Yu00) + Sqr(Yu10) + Sqr(Yu20)))*(0.6*Sqr(g1) + 12*QH2p*
      QQp*Sqr(gN) + 6*(Sqr(Yu00) + Sqr(Yu10) + Sqr(Yu20))) + (0.2*Sqr(g1) + 4*QH2p
      *QQp*Sqr(gN) + 2*(Sqr(Yu01) + Sqr(Yu11) + Sqr(Yu21)))*(0.6*Sqr(g1) + 12*QH2p
      *QQp*Sqr(gN) + 6*(Sqr(Yu01) + Sqr(Yu11) + Sqr(Yu21))) + 24*Sqr(Yu00*Yu01 +
      Yu10*Yu11 + Yu20*Yu21) + (0.2*Sqr(g1) + 4*QH2p*QQp*Sqr(gN) + 2*(Sqr(Yu02) +
      Sqr(Yu12) + Sqr(Yu22)))*(0.6*Sqr(g1) + 12*QH2p*QQp*Sqr(gN) + 6*(Sqr(Yu02) +
      Sqr(Yu12) + Sqr(Yu22))) + (-0.8*Sqr(g1) + 4*QH2p*Qup*Sqr(gN) + 4*(Sqr(Yu20)
      + Sqr(Yu21) + Sqr(Yu22)))*(-1.2*Sqr(g1) + 6*QH2p*Qup*Sqr(gN) + 6*(Sqr(Yu20)
      + Sqr(Yu21) + Sqr(Yu22))) + 4*Lambdax*(4*Power(Lambdax,3) - 0.6*Lambdax*Sqr(
      g1) - 3*Lambdax*Sqr(g2) + 3*Lambdax*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(
      Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(
      Kappa21) + Sqr(Kappa22)) + 2*Lambdax*(Sqr(Lambda1200) + Sqr(Lambda1201) +
      Sqr(Lambda1210) + Sqr(Lambda1211)) - 2*Lambdax*Sqr(gN)*Sqr(QH1p) - 2*Lambdax
      *Sqr(gN)*Sqr(QH2p) - 2*Lambdax*Sqr(gN)*Sqr(QSp) + 3*Lambdax*(Sqr(Yd00) + Sqr
      (Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(
      Yd21) + Sqr(Yd22)) + Lambdax*(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10)
      + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22)) + 3*Lambdax*(
      Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(
      Yu20) + Sqr(Yu21) + Sqr(Yu22))) + 12*Yu00*(Yd00*(Yd00*Yu00 + Yd01*Yu01 +
      Yd02*Yu02) + Yd10*(Yd10*Yu00 + Yd11*Yu01 + Yd12*Yu02) + Yd20*(Yd20*Yu00 +
      Yd21*Yu01 + Yd22*Yu02) + 3*(Yu10*(Yu00*Yu10 + Yu01*Yu11 + Yu02*Yu12) + Yu20*
      (Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu00*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02
      ))) + Yu00*(-0.8666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(
      g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QQp) - 2*Sqr(gN)*
      Sqr(Qup) + 3*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) +
      Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) + 12*Yu01*(Yd01*(Yd00*Yu00
      + Yd01*Yu01 + Yd02*Yu02) + Yd11*(Yd10*Yu00 + Yd11*Yu01 + Yd12*Yu02) + Yd21*(
      Yd20*Yu00 + Yd21*Yu01 + Yd22*Yu02) + 3*(Yu11*(Yu00*Yu10 + Yu01*Yu11 + Yu02*
      Yu12) + Yu21*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu01*(Sqr(Yu00) + Sqr(
      Yu01) + Sqr(Yu02))) + Yu01*(-0.8666666666666667*Sqr(g1) - 3*Sqr(g2) -
      5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*
      Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(
      Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) + 12*
      Yu02*(Yd02*(Yd00*Yu00 + Yd01*Yu01 + Yd02*Yu02) + Yd12*(Yd10*Yu00 + Yd11*Yu01
      + Yd12*Yu02) + Yd22*(Yd20*Yu00 + Yd21*Yu01 + Yd22*Yu02) + 3*(Yu12*(Yu00*
      Yu10 + Yu01*Yu11 + Yu02*Yu12) + Yu22*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) +
      Yu02*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02))) + Yu02*(-0.8666666666666667*Sqr(g1
      ) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(
      QH2p) - 2*Sqr(gN)*Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3*(Sqr(Yu00) + Sqr(Yu01) +
      Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr
      (Yu22)))) + 12*Yu10*(Yd00*(Yd00*Yu10 + Yd01*Yu11 + Yd02*Yu12) + Yd10*(Yd10*
      Yu10 + Yd11*Yu11 + Yd12*Yu12) + Yd20*(Yd20*Yu10 + Yd21*Yu11 + Yd22*Yu12) + 3
      *(Yu00*(Yu00*Yu10 + Yu01*Yu11 + Yu02*Yu12) + Yu20*(Yu10*Yu20 + Yu11*Yu21 +
      Yu12*Yu22) + Yu10*(Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12))) + Yu10*(
      -0.8666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(
      Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3
      *(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) +
      Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) + 12*Yu11*(Yd01*(Yd00*Yu10 + Yd01*Yu11
      + Yd02*Yu12) + Yd11*(Yd10*Yu10 + Yd11*Yu11 + Yd12*Yu12) + Yd21*(Yd20*Yu10 +
      Yd21*Yu11 + Yd22*Yu12) + 3*(Yu01*(Yu00*Yu10 + Yu01*Yu11 + Yu02*Yu12) + Yu21*
      (Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + Yu11*(Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12
      ))) + Yu11*(-0.8666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(
      g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QQp) - 2*Sqr(gN)*
      Sqr(Qup) + 3*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) +
      Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) + 12*Yu12*(Yd02*(Yd00*Yu10
      + Yd01*Yu11 + Yd02*Yu12) + Yd12*(Yd10*Yu10 + Yd11*Yu11 + Yd12*Yu12) + Yd22*(
      Yd20*Yu10 + Yd21*Yu11 + Yd22*Yu12) + 3*(Yu02*(Yu00*Yu10 + Yu01*Yu11 + Yu02*
      Yu12) + Yu22*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + Yu12*(Sqr(Yu10) + Sqr(
      Yu11) + Sqr(Yu12))) + Yu12*(-0.8666666666666667*Sqr(g1) - 3*Sqr(g2) -
      5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*
      Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(
      Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) + 12*
      Yu20*(Yd00*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22) + Yd10*(Yd10*Yu20 + Yd11*Yu21
      + Yd12*Yu22) + Yd20*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22) + 3*(Yu00*(Yu00*
      Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu10*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) +
      Yu20*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + Yu20*(-0.8666666666666667*Sqr(g1
      ) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(
      QH2p) - 2*Sqr(gN)*Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3*(Sqr(Yu00) + Sqr(Yu01) +
      Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr
      (Yu22)))) + 12*Yu21*(Yd01*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22) + Yd11*(Yd10*
      Yu20 + Yd11*Yu21 + Yd12*Yu22) + Yd21*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22) + 3
      *(Yu01*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu11*(Yu10*Yu20 + Yu11*Yu21 +
      Yu12*Yu22) + Yu21*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + Yu21*(
      -0.8666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(
      Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3
      *(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) +
      Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) + 12*Yu22*(Yd02*(Yd00*Yu20 + Yd01*Yu21
      + Yd02*Yu22) + Yd12*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22) + Yd22*(Yd20*Yu20 +
      Yd21*Yu21 + Yd22*Yu22) + 3*(Yu02*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu12*
      (Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + Yu22*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22
      ))) + Yu22*(-0.8666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(
      g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QQp) - 2*Sqr(gN)*
      Sqr(Qup) + 3*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) +
      Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) + 48*Sqr(Yu00*Yu20 + Yu01*
      Yu21 + Yu02*Yu22) + 48*Sqr(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + 24*Sqr(Yu00*
      Yu02 + Yu10*Yu12 + Yu20*Yu22) + 24*Sqr(Yu01*Yu02 + Yu11*Yu12 + Yu21*Yu22) +
      Sqr(0.6*Power(g1,2) + 2*Power(Lambdax,2) + 4*Power(gN,2)*Power(QH2p,2) + 6*(
      Power(Yu00,2) + Power(Yu01,2) + Power(Yu02,2) + Power(Yu10,2) + Power(Yu11,2
      ) + Power(Yu12,2) + Power(Yu20,2) + Power(Yu21,2) + Power(Yu22,2)));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dlead_log_one_loop_mHu2_dmu222() const
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

   double deriv = -23.04*Power(g1,4) + 24*Yu20*Yu21*(Yu00*Yu01 + Yu10*Yu11 +
      Yu20*Yu21) + 24*Yu20*Yu22*(Yu00*Yu02 + Yu10*Yu12 + Yu20*Yu22) + 24*Yu21*Yu22
      *(Yu01*Yu02 + Yu11*Yu12 + Yu21*Yu22) + 3*(0.6*Sqr(g1) + 6*Qdp*QH2p*Sqr(gN))*
      (-0.8*Sqr(g1) + 6*Qdp*Qup*Sqr(gN)) + 3*(0.6*Sqr(g1) + 6*QDxbarp*QH2p*Sqr(gN)
      )*(-0.8*Sqr(g1) + 6*QDxbarp*Qup*Sqr(gN)) + 3*(-0.6*Sqr(g1) + 6*QDxp*QH2p*Sqr
      (gN))*(0.8*Sqr(g1) + 6*QDxp*Qup*Sqr(gN)) + 3*(0.6*Sqr(g1) + 2*Qep*QH2p*Sqr(
      gN))*(-2.4*Sqr(g1) + 6*Qep*Qup*Sqr(gN)) + 2*(-0.6*Sqr(g1) + 4*QH1p*QH2p*Sqr(
      gN))*(1.2*Sqr(g1) + 6*QH1p*Qup*Sqr(gN)) + (0.6*Sqr(g1) + 4*QH2p*QHpbarp*Sqr(
      gN))*(-1.2*Sqr(g1) + 6*QHpbarp*Qup*Sqr(gN)) + (-0.6*Sqr(g1) + 4*QH2p*QHpp*
      Sqr(gN))*(1.2*Sqr(g1) + 6*QHpp*Qup*Sqr(gN)) + 3*(-0.6*Sqr(g1) + 4*QH2p*QLp*
      Sqr(gN))*(1.2*Sqr(g1) + 6*QLp*Qup*Sqr(gN)) + (1.2*Sqr(g1) + 6*QH1p*Qup*Sqr(
      gN))*(-0.6*Sqr(g1) + 4*QH1p*QH2p*Sqr(gN) + 2*Sqr(Lambdax)) + 6*QSp*Qup*Sqr(
      gN)*(2*QH2p*QSp*Sqr(gN) + 2*Sqr(Lambdax)) + 2*(-1.2*Sqr(g1) + 6*QH2p*Qup*Sqr
      (gN))*(0.6*Sqr(g1) + 4*Sqr(gN)*Sqr(QH2p)) + 24*Power(gN,4)*QH2p*Qup*Sqr(QSp)
      + 12*Power(gN,4)*QH2p*Qup*(9*Sqr(Qdp) + 9*Sqr(QDxbarp) + 9*Sqr(QDxp) + 3*
      Sqr(Qep) + 6*Sqr(QH1p) + 6*Sqr(QH2p) + 2*Sqr(QHpbarp) + 2*Sqr(QHpp) + 6*Sqr(
      QLp) + 18*Sqr(QQp) + 3*Sqr(QSp) + 9*Sqr(Qup)) + (1.6*Sqr(g1) + 6*Sqr(gN)*Sqr
      (Qup))*(-1.2*Sqr(g1) + 6*QH2p*Qup*Sqr(gN) + 6*(Sqr(Yu00) + Sqr(Yu01) + Sqr(
      Yu02))) + (1.6*Sqr(g1) + 6*Sqr(gN)*Sqr(Qup))*(-1.2*Sqr(g1) + 6*QH2p*Qup*Sqr(
      gN) + 6*(Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12))) + (-0.4*Sqr(g1) + 6*QQp*Qup*Sqr
      (gN) + 2*Sqr(Yu20))*(0.6*Sqr(g1) + 12*QH2p*QQp*Sqr(gN) + 6*(Sqr(Yu00) + Sqr(
      Yu10) + Sqr(Yu20))) + (-0.4*Sqr(g1) + 6*QQp*Qup*Sqr(gN) + 2*Sqr(Yu21))*(0.6*
      Sqr(g1) + 12*QH2p*QQp*Sqr(gN) + 6*(Sqr(Yu01) + Sqr(Yu11) + Sqr(Yu21))) + (
      -0.4*Sqr(g1) + 6*QQp*Qup*Sqr(gN) + 2*Sqr(Yu22))*(0.6*Sqr(g1) + 12*QH2p*QQp*
      Sqr(gN) + 6*(Sqr(Yu02) + Sqr(Yu12) + Sqr(Yu22))) + (1.6*Sqr(g1) + 6*Sqr(gN)*
      Sqr(Qup) + 4*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))*(-1.2*Sqr(g1) + 6*QH2p*Qup
      *Sqr(gN) + 6*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + (-1.2*Sqr(g1) + 6*QH2p*
      Qup*Sqr(gN) + 6*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))*(0.6*Sqr(g1) + 2*Sqr(
      Lambdax) + 4*Sqr(gN)*Sqr(QH2p) + 6*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(
      Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + 12*
      Yu20*(Yd00*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22) + Yd10*(Yd10*Yu20 + Yd11*Yu21
      + Yd12*Yu22) + Yd20*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22) + 3*(Yu00*(Yu00*
      Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu10*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) +
      Yu20*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + Yu20*(-0.8666666666666667*Sqr(g1
      ) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(
      QH2p) - 2*Sqr(gN)*Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3*(Sqr(Yu00) + Sqr(Yu01) +
      Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr
      (Yu22)))) + 12*Yu21*(Yd01*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22) + Yd11*(Yd10*
      Yu20 + Yd11*Yu21 + Yd12*Yu22) + Yd21*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22) + 3
      *(Yu01*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu11*(Yu10*Yu20 + Yu11*Yu21 +
      Yu12*Yu22) + Yu21*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + Yu21*(
      -0.8666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(
      Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3
      *(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) +
      Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) + 12*Yu22*(Yd02*(Yd00*Yu20 + Yd01*Yu21
      + Yd02*Yu22) + Yd12*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22) + Yd22*(Yd20*Yu20 +
      Yd21*Yu21 + Yd22*Yu22) + 3*(Yu02*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu12*
      (Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + Yu22*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22
      ))) + Yu22*(-0.8666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(
      g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QQp) - 2*Sqr(gN)*
      Sqr(Qup) + 3*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) +
      Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) + 24*Sqr(Yu00*Yu20 + Yu01*
      Yu21 + Yu02*Yu22) + 24*Sqr(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22);


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dlead_log_one_loop_mHu2_dms2() const
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

   double deriv = 8*Power(gN,4)*QH2p*Power(QSp,3) + 6*Qdp*QSp*Sqr(gN)*(0.6*Sqr(
      g1) + 6*Qdp*QH2p*Sqr(gN)) + 6*Qep*QSp*Sqr(gN)*(0.6*Sqr(g1) + 2*Qep*QH2p*Sqr(
      gN)) + 2*QHpbarp*QSp*Sqr(gN)*(0.6*Sqr(g1) + 4*QH2p*QHpbarp*Sqr(gN)) + 2*QHpp
      *QSp*Sqr(gN)*(-0.6*Sqr(g1) + 4*QH2p*QHpp*Sqr(gN)) + 6*QLp*QSp*Sqr(gN)*(-0.6*
      Sqr(g1) + 4*QH2p*QLp*Sqr(gN)) + (-0.6*Sqr(g1) + 6*QDxp*QH2p*Sqr(gN))*(2*QDxp
      *QSp*Sqr(gN) + 2*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02))) + (-0.6*Sqr(
      g1) + 6*QDxp*QH2p*Sqr(gN))*(2*QDxp*QSp*Sqr(gN) + 2*(Sqr(Kappa10) + Sqr(
      Kappa11) + Sqr(Kappa12))) + (0.6*Sqr(g1) + 6*QDxbarp*QH2p*Sqr(gN))*(2*
      QDxbarp*QSp*Sqr(gN) + 2*(Sqr(Kappa00) + Sqr(Kappa10) + Sqr(Kappa20))) + (0.6
      *Sqr(g1) + 6*QDxbarp*QH2p*Sqr(gN))*(2*QDxbarp*QSp*Sqr(gN) + 2*(Sqr(Kappa01)
      + Sqr(Kappa11) + Sqr(Kappa21))) + (0.6*Sqr(g1) + 6*QDxbarp*QH2p*Sqr(gN))*(2*
      QDxbarp*QSp*Sqr(gN) + 2*(Sqr(Kappa02) + Sqr(Kappa12) + Sqr(Kappa22))) + (
      -0.6*Sqr(g1) + 6*QDxp*QH2p*Sqr(gN))*(2*QDxp*QSp*Sqr(gN) + 2*(Sqr(Kappa20) +
      Sqr(Kappa21) + Sqr(Kappa22))) + (-0.6*Sqr(g1) + 4*QH1p*QH2p*Sqr(gN))*(2*QH1p
      *QSp*Sqr(gN) + 2*(Sqr(Lambda1200) + Sqr(Lambda1210))) + (-0.6*Sqr(g1) + 4*
      QH1p*QH2p*Sqr(gN))*(2*QH1p*QSp*Sqr(gN) + 2*(Sqr(Lambda1201) + Sqr(Lambda1211
      ))) + (-0.6*Sqr(g1) + 4*QH1p*QH2p*Sqr(gN) + 2*Sqr(Lambdax))*(2*QH1p*QSp*Sqr(
      gN) + 2*Sqr(Lambdax)) + (2*QH2p*QSp*Sqr(gN) + 2*(Sqr(Lambda1200) + Sqr(
      Lambda1201)))*(0.6*Sqr(g1) + 4*Sqr(gN)*Sqr(QH2p)) + (2*QH2p*QSp*Sqr(gN) + 2*
      (Sqr(Lambda1210) + Sqr(Lambda1211)))*(0.6*Sqr(g1) + 4*Sqr(gN)*Sqr(QH2p)) + (
      2*QH2p*QSp*Sqr(gN) + 2*Sqr(Lambdax))*(6*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(
      Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(
      Kappa21) + Sqr(Kappa22)) + 4*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(
      Lambda1210) + Sqr(Lambda1211)) + 4*Sqr(Lambdax) + 2*Sqr(gN)*Sqr(QSp)) + 4*
      Power(gN,4)*QH2p*QSp*(9*Sqr(Qdp) + 9*Sqr(QDxbarp) + 9*Sqr(QDxp) + 3*Sqr(Qep)
      + 6*Sqr(QH1p) + 6*Sqr(QH2p) + 2*Sqr(QHpbarp) + 2*Sqr(QHpp) + 6*Sqr(QLp) +
      18*Sqr(QQp) + 3*Sqr(QSp) + 9*Sqr(Qup)) + 2*QSp*Qup*Sqr(gN)*(-1.2*Sqr(g1) + 6
      *QH2p*Qup*Sqr(gN) + 6*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02))) + 2*QSp*Qup*Sqr(
      gN)*(-1.2*Sqr(g1) + 6*QH2p*Qup*Sqr(gN) + 6*(Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12
      ))) + 2*QQp*QSp*Sqr(gN)*(0.6*Sqr(g1) + 12*QH2p*QQp*Sqr(gN) + 6*(Sqr(Yu00) +
      Sqr(Yu10) + Sqr(Yu20))) + 2*QQp*QSp*Sqr(gN)*(0.6*Sqr(g1) + 12*QH2p*QQp*Sqr(
      gN) + 6*(Sqr(Yu01) + Sqr(Yu11) + Sqr(Yu21))) + 2*QQp*QSp*Sqr(gN)*(0.6*Sqr(g1
      ) + 12*QH2p*QQp*Sqr(gN) + 6*(Sqr(Yu02) + Sqr(Yu12) + Sqr(Yu22))) + 2*QSp*Qup
      *Sqr(gN)*(-1.2*Sqr(g1) + 6*QH2p*Qup*Sqr(gN) + 6*(Sqr(Yu20) + Sqr(Yu21) + Sqr
      (Yu22))) + (2*QH2p*QSp*Sqr(gN) + 2*Sqr(Lambdax))*(0.6*Sqr(g1) + 2*Sqr(
      Lambdax) + 4*Sqr(gN)*Sqr(QH2p) + 6*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(
      Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + 4*
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

double lowE6SSM_tuning_calculator::deriv_dlead_log_one_loop_mHu2_dMassB() const
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

   double deriv = -138.24*Power(g1,4)*MassB + 4.8*Lambdax*TLambdax*Sqr(g1) +
      20.8*TYu00*Yu00*Sqr(g1) + 20.8*TYu01*Yu01*Sqr(g1) + 20.8*TYu02*Yu02*Sqr(g1)
      + 20.8*TYu10*Yu10*Sqr(g1) + 20.8*TYu11*Yu11*Sqr(g1) + 20.8*TYu12*Yu12*Sqr(g1
      ) + 20.8*TYu20*Yu20*Sqr(g1) + 20.8*TYu21*Yu21*Sqr(g1) + 20.8*TYu22*Yu22*Sqr(
      g1) - 3.2*MassB*Sqr(g1)*(0.6*Sqr(g1) + 6*Qdp*QH2p*Sqr(gN)) - 3.2*MassB*Sqr(
      g1)*(0.6*Sqr(g1) + 6*QDxbarp*QH2p*Sqr(gN)) - 3.2*MassB*Sqr(g1)*(-0.6*Sqr(g1)
      + 6*QDxp*QH2p*Sqr(gN)) - 28.8*MassB*Sqr(g1)*(0.6*Sqr(g1) + 2*Qep*QH2p*Sqr(
      gN)) - 4.8*MassB*Sqr(g1)*(-0.6*Sqr(g1) + 4*QH1p*QH2p*Sqr(gN)) - 2.4*MassB*
      Sqr(g1)*(0.6*Sqr(g1) + 4*QH2p*QHpbarp*Sqr(gN)) - 2.4*MassB*Sqr(g1)*(-0.6*Sqr
      (g1) + 4*QH2p*QHpp*Sqr(gN)) - 7.2*MassB*Sqr(g1)*(-0.6*Sqr(g1) + 4*QH2p*QLp*
      Sqr(gN)) - 2.4*MassB*Sqr(g1)*(-0.6*Sqr(g1) + 4*QH1p*QH2p*Sqr(gN) + 2*Sqr(
      Lambdax)) - 4.8*MassB*Sqr(g1)*(0.6*Sqr(g1) + 4*Sqr(gN)*Sqr(QH2p)) -
      4.266666666666667*MassB*Sqr(g1)*(-1.2*Sqr(g1) + 6*QH2p*Qup*Sqr(gN) + 6*(Sqr(
      Yu00) + Sqr(Yu01) + Sqr(Yu02))) - 4.266666666666667*MassB*Sqr(g1)*(-1.2*Sqr(
      g1) + 6*QH2p*Qup*Sqr(gN) + 6*(Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12))) -
      0.26666666666666666*MassB*Sqr(g1)*(0.6*Sqr(g1) + 12*QH2p*QQp*Sqr(gN) + 6*(
      Sqr(Yu00) + Sqr(Yu10) + Sqr(Yu20))) - 0.26666666666666666*MassB*Sqr(g1)*(0.6
      *Sqr(g1) + 12*QH2p*QQp*Sqr(gN) + 6*(Sqr(Yu01) + Sqr(Yu11) + Sqr(Yu21))) -
      0.26666666666666666*MassB*Sqr(g1)*(0.6*Sqr(g1) + 12*QH2p*QQp*Sqr(gN) + 6*(
      Sqr(Yu02) + Sqr(Yu12) + Sqr(Yu22))) - 4.266666666666667*MassB*Sqr(g1)*(-1.2*
      Sqr(g1) + 6*QH2p*Qup*Sqr(gN) + 6*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) - 2.4*
      MassB*Sqr(g1)*(0.6*Sqr(g1) + 2*Sqr(Lambdax) + 4*Sqr(gN)*Sqr(QH2p) + 6*(Sqr(
      Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20
      ) + Sqr(Yu21) + Sqr(Yu22)));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dlead_log_one_loop_mHu2_dMassWB() const
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

   double deriv = -288*Power(g2,4)*MassWB + 24*Lambdax*TLambdax*Sqr(g2) + 72*
      TYu00*Yu00*Sqr(g2) + 72*TYu01*Yu01*Sqr(g2) + 72*TYu02*Yu02*Sqr(g2) + 72*
      TYu10*Yu10*Sqr(g2) + 72*TYu11*Yu11*Sqr(g2) + 72*TYu12*Yu12*Sqr(g2) + 72*
      TYu20*Yu20*Sqr(g2) + 72*TYu21*Yu21*Sqr(g2) + 72*TYu22*Yu22*Sqr(g2) - 24*
      MassWB*Sqr(g2)*(-0.6*Sqr(g1) + 4*QH1p*QH2p*Sqr(gN)) - 12*MassWB*Sqr(g2)*(0.6
      *Sqr(g1) + 4*QH2p*QHpbarp*Sqr(gN)) - 12*MassWB*Sqr(g2)*(-0.6*Sqr(g1) + 4*
      QH2p*QHpp*Sqr(gN)) - 36*MassWB*Sqr(g2)*(-0.6*Sqr(g1) + 4*QH2p*QLp*Sqr(gN)) -
      12*MassWB*Sqr(g2)*(-0.6*Sqr(g1) + 4*QH1p*QH2p*Sqr(gN) + 2*Sqr(Lambdax)) -
      24*MassWB*Sqr(g2)*(0.6*Sqr(g1) + 4*Sqr(gN)*Sqr(QH2p)) - 12*MassWB*Sqr(g2)*(
      0.6*Sqr(g1) + 12*QH2p*QQp*Sqr(gN) + 6*(Sqr(Yu00) + Sqr(Yu10) + Sqr(Yu20))) -
      12*MassWB*Sqr(g2)*(0.6*Sqr(g1) + 12*QH2p*QQp*Sqr(gN) + 6*(Sqr(Yu01) + Sqr(
      Yu11) + Sqr(Yu21))) - 12*MassWB*Sqr(g2)*(0.6*Sqr(g1) + 12*QH2p*QQp*Sqr(gN) +
      6*(Sqr(Yu02) + Sqr(Yu12) + Sqr(Yu22))) - 12*MassWB*Sqr(g2)*(0.6*Sqr(g1) + 2
      *Sqr(Lambdax) + 4*Sqr(gN)*Sqr(QH2p) + 6*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) +
      Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dlead_log_one_loop_mHu2_dMassG() const
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

   double deriv = 128*TYu00*Yu00*Sqr(g3) + 128*TYu01*Yu01*Sqr(g3) + 128*TYu02*
      Yu02*Sqr(g3) + 128*TYu10*Yu10*Sqr(g3) + 128*TYu11*Yu11*Sqr(g3) + 128*TYu12*
      Yu12*Sqr(g3) + 128*TYu20*Yu20*Sqr(g3) + 128*TYu21*Yu21*Sqr(g3) + 128*TYu22*
      Yu22*Sqr(g3) - 64*MassG*Sqr(g3)*(0.6*Sqr(g1) + 6*Qdp*QH2p*Sqr(gN)) - 64*
      MassG*Sqr(g3)*(0.6*Sqr(g1) + 6*QDxbarp*QH2p*Sqr(gN)) - 64*MassG*Sqr(g3)*(
      -0.6*Sqr(g1) + 6*QDxp*QH2p*Sqr(gN)) - 21.333333333333332*MassG*Sqr(g3)*(-1.2
      *Sqr(g1) + 6*QH2p*Qup*Sqr(gN) + 6*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02))) -
      21.333333333333332*MassG*Sqr(g3)*(-1.2*Sqr(g1) + 6*QH2p*Qup*Sqr(gN) + 6*(Sqr
      (Yu10) + Sqr(Yu11) + Sqr(Yu12))) - 21.333333333333332*MassG*Sqr(g3)*(0.6*Sqr
      (g1) + 12*QH2p*QQp*Sqr(gN) + 6*(Sqr(Yu00) + Sqr(Yu10) + Sqr(Yu20))) -
      21.333333333333332*MassG*Sqr(g3)*(0.6*Sqr(g1) + 12*QH2p*QQp*Sqr(gN) + 6*(Sqr
      (Yu01) + Sqr(Yu11) + Sqr(Yu21))) - 21.333333333333332*MassG*Sqr(g3)*(0.6*Sqr
      (g1) + 12*QH2p*QQp*Sqr(gN) + 6*(Sqr(Yu02) + Sqr(Yu12) + Sqr(Yu22))) -
      21.333333333333332*MassG*Sqr(g3)*(-1.2*Sqr(g1) + 6*QH2p*Qup*Sqr(gN) + 6*(Sqr
      (Yu20) + Sqr(Yu21) + Sqr(Yu22)));


   return twoLoop*deriv;
}

double lowE6SSM_tuning_calculator::deriv_dlead_log_one_loop_mHu2_dMassBp() const
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

   double deriv = -64*Power(gN,4)*MassBp*QH2p*Power(QSp,3) - 48*MassBp*Sqr(gN)*
      (0.6*Sqr(g1) + 6*Qdp*QH2p*Sqr(gN))*Sqr(Qdp) - 48*MassBp*Sqr(gN)*(0.6*Sqr(g1)
      + 6*QDxbarp*QH2p*Sqr(gN))*Sqr(QDxbarp) - 48*MassBp*Sqr(gN)*(-0.6*Sqr(g1) +
      6*QDxp*QH2p*Sqr(gN))*Sqr(QDxp) - 48*MassBp*Sqr(gN)*(0.6*Sqr(g1) + 2*Qep*QH2p
      *Sqr(gN))*Sqr(Qep) - 32*MassBp*Sqr(gN)*(-0.6*Sqr(g1) + 4*QH1p*QH2p*Sqr(gN))*
      Sqr(QH1p) - 16*MassBp*Sqr(gN)*(-0.6*Sqr(g1) + 4*QH1p*QH2p*Sqr(gN) + 2*Sqr(
      Lambdax))*Sqr(QH1p) - 32*MassBp*Sqr(gN)*Sqr(QH2p)*(0.6*Sqr(g1) + 4*Sqr(gN)*
      Sqr(QH2p)) - 16*MassBp*Sqr(gN)*(0.6*Sqr(g1) + 4*QH2p*QHpbarp*Sqr(gN))*Sqr(
      QHpbarp) - 16*MassBp*Sqr(gN)*(-0.6*Sqr(g1) + 4*QH2p*QHpp*Sqr(gN))*Sqr(QHpp)
      - 48*MassBp*Sqr(gN)*(-0.6*Sqr(g1) + 4*QH2p*QLp*Sqr(gN))*Sqr(QLp) - 16*MassBp
      *Sqr(gN)*(2*QH2p*QSp*Sqr(gN) + 2*Sqr(Lambdax))*Sqr(QSp) + 4*TLambdax*(4*
      Lambdax*Sqr(gN)*Sqr(QH1p) + 4*Lambdax*Sqr(gN)*Sqr(QH2p) + 4*Lambdax*Sqr(gN)*
      Sqr(QSp)) - 96*Power(gN,4)*MassBp*Sqr(QH2p)*(9*Sqr(Qdp) + 9*Sqr(QDxbarp) + 9
      *Sqr(QDxp) + 3*Sqr(Qep) + 6*Sqr(QH1p) + 6*Sqr(QH2p) + 2*Sqr(QHpbarp) + 2*Sqr
      (QHpp) + 6*Sqr(QLp) + 18*Sqr(QQp) + 3*Sqr(QSp) + 9*Sqr(Qup)) + 12*TYu00*Yu00
      *(4*Sqr(gN)*Sqr(QH2p) + 4*Sqr(gN)*Sqr(QQp) + 4*Sqr(gN)*Sqr(Qup)) + 12*TYu01*
      Yu01*(4*Sqr(gN)*Sqr(QH2p) + 4*Sqr(gN)*Sqr(QQp) + 4*Sqr(gN)*Sqr(Qup)) + 12*
      TYu02*Yu02*(4*Sqr(gN)*Sqr(QH2p) + 4*Sqr(gN)*Sqr(QQp) + 4*Sqr(gN)*Sqr(Qup)) +
      12*TYu10*Yu10*(4*Sqr(gN)*Sqr(QH2p) + 4*Sqr(gN)*Sqr(QQp) + 4*Sqr(gN)*Sqr(Qup
      )) + 12*TYu11*Yu11*(4*Sqr(gN)*Sqr(QH2p) + 4*Sqr(gN)*Sqr(QQp) + 4*Sqr(gN)*Sqr
      (Qup)) + 12*TYu12*Yu12*(4*Sqr(gN)*Sqr(QH2p) + 4*Sqr(gN)*Sqr(QQp) + 4*Sqr(gN)
      *Sqr(Qup)) + 12*TYu20*Yu20*(4*Sqr(gN)*Sqr(QH2p) + 4*Sqr(gN)*Sqr(QQp) + 4*Sqr
      (gN)*Sqr(Qup)) + 12*TYu21*Yu21*(4*Sqr(gN)*Sqr(QH2p) + 4*Sqr(gN)*Sqr(QQp) + 4
      *Sqr(gN)*Sqr(Qup)) + 12*TYu22*Yu22*(4*Sqr(gN)*Sqr(QH2p) + 4*Sqr(gN)*Sqr(QQp)
      + 4*Sqr(gN)*Sqr(Qup)) - 16*MassBp*Sqr(gN)*Sqr(Qup)*(-1.2*Sqr(g1) + 6*QH2p*
      Qup*Sqr(gN) + 6*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02))) - 16*MassBp*Sqr(gN)*Sqr
      (Qup)*(-1.2*Sqr(g1) + 6*QH2p*Qup*Sqr(gN) + 6*(Sqr(Yu10) + Sqr(Yu11) + Sqr(
      Yu12))) - 16*MassBp*Sqr(gN)*Sqr(QQp)*(0.6*Sqr(g1) + 12*QH2p*QQp*Sqr(gN) + 6*
      (Sqr(Yu00) + Sqr(Yu10) + Sqr(Yu20))) - 16*MassBp*Sqr(gN)*Sqr(QQp)*(0.6*Sqr(
      g1) + 12*QH2p*QQp*Sqr(gN) + 6*(Sqr(Yu01) + Sqr(Yu11) + Sqr(Yu21))) - 16*
      MassBp*Sqr(gN)*Sqr(QQp)*(0.6*Sqr(g1) + 12*QH2p*QQp*Sqr(gN) + 6*(Sqr(Yu02) +
      Sqr(Yu12) + Sqr(Yu22))) - 16*MassBp*Sqr(gN)*Sqr(Qup)*(-1.2*Sqr(g1) + 6*QH2p*
      Qup*Sqr(gN) + 6*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) - 16*MassBp*Sqr(gN)*Sqr
      (QH2p)*(0.6*Sqr(g1) + 2*Sqr(Lambdax) + 4*Sqr(gN)*Sqr(QH2p) + 6*(Sqr(Yu00) +
      Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(
      Yu21) + Sqr(Yu22)));


   return twoLoop*deriv;
}

} // namespace flexiblesusy
