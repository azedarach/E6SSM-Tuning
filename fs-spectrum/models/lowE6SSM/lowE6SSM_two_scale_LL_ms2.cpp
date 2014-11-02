#include "lowE6SSM_two_scale_tuning_calculator.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

double lowE6SSM_tuning_calculator::leading_log_coefficient_ms2() const
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

   double coeff = 6*(Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12)*(
      Kappa10*(Kappa00*mDx200 + Kappa10*mDx201 + Kappa20*mDx202) + Kappa11*(
      Kappa01*mDx200 + Kappa11*mDx201 + Kappa21*mDx202) + Kappa12*(Kappa02*mDx200
      + Kappa12*mDx201 + Kappa22*mDx202) + (Kappa00*Kappa10 + Kappa01*Kappa11 +
      Kappa02*Kappa12)*mDx211 + (Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*
      Kappa22)*mDx221 + 2*(Kappa10*(Kappa00*mDxbar200 + Kappa01*mDxbar210 +
      Kappa02*mDxbar220) + Kappa11*(Kappa00*mDxbar201 + Kappa01*mDxbar211 +
      Kappa02*mDxbar221) + Kappa12*(Kappa00*mDxbar202 + Kappa01*mDxbar212 +
      Kappa02*mDxbar222)) + 2*(Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12
      )*ms2 + 2*(TKappa00*TKappa10 + TKappa01*TKappa11 + TKappa02*TKappa12) +
      mDx201*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02))) + 6*(Kappa00*Kappa20 +
      Kappa01*Kappa21 + Kappa02*Kappa22)*(Kappa20*(Kappa00*mDx200 + Kappa10*mDx201
      + Kappa20*mDx202) + Kappa21*(Kappa01*mDx200 + Kappa11*mDx201 + Kappa21*
      mDx202) + Kappa22*(Kappa02*mDx200 + Kappa12*mDx201 + Kappa22*mDx202) + (
      Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12)*mDx212 + (Kappa00*
      Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)*mDx222 + 2*(Kappa20*(Kappa00*
      mDxbar200 + Kappa01*mDxbar210 + Kappa02*mDxbar220) + Kappa21*(Kappa00*
      mDxbar201 + Kappa01*mDxbar211 + Kappa02*mDxbar221) + Kappa22*(Kappa00*
      mDxbar202 + Kappa01*mDxbar212 + Kappa02*mDxbar222)) + 2*(Kappa00*Kappa20 +
      Kappa01*Kappa21 + Kappa02*Kappa22)*ms2 + 2*(TKappa00*TKappa20 + TKappa01*
      TKappa21 + TKappa02*TKappa22) + mDx202*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(
      Kappa02))) + 6*(Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12)*((
      Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12)*mDx200 + Kappa00*(
      Kappa00*mDx210 + Kappa10*mDx211 + Kappa20*mDx212) + Kappa01*(Kappa01*mDx210
      + Kappa11*mDx211 + Kappa21*mDx212) + Kappa02*(Kappa02*mDx210 + Kappa12*
      mDx211 + Kappa22*mDx212) + (Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*
      Kappa22)*mDx220 + 2*(Kappa00*(Kappa10*mDxbar200 + Kappa11*mDxbar210 +
      Kappa12*mDxbar220) + Kappa01*(Kappa10*mDxbar201 + Kappa11*mDxbar211 +
      Kappa12*mDxbar221) + Kappa02*(Kappa10*mDxbar202 + Kappa11*mDxbar212 +
      Kappa12*mDxbar222)) + 2*(Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12
      )*ms2 + 2*(TKappa00*TKappa10 + TKappa01*TKappa11 + TKappa02*TKappa12) +
      mDx210*(Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12))) + 6*(Kappa10*Kappa20 +
      Kappa11*Kappa21 + Kappa12*Kappa22)*((Kappa00*Kappa10 + Kappa01*Kappa11 +
      Kappa02*Kappa12)*mDx202 + Kappa20*(Kappa00*mDx210 + Kappa10*mDx211 + Kappa20
      *mDx212) + Kappa21*(Kappa01*mDx210 + Kappa11*mDx211 + Kappa21*mDx212) +
      Kappa22*(Kappa02*mDx210 + Kappa12*mDx211 + Kappa22*mDx212) + (Kappa10*
      Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)*mDx222 + 2*(Kappa20*(Kappa10*
      mDxbar200 + Kappa11*mDxbar210 + Kappa12*mDxbar220) + Kappa21*(Kappa10*
      mDxbar201 + Kappa11*mDxbar211 + Kappa12*mDxbar221) + Kappa22*(Kappa10*
      mDxbar202 + Kappa11*mDxbar212 + Kappa12*mDxbar222)) + 2*(Kappa10*Kappa20 +
      Kappa11*Kappa21 + Kappa12*Kappa22)*ms2 + 2*(TKappa10*TKappa20 + TKappa11*
      TKappa21 + TKappa12*TKappa22) + mDx212*(Sqr(Kappa10) + Sqr(Kappa11) + Sqr(
      Kappa12))) + 6*(Kappa00*Kappa01 + Kappa10*Kappa11 + Kappa20*Kappa21)*(2*(
      Kappa01*(Kappa00*mDx200 + Kappa10*mDx210 + Kappa20*mDx220) + Kappa11*(
      Kappa00*mDx201 + Kappa10*mDx211 + Kappa20*mDx221) + Kappa21*(Kappa00*mDx202
      + Kappa10*mDx212 + Kappa20*mDx222)) + Kappa01*(Kappa00*mDxbar200 + Kappa01*
      mDxbar201 + Kappa02*mDxbar202) + Kappa11*(Kappa10*mDxbar200 + Kappa11*
      mDxbar201 + Kappa12*mDxbar202) + Kappa21*(Kappa20*mDxbar200 + Kappa21*
      mDxbar201 + Kappa22*mDxbar202) + (Kappa00*Kappa01 + Kappa10*Kappa11 +
      Kappa20*Kappa21)*mDxbar211 + (Kappa00*Kappa02 + Kappa10*Kappa12 + Kappa20*
      Kappa22)*mDxbar221 + 2*(Kappa00*Kappa01 + Kappa10*Kappa11 + Kappa20*Kappa21)
      *ms2 + 2*(TKappa00*TKappa01 + TKappa10*TKappa11 + TKappa20*TKappa21) +
      mDxbar201*(Sqr(Kappa00) + Sqr(Kappa10) + Sqr(Kappa20))) + 6*(Kappa00*Kappa02
      + Kappa10*Kappa12 + Kappa20*Kappa22)*(2*(Kappa02*(Kappa00*mDx200 + Kappa10*
      mDx210 + Kappa20*mDx220) + Kappa12*(Kappa00*mDx201 + Kappa10*mDx211 +
      Kappa20*mDx221) + Kappa22*(Kappa00*mDx202 + Kappa10*mDx212 + Kappa20*mDx222)
      ) + Kappa02*(Kappa00*mDxbar200 + Kappa01*mDxbar201 + Kappa02*mDxbar202) +
      Kappa12*(Kappa10*mDxbar200 + Kappa11*mDxbar201 + Kappa12*mDxbar202) +
      Kappa22*(Kappa20*mDxbar200 + Kappa21*mDxbar201 + Kappa22*mDxbar202) + (
      Kappa00*Kappa01 + Kappa10*Kappa11 + Kappa20*Kappa21)*mDxbar212 + (Kappa00*
      Kappa02 + Kappa10*Kappa12 + Kappa20*Kappa22)*mDxbar222 + 2*(Kappa00*Kappa02
      + Kappa10*Kappa12 + Kappa20*Kappa22)*ms2 + 2*(TKappa00*TKappa02 + TKappa10*
      TKappa12 + TKappa20*TKappa22) + mDxbar202*(Sqr(Kappa00) + Sqr(Kappa10) + Sqr
      (Kappa20))) + 6*(Kappa00*Kappa01 + Kappa10*Kappa11 + Kappa20*Kappa21)*(2*(
      Kappa00*(Kappa01*mDx200 + Kappa11*mDx210 + Kappa21*mDx220) + Kappa10*(
      Kappa01*mDx201 + Kappa11*mDx211 + Kappa21*mDx221) + Kappa20*(Kappa01*mDx202
      + Kappa11*mDx212 + Kappa21*mDx222)) + (Kappa00*Kappa01 + Kappa10*Kappa11 +
      Kappa20*Kappa21)*mDxbar200 + Kappa00*(Kappa00*mDxbar210 + Kappa01*mDxbar211
      + Kappa02*mDxbar212) + Kappa10*(Kappa10*mDxbar210 + Kappa11*mDxbar211 +
      Kappa12*mDxbar212) + Kappa20*(Kappa20*mDxbar210 + Kappa21*mDxbar211 +
      Kappa22*mDxbar212) + (Kappa01*Kappa02 + Kappa11*Kappa12 + Kappa21*Kappa22)*
      mDxbar220 + 2*(Kappa00*Kappa01 + Kappa10*Kappa11 + Kappa20*Kappa21)*ms2 + 2*
      (TKappa00*TKappa01 + TKappa10*TKappa11 + TKappa20*TKappa21) + mDxbar210*(Sqr
      (Kappa01) + Sqr(Kappa11) + Sqr(Kappa21))) + 6*(Kappa01*Kappa02 + Kappa11*
      Kappa12 + Kappa21*Kappa22)*(2*(Kappa02*(Kappa01*mDx200 + Kappa11*mDx210 +
      Kappa21*mDx220) + Kappa12*(Kappa01*mDx201 + Kappa11*mDx211 + Kappa21*mDx221)
      + Kappa22*(Kappa01*mDx202 + Kappa11*mDx212 + Kappa21*mDx222)) + (Kappa00*
      Kappa01 + Kappa10*Kappa11 + Kappa20*Kappa21)*mDxbar202 + Kappa02*(Kappa00*
      mDxbar210 + Kappa01*mDxbar211 + Kappa02*mDxbar212) + Kappa12*(Kappa10*
      mDxbar210 + Kappa11*mDxbar211 + Kappa12*mDxbar212) + Kappa22*(Kappa20*
      mDxbar210 + Kappa21*mDxbar211 + Kappa22*mDxbar212) + (Kappa01*Kappa02 +
      Kappa11*Kappa12 + Kappa21*Kappa22)*mDxbar222 + 2*(Kappa01*Kappa02 + Kappa11*
      Kappa12 + Kappa21*Kappa22)*ms2 + 2*(TKappa01*TKappa02 + TKappa11*TKappa12 +
      TKappa21*TKappa22) + mDxbar212*(Sqr(Kappa01) + Sqr(Kappa11) + Sqr(Kappa21)))
      + 6*(Kappa00*Kappa02 + Kappa10*Kappa12 + Kappa20*Kappa22)*(2*(Kappa00*(
      Kappa02*mDx200 + Kappa12*mDx210 + Kappa22*mDx220) + Kappa10*(Kappa02*mDx201
      + Kappa12*mDx211 + Kappa22*mDx221) + Kappa20*(Kappa02*mDx202 + Kappa12*
      mDx212 + Kappa22*mDx222)) + (Kappa00*Kappa02 + Kappa10*Kappa12 + Kappa20*
      Kappa22)*mDxbar200 + (Kappa01*Kappa02 + Kappa11*Kappa12 + Kappa21*Kappa22)*
      mDxbar210 + Kappa00*(Kappa00*mDxbar220 + Kappa01*mDxbar221 + Kappa02*
      mDxbar222) + Kappa10*(Kappa10*mDxbar220 + Kappa11*mDxbar221 + Kappa12*
      mDxbar222) + Kappa20*(Kappa20*mDxbar220 + Kappa21*mDxbar221 + Kappa22*
      mDxbar222) + 2*(Kappa00*Kappa02 + Kappa10*Kappa12 + Kappa20*Kappa22)*ms2 + 2
      *(TKappa00*TKappa02 + TKappa10*TKappa12 + TKappa20*TKappa22) + mDxbar220*(
      Sqr(Kappa02) + Sqr(Kappa12) + Sqr(Kappa22))) + 6*(Kappa01*Kappa02 + Kappa11*
      Kappa12 + Kappa21*Kappa22)*(2*(Kappa01*(Kappa02*mDx200 + Kappa12*mDx210 +
      Kappa22*mDx220) + Kappa11*(Kappa02*mDx201 + Kappa12*mDx211 + Kappa22*mDx221)
      + Kappa21*(Kappa02*mDx202 + Kappa12*mDx212 + Kappa22*mDx222)) + (Kappa00*
      Kappa02 + Kappa10*Kappa12 + Kappa20*Kappa22)*mDxbar201 + (Kappa01*Kappa02 +
      Kappa11*Kappa12 + Kappa21*Kappa22)*mDxbar211 + Kappa01*(Kappa00*mDxbar220 +
      Kappa01*mDxbar221 + Kappa02*mDxbar222) + Kappa11*(Kappa10*mDxbar220 +
      Kappa11*mDxbar221 + Kappa12*mDxbar222) + Kappa21*(Kappa20*mDxbar220 +
      Kappa21*mDxbar221 + Kappa22*mDxbar222) + 2*(Kappa01*Kappa02 + Kappa11*
      Kappa12 + Kappa21*Kappa22)*ms2 + 2*(TKappa01*TKappa02 + TKappa11*TKappa12 +
      TKappa21*TKappa22) + mDxbar221*(Sqr(Kappa02) + Sqr(Kappa12) + Sqr(Kappa22)))
      + 6*(Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)*((Kappa00*Kappa20
      + Kappa01*Kappa21 + Kappa02*Kappa22)*mDx200 + (Kappa10*Kappa20 + Kappa11*
      Kappa21 + Kappa12*Kappa22)*mDx210 + Kappa00*(Kappa00*mDx220 + Kappa10*mDx221
      + Kappa20*mDx222) + Kappa01*(Kappa01*mDx220 + Kappa11*mDx221 + Kappa21*
      mDx222) + Kappa02*(Kappa02*mDx220 + Kappa12*mDx221 + Kappa22*mDx222) + 2*(
      Kappa00*(Kappa20*mDxbar200 + Kappa21*mDxbar210 + Kappa22*mDxbar220) +
      Kappa01*(Kappa20*mDxbar201 + Kappa21*mDxbar211 + Kappa22*mDxbar221) +
      Kappa02*(Kappa20*mDxbar202 + Kappa21*mDxbar212 + Kappa22*mDxbar222)) + 2*(
      Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)*ms2 + 2*(TKappa00*
      TKappa20 + TKappa01*TKappa21 + TKappa02*TKappa22) + mDx220*(Sqr(Kappa20) +
      Sqr(Kappa21) + Sqr(Kappa22))) + 6*(Kappa10*Kappa20 + Kappa11*Kappa21 +
      Kappa12*Kappa22)*((Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)*
      mDx201 + (Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)*mDx211 +
      Kappa10*(Kappa00*mDx220 + Kappa10*mDx221 + Kappa20*mDx222) + Kappa11*(
      Kappa01*mDx220 + Kappa11*mDx221 + Kappa21*mDx222) + Kappa12*(Kappa02*mDx220
      + Kappa12*mDx221 + Kappa22*mDx222) + 2*(Kappa10*(Kappa20*mDxbar200 + Kappa21
      *mDxbar210 + Kappa22*mDxbar220) + Kappa11*(Kappa20*mDxbar201 + Kappa21*
      mDxbar211 + Kappa22*mDxbar221) + Kappa12*(Kappa20*mDxbar202 + Kappa21*
      mDxbar212 + Kappa22*mDxbar222)) + 2*(Kappa10*Kappa20 + Kappa11*Kappa21 +
      Kappa12*Kappa22)*ms2 + 2*(TKappa10*TKappa20 + TKappa11*TKappa21 + TKappa12*
      TKappa22) + mDx221*(Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22))) + 4*(
      Lambda1200*Lambda1210 + Lambda1201*Lambda1211)*(2*(Lambda1210*(Lambda1200*
      mH1I200 + Lambda1201*mH1I210) + Lambda1211*(Lambda1200*mH1I201 + Lambda1201*
      mH1I211)) + Lambda1210*(Lambda1200*mH2I200 + Lambda1210*mH2I201) +
      Lambda1211*(Lambda1201*mH2I200 + Lambda1211*mH2I201) + (Lambda1200*
      Lambda1210 + Lambda1201*Lambda1211)*mH2I211 + 2*(Lambda1200*Lambda1210 +
      Lambda1201*Lambda1211)*ms2 + 2*(TLambda1200*TLambda1210 + TLambda1201*
      TLambda1211) + mH2I201*(Sqr(Lambda1200) + Sqr(Lambda1201))) + 4*(Lambda1200*
      Lambda1201 + Lambda1210*Lambda1211)*(Lambda1201*(Lambda1200*mH1I200 +
      Lambda1201*mH1I201) + Lambda1211*(Lambda1210*mH1I200 + Lambda1211*mH1I201) +
      (Lambda1200*Lambda1201 + Lambda1210*Lambda1211)*mH1I211 + 2*(Lambda1201*(
      Lambda1200*mH2I200 + Lambda1210*mH2I210) + Lambda1211*(Lambda1200*mH2I201 +
      Lambda1210*mH2I211)) + 2*(Lambda1200*Lambda1201 + Lambda1210*Lambda1211)*ms2
      + 2*(TLambda1200*TLambda1201 + TLambda1210*TLambda1211) + mH1I201*(Sqr(
      Lambda1200) + Sqr(Lambda1210))) + 4*(Lambda1200*Lambda1201 + Lambda1210*
      Lambda1211)*((Lambda1200*Lambda1201 + Lambda1210*Lambda1211)*mH1I200 +
      Lambda1200*(Lambda1200*mH1I210 + Lambda1201*mH1I211) + Lambda1210*(
      Lambda1210*mH1I210 + Lambda1211*mH1I211) + 2*(Lambda1200*(Lambda1201*mH2I200
      + Lambda1211*mH2I210) + Lambda1210*(Lambda1201*mH2I201 + Lambda1211*mH2I211
      )) + 2*(Lambda1200*Lambda1201 + Lambda1210*Lambda1211)*ms2 + 2*(TLambda1200*
      TLambda1201 + TLambda1210*TLambda1211) + mH1I210*(Sqr(Lambda1201) + Sqr(
      Lambda1211))) + 4*(Lambda1200*Lambda1210 + Lambda1201*Lambda1211)*(2*(
      Lambda1200*(Lambda1210*mH1I200 + Lambda1211*mH1I210) + Lambda1201*(
      Lambda1210*mH1I201 + Lambda1211*mH1I211)) + (Lambda1200*Lambda1210 +
      Lambda1201*Lambda1211)*mH2I200 + Lambda1200*(Lambda1200*mH2I210 + Lambda1210
      *mH2I211) + Lambda1201*(Lambda1201*mH2I210 + Lambda1211*mH2I211) + 2*(
      Lambda1200*Lambda1210 + Lambda1201*Lambda1211)*ms2 + 2*(TLambda1200*
      TLambda1210 + TLambda1201*TLambda1211) + mH2I210*(Sqr(Lambda1210) + Sqr(
      Lambda1211))) + 8*Power(gN,3)*Power(QSp,3)*(gN*(3*(md200 + md211 + md222)*
      Qdp + 3*(mDxbar200 + mDxbar211 + mDxbar222)*QDxbarp + 3*(mDx200 + mDx211 +
      mDx222)*QDxp + (me200 + me211 + me222)*Qep + 2*(mH1I200 + mH1I211)*QH1p + 2*
      mHd2*QH1p + 2*(mH2I200 + mH2I211)*QH2p + 2*mHu2*QH2p + 2*mHpbar2*QHpbarp + 2
      *mHp2*QHpp + 2*(ml200 + ml211 + ml222)*QLp + 6*(mq200 + mq211 + mq222)*QQp +
      ms2*QSp + (msI200 + msI211)*QSp + 3*(mu200 + mu211 + mu222)*Qup) - 4*gN*QSp
      *Sqr(MassBp)) + 4*QHpbarp*QSp*Sqr(gN)*(0.6*(md200 + md211 + md222 - mDx200 -
      mDx211 - mDx222 + mDxbar200 + mDxbar211 + mDxbar222 + me200 + me211 + me222
      - mH1I200 - mH1I211 + mH2I200 + mH2I211 - mHd2 - mHp2 + mHpbar2 + mHu2 -
      ml200 - ml211 - ml222 + mq200 + mq211 + mq222 - 2*(mu200 + mu211 + mu222))*
      Sqr(g1) + 2*QHpbarp*(3*(md200 + md211 + md222)*Qdp + 3*(mDxbar200 +
      mDxbar211 + mDxbar222)*QDxbarp + 3*(mDx200 + mDx211 + mDx222)*QDxp + (me200
      + me211 + me222)*Qep + 2*(mH1I200 + mH1I211)*QH1p + 2*mHd2*QH1p + 2*(mH2I200
      + mH2I211)*QH2p + 2*mHu2*QH2p + 2*mHpbar2*QHpbarp + 2*mHp2*QHpp + 2*(ml200
      + ml211 + ml222)*QLp + 6*(mq200 + mq211 + mq222)*QQp + ms2*QSp + (msI200 +
      msI211)*QSp + 3*(mu200 + mu211 + mu222)*Qup)*Sqr(gN) - 1.2*Sqr(g1)*Sqr(MassB
      ) - 6*Sqr(g2)*Sqr(MassWB) - 8*Sqr(gN)*Sqr(MassBp)*Sqr(QHpbarp)) + 4*QHpp*QSp
      *Sqr(gN)*(-0.6*(md200 + md211 + md222 - mDx200 - mDx211 - mDx222 + mDxbar200
      + mDxbar211 + mDxbar222 + me200 + me211 + me222 - mH1I200 - mH1I211 +
      mH2I200 + mH2I211 - mHd2 - mHp2 + mHpbar2 + mHu2 - ml200 - ml211 - ml222 +
      mq200 + mq211 + mq222 - 2*(mu200 + mu211 + mu222))*Sqr(g1) + 2*QHpp*(3*(
      md200 + md211 + md222)*Qdp + 3*(mDxbar200 + mDxbar211 + mDxbar222)*QDxbarp +
      3*(mDx200 + mDx211 + mDx222)*QDxp + (me200 + me211 + me222)*Qep + 2*(
      mH1I200 + mH1I211)*QH1p + 2*mHd2*QH1p + 2*(mH2I200 + mH2I211)*QH2p + 2*mHu2*
      QH2p + 2*mHpbar2*QHpbarp + 2*mHp2*QHpp + 2*(ml200 + ml211 + ml222)*QLp + 6*(
      mq200 + mq211 + mq222)*QQp + ms2*QSp + (msI200 + msI211)*QSp + 3*(mu200 +
      mu211 + mu222)*Qup)*Sqr(gN) - 1.2*Sqr(g1)*Sqr(MassB) - 6*Sqr(g2)*Sqr(MassWB)
      - 8*Sqr(gN)*Sqr(MassBp)*Sqr(QHpp)) + (6*(2*Kappa00*mDx200 + Kappa10*mDx201
      + Kappa20*mDx202 + Kappa10*mDx210 + Kappa20*mDx220) + 6*(2*Kappa00*mDxbar200
      + Kappa01*mDxbar201 + Kappa02*mDxbar202 + Kappa01*mDxbar210 + Kappa02*
      mDxbar220) + 12*Kappa00*ms2)*(2*(Kappa10*(Kappa00*Kappa10 + Kappa01*Kappa11
      + Kappa02*Kappa12) + Kappa20*(Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*
      Kappa22) + Kappa00*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02))) + Kappa00*(
      -0.26666666666666666*Sqr(g1) - 5.333333333333333*Sqr(g3) + 3*(Sqr(Kappa00) +
      Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) +
      Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*(Sqr(Lambda1200) + Sqr(
      Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 2*Sqr(Lambdax) - 2*Sqr(gN
      )*Sqr(QDxbarp) - 2*Sqr(gN)*Sqr(QDxp) - 2*Sqr(gN)*Sqr(QSp))) + (6*(2*Kappa01*
      mDx200 + Kappa11*mDx201 + Kappa21*mDx202 + Kappa11*mDx210 + Kappa21*mDx220)
      + 6*(Kappa00*mDxbar201 + Kappa00*mDxbar210 + 2*Kappa01*mDxbar211 + Kappa02*
      mDxbar212 + Kappa02*mDxbar221) + 12*Kappa01*ms2)*(2*(Kappa11*(Kappa00*
      Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12) + Kappa21*(Kappa00*Kappa20 +
      Kappa01*Kappa21 + Kappa02*Kappa22) + Kappa01*(Sqr(Kappa00) + Sqr(Kappa01) +
      Sqr(Kappa02))) + Kappa01*(-0.26666666666666666*Sqr(g1) - 5.333333333333333*
      Sqr(g3) + 3*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr
      (Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*(
      Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 2*
      Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QDxbarp) - 2*Sqr(gN)*Sqr(QDxp) - 2*Sqr(gN)*Sqr(
      QSp))) + (6*(2*Kappa02*mDx200 + Kappa12*mDx201 + Kappa22*mDx202 + Kappa12*
      mDx210 + Kappa22*mDx220) + 6*(Kappa00*mDxbar202 + Kappa01*mDxbar212 +
      Kappa00*mDxbar220 + Kappa01*mDxbar221 + 2*Kappa02*mDxbar222) + 12*Kappa02*
      ms2)*(2*(Kappa12*(Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12) +
      Kappa22*(Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22) + Kappa02*(Sqr
      (Kappa00) + Sqr(Kappa01) + Sqr(Kappa02))) + Kappa02*(-0.26666666666666666*
      Sqr(g1) - 5.333333333333333*Sqr(g3) + 3*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(
      Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(
      Kappa21) + Sqr(Kappa22)) + 2*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(
      Lambda1210) + Sqr(Lambda1211)) + 2*Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QDxbarp) - 2
      *Sqr(gN)*Sqr(QDxp) - 2*Sqr(gN)*Sqr(QSp))) + (6*(Kappa00*mDx201 + Kappa00*
      mDx210 + 2*Kappa10*mDx211 + Kappa20*mDx212 + Kappa20*mDx221) + 6*(2*Kappa10*
      mDxbar200 + Kappa11*mDxbar201 + Kappa12*mDxbar202 + Kappa11*mDxbar210 +
      Kappa12*mDxbar220) + 12*Kappa10*ms2)*(2*(Kappa00*(Kappa00*Kappa10 + Kappa01*
      Kappa11 + Kappa02*Kappa12) + Kappa20*(Kappa10*Kappa20 + Kappa11*Kappa21 +
      Kappa12*Kappa22) + Kappa10*(Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12))) +
      Kappa10*(-0.26666666666666666*Sqr(g1) - 5.333333333333333*Sqr(g3) + 3*(Sqr(
      Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(
      Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*(Sqr(Lambda1200)
      + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 2*Sqr(Lambdax) - 2*
      Sqr(gN)*Sqr(QDxbarp) - 2*Sqr(gN)*Sqr(QDxp) - 2*Sqr(gN)*Sqr(QSp))) + (6*(
      Kappa01*mDx201 + Kappa01*mDx210 + 2*Kappa11*mDx211 + Kappa21*mDx212 +
      Kappa21*mDx221) + 6*(Kappa10*mDxbar201 + Kappa10*mDxbar210 + 2*Kappa11*
      mDxbar211 + Kappa12*mDxbar212 + Kappa12*mDxbar221) + 12*Kappa11*ms2)*(2*(
      Kappa01*(Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12) + Kappa21*(
      Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22) + Kappa11*(Sqr(Kappa10)
      + Sqr(Kappa11) + Sqr(Kappa12))) + Kappa11*(-0.26666666666666666*Sqr(g1) -
      5.333333333333333*Sqr(g3) + 3*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) +
      Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) +
      Sqr(Kappa22)) + 2*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr
      (Lambda1211)) + 2*Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QDxbarp) - 2*Sqr(gN)*Sqr(QDxp
      ) - 2*Sqr(gN)*Sqr(QSp))) + (6*(Kappa02*mDx201 + Kappa02*mDx210 + 2*Kappa12*
      mDx211 + Kappa22*mDx212 + Kappa22*mDx221) + 6*(Kappa10*mDxbar202 + Kappa11*
      mDxbar212 + Kappa10*mDxbar220 + Kappa11*mDxbar221 + 2*Kappa12*mDxbar222) +
      12*Kappa12*ms2)*(2*(Kappa02*(Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*
      Kappa12) + Kappa22*(Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22) +
      Kappa12*(Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12))) + Kappa12*(
      -0.26666666666666666*Sqr(g1) - 5.333333333333333*Sqr(g3) + 3*(Sqr(Kappa00) +
      Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) +
      Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*(Sqr(Lambda1200) + Sqr(
      Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 2*Sqr(Lambdax) - 2*Sqr(gN
      )*Sqr(QDxbarp) - 2*Sqr(gN)*Sqr(QDxp) - 2*Sqr(gN)*Sqr(QSp))) + (6*(Kappa00*
      mDx202 + Kappa10*mDx212 + Kappa00*mDx220 + Kappa10*mDx221 + 2*Kappa20*mDx222
      ) + 6*(2*Kappa20*mDxbar200 + Kappa21*mDxbar201 + Kappa22*mDxbar202 + Kappa21
      *mDxbar210 + Kappa22*mDxbar220) + 12*Kappa20*ms2)*(2*(Kappa00*(Kappa00*
      Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22) + Kappa10*(Kappa10*Kappa20 +
      Kappa11*Kappa21 + Kappa12*Kappa22) + Kappa20*(Sqr(Kappa20) + Sqr(Kappa21) +
      Sqr(Kappa22))) + Kappa20*(-0.26666666666666666*Sqr(g1) - 5.333333333333333*
      Sqr(g3) + 3*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr
      (Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*(
      Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 2*
      Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QDxbarp) - 2*Sqr(gN)*Sqr(QDxp) - 2*Sqr(gN)*Sqr(
      QSp))) + (6*(Kappa01*mDx202 + Kappa11*mDx212 + Kappa01*mDx220 + Kappa11*
      mDx221 + 2*Kappa21*mDx222) + 6*(Kappa20*mDxbar201 + Kappa20*mDxbar210 + 2*
      Kappa21*mDxbar211 + Kappa22*mDxbar212 + Kappa22*mDxbar221) + 12*Kappa21*ms2)
      *(2*(Kappa01*(Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22) + Kappa11
      *(Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22) + Kappa21*(Sqr(
      Kappa20) + Sqr(Kappa21) + Sqr(Kappa22))) + Kappa21*(-0.26666666666666666*Sqr
      (g1) - 5.333333333333333*Sqr(g3) + 3*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(
      Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(
      Kappa21) + Sqr(Kappa22)) + 2*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(
      Lambda1210) + Sqr(Lambda1211)) + 2*Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QDxbarp) - 2
      *Sqr(gN)*Sqr(QDxp) - 2*Sqr(gN)*Sqr(QSp))) + (6*(Kappa02*mDx202 + Kappa12*
      mDx212 + Kappa02*mDx220 + Kappa12*mDx221 + 2*Kappa22*mDx222) + 6*(Kappa20*
      mDxbar202 + Kappa21*mDxbar212 + Kappa20*mDxbar220 + Kappa21*mDxbar221 + 2*
      Kappa22*mDxbar222) + 12*Kappa22*ms2)*(2*(Kappa02*(Kappa00*Kappa20 + Kappa01*
      Kappa21 + Kappa02*Kappa22) + Kappa12*(Kappa10*Kappa20 + Kappa11*Kappa21 +
      Kappa12*Kappa22) + Kappa22*(Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22))) +
      Kappa22*(-0.26666666666666666*Sqr(g1) - 5.333333333333333*Sqr(g3) + 3*(Sqr(
      Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(
      Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*(Sqr(Lambda1200)
      + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 2*Sqr(Lambdax) - 2*
      Sqr(gN)*Sqr(QDxbarp) - 2*Sqr(gN)*Sqr(QDxp) - 2*Sqr(gN)*Sqr(QSp))) + (4*(2*
      Lambda1200*mH1I200 + Lambda1201*mH1I201 + Lambda1201*mH1I210) + 4*(2*
      Lambda1200*mH2I200 + Lambda1210*mH2I201 + Lambda1210*mH2I210) + 8*Lambda1200
      *ms2)*(2*(Lambda1210*(Lambda1200*Lambda1210 + Lambda1201*Lambda1211) +
      Lambda1200*(Sqr(Lambda1200) + Sqr(Lambda1201))) + Lambda1200*(-0.6*Sqr(g1) -
      3*Sqr(g2) + 3*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) +
      Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) +
      2*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) +
      2*Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(
      QSp))) + (4*(Lambda1200*mH1I201 + Lambda1200*mH1I210 + 2*Lambda1201*mH1I211)
      + 4*(2*Lambda1201*mH2I200 + Lambda1211*mH2I201 + Lambda1211*mH2I210) + 8*
      Lambda1201*ms2)*(2*(Lambda1211*(Lambda1200*Lambda1210 + Lambda1201*
      Lambda1211) + Lambda1201*(Sqr(Lambda1200) + Sqr(Lambda1201))) + Lambda1201*(
      -0.6*Sqr(g1) - 3*Sqr(g2) + 3*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) +
      Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) +
      Sqr(Kappa22)) + 2*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr
      (Lambda1211)) + 2*Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QH2p) -
      2*Sqr(gN)*Sqr(QSp))) + (4*(2*Lambda1210*mH1I200 + Lambda1211*mH1I201 +
      Lambda1211*mH1I210) + 4*(Lambda1200*mH2I201 + Lambda1200*mH2I210 + 2*
      Lambda1210*mH2I211) + 8*Lambda1210*ms2)*(2*(Lambda1200*(Lambda1200*
      Lambda1210 + Lambda1201*Lambda1211) + Lambda1210*(Sqr(Lambda1210) + Sqr(
      Lambda1211))) + Lambda1210*(-0.6*Sqr(g1) - 3*Sqr(g2) + 3*(Sqr(Kappa00) + Sqr
      (Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(
      Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*(Sqr(Lambda1200) + Sqr(
      Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 2*Sqr(Lambdax) - 2*Sqr(gN
      )*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QSp))) + (4*(Lambda1210*
      mH1I201 + Lambda1210*mH1I210 + 2*Lambda1211*mH1I211) + 4*(Lambda1201*mH2I201
      + Lambda1201*mH2I210 + 2*Lambda1211*mH2I211) + 8*Lambda1211*ms2)*(2*(
      Lambda1201*(Lambda1200*Lambda1210 + Lambda1201*Lambda1211) + Lambda1211*(Sqr
      (Lambda1210) + Sqr(Lambda1211))) + Lambda1211*(-0.6*Sqr(g1) - 3*Sqr(g2) + 3*
      (Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) +
      Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*(Sqr(
      Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 2*Sqr(
      Lambdax) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QSp)))
      + 12*TKappa00*(3*(Kappa00*(Kappa00*TKappa00 + Kappa01*TKappa01 + Kappa02*
      TKappa02) + Kappa10*(Kappa10*TKappa00 + Kappa11*TKappa01 + Kappa12*TKappa02)
      + Kappa20*(Kappa20*TKappa00 + Kappa21*TKappa01 + Kappa22*TKappa02)) -
      0.26666666666666666*TKappa00*Sqr(g1) - 5.333333333333333*TKappa00*Sqr(g3) +
      3*((Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12)*TKappa10 + (Kappa00
      *Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)*TKappa20 + TKappa00*(Sqr(
      Kappa00) + Sqr(Kappa01) + Sqr(Kappa02))) + 3*TKappa00*(Sqr(Kappa00) + Sqr(
      Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(
      Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*TKappa00*(Sqr(Lambda1200) + Sqr(
      Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 2*TKappa00*Sqr(Lambdax) -
      2*TKappa00*Sqr(gN)*Sqr(QDxbarp) - 2*TKappa00*Sqr(gN)*Sqr(QDxp) - 2*TKappa00
      *Sqr(gN)*Sqr(QSp) + Kappa00*(6*(Kappa00*TKappa00 + Kappa01*TKappa01 +
      Kappa02*TKappa02 + Kappa10*TKappa10 + Kappa11*TKappa11 + Kappa12*TKappa12 +
      Kappa20*TKappa20 + Kappa21*TKappa21 + Kappa22*TKappa22) + 4*(Lambda1200*
      TLambda1200 + Lambda1201*TLambda1201 + Lambda1210*TLambda1210 + Lambda1211*
      TLambda1211) + 4*Lambdax*TLambdax + 0.5333333333333333*MassB*Sqr(g1) +
      10.666666666666666*MassG*Sqr(g3) + 4*MassBp*Sqr(gN)*Sqr(QDxbarp) + 4*MassBp*
      Sqr(gN)*Sqr(QDxp) + 4*MassBp*Sqr(gN)*Sqr(QSp))) + 12*TKappa01*(3*(Kappa01*(
      Kappa00*TKappa00 + Kappa01*TKappa01 + Kappa02*TKappa02) + Kappa11*(Kappa10*
      TKappa00 + Kappa11*TKappa01 + Kappa12*TKappa02) + Kappa21*(Kappa20*TKappa00
      + Kappa21*TKappa01 + Kappa22*TKappa02)) - 0.26666666666666666*TKappa01*Sqr(
      g1) - 5.333333333333333*TKappa01*Sqr(g3) + 3*((Kappa00*Kappa10 + Kappa01*
      Kappa11 + Kappa02*Kappa12)*TKappa11 + (Kappa00*Kappa20 + Kappa01*Kappa21 +
      Kappa02*Kappa22)*TKappa21 + TKappa01*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(
      Kappa02))) + 3*TKappa01*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(
      Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(
      Kappa22)) + 2*TKappa01*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210)
      + Sqr(Lambda1211)) + 2*TKappa01*Sqr(Lambdax) - 2*TKappa01*Sqr(gN)*Sqr(
      QDxbarp) - 2*TKappa01*Sqr(gN)*Sqr(QDxp) - 2*TKappa01*Sqr(gN)*Sqr(QSp) +
      Kappa01*(6*(Kappa00*TKappa00 + Kappa01*TKappa01 + Kappa02*TKappa02 + Kappa10
      *TKappa10 + Kappa11*TKappa11 + Kappa12*TKappa12 + Kappa20*TKappa20 + Kappa21
      *TKappa21 + Kappa22*TKappa22) + 4*(Lambda1200*TLambda1200 + Lambda1201*
      TLambda1201 + Lambda1210*TLambda1210 + Lambda1211*TLambda1211) + 4*Lambdax*
      TLambdax + 0.5333333333333333*MassB*Sqr(g1) + 10.666666666666666*MassG*Sqr(
      g3) + 4*MassBp*Sqr(gN)*Sqr(QDxbarp) + 4*MassBp*Sqr(gN)*Sqr(QDxp) + 4*MassBp*
      Sqr(gN)*Sqr(QSp))) + 12*TKappa02*(3*(Kappa02*(Kappa00*TKappa00 + Kappa01*
      TKappa01 + Kappa02*TKappa02) + Kappa12*(Kappa10*TKappa00 + Kappa11*TKappa01
      + Kappa12*TKappa02) + Kappa22*(Kappa20*TKappa00 + Kappa21*TKappa01 + Kappa22
      *TKappa02)) - 0.26666666666666666*TKappa02*Sqr(g1) - 5.333333333333333*
      TKappa02*Sqr(g3) + 3*((Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12)*
      TKappa12 + (Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)*TKappa22 +
      TKappa02*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02))) + 3*TKappa02*(Sqr(
      Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(
      Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*TKappa02*(Sqr(
      Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 2*
      TKappa02*Sqr(Lambdax) - 2*TKappa02*Sqr(gN)*Sqr(QDxbarp) - 2*TKappa02*Sqr(gN)
      *Sqr(QDxp) - 2*TKappa02*Sqr(gN)*Sqr(QSp) + Kappa02*(6*(Kappa00*TKappa00 +
      Kappa01*TKappa01 + Kappa02*TKappa02 + Kappa10*TKappa10 + Kappa11*TKappa11 +
      Kappa12*TKappa12 + Kappa20*TKappa20 + Kappa21*TKappa21 + Kappa22*TKappa22) +
      4*(Lambda1200*TLambda1200 + Lambda1201*TLambda1201 + Lambda1210*TLambda1210
      + Lambda1211*TLambda1211) + 4*Lambdax*TLambdax + 0.5333333333333333*MassB*
      Sqr(g1) + 10.666666666666666*MassG*Sqr(g3) + 4*MassBp*Sqr(gN)*Sqr(QDxbarp) +
      4*MassBp*Sqr(gN)*Sqr(QDxp) + 4*MassBp*Sqr(gN)*Sqr(QSp))) + 12*TKappa10*(3*(
      Kappa00*(Kappa00*TKappa10 + Kappa01*TKappa11 + Kappa02*TKappa12) + Kappa10*(
      Kappa10*TKappa10 + Kappa11*TKappa11 + Kappa12*TKappa12) + Kappa20*(Kappa20*
      TKappa10 + Kappa21*TKappa11 + Kappa22*TKappa12)) - 0.26666666666666666*
      TKappa10*Sqr(g1) - 5.333333333333333*TKappa10*Sqr(g3) + 3*((Kappa00*Kappa10
      + Kappa01*Kappa11 + Kappa02*Kappa12)*TKappa00 + (Kappa10*Kappa20 + Kappa11*
      Kappa21 + Kappa12*Kappa22)*TKappa20 + TKappa10*(Sqr(Kappa10) + Sqr(Kappa11)
      + Sqr(Kappa12))) + 3*TKappa10*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) +
      Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) +
      Sqr(Kappa22)) + 2*TKappa10*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(
      Lambda1210) + Sqr(Lambda1211)) + 2*TKappa10*Sqr(Lambdax) - 2*TKappa10*Sqr(gN
      )*Sqr(QDxbarp) - 2*TKappa10*Sqr(gN)*Sqr(QDxp) - 2*TKappa10*Sqr(gN)*Sqr(QSp)
      + Kappa10*(6*(Kappa00*TKappa00 + Kappa01*TKappa01 + Kappa02*TKappa02 +
      Kappa10*TKappa10 + Kappa11*TKappa11 + Kappa12*TKappa12 + Kappa20*TKappa20 +
      Kappa21*TKappa21 + Kappa22*TKappa22) + 4*(Lambda1200*TLambda1200 +
      Lambda1201*TLambda1201 + Lambda1210*TLambda1210 + Lambda1211*TLambda1211) +
      4*Lambdax*TLambdax + 0.5333333333333333*MassB*Sqr(g1) + 10.666666666666666*
      MassG*Sqr(g3) + 4*MassBp*Sqr(gN)*Sqr(QDxbarp) + 4*MassBp*Sqr(gN)*Sqr(QDxp) +
      4*MassBp*Sqr(gN)*Sqr(QSp))) + 12*TKappa11*(3*(Kappa01*(Kappa00*TKappa10 +
      Kappa01*TKappa11 + Kappa02*TKappa12) + Kappa11*(Kappa10*TKappa10 + Kappa11*
      TKappa11 + Kappa12*TKappa12) + Kappa21*(Kappa20*TKappa10 + Kappa21*TKappa11
      + Kappa22*TKappa12)) - 0.26666666666666666*TKappa11*Sqr(g1) -
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
      Lambda1210*TLambda1210 + Lambda1211*TLambda1211) + 4*Lambdax*TLambdax +
      0.5333333333333333*MassB*Sqr(g1) + 10.666666666666666*MassG*Sqr(g3) + 4*
      MassBp*Sqr(gN)*Sqr(QDxbarp) + 4*MassBp*Sqr(gN)*Sqr(QDxp) + 4*MassBp*Sqr(gN)*
      Sqr(QSp))) + 12*TKappa12*(3*(Kappa02*(Kappa00*TKappa10 + Kappa01*TKappa11 +
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
      Lambda1211*TLambda1211) + 4*Lambdax*TLambdax + 0.5333333333333333*MassB*Sqr(
      g1) + 10.666666666666666*MassG*Sqr(g3) + 4*MassBp*Sqr(gN)*Sqr(QDxbarp) + 4*
      MassBp*Sqr(gN)*Sqr(QDxp) + 4*MassBp*Sqr(gN)*Sqr(QSp))) + 12*TKappa20*(3*(
      Kappa00*(Kappa00*TKappa20 + Kappa01*TKappa21 + Kappa02*TKappa22) + Kappa10*(
      Kappa10*TKappa20 + Kappa11*TKappa21 + Kappa12*TKappa22) + Kappa20*(Kappa20*
      TKappa20 + Kappa21*TKappa21 + Kappa22*TKappa22)) - 0.26666666666666666*
      TKappa20*Sqr(g1) - 5.333333333333333*TKappa20*Sqr(g3) + 3*TKappa20*(Sqr(
      Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(
      Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 3*((Kappa00*Kappa20
      + Kappa01*Kappa21 + Kappa02*Kappa22)*TKappa00 + (Kappa10*Kappa20 + Kappa11*
      Kappa21 + Kappa12*Kappa22)*TKappa10 + TKappa20*(Sqr(Kappa20) + Sqr(Kappa21)
      + Sqr(Kappa22))) + 2*TKappa20*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(
      Lambda1210) + Sqr(Lambda1211)) + 2*TKappa20*Sqr(Lambdax) - 2*TKappa20*Sqr(gN
      )*Sqr(QDxbarp) - 2*TKappa20*Sqr(gN)*Sqr(QDxp) - 2*TKappa20*Sqr(gN)*Sqr(QSp)
      + Kappa20*(6*(Kappa00*TKappa00 + Kappa01*TKappa01 + Kappa02*TKappa02 +
      Kappa10*TKappa10 + Kappa11*TKappa11 + Kappa12*TKappa12 + Kappa20*TKappa20 +
      Kappa21*TKappa21 + Kappa22*TKappa22) + 4*(Lambda1200*TLambda1200 +
      Lambda1201*TLambda1201 + Lambda1210*TLambda1210 + Lambda1211*TLambda1211) +
      4*Lambdax*TLambdax + 0.5333333333333333*MassB*Sqr(g1) + 10.666666666666666*
      MassG*Sqr(g3) + 4*MassBp*Sqr(gN)*Sqr(QDxbarp) + 4*MassBp*Sqr(gN)*Sqr(QDxp) +
      4*MassBp*Sqr(gN)*Sqr(QSp))) + 12*TKappa21*(3*(Kappa01*(Kappa00*TKappa20 +
      Kappa01*TKappa21 + Kappa02*TKappa22) + Kappa11*(Kappa10*TKappa20 + Kappa11*
      TKappa21 + Kappa12*TKappa22) + Kappa21*(Kappa20*TKappa20 + Kappa21*TKappa21
      + Kappa22*TKappa22)) - 0.26666666666666666*TKappa21*Sqr(g1) -
      5.333333333333333*TKappa21*Sqr(g3) + 3*TKappa21*(Sqr(Kappa00) + Sqr(Kappa01)
      + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20)
      + Sqr(Kappa21) + Sqr(Kappa22)) + 3*((Kappa00*Kappa20 + Kappa01*Kappa21 +
      Kappa02*Kappa22)*TKappa01 + (Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*
      Kappa22)*TKappa11 + TKappa21*(Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22))) +
      2*TKappa21*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(
      Lambda1211)) + 2*TKappa21*Sqr(Lambdax) - 2*TKappa21*Sqr(gN)*Sqr(QDxbarp) - 2
      *TKappa21*Sqr(gN)*Sqr(QDxp) - 2*TKappa21*Sqr(gN)*Sqr(QSp) + Kappa21*(6*(
      Kappa00*TKappa00 + Kappa01*TKappa01 + Kappa02*TKappa02 + Kappa10*TKappa10 +
      Kappa11*TKappa11 + Kappa12*TKappa12 + Kappa20*TKappa20 + Kappa21*TKappa21 +
      Kappa22*TKappa22) + 4*(Lambda1200*TLambda1200 + Lambda1201*TLambda1201 +
      Lambda1210*TLambda1210 + Lambda1211*TLambda1211) + 4*Lambdax*TLambdax +
      0.5333333333333333*MassB*Sqr(g1) + 10.666666666666666*MassG*Sqr(g3) + 4*
      MassBp*Sqr(gN)*Sqr(QDxbarp) + 4*MassBp*Sqr(gN)*Sqr(QDxp) + 4*MassBp*Sqr(gN)*
      Sqr(QSp))) + 12*TKappa22*(3*(Kappa02*(Kappa00*TKappa20 + Kappa01*TKappa21 +
      Kappa02*TKappa22) + Kappa12*(Kappa10*TKappa20 + Kappa11*TKappa21 + Kappa12*
      TKappa22) + Kappa22*(Kappa20*TKappa20 + Kappa21*TKappa21 + Kappa22*TKappa22)
      ) - 0.26666666666666666*TKappa22*Sqr(g1) - 5.333333333333333*TKappa22*Sqr(g3
      ) + 3*TKappa22*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) +
      Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) +
      3*((Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)*TKappa02 + (Kappa10
      *Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)*TKappa12 + TKappa22*(Sqr(
      Kappa20) + Sqr(Kappa21) + Sqr(Kappa22))) + 2*TKappa22*(Sqr(Lambda1200) + Sqr
      (Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 2*TKappa22*Sqr(Lambdax)
      - 2*TKappa22*Sqr(gN)*Sqr(QDxbarp) - 2*TKappa22*Sqr(gN)*Sqr(QDxp) - 2*
      TKappa22*Sqr(gN)*Sqr(QSp) + Kappa22*(6*(Kappa00*TKappa00 + Kappa01*TKappa01
      + Kappa02*TKappa02 + Kappa10*TKappa10 + Kappa11*TKappa11 + Kappa12*TKappa12
      + Kappa20*TKappa20 + Kappa21*TKappa21 + Kappa22*TKappa22) + 4*(Lambda1200*
      TLambda1200 + Lambda1201*TLambda1201 + Lambda1210*TLambda1210 + Lambda1211*
      TLambda1211) + 4*Lambdax*TLambdax + 0.5333333333333333*MassB*Sqr(g1) +
      10.666666666666666*MassG*Sqr(g3) + 4*MassBp*Sqr(gN)*Sqr(QDxbarp) + 4*MassBp*
      Sqr(gN)*Sqr(QDxp) + 4*MassBp*Sqr(gN)*Sqr(QSp))) + 8*TLambda1200*(3*(
      Lambda1200*(Lambda1200*TLambda1200 + Lambda1201*TLambda1201) + Lambda1210*(
      Lambda1210*TLambda1200 + Lambda1211*TLambda1201)) - 0.6*TLambda1200*Sqr(g1)
      - 3*TLambda1200*Sqr(g2) + 3*TLambda1200*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(
      Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(
      Kappa21) + Sqr(Kappa22)) + 3*((Lambda1200*Lambda1210 + Lambda1201*Lambda1211
      )*TLambda1210 + TLambda1200*(Sqr(Lambda1200) + Sqr(Lambda1201))) + 2*
      TLambda1200*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(
      Lambda1211)) + 2*TLambda1200*Sqr(Lambdax) - 2*TLambda1200*Sqr(gN)*Sqr(QH1p)
      - 2*TLambda1200*Sqr(gN)*Sqr(QH2p) - 2*TLambda1200*Sqr(gN)*Sqr(QSp) +
      Lambda1200*(6*(Kappa00*TKappa00 + Kappa01*TKappa01 + Kappa02*TKappa02 +
      Kappa10*TKappa10 + Kappa11*TKappa11 + Kappa12*TKappa12 + Kappa20*TKappa20 +
      Kappa21*TKappa21 + Kappa22*TKappa22) + 4*(Lambda1200*TLambda1200 +
      Lambda1201*TLambda1201 + Lambda1210*TLambda1210 + Lambda1211*TLambda1211) +
      4*Lambdax*TLambdax + 1.2*MassB*Sqr(g1) + 6*MassWB*Sqr(g2) + 4*MassBp*Sqr(gN)
      *Sqr(QH1p) + 4*MassBp*Sqr(gN)*Sqr(QH2p) + 4*MassBp*Sqr(gN)*Sqr(QSp))) + 8*
      TLambda1201*(3*(Lambda1201*(Lambda1200*TLambda1200 + Lambda1201*TLambda1201)
      + Lambda1211*(Lambda1210*TLambda1200 + Lambda1211*TLambda1201)) - 0.6*
      TLambda1201*Sqr(g1) - 3*TLambda1201*Sqr(g2) + 3*TLambda1201*(Sqr(Kappa00) +
      Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) +
      Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 3*((Lambda1200*Lambda1210 +
      Lambda1201*Lambda1211)*TLambda1211 + TLambda1201*(Sqr(Lambda1200) + Sqr(
      Lambda1201))) + 2*TLambda1201*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(
      Lambda1210) + Sqr(Lambda1211)) + 2*TLambda1201*Sqr(Lambdax) - 2*TLambda1201*
      Sqr(gN)*Sqr(QH1p) - 2*TLambda1201*Sqr(gN)*Sqr(QH2p) - 2*TLambda1201*Sqr(gN)*
      Sqr(QSp) + Lambda1201*(6*(Kappa00*TKappa00 + Kappa01*TKappa01 + Kappa02*
      TKappa02 + Kappa10*TKappa10 + Kappa11*TKappa11 + Kappa12*TKappa12 + Kappa20*
      TKappa20 + Kappa21*TKappa21 + Kappa22*TKappa22) + 4*(Lambda1200*TLambda1200
      + Lambda1201*TLambda1201 + Lambda1210*TLambda1210 + Lambda1211*TLambda1211)
      + 4*Lambdax*TLambdax + 1.2*MassB*Sqr(g1) + 6*MassWB*Sqr(g2) + 4*MassBp*Sqr(
      gN)*Sqr(QH1p) + 4*MassBp*Sqr(gN)*Sqr(QH2p) + 4*MassBp*Sqr(gN)*Sqr(QSp))) + 8
      *TLambda1210*(3*(Lambda1200*(Lambda1200*TLambda1210 + Lambda1201*TLambda1211
      ) + Lambda1210*(Lambda1210*TLambda1210 + Lambda1211*TLambda1211)) - 0.6*
      TLambda1210*Sqr(g1) - 3*TLambda1210*Sqr(g2) + 3*TLambda1210*(Sqr(Kappa00) +
      Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) +
      Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*TLambda1210*(Sqr(Lambda1200)
      + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 3*((Lambda1200*
      Lambda1210 + Lambda1201*Lambda1211)*TLambda1200 + TLambda1210*(Sqr(
      Lambda1210) + Sqr(Lambda1211))) + 2*TLambda1210*Sqr(Lambdax) - 2*TLambda1210
      *Sqr(gN)*Sqr(QH1p) - 2*TLambda1210*Sqr(gN)*Sqr(QH2p) - 2*TLambda1210*Sqr(gN)
      *Sqr(QSp) + Lambda1210*(6*(Kappa00*TKappa00 + Kappa01*TKappa01 + Kappa02*
      TKappa02 + Kappa10*TKappa10 + Kappa11*TKappa11 + Kappa12*TKappa12 + Kappa20*
      TKappa20 + Kappa21*TKappa21 + Kappa22*TKappa22) + 4*(Lambda1200*TLambda1200
      + Lambda1201*TLambda1201 + Lambda1210*TLambda1210 + Lambda1211*TLambda1211)
      + 4*Lambdax*TLambdax + 1.2*MassB*Sqr(g1) + 6*MassWB*Sqr(g2) + 4*MassBp*Sqr(
      gN)*Sqr(QH1p) + 4*MassBp*Sqr(gN)*Sqr(QH2p) + 4*MassBp*Sqr(gN)*Sqr(QSp))) + 8
      *TLambda1211*(3*(Lambda1201*(Lambda1200*TLambda1210 + Lambda1201*TLambda1211
      ) + Lambda1211*(Lambda1210*TLambda1210 + Lambda1211*TLambda1211)) - 0.6*
      TLambda1211*Sqr(g1) - 3*TLambda1211*Sqr(g2) + 3*TLambda1211*(Sqr(Kappa00) +
      Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) +
      Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*TLambda1211*(Sqr(Lambda1200)
      + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 3*((Lambda1200*
      Lambda1210 + Lambda1201*Lambda1211)*TLambda1201 + TLambda1211*(Sqr(
      Lambda1210) + Sqr(Lambda1211))) + 2*TLambda1211*Sqr(Lambdax) - 2*TLambda1211
      *Sqr(gN)*Sqr(QH1p) - 2*TLambda1211*Sqr(gN)*Sqr(QH2p) - 2*TLambda1211*Sqr(gN)
      *Sqr(QSp) + Lambda1211*(6*(Kappa00*TKappa00 + Kappa01*TKappa01 + Kappa02*
      TKappa02 + Kappa10*TKappa10 + Kappa11*TKappa11 + Kappa12*TKappa12 + Kappa20*
      TKappa20 + Kappa21*TKappa21 + Kappa22*TKappa22) + 4*(Lambda1200*TLambda1200
      + Lambda1201*TLambda1201 + Lambda1210*TLambda1210 + Lambda1211*TLambda1211)
      + 4*Lambdax*TLambdax + 1.2*MassB*Sqr(g1) + 6*MassWB*Sqr(g2) + 4*MassBp*Sqr(
      gN)*Sqr(QH1p) + 4*MassBp*Sqr(gN)*Sqr(QH2p) + 4*MassBp*Sqr(gN)*Sqr(QSp))) -
      32*Power(gN,4)*Sqr(MassBp)*Sqr(QSp)*(9*Sqr(Qdp) + 9*Sqr(QDxbarp) + 9*Sqr(
      QDxp) + 3*Sqr(Qep) + 6*Sqr(QH1p) + 6*Sqr(QH2p) + 2*Sqr(QHpbarp) + 2*Sqr(QHpp
      ) + 6*Sqr(QLp) + 18*Sqr(QQp) + 3*Sqr(QSp) + 9*Sqr(Qup)) + Power(gN,3)*(4*gN*
      QSp*(3*(md200 + md211 + md222)*Qdp + 3*(mDxbar200 + mDxbar211 + mDxbar222)*
      QDxbarp + 3*(mDx200 + mDx211 + mDx222)*QDxp + (me200 + me211 + me222)*Qep +
      2*(mH1I200 + mH1I211)*QH1p + 2*mHd2*QH1p + 2*(mH2I200 + mH2I211)*QH2p + 2*
      mHu2*QH2p + 2*mHpbar2*QHpbarp + 2*mHp2*QHpp + 2*(ml200 + ml211 + ml222)*QLp
      + 6*(mq200 + mq211 + mq222)*QQp + ms2*QSp + (msI200 + msI211)*QSp + 3*(mu200
      + mu211 + mu222)*Qup) - 16*gN*Sqr(MassBp)*Sqr(QSp))*(9*Sqr(Qdp) + 9*Sqr(
      QDxbarp) + 9*Sqr(QDxp) + 3*Sqr(Qep) + 6*Sqr(QH1p) + 6*Sqr(QH2p) + 2*Sqr(
      QHpbarp) + 2*Sqr(QHpp) + 6*Sqr(QLp) + 18*Sqr(QQp) + 3*Sqr(QSp) + 9*Sqr(Qup))
      + (6*QDxp*QSp*Sqr(gN) + 6*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02)))*(
      Kappa00*(Kappa00*mDx200 + Kappa10*mDx201 + Kappa20*mDx202) + Kappa01*(
      Kappa01*mDx200 + Kappa11*mDx201 + Kappa21*mDx202) + Kappa02*(Kappa02*mDx200
      + Kappa12*mDx201 + Kappa22*mDx202) + (Kappa00*Kappa10 + Kappa01*Kappa11 +
      Kappa02*Kappa12)*mDx210 + (Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*
      Kappa22)*mDx220 + 2*(Kappa00*(Kappa00*mDxbar200 + Kappa01*mDxbar210 +
      Kappa02*mDxbar220) + Kappa01*(Kappa00*mDxbar201 + Kappa01*mDxbar211 +
      Kappa02*mDxbar221) + Kappa02*(Kappa00*mDxbar202 + Kappa01*mDxbar212 +
      Kappa02*mDxbar222)) - 0.4*(md200 + md211 + md222 - mDx200 - mDx211 - mDx222
      + mDxbar200 + mDxbar211 + mDxbar222 + me200 + me211 + me222 - mH1I200 -
      mH1I211 + mH2I200 + mH2I211 - mHd2 - mHp2 + mHpbar2 + mHu2 - ml200 - ml211 -
      ml222 + mq200 + mq211 + mq222 - 2*(mu200 + mu211 + mu222))*Sqr(g1) + 2*QDxp
      *(3*(md200 + md211 + md222)*Qdp + 3*(mDxbar200 + mDxbar211 + mDxbar222)*
      QDxbarp + 3*(mDx200 + mDx211 + mDx222)*QDxp + (me200 + me211 + me222)*Qep +
      2*(mH1I200 + mH1I211)*QH1p + 2*mHd2*QH1p + 2*(mH2I200 + mH2I211)*QH2p + 2*
      mHu2*QH2p + 2*mHpbar2*QHpbarp + 2*mHp2*QHpp + 2*(ml200 + ml211 + ml222)*QLp
      + 6*(mq200 + mq211 + mq222)*QQp + ms2*QSp + (msI200 + msI211)*QSp + 3*(mu200
      + mu211 + mu222)*Qup)*Sqr(gN) + mDx200*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(
      Kappa02)) + 2*ms2*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02)) -
      0.5333333333333333*Sqr(g1)*Sqr(MassB) - 10.666666666666666*Sqr(g3)*Sqr(MassG
      ) - 8*Sqr(gN)*Sqr(MassBp)*Sqr(QDxp) + 2*(Sqr(TKappa00) + Sqr(TKappa01) + Sqr
      (TKappa02))) + (6*QDxp*QSp*Sqr(gN) + 6*(Sqr(Kappa10) + Sqr(Kappa11) + Sqr(
      Kappa12)))*((Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12)*mDx201 +
      Kappa10*(Kappa00*mDx210 + Kappa10*mDx211 + Kappa20*mDx212) + Kappa11*(
      Kappa01*mDx210 + Kappa11*mDx211 + Kappa21*mDx212) + Kappa12*(Kappa02*mDx210
      + Kappa12*mDx211 + Kappa22*mDx212) + (Kappa10*Kappa20 + Kappa11*Kappa21 +
      Kappa12*Kappa22)*mDx221 + 2*(Kappa10*(Kappa10*mDxbar200 + Kappa11*mDxbar210
      + Kappa12*mDxbar220) + Kappa11*(Kappa10*mDxbar201 + Kappa11*mDxbar211 +
      Kappa12*mDxbar221) + Kappa12*(Kappa10*mDxbar202 + Kappa11*mDxbar212 +
      Kappa12*mDxbar222)) - 0.4*(md200 + md211 + md222 - mDx200 - mDx211 - mDx222
      + mDxbar200 + mDxbar211 + mDxbar222 + me200 + me211 + me222 - mH1I200 -
      mH1I211 + mH2I200 + mH2I211 - mHd2 - mHp2 + mHpbar2 + mHu2 - ml200 - ml211 -
      ml222 + mq200 + mq211 + mq222 - 2*(mu200 + mu211 + mu222))*Sqr(g1) + 2*QDxp
      *(3*(md200 + md211 + md222)*Qdp + 3*(mDxbar200 + mDxbar211 + mDxbar222)*
      QDxbarp + 3*(mDx200 + mDx211 + mDx222)*QDxp + (me200 + me211 + me222)*Qep +
      2*(mH1I200 + mH1I211)*QH1p + 2*mHd2*QH1p + 2*(mH2I200 + mH2I211)*QH2p + 2*
      mHu2*QH2p + 2*mHpbar2*QHpbarp + 2*mHp2*QHpp + 2*(ml200 + ml211 + ml222)*QLp
      + 6*(mq200 + mq211 + mq222)*QQp + ms2*QSp + (msI200 + msI211)*QSp + 3*(mu200
      + mu211 + mu222)*Qup)*Sqr(gN) + mDx211*(Sqr(Kappa10) + Sqr(Kappa11) + Sqr(
      Kappa12)) + 2*ms2*(Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12)) -
      0.5333333333333333*Sqr(g1)*Sqr(MassB) - 10.666666666666666*Sqr(g3)*Sqr(MassG
      ) - 8*Sqr(gN)*Sqr(MassBp)*Sqr(QDxp) + 2*(Sqr(TKappa10) + Sqr(TKappa11) + Sqr
      (TKappa12))) + (6*QDxbarp*QSp*Sqr(gN) + 6*(Sqr(Kappa00) + Sqr(Kappa10) + Sqr
      (Kappa20)))*(2*(Kappa00*(Kappa00*mDx200 + Kappa10*mDx210 + Kappa20*mDx220) +
      Kappa10*(Kappa00*mDx201 + Kappa10*mDx211 + Kappa20*mDx221) + Kappa20*(
      Kappa00*mDx202 + Kappa10*mDx212 + Kappa20*mDx222)) + Kappa00*(Kappa00*
      mDxbar200 + Kappa01*mDxbar201 + Kappa02*mDxbar202) + Kappa10*(Kappa10*
      mDxbar200 + Kappa11*mDxbar201 + Kappa12*mDxbar202) + Kappa20*(Kappa20*
      mDxbar200 + Kappa21*mDxbar201 + Kappa22*mDxbar202) + (Kappa00*Kappa01 +
      Kappa10*Kappa11 + Kappa20*Kappa21)*mDxbar210 + (Kappa00*Kappa02 + Kappa10*
      Kappa12 + Kappa20*Kappa22)*mDxbar220 + 0.4*(md200 + md211 + md222 - mDx200 -
      mDx211 - mDx222 + mDxbar200 + mDxbar211 + mDxbar222 + me200 + me211 + me222
      - mH1I200 - mH1I211 + mH2I200 + mH2I211 - mHd2 - mHp2 + mHpbar2 + mHu2 -
      ml200 - ml211 - ml222 + mq200 + mq211 + mq222 - 2*(mu200 + mu211 + mu222))*
      Sqr(g1) + 2*QDxbarp*(3*(md200 + md211 + md222)*Qdp + 3*(mDxbar200 +
      mDxbar211 + mDxbar222)*QDxbarp + 3*(mDx200 + mDx211 + mDx222)*QDxp + (me200
      + me211 + me222)*Qep + 2*(mH1I200 + mH1I211)*QH1p + 2*mHd2*QH1p + 2*(mH2I200
      + mH2I211)*QH2p + 2*mHu2*QH2p + 2*mHpbar2*QHpbarp + 2*mHp2*QHpp + 2*(ml200
      + ml211 + ml222)*QLp + 6*(mq200 + mq211 + mq222)*QQp + ms2*QSp + (msI200 +
      msI211)*QSp + 3*(mu200 + mu211 + mu222)*Qup)*Sqr(gN) + mDxbar200*(Sqr(
      Kappa00) + Sqr(Kappa10) + Sqr(Kappa20)) + 2*ms2*(Sqr(Kappa00) + Sqr(Kappa10)
      + Sqr(Kappa20)) - 0.5333333333333333*Sqr(g1)*Sqr(MassB) -
      10.666666666666666*Sqr(g3)*Sqr(MassG) - 8*Sqr(gN)*Sqr(MassBp)*Sqr(QDxbarp) +
      2*(Sqr(TKappa00) + Sqr(TKappa10) + Sqr(TKappa20))) + (6*QDxbarp*QSp*Sqr(gN)
      + 6*(Sqr(Kappa01) + Sqr(Kappa11) + Sqr(Kappa21)))*(2*(Kappa01*(Kappa01*
      mDx200 + Kappa11*mDx210 + Kappa21*mDx220) + Kappa11*(Kappa01*mDx201 +
      Kappa11*mDx211 + Kappa21*mDx221) + Kappa21*(Kappa01*mDx202 + Kappa11*mDx212
      + Kappa21*mDx222)) + (Kappa00*Kappa01 + Kappa10*Kappa11 + Kappa20*Kappa21)*
      mDxbar201 + Kappa01*(Kappa00*mDxbar210 + Kappa01*mDxbar211 + Kappa02*
      mDxbar212) + Kappa11*(Kappa10*mDxbar210 + Kappa11*mDxbar211 + Kappa12*
      mDxbar212) + Kappa21*(Kappa20*mDxbar210 + Kappa21*mDxbar211 + Kappa22*
      mDxbar212) + (Kappa01*Kappa02 + Kappa11*Kappa12 + Kappa21*Kappa22)*mDxbar221
      + 0.4*(md200 + md211 + md222 - mDx200 - mDx211 - mDx222 + mDxbar200 +
      mDxbar211 + mDxbar222 + me200 + me211 + me222 - mH1I200 - mH1I211 + mH2I200
      + mH2I211 - mHd2 - mHp2 + mHpbar2 + mHu2 - ml200 - ml211 - ml222 + mq200 +
      mq211 + mq222 - 2*(mu200 + mu211 + mu222))*Sqr(g1) + 2*QDxbarp*(3*(md200 +
      md211 + md222)*Qdp + 3*(mDxbar200 + mDxbar211 + mDxbar222)*QDxbarp + 3*(
      mDx200 + mDx211 + mDx222)*QDxp + (me200 + me211 + me222)*Qep + 2*(mH1I200 +
      mH1I211)*QH1p + 2*mHd2*QH1p + 2*(mH2I200 + mH2I211)*QH2p + 2*mHu2*QH2p + 2*
      mHpbar2*QHpbarp + 2*mHp2*QHpp + 2*(ml200 + ml211 + ml222)*QLp + 6*(mq200 +
      mq211 + mq222)*QQp + ms2*QSp + (msI200 + msI211)*QSp + 3*(mu200 + mu211 +
      mu222)*Qup)*Sqr(gN) + mDxbar211*(Sqr(Kappa01) + Sqr(Kappa11) + Sqr(Kappa21))
      + 2*ms2*(Sqr(Kappa01) + Sqr(Kappa11) + Sqr(Kappa21)) - 0.5333333333333333*
      Sqr(g1)*Sqr(MassB) - 10.666666666666666*Sqr(g3)*Sqr(MassG) - 8*Sqr(gN)*Sqr(
      MassBp)*Sqr(QDxbarp) + 2*(Sqr(TKappa01) + Sqr(TKappa11) + Sqr(TKappa21))) +
      (6*QDxbarp*QSp*Sqr(gN) + 6*(Sqr(Kappa02) + Sqr(Kappa12) + Sqr(Kappa22)))*(2*
      (Kappa02*(Kappa02*mDx200 + Kappa12*mDx210 + Kappa22*mDx220) + Kappa12*(
      Kappa02*mDx201 + Kappa12*mDx211 + Kappa22*mDx221) + Kappa22*(Kappa02*mDx202
      + Kappa12*mDx212 + Kappa22*mDx222)) + (Kappa00*Kappa02 + Kappa10*Kappa12 +
      Kappa20*Kappa22)*mDxbar202 + (Kappa01*Kappa02 + Kappa11*Kappa12 + Kappa21*
      Kappa22)*mDxbar212 + Kappa02*(Kappa00*mDxbar220 + Kappa01*mDxbar221 +
      Kappa02*mDxbar222) + Kappa12*(Kappa10*mDxbar220 + Kappa11*mDxbar221 +
      Kappa12*mDxbar222) + Kappa22*(Kappa20*mDxbar220 + Kappa21*mDxbar221 +
      Kappa22*mDxbar222) + 0.4*(md200 + md211 + md222 - mDx200 - mDx211 - mDx222 +
      mDxbar200 + mDxbar211 + mDxbar222 + me200 + me211 + me222 - mH1I200 -
      mH1I211 + mH2I200 + mH2I211 - mHd2 - mHp2 + mHpbar2 + mHu2 - ml200 - ml211 -
      ml222 + mq200 + mq211 + mq222 - 2*(mu200 + mu211 + mu222))*Sqr(g1) + 2*
      QDxbarp*(3*(md200 + md211 + md222)*Qdp + 3*(mDxbar200 + mDxbar211 +
      mDxbar222)*QDxbarp + 3*(mDx200 + mDx211 + mDx222)*QDxp + (me200 + me211 +
      me222)*Qep + 2*(mH1I200 + mH1I211)*QH1p + 2*mHd2*QH1p + 2*(mH2I200 + mH2I211
      )*QH2p + 2*mHu2*QH2p + 2*mHpbar2*QHpbarp + 2*mHp2*QHpp + 2*(ml200 + ml211 +
      ml222)*QLp + 6*(mq200 + mq211 + mq222)*QQp + ms2*QSp + (msI200 + msI211)*QSp
      + 3*(mu200 + mu211 + mu222)*Qup)*Sqr(gN) + mDxbar222*(Sqr(Kappa02) + Sqr(
      Kappa12) + Sqr(Kappa22)) + 2*ms2*(Sqr(Kappa02) + Sqr(Kappa12) + Sqr(Kappa22)
      ) - 0.5333333333333333*Sqr(g1)*Sqr(MassB) - 10.666666666666666*Sqr(g3)*Sqr(
      MassG) - 8*Sqr(gN)*Sqr(MassBp)*Sqr(QDxbarp) + 2*(Sqr(TKappa02) + Sqr(
      TKappa12) + Sqr(TKappa22))) + (6*QDxp*QSp*Sqr(gN) + 6*(Sqr(Kappa20) + Sqr(
      Kappa21) + Sqr(Kappa22)))*((Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*
      Kappa22)*mDx202 + (Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)*
      mDx212 + Kappa20*(Kappa00*mDx220 + Kappa10*mDx221 + Kappa20*mDx222) +
      Kappa21*(Kappa01*mDx220 + Kappa11*mDx221 + Kappa21*mDx222) + Kappa22*(
      Kappa02*mDx220 + Kappa12*mDx221 + Kappa22*mDx222) + 2*(Kappa20*(Kappa20*
      mDxbar200 + Kappa21*mDxbar210 + Kappa22*mDxbar220) + Kappa21*(Kappa20*
      mDxbar201 + Kappa21*mDxbar211 + Kappa22*mDxbar221) + Kappa22*(Kappa20*
      mDxbar202 + Kappa21*mDxbar212 + Kappa22*mDxbar222)) - 0.4*(md200 + md211 +
      md222 - mDx200 - mDx211 - mDx222 + mDxbar200 + mDxbar211 + mDxbar222 + me200
      + me211 + me222 - mH1I200 - mH1I211 + mH2I200 + mH2I211 - mHd2 - mHp2 +
      mHpbar2 + mHu2 - ml200 - ml211 - ml222 + mq200 + mq211 + mq222 - 2*(mu200 +
      mu211 + mu222))*Sqr(g1) + 2*QDxp*(3*(md200 + md211 + md222)*Qdp + 3*(
      mDxbar200 + mDxbar211 + mDxbar222)*QDxbarp + 3*(mDx200 + mDx211 + mDx222)*
      QDxp + (me200 + me211 + me222)*Qep + 2*(mH1I200 + mH1I211)*QH1p + 2*mHd2*
      QH1p + 2*(mH2I200 + mH2I211)*QH2p + 2*mHu2*QH2p + 2*mHpbar2*QHpbarp + 2*mHp2
      *QHpp + 2*(ml200 + ml211 + ml222)*QLp + 6*(mq200 + mq211 + mq222)*QQp + ms2*
      QSp + (msI200 + msI211)*QSp + 3*(mu200 + mu211 + mu222)*Qup)*Sqr(gN) +
      mDx222*(Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*ms2*(Sqr(Kappa20) +
      Sqr(Kappa21) + Sqr(Kappa22)) - 0.5333333333333333*Sqr(g1)*Sqr(MassB) -
      10.666666666666666*Sqr(g3)*Sqr(MassG) - 8*Sqr(gN)*Sqr(MassBp)*Sqr(QDxp) + 2*
      (Sqr(TKappa20) + Sqr(TKappa21) + Sqr(TKappa22))) + (4*QH2p*QSp*Sqr(gN) + 4*(
      Sqr(Lambda1200) + Sqr(Lambda1201)))*(2*(Lambda1200*(Lambda1200*mH1I200 +
      Lambda1201*mH1I210) + Lambda1201*(Lambda1200*mH1I201 + Lambda1201*mH1I211))
      + Lambda1200*(Lambda1200*mH2I200 + Lambda1210*mH2I201) + Lambda1201*(
      Lambda1201*mH2I200 + Lambda1211*mH2I201) + (Lambda1200*Lambda1210 +
      Lambda1201*Lambda1211)*mH2I210 + 0.6*(md200 + md211 + md222 - mDx200 -
      mDx211 - mDx222 + mDxbar200 + mDxbar211 + mDxbar222 + me200 + me211 + me222
      - mH1I200 - mH1I211 + mH2I200 + mH2I211 - mHd2 - mHp2 + mHpbar2 + mHu2 -
      ml200 - ml211 - ml222 + mq200 + mq211 + mq222 - 2*(mu200 + mu211 + mu222))*
      Sqr(g1) + 2*QH2p*(3*(md200 + md211 + md222)*Qdp + 3*(mDxbar200 + mDxbar211 +
      mDxbar222)*QDxbarp + 3*(mDx200 + mDx211 + mDx222)*QDxp + (me200 + me211 +
      me222)*Qep + 2*(mH1I200 + mH1I211)*QH1p + 2*mHd2*QH1p + 2*(mH2I200 + mH2I211
      )*QH2p + 2*mHu2*QH2p + 2*mHpbar2*QHpbarp + 2*mHp2*QHpp + 2*(ml200 + ml211 +
      ml222)*QLp + 6*(mq200 + mq211 + mq222)*QQp + ms2*QSp + (msI200 + msI211)*QSp
      + 3*(mu200 + mu211 + mu222)*Qup)*Sqr(gN) + mH2I200*(Sqr(Lambda1200) + Sqr(
      Lambda1201)) + 2*ms2*(Sqr(Lambda1200) + Sqr(Lambda1201)) - 1.2*Sqr(g1)*Sqr(
      MassB) - 6*Sqr(g2)*Sqr(MassWB) - 8*Sqr(gN)*Sqr(MassBp)*Sqr(QH2p) + 2*(Sqr(
      TLambda1200) + Sqr(TLambda1201))) + (4*QH1p*QSp*Sqr(gN) + 4*(Sqr(Lambda1200)
      + Sqr(Lambda1210)))*(Lambda1200*(Lambda1200*mH1I200 + Lambda1201*mH1I201) +
      Lambda1210*(Lambda1210*mH1I200 + Lambda1211*mH1I201) + (Lambda1200*
      Lambda1201 + Lambda1210*Lambda1211)*mH1I210 + 2*(Lambda1200*(Lambda1200*
      mH2I200 + Lambda1210*mH2I210) + Lambda1210*(Lambda1200*mH2I201 + Lambda1210*
      mH2I211)) - 0.6*(md200 + md211 + md222 - mDx200 - mDx211 - mDx222 +
      mDxbar200 + mDxbar211 + mDxbar222 + me200 + me211 + me222 - mH1I200 -
      mH1I211 + mH2I200 + mH2I211 - mHd2 - mHp2 + mHpbar2 + mHu2 - ml200 - ml211 -
      ml222 + mq200 + mq211 + mq222 - 2*(mu200 + mu211 + mu222))*Sqr(g1) + 2*QH1p
      *(3*(md200 + md211 + md222)*Qdp + 3*(mDxbar200 + mDxbar211 + mDxbar222)*
      QDxbarp + 3*(mDx200 + mDx211 + mDx222)*QDxp + (me200 + me211 + me222)*Qep +
      2*(mH1I200 + mH1I211)*QH1p + 2*mHd2*QH1p + 2*(mH2I200 + mH2I211)*QH2p + 2*
      mHu2*QH2p + 2*mHpbar2*QHpbarp + 2*mHp2*QHpp + 2*(ml200 + ml211 + ml222)*QLp
      + 6*(mq200 + mq211 + mq222)*QQp + ms2*QSp + (msI200 + msI211)*QSp + 3*(mu200
      + mu211 + mu222)*Qup)*Sqr(gN) + mH1I200*(Sqr(Lambda1200) + Sqr(Lambda1210))
      + 2*ms2*(Sqr(Lambda1200) + Sqr(Lambda1210)) - 1.2*Sqr(g1)*Sqr(MassB) - 6*
      Sqr(g2)*Sqr(MassWB) - 8*Sqr(gN)*Sqr(MassBp)*Sqr(QH1p) + 2*(Sqr(TLambda1200)
      + Sqr(TLambda1210))) + (4*QH1p*QSp*Sqr(gN) + 4*(Sqr(Lambda1201) + Sqr(
      Lambda1211)))*((Lambda1200*Lambda1201 + Lambda1210*Lambda1211)*mH1I201 +
      Lambda1201*(Lambda1200*mH1I210 + Lambda1201*mH1I211) + Lambda1211*(
      Lambda1210*mH1I210 + Lambda1211*mH1I211) + 2*(Lambda1201*(Lambda1201*mH2I200
      + Lambda1211*mH2I210) + Lambda1211*(Lambda1201*mH2I201 + Lambda1211*mH2I211
      )) - 0.6*(md200 + md211 + md222 - mDx200 - mDx211 - mDx222 + mDxbar200 +
      mDxbar211 + mDxbar222 + me200 + me211 + me222 - mH1I200 - mH1I211 + mH2I200
      + mH2I211 - mHd2 - mHp2 + mHpbar2 + mHu2 - ml200 - ml211 - ml222 + mq200 +
      mq211 + mq222 - 2*(mu200 + mu211 + mu222))*Sqr(g1) + 2*QH1p*(3*(md200 +
      md211 + md222)*Qdp + 3*(mDxbar200 + mDxbar211 + mDxbar222)*QDxbarp + 3*(
      mDx200 + mDx211 + mDx222)*QDxp + (me200 + me211 + me222)*Qep + 2*(mH1I200 +
      mH1I211)*QH1p + 2*mHd2*QH1p + 2*(mH2I200 + mH2I211)*QH2p + 2*mHu2*QH2p + 2*
      mHpbar2*QHpbarp + 2*mHp2*QHpp + 2*(ml200 + ml211 + ml222)*QLp + 6*(mq200 +
      mq211 + mq222)*QQp + ms2*QSp + (msI200 + msI211)*QSp + 3*(mu200 + mu211 +
      mu222)*Qup)*Sqr(gN) + mH1I211*(Sqr(Lambda1201) + Sqr(Lambda1211)) + 2*ms2*(
      Sqr(Lambda1201) + Sqr(Lambda1211)) - 1.2*Sqr(g1)*Sqr(MassB) - 6*Sqr(g2)*Sqr(
      MassWB) - 8*Sqr(gN)*Sqr(MassBp)*Sqr(QH1p) + 2*(Sqr(TLambda1201) + Sqr(
      TLambda1211))) + (4*QH2p*QSp*Sqr(gN) + 4*(Sqr(Lambda1210) + Sqr(Lambda1211))
      )*(2*(Lambda1210*(Lambda1210*mH1I200 + Lambda1211*mH1I210) + Lambda1211*(
      Lambda1210*mH1I201 + Lambda1211*mH1I211)) + (Lambda1200*Lambda1210 +
      Lambda1201*Lambda1211)*mH2I201 + Lambda1210*(Lambda1200*mH2I210 + Lambda1210
      *mH2I211) + Lambda1211*(Lambda1201*mH2I210 + Lambda1211*mH2I211) + 0.6*(
      md200 + md211 + md222 - mDx200 - mDx211 - mDx222 + mDxbar200 + mDxbar211 +
      mDxbar222 + me200 + me211 + me222 - mH1I200 - mH1I211 + mH2I200 + mH2I211 -
      mHd2 - mHp2 + mHpbar2 + mHu2 - ml200 - ml211 - ml222 + mq200 + mq211 + mq222
      - 2*(mu200 + mu211 + mu222))*Sqr(g1) + 2*QH2p*(3*(md200 + md211 + md222)*
      Qdp + 3*(mDxbar200 + mDxbar211 + mDxbar222)*QDxbarp + 3*(mDx200 + mDx211 +
      mDx222)*QDxp + (me200 + me211 + me222)*Qep + 2*(mH1I200 + mH1I211)*QH1p + 2*
      mHd2*QH1p + 2*(mH2I200 + mH2I211)*QH2p + 2*mHu2*QH2p + 2*mHpbar2*QHpbarp + 2
      *mHp2*QHpp + 2*(ml200 + ml211 + ml222)*QLp + 6*(mq200 + mq211 + mq222)*QQp +
      ms2*QSp + (msI200 + msI211)*QSp + 3*(mu200 + mu211 + mu222)*Qup)*Sqr(gN) +
      mH2I211*(Sqr(Lambda1210) + Sqr(Lambda1211)) + 2*ms2*(Sqr(Lambda1210) + Sqr(
      Lambda1211)) - 1.2*Sqr(g1)*Sqr(MassB) - 6*Sqr(g2)*Sqr(MassWB) - 8*Sqr(gN)*
      Sqr(MassBp)*Sqr(QH2p) + 2*(Sqr(TLambda1210) + Sqr(TLambda1211))) + (6*(Sqr(
      Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(
      Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 4*(Sqr(Lambda1200)
      + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 4*Sqr(Lambdax) + 2*
      Sqr(gN)*Sqr(QSp))*(6*(Kappa00*(Kappa00*mDxbar200 + Kappa01*mDxbar210 +
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
      Lambda1211))) + 4*(mHd2 + mHu2 + ms2)*Sqr(Lambdax) - 8*Sqr(gN)*Sqr(MassBp)*
      Sqr(QSp) + 6*(Sqr(TKappa00) + Sqr(TKappa01) + Sqr(TKappa02) + Sqr(TKappa10)
      + Sqr(TKappa11) + Sqr(TKappa12) + Sqr(TKappa20) + Sqr(TKappa21) + Sqr(
      TKappa22)) + 4*(Sqr(TLambda1200) + Sqr(TLambda1201) + Sqr(TLambda1210) + Sqr
      (TLambda1211)) + 4*Sqr(TLambdax)) + 6*Qdp*QSp*Sqr(gN)*(4*(Yd00*(mq200*Yd00 +
      mq210*Yd01 + mq220*Yd02) + Yd01*(mq201*Yd00 + mq211*Yd01 + mq221*Yd02) +
      Yd02*(mq202*Yd00 + mq212*Yd01 + mq222*Yd02)) + 2*(Yd00*(md200*Yd00 + md201*
      Yd10 + md202*Yd20) + Yd01*(md200*Yd01 + md201*Yd11 + md202*Yd21) + Yd02*(
      md200*Yd02 + md201*Yd12 + md202*Yd22)) + 0.4*(md200 + md211 + md222 - mDx200
      - mDx211 - mDx222 + mDxbar200 + mDxbar211 + mDxbar222 + me200 + me211 +
      me222 - mH1I200 - mH1I211 + mH2I200 + mH2I211 - mHd2 - mHp2 + mHpbar2 + mHu2
      - ml200 - ml211 - ml222 + mq200 + mq211 + mq222 - 2*(mu200 + mu211 + mu222)
      )*Sqr(g1) + 2*Qdp*(3*(md200 + md211 + md222)*Qdp + 3*(mDxbar200 + mDxbar211
      + mDxbar222)*QDxbarp + 3*(mDx200 + mDx211 + mDx222)*QDxp + (me200 + me211 +
      me222)*Qep + 2*(mH1I200 + mH1I211)*QH1p + 2*mHd2*QH1p + 2*(mH2I200 + mH2I211
      )*QH2p + 2*mHu2*QH2p + 2*mHpbar2*QHpbarp + 2*mHp2*QHpp + 2*(ml200 + ml211 +
      ml222)*QLp + 6*(mq200 + mq211 + mq222)*QQp + ms2*QSp + (msI200 + msI211)*QSp
      + 3*(mu200 + mu211 + mu222)*Qup)*Sqr(gN) - 0.5333333333333333*Sqr(g1)*Sqr(
      MassB) - 10.666666666666666*Sqr(g3)*Sqr(MassG) - 8*Sqr(gN)*Sqr(MassBp)*Sqr(
      Qdp) + 4*(Sqr(TYd00) + Sqr(TYd01) + Sqr(TYd02)) + 4*mHd2*(Sqr(Yd00) + Sqr(
      Yd01) + Sqr(Yd02)) + 2*(md210*(Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + md220*(
      Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22) + md200*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02
      )))) + 6*Qdp*QSp*Sqr(gN)*(4*(Yd10*(mq200*Yd10 + mq210*Yd11 + mq220*Yd12) +
      Yd11*(mq201*Yd10 + mq211*Yd11 + mq221*Yd12) + Yd12*(mq202*Yd10 + mq212*Yd11
      + mq222*Yd12)) + 2*(Yd10*(md210*Yd00 + md211*Yd10 + md212*Yd20) + Yd11*(
      md210*Yd01 + md211*Yd11 + md212*Yd21) + Yd12*(md210*Yd02 + md211*Yd12 +
      md212*Yd22)) + 0.4*(md200 + md211 + md222 - mDx200 - mDx211 - mDx222 +
      mDxbar200 + mDxbar211 + mDxbar222 + me200 + me211 + me222 - mH1I200 -
      mH1I211 + mH2I200 + mH2I211 - mHd2 - mHp2 + mHpbar2 + mHu2 - ml200 - ml211 -
      ml222 + mq200 + mq211 + mq222 - 2*(mu200 + mu211 + mu222))*Sqr(g1) + 2*Qdp*
      (3*(md200 + md211 + md222)*Qdp + 3*(mDxbar200 + mDxbar211 + mDxbar222)*
      QDxbarp + 3*(mDx200 + mDx211 + mDx222)*QDxp + (me200 + me211 + me222)*Qep +
      2*(mH1I200 + mH1I211)*QH1p + 2*mHd2*QH1p + 2*(mH2I200 + mH2I211)*QH2p + 2*
      mHu2*QH2p + 2*mHpbar2*QHpbarp + 2*mHp2*QHpp + 2*(ml200 + ml211 + ml222)*QLp
      + 6*(mq200 + mq211 + mq222)*QQp + ms2*QSp + (msI200 + msI211)*QSp + 3*(mu200
      + mu211 + mu222)*Qup)*Sqr(gN) - 0.5333333333333333*Sqr(g1)*Sqr(MassB) -
      10.666666666666666*Sqr(g3)*Sqr(MassG) - 8*Sqr(gN)*Sqr(MassBp)*Sqr(Qdp) + 4*(
      Sqr(TYd10) + Sqr(TYd11) + Sqr(TYd12)) + 4*mHd2*(Sqr(Yd10) + Sqr(Yd11) + Sqr(
      Yd12)) + 2*(md201*(Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + md221*(Yd10*Yd20 +
      Yd11*Yd21 + Yd12*Yd22) + md211*(Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12)))) + 6*Qdp
      *QSp*Sqr(gN)*(2*(Yd20*(md220*Yd00 + md221*Yd10 + md222*Yd20) + Yd21*(md220*
      Yd01 + md221*Yd11 + md222*Yd21) + Yd22*(md220*Yd02 + md221*Yd12 + md222*Yd22
      )) + 4*(Yd20*(mq200*Yd20 + mq210*Yd21 + mq220*Yd22) + Yd21*(mq201*Yd20 +
      mq211*Yd21 + mq221*Yd22) + Yd22*(mq202*Yd20 + mq212*Yd21 + mq222*Yd22)) +
      0.4*(md200 + md211 + md222 - mDx200 - mDx211 - mDx222 + mDxbar200 +
      mDxbar211 + mDxbar222 + me200 + me211 + me222 - mH1I200 - mH1I211 + mH2I200
      + mH2I211 - mHd2 - mHp2 + mHpbar2 + mHu2 - ml200 - ml211 - ml222 + mq200 +
      mq211 + mq222 - 2*(mu200 + mu211 + mu222))*Sqr(g1) + 2*Qdp*(3*(md200 + md211
      + md222)*Qdp + 3*(mDxbar200 + mDxbar211 + mDxbar222)*QDxbarp + 3*(mDx200 +
      mDx211 + mDx222)*QDxp + (me200 + me211 + me222)*Qep + 2*(mH1I200 + mH1I211)*
      QH1p + 2*mHd2*QH1p + 2*(mH2I200 + mH2I211)*QH2p + 2*mHu2*QH2p + 2*mHpbar2*
      QHpbarp + 2*mHp2*QHpp + 2*(ml200 + ml211 + ml222)*QLp + 6*(mq200 + mq211 +
      mq222)*QQp + ms2*QSp + (msI200 + msI211)*QSp + 3*(mu200 + mu211 + mu222)*Qup
      )*Sqr(gN) - 0.5333333333333333*Sqr(g1)*Sqr(MassB) - 10.666666666666666*Sqr(
      g3)*Sqr(MassG) - 8*Sqr(gN)*Sqr(MassBp)*Sqr(Qdp) + 4*(Sqr(TYd20) + Sqr(TYd21)
      + Sqr(TYd22)) + 4*mHd2*(Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + 2*(md202*(Yd00
      *Yd20 + Yd01*Yd21 + Yd02*Yd22) + md212*(Yd10*Yd20 + Yd11*Yd21 + Yd12*Yd22) +
      md222*(Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)))) + 2*Qep*QSp*Sqr(gN)*(4*(Ye00*(
      ml200*Ye00 + ml210*Ye01 + ml220*Ye02) + Ye01*(ml201*Ye00 + ml211*Ye01 +
      ml221*Ye02) + Ye02*(ml202*Ye00 + ml212*Ye01 + ml222*Ye02)) + 2*(Ye00*(me200*
      Ye00 + me201*Ye10 + me202*Ye20) + Ye01*(me200*Ye01 + me201*Ye11 + me202*Ye21
      ) + Ye02*(me200*Ye02 + me201*Ye12 + me202*Ye22)) + 1.2*(md200 + md211 +
      md222 - mDx200 - mDx211 - mDx222 + mDxbar200 + mDxbar211 + mDxbar222 + me200
      + me211 + me222 - mH1I200 - mH1I211 + mH2I200 + mH2I211 - mHd2 - mHp2 +
      mHpbar2 + mHu2 - ml200 - ml211 - ml222 + mq200 + mq211 + mq222 - 2*(mu200 +
      mu211 + mu222))*Sqr(g1) + 2*Qep*(3*(md200 + md211 + md222)*Qdp + 3*(
      mDxbar200 + mDxbar211 + mDxbar222)*QDxbarp + 3*(mDx200 + mDx211 + mDx222)*
      QDxp + (me200 + me211 + me222)*Qep + 2*(mH1I200 + mH1I211)*QH1p + 2*mHd2*
      QH1p + 2*(mH2I200 + mH2I211)*QH2p + 2*mHu2*QH2p + 2*mHpbar2*QHpbarp + 2*mHp2
      *QHpp + 2*(ml200 + ml211 + ml222)*QLp + 6*(mq200 + mq211 + mq222)*QQp + ms2*
      QSp + (msI200 + msI211)*QSp + 3*(mu200 + mu211 + mu222)*Qup)*Sqr(gN) - 4.8*
      Sqr(g1)*Sqr(MassB) - 8*Sqr(gN)*Sqr(MassBp)*Sqr(Qep) + 4*(Sqr(TYe00) + Sqr(
      TYe01) + Sqr(TYe02)) + 4*mHd2*(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02)) + 2*(me210
      *(Ye00*Ye10 + Ye01*Ye11 + Ye02*Ye12) + me220*(Ye00*Ye20 + Ye01*Ye21 + Ye02*
      Ye22) + me200*(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02)))) + 2*Qep*QSp*Sqr(gN)*(4*(
      Ye10*(ml200*Ye10 + ml210*Ye11 + ml220*Ye12) + Ye11*(ml201*Ye10 + ml211*Ye11
      + ml221*Ye12) + Ye12*(ml202*Ye10 + ml212*Ye11 + ml222*Ye12)) + 2*(Ye10*(
      me210*Ye00 + me211*Ye10 + me212*Ye20) + Ye11*(me210*Ye01 + me211*Ye11 +
      me212*Ye21) + Ye12*(me210*Ye02 + me211*Ye12 + me212*Ye22)) + 1.2*(md200 +
      md211 + md222 - mDx200 - mDx211 - mDx222 + mDxbar200 + mDxbar211 + mDxbar222
      + me200 + me211 + me222 - mH1I200 - mH1I211 + mH2I200 + mH2I211 - mHd2 -
      mHp2 + mHpbar2 + mHu2 - ml200 - ml211 - ml222 + mq200 + mq211 + mq222 - 2*(
      mu200 + mu211 + mu222))*Sqr(g1) + 2*Qep*(3*(md200 + md211 + md222)*Qdp + 3*(
      mDxbar200 + mDxbar211 + mDxbar222)*QDxbarp + 3*(mDx200 + mDx211 + mDx222)*
      QDxp + (me200 + me211 + me222)*Qep + 2*(mH1I200 + mH1I211)*QH1p + 2*mHd2*
      QH1p + 2*(mH2I200 + mH2I211)*QH2p + 2*mHu2*QH2p + 2*mHpbar2*QHpbarp + 2*mHp2
      *QHpp + 2*(ml200 + ml211 + ml222)*QLp + 6*(mq200 + mq211 + mq222)*QQp + ms2*
      QSp + (msI200 + msI211)*QSp + 3*(mu200 + mu211 + mu222)*Qup)*Sqr(gN) - 4.8*
      Sqr(g1)*Sqr(MassB) - 8*Sqr(gN)*Sqr(MassBp)*Sqr(Qep) + 4*(Sqr(TYe10) + Sqr(
      TYe11) + Sqr(TYe12)) + 4*mHd2*(Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12)) + 2*(me201
      *(Ye00*Ye10 + Ye01*Ye11 + Ye02*Ye12) + me221*(Ye10*Ye20 + Ye11*Ye21 + Ye12*
      Ye22) + me211*(Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12)))) + 4*QLp*QSp*Sqr(gN)*(
      Ye00*(ml200*Ye00 + ml201*Ye01 + ml202*Ye02) + Ye10*(ml200*Ye10 + ml201*Ye11
      + ml202*Ye12) + 2*(Ye00*(me200*Ye00 + me210*Ye10 + me220*Ye20) + Ye10*(me201
      *Ye00 + me211*Ye10 + me221*Ye20) + Ye20*(me202*Ye00 + me212*Ye10 + me222*
      Ye20)) + ml210*(Ye00*Ye01 + Ye10*Ye11 + Ye20*Ye21) + Ye20*(ml200*Ye20 +
      ml201*Ye21 + ml202*Ye22) + ml220*(Ye00*Ye02 + Ye10*Ye12 + Ye20*Ye22) - 0.6*(
      md200 + md211 + md222 - mDx200 - mDx211 - mDx222 + mDxbar200 + mDxbar211 +
      mDxbar222 + me200 + me211 + me222 - mH1I200 - mH1I211 + mH2I200 + mH2I211 -
      mHd2 - mHp2 + mHpbar2 + mHu2 - ml200 - ml211 - ml222 + mq200 + mq211 + mq222
      - 2*(mu200 + mu211 + mu222))*Sqr(g1) + 2*QLp*(3*(md200 + md211 + md222)*Qdp
      + 3*(mDxbar200 + mDxbar211 + mDxbar222)*QDxbarp + 3*(mDx200 + mDx211 +
      mDx222)*QDxp + (me200 + me211 + me222)*Qep + 2*(mH1I200 + mH1I211)*QH1p + 2*
      mHd2*QH1p + 2*(mH2I200 + mH2I211)*QH2p + 2*mHu2*QH2p + 2*mHpbar2*QHpbarp + 2
      *mHp2*QHpp + 2*(ml200 + ml211 + ml222)*QLp + 6*(mq200 + mq211 + mq222)*QQp +
      ms2*QSp + (msI200 + msI211)*QSp + 3*(mu200 + mu211 + mu222)*Qup)*Sqr(gN) -
      1.2*Sqr(g1)*Sqr(MassB) - 6*Sqr(g2)*Sqr(MassWB) - 8*Sqr(gN)*Sqr(MassBp)*Sqr(
      QLp) + 2*(Sqr(TYe00) + Sqr(TYe10) + Sqr(TYe20)) + 2*mHd2*(Sqr(Ye00) + Sqr(
      Ye10) + Sqr(Ye20)) + ml200*(Sqr(Ye00) + Sqr(Ye10) + Sqr(Ye20))) + 4*QLp*QSp*
      Sqr(gN)*(Ye01*(ml210*Ye00 + ml211*Ye01 + ml212*Ye02) + Ye11*(ml210*Ye10 +
      ml211*Ye11 + ml212*Ye12) + ml201*(Ye00*Ye01 + Ye10*Ye11 + Ye20*Ye21) + 2*(
      Ye01*(me200*Ye01 + me210*Ye11 + me220*Ye21) + Ye11*(me201*Ye01 + me211*Ye11
      + me221*Ye21) + Ye21*(me202*Ye01 + me212*Ye11 + me222*Ye21)) + Ye21*(ml210*
      Ye20 + ml211*Ye21 + ml212*Ye22) + ml221*(Ye01*Ye02 + Ye11*Ye12 + Ye21*Ye22)
      - 0.6*(md200 + md211 + md222 - mDx200 - mDx211 - mDx222 + mDxbar200 +
      mDxbar211 + mDxbar222 + me200 + me211 + me222 - mH1I200 - mH1I211 + mH2I200
      + mH2I211 - mHd2 - mHp2 + mHpbar2 + mHu2 - ml200 - ml211 - ml222 + mq200 +
      mq211 + mq222 - 2*(mu200 + mu211 + mu222))*Sqr(g1) + 2*QLp*(3*(md200 + md211
      + md222)*Qdp + 3*(mDxbar200 + mDxbar211 + mDxbar222)*QDxbarp + 3*(mDx200 +
      mDx211 + mDx222)*QDxp + (me200 + me211 + me222)*Qep + 2*(mH1I200 + mH1I211)*
      QH1p + 2*mHd2*QH1p + 2*(mH2I200 + mH2I211)*QH2p + 2*mHu2*QH2p + 2*mHpbar2*
      QHpbarp + 2*mHp2*QHpp + 2*(ml200 + ml211 + ml222)*QLp + 6*(mq200 + mq211 +
      mq222)*QQp + ms2*QSp + (msI200 + msI211)*QSp + 3*(mu200 + mu211 + mu222)*Qup
      )*Sqr(gN) - 1.2*Sqr(g1)*Sqr(MassB) - 6*Sqr(g2)*Sqr(MassWB) - 8*Sqr(gN)*Sqr(
      MassBp)*Sqr(QLp) + 2*(Sqr(TYe01) + Sqr(TYe11) + Sqr(TYe21)) + 2*mHd2*(Sqr(
      Ye01) + Sqr(Ye11) + Sqr(Ye21)) + ml211*(Sqr(Ye01) + Sqr(Ye11) + Sqr(Ye21)))
      + 4*QLp*QSp*Sqr(gN)*(Ye02*(ml220*Ye00 + ml221*Ye01 + ml222*Ye02) + Ye12*(
      ml220*Ye10 + ml221*Ye11 + ml222*Ye12) + Ye22*(ml220*Ye20 + ml221*Ye21 +
      ml222*Ye22) + ml202*(Ye00*Ye02 + Ye10*Ye12 + Ye20*Ye22) + ml212*(Ye01*Ye02 +
      Ye11*Ye12 + Ye21*Ye22) + 2*(Ye02*(me200*Ye02 + me210*Ye12 + me220*Ye22) +
      Ye12*(me201*Ye02 + me211*Ye12 + me221*Ye22) + Ye22*(me202*Ye02 + me212*Ye12
      + me222*Ye22)) - 0.6*(md200 + md211 + md222 - mDx200 - mDx211 - mDx222 +
      mDxbar200 + mDxbar211 + mDxbar222 + me200 + me211 + me222 - mH1I200 -
      mH1I211 + mH2I200 + mH2I211 - mHd2 - mHp2 + mHpbar2 + mHu2 - ml200 - ml211 -
      ml222 + mq200 + mq211 + mq222 - 2*(mu200 + mu211 + mu222))*Sqr(g1) + 2*QLp*
      (3*(md200 + md211 + md222)*Qdp + 3*(mDxbar200 + mDxbar211 + mDxbar222)*
      QDxbarp + 3*(mDx200 + mDx211 + mDx222)*QDxp + (me200 + me211 + me222)*Qep +
      2*(mH1I200 + mH1I211)*QH1p + 2*mHd2*QH1p + 2*(mH2I200 + mH2I211)*QH2p + 2*
      mHu2*QH2p + 2*mHpbar2*QHpbarp + 2*mHp2*QHpp + 2*(ml200 + ml211 + ml222)*QLp
      + 6*(mq200 + mq211 + mq222)*QQp + ms2*QSp + (msI200 + msI211)*QSp + 3*(mu200
      + mu211 + mu222)*Qup)*Sqr(gN) - 1.2*Sqr(g1)*Sqr(MassB) - 6*Sqr(g2)*Sqr(
      MassWB) - 8*Sqr(gN)*Sqr(MassBp)*Sqr(QLp) + 2*(Sqr(TYe02) + Sqr(TYe12) + Sqr(
      TYe22)) + 2*mHd2*(Sqr(Ye02) + Sqr(Ye12) + Sqr(Ye22)) + ml222*(Sqr(Ye02) +
      Sqr(Ye12) + Sqr(Ye22))) + (4*QH1p*QSp*Sqr(gN) + 4*Sqr(Lambdax))*(6*(Yd00*(
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
      2*ms2*Sqr(Lambdax) - 1.2*Sqr(g1)*Sqr(MassB) - 6*Sqr(g2)*Sqr(MassWB) - 8*Sqr
      (gN)*Sqr(MassBp)*Sqr(QH1p) + 2*Sqr(TLambdax) + 6*(Sqr(TYd00) + Sqr(TYd01) +
      Sqr(TYd02) + Sqr(TYd10) + Sqr(TYd11) + Sqr(TYd12) + Sqr(TYd20) + Sqr(TYd21)
      + Sqr(TYd22)) + 2*(Sqr(TYe00) + Sqr(TYe01) + Sqr(TYe02) + Sqr(TYe10) + Sqr(
      TYe11) + Sqr(TYe12) + Sqr(TYe20) + Sqr(TYe21) + Sqr(TYe22)) + 6*mHd2*(Sqr(
      Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20
      ) + Sqr(Yd21) + Sqr(Yd22)) + 2*mHd2*(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr
      (Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22))) + 2*Qep
      *QSp*Sqr(gN)*(2*(Ye20*(me220*Ye00 + me221*Ye10 + me222*Ye20) + Ye21*(me220*
      Ye01 + me221*Ye11 + me222*Ye21) + Ye22*(me220*Ye02 + me221*Ye12 + me222*Ye22
      )) + 4*(Ye20*(ml200*Ye20 + ml210*Ye21 + ml220*Ye22) + Ye21*(ml201*Ye20 +
      ml211*Ye21 + ml221*Ye22) + Ye22*(ml202*Ye20 + ml212*Ye21 + ml222*Ye22)) +
      1.2*(md200 + md211 + md222 - mDx200 - mDx211 - mDx222 + mDxbar200 +
      mDxbar211 + mDxbar222 + me200 + me211 + me222 - mH1I200 - mH1I211 + mH2I200
      + mH2I211 - mHd2 - mHp2 + mHpbar2 + mHu2 - ml200 - ml211 - ml222 + mq200 +
      mq211 + mq222 - 2*(mu200 + mu211 + mu222))*Sqr(g1) + 2*Qep*(3*(md200 + md211
      + md222)*Qdp + 3*(mDxbar200 + mDxbar211 + mDxbar222)*QDxbarp + 3*(mDx200 +
      mDx211 + mDx222)*QDxp + (me200 + me211 + me222)*Qep + 2*(mH1I200 + mH1I211)*
      QH1p + 2*mHd2*QH1p + 2*(mH2I200 + mH2I211)*QH2p + 2*mHu2*QH2p + 2*mHpbar2*
      QHpbarp + 2*mHp2*QHpp + 2*(ml200 + ml211 + ml222)*QLp + 6*(mq200 + mq211 +
      mq222)*QQp + ms2*QSp + (msI200 + msI211)*QSp + 3*(mu200 + mu211 + mu222)*Qup
      )*Sqr(gN) - 4.8*Sqr(g1)*Sqr(MassB) - 8*Sqr(gN)*Sqr(MassBp)*Sqr(Qep) + 4*(Sqr
      (TYe20) + Sqr(TYe21) + Sqr(TYe22)) + 4*mHd2*(Sqr(Ye20) + Sqr(Ye21) + Sqr(
      Ye22)) + 2*(me202*(Ye00*Ye20 + Ye01*Ye21 + Ye02*Ye22) + me212*(Ye10*Ye20 +
      Ye11*Ye21 + Ye12*Ye22) + me222*(Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22)))) + 6*QSp
      *Qup*Sqr(gN)*(4*(Yu00*(mq200*Yu00 + mq210*Yu01 + mq220*Yu02) + Yu01*(mq201*
      Yu00 + mq211*Yu01 + mq221*Yu02) + Yu02*(mq202*Yu00 + mq212*Yu01 + mq222*Yu02
      )) + 2*(Yu00*(mu200*Yu00 + mu201*Yu10 + mu202*Yu20) + Yu01*(mu200*Yu01 +
      mu201*Yu11 + mu202*Yu21) + Yu02*(mu200*Yu02 + mu201*Yu12 + mu202*Yu22)) -
      0.8*(md200 + md211 + md222 - mDx200 - mDx211 - mDx222 + mDxbar200 +
      mDxbar211 + mDxbar222 + me200 + me211 + me222 - mH1I200 - mH1I211 + mH2I200
      + mH2I211 - mHd2 - mHp2 + mHpbar2 + mHu2 - ml200 - ml211 - ml222 + mq200 +
      mq211 + mq222 - 2*(mu200 + mu211 + mu222))*Sqr(g1) + 2*Qup*(3*(md200 + md211
      + md222)*Qdp + 3*(mDxbar200 + mDxbar211 + mDxbar222)*QDxbarp + 3*(mDx200 +
      mDx211 + mDx222)*QDxp + (me200 + me211 + me222)*Qep + 2*(mH1I200 + mH1I211)*
      QH1p + 2*mHd2*QH1p + 2*(mH2I200 + mH2I211)*QH2p + 2*mHu2*QH2p + 2*mHpbar2*
      QHpbarp + 2*mHp2*QHpp + 2*(ml200 + ml211 + ml222)*QLp + 6*(mq200 + mq211 +
      mq222)*QQp + ms2*QSp + (msI200 + msI211)*QSp + 3*(mu200 + mu211 + mu222)*Qup
      )*Sqr(gN) - 2.1333333333333333*Sqr(g1)*Sqr(MassB) - 10.666666666666666*Sqr(
      g3)*Sqr(MassG) - 8*Sqr(gN)*Sqr(MassBp)*Sqr(Qup) + 4*(Sqr(TYu00) + Sqr(TYu01)
      + Sqr(TYu02)) + 4*mHu2*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02)) + 2*(mu210*(Yu00
      *Yu10 + Yu01*Yu11 + Yu02*Yu12) + mu220*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) +
      mu200*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02)))) + 6*QSp*Qup*Sqr(gN)*(4*(Yu10*(
      mq200*Yu10 + mq210*Yu11 + mq220*Yu12) + Yu11*(mq201*Yu10 + mq211*Yu11 +
      mq221*Yu12) + Yu12*(mq202*Yu10 + mq212*Yu11 + mq222*Yu12)) + 2*(Yu10*(mu210*
      Yu00 + mu211*Yu10 + mu212*Yu20) + Yu11*(mu210*Yu01 + mu211*Yu11 + mu212*Yu21
      ) + Yu12*(mu210*Yu02 + mu211*Yu12 + mu212*Yu22)) - 0.8*(md200 + md211 +
      md222 - mDx200 - mDx211 - mDx222 + mDxbar200 + mDxbar211 + mDxbar222 + me200
      + me211 + me222 - mH1I200 - mH1I211 + mH2I200 + mH2I211 - mHd2 - mHp2 +
      mHpbar2 + mHu2 - ml200 - ml211 - ml222 + mq200 + mq211 + mq222 - 2*(mu200 +
      mu211 + mu222))*Sqr(g1) + 2*Qup*(3*(md200 + md211 + md222)*Qdp + 3*(
      mDxbar200 + mDxbar211 + mDxbar222)*QDxbarp + 3*(mDx200 + mDx211 + mDx222)*
      QDxp + (me200 + me211 + me222)*Qep + 2*(mH1I200 + mH1I211)*QH1p + 2*mHd2*
      QH1p + 2*(mH2I200 + mH2I211)*QH2p + 2*mHu2*QH2p + 2*mHpbar2*QHpbarp + 2*mHp2
      *QHpp + 2*(ml200 + ml211 + ml222)*QLp + 6*(mq200 + mq211 + mq222)*QQp + ms2*
      QSp + (msI200 + msI211)*QSp + 3*(mu200 + mu211 + mu222)*Qup)*Sqr(gN) -
      2.1333333333333333*Sqr(g1)*Sqr(MassB) - 10.666666666666666*Sqr(g3)*Sqr(MassG
      ) - 8*Sqr(gN)*Sqr(MassBp)*Sqr(Qup) + 4*(Sqr(TYu10) + Sqr(TYu11) + Sqr(TYu12)
      ) + 4*mHu2*(Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12)) + 2*(mu201*(Yu00*Yu10 + Yu01*
      Yu11 + Yu02*Yu12) + mu221*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + mu211*(Sqr(
      Yu10) + Sqr(Yu11) + Sqr(Yu12)))) + 12*QQp*QSp*Sqr(gN)*(Yd00*(mq200*Yd00 +
      mq201*Yd01 + mq202*Yd02) + Yd10*(mq200*Yd10 + mq201*Yd11 + mq202*Yd12) + 2*(
      Yd00*(md200*Yd00 + md210*Yd10 + md220*Yd20) + Yd10*(md201*Yd00 + md211*Yd10
      + md221*Yd20) + Yd20*(md202*Yd00 + md212*Yd10 + md222*Yd20)) + mq210*(Yd00*
      Yd01 + Yd10*Yd11 + Yd20*Yd21) + Yd20*(mq200*Yd20 + mq201*Yd21 + mq202*Yd22)
      + mq220*(Yd00*Yd02 + Yd10*Yd12 + Yd20*Yd22) + Yu00*(mq200*Yu00 + mq201*Yu01
      + mq202*Yu02) + Yu10*(mq200*Yu10 + mq201*Yu11 + mq202*Yu12) + 2*(Yu00*(mu200
      *Yu00 + mu210*Yu10 + mu220*Yu20) + Yu10*(mu201*Yu00 + mu211*Yu10 + mu221*
      Yu20) + Yu20*(mu202*Yu00 + mu212*Yu10 + mu222*Yu20)) + mq210*(Yu00*Yu01 +
      Yu10*Yu11 + Yu20*Yu21) + Yu20*(mq200*Yu20 + mq201*Yu21 + mq202*Yu22) + mq220
      *(Yu00*Yu02 + Yu10*Yu12 + Yu20*Yu22) + 0.2*(md200 + md211 + md222 - mDx200 -
      mDx211 - mDx222 + mDxbar200 + mDxbar211 + mDxbar222 + me200 + me211 + me222
      - mH1I200 - mH1I211 + mH2I200 + mH2I211 - mHd2 - mHp2 + mHpbar2 + mHu2 -
      ml200 - ml211 - ml222 + mq200 + mq211 + mq222 - 2*(mu200 + mu211 + mu222))*
      Sqr(g1) + 2*QQp*(3*(md200 + md211 + md222)*Qdp + 3*(mDxbar200 + mDxbar211 +
      mDxbar222)*QDxbarp + 3*(mDx200 + mDx211 + mDx222)*QDxp + (me200 + me211 +
      me222)*Qep + 2*(mH1I200 + mH1I211)*QH1p + 2*mHd2*QH1p + 2*(mH2I200 + mH2I211
      )*QH2p + 2*mHu2*QH2p + 2*mHpbar2*QHpbarp + 2*mHp2*QHpp + 2*(ml200 + ml211 +
      ml222)*QLp + 6*(mq200 + mq211 + mq222)*QQp + ms2*QSp + (msI200 + msI211)*QSp
      + 3*(mu200 + mu211 + mu222)*Qup)*Sqr(gN) - 0.13333333333333333*Sqr(g1)*Sqr(
      MassB) - 10.666666666666666*Sqr(g3)*Sqr(MassG) - 6*Sqr(g2)*Sqr(MassWB) - 8*
      Sqr(gN)*Sqr(MassBp)*Sqr(QQp) + 2*(Sqr(TYd00) + Sqr(TYd10) + Sqr(TYd20)) + 2*
      (Sqr(TYu00) + Sqr(TYu10) + Sqr(TYu20)) + 2*mHd2*(Sqr(Yd00) + Sqr(Yd10) + Sqr
      (Yd20)) + mq200*(Sqr(Yd00) + Sqr(Yd10) + Sqr(Yd20)) + 2*mHu2*(Sqr(Yu00) +
      Sqr(Yu10) + Sqr(Yu20)) + mq200*(Sqr(Yu00) + Sqr(Yu10) + Sqr(Yu20))) + 12*QQp
      *QSp*Sqr(gN)*(Yd01*(mq210*Yd00 + mq211*Yd01 + mq212*Yd02) + Yd11*(mq210*Yd10
      + mq211*Yd11 + mq212*Yd12) + mq201*(Yd00*Yd01 + Yd10*Yd11 + Yd20*Yd21) + 2*
      (Yd01*(md200*Yd01 + md210*Yd11 + md220*Yd21) + Yd11*(md201*Yd01 + md211*Yd11
      + md221*Yd21) + Yd21*(md202*Yd01 + md212*Yd11 + md222*Yd21)) + Yd21*(mq210*
      Yd20 + mq211*Yd21 + mq212*Yd22) + mq221*(Yd01*Yd02 + Yd11*Yd12 + Yd21*Yd22)
      + Yu01*(mq210*Yu00 + mq211*Yu01 + mq212*Yu02) + Yu11*(mq210*Yu10 + mq211*
      Yu11 + mq212*Yu12) + mq201*(Yu00*Yu01 + Yu10*Yu11 + Yu20*Yu21) + 2*(Yu01*(
      mu200*Yu01 + mu210*Yu11 + mu220*Yu21) + Yu11*(mu201*Yu01 + mu211*Yu11 +
      mu221*Yu21) + Yu21*(mu202*Yu01 + mu212*Yu11 + mu222*Yu21)) + Yu21*(mq210*
      Yu20 + mq211*Yu21 + mq212*Yu22) + mq221*(Yu01*Yu02 + Yu11*Yu12 + Yu21*Yu22)
      + 0.2*(md200 + md211 + md222 - mDx200 - mDx211 - mDx222 + mDxbar200 +
      mDxbar211 + mDxbar222 + me200 + me211 + me222 - mH1I200 - mH1I211 + mH2I200
      + mH2I211 - mHd2 - mHp2 + mHpbar2 + mHu2 - ml200 - ml211 - ml222 + mq200 +
      mq211 + mq222 - 2*(mu200 + mu211 + mu222))*Sqr(g1) + 2*QQp*(3*(md200 + md211
      + md222)*Qdp + 3*(mDxbar200 + mDxbar211 + mDxbar222)*QDxbarp + 3*(mDx200 +
      mDx211 + mDx222)*QDxp + (me200 + me211 + me222)*Qep + 2*(mH1I200 + mH1I211)*
      QH1p + 2*mHd2*QH1p + 2*(mH2I200 + mH2I211)*QH2p + 2*mHu2*QH2p + 2*mHpbar2*
      QHpbarp + 2*mHp2*QHpp + 2*(ml200 + ml211 + ml222)*QLp + 6*(mq200 + mq211 +
      mq222)*QQp + ms2*QSp + (msI200 + msI211)*QSp + 3*(mu200 + mu211 + mu222)*Qup
      )*Sqr(gN) - 0.13333333333333333*Sqr(g1)*Sqr(MassB) - 10.666666666666666*Sqr(
      g3)*Sqr(MassG) - 6*Sqr(g2)*Sqr(MassWB) - 8*Sqr(gN)*Sqr(MassBp)*Sqr(QQp) + 2*
      (Sqr(TYd01) + Sqr(TYd11) + Sqr(TYd21)) + 2*(Sqr(TYu01) + Sqr(TYu11) + Sqr(
      TYu21)) + 2*mHd2*(Sqr(Yd01) + Sqr(Yd11) + Sqr(Yd21)) + mq211*(Sqr(Yd01) +
      Sqr(Yd11) + Sqr(Yd21)) + 2*mHu2*(Sqr(Yu01) + Sqr(Yu11) + Sqr(Yu21)) + mq211*
      (Sqr(Yu01) + Sqr(Yu11) + Sqr(Yu21))) + 12*QQp*QSp*Sqr(gN)*(Yd02*(mq220*Yd00
      + mq221*Yd01 + mq222*Yd02) + Yd12*(mq220*Yd10 + mq221*Yd11 + mq222*Yd12) +
      Yd22*(mq220*Yd20 + mq221*Yd21 + mq222*Yd22) + mq202*(Yd00*Yd02 + Yd10*Yd12 +
      Yd20*Yd22) + mq212*(Yd01*Yd02 + Yd11*Yd12 + Yd21*Yd22) + 2*(Yd02*(md200*
      Yd02 + md210*Yd12 + md220*Yd22) + Yd12*(md201*Yd02 + md211*Yd12 + md221*Yd22
      ) + Yd22*(md202*Yd02 + md212*Yd12 + md222*Yd22)) + Yu02*(mq220*Yu00 + mq221*
      Yu01 + mq222*Yu02) + Yu12*(mq220*Yu10 + mq221*Yu11 + mq222*Yu12) + Yu22*(
      mq220*Yu20 + mq221*Yu21 + mq222*Yu22) + mq202*(Yu00*Yu02 + Yu10*Yu12 + Yu20*
      Yu22) + mq212*(Yu01*Yu02 + Yu11*Yu12 + Yu21*Yu22) + 2*(Yu02*(mu200*Yu02 +
      mu210*Yu12 + mu220*Yu22) + Yu12*(mu201*Yu02 + mu211*Yu12 + mu221*Yu22) +
      Yu22*(mu202*Yu02 + mu212*Yu12 + mu222*Yu22)) + 0.2*(md200 + md211 + md222 -
      mDx200 - mDx211 - mDx222 + mDxbar200 + mDxbar211 + mDxbar222 + me200 + me211
      + me222 - mH1I200 - mH1I211 + mH2I200 + mH2I211 - mHd2 - mHp2 + mHpbar2 +
      mHu2 - ml200 - ml211 - ml222 + mq200 + mq211 + mq222 - 2*(mu200 + mu211 +
      mu222))*Sqr(g1) + 2*QQp*(3*(md200 + md211 + md222)*Qdp + 3*(mDxbar200 +
      mDxbar211 + mDxbar222)*QDxbarp + 3*(mDx200 + mDx211 + mDx222)*QDxp + (me200
      + me211 + me222)*Qep + 2*(mH1I200 + mH1I211)*QH1p + 2*mHd2*QH1p + 2*(mH2I200
      + mH2I211)*QH2p + 2*mHu2*QH2p + 2*mHpbar2*QHpbarp + 2*mHp2*QHpp + 2*(ml200
      + ml211 + ml222)*QLp + 6*(mq200 + mq211 + mq222)*QQp + ms2*QSp + (msI200 +
      msI211)*QSp + 3*(mu200 + mu211 + mu222)*Qup)*Sqr(gN) - 0.13333333333333333*
      Sqr(g1)*Sqr(MassB) - 10.666666666666666*Sqr(g3)*Sqr(MassG) - 6*Sqr(g2)*Sqr(
      MassWB) - 8*Sqr(gN)*Sqr(MassBp)*Sqr(QQp) + 2*(Sqr(TYd02) + Sqr(TYd12) + Sqr(
      TYd22)) + 2*(Sqr(TYu02) + Sqr(TYu12) + Sqr(TYu22)) + 2*mHd2*(Sqr(Yd02) + Sqr
      (Yd12) + Sqr(Yd22)) + mq222*(Sqr(Yd02) + Sqr(Yd12) + Sqr(Yd22)) + 2*mHu2*(
      Sqr(Yu02) + Sqr(Yu12) + Sqr(Yu22)) + mq222*(Sqr(Yu02) + Sqr(Yu12) + Sqr(Yu22
      ))) + 8*Lambdax*(mHd2 + mHu2 + ms2)*(4*Power(Lambdax,3) - 0.6*Lambdax*Sqr(g1
      ) - 3*Lambdax*Sqr(g2) + 3*Lambdax*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02
      ) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21)
      + Sqr(Kappa22)) + 2*Lambdax*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(
      Lambda1210) + Sqr(Lambda1211)) - 2*Lambdax*Sqr(gN)*Sqr(QH1p) - 2*Lambdax*Sqr
      (gN)*Sqr(QH2p) - 2*Lambdax*Sqr(gN)*Sqr(QSp) + 3*Lambdax*(Sqr(Yd00) + Sqr(
      Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21
      ) + Sqr(Yd22)) + Lambdax*(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) +
      Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22)) + 3*Lambdax*(Sqr(
      Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20
      ) + Sqr(Yu21) + Sqr(Yu22))) + (4*QH2p*QSp*Sqr(gN) + 4*Sqr(Lambdax))*(6*(Yu00
      *(mq200*Yu00 + mq201*Yu01 + mq202*Yu02) + Yu01*(mq210*Yu00 + mq211*Yu01 +
      mq212*Yu02) + Yu02*(mq220*Yu00 + mq221*Yu01 + mq222*Yu02) + Yu10*(mq200*Yu10
      + mq201*Yu11 + mq202*Yu12) + Yu11*(mq210*Yu10 + mq211*Yu11 + mq212*Yu12) +
      Yu12*(mq220*Yu10 + mq221*Yu11 + mq222*Yu12) + Yu20*(mq200*Yu20 + mq201*Yu21
      + mq202*Yu22) + Yu21*(mq210*Yu20 + mq211*Yu21 + mq212*Yu22) + Yu22*(mq220*
      Yu20 + mq221*Yu21 + mq222*Yu22)) + 6*(Yu00*(mu200*Yu00 + mu201*Yu10 + mu202*
      Yu20) + Yu10*(mu210*Yu00 + mu211*Yu10 + mu212*Yu20) + Yu20*(mu220*Yu00 +
      mu221*Yu10 + mu222*Yu20) + Yu01*(mu200*Yu01 + mu201*Yu11 + mu202*Yu21) +
      Yu11*(mu210*Yu01 + mu211*Yu11 + mu212*Yu21) + Yu21*(mu220*Yu01 + mu221*Yu11
      + mu222*Yu21) + Yu02*(mu200*Yu02 + mu201*Yu12 + mu202*Yu22) + Yu12*(mu210*
      Yu02 + mu211*Yu12 + mu212*Yu22) + Yu22*(mu220*Yu02 + mu221*Yu12 + mu222*Yu22
      )) + 0.6*(md200 + md211 + md222 - mDx200 - mDx211 - mDx222 + mDxbar200 +
      mDxbar211 + mDxbar222 + me200 + me211 + me222 - mH1I200 - mH1I211 + mH2I200
      + mH2I211 - mHd2 - mHp2 + mHpbar2 + mHu2 - ml200 - ml211 - ml222 + mq200 +
      mq211 + mq222 - 2*(mu200 + mu211 + mu222))*Sqr(g1) + 2*QH2p*(3*(md200 +
      md211 + md222)*Qdp + 3*(mDxbar200 + mDxbar211 + mDxbar222)*QDxbarp + 3*(
      mDx200 + mDx211 + mDx222)*QDxp + (me200 + me211 + me222)*Qep + 2*(mH1I200 +
      mH1I211)*QH1p + 2*mHd2*QH1p + 2*(mH2I200 + mH2I211)*QH2p + 2*mHu2*QH2p + 2*
      mHpbar2*QHpbarp + 2*mHp2*QHpp + 2*(ml200 + ml211 + ml222)*QLp + 6*(mq200 +
      mq211 + mq222)*QQp + ms2*QSp + (msI200 + msI211)*QSp + 3*(mu200 + mu211 +
      mu222)*Qup)*Sqr(gN) + 2*mHd2*Sqr(Lambdax) + 2*mHu2*Sqr(Lambdax) + 2*ms2*Sqr(
      Lambdax) - 1.2*Sqr(g1)*Sqr(MassB) - 6*Sqr(g2)*Sqr(MassWB) - 8*Sqr(gN)*Sqr(
      MassBp)*Sqr(QH2p) + 2*Sqr(TLambdax) + 6*(Sqr(TYu00) + Sqr(TYu01) + Sqr(TYu02
      ) + Sqr(TYu10) + Sqr(TYu11) + Sqr(TYu12) + Sqr(TYu20) + Sqr(TYu21) + Sqr(
      TYu22)) + 6*mHu2*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11)
      + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + 6*QSp*Qup*Sqr(gN)*(4*(
      Yu20*(mq200*Yu20 + mq210*Yu21 + mq220*Yu22) + Yu21*(mq201*Yu20 + mq211*Yu21
      + mq221*Yu22) + Yu22*(mq202*Yu20 + mq212*Yu21 + mq222*Yu22)) + 2*(Yu20*(
      mu220*Yu00 + mu221*Yu10 + mu222*Yu20) + Yu21*(mu220*Yu01 + mu221*Yu11 +
      mu222*Yu21) + Yu22*(mu220*Yu02 + mu221*Yu12 + mu222*Yu22)) - 0.8*(md200 +
      md211 + md222 - mDx200 - mDx211 - mDx222 + mDxbar200 + mDxbar211 + mDxbar222
      + me200 + me211 + me222 - mH1I200 - mH1I211 + mH2I200 + mH2I211 - mHd2 -
      mHp2 + mHpbar2 + mHu2 - ml200 - ml211 - ml222 + mq200 + mq211 + mq222 - 2*(
      mu200 + mu211 + mu222))*Sqr(g1) + 2*Qup*(3*(md200 + md211 + md222)*Qdp + 3*(
      mDxbar200 + mDxbar211 + mDxbar222)*QDxbarp + 3*(mDx200 + mDx211 + mDx222)*
      QDxp + (me200 + me211 + me222)*Qep + 2*(mH1I200 + mH1I211)*QH1p + 2*mHd2*
      QH1p + 2*(mH2I200 + mH2I211)*QH2p + 2*mHu2*QH2p + 2*mHpbar2*QHpbarp + 2*mHp2
      *QHpp + 2*(ml200 + ml211 + ml222)*QLp + 6*(mq200 + mq211 + mq222)*QQp + ms2*
      QSp + (msI200 + msI211)*QSp + 3*(mu200 + mu211 + mu222)*Qup)*Sqr(gN) -
      2.1333333333333333*Sqr(g1)*Sqr(MassB) - 10.666666666666666*Sqr(g3)*Sqr(MassG
      ) - 8*Sqr(gN)*Sqr(MassBp)*Sqr(Qup) + 4*(Sqr(TYu20) + Sqr(TYu21) + Sqr(TYu22)
      ) + 4*mHu2*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)) + 2*(mu202*(Yu00*Yu20 + Yu01*
      Yu21 + Yu02*Yu22) + mu212*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + mu222*(Sqr(
      Yu20) + Sqr(Yu21) + Sqr(Yu22)))) + 8*TLambdax*(6*Lambdax*(Kappa00*TKappa00 +
      Kappa01*TKappa01 + Kappa02*TKappa02 + Kappa10*TKappa10 + Kappa11*TKappa11 +
      Kappa12*TKappa12 + Kappa20*TKappa20 + Kappa21*TKappa21 + Kappa22*TKappa22)
      + 4*Lambdax*(Lambda1200*TLambda1200 + Lambda1201*TLambda1201 + Lambda1210*
      TLambda1210 + Lambda1211*TLambda1211) + 6*Lambdax*(TYd00*Yd00 + TYd01*Yd01 +
      TYd02*Yd02 + TYd10*Yd10 + TYd11*Yd11 + TYd12*Yd12 + TYd20*Yd20 + TYd21*Yd21
      + TYd22*Yd22) + 2*Lambdax*(TYe00*Ye00 + TYe01*Ye01 + TYe02*Ye02 + TYe10*
      Ye10 + TYe11*Ye11 + TYe12*Ye12 + TYe20*Ye20 + TYe21*Ye21 + TYe22*Ye22) + 6*
      Lambdax*(TYu00*Yu00 + TYu01*Yu01 + TYu02*Yu02 + TYu10*Yu10 + TYu11*Yu11 +
      TYu12*Yu12 + TYu20*Yu20 + TYu21*Yu21 + TYu22*Yu22) + 1.2*Lambdax*MassB*Sqr(
      g1) + 6*Lambdax*MassWB*Sqr(g2) + 4*Lambdax*MassBp*Sqr(gN)*Sqr(QH1p) + 4*
      Lambdax*MassBp*Sqr(gN)*Sqr(QH2p) + 4*Lambdax*MassBp*Sqr(gN)*Sqr(QSp) +
      TLambdax*(-0.6*Sqr(g1) - 3*Sqr(g2) + 3*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(
      Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(
      Kappa21) + Sqr(Kappa22)) + 2*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(
      Lambda1210) + Sqr(Lambda1211)) + 12*Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH1p) - 2*
      Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QSp) + 3*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02
      ) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) +
      Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr
      (Ye20) + Sqr(Ye21) + Sqr(Ye22) + 3*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(
      Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))));


   return twoLoop*coeff;
}

} // namespace flexiblesusy
