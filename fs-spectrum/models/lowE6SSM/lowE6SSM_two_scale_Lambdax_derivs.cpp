#include "lowE6SSM_two_scale_tuning_calculator.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

double lowE6SSM_tuning_calculator::deriv_dbeta_one_loop_Lambdax_dLambdax() const
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

   double deriv = -0.6*Sqr(g1) - 3*Sqr(g2) + 3*(Sqr(Kappa00) + Sqr(Kappa01) +
      Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) +
      Sqr(Kappa21) + Sqr(Kappa22)) + 2*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(
      Lambda1210) + Sqr(Lambda1211)) + 12*Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH1p) - 2*
      Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QSp) + 3*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02
      ) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) +
      Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr
      (Ye20) + Sqr(Ye21) + Sqr(Ye22) + 3*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(
      Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22));


   return oneOver16PiSqr * deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_two_loop_Lambdax_dLambdax() const
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

   double deriv = -0.02*Lambdax*(2000*Power(Lambdax,3) - 20*Lambdax*(6*Sqr(g1)
      + 30*Sqr(g2) - 30*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10)
      + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22))
      - 20*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)
      ) + 20*Sqr(gN)*Sqr(QH1p) + 20*Sqr(gN)*Sqr(QH2p) - 45*(Sqr(Yd00) + Sqr(Yd01)
      + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) +
      Sqr(Yd22)) - 15*(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) +
      Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22)) - 45*(Sqr(Yu00) + Sqr(Yu01)
      + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) +
      Sqr(Yu22)))) + 0.02*(297*Power(g1,4) + 825*Power(g2,4) - 500*Power(Lambdax,4
      ) + 800*Power(gN,4)*Power(QH1p,4) + 800*Power(gN,4)*Power(QH2p,4) + 500*
      Power(gN,4)*Power(QSp,4) - 300*(Yd00*(Yu00*(Yd00*Yu00 + Yd01*Yu01 + Yd02*
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
      + Yd22*Yu12) + Yu22*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22))) + 90*Sqr(g1)*Sqr(
      g2) - 180*Qdp*QH1p*Sqr(g1)*Sqr(gN) - 180*QDxbarp*QH1p*Sqr(g1)*Sqr(gN) + 180*
      QDxp*QH1p*Sqr(g1)*Sqr(gN) - 180*Qep*QH1p*Sqr(g1)*Sqr(gN) + 180*Qdp*QH2p*Sqr(
      g1)*Sqr(gN) + 180*QDxbarp*QH2p*Sqr(g1)*Sqr(gN) - 180*QDxp*QH2p*Sqr(g1)*Sqr(
      gN) + 180*Qep*QH2p*Sqr(g1)*Sqr(gN) - 360*QH1p*QH2p*Sqr(g1)*Sqr(gN) - 60*QH1p
      *QHpbarp*Sqr(g1)*Sqr(gN) + 60*QH2p*QHpbarp*Sqr(g1)*Sqr(gN) + 60*QH1p*QHpp*
      Sqr(g1)*Sqr(gN) - 60*QH2p*QHpp*Sqr(g1)*Sqr(gN) + 180*QH1p*QLp*Sqr(g1)*Sqr(gN
      ) - 180*QH2p*QLp*Sqr(g1)*Sqr(gN) - 180*QH1p*QQp*Sqr(g1)*Sqr(gN) + 180*QH2p*
      QQp*Sqr(g1)*Sqr(gN) + 360*QH1p*Qup*Sqr(g1)*Sqr(gN) - 360*QH2p*Qup*Sqr(g1)*
      Sqr(gN) + 40*Sqr(g1)*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(
      Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(
      Kappa22)) + 800*Sqr(g3)*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(
      Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(
      Kappa22)) - 300*(Kappa00*(Kappa10*(Kappa00*Kappa10 + Kappa01*Kappa11 +
      Kappa02*Kappa12) + Kappa20*(Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*
      Kappa22) + Kappa00*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02))) + Kappa01*(
      Kappa11*(Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12) + Kappa21*(
      Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22) + Kappa01*(Sqr(Kappa00)
      + Sqr(Kappa01) + Sqr(Kappa02))) + Kappa02*(Kappa12*(Kappa00*Kappa10 +
      Kappa01*Kappa11 + Kappa02*Kappa12) + Kappa22*(Kappa00*Kappa20 + Kappa01*
      Kappa21 + Kappa02*Kappa22) + Kappa02*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(
      Kappa02))) + Kappa10*(Kappa00*(Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*
      Kappa12) + Kappa20*(Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22) +
      Kappa10*(Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12))) + Kappa11*(Kappa01*(
      Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12) + Kappa21*(Kappa10*
      Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22) + Kappa11*(Sqr(Kappa10) + Sqr(
      Kappa11) + Sqr(Kappa12))) + Kappa12*(Kappa02*(Kappa00*Kappa10 + Kappa01*
      Kappa11 + Kappa02*Kappa12) + Kappa22*(Kappa10*Kappa20 + Kappa11*Kappa21 +
      Kappa12*Kappa22) + Kappa12*(Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12))) +
      Kappa20*(Kappa00*(Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22) +
      Kappa10*(Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22) + Kappa20*(Sqr
      (Kappa20) + Sqr(Kappa21) + Sqr(Kappa22))) + Kappa21*(Kappa01*(Kappa00*
      Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22) + Kappa11*(Kappa10*Kappa20 +
      Kappa11*Kappa21 + Kappa12*Kappa22) + Kappa21*(Sqr(Kappa20) + Sqr(Kappa21) +
      Sqr(Kappa22))) + Kappa22*(Kappa02*(Kappa00*Kappa20 + Kappa01*Kappa21 +
      Kappa02*Kappa22) + Kappa12*(Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*
      Kappa22) + Kappa22*(Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)))) + 60*Sqr(
      g1)*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211))
      + 300*Sqr(g2)*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(
      Lambda1211)) - 200*(Lambda1200*(Lambda1210*(Lambda1200*Lambda1210 +
      Lambda1201*Lambda1211) + Lambda1200*(Sqr(Lambda1200) + Sqr(Lambda1201))) +
      Lambda1201*(Lambda1211*(Lambda1200*Lambda1210 + Lambda1201*Lambda1211) +
      Lambda1201*(Sqr(Lambda1200) + Sqr(Lambda1201))) + Lambda1210*(Lambda1200*(
      Lambda1200*Lambda1210 + Lambda1201*Lambda1211) + Lambda1210*(Sqr(Lambda1210)
      + Sqr(Lambda1211))) + Lambda1211*(Lambda1201*(Lambda1200*Lambda1210 +
      Lambda1201*Lambda1211) + Lambda1211*(Sqr(Lambda1210) + Sqr(Lambda1211)))) +
      300*Sqr(gN)*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr
      (Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22))*Sqr(
      QDxbarp) + 300*Sqr(gN)*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(
      Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(
      Kappa22))*Sqr(QDxp) + 240*Sqr(g1)*Sqr(gN)*Sqr(QH1p) + 300*Sqr(g2)*Sqr(gN)*
      Sqr(QH1p) + 200*Sqr(gN)*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210)
      + Sqr(Lambda1211))*Sqr(QH1p) + 900*Power(gN,4)*Sqr(Qdp)*Sqr(QH1p) + 900*
      Power(gN,4)*Sqr(QDxbarp)*Sqr(QH1p) + 900*Power(gN,4)*Sqr(QDxp)*Sqr(QH1p) +
      300*Power(gN,4)*Sqr(Qep)*Sqr(QH1p) + 240*Sqr(g1)*Sqr(gN)*Sqr(QH2p) + 300*Sqr
      (g2)*Sqr(gN)*Sqr(QH2p) + 200*Sqr(gN)*(Sqr(Lambda1200) + Sqr(Lambda1201) +
      Sqr(Lambda1210) + Sqr(Lambda1211))*Sqr(QH2p) + 900*Power(gN,4)*Sqr(Qdp)*Sqr(
      QH2p) + 900*Power(gN,4)*Sqr(QDxbarp)*Sqr(QH2p) + 900*Power(gN,4)*Sqr(QDxp)*
      Sqr(QH2p) + 300*Power(gN,4)*Sqr(Qep)*Sqr(QH2p) + 1200*Power(gN,4)*Sqr(QH1p)*
      Sqr(QH2p) + 200*Power(gN,4)*Sqr(QH1p)*Sqr(QHpbarp) + 200*Power(gN,4)*Sqr(
      QH2p)*Sqr(QHpbarp) + 200*Power(gN,4)*Sqr(QH1p)*Sqr(QHpp) + 200*Power(gN,4)*
      Sqr(QH2p)*Sqr(QHpp) + 600*Power(gN,4)*Sqr(QH1p)*Sqr(QLp) + 600*Power(gN,4)*
      Sqr(QH2p)*Sqr(QLp) + 1800*Power(gN,4)*Sqr(QH1p)*Sqr(QQp) + 1800*Power(gN,4)*
      Sqr(QH2p)*Sqr(QQp) - 300*Sqr(gN)*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02)
      + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21)
      + Sqr(Kappa22))*Sqr(QSp) - 200*Sqr(gN)*(Sqr(Lambda1200) + Sqr(Lambda1201) +
      Sqr(Lambda1210) + Sqr(Lambda1211))*Sqr(QSp) + 900*Power(gN,4)*Sqr(Qdp)*Sqr(
      QSp) + 900*Power(gN,4)*Sqr(QDxbarp)*Sqr(QSp) + 900*Power(gN,4)*Sqr(QDxp)*Sqr
      (QSp) + 300*Power(gN,4)*Sqr(Qep)*Sqr(QSp) + 900*Power(gN,4)*Sqr(QH1p)*Sqr(
      QSp) + 900*Power(gN,4)*Sqr(QH2p)*Sqr(QSp) + 200*Power(gN,4)*Sqr(QHpbarp)*Sqr
      (QSp) + 200*Power(gN,4)*Sqr(QHpp)*Sqr(QSp) + 600*Power(gN,4)*Sqr(QLp)*Sqr(
      QSp) + 1800*Power(gN,4)*Sqr(QQp)*Sqr(QSp) + 900*Power(gN,4)*Sqr(QH1p)*Sqr(
      Qup) + 900*Power(gN,4)*Sqr(QH2p)*Sqr(Qup) + 900*Power(gN,4)*Sqr(QSp)*Sqr(Qup
      ) - 20*(Sqr(g1) - 5*(8*Sqr(g3) + 3*Sqr(gN)*(Sqr(Qdp) - Sqr(QH1p) + Sqr(QQp))
      ))*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) +
      Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) - 450*(Yd00*(Yd10*(Yd00*Yd10 + Yd01*Yd11
      + Yd02*Yd12) + Yd20*(Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22) + Yd00*(Sqr(Yd00) +
      Sqr(Yd01) + Sqr(Yd02))) + Yd01*(Yd11*(Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) +
      Yd21*(Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22) + Yd01*(Sqr(Yd00) + Sqr(Yd01) + Sqr
      (Yd02))) + Yd02*(Yd12*(Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + Yd22*(Yd00*Yd20
      + Yd01*Yd21 + Yd02*Yd22) + Yd02*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02))) + Yd10*
      (Yd00*(Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + Yd20*(Yd10*Yd20 + Yd11*Yd21 +
      Yd12*Yd22) + Yd10*(Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12))) + Yd11*(Yd01*(Yd00*
      Yd10 + Yd01*Yd11 + Yd02*Yd12) + Yd21*(Yd10*Yd20 + Yd11*Yd21 + Yd12*Yd22) +
      Yd11*(Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12))) + Yd12*(Yd02*(Yd00*Yd10 + Yd01*
      Yd11 + Yd02*Yd12) + Yd22*(Yd10*Yd20 + Yd11*Yd21 + Yd12*Yd22) + Yd12*(Sqr(
      Yd10) + Sqr(Yd11) + Sqr(Yd12))) + Yd20*(Yd00*(Yd00*Yd20 + Yd01*Yd21 + Yd02*
      Yd22) + Yd10*(Yd10*Yd20 + Yd11*Yd21 + Yd12*Yd22) + Yd20*(Sqr(Yd20) + Sqr(
      Yd21) + Sqr(Yd22))) + Yd21*(Yd01*(Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22) + Yd11*
      (Yd10*Yd20 + Yd11*Yd21 + Yd12*Yd22) + Yd21*(Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22
      ))) + Yd22*(Yd02*(Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22) + Yd12*(Yd10*Yd20 +
      Yd11*Yd21 + Yd12*Yd22) + Yd22*(Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)))) + 60*Sqr
      (g1)*(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12)
      + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22)) + 100*Sqr(gN)*Sqr(Qep)*(Sqr(Ye00) + Sqr
      (Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(
      Ye21) + Sqr(Ye22)) - 100*Sqr(gN)*Sqr(QH1p)*(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02
      ) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22)) +
      100*Sqr(gN)*Sqr(QLp)*(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(
      Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22)) - 150*(Ye00*(Ye10*(
      Ye00*Ye10 + Ye01*Ye11 + Ye02*Ye12) + Ye20*(Ye00*Ye20 + Ye01*Ye21 + Ye02*Ye22
      ) + Ye00*(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02))) + Ye01*(Ye11*(Ye00*Ye10 + Ye01
      *Ye11 + Ye02*Ye12) + Ye21*(Ye00*Ye20 + Ye01*Ye21 + Ye02*Ye22) + Ye01*(Sqr(
      Ye00) + Sqr(Ye01) + Sqr(Ye02))) + Ye02*(Ye12*(Ye00*Ye10 + Ye01*Ye11 + Ye02*
      Ye12) + Ye22*(Ye00*Ye20 + Ye01*Ye21 + Ye02*Ye22) + Ye02*(Sqr(Ye00) + Sqr(
      Ye01) + Sqr(Ye02))) + Ye10*(Ye00*(Ye00*Ye10 + Ye01*Ye11 + Ye02*Ye12) + Ye20*
      (Ye10*Ye20 + Ye11*Ye21 + Ye12*Ye22) + Ye10*(Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12
      ))) + Ye11*(Ye01*(Ye00*Ye10 + Ye01*Ye11 + Ye02*Ye12) + Ye21*(Ye10*Ye20 +
      Ye11*Ye21 + Ye12*Ye22) + Ye11*(Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12))) + Ye12*(
      Ye02*(Ye00*Ye10 + Ye01*Ye11 + Ye02*Ye12) + Ye22*(Ye10*Ye20 + Ye11*Ye21 +
      Ye12*Ye22) + Ye12*(Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12))) + Ye20*(Ye00*(Ye00*
      Ye20 + Ye01*Ye21 + Ye02*Ye22) + Ye10*(Ye10*Ye20 + Ye11*Ye21 + Ye12*Ye22) +
      Ye20*(Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22))) + Ye21*(Ye01*(Ye00*Ye20 + Ye01*
      Ye21 + Ye02*Ye22) + Ye11*(Ye10*Ye20 + Ye11*Ye21 + Ye12*Ye22) + Ye21*(Sqr(
      Ye20) + Sqr(Ye21) + Sqr(Ye22))) + Ye22*(Ye02*(Ye00*Ye20 + Ye01*Ye21 + Ye02*
      Ye22) + Ye12*(Ye10*Ye20 + Ye11*Ye21 + Ye12*Ye22) + Ye22*(Sqr(Ye20) + Sqr(
      Ye21) + Sqr(Ye22)))) + 40*Sqr(g1)*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(
      Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)) + 800*Sqr
      (g3)*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12)
      + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)) - 300*Sqr(gN)*Sqr(QH2p)*(Sqr(Yu00) +
      Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(
      Yu21) + Sqr(Yu22)) + 300*Sqr(gN)*Sqr(QQp)*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02)
      + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)) +
      300*Sqr(gN)*Sqr(Qup)*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(
      Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)) + 10*Sqr(Lambdax)*(6*
      Sqr(g1) + 30*Sqr(g2) - 30*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(
      Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(
      Kappa22)) - 20*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(
      Lambda1211)) + 20*Sqr(gN)*Sqr(QH1p) + 20*Sqr(gN)*Sqr(QH2p) - 45*(Sqr(Yd00) +
      Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr
      (Yd21) + Sqr(Yd22)) - 15*(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) +
      Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22)) - 45*(Sqr(Yu00) +
      Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr
      (Yu21) + Sqr(Yu22))) - 450*(Yu00*(Yu10*(Yu00*Yu10 + Yu01*Yu11 + Yu02*Yu12) +
      Yu20*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu00*(Sqr(Yu00) + Sqr(Yu01) +
      Sqr(Yu02))) + Yu01*(Yu11*(Yu00*Yu10 + Yu01*Yu11 + Yu02*Yu12) + Yu21*(Yu00*
      Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu01*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02))) +
      Yu02*(Yu12*(Yu00*Yu10 + Yu01*Yu11 + Yu02*Yu12) + Yu22*(Yu00*Yu20 + Yu01*Yu21
      + Yu02*Yu22) + Yu02*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02))) + Yu10*(Yu00*(Yu00
      *Yu10 + Yu01*Yu11 + Yu02*Yu12) + Yu20*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) +
      Yu10*(Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12))) + Yu11*(Yu01*(Yu00*Yu10 + Yu01*
      Yu11 + Yu02*Yu12) + Yu21*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + Yu11*(Sqr(
      Yu10) + Sqr(Yu11) + Sqr(Yu12))) + Yu12*(Yu02*(Yu00*Yu10 + Yu01*Yu11 + Yu02*
      Yu12) + Yu22*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + Yu12*(Sqr(Yu10) + Sqr(
      Yu11) + Sqr(Yu12))) + Yu20*(Yu00*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu10*
      (Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + Yu20*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22
      ))) + Yu21*(Yu01*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu11*(Yu10*Yu20 +
      Yu11*Yu21 + Yu12*Yu22) + Yu21*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + Yu22*(
      Yu02*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu12*(Yu10*Yu20 + Yu11*Yu21 +
      Yu12*Yu22) + Yu22*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))));


   return twoLoop * deriv;
}

double lowE6SSM_tuning_calculator::deriv_dlead_log_one_loop_Lambdax_dLambdax() const
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

   double deriv = -11.52*Power(g1,4) - 24*Power(g2,4) + 24*Sqr(Kappa00)*Sqr(
      Lambdax) + 24*Sqr(Kappa01)*Sqr(Lambdax) + 24*Sqr(Kappa02)*Sqr(Lambdax) + 24*
      Sqr(Kappa10)*Sqr(Lambdax) + 24*Sqr(Kappa11)*Sqr(Lambdax) + 24*Sqr(Kappa12)*
      Sqr(Lambdax) + 24*Sqr(Kappa20)*Sqr(Lambdax) + 24*Sqr(Kappa21)*Sqr(Lambdax) +
      24*Sqr(Kappa22)*Sqr(Lambdax) + 16*Sqr(Lambda1200)*Sqr(Lambdax) + 16*Sqr(
      Lambda1201)*Sqr(Lambdax) + 16*Sqr(Lambda1210)*Sqr(Lambdax) + 16*Sqr(
      Lambda1211)*Sqr(Lambdax) + 6*Kappa00*(2*(Kappa10*(Kappa00*Kappa10 + Kappa01*
      Kappa11 + Kappa02*Kappa12) + Kappa20*(Kappa00*Kappa20 + Kappa01*Kappa21 +
      Kappa02*Kappa22) + Kappa00*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02))) +
      Kappa00*(-0.26666666666666666*Sqr(g1) - 5.333333333333333*Sqr(g3) + 3*(Sqr(
      Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(
      Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*(Sqr(Lambda1200)
      + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 2*Sqr(Lambdax) - 2*
      Sqr(gN)*Sqr(QDxbarp) - 2*Sqr(gN)*Sqr(QDxp) - 2*Sqr(gN)*Sqr(QSp))) + 6*
      Kappa01*(2*(Kappa11*(Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12) +
      Kappa21*(Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22) + Kappa01*(Sqr
      (Kappa00) + Sqr(Kappa01) + Sqr(Kappa02))) + Kappa01*(-0.26666666666666666*
      Sqr(g1) - 5.333333333333333*Sqr(g3) + 3*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(
      Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(
      Kappa21) + Sqr(Kappa22)) + 2*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(
      Lambda1210) + Sqr(Lambda1211)) + 2*Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QDxbarp) - 2
      *Sqr(gN)*Sqr(QDxp) - 2*Sqr(gN)*Sqr(QSp))) + 6*Kappa02*(2*(Kappa12*(Kappa00*
      Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12) + Kappa22*(Kappa00*Kappa20 +
      Kappa01*Kappa21 + Kappa02*Kappa22) + Kappa02*(Sqr(Kappa00) + Sqr(Kappa01) +
      Sqr(Kappa02))) + Kappa02*(-0.26666666666666666*Sqr(g1) - 5.333333333333333*
      Sqr(g3) + 3*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr
      (Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*(
      Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 2*
      Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QDxbarp) - 2*Sqr(gN)*Sqr(QDxp) - 2*Sqr(gN)*Sqr(
      QSp))) + 6*Kappa10*(2*(Kappa00*(Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*
      Kappa12) + Kappa20*(Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22) +
      Kappa10*(Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12))) + Kappa10*(
      -0.26666666666666666*Sqr(g1) - 5.333333333333333*Sqr(g3) + 3*(Sqr(Kappa00) +
      Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) +
      Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*(Sqr(Lambda1200) + Sqr(
      Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 2*Sqr(Lambdax) - 2*Sqr(gN
      )*Sqr(QDxbarp) - 2*Sqr(gN)*Sqr(QDxp) - 2*Sqr(gN)*Sqr(QSp))) + 6*Kappa11*(2*(
      Kappa01*(Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12) + Kappa21*(
      Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22) + Kappa11*(Sqr(Kappa10)
      + Sqr(Kappa11) + Sqr(Kappa12))) + Kappa11*(-0.26666666666666666*Sqr(g1) -
      5.333333333333333*Sqr(g3) + 3*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) +
      Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) +
      Sqr(Kappa22)) + 2*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr
      (Lambda1211)) + 2*Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QDxbarp) - 2*Sqr(gN)*Sqr(QDxp
      ) - 2*Sqr(gN)*Sqr(QSp))) + 6*Kappa12*(2*(Kappa02*(Kappa00*Kappa10 + Kappa01*
      Kappa11 + Kappa02*Kappa12) + Kappa22*(Kappa10*Kappa20 + Kappa11*Kappa21 +
      Kappa12*Kappa22) + Kappa12*(Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12))) +
      Kappa12*(-0.26666666666666666*Sqr(g1) - 5.333333333333333*Sqr(g3) + 3*(Sqr(
      Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(
      Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*(Sqr(Lambda1200)
      + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 2*Sqr(Lambdax) - 2*
      Sqr(gN)*Sqr(QDxbarp) - 2*Sqr(gN)*Sqr(QDxp) - 2*Sqr(gN)*Sqr(QSp))) + 6*
      Kappa20*(2*(Kappa00*(Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22) +
      Kappa10*(Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22) + Kappa20*(Sqr
      (Kappa20) + Sqr(Kappa21) + Sqr(Kappa22))) + Kappa20*(-0.26666666666666666*
      Sqr(g1) - 5.333333333333333*Sqr(g3) + 3*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(
      Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(
      Kappa21) + Sqr(Kappa22)) + 2*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(
      Lambda1210) + Sqr(Lambda1211)) + 2*Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QDxbarp) - 2
      *Sqr(gN)*Sqr(QDxp) - 2*Sqr(gN)*Sqr(QSp))) + 6*Kappa21*(2*(Kappa01*(Kappa00*
      Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22) + Kappa11*(Kappa10*Kappa20 +
      Kappa11*Kappa21 + Kappa12*Kappa22) + Kappa21*(Sqr(Kappa20) + Sqr(Kappa21) +
      Sqr(Kappa22))) + Kappa21*(-0.26666666666666666*Sqr(g1) - 5.333333333333333*
      Sqr(g3) + 3*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr
      (Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*(
      Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 2*
      Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QDxbarp) - 2*Sqr(gN)*Sqr(QDxp) - 2*Sqr(gN)*Sqr(
      QSp))) + 6*Kappa22*(2*(Kappa02*(Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*
      Kappa22) + Kappa12*(Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22) +
      Kappa22*(Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22))) + Kappa22*(
      -0.26666666666666666*Sqr(g1) - 5.333333333333333*Sqr(g3) + 3*(Sqr(Kappa00) +
      Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) +
      Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*(Sqr(Lambda1200) + Sqr(
      Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 2*Sqr(Lambdax) - 2*Sqr(gN
      )*Sqr(QDxbarp) - 2*Sqr(gN)*Sqr(QDxp) - 2*Sqr(gN)*Sqr(QSp))) + 4*Lambda1200*(
      2*(Lambda1210*(Lambda1200*Lambda1210 + Lambda1201*Lambda1211) + Lambda1200*(
      Sqr(Lambda1200) + Sqr(Lambda1201))) + Lambda1200*(-0.6*Sqr(g1) - 3*Sqr(g2) +
      3*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11)
      + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*(Sqr(
      Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 2*Sqr(
      Lambdax) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QSp)))
      + 4*Lambda1201*(2*(Lambda1211*(Lambda1200*Lambda1210 + Lambda1201*Lambda1211
      ) + Lambda1201*(Sqr(Lambda1200) + Sqr(Lambda1201))) + Lambda1201*(-0.6*Sqr(
      g1) - 3*Sqr(g2) + 3*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(
      Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(
      Kappa22)) + 2*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(
      Lambda1211)) + 2*Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QH2p) -
      2*Sqr(gN)*Sqr(QSp))) + 4*Lambda1210*(2*(Lambda1200*(Lambda1200*Lambda1210 +
      Lambda1201*Lambda1211) + Lambda1210*(Sqr(Lambda1210) + Sqr(Lambda1211))) +
      Lambda1210*(-0.6*Sqr(g1) - 3*Sqr(g2) + 3*(Sqr(Kappa00) + Sqr(Kappa01) + Sqr(
      Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12) + Sqr(Kappa20) + Sqr(
      Kappa21) + Sqr(Kappa22)) + 2*(Sqr(Lambda1200) + Sqr(Lambda1201) + Sqr(
      Lambda1210) + Sqr(Lambda1211)) + 2*Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH1p) - 2*
      Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QSp))) + 4*Lambda1211*(2*(Lambda1201*(
      Lambda1200*Lambda1210 + Lambda1201*Lambda1211) + Lambda1211*(Sqr(Lambda1210)
      + Sqr(Lambda1211))) + Lambda1211*(-0.6*Sqr(g1) - 3*Sqr(g2) + 3*(Sqr(Kappa00
      ) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(Kappa12)
      + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*(Sqr(Lambda1200) + Sqr(
      Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) + 2*Sqr(Lambdax) - 2*Sqr(gN
      )*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QSp))) + Power(gN,3)*(-4*
      gN*Sqr(QH1p) - 4*gN*Sqr(QH2p) - 4*gN*Sqr(QSp))*(9*Sqr(Qdp) + 9*Sqr(QDxbarp)
      + 9*Sqr(QDxp) + 3*Sqr(Qep) + 6*Sqr(QH1p) + 6*Sqr(QH2p) + 2*Sqr(QHpbarp) + 2*
      Sqr(QHpp) + 6*Sqr(QLp) + 18*Sqr(QQp) + 3*Sqr(QSp) + 9*Sqr(Qup)) + 12*Sqr(
      Lambdax)*Sqr(Yd00) + 12*Sqr(Lambdax)*Sqr(Yd01) + 12*Sqr(Lambdax)*Sqr(Yd02) +
      12*Sqr(Lambdax)*Sqr(Yd10) + 12*Sqr(Lambdax)*Sqr(Yd11) + 12*Sqr(Lambdax)*Sqr
      (Yd12) + 12*Sqr(Lambdax)*Sqr(Yd20) + 12*Sqr(Lambdax)*Sqr(Yd21) + 12*Sqr(
      Lambdax)*Sqr(Yd22) + 4*Sqr(Lambdax)*Sqr(Ye00) + 4*Sqr(Lambdax)*Sqr(Ye01) + 4
      *Sqr(Lambdax)*Sqr(Ye02) + 4*Sqr(Lambdax)*Sqr(Ye10) + 4*Sqr(Lambdax)*Sqr(Ye11
      ) + 4*Sqr(Lambdax)*Sqr(Ye12) + 4*Sqr(Lambdax)*Sqr(Ye20) + 4*Sqr(Lambdax)*Sqr
      (Ye21) + 4*Sqr(Lambdax)*Sqr(Ye22) + 2*Ye00*(3*(Ye10*(Ye00*Ye10 + Ye01*Ye11 +
      Ye02*Ye12) + Ye20*(Ye00*Ye20 + Ye01*Ye21 + Ye02*Ye22) + Ye00*(Sqr(Ye00) +
      Sqr(Ye01) + Sqr(Ye02))) + Ye00*(-1.8*Sqr(g1) - 3*Sqr(g2) + Sqr(Lambdax) - 2*
      Sqr(gN)*Sqr(Qep) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QLp) + 3*(Sqr(Yd00) +
      Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr
      (Yd21) + Sqr(Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(
      Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22))) + 2*Ye01*(3*(Ye11*(
      Ye00*Ye10 + Ye01*Ye11 + Ye02*Ye12) + Ye21*(Ye00*Ye20 + Ye01*Ye21 + Ye02*Ye22
      ) + Ye01*(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02))) + Ye01*(-1.8*Sqr(g1) - 3*Sqr(
      g2) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(Qep) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*
      Sqr(QLp) + 3*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) +
      Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr
      (Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(
      Ye22))) + 2*Ye02*(3*(Ye12*(Ye00*Ye10 + Ye01*Ye11 + Ye02*Ye12) + Ye22*(Ye00*
      Ye20 + Ye01*Ye21 + Ye02*Ye22) + Ye02*(Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02))) +
      Ye02*(-1.8*Sqr(g1) - 3*Sqr(g2) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(Qep) - 2*Sqr(
      gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QLp) + 3*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) +
      Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr
      (Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(
      Ye20) + Sqr(Ye21) + Sqr(Ye22))) + 2*Ye10*(3*(Ye00*(Ye00*Ye10 + Ye01*Ye11 +
      Ye02*Ye12) + Ye20*(Ye10*Ye20 + Ye11*Ye21 + Ye12*Ye22) + Ye10*(Sqr(Ye10) +
      Sqr(Ye11) + Sqr(Ye12))) + Ye10*(-1.8*Sqr(g1) - 3*Sqr(g2) + Sqr(Lambdax) - 2*
      Sqr(gN)*Sqr(Qep) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QLp) + 3*(Sqr(Yd00) +
      Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr
      (Yd21) + Sqr(Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(
      Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22))) + 2*Ye11*(3*(Ye01*(
      Ye00*Ye10 + Ye01*Ye11 + Ye02*Ye12) + Ye21*(Ye10*Ye20 + Ye11*Ye21 + Ye12*Ye22
      ) + Ye11*(Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12))) + Ye11*(-1.8*Sqr(g1) - 3*Sqr(
      g2) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(Qep) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*
      Sqr(QLp) + 3*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) +
      Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr
      (Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(
      Ye22))) + 2*Ye12*(3*(Ye02*(Ye00*Ye10 + Ye01*Ye11 + Ye02*Ye12) + Ye22*(Ye10*
      Ye20 + Ye11*Ye21 + Ye12*Ye22) + Ye12*(Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12))) +
      Ye12*(-1.8*Sqr(g1) - 3*Sqr(g2) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(Qep) - 2*Sqr(
      gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QLp) + 3*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) +
      Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr
      (Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(
      Ye20) + Sqr(Ye21) + Sqr(Ye22))) + 6*Yd00*(Yu00*(Yd00*Yu00 + Yd01*Yu01 + Yd02
      *Yu02) + Yu10*(Yd00*Yu10 + Yd01*Yu11 + Yd02*Yu12) + Yu20*(Yd00*Yu20 + Yd01*
      Yu21 + Yd02*Yu22) + 3*(Yd10*(Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + Yd20*(Yd00
      *Yd20 + Yd01*Yd21 + Yd02*Yd22) + Yd00*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02))) +
      Yd00*(-0.4666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) +
      Sqr(Lambdax) - 2*Sqr(gN)*Sqr(Qdp) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QQp
      ) + 3*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12)
      + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) +
      Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22))) + 6*
      Yd01*(Yu01*(Yd00*Yu00 + Yd01*Yu01 + Yd02*Yu02) + Yu11*(Yd00*Yu10 + Yd01*Yu11
      + Yd02*Yu12) + Yu21*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22) + 3*(Yd11*(Yd00*
      Yd10 + Yd01*Yd11 + Yd02*Yd12) + Yd21*(Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22) +
      Yd01*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02))) + Yd01*(-0.4666666666666667*Sqr(g1
      ) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(Qdp
      ) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QQp) + 3*(Sqr(Yd00) + Sqr(Yd01) +
      Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(
      Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(
      Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22))) + 6*Yd02*(Yu02*(Yd00*Yu00 + Yd01
      *Yu01 + Yd02*Yu02) + Yu12*(Yd00*Yu10 + Yd01*Yu11 + Yd02*Yu12) + Yu22*(Yd00*
      Yu20 + Yd01*Yu21 + Yd02*Yu22) + 3*(Yd12*(Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12)
      + Yd22*(Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22) + Yd02*(Sqr(Yd00) + Sqr(Yd01) +
      Sqr(Yd02))) + Yd02*(-0.4666666666666667*Sqr(g1) - 3*Sqr(g2) -
      5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(Qdp) - 2*Sqr(gN)*
      Sqr(QH1p) - 2*Sqr(gN)*Sqr(QQp) + 3*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(
      Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(
      Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20
      ) + Sqr(Ye21) + Sqr(Ye22))) + 6*Yd10*(Yu00*(Yd10*Yu00 + Yd11*Yu01 + Yd12*
      Yu02) + Yu10*(Yd10*Yu10 + Yd11*Yu11 + Yd12*Yu12) + Yu20*(Yd10*Yu20 + Yd11*
      Yu21 + Yd12*Yu22) + 3*(Yd00*(Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + Yd20*(Yd10
      *Yd20 + Yd11*Yd21 + Yd12*Yd22) + Yd10*(Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12))) +
      Yd10*(-0.4666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) +
      Sqr(Lambdax) - 2*Sqr(gN)*Sqr(Qdp) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QQp
      ) + 3*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12)
      + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) +
      Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22))) + 6*
      Yd11*(Yu01*(Yd10*Yu00 + Yd11*Yu01 + Yd12*Yu02) + Yu11*(Yd10*Yu10 + Yd11*Yu11
      + Yd12*Yu12) + Yu21*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22) + 3*(Yd01*(Yd00*
      Yd10 + Yd01*Yd11 + Yd02*Yd12) + Yd21*(Yd10*Yd20 + Yd11*Yd21 + Yd12*Yd22) +
      Yd11*(Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12))) + Yd11*(-0.4666666666666667*Sqr(g1
      ) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(Qdp
      ) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QQp) + 3*(Sqr(Yd00) + Sqr(Yd01) +
      Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(
      Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(
      Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22))) + 6*Yd12*(Yu02*(Yd10*Yu00 + Yd11
      *Yu01 + Yd12*Yu02) + Yu12*(Yd10*Yu10 + Yd11*Yu11 + Yd12*Yu12) + Yu22*(Yd10*
      Yu20 + Yd11*Yu21 + Yd12*Yu22) + 3*(Yd02*(Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12)
      + Yd22*(Yd10*Yd20 + Yd11*Yd21 + Yd12*Yd22) + Yd12*(Sqr(Yd10) + Sqr(Yd11) +
      Sqr(Yd12))) + Yd12*(-0.4666666666666667*Sqr(g1) - 3*Sqr(g2) -
      5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(Qdp) - 2*Sqr(gN)*
      Sqr(QH1p) - 2*Sqr(gN)*Sqr(QQp) + 3*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(
      Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(
      Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20
      ) + Sqr(Ye21) + Sqr(Ye22))) + 6*Yd20*(Yu00*(Yd20*Yu00 + Yd21*Yu01 + Yd22*
      Yu02) + Yu10*(Yd20*Yu10 + Yd21*Yu11 + Yd22*Yu12) + Yu20*(Yd20*Yu20 + Yd21*
      Yu21 + Yd22*Yu22) + 3*(Yd00*(Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22) + Yd10*(Yd10
      *Yd20 + Yd11*Yd21 + Yd12*Yd22) + Yd20*(Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22))) +
      Yd20*(-0.4666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) +
      Sqr(Lambdax) - 2*Sqr(gN)*Sqr(Qdp) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QQp
      ) + 3*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12)
      + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) +
      Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22))) + 6*
      Yd21*(Yu01*(Yd20*Yu00 + Yd21*Yu01 + Yd22*Yu02) + Yu11*(Yd20*Yu10 + Yd21*Yu11
      + Yd22*Yu12) + Yu21*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22) + 3*(Yd01*(Yd00*
      Yd20 + Yd01*Yd21 + Yd02*Yd22) + Yd11*(Yd10*Yd20 + Yd11*Yd21 + Yd12*Yd22) +
      Yd21*(Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22))) + Yd21*(-0.4666666666666667*Sqr(g1
      ) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(Qdp
      ) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QQp) + 3*(Sqr(Yd00) + Sqr(Yd01) +
      Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(
      Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(
      Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22))) + 6*Yd22*(Yu02*(Yd20*Yu00 + Yd21
      *Yu01 + Yd22*Yu02) + Yu12*(Yd20*Yu10 + Yd21*Yu11 + Yd22*Yu12) + Yu22*(Yd20*
      Yu20 + Yd21*Yu21 + Yd22*Yu22) + 3*(Yd02*(Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22)
      + Yd12*(Yd10*Yd20 + Yd11*Yd21 + Yd12*Yd22) + Yd22*(Sqr(Yd20) + Sqr(Yd21) +
      Sqr(Yd22))) + Yd22*(-0.4666666666666667*Sqr(g1) - 3*Sqr(g2) -
      5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(Qdp) - 2*Sqr(gN)*
      Sqr(QH1p) - 2*Sqr(gN)*Sqr(QQp) + 3*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(
      Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(
      Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20
      ) + Sqr(Ye21) + Sqr(Ye22))) + 2*Ye20*(Ye20*(-1.8*Sqr(g1) - 3*Sqr(g2) + Sqr(
      Lambdax) - 2*Sqr(gN)*Sqr(Qep) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QLp) + 3
      *(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) +
      Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr
      (Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22)) + 3*(
      Ye00*(Ye00*Ye20 + Ye01*Ye21 + Ye02*Ye22) + Ye10*(Ye10*Ye20 + Ye11*Ye21 +
      Ye12*Ye22) + Ye20*(Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22)))) + 2*Ye21*(Ye21*(-1.8
      *Sqr(g1) - 3*Sqr(g2) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(Qep) - 2*Sqr(gN)*Sqr(
      QH1p) - 2*Sqr(gN)*Sqr(QLp) + 3*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10
      ) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Sqr(Ye00) +
      Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) + Sqr
      (Ye21) + Sqr(Ye22)) + 3*(Ye01*(Ye00*Ye20 + Ye01*Ye21 + Ye02*Ye22) + Ye11*(
      Ye10*Ye20 + Ye11*Ye21 + Ye12*Ye22) + Ye21*(Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22)
      ))) + 2*Ye22*(Ye22*(-1.8*Sqr(g1) - 3*Sqr(g2) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(
      Qep) - 2*Sqr(gN)*Sqr(QH1p) - 2*Sqr(gN)*Sqr(QLp) + 3*(Sqr(Yd00) + Sqr(Yd01) +
      Sqr(Yd02) + Sqr(Yd10) + Sqr(Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr
      (Yd22)) + Sqr(Ye00) + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(
      Ye12) + Sqr(Ye20) + Sqr(Ye21) + Sqr(Ye22)) + 3*(Ye02*(Ye00*Ye20 + Ye01*Ye21
      + Ye02*Ye22) + Ye12*(Ye10*Ye20 + Ye11*Ye21 + Ye12*Ye22) + Ye22*(Sqr(Ye20) +
      Sqr(Ye21) + Sqr(Ye22)))) + 12*Sqr(Lambdax)*Sqr(Yu00) + 12*Sqr(Lambdax)*Sqr(
      Yu01) + 12*Sqr(Lambdax)*Sqr(Yu02) + 12*Sqr(Lambdax)*Sqr(Yu10) + 12*Sqr(
      Lambdax)*Sqr(Yu11) + 12*Sqr(Lambdax)*Sqr(Yu12) + 12*Sqr(Lambdax)*Sqr(Yu20) +
      12*Sqr(Lambdax)*Sqr(Yu21) + 12*Sqr(Lambdax)*Sqr(Yu22) + 24*Lambdax*(4*Power
      (Lambdax,3) - 0.6*Lambdax*Sqr(g1) - 3*Lambdax*Sqr(g2) + 3*Lambdax*(Sqr(
      Kappa00) + Sqr(Kappa01) + Sqr(Kappa02) + Sqr(Kappa10) + Sqr(Kappa11) + Sqr(
      Kappa12) + Sqr(Kappa20) + Sqr(Kappa21) + Sqr(Kappa22)) + 2*Lambdax*(Sqr(
      Lambda1200) + Sqr(Lambda1201) + Sqr(Lambda1210) + Sqr(Lambda1211)) - 2*
      Lambdax*Sqr(gN)*Sqr(QH1p) - 2*Lambdax*Sqr(gN)*Sqr(QH2p) - 2*Lambdax*Sqr(gN)*
      Sqr(QSp) + 3*Lambdax*(Sqr(Yd00) + Sqr(Yd01) + Sqr(Yd02) + Sqr(Yd10) + Sqr(
      Yd11) + Sqr(Yd12) + Sqr(Yd20) + Sqr(Yd21) + Sqr(Yd22)) + Lambdax*(Sqr(Ye00)
      + Sqr(Ye01) + Sqr(Ye02) + Sqr(Ye10) + Sqr(Ye11) + Sqr(Ye12) + Sqr(Ye20) +
      Sqr(Ye21) + Sqr(Ye22)) + 3*Lambdax*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(
      Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + 6*Yu00
      *(Yd00*(Yd00*Yu00 + Yd01*Yu01 + Yd02*Yu02) + Yd10*(Yd10*Yu00 + Yd11*Yu01 +
      Yd12*Yu02) + Yd20*(Yd20*Yu00 + Yd21*Yu01 + Yd22*Yu02) + 3*(Yu10*(Yu00*Yu10 +
      Yu01*Yu11 + Yu02*Yu12) + Yu20*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu00*(
      Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02))) + Yu00*(-0.8666666666666667*Sqr(g1) - 3*
      Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2
      *Sqr(gN)*Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02
      ) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))))
      + 6*Yu01*(Yd01*(Yd00*Yu00 + Yd01*Yu01 + Yd02*Yu02) + Yd11*(Yd10*Yu00 + Yd11
      *Yu01 + Yd12*Yu02) + Yd21*(Yd20*Yu00 + Yd21*Yu01 + Yd22*Yu02) + 3*(Yu11*(
      Yu00*Yu10 + Yu01*Yu11 + Yu02*Yu12) + Yu21*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22
      ) + Yu01*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02))) + Yu01*(-0.8666666666666667*
      Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*
      Sqr(QH2p) - 2*Sqr(gN)*Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3*(Sqr(Yu00) + Sqr(
      Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21
      ) + Sqr(Yu22)))) + 6*Yu02*(Yd02*(Yd00*Yu00 + Yd01*Yu01 + Yd02*Yu02) + Yd12*(
      Yd10*Yu00 + Yd11*Yu01 + Yd12*Yu02) + Yd22*(Yd20*Yu00 + Yd21*Yu01 + Yd22*Yu02
      ) + 3*(Yu12*(Yu00*Yu10 + Yu01*Yu11 + Yu02*Yu12) + Yu22*(Yu00*Yu20 + Yu01*
      Yu21 + Yu02*Yu22) + Yu02*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02))) + Yu02*(
      -0.8666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(
      Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3
      *(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) +
      Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) + 6*Yu10*(Yd00*(Yd00*Yu10 + Yd01*Yu11 +
      Yd02*Yu12) + Yd10*(Yd10*Yu10 + Yd11*Yu11 + Yd12*Yu12) + Yd20*(Yd20*Yu10 +
      Yd21*Yu11 + Yd22*Yu12) + 3*(Yu00*(Yu00*Yu10 + Yu01*Yu11 + Yu02*Yu12) + Yu20*
      (Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + Yu10*(Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12
      ))) + Yu10*(-0.8666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(
      g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QQp) - 2*Sqr(gN)*
      Sqr(Qup) + 3*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) +
      Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) + 6*Yu11*(Yd01*(Yd00*Yu10 +
      Yd01*Yu11 + Yd02*Yu12) + Yd11*(Yd10*Yu10 + Yd11*Yu11 + Yd12*Yu12) + Yd21*(
      Yd20*Yu10 + Yd21*Yu11 + Yd22*Yu12) + 3*(Yu01*(Yu00*Yu10 + Yu01*Yu11 + Yu02*
      Yu12) + Yu21*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + Yu11*(Sqr(Yu10) + Sqr(
      Yu11) + Sqr(Yu12))) + Yu11*(-0.8666666666666667*Sqr(g1) - 3*Sqr(g2) -
      5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*
      Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(
      Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) + 6*
      Yu12*(Yd02*(Yd00*Yu10 + Yd01*Yu11 + Yd02*Yu12) + Yd12*(Yd10*Yu10 + Yd11*Yu11
      + Yd12*Yu12) + Yd22*(Yd20*Yu10 + Yd21*Yu11 + Yd22*Yu12) + 3*(Yu02*(Yu00*
      Yu10 + Yu01*Yu11 + Yu02*Yu12) + Yu22*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) +
      Yu12*(Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12))) + Yu12*(-0.8666666666666667*Sqr(g1
      ) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(
      QH2p) - 2*Sqr(gN)*Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3*(Sqr(Yu00) + Sqr(Yu01) +
      Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr
      (Yu22)))) + 6*Yu20*(Yd00*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22) + Yd10*(Yd10*
      Yu20 + Yd11*Yu21 + Yd12*Yu22) + Yd20*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22) + 3
      *(Yu00*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu10*(Yu10*Yu20 + Yu11*Yu21 +
      Yu12*Yu22) + Yu20*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22))) + Yu20*(
      -0.8666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) + Sqr(
      Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3
      *(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) + Sqr(Yu12) +
      Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) + 6*Yu21*(Yd01*(Yd00*Yu20 + Yd01*Yu21 +
      Yd02*Yu22) + Yd11*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22) + Yd21*(Yd20*Yu20 +
      Yd21*Yu21 + Yd22*Yu22) + 3*(Yu01*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + Yu11*
      (Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + Yu21*(Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22
      ))) + Yu21*(-0.8666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(
      g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*Sqr(QQp) - 2*Sqr(gN)*
      Sqr(Qup) + 3*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(Yu10) + Sqr(Yu11) +
      Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) + 6*Yu22*(Yd02*(Yd00*Yu20 +
      Yd01*Yu21 + Yd02*Yu22) + Yd12*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22) + Yd22*(
      Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22) + 3*(Yu02*(Yu00*Yu20 + Yu01*Yu21 + Yu02*
      Yu22) + Yu12*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + Yu22*(Sqr(Yu20) + Sqr(
      Yu21) + Sqr(Yu22))) + Yu22*(-0.8666666666666667*Sqr(g1) - 3*Sqr(g2) -
      5.333333333333333*Sqr(g3) + Sqr(Lambdax) - 2*Sqr(gN)*Sqr(QH2p) - 2*Sqr(gN)*
      Sqr(QQp) - 2*Sqr(gN)*Sqr(Qup) + 3*(Sqr(Yu00) + Sqr(Yu01) + Sqr(Yu02) + Sqr(
      Yu10) + Sqr(Yu11) + Sqr(Yu12) + Sqr(Yu20) + Sqr(Yu21) + Sqr(Yu22)))) + Sqr(
      -0.6*Power(g1,2) - 3*Power(g2,2) + 3*(Power(Kappa00,2) + Power(Kappa01,2) +
      Power(Kappa02,2) + Power(Kappa10,2) + Power(Kappa11,2) + Power(Kappa12,2) +
      Power(Kappa20,2) + Power(Kappa21,2) + Power(Kappa22,2)) + 2*(Power(
      Lambda1200,2) + Power(Lambda1201,2) + Power(Lambda1210,2) + Power(Lambda1211
      ,2)) + 12*Power(Lambdax,2) - 2*Power(gN,2)*Power(QH1p,2) - 2*Power(gN,2)*
      Power(QH2p,2) - 2*Power(gN,2)*Power(QSp,2) + 3*(Power(Yd00,2) + Power(Yd01,2
      ) + Power(Yd02,2) + Power(Yd10,2) + Power(Yd11,2) + Power(Yd12,2) + Power(
      Yd20,2) + Power(Yd21,2) + Power(Yd22,2)) + Power(Ye00,2) + Power(Ye01,2) +
      Power(Ye02,2) + Power(Ye10,2) + Power(Ye11,2) + Power(Ye12,2) + Power(Ye20,2
      ) + Power(Ye21,2) + Power(Ye22,2) + 3*(Power(Yu00,2) + Power(Yu01,2) + Power
      (Yu02,2) + Power(Yu10,2) + Power(Yu11,2) + Power(Yu12,2) + Power(Yu20,2) +
      Power(Yu21,2) + Power(Yu22,2)));


   return twoLoop * deriv;
}

} // namespace flexiblesusy 
