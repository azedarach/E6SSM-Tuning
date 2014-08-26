/**
 * @file lowMSSM_tuning_calculator.cpp
 * @brief implementation of the lowMSSM tuning calculator class
 *
 * Contains the definition of the lowMSSM tuning calculator class
 * methods which compute the fine tuning in the model
 */

#include "lowMSSM_tuning_calculator.hpp"
#include "pv.hpp"

namespace flexiblesusy {
  
double lowMSSM_tuning_calculator::gbar() const
{
   double g1  = model.get_g1();
   double g2  = model.get_g2();
   double thW = model.ThetaW();

   return (g2*Cos(thW) + 0.7745966692414834*g1*Sin(thW));
}

double lowMSSM_tuning_calculator::MQQ2() const
{
   double mq2_22 = model.get_mq2(2,2);
   double mu2_22 = model.get_mu2(2,2);
   double g1 = model.get_g1();
   double g2 = model.get_g2();
   double vd = model.get_vd();
   double vu = model.get_vu();

   return (mq2_22 - mu2_22 + 0.125*(Sqr(vd)-Sqr(vu))*(Sqr(g2)-Sqr(g1)));
}

double lowMSSM_tuning_calculator::RQQ() const
{
   double mqqsq = MQQ2();
   double g1 = model.get_g1();
   double g2 = model.get_g2();

   return mqqsq*(Sqr(g2)-Sqr(g1));
}

double lowMSSM_tuning_calculator::rt() const
{
   double mqq4 = Sqr(MQQ2());
   double vd = model.get_vd();
   double vu = model.get_vu();
   double yt = model.get_Yu(2,2);
   double at = model.get_TYu(2,2);
   double mu = model.get_Mu();

   return (mqq4 + 2.0*Sqr(vu)*Sqr(at-yt*mu*vd/vu));
}

// double deriv_d2DeltaV_dvd_dvd() const
// {
//   double tmp_1;
//   double tmp_2;
//   double tmp_3;
//   double tmp_4;
// }

// double deriv_d2DeltaV_dvu_dvu() const
// {

// }

// double deriv_d2DeltaV_dvu_dvd() const
// {

// }

} // namespace flexiblesusy
