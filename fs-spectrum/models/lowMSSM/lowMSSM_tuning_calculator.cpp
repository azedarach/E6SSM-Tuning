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

      return (g2 * Cos(thW) + 0.7745966692414834 * g1 * Sin(thW));
   }

   double lowMSSM_tuning_calculator::MFtop_DRbar() const
   {
      double yt = model.get_Yu(2, 2);
      double vu = model.get_vu();

      return 0.7071067811865475 * yt * vu;
   }

   double lowMSSM_tuning_calculator::stop_mass_matrix_LL_entry() const
   {
      double mq2_22 = model.get_mq2(2, 2);
      double mt = MFtop_DRbar();
      double vd = model.get_vd();
      double vu = model.get_vu();
      double g1 = model.get_g1();
      double g2 = model.get_g2();
      double delta_ul = 0.125 * (Sqr(g2) - 0.2 * Sqr(g1)) * (Sqr(vd) - Sqr(vu));

      return (mq2_22 + Sqr(mt) + delta_ul);
   }

   double lowMSSM_tuning_calculator::stop_mass_matrix_RR_entry() const
   {
      double mu2_22 = model.get_mu2(2, 2);
      double mt = MFtop_DRbar();
      double vd = model.get_vd();
      double vu = model.get_vu();
      double g1 = model.get_g1();
      double delta_ur = 0.1 * Sqr(g1) * (Sqr(vd) - Sqr(vu));

      return (mu2_22 + Sqr(mt) + delta_ur);
   }

   double lowMSSM_tuning_calculator::stop_mass_matrix_LR_entry() const
   {
      // DH:: currently assumes all parameters are *real*
      double at = model.get_TYu(2, 2);
      double yt = model.get_Yu(2, 2);
      double mu = model.get_Mu();
      double vd = model.get_vd();
      double vu = model.get_vu();

      return (0.7071067811865475 * (at * vu - mu * yt * vd));
   }

   double lowMSSM_tuning_calculator::MQQ2() const
   {
      double mass_matrix_stop_LL = stop_mass_matrix_LL_entry();
      double mass_matrix_stop_RR = stop_mass_matrix_RR_entry();

      return (mass_matrix_stop_LL - mass_matrix_stop_RR);
   }

   double lowMSSM_tuning_calculator::RQQ() const
   {
      double g1 = model.get_g1();
      double g2 = model.get_g2();
      double mqqsq = MQQ2();

      return (Sqr(g2) - Sqr(g1)) * mqqsq;
   }

   double lowMSSM_tuning_calculator::stop_discriminant() const
   {
      double mqqsq = MQQ2();
      double mass_matrix_stop_LR = stop_mass_matrix_LR_entry();

      return (Sqr(mqqsq) + 4.0 * AbsSqr(mass_matrix_stop_LR));
   }

   void lowMSSM_tuning_calculator::calculate_MStop()
   {
      // For now just use the explicit solutions for the mass
      // eigenvalues; avoids any mass ordering issues.
      double stop_mass_matrix_trace = stop_mass_matrix_LL_entry()
                                      + stop_mass_matrix_RR_entry();
      double disc = stop_discriminant();

      MStop(0) = 0.5 * (stop_mass_matrix_trace - Sqrt(disc));
      MStop(1) = 0.5 * (stop_mass_matrix_trace + Sqrt(disc));

      /// DH:: Check that this is appropriate action
      if (MStop.minCoeff() < 0.) model.get_problems().flag_tachyon(lowMSSM_info::Su);

      MStop = AbsSqrt(MStop);
   }

   double lowMSSM_tuning_calculator::deriv_dMFtop2_dparam(lowMSSM_info::Parameters p) const
   {
      switch (p) {
      case lowMSSM_info::vu: {
         return Sqr(model.get_Yu(2, 2)) * model.get_vu();
      }
      default: {
         return 0.;
      }   
      }
   }

   double lowMSSM_tuning_calculator::deriv_dMStop2_dvd(stop_mass which_stop) const
   {
      double mu = model.get_Mu();
      double vd = model.get_vd();
      double g1 = model.get_g1();
      double g2 = model.get_g2();
      double yt = model.get_Yu(2, 2);
      double rqq_val = RQQ();
      double stop_mixing = 1.4142135623730951 * stop_mass_matrix_LR_entry();
      double rt = stop_discriminant();

      double split = (0.25 * rqq_val - 2.0 * yt * mu * stop_mixing / vd) / Sqrt(rt);

      double deriv = 0.25 * Sqr(g2) + 0.15 * Sqr(g1);

      if (which_stop == stop_mass::mstop_1) {
         return 0.5 * vd * (deriv - split);
      } else {
      return 0.5 * vd * (deriv + split);
      }
   }

   double lowMSSM_tuning_calculator::deriv_dMStop2_dvu(stop_mass which_stop) const
   {
      double at = model.get_TYu(2, 2);
      double vu = model.get_vu();
      double g1 = model.get_g1();
      double g2 = model.get_g2();
      double yt = model.get_Yu(2, 2);
      double rqq_val = RQQ();
      double stop_mixing = 1.4142135623730951 * stop_mass_matrix_LR_entry();
      double rt = stop_discriminant();

      double split = (2.0 * at * stop_mixing / vu - 0.25 * rqq_val) / Sqrt(rt);

      double deriv = 2.0 * Sqr(yt) - 0.25 * Sqr(g2) - 0.15 * Sqr(g1);

      if (which_stop == stop_mass::mstop_1) {
         return 0.5 * vu * (deriv - split);
      } else {
         return 0.5 * vu * (deriv + split);
      }
   }

   double lowMSSM_tuning_calculator::deriv_dMStop2_dmq222(stop_mass which_stop) const
   {
      double mqqsq = MQQ2();
      double rt = stop_discriminant();
      if (which_stop == stop_mass::mstop_1) {
         return 0.5 * (1.0 - mqqsq / Sqrt(rt));
      } else {
         return 0.5 * (1.0 + mqqsq / Sqrt(rt));
      }
   }

   double lowMSSM_tuning_calculator::deriv_dMStop2_dmu222(stop_mass which_stop) const
   {
      double mqqsq = MQQ2();
      double rt = stop_discriminant();
      if (which_stop == stop_mass::mstop_1) {
         return 0.5 * (1.0 + mqqsq / Sqrt(rt));
      } else {
         return 0.5 * (1.0 - mqqsq / Sqrt(rt));
      }
   }

   double lowMSSM_tuning_calculator::deriv_dMStop2_dMu(stop_mass which_stop) const
   {
      double yt = model.get_Yu(2, 2);
      double vd = model.get_vd();
      double stop_mixing = 1.4142135623730951 * stop_mass_matrix_LR_entry();
      double rt = stop_discriminant();
      double deriv = yt * vd * stop_mixing / Sqrt(rt);
      if (which_stop == stop_mass::mstop_1) {
         return deriv;
      } else {
         return -deriv; 
      }
   }

   // Note this equals the derivative w.r.t. A_t up to a factor of y_t 
   double lowMSSM_tuning_calculator::deriv_dMStop2_dTYu22(stop_mass which_stop) const
   {
      double vu = model.get_vu();
      double stop_mixing = 1.4142135623730951 * stop_mass_matrix_LR_entry();
      double rt = stop_discriminant();
      double deriv = vu * stop_mixing / Sqrt(rt);
      if (which_stop == stop_mass::mstop_1) {
         return -deriv;
      } else {
         return deriv;
      }
   }

   double lowMSSM_tuning_calculator::deriv_dMStop2_dparam(stop_mass which_stop, lowMSSM_info::Parameters p) const
   {
      switch(p) {
      case lowMSSM_info::mq222: {
         return deriv_dMStop2_dmq222(which_stop);
      }
      case lowMSSM_info::mu222: {
         return deriv_dMStop2_dmu222(which_stop);
      }
      case lowMSSM_info::Mu: {
         return deriv_dMStop2_dMu(which_stop);
      }
      case lowMSSM_info::vd: {
         return deriv_dMStop2_dvd(which_stop);
      }
      case lowMSSM_info::vu: {
         return deriv_dMStop2_dvu(which_stop);
      }
      case lowMSSM_info::TYu22: {
         return deriv_dMStop2_dTYu22(which_stop);
      }
      default: {
         return 0.;
      }
      }
   }

   double lowMSSM_tuning_calculator::deriv_d2MFtop2_dparam_dparam(lowMSSM_info::Parameters p1, 
                                                                  lowMSSM_info::Parameters p2) const
   {
      if (p1 == lowMSSM_info::vu && p2 == lowMSSM_info::vu) {
         return Sqr(model.get_Yu(2, 2));
      } else {
         return 0.;
      }
   }

   double lowMSSM_tuning_calculator::deriv_d2MStop2_dvd_dvd(stop_mass which_stop) const
   {
      double mu = model.get_Mu();
      double vd = model.get_vd();
      double yt = model.get_Yu(2, 2);
      double g1 = model.get_g1();
      double g2 = model.get_g2();

      double gbar_val = gbar();
      double stop_mixing = 1.4142135623730951 * stop_mass_matrix_LR_entry();
      double RQQ_val = RQQ();
      double rt = stop_discriminant();

      double split = 0.03125 * ((4.0 * RQQ_val + Sqr((g2 * g2 - g1 * g1) * vd) + 32.0 * Sqr(yt * mu)) / Sqrt(rt)
                                - Sqr(vd * RQQ_val - 8.0 * yt * mu * stop_mixing) / (rt * Sqrt(rt)));

      if (which_stop == stop_mass::mstop_1) {
         return 0.125 * Sqr(gbar_val) - split;
      } else {
         return 0.125 * Sqr(gbar_val) + split;
      }
   }

   double lowMSSM_tuning_calculator::deriv_d2MStop2_dvu_dvu(stop_mass which_stop) const
   {
      double vu = model.get_vu();
      double yt = model.get_Yu(2, 2);
      double at = model.get_TYu(2, 2);
      double g1 = model.get_g1();
      double g2 = model.get_g2();

      double gbar_val = gbar();
      double stop_mixing = 1.4142135623730951 * stop_mass_matrix_LR_entry();
      double RQQ_val = RQQ();
      double rt = stop_discriminant();

      double split = (0.03125 * Sqr(vu * (g2 * g2 - g1 * g1)) - 0.125 * RQQ_val + Sqr(at)
                      - 0.03125 * Sqr(vu) * Sqr(8.0 * at * stop_mixing / vu - RQQ_val) / rt) / Sqrt(rt);

      if (which_stop == stop_mass::mstop_1) {
         return Sqr(yt) - 0.125 * Sqr(gbar_val) - split;
      } else {
         return Sqr(yt) - 0.125 * Sqr(gbar_val) + split;
      }
   }

   double lowMSSM_tuning_calculator::deriv_d2MStop2_dvd_dvu(stop_mass which_stop) const
   {
      double mu = model.get_Mu();
      double vd = model.get_vd();
      double vu = model.get_vu();
      double yt = model.get_Yu(2, 2);
      double at = model.get_TYu(2, 2);
      double g1 = model.get_g1();
      double g2 = model.get_g2();

      double stop_mixing = 1.4142135623730951 * stop_mass_matrix_LR_entry();
      double RQQ_val = RQQ();
      double rt = stop_discriminant();

      double deriv = 0.03125 * vu * vd * Sqr(g2 * g2 - g1 * g1) + yt * mu * at
         + 2.0 * vu * vd * ((0.125 * RQQ_val - yt * mu * stop_mixing / vd)
                            * (at * stop_mixing / vu - 0.125 * RQQ_val)) / rt;

      if (which_stop == stop_mass::mstop_1) {
         return deriv / Sqrt(rt);
      } else {
         return -deriv / Sqrt(rt);
      }
   }

   double lowMSSM_tuning_calculator::deriv_d2MStop2_dmq222_dvd(stop_mass which_stop) const
   {
      double mu = model.get_Mu();
      double vd = model.get_vd();
      double yt = model.get_Yu(2, 2);
      double g1 = model.get_g1();
      double g2 = model.get_g2();

      double stop_mixing = 1.4142135623730951 * stop_mass_matrix_LR_entry();
      double mqqsq = MQQ2();
      double RQQ_val = RQQ();
      double rt = stop_discriminant();

      double deriv = mqqsq * (0.125 * RQQ_val - yt * mu * stop_mixing / vd) / rt
         - 0.125 * ( g2 * g2 - g1 * g1);

      if (which_stop == stop_mass::mstop_1) {
         return vd * deriv / Sqrt(rt);
      } else {
         return -vd * deriv / Sqrt(rt);
      }
   }

   double lowMSSM_tuning_calculator::deriv_d2MStop2_dmu222_dvd(stop_mass which_stop) const
   {
      double mu = model.get_Mu();
      double vd = model.get_vd();
      double yt = model.get_Yu(2, 2);
      double g1 = model.get_g1();
      double g2 = model.get_g2();

      double stop_mixing = 1.4142135623730951 * stop_mass_matrix_LR_entry();
      double mqqsq = MQQ2();
      double RQQ_val = RQQ();
      double rt = stop_discriminant();

      double deriv = mqqsq * (0.125 * RQQ_val - yt * mu * stop_mixing / vd) / rt
         - 0.125 * ( g2 * g2 - g1 * g1);

      if (which_stop == stop_mass::mstop_1) {
         return -vd * deriv / Sqrt(rt);
      } else {
         return vd * deriv / Sqrt(rt);
      }
   }

   double lowMSSM_tuning_calculator::deriv_d2MStop2_dMu_dvd(stop_mass which_stop) const
   {
      double mu = model.get_Mu();
      double vd = model.get_vd();
      double yt = model.get_Yu(2, 2);

      double stop_mixing = 1.4142135623730951 * stop_mass_matrix_LR_entry();
      double RQQ_val = RQQ();
      double rt = stop_discriminant();

      double deriv = 2.0 * yt * vd * stop_mixing * (0.125 * RQQ_val - yt * mu * stop_mixing / vd) / rt
         - yt * stop_mixing / vd + yt * yt * mu;

      if (which_stop == stop_mass::mstop_1) {
         return -vd * deriv / Sqrt(rt);
      } else {
         return vd * deriv / Sqrt(rt);
      }
   }

   double lowMSSM_tuning_calculator::deriv_d2MStop2_dTYu22_dvd(stop_mass which_stop) const
   {
      double mu = model.get_Mu();
      double vd = model.get_vd();
      double vu = model.get_vu();
      double yt = model.get_Yu(2, 2);

      double stop_mixing = 1.4142135623730951 * stop_mass_matrix_LR_entry();
      double RQQ_val = RQQ();
      double rt = stop_discriminant();

      double deriv = 2.0 * vu * stop_mixing * (0.125 * RQQ_val - yt * mu * stop_mixing / vd) / rt
         + yt * vu * mu / vd;

      if (which_stop == stop_mass::mstop_1) {
         return vd * deriv / Sqrt(rt);
      } else {
         return -vd * deriv / Sqrt(rt);
      }
   }

   double lowMSSM_tuning_calculator::deriv_d2MStop2_dmq222_dvu(stop_mass which_stop) const
   {
      double at = model.get_TYu(2, 2);
      double vu = model.get_vu();
      double g1 = model.get_g1();
      double g2 = model.get_g2();

      double stop_mixing = 1.4142135623730951 * stop_mass_matrix_LR_entry();
      double mqqsq = MQQ2();
      double RQQ_val = RQQ();
      double rt = stop_discriminant();

      double deriv = mqqsq * (at * stop_mixing / vu - 0.125 * RQQ_val) / rt
         + 0.125 * ( g2 * g2 - g1 * g1);

      if (which_stop == stop_mass::mstop_1) {
         return vu * deriv / Sqrt(rt);
      } else {
         return -vu * deriv / Sqrt(rt);
      }
   }

   double lowMSSM_tuning_calculator::deriv_d2MStop2_dmu222_dvu(stop_mass which_stop) const
   {
      double at = model.get_TYu(2, 2);
      double vu = model.get_vu();
      double g1 = model.get_g1();
      double g2 = model.get_g2();

      double stop_mixing = 1.4142135623730951 * stop_mass_matrix_LR_entry();
      double mqqsq = MQQ2();
      double RQQ_val = RQQ();
      double rt = stop_discriminant();

      double deriv = mqqsq * (at * stop_mixing / vu - 0.125 * RQQ_val) / rt
         + 0.125 * ( g2 * g2 - g1 * g1);

      if (which_stop == stop_mass::mstop_1) {
         return -vu * deriv / Sqrt(rt);
      } else {
         return vu * deriv / Sqrt(rt);
      }
   }

   double lowMSSM_tuning_calculator::deriv_d2MStop2_dMu_dvu(stop_mass which_stop) const
   {
      double at = model.get_TYu(2, 2);
      double vd = model.get_vd();
      double vu = model.get_vu();
      double yt = model.get_Yu(2, 2);

      double stop_mixing = 1.4142135623730951 * stop_mass_matrix_LR_entry();
      double RQQ_val = RQQ();
      double rt = stop_discriminant();

      double deriv = 2.0 * yt * vd * stop_mixing * (at * stop_mixing / vu - 0.125 * RQQ_val) / rt
         - yt * at * vd / vu;

      if (which_stop == stop_mass::mstop_1) {
         return -vu * deriv / Sqrt(rt);
      } else {
         return vu * deriv / Sqrt(rt);
      }
   }

   double lowMSSM_tuning_calculator::deriv_d2MStop2_dTYu22_dvu(stop_mass which_stop) const
   {
      double at = model.get_TYu(2, 2);
      double vu = model.get_vu();

      double stop_mixing = 1.4142135623730951 * stop_mass_matrix_LR_entry();
      double RQQ_val = RQQ();
      double rt = stop_discriminant();

      double deriv = 2.0 * vu * stop_mixing * (at * stop_mixing / vu - 0.125 * RQQ_val) / rt
         - stop_mixing / vu - at;

      if (which_stop == stop_mass::mstop_1) {
         return vu * deriv / Sqrt(rt);
      } else {
         return -vu * deriv / Sqrt(rt);
      }
   }

   double lowMSSM_tuning_calculator::deriv_d2MStop2_dparam_dparam(stop_mass which_stop, lowMSSM_info::Parameters p1, 
                                                                 lowMSSM_info::Parameters p2) const
   {
      if (p1 == lowMSSM_info::vd && p2 == lowMSSM_info::vd) {
         return deriv_d2MStop2_dvd_dvd(which_stop);
      } else if (p1 == lowMSSM_info::vu && p2 == lowMSSM_info::vu) {
         return deriv_d2MStop2_dvu_dvu(which_stop);
      } else if ((p1 == lowMSSM_info::vu && p2 == lowMSSM_info::vd) ||
                 (p1 == lowMSSM_info::vd && p2 == lowMSSM_info::vu)) {
         return deriv_d2MStop2_dvd_dvu(which_stop);
      } else if ((p1 == lowMSSM_info::mq222 && p2 == lowMSSM_info::vd) ||
                 (p1 == lowMSSM_info::vd && p2 == lowMSSM_info::mq222)) {
         return deriv_d2MStop2_dmq222_dvd(which_stop);
      } else if ((p1 == lowMSSM_info::mu222 && p2 == lowMSSM_info::vd) ||
                 (p1 == lowMSSM_info::vd && p2 == lowMSSM_info::mu222)) {
         return deriv_d2MStop2_dmu222_dvd(which_stop);
      } else if ((p1 == lowMSSM_info::Mu && p2 == lowMSSM_info::vd) ||
                 (p1 == lowMSSM_info::vd && p2 == lowMSSM_info::Mu)) {
         return deriv_d2MStop2_dMu_dvd(which_stop);
      } else if ((p1 == lowMSSM_info::TYu22 && p2 == lowMSSM_info::vd) ||
                 (p1 == lowMSSM_info::vd && p2 == lowMSSM_info::TYu22)) {
         return deriv_d2MStop2_dTYu22_dvd(which_stop);
      } else if ((p1 == lowMSSM_info::mq222 && p2 == lowMSSM_info::vu) ||
                 (p1 == lowMSSM_info::vu && p2 == lowMSSM_info::mq222)) {
         return deriv_d2MStop2_dmq222_dvu(which_stop);
      } else if ((p1 == lowMSSM_info::mu222 && p2 == lowMSSM_info::vu) ||
                 (p1 == lowMSSM_info::vu && p2 == lowMSSM_info::mu222)) {
         return deriv_d2MStop2_dmu222_dvu(which_stop);
      } else if ((p1 == lowMSSM_info::Mu && p2 == lowMSSM_info::vu) ||
                 (p1 == lowMSSM_info::vu && p2 == lowMSSM_info::Mu)) {
         return deriv_d2MStop2_dMu_dvu(which_stop);
      } else if ((p1 == lowMSSM_info::TYu22 && p2 == lowMSSM_info::vu) ||
                 (p1 == lowMSSM_info::vu && p2 == lowMSSM_info::TYu22)) {
         return deriv_d2MStop2_dTYu22_dvu(which_stop);
      } else {
         return 0.; //< while not actually zero, all other derivatives are not needed for now
      }
   }

   double lowMSSM_tuning_calculator::deriv_d2DeltaV_dvd_dvd() const
   {
      double scale = model.get_scale();

      double dMStop20_dvd = deriv_dMStop2_dparam(stop_mass::mstop_1, lowMSSM_info::vd);
      double d2MStop20_dvd_dvd = deriv_d2MStop2_dparam_dparam(stop_mass::mstop_1, lowMSSM_info::vd, lowMSSM_info::vd);
      double dMStop21_dvd = deriv_dMStop2_dparam(stop_mass::mstop_2, lowMSSM_info::vd);
      double d2MStop21_dvd_dvd = deriv_d2MStop2_dparam_dparam(stop_mass::mstop_2, lowMSSM_info::vd, lowMSSM_info::vd);

      double a0_mstop1 = passarino_veltman::ReA0(Sqr(MStop(0)), Sqr(scale));
      double a0_mstop2 = passarino_veltman::ReA0(Sqr(MStop(1)), Sqr(scale));

      double tmp_1 = Sqr(dMStop20_dvd) * Log(Sqr(MStop(0) / scale));
      double tmp_2 = -a0_mstop1 * d2MStop20_dvd_dvd;
      double tmp_3 = Sqr(dMStop21_dvd) * Log(Sqr(MStop(1) / scale));
      double tmp_4 = -a0_mstop2 * d2MStop21_dvd_dvd;

      return (3.0 * oneOver16PiSqr * (tmp_1 + tmp_2 + tmp_3 + tmp_4));
   }

   double lowMSSM_tuning_calculator::deriv_d2DeltaV_dvu_dvu() const
   {
      double scale = model.get_scale();

      double dMStop20_dvu = deriv_dMStop2_dparam(stop_mass::mstop_1, lowMSSM_info::vu);
      double d2MStop20_dvu_dvu = deriv_d2MStop2_dparam_dparam(stop_mass::mstop_1, lowMSSM_info::vu, lowMSSM_info::vu);
      double dMStop21_dvu = deriv_dMStop2_dparam(stop_mass::mstop_2, lowMSSM_info::vu);
      double d2MStop21_dvu_dvu = deriv_d2MStop2_dparam_dparam(stop_mass::mstop_2, lowMSSM_info::vu, lowMSSM_info::vu);
      double dMFtop2_dvu = deriv_dMFtop2_dparam(lowMSSM_info::vu);
      double d2MFtop2_dvu_dvu = deriv_d2MFtop2_dparam_dparam(lowMSSM_info::vu, lowMSSM_info::vu); 

      double mt = MFtop_DRbar();

      double a0_mstop1 = passarino_veltman::ReA0(Sqr(MStop(0)), Sqr(scale));
      double a0_mstop2 = passarino_veltman::ReA0(Sqr(MStop(1)), Sqr(scale));
      double a0_mtop = passarino_veltman::ReA0(Sqr(mt), Sqr(scale));

      double tmp_1 = Sqr(dMStop20_dvu) * Log(Sqr(MStop(0) / scale));
      double tmp_2 = -a0_mstop1 * d2MStop20_dvu_dvu;
      double tmp_3 = Sqr(dMStop21_dvu) * Log(Sqr(MStop(1) / scale));
      double tmp_4 = -a0_mstop2 * d2MStop21_dvu_dvu;
      double tmp_5 = -2.0 * Sqr(dMFtop2_dvu) * Log(Sqr(mt / scale));
      double tmp_6 = 2.0 * a0_mtop * d2MFtop2_dvu_dvu;

      return (3.0 * oneOver16PiSqr * (tmp_1 + tmp_2 + tmp_3 + tmp_4 + tmp_5 + tmp_6));
   }

   double lowMSSM_tuning_calculator::deriv_d2DeltaV_dvu_dvd() const
   {
      double scale = model.get_scale();

      double dMStop20_dvd = deriv_dMStop2_dparam(stop_mass::mstop_1, lowMSSM_info::vd);
      double dMStop20_dvu = deriv_dMStop2_dparam(stop_mass::mstop_1, lowMSSM_info::vu);
      double d2MStop20_dvd_dvu = deriv_d2MStop2_dparam_dparam(stop_mass::mstop_1, lowMSSM_info::vd, lowMSSM_info::vu);
      double dMStop21_dvd = deriv_dMStop2_dparam(stop_mass::mstop_2, lowMSSM_info::vd);
      double dMStop21_dvu = deriv_dMStop2_dparam(stop_mass::mstop_2, lowMSSM_info::vu);
      double d2MStop21_dvd_dvu = deriv_d2MStop2_dparam_dparam(stop_mass::mstop_2, lowMSSM_info::vd, lowMSSM_info::vu);

      double a0_mstop1 = passarino_veltman::ReA0(Sqr(MStop(0)), Sqr(scale));
      double a0_mstop2 = passarino_veltman::ReA0(Sqr(MStop(1)), Sqr(scale));

      double tmp_1 = dMStop20_dvd * dMStop20_dvu * Log(Sqr(MStop(0) / scale));
      double tmp_2 = -a0_mstop1 * d2MStop20_dvd_dvu;
      double tmp_3 = dMStop21_dvd * dMStop21_dvu * Log(Sqr(MStop(1) / scale));
      double tmp_4 = -a0_mstop2 * d2MStop21_dvd_dvu;

      return (3.0 * oneOver16PiSqr * (tmp_1 + tmp_2 + tmp_3 + tmp_4));
   }

   double lowMSSM_tuning_calculator::deriv_d2DeltaV_dmq222_dvd() const
   {
      double scale = model.get_scale();

      double dMStop20_dvd = deriv_dMStop2_dparam(stop_mass::mstop_1, lowMSSM_info::vd);
      double dMStop20_dmq222 = deriv_dMStop2_dparam(stop_mass::mstop_1, lowMSSM_info::mq222);
      double d2MStop20_dmq222_dvd = deriv_d2MStop2_dparam_dparam(stop_mass::mstop_1, lowMSSM_info::mq222, lowMSSM_info::vd);
      double dMStop21_dvd = deriv_dMStop2_dparam(stop_mass::mstop_2, lowMSSM_info::vd);
      double dMStop21_dmq222 = deriv_dMStop2_dparam(stop_mass::mstop_2, lowMSSM_info::mq222);
      double d2MStop21_dmq222_dvd = deriv_d2MStop2_dparam_dparam(stop_mass::mstop_2, lowMSSM_info::mq222, lowMSSM_info::vd);

      double a0_mstop1 = passarino_veltman::ReA0(Sqr(MStop(0)), Sqr(scale));
      double a0_mstop2 = passarino_veltman::ReA0(Sqr(MStop(1)), Sqr(scale));

      double tmp_1 = dMStop20_dmq222 * dMStop20_dvd * Log(Sqr(MStop(0) / scale));
      double tmp_2 = -a0_mstop1 * d2MStop20_dmq222_dvd;
      double tmp_3 = dMStop21_dmq222 * dMStop21_dvd * Log(Sqr(MStop(1) / scale));
      double tmp_4 = -a0_mstop2 * d2MStop21_dmq222_dvd;

      return (3.0 * oneOver16PiSqr * (tmp_1 + tmp_2 + tmp_3 + tmp_4));
   }

   double lowMSSM_tuning_calculator::deriv_d2DeltaV_dmu222_dvd() const
   {
      double scale = model.get_scale();

      double dMStop20_dvd = deriv_dMStop2_dparam(stop_mass::mstop_1, lowMSSM_info::vd);
      double dMStop20_dmu222 = deriv_dMStop2_dparam(stop_mass::mstop_1, lowMSSM_info::mu222);
      double d2MStop20_dmu222_dvd = deriv_d2MStop2_dparam_dparam(stop_mass::mstop_1, lowMSSM_info::mu222, lowMSSM_info::vd);
      double dMStop21_dvd = deriv_dMStop2_dparam(stop_mass::mstop_2, lowMSSM_info::vd);
      double dMStop21_dmu222 = deriv_dMStop2_dparam(stop_mass::mstop_2, lowMSSM_info::mu222);
      double d2MStop21_dmu222_dvd = deriv_d2MStop2_dparam_dparam(stop_mass::mstop_2, lowMSSM_info::mu222, lowMSSM_info::vd);

      double a0_mstop1 = passarino_veltman::ReA0(Sqr(MStop(0)), Sqr(scale));
      double a0_mstop2 = passarino_veltman::ReA0(Sqr(MStop(1)), Sqr(scale));

      double tmp_1 = dMStop20_dmu222 * dMStop20_dvd * Log(Sqr(MStop(0) / scale));
      double tmp_2 = -a0_mstop1 * d2MStop20_dmu222_dvd;
      double tmp_3 = dMStop21_dmu222 * dMStop21_dvd * Log(Sqr(MStop(1) / scale));
      double tmp_4 = -a0_mstop2 * d2MStop21_dmu222_dvd;

      return (3.0 * oneOver16PiSqr * (tmp_1 + tmp_2 + tmp_3 + tmp_4));
   }

   double lowMSSM_tuning_calculator::deriv_d2DeltaV_dMu_dvd() const
   {
      double scale = model.get_scale();

      double dMStop20_dvd = deriv_dMStop2_dparam(stop_mass::mstop_1, lowMSSM_info::vd);
      double dMStop20_dMu = deriv_dMStop2_dparam(stop_mass::mstop_1, lowMSSM_info::Mu);
      double d2MStop20_dMu_dvd = deriv_d2MStop2_dparam_dparam(stop_mass::mstop_1, lowMSSM_info::Mu, lowMSSM_info::vd);
      double dMStop21_dvd = deriv_dMStop2_dparam(stop_mass::mstop_2, lowMSSM_info::vd);
      double dMStop21_dMu = deriv_dMStop2_dparam(stop_mass::mstop_2, lowMSSM_info::Mu);
      double d2MStop21_dMu_dvd = deriv_d2MStop2_dparam_dparam(stop_mass::mstop_2, lowMSSM_info::Mu, lowMSSM_info::vd);

      double a0_mstop1 = passarino_veltman::ReA0(Sqr(MStop(0)), Sqr(scale));
      double a0_mstop2 = passarino_veltman::ReA0(Sqr(MStop(1)), Sqr(scale));

      double tmp_1 = dMStop20_dMu * dMStop20_dvd * Log(Sqr(MStop(0) / scale));
      double tmp_2 = -a0_mstop1 * d2MStop20_dMu_dvd;
      double tmp_3 = dMStop21_dMu * dMStop21_dvd * Log(Sqr(MStop(1) / scale));
      double tmp_4 = -a0_mstop2 * d2MStop21_dMu_dvd;

      return (3.0 * oneOver16PiSqr * (tmp_1 + tmp_2 + tmp_3 + tmp_4));
   }

   double lowMSSM_tuning_calculator::deriv_d2DeltaV_dTYu22_dvd() const
   {
      double scale = model.get_scale();

      double dMStop20_dvd = deriv_dMStop2_dparam(stop_mass::mstop_1, lowMSSM_info::vd);
      double dMStop20_dTYu22 = deriv_dMStop2_dparam(stop_mass::mstop_1, lowMSSM_info::TYu22);
      double d2MStop20_dTYu22_dvd = deriv_d2MStop2_dparam_dparam(stop_mass::mstop_1, lowMSSM_info::TYu22, lowMSSM_info::vd);
      double dMStop21_dvd = deriv_dMStop2_dparam(stop_mass::mstop_2, lowMSSM_info::vd);
      double dMStop21_dTYu22 = deriv_dMStop2_dparam(stop_mass::mstop_2, lowMSSM_info::TYu22);
      double d2MStop21_dTYu22_dvd = deriv_d2MStop2_dparam_dparam(stop_mass::mstop_2, lowMSSM_info::TYu22, lowMSSM_info::vd);

      double a0_mstop1 = passarino_veltman::ReA0(Sqr(MStop(0)), Sqr(scale));
      double a0_mstop2 = passarino_veltman::ReA0(Sqr(MStop(1)), Sqr(scale));

      double tmp_1 = dMStop20_dTYu22 * dMStop20_dvd * Log(Sqr(MStop(0) / scale));
      double tmp_2 = -a0_mstop1 * d2MStop20_dTYu22_dvd;
      double tmp_3 = dMStop21_dTYu22 * dMStop21_dvd * Log(Sqr(MStop(1) / scale));
      double tmp_4 = -a0_mstop2 * d2MStop21_dTYu22_dvd;

      return (3.0 * oneOver16PiSqr * (tmp_1 + tmp_2 + tmp_3 + tmp_4));
   }

   double lowMSSM_tuning_calculator::deriv_d2DeltaV_dmq222_dvu() const
   {
      double scale = model.get_scale();

      double dMStop20_dvu = deriv_dMStop2_dparam(stop_mass::mstop_1, lowMSSM_info::vu);
      double dMStop20_dmq222 = deriv_dMStop2_dparam(stop_mass::mstop_1, lowMSSM_info::mq222);
      double d2MStop20_dmq222_dvu = deriv_d2MStop2_dparam_dparam(stop_mass::mstop_1, lowMSSM_info::mq222, lowMSSM_info::vu);
      double dMStop21_dvu = deriv_dMStop2_dparam(stop_mass::mstop_2, lowMSSM_info::vu);
      double dMStop21_dmq222 = deriv_dMStop2_dparam(stop_mass::mstop_2, lowMSSM_info::mq222);
      double d2MStop21_dmq222_dvu = deriv_d2MStop2_dparam_dparam(stop_mass::mstop_2, lowMSSM_info::mq222, lowMSSM_info::vu);

      double a0_mstop1 = passarino_veltman::ReA0(Sqr(MStop(0)), Sqr(scale));
      double a0_mstop2 = passarino_veltman::ReA0(Sqr(MStop(1)), Sqr(scale));

      double tmp_1 = dMStop20_dmq222 * dMStop20_dvu * Log(Sqr(MStop(0) / scale));
      double tmp_2 = -a0_mstop1 * d2MStop20_dmq222_dvu;
      double tmp_3 = dMStop21_dmq222 * dMStop21_dvu * Log(Sqr(MStop(1) / scale));
      double tmp_4 = -a0_mstop2 * d2MStop21_dmq222_dvu;

      return (3.0 * oneOver16PiSqr * (tmp_1 + tmp_2 + tmp_3 + tmp_4));
   }

   double lowMSSM_tuning_calculator::deriv_d2DeltaV_dmu222_dvu() const
   {
      double scale = model.get_scale();

      double dMStop20_dvu = deriv_dMStop2_dparam(stop_mass::mstop_1, lowMSSM_info::vu);
      double dMStop20_dmu222 = deriv_dMStop2_dparam(stop_mass::mstop_1, lowMSSM_info::mu222);
      double d2MStop20_dmu222_dvu = deriv_d2MStop2_dparam_dparam(stop_mass::mstop_1, lowMSSM_info::mu222, lowMSSM_info::vu);
      double dMStop21_dvu = deriv_dMStop2_dparam(stop_mass::mstop_2, lowMSSM_info::vu);
      double dMStop21_dmu222 = deriv_dMStop2_dparam(stop_mass::mstop_2, lowMSSM_info::mu222);
      double d2MStop21_dmu222_dvu = deriv_d2MStop2_dparam_dparam(stop_mass::mstop_2, lowMSSM_info::mu222, lowMSSM_info::vu);

      double a0_mstop1 = passarino_veltman::ReA0(Sqr(MStop(0)), Sqr(scale));
      double a0_mstop2 = passarino_veltman::ReA0(Sqr(MStop(1)), Sqr(scale));

      double tmp_1 = dMStop20_dmu222 * dMStop20_dvu * Log(Sqr(MStop(0) / scale));
      double tmp_2 = -a0_mstop1 * d2MStop20_dmu222_dvu;
      double tmp_3 = dMStop21_dmu222 * dMStop21_dvu * Log(Sqr(MStop(1) / scale));
      double tmp_4 = -a0_mstop2 * d2MStop21_dmu222_dvu;

      return (3.0 * oneOver16PiSqr * (tmp_1 + tmp_2 + tmp_3 + tmp_4));
   }

   double lowMSSM_tuning_calculator::deriv_d2DeltaV_dMu_dvu() const
   {
      double scale = model.get_scale();

      double dMStop20_dvu = deriv_dMStop2_dparam(stop_mass::mstop_1, lowMSSM_info::vu);
      double dMStop20_dMu = deriv_dMStop2_dparam(stop_mass::mstop_1, lowMSSM_info::Mu);
      double d2MStop20_dMu_dvu = deriv_d2MStop2_dparam_dparam(stop_mass::mstop_1, lowMSSM_info::Mu, lowMSSM_info::vu);
      double dMStop21_dvu = deriv_dMStop2_dparam(stop_mass::mstop_2, lowMSSM_info::vu);
      double dMStop21_dMu = deriv_dMStop2_dparam(stop_mass::mstop_2, lowMSSM_info::Mu);
      double d2MStop21_dMu_dvu = deriv_d2MStop2_dparam_dparam(stop_mass::mstop_2, lowMSSM_info::Mu, lowMSSM_info::vu);

      double a0_mstop1 = passarino_veltman::ReA0(Sqr(MStop(0)), Sqr(scale));
      double a0_mstop2 = passarino_veltman::ReA0(Sqr(MStop(1)), Sqr(scale));

      double tmp_1 = dMStop20_dMu * dMStop20_dvu * Log(Sqr(MStop(0) / scale));
      double tmp_2 = -a0_mstop1 * d2MStop20_dMu_dvu;
      double tmp_3 = dMStop21_dMu * dMStop21_dvu * Log(Sqr(MStop(1) / scale));
      double tmp_4 = -a0_mstop2 * d2MStop21_dMu_dvu;

      return (3.0 * oneOver16PiSqr * (tmp_1 + tmp_2 + tmp_3 + tmp_4));
   }

   double lowMSSM_tuning_calculator::deriv_d2DeltaV_dTYu22_dvu() const
   {
      double scale = model.get_scale();

      double dMStop20_dvu = deriv_dMStop2_dparam(stop_mass::mstop_1, lowMSSM_info::vu);
      double dMStop20_dTYu22 = deriv_dMStop2_dparam(stop_mass::mstop_1, lowMSSM_info::TYu22);
      double d2MStop20_dTYu22_dvu = deriv_d2MStop2_dparam_dparam(stop_mass::mstop_1, lowMSSM_info::TYu22, lowMSSM_info::vu);
      double dMStop21_dvu = deriv_dMStop2_dparam(stop_mass::mstop_2, lowMSSM_info::vu);
      double dMStop21_dTYu22 = deriv_dMStop2_dparam(stop_mass::mstop_2, lowMSSM_info::TYu22);
      double d2MStop21_dTYu22_dvu = deriv_d2MStop2_dparam_dparam(stop_mass::mstop_2, lowMSSM_info::TYu22, lowMSSM_info::vu);

      double a0_mstop1 = passarino_veltman::ReA0(Sqr(MStop(0)), Sqr(scale));
      double a0_mstop2 = passarino_veltman::ReA0(Sqr(MStop(1)), Sqr(scale));

      double tmp_1 = dMStop20_dTYu22 * dMStop20_dvu * Log(Sqr(MStop(0) / scale));
      double tmp_2 = -a0_mstop1 * d2MStop20_dTYu22_dvu;
      double tmp_3 = dMStop21_dTYu22 * dMStop21_dvu * Log(Sqr(MStop(1) / scale));
      double tmp_4 = -a0_mstop2 * d2MStop21_dTYu22_dvu;

      return (3.0 * oneOver16PiSqr * (tmp_1 + tmp_2 + tmp_3 + tmp_4));
   }

} // namespace flexiblesusy