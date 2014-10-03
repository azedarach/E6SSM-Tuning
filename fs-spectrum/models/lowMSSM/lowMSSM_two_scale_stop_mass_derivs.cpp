#include "lowMSSM_ewsb_conditions.hpp"

namespace flexiblesusy {

   double lowMSSM_ewsb_conditions::deriv_dMStop2_dvd(stop_mass which_stop) const
   {
      double mu = model.get_Mu();
      double vd = model.get_vd();
      double g1 = model.get_g1();
      double g2 = model.get_g2();
      double yt = model.get_Yu(2, 2);
      double rqq_val = stop_RQQ();
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

   double lowMSSM_ewsb_conditions::deriv_dMStop2_dvu(stop_mass which_stop) const
   {
      double at = model.get_TYu(2, 2);
      double vu = model.get_vu();
      double g1 = model.get_g1();
      double g2 = model.get_g2();
      double yt = model.get_Yu(2, 2);
      double rqq_val = stop_RQQ();
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
//TODO
   double lowMSSM_ewsb_conditions::deriv_dMStop2_dg1(stop_mass which_stop) const;
//TODO
   double lowMSSM_ewsb_conditions::deriv_dMStop2_dg2(stop_mass which_stop) const;
//TODO
   double lowMSSM_ewsb_conditions::deriv_dMStop2_dYu22(stop_mass which_stop) const;

   double lowMSSM_ewsb_conditions::deriv_dMStop2_dmq222(stop_mass which_stop) const
   {
      double mqqsq = stop_MQQ2();
      double rt = stop_discriminant();
      if (which_stop == stop_mass::mstop_1) {
         return 0.5 * (1.0 - mqqsq / Sqrt(rt));
      } else {
         return 0.5 * (1.0 + mqqsq / Sqrt(rt));
      }
   }

   double lowMSSM_ewsb_conditions::deriv_dMStop2_dmu222(stop_mass which_stop) const
   {
      double mqqsq = stop_MQQ2();
      double rt = stop_discriminant();
      if (which_stop == stop_mass::mstop_1) {
         return 0.5 * (1.0 + mqqsq / Sqrt(rt));
      } else {
         return 0.5 * (1.0 - mqqsq / Sqrt(rt));
      }
   }

   double lowMSSM_ewsb_conditions::deriv_dMStop2_dMu(stop_mass which_stop) const
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

   double lowMSSM_ewsb_conditions::deriv_dMStop2_dTYu22(stop_mass which_stop) const
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

   double lowMSSM_ewsb_conditions::deriv_dMStop2_dparam(stop_mass which_stop, lowMSSM_info::Parameters p) const
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
      case lowMSSM_info::g1 {
         return deriv_dMStop2_dg1(which_stop);
      }
      case lowMSSM_info::g2 {
         return deriv_dMStop2_dg2(which_stop);
      }
      case lowMSSM_info::Yu22 {
         return deriv_dMStop2_dYu22(which_stop);
      }
      default: {
         return 0.;
      }
      }
   }

   double lowMSSM_ewsb_conditions::deriv_d2MStop2_dvd_dvd(stop_mass which_stop) const
   {
      double mu = model.get_Mu();
      double vd = model.get_vd();
      double yt = model.get_Yu(2, 2);
      double g1 = model.get_g1();
      double g2 = model.get_g2();

      double gbar_val = gbar();
      double stop_mixing = 1.4142135623730951 * stop_mass_matrix_LR_entry();
      double RQQ_val = stop_RQQ();
      double rt = stop_discriminant();

      double split = 0.03125 * ((4.0 * RQQ_val + Sqr((g2 * g2 - g1 * g1) * vd) + 32.0 * Sqr(yt * mu)) / Sqrt(rt)
                                - Sqr(vd * RQQ_val - 8.0 * yt * mu * stop_mixing) / (rt * Sqrt(rt)));

      if (which_stop == stop_mass::mstop_1) {
         return 0.125 * Sqr(gbar_val) - split;
      } else {
         return 0.125 * Sqr(gbar_val) + split;
      }
   }

   double lowMSSM_ewsb_conditions::deriv_d2MStop2_dvu_dvu(stop_mass which_stop) const
   {
      double vu = model.get_vu();
      double yt = model.get_Yu(2, 2);
      double at = model.get_TYu(2, 2);
      double g1 = model.get_g1();
      double g2 = model.get_g2();

      double gbar_val = gbar();
      double stop_mixing = 1.4142135623730951 * stop_mass_matrix_LR_entry();
      double RQQ_val = stop_RQQ();
      double rt = stop_discriminant();

      double split = (0.03125 * Sqr(vu * (g2 * g2 - g1 * g1)) - 0.125 * RQQ_val + Sqr(at)
                      - 0.03125 * Sqr(vu) * Sqr(8.0 * at * stop_mixing / vu - RQQ_val) / rt) / Sqrt(rt);

      if (which_stop == stop_mass::mstop_1) {
         return Sqr(yt) - 0.125 * Sqr(gbar_val) - split;
      } else {
         return Sqr(yt) - 0.125 * Sqr(gbar_val) + split;
      }
   }

   double lowMSSM_ewsb_conditions::deriv_d2MStop2_dvd_dvu(stop_mass which_stop) const
   {
      double mu = model.get_Mu();
      double vd = model.get_vd();
      double vu = model.get_vu();
      double yt = model.get_Yu(2, 2);
      double at = model.get_TYu(2, 2);
      double g1 = model.get_g1();
      double g2 = model.get_g2();

      double stop_mixing = 1.4142135623730951 * stop_mass_matrix_LR_entry();
      double RQQ_val = stop_RQQ();
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
//TODO
   double lowMSSM_ewsb_conditions::deriv_d2MStop2_dg1_dvd(stop_mass which_stop) const;
//TODO
   double lowMSSM_ewsb_conditions::deriv_d2MStop2_dg2_dvd(stop_mass which_stop) const;
//TODO
   double lowMSSM_ewsb_conditions::deriv_d2MStop2_dYu22_dvd(stop_mass which_stop) const;

   double lowMSSM_ewsb_conditions::deriv_d2MStop2_dmq222_dvd(stop_mass which_stop) const
   {
      double mu = model.get_Mu();
      double vd = model.get_vd();
      double yt = model.get_Yu(2, 2);
      double g1 = model.get_g1();
      double g2 = model.get_g2();

      double stop_mixing = 1.4142135623730951 * stop_mass_matrix_LR_entry();
      double mqqsq = stop_MQQ2();
      double RQQ_val = stop_RQQ();
      double rt = stop_discriminant();

      double deriv = mqqsq * (0.125 * RQQ_val - yt * mu * stop_mixing / vd) / rt
         - 0.125 * ( g2 * g2 - g1 * g1);

      if (which_stop == stop_mass::mstop_1) {
         return vd * deriv / Sqrt(rt);
      } else {
         return -vd * deriv / Sqrt(rt);
      }
   }

   double lowMSSM_ewsb_conditions::deriv_d2MStop2_dmu222_dvd(stop_mass which_stop) const
   {
      double mu = model.get_Mu();
      double vd = model.get_vd();
      double yt = model.get_Yu(2, 2);
      double g1 = model.get_g1();
      double g2 = model.get_g2();

      double stop_mixing = 1.4142135623730951 * stop_mass_matrix_LR_entry();
      double mqqsq = stop_MQQ2();
      double RQQ_val = stop_RQQ();
      double rt = stop_discriminant();

      double deriv = mqqsq * (0.125 * RQQ_val - yt * mu * stop_mixing / vd) / rt
         - 0.125 * ( g2 * g2 - g1 * g1);

      if (which_stop == stop_mass::mstop_1) {
         return -vd * deriv / Sqrt(rt);
      } else {
         return vd * deriv / Sqrt(rt);
      }
   }

   double lowMSSM_ewsb_conditions::deriv_d2MStop2_dMu_dvd(stop_mass which_stop) const
   {
      double mu = model.get_Mu();
      double vd = model.get_vd();
      double yt = model.get_Yu(2, 2);

      double stop_mixing = 1.4142135623730951 * stop_mass_matrix_LR_entry();
      double RQQ_val = stop_RQQ();
      double rt = stop_discriminant();

      double deriv = 2.0 * yt * vd * stop_mixing * (0.125 * RQQ_val - yt * mu * stop_mixing / vd) / rt
         - yt * stop_mixing / vd + yt * yt * mu;

      if (which_stop == stop_mass::mstop_1) {
         return -vd * deriv / Sqrt(rt);
      } else {
         return vd * deriv / Sqrt(rt);
      }
   }

   double lowMSSM_ewsb_conditions::deriv_d2MStop2_dTYu22_dvd(stop_mass which_stop) const
   {
      double mu = model.get_Mu();
      double vd = model.get_vd();
      double vu = model.get_vu();
      double yt = model.get_Yu(2, 2);

      double stop_mixing = 1.4142135623730951 * stop_mass_matrix_LR_entry();
      double RQQ_val = stop_RQQ();
      double rt = stop_discriminant();

      double deriv = 2.0 * vu * stop_mixing * (0.125 * RQQ_val - yt * mu * stop_mixing / vd) / rt
         + yt * vu * mu / vd;

      if (which_stop == stop_mass::mstop_1) {
         return vd * deriv / Sqrt(rt);
      } else {
         return -vd * deriv / Sqrt(rt);
      }
   }
//TODO
   double lowMSSM_ewsb_conditions::deriv_d2MStop2_dg1_dvu(stop_mass which_stop) const;
//TODO
   double lowMSSM_ewsb_conditions::deriv_d2MStop2_dg2_dvu(stop_mass which_stop) const;
//TODO
   double lowMSSM_ewsb_conditions::deriv_d2MStop2_dYu22_dvu(stop_mass which_stop) const;

   double lowMSSM_ewsb_conditions::deriv_d2MStop2_dmq222_dvu(stop_mass which_stop) const
   {
      double at = model.get_TYu(2, 2);
      double vu = model.get_vu();
      double g1 = model.get_g1();
      double g2 = model.get_g2();

      double stop_mixing = 1.4142135623730951 * stop_mass_matrix_LR_entry();
      double mqqsq = stop_MQQ2();
      double RQQ_val = stop_RQQ();
      double rt = stop_discriminant();

      double deriv = mqqsq * (at * stop_mixing / vu - 0.125 * RQQ_val) / rt
         + 0.125 * ( g2 * g2 - g1 * g1);

      if (which_stop == stop_mass::mstop_1) {
         return vu * deriv / Sqrt(rt);
      } else {
         return -vu * deriv / Sqrt(rt);
      }
   }

   double lowMSSM_ewsb_conditions::deriv_d2MStop2_dmu222_dvu(stop_mass which_stop) const
   {
      double at = model.get_TYu(2, 2);
      double vu = model.get_vu();
      double g1 = model.get_g1();
      double g2 = model.get_g2();

      double stop_mixing = 1.4142135623730951 * stop_mass_matrix_LR_entry();
      double mqqsq = stop_MQQ2();
      double RQQ_val = stop_RQQ();
      double rt = stop_discriminant();

      double deriv = mqqsq * (at * stop_mixing / vu - 0.125 * RQQ_val) / rt
         + 0.125 * ( g2 * g2 - g1 * g1);

      if (which_stop == stop_mass::mstop_1) {
         return -vu * deriv / Sqrt(rt);
      } else {
         return vu * deriv / Sqrt(rt);
      }
   }

   double lowMSSM_ewsb_conditions::deriv_d2MStop2_dMu_dvu(stop_mass which_stop) const
   {
      double at = model.get_TYu(2, 2);
      double vd = model.get_vd();
      double vu = model.get_vu();
      double yt = model.get_Yu(2, 2);

      double stop_mixing = 1.4142135623730951 * stop_mass_matrix_LR_entry();
      double RQQ_val = stop_RQQ();
      double rt = stop_discriminant();

      double deriv = 2.0 * yt * vd * stop_mixing * (at * stop_mixing / vu - 0.125 * RQQ_val) / rt
         - yt * at * vd / vu;

      if (which_stop == stop_mass::mstop_1) {
         return -vu * deriv / Sqrt(rt);
      } else {
         return vu * deriv / Sqrt(rt);
      }
   }

   double lowMSSM_ewsb_conditions::deriv_d2MStop2_dTYu22_dvu(stop_mass which_stop) const
   {
      double at = model.get_TYu(2, 2);
      double vu = model.get_vu();

      double stop_mixing = 1.4142135623730951 * stop_mass_matrix_LR_entry();
      double RQQ_val = stop_RQQ();
      double rt = stop_discriminant();

      double deriv = 2.0 * vu * stop_mixing * (at * stop_mixing / vu - 0.125 * RQQ_val) / rt
         - stop_mixing / vu - at;

      if (which_stop == stop_mass::mstop_1) {
         return vu * deriv / Sqrt(rt);
      } else {
         return -vu * deriv / Sqrt(rt);
      }
   }

   double lowMSSM_ewsb_conditions::deriv_d2MStop2_dparam_dparam(stop_mass which_stop, lowMSSM_info::Parameters p1, 
                                                                lowMSSM_info::Parameters p2) const
   {
      if (equal_as_unordered_pairs(p1, p2, lowMSSM_info::vd, lowMSSM_info::vd)) {
         return deriv_d2MStop2_dvd_dvd(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowMSSM_info::vu, lowMSSM_info::vu)) {
         return deriv_d2MStop2_dvu_dvu(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowMSSM_info::vd, lowMSSM_info::vu)) {
         return deriv_d2MStop2_dvd_dvu(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowMSSM_info::g1, lowMSSM_info::vd)) {
         return deriv_d2MStop2_dg1_dvd(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowMSSM_info::g2, lowMSSM_info::vd)) {
         return deriv_d2MStop2_dg2_dvd(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowMSSM_info::Yu22, lowMSSM_info::vd)) {
         return deriv_d2MStop2_dYu22_dvd(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowMSSM_info::mq222, lowMSSM_info::vd)) {
         return deriv_d2MStop2_dmq222_dvd(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowMSSM_info::mu222, lowMSSM_info::vd)) {
         return deriv_d2MStop2_dmu222_dvd(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowMSSM_info::Mu, lowMSSM_info::vd)) {
         return deriv_d2MStop2_dMu_dvd(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowMSSM_info::TYu22, lowMSSM_info::vd)) {
         return deriv_d2MStop2_dTYu22_dvd(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowMSSM_info::g1, lowMSSM_info::vu)) {
         return deriv_d2MStop2_dg1_dvu(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowMSSM_info::g2, lowMSSM_info::vu)) {
         return deriv_d2MStop2_dg2_dvu(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowMSSM_info::Yu22, lowMSSM_info::vu)) {
         return deriv_d2MStop2_dYu22_dvu(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowMSSM_info::mq222, lowMSSM_info::vu)) {
         return deriv_d2MStop2_dmq222_dvu(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowMSSM_info::mu222, lowMSSM_info::vu)) {
         return deriv_d2MStop2_dmu222_dvu(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowMSSM_info::Mu, lowMSSM_info::vu)) {
         return deriv_d2MStop2_dMu_dvu(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowMSSM_info::TYu22, lowMSSM_info::vu)) {
         return deriv_d2MStop2_dTYu22_dvu(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowMSSM_info::g1, lowMSSM_info::g1)) {
         return deriv_d2MStop2_dg1_dg1(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowMSSM_info::g2, lowMSSM_info::g2)) {
         return deriv_d2MStop2_dg2_dg2(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowMSSM_info::Yu22, lowMSSM_info::Yu22)) {
         return deriv_d2MStop2_dYu22_dYu22(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowMSSM_info::mq222, lowMSSM_info::mq222)) {
         return deriv_d2MStop2_dmq222_dmq222(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowMSSM_info::mu222, lowMSSM_info::mu222)) {
         return deriv_d2MStop2_dmu222_dmu222(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowMSSM_info::Mu, lowMSSM_info::Mu)) {
         return deriv_d2MStop2_dMu_dMu(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowMSSM_info::TYu22, lowMSSM_info::TYu22)) {
         return deriv_d2MStop2_dTYu22_dTYu22(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowMSSM_info::g2, lowMSSM_info::g1)) {
         return deriv_d2MStop2_dg1_dg2(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowMSSM_info::Yu22, lowMSSM_info::g1)) {
         return deriv_d2MStop2_dYu22_dg1(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowMSSM_info::mq222, lowMSSM_info::g1)) {
         return deriv_d2MStop2_dmq222_dg1(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowMSSM_info::mu222, lowMSSM_info::g1)) {
         return deriv_d2MStop2_dmu222_dg1(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowMSSM_info::Mu, lowMSSM_info::g1)) {
         return deriv_d2MStop2_dMu_dg1(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowMSSM_info::TYu22, lowMSSM_info::g1)) {
         return deriv_d2MStop2_dTYu22_dg1(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowMSSM_info::Yu22, lowMSSM_info::g2)) {
         return deriv_d2MStop2_dYu22_dg2(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowMSSM_info::mq222, lowMSSM_info::g2)) {
         return deriv_d2MStop2_dmq222_dg2(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowMSSM_info::mu222, lowMSSM_info::g2)) {
         return deriv_d2MStop2_dmu222_dg2(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowMSSM_info::Mu, lowMSSM_info::g2)) {
         return deriv_d2MStop2_dMu_dg2(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowMSSM_info::TYu22, lowMSSM_info::g2)) {
         return deriv_d2MStop2_dTYu22_dg2(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowMSSM_info::mq222, lowMSSM_info::Yu22)) {
         return deriv_d2MStop2_dmq222_dYu22(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowMSSM_info::mu222, lowMSSM_info::Yu22)) {
         return deriv_d2MStop2_dmu222_dYu22(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowMSSM_info::Mu, lowMSSM_info::Yu22)) {
         return deriv_d2MStop2_dMu_dYu22(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowMSSM_info::TYu22, lowMSSM_info::Yu22)) {
         return deriv_d2MStop2_dTYu22_dYu22(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowMSSM_info::mu222, lowMSSM_info::mq222)) {
         return deriv_d2MStop2_dmu222_dmq222(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowMSSM_info::Mu, lowMSSM_info::mq222)) {
         return deriv_d2MStop2_dMu_dmq222(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowMSSM_info::TYu22, lowMSSM_info::mq222)) {
         return deriv_d2MStop2_dTYu22_dmq222(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowMSSM_info::Mu, lowMSSM_info::mu222)) {
         return deriv_d2MStop2_dMu_dmu222(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowMSSM_info::TYu22, lowMSSM_info::mu222)) {
         return deriv_d2MStop2_dTYu22_dmu222(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowMSSM_info::TYu22, lowMSSM_info::Mu)) {
         return deriv_d2MStop2_dTYu22_dMu(which_stop);
      } else {
         return 0.;
      }
   }
 

} // namespace flexiblesusy
