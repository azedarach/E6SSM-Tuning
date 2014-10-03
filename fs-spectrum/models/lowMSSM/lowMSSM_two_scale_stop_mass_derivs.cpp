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

   double lowMSSM_ewsb_conditions::deriv_d2MStop2_dvd_dvd(stop_mass which_stop) const;
   double lowMSSM_ewsb_conditions::deriv_d2MStop2_dvu_dvu(stop_mass which_stop) const;
   double lowMSSM_ewsb_conditions::deriv_d2MStop2_dvd_dvu(stop_mass which_stop) const;
   double lowMSSM_ewsb_conditions::deriv_d2MStop2_dg1_dvd(stop_mass which_stop) const;
   double lowMSSM_ewsb_conditions::deriv_d2MStop2_dg2_dvd(stop_mass which_stop) const;
   double lowMSSM_ewsb_conditions::deriv_d2MStop2_dYu22_dvd(stop_mass which_stop) const;
   double lowMSSM_ewsb_conditions::deriv_d2MStop2_dmq222_dvd(stop_mass which_stop) const;
   double lowMSSM_ewsb_conditions::deriv_d2MStop2_dmu222_dvd(stop_mass which_stop) const;
   double lowMSSM_ewsb_conditions::deriv_d2MStop2_dMu_dvd(stop_mass which_stop) const;
   double lowMSSM_ewsb_conditions::deriv_d2MStop2_dTYu22_dvd(stop_mass which_stop) const;
   double lowMSSM_ewsb_conditions::deriv_d2MStop2_dg1_dvu(stop_mass which_stop) const;
   double lowMSSM_ewsb_conditions::deriv_d2MStop2_dg2_dvu(stop_mass which_stop) const;
   double lowMSSM_ewsb_conditions::deriv_d2MStop2_dYu22_dvu(stop_mass which_stop) const;
   double lowMSSM_ewsb_conditions::deriv_d2MStop2_dmq222_dvu(stop_mass which_stop) const;
   double lowMSSM_ewsb_conditions::deriv_d2MStop2_dmu222_dvu(stop_mass which_stop) const;
   double lowMSSM_ewsb_conditions::deriv_d2MStop2_dMu_dvu(stop_mass which_stop) const;
   double lowMSSM_ewsb_conditions::deriv_d2MStop2_dTYu22_dvu(stop_mass which_stop) const;
   double lowMSSM_ewsb_conditions::deriv_d2MStop2_dparam_dparam(stop_mass which_stop, lowMSSM_info::Parameters p1, 
                                                                lowMSSM_info::Parameters p2) const
   {
      if (p1 == lowMSSM_info::vd && p2 == lowMSSM_info::vd) {
         return deriv_d2MStop2_dvd_dvd(which_stop);
      } else if (p1 == lowMSSM_info::vu && p2 == lowMSSM_info::vu) {
         return deriv_d2MStop2_dvu_dvu(which_stop);
      } else if ((p1 == lowMSSM_info::vu && p2 == lowMSSM_info::vd) ||
                 (p1 == lowMSSM_info::vd && p2 == lowMSSM_info::vu)) {
         return deriv_d2MStop2_dvd_dvu(which_stop);
      } else if ((p1 == lowMSSM_info::g1 && p2 == lowMSSM_info::vd) || 
                 (p1 == lowMSSM_info::vd && p2 == lowMSSM_info::g1)) {
         return deriv_d2MStop2_dg1_dvd(which_stop);
      } else if ((p1 == lowMSSM_info::g2 && p2 == lowMSSM_info::vd) || 
                 (p1 == lowMSSM_info::vd && p2 == lowMSSM_info::g2)) {
         return deriv_d2MStop2_dg2_dvd(which_stop);
      } else if ((p1 == lowMSSM_info::Yu22 && p2 == lowMSSM_info::vd) ||
                 (p1 == lowMSSM_info::vd && p2 == lowMSSM_info::Yu22)) {
         return deriv_d2MStop2_dYu22_dvd(which_stop);
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
      } else if ((p1 == lowMSSM_info::g1 && p2 == lowMSSM_info::vu) || 
                 (p1 == lowMSSM_info::vu && p2 == lowMSSM_info::g1)) {
         return deriv_d2MStop2_dg1_dvu(which_stop);
      } else if ((p1 == lowMSSM_info::g2 && p2 == lowMSSM_info::vu) || 
                 (p1 == lowMSSM_info::vu && p2 == lowMSSM_info::g2)) {
         return deriv_d2MStop2_dg2_dvu(which_stop);
      } else if ((p1 == lowMSSM_info::Yu22 && p2 == lowMSSM_info::vu) ||
                 (p1 == lowMSSM_info::vu && p2 == lowMSSM_info::Yu22)) {
         return deriv_d2MStop2_dYu22_dvu(which_stop);
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
 

} // namespace flexiblesusy
