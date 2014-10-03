#include "lowMSSM_ewsb_conditions.hpp"

namespace flexiblesusy {

   //TODO
   double lowMSSM_ewsb_conditions::deriv_dMSbottom2_dvd(sbottom_mass which_sbottom) const;
//TODO
   double lowMSSM_ewsb_conditions::deriv_dMSbottom2_dvu(sbottom_mass which_sbottom) const;
//TODO
   double lowMSSM_ewsb_conditions::deriv_dMSbottom2_dg1(sbottom_mass which_sbottom) const;
//TODO
   double lowMSSM_ewsb_conditions::deriv_dMSbottom2_dg2(sbottom_mass which_sbottom) const;
//TODO
   double lowMSSM_ewsb_conditions::deriv_dMSbottom2_dYd22(sbottom_mass which_sbottom) const;
//TODO
   double lowMSSM_ewsb_conditions::deriv_dMSbottom2_dmq222(sbottom_mass which_sbottom) const;
//TODO
   double lowMSSM_ewsb_conditions::deriv_dMSbottom2_dmd222(sbottom_mass which_sbottom) const;
//TODO
   double lowMSSM_ewsb_conditions::deriv_dMSbottom2_dMu(sbottom_mass which_sbottom) const;
//TODO
   double lowMSSM_ewsb_conditions::deriv_dMSbottom2_dTYd22(sbottom_mass which_sbottom) const;

   double lowMSSM_ewsb_conditions::deriv_dMSbottom2_dparam(sbottom_mass which_sbottom, lowMSSM_info::Parameters p) const
   {
      switch(p) {
      case lowMSSM_info::mq222: {
         return deriv_dMSbottom2_dmq222(which_sbottom);
      }
      case lowMSSM_info::md222: {
         return deriv_dMSbottom2_dmd222(which_sbottom);
      }
      case lowMSSM_info::Mu: {
         return deriv_dMSbottom2_dMu(which_sbottom);
      }
      case lowMSSM_info::vd: {
         return deriv_dMSbottom2_dvd(which_sbottom);
      }
      case lowMSSM_info::vu: {
         return deriv_dMSbottom2_dvu(which_sbottom);
      }
      case lowMSSM_info::TYd22: {
         return deriv_dMSbottom2_dTYd22(which_sbottom);
      }
      case lowMSSM_info::g1 {
         return deriv_dMSbottom2_dg1(which_sbottom);
      }
      case lowMSSM_info::g2 {
         return deriv_dMSbottom2_dg2(which_sbottom);
      }
      case lowMSSM_info::Yd22 {
         return deriv_dMSbottom2_dYd22(which_sbottom);
      }
      default: {
         return 0.;
      }
      }
   }

   double lowMSSM_ewsb_conditions::deriv_d2MSbottom2_dvd_dvd(sbottom_mass which_sbottom) const;
   double lowMSSM_ewsb_conditions::deriv_d2MSbottom2_dvu_dvu(sbottom_mass which_sbottom) const;
   double lowMSSM_ewsb_conditions::deriv_d2MSbottom2_dvd_dvu(sbottom_mass which_sbottom) const;
   double lowMSSM_ewsb_conditions::deriv_d2MSbottom2_dg1_dvd(sbottom_mass which_sbottom) const;
   double lowMSSM_ewsb_conditions::deriv_d2MSbottom2_dg2_dvd(sbottom_mass which_sbottom) const;
   double lowMSSM_ewsb_conditions::deriv_d2MSbottom2_dYd22_dvd(sbottom_mass which_sbottom) const;
   double lowMSSM_ewsb_conditions::deriv_d2MSbottom2_dmq222_dvd(sbottom_mass which_sbottom) const;
   double lowMSSM_ewsb_conditions::deriv_d2MSbottom2_dmd222_dvd(sbottom_mass which_sbottom) const;
   double lowMSSM_ewsb_conditions::deriv_d2MSbottom2_dMu_dvd(sbottom_mass which_sbottom) const;
   double lowMSSM_ewsb_conditions::deriv_d2MSbottom2_dTYd22_dvd(sbottom_mass which_sbottom) const;
   double lowMSSM_ewsb_conditions::deriv_d2MSbottom2_dg1_dvu(sbottom_mass which_sbottom) const;
   double lowMSSM_ewsb_conditions::deriv_d2MSbottom2_dg2_dvu(sbottom_mass which_sbottom) const;
   double lowMSSM_ewsb_conditions::deriv_d2MSbottom2_dYd22_dvu(sbottom_mass which_sbottom) const;
   double lowMSSM_ewsb_conditions::deriv_d2MSbottom2_dmq222_dvu(sbottom_mass which_sbottom) const;
   double lowMSSM_ewsb_conditions::deriv_d2MSbottom2_dmd222_dvu(sbottom_mass which_sbottom) const;
   double lowMSSM_ewsb_conditions::deriv_d2MSbottom2_dMu_dvu(sbottom_mass which_sbottom) const;
   double lowMSSM_ewsb_conditions::deriv_d2MSbottom2_dTYd22_dvu(sbottom_mass which_sbottom) const;
   double lowMSSM_ewsb_conditions::deriv_d2MSbottom2_dparam_dparam(sbottom_mass which_sbottom, lowMSSM_info::Parameters p1, 
                                          lowMSSM_info::Parameters p2) const
   {
      if (p1 == lowMSSM_info::vd && p2 == lowMSSM_info::vd) {
         return deriv_d2MSbottom2_dvd_dvd(which_sbottom);
      } else if (p1 == lowMSSM_info::vu && p2 == lowMSSM_info::vu) {
         return deriv_d2MSbottom2_dvu_dvu(which_sbottom);
      } else if ((p1 == lowMSSM_info::vu && p2 == lowMSSM_info::vd) ||
                 (p1 == lowMSSM_info::vd && p2 == lowMSSM_info::vu)) {
         return deriv_d2MSbottom2_dvd_dvu(which_sbottom);
      } else if ((p1 == lowMSSM_info::g1 && p2 == lowMSSM_info::vd) || 
                 (p1 == lowMSSM_info::vd && p2 == lowMSSM_info::g1)) {
         return deriv_d2MSbottom2_dg1_dvd(which_sbottom);
      } else if ((p1 == lowMSSM_info::g2 && p2 == lowMSSM_info::vd) || 
                 (p1 == lowMSSM_info::vd && p2 == lowMSSM_info::g2)) {
         return deriv_d2MSbottom2_dg2_dvd(which_sbottom);
      } else if ((p1 == lowMSSM_info::Yd22 && p2 == lowMSSM_info::vd) ||
                 (p1 == lowMSSM_info::vd && p2 == lowMSSM_info::Yd22)) {
         return deriv_d2MSbottom2_dYd22_dvd(which_sbottom);
      } else if ((p1 == lowMSSM_info::mq222 && p2 == lowMSSM_info::vd) ||
                 (p1 == lowMSSM_info::vd && p2 == lowMSSM_info::mq222)) {
         return deriv_d2MSbottom2_dmq222_dvd(which_sbottom);
      } else if ((p1 == lowMSSM_info::md222 && p2 == lowMSSM_info::vd) ||
                 (p1 == lowMSSM_info::vd && p2 == lowMSSM_info::md222)) {
         return deriv_d2MSbottom2_dmd222_dvd(which_sbottom);
      } else if ((p1 == lowMSSM_info::Mu && p2 == lowMSSM_info::vd) ||
                 (p1 == lowMSSM_info::vd && p2 == lowMSSM_info::Mu)) {
         return deriv_d2MSbottom2_dMu_dvd(which_sbottom);
      } else if ((p1 == lowMSSM_info::TYd22 && p2 == lowMSSM_info::vd) ||
                 (p1 == lowMSSM_info::vd && p2 == lowMSSM_info::TYd22)) {
         return deriv_d2MSbottom2_dTYd22_dvd(which_sbottom);
      } else if ((p1 == lowMSSM_info::g1 && p2 == lowMSSM_info::vu) || 
                 (p1 == lowMSSM_info::vu && p2 == lowMSSM_info::g1)) {
         return deriv_d2MSbottom2_dg1_dvu(which_sbottom);
      } else if ((p1 == lowMSSM_info::g2 && p2 == lowMSSM_info::vu) || 
                 (p1 == lowMSSM_info::vu && p2 == lowMSSM_info::g2)) {
         return deriv_d2MSbottom2_dg2_dvu(which_sbottom);
      } else if ((p1 == lowMSSM_info::Yd22 && p2 == lowMSSM_info::vu) ||
                 (p1 == lowMSSM_info::vu && p2 == lowMSSM_info::Yd22)) {
         return deriv_d2MSbottom2_dYd22_dvu(which_sbottom);
      } else if ((p1 == lowMSSM_info::mq222 && p2 == lowMSSM_info::vu) ||
                 (p1 == lowMSSM_info::vu && p2 == lowMSSM_info::mq222)) {
         return deriv_d2MSbottom2_dmq222_dvu(which_sbottom);
      } else if ((p1 == lowMSSM_info::md222 && p2 == lowMSSM_info::vu) ||
                 (p1 == lowMSSM_info::vu && p2 == lowMSSM_info::md222)) {
         return deriv_d2MSbottom2_dmd222_dvu(which_sbottom);
      } else if ((p1 == lowMSSM_info::Mu && p2 == lowMSSM_info::vu) ||
                 (p1 == lowMSSM_info::vu && p2 == lowMSSM_info::Mu)) {
         return deriv_d2MSbottom2_dMu_dvu(which_sbottom);
      } else if ((p1 == lowMSSM_info::TYd22 && p2 == lowMSSM_info::vu) ||
                 (p1 == lowMSSM_info::vu && p2 == lowMSSM_info::TYd22)) {
         return deriv_d2MSbottom2_dTYd22_dvu(which_sbottom);
      } else {
         return 0.; //< while not actually zero, all other derivatives are not needed for now
      }
   }

} // namespace flexiblesusy
