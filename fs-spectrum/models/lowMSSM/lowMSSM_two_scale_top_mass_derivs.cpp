#include "lowMSSM_ewsb_conditions.hpp"

namespace flexiblesusy {

   double lowMSSM_ewsb_conditions::deriv_dMFtop2_dvu() const
   {
      double vu = model.get_vu();
      double yt = model.get_Yu(2,2);

      return Sqr(yt) * vu;
   }

   double lowMSSM_ewsb_conditions::deriv_dMFtop2_dYu22() const
   {
      double vu = model.get_vu();
      double yt = model.get_Yu(2,2);

      return yt * Sqr(vu);
   }

   double lowMSSM_ewsb_conditions::deriv_dMFtop2_dparam(lowMSSM_info::Parameters p) const
   {
      switch (p) {
      case lowMSSM_info::vu: {
         return deriv_dMFtop2_dvu();
      }
      case lowMSSM_info::Yu22: {
         return deriv_dMFtop2_dYu22();
      }
      default: {
         return 0.;
      }
      }
   }
  
   double lowMSSM_ewsb_conditions::deriv_d2MFtop2_dvu_dvu() const
   {
      double yt = model.get_Yu(2,2);

      return Sqr(yt);
   }

   double lowMSSM_ewsb_conditions::deriv_d2MFtop2_dYu22_dvu() const
   {
      double yt = model.get_Yu(2,2);
      double vu = model.get_vu();

      return 2.0 * yt * vu;
   }

   double lowMSSM_ewsb_conditions::deriv_d2MFtop2_dYu22_dYu22() const
   {
      double vu = model.get_vu();

      return Sqr(vu);
   }

   double lowMSSM_ewsb_conditions::deriv_d2MFtop2_dparam_dparam(lowMSSM_info::Parameters p1, lowMSSM_info::Parameters p2) const
   {
      if (p1 == lowMSSM_info::vu && p2 == lowMSSM_info::vu) {
         return deriv_d2MFtop2_dvu_dvu();
      } else if (p1 == lowMSSM_info::Yu22 && p2 == lowMSSM_info::Yu22) {
         return deriv_d2MFtop2_dYu22_dYu22();
      } else if ((p1 == lowMSSM_info::Yu22 && p2 == lowMSSM_info::vu) ||
                 (p1 == lowMSSM_info::vu && p2 == lowMSSM_info::Yu22)) {
         return deriv_d2MFtop2_dYu22_dvu();
      } else {
         return 0.;
      }
   }
   
} // namespace flexiblesusy
