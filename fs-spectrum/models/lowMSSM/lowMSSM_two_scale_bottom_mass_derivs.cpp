#include "lowMSSM_ewsb_conditions.hpp"

namespace flexiblesusy {

  double lowMSSM_ewsb_conditions::deriv_dMFbottom2_dvd() const
   {
      double vd = model.get_vd();
      double yb = model.get_Yd(2,2);

      return Sqr(yb) * vd;
   }

   double lowMSSM_ewsb_conditions::deriv_dMFbottom2_dYd22() const
   {
      double vd = model.get_vd();
      double yb = model.get_Yd(2,2);

      return yb * Sqr(vd);
   }

   double lowMSSM_ewsb_conditions::deriv_dMFbottom2_dparam(lowMSSM_info::Parameters p) const
   {
      switch (p) {
      case lowMSSM_info::vd: {
         return deriv_dMFbottom2_dvd();
      }
      case lowMSSM_info::Yd22: {
         return deriv_dMFbottom2_dYd22();
      }
      default: {
         return 0.;
      }
      }
   }
 
   double lowMSSM_ewsb_conditions::deriv_d2MFbottom2_dvd_dvd() const
   {
      double yb = model.get_Yd(2,2);

      return Sqr(yb);
   }

   double lowMSSM_ewsb_conditions::deriv_d2MFbottom2_dYd22_dvd() const
   {
      double vd = model.get_vd();

      return Sqr(vd);
   }

   double lowMSSM_ewsb_conditions::deriv_d2MFbottom2_dYd22_dYd22() const
   {
      double yb = model.get_Yd(2,2);
      double vd = model.get_vd();

      return 2.0 * yb * vd;
   }

   double lowMSSM_ewsb_conditions::deriv_d2MFbottom2_dparam_dparam(lowMSSM_info::Parameters p1, lowMSSM_info::Parameters p2) const
   {
      if (p1 == lowMSSM_info::vd && p2 == lowMSSM_info::vd) {
         return deriv_d2MFbottom2_dvd_dvd();
      } else if (p1 == lowMSSM_info::Yd22 && p2 == lowMSSM_info::Yd22) {
         return deriv_d2MFbottom2_dYd22_dYd22();
      } else if ((p1 == lowMSSM_info::Yd22 && p2 == lowMSSM_info::vd) ||
                 (p1 == lowMSSM_info::vd && p2 == lowMSSM_info::Yd22)) {
         return deriv_d2MFbottom2_dYd22_dvd();
      } else {
         return 0.
      }
   }

} // namespace flexiblesusy
