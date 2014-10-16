#include "lowE6SSM_two_scale_ew_derivs.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

   double lowE6SSM_ew_derivs::deriv_dMFtop2_dvu() const
   {
      double yt = model.get_Yu(2,2);
      double vu = model.get_vu();

      return Sqr(yt) * vu;
   }

   double lowE6SSM_ew_derivs::deriv_dMFtop2_dYu22() const
   {
      double yt = model.get_Yu(2,2);
      double vu = model.get_vu();

      return yt * Sqr(vu);
   }

   double lowE6SSM_ew_derivs::deriv_dMFtop2_dparam(lowE6SSM_info::Parameters p) const
   {
      switch (p) {
      case lowE6SSM_info::vu : {
         return lowE6SSM_ew_derivs::deriv_dMFtop2_dvu();
      }
      case lowE6SSM_info::Yu22 : {
         return lowE6SSM_ew_derivs::deriv_dMFtop2_dYu22();
      }
      default : {
         return 0.;
      }
      }
   }
   
   double lowE6SSM_ew_derivs::deriv_d2MFtop2_dvu_dvu() const
   {
      double yt = model.get_Yu(2,2);

      return Sqr(yt);
   }

   double lowE6SSM_ew_derivs::deriv_d2MFtop2_dYu22_dvu() const
   {
      double yt = model.get_Yu(2,2);
      double vu = model.get_vu();

      return 2.0 * yt * vu;
   }

   double lowE6SSM_ew_derivs::deriv_d2MFtop2_dYu22_dYu22() const
   {
      double vu = model.get_vu();

      return Sqr(vu);
   }

   double lowE6SSM_ew_derivs::deriv_d2MFtop2_dparam_dparam(lowE6SSM_info::Parameters p1, 
                                       lowE6SSM_info::Parameters p2) const
   {
      if (equal_as_unordered_pairs(p1, p2, lowE6SSM_info::vu, lowE6SSM_info::vu)) {
         return deriv_d2MFtop2_dvu_dvu();
      } else if (equal_as_unordered_pairs(p1, p2, lowE6SSM_info::vu, lowE6SSM_info::Yu22)) {
         return deriv_d2MFtop2_dYu22_dvu();
      } else if (equal_as_unordered_pairs(p1, p2, lowE6SSM_info::Yu22, lowE6SSM_info::Yu22)) {
         return deriv_d2MFtop2_dYu22_dYu22();
      } else {
         return 0.;
      }
   }

} // namespace flexiblesusy
