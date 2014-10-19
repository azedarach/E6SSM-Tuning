#include "lowE6SSM_two_scale_ew_derivs.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

   double lowE6SSM_ew_derivs::stop_discriminant() const
   {
      const double QH1p = model.get_input().QH1p;
      const double QH2p = model.get_input().QH2p;
      const double QSp = model.get_input().QSp;
      const double QQp = model.get_input().QQp;
      const double Qup = model.get_input().Qup;
      
      const double g1 = model.get_g1();
      const double g2 = model.get_g2();
      const double gN = model.get_gN();
      const double vd = model.get_vd();
      const double vu = model.get_vu();
      const double vs = model.get_vs();
      const double yt = model.get_Yu(2,2);
      const double Lambdax = model.get_Lambdax();
      const double mq222 = model.get_mq2(2,2);
      const double mu222 = model.get_mu2(2,2);
      const double TYu22 = model.get_TYu(2,2);

      double result = Sqr(mq222 - mu222 + 0.125 * Sqr(g2) * Sqr(vd) - 
                          0.125 * Sqr(g1) * Sqr(vd) - 0.125 * Sqr(g2) *
                          Sqr(vu) + 0.125 * Sqr(g1) * Sqr(vu) + 0.5 * QQp *
                          Sqr(gN) * (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp
                                     * Sqr(vs)) - 0.5 * Qup * Sqr(gN) *
                          (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs)))
         + 2.0 * Sqr(TYu22 * vu - 0.7071067811865475 * vd * vs * yt
                     * Lambdax);

      return result;
   }

   double lowE6SSM_ew_derivs::deriv_dMStop2_dvd(stop_mass which_stop) const
   {
      const double QH1p = model.get_input().QH1p;
      const double QH2p = model.get_input().QH2p;
      const double QSp = model.get_input().QSp;
      const double QQp = model.get_input().QQp;
      const double Qup = model.get_input().Qup;

      const double g1 = model.get_g1();
      const double g2 = model.get_g2();
      const double gN = model.get_gN();
      const double vd = model.get_vd();
      const double vu = model.get_vu();
      const double vs = model.get_vs();
      const double yt = model.get_Yu(2,2);
      const double Lambdax = model.get_Lambdax();
      const double mq222 = model.get_mq2(2,2);
      const double mu222 = model.get_mu2(2,2);
      const double TYu22 = model.get_TYu(2,2);

      const double rt = stop_discriminant();

      double tmp_1 = 0.25 * Sqr(g2) * vd + 0.15 * Sqr(g1) * vd + QH1p * QQp * 
         Sqr(gN) * vd + QH1p * Qup * Sqr(gN) * vd;
      double tmp_2 = (2.0 * (0.25 * Sqr(g2) * vd - 0.25 * Sqr(g1) * vd + QH1p 
                             * QQp * Sqr(gN) * vd - QH1p * Qup * Sqr(gN) * vd)
                      * (mq222 - mu222 + 0.125 * Sqr(g2) * Sqr(vd) - 0.125 * 
                         Sqr(g1) * Sqr(vd) - 0.125 * Sqr(g2) * Sqr(vu) + 0.125 *
                         Sqr(g1) * Sqr(vu) + 0.5 * QQp * Sqr(gN) *
                         (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
                         - 0.5 * Qup * Sqr(gN) * (QH1p * Sqr(vd) + QH2p * Sqr(vu)
                                                   + QSp * Sqr(vs))) 
                      - 2.8284271247461903 * vs * yt * Lambdax * 
                      (TYu22 * vu - 0.7071067811865475 * vd * vs * yt * Lambdax))
         / (2.0 * Sqrt(rt));

      if (which_stop == stop_mass::mstop_1) {
         return 0.5 * (tmp_1 - tmp_2);
      } else {
         return 0.5 * (tmp_1 + tmp_2);
      }
      
   }
   
   double lowE6SSM_ew_derivs::deriv_dMStop2_dvu(stop_mass which_stop) const
   {
      const double QH1p = model.get_input().QH1p;
      const double QH2p = model.get_input().QH2p;
      const double QSp = model.get_input().QSp;
      const double QQp = model.get_input().QQp;
      const double Qup = model.get_input().Qup;

      const double g1 = model.get_g1();
      const double g2 = model.get_g2();
      const double gN = model.get_gN();
      const double vd = model.get_vd();
      const double vu = model.get_vu();
      const double vs = model.get_vs();
      const double yt = model.get_Yu(2,2);
      const double Lambdax = model.get_Lambdax();
      const double mq222 = model.get_mq2(2,2);
      const double mu222 = model.get_mu2(2,2);
      const double TYu22 = model.get_TYu(2,2);

      const double rt = stop_discriminant();

      double tmp_1 = -0.25 * Sqr(g2) * vu - 0.15 * Sqr(g1) * vu + QH2p * QQp * 
         Sqr(gN) * vu + QH2p * Qup * Sqr(gN) * vu + 2.0 * Sqr(yt) * vu;
      double tmp_2 = (2.0 * (-0.25 * Sqr(g2) * vu + 0.25 * Sqr(g1) * vu + QH2p 
                             * QQp * Sqr(gN) * vu - QH2p * Qup * Sqr(gN) * vu)
                      * (mq222 - mu222 + 0.125 * Sqr(g2) * Sqr(vd) - 0.125 * 
                         Sqr(g1) * Sqr(vd) - 0.125 * Sqr(g2) * Sqr(vu) + 0.125 *
                         Sqr(g1) * Sqr(vu) + 0.5 * QQp * Sqr(gN) *
                         (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
                         - 0.5 * Qup * Sqr(gN) * (QH1p * Sqr(vd) + QH2p * Sqr(vu)
                                                   + QSp * Sqr(vs))) 
                      + 4.0 * TYu22 * (TYu22 * vu - 0.7071067811865475 * vd * vs * yt
                                       * Lambdax))
         / (2.0 * Sqrt(rt));

      if (which_stop == stop_mass::mstop_1) {
         return 0.5 * (tmp_1 - tmp_2);
      } else {
         return 0.5 * (tmp_1 + tmp_2);
      }
   }
   
   double lowE6SSM_ew_derivs::deriv_dMStop2_dvs(stop_mass which_stop) const
   {
      const double QH1p = model.get_input().QH1p;
      const double QH2p = model.get_input().QH2p;
      const double QSp = model.get_input().QSp;
      const double QQp = model.get_input().QQp;
      const double Qup = model.get_input().Qup;

      const double g1 = model.get_g1();
      const double g2 = model.get_g2();
      const double gN = model.get_gN();
      const double vd = model.get_vd();
      const double vu = model.get_vu();
      const double vs = model.get_vs();
      const double yt = model.get_Yu(2,2);
      const double Lambdax = model.get_Lambdax();
      const double mq222 = model.get_mq2(2,2);
      const double mu222 = model.get_mu2(2,2);
      const double TYu22 = model.get_TYu(2,2);

      const double rt = stop_discriminant();

      double tmp_1 = QSp * QQp * Sqr(gN) * vs + QSp * Qup * Sqr(gN) * vs;
      double tmp_2 = (2.0 * (QSp * QQp * Sqr(gN) * vs - QSp * Qup * Sqr(gN) * vs)
                      * (mq222 - mu222 + 0.125 * Sqr(g2) * Sqr(vd) - 0.125 * 
                         Sqr(g1) * Sqr(vd) - 0.125 * Sqr(g2) * Sqr(vu) + 0.125 *
                         Sqr(g1) * Sqr(vu) + 0.5 * QQp * Sqr(gN) *
                         (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
                         - 0.5 * Qup * Sqr(gN) * (QH1p * Sqr(vd) + QH2p * Sqr(vu)
                                                  + QSp * Sqr(vs))) 
                      - 2.8284271247461903 * vd * yt * Lambdax * 
                      (TYu22 * vu - 0.7071067811865475 * vd * vs * yt * Lambdax))
         / (2.0 * Sqrt(rt));

      if (which_stop == stop_mass::mstop_1) {
         return 0.5 * (tmp_1 - tmp_2);
      } else {
         return 0.5 * (tmp_1 + tmp_2);
      }
   }
   
   double lowE6SSM_ew_derivs::deriv_dMStop2_dg1(stop_mass which_stop) const
   {
      const double QH1p = model.get_input().QH1p;
      const double QH2p = model.get_input().QH2p;
      const double QSp = model.get_input().QSp;
      const double QQp = model.get_input().QQp;
      const double Qup = model.get_input().Qup;

      const double g1 = model.get_g1();
      const double g2 = model.get_g2();
      const double gN = model.get_gN();
      const double vd = model.get_vd();
      const double vu = model.get_vu();
      const double vs = model.get_vs();
      const double mq222 = model.get_mq2(2,2);
      const double mu222 = model.get_mu2(2,2);

      const double rt = stop_discriminant();

      double tmp_1 = 0.15 * g1 * Sqr(vd) - 0.15 * g1 * Sqr(vu);
      double tmp_2 = (-0.25 * g1 * (Sqr(vd) - Sqr(vu))
                      * (mq222 - mu222 + 0.125 * Sqr(g2) * Sqr(vd) - 0.125 * 
                         Sqr(g1) * Sqr(vd) - 0.125 * Sqr(g2) * Sqr(vu) + 0.125 *
                         Sqr(g1) * Sqr(vu) + 0.5 * QQp * Sqr(gN) *
                         (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
                         - 0.5 * Qup * Sqr(gN) * (QH1p * Sqr(vd) + QH2p * Sqr(vu)
                                                  + QSp * Sqr(vs)))) / Sqrt(rt);
      
      if (which_stop == stop_mass::mstop_1) {
         return 0.5 * (tmp_1 - tmp_2);
      } else {
         return 0.5 * (tmp_1 + tmp_2);
      }
   }
   
   double lowE6SSM_ew_derivs::deriv_dMStop2_dg2(stop_mass which_stop) const
   {
      const double QH1p = model.get_input().QH1p;
      const double QH2p = model.get_input().QH2p;
      const double QSp = model.get_input().QSp;
      const double QQp = model.get_input().QQp;
      const double Qup = model.get_input().Qup;

      const double g1 = model.get_g1();
      const double g2 = model.get_g2();
      const double gN = model.get_gN();
      const double vd = model.get_vd();
      const double vu = model.get_vu();
      const double vs = model.get_vs();
      const double mq222 = model.get_mq2(2,2);
      const double mu222 = model.get_mu2(2,2);

      const double rt = stop_discriminant();

      double tmp_1 = 0.25 * g2 * Sqr(vd) - 0.25 * g2 * Sqr(vu);
      double tmp_2 = (0.25 * g2 * (Sqr(vd) - Sqr(vu))
                      * (mq222 - mu222 + 0.125 * Sqr(g2) * Sqr(vd) - 0.125 * 
                         Sqr(g1) * Sqr(vd) - 0.125 * Sqr(g2) * Sqr(vu) + 0.125 *
                         Sqr(g1) * Sqr(vu) + 0.5 * QQp * Sqr(gN) *
                         (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
                         - 0.5 * Qup * Sqr(gN) * (QH1p * Sqr(vd) + QH2p * Sqr(vu)
                                                  + QSp * Sqr(vs)))) / Sqrt(rt);
      
      if (which_stop == stop_mass::mstop_1) {
         return 0.5 * (tmp_1 - tmp_2);
      } else {
         return 0.5 * (tmp_1 + tmp_2);
      }   
   }
   
   double lowE6SSM_ew_derivs::deriv_dMStop2_dgN(stop_mass which_stop) const
   {
      const double QH1p = model.get_input().QH1p;
      const double QH2p = model.get_input().QH2p;
      const double QSp = model.get_input().QSp;
      const double QQp = model.get_input().QQp;
      const double Qup = model.get_input().Qup;

      const double g1 = model.get_g1();
      const double g2 = model.get_g2();
      const double gN = model.get_gN();
      const double vd = model.get_vd();
      const double vu = model.get_vu();
      const double vs = model.get_vs();
      const double mq222 = model.get_mq2(2,2);
      const double mu222 = model.get_mu2(2,2);

      const double rt = stop_discriminant();

      double tmp_1 = QQp * gN * (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
         + Qup * gN * (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs));
      double tmp_2 = ((QQp * gN * (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
                       - Qup * gN * (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs)))
                      * (mq222 - mu222 + 0.125 * Sqr(g2) * Sqr(vd) - 0.125 * 
                         Sqr(g1) * Sqr(vd) - 0.125 * Sqr(g2) * Sqr(vu) + 0.125 *
                         Sqr(g1) * Sqr(vu) + 0.5 * QQp * Sqr(gN) *
                         (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
                         - 0.5 * Qup * Sqr(gN) * (QH1p * Sqr(vd) + QH2p * Sqr(vu)
                                                  + QSp * Sqr(vs)))) / Sqrt(rt);
      
      if (which_stop == stop_mass::mstop_1) {
         return 0.5 * (tmp_1 - tmp_2);
      } else {
         return 0.5 * (tmp_1 + tmp_2);
      }      
   }
   
   double lowE6SSM_ew_derivs::deriv_dMStop2_dYu22(stop_mass which_stop) const
   {
      const double vd = model.get_vd();
      const double vu = model.get_vu();
      const double vs = model.get_vs();
      const double yt = model.get_Yu(2,2);
      const double Lambdax = model.get_Lambdax();
      const double TYu22 = model.get_TYu(2,2);

      const double rt = stop_discriminant();

      double tmp_1 = 2.0 * Sqr(vu) * yt;
      double tmp_2 = (-1.4142135623730951 * vd * vs * Lambdax * 
                      (TYu22 * vu - 0.7071067811865475 * vd * vs * yt 
                       * Lambdax)) / Sqrt(rt);

      if (which_stop == stop_mass::mstop_1) {
         return 0.5 * (tmp_1 - tmp_2);
      } else {
         return 0.5 * (tmp_1 + tmp_2);
      }
   }
   
   double lowE6SSM_ew_derivs::deriv_dMStop2_dmq222(stop_mass which_stop) const
   {
      const double QH1p = model.get_input().QH1p;
      const double QH2p = model.get_input().QH2p;
      const double QSp = model.get_input().QSp;
      const double QQp = model.get_input().QQp;
      const double Qup = model.get_input().Qup;

      const double g1 = model.get_g1();
      const double g2 = model.get_g2();
      const double gN = model.get_gN();
      const double vd = model.get_vd();
      const double vu = model.get_vu();
      const double vs = model.get_vs();
      const double mq222 = model.get_mq2(2,2);
      const double mu222 = model.get_mu2(2,2);

      const double rt = stop_discriminant();

      double tmp_1 = 1.;
      double tmp_2 = (mq222 - mu222 + 0.125 * Sqr(g2) * Sqr(vd) - 0.125 * 
                      Sqr(g1) * Sqr(vd) - 0.125 * Sqr(g2) * Sqr(vu) + 0.125 *
                      Sqr(g1) * Sqr(vu) + 0.5 * QQp * Sqr(gN) *
                      (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
                      - 0.5 * Qup * Sqr(gN) * (QH1p * Sqr(vd) + QH2p * Sqr(vu)
                                               + QSp * Sqr(vs))) / Sqrt(rt);

      if (which_stop == stop_mass::mstop_1) {
         return 0.5 * (tmp_1 - tmp_2);
      } else {
         return 0.5 * (tmp_1 + tmp_2);
      } 
   }
   
   double lowE6SSM_ew_derivs::deriv_dMStop2_dmu222(stop_mass which_stop) const
   {
      const double QH1p = model.get_input().QH1p;
      const double QH2p = model.get_input().QH2p;
      const double QSp = model.get_input().QSp;
      const double QQp = model.get_input().QQp;
      const double Qup = model.get_input().Qup;

      const double g1 = model.get_g1();
      const double g2 = model.get_g2();
      const double gN = model.get_gN();
      const double vd = model.get_vd();
      const double vu = model.get_vu();
      const double vs = model.get_vs();
      const double mq222 = model.get_mq2(2,2);
      const double mu222 = model.get_mu2(2,2);

      const double rt = stop_discriminant();

      double tmp_1 = 1.;
      double tmp_2 = -(mq222 - mu222 + 0.125 * Sqr(g2) * Sqr(vd) - 0.125 * 
                       Sqr(g1) * Sqr(vd) - 0.125 * Sqr(g2) * Sqr(vu) + 0.125 *
                       Sqr(g1) * Sqr(vu) + 0.5 * QQp * Sqr(gN) *
                       (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
                       - 0.5 * Qup * Sqr(gN) * (QH1p * Sqr(vd) + QH2p * Sqr(vu)
                                                + QSp * Sqr(vs))) / Sqrt(rt);

      if (which_stop == stop_mass::mstop_1) {
         return 0.5 * (tmp_1 - tmp_2);
      } else {
         return 0.5 * (tmp_1 + tmp_2);
      }
   }
   
   double lowE6SSM_ew_derivs::deriv_dMStop2_dLambdax(stop_mass which_stop) const
   {
      const double vd = model.get_vd();
      const double vu = model.get_vu();
      const double vs = model.get_vs();
      const double yt = model.get_Yu(2,2);
      const double Lambdax = model.get_Lambdax();
      const double TYu22 = model.get_TYu(2,2);

      const double rt = stop_discriminant();

      double tmp_1 = 0.;
      double tmp_2 = -(1.4142135623730951 * vd * vs * yt *
                       (TYu22 * vu - 0.7071067811865475 * vd * vs * yt
                        * Lambdax)) / Sqrt(rt);

      if (which_stop == stop_mass::mstop_1) {
         return 0.5 * (tmp_1 - tmp_2);
      } else {
         return 0.5 * (tmp_1 + tmp_2);
      }
   }
   
   double lowE6SSM_ew_derivs::deriv_dMStop2_dTYu22(stop_mass which_stop) const
   {
      const double vd = model.get_vd();
      const double vu = model.get_vu();
      const double vs = model.get_vs();
      const double yt = model.get_Yu(2,2);
      const double Lambdax = model.get_Lambdax();
      const double TYu22 = model.get_TYu(2,2);

      const double rt = stop_discriminant();

      double tmp_1 = 0.;
      double tmp_2 = (2.0 * vu * (TYu22 * vu - 0.7071067811865475 * vd * vs * yt
                       * Lambdax)) / Sqrt(rt);

      if (which_stop == stop_mass::mstop_1) {
         return 0.5 * (tmp_1 - tmp_2);
      } else {
         return 0.5 * (tmp_1 + tmp_2);
      }
   }

   double lowE6SSM_ew_derivs::deriv_dMStop2_dparam(stop_mass which_stop, lowE6SSM_info::Parameters p) const
   {
      switch (p) {
      case lowE6SSM_info::vd : {
         return deriv_dMStop2_dvd(which_stop);
      }
      case lowE6SSM_info::vu : {
         return deriv_dMStop2_dvu(which_stop);
      }
      case lowE6SSM_info::vs : {
         return deriv_dMStop2_dvs(which_stop);
      }
      case lowE6SSM_info::g1 : {
         return deriv_dMStop2_dg1(which_stop);
      }
      case lowE6SSM_info::g2 : {
         return deriv_dMStop2_dg2(which_stop);
      }
      case lowE6SSM_info::gN : {
         return deriv_dMStop2_dgN(which_stop);
      }
      case lowE6SSM_info::Yu22 : {
         return deriv_dMStop2_dYu22(which_stop);
      }
      case lowE6SSM_info::mq222 : {
         return deriv_dMStop2_dmq222(which_stop);
      }
      case lowE6SSM_info::mu222 : {
         return deriv_dMStop2_dmu222(which_stop);
      }
      case lowE6SSM_info::Lambdax : {
         return deriv_dMStop2_dLambdax(which_stop);
      }
      case lowE6SSM_info::TYu22 : {
         return deriv_dMStop2_dTYu22(which_stop);
      }
      default : {
         return 0.;
      }

      }
   }

   double lowE6SSM_ew_derivs::deriv_d2MStop2_dvd_dvd(stop_mass which_stop) const
   {
      const double QH1p = model.get_input().QH1p;
      const double QH2p = model.get_input().QH2p;
      const double QSp = model.get_input().QSp;
      const double QQp = model.get_input().QQp;
      const double Qup = model.get_input().Qup;

      const double g1 = model.get_g1();
      const double g2 = model.get_g2();
      const double gN = model.get_gN();
      const double vd = model.get_vd();
      const double vu = model.get_vu();
      const double vs = model.get_vs();
      const double yt = model.get_Yu(2,2);
      const double Lambdax = model.get_Lambdax();
      const double mq222 = model.get_mq2(2,2);
      const double mu222 = model.get_mu2(2,2);
      const double TYu22 = model.get_TYu(2,2);

      const double rt = stop_discriminant();

      double tmp_1 = 0.25 * Sqr(g2) + 0.15 * Sqr(g1) + QH1p * QQp * Sqr(gN)
         + QH1p * Qup * Sqr(gN);
      double tmp_2 = (2.0 * Sqr(0.25 * Sqr(g2) * vd - 0.25 * Sqr(g1) * vd + QH1p * QQp * Sqr(gN) * vd
                                - QH1p * Qup * Sqr(gN) * vd) + 
                      2.0 * (0.25 * Sqr(g2) - 0.25 * Sqr(g1) + QH1p * QQp * Sqr(gN) - QH1p * Qup * Sqr(gN))
                      * (mq222 - mu222 + 0.125 * Sqr(g2) * Sqr(vd) - 0.125 * 
                         Sqr(g1) * Sqr(vd) - 0.125 * Sqr(g2) * Sqr(vu) + 0.125 *
                         Sqr(g1) * Sqr(vu) + 0.5 * QQp * Sqr(gN) *
                         (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
                         - 0.5 * Qup * Sqr(gN) * (QH1p * Sqr(vd) + QH2p * Sqr(vu)
                                                  + QSp * Sqr(vs))) + 2.0 * Sqr(vs) * Sqr(yt)
                      * Sqr(Lambdax))/ (2.0 * Sqrt(rt));
      double tmp_3 = Sqr(2.0 * (0.25 * Sqr(g2) * vd - 0.25 * Sqr(g1) * vd 
                                + QH1p * QQp * Sqr(gN) * vd - QH1p * Qup * Sqr(gN)
                                * vd) * (mq222 - mu222 + 0.125 * Sqr(g2) * Sqr(vd) - 0.125 * 
                                         Sqr(g1) * Sqr(vd) - 0.125 * Sqr(g2) * Sqr(vu) + 0.125 *
                                         Sqr(g1) * Sqr(vu) + 0.5 * QQp * Sqr(gN) *
                                         (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
                                         - 0.5 * Qup * Sqr(gN) * (QH1p * Sqr(vd) + QH2p * Sqr(vu)
                                                                  + QSp * Sqr(vs)))
                         - 2.8284271247461903 * vs * yt * Lambdax * 
                         (TYu22 * vu - 0.7071067811865475 * vd * vs * yt * Lambdax)) 
         / (4.0 * rt * Sqrt(rt));

      if (which_stop == stop_mass::mstop_1) {
         return 0.5 * (tmp_1 - (tmp_2 - tmp_3));
      } else {
         return 0.5 * (tmp_1 + (tmp_2 - tmp_3));
      }
   }

   double lowE6SSM_ew_derivs::deriv_d2MStop2_dvd_dvu(stop_mass which_stop) const
   {
      const double QH1p = model.get_input().QH1p;
      const double QH2p = model.get_input().QH2p;
      const double QSp = model.get_input().QSp;
      const double QQp = model.get_input().QQp;
      const double Qup = model.get_input().Qup;

      const double g1 = model.get_g1();
      const double g2 = model.get_g2();
      const double gN = model.get_gN();
      const double vd = model.get_vd();
      const double vu = model.get_vu();
      const double vs = model.get_vs();
      const double yt = model.get_Yu(2,2);
      const double Lambdax = model.get_Lambdax();
      const double mq222 = model.get_mq2(2,2);
      const double mu222 = model.get_mu2(2,2);
      const double TYu22 = model.get_TYu(2,2);

      const double rt = stop_discriminant();

      double tmp_1 = 0.;
      double tmp_2 = (2.0 * (0.25 * Sqr(g2) * vd - 0.25 * Sqr(g1) * vd + QH1p * QQp *
                             Sqr(gN) * vd - QH1p * Qup * Sqr(gN) * vd) * 
                      (-0.25 * Sqr(g2) * vu + 0.25 * Sqr(g1) * vu + QH2p * QQp * Sqr(gN)
                       * vu - QH2p * Qup * Sqr(gN) * vu) - 2.8284271247461903 * TYu22 * vs
                      * yt * Lambdax)/ (2.0 * Sqrt(rt));
      double tmp_3 = ((2.0 * (-0.25 * Sqr(g2) * vu + 0.25 * Sqr(g1) * vu + QH2p * QQp * Sqr(gN)
                             * vu - QH2p * Qup * Sqr(gN) * vu) * 
                      (mq222 - mu222 + 0.125 * Sqr(g2) * Sqr(vd) - 0.125 * 
                       Sqr(g1) * Sqr(vd) - 0.125 * Sqr(g2) * Sqr(vu) + 0.125 *
                       Sqr(g1) * Sqr(vu) + 0.5 * QQp * Sqr(gN) *
                       (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
                       - 0.5 * Qup * Sqr(gN) * (QH1p * Sqr(vd) + QH2p * Sqr(vu)
                                                + QSp * Sqr(vs))) + 4.0 * TYu22 * 
                      (TYu22 * vu - 0.7071067811865475 * vd * vs * yt * Lambdax))
         * (2.0 * (0.25 * Sqr(g2) * vd - 0.25 * Sqr(g1) * vd + QH1p * QQp * Sqr(gN) * vd
                   - QH1p * Qup * Sqr(gN) * vd) * 
            (mq222 - mu222 + 0.125 * Sqr(g2) * Sqr(vd) - 0.125 * 
             Sqr(g1) * Sqr(vd) - 0.125 * Sqr(g2) * Sqr(vu) + 0.125 *
             Sqr(g1) * Sqr(vu) + 0.5 * QQp * Sqr(gN) *
             (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
                       - 0.5 * Qup * Sqr(gN) * (QH1p * Sqr(vd) + QH2p * Sqr(vu)
                                                + QSp * Sqr(vs))) 
            - 2.8284271247461903 * vs * yt * Lambdax * 
            (TYu22 * vu - 0.7071067811865475 * vd * vs * yt * Lambdax)))
         / (4.0 * rt * Sqrt(rt));

      if (which_stop == stop_mass::mstop_1) {
         return 0.5 * (tmp_1 - (tmp_2 - tmp_3));
      } else {
         return 0.5 * (tmp_1 + (tmp_2 - tmp_3));
      }
   }

   double lowE6SSM_ew_derivs::deriv_d2MStop2_dvd_dvs(stop_mass which_stop) const
   {
      const double QH1p = model.get_input().QH1p;
      const double QH2p = model.get_input().QH2p;
      const double QSp = model.get_input().QSp;
      const double QQp = model.get_input().QQp;
      const double Qup = model.get_input().Qup;

      const double g1 = model.get_g1();
      const double g2 = model.get_g2();
      const double gN = model.get_gN();
      const double vd = model.get_vd();
      const double vu = model.get_vu();
      const double vs = model.get_vs();
      const double yt = model.get_Yu(2,2);
      const double Lambdax = model.get_Lambdax();
      const double mq222 = model.get_mq2(2,2);
      const double mu222 = model.get_mu2(2,2);
      const double TYu22 = model.get_TYu(2,2);

      const double rt = stop_discriminant();

      double tmp_1 = 0.;
      double tmp_2 = (2.0 * (0.25 * Sqr(g2) * vd - 0.25 * Sqr(g1) * vd + QH1p * QQp * Sqr(gN)
                             * vd - QH1p * Qup * Sqr(gN) * vd) * 
                      (QSp * QQp * Sqr(gN) * vs - QSp * Qup * Sqr(gN) * vs) 
                      + 2.0 * vd * vs * Sqr(yt) * Sqr(Lambdax) 
                      - 2.8284271247461903 * yt * Lambdax * 
                      (TYu22 * vu - 0.7071067811865475 * vd * vs * yt * Lambdax)) / (2.0 * Sqrt(rt));
      double tmp_3 = ((2.0 * (QSp * QQp * Sqr(gN) * vs - QSp * Qup * Sqr(gN) * vs) *
                       (mq222 - mu222 + 0.125 * Sqr(g2) * Sqr(vd) - 0.125 * 
                        Sqr(g1) * Sqr(vd) - 0.125 * Sqr(g2) * Sqr(vu) + 0.125 *
                        Sqr(g1) * Sqr(vu) + 0.5 * QQp * Sqr(gN) *
                        (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
                        - 0.5 * Qup * Sqr(gN) * (QH1p * Sqr(vd) + QH2p * Sqr(vu)
                                                 + QSp * Sqr(vs)))
                       - 2.8284271247461903 * vd * yt * Lambdax *
                       (TYu22 * vu - 0.7071067811865475 * vd * vs * yt * Lambdax)) *
                      (2.0 * (0.25 * Sqr(g2) * vd - 0.25 * Sqr(g1) * vd + QH1p * QQp * Sqr(gN) * vd
                              - QH1p * Qup * Sqr(gN) * vd) * 
                       (mq222 - mu222 + 0.125 * Sqr(g2) * Sqr(vd) - 0.125 * 
                        Sqr(g1) * Sqr(vd) - 0.125 * Sqr(g2) * Sqr(vu) + 0.125 *
                        Sqr(g1) * Sqr(vu) + 0.5 * QQp * Sqr(gN) *
                        (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
                        - 0.5 * Qup * Sqr(gN) * (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs)))
                       - 2.8284271247461903 * vs * yt * Lambdax *
                       (TYu22 * vu - 0.7071067811865475 * vd * vs * yt * Lambdax)))
         / (4.0 * rt * Sqrt(rt));

      if (which_stop == stop_mass::mstop_1) {
         return 0.5 * (tmp_1 - (tmp_2 - tmp_3));
      } else {
         return 0.5 * (tmp_1 + (tmp_2 - tmp_3));
      }
   }

   double lowE6SSM_ew_derivs::deriv_d2MStop2_dg1_dvd(stop_mass which_stop) const
   {
      const double QH1p = model.get_input().QH1p;
      const double QH2p = model.get_input().QH2p;
      const double QSp = model.get_input().QSp;
      const double QQp = model.get_input().QQp;
      const double Qup = model.get_input().Qup;

      const double g1 = model.get_g1();
      const double g2 = model.get_g2();
      const double gN = model.get_gN();
      const double vd = model.get_vd();
      const double vu = model.get_vu();
      const double vs = model.get_vs();
      const double yt = model.get_Yu(2,2);
      const double Lambdax = model.get_Lambdax();
      const double mq222 = model.get_mq2(2,2);
      const double mu222 = model.get_mu2(2,2);
      const double TYu22 = model.get_TYu(2,2);

      const double rt = stop_discriminant();

      double tmp_1 = 0.3 * g1 * vd;
      double tmp_2 = (-0.5 * g1 * (0.25 * Sqr(g2) * vd - 0.25 * Sqr(g1) * vd + QH1p * QQp * Sqr(gN)
                                   * vd - QH1p * Qup * Sqr(gN) * vd) * (Sqr(vd) - Sqr(vu)) - g1 * vd
                      * (mq222 - mu222 + 0.125 * Sqr(g2) * Sqr(vd) - 0.125 * 
                         Sqr(g1) * Sqr(vd) - 0.125 * Sqr(g2) * Sqr(vu) + 0.125 *
                         Sqr(g1) * Sqr(vu) + 0.5 * QQp * Sqr(gN) *
                         (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
                         - 0.5 * Qup * Sqr(gN) * (QH1p * Sqr(vd) + QH2p * Sqr(vu)
                                                  + QSp * Sqr(vs)))) / (2.0 * Sqrt(rt));
      double tmp_3 = -(0.125 * g1 * (Sqr(vd) - Sqr(vu)) * (mq222 - mu222 + 0.125 * Sqr(g2) * Sqr(vd) - 0.125 * 
                                                           Sqr(g1) * Sqr(vd) - 0.125 * Sqr(g2) * Sqr(vu) + 0.125 *
                                                           Sqr(g1) * Sqr(vu) + 0.5 * QQp * Sqr(gN) *
                                                           (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
                                                           - 0.5 * Qup * Sqr(gN) * (QH1p * Sqr(vd) + QH2p * Sqr(vu)
                                                                                    + QSp * Sqr(vs)))
                       * (2.0 * (0.25 * Sqr(g2) * vd - 0.25 * Sqr(g1) * vd + QH1p * QQp * Sqr(gN) * vd 
                                 - QH1p * Qup * Sqr(gN) * vd) * (mq222 - mu222 + 0.125 * Sqr(g2) * Sqr(vd) - 0.125 * 
                                                                 Sqr(g1) * Sqr(vd) - 0.125 * Sqr(g2) * Sqr(vu) + 0.125 *
                                                                 Sqr(g1) * Sqr(vu) + 0.5 * QQp * Sqr(gN) *
                                                                 (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
                                                                 - 0.5 * Qup * Sqr(gN) * (QH1p * Sqr(vd) + QH2p * 
                                                                                          Sqr(vu) + QSp * Sqr(vs)))
                          - 2.8284271247461903 * vs * yt * Lambdax *
                          (TYu22 * vu - 0.7071067811865475 * vd * vs * yt * Lambdax)))
         / (rt * Sqrt(rt));

      if (which_stop == stop_mass::mstop_1) {
         return 0.5 * (tmp_1 - (tmp_2 - tmp_3));
      } else {
         return 0.5 * (tmp_1 + (tmp_2 - tmp_3));
      }
   }

   double lowE6SSM_ew_derivs::deriv_d2MStop2_dg2_dvd(stop_mass which_stop) const
   {
      const double QH1p = model.get_input().QH1p;
      const double QH2p = model.get_input().QH2p;
      const double QSp = model.get_input().QSp;
      const double QQp = model.get_input().QQp;
      const double Qup = model.get_input().Qup;

      const double g1 = model.get_g1();
      const double g2 = model.get_g2();
      const double gN = model.get_gN();
      const double vd = model.get_vd();
      const double vu = model.get_vu();
      const double vs = model.get_vs();
      const double yt = model.get_Yu(2,2);
      const double Lambdax = model.get_Lambdax();
      const double mq222 = model.get_mq2(2,2);
      const double mu222 = model.get_mu2(2,2);
      const double TYu22 = model.get_TYu(2,2);

      const double rt = stop_discriminant();

      double tmp_1 = 0.5 * g2 * vd;
      double tmp_2 = (0.5 * g2 * (0.25 * Sqr(g2) * vd - 0.25 * Sqr(g1) * vd + QH1p * QQp * Sqr(gN)
                                   * vd - QH1p * Qup * Sqr(gN) * vd) * (Sqr(vd) - Sqr(vu)) + g2 * vd
                      * (mq222 - mu222 + 0.125 * Sqr(g2) * Sqr(vd) - 0.125 * 
                         Sqr(g1) * Sqr(vd) - 0.125 * Sqr(g2) * Sqr(vu) + 0.125 *
                         Sqr(g1) * Sqr(vu) + 0.5 * QQp * Sqr(gN) *
                         (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
                         - 0.5 * Qup * Sqr(gN) * (QH1p * Sqr(vd) + QH2p * Sqr(vu)
                                                  + QSp * Sqr(vs)))) / (2.0 * Sqrt(rt));
      double tmp_3 = (0.125 * g2 * (Sqr(vd) - Sqr(vu)) * (mq222 - mu222 + 0.125 * Sqr(g2) * Sqr(vd) - 0.125 * 
                                                           Sqr(g1) * Sqr(vd) - 0.125 * Sqr(g2) * Sqr(vu) + 0.125 *
                                                           Sqr(g1) * Sqr(vu) + 0.5 * QQp * Sqr(gN) *
                                                           (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
                                                           - 0.5 * Qup * Sqr(gN) * (QH1p * Sqr(vd) + QH2p * Sqr(vu)
                                                                                    + QSp * Sqr(vs)))
                       * (2.0 * (0.25 * Sqr(g2) * vd - 0.25 * Sqr(g1) * vd + QH1p * QQp * Sqr(gN) * vd 
                                 - QH1p * Qup * Sqr(gN) * vd) * (mq222 - mu222 + 0.125 * Sqr(g2) * Sqr(vd) - 0.125 * 
                                                                 Sqr(g1) * Sqr(vd) - 0.125 * Sqr(g2) * Sqr(vu) + 0.125 *
                                                                 Sqr(g1) * Sqr(vu) + 0.5 * QQp * Sqr(gN) *
                                                                 (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
                                                                 - 0.5 * Qup * Sqr(gN) * (QH1p * Sqr(vd) + QH2p * 
                                                                                          Sqr(vu) + QSp * Sqr(vs)))
                          - 2.8284271247461903 * vs * yt * Lambdax *
                          (TYu22 * vu - 0.7071067811865475 * vd * vs * yt * Lambdax)))
         / (rt * Sqrt(rt));

      if (which_stop == stop_mass::mstop_1) {
         return 0.5 * (tmp_1 - (tmp_2 - tmp_3));
      } else {
         return 0.5 * (tmp_1 + (tmp_2 - tmp_3));
      }
   }

   double lowE6SSM_ew_derivs::deriv_d2MStop2_dgN_dvd(stop_mass which_stop) const
   {
      const double QH1p = model.get_input().QH1p;
      const double QH2p = model.get_input().QH2p;
      const double QSp = model.get_input().QSp;
      const double QQp = model.get_input().QQp;
      const double Qup = model.get_input().Qup;

      const double g1 = model.get_g1();
      const double g2 = model.get_g2();
      const double gN = model.get_gN();
      const double vd = model.get_vd();
      const double vu = model.get_vu();
      const double vs = model.get_vs();
      const double yt = model.get_Yu(2,2);
      const double Lambdax = model.get_Lambdax();
      const double mq222 = model.get_mq2(2,2);
      const double mu222 = model.get_mu2(2,2);
      const double TYu22 = model.get_TYu(2,2);

      const double rt = stop_discriminant();

      double tmp_1 = 2.0 * QH1p * QQp * gN * vd + 2.0 * QH1p * Qup * gN * vd;
      double tmp_2 = (2.0 * (0.25 * Sqr(g2) * vd - 0.25 * Sqr(g1) * vd + QH1p * QQp * Sqr(gN)
                             * vd - QH1p * Qup * Sqr(gN) * vd) * 
                      (QQp * gN * (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
                       - Qup * gN * (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))) 
                      + 2.0 * (2.0 * QH1p * QQp * gN * vd - 2.0 * QH1p * Qup * gN * vd) *
                      (mq222 - mu222 + 0.125 * Sqr(g2) * Sqr(vd) - 0.125 * 
                       Sqr(g1) * Sqr(vd) - 0.125 * Sqr(g2) * Sqr(vu) + 0.125 *
                       Sqr(g1) * Sqr(vu) + 0.5 * QQp * Sqr(gN) *
                       (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
                       - 0.5 * Qup * Sqr(gN) * (QH1p * Sqr(vd) + QH2p * Sqr(vu)
                                                + QSp * Sqr(vs)))) / (2.0 * Sqrt(rt));
      double tmp_3 = ((QQp * gN * (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
                       - Qup * gN * (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))) * 
                      (mq222 - mu222 + 0.125 * Sqr(g2) * Sqr(vd) - 0.125 * 
                       Sqr(g1) * Sqr(vd) - 0.125 * Sqr(g2) * Sqr(vu) + 0.125 *
                       Sqr(g1) * Sqr(vu) + 0.5 * QQp * Sqr(gN) *
                       (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
                       - 0.5 * Qup * Sqr(gN) * (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs)))
                      * (2.0 * (0.25 * Sqr(g2) * vd - 0.25 * Sqr(g1) * vd + QH1p * QQp * Sqr(gN) * vd 
                                - QH1p * Qup * Sqr(gN) * vd) * (mq222 - mu222 + 0.125 * Sqr(g2) * Sqr(vd) - 0.125 * 
                                                                Sqr(g1) * Sqr(vd) - 0.125 * Sqr(g2) * Sqr(vu) + 0.125 *
                                                                Sqr(g1) * Sqr(vu) + 0.5 * QQp * Sqr(gN) *
                                                                (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
                                                                - 0.5 * Qup * Sqr(gN) * (QH1p * Sqr(vd) + QH2p * 
                                                                                         Sqr(vu) + QSp * Sqr(vs)))
                         - 2.8284271247461903 * vs * yt * Lambdax *
                         (TYu22 * vu - 0.7071067811865475 * vd * vs * yt * Lambdax)))
         / (2.0 * rt * Sqrt(rt));

      if (which_stop == stop_mass::mstop_1) {
         return 0.5 * (tmp_1 - (tmp_2 - tmp_3));
      } else {
         return 0.5 * (tmp_1 + (tmp_2 - tmp_3));
      }
   }

   double lowE6SSM_ew_derivs::deriv_d2MStop2_dYu22_dvd(stop_mass which_stop) const
   {
      const double QH1p = model.get_input().QH1p;
      const double QH2p = model.get_input().QH2p;
      const double QSp = model.get_input().QSp;
      const double QQp = model.get_input().QQp;
      const double Qup = model.get_input().Qup;

      const double g1 = model.get_g1();
      const double g2 = model.get_g2();
      const double gN = model.get_gN();
      const double vd = model.get_vd();
      const double vu = model.get_vu();
      const double vs = model.get_vs();
      const double yt = model.get_Yu(2,2);
      const double Lambdax = model.get_Lambdax();
      const double mq222 = model.get_mq2(2,2);
      const double mu222 = model.get_mu2(2,2);
      const double TYu22 = model.get_TYu(2,2);

      const double rt = stop_discriminant();

      double tmp_1 = 0.;
      double tmp_2 = (2.0 * vd * Sqr(vs) * yt * Sqr(Lambdax) - 
                      2.8284271247461903 * vs * Lambdax * 
                      (TYu22 * vu - 0.7071067811865475 * vd * vs * yt * Lambdax)) / (2.0 * Sqrt(rt));
      double tmp_3 = -0.7071067811865475 * (vd * vs * Lambdax * 
                                            (TYu22 * vu - 0.7071067811865475 * vd * vs * yt * Lambdax)
                                            * (2.0 * (0.25 * Sqr(g2) * vd - 0.25 * Sqr(g1) * vd + QH1p * QQp * Sqr(gN) * vd 
                                                      - QH1p * Qup * Sqr(gN) * vd) * 
                                               (mq222 - mu222 + 0.125 * Sqr(g2) * Sqr(vd) - 0.125 * 
                                                Sqr(g1) * Sqr(vd) - 0.125 * Sqr(g2) * Sqr(vu) + 0.125 *
                                                Sqr(g1) * Sqr(vu) + 0.5 * QQp * Sqr(gN) *
                                                (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
                                                - 0.5 * Qup * Sqr(gN) * (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs)))
                                               - 2.8284271247461903 * vs * yt * Lambdax *
                                               (TYu22 * vu - 0.7071067811865475 * vd * vs * yt * Lambdax)))
         / (rt * Sqrt(rt));

      if (which_stop == stop_mass::mstop_1) {
         return 0.5 * (tmp_1 - (tmp_2 - tmp_3));
      } else {
         return 0.5 * (tmp_1 + (tmp_2 - tmp_3));
      }
   }

   double lowE6SSM_ew_derivs::deriv_d2MStop2_dmq222_dvd(stop_mass which_stop) const
   {
      const double QH1p = model.get_input().QH1p;
      const double QH2p = model.get_input().QH2p;
      const double QSp = model.get_input().QSp;
      const double QQp = model.get_input().QQp;
      const double Qup = model.get_input().Qup;

      const double g1 = model.get_g1();
      const double g2 = model.get_g2();
      const double gN = model.get_gN();
      const double vd = model.get_vd();
      const double vu = model.get_vu();
      const double vs = model.get_vs();
      const double yt = model.get_Yu(2,2);
      const double Lambdax = model.get_Lambdax();
      const double mq222 = model.get_mq2(2,2);
      const double mu222 = model.get_mu2(2,2);
      const double TYu22 = model.get_TYu(2,2);

      const double rt = stop_discriminant();

      double tmp_1 = 0.;
      double tmp_2 = (0.25 * Sqr(g2) * vd - 0.25 * Sqr(g1) * vd + QH1p * QQp * Sqr(gN) * vd
                      - QH1p * Qup * Sqr(gN) * vd) / Sqrt(rt);
      double tmp_3 = ((mq222 - mu222 + 0.125 * Sqr(g2) * Sqr(vd) - 0.125 * 
                       Sqr(g1) * Sqr(vd) - 0.125 * Sqr(g2) * Sqr(vu) + 0.125 *
                       Sqr(g1) * Sqr(vu) + 0.5 * QQp * Sqr(gN) *
                       (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
                       - 0.5 * Qup * Sqr(gN) * (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs)))
                      * (2.0 * (0.25 * Sqr(g2) * vd - 0.25 * Sqr(g1) * vd + QH1p * QQp * Sqr(gN) * vd 
                                - QH1p * Qup * Sqr(gN) * vd) * (mq222 - mu222 + 0.125 * Sqr(g2) * Sqr(vd) - 0.125 * 
                                                                Sqr(g1) * Sqr(vd) - 0.125 * Sqr(g2) * Sqr(vu) + 0.125 *
                                                                Sqr(g1) * Sqr(vu) + 0.5 * QQp * Sqr(gN) *
                                                                (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
                                                                - 0.5 * Qup * Sqr(gN) * (QH1p * Sqr(vd) + QH2p * 
                                                                                         Sqr(vu) + QSp * Sqr(vs)))
                         - 2.8284271247461903 * vs * yt * Lambdax *
                         (TYu22 * vu - 0.7071067811865475 * vd * vs * yt * Lambdax)))
         / (2.0 * rt * Sqrt(rt));

      if (which_stop == stop_mass::mstop_1) {
         return 0.5 * (tmp_1 - (tmp_2 - tmp_3));
      } else {
         return 0.5 * (tmp_1 + (tmp_2 - tmp_3));
      }
   }

   double lowE6SSM_ew_derivs::deriv_d2MStop2_dmu222_dvd(stop_mass which_stop) const
   {
      const double QH1p = model.get_input().QH1p;
      const double QH2p = model.get_input().QH2p;
      const double QSp = model.get_input().QSp;
      const double QQp = model.get_input().QQp;
      const double Qup = model.get_input().Qup;

      const double g1 = model.get_g1();
      const double g2 = model.get_g2();
      const double gN = model.get_gN();
      const double vd = model.get_vd();
      const double vu = model.get_vu();
      const double vs = model.get_vs();
      const double yt = model.get_Yu(2,2);
      const double Lambdax = model.get_Lambdax();
      const double mq222 = model.get_mq2(2,2);
      const double mu222 = model.get_mu2(2,2);
      const double TYu22 = model.get_TYu(2,2);

      const double rt = stop_discriminant();

      double tmp_1 = 0.;
      double tmp_2 = -(0.25 * Sqr(g2) * vd - 0.25 * Sqr(g1) * vd + QH1p * QQp * Sqr(gN) * vd
                       - QH1p * Qup * Sqr(gN) * vd) / Sqrt(rt);
      double tmp_3 = -((mq222 - mu222 + 0.125 * Sqr(g2) * Sqr(vd) - 0.125 * 
                        Sqr(g1) * Sqr(vd) - 0.125 * Sqr(g2) * Sqr(vu) + 0.125 *
                        Sqr(g1) * Sqr(vu) + 0.5 * QQp * Sqr(gN) *
                        (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
                        - 0.5 * Qup * Sqr(gN) * (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs)))
                       * (2.0 * (0.25 * Sqr(g2) * vd - 0.25 * Sqr(g1) * vd + QH1p * QQp * Sqr(gN) * vd 
                                 - QH1p * Qup * Sqr(gN) * vd) * (mq222 - mu222 + 0.125 * Sqr(g2) * Sqr(vd) - 0.125 * 
                                                                 Sqr(g1) * Sqr(vd) - 0.125 * Sqr(g2) * Sqr(vu) + 0.125 *
                                                                 Sqr(g1) * Sqr(vu) + 0.5 * QQp * Sqr(gN) *
                                                                 (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
                                                                 - 0.5 * Qup * Sqr(gN) * (QH1p * Sqr(vd) + QH2p * 
                                                                                          Sqr(vu) + QSp * Sqr(vs)))
                          - 2.8284271247461903 * vs * yt * Lambdax *
                          (TYu22 * vu - 0.7071067811865475 * vd * vs * yt * Lambdax)))
      / (2.0 * rt * Sqrt(rt));

      if (which_stop == stop_mass::mstop_1) {
         return 0.5 * (tmp_1 - (tmp_2 - tmp_3));
      } else {
         return 0.5 * (tmp_1 + (tmp_2 - tmp_3));
      }
   }

   double lowE6SSM_ew_derivs::deriv_d2MStop2_dLambdax_dvd(stop_mass which_stop) const
   {
      const double QH1p = model.get_input().QH1p;
      const double QH2p = model.get_input().QH2p;
      const double QSp = model.get_input().QSp;
      const double QQp = model.get_input().QQp;
      const double Qup = model.get_input().Qup;

      const double g1 = model.get_g1();
      const double g2 = model.get_g2();
      const double gN = model.get_gN();
      const double vd = model.get_vd();
      const double vu = model.get_vu();
      const double vs = model.get_vs();
      const double yt = model.get_Yu(2,2);
      const double Lambdax = model.get_Lambdax();
      const double mq222 = model.get_mq2(2,2);
      const double mu222 = model.get_mu2(2,2);
      const double TYu22 = model.get_TYu(2,2);

      const double rt = stop_discriminant();

      double tmp_1 = 0.;
      double tmp_2 = (2.0 * vd * Sqr(vs) * Sqr(yt) * Lambdax 
                      - 2.8284271247461903 * vs * yt *
                      (TYu22 * vu - 0.7071067811865475 * vd * vs * yt * Lambdax)) / (2.0 * Sqrt(rt));
      double tmp_3 = -0.7071067811865475 * 
      (vd * vs * yt * (TYu22 * vu - 0.7071067811865475 * vd * vs * yt * Lambdax) 
       * (2.0 * (0.25 * Sqr(g2) * vd - 0.25 * Sqr(g1) * vd + QH1p * QQp * Sqr(gN) * vd 
                 - QH1p * Qup * Sqr(gN) * vd) * (mq222 - mu222 + 0.125 * Sqr(g2) * Sqr(vd) - 0.125 * 
                                                 Sqr(g1) * Sqr(vd) - 0.125 * Sqr(g2) * Sqr(vu) + 0.125 *
                                                 Sqr(g1) * Sqr(vu) + 0.5 * QQp * Sqr(gN) *
                                                 (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
                                                 - 0.5 * Qup * Sqr(gN) * (QH1p * Sqr(vd) + QH2p * 
                                                                          Sqr(vu) + QSp * Sqr(vs)))
          - 2.8284271247461903 * vs * yt * Lambdax *
          (TYu22 * vu - 0.7071067811865475 * vd * vs * yt * Lambdax)))
      / (rt * Sqrt(rt));

      if (which_stop == stop_mass::mstop_1) {
         return 0.5 * (tmp_1 - (tmp_2 - tmp_3));
      } else {
         return 0.5 * (tmp_1 + (tmp_2 - tmp_3));
      }
   }

   double lowE6SSM_ew_derivs::deriv_d2MStop2_dTYu22_dvd(stop_mass which_stop) const
   {
      const double QH1p = model.get_input().QH1p;
      const double QH2p = model.get_input().QH2p;
      const double QSp = model.get_input().QSp;
      const double QQp = model.get_input().QQp;
      const double Qup = model.get_input().Qup;

      const double g1 = model.get_g1();
      const double g2 = model.get_g2();
      const double gN = model.get_gN();
      const double vd = model.get_vd();
      const double vu = model.get_vu();
      const double vs = model.get_vs();
      const double yt = model.get_Yu(2,2);
      const double Lambdax = model.get_Lambdax();
      const double mq222 = model.get_mq2(2,2);
      const double mu222 = model.get_mu2(2,2);
      const double TYu22 = model.get_TYu(2,2);

      const double rt = stop_discriminant();

      double tmp_1 = 0.;
      double tmp_2 = -(1.4142135623730951 * vs * vu * yt * Lambdax) / Sqrt(rt);
      double tmp_3 = (vu * (TYu22 * vu - 0.7071067811865475 * vd * vs * yt * Lambdax) 
       * (2.0 * (0.25 * Sqr(g2) * vd - 0.25 * Sqr(g1) * vd + QH1p * QQp * Sqr(gN) * vd 
                 - QH1p * Qup * Sqr(gN) * vd) * (mq222 - mu222 + 0.125 * Sqr(g2) * Sqr(vd) - 0.125 * 
                                                 Sqr(g1) * Sqr(vd) - 0.125 * Sqr(g2) * Sqr(vu) + 0.125 *
                                                 Sqr(g1) * Sqr(vu) + 0.5 * QQp * Sqr(gN) *
                                                 (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
                                                 - 0.5 * Qup * Sqr(gN) * (QH1p * Sqr(vd) + QH2p * 
                                                                          Sqr(vu) + QSp * Sqr(vs)))
          - 2.8284271247461903 * vs * yt * Lambdax *
          (TYu22 * vu - 0.7071067811865475 * vd * vs * yt * Lambdax)))
      / (rt * Sqrt(rt));

      if (which_stop == stop_mass::mstop_1) {
         return 0.5 * (tmp_1 - (tmp_2 - tmp_3));
      } else {
         return 0.5 * (tmp_1 + (tmp_2 - tmp_3));
      }
   }

   double lowE6SSM_ew_derivs::deriv_d2MStop2_dvu_dvu(stop_mass which_stop) const
   {
      const double QH1p = model.get_input().QH1p;
      const double QH2p = model.get_input().QH2p;
      const double QSp = model.get_input().QSp;
      const double QQp = model.get_input().QQp;
      const double Qup = model.get_input().Qup;

      const double g1 = model.get_g1();
      const double g2 = model.get_g2();
      const double gN = model.get_gN();
      const double vd = model.get_vd();
      const double vu = model.get_vu();
      const double vs = model.get_vs();
      const double yt = model.get_Yu(2,2);
      const double Lambdax = model.get_Lambdax();
      const double mq222 = model.get_mq2(2,2);
      const double mu222 = model.get_mu2(2,2);
      const double TYu22 = model.get_TYu(2,2);

      const double rt = stop_discriminant();

      double tmp_1 = -0.25 * Sqr(g2) - 0.15 * Sqr(g1) + QH2p * QQp * Sqr(gN)
         + QH2p * Qup * Sqr(gN) + 2.0 * Sqr(yt);
      double tmp_2 = (4.0 * Sqr(TYu22) + 2.0 * Sqr(-0.25 * Sqr(g2) * vu + 0.25 * Sqr(g1) * vu + QH2p * QQp * Sqr(gN) * vu
                                - QH2p * Qup * Sqr(gN) * vu) + 
                      2.0 * (-0.25 * Sqr(g2) + 0.25 * Sqr(g1) + QH2p * QQp * Sqr(gN) - QH2p * Qup * Sqr(gN))
                      * (mq222 - mu222 + 0.125 * Sqr(g2) * Sqr(vd) - 0.125 * 
                         Sqr(g1) * Sqr(vd) - 0.125 * Sqr(g2) * Sqr(vu) + 0.125 *
                         Sqr(g1) * Sqr(vu) + 0.5 * QQp * Sqr(gN) *
                         (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
                         - 0.5 * Qup * Sqr(gN) * (QH1p * Sqr(vd) + QH2p * Sqr(vu)
                                                  + QSp * Sqr(vs))))/ (2.0 * Sqrt(rt));
      double tmp_3 = Sqr(2.0 * (-0.25 * Sqr(g2) * vu + 0.25 * Sqr(g1) * vu 
                                + QH2p * QQp * Sqr(gN) * vu - QH2p * Qup * Sqr(gN)
                                * vu) * (mq222 - mu222 + 0.125 * Sqr(g2) * Sqr(vd) - 0.125 * 
                                         Sqr(g1) * Sqr(vd) - 0.125 * Sqr(g2) * Sqr(vu) + 0.125 *
                                         Sqr(g1) * Sqr(vu) + 0.5 * QQp * Sqr(gN) *
                                         (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
                                         - 0.5 * Qup * Sqr(gN) * (QH1p * Sqr(vd) + QH2p * Sqr(vu)
                                                                  + QSp * Sqr(vs)))
                         + 4.0 * TYu22 * 
                         (TYu22 * vu - 0.7071067811865475 * vd * vs * yt * Lambdax)) 
         / (4.0 * rt * Sqrt(rt));

      if (which_stop == stop_mass::mstop_1) {
         return 0.5 * (tmp_1 - (tmp_2 - tmp_3));
      } else {
         return 0.5 * (tmp_1 + (tmp_2 - tmp_3));
      }
   }

   double lowE6SSM_ew_derivs::deriv_d2MStop2_dvu_dvs(stop_mass which_stop) const
   {
      const double QH1p = model.get_input().QH1p;
      const double QH2p = model.get_input().QH2p;
      const double QSp = model.get_input().QSp;
      const double QQp = model.get_input().QQp;
      const double Qup = model.get_input().Qup;

      const double g1 = model.get_g1();
      const double g2 = model.get_g2();
      const double gN = model.get_gN();
      const double vd = model.get_vd();
      const double vu = model.get_vu();
      const double vs = model.get_vs();
      const double yt = model.get_Yu(2,2);
      const double Lambdax = model.get_Lambdax();
      const double mq222 = model.get_mq2(2,2);
      const double mu222 = model.get_mu2(2,2);
      const double TYu22 = model.get_TYu(2,2);

      const double rt = stop_discriminant();

      double tmp_1 = 0.;
      double tmp_2 = (2.0 * (-0.25 * Sqr(g2) * vu + 0.25 * Sqr(g1) * vu + QH2p * QQp * Sqr(gN)
                             * vu - QH2p * Qup * Sqr(gN) * vu) * 
                      (QSp * QQp * Sqr(gN) * vs - QSp * Qup * Sqr(gN) * vs) 
                      - 2.8284271247461903 * TYu22 * vd * yt * Lambdax) / (2.0 * Sqrt(rt));
      double tmp_3 = ((2.0 * (QSp * QQp * Sqr(gN) * vs - QSp * Qup * Sqr(gN) * vs) *
                       (mq222 - mu222 + 0.125 * Sqr(g2) * Sqr(vd) - 0.125 * 
                        Sqr(g1) * Sqr(vd) - 0.125 * Sqr(g2) * Sqr(vu) + 0.125 *
                        Sqr(g1) * Sqr(vu) + 0.5 * QQp * Sqr(gN) *
                        (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
                        - 0.5 * Qup * Sqr(gN) * (QH1p * Sqr(vd) + QH2p * Sqr(vu)
                                                 + QSp * Sqr(vs)))
                       - 2.8284271247461903 * vd * yt * Lambdax *
                       (TYu22 * vu - 0.7071067811865475 * vd * vs * yt * Lambdax)) *
                      (2.0 * (-0.25 * Sqr(g2) * vu + 0.25 * Sqr(g1) * vu + QH2p * QQp * Sqr(gN) * vu
                              - QH2p * Qup * Sqr(gN) * vu) * 
                       (mq222 - mu222 + 0.125 * Sqr(g2) * Sqr(vd) - 0.125 * 
                        Sqr(g1) * Sqr(vd) - 0.125 * Sqr(g2) * Sqr(vu) + 0.125 *
                        Sqr(g1) * Sqr(vu) + 0.5 * QQp * Sqr(gN) *
                        (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
                        - 0.5 * Qup * Sqr(gN) * (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs)))
                       + 4.0 * TYu22 *
                       (TYu22 * vu - 0.7071067811865475 * vd * vs * yt * Lambdax)))
         / (4.0 * rt * Sqrt(rt));

      if (which_stop == stop_mass::mstop_1) {
         return 0.5 * (tmp_1 - (tmp_2 - tmp_3));
      } else {
         return 0.5 * (tmp_1 + (tmp_2 - tmp_3));
      }
   }

   double lowE6SSM_ew_derivs::deriv_d2MStop2_dg1_dvu(stop_mass which_stop) const
   {
      const double QH1p = model.get_input().QH1p;
      const double QH2p = model.get_input().QH2p;
      const double QSp = model.get_input().QSp;
      const double QQp = model.get_input().QQp;
      const double Qup = model.get_input().Qup;

      const double g1 = model.get_g1();
      const double g2 = model.get_g2();
      const double gN = model.get_gN();
      const double vd = model.get_vd();
      const double vu = model.get_vu();
      const double vs = model.get_vs();
      const double yt = model.get_Yu(2,2);
      const double Lambdax = model.get_Lambdax();
      const double mq222 = model.get_mq2(2,2);
      const double mu222 = model.get_mu2(2,2);
      const double TYu22 = model.get_TYu(2,2);

      const double rt = stop_discriminant();

      double tmp_1 = -0.3 * g1 * vu;
      double tmp_2 = (-0.5 * g1 * (-0.25 * Sqr(g2) * vu + 0.25 * Sqr(g1) * vu + QH2p * QQp * Sqr(gN)
                                   * vu - QH2p * Qup * Sqr(gN) * vu) * (Sqr(vd) - Sqr(vu)) + g1 * vu
                      * (mq222 - mu222 + 0.125 * Sqr(g2) * Sqr(vd) - 0.125 * 
                         Sqr(g1) * Sqr(vd) - 0.125 * Sqr(g2) * Sqr(vu) + 0.125 *
                         Sqr(g1) * Sqr(vu) + 0.5 * QQp * Sqr(gN) *
                         (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
                         - 0.5 * Qup * Sqr(gN) * (QH1p * Sqr(vd) + QH2p * Sqr(vu)
                                                  + QSp * Sqr(vs)))) / (2.0 * Sqrt(rt));
      double tmp_3 = -(0.125 * g1 * (Sqr(vd) - Sqr(vu)) * (mq222 - mu222 + 0.125 * Sqr(g2) * Sqr(vd) - 0.125 * 
                                                           Sqr(g1) * Sqr(vd) - 0.125 * Sqr(g2) * Sqr(vu) + 0.125 *
                                                           Sqr(g1) * Sqr(vu) + 0.5 * QQp * Sqr(gN) *
                                                           (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
                                                           - 0.5 * Qup * Sqr(gN) * (QH1p * Sqr(vd) + QH2p * Sqr(vu)
                                                                                    + QSp * Sqr(vs)))
                       * (2.0 * (-0.25 * Sqr(g2) * vu + 0.25 * Sqr(g1) * vu + QH2p * QQp * Sqr(gN) * vu 
                                 - QH2p * Qup * Sqr(gN) * vu) * (mq222 - mu222 + 0.125 * Sqr(g2) * Sqr(vd) - 0.125 * 
                                                                 Sqr(g1) * Sqr(vd) - 0.125 * Sqr(g2) * Sqr(vu) + 0.125 *
                                                                 Sqr(g1) * Sqr(vu) + 0.5 * QQp * Sqr(gN) *
                                                                 (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
                                                                 - 0.5 * Qup * Sqr(gN) * (QH1p * Sqr(vd) + QH2p * 
                                                                                          Sqr(vu) + QSp * Sqr(vs)))
                          + 4.0 * TYu22 *
                          (TYu22 * vu - 0.7071067811865475 * vd * vs * yt * Lambdax)))
         / (rt * Sqrt(rt));

      if (which_stop == stop_mass::mstop_1) {
         return 0.5 * (tmp_1 - (tmp_2 - tmp_3));
      } else {
         return 0.5 * (tmp_1 + (tmp_2 - tmp_3));
      }
   }

   double lowE6SSM_ew_derivs::deriv_d2MStop2_dg2_dvu(stop_mass which_stop) const
   {
      const double QH1p = model.get_input().QH1p;
      const double QH2p = model.get_input().QH2p;
      const double QSp = model.get_input().QSp;
      const double QQp = model.get_input().QQp;
      const double Qup = model.get_input().Qup;

      const double g1 = model.get_g1();
      const double g2 = model.get_g2();
      const double gN = model.get_gN();
      const double vd = model.get_vd();
      const double vu = model.get_vu();
      const double vs = model.get_vs();
      const double yt = model.get_Yu(2,2);
      const double Lambdax = model.get_Lambdax();
      const double mq222 = model.get_mq2(2,2);
      const double mu222 = model.get_mu2(2,2);
      const double TYu22 = model.get_TYu(2,2);

      const double rt = stop_discriminant();

      double tmp_1 = -0.5 * g2 * vu;
      double tmp_2 = (0.5 * g2 * (-0.25 * Sqr(g2) * vu + 0.25 * Sqr(g1) * vu + QH2p * QQp * Sqr(gN)
                                   * vu - QH2p * Qup * Sqr(gN) * vu) * (Sqr(vd) - Sqr(vu)) - g2 * vu
                      * (mq222 - mu222 + 0.125 * Sqr(g2) * Sqr(vd) - 0.125 * 
                         Sqr(g1) * Sqr(vd) - 0.125 * Sqr(g2) * Sqr(vu) + 0.125 *
                         Sqr(g1) * Sqr(vu) + 0.5 * QQp * Sqr(gN) *
                         (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
                         - 0.5 * Qup * Sqr(gN) * (QH1p * Sqr(vd) + QH2p * Sqr(vu)
                                                  + QSp * Sqr(vs)))) / (2.0 * Sqrt(rt));
      double tmp_3 = (0.125 * g2 * (Sqr(vd) - Sqr(vu)) * (mq222 - mu222 + 0.125 * Sqr(g2) * Sqr(vd) - 0.125 * 
                                                           Sqr(g1) * Sqr(vd) - 0.125 * Sqr(g2) * Sqr(vu) + 0.125 *
                                                           Sqr(g1) * Sqr(vu) + 0.5 * QQp * Sqr(gN) *
                                                           (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
                                                           - 0.5 * Qup * Sqr(gN) * (QH1p * Sqr(vd) + QH2p * Sqr(vu)
                                                                                    + QSp * Sqr(vs)))
                       * (2.0 * (-0.25 * Sqr(g2) * vu + 0.25 * Sqr(g1) * vu + QH2p * QQp * Sqr(gN) * vu 
                                 - QH2p * Qup * Sqr(gN) * vu) * (mq222 - mu222 + 0.125 * Sqr(g2) * Sqr(vd) - 0.125 * 
                                                                 Sqr(g1) * Sqr(vd) - 0.125 * Sqr(g2) * Sqr(vu) + 0.125 *
                                                                 Sqr(g1) * Sqr(vu) + 0.5 * QQp * Sqr(gN) *
                                                                 (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
                                                                 - 0.5 * Qup * Sqr(gN) * (QH1p * Sqr(vd) + QH2p * 
                                                                                          Sqr(vu) + QSp * Sqr(vs)))
                          + 4.0 * TYu22 *
                          (TYu22 * vu - 0.7071067811865475 * vd * vs * yt * Lambdax)))
         / (rt * Sqrt(rt));

      if (which_stop == stop_mass::mstop_1) {
         return 0.5 * (tmp_1 - (tmp_2 - tmp_3));
      } else {
         return 0.5 * (tmp_1 + (tmp_2 - tmp_3));
      }
   }

   double lowE6SSM_ew_derivs::deriv_d2MStop2_dgN_dvu(stop_mass which_stop) const
   {
      const double QH1p = model.get_input().QH1p;
      const double QH2p = model.get_input().QH2p;
      const double QSp = model.get_input().QSp;
      const double QQp = model.get_input().QQp;
      const double Qup = model.get_input().Qup;

      const double g1 = model.get_g1();
      const double g2 = model.get_g2();
      const double gN = model.get_gN();
      const double vd = model.get_vd();
      const double vu = model.get_vu();
      const double vs = model.get_vs();
      const double yt = model.get_Yu(2,2);
      const double Lambdax = model.get_Lambdax();
      const double mq222 = model.get_mq2(2,2);
      const double mu222 = model.get_mu2(2,2);
      const double TYu22 = model.get_TYu(2,2);

      const double rt = stop_discriminant();

      double tmp_1 = 2.0 * QH2p * QQp * gN * vu + 2.0 * QH2p * Qup * gN * vu;
      double tmp_2 = (2.0 * (-0.25 * Sqr(g2) * vu + 0.25 * Sqr(g1) * vu + QH2p * QQp * Sqr(gN)
                             * vu - QH2p * Qup * Sqr(gN) * vu) * 
                      (QQp * gN * (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
                       - Qup * gN * (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))) 
                      + 2.0 * (2.0 * QH1p * QQp * gN * vd - 2.0 * QH1p * Qup * gN * vd) *
                      (mq222 - mu222 + 0.125 * Sqr(g2) * Sqr(vd) - 0.125 * 
                       Sqr(g1) * Sqr(vd) - 0.125 * Sqr(g2) * Sqr(vu) + 0.125 *
                       Sqr(g1) * Sqr(vu) + 0.5 * QQp * Sqr(gN) *
                       (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
                       - 0.5 * Qup * Sqr(gN) * (QH1p * Sqr(vd) + QH2p * Sqr(vu)
                                                + QSp * Sqr(vs)))) / (2.0 * Sqrt(rt));
      double tmp_3 = ((QQp * gN * (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
                       - Qup * gN * (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))) * 
                      (mq222 - mu222 + 0.125 * Sqr(g2) * Sqr(vd) - 0.125 * 
                       Sqr(g1) * Sqr(vd) - 0.125 * Sqr(g2) * Sqr(vu) + 0.125 *
                       Sqr(g1) * Sqr(vu) + 0.5 * QQp * Sqr(gN) *
                       (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
                       - 0.5 * Qup * Sqr(gN) * (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs)))
                      * (2.0 * (-0.25 * Sqr(g2) * vu + 0.25 * Sqr(g1) * vu + QH2p * QQp * Sqr(gN) * vu 
                                - QH2p * Qup * Sqr(gN) * vu) * (mq222 - mu222 + 0.125 * Sqr(g2) * Sqr(vd) - 0.125 * 
                                                                Sqr(g1) * Sqr(vd) - 0.125 * Sqr(g2) * Sqr(vu) + 0.125 *
                                                                Sqr(g1) * Sqr(vu) + 0.5 * QQp * Sqr(gN) *
                                                                (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
                                                                - 0.5 * Qup * Sqr(gN) * (QH1p * Sqr(vd) + QH2p * 
                                                                                         Sqr(vu) + QSp * Sqr(vs)))
                         + 4.0 * TYu22 *
                         (TYu22 * vu - 0.7071067811865475 * vd * vs * yt * Lambdax)))
         / (2.0 * rt * Sqrt(rt));

      if (which_stop == stop_mass::mstop_1) {
         return 0.5 * (tmp_1 - (tmp_2 - tmp_3));
      } else {
         return 0.5 * (tmp_1 + (tmp_2 - tmp_3));
      }
   }

   double lowE6SSM_ew_derivs::deriv_d2MStop2_dYu22_dvu(stop_mass which_stop) const
   {
      const double QH1p = model.get_input().QH1p;
      const double QH2p = model.get_input().QH2p;
      const double QSp = model.get_input().QSp;
      const double QQp = model.get_input().QQp;
      const double Qup = model.get_input().Qup;

      const double g1 = model.get_g1();
      const double g2 = model.get_g2();
      const double gN = model.get_gN();
      const double vd = model.get_vd();
      const double vu = model.get_vu();
      const double vs = model.get_vs();
      const double yt = model.get_Yu(2,2);
      const double Lambdax = model.get_Lambdax();
      const double mq222 = model.get_mq2(2,2);
      const double mu222 = model.get_mu2(2,2);
      const double TYu22 = model.get_TYu(2,2);

      const double rt = stop_discriminant();

      double tmp_1 = 4.0 * yt * vu;
      double tmp_2 = (-1.4142135623730951 * TYu22 * vd * vs * Lambdax) / Sqrt(rt);
      double tmp_3 = -0.7071067811865475 * (vd * vs * Lambdax * 
                                            (TYu22 * vu - 0.7071067811865475 * vd * vs * yt * Lambdax)
                                            * (2.0 * (-0.25 * Sqr(g2) * vu + 0.25 * Sqr(g1) * vu + QH2p * QQp * Sqr(gN) * vu 
                                                      - QH2p * Qup * Sqr(gN) * vu) * 
                                               (mq222 - mu222 + 0.125 * Sqr(g2) * Sqr(vd) - 0.125 * 
                                                Sqr(g1) * Sqr(vd) - 0.125 * Sqr(g2) * Sqr(vu) + 0.125 *
                                                Sqr(g1) * Sqr(vu) + 0.5 * QQp * Sqr(gN) *
                                                (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
                                                - 0.5 * Qup * Sqr(gN) * (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs)))
                                               + 4.0 * TYu22 *
                                               (TYu22 * vu - 0.7071067811865475 * vd * vs * yt * Lambdax)))
         / (rt * Sqrt(rt));

      if (which_stop == stop_mass::mstop_1) {
         return 0.5 * (tmp_1 - (tmp_2 - tmp_3));
      } else {
         return 0.5 * (tmp_1 + (tmp_2 - tmp_3));
      }
   }

   double lowE6SSM_ew_derivs::deriv_d2MStop2_dmq222_dvu(stop_mass which_stop) const
   {
      const double QH1p = model.get_input().QH1p;
      const double QH2p = model.get_input().QH2p;
      const double QSp = model.get_input().QSp;
      const double QQp = model.get_input().QQp;
      const double Qup = model.get_input().Qup;

      const double g1 = model.get_g1();
      const double g2 = model.get_g2();
      const double gN = model.get_gN();
      const double vd = model.get_vd();
      const double vu = model.get_vu();
      const double vs = model.get_vs();
      const double yt = model.get_Yu(2,2);
      const double Lambdax = model.get_Lambdax();
      const double mq222 = model.get_mq2(2,2);
      const double mu222 = model.get_mu2(2,2);
      const double TYu22 = model.get_TYu(2,2);

      const double rt = stop_discriminant();

      double tmp_1 = 0.;
      double tmp_2 = (-0.25 * Sqr(g2) * vu + 0.25 * Sqr(g1) * vu + QH2p * QQp * Sqr(gN) * vu
                      - QH2p * Qup * Sqr(gN) * vu) / Sqrt(rt);
      double tmp_3 = ((mq222 - mu222 + 0.125 * Sqr(g2) * Sqr(vd) - 0.125 * 
                       Sqr(g1) * Sqr(vd) - 0.125 * Sqr(g2) * Sqr(vu) + 0.125 *
                       Sqr(g1) * Sqr(vu) + 0.5 * QQp * Sqr(gN) *
                       (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
                       - 0.5 * Qup * Sqr(gN) * (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs)))
                      * (2.0 * (-0.25 * Sqr(g2) * vu + 0.25 * Sqr(g1) * vu + QH2p * QQp * Sqr(gN) * vu 
                                - QH2p * Qup * Sqr(gN) * vu) * (mq222 - mu222 + 0.125 * Sqr(g2) * Sqr(vd) - 0.125 * 
                                                                Sqr(g1) * Sqr(vd) - 0.125 * Sqr(g2) * Sqr(vu) + 0.125 *
                                                                Sqr(g1) * Sqr(vu) + 0.5 * QQp * Sqr(gN) *
                                                                (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
                                                                - 0.5 * Qup * Sqr(gN) * (QH1p * Sqr(vd) + QH2p * 
                                                                                         Sqr(vu) + QSp * Sqr(vs)))
                         + 4.0 * TYu22 *
                         (TYu22 * vu - 0.7071067811865475 * vd * vs * yt * Lambdax)))
         / (2.0 * rt * Sqrt(rt));

      if (which_stop == stop_mass::mstop_1) {
         return 0.5 * (tmp_1 - (tmp_2 - tmp_3));
      } else {
         return 0.5 * (tmp_1 + (tmp_2 - tmp_3));
      }
   }

   double lowE6SSM_ew_derivs::deriv_d2MStop2_dmu222_dvu(stop_mass which_stop) const
   {
      const double QH1p = model.get_input().QH1p;
      const double QH2p = model.get_input().QH2p;
      const double QSp = model.get_input().QSp;
      const double QQp = model.get_input().QQp;
      const double Qup = model.get_input().Qup;

      const double g1 = model.get_g1();
      const double g2 = model.get_g2();
      const double gN = model.get_gN();
      const double vd = model.get_vd();
      const double vu = model.get_vu();
      const double vs = model.get_vs();
      const double yt = model.get_Yu(2,2);
      const double Lambdax = model.get_Lambdax();
      const double mq222 = model.get_mq2(2,2);
      const double mu222 = model.get_mu2(2,2);
      const double TYu22 = model.get_TYu(2,2);

      const double rt = stop_discriminant();

      double tmp_1 = 0.;
      double tmp_2 = -(-0.25 * Sqr(g2) * vu + 0.25 * Sqr(g1) * vu + QH2p * QQp * Sqr(gN) * vu
                       - QH2p * Qup * Sqr(gN) * vu) / Sqrt(rt);
      double tmp_3 = -((mq222 - mu222 + 0.125 * Sqr(g2) * Sqr(vd) - 0.125 * 
                        Sqr(g1) * Sqr(vd) - 0.125 * Sqr(g2) * Sqr(vu) + 0.125 *
                        Sqr(g1) * Sqr(vu) + 0.5 * QQp * Sqr(gN) *
                        (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
                        - 0.5 * Qup * Sqr(gN) * (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs)))
                       * (2.0 * (-0.25 * Sqr(g2) * vu + 0.25 * Sqr(g1) * vu + QH2p * QQp * Sqr(gN) * vu 
                                 - QH2p * Qup * Sqr(gN) * vu) * (mq222 - mu222 + 0.125 * Sqr(g2) * Sqr(vd) - 0.125 * 
                                                                 Sqr(g1) * Sqr(vd) - 0.125 * Sqr(g2) * Sqr(vu) + 0.125 *
                                                                 Sqr(g1) * Sqr(vu) + 0.5 * QQp * Sqr(gN) *
                                                                 (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
                                                                 - 0.5 * Qup * Sqr(gN) * (QH1p * Sqr(vd) + QH2p * 
                                                                                          Sqr(vu) + QSp * Sqr(vs)))
                          + 4.0 * TYu22 *
                          (TYu22 * vu - 0.7071067811865475 * vd * vs * yt * Lambdax)))
      / (2.0 * rt * Sqrt(rt));

      if (which_stop == stop_mass::mstop_1) {
         return 0.5 * (tmp_1 - (tmp_2 - tmp_3));
      } else {
         return 0.5 * (tmp_1 + (tmp_2 - tmp_3));
      }
   }

   double lowE6SSM_ew_derivs::deriv_d2MStop2_dLambdax_dvu(stop_mass which_stop) const
   {
      const double QH1p = model.get_input().QH1p;
      const double QH2p = model.get_input().QH2p;
      const double QSp = model.get_input().QSp;
      const double QQp = model.get_input().QQp;
      const double Qup = model.get_input().Qup;

      const double g1 = model.get_g1();
      const double g2 = model.get_g2();
      const double gN = model.get_gN();
      const double vd = model.get_vd();
      const double vu = model.get_vu();
      const double vs = model.get_vs();
      const double yt = model.get_Yu(2,2);
      const double Lambdax = model.get_Lambdax();
      const double mq222 = model.get_mq2(2,2);
      const double mu222 = model.get_mu2(2,2);
      const double TYu22 = model.get_TYu(2,2);

      const double rt = stop_discriminant();

      double tmp_1 = 0.;
      double tmp_2 = (-1.4142135623730951 * TYu22 * vd * vs * yt) / Sqrt(rt);
      double tmp_3 = -0.7071067811865475 * 
      (vd * vs * yt * (TYu22 * vu - 0.7071067811865475 * vd * vs * yt * Lambdax) 
       * (2.0 * (-0.25 * Sqr(g2) * vu + 0.25 * Sqr(g1) * vu + QH2p * QQp * Sqr(gN) * vu 
                 - QH2p * Qup * Sqr(gN) * vu) * (mq222 - mu222 + 0.125 * Sqr(g2) * Sqr(vd) - 0.125 * 
                                                 Sqr(g1) * Sqr(vd) - 0.125 * Sqr(g2) * Sqr(vu) + 0.125 *
                                                 Sqr(g1) * Sqr(vu) + 0.5 * QQp * Sqr(gN) *
                                                 (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
                                                 - 0.5 * Qup * Sqr(gN) * (QH1p * Sqr(vd) + QH2p * 
                                                                          Sqr(vu) + QSp * Sqr(vs)))
          + 4.0 * TYu22 *
          (TYu22 * vu - 0.7071067811865475 * vd * vs * yt * Lambdax)))
      / (rt * Sqrt(rt));

      if (which_stop == stop_mass::mstop_1) {
         return 0.5 * (tmp_1 - (tmp_2 - tmp_3));
      } else {
         return 0.5 * (tmp_1 + (tmp_2 - tmp_3));
      }
   }

   double lowE6SSM_ew_derivs::deriv_d2MStop2_dTYu22_dvu(stop_mass which_stop) const
   {
      const double QH1p = model.get_input().QH1p;
      const double QH2p = model.get_input().QH2p;
      const double QSp = model.get_input().QSp;
      const double QQp = model.get_input().QQp;
      const double Qup = model.get_input().Qup;

      const double g1 = model.get_g1();
      const double g2 = model.get_g2();
      const double gN = model.get_gN();
      const double vd = model.get_vd();
      const double vu = model.get_vu();
      const double vs = model.get_vs();
      const double yt = model.get_Yu(2,2);
      const double Lambdax = model.get_Lambdax();
      const double mq222 = model.get_mq2(2,2);
      const double mu222 = model.get_mu2(2,2);
      const double TYu22 = model.get_TYu(2,2);

      const double rt = stop_discriminant();

      double tmp_1 = 0.;
      double tmp_2 = (2.0 * TYu22 * vu + 2.0 * 
                      (TYu22 * vu - 0.7071067811865475 * vd * vs * yt * Lambdax)) / Sqrt(rt);
      double tmp_3 = (vu * (TYu22 * vu - 0.7071067811865475 * vd * vs * yt * Lambdax) 
       * (2.0 * (-0.25 * Sqr(g2) * vu + 0.25 * Sqr(g1) * vu + QH2p * QQp * Sqr(gN) * vu 
                 - QH2p * Qup * Sqr(gN) * vu) * (mq222 - mu222 + 0.125 * Sqr(g2) * Sqr(vd) - 0.125 * 
                                                 Sqr(g1) * Sqr(vd) - 0.125 * Sqr(g2) * Sqr(vu) + 0.125 *
                                                 Sqr(g1) * Sqr(vu) + 0.5 * QQp * Sqr(gN) *
                                                 (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
                                                 - 0.5 * Qup * Sqr(gN) * (QH1p * Sqr(vd) + QH2p * 
                                                                          Sqr(vu) + QSp * Sqr(vs)))
          + 4.0 * TYu22 *
          (TYu22 * vu - 0.7071067811865475 * vd * vs * yt * Lambdax)))
      / (rt * Sqrt(rt));

      if (which_stop == stop_mass::mstop_1) {
         return 0.5 * (tmp_1 - (tmp_2 - tmp_3));
      } else {
         return 0.5 * (tmp_1 + (tmp_2 - tmp_3));
      }
   }

   double lowE6SSM_ew_derivs::deriv_d2MStop2_dvs_dvs(stop_mass which_stop) const
   {
      const double QH1p = model.get_input().QH1p;
      const double QH2p = model.get_input().QH2p;
      const double QSp = model.get_input().QSp;
      const double QQp = model.get_input().QQp;
      const double Qup = model.get_input().Qup;

      const double g1 = model.get_g1();
      const double g2 = model.get_g2();
      const double gN = model.get_gN();
      const double vd = model.get_vd();
      const double vu = model.get_vu();
      const double vs = model.get_vs();
      const double yt = model.get_Yu(2,2);
      const double Lambdax = model.get_Lambdax();
      const double mq222 = model.get_mq2(2,2);
      const double mu222 = model.get_mu2(2,2);
      const double TYu22 = model.get_TYu(2,2);

      const double rt = stop_discriminant();

      double tmp_1 = QSp * QQp * Sqr(gN) + QSp * Qup * Sqr(gN);
      double tmp_2 = (2.0 * Sqr(vd) * Sqr(yt) * Sqr(Lambdax) 
                      + 2.0 * Sqr(QSp * QQp * Sqr(gN) * vs - QSp * Qup * Sqr(gN) * vs) + 
                      2.0 * (QSp * QQp * Sqr(gN) - QSp * Qup * Sqr(gN))
                      * (mq222 - mu222 + 0.125 * Sqr(g2) * Sqr(vd) - 0.125 * 
                         Sqr(g1) * Sqr(vd) - 0.125 * Sqr(g2) * Sqr(vu) + 0.125 *
                         Sqr(g1) * Sqr(vu) + 0.5 * QQp * Sqr(gN) *
                         (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
                         - 0.5 * Qup * Sqr(gN) * (QH1p * Sqr(vd) + QH2p * Sqr(vu)
                                                  + QSp * Sqr(vs))))/ (2.0 * Sqrt(rt));
      double tmp_3 = Sqr(2.0 * (QSp * QQp * Sqr(gN) * vs - QSp * Qup * Sqr(gN) * vs) 
                         * (mq222 - mu222 + 0.125 * Sqr(g2) * Sqr(vd) - 0.125 * 
                            Sqr(g1) * Sqr(vd) - 0.125 * Sqr(g2) * Sqr(vu) + 0.125 *
                            Sqr(g1) * Sqr(vu) + 0.5 * QQp * Sqr(gN) *
                            (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
                            - 0.5 * Qup * Sqr(gN) * (QH1p * Sqr(vd) + QH2p * Sqr(vu)
                                                     + QSp * Sqr(vs)))
                         - 2.8284271247461903 * vd * yt * Lambdax * 
                         (TYu22 * vu - 0.7071067811865475 * vd * vs * yt * Lambdax)) 
         / (4.0 * rt * Sqrt(rt));

      if (which_stop == stop_mass::mstop_1) {
         return 0.5 * (tmp_1 - (tmp_2 - tmp_3));
      } else {
         return 0.5 * (tmp_1 + (tmp_2 - tmp_3));
      }
   }

   double lowE6SSM_ew_derivs::deriv_d2MStop2_dg1_dvs(stop_mass which_stop) const
   {
      const double QH1p = model.get_input().QH1p;
      const double QH2p = model.get_input().QH2p;
      const double QSp = model.get_input().QSp;
      const double QQp = model.get_input().QQp;
      const double Qup = model.get_input().Qup;

      const double g1 = model.get_g1();
      const double g2 = model.get_g2();
      const double gN = model.get_gN();
      const double vd = model.get_vd();
      const double vu = model.get_vu();
      const double vs = model.get_vs();
      const double yt = model.get_Yu(2,2);
      const double Lambdax = model.get_Lambdax();
      const double mq222 = model.get_mq2(2,2);
      const double mu222 = model.get_mu2(2,2);
      const double TYu22 = model.get_TYu(2,2);

      const double rt = stop_discriminant();

      double tmp_1 = 0;
      double tmp_2 = (-0.25 * g1 * (QSp * QQp * Sqr(gN) * vs - QSp * Qup * Sqr(gN) * vs) 
                       * (Sqr(vd) - Sqr(vu))) / Sqrt(rt);
      double tmp_3 = -(0.125 * g1 * (Sqr(vd) - Sqr(vu)) * (mq222 - mu222 + 0.125 * Sqr(g2) * Sqr(vd) - 0.125 * 
                                                           Sqr(g1) * Sqr(vd) - 0.125 * Sqr(g2) * Sqr(vu) + 0.125 *
                                                           Sqr(g1) * Sqr(vu) + 0.5 * QQp * Sqr(gN) *
                                                           (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
                                                           - 0.5 * Qup * Sqr(gN) * (QH1p * Sqr(vd) + QH2p * Sqr(vu)
                                                                                    + QSp * Sqr(vs)))
                       * (2.0 * (QSp * QQp * Sqr(gN) * vs - QSp * Qup * Sqr(gN) * vs) 
                          * (mq222 - mu222 + 0.125 * Sqr(g2) * Sqr(vd) - 0.125 * 
                             Sqr(g1) * Sqr(vd) - 0.125 * Sqr(g2) * Sqr(vu) + 0.125 *
                             Sqr(g1) * Sqr(vu) + 0.5 * QQp * Sqr(gN) *
                             (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
                             - 0.5 * Qup * Sqr(gN) * (QH1p * Sqr(vd) + QH2p * 
                                                      Sqr(vu) + QSp * Sqr(vs)))
                          - 2.8284271247461903 * vd * yt * Lambdax *
                          (TYu22 * vu - 0.7071067811865475 * vd * vs * yt * Lambdax)))
         / (rt * Sqrt(rt));

      if (which_stop == stop_mass::mstop_1) {
         return 0.5 * (tmp_1 - (tmp_2 - tmp_3));
      } else {
         return 0.5 * (tmp_1 + (tmp_2 - tmp_3));
      }
   }

   double lowE6SSM_ew_derivs::deriv_d2MStop2_dg2_dvs(stop_mass which_stop) const
   {
      const double QH1p = model.get_input().QH1p;
      const double QH2p = model.get_input().QH2p;
      const double QSp = model.get_input().QSp;
      const double QQp = model.get_input().QQp;
      const double Qup = model.get_input().Qup;

      const double g1 = model.get_g1();
      const double g2 = model.get_g2();
      const double gN = model.get_gN();
      const double vd = model.get_vd();
      const double vu = model.get_vu();
      const double vs = model.get_vs();
      const double yt = model.get_Yu(2,2);
      const double Lambdax = model.get_Lambdax();
      const double mq222 = model.get_mq2(2,2);
      const double mu222 = model.get_mu2(2,2);
      const double TYu22 = model.get_TYu(2,2);

      const double rt = stop_discriminant();

      double tmp_1 = 0.;
      double tmp_2 = (0.25 * g2 * (QSp * QQp * Sqr(gN) * vs - QSp * Qup * Sqr(gN) * vs) 
                      * (Sqr(vd) - Sqr(vu))) / Sqrt(rt);
      double tmp_3 = (0.125 * g2 * (Sqr(vd) - Sqr(vu)) * (mq222 - mu222 + 0.125 * Sqr(g2) * Sqr(vd) - 0.125 * 
                                                           Sqr(g1) * Sqr(vd) - 0.125 * Sqr(g2) * Sqr(vu) + 0.125 *
                                                           Sqr(g1) * Sqr(vu) + 0.5 * QQp * Sqr(gN) *
                                                           (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
                                                           - 0.5 * Qup * Sqr(gN) * (QH1p * Sqr(vd) + QH2p * Sqr(vu)
                                                                                    + QSp * Sqr(vs)))
                       * (2.0 * (QSp * QQp * Sqr(gN) * vs - QSp * Qup * Sqr(gN) * vs) 
                          * (mq222 - mu222 + 0.125 * Sqr(g2) * Sqr(vd) - 0.125 * 
                             Sqr(g1) * Sqr(vd) - 0.125 * Sqr(g2) * Sqr(vu) + 0.125 *
                             Sqr(g1) * Sqr(vu) + 0.5 * QQp * Sqr(gN) *
                             (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
                             - 0.5 * Qup * Sqr(gN) * (QH1p * Sqr(vd) + QH2p * 
                                                      Sqr(vu) + QSp * Sqr(vs)))
                          - 2.8284271247461903 * vd * yt * Lambdax *
                          (TYu22 * vu - 0.7071067811865475 * vd * vs * yt * Lambdax)))
         / (rt * Sqrt(rt));

      if (which_stop == stop_mass::mstop_1) {
         return 0.5 * (tmp_1 - (tmp_2 - tmp_3));
      } else {
         return 0.5 * (tmp_1 + (tmp_2 - tmp_3));
      }
   }

   double lowE6SSM_ew_derivs::deriv_d2MStop2_dgN_dvs(stop_mass which_stop) const
   {
      const double QH1p = model.get_input().QH1p;
      const double QH2p = model.get_input().QH2p;
      const double QSp = model.get_input().QSp;
      const double QQp = model.get_input().QQp;
      const double Qup = model.get_input().Qup;

      const double g1 = model.get_g1();
      const double g2 = model.get_g2();
      const double gN = model.get_gN();
      const double vd = model.get_vd();
      const double vu = model.get_vu();
      const double vs = model.get_vs();
      const double yt = model.get_Yu(2,2);
      const double Lambdax = model.get_Lambdax();
      const double mq222 = model.get_mq2(2,2);
      const double mu222 = model.get_mu2(2,2);
      const double TYu22 = model.get_TYu(2,2);

      const double rt = stop_discriminant();

      double tmp_1 = 2.0 * QSp * QQp * gN * vs + 2.0 * QSp * Qup * gN * vs;
      double tmp_2 = ((QSp * QQp * Sqr(gN) * vs - QSp * Qup * Sqr(gN) * vs) * 
                      (QQp * gN * (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
                       - Qup * gN * (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))) 
                      + 2.0 * (2.0 * QSp * QQp * gN * vs - 2.0 * QSp * Qup * gN * vs) *
                      (mq222 - mu222 + 0.125 * Sqr(g2) * Sqr(vd) - 0.125 * 
                       Sqr(g1) * Sqr(vd) - 0.125 * Sqr(g2) * Sqr(vu) + 0.125 *
                       Sqr(g1) * Sqr(vu) + 0.5 * QQp * Sqr(gN) *
                       (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
                       - 0.5 * Qup * Sqr(gN) * (QH1p * Sqr(vd) + QH2p * Sqr(vu)
                                                + QSp * Sqr(vs)))) / (2.0 * Sqrt(rt));
      double tmp_3 = ((QQp * gN * (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
                       - Qup * gN * (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))) * 
                      (mq222 - mu222 + 0.125 * Sqr(g2) * Sqr(vd) - 0.125 * 
                       Sqr(g1) * Sqr(vd) - 0.125 * Sqr(g2) * Sqr(vu) + 0.125 *
                       Sqr(g1) * Sqr(vu) + 0.5 * QQp * Sqr(gN) *
                       (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
                       - 0.5 * Qup * Sqr(gN) * (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs)))
                      * (2.0 * (QSp * QQp * Sqr(gN) * vs - QSp * Qup * Sqr(gN) * vs) 
                         * (mq222 - mu222 + 0.125 * Sqr(g2) * Sqr(vd) - 0.125 * 
                            Sqr(g1) * Sqr(vd) - 0.125 * Sqr(g2) * Sqr(vu) + 0.125 *
                            Sqr(g1) * Sqr(vu) + 0.5 * QQp * Sqr(gN) *
                            (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
                            - 0.5 * Qup * Sqr(gN) * (QH1p * Sqr(vd) + QH2p * 
                                                     Sqr(vu) + QSp * Sqr(vs)))
                         - 2.8284271247461903 * vd * yt * Lambdax *
                         (TYu22 * vu - 0.7071067811865475 * vd * vs * yt * Lambdax)))
         / (2.0 * rt * Sqrt(rt));

      if (which_stop == stop_mass::mstop_1) {
         return 0.5 * (tmp_1 - (tmp_2 - tmp_3));
      } else {
         return 0.5 * (tmp_1 + (tmp_2 - tmp_3));
      }
   }

   double lowE6SSM_ew_derivs::deriv_d2MStop2_dYu22_dvs(stop_mass which_stop) const
   {
      const double QH1p = model.get_input().QH1p;
      const double QH2p = model.get_input().QH2p;
      const double QSp = model.get_input().QSp;
      const double QQp = model.get_input().QQp;
      const double Qup = model.get_input().Qup;

      const double g1 = model.get_g1();
      const double g2 = model.get_g2();
      const double gN = model.get_gN();
      const double vd = model.get_vd();
      const double vu = model.get_vu();
      const double vs = model.get_vs();
      const double yt = model.get_Yu(2,2);
      const double Lambdax = model.get_Lambdax();
      const double mq222 = model.get_mq2(2,2);
      const double mu222 = model.get_mu2(2,2);
      const double TYu22 = model.get_TYu(2,2);

      const double rt = stop_discriminant();

      double tmp_1 = 0.;
      double tmp_2 = (2.0 * vs * Sqr(vd) * yt * Sqr(Lambdax) - 
                      2.8284271247461903 * vd * Lambdax * 
                      (TYu22 * vu - 0.7071067811865475 * vd * vs * yt * Lambdax)) / (2.0 * Sqrt(rt));
      double tmp_3 = -0.7071067811865475 * (vd * vs * Lambdax * 
                                            (TYu22 * vu - 0.7071067811865475 * vd * vs * yt * Lambdax)
                                            * (2.0 * (QSp * QQp * Sqr(gN) * vs - QSp * Qup * Sqr(gN) * vs) * 
                                               (mq222 - mu222 + 0.125 * Sqr(g2) * Sqr(vd) - 0.125 * 
                                                Sqr(g1) * Sqr(vd) - 0.125 * Sqr(g2) * Sqr(vu) + 0.125 *
                                                Sqr(g1) * Sqr(vu) + 0.5 * QQp * Sqr(gN) *
                                                (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
                                                - 0.5 * Qup * Sqr(gN) * (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs)))
                                               - 2.8284271247461903 * vd * yt * Lambdax *
                                               (TYu22 * vu - 0.7071067811865475 * vd * vs * yt * Lambdax)))
         / (rt * Sqrt(rt));

      if (which_stop == stop_mass::mstop_1) {
         return 0.5 * (tmp_1 - (tmp_2 - tmp_3));
      } else {
         return 0.5 * (tmp_1 + (tmp_2 - tmp_3));
      }
   }

   double lowE6SSM_ew_derivs::deriv_d2MStop2_dmq222_dvs(stop_mass which_stop) const
   {
      const double QH1p = model.get_input().QH1p;
      const double QH2p = model.get_input().QH2p;
      const double QSp = model.get_input().QSp;
      const double QQp = model.get_input().QQp;
      const double Qup = model.get_input().Qup;

      const double g1 = model.get_g1();
      const double g2 = model.get_g2();
      const double gN = model.get_gN();
      const double vd = model.get_vd();
      const double vu = model.get_vu();
      const double vs = model.get_vs();
      const double yt = model.get_Yu(2,2);
      const double Lambdax = model.get_Lambdax();
      const double mq222 = model.get_mq2(2,2);
      const double mu222 = model.get_mu2(2,2);
      const double TYu22 = model.get_TYu(2,2);

      const double rt = stop_discriminant();

      double tmp_1 = 0.;
      double tmp_2 = (QSp * QQp * Sqr(gN) * vs - QSp * Qup * Sqr(gN) * vs) / Sqrt(rt);
      double tmp_3 = ((mq222 - mu222 + 0.125 * Sqr(g2) * Sqr(vd) - 0.125 * 
                       Sqr(g1) * Sqr(vd) - 0.125 * Sqr(g2) * Sqr(vu) + 0.125 *
                       Sqr(g1) * Sqr(vu) + 0.5 * QQp * Sqr(gN) *
                       (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
                       - 0.5 * Qup * Sqr(gN) * (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs)))
                      * (2.0 * (QSp * QQp * Sqr(gN) * vs - QSp * Qup * Sqr(gN) * vs) 
                         * (mq222 - mu222 + 0.125 * Sqr(g2) * Sqr(vd) - 0.125 * 
                            Sqr(g1) * Sqr(vd) - 0.125 * Sqr(g2) * Sqr(vu) + 0.125 *
                            Sqr(g1) * Sqr(vu) + 0.5 * QQp * Sqr(gN) *
                            (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
                            - 0.5 * Qup * Sqr(gN) * (QH1p * Sqr(vd) + QH2p * 
                                                     Sqr(vu) + QSp * Sqr(vs)))
                         - 2.8284271247461903 * vd * yt * Lambdax *
                         (TYu22 * vu - 0.7071067811865475 * vd * vs * yt * Lambdax)))
         / (2.0 * rt * Sqrt(rt));

      if (which_stop == stop_mass::mstop_1) {
         return 0.5 * (tmp_1 - (tmp_2 - tmp_3));
      } else {
         return 0.5 * (tmp_1 + (tmp_2 - tmp_3));
      }
   }

   double lowE6SSM_ew_derivs::deriv_d2MStop2_dmu222_dvs(stop_mass which_stop) const
   {
      const double QH1p = model.get_input().QH1p;
      const double QH2p = model.get_input().QH2p;
      const double QSp = model.get_input().QSp;
      const double QQp = model.get_input().QQp;
      const double Qup = model.get_input().Qup;

      const double g1 = model.get_g1();
      const double g2 = model.get_g2();
      const double gN = model.get_gN();
      const double vd = model.get_vd();
      const double vu = model.get_vu();
      const double vs = model.get_vs();
      const double yt = model.get_Yu(2,2);
      const double Lambdax = model.get_Lambdax();
      const double mq222 = model.get_mq2(2,2);
      const double mu222 = model.get_mu2(2,2);
      const double TYu22 = model.get_TYu(2,2);

      const double rt = stop_discriminant();

      double tmp_1 = 0.;
      double tmp_2 = -(QSp * QQp * Sqr(gN) * vs - QSp * Qup * Sqr(gN) * vs) / Sqrt(rt);
      double tmp_3 = -((mq222 - mu222 + 0.125 * Sqr(g2) * Sqr(vd) - 0.125 * 
                        Sqr(g1) * Sqr(vd) - 0.125 * Sqr(g2) * Sqr(vu) + 0.125 *
                        Sqr(g1) * Sqr(vu) + 0.5 * QQp * Sqr(gN) *
                        (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
                        - 0.5 * Qup * Sqr(gN) * (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs)))
                       * (2.0 * (QSp * QQp * Sqr(gN) * vs - QSp * Qup * Sqr(gN) * vs) 
                          * (mq222 - mu222 + 0.125 * Sqr(g2) * Sqr(vd) - 0.125 * 
                             Sqr(g1) * Sqr(vd) - 0.125 * Sqr(g2) * Sqr(vu) + 0.125 *
                             Sqr(g1) * Sqr(vu) + 0.5 * QQp * Sqr(gN) *
                             (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
                             - 0.5 * Qup * Sqr(gN) * (QH1p * Sqr(vd) + QH2p * 
                                                      Sqr(vu) + QSp * Sqr(vs)))
                          - 2.8284271247461903 * vd * yt * Lambdax *
                          (TYu22 * vu - 0.7071067811865475 * vd * vs * yt * Lambdax)))
      / (2.0 * rt * Sqrt(rt));

      if (which_stop == stop_mass::mstop_1) {
         return 0.5 * (tmp_1 - (tmp_2 - tmp_3));
      } else {
         return 0.5 * (tmp_1 + (tmp_2 - tmp_3));
      }
   }

   double lowE6SSM_ew_derivs::deriv_d2MStop2_dLambdax_dvs(stop_mass which_stop) const
   {
      const double QH1p = model.get_input().QH1p;
      const double QH2p = model.get_input().QH2p;
      const double QSp = model.get_input().QSp;
      const double QQp = model.get_input().QQp;
      const double Qup = model.get_input().Qup;

      const double g1 = model.get_g1();
      const double g2 = model.get_g2();
      const double gN = model.get_gN();
      const double vd = model.get_vd();
      const double vu = model.get_vu();
      const double vs = model.get_vs();
      const double yt = model.get_Yu(2,2);
      const double Lambdax = model.get_Lambdax();
      const double mq222 = model.get_mq2(2,2);
      const double mu222 = model.get_mu2(2,2);
      const double TYu22 = model.get_TYu(2,2);

      const double rt = stop_discriminant();

      double tmp_1 = 0.;
      double tmp_2 = (2.0 * vs * Sqr(vd) * Sqr(yt) * Lambdax 
                      - 2.8284271247461903 * vd * yt *
                      (TYu22 * vu - 0.7071067811865475 * vd * vs * yt * Lambdax)) / (2.0 * Sqrt(rt));
      double tmp_3 = -0.7071067811865475 * 
      (vd * vs * yt * (TYu22 * vu - 0.7071067811865475 * vd * vs * yt * Lambdax) 
       * (2.0 * (QSp * QQp * Sqr(gN) * vs - QSp * Qup * Sqr(gN) * vs) 
          * (mq222 - mu222 + 0.125 * Sqr(g2) * Sqr(vd) - 0.125 * 
             Sqr(g1) * Sqr(vd) - 0.125 * Sqr(g2) * Sqr(vu) + 0.125 *
             Sqr(g1) * Sqr(vu) + 0.5 * QQp * Sqr(gN) *
             (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
             - 0.5 * Qup * Sqr(gN) * (QH1p * Sqr(vd) + QH2p * 
                                      Sqr(vu) + QSp * Sqr(vs)))
          - 2.8284271247461903 * vd * yt * Lambdax *
          (TYu22 * vu - 0.7071067811865475 * vd * vs * yt * Lambdax)))
      / (rt * Sqrt(rt));

      if (which_stop == stop_mass::mstop_1) {
         return 0.5 * (tmp_1 - (tmp_2 - tmp_3));
      } else {
         return 0.5 * (tmp_1 + (tmp_2 - tmp_3));
      }
   }

   double lowE6SSM_ew_derivs::deriv_d2MStop2_dTYu22_dvs(stop_mass which_stop) const
   {
      const double QH1p = model.get_input().QH1p;
      const double QH2p = model.get_input().QH2p;
      const double QSp = model.get_input().QSp;
      const double QQp = model.get_input().QQp;
      const double Qup = model.get_input().Qup;

      const double g1 = model.get_g1();
      const double g2 = model.get_g2();
      const double gN = model.get_gN();
      const double vd = model.get_vd();
      const double vu = model.get_vu();
      const double vs = model.get_vs();
      const double yt = model.get_Yu(2,2);
      const double Lambdax = model.get_Lambdax();
      const double mq222 = model.get_mq2(2,2);
      const double mu222 = model.get_mu2(2,2);
      const double TYu22 = model.get_TYu(2,2);

      const double rt = stop_discriminant();

      double tmp_1 = 0.;
      double tmp_2 = -(1.4142135623730951 * vd * vu * yt * Lambdax) / Sqrt(rt);
      double tmp_3 = (vu * (TYu22 * vu - 0.7071067811865475 * vd * vs * yt * Lambdax) 
       * (2.0 * (QSp * QQp * Sqr(gN) * vs - QSp * Qup * Sqr(gN) * vs) 
          * (mq222 - mu222 + 0.125 * Sqr(g2) * Sqr(vd) - 0.125 * 
             Sqr(g1) * Sqr(vd) - 0.125 * Sqr(g2) * Sqr(vu) + 0.125 *
             Sqr(g1) * Sqr(vu) + 0.5 * QQp * Sqr(gN) *
             (QH1p * Sqr(vd) + QH2p * Sqr(vu) + QSp * Sqr(vs))
             - 0.5 * Qup * Sqr(gN) * (QH1p * Sqr(vd) + QH2p * 
                                      Sqr(vu) + QSp * Sqr(vs)))
          - 2.8284271247461903 * vd * yt * Lambdax *
          (TYu22 * vu - 0.7071067811865475 * vd * vs * yt * Lambdax)))
      / (rt * Sqrt(rt));

      if (which_stop == stop_mass::mstop_1) {
         return 0.5 * (tmp_1 - (tmp_2 - tmp_3));
      } else {
         return 0.5 * (tmp_1 + (tmp_2 - tmp_3));
      }
   }

   double lowE6SSM_ew_derivs::deriv_d2MStop2_dg1_dg1(stop_mass which_stop) const
   {
      double result = (which_stop == stop_mass::mstop_1 ? 1.0 : 2.0);
      return result;
   }

   double lowE6SSM_ew_derivs::deriv_d2MStop2_dg2_dg1(stop_mass which_stop) const
   {
      double result = (which_stop == stop_mass::mstop_1 ? 1.0 : 2.0);
      return result;
   }

   double lowE6SSM_ew_derivs::deriv_d2MStop2_dgN_dg1(stop_mass which_stop) const
   {
      double result = (which_stop == stop_mass::mstop_1 ? 1.0 : 2.0);
      return result;
   }

   double lowE6SSM_ew_derivs::deriv_d2MStop2_dYu22_dg1(stop_mass which_stop) const
   {
      double result = (which_stop == stop_mass::mstop_1 ? 1.0 : 2.0);
      return result;
   }

   double lowE6SSM_ew_derivs::deriv_d2MStop2_dmq222_dg1(stop_mass which_stop) const
   {
      double result = (which_stop == stop_mass::mstop_1 ? 1.0 : 2.0);
      return result;
   }

   double lowE6SSM_ew_derivs::deriv_d2MStop2_dmu222_dg1(stop_mass which_stop) const
   {
      double result = (which_stop == stop_mass::mstop_1 ? 1.0 : 2.0);
      return result;
   }

   double lowE6SSM_ew_derivs::deriv_d2MStop2_dLambdax_dg1(stop_mass which_stop) const
   {
      double result = (which_stop == stop_mass::mstop_1 ? 1.0 : 2.0);
      return result;
   }

   double lowE6SSM_ew_derivs::deriv_d2MStop2_dTYu22_dg1(stop_mass which_stop) const
   {
      double result = (which_stop == stop_mass::mstop_1 ? 1.0 : 2.0);
      return result;
   }

   double lowE6SSM_ew_derivs::deriv_d2MStop2_dg2_dg2(stop_mass which_stop) const
   {
      double result = (which_stop == stop_mass::mstop_1 ? 1.0 : 2.0);
      return result;
   }

   double lowE6SSM_ew_derivs::deriv_d2MStop2_dgN_dg2(stop_mass which_stop) const
   {
      double result = (which_stop == stop_mass::mstop_1 ? 1.0 : 2.0);
      return result;
   }

   double lowE6SSM_ew_derivs::deriv_d2MStop2_dYu22_dg2(stop_mass which_stop) const
   {
      double result = (which_stop == stop_mass::mstop_1 ? 1.0 : 2.0);
      return result;
   }

   double lowE6SSM_ew_derivs::deriv_d2MStop2_dmq222_dg2(stop_mass which_stop) const
   {
      double result = (which_stop == stop_mass::mstop_1 ? 1.0 : 2.0);
      return result;
   }

   double lowE6SSM_ew_derivs::deriv_d2MStop2_dmu222_dg2(stop_mass which_stop) const
   {
      double result = (which_stop == stop_mass::mstop_1 ? 1.0 : 2.0);
      return result;
   }

   double lowE6SSM_ew_derivs::deriv_d2MStop2_dLambdax_dg2(stop_mass which_stop) const
   {
      double result = (which_stop == stop_mass::mstop_1 ? 1.0 : 2.0);
      return result;
   }

   double lowE6SSM_ew_derivs::deriv_d2MStop2_dTYu22_dg2(stop_mass which_stop) const
   {
      double result = (which_stop == stop_mass::mstop_1 ? 1.0 : 2.0);
      return result;
   }

   double lowE6SSM_ew_derivs::deriv_d2MStop2_dgN_dgN(stop_mass which_stop) const
   {
      double result = (which_stop == stop_mass::mstop_1 ? 1.0 : 2.0);
      return result;
   }

   double lowE6SSM_ew_derivs::deriv_d2MStop2_dYu22_dgN(stop_mass which_stop) const
   {
      double result = (which_stop == stop_mass::mstop_1 ? 1.0 : 2.0);
      return result;
   }

   double lowE6SSM_ew_derivs::deriv_d2MStop2_dmq222_dgN(stop_mass which_stop) const
   {
      double result = (which_stop == stop_mass::mstop_1 ? 1.0 : 2.0);
      return result;
   }

   double lowE6SSM_ew_derivs::deriv_d2MStop2_dmu222_dgN(stop_mass which_stop) const
   {
      double result = (which_stop == stop_mass::mstop_1 ? 1.0 : 2.0);
      return result;
   }

   double lowE6SSM_ew_derivs::deriv_d2MStop2_dLambdax_dgN(stop_mass which_stop) const
   {
      double result = (which_stop == stop_mass::mstop_1 ? 1.0 : 2.0);
      return result;
   }

   double lowE6SSM_ew_derivs::deriv_d2MStop2_dTYu22_dgN(stop_mass which_stop) const
   {
      double result = (which_stop == stop_mass::mstop_1 ? 1.0 : 2.0);
      return result;
   }

   double lowE6SSM_ew_derivs::deriv_d2MStop2_dYu22_dYu22(stop_mass which_stop) const
   {
      double result = (which_stop == stop_mass::mstop_1 ? 1.0 : 2.0);
      return result;
   }

   double lowE6SSM_ew_derivs::deriv_d2MStop2_dmq222_dYu22(stop_mass which_stop) const
   {
      double result = (which_stop == stop_mass::mstop_1 ? 1.0 : 2.0);
      return result;
   }

   double lowE6SSM_ew_derivs::deriv_d2MStop2_dmu222_dYu22(stop_mass which_stop) const
   {
      double result = (which_stop == stop_mass::mstop_1 ? 1.0 : 2.0);
      return result;
   }

   double lowE6SSM_ew_derivs::deriv_d2MStop2_dLambdax_dYu22(stop_mass which_stop) const
   {
      double result = (which_stop == stop_mass::mstop_1 ? 1.0 : 2.0);
      return result;
   }

   double lowE6SSM_ew_derivs::deriv_d2MStop2_dTYu22_dYu22(stop_mass which_stop) const
   {
      double result = (which_stop == stop_mass::mstop_1 ? 1.0 : 2.0);
      return result;
   }

   double lowE6SSM_ew_derivs::deriv_d2MStop2_dmq222_dmq222(stop_mass which_stop) const
   {
      double result = (which_stop == stop_mass::mstop_1 ? 1.0 : 2.0);
      return result;
   }

   double lowE6SSM_ew_derivs::deriv_d2MStop2_dmu222_dmq222(stop_mass which_stop) const
   {
      double result = (which_stop == stop_mass::mstop_1 ? 1.0 : 2.0);
      return result;
   }

   double lowE6SSM_ew_derivs::deriv_d2MStop2_dLambdax_dmq222(stop_mass which_stop) const
   {
      double result = (which_stop == stop_mass::mstop_1 ? 1.0 : 2.0);
      return result;
   }

   double lowE6SSM_ew_derivs::deriv_d2MStop2_dTYu22_dmq222(stop_mass which_stop) const
   {
      double result = (which_stop == stop_mass::mstop_1 ? 1.0 : 2.0);
      return result;
   }

   double lowE6SSM_ew_derivs::deriv_d2MStop2_dmu222_dmu222(stop_mass which_stop) const
   {
      double result = (which_stop == stop_mass::mstop_1 ? 1.0 : 2.0);
      return result;
   }

   double lowE6SSM_ew_derivs::deriv_d2MStop2_dLambdax_dmu222(stop_mass which_stop) const
   {
      double result = (which_stop == stop_mass::mstop_1 ? 1.0 : 2.0);
      return result;
   }

   double lowE6SSM_ew_derivs::deriv_d2MStop2_dTYu22_dmu222(stop_mass which_stop) const
   {
      double result = (which_stop == stop_mass::mstop_1 ? 1.0 : 2.0);
      return result;
   }

   double lowE6SSM_ew_derivs::deriv_d2MStop2_dLambdax_dLambdax(stop_mass which_stop) const
   {
      double result = (which_stop == stop_mass::mstop_1 ? 1.0 : 2.0);
      return result;
   }

   double lowE6SSM_ew_derivs::deriv_d2MStop2_dTYu22_dLambdax(stop_mass which_stop) const
   {
      double result = (which_stop == stop_mass::mstop_1 ? 1.0 : 2.0);
      return result;
   }

   double lowE6SSM_ew_derivs::deriv_d2MStop2_dTYu22_dTYu22(stop_mass which_stop) const
   {
      double result = (which_stop == stop_mass::mstop_1 ? 1.0 : 2.0);
      return result;
   }

   double lowE6SSM_ew_derivs::deriv_d2MStop2_dparam_dparam(stop_mass which_stop, lowE6SSM_info::Parameters p1, lowE6SSM_info::Parameters p2) const
   {
      if (equal_as_unordered_pairs(p1, p2, lowE6SSM_info::vd, lowE6SSM_info::vd)) {
         return deriv_d2MStop2_dvd_dvd(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowE6SSM_info::vd, lowE6SSM_info::vu)) {
         return deriv_d2MStop2_dvd_dvu(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowE6SSM_info::vd, lowE6SSM_info::vs)) {
         return deriv_d2MStop2_dvd_dvs(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowE6SSM_info::vd, lowE6SSM_info::g1)) {
         return deriv_d2MStop2_dg1_dvd(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowE6SSM_info::vd, lowE6SSM_info::g2)) {
         return deriv_d2MStop2_dg2_dvd(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowE6SSM_info::vd, lowE6SSM_info::gN)) {
         return deriv_d2MStop2_dgN_dvd(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowE6SSM_info::vd, lowE6SSM_info::Yu22)) {
         return deriv_d2MStop2_dYu22_dvd(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowE6SSM_info::vd, lowE6SSM_info::mq222)) {
         return deriv_d2MStop2_dmq222_dvd(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowE6SSM_info::vd, lowE6SSM_info::mu222)) {
         return deriv_d2MStop2_dmu222_dvd(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowE6SSM_info::vd, lowE6SSM_info::Lambdax)) {
         return deriv_d2MStop2_dLambdax_dvd(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowE6SSM_info::vd, lowE6SSM_info::TYu22)) {
         return deriv_d2MStop2_dTYu22_dvd(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowE6SSM_info::vu, lowE6SSM_info::vu)) {
         return deriv_d2MStop2_dvu_dvu(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowE6SSM_info::vu, lowE6SSM_info::vs)) {
         return deriv_d2MStop2_dvu_dvs(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowE6SSM_info::vu, lowE6SSM_info::g1)) {
         return deriv_d2MStop2_dg1_dvu(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowE6SSM_info::vu, lowE6SSM_info::g2)) {
         return deriv_d2MStop2_dg2_dvu(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowE6SSM_info::vu, lowE6SSM_info::gN)) {
         return deriv_d2MStop2_dgN_dvu(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowE6SSM_info::vu, lowE6SSM_info::Yu22)) {
         return deriv_d2MStop2_dYu22_dvu(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowE6SSM_info::vu, lowE6SSM_info::mq222)) {
         return deriv_d2MStop2_dmq222_dvu(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowE6SSM_info::vu, lowE6SSM_info::mu222)) {
         return deriv_d2MStop2_dmu222_dvu(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowE6SSM_info::vu, lowE6SSM_info::Lambdax)) {
         return deriv_d2MStop2_dLambdax_dvu(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowE6SSM_info::vu, lowE6SSM_info::TYu22)) {
         return deriv_d2MStop2_dTYu22_dvu(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowE6SSM_info::vs, lowE6SSM_info::vs)) {
         return deriv_d2MStop2_dvs_dvs(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowE6SSM_info::vs, lowE6SSM_info::g1)) {
         return deriv_d2MStop2_dg1_dvs(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowE6SSM_info::vs, lowE6SSM_info::g2)) {
         return deriv_d2MStop2_dg2_dvs(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowE6SSM_info::vs, lowE6SSM_info::gN)) {
         return deriv_d2MStop2_dgN_dvs(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowE6SSM_info::vs, lowE6SSM_info::Yu22)) {
         return deriv_d2MStop2_dYu22_dvs(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowE6SSM_info::vs, lowE6SSM_info::mq222)) {
         return deriv_d2MStop2_dmq222_dvs(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowE6SSM_info::vs, lowE6SSM_info::mu222)) {
         return deriv_d2MStop2_dmu222_dvs(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowE6SSM_info::vs, lowE6SSM_info::Lambdax)) {
         return deriv_d2MStop2_dLambdax_dvs(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowE6SSM_info::vs, lowE6SSM_info::TYu22)) {
         return deriv_d2MStop2_dTYu22_dvs(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowE6SSM_info::g1, lowE6SSM_info::g1)) {
         return deriv_d2MStop2_dg1_dg1(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowE6SSM_info::g1, lowE6SSM_info::g2)) {
         return deriv_d2MStop2_dg2_dg1(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowE6SSM_info::g1, lowE6SSM_info::gN)) {
         return deriv_d2MStop2_dgN_dg1(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowE6SSM_info::g1, lowE6SSM_info::Yu22)) {
         return deriv_d2MStop2_dYu22_dg1(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowE6SSM_info::g1, lowE6SSM_info::mq222)) {
         return deriv_d2MStop2_dmq222_dg1(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowE6SSM_info::g1, lowE6SSM_info::mu222)) {
         return deriv_d2MStop2_dmu222_dg1(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowE6SSM_info::g1, lowE6SSM_info::Lambdax)) {
         return deriv_d2MStop2_dLambdax_dg1(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowE6SSM_info::g1, lowE6SSM_info::TYu22)) {
         return deriv_d2MStop2_dTYu22_dg1(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowE6SSM_info::g2, lowE6SSM_info::g2)) {
         return deriv_d2MStop2_dg2_dg2(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowE6SSM_info::g2, lowE6SSM_info::gN)) {
         return deriv_d2MStop2_dgN_dg2(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowE6SSM_info::g2, lowE6SSM_info::Yu22)) {
         return deriv_d2MStop2_dYu22_dg2(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowE6SSM_info::g2, lowE6SSM_info::mq222)) {
         return deriv_d2MStop2_dmq222_dg2(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowE6SSM_info::g2, lowE6SSM_info::mu222)) {
         return deriv_d2MStop2_dmu222_dg2(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowE6SSM_info::g2, lowE6SSM_info::Lambdax)) {
         return deriv_d2MStop2_dLambdax_dg2(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowE6SSM_info::g2, lowE6SSM_info::TYu22)) {
         return deriv_d2MStop2_dTYu22_dg2(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowE6SSM_info::gN, lowE6SSM_info::gN)) {
         return deriv_d2MStop2_dgN_dgN(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowE6SSM_info::gN, lowE6SSM_info::Yu22)) {
         return deriv_d2MStop2_dYu22_dgN(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowE6SSM_info::gN, lowE6SSM_info::mq222)) {
         return deriv_d2MStop2_dmq222_dgN(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowE6SSM_info::gN, lowE6SSM_info::mu222)) {
         return deriv_d2MStop2_dmu222_dgN(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowE6SSM_info::gN, lowE6SSM_info::Lambdax)) {
         return deriv_d2MStop2_dLambdax_dgN(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowE6SSM_info::gN, lowE6SSM_info::TYu22)) {
         return deriv_d2MStop2_dTYu22_dgN(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowE6SSM_info::Yu22, lowE6SSM_info::Yu22)) {
         return deriv_d2MStop2_dYu22_dYu22(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowE6SSM_info::Yu22, lowE6SSM_info::mq222)) {
         return deriv_d2MStop2_dmq222_dYu22(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowE6SSM_info::Yu22, lowE6SSM_info::mu222)) {
         return deriv_d2MStop2_dmu222_dYu22(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowE6SSM_info::Yu22, lowE6SSM_info::Lambdax)) {
         return deriv_d2MStop2_dLambdax_dYu22(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowE6SSM_info::Yu22, lowE6SSM_info::TYu22)) {
         return deriv_d2MStop2_dTYu22_dYu22(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowE6SSM_info::mq222, lowE6SSM_info::mq222)) {
         return deriv_d2MStop2_dmq222_dmq222(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowE6SSM_info::mq222, lowE6SSM_info::mu222)) {
         return deriv_d2MStop2_dmu222_dmq222(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowE6SSM_info::mq222, lowE6SSM_info::Lambdax)) {
         return deriv_d2MStop2_dLambdax_dmq222(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowE6SSM_info::mq222, lowE6SSM_info::TYu22)) {
         return deriv_d2MStop2_dTYu22_dmq222(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowE6SSM_info::mu222, lowE6SSM_info::mu222)) {
         return deriv_d2MStop2_dmu222_dmu222(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowE6SSM_info::mu222, lowE6SSM_info::Lambdax)) {
         return deriv_d2MStop2_dLambdax_dmu222(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowE6SSM_info::mu222, lowE6SSM_info::TYu22)) {
         return deriv_d2MStop2_dTYu22_dmu222(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowE6SSM_info::Lambdax, lowE6SSM_info::Lambdax)) {
         return deriv_d2MStop2_dLambdax_dLambdax(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowE6SSM_info::Lambdax, lowE6SSM_info::TYu22)) {
         return deriv_d2MStop2_dTYu22_dLambdax(which_stop);
      } else if (equal_as_unordered_pairs(p1, p2, lowE6SSM_info::TYu22, lowE6SSM_info::TYu22)) {
         return deriv_d2MStop2_dTYu22_dTYu22(which_stop);
      } else {
         return 0.;
      }
   }
   
} // namespace flexiblesusy
