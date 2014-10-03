#include "lowMSSM_ewsb_conditions.hpp"

namespace flexiblesusy {
   
   lowMSSM_ewsb_conditions::lowMSSM_ewsb_conditions(const lowMSSM<Two_scale>& m)
      : model(m)
   {

   }
   
   Eigen::Matrix<double,2,2> lowMSSM_ewsb_conditions::calculate_unrotated_mass_matrix_hh() const
   {
      Eigen::Matrix<double,2,2> mass_matrix(model.get_mass_matrix_hh());
   
      // Add in tadpole corrections
      if (model.get_ewsb_loop_order() > 0) {
         mass_matrix(0,0) += deriv_d2DeltaV_dvd_dvd();
         mass_matrix(0,1) += deriv_d2DeltaV_dvd_dvu();
         mass_matrix(1,0) += deriv_d2DeltaV_dvd_dvu();
         mass_matrix(1,1) += deriv_d2DeltaV_dvu_dvu();
      }
      
      return mass_matrix;
   }
//TODO
   Eigen::Matrix<double,2,1> lowMSSM_ewsb_conditions::calculate_derivs_dewsb_eqs_dparam(lowMSSM_info::Parameters p) const
   {
      Eigen::Matrix<double,number_of_ewsb_eqs,1> derivs;

      return derivs;
   }
//TODO
   int lowMSSM_ewsb_conditions::solve_ewsb_conditions_for_vevs()
   {
      return 0;
   }

   double lowMSSM_ewsb_conditions::gbar() const
   {
      double g1  = model.get_g1();
      double g2  = model.get_g2();
      double thW = model.ThetaW();

      return (g2 * Cos(thW) + 0.7745966692414834 * g1 * Sin(thW));
   }

   double lowMSSM_ewsb_conditions::MFtop_DRbar() const
   {
      double yt = model.get_Yu(2, 2);
      double vu = model.get_vu();

      return 0.7071067811865475 * yt * vu;
   }
   
   double lowMSSM_ewsb_conditions::MFbottom_DRbar() const
   {
      double yb = model.get_Yd(2,2);
      double vd = model.get_vd();

      return 0.7071067811865475 * yb * vd;
   }

   double lowMSSM_ewsb_conditions::stop_mass_matrix_LL_entry() const
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

   double lowMSSM_ewsb_conditions::stop_mass_matrix_RR_entry() const
   {
      double mu2_22 = model.get_mu2(2, 2);
      double mt = MFtop_DRbar();
      double vd = model.get_vd();
      double vu = model.get_vu();
      double g1 = model.get_g1();
      double delta_ur = 0.1 * Sqr(g1) * (Sqr(vd) - Sqr(vu));

      return (mu2_22 + Sqr(mt) + delta_ur);
   }

   double lowMSSM_ewsb_conditions::stop_mass_matrix_LR_entry() const
   {
      // Assumes all parameters are real
      double at = model.get_TYu(2, 2);
      double yt = model.get_Yu(2, 2);
      double mu = model.get_Mu();
      double vd = model.get_vd();
      double vu = model.get_vu();

      return (0.7071067811865475 * (at * vu - mu * yt * vd));
   }
   
   double lowMSSM_ewsb_conditions::stop_MQQ2() const
   {
      double mass_matrix_stop_LL = stop_mass_matrix_LL_entry();
      double mass_matrix_stop_RR = stop_mass_matrix_RR_entry();

      return (mass_matrix_stop_LL - mass_matrix_stop_RR);
   }

   double lowMSSM_ewsb_conditions::stop_RQQ() const
   {
      double g1 = model.get_g1();
      double g2 = model.get_g2();
      double mqqsq = stop_MQQ2();

      return (Sqr(g2) - Sqr(g1)) * mqqsq;
   }

   double lowMSSM_ewsb_conditions::stop_discriminant() const
   {
      double mqqsq = stop_MQQ2();
      double mass_matrix_stop_LR = stop_mass_matrix_LR_entry();

      return (Sqr(mqqsq) + 4.0 * AbsSqr(mass_matrix_stop_LR));
   }
   
   void lowMSSM_ewsb_conditions::calculate_MStop()
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

   double lowMSSM_ewsb_conditions::sbottom_mass_matrix_LL_entry() const
   {
      double mq2_22 = model.get_mq2(2,2);
      double mb = MFbottom_DRbar();
      double vd = model.get_vd();
      double vu = model.get_vu();
      double g1 = model.get_g1();
      double g2 = model.get_g2();
      double delta_dl = 0.125 * (Sqr(g2) + 0.2 * Sqr(g1)) * (Sqr(vu) - Sqr(vd));

      return (mq2_22 + Sqr(mb) + delta_dl);
   }

   double lowMSSM_ewsb_conditions::sbottom_mass_matrix_RR_entry() const
   {
      double md2_22 = model.get_md2(2,2);
      double mb = MFbottom_DRbar();
      double vd = model.get_vd();
      double vu = model.get_vu();
      double g1 = model.get_g1();
      double delta_dr = 0.05 * Sqr(g1) * (Sqr(vu) - Sqr(vd));

      return (md2_22 + Sqr(mb) + delta_dr);
   }

   double lowMSSM_ewsb_conditions::sbottom_mass_matrix_LR_entry() const
   {
      double ab = model.get_TYd(2,2);
      double yb = model.get_Yd(2,2);
      double mu = model.get_Mu();
      double vd = model.get_vd();
      double vu = model.get_vu();

      return (0.7071067811865475 * (ab * vd - mu * yb * vu))
   }
   
   double lowMSSM_ewsb_conditions::sbottom_MQQ2() const
   {
      double mass_matrix_sbottom_LL = sbottom_mass_matrix_LL_entry();
      double mass_matrix_sbottom_RR = sbottom_mass_matrix_RR_entry();

      return (mass_matrix_sbottom_LL - mass_matrix_sbottom_RR);
   }

   double lowMSSM_ewsb_conditions::sbottom_discriminant() const
   {
      double mqqsq = sbottom_MQQ2();
      double mass_matrix_sbottom_LR = sbottom_mass_matrix_LR_entry();

      return (Sqr(mqqsq) + 4.0 * AbsSqr(mass_matrix_sbottom_LR));
   }

   void lowMSSM_ewsb_conditions::calculate_MSbottom()
   {
      // For now just use the explicit solutions for the mass
      // eigenvalues; avoids any mass ordering issues.
      double sbottom_mass_matrix_trace = sbottom_mass_matrix_LL_entry()
         + sbottom_mass_matrix_RR_entry();
      double disc = sbottom_discriminant();

      MSbottom(0) = 0.5 * (sbottom_mass_matrix_trace - Sqrt(disc));
      MSbottom(1) = 0.5 * (sbottom_mass_matrix_trace + Sqrt(disc));

      if (MSbottom.minCoeff() < 0.) model.get_problems().flag_tachyon(lowMSSM_info::Sd);

      MSbottom = AbsSqrt(MSbottom);
   }
   
   double lowMSSM_ewsb_conditions::deriv_d2DeltaV_dparam_dparam(lowMSSM_info::Parameters p1, lowMSSM_info::Parameters p2) const
   {
      double scale = model.get_scale();

      double dMStop20_dp1 = deriv_dMStop2_dparam(stop_mass::mstop_1, p1);
      double dMStop20_dp2 = deriv_dMStop2_dparam(stop_mass::mstop_1, p2);
      double d2MStop20_dp1_dp2 = deriv_d2MStop2_dparam_dparam(stop_mass::mstop_1, p1, p2);

      double dMStop21_dp1 = deriv_dMStop2_dparam(stop_mass::mstop_2, p1);
      double dMStop21_dp2 = deriv_dMStop2_dparam(stop_mass::mstop_2, p2);
      double d2MStop21_dp1_dp2 = deriv_d2MStop2_dparam_dparam(stop_mass::mstop_2, p1, p2);

      double dMFtop2_dp1 = deriv_dMFtop2_dparam(p1);
      double dMFtop2_dp2 = deriv_dMFtop2_dparam(p2);
      double d2MFtop2_dp1_dp2 = deriv_d2MFtop2_dparam_dparam(p1, p2);

      double mt = MFtop_DRbar();

      double a0_mstop1 = passarino_veltman::ReA0(Sqr(MStop(0)), Sqr(scale));
      double a0_mstop2 = passarino_veltman::ReA0(Sqr(MStop(1)), Sqr(scale));
      double a0_mtop = passarino_veltman::ReA0(Sqr(mt), Sqr(scale));

      double tmp_1 = dMStop20_dp1 * dMStop20_dp2 * Log(Sqr(MStop(0) / scale));
      double tmp_2 = -a0_mstop1 * d2MStop20_dp1_dp2;
      double tmp_3 = dMStop21_dp1 * dMStop21_dp2 * Log(Sqr(MStop(1) / scale));
      double tmp_4 = -a0_mstop2 * d2MStop21_dp1_dp2;
      double tmp_5 = -2.0 * dMFtop2_dp1 * dMFtop2_dp2 * Log(Sqr(mt / scale));
      double tmp_6 = 2.0 * a0_mtop * d2MFtop2_dp1_dp2;

      double result = 3.0 * oneOver16PiSqr * (tmp_1 + tmp_2 + tmp_3 + tmp_4 + tmp_5 + tmp_6);

      if (include_sbottom_loops) {
         double dMSbottom20_dp1 = deriv_dMSbottom2_dparam(sbottom_mass::msbottom_1, p1);
         double dMSbottom20_dp2 = deriv_dMSbottom2_dparam(sbottom_mass::msbottom_1, p2);
         double d2MSbottom20_dp1_dp2 = deriv_d2MSbottom2_dparam_dparam(sbottom_mass::msbottom_1, p1, p2);

         double dMSbottom21_dp1 = deriv_dMSbottom2_dparam(sbottom_mass::msbottom_2, p1);
         double dMSbottom21_dp2 = deriv_dMSbottom2_dparam(sbottom_mass::msbottom_2, p2);
         double d2MSbottom21_dp1_dp2 = deriv_d2MSbottom2_dparam_dparam(sbottom_mass::msbottom_2, p1, p2);

         double dMFbottom2_dp1 = deriv_dMFbottom2_dparam(p1);
         double dMFbottom2_dp2 = deriv_dMFbottom2_dparam(p2);
         double d2MFbottom2_dp1_dp2 = deriv_d2MFbottom2_dparam_dparam(p1, p2);
         
         double mb = MFbottom_DRbar();
         
         double a0_msbottom1 = passarino_veltman::ReA0(Sqr(MSbottom(0)), Sqr(scale));
         double a0_msbottom2 = passarino_veltman::ReA0(Sqr(MSbottom(1)), Sqr(scale));
         double a0_mbottom = passarino_veltman::ReA0(Sqr(mb), Sqr(scale));

         double sb_tmp_1 = dMSbottom20_dp1 * dMSbottom20_dp2 * Log(Sqr(MSbottom(0) / scale));
         double sb_tmp_2 = -a0_msbottom1 * d2MSbottom20_dp1_dp2;
         double sb_tmp_3 = dMSbottom21_dp1 * dMSbottom21_dp2 * Log(Sqr(MSbottom(1) / scale));
         double sb_tmp_4 = -a0_msbottom2 * d2MSbottom21_dp1_dp2;
         double sb_tmp_5 = -2.0 * dMFbottom2_dp1 * dMFbottom2_dp2 * Log(Sqr(mt / scale));
         double sb_tmp_6 = 2.0 * a0_mbottom * d2MFbottom2_dp1_dp2;

         result += 3.0 * oneOver16PiSqr * (sb_tmp_1 + sb_tmp_2 + sb_tmp_3 + sb_tmp_4 + sb_tmp_5 + sb_tmp_6);

      }

      return result;
   }

} // namespace flexiblesusy
