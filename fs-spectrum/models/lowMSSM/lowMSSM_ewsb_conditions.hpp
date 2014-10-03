// ====================================================================
// Helper class providing EWSB conditions and derivatives for use
// in the fine tuning calculation.
// Notes:
//   - currently we only include stop contributions to the 
//     Coleman-Weinberg potential
// ====================================================================

#ifndef lowMSSM_EWSB_CONDITIONS_H
#define lowMSSM_EWSB_CONDITIONS_H

#include "lowMSSM_two_scale_model.hpp"

#include <Eigen/Core>

namespace flexiblesusy {

   class lowMSSM_ewsb_conditions {
   public:

      enum class stop_mass : char {mstop_1, mstop_2};
      enum class sbottom_mass : char {msbottom_1, msbottom_2};

      explicit lowMSSM_ewsb_conditions(const lowMSSM<Two_scale>&);

      double get_ewsb_loop_order() const { return model.get_ewsb_loop_order(); }

      void set_ewsb_loop_order(unsigned l) { model.set_ewsb_loop_order(l); }

      bool get_include_sbottom_loops() const { return include_sbottom_loops; }

      void set_include_sbottom_loops(bool b) { include_sbottom_loops = b; }

      /// Values of EWSB conditions
      double get_ewsb_condition_1() const;
      double get_ewsb_condition_2() const;

      /// Derivatives of EWSB conditions w.r.t. VEVs
      Eigen::Matrix<double,2,2> calculate_unrotated_mass_matrix_hh() const;

      /// Derivatives of EWSB conditions w.r.t. model parameters
      Eigen::Matrix<double,2,1> calculate_derivs_dewsb_eqs_dparam(lowMSSM_info::Parameters) const;

      /// Function to solve for VEVs given parameters
      int solve_ewsb_conditions_for_vevs();

      /// Derivatives of DR bar masses appearing in the Coleman-Weinberg potential
      double deriv_dMFtop2_dparam(lowMSSM_info::Parameters p) const;
      double deriv_dMFbottom2_dparam(lowMSSM_info::Parameters p) const;
      double deriv_dMStop2_dparam(stop_mass which_stop, lowMSSM_info::Parameters p) const;
      double deriv_dMSbottom2_dparam(sbottom_mass which_sbottom, lowMSSM_info::Parameters p) const;
      double deriv_d2MFtop2_dparam_dparam(lowMSSM_info::Parameters p1, lowMSSM_info::Parameters p2) const;
      double deriv_d2MFbottom2_dparam_dparam(lowMSSM_info::Parameters p1, lowMSSM_info::Parameters p2) const;
      double deriv_d2MStop2_dparam_dparam(stop_mass which_stop, lowMSSM_info::Parameters p1, 
                                         lowMSSM_info::Parameters p2) const;
      double deriv_d2MSbottom2_dparam_dparam(sbottom_mass which_sbottom, lowMSSM_info::Parameters p1, 
                                         lowMSSM_info::Parameters p2) const;

      /// Derivatives of the Coleman-Weinberg loop contributions. Note that by default
      /// only top and stop loops are included.
      double deriv_d2DeltaV_dparam_dparam(lowMSSM_info::Parameters p1, lowMSSM_info::Parameters p2) const;

   private:
      static const std::size_t number_of_ewsb_eqs = 2;

      lowMSSM<Two_scale> model;
      bool include_sbottom_loops;

      Eigen::Matrix<double,2,1> MStop;
      Eigen::Matrix<double,2,1> MSbottom;

      double gbar() const;
      double MFtop_DRbar() const;
      double MFbottom_DRbar() const;

      double stop_mass_matrix_LL_entry() const;
      double stop_mass_matrix_RR_entry() const;
      double stop_mass_matrix_LR_entry() const;

      double stop_MQQ2() const;
      double stop_RQQ() const;
      double stop_discriminant() const;

      double sbottom_mass_matrix_LL_entry() const;
      double sbottom_mass_matrix_RR_entry() const;
      double sbottom_mass_matrix_LR_entry() const;

      double sbottom_MQQ2() const;
      double sbottom_discriminant() const;

      void calculate_MStop();
      void calculate_MSbottom();

      double deriv_dMFtop2_dvu() const;
      double deriv_dMFtop2_dYu22() const;
      double deriv_d2MFtop2_dvu_dvu() const;
      double deriv_d2MFtop2_dYu22_dvu() const;
      double deriv_d2MFtop2_dYu22_dYu22() const;

      double deriv_dMFbottom2_dvd() const;
      double deriv_dMFbottom2_dYd22() const;
      double deriv_d2MFbottom2_dvd_dvd() const;
      double deriv_d2MFbottom2_dYd22_dvd() const;
      double deriv_d2MFbottom2_dYd22_dYd22() const;

      double deriv_dMStop2_dvd(stop_mass which_stop) const;
      double deriv_dMStop2_dvu(stop_mass which_stop) const;
      double deriv_dMStop2_dg1(stop_mass which_stop) const;
      double deriv_dMStop2_dg2(stop_mass which_stop) const;
      double deriv_dMStop2_dYu22(stop_mass which_stop) const;
      double deriv_dMStop2_dmq222(stop_mass which_stop) const;
      double deriv_dMStop2_dmu222(stop_mass which_stop) const;
      double deriv_dMStop2_dMu(stop_mass which_stop) const;
      double deriv_dMStop2_dTYu22(stop_mass which_stop) const;

      double deriv_d2MStop2_dvd_dvd(stop_mass which_stop) const;
      double deriv_d2MStop2_dvu_dvu(stop_mass which_stop) const;
      double deriv_d2MStop2_dvd_dvu(stop_mass which_stop) const;
      double deriv_d2MStop2_dmq222_dvd(stop_mass which_stop) const;
      double deriv_d2MStop2_dmu222_dvd(stop_mass which_stop) const;
      double deriv_d2MStop2_dMu_dvd(stop_mass which_stop) const;
      double deriv_d2MStop2_dTYu22_dvd(stop_mass which_stop) const;
      double deriv_d2MStop2_dmq222_dvu(stop_mass which_stop) const;
      double deriv_d2MStop2_dmu222_dvu(stop_mass which_stop) const;
      double deriv_d2MStop2_dMu_dvu(stop_mass which_stop) const;
      double deriv_d2MStop2_dTYu22_dvu(stop_mass which_stop) const;

      double deriv_dMSbottom2_dvd(sbottom_mass which_sbottom) const;
      double deriv_dMSbottom2_dvu(sbottom_mass which_sbottom) const;
      double deriv_dMSbottom2_dg1(sbottom_mass which_sbottom) const;
      double deriv_dMSbottom2_dg2(sbottom_mass which_sbottom) const;
      double deriv_dMSbottom2_dYd22(sbottom_mass which_sbottom) const;
      double deriv_dMSbottom2_dmq222(sbottom_mass which_sbottom) const;
      double deriv_dMSbottom2_dmd222(sbottom_mass which_sbottom) const;
      double deriv_dMSbottom2_dMu(sbottom_mass which_sbottom) const;
      double deriv_dMSbottom2_dTYd22(sbottom_mass which_sbottom) const;
      double deriv_d2MSbottom2_dvd_dvd(sbottom_mass which_sbottom) const;
      double deriv_d2MSbottom2_dvu_dvu(sbottom_mass which_sbottom) const;
      double deriv_d2MSbottom2_dvd_dvu(sbottom_mass which_sbottom) const;
      double deriv_d2MSbottom2_dmq222_dvd(sbottom_mass which_sbottom) const;
      double deriv_d2MSbottom2_dmd222_dvd(sbottom_mass which_sbottom) const;
      double deriv_d2MSbottom2_dMu_dvd(sbottom_mass which_sbottom) const;
      double deriv_d2MSbottom2_dTYd22_dvd(sbottom_mass which_sbottom) const;
      double deriv_d2MSbottom2_dmq222_dvu(sbottom_mass which_sbottom) const;
      double deriv_d2MSbottom2_dmd222_dvu(sbottom_mass which_sbottom) const;
      double deriv_d2MSbottom2_dMu_dvu(sbottom_mass which_sbottom) const;
      double deriv_d2MSbottom2_dTYd22_dvu(sbottom_mass which_sbottom) const;
   };

} // namespace flexiblesusy

#endif
