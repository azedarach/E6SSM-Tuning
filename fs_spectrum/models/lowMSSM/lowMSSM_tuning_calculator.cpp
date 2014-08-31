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

   return (g2*Cos(thW) + 0.7745966692414834*g1*Sin(thW));
}

double lowMSSM_tuning_calculator::MFtop_DRbar() const
{
   double yt = model.get_Yu(2,2);
   double vu = model.get_vu();

   return 0.7071067811865475*yt*vu;
}
 
double lowMSSM_tuning_calculator::stop_mass_matrix_LL_entry() const
{
   double mq2_22 = model.get_mq2(2,2);
   double mt = MFtop_DRbar();
   double vd = model.get_vd();
   double vu = model.get_vu();
   double g1 = model.get_g1();
   double g2 = model.get_g2();
   double delta_ul = 0.125*(Sqr(g2)-0.2*Sqr(g1))*(Sqr(vd)-Sqr(vu));

   return (mq2_22 + Sqr(mt) + delta_ul);
}

double lowMSSM_tuning_calculator::stop_mass_matrix_RR_entry() const
{
   double mu2_22 = model.get_mu2(2,2);
   double mt = MFtop_DRbar();
   double vd = model.get_vd();
   double vu = model.get_vu();
   double g1 = model.get_g1();
   double delta_ur = 0.1*Sqr(g1)*(Sqr(vd)-Sqr(vu));

   return (mu2_22 + Sqr(mt) + delta_ur);
}

double lowMSSM_tuning_calculator::stop_mass_matrix_LR_entry() const
{
   // DH:: currently assumes all parameters are *real*
   double at = model.get_TYu(2,2);
   double yt = model.get_Yu(2,2);
   double mu = model.get_Mu();   
   double vd = model.get_vd();
   double vu = model.get_vu();

   return (0.7071067811865475*(at*vu-mu*yt*vd));
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

   return (Sqr(g2)-Sqr(g1))*mqqsq;
}

double lowMSSM_tuning_calculator::stop_discriminant() const
{
   double mqqsq = MQQ2();
   double mass_matrix_stop_LR = stop_mass_matrix_LR_entry();

   return (Sqr(mqqsq) + 4.0*AbsSqr(mass_matrix_stop_LR));
}

void lowMSSM_tuning_calculator::calculate_MStop() 
{
   // For now just use the explicit solutions for the mass
   // eigenvalues; avoids any mass ordering issues.
   double stop_mass_matrix_trace = stop_mass_matrix_LL_entry() 
     + stop_mass_matrix_RR_entry();
   double disc = stop_discriminant();

   MStop(0) = 0.5*(stop_mass_matrix_trace - Sqrt(disc));
   MStop(1) = 0.5*(stop_mass_matrix_trace + Sqrt(disc));

   /// DH:: Check that this is appropriate action
   if (MStop.minCoeff() < 0.) model.get_problems().flag_tachyon(lowMSSM_info::Su);

   MStop = AbsSqrt(MStop);
}

double lowMSSM_tuning_calculator::deriv_dMStop2_dv1(stop_mass which_stop) const
{
   double mu = model.get_Mu();
   double vd = model.get_vd();
   double g1 = model.get_g1();
   double g2 = model.get_g2();
   double yt = model.get_Yu(2,2);
   double rqq_val = RQQ();
   double stop_mixing = 1.4142135623730951*stop_mass_matrix_LR_entry();
   double rt = stop_discriminant();

   double split = (0.25*rqq_val - 2.0*yt*mu*stop_mixing/vd)/Sqrt(rt);

   double deriv = 0.25*Sqr(g2) + 0.15 * Sqr(g1);

   if (which_stop == stop_mass::mstop_1)
     {
       return 0.5*vd*(deriv - split);
     }

   return 0.5*vd*(deriv + split);
}

double lowMSSM_tuning_calculator::deriv_dMStop2_dv2(stop_mass which_stop) const
{
   double at = model.get_TYu(2,2);
   double vu = model.get_vu();
   double g1 = model.get_g1();
   double g2 = model.get_g2();
   double yt = model.get_Yu(2,2);
   double rqq_val = RQQ();
   double stop_mixing = 1.4142135623730951*stop_mass_matrix_LR_entry();
   double rt = stop_discriminant();

   double split = (2.0*at*stop_mixing/vu - 0.25*rqq_val)/Sqrt(rt);

   double deriv = 2.0*Sqr(yt) - 0.25*Sqr(g2) - 0.15 * Sqr(g1);

   if (which_stop == stop_mass::mstop_1)
     {
       return 0.5*vu*(deriv - split);
     }

   return 0.5*vu*(deriv + split);
}

double lowMSSM_tuning_calculator::deriv_d2DeltaV_dvd_dvd() const
{
   double scale = model.get_scale();

   double mu = model.get_Mu();
   double vd = model.get_vd();
   double yt = model.get_Yu(2,2);
   double g1 = model.get_g1();
   double g2 = model.get_g2();

   double gbar_val = gbar();
   double stop_mixing = 1.4142135623730951*stop_mass_matrix_LR_entry();
   double RQQ_val = RQQ();
   double rt = stop_discriminant();

   double a0_mstop1 = passarino_veltman::ReA0(Sqr(MStop(0)), Sqr(scale));
   double a0_mstop2 = passarino_veltman::ReA0(Sqr(MStop(1)), Sqr(scale));

   double tmp_1 = 0.;
   tmp_1 += Sqr(0.125*Sqr(gbar_val)*vd) 
     + Sqr(0.125*vd*RQQ_val-yt*mu*stop_mixing)/rt;
   tmp_1 *= Log(Sqr(MStop(0)*MStop(1))/Sqr(scale*scale));

   double tmp_2 = 0.;
   tmp_2 += (0.03125*vd*Sqr(gbar_val)/Sqrt(rt))
     *(vd*RQQ_val-8.0*yt*mu*stop_mixing);
   tmp_2 *= Log(Sqr(MStop(1))/Sqr(MStop(0)));

   double tmp_3 = -0.125*Sqr(gbar_val)*(a0_mstop1+a0_mstop2);

   double tmp_4 = 0.;
   tmp_4 += 0.03125*((4.0*RQQ_val+Sqr((g2*g2-g1*g1)*vd)+32.0*Sqr(yt*mu))/Sqrt(rt)
		     -Sqr(vd*RQQ_val-8.0*yt*mu*stop_mixing)/(rt*Sqrt(rt)));
   tmp_4 *= (a0_mstop1-a0_mstop2);

  return (3.0*oneOver16PiSqr*(tmp_1+tmp_2+tmp_3+tmp_4));
}

double lowMSSM_tuning_calculator::deriv_d2DeltaV_dvu_dvu() const
{
   double scale = model.get_scale();

   double vu = model.get_vu();
   double yt = model.get_Yu(2,2);
   double at = model.get_TYu(2,2);
   double g1 = model.get_g1();
   double g2 = model.get_g2();

   double gbar_val = gbar();
   double stop_mixing = 1.4142135623730951*stop_mass_matrix_LR_entry();
   double RQQ_val = RQQ();
   double rt = stop_discriminant();
   double mt = MFtop_DRbar();

   double a0_mstop1 = passarino_veltman::ReA0(Sqr(MStop(0)), Sqr(scale));
   double a0_mstop2 = passarino_veltman::ReA0(Sqr(MStop(1)), Sqr(scale));
   double a0_mtop = passarino_veltman::ReA0(Sqr(mt), Sqr(scale));

   double tmp_1 = 0.;
   tmp_1 += Sqr(yt*yt-0.125*Sqr(gbar_val)) 
     + Sqr(8.0*at*stop_mixing/vu-RQQ_val)/(64.0*rt);
   tmp_1 *= Sqr(vu)*Log(Sqr(MStop(0)*MStop(1))/Sqr(scale*scale));

   double tmp_2 = 0.;
   tmp_2 += Sqr(vu)*(yt*yt-0.125*Sqr(gbar_val))
     *(8.0*at*stop_mixing/vu-RQQ_val)/(4.0*Sqrt(rt));
   tmp_2 *= Log(Sqr(MStop(1))/Sqr(MStop(0)));

   double tmp_3 = -(yt*yt-0.125*Sqr(gbar_val))*(a0_mstop1+a0_mstop2);

   double tmp_4 = 0.;
   tmp_4 += (0.03125*Sqr(vu*(g2*g2-g1*g1)) - 0.125*RQQ_val+Sqr(at)
	     - 0.03125*Sqr(vu)*Sqr(8.0*at*stop_mixing/vu-RQQ_val)/rt)/Sqrt(rt);
   tmp_4 *= (a0_mstop1 - a0_mstop2);

   double tmp_5 = -2.0*Sqr(yt*yt*vu)*Log(Sqr(mt)/Sqr(scale)) 
     + 2.0*Sqr(yt)*a0_mtop;

   return (3.0*oneOver16PiSqr*(tmp_1+tmp_2+tmp_3+tmp_4+tmp_5));
}

double lowMSSM_tuning_calculator::deriv_d2DeltaV_dvu_dvd() const
{
   double scale = model.get_scale();

   double mu = model.get_Mu();
   double vd = model.get_vd();
   double vu = model.get_vu();
   double yt = model.get_Yu(2,2);
   double at = model.get_TYu(2,2);
   double g1 = model.get_g1();
   double g2 = model.get_g2();

   double gbar_val = gbar();
   double stop_mixing = 1.4142135623730951*stop_mass_matrix_LR_entry();
   double RQQ_val = RQQ();
   double rt = stop_discriminant();

   double a0_mstop1 = passarino_veltman::ReA0(Sqr(MStop(0)), Sqr(scale));
   double a0_mstop2 = passarino_veltman::ReA0(Sqr(MStop(1)), Sqr(scale));

   double tmp_1 = 0.;
   tmp_1 += 0.125*Sqr(gbar_val)*(yt*yt-0.125*Sqr(gbar_val))
     + ((0.125*RQQ_val-yt*mu*stop_mixing/vd)
	*(at*stop_mixing/vu-0.125*RQQ_val))/rt;
   tmp_1 *= vu*vd*Log(Sqr(MStop(0)*MStop(1))/Sqr(scale*scale));

   double tmp_2 = 0.;
   tmp_2 += 0.125*Sqr(gbar_val)*(at*stop_mixing/vu-0.125*RQQ_val)
     + (yt*yt-0.125*Sqr(gbar_val))*(0.125*RQQ_val-yt*mu*stop_mixing/vd);
   tmp_2 *= vu*vd*Log(Sqr(MStop(1))/Sqr(MStop(0)))/Sqrt(rt);

   double tmp_3 = 0.;
   tmp_3 += 0.03125*vu*vd*Sqr(g2*g2-g1*g1) + yt*mu*at
     +2.0*vu*vd*((0.125*RQQ_val-yt*mu*stop_mixing/vd)
		 *(at*stop_mixing/vu-0.125*RQQ_val))/rt;
   tmp_3 *= (a0_mstop2-a0_mstop1)/Sqrt(rt);

   return (3.0*oneOver16PiSqr*(tmp_1+tmp_2+tmp_3));
}

} // namespace flexiblesusy
