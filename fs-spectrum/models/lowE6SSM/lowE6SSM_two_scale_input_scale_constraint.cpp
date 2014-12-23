// ====================================================================
// Additional constraint for fixing parameters at an input scale
// other than that defined by g1 == g2 or M_{SUSY}
// ====================================================================

#include "lowE6SSM_two_scale_input_scale_constraint.hpp"
#include "lowE6SSM_two_scale_model.hpp"
#include "wrappers.hpp"
#include "logger.hpp"
#include "ew_input.hpp"
#include "gsl_utils.hpp"
#include "minimizer.hpp"
#include "root_finder.hpp"
#include "numerics.hpp"

#include <cassert>
#include <cmath>
#include <cerrno>
#include <cstring>

namespace flexiblesusy {

#define INPUTPARAMETER(p) inputPars.p
#define MODELPARAMETER(p) model->get_##p()
#define BETAPARAMETER(p) beta_functions.get_##p()
#define BETA(p) beta_##p
#define SM(p) Electroweak_constrants::p
#define STANDARDDEVIATION(p) Electroweak_constants::Error_##p
#define Pole(p) model->get_physical().p
#define MODEL model
#define MODELCLASSNAME lowE6SSM<Two_scale>

   lowE6SSM_input_scale_constraint<Two_scale>::lowE6SSM_input_scale_constraint()
      : Constraint<Two_scale>()
      , initial_scale_guess(0.)
      , fixed_scale(0.)
      , model(0)
      , inputPars()
   {
   }

   lowE6SSM_input_scale_constraint<Two_scale>::lowE6SSM_input_scale_constraint(const lowE6SSM_input_parameters& inputPars_)
      : Constraint<Two_scale>()
      , fixed_scale(0.)
      , model(0)
      , inputPars(inputPars_)
   {
      initialize();
   }

   lowE6SSM_input_scale_constraint<Two_scale>::~lowE6SSM_input_scale_constraint()
   {
   }

   void lowE6SSM_input_scale_constraint<Two_scale>::apply()
   {
      assert(model && "Error: lowE6SSM_input_scale_constraint:"
             " model pointer must not be zero");

      // apply user defined input scale constraints
      const auto KappaInput = INPUTPARAMETER(KappaInput);
      const auto Lambda12Input = INPUTPARAMETER(Lambda12Input);
      const auto LambdaxInput = INPUTPARAMETER(LambdaxInput);
      const auto MuPrInput = INPUTPARAMETER(MuPrInput);
      const auto AYdInput = INPUTPARAMETER(AYdInput);
      const auto AYeInput = INPUTPARAMETER(AYeInput);
      const auto TKappaInput = INPUTPARAMETER(TKappaInput);
      const auto TLambda12Input = INPUTPARAMETER(TLambda12Input);
      const auto TLambdaxInput = INPUTPARAMETER(TLambdaxInput);
      const auto AYuInput = INPUTPARAMETER(AYuInput);
      const auto BMuPrInput = INPUTPARAMETER(BMuPrInput);
      const auto mq2Input = INPUTPARAMETER(mq2Input);
      const auto ml2Input = INPUTPARAMETER(ml2Input);
      const auto md2Input = INPUTPARAMETER(md2Input);
      const auto mu2Input = INPUTPARAMETER(mu2Input);
      const auto me2Input = INPUTPARAMETER(me2Input);
      const auto mH1I2Input = INPUTPARAMETER(mH1I2Input);
      const auto mH2I2Input = INPUTPARAMETER(mH2I2Input);
      const auto msI2Input = INPUTPARAMETER(msI2Input);
      const auto mDx2Input = INPUTPARAMETER(mDx2Input);
      const auto mDxbar2Input = INPUTPARAMETER(mDxbar2Input);
      const auto mHp2Input = INPUTPARAMETER(mHp2Input);
      const auto mHpbar2Input = INPUTPARAMETER(mHpbar2Input);
      const auto MassBInput = INPUTPARAMETER(MassBInput);
      const auto MassWBInput = INPUTPARAMETER(MassWBInput);
      const auto MassGInput = INPUTPARAMETER(MassGInput);
      const auto MassBpInput = INPUTPARAMETER(MassBpInput);

      const auto Yd = MODELPARAMETER(Yd);
      const auto Ye = MODELPARAMETER(Ye);
      const auto Yu = MODELPARAMETER(Yu);

      MODEL->set_Kappa(KappaInput);
      MODEL->set_Lambda12(Lambda12Input);
      MODEL->set_Lambdax(LambdaxInput);
      MODEL->set_MuPr(MuPrInput);
      MODEL->set_TYd((Yd.array() * AYdInput.array()).matrix());
      MODEL->set_TYe((Ye.array() * AYeInput.array()).matrix());
      MODEL->set_TKappa(TKappaInput);
      MODEL->set_TLambda12(TLambda12Input);
      MODEL->set_TLambdax(TLambdaxInput);
      MODEL->set_TYu((Yu.array() * AYuInput.array()).matrix());
      MODEL->set_BMuPr(BMuPrInput);
      MODEL->set_mq2(mq2Input);
      MODEL->set_ml2(ml2Input);
      MODEL->set_md2(md2Input);
      MODEL->set_mu2(mu2Input);
      MODEL->set_me2(me2Input);
      MODEL->set_mH1I2(mH1I2Input);
      MODEL->set_mH2I2(mH2I2Input);
      MODEL->set_msI2(msI2Input);
      MODEL->set_mDx2(mDx2Input);
      MODEL->set_mDxbar2(mDxbar2Input);
      MODEL->set_mHp2(mHp2Input);
      MODEL->set_mHpbar2(mHpbar2Input);
      MODEL->set_MassB(MassBInput);
      MODEL->set_MassWB(MassWBInput);
      MODEL->set_MassG(MassGInput);
      MODEL->set_MassBp(MassBpInput);

   }

   double lowE6SSM_input_scale_constraint<Two_scale>::get_scale() const
   {
      return fixed_scale;
   }

   double lowE6SSM_input_scale_constraint<Two_scale>::get_initial_scale_guess() const
   {
      return initial_scale_guess;
   }

   void lowE6SSM_input_scale_constraint<Two_scale>::set_model(Two_scale_model* model_)
   {
      model = cast_model<lowE6SSM<Two_scale> >(model_);
   }

   void lowE6SSM_input_scale_constraint<Two_scale>::set_input_parameters(const lowE6SSM_input_parameters& inputPars_)
   {
      inputPars = inputPars_;
   }

   void lowE6SSM_input_scale_constraint<Two_scale>::set_scale(double s)
   {
      fixed_scale = s;
   }

   void lowE6SSM_input_scale_constraint<Two_scale>::clear()
   {
      initial_scale_guess = 0.;
      fixed_scale = 0.;
      model = NULL;
   }

   void lowE6SSM_input_scale_constraint<Two_scale>::initialize()
   {
      initial_scale_guess = 2.e16;

      fixed_scale = initial_scale_guess;
   }

} // namespace flexiblesusy
