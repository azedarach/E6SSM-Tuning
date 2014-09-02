// ====================================================================
// This file is part of FlexibleSUSY.
//
// FlexibleSUSY is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// FlexibleSUSY is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with FlexibleSUSY.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

// File generated at Sun 24 Aug 2014 16:16:19

#include "lowMSSM_two_scale_susy_scale_constraint.hpp"
#include "lowMSSM_two_scale_model.hpp"
#include "wrappers.hpp"
#include "logger.hpp"
#include "ew_input.hpp"
#include "gsl_utils.hpp"
#include "minimizer.hpp"
#include "root_finder.hpp"

#include <cassert>
#include <cmath>

namespace flexiblesusy {

#define INPUTPARAMETER(p) inputPars.p
#define MODELPARAMETER(p) model->get_##p()
#define BETAPARAMETER(p) beta_functions.get_##p()
#define BETA(p) beta_##p
#define SM(p) Electroweak_constants::p
#define STANDARDDEVIATION(p) Electroweak_constants::Error_##p
#define Pole(p) model->get_physical().p
#define MODEL model
#define MODELCLASSNAME lowMSSM<Two_scale>

lowMSSM_susy_scale_constraint<Two_scale>::lowMSSM_susy_scale_constraint()
   : Constraint<Two_scale>()
   , scale(0.)
   , initial_scale_guess(0.)
   , model(0)
   , inputPars()
{
}

lowMSSM_susy_scale_constraint<Two_scale>::lowMSSM_susy_scale_constraint(const lowMSSM_input_parameters& inputPars_)
   : Constraint<Two_scale>()
   , model(0)
   , inputPars(inputPars_)
{
   initialize();
}

lowMSSM_susy_scale_constraint<Two_scale>::~lowMSSM_susy_scale_constraint()
{
}

void lowMSSM_susy_scale_constraint<Two_scale>::apply()
{
   assert(model && "Error: lowMSSM_susy_scale_constraint:"
          " model pointer must not be zero");

   model->calculate_DRbar_parameters();
   update_scale();

   // apply user-defined susy scale constraints
   const auto MuInput = INPUTPARAMETER(MuInput);
   const auto TYdInput = INPUTPARAMETER(TYdInput);
   const auto TYeInput = INPUTPARAMETER(TYeInput);
   const auto TYuInput = INPUTPARAMETER(TYuInput);
   const auto BMuInput = INPUTPARAMETER(BMuInput);
   const auto mq2Input = INPUTPARAMETER(mq2Input);
   const auto ml2Input = INPUTPARAMETER(ml2Input);
   const auto md2Input = INPUTPARAMETER(md2Input);
   const auto mu2Input = INPUTPARAMETER(mu2Input);
   const auto me2Input = INPUTPARAMETER(me2Input);
   const auto MassBInput = INPUTPARAMETER(MassBInput);
   const auto MassWBInput = INPUTPARAMETER(MassWBInput);
   const auto MassGInput = INPUTPARAMETER(MassGInput);

   MODEL->set_Mu(MuInput);
   MODEL->set_TYd(TYdInput);
   MODEL->set_TYe(TYeInput);
   MODEL->set_TYu(TYuInput);
   MODEL->set_BMu(BMuInput);
   MODEL->set_mq2(mq2Input);
   MODEL->set_ml2(ml2Input);
   MODEL->set_md2(md2Input);
   MODEL->set_mu2(mu2Input);
   MODEL->set_me2(me2Input);
   MODEL->set_MassB(MassBInput);
   MODEL->set_MassWB(MassWBInput);
   MODEL->set_MassG(MassGInput);


   // the parameters, which are fixed by the EWSB eqs., will now be
   // defined at this scale (at the EWSB loop level defined in the
   // model)
   model->solve_ewsb();
}

double lowMSSM_susy_scale_constraint<Two_scale>::get_scale() const
{
   return scale;
}

double lowMSSM_susy_scale_constraint<Two_scale>::get_initial_scale_guess() const
{
   return initial_scale_guess;
}

void lowMSSM_susy_scale_constraint<Two_scale>::set_model(Two_scale_model* model_)
{
   model = cast_model<lowMSSM<Two_scale> >(model_);
}

void lowMSSM_susy_scale_constraint<Two_scale>::set_input_parameters(const lowMSSM_input_parameters& inputPars_)
{
   inputPars = inputPars_;
}

void lowMSSM_susy_scale_constraint<Two_scale>::clear()
{
   scale = 0.;
   initial_scale_guess = 0.;
   model = NULL;
}

void lowMSSM_susy_scale_constraint<Two_scale>::initialize()
{
   initial_scale_guess = 1000;

   scale = initial_scale_guess;
}

void lowMSSM_susy_scale_constraint<Two_scale>::update_scale()
{
   const auto MSu = MODELPARAMETER(MSu);

   scale = Sqrt(MSu(0)*MSu(5));


}

} // namespace flexiblesusy
