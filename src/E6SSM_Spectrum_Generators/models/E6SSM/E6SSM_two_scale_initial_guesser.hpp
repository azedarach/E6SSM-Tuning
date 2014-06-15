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

// File generated at Sun 15 Jun 2014 19:16:28

#ifndef E6SSM_TWO_SCALE_INITIAL_GUESSER_H
#define E6SSM_TWO_SCALE_INITIAL_GUESSER_H

#include "E6SSM_initial_guesser.hpp"
#include "E6SSM_input_parameters.hpp"
#include "E6SSM_two_scale_low_scale_constraint.hpp"
#include "E6SSM_two_scale_susy_scale_constraint.hpp"
#include "E6SSM_two_scale_high_scale_constraint.hpp"
#include "two_scale_initial_guesser.hpp"

#include <sstream>

namespace flexiblesusy {

template <class T>
class E6SSM;

class Two_scale;

template<>
class E6SSM_initial_guesser<Two_scale> : public Initial_guesser<Two_scale> {
public:
   E6SSM_initial_guesser(E6SSM<Two_scale>*,
                               const E6SSM_input_parameters&,
                               const QedQcd&,
                               const E6SSM_low_scale_constraint<Two_scale>&,
                               const E6SSM_susy_scale_constraint<Two_scale>&,
                               const E6SSM_high_scale_constraint<Two_scale>&);
   virtual ~E6SSM_initial_guesser();
   virtual void guess();

   void set_running_precision(double p) { running_precision = p; }

private:
   E6SSM<Two_scale>* model;
   E6SSM_input_parameters input_pars;
   QedQcd oneset;
   double mu_guess;
   double mc_guess;
   double mt_guess;
   double md_guess;
   double ms_guess;
   double mb_guess;
   double me_guess;
   double mm_guess;
   double mtau_guess;
   double running_precision;
   E6SSM_low_scale_constraint<Two_scale> low_constraint;
   E6SSM_susy_scale_constraint<Two_scale> susy_constraint;
   E6SSM_high_scale_constraint<Two_scale> high_constraint;

   void guess_susy_parameters();
   void guess_soft_parameters();
   void calculate_DRbar_yukawa_couplings();
   void calculate_Yu_DRbar();
   void calculate_Yd_DRbar();
   void calculate_Ye_DRbar();
};

} // namespace flexiblesusy

#endif
