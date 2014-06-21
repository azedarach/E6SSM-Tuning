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

// File generated at Sun 15 Jun 2014 19:16:44

#ifndef genericE6SSM_TWO_SCALE_CONVERGENCE_TESTER_H
#define genericE6SSM_TWO_SCALE_CONVERGENCE_TESTER_H

#include "genericE6SSM_convergence_tester.hpp"
#include "genericE6SSM_two_scale_model.hpp"
#include "two_scale_convergence_tester_drbar.hpp"

namespace flexiblesusy {

class Two_scale;

template<>
class genericE6SSM_convergence_tester<Two_scale> : public Convergence_tester_DRbar<genericE6SSM<Two_scale> > {
public:
   genericE6SSM_convergence_tester(genericE6SSM<Two_scale>*, double);
   virtual ~genericE6SSM_convergence_tester();

protected:
   virtual double max_rel_diff() const;
};

} // namespace flexiblesusy

#endif
