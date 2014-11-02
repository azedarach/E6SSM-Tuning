#include "lowE6SSM_two_scale_tuning_calculator.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

double lowE6SSM_tuning_calculator::leading_log_coefficient_g1() const
{

   const double g1 = model.get_g1();

   double coeff = 276.48*Power(g1,5);


   return twoLoop*coeff;
}

} // namespace flexiblesusy
