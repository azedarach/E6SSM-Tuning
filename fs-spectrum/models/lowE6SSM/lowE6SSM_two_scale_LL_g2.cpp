#include "lowE6SSM_two_scale_tuning_calculator.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

double lowE6SSM_tuning_calculator::leading_log_coefficient_g2() const
{

   const double g2 = model.get_g2();

   double coeff = 48*Power(g2,5);


   return twoLoop*coeff;
}

} // namespace flexiblesusy
