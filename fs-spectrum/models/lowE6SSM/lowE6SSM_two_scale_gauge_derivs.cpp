#include "lowE6SSM_two_scale_tuning_calculator.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

double lowE6SSM_tuning_calculator::deriv_dbeta_two_loop_g1_dLambdax() const
{

  const double Lambdax = model.get_Lambdax();
  const double g1 = model.get_g1();

  double deriv = -2.4*Power(g1,3)*Lambdax;

  
  return twoLoop * deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_two_loop_g2_dLambdax() const
{

  const double Lambdax = model.get_Lambdax();
  const double g2 = model.get_g2();
  
  double deriv = -4*Power(g2,3)*Lambdax;
  
  
  return twoLoop * deriv;
}

double lowE6SSM_tuning_calculator::deriv_dbeta_two_loop_gN_dLambdax() const
{
  const lowE6SSM_input_parameters inputs = model.get_input();

  const auto QH1p = inputs.QH1p;
  const auto QH2p = inputs.QH2p;
  const auto QSp = inputs.QSp;

  const double Lambdax = model.get_Lambdax();
  const double gN = model.get_gN();
  
  double deriv = -8*Power(gN,3)*Lambdax*(Sqr(QH1p) + Sqr(QH2p) + Sqr(QSp));
  
  
  return twoLoop * deriv;
}

} // namespace flexiblesusy
