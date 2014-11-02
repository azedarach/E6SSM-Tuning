#include "lowE6SSM_two_scale_tuning_calculator.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

double lowE6SSM_tuning_calculator::leading_log_coefficient_gN() const
{
   const lowE6SSM_input_parameters inputs = model.get_input();

   const auto QQp = inputs.QQp;
   const auto QLp = inputs.QLp;
   const auto QH1p = inputs.QH1p;
   const auto QH2p = inputs.QH2p;
   const auto Qdp = inputs.Qdp;
   const auto Qup = inputs.Qup;
   const auto Qep = inputs.Qep;
   const auto QSp = inputs.QSp;
   const auto QDxp = inputs.QDxp;
   const auto QDxbarp = inputs.QDxbarp;
   const auto QHpp = inputs.QHpp;
   const auto QHpbarp = inputs.QHpbarp;
   
   const double gN = model.get_gN();

   double coeff = 3*Power(gN,5)*Sqr(9*Power(Qdp,2) + 9*Power(QDxbarp,2) + 9*
      Power(QDxp,2) + 3*Power(Qep,2) + 6*Power(QH1p,2) + 6*Power(QH2p,2) + 2*Power
      (QHpbarp,2) + 2*Power(QHpp,2) + 6*Power(QLp,2) + 18*Power(QQp,2) + 3*Power(
      QSp,2) + 9*Power(Qup,2));


   return twoLoop*coeff;
}

} // namespace flexiblesusy
