// ====================================================================
// Additional constraint for fixing parameters at an input scale
// other than that defined by g1 == g2 or M_{SUSY}
// ====================================================================

#ifndef lowE6SSM_TWO_SCALE_INPUT_SCALE_CONSTRAINT_H
#define lowE6SSM_TWO_SCALE_INPUT_SCALE_CONSTRAINT_H

#include "lowE6SSM_input_scale_constraint.hpp"
#include "lowE6SSM_input_parameters.hpp"
#include "two_scale_constraint.hpp"

namespace flexiblesusy {

   template <class T>
   class lowE6SSM;

   class Two_scale;

   template<>
   class lowE6SSM_input_scale_constraint<Two_scale> : public Constraint<Two_scale> {
   public:
      lowE6SSM_input_scale_constraint();
      lowE6SSM_input_scale_constraint(const lowE6SSM_input_parameters&);
      virtual ~lowE6SSM_input_scale_constraint();
      virtual void apply();
      virtual double get_scale() const;
      virtual void set_model(Two_scale_model*);

      void clear();
      double get_initial_scale_guess() const;
      void initialize();
      void set_input_parameters(const lowE6SSM_input_parameters&);
      void set_scale(double); ///< fixed input scale

   private:
      double initial_scale_guess;
      double fixed_scale; ///< fixed input scale
      lowE6SSM<Two_scale>* model;
      lowE6SSM_input_parameters inputPars;

   };

} // namespace flexiblesusy

#endif
