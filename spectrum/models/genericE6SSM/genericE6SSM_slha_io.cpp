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

// File generated at Sun 15 Jun 2014 19:16:45

#include "genericE6SSM_slha_io.hpp"
#include "genericE6SSM_input_parameters.hpp"
#include "logger.hpp"
#include "wrappers.hpp"
#include "spectrum_generator_settings.hpp"
#include "lowe.h"
#include "config.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <boost/bind.hpp>

using namespace softsusy;

namespace flexiblesusy {

genericE6SSM_slha_io::genericE6SSM_slha_io()
   : slha_io()
{
}

void genericE6SSM_slha_io::clear()
{
   slha_io.clear();
}

void genericE6SSM_slha_io::set_extpar(const genericE6SSM_input_parameters& input)
{
   std::ostringstream extpar;

   extpar << "Block EXTPAR\n";
   extpar << FORMAT_ELEMENT(61, input.LambdaInput, "LambdaInput");
   extpar << FORMAT_ELEMENT(62, input.KappaInput, "KappaInput");
   extpar << FORMAT_ELEMENT(63, input.muPrimeInput, "muPrimeInput");
   extpar << FORMAT_ELEMENT(64, input.BmuPrimeInput, "BmuPrimeInput");
   extpar << FORMAT_ELEMENT(65, input.vSInput, "vSInput");
   extpar << FORMAT_ELEMENT(66, input.Lambda12Input, "Lambda12Input");
   extpar << FORMAT_ELEMENT(200, input.QQp, "QQp");
   extpar << FORMAT_ELEMENT(201, input.QLp, "QLp");
   extpar << FORMAT_ELEMENT(202, input.QH1p, "QH1p");
   extpar << FORMAT_ELEMENT(203, input.QH2p, "QH2p");
   extpar << FORMAT_ELEMENT(204, input.Qdp, "Qdp");
   extpar << FORMAT_ELEMENT(205, input.Qup, "Qup");
   extpar << FORMAT_ELEMENT(206, input.Qep, "Qep");
   extpar << FORMAT_ELEMENT(207, input.QSp, "QSp");
   extpar << FORMAT_ELEMENT(208, input.QDxp, "QDxp");
   extpar << FORMAT_ELEMENT(209, input.QDxbarp, "QDxbarp");
   extpar << FORMAT_ELEMENT(210, input.QHpp, "QHpp");
   extpar << FORMAT_ELEMENT(211, input.QHpbarp, "QHpbarp");
   slha_io.set_block(extpar);

}

void genericE6SSM_slha_io::set_minpar(const genericE6SSM_input_parameters& input)
{
   std::ostringstream minpar;

   minpar << "Block MINPAR\n";
   minpar << FORMAT_ELEMENT(1, input.m0, "m0");
   minpar << FORMAT_ELEMENT(2, input.m12, "m12");
   minpar << FORMAT_ELEMENT(3, input.TanBeta, "TanBeta");
   minpar << FORMAT_ELEMENT(5, input.Azero, "Azero");
   slha_io.set_block(minpar);

}

void genericE6SSM_slha_io::set_sminputs(const softsusy::QedQcd& qedqcd)
{
   slha_io.set_sminputs(qedqcd);
}

void genericE6SSM_slha_io::set_spinfo(const Problems<genericE6SSM_info::NUMBER_OF_PARTICLES>& problems)
{
   std::ostringstream spinfo;
   spinfo << "# FlexibleSUSY " FLEXIBLESUSY_VERSION " SLHA compliant output\n"
             "# P. Athron, Jae-hyeon Park, D. Stöckinger, A. Voigt\n"
             "Block SPINFO\n"
          << FORMAT_SPINFO(1, PKGNAME)
          << FORMAT_SPINFO(2, FLEXIBLESUSY_VERSION);

   if (problems.have_serious_problem()) {
      std::ostringstream serious_problems;
      problems.print(serious_problems);
      spinfo << FORMAT_SPINFO(4, serious_problems.str());
   }

   slha_io.set_block(spinfo, SLHA_io::front);
}

void genericE6SSM_slha_io::write_to_file(const std::string& file_name)
{
   slha_io.write_to_file(file_name);
}

double genericE6SSM_slha_io::get_input_scale() const
{
   return slha_io.get_extpar().input_scale;
}

double genericE6SSM_slha_io::get_parameter_output_scale() const
{
   return slha_io.get_modsel().parameter_output_scale;
}

void genericE6SSM_slha_io::read_from_file(const std::string& file_name)
{
   slha_io.read_from_file(file_name);
   slha_io.read_modsel();
   slha_io.read_extpar();
}

void genericE6SSM_slha_io::fill(genericE6SSM_input_parameters& input) const
{
   SLHA_io::Tuple_processor minpar_processor
      = boost::bind(&genericE6SSM_slha_io::fill_minpar_tuple, boost::ref(input), _1, _2);
   SLHA_io::Tuple_processor extpar_processor
      = boost::bind(&genericE6SSM_slha_io::fill_extpar_tuple, boost::ref(input), _1, _2);

   slha_io.read_block("MINPAR", minpar_processor);
   slha_io.read_block("EXTPAR", extpar_processor);


}

void genericE6SSM_slha_io::fill(Spectrum_generator_settings& settings) const
{
   SLHA_io::Tuple_processor flexiblesusy_processor
      = boost::bind(&genericE6SSM_slha_io::fill_flexiblesusy_tuple, boost::ref(settings), _1, _2);

   slha_io.read_block("FlexibleSUSY", flexiblesusy_processor);
}

void genericE6SSM_slha_io::fill_minpar_tuple(genericE6SSM_input_parameters& input,
                                                int key, double value)
{
   switch (key) {
   case 1: input.m0 = value; break;
   case 2: input.m12 = value; break;
   case 3: input.TanBeta = value; break;
   case 5: input.Azero = value; break;
   default: WARNING("Unrecognized key: " << key); break;
   }

}

void genericE6SSM_slha_io::fill_extpar_tuple(genericE6SSM_input_parameters& input,
                                                int key, double value)
{
   // key 0 is the model parameter input scale, which is read in
   // slha_io.{hpp,cpp}
   if (key == 0)
      return;

   switch (key) {
   case 61: input.LambdaInput = value; break;
   case 62: input.KappaInput = value; break;
   case 63: input.muPrimeInput = value; break;
   case 64: input.BmuPrimeInput = value; break;
   case 65: input.vSInput = value; break;
   case 66: input.Lambda12Input = value; break;
   case 200: input.QQp = value; break;
   case 201: input.QLp = value; break;
   case 202: input.QH1p = value; break;
   case 203: input.QH2p = value; break;
   case 204: input.Qdp = value; break;
   case 205: input.Qup = value; break;
   case 206: input.Qep = value; break;
   case 207: input.QSp = value; break;
   case 208: input.QDxp = value; break;
   case 209: input.QDxbarp = value; break;
   case 210: input.QHpp = value; break;
   case 211: input.QHpbarp = value; break;
   default: WARNING("Unrecognized key: " << key); break;
   }

}

void genericE6SSM_slha_io::fill_flexiblesusy_tuple(Spectrum_generator_settings& settings,
                                                  int key, double value)
{
   if (0 <= key && key < static_cast<int>(Spectrum_generator_settings::NUMBER_OF_OPTIONS)) {
      settings.set((Spectrum_generator_settings::Settings)key, value);
   } else {
      WARNING("Unrecognized key in block FlexibleSUSY: " << key);
   }
}

} // namespace flexiblesusy
