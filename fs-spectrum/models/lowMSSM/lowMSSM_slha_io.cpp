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

// File generated at Sun 24 Aug 2014 16:16:18

#include "lowMSSM_slha_io.hpp"
#include "lowMSSM_input_parameters.hpp"
#include "logger.hpp"
#include "wrappers.hpp"
#include "numerics.hpp"
#include "spectrum_generator_settings.hpp"
#include "lowe.h"
#include "config.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <boost/bind.hpp>

using namespace softsusy;

namespace flexiblesusy {

char const * const lowMSSM_slha_io::drbar_blocks[NUMBER_OF_DRBAR_BLOCKS] =
   { "gauge", "Yu", "Yd", "Ye", "Te", "Td", "Tu", "HMIX", "MSQ2", "MSE2",
   "MSL2", "MSU2", "MSD2", "MSOFT" }
;

lowMSSM_slha_io::lowMSSM_slha_io()
   : slha_io()
{
}

void lowMSSM_slha_io::clear()
{
   slha_io.clear();
}

/**
 * Stores the EXTPAR input parameters in the SLHA object.
 *
 * @param input struct of input parameters
 */
void lowMSSM_slha_io::set_extpar(const lowMSSM_input_parameters& input)
{

}

/**
 * Stores the MINPAR input parameters in the SLHA object.
 *
 * @param input struct of input parameters
 */
void lowMSSM_slha_io::set_minpar(const lowMSSM_input_parameters& input)
{
   std::ostringstream minpar;

   minpar << "Block MINPAR\n";
   minpar << FORMAT_ELEMENT(3, input.TanBeta, "TanBeta");
   slha_io.set_block(minpar);

}

/**
 * Stores the SMINPUTS input parameters in the SLHA object.
 *
 * @param qedqcd class of Standard Model parameters
 */
void lowMSSM_slha_io::set_sminputs(const softsusy::QedQcd& qedqcd)
{
   slha_io.set_sminputs(qedqcd);
}

/**
 * Stores the spectrum generator information in the SPINFO block in
 * the SLHA object.
 *
 * @param problems struct with parameter point problems
 */
void lowMSSM_slha_io::set_spinfo(const Problems<lowMSSM_info::NUMBER_OF_PARTICLES>& problems)
{
   std::ostringstream spinfo;
   spinfo << "Block SPINFO\n"
          << FORMAT_SPINFO(1, PKGNAME)
          << FORMAT_SPINFO(2, FLEXIBLESUSY_VERSION);

   if (problems.have_serious_problem()) {
      std::ostringstream serious_problems;
      problems.print(serious_problems);
      spinfo << FORMAT_SPINFO(4, serious_problems.str());
   }

   slha_io.set_block(spinfo, SLHA_io::front);
}

/**
 * Stores the particle masses in the SLHA object.
 *
 * @param physical struct of physical parameters
 *
 * @param write_sm_masses flag to indicate if Standard Model
 *    particle masses should be written as well
 */
void lowMSSM_slha_io::set_mass(const lowMSSM_physical& physical,
                                   bool write_sm_masses)
{
   std::ostringstream mass;

   mass << "Block MASS\n"
      << FORMAT_MASS(1000021, LOCALPHYSICAL(MGlu), "Glu")
      << FORMAT_MASS(1000024, LOCALPHYSICAL(MCha(0)), "Cha_1")
      << FORMAT_MASS(1000037, LOCALPHYSICAL(MCha(1)), "Cha_2")
      << FORMAT_MASS(25, LOCALPHYSICAL(Mhh(0)), "hh_1")
      << FORMAT_MASS(35, LOCALPHYSICAL(Mhh(1)), "hh_2")
      << FORMAT_MASS(37, LOCALPHYSICAL(MHpm(1)), "Hpm_2")
      << FORMAT_MASS(36, LOCALPHYSICAL(MAh(1)), "Ah_2")
      << FORMAT_MASS(1000012, LOCALPHYSICAL(MSv(0)), "Sv_1")
      << FORMAT_MASS(1000014, LOCALPHYSICAL(MSv(1)), "Sv_2")
      << FORMAT_MASS(1000016, LOCALPHYSICAL(MSv(2)), "Sv_3")
      << FORMAT_MASS(1000022, LOCALPHYSICAL(MChi(0)), "Chi_1")
      << FORMAT_MASS(1000023, LOCALPHYSICAL(MChi(1)), "Chi_2")
      << FORMAT_MASS(1000025, LOCALPHYSICAL(MChi(2)), "Chi_3")
      << FORMAT_MASS(1000035, LOCALPHYSICAL(MChi(3)), "Chi_4")
      << FORMAT_MASS(1000001, LOCALPHYSICAL(MSd(0)), "Sd_1")
      << FORMAT_MASS(1000003, LOCALPHYSICAL(MSd(1)), "Sd_2")
      << FORMAT_MASS(1000005, LOCALPHYSICAL(MSd(2)), "Sd_3")
      << FORMAT_MASS(2000001, LOCALPHYSICAL(MSd(3)), "Sd_4")
      << FORMAT_MASS(2000003, LOCALPHYSICAL(MSd(4)), "Sd_5")
      << FORMAT_MASS(2000005, LOCALPHYSICAL(MSd(5)), "Sd_6")
      << FORMAT_MASS(1000011, LOCALPHYSICAL(MSe(0)), "Se_1")
      << FORMAT_MASS(1000013, LOCALPHYSICAL(MSe(1)), "Se_2")
      << FORMAT_MASS(1000015, LOCALPHYSICAL(MSe(2)), "Se_3")
      << FORMAT_MASS(2000011, LOCALPHYSICAL(MSe(3)), "Se_4")
      << FORMAT_MASS(2000013, LOCALPHYSICAL(MSe(4)), "Se_5")
      << FORMAT_MASS(2000015, LOCALPHYSICAL(MSe(5)), "Se_6")
      << FORMAT_MASS(1000002, LOCALPHYSICAL(MSu(0)), "Su_1")
      << FORMAT_MASS(1000004, LOCALPHYSICAL(MSu(1)), "Su_2")
      << FORMAT_MASS(1000006, LOCALPHYSICAL(MSu(2)), "Su_3")
      << FORMAT_MASS(2000002, LOCALPHYSICAL(MSu(3)), "Su_4")
      << FORMAT_MASS(2000004, LOCALPHYSICAL(MSu(4)), "Su_5")
      << FORMAT_MASS(2000006, LOCALPHYSICAL(MSu(5)), "Su_6")
   ;

   if (write_sm_masses) {
      mass
         << FORMAT_MASS(21, LOCALPHYSICAL(MVG), "VG")
         << FORMAT_MASS(12, LOCALPHYSICAL(MFv(0)), "Fv_1")
         << FORMAT_MASS(14, LOCALPHYSICAL(MFv(1)), "Fv_2")
         << FORMAT_MASS(16, LOCALPHYSICAL(MFv(2)), "Fv_3")
         << FORMAT_MASS(22, LOCALPHYSICAL(MVP), "VP")
         << FORMAT_MASS(23, LOCALPHYSICAL(MVZ), "VZ")
         << FORMAT_MASS(11, LOCALPHYSICAL(MFe(0)), "Fe_1")
         << FORMAT_MASS(13, LOCALPHYSICAL(MFe(1)), "Fe_2")
         << FORMAT_MASS(15, LOCALPHYSICAL(MFe(2)), "Fe_3")
         << FORMAT_MASS(1, LOCALPHYSICAL(MFd(0)), "Fd_1")
         << FORMAT_MASS(3, LOCALPHYSICAL(MFd(1)), "Fd_2")
         << FORMAT_MASS(5, LOCALPHYSICAL(MFd(2)), "Fd_3")
         << FORMAT_MASS(2, LOCALPHYSICAL(MFu(0)), "Fu_1")
         << FORMAT_MASS(4, LOCALPHYSICAL(MFu(1)), "Fu_2")
         << FORMAT_MASS(6, LOCALPHYSICAL(MFu(2)), "Fu_3")
         << FORMAT_MASS(24, LOCALPHYSICAL(MVWm), "VWm")
      ;
   }

   slha_io.set_block(mass);

}

/**
 * Stores the mixing matrices in the SLHA object.
 *
 * @param physical struct of physical parameters
 *
 * @param write_sm_mixing_matrics flag to indicate if Standard Model
 *    particle mixing matrices should be written as well
 */
void lowMSSM_slha_io::set_mixing_matrices(const lowMSSM_physical& physical,
                                              bool write_sm_mixing_matrics)
{
   slha_io.set_block("UMIX", LOCALPHYSICAL(UM), "UM");
   slha_io.set_block("VMIX", LOCALPHYSICAL(UP), "UP");
   slha_io.set_block("PSEUDOSCALARMIX", LOCALPHYSICAL(ZA), "ZA");
   slha_io.set_block("DSQMIX", LOCALPHYSICAL(ZD), "ZD");
   slha_io.set_block("SELMIX", LOCALPHYSICAL(ZE), "ZE");
   slha_io.set_block("SCALARMIX", LOCALPHYSICAL(ZH), "ZH");
   slha_io.set_block("NMIX", LOCALPHYSICAL(ZN), "ZN");
   slha_io.set_block("CHARGEMIX", LOCALPHYSICAL(ZP), "ZP");
   slha_io.set_block("USQMIX", LOCALPHYSICAL(ZU), "ZU");
   slha_io.set_block("SNUMIX", LOCALPHYSICAL(ZV), "ZV");

   if (write_sm_mixing_matrics) {
      slha_io.set_block("UELMIX", LOCALPHYSICAL(ZEL), "ZEL");
      slha_io.set_block("UERMIX", LOCALPHYSICAL(ZER), "ZER");
      slha_io.set_block("UDLMIX", LOCALPHYSICAL(ZDL), "ZDL");
      slha_io.set_block("UDRMIX", LOCALPHYSICAL(ZDR), "ZDR");
      slha_io.set_block("UULMIX", LOCALPHYSICAL(ZUL), "ZUL");
      slha_io.set_block("UURMIX", LOCALPHYSICAL(ZUR), "ZUR");
   }

}

/**
 * Write SLHA object to file.
 *
 * @param file_name file name
 */
void lowMSSM_slha_io::write_to_file(const std::string& file_name)
{
   slha_io.write_to_file(file_name);
}

/**
 * Read (DR-bar) model parameter input scale from EXTPAR entry 0
 */
double lowMSSM_slha_io::get_input_scale() const
{
   return slha_io.get_extpar().input_scale;
}

/**
 * Read (DR-bar) model parameter output scale from MODSEL entry 12
 */
double lowMSSM_slha_io::get_parameter_output_scale() const
{
   return slha_io.get_modsel().parameter_output_scale;
}

/**
 * Read SLHA object from file
 *
 * @param file_name file name
 */
void lowMSSM_slha_io::read_from_file(const std::string& file_name)
{
   slha_io.read_from_file(file_name);
   slha_io.read_modsel();
   slha_io.read_extpar();
}

/**
 * Fill struct of model input parameters from SLHA object (MINPAR and
 * EXTPAR blocks)
 *
 * @param input struct of model input parameters
 */
void lowMSSM_slha_io::fill(lowMSSM_input_parameters& input) const
{
   SLHA_io::Tuple_processor minpar_processor
      = boost::bind(&lowMSSM_slha_io::fill_minpar_tuple, boost::ref(input), _1, _2);
   SLHA_io::Tuple_processor extpar_processor
      = boost::bind(&lowMSSM_slha_io::fill_extpar_tuple, boost::ref(input), _1, _2);

   slha_io.read_block("MINPAR", minpar_processor);
   slha_io.read_block("EXTPAR", extpar_processor);

   slha_io.read_block("TeIN", input.TYeInput);
   slha_io.read_block("TdIN", input.TYdInput);
   slha_io.read_block("TuIN", input.TYuInput);
   input.MuInput = slha_io.read_entry("HMIXIN", 1);
   input.BMuInput = slha_io.read_entry("HMIXIN", 101);
   slha_io.read_block("MSQ2IN", input.mq2Input);
   slha_io.read_block("MSE2IN", input.me2Input);
   slha_io.read_block("MSL2IN", input.ml2Input);
   slha_io.read_block("MSU2IN", input.mu2Input);
   slha_io.read_block("MSD2IN", input.md2Input);
   input.MassBInput = slha_io.read_entry("MSOFTIN", 1);
   input.MassWBInput = slha_io.read_entry("MSOFTIN", 2);
   input.MassGInput = slha_io.read_entry("MSOFTIN", 3);

}

/**
 * Fill struct of spectrum generator settings from SLHA object
 * (FlexibleSUSY block)
 *
 * @param settings struct of spectrum generator settings
 */
void lowMSSM_slha_io::fill(Spectrum_generator_settings& settings) const
{
   SLHA_io::Tuple_processor flexiblesusy_processor
      = boost::bind(&lowMSSM_slha_io::fill_flexiblesusy_tuple, boost::ref(settings), _1, _2);

   slha_io.read_block("FlexibleSUSY", flexiblesusy_processor);
}

void lowMSSM_slha_io::fill_minpar_tuple(lowMSSM_input_parameters& input,
                                                int key, double value)
{
   switch (key) {
   case 3: input.TanBeta = value; break;
   default: WARNING("Unrecognized key: " << key); break;
   }

}

void lowMSSM_slha_io::fill_extpar_tuple(lowMSSM_input_parameters& input,
                                                int key, double value)
{
   // key 0 is the model parameter input scale, which is read in
   // slha_io.{hpp,cpp}
   if (key == 0)
      return;

   switch (key) {
   default: WARNING("Unrecognized key: " << key); break;
   }

}

void lowMSSM_slha_io::fill_flexiblesusy_tuple(Spectrum_generator_settings& settings,
                                                  int key, double value)
{
   if (0 <= key && key < static_cast<int>(Spectrum_generator_settings::NUMBER_OF_OPTIONS)) {
      settings.set((Spectrum_generator_settings::Settings)key, value);
   } else {
      WARNING("Unrecognized key in block FlexibleSUSY: " << key);
   }
}

/**
 * Reads the renormalization scales from all DR-bar parameter blocks.
 * If blocks with different scales are found the last scale is
 * returned and a warning is printed.
 *
 * @return common renormalization scale
 */
double lowMSSM_slha_io::read_scale() const
{
   double scale = 0.;

   for (unsigned i = 0; i < NUMBER_OF_DRBAR_BLOCKS; i++) {
      const double block_scale = slha_io.read_scale(drbar_blocks[i]);
      if (!is_zero(block_scale)) {
         if (!is_zero(scale) && !is_equal(scale, block_scale))
            WARNING("DR-bar parameters defined at different scales");
         scale = block_scale;
      }
   }

   return scale;
}

/**
 * Convert masses and mixing matrices to SLHA convention: Fermion
 * mixing matrices are always real and fermion masses are allowed to
 * be negative.
 *
 * @param physical struct of physical parameters to convert
 */
void lowMSSM_slha_io::convert_to_slha_convention(lowMSSM_physical& physical)
{
   SLHA_io::convert_symmetric_fermion_mixings_to_slha(LOCALPHYSICAL(MChi), LOCALPHYSICAL(ZN));

}

} // namespace flexiblesusy
