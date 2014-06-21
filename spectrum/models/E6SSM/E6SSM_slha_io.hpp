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

// File generated at Sun 15 Jun 2014 19:16:26

#ifndef E6SSM_SLHA_IO_H
#define E6SSM_SLHA_IO_H

#include "E6SSM_two_scale_model.hpp"
#include "E6SSM_info.hpp"
#include "slha_io.hpp"

#include <Eigen/Core>
#include <string>
#include <utility>

#define PHYSICAL(p) model.get_physical().p
#define MODELPARAMETER(p) model.get_##p()
#define DEFINE_PARAMETER(p)                                            \
   typename std::remove_const<typename std::remove_reference<decltype(MODELPARAMETER(p))>::type>::type p;

namespace flexiblesusy {

struct E6SSM_input_parameters;
class Spectrum_generator_settings;

class E6SSM_slha_io {
public:
   E6SSM_slha_io();
   ~E6SSM_slha_io() {}

   void clear();

   void fill(QedQcd& qedqcd) const { slha_io.fill(qedqcd); }
   void fill(E6SSM_input_parameters&) const;
   template <class T> void fill(E6SSM<T>&) const;
   void fill(Spectrum_generator_settings&) const;
   double get_input_scale() const;
   double get_parameter_output_scale() const;
   void read_from_file(const std::string&);
   void set_extpar(const E6SSM_input_parameters&);
   void set_minpar(const E6SSM_input_parameters&);
   void set_sminputs(const softsusy::QedQcd&);
   template <class T> void set_spectrum(const E6SSM<T>&);
   void set_spinfo(const Problems<E6SSM_info::NUMBER_OF_PARTICLES>&);
   void write_to_file(const std::string&);
   void write_to_stream(std::ostream& ostr = std::cout) { slha_io.write_to_stream(ostr); }

   static void fill_minpar_tuple(E6SSM_input_parameters&, int, double);
   static void fill_extpar_tuple(E6SSM_input_parameters&, int, double);
   static void fill_flexiblesusy_tuple(Spectrum_generator_settings&, int, double);

private:
   SLHA_io slha_io; ///< SLHA io class

   template <class T> void set_mass(const E6SSM<T>&);
   template <class T> void set_mixing_matrices(const E6SSM<T>&);
   template <class T> void set_model_parameters(const E6SSM<T>&);
};

/**
 * Reads all model parameters from a SLHA output file.
 */
template <class T>
void E6SSM_slha_io::fill(E6SSM<T>& model) const
{
   model.set_g1(slha_io.read_entry("gauge", 1) * 1.2909944487358056);
   model.set_g2(slha_io.read_entry("gauge", 2));
   model.set_g3(slha_io.read_entry("gauge", 3));
   model.set_gN(slha_io.read_entry("gauge", 4));
   {
      DEFINE_PARAMETER(Yu);
      slha_io.read_block("Yu", Yu);
      model.set_Yu(Yu);
   }
   {
      DEFINE_PARAMETER(Yd);
      slha_io.read_block("Yd", Yd);
      model.set_Yd(Yd);
   }
   {
      DEFINE_PARAMETER(Ye);
      slha_io.read_block("Ye", Ye);
      model.set_Ye(Ye);
   }
   {
      DEFINE_PARAMETER(TYe);
      slha_io.read_block("Te", TYe);
      model.set_TYe(TYe);
   }
   {
      DEFINE_PARAMETER(TYd);
      slha_io.read_block("Td", TYd);
      model.set_TYd(TYd);
   }
   {
      DEFINE_PARAMETER(TYu);
      slha_io.read_block("Tu", TYu);
      model.set_TYu(TYu);
   }
   {
      DEFINE_PARAMETER(mq2);
      slha_io.read_block("MSQ2", mq2);
      model.set_mq2(mq2);
   }
   {
      DEFINE_PARAMETER(me2);
      slha_io.read_block("MSE2", me2);
      model.set_me2(me2);
   }
   {
      DEFINE_PARAMETER(ml2);
      slha_io.read_block("MSL2", ml2);
      model.set_ml2(ml2);
   }
   {
      DEFINE_PARAMETER(mu2);
      slha_io.read_block("MSU2", mu2);
      model.set_mu2(mu2);
   }
   {
      DEFINE_PARAMETER(md2);
      slha_io.read_block("MSD2", md2);
      model.set_md2(md2);
   }
   model.set_mHd2(slha_io.read_entry("MSOFT", 21));
   model.set_mHu2(slha_io.read_entry("MSOFT", 22));
   {
      DEFINE_PARAMETER(mH1I2);
      slha_io.read_block("mHdInert2", mH1I2);
      model.set_mH1I2(mH1I2);
   }
   {
      DEFINE_PARAMETER(mH2I2);
      slha_io.read_block("mHuInert2", mH2I2);
      model.set_mH2I2(mH2I2);
   }
   {
      DEFINE_PARAMETER(mDx2);
      slha_io.read_block("mX2", mDx2);
      model.set_mDx2(mDx2);
   }
   {
      DEFINE_PARAMETER(mDxbar2);
      slha_io.read_block("mXBar2", mDxbar2);
      model.set_mDxbar2(mDxbar2);
   }
   model.set_ms2(slha_io.read_entry("MSOFT", 23));
   {
      DEFINE_PARAMETER(msI2);
      slha_io.read_block("msInert2", msI2);
      model.set_msI2(msI2);
   }
   model.set_mHp2(slha_io.read_entry("MSOFT", 24));
   model.set_mHpbar2(slha_io.read_entry("MSOFT", 25));
   model.set_MassB(slha_io.read_entry("MSOFT", 1));
   model.set_MassWB(slha_io.read_entry("MSOFT", 2));
   model.set_MassG(slha_io.read_entry("MSOFT", 3));
   model.set_MassBp(slha_io.read_entry("MSOFT", 4));
   model.set_vd(slha_io.read_entry("HMIX", 102));
   model.set_vu(slha_io.read_entry("HMIX", 103));
   model.set_vs(slha_io.read_entry("ESIXRUN", 11));
   {
      DEFINE_PARAMETER(Kappa);
      slha_io.read_block("ESIXKAPPA", Kappa);
      model.set_Kappa(Kappa);
   }
   {
      DEFINE_PARAMETER(TKappa);
      slha_io.read_block("ESIXTKAPPA", TKappa);
      model.set_TKappa(TKappa);
   }
   model.set_Lambdax(slha_io.read_entry("ESIXRUN", 1));
   model.set_TLambdax(slha_io.read_entry("ESIXRUN", 2));
   {
      DEFINE_PARAMETER(Lambda12);
      slha_io.read_block("ESIXLAMBDA", Lambda12);
      model.set_Lambda12(Lambda12);
   }
   {
      DEFINE_PARAMETER(TLambda12);
      slha_io.read_block("ESIXTLAMBDA", TLambda12);
      model.set_TLambda12(TLambda12);
   }
   model.set_MuPr(slha_io.read_entry("ESIXRUN", 0));
   model.set_BMuPr(slha_io.read_entry("ESIXRUN", 101));

}

template <class T>
void E6SSM_slha_io::set_mass(const E6SSM<T>& model)
{
   const auto MVG = PHYSICAL(MVG);
   const auto MGlu = PHYSICAL(MGlu);
   const auto MFv = PHYSICAL(MFv);
   const auto MChaP = PHYSICAL(MChaP);
   const auto MVP = PHYSICAL(MVP);
   const auto MVZ = PHYSICAL(MVZ);
   const auto MVZp = PHYSICAL(MVZp);
   const auto MSd = PHYSICAL(MSd);
   const auto MSv = PHYSICAL(MSv);
   const auto MSu = PHYSICAL(MSu);
   const auto MSe = PHYSICAL(MSe);
   const auto MSDX = PHYSICAL(MSDX);
   const auto Mhh = PHYSICAL(Mhh);
   const auto MAh = PHYSICAL(MAh);
   const auto MHpm = PHYSICAL(MHpm);
   const auto MChi = PHYSICAL(MChi);
   const auto MCha = PHYSICAL(MCha);
   const auto MFe = PHYSICAL(MFe);
   const auto MFd = PHYSICAL(MFd);
   const auto MFu = PHYSICAL(MFu);
   const auto MFDX = PHYSICAL(MFDX);
   const auto MSHI0 = PHYSICAL(MSHI0);
   const auto MSHIp = PHYSICAL(MSHIp);
   const auto MChaI = PHYSICAL(MChaI);
   const auto MChiI = PHYSICAL(MChiI);
   const auto MSSI0 = PHYSICAL(MSSI0);
   const auto MFSI = PHYSICAL(MFSI);
   const auto MSHp0 = PHYSICAL(MSHp0);
   const auto MSHpp = PHYSICAL(MSHpp);
   const auto MChiP = PHYSICAL(MChiP);
   const auto MVWm = PHYSICAL(MVWm);

   std::ostringstream mass;

   mass << "Block MASS\n"
      << FORMAT_MASS(1000021, MGlu, "Glu")
      << FORMAT_MASS(1000091, MChaP, "ChaP")
      << FORMAT_MASS(31, MVZp, "VZp")
      << FORMAT_MASS(1000089, MFSI(0), "FSI_1")
      << FORMAT_MASS(1000090, MFSI(1), "FSI_2")
      << FORMAT_MASS(1000092, MChiP(0), "ChiP_1")
      << FORMAT_MASS(1000094, MChiP(1), "ChiP_2")
      << FORMAT_MASS(1000024, MCha(0), "Cha_1")
      << FORMAT_MASS(1000037, MCha(1), "Cha_2")
      << FORMAT_MASS(37, MHpm(1), "Hpm_2")
      << FORMAT_MASS(92, MSHp0(0), "SHp0_1")
      << FORMAT_MASS(94, MSHp0(1), "SHp0_2")
      << FORMAT_MASS(91, MSHpp(0), "SHpp_1")
      << FORMAT_MASS(93, MSHpp(1), "SHpp_2")
      << FORMAT_MASS(89, MSSI0(0), "SSI0_1")
      << FORMAT_MASS(90, MSSI0(1), "SSI0_2")
      << FORMAT_MASS(1000085, MChaI(0), "ChaI_1")
      << FORMAT_MASS(1000086, MChaI(1), "ChaI_2")
      << FORMAT_MASS(25, Mhh(0), "hh_1")
      << FORMAT_MASS(35, Mhh(1), "hh_2")
      << FORMAT_MASS(45, Mhh(2), "hh_3")
      << FORMAT_MASS(1000012, MSv(0), "Sv_1")
      << FORMAT_MASS(1000014, MSv(1), "Sv_2")
      << FORMAT_MASS(1000016, MSv(2), "Sv_3")
      << FORMAT_MASS(36, MAh(2), "Ah_3")
      << FORMAT_MASS(51, MFDX(0), "FDX_1")
      << FORMAT_MASS(52, MFDX(1), "FDX_2")
      << FORMAT_MASS(53, MFDX(2), "FDX_3")
      << FORMAT_MASS(1000081, MChiI(0), "ChiI_1")
      << FORMAT_MASS(1000082, MChiI(1), "ChiI_2")
      << FORMAT_MASS(1000083, MChiI(2), "ChiI_3")
      << FORMAT_MASS(1000084, MChiI(3), "ChiI_4")
      << FORMAT_MASS(82, MSHI0(0), "SHI0_1")
      << FORMAT_MASS(86, MSHI0(1), "SHI0_2")
      << FORMAT_MASS(84, MSHI0(2), "SHI0_3")
      << FORMAT_MASS(88, MSHI0(3), "SHI0_4")
      << FORMAT_MASS(81, MSHIp(0), "SHIp_1")
      << FORMAT_MASS(85, MSHIp(1), "SHIp_2")
      << FORMAT_MASS(83, MSHIp(2), "SHIp_3")
      << FORMAT_MASS(87, MSHIp(3), "SHIp_4")
      << FORMAT_MASS(1000022, MChi(0), "Chi_1")
      << FORMAT_MASS(1000023, MChi(1), "Chi_2")
      << FORMAT_MASS(1000025, MChi(2), "Chi_3")
      << FORMAT_MASS(1000035, MChi(3), "Chi_4")
      << FORMAT_MASS(1000045, MChi(4), "Chi_5")
      << FORMAT_MASS(1000055, MChi(5), "Chi_6")
      << FORMAT_MASS(1000001, MSd(0), "Sd_1")
      << FORMAT_MASS(1000003, MSd(1), "Sd_2")
      << FORMAT_MASS(1000005, MSd(2), "Sd_3")
      << FORMAT_MASS(2000001, MSd(3), "Sd_4")
      << FORMAT_MASS(2000003, MSd(4), "Sd_5")
      << FORMAT_MASS(2000005, MSd(5), "Sd_6")
      << FORMAT_MASS(1000011, MSe(0), "Se_1")
      << FORMAT_MASS(1000013, MSe(1), "Se_2")
      << FORMAT_MASS(1000015, MSe(2), "Se_3")
      << FORMAT_MASS(2000011, MSe(3), "Se_4")
      << FORMAT_MASS(2000013, MSe(4), "Se_5")
      << FORMAT_MASS(2000015, MSe(5), "Se_6")
      << FORMAT_MASS(1000002, MSu(0), "Su_1")
      << FORMAT_MASS(1000004, MSu(1), "Su_2")
      << FORMAT_MASS(1000006, MSu(2), "Su_3")
      << FORMAT_MASS(2000002, MSu(3), "Su_4")
      << FORMAT_MASS(2000004, MSu(4), "Su_5")
      << FORMAT_MASS(2000006, MSu(5), "Su_6")
      << FORMAT_MASS(1000051, MSDX(0), "SDX_1")
      << FORMAT_MASS(2000051, MSDX(1), "SDX_2")
      << FORMAT_MASS(1000052, MSDX(2), "SDX_3")
      << FORMAT_MASS(2000052, MSDX(3), "SDX_4")
      << FORMAT_MASS(1000053, MSDX(4), "SDX_5")
      << FORMAT_MASS(2000053, MSDX(5), "SDX_6")
   ;

   if (model.do_calculate_sm_pole_masses()) {
      mass
         << FORMAT_MASS(21, MVG, "VG")
         << FORMAT_MASS(12, MFv(0), "Fv_1")
         << FORMAT_MASS(14, MFv(1), "Fv_2")
         << FORMAT_MASS(16, MFv(2), "Fv_3")
         << FORMAT_MASS(22, MVP, "VP")
         << FORMAT_MASS(23, MVZ, "VZ")
         << FORMAT_MASS(11, MFe(0), "Fe_1")
         << FORMAT_MASS(13, MFe(1), "Fe_2")
         << FORMAT_MASS(15, MFe(2), "Fe_3")
         << FORMAT_MASS(1, MFd(0), "Fd_1")
         << FORMAT_MASS(3, MFd(1), "Fd_2")
         << FORMAT_MASS(5, MFd(2), "Fd_3")
         << FORMAT_MASS(2, MFu(0), "Fu_1")
         << FORMAT_MASS(4, MFu(1), "Fu_2")
         << FORMAT_MASS(6, MFu(2), "Fu_3")
         << FORMAT_MASS(24, MVWm, "VWm")
      ;
   }

   slha_io.set_block(mass);

}

template <class T>
void E6SSM_slha_io::set_mixing_matrices(const E6SSM<T>& model)
{
   slha_io.set_block("INHMIX", PHYSICAL(UHI0), "UHI0");
   slha_io.set_block("ICHMIX", PHYSICAL(UHIp), "UHIp");
   slha_io.set_block("UHNPMIX", PHYSICAL(UHp0), "UHp0");
   slha_io.set_block("UHPPMIX", PHYSICAL(UHpp), "UHpp");
   slha_io.set_block("UMIX", PHYSICAL(UM), "UM");
   slha_io.set_block("VMIX", PHYSICAL(UP), "UP");
   slha_io.set_block("NMAMIX", PHYSICAL(ZA), "ZA");
   slha_io.set_block("DSQMIX", PHYSICAL(ZD), "ZD");
   slha_io.set_block("ESIXZDX", PHYSICAL(ZDX), "ZDX");
   slha_io.set_block("ESIXZXL", PHYSICAL(ZDXL), "ZDXL");
   slha_io.set_block("ESIXZXR", PHYSICAL(ZDXR), "ZDXR");
   slha_io.set_block("SELMIX", PHYSICAL(ZE), "ZE");
   slha_io.set_block("ESIXZSI", PHYSICAL(ZFSI), "ZFSI");
   slha_io.set_block("NMHMIX", PHYSICAL(ZH), "ZH");
   slha_io.set_block("ESIXZMI", PHYSICAL(ZMI), "ZMI");
   slha_io.set_block("NMNMIX", PHYSICAL(ZN), "ZN");
   slha_io.set_block("ESIXZNI", PHYSICAL(ZNI), "ZNI");
   slha_io.set_block("ZNPMIX", PHYSICAL(ZNp), "ZNp");
   slha_io.set_block("CHARGEMIX", PHYSICAL(ZP), "ZP");
   slha_io.set_block("ESIXZPI", PHYSICAL(ZPI), "ZPI");
   slha_io.set_block("ZSSI", PHYSICAL(ZSSI), "ZSSI");
   slha_io.set_block("USQMIX", PHYSICAL(ZU), "ZU");
   slha_io.set_block("SNUMIX", PHYSICAL(ZV), "ZV");

   if (model.do_calculate_sm_pole_masses()) {
      slha_io.set_block("UELMIX", PHYSICAL(ZEL), "ZEL");
      slha_io.set_block("UERMIX", PHYSICAL(ZER), "ZER");
      slha_io.set_block("UDLMIX", PHYSICAL(ZDL), "ZDL");
      slha_io.set_block("UDRMIX", PHYSICAL(ZDR), "ZDR");
      slha_io.set_block("UULMIX", PHYSICAL(ZUL), "ZUL");
      slha_io.set_block("UURMIX", PHYSICAL(ZUR), "ZUR");
   }

}

template <class T>
void E6SSM_slha_io::set_model_parameters(const E6SSM<T>& model)
{
   {
      std::ostringstream block;
      block << "Block gauge Q= " << FORMAT_NUMBER(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(1, (MODELPARAMETER(g1) * 0.7745966692414834), "gY")
            << FORMAT_ELEMENT(2, MODELPARAMETER(g2), "g2")
            << FORMAT_ELEMENT(3, MODELPARAMETER(g3), "g3")
            << FORMAT_ELEMENT(4, MODELPARAMETER(gN), "gN");
      slha_io.set_block(block);
   }
   slha_io.set_block("Yu", MODELPARAMETER(Yu), "Yu", model.get_scale());
   slha_io.set_block("Yd", MODELPARAMETER(Yd), "Yd", model.get_scale());
   slha_io.set_block("Ye", MODELPARAMETER(Ye), "Ye", model.get_scale());
   slha_io.set_block("Te", MODELPARAMETER(TYe), "TYe", model.get_scale());
   slha_io.set_block("Td", MODELPARAMETER(TYd), "TYd", model.get_scale());
   slha_io.set_block("Tu", MODELPARAMETER(TYu), "TYu", model.get_scale());
   slha_io.set_block("MSQ2", MODELPARAMETER(mq2), "mq2", model.get_scale());
   slha_io.set_block("MSE2", MODELPARAMETER(me2), "me2", model.get_scale());
   slha_io.set_block("MSL2", MODELPARAMETER(ml2), "ml2", model.get_scale());
   slha_io.set_block("MSU2", MODELPARAMETER(mu2), "mu2", model.get_scale());
   slha_io.set_block("MSD2", MODELPARAMETER(md2), "md2", model.get_scale());
   {
      std::ostringstream block;
      block << "Block MSOFT Q= " << FORMAT_NUMBER(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(21, MODELPARAMETER(mHd2), "mHd2")
            << FORMAT_ELEMENT(22, MODELPARAMETER(mHu2), "mHu2")
            << FORMAT_ELEMENT(23, MODELPARAMETER(ms2), "ms2")
            << FORMAT_ELEMENT(24, MODELPARAMETER(mHp2), "mHp2")
            << FORMAT_ELEMENT(25, MODELPARAMETER(mHpbar2), "mHpbar2")
            << FORMAT_ELEMENT(1, MODELPARAMETER(MassB), "MassB")
            << FORMAT_ELEMENT(2, MODELPARAMETER(MassWB), "MassWB")
            << FORMAT_ELEMENT(3, MODELPARAMETER(MassG), "MassG")
            << FORMAT_ELEMENT(4, MODELPARAMETER(MassBp), "MassBp");
      slha_io.set_block(block);
   }
   slha_io.set_block("mHdInert2", MODELPARAMETER(mH1I2), "mH1I2", model.get_scale());
   slha_io.set_block("mHuInert2", MODELPARAMETER(mH2I2), "mH2I2", model.get_scale());
   slha_io.set_block("mX2", MODELPARAMETER(mDx2), "mDx2", model.get_scale());
   slha_io.set_block("mXBar2", MODELPARAMETER(mDxbar2), "mDxbar2", model.get_scale());
   slha_io.set_block("msInert2", MODELPARAMETER(msI2), "msI2", model.get_scale());
   {
      std::ostringstream block;
      block << "Block HMIX Q= " << FORMAT_NUMBER(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(102, MODELPARAMETER(vd), "vd")
            << FORMAT_ELEMENT(103, MODELPARAMETER(vu), "vu");
      slha_io.set_block(block);
   }
   {
      std::ostringstream block;
      block << "Block ESIXRUN Q= " << FORMAT_NUMBER(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(11, MODELPARAMETER(vs), "vs")
            << FORMAT_ELEMENT(1, MODELPARAMETER(Lambdax), "Lambdax")
            << FORMAT_ELEMENT(2, MODELPARAMETER(TLambdax), "TLambdax")
            << FORMAT_ELEMENT(0, MODELPARAMETER(MuPr), "MuPr")
            << FORMAT_ELEMENT(101, MODELPARAMETER(BMuPr), "BMuPr");
      slha_io.set_block(block);
   }
   slha_io.set_block("ESIXKAPPA", MODELPARAMETER(Kappa), "Kappa", model.get_scale());
   slha_io.set_block("ESIXTKAPPA", MODELPARAMETER(TKappa), "TKappa", model.get_scale());
   slha_io.set_block("ESIXLAMBDA", MODELPARAMETER(Lambda12), "Lambda12", model.get_scale());
   slha_io.set_block("ESIXTLAMBDA", MODELPARAMETER(TLambda12), "TLambda12", model.get_scale());

}

template <class T>
void E6SSM_slha_io::set_spectrum(const E6SSM<T>& model)
{
   set_model_parameters(model);
   set_mass(model);
   set_mixing_matrices(model);
}

} // namespace flexiblesusy

#endif
