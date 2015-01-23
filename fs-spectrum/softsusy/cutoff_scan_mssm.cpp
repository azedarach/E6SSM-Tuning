// ====================================================================
// Scanning code for calculating fine tuning in the MSSM, 
// keeping the parameters fixed and varying the cut-off (i.e.
// input) scale
// ====================================================================

#include "flags.h"
#include "tuning_utils.h"
#include "tuning_numerics.h"
#include "mssm_tuning_utils.h"

#include "softsusy_mycomplex.h"
#include "softsusy_def.h"
#include "softsusy_linalg.h"
#include "softsusy_lowe.h"
#include "softsusy_rge.h"
#include "softsusy_softsusy.h"
#include "softsusy_softpars.h"
#include "softsusy_susy.h"
#include "softsusy_utils.h"
#include "softsusy_numerics.h"

#include "error.hpp"
#include "grid_scanner.hpp"
#include "logger.hpp"
#include "scan.hpp"
#include "wrappers.hpp"

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include <iostream>
#include <iomanip>

using namespace essmsoftsusy;

std::string ToUpper(const std::string & s)
{
   std::string result;
   for (std::size_t index = 0; index < s.length(); ++index) {
      char a = s[index];
      a = std::toupper(a);
      result = result + a;
   }

   return result;
}

Eigen::VectorXd fill_linear_values(std::size_t num_points, double lower, double upper)
{
   Eigen::VectorXd vec;
   vec.resize(num_points);

   double incr = 0.;

   if (num_points > 1)
      incr = (upper - lower) / (num_points - 1.);

   for (std::size_t i = 0; i < num_points; ++i) {
      vec(i) = lower + incr * i;
   }

   return vec;
}

Eigen::VectorXd fill_log_values(std::size_t num_points, double lower, double upper)
{
   Eigen::VectorXd vec;
   vec.resize(num_points);

   double incr = 0.;

   if (lower <= 0. || upper <= 0.) {
      ERROR("Bounds of logarithmic scan must be positive.");
      return vec;
   }

   if (num_points > 1)
      incr = (flexiblesusy::Log(upper) - flexiblesusy::Log(lower)) / (num_points - 1.);

   for (std::size_t i = 0; i < num_points; ++i) {
      vec(i) = std::exp(flexiblesusy::Log(lower) + incr * i);
   }

   return vec;
}

void fill_default_input_pars(DoubleVector & pars)
{
   const bool use_default_pars = false;

   if (use_default_pars) {
      // default values used in scan
      pars(1) = 300.; //< M1
      pars(3) = 1170.0;//2000.; //< M3
      pars(12) = 0.; //< AYd(2,2)
      pars(13) = 0.; //< AYe(2,2)
      pars(31) = 5000.; //< ml(0,0)
      pars(32) = 5000.; //< ml(1,1)
      pars(33) = 5000.; //< ml(2,2)
      pars(34) = 5000.; //< me(0,0)
      pars(35) = 5000.; //< me(1,1)
      pars(36) = 5000.; //< me(2,2)
      pars(41) = 5000.; //< mq(0,0)
      pars(42) = 5000.; //< mq(1,1)
      pars(44) = 5000.; //< mu(0,0)
      pars(45) = 5000.; //< mu(1,1)
      pars(46) = 5000.; //< mu(2,2)
      pars(47) = 5000.; //< md(0,0)
      pars(48) = 5000.; //< md(1,1)
      pars(49) = 5000.; //< md(2,2)
   } else {
      pars(1) = 4.35809794e+02; //< M1
      pars(3) = 1160.0;//2.09996302e+03; //< M3
      pars(12) = -1.92476064e+03; //< AYd(2,2)
      pars(13) = -8.31473226e+01; //< AYe(2,2)
      pars(31) = 2.28409908e+03; //< ml(0,0)
      pars(32) = 2.28406934e+03; //< ml(1,1)
      pars(33) = 2.27516902e+03; //< ml(2,2)
      pars(34) = 2.22578992e+03; //< me(0,0)
      pars(35) = 2.22572861e+03; //< me(1,1)
      pars(36) = 2.20733853e+03; //< me(2,2)
      pars(41) = 2.86004800e+03; //< mq(0,0)
      pars(42) = 2.86003809e+03; //< mq(1,1)
      pars(44) = 2.81426818e+03; //< mu(0,0)
      pars(45) = 2.81425727e+03; //< mu(1,1)
      pars(46) = 1.99865480e+03; //< mu(2,2)
      pars(47) = 2.80871491e+03; //< md(0,0)
      pars(48) = 2.80870549e+03; //< md(1,1)
      pars(49) = 2.79059537e+03; //< md(2,2)
   }
   
}

double maximum_tuning(const DoubleVector& tunings)
{
   double max = -1.0e6;

   for (int i = 1; i <= tunings.displayEnd(); ++i) {
      if (tunings(i) > max) max = tunings(i);
   }

   return max;
}

int main() {

   tryToConvergeHard = true;

   /// Sets up exception handling 
   signal(SIGFPE, FPE_ExceptionHandler);

   /// Sets format of output: 8 decimal places
   outputCharacteristics(8);

   /// SM parameters
   const double poleM_t = 173.2; // top quark pole mass in GeV
   const double mbmb = 4.16; // bottom quark running mass m_b(m_b)^MSbar in GeV
   const double poleM_tau = 1.777;// tau pole mass in GeV
   const double poleM_Z = 91.1876; // Z boson pole mass in GeV
   const double G_F = 1.16637900e-5; // G_F^MSbar in GeV^-2.
   const double alphasmz = 0.1193; // alpha_s(mz)^MSbar (strong coupling constant = g_3^2/(4*PI))
   const double alphaemmz = 127.9568; // alpha_em(mz)^MSbar (electromagnetic coupling = e^2/(4*PI))

   QedQcd oneset; 
   
   oneset.setAlpha(ALPHA, 1.0 / alphaemmz);			 
   oneset.setAlpha(ALPHAS, alphasmz);
   oneset.setMu(poleM_Z); MZ = poleM_Z;
   oneset.setMass(mBottom, mbmb);
   oneset.setPoleMt(poleM_t);
   oneset.setMass(mTau, poleM_tau); oneset.setPoleMtau(poleM_tau);

   QedQcd dataset = oneset;

   oneset.toMz();

   // default number of cut-off scales
   int init_mx_npts = 1; //< for dealing with user input
   std::size_t mx_npts = 1; //< actual value
   
   // default lower bound
   double mx_lower = 20000.0;
   
   // default upper bound
   double mx_upper = 1.0e16;
   
   // default values for input parameters
   double TanBeta_value = 10.0;
   double Mu_value = 0.1; // GeV
   double MAPole_value = 500.0; // GeV
   double AYu22_value = 5000.0; // GeV
   double mq222_value = 4.0e4; // GeV^2
   double mu222_value = 4.0e4; // GeV^2
   double MassWB_value = 1050.0; // GeV
   
   bool has_mx_npts = false;
   bool has_mx_lower = false;
   bool has_mx_upper = false;
   bool has_TanBeta_value = false;
   bool has_Mu_value = false;
   bool has_MAPole_value = false;
   bool has_AYu22_value = false;
   bool has_mq222_value = false;
   bool has_mu222_value = false;
   bool has_MassWB_value = false;

   bool fix_at_msusy = false;

   // read in scan parameters
   try {
      std::string line;
      std::size_t line_num = 0;

      while (getline(cin, line)) {
         line_num++;

         std::istringstream input(line);
         std::string word1, word2;

         input >> word1;

         if (word1.find("#") != 0) {
            // remove anything after the first comment character
            if (word1.find("#") != string::npos) {
               word1 = word1.substr(0, word1.find("#"));
            }

            // all valid options have an = in them, so check for that
            if (word1.find("=") == string::npos) {
               WARNING("Invalid option " + word1 + ": ignoring it.");
            } else {
               // divide into two strings: option and value
               word2 = word1.substr(word1.find("=") + 1, string::npos); //< value
               word1 = word1.substr(0, word1.find("=")); //< option
               word1 = ToUpper(word1);

               std::istringstream kk(word2);

               if (word1 == "MXLL") {
                  kk >> mx_lower;
                  has_mx_lower = true;
               } else if (word1 == "MXUL") {
                  kk >> mx_upper;
                  has_mx_upper = true;
               } else if (word1 == "MXNPTS") {
                  kk >> init_mx_npts;
                  has_mx_npts = true;
               } else if (word1 == "TBVAL") {
                  kk >> TanBeta_value;
                  has_TanBeta_value = true;
               } else if (word1 == "MUVAL") {
                  kk >> Mu_value;
                  has_Mu_value = true;
               } else if (word1 == "MAPOLEVAL") {
                  kk >> MAPole_value;
                  has_MAPole_value = true;
               } else if (word1 == "ATVAL") {
                  kk >> AYu22_value;
                  has_AYu22_value = true;
               } else if (word1 == "MQLSQVAL") {
                  kk >> mq222_value;
                  has_mq222_value = true;
               } else if (word1 == "MURSQVAL") {
                  kk >> mu222_value;
                  has_mu222_value = true;
               } else if (word1 == "M2VAL") {
                  kk >> MassWB_value;
                  has_MassWB_value = true;
               } else if (word1 == "INPUTSCALE") {
                  if (ToUpper(word2) == "MSUSY") {
                     fix_at_msusy = true;
                  } else if (ToUpper(word2) == "MX") {
                     fix_at_msusy = false;
                  } else {
                     WARNING("Unrecognised value '" + word2 + "' for option " + word1 + ": ignoring it");
                  }
               } else {
                  WARNING("Unrecognised option '" + word1 + "': ignoring it");
               }

            }
         } // if (word1.find("#") != 0)
      } // while (getline(cin, line))

      // check inputs are valid
      if (!has_mx_lower) {
         WARNING("MX lower limit not found: using default value");
      }

      if (!has_mx_upper) {
         WARNING("MX upper limit not found: using default value");
      }

      if (!has_mx_npts) {
         WARNING("MX number of points not found: using default value");
      }

      if (!has_TanBeta_value) {
         WARNING("TanBeta value not found: using default value");
      }

      if (!has_Mu_value) {
         WARNING("Mu value not found: using default value");
      }

      if (!has_MAPole_value) {
         WARNING("M_A(Pole) value not found: using default value");
      }

      if (!has_AYu22_value) {
         WARNING("AYu(2,2) value not found: using default value");
      }

      if (!has_mq222_value) {
         WARNING("mq2(2,2) value not found: using default value");
      }

      if (!has_mu222_value) {
         WARNING("mu2(2,2) value not found: using default value");
      }

      if (!has_MassWB_value) {
         WARNING("MassWB value not found: using default value");
      }

      double temp;
      if (mx_lower > mx_upper) {
         temp = mx_lower;
         mx_lower = mx_upper;
         mx_upper = temp;
      }

      // check for valid values of TanBeta
      if (TanBeta_value < 1.0 || TanBeta_value > 500.0) {
         WARNING("Invalid TanBeta value of TanBeta: using default value");
         TanBeta_value = 10.0;
      }

      // check for valid number of points for cut-off scale
      if (init_mx_npts < 1) {
         WARNING("Invalid number of MX points: using default 1");
         mx_npts = 1;
      } else {
         mx_npts = init_mx_npts;
      }

      // get scanned cut-off scale values
      Eigen::VectorXd mx_vals = fill_log_values(mx_npts, mx_lower, mx_upper);

      std::vector<std::size_t> scan_dimensions
         = {mx_npts};

      flexiblesusy::Grid_scanner scan(scan_dimensions);

      // print header line
      std::cout << "# "
                << std::setw(12) << std::left << "TanBeta" << ' '
                << std::setw(12) << std::left << "Mu(MX)" << ' '
                << std::setw(12) << std::left << "B(MX)/GeV" << ' '
                << std::setw(12) << std::left << "AYu22(MX)/GeV" << ' '
                << std::setw(12) << std::left << "mq222(MX)/GeV^2" << ' '
                << std::setw(12) << std::left << "mu222(MX)/GeV^2" << ' '
                << std::setw(12) << std::left << "MassWB(MX)/GeV" << ' '
                << std::setw(12) << std::left << "Mhh(1)/GeV" << ' '
                << std::setw(12) << std::left << "Mhh(2)/GeV" << ' '
                << std::setw(12) << std::left << "mHd2/GeV^2" << ' '
                << std::setw(12) << std::left << "mHu2/GeV^2" << ' '
                << std::setw(12) << std::left << "f1" << ' '
                << std::setw(12) << std::left << "f2" << ' '
                << std::setw(12) << std::left << "MS/Gev" << ' '
                << std::setw(12) << std::left << "MX/GeV" << ' '
                << std::setw(12) << std::left << "MVZ/GeV" << ' '
                << std::setw(12) << std::left << "MSu(1)/GeV" << ' '
                << std::setw(12) << std::left << "MSu(2)/GeV" << ' '
                << std::setw(12) << std::left << "MSb(1)/GeV" << ' '
                << std::setw(12) << std::left << "MSb(2)/GeV" << ' '
                << std::setw(12) << std::left << "MCha(1)/GeV" << ' '
                << std::setw(12) << std::left << "MCha(2)/GeV" << ' '
                << std::setw(12) << std::left << "MChi(1)/GeV" << ' '
                << std::setw(12) << std::left << "MChi(2)/GeV" << ' '
                << std::setw(12) << std::left << "MChi(3)/GeV" << ' '
                << std::setw(12) << std::left << "MChi(4)/GeV" << ' '
                << std::setw(12) << std::left << "Mhh(1)DRbar/GeV" << ' '
                << std::setw(12) << std::left << "Mhh(2)DRbar/GeV" << ' '
                << std::setw(12) << std::left << "MAh(1)DRbar/GeV" << ' '
                << std::setw(12) << std::left << "D(Max)" << ' '
                << std::setw(12) << std::left << "D(MassB)" << ' '
                << std::setw(12) << std::left << "D(MassWB)" << ' '
                << std::setw(12) << std::left << "D(MassG)" << ' '
                << std::setw(12) << std::left << "D(AYu22)" << ' '
                << std::setw(12) << std::left << "D(AYd22)" << ' '
                << std::setw(12) << std::left << "D(AYe22)" << ' '
                << std::setw(12) << std::left << "D(mHd2)" << ' '
                << std::setw(12) << std::left << "D(mHu2)" << ' '
                << std::setw(12) << std::left << "D(Mu)" << ' '
                << std::setw(12) << std::left << "D(B)" << ' '
                << std::setw(12) << std::left << "D(ml002)" << ' '
                << std::setw(12) << std::left << "D(me002)" << ' '
                << std::setw(12) << std::left << "D(mq002)" << ' '
                << std::setw(12) << std::left << "D(mu002)" << ' '
                << std::setw(12) << std::left << "D(md002)" << ' '
                << std::setw(12) << std::left << "D(ml222)" << ' '
                << std::setw(12) << std::left << "D(me222)" << ' '
                << std::setw(12) << std::left << "D(mq222)" << ' '
                << std::setw(12) << std::left << "D(mu222)" << ' '
                << std::setw(12) << std::left << "D(md222)" << ' '
                << std::setw(12) << std::left << "tuning_err" << ' '
                << std::setw(12) << std::left << "error"
                << '\n';

      // DH::note
      std::cerr << "# "
                << std::setw(12) << std::left << "TanBeta" << ' '
                << std::setw(12) << std::left << "Mu(MX)" << ' '
                << std::setw(12) << std::left << "B(MX)/GeV" << ' '
                << std::setw(12) << std::left << "AYu22(MX)/GeV" << ' '
                << std::setw(12) << std::left << "mq222(MX)/GeV^2" << ' '
                << std::setw(12) << std::left << "mu222(MX)/GeV^2" << ' '
                << std::setw(12) << std::left << "MassWB(MX)/GeV" << ' '
                << std::setw(12) << std::left << "Mhh(1)/GeV" << ' '
                << std::setw(12) << std::left << "Mhh(2)/GeV" << ' '
                << std::setw(12) << std::left << "mHd2/GeV^2" << ' '
                << std::setw(12) << std::left << "mHu2/GeV^2" << ' '
                << std::setw(12) << std::left << "MS/Gev" << ' '
                << std::setw(12) << std::left << "MX/GeV" << ' '
                << std::setw(12) << std::left << "Yu(2,2)" << ' '
                << std::setw(12) << std::left << "Yd(2,2)" << ' '
                << std::setw(12) << std::left << "Ye(2,2)" << ' '
                << std::setw(12) << std::left << "g1" << ' '
                << std::setw(12) << std::left << "g2" << ' '
                << std::setw(12) << std::left << "g3" << ' '
                << std::setw(12) << std::left << "Mu/GeV" << ' '
                << std::setw(12) << std::left << "MassB/GeV" << ' '
                << std::setw(12) << std::left << "MassWB/GeV" << ' '
                << std::setw(12) << std::left << "MassG/GeV" << ' '
                << std::setw(12) << std::left << "TYu(2,2)/GeV" << ' '
                << std::setw(12) << std::left << "TYd(2,2)/GeV" << ' '
                << std::setw(12) << std::left << "TYe(2,2)/GeV" << ' '
                << std::setw(12) << std::left << "mq2(2,2)/GeV^2" << ' '
                << std::setw(12) << std::left << "ml2(2,2)/GeV^2" << ' '
                << std::setw(12) << std::left << "me2(2,2)/GeV^2" << ' '
                << std::setw(12) << std::left << "mu2(2,2)/GeV^2" << ' '
                << std::setw(12) << std::left << "md2(2,2)/GeV^2" << ' '
                << std::setw(12) << std::left << "BMu/GeV^2" << ' '
                << std::setw(12) << std::left << "MSu(1)/GeV" << ' '
                << std::setw(12) << std::left << "MSu(2)/GeV" << ' '
                << std::setw(12) << std::left << "tuning_err" << ' '
                << std::setw(12) << std::left << "error"
                << '\n';

      // initialise model and inputs
      const int sgnMu = 0;
      double mgut_guess = 1.0e3; //< default for fixing at M_SUSY
      bool gauge_unification = true; //< actually ignored
      bool ewsb_BC_scale = false;

      void (*boundary_condition)(MssmSoftsusy &, const DoubleVector &)
         = &extendedSugraBcs;

      if (!fix_at_msusy) {
         gauge_unification = false;
      } else {
         ewsb_BC_scale = true;
         gauge_unification = false;
      }

      DoubleVector pars(49);
      
      // note TanBeta is set at MZ
      fill_default_input_pars(pars);

      pars(2) = MassWB_value;
      pars(11) = AYu22_value;
      pars(23) = Mu_value;
      pars(26) = MAPole_value;
      pars(43) = flexiblesusy::Sqrt(mq222_value);
      pars(46) = flexiblesusy::Sqrt(mu222_value);

      std::vector<std::size_t> position;
      while (!scan.has_finished()) {
         position = scan.get_position();
         double mx_value = mx_vals(position.at(0));

         MssmSoftsusy model;
         model.useAlternativeEwsb();
         model.setMuCond(Mu_value);
         model.setSusyMu(Mu_value);
         model.setMaCond(MAPole_value);

         if (!fix_at_msusy) {
            mgut_guess = mx_value;
         }

         model.lowOrg(boundary_condition, mgut_guess, pars, sgnMu, 
                      TanBeta_value, oneset, gauge_unification,
                      ewsb_BC_scale);

         const bool error = model.displayProblem().test();

         // model output is at M_SUSY
         model.runto(model.displayMsusy());

         //model.lesHouchesAccordOutput(cout, "nonUniversal", pars, sgnMu, TanBeta_value, 0.0,  
         //                             1, ewsb_BC_scale);

         // DH::note
         MssmSoftsusy high_scale_model(model);
         high_scale_model.runto(mx_value);
         high_scale_model.calcDrBarPars();
         double B = 0.;
         double f1 = 0.;
         double f2 = 0.;
         DoubleVector tunings(NUMPMSSMPARS);
         double max_tuning = 0.;

         bool tuning_error = false;
         // if no problems, calculate tuning
         if (!error) {
            try {
               MssmSoftsusy tuning_model(model);
               const bool use_mx_equals_ms = true;
               DoubleVector soft_masses_1lp(2);
               
               // set the soft masses correctly for the tuning measures
               bool ewsb_error = MSSM_ImplementEWSBConstraints_SoftMasses(tuning_model, tuning_model.displayMu(), tuning_model.displayMu(),
                                                                          use_mx_equals_ms, soft_masses_1lp); 
               
               tuning_model.setMh1Squared(soft_masses_1lp(1));
               tuning_model.setMh2Squared(soft_masses_1lp(2));
               
               f1 = MSSM_EWSBCondition1(tuning_model);
               f2 = MSSM_EWSBCondition2(tuning_model);

               // calculate the parameter B
               if (fix_at_msusy) {
                  B = model.displayM3Squared() / model.displaySusyMu();
               }
               
               tuning_model.runto(mx_value);
               
               if (!fix_at_msusy) {
                  B = model.displayM3Squared() / model.displaySusyMu();
               }
               
               // calculate the fine tuning
               DoubleVector tuningPars(NUMPMSSMPARS);
               
               tuningPars(1) = tuning_model.displayGaugino(1);
               tuningPars(2) = tuning_model.displayGaugino(2);
               tuningPars(3) = tuning_model.displayGaugino(3);
               tuningPars(4) = tuning_model.displaySoftA(UA, 3, 3);
               tuningPars(5) = tuning_model.displaySoftA(DA, 3, 3);
               tuningPars(6) = tuning_model.displaySoftA(EA, 3, 3);
               tuningPars(7) = tuning_model.displayMh1Squared();
               tuningPars(8) = tuning_model.displayMh2Squared();
               tuningPars(9) = tuning_model.displaySusyMu();
               tuningPars(10) = tuning_model.displayM3Squared()/tuning_model.displaySusyMu();
               tuningPars(11) = tuning_model.displaySoftMassSquared(mLl, 1, 1);
               tuningPars(12) = tuning_model.displaySoftMassSquared(mEr, 1, 1);
               tuningPars(13) = tuning_model.displaySoftMassSquared(mQl, 1, 1);
               tuningPars(14) = tuning_model.displaySoftMassSquared(mUr, 1, 1);
               tuningPars(15) = tuning_model.displaySoftMassSquared(mDr, 1, 1);
               tuningPars(16) = tuning_model.displaySoftMassSquared(mLl, 3, 3);
               tuningPars(17) = tuning_model.displaySoftMassSquared(mEr, 3, 3);
               tuningPars(18) = tuning_model.displaySoftMassSquared(mQl, 3, 3);
               tuningPars(19) = tuning_model.displaySoftMassSquared(mUr, 3, 3);
               tuningPars(20) = tuning_model.displaySoftMassSquared(mDr, 3, 3);
               
               // This is faster...
               //tunings = doCalcMSSMTuningNumerically(tuning_model, tuning_model.displayMsusy(), 
               //                                      mx_value, tuningPars, pMSSMftBCs); 
               tunings = doCalcpMSSMFineTuning(tuning_model, tuning_model.displayMsusy(), ewsb_error, tuning_error, false, 10.);
               
               for (int i = 1; i <= tunings.displayEnd(); ++i) {
                  // Flag errors with negative values
                  if (tunings(i) < 0.) {
                     tuning_error = true;
                     break;
                  }
               }
               
               if (tuning_error || ewsb_error) {
                  tuning_error = true;
               }

               max_tuning = maximum_tuning(tunings);
            }
            catch(const string & a) 
	    { 
               WARNING("Serious numerical problem encountered in fine tuning calculation");
               tuning_error = true; 
	    }
            catch(const char * a) 
	    { 
               WARNING("Serious numerical problem encountered in fine tuning calculation");
               tuning_error = true; 
	    }
            catch(...) 
	    { 
               WARNING("Serious numerical problem encountered in fine tuning calculation");
               tuning_error = true; 
	    }
         }

         // print output
         //std::cout << "MGlu = " << model.displayPhys().mGluino << "\n";
         std::cout << std::setw(12) << std::left << TanBeta_value << ' '
                   << std::setw(12) << std::left << Mu_value << ' '
                   << std::setw(12) << std::left << B << ' '
                   << std::setw(12) << std::left << AYu22_value << ' '
                   << std::setw(12) << std::left << mq222_value << ' '
                   << std::setw(12) << std::left << mu222_value << ' '
                   << std::setw(12) << std::left << MassWB_value << ' '
                   << std::setw(12) << std::left << model.displayPhys().mh0 << ' '
                   << std::setw(12) << std::left << model.displayPhys().mH0 << ' '
                   << std::setw(12) << std::left << model.displayMh1Squared() << ' '
                   << std::setw(12) << std::left << model.displayMh2Squared() << ' '
                   << std::setw(12) << std::left << f1 << ' '
                   << std::setw(12) << std::left << f2 << ' '
                   << std::setw(12) << std::left << model.displayMu() << ' '
                   << std::setw(12) << std::left << mx_value << ' '
                   << std::setw(12) << std::left << model.displayMzRun() << ' '
                   << std::setw(12) << std::left << model.displayDrBarPars().mu(1,3) << ' '
                   << std::setw(12) << std::left << model.displayDrBarPars().mu(2,3) << ' '
                   << std::setw(12) << std::left << model.displayDrBarPars().md(1,3) << ' '
                   << std::setw(12) << std::left << model.displayDrBarPars().md(2,3) << ' '
                   << std::setw(12) << std::left << model.displayDrBarPars().mch(1) << ' '
                   << std::setw(12) << std::left << model.displayDrBarPars().mch(2) << ' '
                   << std::setw(12) << std::left << model.displayDrBarPars().mneut(1) << ' '
                   << std::setw(12) << std::left << model.displayDrBarPars().mneut(2) << ' '
                   << std::setw(12) << std::left << model.displayDrBarPars().mneut(3) << ' '
                   << std::setw(12) << std::left << model.displayDrBarPars().mneut(4) << ' '
                   << std::setw(12) << std::left << model.displayDrBarPars().mh0 << ' '
                   << std::setw(12) << std::left << model.displayDrBarPars().mH0 << ' '
                   << std::setw(12) << std::left << model.displayDrBarPars().mA0 << ' '
                   << std::setw(12) << std::left << max_tuning << ' '
                   << std::setw(12) << std::left << tunings(1) << ' '
                   << std::setw(12) << std::left << tunings(2) << ' '
                   << std::setw(12) << std::left << tunings(3) << ' '
                   << std::setw(12) << std::left << tunings(4) << ' '
                   << std::setw(12) << std::left << tunings(5) << ' '
                   << std::setw(12) << std::left << tunings(6) << ' '
                   << std::setw(12) << std::left << tunings(7) << ' '
                   << std::setw(12) << std::left << tunings(8) << ' '
                   << std::setw(12) << std::left << tunings(9) << ' '
                   << std::setw(12) << std::left << tunings(10) << ' '
                   << std::setw(12) << std::left << tunings(11) << ' '
                   << std::setw(12) << std::left << tunings(12) << ' '
                   << std::setw(12) << std::left << tunings(13) << ' '
                   << std::setw(12) << std::left << tunings(14) << ' '
                   << std::setw(12) << std::left << tunings(15) << ' '
                   << std::setw(12) << std::left << tunings(16) << ' '
                   << std::setw(12) << std::left << tunings(17) << ' '
                   << std::setw(12) << std::left << tunings(18) << ' '
                   << std::setw(12) << std::left << tunings(19) << ' '
                   << std::setw(12) << std::left << tunings(20) << ' '
                   << std::setw(12) << std::left << tuning_error << ' '
                   << std::setw(12) << std::left << error;

         if (error || tuning_error) {
            std::cout << "# " << model.displayProblem();
            if (tuning_error) {
               std::cout << ", tuning error\n"; 
            } else {
               std::cout << '\n';
            }
         } else {
            std::cout << '\n';
         }
         

         // DH::note
         std::cerr << std::setw(12) << std::left << TanBeta_value << ' '
                   << std::setw(12) << std::left << Mu_value << ' '
                   << std::setw(12) << std::left << B << ' '
                   << std::setw(12) << std::left << AYu22_value << ' '
                   << std::setw(12) << std::left << mq222_value << ' '
                   << std::setw(12) << std::left << mu222_value << ' '
                   << std::setw(12) << std::left << MassWB_value << ' '
                   << std::setw(12) << std::left << model.displayPhys().mh0 << ' '
                   << std::setw(12) << std::left << model.displayPhys().mH0 << ' '
                   << std::setw(12) << std::left << high_scale_model.displayMh1Squared() << ' '
                   << std::setw(12) << std::left << high_scale_model.displayMh2Squared() << ' '
                   << std::setw(12) << std::left << model.displayMu() << ' '
                   << std::setw(12) << std::left << mx_value << ' '
                   << std::setw(12) << std::left << high_scale_model.displayYukawaElement(YU, 3, 3) << ' '
                   << std::setw(12) << std::left << high_scale_model.displayYukawaElement(YD, 3, 3) << ' '
                   << std::setw(12) << std::left << high_scale_model.displayYukawaElement(YE, 3, 3) << ' '
                   << std::setw(12) << std::left << high_scale_model.displayGaugeCoupling(1) << ' '
                   << std::setw(12) << std::left << high_scale_model.displayGaugeCoupling(2) << ' '
                   << std::setw(12) << std::left << high_scale_model.displayGaugeCoupling(3) << ' '
                   << std::setw(12) << std::left << high_scale_model.displaySusyMu() << ' '
                   << std::setw(12) << std::left << high_scale_model.displayGaugino(1) << ' '
                   << std::setw(12) << std::left << high_scale_model.displayGaugino(2) << ' '
                   << std::setw(12) << std::left << high_scale_model.displayGaugino(3) << ' '
                   << std::setw(12) << std::left << high_scale_model.displayTrilinear(UA, 3, 3) << ' '
                   << std::setw(12) << std::left << high_scale_model.displayTrilinear(DA, 3, 3) << ' '
                   << std::setw(12) << std::left << high_scale_model.displayTrilinear(EA, 3, 3) << ' '
                   << std::setw(12) << std::left << high_scale_model.displaySoftMassSquared(mQl, 3, 3) << ' '
                   << std::setw(12) << std::left << high_scale_model.displaySoftMassSquared(mLl, 3, 3) << ' '
                   << std::setw(12) << std::left << high_scale_model.displaySoftMassSquared(mEr, 3, 3) << ' '
                   << std::setw(12) << std::left << high_scale_model.displaySoftMassSquared(mUr, 3, 3) << ' '
                   << std::setw(12) << std::left << high_scale_model.displaySoftMassSquared(mDr, 3, 3) << ' '
                   << std::setw(12) << std::left << high_scale_model.displayM3Squared() << ' '
                   << std::setw(12) << std::left << high_scale_model.displayDrBarPars().mu(1,3) << ' '
                   << std::setw(12) << std::left << high_scale_model.displayDrBarPars().mu(2,3) << ' '
                   << std::setw(12) << std::left << tuning_error << ' '
                   << std::setw(12) << std::left << error;

         if (error || tuning_error) {
            std::cerr << "# " << model.displayProblem();
            if (tuning_error) {
               std::cerr << ", tuning error\n"; 
            } else {
               std::cerr << '\n';
            }
         } else {
            std::cerr << '\n';
         }

         scan.step_forward();
      } // while (!scan.has_finished())
   }
   catch(const string & a) { cout << a; return -1; }
   catch(const char * a) { cout << a; return -1; }
   catch(...) { cout << "Unknown type of exception caught.\n"; return -1; }
   
   return 0;
}
