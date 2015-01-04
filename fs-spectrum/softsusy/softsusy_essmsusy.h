
/** \file susy.h
   - Project:     SOFTSUSY 
   - Author:      Ben Allanach 
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   - Webpage:     http://allanach.home.cern.ch/allanach/softsusy.html
   - Description: All display and set commands apart from (DoubleVector &) 
                involved with running to different energy scales have been
                explicitly named.

   $Log: susy.h,v $
   Revision 1.3  2005/11/09 14:12:25  allanach
   Updated for SOFTSUSY 2.0.1 - cleaned everything up etc

   Revision 1.1.1.1  2004/11/19 16:18:31  allanach


   Revision 1.11  2004/01/15 13:54:55  allanach
   New heaer style implemented

   Revision 1.10  2003/10/24 16:09:04  allanach
   Implemented running Higgs DRbar vev

   Revision 1.8  2003/07/28 12:11:37  allanach
   More error trapping, and rearranging rpvsoftsusy to use correct Higgs VEV
   (which is sometimes called at MZ)

   Revision 1.7  2003/05/20 15:19:40  allanach
   doxygen comment style implemented

   Revision 1.6  2003/02/21 17:59:36  allanach
   Added drbar parameter class and calculation, starting to move to DRbar
   parameters in the 1-loop corrections

   Revision 1.5  2003/02/21 13:02:07  allanach
   Changed headings to new conventions

   Revision 1.4  2002/11/19 16:59:30  allanach
   Added routine to get quark Yukawa diagonalising matrices

   Revision 1.3  2002/04/12 06:24:50  allanach
   Code maintenance - returning a subobject made simpler

   Revision 1.2  2001/07/18 14:42:51  allanach
   Added proper header info
*/

#ifndef SOFTSUSY_ESSMSUSY_H
#define SOFTSUSY_ESSMSUSY_H

#include <iostream>
#include <cmath>
#include <fstream>
#include "softsusy_lowe.h"
#include "softsusy_utils.h"
#include "softsusy_susy.h"
#include "softsusy_linalg.h"
using namespace std;
//changed to add on extra Essm susy parameters
//const static int numSusyESSMPars = 33;

const static int numSusyESSMPars = 46;

/// For accessing up, down and charged lepton Yukawa matrices respectively
//typedef enum {YU=1, YD, YE} yukawa;

/// Contains data needed in beta function calculation to make it faster
/*struct sBrevity {
  /// dt=\f$Y_D^T\f$, ut=\f$Y_U^T\f$, et=\f$Y_E^T\f$, u2=\f$Y_U Y_U^T\f$,
  /// d2=\f$Y_D Y_D^T\f$, e2=\f$Y_E Y_E^T\f$, 
  /// u2t=\f$Y_U^T Y_U\f$, d2t=\f$Y_D^T Y_D\f$, e2t=\f$Y_E^T Y_E\f$
  DoubleMatrix dt, ut, et, u2, d2, e2, u2t, e2t, d2t;
  /// gsq=\f$g_i^2\f$, g3=\f$g_i^3\f$, g4=\f$g_i^4\f$
  DoubleVector  gsq, g3, g4;
  /// uuT=\f$Tr(Y_U^T Y_U)\f$, ddt=\f$Tr(Y_D^T Y_D)\f$, eet=\f$Tr(Y_E^T Y_E)\f$
  double uuT, ddT, eeT;
  /// d1=\f$Y_D\f$, u1=\f$Y_U\f$, e1=\f$Y_E\f$
  DoubleMatrix u1, d1, e1;
  
  /// Constructor fills sruct with zeroes by default
  sBrevity();
  /// Constructor sets struct to be equal to another
  sBrevity(const sBrevity &);
  /// Sets struct to be equal to another
  const sBrevity & operator=(const sBrevity &);
  
  /// Calculates all the data given Yukawa matrices yu, yd and ye and gauge
  /// couplings g
  void calculate(const DoubleMatrix & yu, const DoubleMatrix & yd, const
		 DoubleMatrix & ye, const DoubleVector  & g);
};

inline sBrevity::sBrevity()
  : dt(3, 3), ut(3, 3), et(3, 3), u2(3, 3), 
     d2(3, 3), e2(3, 3), u2t(3, 3), 
    e2t(3, 3), d2t(3, 3), gsq(1, 3), g3(1, 3), g4(1, 3), uuT(0.0),
     ddT(0.0), eeT(0.0), u1(3, 3), d1(3, 3), e1(3, 3)
{}

inline sBrevity::sBrevity(const sBrevity &s)
  : dt(s.dt), ut(s.ut), et(s.ut), u2(s.u2), d2(s.d2),
     e2(s.e2), u2t(s.u2t),
     e2t(s.e2t), d2t(s.d2t), gsq(s.gsq), g3(s.g3), g4(s.g4),
     uuT(s.uuT), ddT(s.ddT), eeT(s.eeT), u1(s.u1), d1(s.d1), e1(s.e1) 
{}

*/

// Dylan:: I have added in the GUT normalised U(1)_Y and U(1)_N charges
// of the third generation Higgses and singlets as named
// constants for convenience.
const double QY_H1=-0.5*sqrt(3.0/5.0), QY_H2=0.5*sqrt(3.0/5.0), QY_S=0;
const double QN_H1=-3.0/sqrt(40), QN_H2=-2.0/sqrt(40), QN_S=5.0/sqrt(40);

/// Contains all supersymmetric RPCMSSM parameters and RGEs for interpolating
/// them (In MSSM domain) 
//Peter:: In MSSM version, this is an inherits from RGE, I am making this inherit from MssmSusy instead
class EssmSusy: public MssmSusy //< Dylan:: changed RGE to MssmSusy
{
 private:
  DoubleMatrix u, d, e; ///< Yukawa matrices for ups, downs and leptons
  DoubleVector g; ///< Gauge couplings
  /// Bilinear Higgs superpotential parameter and ratio of Higgs VEVs,
  /// \f$ v_1/v_2 \f$ 
  double smu, tanb, hVev; 
  
  //Peter:: New ESSM susy parameters
  
  DoubleVector lambda, kappa, MN_SUSY; 
  
  double h_N, mu_0, gdash_1, g_11;
  
  
 public:
  EssmSusy(); ///< Constructor fills object with zeroes by default
  /// Constructor sets object to be equal to another
  EssmSusy(const EssmSusy &); 
  /// Constructor given Yukawa matrices u,d,e, gauge couplings v, mu
  /// parameter=m, tan beta=tb, renormalisation scale MU, number of loops in
  /// RG evolution l and thresholds parameter t
  //  MssmSusy(const DoubleMatrix & u, const DoubleMatrix & d, const
  //DoubleMatrix & e, const DoubleVector & v, double m,
  //   double tb, double MU, int l, int t, double h);
  
  //Peter:: New constructor which takes new ESSM parameters


// u = up type Yukawa couplings
// d = down type Yukawa couplings
// e = lepton Yukawa couplings
// v = SM/MSSM gauge couplings
// m = MSSM mu parameter
// tb = tan(beta)
// MU = renormalisation scale
// l = number of loops
// t = threshold accuracy
// h = Higgs vev (no class member for S vev)
// lam = ESSM SH_1iH_2i Yukawa couplings (only lambda(3) appears in tuning measures)
// kap = ESSM SD_iDbar_i couplings (don't appear in my tuning measures)
// MN_SUSY = right-handed neutrino masses (not needed)
// h_N = Yukawa coupling to right-handed neutrino (not needed)
// mu_0 = ?
// gdash_1 = ESSM U(1)_N gauge coupling g_1'
// g_11 = ESSM gauge coupling mixing g_11
  EssmSusy(const DoubleMatrix & u, const DoubleMatrix & d, const
	   DoubleMatrix & e, const DoubleVector & v, double m,
	   double tb, double MU, int l, int t, double h, DoubleVector lambda, DoubleVector kappa, DoubleVector MN_SUSY, double h_N, double mu_0, double gdash_1, double g_11);
  
  virtual ~EssmSusy() {}; ///< Default destructor
  
  /// sets object to be equal to another
  const EssmSusy & operator=(const EssmSusy & s);
  /// sets object to be equal to another
  void setSusy(const EssmSusy &s);
  
  /// Sets DRbar running Higgs vev
  void setHvev(double h);
  /// Copies Yukawa matrices and gauge couplings from s only
  void setSomePars(const EssmSusy & s);
  /// Sets one element of a Yukawa matrix
  void setYukawaElement(yukawa, int, int, double);
  /// Sets whole Yukawa matrix 
  void setYukawaMatrix(yukawa, const DoubleMatrix &);
  /// Set a single gauge coupling
  void setGaugeCoupling(int, double);
  /// Set all gauge couplings
  void setAllGauge(const DoubleVector  &);
  /// Sets superpotential mu parameter
  void setSusyMu(double);
  /// Sets all RGE parameters to elements of vector
  void set(const DoubleVector &);
  /// Sets tan beta
  void setTanb(double);
  
  //Peter:: setting new ESSM parameters.
  void seth_N(double);
  void setgash_1(double);
  void setg_11(double);
  void setmu_0(double);
  void setkappa(DoubleVector);
  void setlambda(DoubleVector);
  void setMN_SUSY(DoubleVector);
  /// Returns DRbar running Higgs vev
  double displayHvev() const;
  /// Returns whole object as a const
  inline EssmSusy displaySusy() const;
  /// Returns a single Yukawa matrix element
  double displayYukawaElement(yukawa, int, int) const;
  /// Returns a whole Yukawa matrix
  DoubleMatrix displayYukawaMatrix(yukawa) const;
  /// Returns a single gauge coupling
  double displayGaugeCoupling(int) const;
  /// Returns all gauge couplings
  DoubleVector  displayGauge() const;
  /// Returns superpotential mu parameter
  double displaySusyMu() const;
  /// Returns all parameters as elements of a vector
  const  DoubleVector display() const; // Dylan:: Made this return a const reference
  /// Returns tan beta
  double displayTanb() const;
  double displayh_N() const;
  double displaymu_0() const;
  double displaygdash_1() const;
  double displayg_11() const;
  double displayMN_SUSY(int i) const;
  DoubleVector displaykappa() const;
  DoubleVector displaylambda() const;

  // Dylan:: here I have added some methods to get the value of s
  // given M_Z'.
  double displaySvev(double) const;

  // Dylan:: I am also adding methods to display the effective
  // U(1)_N charges at the current scale (i.e. using the current
  // value of g_11 and g_dash_1), but note the effective charges
  // are not stored as variables.
  double displayQH1tilde() const;
  double displayQH2tilde() const;
  double displayQStilde() const;
  DoubleVector displayQtilde() const;

/// outputs object QedQcd & r valid at 1 GeV from SUSY data at mt, from
  /// diagonal elements of Yukawa couplings and Higgs VEV vev. 
  void getMasses(QedQcd & r, double vev) const;
  /// This turns diagonal Yukawa couplings at MZ into CKM mixed ones
  /// Takes diagonal quark Yukawa matrices and mixes them up
  /// according to the CKM matrix assuming:
  /// mix=2, all mixing is in down sector
  /// mix=1, all mixing is in up sector
  void quarkMixing(const DoubleMatrix & CKM, int mix);
  /// Sets diagonal Yukawa couplings according to data in QedQcd input and
  /// Higgs VEV parameter vev=\f$v_1^2+v_2^2\f$
  void setDiagYukawas(const QedQcd &, double vev);
  /// Defines mixed Yukawa matrices from data input in form of CKM matrix and
  /// r, vev. If mix=2, all mixing is in down sector
  /// mix=1, all mixing is in up sector
  void getQuarkMixedYukawas(const QedQcd & r, const DoubleMatrix &
			    CKM, int mix, double vev);
  /// Calculate beta functions of SUSY preserving parameters of RPC MSSM
  DoubleVector beta() const;
  /// Calculate beta functions of SUSY preserving parameters of RPC MSSM
  EssmSusy beta(sBrevity &) const;
  /// Outputs one-loop anomlous dimensions gii given matrix inputs.
  /// for RH leptons, LH leptons, LH quarks, RH downs, RH ups, H1 and H2
  /// respectively. Note that we use the convention (for matrices in terms of
  /// gamma's): gamma^Li_Lj = M_ij for LH fields and
  /// gamma^Rj_Ri = M_ij for RH fields (since they are really the complex
  /// conjugates of the RH fields). a should already be defined.
  void getOneLpAnom(DoubleMatrix & gEE, DoubleMatrix & gLL,
		    DoubleMatrix & gQQ, DoubleMatrix & gDD,
		    DoubleMatrix & gUU, double & gH1H1, double &
		    gH2H2, sBrevity & a) const; 
  /// Outputs two-loop anomlous dimensions gii given matrix inputs.
  /// for RH leptons, LH leptons, LH quarks, RH downs, RH ups, H1 and H2
  /// respectively. Note that we use the convention (for matrices in terms of
  /// gamma's): gamma^Li_Lj = M_ij for LH fields and
  /// gamma^Rj_Ri = M_ij for RH fields (since they are really the complex
  /// conjugates of the RH fields). a should already be defined.
  void getTwoLpAnom(DoubleMatrix & gEE, DoubleMatrix & gLL,
		    DoubleMatrix & gQQ, DoubleMatrix & gDD,
		    DoubleMatrix & gUU, double & gH1H1, double &
		    gH2H2, sBrevity & a) const; 
  /// Outputs wave function renormalisation for SUSY parameters and gauge beta
  /// functions up to 2 loops. Also calculates and outputs a.
  /// IO parameters: RH leptons, LH leptons, LH quarks, RH downs, RH ups, H1
  /// and H2 respectively. 
  /// g^Li_Lj = m_{ij} for LH fields
  /// g^Ei_Ej = m_{ji} for RH fields
  //Peter:: edited to include new ESSM gage coupling
  void anomalousDimension(DoubleMatrix & gEE, DoubleMatrix & gLL,
			  DoubleMatrix & gQQ, DoubleMatrix & gUU,
			  DoubleMatrix & gDD, DoubleVector & dg, 
			  double & dgdash_1, double & dg_11, double & gH1H1, 
			  double & gH2H2, sBrevity & a) const; 
  /// Rotates to quark mass basis, returning the mixing matrices defined as 
  /// \f$(Y_U)_{diag}\f$ = vul \f$Y_U\f$ vur^+,  
  /// \f$(Y_D)_{diag}\f$ = vdl \f$Y_D\f$ vdr^+  
  /// All matrices should be 3 by 3
  void diagQuarkBasis(DoubleMatrix & vdl, DoubleMatrix & vdr, 
		      DoubleMatrix & vul, DoubleMatrix & vur) const;
};
/// Formatted output
ostream & operator <<(ostream &, const EssmSusy &);
/// Formatted input
istream & operator >>(istream &left, EssmSusy &s);
/// Outputs beta function coefficients for MSSM gauge coupling evolution in
/// arguments. 
void setESSMBetas(DoubleMatrix &, DoubleVector  &, DoubleVector  &, DoubleVector
	       &, DoubleVector  &);

inline EssmSusy EssmSusy::displaySusy() const { return *this; }

inline void EssmSusy::setGaugeCoupling(int i, double f) { g(i) = f; }

inline void EssmSusy::setAllGauge(const DoubleVector & v) { 
  if (v.displayStart() != 1 || v.displayEnd() !=3) {
    ostringstream ii;
    ii << 
      "Initialising SUSY params gauge function with vector NOT 1..3\n" <<
      v;
    throw ii.str();
  }
  g = v; 
}

inline double EssmSusy::displayHvev() const { return hVev; } 

inline void EssmSusy::setHvev(double h) { hVev = h; }
inline void EssmSusy::setSusyMu(double f) { smu = f; }
inline void EssmSusy::setTanb(double f) { tanb = f; }
inline void EssmSusy::seth_N(double f ){h_N = f; }
 inline  void  EssmSusy::setgash_1(double f){gdash_1 = f; }
  inline  void  EssmSusy::setg_11(double f){g_11 = f; }
inline void  EssmSusy::setmu_0(double f){mu_0 = f;} 
 inline void  EssmSusy::setkappa(DoubleVector v){kappa = v;}
 inline void  EssmSusy::setlambda(DoubleVector v){lambda = v;}
inline void  EssmSusy::setMN_SUSY(DoubleVector v){MN_SUSY = v;}

inline DoubleVector EssmSusy::displayGauge() const { return g; }
inline double EssmSusy::displayGaugeCoupling(int i) const { 
  return g.display(i); 
}
inline double EssmSusy::displaySusyMu() const { return smu; }

 

#endif



