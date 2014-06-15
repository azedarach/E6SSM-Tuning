/*
  E6SSM version of SOFTSUSY's susy.h
  Currently neglects U(1) mixing in RGEs
  Search for DH::TODO to find incomplete sections...
 */

#ifndef ESSMSUSY_H
#define ESSMSUSY_H

#include <iostream>
#include <cmath>
#include <fstream>
#include <lowe.h>
#include <utils.h>
#include <susy.h> ///< DH:: for Yukawa type-defs etc.

using namespace softsusy; ///< DH:: while using SOFTSUSY-3.3.10, might consider 
                          ///  moving into softsusy namespace if using with 3.4.x

const static int numEssmSusyPars = 44;
const static int numEssmCharges = 13;

/// DH:: we want this class to be able to represent a 
/// more generic E6 model, with the surviving U(1)'
/// not necessarily U(1)_N. Define the U(1)_\chi and
/// U(1)_\psi charges as constants, in the usual
/// GUT normalisation.
const static double QQPsi = 1.0/(2.0*sqrt(6.0));
const static double QuPsi = 1.0/(2.0*sqrt(6.0));
const static double QdPsi = 1.0/(2.0*sqrt(6.0));
const static double QLPsi = 1.0/(2.0*sqrt(6.0));
const static double QePsi = 1.0/(2.0*sqrt(6.0));
const static double QNPsi = 1.0/(2.0*sqrt(6.0));
const static double QSPsi = 4.0/(2.0*sqrt(6.0));
const static double QH1Psi = -1.0/sqrt(6.0);
const static double QH2Psi = -1.0/sqrt(6.0);
const static double QXPsi = -1.0/sqrt(6.0);
const static double QXbarPsi = -1.0/sqrt(6.0);
const static double QHPrPsi = 1.0/(2.0*sqrt(6.0));
const static double QHbarPrPsi = -1.0/(2.0*sqrt(6.0));

const static double QQChi = -1.0/(2.0*sqrt(10.0));
const static double QuChi = -1.0/(2.0*sqrt(10.0));
const static double QdChi = 3.0/(2.0*sqrt(10.0));
const static double QLChi = 3.0/(2.0*sqrt(10.0));
const static double QeChi = -1.0/(2.0*sqrt(10.0));
const static double QNChi = -5.0/(2.0*sqrt(10.0));
const static double QSChi = 0.0;
const static double QH1Chi = -1.0/sqrt(10.0);
const static double QH2Chi = 1.0/sqrt(10.0);
const static double QXChi = 1.0/sqrt(10.0);
const static double QXbarChi = -1.0/sqrt(10.0);
const static double QHPrChi = 3.0/(2.0*sqrt(10.0));
const static double QHbarPrChi = -3.0/(2.0*sqrt(10.0));

/// DH:: For convenience we will also define the GUT normalised
/// U(1)_Y charges as constants here as well.
const static double QQY = sqrt(0.6)*(1.0/6.0);
const static double QuY = sqrt(0.6)*(-2.0/3.0);
const static double QdY = sqrt(0.6)*(1.0/3.0);
const static double QLY = sqrt(0.6)*(-1.0/2.0);
const static double QeY = sqrt(0.6);
const static double QNY = 0.0;
const static double QSY = 0.0;
const static double QH1Y = sqrt(0.6)*(-1.0/2.0);
const static double QH2Y = sqrt(0.6)*(1.0/2.0);
const static double QXY = sqrt(0.6)*(-1.0/3.0);
const static double QXbarY = sqrt(0.6)*(1.0/3.0);
const static double QHPrY = sqrt(0.6)*(-1.0/2.0);
const static double QHbarPrY = sqrt(0.6)*(1.0/2.0);

/// DH:: For accessing the U(1)' charges
typedef enum {QQ=1, Qu, Qd, QL, Qe, QN, QS, QH1, QH2, QX, QXbar, QHPr, QHbarPr} U1Charge;

/// DH::TODO In SOFTSUSY-3.4.1, the corresponding construction for the NMSSM
/// inherits from the version for the MSSM. For now we will just create
/// a separate structure, consider merging with everything else later.
struct essBrevity {
  /// dt=\f$Y_D^T\f$, ut=\f$Y_U^T\f$, et=\f$Y_E^T\f$, u2=\f$Y_U Y_U^T\f$,
  /// d2=\f$Y_D Y_D^T\f$, e2=\f$Y_E Y_E^T\f$, 
  /// u2t=\f$Y_U^T Y_U\f$, d2t=\f$Y_D^T Y_D\f$, e2t=\f$Y_E^T Y_E\f$
  DoubleMatrix dt, ut, et, u2, d2, e2, u2t, e2t, d2t;
  /// DH:: Note that this only contains the non-Abelian gauge couplings
  /// gsq=\f$g_i^2\f$, g3=\f$g_i^3\f$, g4=\f$g_i^4\f$
  DoubleVector gsq, g3, g4;

  DoubleMatrix Gasq, Ga3, Ga4;


  /// uuT=\f$Tr(Y_U^T Y_U)\f$, ddt=\f$Tr(Y_D^T Y_D)\f$, eet=\f$Tr(Y_E^T Y_E)\f$
  double uuT, ddT, eeT;
  /// d1=\f$Y_D\f$, u1=\f$Y_U\f$, e1=\f$Y_E\f$
  DoubleMatrix u1, d1, e1;

  /// DH:: Written in matrix notation for later, even though currently flavour diagonal.
  /// DH:: l2=\f$\lambda \lambda^T\f$, k2=\f$\kappa \kappa^T\f$
  DoubleMatrix lt, kt, l2, k2, l2t, k2t;
  /// DH:: l1 = \f$\lambda\f$, k1=\f$\kappa\f$
  DoubleMatrix l1, k1;
  /// DH:: lllT=\f$Tr(lambda^T lambda)\f$, kkT=\f$Tr(kappa^T kappa)\f$
  double llT, kkT;

  /// Constructor fills sruct with zeroes by default
  essBrevity();
  /// Constructor sets struct to be equal to another
  essBrevity(const essBrevity &);
  /// Sets struct to be equal to another
  const essBrevity & operator=(const essBrevity &);
  
  /// Calculates all the data given Yukawa matrices yu, yd and ye and gauge
  /// couplings g
  void calculate(const DoubleMatrix & yu, const DoubleMatrix & yd, const
		 DoubleMatrix & ye, const DoubleVector & g, 
		 const DoubleMatrix & Ga, const DoubleVector & lam, 
		 const DoubleVector & kap);
};

inline essBrevity::essBrevity()
		  : dt(3, 3), ut(3, 3), et(3, 3), u2(3, 3), 
		  d2(3, 3), e2(3, 3), u2t(3, 3), 
		  e2t(3, 3), d2t(3, 3), gsq(1, 2), g3(1, 2), g4(1, 2), uuT(0.0),
		  ddT(0.0), eeT(0.0), u1(3, 3), d1(3, 3), e1(3, 3),
		  Gasq(2,2), Ga3(2,2), Ga4(2,2), lt(3,3), kt(3,3), l2(3,3),
		  k2(3,3), l2t(3,3), k2t(3,3), l1(3,3), k1(3,3), llT(0.0), kkT(0.0)
{}

inline essBrevity::essBrevity(const essBrevity &s)
		  : dt(s.dt), ut(s.ut), et(s.ut), u2(s.u2), d2(s.d2),
		  e2(s.e2), u2t(s.u2t),
		  e2t(s.e2t), d2t(s.d2t), gsq(s.gsq), g3(s.g3), g4(s.g4),
		  uuT(s.uuT), ddT(s.ddT), eeT(s.eeT), u1(s.u1), d1(s.d1), e1(s.e1),
		  Gasq(s.Gasq), Ga3(s.Ga3), Ga4(s.Ga4), lt(s.lt), kt(s.kt), l2(s.l2),
		  k2(s.k2), l2t(s.l2t), k2t(s.k2t), l1(s.l1), k1(s.k1), llT(s.llT), kkT(s.kkT) 
{}


/// DH:: Contains all supersymmetric ESSM parameters and RGEs 
/// Unlike in the NMSSM, the ESSM has a different gauge group and
/// so parameters that make sense in the MSSM (e.g. \mu) don't
/// in the ESSM. Obviously we therefore don't want to inherit 
/// any methods that allow the user to access or set any of those
/// parameters. For that reason I going to change from Peter's code
/// and inherit directly from RGE instead (at least that is the
/// simplest way to do it).
class EssmSusy: public RGE {
private:
  DoubleMatrix u, d, e; ///< Yukawa matrices for ups, downs and leptons
  DoubleVector g; ///< Non-Abelian gauge couplings in GUT normalisation, order is \{ g_2, g_3 \}.
  DoubleMatrix Gmat; ///< DH:: Abelian gauge couplings in GUT normalisation

  /// DH:: Note the MSSM \mu parameter is forbidden by U(1)' gauge invariance.
  /// DH:: Ratio of third generation Higgs VEVs, tan(beta) = v2/v1 and
  /// v = \sqrt{v1^2+v_2^2}.
  double tanb, hVev;

  /// DH:: Third generation singlet VEV
  double sVev; 

  /// DH:: ESSM superpotential Yukawa couplings. Assumes approximate
  /// Z_2^H symmetry and defined in the basis in which the singlet
  /// couplings are flavour diagonal.
  DoubleVector lambda, kappa;
  
  /// DH:: E6 mixing angle for U(1)' charges
  double thetaE6;

  /// DH:: SUSY \mu' bilinear mass parameter for survival Higgs fields
  double smuPr;

  /// DH:: We are currently neglecting a large number of exotic couplings,
  /// that have to be small anyway. I am also going to leave out the 
  /// heavy right-handed neutrinos - if U(1)' != U(1)_N then there will
  /// be problems with e.g. including mass terms for them. 

public:
  EssmSusy(); ///< Constructor fills object with zeroes by default
  /// Constructor sets object to be equal to another
  EssmSusy(const EssmSusy &); 

  /// DH:: Constructor given Yukawa matrices u, d, e, non-Abelian gauge 
  /// couplings v, Abelian gauge couplings G, singlet-Higgs couplings
  /// lam, singlet-exotic couplings kap, mu' parameter mp, 
  /// tan(beta) tb, Higgs VEV h, singlet VEV s, E6 mixing angle th, 
  /// renormalisation scale Q, number of loops in RG evolution l
  /// and thresholds parameter t
  EssmSusy(const DoubleMatrix & u, const DoubleMatrix & d, const
	   DoubleMatrix & e, const DoubleVector & v, 
	   const DoubleMatrix & G, const DoubleVector & lam, 
	   const DoubleVector & kap, double mp,
	   double tb, double h, double s, double th, 
	   double Q, int l, int t);
  virtual ~EssmSusy() {}; ///< Default destructor

  /// sets object to be equal to another
  const EssmSusy & operator=(const EssmSusy & s);
  /// sets object to be equal to another
  void setSusy(const EssmSusy &s); 
  /// Sets DRbar running Higgs vev
  void setHvev(double h);
  /// DH:: Sets DRbar running singlet vev
  void setSvev(double s);
  /// DH:: Copies SM Yukawas, E6 Yukawas and all gauge couplings
  /// from s only (includes both non-Abelian and Abelian gauge 
  /// couplings, excludes VEVs, tan(beta) and \mu').
  void setSomePars(const EssmSusy & s);
  /// Sets one element of a Yukawa matrix
  void setYukawaElement(yukawa, int, int, double);
  /// Sets whole Yukawa matrix
  void setYukawaMatrix(yukawa, const DoubleMatrix &);
  /// DH:: Set a single non-Abelian gauge coupling
  void setNonAbelianGaugeCoupling(int, double);
  /// DH:: Set all non-Abelian gauge couplings
  void setAllNonAbelianGauge(const DoubleVector  &);
  /// DH:: Set a single Abelian gauge coupling
  void setAbelianGaugeCoupling(int, int, double);
  /// DH:: Set all non-Abelian gauge couplings
  void setAllAbelianGauge(const DoubleMatrix &);
  /// DH:: Sets superpotential mu' parameter
  void setSusyMuPrime(double);
  /// DH:: Sets a single E6SSM \lambda coupling
  void setLambdaElement(int, double);
  /// DH:: Sets all E6SSM \lambda couplings
  void setLambda(const DoubleVector &);
  /// DH:: Sets a single E6SSM \kappa coupling
  void setKappaElement(int, double);
  /// DH:: Sets all E6SSM \kappa couplings
  void setKappa(const DoubleVector &);
  /// Sets all RGE parameters to elements of vector
  void set(const DoubleVector &);
  /// Sets tan beta
  void setTanb(double);

  /// Returns DRbar running Higgs vev
  double displayHvev() const;
  /// DH:: Returns DRbar running singlet vev
  double displaySvev() const;
  /// Returns whole object as a const
  inline const EssmSusy & displaySusy() const;
  /// Returns a single Yukawa matrix element
  double displayYukawaElement(yukawa, int, int) const;
  /// Returns a whole Yukawa matrix
  const DoubleMatrix & displayYukawaMatrix(yukawa) const;
  /// DH:: Returns a single non-Abelian gauge coupling
  double displayNonAbelianGaugeCoupling(int) const;
  /// Returns all non-Abelian gauge couplings
  DoubleVector displayNonAbelianGauge() const;
  /// DH:: Returns a single U(1) coupling
  double displayAbelianGaugeCoupling(int, int) const;
  /// DH:: Returns matrix of U(1) couplings
  DoubleMatrix displayAbelianGauge() const;
  /// DH:: Display a single E6SSM \lambda coupling
  double displayLambdaElement(int) const;
  /// DH:: Display all \lambda couplings
  DoubleVector displayLambda() const;
  /// DH:: Display a single E6SSM \kappa coupling
  double displayKappaElement(int) const;
  /// DH:: Display all \kappa couplings
  DoubleVector displayKappa() const;
  /// DH:: Returns superpotential mu' parameter
  double displaySusyMuPrime() const;
  /// DH:: Returns E6 mixing angle
  double displayThetaE6() const;
  /// Returns all parameters as elements of a vector
  const DoubleVector display() const;
  /// Returns tan beta
  double displayTanb() const;
  /// DH:: Returns U(1)_Y charges for the superfields
  /// as a vector. Ordering is {Q, u, d, L, e, N, S, H1, H2, D, Dbar, H', H'bar}
  DoubleVector displayU1Y() const;
  /// DH:: Returns a single U(1)_Y charge
  double displayU1YCharge(U1Charge) const;
  /// DH:: Returns U(1)' charges for the superfields
  /// as a vector. Ordering is {Q, u, d, L, e, N, S, H1, H2, D, Dbar, H', H'bar}
  DoubleVector displayU1Prime() const;
  /// DH:: Returns a single U(1)' charge
  double displayU1PrimeCharge(U1Charge) const;
  /// DH:: Returns scale-dependent effective U(1)' charges for the superfields
  /// as a vector. Ordering is {Q, u, d, L, e, N, S, H1, H2, D, Dbar, H', H'bar}
  DoubleVector displayU1PrimeEff() const;
  /// DH:: Returns a single scale-dependent effective U(1)' charge
  double displayU1PrimeEffCharge(U1Charge) const;
  
  /// outputs object QedQcd & r valid at 1 GeV from SUSY data at mt, from
  /// diagonal elements of Yukawa couplings and Higgs VEV vev.
  /// DH:: NB vev = \frac{v}{\sqrt{2}}
  void getMasses(QedQcd & r, double vev) const;

  /// This turns diagonal Yukawa couplings at MZ into CKM mixed ones
  /// Takes diagonal quark Yukawa matrices and mixes them up
  /// according to the CKM matrix assuming:
  /// mix=2, all mixing is in down sector
  /// mix=1, all mixing is in up sector
  /// DH:: I will neglect any possible changes that need to be done here
  void quarkMixing(const DoubleMatrix & CKM, int mix);

  /// Sets diagonal Yukawa couplings according to data in QedQcd input and
  /// Higgs VEV parameter vev=\f$v_1^2+v_2^2\f$
  /// DH:: NB vev = \frac{v}{\sqrt{2}}
  void setDiagYukawas(const QedQcd &, double vev);

  /// Defines mixed Yukawa matrices from data input in form of CKM matrix and
  /// r, vev. If mix=2, all mixing is in down sector
  /// mix=1, all mixing is in up sector
  /// DH:: NB vev = \frac{v}{\sqrt{2}}
  void getQuarkMixedYukawas(const QedQcd & r, const DoubleMatrix &
				    CKM, int mix, double vev);

  /// DH:: Temporary dominant third-family versions DH::TODO modify to 
  /// include full 3 family effects
  /// DH:: Calculate beta functions of SUSY preserving parameters of ESSM DH::TODO
  DoubleVector beta() const;
  /// DH:: Calculate beta functions of SUSY preserving parameters of ESSM DH::TODO
  EssmSusy beta(essBrevity &) const;

  

  
  /// Outputs one-loop anomalous dimensions gii given matrix inputs.
  /// for RH leptons, LH leptons, LH quarks, RH downs, RH ups, H1 and H2
  /// respectively. Note that we use the convention (for matrices in terms of
  /// gamma's): gamma^Li_Lj = M_ij for LH fields and
  /// gamma^Rj_Ri = M_ij for RH fields (since they are really the complex
  /// conjugates of the RH fields). a should already be defined. DH::TODO
  void getOneLpAnom(DoubleMatrix & gEE, DoubleMatrix & gLL,
		    DoubleMatrix & gQQ, DoubleMatrix & gDD,
		    DoubleMatrix & gUU, DoubleMatrix & gXX,
		    DoubleMatrix & gXbarXbar, DoubleVector & gH1H1, 
		    DoubleVector & gH2H2, DoubleVector & gSS,
		    double & gHPrHPr, double & gHbarPrHbarPr, 
		    essBrevity & a) const;
		    /*
  /// Outputs two-loop anomalous dimensions gii given matrix inputs.
  /// for RH leptons, LH leptons, LH quarks, RH downs, RH ups, H1 and H2
  /// respectively. Note that we use the convention (for matrices in terms of
  /// gamma's): gamma^Li_Lj = M_ij for LH fields and
  /// gamma^Rj_Ri = M_ij for RH fields (since they are really the complex
  /// conjugates of the RH fields). a should already be defined. DH::TODO
  void getTwoLpAnom(DoubleMatrix & gEE, DoubleMatrix & gLL,
		    DoubleMatrix & gQQ, DoubleMatrix & gDD,
		    DoubleMatrix & gUU, double & gH1H1, double &
		    gH2H2, essBrevity & a) const; 
  */
  /// DH:: Three-loop anomalous dimensions for the E6SSM are not available
  /// presently.

  /// Outputs wave function renormalisation for SUSY parameters and gauge beta
  /// functions up to 2 loops. Also calculates and outputs a.
  /// IO parameters: RH leptons, LH leptons, LH quarks, RH downs, RH ups, H1
  /// and H2 respectively. 
  /// g^Li_Lj = m_{ij} for LH fields
  /// g^Ei_Ej = m_{ji} for RH fields DH::TODO
  void anomalousDimension(DoubleMatrix & gEE, DoubleMatrix & gLL,
			  DoubleMatrix & gQQ, DoubleMatrix & gUU,
			  DoubleMatrix & gDD, DoubleMatrix & gXX, 
			  DoubleMatrix & gXbarXbar, DoubleVector & gH1H1, 
			  DoubleVector & gH2H2, DoubleVector & gSS, 
			  double & gHPrHPr, double & gHbarPrHbarPr,
			  DoubleVector & dg, DoubleMatrix & dGa,
			  essBrevity & a)  const;

  

  /// Rotates to quark mass basis, returning the mixing matrices defined as 
  /// \f$(Y_U)_{diag}\f$ = vul \f$Y_U\f$ vur^+,  
  /// \f$(Y_D)_{diag}\f$ = vdl \f$Y_D\f$ vdr^+  
  /// All matrices should be 3 by 3
  void diagQuarkBasis(DoubleMatrix & vdl, DoubleMatrix & vdr, 
			  DoubleMatrix & vul, DoubleMatrix & vur) const;

  /// Outputs beta function coefficients for ESSM gauge coupling evolution in
  /// arguments. These depend on the choice of theta_E6, which is why it
  /// has been made a member function instead. 
  void setEssmBetas(DoubleMatrix &, DoubleVector  &, DoubleVector  &, DoubleVector &,
		    DoubleVector &, DoubleVector &, DoubleVector &, double) const;

};
/// Formatted output
ostream & operator <<(ostream &, const EssmSusy &);
/// Formatted input
istream & operator >>(istream &left, EssmSusy &s);

inline const EssmSusy & EssmSusy::displaySusy() const { return *this; }

inline void EssmSusy::setNonAbelianGaugeCoupling(int i, double f) { g(i) = f; }

inline void EssmSusy::setAbelianGaugeCoupling(int i, int j, double f) { Gmat(i, j) = f; }

// DH:: Note there are only two non-Abelian gauge couplings, g_2 = g(1) and g_3 = g(2)
inline void EssmSusy::setAllNonAbelianGauge(const DoubleVector & v) { 
  if (v.displayStart() != 1 || v.displayEnd() != 2) {
    ostringstream ii;
    ii << 
      "Initialising ESSM SUSY params non-Abelian gauge function with vector NOT 1..2\n" <<
      v;
    throw ii.str();
  }
  g = v; 

}

inline void EssmSusy::setAllAbelianGauge(const DoubleMatrix & G) {
  if (G.displayRows() != 2 || G.displayCols() != 2) {
    ostringstream ii;
    ii << "Initialising ESSM SUSY params Abelian gauge function with matrix NOT 2x2\n"
       << G;
    throw ii.str();
  }
  Gmat = G;
}

inline double EssmSusy::displayHvev() const { return hVev; } 
inline double EssmSusy::displaySvev() const { return sVev; }

inline void EssmSusy::setHvev(double h) { hVev = h; }
inline void EssmSusy::setSvev(double s) { sVev = s; }

// DH:: Note no \mu parameter
inline void EssmSusy::setSusyMuPrime(double f) { smuPr = f; }

inline void EssmSusy::setTanb(double f) { tanb = f; }

inline void EssmSusy::setLambdaElement(int i, double f) { lambda(i) = f; }
inline void EssmSusy::setLambda(const DoubleVector & lam) { 
  if (lam.displayStart() != 1 || lam.displayEnd() != 3) {
    ostringstream ii;
    ii << 
      "Initialising ESSM SUSY params lambda function with vector NOT 1..3\n" <<
      lam;
    throw ii.str();
  }
  lambda = lam; 

}
inline void EssmSusy::setKappaElement(int i, double f) { kappa(i) = f; }
inline void EssmSusy::setKappa(const DoubleVector & kap) { 
  if (kap.displayStart() != 1 || kap.displayEnd() != 3) {
    ostringstream ii;
    ii << 
      "Initialising ESSM SUSY params kappa function with vector NOT 1..3\n" <<
      kap;
    throw ii.str();
  }
  kappa = kap; 

}

// DH:: Note only returns non-Abelian gauge couplings
inline DoubleVector EssmSusy::displayNonAbelianGauge() const { return g; }

// DH:: Returns matrix of Abelian gauge couplings
inline DoubleMatrix EssmSusy::displayAbelianGauge() const { return Gmat; }

/// DH:: Returns a single non-Abelian gauge coupling
inline double EssmSusy::displayNonAbelianGaugeCoupling(int i) const { 
  return g.display(i); 
}

/// DH:: Returns a single Abelian gauge coupling
inline double EssmSusy::displayAbelianGaugeCoupling(int i, int j) const {
  return Gmat(i, j);
}

/// DH:: Returns a single E6SSM \lambda coupling
inline double EssmSusy::displayLambdaElement(int i) const {
  return lambda.display(i);
}

/// DH:: Returns vector of \lambda couplings
inline DoubleVector EssmSusy::displayLambda() const { return lambda; }

/// DH:: Returns a single E6SSM \kappa coupling
inline double EssmSusy::displayKappaElement(int i) const {
  return kappa.display(i);
}

/// DH:: Returns vector of \kappa couplings
inline DoubleVector EssmSusy::displayKappa() const { return kappa; }

/// DH:: Returns superpotential \mu' parameter, no \mu parameter
inline double EssmSusy::displaySusyMuPrime() const { return smuPr; }

/// DH:: Returns E6 mixing angle
inline double EssmSusy::displayThetaE6() const { return thetaE6; }

#endif




