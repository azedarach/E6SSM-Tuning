/*
  E6SSM version of softpars.cpp
  Neglects U(1) mixing for now
 */

#ifndef ESSMSOFTPARS_H
#define ESSMSOFTPARS_H

#include <cmath>
#include <susy.h>
#include <def.h>
#include <linalg.h>
#include <utils.h>
#include <numerics.h>
#include <softpars.h> ///< DH:: For soft masses typedefs etc.
#include "essmsusy.h"
using namespace softsusy;


/// Number of parameters contained in RGEs
const static int numSoftParsEssm = 99 + numEssmSusyPars;

/// Soft SUSY breaking parameters and beta functions. 
class SoftParsEssm: public EssmSusy {
private:
  DoubleVector mGaugino; ///< Gaugino masses, see ::beta for definitions
  DoubleMatrix ua, da, ea; ///< MSSM trilinear soft terms..
  /// soft mass squared matrices of \f$ m_Q^2, m_U^2, m_D^2, m_L^2, m_E^2 \f$
  /// respectively.
  DoubleMatrix mQLsq, mURsq, mDRsq, mLLsq, mSEsq; 
  /// Bilinear Higgs soft parameters: \f$ m_3^2, m_{H_1}^2, m_{H_2}^2 \f$
  /// respectively. 
  /// DH:: m_3^2 doesn't exist in the E6SSM, but m_3'^2 (= B\mu') does.
  double m3Prsq;
  double m32;         ///< Gravitino mass

  /// DH:: E6SSM soft parameters
  /// DH:: Soft squared masses for Higgs, singlets and exotics.
  /// For now we are assuming that they are diagonal (hence double vectors)
  DoubleVector mH1sq, mH2sq, mSsq, mDsq, mDbarsq;

  /// DH:: Diagonal E6SSM soft trilinears (la(i) = lambda(i)*A_lambda(i) etc.)
  DoubleVector la, ka;

  /// DH:: Soft squared masses for survival Higgs fields
  double mHPrsq, mHbarPrsq;

public:
  /// Default constructor fills object with zeroes
  SoftParsEssm();
  /// Constructor fills SUSY conserving parts with another object, all
  /// SUSY breaking parameters set to zero
  SoftParsEssm(const EssmSusy &);
  /// Constructor sets all parameters equal to those in another object
  SoftParsEssm(const SoftParsEssm &);
  /// Sets all parameters equal to those in another object
  const SoftParsEssm & operator=(const SoftParsEssm & s);
  /// Constructor sets RPC SUSY parameters to s, gaugino masses to mG,
  /// trilinears to aU, aD, aE for au, ad, ae
  /// trilnears respectively,  \f$m_Q^2\f$=mQl, \f$m_U^2\f$=mUr,
  /// \f$m_D^2\f$=mDr, \f$m_L^2\f$=mLl, \f$m_E^2\f$=mEr, \f$ m_3^2\f$=m3sq,
  /// \f$m_{H_1}^2\f$=mH1sq, \f$m_{H_2}^2\f$=mH2sq, mu parameter, number of
  /// loops=l, and threshold parameter=t
  SoftParsEssm(const EssmSusy & s, const DoubleVector & mG, const
	       DoubleMatrix & aU, const DoubleMatrix & aD, const DoubleMatrix
	       & aE, const DoubleMatrix & mQl, const DoubleMatrix & mUr, const
	       DoubleMatrix & mDr, const DoubleMatrix & mLl, const
	       DoubleMatrix & mEr, double m3sq, double mH1sq, double mH2sq, 
	       double mGravitino, double mu, int l, int t);

  /// Returns whole object as a const
  inline const SoftParsEssm & displaySoftPars() const;

  /// Return a trilinear coupling matrix
  const DoubleMatrix & displayTrilinear(trilinears) const;
  /// Return a trilinear element
  double displayTrilinear(trilinears, int i, int j) const;
  /// Return a trilinear element in "SUGRA style"
  double displaySoftA(trilinears, int, int) const;
  /// Return a soft mass squared matrix
  const DoubleMatrix & displaySoftMassSquared(softMasses) const;
  /// Return a soft mass squared element
  double displaySoftMassSquared(softMasses, int i, int j) const;

  double displayGravitino() const; ///< Returns the gravitino mass
  inline double displayM3Squared() const;     ///< Return \f$ m_3^2\f$
  inline double displayMh1Squared() const;    ///< Return \f$m_{H_1}^2\f$
  inline double displayMh2Squared() const;    ///< Return \f$m_{H_2}^2\f$=mH2sq
  inline DoubleVector displayGaugino() const; ///< Return \f$M_{G_i}\f$
  inline double displayGaugino(int i) const;  ///< Return \f$M_{G_i}\f$
  /// Return contents of object in a vector: for RG evolution
  virtual const DoubleVector display() const; 
  
  /// Sets gravitino mass
  void setM32(double);
  /// Sets whole thing equal to another object
  void setSoftPars(SoftParsEssm const &);
  /// Set one element of a soft mass squared matrix 
  void setSoftMassElement(softMasses, int, int, double);
  /// Set whole of a soft mass squared matrix
  void setSoftMassMatrix(softMasses, const DoubleMatrix &);
  /// Set whole of a trilinear SUSY breaking parameter matrix
  void setTrilinearMatrix(trilinears, const DoubleMatrix &);
  /// Set one element of a trilinear SUSY breaking parameter matrix
  void setTrilinearElement(trilinears k, int i, int j, double a);
  /// Set one gaugino mass
  void setGauginoMass(int, double);
  /// Set all gaugino masses
  void setAllGauginos(const DoubleVector &);
  void setM3Squared(double);  ///< Sets \f$ m_3^2\f$
  void setMh1Squared(double); ///< Sets \f$ m_{H_1}^2\f$
  void setMh2Squared(double); ///< Sets \f$ m_{H_2}^2\f$
  /// Sets total set of RGE parameters equal to elements of a vector
  void set(const DoubleVector &);
  //  void setSusy(const EssmSusy &);

  /// Returns double vector containing numerical beta functions of parameters
  DoubleVector beta() const; 
  /// Returns numerical beta functions of parameters  
  SoftParsEssm beta2() const; 
  /// Returns derivatives of anomalous dimensions of fields with respect to
  /// renormalisation scale in MSSM for: RH leptons, LH leptons, LH quarks, RH
  /// up quarks, RH down quarks, H1 and H2 respectively
  void anomalousDeriv(DoubleMatrix & gEE, DoubleMatrix & gLL,
		      DoubleMatrix & gQQ, DoubleMatrix & gUU,
		      DoubleMatrix & gDD, 
		      double & gH1H1, double & gH2H2) const;
  /// Ytilde quantities are for calculational brevity in beta functions.
  void yTildes(DoubleMatrix & yu, DoubleMatrix & yd, DoubleMatrix &ye) const;
  
  /// Reads in universal cE6SSM boundary conditions at the current scale:
  /// m0, M1/2, A0, B-parameter and mu
  void universal(double m0,  double m12,  double a0,  double mu,
		 double m3sq);
  /// Reads in universal boundary conditions at the current scale: m0, M1/2, A0
  void standardSugra(double m0,  double m12, double a0);
  /// Sets all flavour-diagonal SUSY breaking scalar masses to m0
  void universalScalars(double m0);
  /// Sets all flavour-diagonal SUSY breaking gaugino masses to m12
  void universalGauginos(double m12);
  /// Sets all SUSY breaking trilinear couplings to a0
  void universalTrilinears(double a0);

  /// Reads in soft SUSY breaking parameters from a file
  void inputSoftParsOnly();
};

/// Formatted ouput of whole object
ostream & operator <<(ostream &left, const SoftParsEssm &s);
/// Formatted input of whole object
istream & operator >>(istream &left, SoftParsEssm &s);

inline SoftParsEssm::SoftParsEssm()
  : EssmSusy(), mGaugino(3), ua(3, 3), da(3, 3), ea(3, 3), 
  mQLsq(3, 3), mURsq(3, 3), mDRsq(3, 3), mLLsq(3, 3), mSEsq(3, 3), m3sq(0.0),
  mH1sq(0.0), mH2sq(0.0), m32(1.e19) {      
  
  setPars(numSoftParsEssm);
  setMu(0.0);
  setLoops(0);
  setThresholds(0);
}

inline SoftParsEssm::SoftParsEssm(const SoftParsEssm & s)
  : EssmSusy(s.displaySusy()), 
  mGaugino(s.displayGaugino()), ua(s.displayTrilinear(UA)),
  da(s.displayTrilinear(DA)), ea(s.displayTrilinear(EA)),
  mQLsq(s.displaySoftMassSquared(mQl)), 
  mURsq(s.displaySoftMassSquared(mUr)), 
  mDRsq(s.displaySoftMassSquared(mDr)),
  mLLsq(s.displaySoftMassSquared(mLl)),
  mSEsq(s.displaySoftMassSquared(mEr)), 
  m3sq(s.displayM3Squared()), mH1sq(s.displayMh1Squared()),
  mH2sq(s.displayMh2Squared()), m32(s.displayGravitino()) {
  
  setPars(numSoftParsEssm);
  setMu(s.displayMu()); 
  setLoops(s.displayLoops());
  setThresholds(s.displayThresholds());
}

inline SoftParsEssm::SoftParsEssm(const EssmSusy &s)
  : EssmSusy(s), mGaugino(3), ua(3, 3), da(3, 3), ea(3, 3), 
    mQLsq(3, 3), mURsq(3, 3), mDRsq(3, 3), mLLsq(3, 3), mSEsq(3, 3), m3sq(0.0),
    mH1sq(0.0),  mH2sq(0.0), m32(1.e19) { 
      setPars(numSoftParsEssm);
      setMu(s.displayMu()); 
      setLoops(s.displayLoops());
      setThresholds(s.displayThresholds());
}

inline SoftParsEssm::SoftParsEssm
(const EssmSusy & s, const DoubleVector & mG, const
 DoubleMatrix & aU, const DoubleMatrix & aD, const DoubleMatrix & aE, const
 DoubleMatrix & mQl, const DoubleMatrix & mUr, const DoubleMatrix & mDr, const
 DoubleMatrix & mLl, const DoubleMatrix & mEr, double m3sqn, double mH1sq,
 double mH2sq, double mg, double mu, int l, int t)
  : EssmSusy(s), mGaugino(mG), ua(aU), da(aD), ea(aE),
    mQLsq(mQl), mURsq(mUr), mDRsq(mDr), mLLsq(mLl), mSEsq(mEr), m3sq(m3sqn),
    mH1sq(mH1sq), mH2sq(mH2sq), m32(mg) {
      setPars(numSoftParsEssm);
      setMu(mu);
      setLoops(l);
      setThresholds(t);
}

inline const SoftParsEssm & SoftParsEssm::displaySoftPars() const { return *this; }

inline double SoftParsEssm::displayM3Squared() const { return m3sq; }

inline double SoftParsEssm::displayMh1Squared() const { return mH1sq; }

inline double SoftParsEssm::displayMh2Squared() const { return mH2sq; }

inline DoubleVector SoftParsEssm::displayGaugino() const { return mGaugino; }

inline double SoftParsEssm::displayGaugino(int i) const { 
  return mGaugino.display(i); 
}

inline double SoftParsEssm::displayGravitino() const { return m32; }

inline void SoftParsEssm::setGauginoMass(int i, double f) { 
  mGaugino(i) = f; 
}

inline void SoftParsEssm::setM3Squared(double f) { m3sq = f; }
inline void SoftParsEssm::setMh1Squared(double f) { mH1sq = f; }
inline void SoftParsEssm::setMh2Squared(double f) { mH2sq = f; }
inline void SoftParsEssm::setSoftPars(SoftParsEssm const & s) { *this = s; }
inline void SoftParsEssm::setM32(double a) { m32 = a; }
#endif

