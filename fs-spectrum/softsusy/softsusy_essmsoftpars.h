
/** \file softpars.h
   - Project:     SOFTSUSY
   - Author:      Ben Allanach 
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   - Webpage:     http://allanach.home.cern.ch/allanach/softsusy.html
   - Description: Soft SUSY breaking parameters

   $Log: softpars.h,v $
   Revision 1.3  2005/11/09 14:12:24  allanach
   Updated for SOFTSUSY 2.0.1 - cleaned everything up etc

   Revision 1.1.1.1  2004/11/19 16:18:31  allanach


   Revision 1.17  2004/01/15 13:54:55  allanach
   New heaer style implemented

   Revision 1.16  2003/10/24 16:09:04  allanach
   Implemented running Higgs DRbar vev

   Revision 1.14  2003/07/21 14:00:18  allanach
   MZ fully implemented as an input now. Kept MZ as the central PDG 2002 value,
   for defaults etc

   Revision 1.12  2003/05/27 15:05:52  allanach
   Purely efficiency corrections: used variable rather than display() methods
   whenever possible

   Revision 1.11  2003/05/20 15:19:40  allanach
   doxygen comment style implemented

   Revision 1.8  2002/07/30 12:57:32  allanach
   SOFTSUSY1.5

   Revision 1.7  2002/06/14 16:26:30  allanach
   Switches included for 2-loop running of scalar masses, and calulating 
   mt at mt.

   Revision 1.6  2002/04/26 15:14:44  allanach
   Deleted all translation routines and defined boundary conditions within
   softsusy.h and softsusy.cpp

   Revision 1.5  2002/04/14 13:50:41  allanach
   Now use V=m3^2 H1 H2 instead of V=mu B H1 H2. It's more natural!

   Revision 1.4  2002/04/12 16:51:27  allanach
   Added display/set functions to work automatically

   Revision 1.3  2002/04/12 06:24:50  allanach
   Code maintenance - returning a subobject made simpler

   Revision 1.3  2001/09/28 14:49:02  allanach
   GMSB boundary conditions added

   Revision 1.2  2001/07/18 14:42:51  allanach
   Added proper header info
*/

#ifndef SOFTSUSY_ESSMSOFTPARS_H
#define SOFTSUSY_ESSMSOFTPARS_H

#include <cmath>
#include "softsusy_essmsusy.h"
#include "softsusy_def.h"
#include "softsusy_linalg.h"
#include "softsusy_utils.h"
#include "softsusy_numerics.h"
#include "softsusy_softpars.h"

/// SUSY breaking soft mass squared parameters
//typedef enum {mQl=1, mUr, mDr, mLl, mEr} softMasses;
/// SUSY breaking trilinear parameters
//typedef enum {UA=1, DA, EA} trilinears;

/// Number of parameters contained in RGEs
//const static int numSoftParsMssm = 78 + numSusyPars;
//Peter:: I think the 78 actually gives the number -1.

//Peter:: edited to include new softssusy parameters.
const static int numSoftParsEssm = 111 + numSusyESSMPars;
/// Soft SUSY breaking para(meters and beta functions. 
class SoftParsEssm: public EssmSusy
{
private:
  DoubleVector mGaugino; ///< Gaugino masses, see ::beta for definitions
  DoubleMatrix ua, da, ea; ///< Trilinear soft terms..
  /// soft mass squared matrices of \f$ m_Q^2, m_U^2, m_D^2, m_L^2, m_E^2 \f$
  /// respectively.
  DoubleMatrix mQLsq, mURsq, mDRsq, mLLsq, mSEsq; 
  /// Bilinear Higgs soft parameters: \f$ m_3^2, m_{H_1}^2, m_{H_2}^2 \f$
  /// respectively. 
  //Peter:: i have commented out the line below since i need to use mH1sq as a double vector
//double m3sq, mH1sq, mH2sq; 
  double m32;         ///< Gravitino mass

  //Peter:: So i should add the new soft terms here I think?
  DoubleVector  mHdashsq, mHbardashsq,A_lambda, A_kappa, mSsq, mDsq, mDbarsq, MNcsq, mH1sq, mH2sq  ;
  double A_top, A_b, A_tau, A_N, B_0, mGdash_1;


  //Peter:: mH1/2sq may conflict with above possibly will have to change names or delete old versions and assosiate them with third generATION.
  //DoubleVector mH1sq, mH2sq, mDbarsq, mDsq;
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
  // s = SUSY part of ESSM model.
  // mG = gaugino masses (ASSUME THIS IS SO FOR NOW) - mG(3)~ gluino mass, mG(2)~ wino mass, mG(1)~bino mass
  // aU = SUSY breaking up Yukawas
  // aD = SUSY breaking down Yukawas
  // aE = SUSY breaking lepton Yukawas
  // mQl = squared soft masses for left-handed squarks
  // mUr = squared soft masses for right-handed up-type squarks
  // mDr = squared soft masses for right-handed down-type squarks
  // mLl = squared soft masses for left-handed sleptons
  // mEr = squared soft masses for right-handed sleptons
  // mH1sq = ESSM squared Higgs soft masses for H_1
  // mH2sq = ESSM squared Higgs soft masses for H_2
  // mg = gravitino soft mass
  // mu = renormalisation scale
  // l = # loops
  // t = threshold accuracy
  // mDbarsq = ESSM squared soft masses for exotic D_bar quarks
  // mDsq = ESSM squared soft masses for exotic D quarks
  // mGd = additional U(1)' gaugino mass (ASSUME THIS FOR NOW)
  // mHd = ESSM squared soft masses for H'
  // mHdb = ESSM squared soft masses for H_bar'
  // A_lambda = ESSM SUSY breaking SHH trilinear coupling (interested in A_lambda(3))
  // A_kap = ESSM SUSY breaking SDD_bar trilinear coupling
  // B_0 = ESSM SUSY breaking mass term corresponding to mu'H'H_bar' at GUT scale? CHECK THIS
  // mS = ESSM squared soft masses for singlet fields (m_s^2 = mS(3))
  // A_N = ESSM SUSY breaking right-handed neutrino coupling?
  // A_top = ESSM SUSY breaking t coupling
  // A_tau = ESSM SUSY breaking tau coupling
  // A_b = ESSM SUSY breaking b coupling
  // MNc = ESSM SUSY breaking right-handed neutrino masses? CHECK THIS
  SoftParsEssm(const EssmSusy & s, const DoubleVector & mG, const
 DoubleMatrix & aU, const DoubleMatrix & aD, const DoubleMatrix & aE, const
 DoubleMatrix & mQl, const DoubleMatrix & mUr, const DoubleMatrix & mDr, const
 DoubleMatrix & mLl, const DoubleMatrix & mEr, DoubleVector mH1sq,
	       DoubleVector mH2sq, double mg, double mu, int l, int t, DoubleVector mDbarsq, DoubleVector mDsq, double mGd,DoubleVector mHd,DoubleVector mHdb,DoubleVector A_lambda , DoubleVector A_kap, double B_0, DoubleVector mS, double A_N,double A_top, double A_tau, double A_b, DoubleVector MNc );



    //(const MssSusy & s, const DoubleVector & mG, const
    //       DoubleMatrix & aU, const DoubleMatrix & aD, const DoubleMatrix
    //       & aE, const DoubleMatrix & mQl, const DoubleMatrix & mUr, const
    //       DoubleMatrix & mDr, const DoubleMatrix & mLl, const
    //       DoubleMatrix & mEr, double m3sq, 
    //       double mGravitino, double mu, int l, int t, DoubleVector mDbarsq, DoubleVector mDsq, DoubleVector mGd,DoubleVector mHd,DoubleVector mHdb,DoubleVector A_lambda , DoubleVector A_kap,DoubleVector mHd,DoubleVector mHdb, double B_0, DoubleVector mS, double A_N,double A_top, double A_tau, double A_b, DoubleVector MNc );


	       //DoubleVector mHdashsq, DoubleVector mHbardashsq, DoubleVector A_lambda, DoubleVector A_kappa, DoubleVector mSsq, DoubleVector mDsq, DoubleVector MNcsq, DoubleVector mH1sq, DoubleVector mH2sq, double A_top, double A_b, double A_tau, double A_N, double B_0, double mGdash_1);

  /// Returns whole object as a const
  inline SoftParsEssm displaySoftPars() const;

  /// Return a trilinear coupling matrix
  DoubleMatrix displayTrilinear(trilinears) const;
  /// Return a trilinear element
  double displayTrilinear(trilinears, int i, int j) const;
  /// Return a trilinear element in "SUGRA style"
  double displaySoftA(trilinears, int, int) const;
  /// Return a soft mass squared matrix
  DoubleMatrix displaySoftMassSquared(softMasses) const;
  /// Return a soft mass squared element
  double displaySoftMassSquared(softMasses, int i, int j) const;

  double displayGravitino() const; ///< Returns the gravitino mass
  //  inline double displayM3Squared() const;     ///< Return \f$ m_3^2\f$
  inline DoubleVector displayMh1Squared() const;    ///< Return \f$m_{H_1}^2\f$
  inline DoubleVector displayMh2Squared() const;    ///< Return \f$m_{H_2}^2\f$=mH2sq
  inline DoubleVector displayGaugino() const; ///< Return \f$M_{G_i}\f$
  inline double displayGaugino(int i) const;  ///< Return \f$M_{G_i}\f$
  /// Return contents of object in a vector: for RG evolution
  virtual const DoubleVector display() const; 

  inline DoubleVector displaymDbarsq() const; 
  inline DoubleVector displaymDSq() const;

  inline double displaymGdash_1() const;
  inline DoubleVector displayA_lambda() const;
  inline DoubleVector displayA_kappa() const;



  inline DoubleVector displaymHdashsq() const;
  inline DoubleVector displaymHbardashsq() const;

  inline double displayB_0() const;
  inline DoubleVector displaymSsq() const;
  inline DoubleVector displayMNcsq() const;
  inline double displayA_top() const;
  inline double displayA_b() const;
  inline double displayA_N() const;
  inline double displayA_tau() const;

    


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
  
  //Peter:: ESSM changes
//void setM3Squared(double);  ///< Sets \f$ m_3^2\f$
  void setMh1Squared(DoubleVector); ///< Sets \f$ m_{H_1}^2\f$
  void setMh2Squared(DoubleVector); ///< Sets \f$ m_{H_2}^2\f$
  //Peter:: New functions to set new ESSm parameters:
  //double mGdash_1,A_lambda, A_kappa, mHdashsq, mHbardashsq,  B_0,mSsq;

  //Peter:: mH1/2sq may conflict with above possibly will have to change names or delete old versions and assosiate them with third generATION.
  //DoubleVector mH1sq, mH2sq, mDbarsq, mDsq;

 void setMGdash_1(double);
 void setA_lambda(DoubleVector);
 void setA_kappa(DoubleVector);
 void setmHdashsq(DoubleVector);
void setmHbardashsq(DoubleVector);

void setB_0(double);
//void setB0_tilde(double);
void setmSsq(DoubleVector);
void setmDsq(DoubleVector);

void setmDbarsq(DoubleVector);
void setA_top(double);
void setA_tau(double);
void setA_b(double);
void setA_N(double);
void setMNcsq(DoubleVector);

  /// Sets total set of RGE parameters equal to elements of a vector
  void set(const DoubleVector &);
  //  void setSusy(const MssmSusy &);

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
  
  /// Reads in universal boundary conditions at the current scale:
  /// m0, M1/2, A0, B-parameter and mu
  void universal(double m0,  double m12,  double a0,  double mu);
  /// Give it a SUSY object and a value of M3/2, and it will return a soft
  /// object with AMSB soft breaking terms. Note that the sleptons will be
  /// tachyonic, ie nothing has been done to fix that problem.
  /// Note that in the following, we are neglecting all Yukawa couplings except
  /// that of the third family.
  void addAmsb(double m32);
  /// Reads in universal boundary conditions at the current scale: m0, M1/2, A0
  void standardSugra(double m0,  double m12, double a0);
  /// Sets all flavour-diagonal SUSY breaking scalar masses to m0
  void universalScalars(double m0);
  /// Sets all flavour-diagonal SUSY breaking gaugino masses to m12
  void universalGauginos(double m12);
  /// Sets all SUSY breaking trilinear couplings to a0
  void universalTrilinears(double a0);

  /// Sets all flavour-diagonal SUSY breaking scalar masses to m0
  void universalScalars_with_new_ESSM_pars(double m0);
  /// Sets all flavour-diagonal SUSY breaking gaugino masses to m12
  void universalGauginos_with_new_ESSM_pars(double m12);
  /// Sets all SUSY breaking trilinear couplings to a0
  void universalTrilinears_with_new_ESSM_pars(double a0);
  /// Boundary conditions to be applied at messenger scale for Gauge mediated
  /// SUSY breaking (see hep-ph/9703211 for example), n5 is the number of
  /// 5-plets, mMess is the messenger scale and lambda is the GMSB scale
  void minimalGmsb(int n5, double lambda, double mMess, double cgrav);  

  /// Reads in soft SUSY breaking parameters from a file
  
  //Peter:: if i want to make this have ESSM code to have the same user input options as the MSSM origional I should edit this to take new params but for now hopefully this won't casue any compliation problems as simplys sets each soft parameter individually as entered by user
  void inputSoftParsOnly();
};

//Peter:: i think this is used when an error comes. just calls individual diplay functions so hopefully won't casue problems.  Edit to get more info on erros at some point
/// Formatted ouput of whole object
ostream & operator <<(ostream &left, const SoftParsEssm &s);
//Peter: i don't undersand how this one works fully, but may be ok without editing it??
/// Formatted input of whole object
istream & operator >>(istream &left, SoftParsEssm &s);


inline SoftParsEssm::SoftParsEssm()
  : EssmSusy(), mGaugino(3), ua(3, 3), da(3, 3), ea(3, 3), 
     mQLsq(3, 3), mURsq(3, 3), mDRsq(3, 3), mLLsq(3, 3), mSEsq(3, 3), mDbarsq(3), mDsq(3),mH1sq(3), mH2sq(3), m32(0.0), mGdash_1(0.0),A_lambda(3), A_kappa(3), mHdashsq(3), mHbardashsq(3), B_0(0.0),mSsq(3), A_N(0.0), A_top(0.0), A_b(0.0), A_tau(0.0),MNcsq(3)  {      
  
  setPars(numSoftParsEssm);
  setMu(0.0);
  setLoops(0);
  setThresholds(0);
}


//Peter:: editing to also transfer values on new variables
//Peter:: particularly unsure about the vectors
inline SoftParsEssm::SoftParsEssm(const SoftParsEssm & s)
  : EssmSusy(s.displaySusy()), 
  mGaugino(s.displayGaugino()), ua(s.displayTrilinear(UA)),
  da(s.displayTrilinear(DA)), ea(s.displayTrilinear(EA)),
  mQLsq(s.displaySoftMassSquared(mQl)), 
  mURsq(s.displaySoftMassSquared(mUr)), 
  mDRsq(s.displaySoftMassSquared(mDr)),
  mLLsq(s.displaySoftMassSquared(mLl)),
  mSEsq(s.displaySoftMassSquared(mEr)), 
  mH1sq(s.displayMh1Squared()),
     mH2sq(s.displayMh2Squared()), m32(s.displayGravitino()), 
     mDbarsq(s.displaymDbarsq()), mDsq(s.displaymDSq()), mGdash_1(s.displaymGdash_1()),A_lambda(s.displayA_lambda()), A_kappa(s.displayA_kappa()), mHdashsq(s.displaymHdashsq()), mHbardashsq(s.displaymHbardashsq()),  B_0(s.displayB_0()),mSsq(s.displaymSsq()), MNcsq(s.displayMNcsq()), A_top(s.displayA_top()), A_tau(s.displayA_tau()), A_b(s.displayA_b()), A_N(s.displayA_N()){
  setPars(numSoftParsEssm);
  setMu(s.displayMu()); 
  setLoops(s.displayLoops());
  setThresholds(s.displayThresholds());
}

inline SoftParsEssm::SoftParsEssm(const EssmSusy &s)
  : EssmSusy(s), mGaugino(3), ua(3, 3), da(3, 3), ea(3, 3), 
     mQLsq(3, 3), mURsq(3, 3), mDRsq(3, 3), mLLsq(3, 3), mSEsq(3, 3),
     m32(0.0), mDbarsq(3), mDsq(3),mH1sq(3), mH2sq(3), mGdash_1(0.0),A_lambda(3), A_kappa(3), mHdashsq(3), mHbardashsq(3), B_0(0.0),mSsq(3), A_N(0.0), A_top(0.0), A_b(0.0), A_tau(0.0),MNcsq(3)  {      
   
    setPars(numSoftParsEssm);
    setMu(s.displayMu()); 
    setLoops(s.displayLoops());
    setThresholds(s.displayThresholds());
}
  
inline SoftParsEssm::SoftParsEssm
(const EssmSusy & s, const DoubleVector & mG, const
 DoubleMatrix & aU, const DoubleMatrix & aD, const DoubleMatrix & aE, const
 DoubleMatrix & mQl, const DoubleMatrix & mUr, const DoubleMatrix & mDr, const
 DoubleMatrix & mLl, const DoubleMatrix & mEr, DoubleVector mH1sq,
 DoubleVector mH2sq, double mg, double mu, int l, int t, DoubleVector mDbarsq, DoubleVector mDsq, double mGd,DoubleVector mHd,DoubleVector mHdb,DoubleVector A_lambda , DoubleVector A_kap, double B_0, DoubleVector mS, double A_N,double A_top, double A_tau, double A_b, DoubleVector MNc )
  : EssmSusy(s), mGaugino(mG), ua(aU), da(aD), ea(aE),
    mQLsq(mQl), mURsq(mUr), mDRsq(mDr), mLLsq(mLl), mSEsq(mEr), 
     mH1sq(mH1sq), mH2sq(mH2sq), m32(mg), mDbarsq(mDbarsq), mDsq(mDsq), mGdash_1(mGd),A_lambda(A_lambda), A_kappa(A_kap), mHdashsq(mHd), mHbardashsq(mHdb), B_0(B_0),mSsq(mS), A_N(A_N), A_top(A_top), A_b(A_b), A_tau(A_tau),MNcsq(MNc)  {      
  
  setPars(numSoftParsEssm);
  setMu(mu);
  setLoops(l);
  setThresholds(t);
}

inline SoftParsEssm SoftParsEssm::displaySoftPars() const { return *this; }

//inline double SoftParsEssm::displayM3Squared() const { return m3sq; }

//inline double SoftParsEssm::displayMh1Squared() const { return mH1sq; }

//inline double SoftParsEssm::displayMh2Squared() const { return mH2sq; }

//Peter:: replaced with ESSM vector versions


inline DoubleVector SoftParsEssm::displayMh1Squared() const { return mH1sq; }
inline DoubleVector SoftParsEssm::displayMh2Squared() const { return mH2sq; }

//Peter:: new ESSM methods.

inline DoubleVector SoftParsEssm::displaymDbarsq() const { return mDbarsq; }
inline DoubleVector SoftParsEssm::displaymDSq() const { return mDsq; }

inline double SoftParsEssm::displaymGdash_1() const { return mGdash_1; }
inline DoubleVector SoftParsEssm::displayA_lambda() const { return A_lambda; }
inline DoubleVector SoftParsEssm::displayA_kappa() const { return A_kappa; }


inline DoubleVector SoftParsEssm::displaymHdashsq() const { return mHdashsq; }
inline DoubleVector SoftParsEssm::displaymHbardashsq() const { return mHbardashsq; }

inline double SoftParsEssm::displayB_0() const { return B_0; }
inline DoubleVector SoftParsEssm::displaymSsq() const { return mSsq; }
inline DoubleVector SoftParsEssm::displayMNcsq() const { return MNcsq; }
inline double SoftParsEssm::displayA_top() const { return A_top; }
inline double SoftParsEssm::displayA_b() const { return A_b; }
inline double SoftParsEssm::displayA_N() const { return A_N; }
inline double SoftParsEssm::displayA_tau() const { return A_tau; }

    

inline DoubleVector SoftParsEssm::displayGaugino() const { return mGaugino; }

inline double SoftParsEssm::displayGaugino(int i) const { 
  return mGaugino.display(i); 
}

inline double SoftParsEssm::displayGravitino() const { return m32; }

inline void SoftParsEssm::setGauginoMass(int i, double f) { 
  mGaugino(i) = f; 
}

//inline void SoftParsEssm::setM3Squared(double f) { m3sq = f; }
//Peter:: ESSM Changed
inline void SoftParsEssm::setMh1Squared(DoubleVector f) { mH1sq = f; }
inline void SoftParsEssm::setMh2Squared(DoubleVector f) { mH2sq = f; }

inline void SoftParsEssm::setSoftPars(SoftParsEssm const & s) { *this = s; }
inline void SoftParsEssm::setM32(double a) { m32 = a; }
inline void SoftParsEssm::setMGdash_1(double f) { mGdash_1 = f; }
inline void SoftParsEssm::setA_lambda(DoubleVector f) { A_lambda = f; }
inline void SoftParsEssm::setA_kappa(DoubleVector f) { A_kappa = f; }
 

inline void SoftParsEssm::setA_top(double f) { A_top = f; }
inline void SoftParsEssm::setA_b(double f) { A_b = f; }
inline void SoftParsEssm::setA_tau(double f) { A_tau = f; }
inline void SoftParsEssm::setA_N(double f) { A_N = f; }



inline void SoftParsEssm::setmHdashsq(DoubleVector f) { mHdashsq = f; }
 inline void SoftParsEssm::setmHbardashsq(DoubleVector f) { mHbardashsq = f; }
 inline void SoftParsEssm::setB_0(double f) { B_0 = f; }
 
 //inline void SoftParsEssm::setB0_tilde(double f) { B0_tilde = f; }
inline void SoftParsEssm::setmSsq(DoubleVector f) { mSsq = f; }
inline void SoftParsEssm::setmDsq(DoubleVector f) { mDsq = f; }
inline void SoftParsEssm::setmDbarsq(DoubleVector f) {mDbarsq = f; }

inline void SoftParsEssm::setMNcsq(DoubleVector f) { MNcsq = f; }

#endif
