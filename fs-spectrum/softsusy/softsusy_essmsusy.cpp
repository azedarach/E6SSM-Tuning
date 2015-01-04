
/** \file susy.cpp
   - Project:     SOFTSUSY 
   - Author:      Ben Allanach 
   - Manual:      hep-ph/0104145 
   - Webpage:     http://allanach.home.cern.ch/allanach/softsusy.html

   $Log: susy.cpp,v $
   Revision 1.3  2005/11/09 14:12:25  allanach
   Updated for SOFTSUSY 2.0.1 - cleaned everything up etc

   Revision 1.3  2005/08/16 13:52:27  allanach
   Corrected comment on hep-ph number

   Revision 1.20  2003/11/19 17:07:22  allanach
   Higgs vev included in reading in SUSY object

   Revision 1.19  2003/10/24 16:09:04  allanach
   Implemented running Higgs DRbar vev

   Revision 1.17  2003/10/24 10:08:01  allanach
   Changed signs in superpotential

   Revision 1.16  2003/09/12 14:13:24  allanach
   Bug-fixed 2-loop contribution to wave-function renormalisation of d_R

   Revision 1.15  2003/07/25 13:39:15  allanach
   Trapped errors properly rather than exiting

   Revision 1.14  2003/07/16 11:10:23  allanach
   Cosmetic changes and removed a static (useless one)

   Revision 1.13  2003/05/27 15:05:53  allanach
   Purely efficiency corrections: used variable rather than display() methods
   whenever possible

   Revision 1.12  2003/05/20 15:19:40  allanach
   doxygen comment style implemented

   Revision 1.11  2003/02/21 17:59:36  allanach
   Added drbar parameter class and calculation, starting to move to DRbar
   parameters in the 1-loop corrections

   Revision 1.10  2002/11/19 16:59:30  allanach
   Added routine to get quark Yukawa diagonalising matrices

   Revision 1.9  2002/09/20 15:37:10  allanach
   Adding quark-mixing routines

   Revision 1.8  2002/08/30 12:31:56  allanach
   Test RPV version that works!

   Revision 1.7  2002/07/19 15:54:07  allanach
   SOFTSUSY1.5 version

   Revision 1.5  2002/04/18 14:32:05  allanach
   Changed RGEs and anomalous dimensions to be compatible with new notation;
   started implementation of rewsb in R-parity violation

   Revision 1.3  2002/04/12 16:51:27  allanach
   Added display/set functions to work automatically

   Revision 1.4  2001/07/26 14:19:08  allanach
   Changed output format in Mssmsusy::operator<<

   Revision 1.3  2001/07/18 14:42:51  allanach
   Added proper header info
*/

#include "softsusy_essmsusy.h"
#include "softsusy_susy.h"
#define HR "---------------------------------------------------------------\n"


/*
const sBrevity & sBrevity::operator=(const sBrevity &s) {
  if (this == &s) return *this;
  dt = s.dt; ut = s.ut; et = s.et; 
  u2 = s.u2; d2 = s.d2; e2 = s.e2; 
  u2t = s.u2t; e2t = s.e2t; d2t = s.d2t; 
  gsq = s.gsq; g3 = s.g3; g4 = s.g4;
  uuT = s.uuT; ddT = s.ddT; eeT = s.eeT;
  d1 = s.d1; e1 = s.e1; u1 = s.u1;
  return *this;
}

void sBrevity::calculate(const DoubleMatrix & yu, const DoubleMatrix & yd,
			  const DoubleMatrix & ye, const DoubleVector & g) {
  static DoubleVector g1(1, 3);
  g1 = g.display();
  u1 = yu.display(); 
  d1 = yd.display(); 
  e1 = ye.display();
  dt = d1.transpose(); ut = u1.transpose(); et = e1.transpose();
  u2 = u1 * ut; d2 = d1 * dt; e2 = e1 * et; 
  u2t = ut * u1; d2t = dt * d1; e2t = et * e1;
  uuT = u2.trace(), ddT = d2.trace(), eeT = e2.trace();
  gsq = g1 * g1; g3 = gsq * g1; g4 = g3 * g1;
  }

*/


//MssmSusy::MssmSusy()
//: u(3, 3), d(3, 3), e(3, 3), g(3), smu(0.0), tanb(0.0), hVev(0.0) {
//  setPars(numSusyESSMPars);
//  setMu(0.0);
//  setLoops(2);
//  setThresholds(0);
//}



//Peter:: New ESSM versions of constructors
EssmSusy::EssmSusy()
  : u(3, 3), d(3, 3), e(3, 3), g(3), smu(0.0), tanb(0.0), hVev(0.0), lambda(3), kappa(3),MN_SUSY(3), h_N(0.0),mu_0(0.0), gdash_1(0.0), g_11(0.0) {
    setPars(numSusyESSMPars);
    setMu(0.0);
    setLoops(2);
    setThresholds(0);
}

EssmSusy::EssmSusy(const EssmSusy &s)
  : u(s.u), d(s.d), e(s.e), g(s.g), smu(s.smu), tanb(s.tanb), hVev(s.hVev),  lambda(s.lambda), kappa(s.kappa),MN_SUSY(s.MN_SUSY), h_N(s.h_N),mu_0(s.mu_0), gdash_1(s.gdash_1), g_11(s.g_11) { 
    setPars(numSusyESSMPars);
    setMu(s.displayMu()); 
    setLoops(s.displayLoops());
    setThresholds(s.displayThresholds());
}

// u = up type Yukawa couplings
// d = down type Yukawa couplings
// e = lepton Yukawa couplings
// v = SM/MSSM gauge couplings
// m = MSSM mu parameter
// tb = tan(beta)
// MU = renormalisation scale
// l = number of loops
// t = threshold accuracy
// hv = Higgs vev (no class member for S vev)
// lam = ESSM SH_1iH_2i Yukawa couplings (only lambda(3) appears in tuning measures)
// kap = ESSM SD_iDbar_i couplings (don't appear in my tuning measures)
// MN_SUSY = right-handed neutrino masses (not needed)
// h_N = Yukawa coupling to right-handed neutrino (not needed)
// mu_0 = ?
// gd = ESSM U(1) gauge coupling g_1'
// g_11 = ESSM gauge coupling mixing g_11
EssmSusy::EssmSusy(const DoubleMatrix & u, const DoubleMatrix & d, const
		     DoubleMatrix & e, const DoubleVector & v, double m,
		     double tb, double MU, int l, int t, double hv, DoubleVector lam, DoubleVector kap, DoubleVector MN_SUSY, double h_N, double mu_0, double gd, double g_11)
  : u(u), d(d), e(e), g(v), smu(m), tanb(tb), hVev(hv), lambda(lam), kappa(kap),MN_SUSY(MN_SUSY), h_N(h_N), mu_0(mu_0), gdash_1(gd), g_11(g_11) { 
    setPars(numSusyESSMPars);
    setMu(MU); 
    setLoops(l);
    setThresholds(t);
}

const EssmSusy & EssmSusy::operator=(const EssmSusy & s) {
  if (this == &s) return *this;
  u = s.u;
  d = s.d;
  e = s.e;
  smu = s.smu;
  tanb = s.tanb;
  g = s.g;
  setMu(s.displayMu());
  setLoops(s.displayLoops());
  setThresholds(s.displayThresholds());
  hVev = s.hVev;
  kappa = s.kappa;
  lambda = s.lambda;
  h_N = s.h_N;
  mu_0 = s.mu_0;
  gdash_1= s.gdash_1;
  g_11 = s.g_11;
  MN_SUSY = s.MN_SUSY;
  return *this;
}

void EssmSusy::setSomePars(const EssmSusy & s) {
  u = s.u;
  d = s.d;
  e = s.e;
  g = s.g;
}

void EssmSusy::setYukawaElement(yukawa k, int i, int j, double f) { 
  switch(k) {
  case YU: u(i, j) = f; break;
  case YD: d(i, j) = f; break;
  case YE: e(i, j) = f; break;
  default: 
    ostringstream ii;
    ii << "MssmSusy::set called with illegal " << int(k) << "\n";
    throw ii.str(); break;
  }
}

void EssmSusy::setYukawaMatrix(yukawa k, const DoubleMatrix & m) { 
  switch(k) {
  case YU: u = m; break;
  case YD: d = m; break;
  case YE: e = m; break;
  default: 
    ostringstream ii;
    ii << "MssmSusy::set called with illegal " << int(k) << "\n";
    throw ii.str(); break;
  }
}

double EssmSusy::displayYukawaElement(yukawa k, int i, int j) const {
  switch(k) {
  case YU: return u.display(i, j); break;
  case YD: return d.display(i, j); break;
  case YE: return e.display(i, j); break;
  default: 
    ostringstream ii;
    ii << "MssmSusy::display called with illegal " << int(k) << "\n";
    throw ii.str(); break;
  }
  return 0.0;
}

DoubleMatrix EssmSusy::displayYukawaMatrix(yukawa k) const {
  switch(k) {
  case YU: return u; break;
  case YD: return d; break;
  case YE: return e; break;
  default: 
    ostringstream ii;    
    ii << "MssmSusy::display called with illegal " << int(k) << "\n";
    throw ii.str(); break;
  }
}

//Peter:: edited to include new Essm parameters
//Dylan:: changed return type to const DoubleVector
const DoubleVector EssmSusy::display() const {
  DoubleVector y(numSusyESSMPars);
  int i, j, k=0;
  for (i=1; i<=3; i++)    
    for (j=1; j<=3; j++) {
      k++;
      y(k) = u.display(i, j);
      y(k+9) = d.display(i, j);
      y(k+18) = e.display(i, j);
    }
  k=27;
  for (i=1; i<=3; i++) {
    k++;
    y(k) = g.display(i);
  }
  y(31) = smu;
  y(32) = tanb;
  y(33) = hVev;
  
  //Peter:: New Essm parameters
  k= 33;
for (i=1; i<=3; i++) {
    k++;
    y(k) = lambda.display(i);
  }  

for (i=1; i<=3; i++) {
    k++;
    y(k) = kappa.display(i);
  }

 y(40) = h_N;
 y(41) = mu_0;
 y(42) = gdash_1;
 y(43) = g_11;
 y(44) = MN_SUSY.display(1);
 y(45) = MN_SUSY.display(2);
 y(46) = MN_SUSY.display(3);
  return y;
}
//Peter:: edited to include new Essm parameters

void EssmSusy::set(const DoubleVector & y) {
  int i, j, k=0;
  for (i=1; i<=3; i++)    
    for (j=1; j<=3; j++){
      k++;
      u(i, j) = y.display(k);
      d(i, j) = y.display(k+9);
      e(i, j) = y.display(k+18);
    }
  k=27;
  for (i=1; i<=3; i++) {
    k++;
    g(i) = y.display(k);
  }
  smu = y.display(31);
  tanb = y.display(32);
  hVev = y.display(33);
  k = 33;
  for (i=1; i<=3; i++) {
    k++;
    lambda(i) = y.display(k);
    }
for (i=1; i<=3; i++) {
    k++;
    kappa(i) = y.display(k);
    }

 h_N = y.display(40);
 mu_0 = y.display(41);
 gdash_1 = y.display(42);
 g_11 = y.display(43);
 k= 43;
for (i=1; i<=3; i++) {
    k++;
    MN_SUSY(i) = y.display(k);
    }

}






double EssmSusy::displayTanb() const { return tanb; }

//Peter:: new Essm display functions

double EssmSusy::displayh_N() const { return h_N; }
double EssmSusy::displaymu_0() const { return mu_0; }
double EssmSusy::displaygdash_1() const { return gdash_1; }
double EssmSusy::displayg_11() const { return g_11; }
double EssmSusy::displayMN_SUSY(int i) const { return MN_SUSY.display(i); }
DoubleVector EssmSusy::displaykappa() const{ return kappa;}
DoubleVector EssmSusy::displaylambda() const{ return lambda;}

//Dylan:: additional display function. Note neglects
// Z-Z' mixing.
double EssmSusy::displaySvev(double M_Zprime) const
{
  double delta = g_11/gdash_1;
  double Q1eff = QN_H1 + QY_H1*delta;
  double Q2eff = QN_H2 + QY_H2*delta;
  double Qseff = QN_S + QY_S*delta;
  double cSqb = 1.0/(1.0+tanb*tanb);
  double sSqb = (tanb*tanb)/(1.0+tanb*tanb);
  double s = (1.0/(gdash_1*Qseff))*sqrt(M_Zprime*M_Zprime-gdash_1*gdash_1*hVev*hVev*
					(Q1eff*Q1eff*cSqb+Q2eff*Q2eff*sSqb));
  return s;
}

//Dylan:: returns effective U(1)_N charges for Higgs and S
//fields.
double EssmSusy::displayQH1tilde() const
{
  double delta = g_11/gdash_1;
  double Q1eff = QN_H1 + QY_H1*delta;
  return Q1eff;
}

double EssmSusy::displayQH2tilde() const
{
  double delta = g_11/gdash_1;
  double Q2eff = QN_H2 + QY_H2*delta;
  return Q2eff;
}

double EssmSusy::displayQStilde() const
{
  double delta = g_11/gdash_1;
  double QSeff = QN_S + QY_S*delta;
  return QSeff;
}

DoubleVector EssmSusy::displayQtilde() const
{
  DoubleVector Qtildes(3);
  double delta = g_11/gdash_1;
  double Q1eff = QN_H1 + QY_H1*delta;
  double Q2eff = QN_H2 + QY_H2*delta;
  double QSeff = QN_S + QY_S*delta;

  Qtildes.set(1, Q1eff);
  Qtildes.set(2, Q2eff);
  Qtildes.set(3, QSeff);

  return Qtildes;
}

ostream & operator <<(ostream &left, const EssmSusy &s) {
  left << "Supersymmetric parameters at Q: " << s.displayMu() << endl;
  left << " Y^U" << s.displayYukawaMatrix(YU) << " Y^D" <<
    s.displayYukawaMatrix(YD) << " Y^E" << s.displayYukawaMatrix(YE);
  left << "higgs VEV: " << s.displayHvev() 
       << " tan beta: " << s.displayTanb() << " smu: " << s.displaySusyMu() << 
    "\n";
  left << "g1: " << s.displayGaugeCoupling(1) << " g2: " <<
    s.displayGaugeCoupling(2) << " g3: " << 
    s.displayGaugeCoupling(3) << endl; 
  left << "thresholds: " << s.displayThresholds() 
       << " #loops: " << s.displayLoops() << '\n';
  return left;
}

void EssmSusy::setSusy(const EssmSusy & s) {
  setLoops(s.displayLoops());
  setThresholds(s.displayThresholds());
  setMu(s.displayMu());
  setYukawaMatrix(YU, s.displayYukawaMatrix(YU)); 
  setYukawaMatrix(YD, s.displayYukawaMatrix(YD)); 
  setYukawaMatrix(YE, s.displayYukawaMatrix(YE)); 
  setHvev(s.displayHvev());
  setTanb(s.displayTanb());
  setSusyMu(s.displaySusyMu());
  setAllGauge(s.displayGauge());
 seth_N(s.displayh_N());
   setgash_1(s.displaygdash_1());
   setg_11(s.displayg_11());
   setmu_0(s.displaymu_0());
  setkappa(s.displaykappa());
  setlambda(s.displaylambda());
 DoubleVector MN(3);
 
 int gen;
for(gen=1; gen<4; gen++) MN(gen) = s.displayMN_SUSY(gen);
  setMN_SUSY(MN);

}

istream & operator >>(istream &left, EssmSusy &s) {
  char c[70];
  DoubleMatrix u(3, 3), d(3, 3), e(3, 3);
  double g1, g2, g3, smu, mu, tanb, hv;
  int loops, thresh;
  left >> c >> c >> c >> c >> mu;
  left >> c >> u >> c >> d >> c >> e >> c >> c >> hv;
  left >> c >> c >> tanb >> c >> smu;
  left >> c >> g1 >> c >> g2 >> c >> g3;
  left >> c >> thresh >> c >> loops;
  s.setYukawaMatrix(YU, u);
  s.setYukawaMatrix(YD, d);
  s.setYukawaMatrix(YE, e);
  s.setHvev(hv);
  s.setTanb(tanb);
  s.setGaugeCoupling(1, g1);
  s.setGaugeCoupling(2, g2);
  s.setGaugeCoupling(3, g3);
  s.setThresholds(thresh);
  s.setSusyMu(smu);
  s.setMu(mu);
  s.setLoops(loops);
  return left;
}



// Outputs derivatives (DRbar scheme) in the form of ds. a contains the
// matrices calculated that are handy for computation.
// W=  LL Y^E H1 ER + QL Y^D H1 DR + QL Y^U H2 UR + smu H2 H1
// is the superpotential. Consistent with Allanach, Dedes, Dreiner
// hep-ph/9902251 and Barger, Berger and Ohmann hep-ph/9209232, 9311269
// EXCEPT for the sign of smu, which is opposite. These equations are also
// valid for W=  - LL Y^E H1 ER - QL Y^D H1 DR + QL Y^U H2 UR + smu H2 H1, the
// new SOFTSUSY convention
EssmSusy EssmSusy::beta(sBrevity & a) const {
  // Wave function renormalisations: convention for g**(i, j) is that i is the
  // LOWER index and j the upper in our paper hep-ph/9902251
  //Peter:  In third family approximation we are only using these as numbers not matrices
   static DoubleMatrix gEE(3, 3), gLL(3, 3), gQQ(3, 3), gDD(3, 3), gUU(3, 3);
   double gH1H1, gH2H2;
 static double gEE_3, gLL_3, gQQ_3, gDD_3, gUU_3;
  double gH1H1_3=0.0, gH2H2_3=0.0;
  static DoubleVector dg(1,3);
  static double dgdash_1, dg_11;
  // keep this option in order to interface with RPVSUSY  
  //Peter:: keep this to do the gauge couplings but then ignore the matrices gEE etc and replace with  3rd gen numbers 
  anomalousDimension(gEE, gLL, gQQ, gUU, gDD, dg,dgdash_1, dg_11, gH1H1, gH2H2, a);
  
  // To keep this a const function
  const DoubleMatrix &u1 = u.display(), &d1 = d.display(), &e1 = e.display();
  
  // contain derivatives of up, down quarks and leptons
  static DoubleMatrix du(3, 3), dd(3, 3), de(3, 3); 
  //Peter:: initialise all entries to zero since we are going to only use the third generation i don't want any funny preset vvalues for the other betas screwing up the evolution
  int just_mak, ing_sure;


  for(just_mak = 1; just_mak < 4; just_mak++){
    for(ing_sure = 1; ing_sure < 4;ing_sure++){
        du(just_mak,ing_sure) = 0;
	   dd(just_mak,ing_sure) = 0; 
	   de(just_mak,ing_sure) = 0;}}

// mu parameter derivatives
  double dmu;
  // DoubleMatrix mu_tilde(3,3), dmu_tilde(3,3);

  // RGEs of SUSY parameters//Peter:: for MSSM
  /* du = u1 * (gUU + gH2H2) + gQQ * u1;
     dd = d1 * (gDD + gH1H1) + gQQ * d1;
     de = e1 * (gEE + gH1H1) + gLL * e1;*/
  //Peter:: Yukawas for ESSM
  //Peter:: note a 1/16pi^2 factor is contained in gUU etc
  DoubleVector lambda = displaylambda(); DoubleVector kappa = displaykappa();
  //double lambda = lambdan(3) ; double kappa = kappan(3);
double h_N = displayh_N();
 DoubleVector MN_SUSY(3) ;
   MN_SUSY(1)  = displayMN_SUSY(1);
  MN_SUSY(2)  = displayMN_SUSY(2);
  MN_SUSY(3)  = displayMN_SUSY(3);
  bool eta_N = 0;
  //Peter:: to implement this threshold codition for the nuetrino coupling in the rges we need to be able to access the the current scale.  To a firts approximation we may simlt ignore these terms.  
  if(displayMu() > displayMN_SUSY(3)) eta_N=1;
  else eta_N = 0;
  double dmu_0, dh_N;
  DoubleVector dlambda(3), dkappa(3), dMN_SUSY(3); 

 //Peter:: since we are using the third family approximation the matrices gEE ect need to be changed appropriately  
  double  oneO16Pisq = 1.0/(16.0*PI*PI);

 gEE_3 = oneO16Pisq * (2.0 * sqr(e1.display(3,3)) - 1.2 * sqr(displayGaugeCoupling(1)));
     gLL_3 = oneO16Pisq * (sqr(e1.display(3,3)) - (0.3 *sqr(displayGaugeCoupling (1)) + 1.5 * sqr(displayGaugeCoupling(2))));
    
    gQQ_3 =  oneO16Pisq *(sqr(d1.display(3,3)) + sqr(u1.display(3,3)) - (sqr(displayGaugeCoupling(1)) / 30.0 + 1.5 * sqr(displayGaugeCoupling(2)) + 8 *sqr(displayGaugeCoupling(3)) / 3.0));
  gUU_3 =   oneO16Pisq *(2.0 *u1.display(3,3)*u1.display(3,3)  - (8 *sqr(displayGaugeCoupling(1))  / 15.0 + 8 *sqr(displayGaugeCoupling(3))  /3.0)); 
gDD_3 = oneO16Pisq * (2.0 * sqr(d1.display(3,3)) - (2 * sqr(displayGaugeCoupling(1)) / 15.0 + 8 * sqr(displayGaugeCoupling(3)) / 3.0));
gH1H1_3 = oneO16Pisq * (3.0 * sqr(d1.display(3,3)) + sqr(e1.display(3,3)) - (0.3 * sqr(displayGaugeCoupling(1)) + 1.5 *sqr(displayGaugeCoupling(2))));

  gH2H2_3 = oneO16Pisq * (3.0 * sqr(u1.display(3,3)) - (0.3 * sqr(displayGaugeCoupling(1)) + 1.5 * sqr(displayGaugeCoupling(2))));


du(3,3) = u1.display(3,3) * (gUU_3 + gH2H2_3) + gQQ_3 * u1.display(3,3) + 1.0 / (16.0 * sqr(PI))*u1.display(3,3)*(lambda(3)*lambda(3) - 0.3*displaygdash_1()*displaygdash_1() + h_N*h_N*eta_N);

//cout << "du(3,3) =  " <<du(3,3) << endl;
//cout << "gUU_3 = " << gUU_3 << endl;
// cout << "gH2H2_3 = " << gH2H2_3 << endl;  

// cout << "sqr(u1.display(3,3)) = " << sqr(u1.display(3,3)) << endl;

// cout << "u1.display(3,3) = " << u1.display(3,3) << endl;
// cout << "gUU_3 + gH2H2_3 = " << gUU_3 + gH2H2_3 << endl;
// cout << "gQQ_3 = " << gQQ_3 << endl;
// cout << "eta_N = " << eta_N << endl;
// cout << " h_N = " << h_N << endl;
// cout << " lambda(3)    = "   << lambda(3)    << endl;  
//cout << "displaygdash_1() = "  << displaygdash_1() << endl;


  dd(3,3) = d1.display(3,3) * (gDD_3 + gH1H1_3) + gQQ_3 * d1.display(3,3) + 1.0 / (16.0 * sqr(PI))*d1.display(3,3)*(lambda(3)*lambda(3) -0.7*displaygdash_1()*displaygdash_1());
  de(3,3) = e1.display(3,3) * (gEE_3 + gH1H1_3) + gLL_3 * e1.display(3,3) + 1.0 / (16.0 * sqr(PI))*e1.display(3,3)*(lambda(3)*lambda(3) - 0.7*displaygdash_1()*displaygdash_1() + h_N*h_N*eta_N);
  dmu = smu * (gH1H1 + gH2H2);
  
  //cout << "de = " << de << endl;

  dh_N =   1.0 / (16.0 * sqr(PI))*h_N*(lambda(3)*lambda(3) + 3* u1.display(3,3)*u1.display(3,3) + sqr(e1.display(3,3)) + 2*h_N*h_N*(1 + eta_N) - 3*displayGaugeCoupling(2)*displayGaugeCoupling(2) - 0.6*displayGaugeCoupling(1)*displayGaugeCoupling(1) - 0.4*displaygdash_1()*displaygdash_1());

  /* int gen_LD, gen_rd;
  for(gen_LD = 1; gen_LD< 4; gen_LD++){
     for(gen_rd = 1; gen_rd< 4; gen_rd++){
       
       dmu_tilde(gen_LD,genrd) = mu_tilde(gen_LD, gen_rd)*(kappa(gen_LD)*kappa(gen_LD)  - (float)16/(float)3 * displayGaugeCoupling(3)*displaygaugeCoupling(3) - (float)4/(float)15*displayGaugeCoupling(1)*displayGaugeCoupling(1) -0.4*gdash_1*gdash_1);
       if(gen_rd == 3){   dmu_tilde(gen_LD,genrd) =   dmu_tilde(gen_LD,genrd) + + 2*d1(3)}
  
     }
     }*/

  dmu_0 = 1.0 / (16.0 * sqr(PI))* mu_0*(-3*displayGaugeCoupling(2)*displayGaugeCoupling(2) - 0.6*displayGaugeCoupling(1)*displayGaugeCoupling(1) - 0.4*displaygdash_1()*displaygdash_1());
  int gen;
  
//Peter:: origionally wrote this out from rges in ESSM theory and Phenomonology paper. There the charges, Q_tilde etc apear explicitly.  In the rge copy Roman has put in numerical value for them so i have updated the beta correspondingly   
  
  //for(gen=1; gen < 3; gen++){
  //dlambda(gen) = 1.0 / (16.0 * sqr(PI))*lambda(gen)*(2*lambda(gen)*lambda(gen) + 2*(lambda(1)*lambda(1) + lambda(2)*lambda(2) + lambda(3)*lambda(3)) + 3*(kappa(1)*kappa(1) + kappa(2)*kappa(2) = kappa(3)*kappa(3))  - 3*dispalyGaugeCoupling(2)*displaygaugeCoupling(2) - 0.6*displayGaugeCoupling(1)*displayGaugeCoupling(1) - 2*gdash_1*gdash_1(Qtilde_S*Qtilde_S + Qtilde_2*Qtilde_2 +  Qtilde_1*Qtilde_1));}
  
  for(gen=1; gen < 3; gen++){
    dlambda(gen) = 1.0 / (16.0 * sqr(PI))*lambda(gen)*(2*lambda(gen)*lambda(gen) + 2*(lambda(1)*lambda(1) + lambda(2)*lambda(2) + lambda(3)*lambda(3)) + 3*(kappa(1)*kappa(1) + kappa(2)*kappa(2) + kappa(3)*kappa(3))  - 3*displayGaugeCoupling(2)*displayGaugeCoupling(2) - 0.6*displayGaugeCoupling(1)*displayGaugeCoupling(1) - 1.9*displaygdash_1()*displaygdash_1()) ;
  }
  
  //dlambda(3) = 1.0 / (16.0 * sqr(PI))*lambda(3)*(2*lambda(3)*lambda(3) + 2*(lambda(1)*lambda(1) + lambda(2)*lambda(2) + lambda(3)*lambda(3)) + 3*(kappa(1)*kappa(1) + kappa(2)*kappa(2) = kappa(3)*kappa(3))  - 3*dispalyGaugeCoupling(2)*displaygaugeCoupling(2) - 0.6*displayGaugeCoupling(1)*displayGaugeCoupling(1) - 2*gdash_1*gdash_1(Qtilde_S*Qtilde_S + Qtilde_2*Qtilde_2 +  Qtilde_1*Qtilde_1) + 3* u1.transpose()*u1 + 3*d1.transpose()*d1 + .e1.transpose()*e1 );
  
  dlambda(3) = 1.0 / (16.0 * sqr(PI))*lambda(3)*(2*lambda(3)*lambda(3) + 2*(lambda(1)*lambda(1) + lambda(2)*lambda(2) + lambda(3)*lambda(3)) + 3*(kappa(1)*kappa(1) + kappa(2)*kappa(2) + kappa(3)*kappa(3))  - 3*displayGaugeCoupling(2)*displayGaugeCoupling(2) - 0.6*displayGaugeCoupling(1)*displayGaugeCoupling(1) - 1.9*displaygdash_1()*displaygdash_1() + 3* u1.display(3,3)*u1.display(3,3) + 3*d1.display(3,3)*d1.display(3,3) + e1.display(3,3)*e1.display(3,3) + h_N*h_N*eta_N);
  
  
  for(gen=1; gen < 4; gen++){
    dkappa(gen)= 1.0 / (16.0 * sqr(PI))*kappa(gen)*( 2*kappa(gen)*kappa(gen) + 2*(lambda(1)*lambda(1) + lambda(2)*lambda(2) + lambda(3)*lambda(3)) + 3*(kappa(1)*kappa(1) + kappa(2)*kappa(2) + kappa(3)*kappa(3)) - (float)16/(float)3* displayGaugeCoupling(3)*displayGaugeCoupling(3) - (float)4/(float)15*displayGaugeCoupling(1)*displayGaugeCoupling(1) - 1.9*displaygdash_1()*displaygdash_1());
  }
  
  
  //Peter:: below is unedited thus far.  

  dMN_SUSY(1) = 0;
  dMN_SUSY(2) = 0;
  dMN_SUSY(3) =  oneO16Pisq * 4*h_N*h_N*MN_SUSY(3);

  //cout << "dMN_SUSY(3) = " <<  dMN_SUSY(3) << endl;

  double cosb2 = sqr(cos(atan(tanb))), sinb2 = 1.0 - cosb2;
  
  // Following is from hep-ph/9308335
  double dt = displayTanb() * (gH1H1 - gH2H2);
  
  double dHvev = -hVev * (cosb2 * gH1H1 + sinb2 * gH2H2);
  
  // Contains all susy derivatives:
  // cout << "dgdash_1 = " << dgdash_1 << endl;
  /*
for (gen =1; gen<4; gen++){
    //cout << "dlambda(" <<gen << ") = " << dlambda(gen) << endl;

  } 
  
 for (gen =1; gen<4; gen++){
   cout << "dkappa(" <<gen << ") = " << dlambda(gen) << endl;
    
 } 
 for (gen =1; gen<4; gen++){
   cout << "dg(" <<gen << ") = " << dg(gen) << endl;
    
 } 
 

 cout << "displayLoops() = " <<  displayLoops() << endl;
  */

 if (displayLoops() > 1) { 
   //cout << "In two loops" << endl;
bool  yukawas_2 = true;
 if(yukawas_2){    
for(gen=1;gen<4; gen++){
     
    dlambda(gen) = dlambda(gen) + sqr(oneO16Pisq) *lambda(gen)*( - 2* lambda(gen)*lambda(gen)*(lambda(gen)*lambda(gen) + 2*(lambda(1)*lambda(1) + lambda(2)*lambda(2) + lambda(3)*lambda(3)) +3*(kappa(1)*kappa(1) + kappa(2)*kappa(2) + kappa(3)*kappa(3))) - 4*(sqr(sqr(lambda(1))) + sqr(sqr(lambda(2))) + sqr(sqr(lambda(3)))) - 6*(sqr(sqr(kappa(1)))+ sqr(sqr(kappa (2)))+sqr(sqr(kappa (3)))) - 2*sqr(lambda(3))*(3*sqr(u1.display(3,3)) + 3*sqr(d1.display(3,3)) + displayh_N()*displayh_N() + sqr(e1.display(3,3)))  + 16.0*sqr(displayGaugeCoupling(3))*( kappa(1)*kappa(1) + kappa(2)*kappa(2) + kappa(3)*kappa(3)) + 6.0*sqr(displayGaugeCoupling(2)  )*( lambda(1)*lambda(1) + lambda(2)*lambda(2) + lambda(3)*lambda(3))   + sqr(displayGaugeCoupling(1))*(0.8*(kappa(1)*kappa(1) + kappa(2)*kappa(2) + kappa(3)*kappa(3)) + 1.2*(lambda(1)*lambda(1) + lambda(2)*lambda(2) + lambda(3)*lambda(3))  ) + sqr(displaygdash_1())*(   2.5*sqr(lambda(gen)) - 1.8* (kappa(1)*kappa(1) + kappa(2)*kappa(2) + kappa(3)*kappa(3)) - 1.2*(lambda(1)*lambda(1) + lambda(2)*lambda(2) + lambda(3)*lambda(3))) +3*sqr(sqr(displayGaugeCoupling(2)))*5.5 + 0.6*9.9* sqr(sqr(displayGaugeCoupling(1))) + 1.9*sqr(sqr(displaygdash_1()))*10.35 + 1.8* sqr(displayGaugeCoupling(2))*sqr(displayGaugeCoupling(1)) + 1.95* sqr(displayGaugeCoupling(2))*sqr(displaygdash_1()) + 0.39*sqr(displayGaugeCoupling(1))*sqr(displaygdash_1()));
     
     dkappa(gen) = dkappa(gen) +  sqr(oneO16Pisq) *kappa(gen)*( 

- 2* kappa(gen)*kappa(gen)*(kappa(gen)*kappa(gen) + 2*(lambda(1)*lambda(1) + lambda(2)*lambda(2) + lambda(3)*lambda(3)) +3*(kappa(1)*kappa(1) + kappa(2)*kappa(2) + kappa(3)*kappa(3))) - 4*(sqr(sqr(lambda(1))) + sqr(sqr(lambda(2))) + sqr(sqr(lambda(3)))) - 6*(sqr(sqr(kappa (1)))+ sqr(sqr(kappa (2)))+sqr(sqr(kappa (3)))) - 2*sqr(lambda(3))*(3*sqr(u1.display(3,3)) + 3*sqr(d1.display(3,3)) + displayh_N()*displayh_N() + sqr(e1.display(3,3))) + 16*sqr(displayGaugeCoupling(3))*(kappa(1)*kappa(1) + kappa(2)*kappa(2) + kappa(3)*kappa(3)) + 6* sqr(displayGaugeCoupling(2))*(lambda(1)*lambda(1) + lambda(2)*lambda(2) + lambda(3)*lambda(3)) +   sqr(displayGaugeCoupling(1))*(0.8*(kappa(1)*kappa(1) + kappa(2)*kappa(2) + kappa(3)*kappa(3)) +1.2*(lambda(1)*lambda(1) + lambda(2)*lambda(2) + lambda(3)*lambda(3))) + sqr(displaygdash_1())*(2.5*sqr(kappa(gen))- 1.8*(kappa(1)*kappa(1) + kappa(2)*kappa(2) + kappa(3)*kappa(3)) - 1.2* (lambda(1)*lambda(1) + lambda(2)*lambda(2) + lambda(3)*lambda(3))) +sqr(sqr(displayGaugeCoupling(3)))*(float)128.0/(float)9.0 +sqr(sqr(displayGaugeCoupling(1)))*584.0/225.0 + sqr(sqr(displaygdash_1()))*19.665 + sqr(displayGaugeCoupling(3))*sqr(displayGaugeCoupling(1))*64.0/45.0 + sqr(displayGaugeCoupling(3))*sqr(displaygdash_1())*52.0/15.0 +13.0/75.0*sqr(displayGaugeCoupling(1))*sqr(displaygdash_1()));

}

   dlambda(3) = dlambda(3) +    sqr(oneO16Pisq) *lambda(3)*(-sqr(lambda(3))*(3*sqr(u1.display(3,3)) + 3*sqr(d1.display(3,3)) + displayh_N()*displayh_N() + sqr(e1.display(3,3))) - 9*sqr(sqr( u1.display(3,3))) - 9*sqr(sqr( d1.display(3,3))) - 3*sqr(sqr( e1.display(3,3)))  - 3* sqr(sqr(displayh_N())) - 6*sqr( u1.display(3,3))*sqr( d1.display(3,3)) - 2*sqr(e1.display(3,3))*sqr(displayh_N()) +   16*sqr(displayGaugeCoupling(3))*(sqr(u1.display(3,3)) +sqr( d1.display(3,3))) + sqr(displayGaugeCoupling(1))*(0.8*sqr(u1.display(3,3)) -0.4*sqr( d1.display(3,3)) + 1.2*sqr( e1.display(3,3))) + sqr(displaygdash_1())*(-0.3*sqr(u1.display(3,3)) -0.2*sqr( d1.display(3,3)) - 0.2* sqr( e1.display(3,3))) );
   //cout << "changed AND IN RIGHT CODE!" << endl;

   //cout << "du.display(3,3) = " << du.display(3,3) << endl;

   du(3,3) = du.display(3,3) +  sqr(oneO16Pisq)*u1.display(3,3) *( 
					     -22.0*sqr((sqr(u1.display(3,3)))) - 5* sqr((sqr(d1.display(3,3)))) - 5*sqr(u1.display(3,3))*sqr(d1.display(3,3)) - 3*sqr(u1.display(3,3))*sqr(displayh_N()) -3*sqr(d1.display(3,3))*sqr(e1.display(3,3)) -sqr(e1.display(3,3))*sqr(displayh_N()) - 3*(sqr(sqr(displayh_N()))) - sqr(lambda(3))*(sqr(lambda(3)) + 3*sqr(u1.display(3,3)) + 4* sqr(d1.display(3,3)) + sqr(e1.display(3,3)) + 2*(lambda(1)*lambda(1) + lambda(2)*lambda(2) + lambda(3)*lambda(3)) +3*(kappa(1)*kappa(1) + kappa(2)*kappa(2) + kappa(3)*kappa(3))) + 16*sqr(displayGaugeCoupling(3))*sqr(u1.display(3,3)) + 6*sqr(displayGaugeCoupling(2))*sqr(u1.display(3,3)) +sqr(displayGaugeCoupling(1))*(1.2* sqr(u1.display(3,3)) + 0.4*sqr(d1.display(3,3))) + sqr(displaygdash_1())*(1.5*lambda(3)*lambda(3) + 0.3*sqr(u1.display(3,3)) + 0.6*sqr(d1.display(3,3))) + sqr(sqr(displayGaugeCoupling(3)))*128.0/9.0 +  sqr(sqr(displayGaugeCoupling(2)))*16.5 + sqr(sqr(displayGaugeCoupling(1)))*3913.0/450.0 +  sqr(sqr(displaygdash_1()))*2.865 + 8*sqr(displayGaugeCoupling(3))*sqr(displayGaugeCoupling(2)) + 136.0/45.0*sqr(displayGaugeCoupling(3))*sqr(displayGaugeCoupling(1)) + 8.0/15.0*sqr(displayGaugeCoupling(3))*sqr(displaygdash_1()) + sqr(displayGaugeCoupling(2))*sqr(displayGaugeCoupling(1)) + 0.75*sqr(displayGaugeCoupling(2))*sqr(displaygdash_1()) + 53.0/300.0* sqr(displayGaugeCoupling(1))*sqr(displaygdash_1()) );
  

   //cout << " sqr(oneO16Pisq) *( -22.0*sqr((sqr(u1.display(3,3)))) - 5* sqr((sqr(d1.display(3,3)))) - 5*sqr(u1.display(3,3))*sqr(d1.display(3,3)) - 3*sqr(u1.display(3,3))*sqr(displayh_N()) -3*sqr(d1.display(3,3))*sqr(e1.display(3,3)) -sqr(e1.display(3,3))*sqr(displayh_N()) - 3*(sqr(sqr(displayh_N()))) = "	<< sqr(oneO16Pisq) *( -22.0*sqr((sqr(u1.display(3,3)))) - 5* sqr((sqr(d1.display(3,3)))) - 5*sqr(u1.display(3,3))*sqr(d1.display(3,3)) - 3*sqr(u1.display(3,3))*sqr(displayh_N()) -3*sqr(d1.display(3,3))*sqr(e1.display(3,3)) -sqr(e1.display(3,3))*sqr(displayh_N()) - 3*(sqr(sqr(displayh_N())))) << endl;


   //cout << "After 2lp du(3,3) = " << du.display(3,3) << endl;
   //Peter:: seems very similar to  du, in particular many of the coefficients on the g^4 terms are exactly the same check this against MSSM then ask Roman

dd(3,3) =  dd.display(3,3) +  sqr(oneO16Pisq) *d1.display(3,3)*( 
								-5.0*sqr((sqr(u1.display(3,3)))) - 22.0* sqr((sqr(d1.display(3,3)))) - 5*sqr(u1.display(3,3))*sqr(d1.display(3,3)-3*sqr(d1.display(3,3))*sqr(e1.display(3,3))) - sqr(u1.display(3,3))*sqr(displayh_N())  -sqr(e1.display(3,3))*sqr(displayh_N()) - 3*(sqr(sqr(e1.display(3,3)))) - sqr(lambda(3))*(sqr(lambda(3)) + 4*sqr(u1.display(3,3)) + 3* sqr(d1.display(3,3)) + sqr(displayh_N()) + 2*(lambda(1)*lambda(1) + lambda(2)*lambda(2) + lambda(3)*lambda(3)) +3*(kappa(1)*kappa(1) + kappa(2)*kappa(2) + kappa(3)*kappa(3))) 

+ 16*sqr(displayGaugeCoupling(3))*sqr(d1.display(3,3)) + 6*sqr(displayGaugeCoupling(2))*sqr(d1.display(3,3)) +sqr(displayGaugeCoupling(1))*(0.8*sqr(u1.display(3,3)) + 0.4*sqr(d1.display(3,3))+1.2*sqr(e1.display(3,3)) ) + sqr(displaygdash_1())*(lambda(3)*lambda(3) + 0.2*sqr(u1.display(3,3)) + sqr(d1.display(3,3)) - 0.2*sqr(e1.display(3,3))) 

+ sqr(sqr(displayGaugeCoupling(3)))*128.0/9.0 +  sqr(sqr(displayGaugeCoupling(2)))*16.5 + sqr(sqr(displayGaugeCoupling(1)))*413.0/90.0 +  sqr(sqr(displaygdash_1()))*6.825 + 8*sqr(displayGaugeCoupling(3))*sqr(displayGaugeCoupling(2)) + 8.0/9.0*sqr(displayGaugeCoupling(3))*sqr(displayGaugeCoupling(1)) 

+ 4.0/3.0*sqr(displayGaugeCoupling(3))*sqr(displaygdash_1()) + sqr(displayGaugeCoupling(2))*sqr(displayGaugeCoupling(1)) + 1.5*sqr(displayGaugeCoupling(2))*sqr(displaygdash_1()) + 49.0/150.0* sqr(displayGaugeCoupling(1))*sqr(displaygdash_1()) );

		      
 de(3,3) =  de.display(3,3) +  sqr(oneO16Pisq) *e1.display(3,3)*( 
			     - 9.0* sqr((sqr(d1.display(3,3)))) - 3*sqr(u1.display(3,3))*sqr(d1.display(3,3)-9*sqr(d1.display(3,3))*sqr(e1.display(3,3))) - 3*sqr(u1.display(3,3))*sqr(displayh_N())  -3* sqr(e1.display(3,3))*sqr(displayh_N()) - 10.0*(sqr(sqr(e1.display(3,3)))) -3*sqr(displayh_N())*sqr((e1.display(3,3))) - sqr(lambda(3))*(sqr(lambda(3)) + 3*sqr(u1.display(3,3)) + 3* sqr(e1.display(3,3)) + 2*sqr(displayh_N()) + 2*(lambda(1)*lambda(1) + lambda(2)*lambda(2) + lambda(3)*lambda(3)) +3*(kappa(1)*kappa(1) + kappa(2)*kappa(2) + kappa(3)*kappa(3))) + 16*sqr(displayGaugeCoupling(3))*sqr(d1.display(3,3)) + 6*sqr(displayGaugeCoupling(2))*sqr(e1.display(3,3)) +sqr(displayGaugeCoupling(1))*(- 0.4*sqr(d1.display(3,3))+1.2*sqr(e1.display(3,3)) ) + sqr(displaygdash_1())*(lambda(3)*lambda(3) - 0.2*sqr(d1.display(3,3)) + 1.3*sqr(e1.display(3,3)) )  +  sqr(sqr(displayGaugeCoupling(2)))*16.5 + sqr(sqr(displayGaugeCoupling(1)))*18.9 +  sqr(sqr(displaygdash_1()))*6.825  + 1.8*sqr(displayGaugeCoupling(2))*sqr(displayGaugeCoupling(1)) + 39.0/20.0*sqr(displayGaugeCoupling(2))*sqr(displaygdash_1()) + 51.0/100.0* sqr(displayGaugeCoupling(1))*sqr(displaygdash_1()) );
 //cout << "I changed them. I chaneged them. I really really changed them!!"   <<endl;
 
 dh_N =  dh_N +  sqr(oneO16Pisq) *displayh_N()*( - 9.0* sqr((sqr(u1.display(3,3)))) - 3*sqr(u1.display(3,3))*sqr(d1.display(3,3)-3*sqr(d1.display(3,3))*sqr(e1.display(3,3))) - 9*sqr(u1.display(3,3))*sqr(displayh_N())  -3* sqr(e1.display(3,3))*sqr(displayh_N()) - 3.0*(sqr(sqr(e1.display(3,3))))-10*sqr(sqr(displayh_N()))  -3*sqr(displayh_N())*sqr((e1.display(3,3))) - sqr(lambda(3))*(sqr(lambda(3)) + 3*sqr(u1.display(3,3)) + 2* sqr(e1.display(3,3)) + 3*sqr(displayh_N()) + 2*(lambda(1)*lambda(1) + lambda(2)*lambda(2) + lambda(3)*lambda(3)) +3*(kappa(1)*kappa(1) + kappa(2)*kappa(2) + kappa(3)*kappa(3))) + 16*sqr(displayGaugeCoupling(3))*sqr(u1.display(3,3)) + 6*sqr(displayGaugeCoupling(2))*sqr(displayh_N()) +sqr(displayGaugeCoupling(1))*(0.8*sqr(u1.display(3,3))+1.2*sqr(e1.display(3,3))+ 1.2*sqr(displayh_N()) ) + sqr(displaygdash_1())*(1.5*lambda(3)*lambda(3) - 0.3*sqr(d1.display(3,3)) + 0.3*sqr(e1.display(3,3)) + 0.8*sqr(displayh_N())  )  +  sqr(sqr(displayGaugeCoupling(2)))*16.5 + sqr(sqr(displayGaugeCoupling(1)))*5.94 +  sqr(sqr(displaygdash_1()))*3.84  + 1.8*sqr(displayGaugeCoupling(2))*sqr(displayGaugeCoupling(1)) + 1.2*sqr(displayGaugeCoupling(2))*sqr(displaygdash_1()) + 0.24* sqr(displayGaugeCoupling(1))*sqr(displaygdash_1()) );
 }
 bool TwoloopGauges = true;
 if(TwoloopGauges){
   //cout << " dg(1) = " << dg(1) << endl;
 dg(1) = dg(1) + sqr(oneO16Pisq)*sqr(displayGaugeCoupling(1))*displayGaugeCoupling(1)*( 24*sqr(displayGaugeCoupling(3)) + 10.8* sqr(displayGaugeCoupling(2)) + 234.0/25.0*sqr(displayGaugeCoupling(1)) + 81.0/25.0*sqr(displaygdash_1()) - 26.0/5.0*sqr(u1.display(3,3)) - 14.0/5.0*sqr(d1.display(3,3)) - 18.0/5.0*sqr(e1.display(3,3)) -  6.0/5.0 *((lambda(1)*lambda(1) + lambda(2)*lambda(2) + lambda(3)*lambda(3))) - 4.0/5.0*(kappa(1)*kappa(1) + kappa(2)*kappa(2) + kappa(3)*kappa(3) )- 6.0/5.0*sqr(displayh_N())*eta_N );  
				  

 //cout << " dg(1) = " << dg(1) << endl;
 dg(2) = dg(2) + sqr(oneO16Pisq)*sqr(displayGaugeCoupling(2))*displayGaugeCoupling(2)*(24*sqr(displayGaugeCoupling(3)) + 46.0* sqr(displayGaugeCoupling(2)) + 18.0/5.0*sqr(displayGaugeCoupling(1)) + 17.0/5.0*sqr(displaygdash_1()) - 6.0*sqr(u1.display(3,3)) - 6.0*sqr(d1.display(3,3)) - 2.0*sqr(e1.display(3,3)) -  2.0*((lambda(1)*lambda(1) + lambda(2)*lambda(2) + lambda(3)*lambda(3))) - 2.0*sqr(displayh_N())*eta_N );
				  

 dg(3) = dg(3) + sqr(oneO16Pisq)*sqr(displayGaugeCoupling(3))*displayGaugeCoupling(3)*((34*3 - 54)*sqr(displayGaugeCoupling(3)) + 9.0* sqr(displayGaugeCoupling(2)) + 3.0*sqr(displayGaugeCoupling(1)) + 3.0*sqr(displaygdash_1()) - 4.0*sqr(u1.display(3,3)) - 4.0*sqr(d1.display(3,3))-2.0*(kappa(1)*kappa(1) + kappa(2)*kappa(2) + kappa(3)*kappa(3) ) );  				  

dgdash_1 = dgdash_1 + sqr(oneO16Pisq)*sqr(displaygdash_1())*displaygdash_1()*(24*sqr(displayGaugeCoupling(3)) + 51.0/5.0* sqr(displayGaugeCoupling(2)) + 81.0/25.0*sqr(displayGaugeCoupling(1)) + 229.0/25.0*sqr(displaygdash_1()) - 9.0/5.0*sqr(u1.display(3,3)) - 21.0/5.0*sqr(d1.display(3,3)) - 7.0/5.0*sqr(e1.display(3,3)) -  19.0/5.0 *((lambda(1)*lambda(1) + lambda(2)*lambda(2) + lambda(3)*lambda(3))) - 5.7*(kappa(1)*kappa(1) + kappa(2)*kappa(2) + kappa(3)*kappa(3) )- 4.0/5.0*sqr(displayh_N())*eta_N ) +sqr(oneO16Pisq)* displaygdash_1()*displayg_11()*displayg_11()*(24*sqr(displayGaugeCoupling(3)) + 10.8* sqr(displayGaugeCoupling(2)) + 234.0/25.0*sqr(displayGaugeCoupling(1)) + 81.0/25.0*sqr(displaygdash_1()) - 26.0/5.0*sqr(u1.display(3,3)) - 14.0/5.0*sqr(d1.display(3,3)) - 18.0/5.0*sqr(e1.display(3,3)) -  6.0/5.0 *((lambda(1)*lambda(1) + lambda(2)*lambda(2) + lambda(3)*lambda(3))) - 4.0/5.0*(kappa(1)*kappa(1) + kappa(2)*kappa(2) + kappa(3)*kappa(3) )- 6.0/5.0*sqr(displayh_N())*eta_N ) ;  				  


 dg_11 = dg_11 +   sqr(oneO16Pisq)*2*sqr(displayGaugeCoupling(1))*displayg_11()*( 24*sqr(displayGaugeCoupling(3)) + 10.8* sqr(displayGaugeCoupling(2)) + 234.0/25.0*sqr(displayGaugeCoupling(1)) + 81.0/25.0*sqr(displaygdash_1()) - 26.0/5.0*sqr(u1.display(3,3)) - 14.0/5.0*sqr(d1.display(3,3)) - 18.0/5.0*sqr(e1.display(3,3)) -  6.0/5.0 *((lambda(1)*lambda(1) + lambda(2)*lambda(2) + lambda(3)*lambda(3))) - 4.0/5.0*(kappa(1)*kappa(1) + kappa(2)*kappa(2) + kappa(3)*kappa(3) )- 6.0/5.0*sqr(displayh_N())*eta_N ) +   sqr(oneO16Pisq)*displayg_11()*sqr(displaygdash_1())* (24*sqr(displayGaugeCoupling(3)) + 51.0/5.0* sqr(displayGaugeCoupling(2)) + 81.0/25.0*sqr(displayGaugeCoupling(1)) + 229.0/25.0*sqr(displaygdash_1()) - 9.0/5.0*sqr(u1.display(3,3)) - 21.0/5.0*sqr(d1.display(3,3)) - 7.0/5.0*sqr(e1.display(3,3)) -  19.0/5.0 *((lambda(1)*lambda(1) + lambda(2)*lambda(2) + lambda(3)*lambda(3))) - 5.7*(kappa(1)*kappa(1) + kappa(2)*kappa(2) + kappa(3)*kappa(3) )- 4.0/5.0*sqr(displayh_N())*eta_N ) +  sqr(oneO16Pisq)*sqr(displayg_11())*displayg_11()*(24*sqr(displayGaugeCoupling(3)) + 10.8* sqr(displayGaugeCoupling(2)) + 234.0/25.0*sqr(displayGaugeCoupling(1)) + 81.0/25.0*sqr(displaygdash_1()) - 26.0/5.0*sqr(u1.display(3,3)) - 14.0/5.0*sqr(d1.display(3,3)) - 18.0/5.0*sqr(e1.display(3,3)) -  6.0/5.0 *((lambda(1)*lambda(1) + lambda(2)*lambda(2) + lambda(3)*lambda(3))) - 4.0/5.0*(kappa(1)*kappa(1) + kappa(2)*kappa(2) + kappa(3)*kappa(3) )- 6.0/5.0*sqr(displayh_N())*eta_N );


}
 }


  EssmSusy ds(du, dd, de, dg, dmu, dt, displayMu(), displayLoops(),
	      displayThresholds(), dHvev, dlambda, dkappa, dMN_SUSY, dh_N, dmu_0, dgdash_1, dg_11 ); 
  
  // cout << "ds = " << ds.display() << endl; 
  return ds;
}

void setESSMBetas(DoubleMatrix & babBeta, DoubleVector &cuBeta, DoubleVector
	       & cdBeta, DoubleVector & ceBeta, DoubleVector & bBeta) {
  // 1 loop gauge beta fns //Peter:: for MSSM
  //  bBeta(1) = 33.0 / 5.0; bBeta(2) = 1.0; bBeta(3) = -3.0; 
  
 //Peter:: 1 loop gauge beta fns  for ESSM
  bBeta(1) = 48 / 5.0; bBeta(2) = 4.0; bBeta(3) = 0; 
  

  // Extra sleptons included in vectorlike rep.s: 3 ER + 2 LL
  //#ifdef SLEPTONS
  //bBeta(1) = bBeta(1) + (3.0 * 1.2 + 0.6 * 2.0);
  //bBeta(2) = bBeta(2) + 2.0;
  //#endif
  
  //Peter:: ALL of below will need to be changed for ESSMM

  // Next come the two loop MSSM constants for gauge beta fns
  
  //Peter:: MssM ones below:
  /*
babBeta(1, 1) = 199.0 / 25.0; babBeta(1, 2) = 27.0 / 5.0; 
  babBeta(1, 3) = 88.0 / 5.0; 
  babBeta(2, 1) = 9.0 / 5.0;    babBeta(2, 2) = 25.0;       
  babBeta(2, 3) = 24.0;
  babBeta(3, 1) = 11.0 / 5.0;   babBeta(3, 2) = 9.0;        
  babBeta(3, 3) = 14.0;
  cuBeta(1) = 26.0 / 5.0; cuBeta(2) = 6.0; cuBeta(3) = 4.0;
  cdBeta(1) = 14.0 / 5.0; cdBeta(2) = 6.0; cdBeta(3) = 4.0;
  ceBeta(1) = 18.0 / 5.0; ceBeta(2) = 2.0; ceBeta(3) = 0.0;
  */

  //ESSM ones
babBeta(1, 1) = 199.0 / 25.0; babBeta(1, 2) = 27.0 / 5.0; 
  babBeta(1, 3) = 88.0 / 5.0; 
  babBeta(2, 1) = 9.0 / 5.0;    babBeta(2, 2) = 25.0;       
  babBeta(2, 3) = 24.0;
  babBeta(3, 1) = 11.0 / 5.0;   babBeta(3, 2) = 9.0;        
  babBeta(3, 3) = 14.0;
  cuBeta(1) = 26.0 / 5.0; cuBeta(2) = 6.0; cuBeta(3) = 4.0;
  cdBeta(1) = 14.0 / 5.0; cdBeta(2) = 6.0; cdBeta(3) = 4.0;
  ceBeta(1) = 18.0 / 5.0; ceBeta(2) = 2.0; ceBeta(3) = 0.0;

}

// outputs one-loop anomlous dimensions gii given matrix inputs
// Note that we use the convention (for matrices in terms of gamma's)
// gamma^Li_Lj = M_ij for LH fields and
// gamma^Rj_Ri = M_ij for RH fields (since they are really the complex
// conjugates of the RH fields): CHECKED 23/5/02

void EssmSusy::getOneLpAnom(DoubleMatrix & gEE, DoubleMatrix & gLL,
				DoubleMatrix & gQQ, DoubleMatrix & gDD,
				DoubleMatrix & gUU, double & gH1H1, double &
				gH2H2, sBrevity & a) const {
  // For calculational brevity
  DoubleMatrix &u2=a.u2, &d2=a.d2, &e2=a.e2, &u2t=a.u2t, &d2t=a.d2t,
    &e2t=a.e2t;      
  double &uuT = a.uuT, &ddT = a.ddT, &eeT = a.eeT;
  DoubleVector &gsq=a.gsq;
  
  static const double oneO16Pisq = 1.0 / (16.0 * sqr(PI));
  
  

  gEE = oneO16Pisq * (2.0 * e2t - 1.2 * gsq(1));
  gLL = oneO16Pisq * (e2 - (0.3 * gsq(1) + 1.5 * gsq(2)));
  gQQ = oneO16Pisq * (d2 + u2 - (gsq(1) / 30.0 + 1.5 * gsq(2) + 8 *
				    gsq(3) / 3.0));
  gUU = oneO16Pisq * (2.0 * u2t - (8 * gsq(1) / 15.0 + 8 * gsq(3) /
				      3.0)); 
  gDD = oneO16Pisq * (2.0 * d2t - 
			 (2 * gsq(1) / 15.0 + 8 * gsq(3) / 3.0));
  gH1H1 = oneO16Pisq * (3.0 * ddT + eeT - (0.3 * gsq(1) + 1.5 *
					      gsq(2)));
  gH2H2 = oneO16Pisq * (3.0 * uuT - (0.3 * gsq(1) + 1.5 * gsq(2)));
}

// adds two-loop anomalous dimension contribution to gii given matrix inputs
// g^Li_Lj = m_{ij} for LH fields
// g^Ei_Ej = m_{ji} for RH fields CHECKED: 23/5/02
void EssmSusy::getTwoLpAnom(DoubleMatrix & gEE, DoubleMatrix & gLL,
				DoubleMatrix & gQQ, DoubleMatrix & gDD,
				DoubleMatrix & gUU, double & gH1H1, double &
				gH2H2, sBrevity & a) const {
  // For calculational brevity
  DoubleMatrix &dt=a.dt, &ut=a.ut, &u2=a.u2, &d2=a.d2, &e2=a.e2,
    &u2t=a.u2t, &d2t=a.d2t, &e2t=a.e2t, &u1=a.u1, &d1=a.d1;      
  double &uuT = a.uuT, &ddT = a.ddT, &eeT = a.eeT;
  DoubleVector &gsq=a.gsq, &g4=a.g4;
  
  // Everything gets the (1/16pi^2)^2 factor at the bottom
  DoubleMatrix ee(3, 3), ll(3, 3), qq(3, 3), dd(3, 3), uu(3, 3); 
  
  // Two-loop pure gauge anom dimensions
  double h1h1 = (3.75 * g4(2) + 2.07 * g4(1) + 0.9 * gsq(2) * gsq(1));
  double h2h2 = h1h1;
  ll = h1h1;
  ee = (234. * g4(1) / 25.0);
  qq = (-8.0 * g4(3) / 9.0 + 3.75 * g4(2) + 199.0 * g4(1) / 900.0 + 8.0 *
	gsq(3) * gsq(2) + 8 * gsq(3) * gsq(1) / 45.0 + 0.1 * gsq(1) *
	gsq(2));
  dd = (-8.0 * g4(3) / 9.0 + 202.0 / 225.0 * g4(1) + 32.0 / 45.0 *
	gsq(3) * gsq(1));
  uu = (-8.0 * g4(3) / 9.0 + 856.0 / 225.0 * g4(1) + 128.0 / 45.0 *
	gsq(3) * gsq(1));

  ll = ll + 1.2 * gsq(1) * e2;
  ee = ee + (6.0 * gsq(2) - 1.2 * gsq(1)) * e2t;
  qq = qq + 0.4 * gsq(1) * (d2 + 2.0 * u2);
  dd = dd + (6.0 * gsq(2) + 0.4 * gsq(1)) * d2t;
  uu = uu + (6.0 * gsq(2) - 0.4 * gsq(1)) * u2t;

  h1h1 = h1h1 + (16 * gsq(3) - 0.4 * gsq(1)) * ddT + 1.2 * gsq(1) *
    eeT; 
  h2h2 = h2h2 + (16 * gsq(3) + 0.8 * gsq(1)) * uuT;

  // Two-loop pure Yukawa contributions
  double s = (eeT + 3.0 * ddT), t = (d2 * u2).trace();

  ll = ll - (2.0 * e2 * e2 + s * e2);
  ee = ee - (2.0 * e2t * e2t + 2.0 * s * e2t);
  qq = qq - (2.0 * d2 * d2 + d2 * s + 2.0 * u2 * u2 + 3.0 * uuT * u2);
  dd = dd - (2.0 * d2t * d2t + 2.0 * (dt * u2 * d1) + 2 * s * d2t);
  uu = uu - (2.0 * u2t * u2t + 2.0 * (ut * d2 * u1) + 6.0 * uuT * u2t);
  h1h1 = h1h1 - (3.0 * (e2 * e2).trace() + 9.0 * (d2t * d2t).trace() +
		 3.0 * t);
  h2h2 = h2h2 - (9.0 * (u2 * u2).trace() + 3.0 * t);

  const static double twolp = 4.010149318236068e-5; // 1/(16 pi^2)^2
  
  gLL = gLL + twolp * ll;
  gEE = gEE + twolp * ee;
  gQQ = gQQ + twolp * qq;
  gDD = gDD + twolp * dd;
  gUU = gUU + twolp * uu;
  gH1H1 = gH1H1 + twolp * h1h1;
  gH2H2 = gH2H2 + twolp * h2h2;
}

// Outputs wave function renormalisation for SUSY parameters and gauge beta
// functions up to 2 loops. 
void EssmSusy::anomalousDimension(DoubleMatrix & gEE, DoubleMatrix & gLL,
				    DoubleMatrix & gQQ, DoubleMatrix & gUU,
				    DoubleMatrix & gDD, DoubleVector & dg, 
				    double & dgdash_1, double & dg_11, double & gH1H1, double & gH2H2, 
				    sBrevity & a)  const {
  // Constants for gauge running
  static DoubleVector bBeta(3), cuBeta(3), cdBeta(3), ceBeta(3);
  static DoubleMatrix babBeta(3, 3);
  if (bBeta(1) < 1.0e-5) // Constants not set yet
    setESSMBetas(babBeta, cuBeta, cdBeta, ceBeta, bBeta);
  
  //  sBrevity a contains all of the shortcutted matrices etc;
  a.calculate(u.display(), d.display(), e.display(), g.display());
  
  // For calculational brevity
double &uuT = a.uuT, &ddT = a.ddT, &eeT = a.eeT;
DoubleVector &gsq = a.gsq, &g3 = a.g3;

// 1 loop contributions: 
if (displayLoops() > 0) {
static const double oneO16Pisq = 1.0 / (16.0 * sqr(PI)); 

getOneLpAnom(gEE, gLL, gQQ, gDD, gUU, gH1H1, gH2H2, a);
dg = oneO16Pisq * g3 * bBeta;  

if (displayLoops() > 1) {
dg_11 =  oneO16Pisq *displayGaugeCoupling(1)*(2*displayGaugeCoupling(1)*displaygdash_1()*(- sqrt(6.0)/5.0 ) + 2*displayGaugeCoupling(1)*displayg_11()*( 48.0/5.0 ))  +  oneO16Pisq *displayg_11()*(sqr(displaygdash_1())*(47.0/5.0  ) 
+ 2*displaygdash_1()*displayg_11()*(-sqrt(6.0)/5.0  )  + sqr(displayg_11())*(48.0/5.0  ));

dgdash_1 = oneO16Pisq*(float)47/(float)5*pow(displaygdash_1(), 3) + oneO16Pisq*sqr(displaygdash_1())*displayg_11()*( - 2*sqrt(6.0)/5.0) +oneO16Pisq* displaygdash_1()*sqr(displayg_11())*( 48.0/5.0 ) ;


} }
  
 bool Petesays = false;
 if(Petesays){  
if (displayLoops() > 1) {
    cout << "WARNING SHOULDN'T BE IN HERE THIS IS TWO LOOPS! " << endl;
    getTwoLpAnom(gEE, gLL, gQQ, gDD, gUU, gH1H1, gH2H2, a);

    const static double twolp = 4.010149318236068e-5; 


    //Peter:: MSSM
    //    dg = dg + g3 * (babBeta * gsq - cuBeta * uuT - cdBeta *   ddT - ceBeta * eeT) * twolp;  
    //ESSM
    //    dg(1) = dg(1) + twolp*(24*sqr(displayGaugeCoupling(3)) + 54.0/5.0* sqr(displayGaugeCoupling(2)) + 234.0/25.0*sqr(displayGaugeCoupling(1)) + 231.0/25.0*sqr(displaygdash_1()) - 26.0/5.0*u1.display(3,3) - 14.0/5.0*d1.display(3,3) - 18.0/5.0*e1.display(3,3) -  6.0/5.0 ((lambda(1)*lambda(1) + lambda(2)*lambda(2) + lambda(3)*lambda(3)) - 4.0/5.0*(kappa(1)*kappa(1) + kappa(2)*kappa(2) + kappa(3)*kappa(3) ) ); 

}
 }
}

// Outputs derivatives vector y[n] for SUSY parameters: interfaces to
// integration routines
DoubleVector EssmSusy::beta() const {
  static sBrevity a;
  // cout <<"In susy beta(() " << endl;
  // calculate the derivatives
  static EssmSusy ds;
  
  ds = beta(a);

  return ds.display(); // convert to a long vector
}



// r should be valid AT mt
void EssmSusy::setDiagYukawas(const QedQcd & r, double vev) {

  double v1, v2; // Higgs VEVs

  v1 = vev * cos(atan(displayTanb()));
  v2 = vev * sin(atan(displayTanb()));

  DoubleMatrix u1(3, 3), d1(3, 3), e1(3, 3);
  
  double invv2, invv1; 
  invv2 = 1.0 / v2; invv1 = 1.0 / v1;
  u1(1, 1) = r.displayMass(mUp) * invv2;
  u1(2, 2) = r.displayMass(mCharm) * invv2;
  u1(3, 3) = r.displayMass(mTop) * invv2;
  
  d1(1, 1) = r.displayMass(mDown) * invv1;
  d1(2, 2) = r.displayMass(mStrange) * invv1;
  d1(3, 3) = r.displayMass(mBottom) * invv1;
  e1(1, 1) = r.displayMass(mElectron) * invv1;
  e1(2, 2) = r.displayMass(mMuon) * invv1;
  e1(3, 3) = r.displayMass(mTau) * invv1;
  
  setYukawaMatrix(YU, u1); 
  setYukawaMatrix(YD, d1);
  setYukawaMatrix(YE, e1);
}
  
// mix = 0 for all mixing in downs...at present this is the only possibility.
// Takes diagonal quark Yukawa matrices and mixes them up according to the CKM
// matrix assuming:
// mix=2, all mixing is in down sector
// mix=1, all mixing is in up sector
void EssmSusy::quarkMixing(const DoubleMatrix & CKM, int mix) {
  switch(mix) {
    case 1:       
      setYukawaMatrix(YU, CKM.transpose() * displayYukawaMatrix(YU) * CKM);
      break;
    case 2: 
      setYukawaMatrix(YD, CKM * displayYukawaMatrix(YD) * CKM.transpose()); 
      break;
     
    default:
    ostringstream ii;
    ii << "Error. EssmSusy::quarkMixing called with mix=" << mix;
    throw ii.str();
  }
}

void EssmSusy::getQuarkMixedYukawas(const QedQcd & r, const DoubleMatrix &
				    CKM, int mix, double vev) { 
  setDiagYukawas(r, vev);
  quarkMixing(CKM, mix);
}

// outputs object QedQcd & r valid at 1 GeV from SUSY data at mt, at
// present, only diagonal masses are handled. 
void EssmSusy::getMasses(QedQcd & r, double vev) const {
  double v1, v2;
  v1 = vev * cos(atan(displayTanb()));
  v2 = vev * sin(atan(displayTanb()));
  
  DoubleMatrix u1(displayYukawaMatrix(YU)), d1(displayYukawaMatrix(YD)),
    e1(displayYukawaMatrix(YE));
  r.setMass(mUp, v2 * u1(1, 1));
  r.setMass(mCharm, v2 * u1(2, 2));
  r.setMass(mTop, v2 * u1(3, 3));
  r.setMass(mDown, v1 * d1(1, 1));
  r.setMass(mStrange, v1 * d1(2, 2));
  r.setMass(mBottom, v1 * d1(3, 3));
  r.setMass(mElectron, v1 * e1(1, 1));
  r.setMass(mMuon, v1 * e1(2, 2));
  r.setMass(mTau, v1 * e1(3, 3));
}

#undef HR

// Rotates to quark mass basis, returning the mixing matrices defined as 
// yu_diag = vul yu vur^+  
// yd_diag = vdl yd vdr^+ 
// All matrices should be 3 by 3
void EssmSusy::diagQuarkBasis(DoubleMatrix & vdl, DoubleMatrix & vdr, 
			DoubleMatrix & vul, DoubleMatrix & vur) const {
  DoubleMatrix u(3, 3), v(3, 3);
  DoubleVector ydDiag(3), yuDiag(3);
  displayYukawaMatrix(YU).diagonalise(u, v, yuDiag);
  vul = u.transpose(); vur = v.transpose();

  displayYukawaMatrix(YD).diagonalise(u, v, ydDiag);
  vdl = u.transpose(); vdr = v.transpose();
}

