/*
  E6SSM version of susy.cpp
*/

#include "essmsusy.h"

#define HR "---------------------------------------------------------------\n"


const essBrevity & essBrevity::operator=(const essBrevity &s) {
  if (this == &s) return *this;
  dt = s.dt; ut = s.ut; et = s.et; 
  u2 = s.u2; d2 = s.d2; e2 = s.e2; 
  u2t = s.u2t; e2t = s.e2t; d2t = s.d2t; 
  gsq = s.gsq; g3 = s.g3; g4 = s.g4;
  uuT = s.uuT; ddT = s.ddT; eeT = s.eeT;
  d1 = s.d1; e1 = s.e1; u1 = s.u1;
  Gasq = s.Gasq; Ga3 = s.Ga3; Ga4 = s.Ga4;
  lt = s.lt; kt = s.kt;
  l2 = s.l2; k2 = s.k2;
  l2t = s.l2t; k2t = s.k2t;
  l1 = s.l1; k1 = s.k1;
  llT = s.llT; kkT = s.kkT; 
  return *this;
}

void essBrevity::calculate(const DoubleMatrix & yu, const DoubleMatrix & yd,
			  const DoubleMatrix & ye, const DoubleVector & g, 
			   const DoubleMatrix & Ga, const DoubleVector & lam, 
			   const DoubleVector & kap) {
  static DoubleVector g1(1, 2);
  static DoubleMatrix Ga1(2,2);

  g1 = g.display();
  Ga1 = Ga.display();
  u1 = yu.display(); 
  d1 = yd.display(); 
  e1 = ye.display();
  dt = d1.transpose(); ut = u1.transpose(); et = e1.transpose();
  u2 = u1 * ut; d2 = d1 * dt; e2 = e1 * et; 
  u2t = ut * u1; d2t = dt * d1; e2t = et * e1;
  uuT = u2.trace(), ddT = d2.trace(), eeT = e2.trace();
  gsq = g1 * g1; g3 = gsq * g1; g4 = g3 * g1;
  Gasq = Ga1.apply(sqr);
  for (int i = 1; i <= 2; i++)
    {
      for (int j = 1; j <= 2; j++)
	{
	  Ga3(i,j) = Gasq(i,j) * Ga1(i,j);
	  Ga4(i,j) = Gasq(i,j) * Gasq(i,j);
	}
    }

  // DH:: Even though \kappa and \lambda are vectors,
  // we treat them as matrices to make the structure of
  // the anomalous dimension matrices clearer
  l1 = DoubleMatrix(lam.display());
  k1 = DoubleMatrix(kap.display());
  lt = l1.transpose(); kt = k1.transpose();
  l2 = l1 * lt; k2 = k1 * kt;
  l2t = lt * l1; k2t = kt * k1;
  llT = l2.trace(); kkT = k2.trace();

}


/// DH:: Note that vector of gauge couplings is length 2, containing
/// only the non-Abelian couplings
EssmSusy::EssmSusy()
  : u(3, 3), d(3, 3), e(3, 3), g(2), Gmat(2,2), lambda(3), kappa(3),
    thetaE6(0.0), smuPr(0.0), tanb(0.0), hVev(0.0), sVev(0.0) 
{
  setPars(numEssmSusyPars);
  setMu(0.0);
  setLoops(2);
  setThresholds(0);
}

EssmSusy::EssmSusy(const EssmSusy &s)
  : u(s.u), d(s.d), e(s.e), g(s.g), Gmat(s.Gmat), lambda(s.lambda), kappa(s.kappa), 
    thetaE6(s.thetaE6), smuPr(s.smuPr), tanb(s.tanb), hVev(s.hVev), sVev(s.sVev) 
{ 
  setPars(numEssmSusyPars);
  setMu(s.displayMu()); 
  setLoops(s.displayLoops());
  setThresholds(s.displayThresholds());
}

EssmSusy::EssmSusy(const DoubleMatrix & u, const DoubleMatrix & d, const
		   DoubleMatrix & e, const DoubleVector & v, 
		   const DoubleMatrix & G, const DoubleVector & lam, 
		   const DoubleVector & kap, double mp, double tb, 
		   double h, double s, double th,
		   double Q, int l, int t)
  : u(u), d(d), e(e), g(v), Gmat(G), lambda(lam), kappa(kap), smuPr(mp), tanb(tb), hVev(h),
    sVev(s), thetaE6(th) 
{ 
  setPars(numEssmSusyPars);
  setMu(Q); 
  setLoops(l);
  setThresholds(t);
}

const EssmSusy & EssmSusy::operator=(const EssmSusy & s) {
  if (this == &s) return *this;
  u = s.u;
  d = s.d;
  e = s.e;
  smuPr = s.smuPr;
  tanb = s.tanb;
  hVev = s.hVev;
  sVev = s.sVev;
  g = s.g;
  Gmat = s.Gmat;
  lambda = s.lambda;
  kappa = s.kappa;
  thetaE6 = s.thetaE6;


  setMu(s.displayMu());
  setLoops(s.displayLoops());
  setThresholds(s.displayThresholds());

  return *this;
}

void EssmSusy::setSomePars(const EssmSusy & s) {
  u = s.u;
  d = s.d;
  e = s.e;
  g = s.g;
  Gmat = s.Gmat;
  lambda = s.lambda;
  kappa = s.kappa;
}

void EssmSusy::setYukawaElement(yukawa k, int i, int j, double f) { 
  switch(k) {
  case YU: u(i, j) = f; break;
  case YD: d(i, j) = f; break;
  case YE: e(i, j) = f; break;
  default: 
    ostringstream ii;
    ii << "EssmSusy::set called with illegal " << int(k) << "\n";
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
    ii << "EssmSusy::set called with illegal " << int(k) << "\n";
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
    ii << "EssmSusy::display called with illegal " << int(k) << "\n";
    throw ii.str(); break;
  }
  return 0.0;
}

const DoubleMatrix & EssmSusy::displayYukawaMatrix(yukawa k) const {
  switch(k) {
  case YU: return u; break;
  case YD: return d; break;
  case YE: return e; break;
  default: 
    ostringstream ii;    
    ii << "EssmSusy::display called with illegal " << int(k) << "\n";
    throw ii.str(); break;
  }
}

  /// DH:: Returns U(1)_Y charges for the superfields
  /// as a vector. Ordering is {Q, u, d, L, e, N, S, H1, H2, D, Dbar, H', H'bar}
DoubleVector EssmSusy::displayU1Y() const
{
  DoubleVector y(numEssmCharges);
  y(1) = QQY;
  y(2) = QuY;
  y(3) = QdY;
  y(4) = QLY;
  y(5) = QeY;
  y(6) = QNY;
  y(7) = QSY;
  y(8) = QH1Y;
  y(9) = QH2Y;
  y(10) = QXY;
  y(11) = QXbarY;
  y(12) = QHPrY;
  y(13) = QHbarPrY;

  return y;
}
/// DH:: Returns a single U(1)_Y charge
double EssmSusy::displayU1YCharge(U1Charge k) const
{
  switch(k)
    {
    case QQ: return QQY; break;
    case Qu: return QuY; break;
    case Qd: return QdY; break;
    case QL: return QLY; break;
    case Qe: return QeY; break;
    case QN: return QNY; break;
    case QS: return QSY; break;
    case QH1: return QH1Y; break;
    case QH2: return QH2Y; break;
    case QX: return QXY; break;
    case QXbar: return QXbarY; break;
    case QHPr: return QHPrY; break;
    case QHbarPr: return QHbarPrY; break;
    default: 
      ostringstream ii;    
      ii << "EssmSusy::display called with illegal " << int(k) << "\n";
      throw ii.str(); break;
    }
}

/// DH:: Returns U(1)' charges for the superfields
/// as a vector. Ordering is {Q, u, d, L, e, N, S, H1, H2, D, Dbar, H', H'bar}
DoubleVector EssmSusy::displayU1Prime() const
{
  double cE6 = cos(thetaE6);
  double sE6 = sin(thetaE6);

  DoubleVector y(numEssmCharges);
  y(1) = QQChi*cE6+QQPsi*sE6;
  y(2) = QuChi*cE6+QuPsi*sE6;
  y(3) = QdChi*cE6+QdPsi*sE6;
  y(4) = QLChi*cE6+QLPsi*sE6;
  y(5) = QeChi*cE6+QePsi*sE6;
  y(6) = QNChi*cE6+QNPsi*sE6;
  y(7) = QSChi*cE6+QSPsi*sE6;
  y(8) = QH1Chi*cE6+QH1Psi*sE6;
  y(9) = QH2Chi*cE6+QH2Psi*sE6;
  y(10) = QXChi*cE6+QXPsi*sE6;
  y(11) = QXbarChi*cE6+QXbarPsi*sE6;
  y(12) = QHPrChi*cE6+QHPrPsi*sE6;
  y(13) = QHbarPrChi*cE6+QHbarPrPsi*sE6;

  return y;
}

/// DH:: Returns a single U(1)' charge 
double EssmSusy::displayU1PrimeCharge(U1Charge k) const
{
  double cE6 = cos(thetaE6);
  double sE6 = sin(thetaE6);

  switch (k)
    {
    case QQ: return QQChi*cE6+QQPsi*sE6; break;
    case Qu: return QuChi*cE6+QuPsi*sE6; break;
    case Qd: return QdChi*cE6+QdPsi*sE6; break;
    case QL: return QLChi*cE6+QLPsi*sE6; break;
    case Qe: return QeChi*cE6+QePsi*sE6; break;
    case QN: return QNChi*cE6+QNPsi*sE6; break;
    case QS: return QSChi*cE6+QSPsi*sE6; break;
    case QH1: return QH1Chi*cE6+QH1Psi*sE6; break;
    case QH2: return QH2Chi*cE6+QH2Psi*sE6; break;
    case QX: return QXChi*cE6+QXPsi*sE6; break;
    case QXbar: return QXbarChi*cE6+QXbarPsi*sE6; break;
    case QHPr: return QHPrChi*cE6+QHPrPsi*sE6; break;
    case QHbarPr: return QHbarPrChi*cE6+QHbarPrPsi*sE6; break;
    default: 
      ostringstream ii;    
      ii << "EssmSusy::display called with illegal " << int(k) << "\n";
      throw ii.str(); break;
    }
}

/// DH:: Returns scale-dependent effective U(1)' charges for the superfields
/// as a vector. Ordering is {Q, u, d, L, e, N, S, H1, H2, D, Dbar, H', H'bar}
DoubleVector EssmSusy::displayU1PrimeEff() const
{
  DoubleVector QY = displayU1Y(); 
  DoubleVector QPr = displayU1Prime();

  double g11 = Gmat(1,2);
  double g1p = Gmat(2,2);

  if (fabs(g1p) < 1.0e-100) {
    ostringstream ii;
    ii << "WARNING: asking for EssmSusy::displayU1PrimeEff(), where gauge coupling g_1' is " <<
      g1p << endl;
    throw ii.str();
  }

  double delta;

  if (fabs(g11) < EPSTOL)
    {
      delta = 0.0;
    }
  else
    {
      delta = g11/g1p;
    }

  DoubleVector y(numEssmCharges);

  for (int i = 1; i <= numEssmCharges; i++)
    {
      y(i) = QPr.display(i)+delta*QY.display(i); 
    }

  return y;

}

/// DH:: Returns a single scale-dependent effective U(1)' charge
double EssmSusy::displayU1PrimeEffCharge(U1Charge k) const
{
  double g11 = Gmat(1,2);
  double g1p = Gmat(2,2);

  if (fabs(g1p) < 1.0e-100) {
    ostringstream ii;
    ii << "WARNING: asking for EssmSusy::displayU1PrimeEff(), where gauge coupling g_1' is " <<
      g1p << endl;
    throw ii.str();
  }

  double delta;

  if (fabs(g11) < EPSTOL)
    {
      delta = 0.0;
    }
  else
    {
      delta = g11/g1p;
    }

  double temp = displayU1PrimeCharge(k)+delta*displayU1YCharge(k);

  return temp;
}

const DoubleVector EssmSusy::display() const {
  DoubleVector y(numEssmSusyPars);
  int i, j, k=0;

  /// DH:: MSSM Yukawas (27 of them)
  for (i=1; i<=3; i++) 
    {   
      for (j=1; j<=3; j++) 
	{
	  k++;
	  y(k) = u.display(i, j);
	  y(k+9) = d.display(i, j);
	  y(k+18) = e.display(i, j);
	}
    }
  k=27;
  /// DH:: Abelian gauge couplings (4 of them), in the
  /// order G(1,1), G(1, 2), G(2, 1), G(2, 2)
  for (i = 1; i <= 2; i++)
    {
      for (j = 1; j<=2; j++)
	{
	  k++;
	  y(k) = Gmat.display(i, j);
	}
    }

  k=31;
  /// DH:: Non-Abelian gauge couplings (2 of them) \{ g_2, g_3 \}
  for (i=1; i<=2; i++) {
    k++;
    y(k) = g.display(i);
  }

  k=33;
  /// DH:: E6SSM \lambda couplings (3 of them)
  for (i = 1; i <= 3; i++)
    {
      k++;
      y(k) = lambda.display(i);
    }

  k=36;
  /// DH:: E6SSM \kappa couplings (3 of them)
  for (i = 1; i <= 3; i++)
    {
      k++;
      y(k) = kappa.display(i);
    }

  y(40) = smuPr;
  y(41) = tanb;
  y(42) = hVev;
  y(43) = sVev;
  y(44) = thetaE6;

  return y;
}

void EssmSusy::set(const DoubleVector & y) {
  int i, j, k=0;

  /// DH:: MSSM Yukawa couplings (27 of them)
  for (i=1; i<=3; i++) 
    {   
      for (j=1; j<=3; j++)
	{
	  k++;
	  u(i, j) = y.display(k);
	  d(i, j) = y.display(k+9);
	  e(i, j) = y.display(k+18);
	}
    }

  k=27;
  /// DH:: Abelian gauge couplings (4 of them), 
  /// ordering is G(1,1), G(1, 2), G(2, 1), G(2, 2)
  for (i = 1; i <= 2; i++)
    {
      for (j = 1; j<= 2; j++)
	{
	  k++;
	  Gmat(i, j) = y.display(k);
	}
    }

  k=31;
  /// DH:: Non-Abelian gauge couplings (2 of them) \{g_2, g_3 \}
  for (i=1; i<=2; i++) {
    k++;
    g(i) = y.display(k);
  }

  k=33;
  /// DH:: E6SSM \lambda couplings (3 of them)
  for (i = 1; i <= 3; i++)
    {
      k++;
      lambda(i) = y.display(k);
    }

  k=36;
  /// DH:: E6SSM \kappa couplings (3 of them)
  for (i = 1; i <= 3; i++)
    {
      k++;
      kappa(i) = y.display(k);
    }

  smuPr = y.display(40);
  tanb = y.display(41);
  hVev = y.display(42);
  sVev = y.display(43);
  thetaE6 = y.display(44);
}

double EssmSusy::displayTanb() const { return tanb; }

ostream & operator <<(ostream &left, const EssmSusy &s) {
  left << "Supersymmetric parameters at Q: " << s.displayMu() << endl;
  left << " Y^U" << s.displayYukawaMatrix(YU) << " Y^D" <<
    s.displayYukawaMatrix(YD) << " Y^E" << s.displayYukawaMatrix(YE);
  left << " lambda " << s.displayLambda() << " kappa " << s.displayKappa();
  left << "higgs VEV: " << s.displayHvev() 
       << " tan beta: " << s.displayTanb() << " singlet VEV: " << s.displaySvev() 
       << " smuPr: " << s.displaySusyMuPrime() << "\n";
  left << "g1: " << s.displayAbelianGaugeCoupling(1,1)
       << " g1Pr: "<< s.displayAbelianGaugeCoupling(2,2) 
       << " g11: " << s.displayAbelianGaugeCoupling(1,2) 
       << " g2: " <<
    s.displayNonAbelianGaugeCoupling(1) << " g3: " << 
    s.displayNonAbelianGaugeCoupling(2) << endl; 
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
  setSvev(s.displaySvev());

  setSusyMuPrime(s.displaySusyMuPrime());

  setAllNonAbelianGauge(s.displayNonAbelianGauge());
  setAllAbelianGauge(s.displayAbelianGauge());

  setLambda(s.displayLambda());
  setKappa(s.displayKappa());

  thetaE6 = s.displayThetaE6();

}

istream & operator >>(istream &left, EssmSusy &s) {
  char c[70];
  DoubleMatrix u(3, 3), d(3, 3), e(3, 3);
  DoubleVector lam(3), kap(3);
  double g1, g1p, g11, g2, g3, thE6, muP, mu, tanb, hv, sv;
  int loops, thresh;
  left >> c >> c >> c >> c >> mu;
  left >> c >> u >> c >> d >> c >> e >> c >> lam >> c >> kap;
  left >> c >> c >> hv;
  left >> c >> c >> tanb >> c >> c >> sv >> c >> muP;
  left >> c >> g1 >> c >> g1p >> c >> g11 >> c >> g2 >> c >> g3;
  left >> c >> thresh >> c >> loops;

  s.setYukawaMatrix(YU, u);
  s.setYukawaMatrix(YD, d);
  s.setYukawaMatrix(YE, e);

  s.setHvev(hv);
  s.setTanb(tanb);
  s.setSvev(sv);

  s.setAbelianGaugeCoupling(1, 1, g1);
  s.setAbelianGaugeCoupling(2, 2, g1p);
  s.setAbelianGaugeCoupling(1, 2, g11);
  s.setAbelianGaugeCoupling(2, 1, 0.0);
  s.setNonAbelianGaugeCoupling(1, g2);
  s.setNonAbelianGaugeCoupling(2, g3);

  s.setLambda(lam);
  s.setKappa(kap);

  s.setThresholds(thresh);
  s.setSusyMuPrime(muP);
  s.setMu(mu);
  s.setLoops(loops);
  return left;
}

/*
DH:: For now we are just coding the beta functions directly
while I don't have time to work out the anomalous dimensions.
Later come back and rewrite in this form DH::TODO write in MSSM
like form

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
  static DoubleMatrix gEE(3, 3), gLL(3, 3), gQQ(3, 3), gDD(3, 3), 
    gUU(3, 3);
  
  double gH1H1=0.0, gH2H2=0.0;
  static DoubleVector dg(1,3);
  
  // keep this option in order to interface with RPVSUSY  
  anomalousDimension(gEE, gLL, gQQ, gUU, gDD, dg, gH1H1, gH2H2, a);
  
  // To keep this a const function
  const DoubleMatrix &u1 = u.display(), &d1 = d.display(), &e1 = e.display();
  
  // contain derivatives of up, down quarks and leptons
  static DoubleMatrix du(3, 3), dd(3, 3), de(3, 3); 
  // mu parameter derivatives
  double dmu;
  
  // RGEs of SUSY parameters
  du = u1 * (gUU + gH2H2) + gQQ * u1;
  dd = d1 * (gDD + gH1H1) + gQQ * d1;
  de = e1 * (gEE + gH1H1) + gLL * e1;
  
  dmu = smu * (gH1H1 + gH2H2);

  // Following is from hep-ph/9308335: scalar H anomalous dimensions (as
  // opposed to the chiral superfield one - see hep-ph/0111209).
  // Additional contribution from Feynman gauge running at two-loops of tan
  // beta: we need this to link up with BPMZ: hep-ph/0112251
  double &uuT = a.uuT, &ddT = a.ddT, &eeT = a.eeT;
  DoubleVector &gsq=a.gsq;
  DoubleMatrix &u2=a.u2, &d2=a.d2, &e2=a.e2, &d2t=a.d2t;
  double t = (d2 * u2).trace();
  static const double oneLoop = 1.0 / (16.0 * sqr(PI));
  double sH1H1 = oneLoop * (3.0 * ddT + eeT);
  double sH2H2 = oneLoop * 3.0 * uuT;

  const static double twolp = 4.010149318236068e-5; // 1/(16 pi^2)^2
  if (displayLoops() > 1) {
    // I don't posess the O(g^4) terms for these RGEs in the Feynman gauge
    // and consequently have neglected. They CANCEL in the RGE for tan
    // beta, but not in the RGE of the Higgs vev. 
    sH1H1 = sH1H1 + twolp * 
      (-(3.0 * (e2 * e2).trace() + 9.0 * (d2t * d2t).trace() + 3.0 * t) + 
       (16 * gsq(3) - 0.4 * gsq(1)) * ddT + 1.2 * gsq(1) * eeT);
    sH2H2 = sH2H2 + twolp *
      (- (9.0 * (u2 * u2).trace() + 3.0 * t) +
       (16 * gsq(3) + 0.8 * gsq(1)) * uuT);
  }

  double cosb2 = sqr(cos(atan(tanb))), sinb2 = 1.0 - cosb2;
  double feynman = 1.5 * gsq(2) + 0.3 * gsq(1);
  /// One-loop RGEs in Feynman gauge
  double dt = displayTanb() * (sH1H1 - sH2H2);
  double dHvev = hVev * 
    (cosb2 * (-sH1H1 + feynman * oneLoop) + 
     sinb2 * (-sH2H2 + feynman * oneLoop)); 

  if (displayLoops() > 1) {
    /// Two-loop pieces
    dt = dt + displayTanb() * twolp * (3.0 * ddT + eeT - 3.0 * uuT) * feynman;
    dHvev = dHvev - hVev * twolp * (cosb2 * (3.0 * ddT + eeT) +
				    sinb2 * 3.0 * uuT) * feynman;
  }

  // Contains all susy derivatives:
  EssmSusy ds(du, dd, de, dg, dmu, dt, displayMu(), displayLoops(),
	       displayThresholds(), dHvev); 

  return ds;
}
*/

// DH:: temporary, third family only version. Neglects
// most U(1) mixing, though does have 1-loop running for
// off-diagonal g_{11}.
EssmSusy EssmSusy::beta(essBrevity & a) const {

  // Wave function renormalisations: convention for g**(i, j) is that i is the
  // LOWER index and j the upper in our paper hep-ph/9902251
  static DoubleMatrix gEE(3, 3), gLL(3, 3), gQQ(3, 3), gDD(3, 3), 
    gUU(3, 3), gXX(3,3), gXbarXbar(3,3);
  static DoubleVector gH1H1(1,3), gH2H2(1,3), gSS(1,3);
  
  double gHPrHPr=0.0, gHbarPrHbarPr=0.0;

  static DoubleVector dg(1,2);//< DH:: non-Abelian couplings
  static DoubleMatrix dGa(2,2); //< DH:: Abelian couplings

  anomalousDimension(gEE, gLL, gQQ, gUU, gDD, gXX, gXbarXbar, gH1H1, 
		     gH2H2, gSS, gHPrHPr, gHbarPrHbarPr, dg, dGa, a);

  // To keep this a const function
  const DoubleMatrix &u1 = u.display(), &d1 = d.display(), &e1 = e.display();
  const DoubleVector &l1 = lambda.display(), &k1 = kappa.display();

  // contain derivatives of up, down quarks and leptons
  static DoubleMatrix du(3, 3), dd(3, 3), de(3, 3); 
  // DH:: mu' parameter derivatives
  double dmuPr;

  // DH:: E6SSM Yukawa derivatives
  static DoubleVector dlam(1,3), dkap(1,3);

  du = u1 * (gUU + gH2H2(3)) + gQQ * u1;
  dd = d1 * (gDD + gH1H1(3)) + gQQ * d1;
  de = e1 * (gEE + gH1H1(3)) + gLL * e1;
  for (int i = 1; i <= 3; i++)
    {
      dlam(i) = l1(i) * (gH1H1(i) + gSS(3)) + gH2H2(i) * l1(i);
      dkap(i) = k1(i) * (gXbarXbar(i,i) + gSS(3)) + gXX(i,i) * k1(i);
    }
   
  dmuPr = smuPr * (gHPrHPr + gHbarPrHbarPr);

  // RGEs for VEVs in the Feynman gauge, from hep-ph/1310.7629, 
  // but neglecting the U(1) mixing effects that appear there, 
  // because at the moment that is a bit messy...
  const static int xi = 1; //< DH:: select Feynman gauge
  double &uuT = a.uuT, &ddT = a.ddT, &eeT = a.eeT;
  double &llT = a.llT, &kkT = a.kkT;
  DoubleVector &gsq=a.gsq;
  DoubleMatrix &Gasq = a.Gasq;
  DoubleMatrix &u2=a.u2, &d2=a.d2, &e2=a.e2, &d2t=a.d2t;
  DoubleMatrix &l2=a.l2, &k2=a.k2;
  double t = (d2 * u2).trace();
  
  static const double oneLoop = 1.0 / (16.0 * sqr(PI));

  // 1-loop betas for v_u, v_d and v_s
  // DH::TODO express in terms of proper scalar anomalous dimensions instead
  double sH1H1 = oneLoop * (sqr(l1(3)) + eeT + 3.0 * ddT - (1+xi)*(0.75 * gsq(1) + 0.15 * Gasq(1,1) 
								   + Gasq(2,2) * sqr(displayU1PrimeCharge(QH1))));
  double sH2H2 = oneLoop * (sqr(l1(3)) + 3.0 * uuT - (1+xi)*(0.75 * gsq(1) + 0.15 * Gasq(1,1)
							     + Gasq(2,2) * sqr(displayU1PrimeCharge(QH2))));
  double sSS = oneLoop * (2.0 * llT + 3.0 * kkT - (1+xi)*Gasq(2,2)*sqr(displayU1PrimeCharge(QS))); 

  const static double twolp = 4.010149318236068e-5; // 1/(16 pi^2)^2
  if (displayLoops() > 1) {
    // DH:: add in 2-loop pieces
    // DH:: First do pure gauge parts...
    sH1H1 = sH1H1 + twolp * 
      (sqr(Gasq(1,1)) * (1.4875 + 0.0225 * xi * (1-xi)) //< g_1^4
       + sqr(gsq(1)) * (5.25 - 2.1875 * xi + 0.5625 * xi * xi) //< g_2^4
       + sqr(Gasq(2,2)) * sqr(displayU1PrimeCharge(QH1)) * (9.0 * sqr(displayU1PrimeCharge(Qd))
							    + 3.0 * sqr(displayU1PrimeCharge(Qe))
							    + 8.0 * sqr(displayU1PrimeCharge(QH1))
							    + 6.0 * sqr(displayU1PrimeCharge(QH2))
							    + 2.0 * sqr(displayU1PrimeCharge(QHbarPr))
							    + 2.0 * sqr(displayU1PrimeCharge(QHPr))
							    + 6.0 * sqr(displayU1PrimeCharge(QL))
							    + 18.0 * sqr(displayU1PrimeCharge(QQ))
							    + 3.0 * sqr(displayU1PrimeCharge(QS))
							    + 9.0 * sqr(displayU1PrimeCharge(Qu))
							    + 9.0 * sqr(displayU1PrimeCharge(QXbar))
							    + 9.0 * sqr(displayU1PrimeCharge(QX))
							    + sqr(displayU1PrimeCharge(QH1)) * xi * (1-xi))//< g_1'^4
       + gsq(1) * Gasq(1,1) * (0.45 + 0.225 * xi * (1-xi)) //< g_1^2g_2^2
       + gsq(1) * Gasq(2,2) * (1.5 * sqr(displayU1PrimeCharge(QH1)) * (2+xi*(1-xi)) ) //< g_2^2g_1'^2
       + Gasq(1,1) * Gasq(2,2) * displayU1PrimeCharge(QH1) * (-1.8 * displayU1PrimeCharge(Qd)
							      -1.8 * displayU1PrimeCharge(Qe)
							      +2.4 * displayU1PrimeCharge(QH1)
							      -1.8 * displayU1PrimeCharge(QH2)
							      -0.6 * displayU1PrimeCharge(QHbarPr)
							      +0.6 * displayU1PrimeCharge(QHPr)
							      +1.8 * displayU1PrimeCharge(QL)
							      -1.8 * displayU1PrimeCharge(QQ)
							      +3.6 * displayU1PrimeCharge(Qu)
							      -1.8 * displayU1PrimeCharge(QXbar)
							      +1.8 * displayU1PrimeCharge(QX)
							      +0.3 * xi * (1-xi) * displayU1PrimeCharge(QH1))); //<g_1^2g_1'^2
    sH2H2 = sH2H2 + twolp * 
      (sqr(Gasq(1,1)) * (1.485 + 0.0225 * xi * (1-xi)) //< g_1^4
       + sqr(gsq(1)) * (5.25 - 2.1875 * xi + 0.5625 * xi * xi) //< g_2^4
       + sqr(Gasq(2,2)) * sqr(displayU1PrimeCharge(QH2)) * (9.0 * sqr(displayU1PrimeCharge(Qd))
							    + 3.0 * sqr(displayU1PrimeCharge(Qe))
							    + 6.0 * sqr(displayU1PrimeCharge(QH1))
							    + 8.0 * sqr(displayU1PrimeCharge(QH2))
							    + 2.0 * sqr(displayU1PrimeCharge(QHbarPr))
							    + 2.0 * sqr(displayU1PrimeCharge(QHPr))
							    + 6.0 * sqr(displayU1PrimeCharge(QL))
							    + 18.0 * sqr(displayU1PrimeCharge(QQ))
							    + 3.0 * sqr(displayU1PrimeCharge(QS))
							    + 9.0 * sqr(displayU1PrimeCharge(Qu))
							    + 9.0 * sqr(displayU1PrimeCharge(QXbar))
							    + 9.0 * sqr(displayU1PrimeCharge(QX))
							    + xi * (1-xi) * sqr(displayU1PrimeCharge(QH2)))//< g_1'^4
       + Gasq(1,1) * gsq(1) * (0.45 + 0.225 * xi * (1-xi)) //< g_1^2g_2^2
       + Gasq(2,2) * gsq(1) * 1.5 * sqr(displayU1PrimeCharge(QH2)) * ( 2 + xi * (1-xi)) //< g_1'^2g_2^2
       + Gasq(1,1) * Gasq(2,2) * displayU1PrimeCharge(QH2) * (1.8 * displayU1PrimeCharge(Qd)
							      + 1.8 * displayU1PrimeCharge(Qe)
							      - 1.8 * displayU1PrimeCharge(QH1)
							      + 2.4 * displayU1PrimeCharge(QH2)
							      + 0.6 * displayU1PrimeCharge(QHbarPr)
							      - 0.6 * displayU1PrimeCharge(QHPr)
							      - 1.8 * displayU1PrimeCharge(QL)
							      + 1.8 * displayU1PrimeCharge(QQ)
							      - 3.6 * displayU1PrimeCharge(Qu)
							      + 1.8 * displayU1PrimeCharge(QXbar)
							      - 1.8 * displayU1PrimeCharge(QX)
							      + 0.3 * xi * (1-xi) * displayU1PrimeCharge(QH2))); //< g_1^2g_1'^2
    sSS = sSS + twolp * (sqr(Gasq(2,2)) * sqr(displayU1PrimeCharge(QS)) * (9.0 * sqr(displayU1PrimeCharge(Qd))
									   + 3.0 * sqr(displayU1PrimeCharge(Qe))
									   + 6.0 * sqr(displayU1PrimeCharge(QH1))
									   + 6.0 * sqr(displayU1PrimeCharge(QH2))
									   + 2.0 * sqr(displayU1PrimeCharge(QHbarPr))
									   + 2.0 * sqr(displayU1PrimeCharge(QHPr))
									   + 6.0 * sqr(displayU1PrimeCharge(QL))
									   + 18.0 * sqr(displayU1PrimeCharge(QQ))
									   + 5.0 * sqr(displayU1PrimeCharge(QS)) * (1 + xi*(1-xi))
									   + 9.0 * sqr(displayU1PrimeCharge(Qu))
									   + 9.0 * sqr(displayU1PrimeCharge(QXbar))
									   + 9.0 * sqr(displayU1PrimeCharge(QX))));

    // DH:: ...then pure Yukawa parts...
    sH1H1 = sH1H1 + twolp * (-sqr(l1(3)) * (sqr(l1(3)) + 3.0 * uuT + 3.0 * kkT + 2.0 * llT)
			     - 9.0 * (d2t * d2t).trace() - 3.0 * (e2 * e2).trace() - 3.0 * (d2 * u2).trace());
    sH2H2 = sH2H2 + twolp * (-sqr(l1(3)) * (sqr(l1(3)) + 3.0 * ddT + eeT + 3.0 * kkT + 2.0 * llT)
			     - 3.0 * (d2 * u2).trace() - 9.0 * (u2 * u2).trace());
    sSS = sSS + twolp * ( -sqr(l1(3)) * (6.0 * ddT + 6.0 * uuT + 2.0 * eeT) 
			  - 6.0 * (k2 * k2).trace() - 4.0 * (l2 * l2).trace());

    // DH:: ...and finally mixed gauge and Yukawa terms
    sH1H1 = sH1H1 + twolp * (ddT * (16.0 * gsq(2) + 4.5 * gsq(1) * xi 
				    + 6.0 * Gasq(2,2) * ((xi-1)*sqr(displayU1PrimeCharge(QH1))
							 +sqr(displayU1PrimeCharge(Qd))
							 +sqr(displayU1PrimeCharge(QQ)))
				    + 0.1 * Gasq(1,1) * (9 * xi - 4) )
			     + eeT * (0.3 * Gasq(1,1) * (4 + xi) + 2.0 * Gasq(2,2) * 
				      (sqr(displayU1PrimeCharge(Qe))+sqr(displayU1PrimeCharge(QL))
				       +sqr(displayU1PrimeCharge(QH1))*(xi-1))
				      + 1.5 * gsq(1) * xi)
			     + sqr(l1(3)) * (0.3 * Gasq(1,1) * xi + 1.5 * gsq(1) * xi 
					     + 2.0 * Gasq(2,2) * (sqr(displayU1PrimeCharge(QS))
								  +sqr(displayU1PrimeCharge(QH2))
								  +sqr(displayU1PrimeCharge(QH1))*(xi-1))));
    sH2H2 = sH2H2 + twolp * (uuT * (16.0 * gsq(2) + 4.5 * xi * gsq(1) 
				    + 6.0 * Gasq(2,2) * ((xi-1)*sqr(displayU1PrimeCharge(QH2))
							 +sqr(displayU1PrimeCharge(QQ))
							 +sqr(displayU1PrimeCharge(Qu)))
				    + 0.1 * Gasq(1,1) * (9 * xi + 8 ))
			     + sqr(l1(3)) * (0.3 * Gasq(1,1) * xi + 1.5 * gsq(1) * xi
					     + 2.0 * Gasq(2,2) * (sqr(displayU1PrimeCharge(QS))
								  +sqr(displayU1PrimeCharge(QH1))
								  +sqr(displayU1PrimeCharge(QH2))*(xi-1))));
    sSS = sSS + twolp * (llT * (1.2 * Gasq(1,1) + 6.0 * gsq(1) 
				+ 4.0 * Gasq(2,2) *(sqr(displayU1PrimeCharge(QH1)) + sqr(displayU1PrimeCharge(QH2))
						    + sqr(displayU1PrimeCharge(QS))*(xi-1)))
			 + kkT * (0.8 * Gasq(1,1) + 16.0 * gsq(2) 
				  + 6.0 * Gasq(2,2) * (sqr(displayU1PrimeCharge(QX))+sqr(displayU1PrimeCharge(QXbar))
						       +sqr(displayU1PrimeCharge(QS))*(xi-1))));

  }
  
  double cosb2 = sqr(cos(atan(tanb))), sinb2 = 1.0 - cosb2;
 
  double dtb = displayTanb() * (sH1H1 - sH2H2);
  double dHvev = hVev * (-cosb2 * sH1H1 -sinb2 * sH2H2);
  double dSvev = -sVev * sSS; 
  
  
  // Contains all susy derivatives:
  // DH:: note theta_E6 doesn't vary
  EssmSusy ds(du, dd, de, dg, dGa, dlam, dkap, dmuPr, dtb, dHvev, dSvev, 0.0, displayMu(), displayLoops(),
  	       displayThresholds()); 
  
  return *this;
}

// DH:: moved to be a member function since the constants depend on the values of the
// U(1)' charges, and so on the E6 mixing angle
void EssmSusy::setEssmBetas(DoubleMatrix & babBeta, DoubleVector &cuBeta, DoubleVector
			    & cdBeta, DoubleVector & ceBeta, DoubleVector & clBeta, DoubleVector & ckBeta,
			    DoubleVector & bBeta, double bBeta11) const
{

  // DH:: traces of U(1) anomalous dimension matrices
  bBeta(1) = 18.0*sqr(QQY)+2.0*sqr(QHPrY)+2.0*sqr(QHbarPrY)+3.0*sqr(QeY)
    +3.0*sqr(QSY)+6.0*sqr(QH1Y)+6.0*sqr(QH2Y)+6.0*sqr(QLY)+9.0*sqr(QXY)+9.0*sqr(QXbarY)
    +9.0*sqr(QdY)+9.0*sqr(QuY);
  bBeta(2) = 18.0*sqr(displayU1PrimeCharge(QQ))+2.0*sqr(displayU1PrimeCharge(QHPr))+2.0*sqr(displayU1PrimeCharge(QHbarPr))
    +3.0*sqr(displayU1PrimeCharge(Qe))+3.0*sqr(displayU1PrimeCharge(QS))+6.0*sqr(displayU1PrimeCharge(QH1))
    +6.0*sqr(displayU1PrimeCharge(QH2))+6.0*sqr(displayU1PrimeCharge(QL))+9.0*sqr(displayU1PrimeCharge(QX))
    +9.0*sqr(displayU1PrimeCharge(QXbar))+9.0*sqr(displayU1PrimeCharge(Qd))+9.0*sqr(displayU1PrimeCharge(Qu));
  bBeta11 = 18.0*QQY*displayU1PrimeCharge(QQ)+2.0*QHPrY*displayU1PrimeCharge(QHPr)+2.0*QHbarPrY*displayU1PrimeCharge(QHbarPr)
    +3.0*QeY*displayU1PrimeCharge(Qe)+3.0*QSY*displayU1PrimeCharge(QS)+6.0*QH1Y*displayU1PrimeCharge(QH1)
    +6.0*QH2Y*displayU1PrimeCharge(QH2)+6.0*QLY*displayU1PrimeCharge(QL)+9.0*QXY*displayU1PrimeCharge(QX)
    +9.0*QXbar*displayU1PrimeCharge(QXbar)+9.0*QdY*displayU1PrimeCharge(Qd)+9.0*QuY*displayU1PrimeCharge(Qu);

  // 1 loop gauge beta fns
  bBeta(3) = 4.0; bBeta(4) = 0.0; 
  
  // Next come the two loop constants for gauge beta fns, 
  // neglects U(1) mixing
  babBeta(1, 1) = 234.0 / 25.0; 
  babBeta(1, 2) = 2.4*sqr(displayU1PrimeCharge(Qd))+7.2*sqr(displayU1PrimeCharge(Qe))+3.6*sqr(displayU1PrimeCharge(QH1))
    +3.6*sqr(displayU1PrimeCharge(QH2))+1.2*sqr(displayU1PrimeCharge(QHbarPr))+1.2*sqr(displayU1PrimeCharge(QHPr))
    +3.6*sqr(displayU1PrimeCharge(QL))+1.2*sqr(displayU1PrimeCharge(QQ))+9.6*sqr(displayU1PrimeCharge(Qu))
    +2.4*sqr(displayU1PrimeCharge(QXbar))+2.4*sqr(displayU1PrimeCharge(QX));
  babBeta(1, 3) = 54.0 / 5.0;
  babBeta(1, 4) = 24.0;
  babBeta(2, 1) = 2.4*sqr(displayU1PrimeCharge(Qd))+7.2*sqr(displayU1PrimeCharge(Qe))+3.6*sqr(displayU1PrimeCharge(QH1))
    +3.6*sqr(displayU1PrimeCharge(QH2))+1.2*sqr(displayU1PrimeCharge(QHbarPr))+1.2*sqr(displayU1PrimeCharge(QHPr))
    +3.6*sqr(displayU1PrimeCharge(QL))+1.2*sqr(displayU1PrimeCharge(QQ))+9.6*sqr(displayU1PrimeCharge(Qu))
    +2.4*sqr(displayU1PrimeCharge(QXbar))+2.4*sqr(displayU1PrimeCharge(QX));
  babBeta(2, 2) = 36.0*sqr(sqr(displayU1PrimeCharge(Qd)))+12.0*sqr(sqr(displayU1PrimeCharge(Qe)))
    +24.0*sqr(sqr(displayU1PrimeCharge(QH1)))+24.0*sqr(sqr(displayU1PrimeCharge(QH2)))
    +8.0*sqr(sqr(displayU1PrimeCharge(QHbarPr)))+8.0*sqr(sqr(displayU1PrimeCharge(QHPr)))
    +24.0*sqr(sqr(displayU1PrimeCharge(QL)))+72.0*sqr(sqr(displayU1PrimeCharge(QQ)))+12.0*sqr(sqr(displayU1PrimeCharge(QS)))
    +36.0*sqr(sqr(displayU1PrimeCharge(Qu)))+36.0*sqr(sqr(displayU1PrimeCharge(QXbar)))+36.0*sqr(sqr(displayU1PrimeCharge(QX)));
  babBeta(2, 3) = 18.0*sqr(displayU1PrimeCharge(QH1))+18.0*sqr(displayU1PrimeCharge(QH2))+6.0*sqr(displayU1PrimeCharge(QHbarPr))
    +6.0*sqr(displayU1PrimeCharge(QHPr))+18.0*sqr(displayU1PrimeCharge(QL))+54.0*sqr(displayU1PrimeCharge(QQ));
  babBeta(2, 4) = 48.0*sqr(displayU1PrimeCharge(Qd))+96.0*sqr(displayU1PrimeCharge(QQ))+48.0*sqr(displayU1PrimeCharge(Qu))
    +48.0*sqr(displayU1PrimeCharge(QXbar))+48.0*sqr(displayU1PrimeCharge(QX));
  babBeta(3, 1) = 18.0 / 5.0;   
  babBeta(3, 2) = 6.0*sqr(displayU1PrimeCharge(QH1))+6.0*sqr(displayU1PrimeCharge(QH2))+2.0*sqr(displayU1PrimeCharge(QHbarPr))
    +2.0*sqr(displayU1PrimeCharge(QHPr))+6.0*sqr(displayU1PrimeCharge(QL))+18.0*sqr(displayU1PrimeCharge(QQ));
  babBeta(3, 3) = 46.0;
  babBeta(3, 4) = 24.0;
  babBeta(4, 1) = 3.0;
  babBeta(4, 2) = 6.0*sqr(displayU1PrimeCharge(Qd))+6.0*sqr(displayU1PrimeCharge(Qu))+12.0*sqr(displayU1PrimeCharge(QQ))
    +6.0*sqr(displayU1PrimeCharge(QXbar))+6.0*sqr(displayU1PrimeCharge(QX));
  babBeta(4, 3) = 9.0;
  babBeta(4, 4) = 48.0;
  cuBeta(1) = 26.0 / 5.0; 
  cuBeta(2) = 12.0*(sqr(displayU1PrimeCharge(QH2))+sqr(displayU1PrimeCharge(QQ))+sqr(displayU1PrimeCharge(Qu)));
  cuBeta(3) = 6.0; 
  cuBeta(4) = 4.0;
  cdBeta(1) = 14.0 / 5.0; 
  cdBeta(2) = 12.0*(sqr(displayU1PrimeCharge(QH1))+sqr(displayU1PrimeCharge(Qd))+sqr(displayU1PrimeCharge(QQ)));
  cdBeta(3) = 6.0;
  cdBeta(4) = 4.0;
  ceBeta(1) = 18.0 / 5.0; 
  ceBeta(2) = 4.0*(sqr(displayU1PrimeCharge(Qe))+sqr(displayU1PrimeCharge(QH1))+sqr(displayU1PrimeCharge(QL)));
  ceBeta(3) = 2.0;
  ceBeta(4) = 0.0;
  clBeta(1) = 6.0 / 5.0; 
  clBeta(2) = 4.0*(sqr(displayU1PrimeCharge(QH1))+sqr(displayU1PrimeCharge(QH2))+sqr(displayU1PrimeCharge(QS)));
  clBeta(3) = 2.0;
  clBeta(4) = 0.0;
  ckBeta(1) = 4.0 / 5.0; 
  ckBeta(2) = 6.0*(sqr(displayU1PrimeCharge(QS))+sqr(displayU1PrimeCharge(QXbar))+sqr(displayU1PrimeCharge(QX)));
  ckBeta(3) = 0.0;
  ckBeta(4) = 2.0;
}

// DH::TODO
// outputs one-loop anomalous dimensions gii given matrix inputs
// Note that we use the convention (for matrices in terms of gamma's)
// gamma^Li_Lj = M_ij for LH fields and
// gamma^Rj_Ri = M_ij for RH fields (since they are really the complex
// conjugates of the RH fields): CHECKED 23/5/02
void EssmSusy::getOneLpAnom(DoubleMatrix & gEE, DoubleMatrix & gLL,
			    DoubleMatrix & gQQ, DoubleMatrix & gDD,
			    DoubleMatrix & gUU, DoubleMatrix & gXX,
			    DoubleMatrix & gXbarXbar, DoubleVector & gH1H1, 
			    DoubleVector & gH2H2, DoubleVector & gSS,
			    double & gHPrHPr, double & gHbarPrHbarPr, 
			    essBrevity & a) const {
  // For calculational brevity
  DoubleMatrix &u2=a.u2, &d2=a.d2, &e2=a.e2, &u2t=a.u2t, &d2t=a.d2t,
    &e2t=a.e2t;      
  double &uuT = a.uuT, &ddT = a.ddT, &eeT = a.eeT;
  DoubleVector &gsq=a.gsq;
  DoubleMatrix &Gasq=a.Gasq;  

  DoubleMatrix &l2=a.l2, &k2=a.k2, &l2t=a.l2t, &k2t=a.k2t;
  double &llT = a.llT, &kkT = a.kkT;

  static const double oneO16Pisq = 1.0 / (16.0 * sqr(PI));

  gQQ = oneO16Pisq * ( u2 + d2 - (1.5 * gsq(1) + 2.0 * Gasq(2,2)*sqr(displayU1PrimeCharge(QQ)) 
				  + (8.0/3.0) * gsq(2) + Gasq(1,1) / 30.0 ));
  gLL = oneO16Pisq * ( e2 - (1.5 * gsq(1) + 2.0 * Gasq(2,2) * sqr(displayU1PrimeCharge(QL))
			     + 0.3 * Gasq(1,1) )); 
  gDD = oneO16Pisq * (2.0 * d2t - (2.0 * Gasq(2,2) * sqr(displayU1PrimeCharge(Qd))
				   + (8.0/3.0) * gsq(2) + (2.0/15.0) * Gasq(1,1)));
  gUU = oneO16Pisq * (2.0 * u2t - (2.0 * Gasq(2,2) * sqr(displayU1PrimeCharge(Qu))
				   + (8.0/3.0) * gsq(2) + (8.0/15.0) * Gasq(1,1)));
  gEE = oneO16Pisq * (2.0 * e2t - (1.2 * Gasq(1,1) + 2.0 * Gasq(2,2) * sqr(displayU1PrimeCharge(Qe))));
  gXX = oneO16Pisq * ( k2 - (2.0 * Gasq(2,2) * sqr(displayU1PrimeCharge(QX)) + (8.0/3.0) * gsq(2)
			     + (2.0/15.0) * Gasq(1,1)));
  gXbarXbar = oneO16Pisq * ( k2t - ( 2.0 * Gasq(2,2) * sqr(displayU1PrimeCharge(QXbar)) 
				     + (8.0/3.0) * gsq(2) + (2.0/15.0) * Gasq(1,1)));

  for (int i = 1; i <= 2; i++)
    {
      gH1H1(i) = oneO16Pisq * ( l2t(i,i) - (1.5 * gsq(1) + 2.0 * Gasq(2,2) * sqr(displayU1PrimeCharge(QH1)) + 0.3 * Gasq(1,1)));
      gH2H2(i) = oneO16Pisq * ( l2(i,i) - (1.5 * gsq(1) + 2.0 * Gasq(2,2) * sqr(displayU1PrimeCharge(QH2)) + 0.3 * Gasq(1,1)));
      gSS(i) = oneO16Pisq * ( -2.0 * Gasq(2,2) * sqr(displayU1PrimeCharge(QS)));
    }

  gH1H1(3) = oneO16Pisq * ( l2t(3,3) + eeT + 3.0 * ddT - (0.3 * Gasq(1,1) + 1.5 * gsq(1) 
							  + 2.0 * Gasq(2,2) * sqr(displayU1PrimeCharge(QH1))));
  gH2H2(3) = oneO16Pisq * ( l2(3,3) + 3.0 * uuT - (0.3 * Gasq(1,1) + 1.5 * gsq(1) 
						   + 2.0 * Gasq(2,2) * sqr(displayU1PrimeCharge(QH2))));
  gSS(3) = oneO16Pisq * ( 2.0 * llT + 3.0 * kkT - 2.0 * Gasq(2,2) * sqr(displayU1PrimeCharge(QS)));

  gHPrHPr = oneO16Pisq * ( -2.0 * Gasq(2,2) * sqr(displayU1PrimeCharge(QHPr)) -0.3 * Gasq(1,1) - 1.5 * gsq(1));
  gHbarPrHbarPr = oneO16Pisq * ( -2.0 * Gasq(2,2) * sqr(displayU1PrimeCharge(QHbarPr)) - 0.3 * Gasq(1,1) - 1.5 * gsq(1));

}
/*
//DH::TODO
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
*/
// Outputs wave function renormalisation for SUSY parameters and gauge beta
// functions up to 2 loops.
// DH:: note modification so dg includes only non-Abelian couplings, 
// and neglects U(1) mixing
void EssmSusy::anomalousDimension(DoubleMatrix & gEE, DoubleMatrix & gLL,
				  DoubleMatrix & gQQ, DoubleMatrix & gUU,
				  DoubleMatrix & gDD, DoubleMatrix & gXX, 
				  DoubleMatrix & gXbarXbar, DoubleVector & gH1H1, 
				  DoubleVector & gH2H2, DoubleVector & gSS, 
				  double & gHPrHPr, double & gHbarPrHbarPr,
				  DoubleVector & dg, DoubleMatrix & dGa,
				  essBrevity & a)  const {
  // Constants for gauge running: order is \{ g_1, g_1', g_2, g_3 \}
  static DoubleVector bBeta(4), cuBeta(4), cdBeta(4), ceBeta(4), clBeta(4), ckBeta(4);
  static DoubleMatrix babBeta(4,4);

  // DH:: while we are neglecting U(1) mixing, we can just have this as a constant...
  static double bBeta11;

  if (bBeta(1) < 1.0e-5) // Constants not set yet
    setEssmBetas(babBeta, cuBeta, cdBeta, ceBeta, clBeta, ckBeta, bBeta, bBeta11);
  
  //  sBrevity a contains all of the shortcutted matrices etc;
  a.calculate(u.display(), d.display(), e.display(), g.display(), Gmat.display(),
	      lambda.display(), kappa.display());
  
  // For calculational brevity
  double &uuT = a.uuT, &ddT = a.ddT, &eeT = a.eeT;
  DoubleVector &gsq=a.gsq, &g3=a.g3;
  const DoubleMatrix Ga1 = Gmat.display();
  DoubleMatrix &Gasq=a.Gasq, &Ga3=a.Ga3;
  double &llT = a.llT, &kkT = a.kkT;

  static DoubleMatrix B(2,2);  

  // 1 loop contributions: 
  if (displayLoops() > 0) {

    static const double oneO16Pisq = 1.0 / (16.0 * sqr(PI)); 
    
    getOneLpAnom(gEE, gLL, gQQ, gDD, gUU, gXX, gXbarXbar, gH1H1, gH2H2, gSS, gHPrHPr, gHbarPrHbarPr, a);

    // DH:: 1-loop contributions to non-Abelian couplings
    dg(1) = oneO16Pisq * g3(1) * bBeta(3);  
    dg(2) = oneO16Pisq * g3(2) * bBeta(4);

    // DH:: 1-loop contributions to Abelian couplings
    B(1,1) = bBeta(1)*Gasq(1,1);
    B(1,2) = 2.0*Ga1(1,1)*(Ga1(2,2)*bBeta11+Ga1(1,2)*bBeta(1));
    B(2,1) = 0.0;
    B(2,2) = Gasq(2,2)*bBeta(2)+2.0*Ga1(2,2)*Ga1(1,2)*bBeta11+Gasq(1,2)*bBeta(1);

    dGa = oneO16Pisq * Ga1 * B;
  } 
  
  if (displayLoops() > 1) { 
    // DH::TODO
    //getTwoLpAnom(gEE, gLL, gQQ, gDD, gUU, gH1H1, gH2H2, a); 
			
    const static double twolp = 4.010149318236068e-5; 				      
    
    // DH:: 2-loop contributions to non-Abelian couplings
    for (int i = 1; i <= 2; i++)
      {
	dg(i) = dg(i) + g3(i) * (babBeta(i+2,1) * Gasq(1,1) + babBeta(i+2,2) * Gasq(2,2)
				 + babBeta(i+2,3) * gsq(1) + babBeta(i+2,4) * gsq(2)
				 - cuBeta(i+2) * uuT - cdBeta(i+2) *
				 ddT - ceBeta(i+2) * eeT - clBeta(i+2) * llT - ckBeta(i+2) * kkT) * twolp; 
      }

    // DH:: 2-loop contributions to Abelian couplings, neglecting
    // U(1) mixing effects (i.e. \beta_{11}^{(2)} = 0 )
    static double beta1_2lp = babBeta(1,1) * Gasq(1,1) + babBeta(1,2) * Gasq(2,2)
      + babBeta(1,3) * gsq(1) + babBeta(1,4) * gsq(2)
      - cuBeta(1) * uuT - cdBeta(1) *
      ddT - ceBeta(1) * eeT - clBeta(1) * llT - ckBeta(1) * kkT;

    static double beta1Pr_2lp = babBeta(2,1) * Gasq(1,1) + babBeta(2,2) * Gasq(2,2)
      + babBeta(2,3) * gsq(1) + babBeta(2,4) * gsq(2)
      - cuBeta(2) * uuT - cdBeta(2) *
      ddT - ceBeta(2) * eeT - clBeta(2) * llT - ckBeta(2) * kkT;    

    dGa(1,1) = dGa(1,1) + Ga3(1,1) * beta1_2lp * twolp;
    dGa(1,2) = dGa(1,2) + ( 2.0 * Gasq(1,1) * Ga1(1,2) * beta1_2lp 
			    +Gasq(2,2) * Ga1(1,2) * beta1Pr_2lp
			    +Ga3(1,2) * beta1_2lp ) * twolp;

    dGa(2,2) = dGa(2,2) + Ga1(2,2) * ( Gasq(2,2) * beta1Pr_2lp
				       +Gasq(1,2) * beta1_2lp ) * twolp;

  }

  // DH:: in this basis, other off-diagonal should always be zero, 
  // so set explicitly (might reduce numerical errors)
  dGa(2,1) = 0.0;

  // DH:: We don't have three loop expressions for the E6SSM...

}


// Outputs derivatives vector y[n] for SUSY parameters: interfaces to
// integration routines
DoubleVector EssmSusy::beta() const {
  static essBrevity a;
  
  // calculate the derivatives
  static EssmSusy ds;
  
  ds = beta(a);

  return ds.display(); // convert to a long vector
}



// r should be valid AT mt
// DH:: NB vev = \frac{v}{\sqrt{2}}
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

// DH:: NB vev = \frac{v}{\sqrt{2}}
void EssmSusy::getQuarkMixedYukawas(const QedQcd & r, const DoubleMatrix &
				    CKM, int mix, double vev) { 
  setDiagYukawas(r, vev);
  quarkMixing(CKM, mix);
}

// outputs object QedQcd & r valid at 1 GeV from SUSY data at mt, at
// present, only diagonal masses are handled.
// DH:: NB vev = \frac{v}{\sqrt{2}} 
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
// yu_diag = vul yu vur^T  
// yd_diag = vdl yd vdr^T 
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

