/*
  pmssmTuningMeasures.cpp contains the function implementations
  to calculate the fine tuning in the pMSSM.
  Author: Dylan Harries
  Date created: 5/2/2014
 */

#include "pmssmTuningMeasures.h"

// pMSSM_EWSBConditioni_TreeLevel, i = 1, 2, calculates the value of the 
// EWSB condition for m_H_i^2.
// Inputs:
//     SoftParsMssm r = pMSSM model object
//     double tb = value of tan(beta)
//     int l = # loops used. Can be either 0 (tree level EWSB conditions)
//             or 1 (leading one-loop tadpole corrections included)
// Tree level checked against SOFTSUSY with model 2403883, seems 
// correct.
double pMSSM_EWSBCondition1(SoftParsMssm r, double tb, int l)
{
  double mu = r.displaySusyMu();
  // NB SOFTSUSY convention B\mu = m_3^2
  double Bmu = r.displayM3Squared();
  double mH1Sq = r.displayMh1Squared();
  double mH2Sq = r.displayMh2Squared();
  double g1 = r.displayGaugeCoupling(1);
  double g2 = r.displayGaugeCoupling(2);
  double gbar = sqrt(g2*g2+0.6*g1*g1);
  double v = r.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;

  double f1 = 0.0;

  if (l == 0) // Tree level
    {
      f1 = mu*mu+mH1Sq-Bmu*tb+0.125*gbar*gbar*(v1*v1-v2*v2);
    }
  else if (l == 1) // One-loop tadpoles included
    {
      f1 = mu*mu+mH1Sq-Bmu*tb+0.125*gbar*gbar*(v1*v1-v2*v2)+doCalcTadpolepMSSMH1(r, tb);
    }
  else
    {
      cout << "Warning: only l = 0 or l = 1 loops currently supported." << endl;
      cout << "Skipping calculating EWSB condition." << endl;
    }

  return f1;
}

// Tree level checked against SOFTSUSY with model 2403883, seems 
// correct.
double pMSSM_EWSBCondition2(SoftParsMssm r, double tb, int l)
{
  double mu = r.displaySusyMu();
  // NB SOFTSUSY convention B\mu = m_3^2
  double Bmu = r.displayM3Squared();
  double mH1Sq = r.displayMh1Squared();
  double mH2Sq = r.displayMh2Squared();
  double g1 = r.displayGaugeCoupling(1);
  double g2 = r.displayGaugeCoupling(2);
  double gbar = sqrt(g2*g2+0.6*g1*g1);
  double v = r.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;

  double f2 = 0.0;

  if (l == 0) // Tree level
    {
      f2 = mu*mu+mH2Sq-Bmu/tb-0.125*gbar*gbar*(v1*v1-v2*v2);
    }
  else if (l == 1) // One-loop tadpoles included
    {
      f2 = mu*mu+mH2Sq-Bmu/tb-0.125*gbar*gbar*(v1*v1-v2*v2)+doCalcTadpolepMSSMH2(r, tb);
    }
  else
    {
      cout << "Warning: only l = 0 or l = 1 loops currently supported." << endl;
      cout << "Skipping calculating EWSB condition." << endl;
    }
  return f2;

}

// pMSSM_ImplementEWSBConstraints fixes the soft masses m_H_1^2 and
// m_H_2^2 to satisfy the EWSB conditions at the requested order in
// perturbation theory.
// Inputs:
//     SoftParsMssm & r = pMSSM model object
//     double tb = value of tan(beta)
//     int l = # loops used. Can be either 0 (tree level EWSB conditions)
//             or 1 (leading one-loop tadpole corrections included)
// Tree level seems correct by testing on 2403883, as does one-loop
// (assuming tadpoles correct)
void pMSSM_ImplementEWSBConstraints(SoftParsMssm & r, double tb, int l)
{
  double mu = r.displaySusyMu();
  // NB SOFTSUSY convention B\mu = m_3^2
  double Bmu = r.displayM3Squared();
  double mH1Sq, mH2Sq;
  double g1 = r.displayGaugeCoupling(1);
  double g2 = r.displayGaugeCoupling(2);
  double gbar = sqrt(g2*g2+0.6*g1*g1);
  double v = r.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;

  if (l == 0)
    {
      mH1Sq = -(mu*mu-Bmu*tb+0.125*gbar*gbar*(v1*v1-v2*v2));
      mH2Sq = -(mu*mu-Bmu/tb-0.125*gbar*gbar*(v1*v1-v2*v2));
    }
  else if (l >= 1)
    {
      if (l != 1)
	{
	  cout << "Warning: only l = 0 or l = 1 loops currently supported." << endl;
	  cout << "Using l = 1 loops." << endl;
	}
      mH1Sq = -(mu*mu-Bmu*tb+0.125*gbar*gbar*(v1*v1-v2*v2)+doCalcTadpolepMSSMH1(r, tb));
      mH2Sq = -(mu*mu-Bmu/tb-0.125*gbar*gbar*(v1*v1-v2*v2)+doCalcTadpolepMSSMH2(r, tb));
    }

  r.setMh1Squared(mH1Sq);
  r.setMh2Squared(mH2Sq);

  return;

}

// pMSSM_FineTunings calculates the fine tuning sensitivities 
// \Delta_a for the parameters a = \mu, B, m_1^2, m_2^2, m_Q^2, m_u^2, A_t at
// the requested order in perturbation theory.
// Inputs:
//     SoftParsMssm const & r = pMSSM model object
//     double tb = the value of tan(beta) to use
//     int l = # loops used. Can be either 0 (tree level EWSB conditions)
//             or 1 (leading one-loop tadpole corrections included)
//     DoubleVector & deltas = a vector (of length 7) to store the sensitivies,
//                             in the order \mu, B, m_1^2, m_2^2, m_Q^2, m_u^2, 
//                             A_t.
//     int & sing = integer flag to indicate if QR solution is singular
//     double & detF = double to store the determinant of the matrix used in 
//                     calculating the tuning
void pMSSM_FineTunings(SoftParsMssm const & r, double tb, int l, DoubleVector & deltas, 
		       int & sing, double & detF)
{
  double mu = r.displaySusyMu();
  double Bmu = r.displayM3Squared();
  // Calculate B
  double B = Bmu/mu;

  //cout << "m_A^2 = " << 2*Bmu/((2*tb/(1.0+tb*tb))) << endl;

  double mH1Sq = r.displayMh1Squared();
  double mH2Sq = r.displayMh2Squared();

  double mQlSq = r.displaySoftMassSquared(mQl, 3, 3);
  double mUrSq = r.displaySoftMassSquared(mUr, 3, 3);

  double v = r.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;

  double g1 = r.displayGaugeCoupling(1);
  double g2 = r.displayGaugeCoupling(2);
  double gbar = sqrt(g2*g2+0.6*g1*g1);

  double yt = r.displayYukawaElement(YU, 3, 3);
  double mt = yt*v2/sqrt(2.0);

  double At = r.displaySoftA(UA, 3, 3);


  // Construct relevant partial derivatives
  double df1dv1 = B*mu*v2/(v1*v1)+0.25*gbar*gbar*v1;
  double df1dv2 = -B*mu/v1-0.25*gbar*gbar*v2;
  double df2dv1 = -B*mu/v2-0.25*gbar*gbar*v1;
  double df2dv2 = B*mu*v1/(v2*v2)+0.25*gbar*gbar*v2;

  double df1dmu =2*mu-B*tb;
  double df2dmu = 2*mu-B/tb;

  double df1dmH1Sq = 1.0;
  double df2dmH1Sq = 0.0;

  double df1dmH2Sq = 0.0;
  double df2dmH2Sq = 1.0;

  double df1dB = -mu*tb;
  double df2dB = -mu/tb;

  double df1dmQlSq = 0.0;
  double df2dmQlSq = 0.0;

  double df1dmUrSq = 0.0;
  double df2dmUrSq = 0.0;

  double df1dAt = 0.0;
  double df2dAt = 0.0;

  double dv1dmu, dv2dmu, dv1dB, dv2dB, dv1dmH1Sq, dv2dmH1Sq, dv1dmH2Sq, dv2dmH2Sq;
  double dv1dmQlSq, dv2dmQlSq, dv1dmUrSq, dv2dmUrSq, dv1dAt, dv2dAt;

  // Add in loop corrections if requested
  if (l >= 1)
    {
      if (l != 1)
	{
	  cout << "Warning: only l = 0 or l = 1 loops currently supported." << endl;
	  cout << "Using l = 1 loops." << endl;
	}
      df1dv1 = df1dv1+(1.0/v1)*(doCalcpMSSMDeltaPrime11(r, tb)-doCalcTadpolepMSSMH1(r,tb));
      df1dv2 = df1dv2+doCalcpMSSMDeltaPrime12(r, tb)/v1;
      df2dv1 = df2dv1+doCalcpMSSMDeltaPrime12(r, tb)/v2;
      df2dv2 = df2dv2+(1.0/v2)*(doCalcpMSSMDeltaPrime22(r, tb)-doCalcTadpolepMSSMH2(r,tb));

      df1dmu = df1dmu+doCalcpMSSMd2DeltaVdMudv1(r, tb);
      df2dmu = df2dmu+doCalcpMSSMd2DeltaVdMudv2(r, tb);

      df1dmQlSq = df1dmQlSq+doCalcpMSSMd2DeltaVdmQlSqdv1(r, tb);
      df2dmQlSq = df2dmQlSq+doCalcpMSSMd2DeltaVdmQlSqdv2(r, tb);

      df1dmUrSq = df1dmUrSq+doCalcpMSSMd2DeltaVdmUrSqdv1(r, tb);
      df2dmUrSq = df2dmUrSq+doCalcpMSSMd2DeltaVdmUrSqdv2(r, tb);

      df1dAt = df1dAt+doCalcpMSSMd2DeltaVdAtdv1(r, tb);
      df2dAt = df2dAt+doCalcpMSSMd2DeltaVdAtdv2(r, tb);
    }

  // Construct determinant of lhs
  detF = df1dv1*df2dv2-df1dv2*df2dv1;
  if (fabs(detF) < EPSTOL) //detF = 0
    {
      cout << "Warning: accuracy bad in fine tuning calculation." << endl;
      cout << "Returned values are incorrect." << endl;
      deltas(1) = -1.0;
      deltas(2) = -1.0;
      deltas(3) = -1.0;
      deltas(4) = -1.0;
      deltas(5) = -1.0;
      deltas(6) = -1.0;
      deltas(7) = -1.0;
      sing = 1;
    }
  else
    {
      // Otherwise we are alright to calculate the fine tuning
      dv1dmu = (df1dv2*df2dmu-df2dv2*df1dmu)/detF;
      dv2dmu = (df2dv1*df1dmu-df1dv1*df2dmu)/detF;
      deltas(1) = fabs((2.0*mu/(v*v))*(v1*dv1dmu+v2*dv2dmu));

      dv1dB = (df1dv2*df2dB-df2dv2*df1dB)/detF;
      dv2dB = (df2dv1*df1dB-df1dv1*df2dB)/detF;
      deltas(2) = fabs((2.0*B/(v*v))*(v1*dv1dB+v2*dv2dB));

      dv1dmH1Sq = (df1dv2*df2dmH1Sq-df2dv2*df1dmH1Sq)/detF;
      dv2dmH1Sq = (df2dv1*df1dmH1Sq-df1dv1*df2dmH1Sq)/detF;
      deltas(3) = fabs((2.0*mH1Sq/(v*v))*(v1*dv1dmH1Sq+v2*dv2dmH1Sq));

      dv1dmH2Sq = (df1dv2*df2dmH2Sq-df2dv2*df1dmH2Sq)/detF;
      dv2dmH2Sq = (df2dv1*df1dmH2Sq-df1dv1*df2dmH2Sq)/detF;
      deltas(4) = fabs((2.0*mH2Sq/(v*v))*(v1*dv1dmH2Sq+v2*dv2dmH2Sq));

      dv1dmQlSq = (df1dv2*df2dmQlSq-df2dv2*df1dmQlSq)/detF;
      dv2dmQlSq = (df2dv1*df1dmQlSq-df1dv1*df2dmQlSq)/detF;
      deltas(5) = fabs((2.0*mQlSq/(v*v))*(v1*dv1dmQlSq+v2*dv2dmQlSq));

      dv1dmUrSq = (df1dv2*df2dmUrSq-df2dv2*df1dmUrSq)/detF;
      dv2dmUrSq = (df2dv1*df1dmUrSq-df1dv1*df2dmUrSq)/detF;
      deltas(6) = fabs((2.0*mUrSq/(v*v))*(v1*dv1dmUrSq+v2*dv2dmUrSq));

      dv1dAt = (df1dv2*df2dAt-df2dv2*df1dAt)/detF;
      dv2dAt = (df2dv1*df1dAt-df1dv1*df2dAt)/detF;
      deltas(7) = fabs((2.0*At/(v*v))*(v1*dv1dAt+v2*dv2dAt));
    }
  return;
}

// Calculates the DRbar' masses of the stops in the pMSSM.
// Inputs:
//     SoftParsMssm r = pMSSM model
//     DoubleVector & mstop = masses of stops
//     DoubleVector & mstop sq = squared masses of stops
//     double tb = value of tan(beta)
// Checked against SOFTSUSY results for model 2403883, seem correct.
void physical_pMSSM(SoftParsMssm r, DoubleVector & mstop, DoubleVector & mstopsq, double tb)
{
  bool speak = false;

  double mu = r.displaySusyMu();
  double mQlSq = r.displaySoftMassSquared(mQl, 3, 3);
  double mUrSq = r.displaySoftMassSquared(mUr, 3, 3);

  double v = r.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;

  double g1 = r.displayGaugeCoupling(1);
  double g2 = r.displayGaugeCoupling(2);
  double gbar = sqrt(g2*g2+0.6*g1*g1);

  double yt = r.displayYukawaElement(YU, 3, 3);
  double mt = yt*v2/sqrt(2.0);

  double At = r.displaySoftA(UA, 3, 3);

  double Xt = At-mu/tb;

  mstopsq(1) = 0.5*(mQlSq+mUrSq+0.125*gbar*gbar*(v1*v1-v2*v2)+2.0*mt*mt
		    -sqrt(sqr(mQlSq-mUrSq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2))
			  +4.0*mt*mt*Xt*Xt));
  mstopsq(2) = 0.5*(mQlSq+mUrSq+0.125*gbar*gbar*(v1*v1-v2*v2)+2.0*mt*mt
		    +sqrt(sqr(mQlSq-mUrSq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2))
			  +4.0*mt*mt*Xt*Xt));


  if (mstopsq(1) < 0.0 || mstopsq(2) < 0.0)
    {
      cout << "Warning: tachyonic stop masses." << endl;
      cout << "m_stop_1^2 = " << mstopsq(1) << " GeV^2." << endl;
      cout << "m_stop_2^2 = " << mstopsq(1) << " GeV^2." << endl;
    }

  mstop(1) = sqrt(fabs(mstopsq(1)));
  mstop(2) = sqrt(fabs(mstopsq(2)));

  if (speak)
    {
      cout << "m_stop_1 = " << mstop(1) << " GeV." << endl;
      cout << "m_stop_2 = " << mstop(2) << " GeV." << endl;
    }
}

// Function needed to evaluate tadpoles.
// Inputs:
//     double mSq = squared mass m^2
//     double Q = renormalisation scale Q
double a0pMSSM(double mSq, double Q)
{
  return mSq*(log(mSq/(Q*Q))-1.0);
}

// Calculates the tadpole corrections to the EWSB conditions. Note
// returns (1/v_i)\partial \Delta V/\partial v_i (i.e. opposite sign
// to equivalent E_6SSM methods).
// Inputs:
//     SoftParsMssm r = pMSSM model
//     double tb = tan(beta)
double doCalcTadpolepMSSMH1(SoftParsMssm r, double tb)
{
  bool speak = false;

  double v = r.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);

  DoubleVector mstop(2), mstopsq(2);

  physical_pMSSM(r, mstop, mstopsq, tb);

  if (mstopsq(1) < 0.0 || mstopsq(2) < 0.0)
    {
      mstopsq(1) = fabs(mstopsq(1));
      mstopsq(2) = fabs(mstopsq(2));
    }

  double q = r.displayMu();

  double a0mstop1 = a0pMSSM(mstopsq(1), q);
  double dmstop1sqdv1 = doCalcpMSSMdmstop1sqdv1(r, tb);
  double a0mstop2 = a0pMSSM(mstopsq(2), q);
  double dmstop2sqdv1 = doCalcpMSSMdmstop2sqdv1(r, tb);

  double delta1tp = (1.0/v1)*(3.0/(16.0*PI*PI))*(a0mstop1*dmstop1sqdv1+a0mstop2*dmstop2sqdv1);

  if (speak)
    {
      cout << "Delta_1^tp = " << delta1tp << endl;
    }

  return delta1tp;
}

double doCalcTadpolepMSSMH2(SoftParsMssm r, double tb)
{
  bool speak = false;

  double v = r.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;

  DoubleVector mstop(2), mstopsq(2);

  physical_pMSSM(r, mstop, mstopsq, tb);

  if (mstopsq(1) < 0.0 || mstopsq(2) < 0.0)
    {
      mstopsq(1) = fabs(mstopsq(1));
      mstopsq(2) = fabs(mstopsq(2));
    }

  double q = r.displayMu();

  double a0mstop1 = a0pMSSM(mstopsq(1), q);
  double dmstop1sqdv2 = doCalcpMSSMdmstop1sqdv2(r, tb);
  double a0mstop2 = a0pMSSM(mstopsq(2), q);
  double dmstop2sqdv2 = doCalcpMSSMdmstop2sqdv2(r, tb);

  double yt = r.displayYukawaElement(YU, 3, 3);
  double mt = yt*v2/sqrt(2.0);
  double a0mtop = a0pMSSM(mt*mt, q);
  double dmtopsqdv2 = yt*yt*v2;

  double delta2tp = (1.0/v2)*(3.0/(16.0*PI*PI))*(a0mstop1*dmstop1sqdv2+a0mstop2*dmstop2sqdv2
						 -2.0*a0mtop*dmtopsqdv2);

  if (speak)
    {
      cout << "Delta_2^tp = " << delta2tp << endl;
    }

  return delta2tp;
}

// The following are helper functions useful for constructing the 
// tuning measures at one-loop order.

// pMSSMdmstop1sqdvi calculates the derivative m_stop_1^2 wrt v_i.
// Note our convention is m_stop_1 is the lighter stop.
// Inputs:
//     SoftParsMssm r = pMSSM model
//     double tb = tan(beta)
double doCalcpMSSMdmstop1sqdv1(SoftParsMssm r, double tb)
{
  double mu = r.displaySusyMu();
  double mQlSq = r.displaySoftMassSquared(mQl, 3, 3);
  double mUrSq = r.displaySoftMassSquared(mUr, 3, 3);

  double v = r.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;

  double g1 = r.displayGaugeCoupling(1);
  double g2 = r.displayGaugeCoupling(2);
  double gbar = sqrt(g2*g2+0.6*g1*g1);

  double yt = r.displayYukawaElement(YU, 3, 3);
  double mt = yt*v2/sqrt(2.0);

  double At = r.displaySoftA(UA, 3, 3);

  double Xt = At-mu/tb;
  double MQQSq = mQlSq-mUrSq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2);
  double rt = MQQSq*MQQSq+4.0*mt*mt*Xt*Xt;
  double RQQ = MQQSq*(g2*g2-g1*g1);

  double deriv = 0.5*v1*(0.25*gbar*gbar-(0.25*RQQ-4.0*mt*mt*Xt*mu/(v1*v2))/sqrt(rt));

  return deriv;
}

double doCalcpMSSMdmstop1sqdv2(SoftParsMssm r, double tb)
{
  double mu = r.displaySusyMu();
  double mQlSq = r.displaySoftMassSquared(mQl, 3, 3);
  double mUrSq = r.displaySoftMassSquared(mUr, 3, 3);

  double v = r.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;

  double g1 = r.displayGaugeCoupling(1);
  double g2 = r.displayGaugeCoupling(2);
  double gbar = sqrt(g2*g2+0.6*g1*g1);

  double yt = r.displayYukawaElement(YU, 3, 3);
  double mt = yt*v2/sqrt(2.0);

  double At = r.displaySoftA(UA, 3, 3);

  double Xt = At-mu/tb;
  double MQQSq = mQlSq-mUrSq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2);
  double rt = MQQSq*MQQSq+4.0*mt*mt*Xt*Xt;
  double RQQ = MQQSq*(g2*g2-g1*g1);

  double deriv = 0.5*v2*(2.0*yt*yt-0.25*gbar*gbar-(2.0*yt*yt*At*Xt-0.25*RQQ)/sqrt(rt));

  return deriv;
}

// As above but for m_stop_2^2.
double doCalcpMSSMdmstop2sqdv1(SoftParsMssm r, double tb)
{
  double mu = r.displaySusyMu();
  double mQlSq = r.displaySoftMassSquared(mQl, 3, 3);
  double mUrSq = r.displaySoftMassSquared(mUr, 3, 3);

  double v = r.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;

  double g1 = r.displayGaugeCoupling(1);
  double g2 = r.displayGaugeCoupling(2);
  double gbar = sqrt(g2*g2+0.6*g1*g1);

  double yt = r.displayYukawaElement(YU, 3, 3);
  double mt = yt*v2/sqrt(2.0);

  double At = r.displaySoftA(UA, 3, 3);

  double Xt = At-mu/tb;
  double MQQSq = mQlSq-mUrSq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2);
  double rt = MQQSq*MQQSq+4.0*mt*mt*Xt*Xt;
  double RQQ = MQQSq*(g2*g2-g1*g1);

  double deriv = 0.5*v1*(0.25*gbar*gbar+(0.25*RQQ-4.0*mt*mt*Xt*mu/(v1*v2))/sqrt(rt));

  return deriv;
}

double doCalcpMSSMdmstop2sqdv2(SoftParsMssm r, double tb)
{
  double mu = r.displaySusyMu();
  double mQlSq = r.displaySoftMassSquared(mQl, 3, 3);
  double mUrSq = r.displaySoftMassSquared(mUr, 3, 3);

  double v = r.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;

  double g1 = r.displayGaugeCoupling(1);
  double g2 = r.displayGaugeCoupling(2);
  double gbar = sqrt(g2*g2+0.6*g1*g1);

  double yt = r.displayYukawaElement(YU, 3, 3);
  double mt = yt*v2/sqrt(2.0);

  double At = r.displaySoftA(UA, 3, 3);

  double Xt = At-mu/tb;
  double MQQSq = mQlSq-mUrSq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2);
  double rt = MQQSq*MQQSq+4.0*mt*mt*Xt*Xt;
  double RQQ = MQQSq*(g2*g2-g1*g1);

  double deriv = 0.5*v2*(2.0*yt*yt-0.25*gbar*gbar+(2.0*yt*yt*At*Xt-0.25*RQQ)/sqrt(rt));

  return deriv;
}

// These functions calculate the second derivatives of the one-loop
// corrections to the effective potential. Returns
// \Delta_{ij}'=\partial^2\Delta V/\partial v_i\partial v_j.
// Inputs:
//     SoftParsMssm r = pMSSM model
//     double tb = tan(beta)
double doCalcpMSSMDeltaPrime11(SoftParsMssm r, double tb)
{
  double mu = r.displaySusyMu();
  double mQlSq = r.displaySoftMassSquared(mQl, 3, 3);
  double mUrSq = r.displaySoftMassSquared(mUr, 3, 3);

  double v = r.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;

  double g1 = r.displayGaugeCoupling(1);
  double g2 = r.displayGaugeCoupling(2);
  double gbar = sqrt(g2*g2+0.6*g1*g1);

  double yt = r.displayYukawaElement(YU, 3, 3);
  double mt = yt*v2/sqrt(2.0);

  double At = r.displaySoftA(UA, 3, 3);

  double Xt = At-mu/tb;
  double MQQSq = mQlSq-mUrSq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2);
  double rt = MQQSq*MQQSq+4.0*mt*mt*Xt*Xt;
  double RQQ = MQQSq*(g2*g2-g1*g1);

  double q = r.displayMu();

  // Get stop masses
  DoubleVector mstop(2), mstopsq(2);
  physical_pMSSM(r, mstop, mstopsq, tb);

  if (mstopsq(1) < 0.0 || mstopsq(2) < 0.0)
    {
      mstopsq(1) = fabs(mstopsq(1)); // in case tachyonic
      mstopsq(2) = fabs(mstopsq(2));
    }

  double term1 = sqr(0.125*gbar*gbar*v1)+(1.0/rt)*sqr(0.125*v1*RQQ-2.0*mt*mt*Xt*mu/v2);
  term1 = term1*log(mstopsq(1)*mstopsq(2)/(q*q*q*q));

  double term2 = (gbar*gbar*v1/(32.0*sqrt(rt)))*(v1*RQQ-16.0*mt*mt*Xt*mu/v2);
  term2 = term2*log(mstopsq(2)/mstopsq(1));

  double term3 = 0.125*gbar*gbar*(a0pMSSM(mstopsq(1), q)+a0pMSSM(mstopsq(2), q));

  double term4 = (1/32.0)*((1.0/sqrt(rt))*(4.0*RQQ+sqr(g2*g2-g1*g1)*v1*v1+32.0*yt*yt*mu*mu)
			   -(1.0/(rt*sqrt(rt)))*sqr(v1*RQQ-16.0*mt*mt*Xt*mu/v2));
  term4 = term4*(a0pMSSM(mstopsq(2), q)-a0pMSSM(mstopsq(1), q));

  double delta11p = (3.0/(16.0*PI*PI))*(term1+term2+term3+term4);

  return delta11p;
}

double doCalcpMSSMDeltaPrime12(SoftParsMssm r, double tb)
{

  double mu = r.displaySusyMu();
  double mQlSq = r.displaySoftMassSquared(mQl, 3, 3);
  double mUrSq = r.displaySoftMassSquared(mUr, 3, 3);

  double v = r.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;

  double g1 = r.displayGaugeCoupling(1);
  double g2 = r.displayGaugeCoupling(2);
  double gbar = sqrt(g2*g2+0.6*g1*g1);

  double yt = r.displayYukawaElement(YU, 3, 3);
  double mt = yt*v2/sqrt(2.0);

  double At = r.displaySoftA(UA, 3, 3);

  double Xt = At-mu/tb;
  double MQQSq = mQlSq-mUrSq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2);
  double rt = MQQSq*MQQSq+4.0*mt*mt*Xt*Xt;
  double RQQ = MQQSq*(g2*g2-g1*g1);

  double q = r.displayMu();

  // Get stop masses
  DoubleVector mstop(2), mstopsq(2);
  physical_pMSSM(r, mstop, mstopsq, tb);

  if (mstopsq(1) < 0.0 || mstopsq(2) < 0.0)
    {
      mstopsq(1) = fabs(mstopsq(1)); // in case tachyonic
      mstopsq(2) = fabs(mstopsq(2));
    }

  double term1 = 0.125*gbar*gbar*(yt*yt-0.125*gbar*gbar)+(1.0/rt)*(0.125*RQQ-yt*yt*Xt*mu*tb)*(yt*yt*Xt*At-0.125*RQQ);
  term1 = term1*v1*v2*log(mstopsq(1)*mstopsq(2)/(q*q*q*q));

  double term2 = 0.125*gbar*gbar*(yt*yt*Xt*At-0.125*RQQ)+(yt*yt-0.125*gbar*gbar)*(0.125*RQQ-yt*yt*Xt*mu*tb);
  term2 = term2*(v1*v2/sqrt(rt))*log(mstopsq(2)/mstopsq(1));

  double term3 = sqr(g2*g2-g1*g1)*v1*v2/32.0+yt*yt*mu*At+(v1*v2/rt)*(0.125*RQQ-yt*yt*Xt*mu*tb)*(2.0*yt*yt*Xt*At-0.25*RQQ);
  term3 = term3*((a0pMSSM(mstopsq(2), q)-a0pMSSM(mstopsq(1), q)))/sqrt(rt);

  double delta12p = (3.0/(16.0*PI*PI))*(term1+term2-term3);

  return delta12p;
}

double doCalcpMSSMDeltaPrime22(SoftParsMssm r, double tb)
{
  double mu = r.displaySusyMu();
  double mQlSq = r.displaySoftMassSquared(mQl, 3, 3);
  double mUrSq = r.displaySoftMassSquared(mUr, 3, 3);

  double v = r.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;

  double g1 = r.displayGaugeCoupling(1);
  double g2 = r.displayGaugeCoupling(2);
  double gbar = sqrt(g2*g2+0.6*g1*g1);

  double yt = r.displayYukawaElement(YU, 3, 3);
  double mt = yt*v2/sqrt(2.0);

  double At = r.displaySoftA(UA, 3, 3);

  double Xt = At-mu/tb;
  double MQQSq = mQlSq-mUrSq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2);
  double rt = MQQSq*MQQSq+4.0*mt*mt*Xt*Xt;
  double RQQ = MQQSq*(g2*g2-g1*g1);

  double q = r.displayMu();

  // Get stop masses
  DoubleVector mstop(2), mstopsq(2);
  physical_pMSSM(r, mstop, mstopsq, tb);

  if (mstopsq(1) < 0.0 || mstopsq(2) < 0.0)
    {
      mstopsq(1) = fabs(mstopsq(1)); // in case tachyonic
      mstopsq(2) = fabs(mstopsq(2));
    }

  double term1 = sqr(yt*yt-0.125*gbar*gbar)+sqr(8.0*yt*yt*At*Xt-RQQ)/(64.0*rt);
  term1 = term1*v2*v2*log(mstopsq(1)*mstopsq(2)/(q*q*q*q));

  double term2 =(v2*v2/(4.0*sqrt(rt)))*(yt*yt-0.125*gbar*gbar)*(8.0*yt*yt*At*Xt-RQQ);
  term2 = term2*log(mstopsq(2)/mstopsq(1));

  double term3 = (yt*yt-0.125*gbar*gbar)*(a0pMSSM(mstopsq(1), q)+a0pMSSM(mstopsq(2), q));

  double term4 = (1.0/sqrt(rt))*(sqr(g2*g2-g1*g1)*v2*v2/32.0-0.125*RQQ+yt*yt*At*At-(v2*v2/(32.0*rt))*sqr(8.0*Xt*At*yt*yt-RQQ));
  term4 = term4*(a0pMSSM(mstopsq(2), q)-a0pMSSM(mstopsq(1), q));

  double delta22p = (3.0/(16.0*PI*PI))*(term1+term2+term3+term4-2.0*yt*yt*yt*yt*v2*v2*log(mt*mt/(q*q))-2.0*yt*yt*a0pMSSM(mt*mt,q));

  return delta22p;
}

// These helper functions calculate the second derivatives of 
// the one-loop corrections wrt the input parameters. They 
// return (1/v_i)\partial^2\Delta V/\partial p_j\partial v_i.
// Inputs:
//     SoftParsMssm r = pMSSM model
//     double tb = tan(beta)
double doCalcpMSSMd2DeltaVdMudv1(SoftParsMssm r, double tb)
{
  double mu = r.displaySusyMu();
  double mQlSq = r.displaySoftMassSquared(mQl, 3, 3);
  double mUrSq = r.displaySoftMassSquared(mUr, 3, 3);

  double v = r.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;

  double g1 = r.displayGaugeCoupling(1);
  double g2 = r.displayGaugeCoupling(2);
  double gbar = sqrt(g2*g2+0.6*g1*g1);

  double yt = r.displayYukawaElement(YU, 3, 3);
  double mt = yt*v2/sqrt(2.0);

  double At = r.displaySoftA(UA, 3, 3);

  double Xt = At-mu/tb;
  double MQQSq = mQlSq-mUrSq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2);
  double rt = MQQSq*MQQSq+4.0*mt*mt*Xt*Xt;
  double RQQ = MQQSq*(g2*g2-g1*g1);

  double q = r.displayMu();

  // Get stop masses
  DoubleVector mstop(2), mstopsq(2);
  physical_pMSSM(r, mstop, mstopsq, tb);

  if (mstopsq(1) < 0.0 || mstopsq(2) < 0.0)
    {
      mstopsq(1) = fabs(mstopsq(1)); // in case tachyonic
      mstopsq(2) = fabs(mstopsq(2));
    }

  // Calculate required first derivatives of stops
  double dmstop1sqdv1 =  doCalcpMSSMdmstop1sqdv1(r, tb);
  double dmstop2sqdv1 =  doCalcpMSSMdmstop2sqdv1(r, tb);
  double dmstop1sqdmu = 2.0*mt*mt*Xt/(sqrt(rt)*tb);
  double dmstop2sqdmu = -2.0*mt*mt*Xt/(sqrt(rt)*tb);

  // Calculate required second derivatives of stops
  double d2mstop1sqdmudv1 = -(v1/(2.0*sqrt(rt)))*((4.0*mt*mt*Xt/(rt*tb))*(0.25*RQQ-4.0*mt*mt*Xt*mu/(v1*v2))
						  -4.0*mt*mt*Xt/(v1*v2)+4.0*mt*mt*mu/(v2*v2));
  double d2mstop2sqdmudv1 = (v1/(2.0*sqrt(rt)))*((4.0*mt*mt*Xt/(rt*tb))*(0.25*RQQ-4.0*mt*mt*Xt*mu/(v1*v2))
						 -4.0*mt*mt*Xt/(v1*v2)+4.0*mt*mt*mu/(v2*v2));

  // Calculate total
  double d2DeltaVdmudv1 = (3.0/(16.0*PI*PI))*(1.0/v1)*(dmstop1sqdmu*dmstop1sqdv1*log(mstopsq(1)/(q*q))
					      +a0pMSSM(mstopsq(1), q)*d2mstop1sqdmudv1
					      +dmstop2sqdmu*dmstop2sqdv1*log(mstopsq(2)/(q*q))
					      +a0pMSSM(mstopsq(2), q)*d2mstop2sqdmudv1);
  return d2DeltaVdmudv1;
}

double doCalcpMSSMd2DeltaVdMudv2(SoftParsMssm r, double tb)
{
  double mu = r.displaySusyMu();
  double mQlSq = r.displaySoftMassSquared(mQl, 3, 3);
  double mUrSq = r.displaySoftMassSquared(mUr, 3, 3);

  double v = r.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;

  double g1 = r.displayGaugeCoupling(1);
  double g2 = r.displayGaugeCoupling(2);
  double gbar = sqrt(g2*g2+0.6*g1*g1);

  double yt = r.displayYukawaElement(YU, 3, 3);
  double mt = yt*v2/sqrt(2.0);

  double At = r.displaySoftA(UA, 3, 3);

  double Xt = At-mu/tb;
  double MQQSq = mQlSq-mUrSq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2);
  double rt = MQQSq*MQQSq+4.0*mt*mt*Xt*Xt;
  double RQQ = MQQSq*(g2*g2-g1*g1);

  double q = r.displayMu();

  // Get stop masses
  DoubleVector mstop(2), mstopsq(2);
  physical_pMSSM(r, mstop, mstopsq, tb);

  if (mstopsq(1) < 0.0 || mstopsq(2) < 0.0)
    {
      mstopsq(1) = fabs(mstopsq(1)); // in case tachyonic
      mstopsq(2) = fabs(mstopsq(2));
    }

  // Calculate required first derivatives of stops
  double dmstop1sqdv2 =  doCalcpMSSMdmstop1sqdv2(r, tb);
  double dmstop2sqdv2 =  doCalcpMSSMdmstop2sqdv2(r, tb);
  double dmstop1sqdmu = 2.0*mt*mt*Xt/(sqrt(rt)*tb);
  double dmstop2sqdmu = -2.0*mt*mt*Xt/(sqrt(rt)*tb);

  // Calculate required second derivatives of stops
  double d2mstop1sqdmudv2 = -(v2/(2.0*sqrt(rt)))*((4.0*mt*mt*Xt/(rt*tb))*(2.0*yt*yt*At*Xt-0.25*RQQ)-2.0*yt*yt*At/tb);
  double d2mstop2sqdmudv2 = (v2/(2.0*sqrt(rt)))*((4.0*mt*mt*Xt/(rt*tb))*(2.0*yt*yt*At*Xt-0.25*RQQ)-2.0*yt*yt*At/tb);
			

  // Calculate total
  double d2DeltaVdmudv2 = (3.0/(16.0*PI*PI))*(1.0/v2)*(dmstop1sqdmu*dmstop1sqdv2*log(mstopsq(1)/(q*q))
					      +a0pMSSM(mstopsq(1), q)*d2mstop1sqdmudv2
					      +dmstop2sqdmu*dmstop2sqdv2*log(mstopsq(2)/(q*q))
					      +a0pMSSM(mstopsq(2), q)*d2mstop2sqdmudv2);
  return d2DeltaVdmudv2;
}

double doCalcpMSSMd2DeltaVdAtdv1(SoftParsMssm r, double tb)
{
  double mu = r.displaySusyMu();
  double mQlSq = r.displaySoftMassSquared(mQl, 3, 3);
  double mUrSq = r.displaySoftMassSquared(mUr, 3, 3);

  double v = r.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;

  double g1 = r.displayGaugeCoupling(1);
  double g2 = r.displayGaugeCoupling(2);
  double gbar = sqrt(g2*g2+0.6*g1*g1);

  double yt = r.displayYukawaElement(YU, 3, 3);
  double mt = yt*v2/sqrt(2.0);

  double At = r.displaySoftA(UA, 3, 3);

  double Xt = At-mu/tb;
  double MQQSq = mQlSq-mUrSq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2);
  double rt = MQQSq*MQQSq+4.0*mt*mt*Xt*Xt;
  double RQQ = MQQSq*(g2*g2-g1*g1);

  double q = r.displayMu();

  // Get stop masses
  DoubleVector mstop(2), mstopsq(2);
  physical_pMSSM(r, mstop, mstopsq, tb);

  if (mstopsq(1) < 0.0 || mstopsq(2) < 0.0)
    {
      mstopsq(1) = fabs(mstopsq(1)); // in case tachyonic
      mstopsq(2) = fabs(mstopsq(2));
    }

  // Calculate required first derivatives of stops
  double dmstop1sqdv1 =  doCalcpMSSMdmstop1sqdv1(r, tb);
  double dmstop2sqdv1 =  doCalcpMSSMdmstop2sqdv1(r, tb);
  double dmstop1sqdAt = -2.0*mt*mt*Xt/sqrt(rt);
  double dmstop2sqdAt = 2.0*mt*mt*Xt/sqrt(rt);

  // Calculate required second derivatives of stops
  double d2mstop1sqdAtdv1 = (v1/(2.0*sqrt(rt)))*((4.0*mt*mt*Xt/rt)*(0.25*RQQ-4.0*mt*mt*Xt*mu/(v1*v2))+4.0*mt*mt*mu/(v1*v2));
  double d2mstop2sqdAtdv1 = -(v1/(2.0*sqrt(rt)))*((4.0*mt*mt*Xt/rt)*(0.25*RQQ-4.0*mt*mt*Xt*mu/(v1*v2))+4.0*mt*mt*mu/(v1*v2));
			

  // Calculate total
  double d2DeltaVdAtdv1 = (3.0/(16.0*PI*PI))*(1.0/v1)*(dmstop1sqdAt*dmstop1sqdv1*log(mstopsq(1)/(q*q))
					      +a0pMSSM(mstopsq(1), q)*d2mstop1sqdAtdv1
					      +dmstop2sqdAt*dmstop2sqdv1*log(mstopsq(2)/(q*q))
					      +a0pMSSM(mstopsq(2), q)*d2mstop2sqdAtdv1);
  return d2DeltaVdAtdv1;
}

double doCalcpMSSMd2DeltaVdAtdv2(SoftParsMssm r, double tb)
{
  double mu = r.displaySusyMu();
  double mQlSq = r.displaySoftMassSquared(mQl, 3, 3);
  double mUrSq = r.displaySoftMassSquared(mUr, 3, 3);

  double v = r.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;

  double g1 = r.displayGaugeCoupling(1);
  double g2 = r.displayGaugeCoupling(2);
  double gbar = sqrt(g2*g2+0.6*g1*g1);

  double yt = r.displayYukawaElement(YU, 3, 3);
  double mt = yt*v2/sqrt(2.0);

  double At = r.displaySoftA(UA, 3, 3);

  double Xt = At-mu/tb;
  double MQQSq = mQlSq-mUrSq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2);
  double rt = MQQSq*MQQSq+4.0*mt*mt*Xt*Xt;
  double RQQ = MQQSq*(g2*g2-g1*g1);

  double q = r.displayMu();

  // Get stop masses
  DoubleVector mstop(2), mstopsq(2);
  physical_pMSSM(r, mstop, mstopsq, tb);

  if (mstopsq(1) < 0.0 || mstopsq(2) < 0.0)
    {
      mstopsq(1) = fabs(mstopsq(1)); // in case tachyonic
      mstopsq(2) = fabs(mstopsq(2));
    }

  // Calculate required first derivatives of stops
  double dmstop1sqdv2 =  doCalcpMSSMdmstop1sqdv2(r, tb);
  double dmstop2sqdv2 =  doCalcpMSSMdmstop2sqdv2(r, tb);
  double dmstop1sqdAt = -2.0*mt*mt*Xt/sqrt(rt);
  double dmstop2sqdAt = 2.0*mt*mt*Xt/sqrt(rt);

  // Calculate required second derivatives of stops
  double d2mstop1sqdAtdv2 = (v2/(2.0*sqrt(rt)))*((4.0*mt*mt*Xt/rt)*(2.0*yt*yt*At*Xt-0.25*RQQ)-2.0*yt*yt*(Xt+At));
  double d2mstop2sqdAtdv2 = -(v2/(2.0*sqrt(rt)))*((4.0*mt*mt*Xt/rt)*(2.0*yt*yt*At*Xt-0.25*RQQ)-2.0*yt*yt*(Xt+At));
			

  // Calculate total
  double d2DeltaVdAtdv2 = (3.0/(16.0*PI*PI))*(1.0/v2)*(dmstop1sqdAt*dmstop1sqdv2*log(mstopsq(1)/(q*q))
					      +a0pMSSM(mstopsq(1), q)*d2mstop1sqdAtdv2
					      +dmstop2sqdAt*dmstop2sqdv2*log(mstopsq(2)/(q*q))
					      +a0pMSSM(mstopsq(2), q)*d2mstop2sqdAtdv2);
  return d2DeltaVdAtdv2;
}

double doCalcpMSSMd2DeltaVdmQlSqdv1(SoftParsMssm r, double tb)
{
  double mu = r.displaySusyMu();
  double mQlSq = r.displaySoftMassSquared(mQl, 3, 3);
  double mUrSq = r.displaySoftMassSquared(mUr, 3, 3);

  double v = r.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;

  double g1 = r.displayGaugeCoupling(1);
  double g2 = r.displayGaugeCoupling(2);
  double gbar = sqrt(g2*g2+0.6*g1*g1);

  double yt = r.displayYukawaElement(YU, 3, 3);
  double mt = yt*v2/sqrt(2.0);

  double At = r.displaySoftA(UA, 3, 3);

  double Xt = At-mu/tb;
  double MQQSq = mQlSq-mUrSq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2);
  double rt = MQQSq*MQQSq+4.0*mt*mt*Xt*Xt;
  double RQQ = MQQSq*(g2*g2-g1*g1);

  double q = r.displayMu();

  // Get stop masses
  DoubleVector mstop(2), mstopsq(2);
  physical_pMSSM(r, mstop, mstopsq, tb);

  if (mstopsq(1) < 0.0 || mstopsq(2) < 0.0)
    {
      mstopsq(1) = fabs(mstopsq(1)); // in case tachyonic
      mstopsq(2) = fabs(mstopsq(2));
    }

  // Calculate required first derivatives of stops
  double dmstop1sqdv1 =  doCalcpMSSMdmstop1sqdv1(r, tb);
  double dmstop2sqdv1 =  doCalcpMSSMdmstop2sqdv1(r, tb);
  double dmstop1sqdmQlSq = 0.5*(1.0-MQQSq/sqrt(rt));
  double dmstop2sqdmQlSq = 0.5*(1.0+MQQSq/sqrt(rt));

  // Calculate required second derivatives of stops
  double d2mstop1sqdmQlSqdv1 = (v1/(2.0*sqrt(rt)))*((MQQSq/rt)*(0.25*RQQ-4.0*mt*mt*Xt*mu/(v1*v2))-0.25*(g2*g2-g1*g1));
  double d2mstop2sqdmQlSqdv1 = -(v1/(2.0*sqrt(rt)))*((MQQSq/rt)*(0.25*RQQ-4.0*mt*mt*Xt*mu/(v1*v2))-0.25*(g2*g2-g1*g1));
			

  // Calculate total
  double d2DeltaVdmQlSqdv1 = (3.0/(16.0*PI*PI))*(1.0/v1)*(dmstop1sqdmQlSq*dmstop1sqdv1*log(mstopsq(1)/(q*q))
					      +a0pMSSM(mstopsq(1), q)*d2mstop1sqdmQlSqdv1
					      +dmstop2sqdmQlSq*dmstop2sqdv1*log(mstopsq(2)/(q*q))
					      +a0pMSSM(mstopsq(2), q)*d2mstop2sqdmQlSqdv1);
  return d2DeltaVdmQlSqdv1;
}

double doCalcpMSSMd2DeltaVdmQlSqdv2(SoftParsMssm r, double tb)
{
  double mu = r.displaySusyMu();
  double mQlSq = r.displaySoftMassSquared(mQl, 3, 3);
  double mUrSq = r.displaySoftMassSquared(mUr, 3, 3);

  double v = r.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;

  double g1 = r.displayGaugeCoupling(1);
  double g2 = r.displayGaugeCoupling(2);
  double gbar = sqrt(g2*g2+0.6*g1*g1);

  double yt = r.displayYukawaElement(YU, 3, 3);
  double mt = yt*v2/sqrt(2.0);

  double At = r.displaySoftA(UA, 3, 3);

  double Xt = At-mu/tb;
  double MQQSq = mQlSq-mUrSq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2);
  double rt = MQQSq*MQQSq+4.0*mt*mt*Xt*Xt;
  double RQQ = MQQSq*(g2*g2-g1*g1);

  double q = r.displayMu();

  // Get stop masses
  DoubleVector mstop(2), mstopsq(2);
  physical_pMSSM(r, mstop, mstopsq, tb);

  if (mstopsq(1) < 0.0 || mstopsq(2) < 0.0)
    {
      mstopsq(1) = fabs(mstopsq(1)); // in case tachyonic
      mstopsq(2) = fabs(mstopsq(2));
    } 
  // Calculate required first derivatives of stops
  double dmstop1sqdv2 =  doCalcpMSSMdmstop1sqdv2(r, tb);
  double dmstop2sqdv2 =  doCalcpMSSMdmstop2sqdv2(r, tb);
  double dmstop1sqdmQlSq = 0.5*(1.0-MQQSq/sqrt(rt));
  double dmstop2sqdmQlSq = 0.5*(1.0+MQQSq/sqrt(rt));

  // Calculate required second derivatives of stops
  double d2mstop1sqdmQlSqdv2 = (v2/(2.0*sqrt(rt)))*((MQQSq/rt)*(2.0*yt*yt*At*Xt-0.25*RQQ)+0.25*(g2*g2-g1*g1));
  double d2mstop2sqdmQlSqdv2 = -(v2/(2.0*sqrt(rt)))*((MQQSq/rt)*(2.0*yt*yt*At*Xt-0.25*RQQ)+0.25*(g2*g2-g1*g1));
			

  // Calculate total
  double d2DeltaVdmQlSqdv2 = (3.0/(16.0*PI*PI))*(1.0/v2)*(dmstop1sqdmQlSq*dmstop1sqdv2*log(mstopsq(1)/(q*q))
					      +a0pMSSM(mstopsq(1), q)*d2mstop1sqdmQlSqdv2
					      +dmstop2sqdmQlSq*dmstop2sqdv2*log(mstopsq(2)/(q*q))
					      +a0pMSSM(mstopsq(2), q)*d2mstop2sqdmQlSqdv2);
  return d2DeltaVdmQlSqdv2;
}

double doCalcpMSSMd2DeltaVdmUrSqdv1(SoftParsMssm r, double tb)
{
  double mu = r.displaySusyMu();
  double mQlSq = r.displaySoftMassSquared(mQl, 3, 3);
  double mUrSq = r.displaySoftMassSquared(mUr, 3, 3);

  double v = r.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;

  double g1 = r.displayGaugeCoupling(1);
  double g2 = r.displayGaugeCoupling(2);
  double gbar = sqrt(g2*g2+0.6*g1*g1);

  double yt = r.displayYukawaElement(YU, 3, 3);
  double mt = yt*v2/sqrt(2.0);

  double At = r.displaySoftA(UA, 3, 3);

  double Xt = At-mu/tb;
  double MQQSq = mQlSq-mUrSq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2);
  double rt = MQQSq*MQQSq+4.0*mt*mt*Xt*Xt;
  double RQQ = MQQSq*(g2*g2-g1*g1);

  double q = r.displayMu();

  // Get stop masses
  DoubleVector mstop(2), mstopsq(2);
  physical_pMSSM(r, mstop, mstopsq, tb);

  if (mstopsq(1) < 0.0 || mstopsq(2) < 0.0)
    {
      mstopsq(1) = fabs(mstopsq(1)); // in case tachyonic
      mstopsq(2) = fabs(mstopsq(2));
    } 
  // Calculate required first derivatives of stops
  double dmstop1sqdv1 =  doCalcpMSSMdmstop1sqdv1(r, tb);
  double dmstop2sqdv1 =  doCalcpMSSMdmstop2sqdv1(r, tb);
  double dmstop1sqdmUrSq = 0.5*(1.0+MQQSq/sqrt(rt));
  double dmstop2sqdmUrSq = 0.5*(1.0-MQQSq/sqrt(rt));

  // Calculate required second derivatives of stops
  double d2mstop1sqdmUrSqdv1 = -(v1/(2.0*sqrt(rt)))*((MQQSq/rt)*(0.25*RQQ-4.0*mt*mt*Xt*mu/(v1*v2))-0.25*(g2*g2-g1*g1));
  double d2mstop2sqdmUrSqdv1 = (v1/(2.0*sqrt(rt)))*((MQQSq/rt)*(0.25*RQQ-4.0*mt*mt*Xt*mu/(v1*v2))-0.25*(g2*g2-g1*g1));
			

  // Calculate total
  double d2DeltaVdmUrSqdv1 = (3.0/(16.0*PI*PI))*(1.0/v1)*(dmstop1sqdmUrSq*dmstop1sqdv1*log(mstopsq(1)/(q*q))
					      +a0pMSSM(mstopsq(1), q)*d2mstop1sqdmUrSqdv1
					      +dmstop2sqdmUrSq*dmstop2sqdv1*log(mstopsq(2)/(q*q))
					      +a0pMSSM(mstopsq(2), q)*d2mstop2sqdmUrSqdv1);
  return d2DeltaVdmUrSqdv1;
}

double doCalcpMSSMd2DeltaVdmUrSqdv2(SoftParsMssm r, double tb)
{
  double mu = r.displaySusyMu();
  double mQlSq = r.displaySoftMassSquared(mQl, 3, 3);
  double mUrSq = r.displaySoftMassSquared(mUr, 3, 3);

  double v = r.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;

  double g1 = r.displayGaugeCoupling(1);
  double g2 = r.displayGaugeCoupling(2);
  double gbar = sqrt(g2*g2+0.6*g1*g1);

  double yt = r.displayYukawaElement(YU, 3, 3);
  double mt = yt*v2/sqrt(2.0);

  double At = r.displaySoftA(UA, 3, 3);

  double Xt = At-mu/tb;
  double MQQSq = mQlSq-mUrSq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2);
  double rt = MQQSq*MQQSq+4.0*mt*mt*Xt*Xt;
  double RQQ = MQQSq*(g2*g2-g1*g1);

  double q = r.displayMu();

  // Get stop masses
  DoubleVector mstop(2), mstopsq(2);
  physical_pMSSM(r, mstop, mstopsq, tb);

  if (mstopsq(1) < 0.0 || mstopsq(2) < 0.0)
    {
      mstopsq(1) = fabs(mstopsq(1)); // in case tachyonic
      mstopsq(2) = fabs(mstopsq(2));
    } 
  // Calculate required first derivatives of stops
  double dmstop1sqdv2 =  doCalcpMSSMdmstop1sqdv2(r, tb);
  double dmstop2sqdv2 =  doCalcpMSSMdmstop2sqdv2(r, tb);
  double dmstop1sqdmUrSq = 0.5*(1.0+MQQSq/sqrt(rt));
  double dmstop2sqdmUrSq = 0.5*(1.0-MQQSq/sqrt(rt));

  // Calculate required second derivatives of stops
  double d2mstop1sqdmUrSqdv2 = -(v2/(2.0*sqrt(rt)))*((MQQSq/rt)*(2.0*yt*yt*At*Xt-0.25*RQQ)+0.25*(g2*g2-g1*g1));
  double d2mstop2sqdmUrSqdv2 = (v2/(2.0*sqrt(rt)))*((MQQSq/rt)*(2.0*yt*yt*At*Xt-0.25*RQQ)+0.25*(g2*g2-g1*g1));
			

  // Calculate total
  double d2DeltaVdmUrSqdv2 = (3.0/(16.0*PI*PI))*(1.0/v2)*(dmstop1sqdmUrSq*dmstop1sqdv2*log(mstopsq(1)/(q*q))
					      +a0pMSSM(mstopsq(1), q)*d2mstop1sqdmUrSqdv2
					      +dmstop2sqdmUrSq*dmstop2sqdv2*log(mstopsq(2)/(q*q))
					      +a0pMSSM(mstopsq(2), q)*d2mstop2sqdmUrSqdv2);
  return d2DeltaVdmUrSqdv2;
}

/*
  --------------------------------------------------------------------------------------
  Following are methods to calculate the fine tuning numerically. They are used only
  for testing, as the analytic methods when properly checked should be reliable and
  fast enough.
  --------------------------------------------------------------------------------------
 */

// pMSSM_NumericalTadpoles calculates the tadpole corrections
// numerically to verify the analytic expressions. Uses a 5-pt
// stencil to approximate the derivatives.
// Inputs:
//     SoftParsMssm r = pMSSM model
//     double tb = tan(beta)
//     double epsilon = value of epsilon used to estimate partial derivatives
DoubleVector pMSSM_NumericalTadpoles(SoftParsMssm r, double tb, double epsilon)
{
  double q = r.displayMu();

  double yt = r.displayYukawaElement(YU, 3, 3);

  DoubleVector mstop(2), mstopsq(2);

  double v = r.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;

  double mt = yt*v2/sqrt(2.0);

  DoubleVector ti(2);

  double vnew, tbnew;

  double deltaV2stpFor, deltaV1stpFor, deltaV1stpBck, deltaV2stpBck;

  // Get 2 steps forward v_1 and value of \Delta V
  vnew = sqrt(v1*v1+v2*v2+4.0*epsilon*epsilon+4.0*v1*epsilon);
  tbnew = v2/(v1+2.0*epsilon);
  r.setHvev(vnew);
  r.setTanb(tbnew);
 
  physical_pMSSM(r, mstop, mstopsq, tbnew);

  deltaV2stpFor = (3.0/(32.0*PI*PI))*(sqr(mstopsq(1))*(log(mstopsq(1)/(q*q))-1.5)
				      +sqr(mstopsq(2))*(log(mstopsq(2)/(q*q))-1.5)
				      -2.0*sqr(mt*mt)*(log(mt*mt/(q*q))-1.5));

  // Get 1 step forward v_1
  vnew = sqrt(v1*v1+v2*v2+epsilon*epsilon+2.0*v1*epsilon);
  tbnew = v2/(v1+epsilon);
  r.setHvev(vnew);
  r.setTanb(tbnew);
 
  physical_pMSSM(r, mstop, mstopsq, tbnew);

  deltaV1stpFor = (3.0/(32.0*PI*PI))*(sqr(mstopsq(1))*(log(mstopsq(1)/(q*q))-1.5)
				      +sqr(mstopsq(2))*(log(mstopsq(2)/(q*q))-1.5)
				      -2.0*sqr(mt*mt)*(log(mt*mt/(q*q))-1.5));

  // Get 1 step backway v_1
  vnew = sqrt(v1*v1+v2*v2+epsilon*epsilon-2.0*v1*epsilon);
  tbnew = v2/(v1-epsilon);
  r.setHvev(vnew);
  r.setTanb(tbnew);
 
  physical_pMSSM(r, mstop, mstopsq, tbnew);

  deltaV1stpBck = (3.0/(32.0*PI*PI))*(sqr(mstopsq(1))*(log(mstopsq(1)/(q*q))-1.5)
				      +sqr(mstopsq(2))*(log(mstopsq(2)/(q*q))-1.5)
				      -2.0*sqr(mt*mt)*(log(mt*mt/(q*q))-1.5));

  // Get 2 steps backward v_1
  vnew = sqrt(v1*v1+v2*v2+4.0*epsilon*epsilon-4.0*v1*epsilon);
  tbnew = v2/(v1-2.0*epsilon);
  r.setHvev(vnew);
  r.setTanb(tbnew);
 
  physical_pMSSM(r, mstop, mstopsq, tbnew);

  deltaV2stpBck = (3.0/(32.0*PI*PI))*(sqr(mstopsq(1))*(log(mstopsq(1)/(q*q))-1.5)
				      +sqr(mstopsq(2))*(log(mstopsq(2)/(q*q))-1.5)
				      -2.0*sqr(mt*mt)*(log(mt*mt/(q*q))-1.5));

  // Estimate partial derivative
  ti(1)=(1.0/(12.0*epsilon))*(-1.0*deltaV2stpFor+8.0*deltaV1stpFor-8.0*deltaV1stpBck+deltaV2stpBck);
  ti(1)=ti(1)/v1;

  // Repeat for v_2
  vnew = sqrt(v1*v1+v2*v2+4.0*epsilon*epsilon+4.0*v2*epsilon);
  tbnew = (v2+2.0*epsilon)/v1;
  r.setHvev(vnew);
  r.setTanb(tbnew);

  mt = yt*(v2+2.0*epsilon)/sqrt(2.0); 

  physical_pMSSM(r, mstop, mstopsq, tbnew);

  deltaV2stpFor = (3.0/(32.0*PI*PI))*(sqr(mstopsq(1))*(log(mstopsq(1)/(q*q))-1.5)
				      +sqr(mstopsq(2))*(log(mstopsq(2)/(q*q))-1.5)
				      -2.0*sqr(mt*mt)*(log(mt*mt/(q*q))-1.5));

  vnew = sqrt(v1*v1+v2*v2+epsilon*epsilon+2.0*v2*epsilon);
  tbnew = (v2+epsilon)/v1;
  r.setHvev(vnew);
  r.setTanb(tbnew);

  mt = yt*(v2+epsilon)/sqrt(2.0); 

  physical_pMSSM(r, mstop, mstopsq, tbnew);

  deltaV1stpFor = (3.0/(32.0*PI*PI))*(sqr(mstopsq(1))*(log(mstopsq(1)/(q*q))-1.5)
				      +sqr(mstopsq(2))*(log(mstopsq(2)/(q*q))-1.5)
				      -2.0*sqr(mt*mt)*(log(mt*mt/(q*q))-1.5));

  vnew = sqrt(v1*v1+v2*v2+epsilon*epsilon-2.0*v2*epsilon);
  tbnew = (v2-epsilon)/v1;
  r.setHvev(vnew);
  r.setTanb(tbnew);

  mt = yt*(v2-epsilon)/sqrt(2.0); 

  physical_pMSSM(r, mstop, mstopsq, tbnew);

  deltaV1stpBck = (3.0/(32.0*PI*PI))*(sqr(mstopsq(1))*(log(mstopsq(1)/(q*q))-1.5)
				      +sqr(mstopsq(2))*(log(mstopsq(2)/(q*q))-1.5)
				      -2.0*sqr(mt*mt)*(log(mt*mt/(q*q))-1.5));

  vnew = sqrt(v1*v1+v2*v2+4.0*epsilon*epsilon-4.0*v2*epsilon);
  tbnew = (v2-2.0*epsilon)/v1;
  r.setHvev(vnew);
  r.setTanb(tbnew);

  mt = yt*(v2-2.0*epsilon)/sqrt(2.0); 

  physical_pMSSM(r, mstop, mstopsq, tbnew);

  deltaV2stpBck = (3.0/(32.0*PI*PI))*(sqr(mstopsq(1))*(log(mstopsq(1)/(q*q))-1.5)
				      +sqr(mstopsq(2))*(log(mstopsq(2)/(q*q))-1.5)
				      -2.0*sqr(mt*mt)*(log(mt*mt/(q*q))-1.5));

  ti(2)=(1.0/(12.0*epsilon))*(-1.0*deltaV2stpFor+8.0*deltaV1stpFor-8.0*deltaV1stpBck+deltaV2stpBck);
  ti(2)=ti(2)/v2;

  return ti;

}

// pMSSM_NumericalDeltaPrime calculates the three second derivatives of the 
// loop corrections to the effective potential, using the analytic
// expressions for the tadpole corrections (i.e. the first derivatives).
// These analytic expressions have been checked against the numerical
// methods for calculating the first derivatives. Outputs the 4 second
// derivatives as a matrix with elements M_{ij}=\Delta_{ij}'.
// Inputs:
//     SoftParsMssm r = pMSSM model
//     double tb = tan(beta)
//     double epsilon = value of epsilon used in estimating the derivative
DoubleMatrix pMSSM_NumericalDeltaPrime(SoftParsMssm r, double tb, double epsilon)
{
  double q = r.displayMu();
  double v = r.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;
  double yt = r.displayYukawaElement(YU, 3, 3);
  double mt = yt*v2/sqrt(2.0);

  double vnew, tbnew;

  DoubleMatrix deltaij(2,2);

  double ti2stpFor, ti1stpFor, ti1stpBck, ti2stpBck;

  // Estimate \Delta_{11}'
  vnew = sqrt(v1*v1+v2*v2+4.0*epsilon*epsilon+4.0*epsilon*v1);
  tbnew = v2/(v1+2.0*epsilon);
  r.setHvev(vnew);
  r.setTanb(tbnew);

  ti2stpFor = (v1+2.0*epsilon)*doCalcTadpolepMSSMH1(r, tbnew);

  vnew = sqrt(v1*v1+v2*v2+epsilon*epsilon+2.0*epsilon*v1);
  tbnew = v2/(v1+epsilon);
  r.setHvev(vnew);
  r.setTanb(tbnew);

  ti1stpFor = (v1+epsilon)*doCalcTadpolepMSSMH1(r, tbnew);

  vnew = sqrt(v1*v1+v2*v2+epsilon*epsilon-2.0*epsilon*v1);
  tbnew = v2/(v1-epsilon);
  r.setHvev(vnew);
  r.setTanb(tbnew);

  ti1stpBck = (v1-epsilon)*doCalcTadpolepMSSMH1(r, tbnew);

  vnew = sqrt(v1*v1+v2*v2+4.0*epsilon*epsilon-4.0*epsilon*v1);
  tbnew = v2/(v1-2.0*epsilon);
  r.setHvev(vnew);
  r.setTanb(tbnew);

  ti2stpBck = (v1-2.0*epsilon)*doCalcTadpolepMSSMH1(r, tbnew);

  deltaij(1,1) = (1.0/(12.0*epsilon))*(-1.0*ti2stpFor+8.0*ti1stpFor-8.0*ti1stpBck+ti2stpBck);

  // Approximate \Delta_{12}
  vnew = sqrt(v1*v1+v2*v2+4.0*epsilon*epsilon+4.0*epsilon*v1);
  tbnew = v2/(v1+2.0*epsilon);
  r.setHvev(vnew);
  r.setTanb(tbnew);

  ti2stpFor = v2*doCalcTadpolepMSSMH2(r, tbnew);

  vnew = sqrt(v1*v1+v2*v2+epsilon*epsilon+2.0*epsilon*v1);
  tbnew = v2/(v1+epsilon);
  r.setHvev(vnew);
  r.setTanb(tbnew);

  ti1stpFor = v2*doCalcTadpolepMSSMH2(r, tbnew);

  vnew = sqrt(v1*v1+v2*v2+epsilon*epsilon-2.0*epsilon*v1);
  tbnew = v2/(v1-epsilon);
  r.setHvev(vnew);
  r.setTanb(tbnew);

  ti1stpBck = v2*doCalcTadpolepMSSMH2(r, tbnew);

  vnew = sqrt(v1*v1+v2*v2+4.0*epsilon*epsilon-4.0*epsilon*v1);
  tbnew = v2/(v1-2.0*epsilon);
  r.setHvev(vnew);
  r.setTanb(tbnew);

  ti2stpBck = v2*doCalcTadpolepMSSMH2(r, tbnew);

  deltaij(1,2) = (1.0/(12.0*epsilon))*(-1.0*ti2stpFor+8.0*ti1stpFor-8.0*ti1stpBck+ti2stpBck);

  // Approximate \Delta_{21}
  vnew = sqrt(v1*v1+v2*v2+4.0*epsilon*epsilon+4.0*epsilon*v2);
  tbnew = (v2+2.0*epsilon)/v1;
  r.setHvev(vnew);
  r.setTanb(tbnew);

  ti2stpFor = v1*doCalcTadpolepMSSMH1(r, tbnew);

  vnew = sqrt(v1*v1+v2*v2+epsilon*epsilon+2.0*epsilon*v2);
  tbnew = (v2+epsilon)/v1;
  r.setHvev(vnew);
  r.setTanb(tbnew);

  ti1stpFor = v1*doCalcTadpolepMSSMH1(r, tbnew);

  vnew = sqrt(v1*v1+v2*v2+epsilon*epsilon-2.0*epsilon*v2);
  tbnew = (v2-epsilon)/v1;
  r.setHvev(vnew);
  r.setTanb(tbnew);

  ti1stpBck = v1*doCalcTadpolepMSSMH1(r, tbnew);

  vnew = sqrt(v1*v1+v2*v2+4.0*epsilon*epsilon-4.0*epsilon*v2);
  tbnew = (v2-2.0*epsilon)/v1;
  r.setHvev(vnew);
  r.setTanb(tbnew);

  ti2stpBck = v1*doCalcTadpolepMSSMH1(r, tbnew);

  deltaij(2,1) = (1.0/(12.0*epsilon))*(-1.0*ti2stpFor+8.0*ti1stpFor-8.0*ti1stpBck+ti2stpBck);

  // Approximate \Delta_{22}
  vnew = sqrt(v1*v1+v2*v2+4.0*epsilon*epsilon+4.0*epsilon*v2);
  tbnew = (v2+2.0*epsilon)/v1;
  r.setHvev(vnew);
  r.setTanb(tbnew);

  ti2stpFor = (v2+2.0*epsilon)*doCalcTadpolepMSSMH2(r, tbnew);

  vnew = sqrt(v1*v1+v2*v2+epsilon*epsilon+2.0*epsilon*v2);
  tbnew = (v2+epsilon)/v1;
  r.setHvev(vnew);
  r.setTanb(tbnew);

  ti1stpFor = (v2+epsilon)*doCalcTadpolepMSSMH2(r, tbnew);

  vnew = sqrt(v1*v1+v2*v2+epsilon*epsilon-2.0*epsilon*v2);
  tbnew = (v2-epsilon)/v1;
  r.setHvev(vnew);
  r.setTanb(tbnew);

  ti1stpBck = (v2-epsilon)*doCalcTadpolepMSSMH2(r, tbnew);

  vnew = sqrt(v1*v1+v2*v2+4.0*epsilon*epsilon-4.0*epsilon*v2);
  tbnew = (v2-2.0*epsilon)/v1;
  r.setHvev(vnew);
  r.setTanb(tbnew);

  ti2stpBck = (v2-2.0*epsilon)*doCalcTadpolepMSSMH2(r, tbnew);

  deltaij(2,2) = (1.0/(12.0*epsilon))*(-1.0*ti2stpFor+8.0*ti1stpFor-8.0*ti1stpBck+ti2stpBck);


  return deltaij;
}

// pMSSM_EWSBNewtonSolver computes the solution to the 
// EWSB conditions at either tree level or one loop order using Newton's method.
// Returns the solution for the VEVs v_1, v_2 as a DoubleVector.
// Note that the provided object is taken to contain both the values of
// the input parameters and the initial guess for the VEVs 
// Inputs:
//     SoftParsMssm r = the object to calculate the solution for
//     double tb = the initial guess for the value of tan(beta)
//     int l = # loops to use (0 or 1)
//     double tol = the required tolerance to test for convergence
//     int maxIters = the maximum number of allowed iterations
//     double epsilon = value of epsilon used to numerically calculate Jacobian
//     int & sing = integer flag to indicate if singular (in which case sing is non-zero)
DoubleVector pMSSM_EWSBNewtonSolver(SoftParsMssm r, double tb, int l, double tol, int maxIters, double epsilon, int & sing)
{

  // Extract initial guess for soln
  double v = r.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;

  // Iterate using Newton's method to obtain new solution
  double v1old, v2old;
  double v1new, v2new;
  double f1old, f2old;
  double fi2stpFor, fi1stpFor, fi1stpBck, fi2stpBck;
  double df1dv1, df1dv2, df2dv1, df2dv2;
  double vold, tbold, vtmp, tbtmp;
  double detF;

  bool hasSoln = false;

  DoubleVector soln(2);

  v1old = v1;
  v2old = v2;
  vold = v;
  tbold = tb;

  double Bmu = r.displayM3Squared();
  double g1 = r.displayGaugeCoupling(1);
  double g2 = r.displayGaugeCoupling(2);
  double gbar =sqrt(g2*g2+0.6*g1*g1);
  if (l >= 1)
    {
      if (l != 1)
	{
	  cout << "Warning: only l = 0 or l = 1 loops currently supported." << endl;
	  cout << "Using l = 1 loops." << endl;
	}
    }

  for (int i = 0; i < maxIters; i++)
    {
      r.setHvev(vold);
      r.setTanb(tbold);

      // Evaluate elements of Jacobian and residual
      f1old = pMSSM_EWSBCondition1(r, tbold, l);
      f2old = pMSSM_EWSBCondition2(r, tbold, l);

      // Jacobian elements have been checked numerically elsewhere;
      // here we use the analytic formulae.
      df1dv1 = Bmu*v2old/(v1old*v1old)+0.25*gbar*gbar*v1old;
      df1dv2 = -Bmu/v1old-0.25*gbar*gbar*v2old;
      df2dv1 = -Bmu/v2old-0.25*gbar*gbar*v1old;
      df2dv2 = Bmu*v1old/(v2old*v2old)+0.25*gbar*gbar*v2old;
      if (l >= 1)
	{
	  df1dv1 = df1dv1+(1.0/v1old)*(doCalcpMSSMDeltaPrime11(r, tbold)-doCalcTadpolepMSSMH1(r,tbold));
	  df1dv2 = df1dv2+doCalcpMSSMDeltaPrime12(r, tbold)/v1old;
	  df2dv1 = df2dv1+doCalcpMSSMDeltaPrime12(r, tbold)/v2old;
	  df2dv2 = df2dv2+(1.0/v2old)*(doCalcpMSSMDeltaPrime22(r, tbold)-doCalcTadpolepMSSMH2(r,tbold));
	}
      detF = df1dv1*df2dv2-df1dv2*df2dv1;

      if (fabs(detF) < EPSTOL)
	{
	  cout << "Warning: singularity encountered in estimating EWSB conditions solution." << endl;
	  v1new = v1old;
	  v2new = v2old;
	  sing = 1;
	}
      else
	{
	  v1new = v1old + (df1dv2*f2old-df2dv2*f1old)/detF;
	  v2new = v2old + (df2dv1*f1old-df1dv1*f2old)/detF;
	}
      if (sqrt(sqr(v1new-v1old)+sqr(v2new-v2old)) < tol)
	{
	  soln(1) = v1new;
	  soln(2) = v2new;
	  hasSoln = true;
	  break;
	}
      v1old = v1new;
      v2old = v2new;

      vold = sqrt(v1old*v1old+v2old*v2old);
      tbold = v2old/v1old;
    }

  if (hasSoln == false)
    {
      cout << "Warning: max. iterations reached before solution converged to required accuracy." << endl;
      soln(1) = v1new;
      soln(2) = v2new;
    }

  return soln;

}

// pMSSM_vevNumericalDerivatives computes estimates for the partial
// derivatives of the Higgs VEVs wrt the low energy input parameters
// \mu, B,... at either tree level or one loop order. This is done
// using Newton's method and a 5-pt stencil approximation to the derivative. 
// The parameter point at which the derivative is to be evaluated, together
// with the initial estimate for the VEVs, is assumed to be contained in the provided
// MSSM model object.
// Inputs:
//     SoftParsMssm r = the object to calculate the derivatives for
//     double tb = the value of tan(beta) to use
//     int l = # loops to use(0 or 1)
//     DoubleVector epsVec = a vector indicating the direction in the low-energy
//                           parameter space in which to calculate the directional
//                           derivative. For the tree level method, epsVec must be
//                           of length 5, corresponding to the vector in parameter
//                           space (\mu, B, m_1^2, m_2^2). For
//                           the one loop method epsVec must be of length 7, containing
//                           the above 4 values followed by (m_Q^2, m_U^2, A_t). For
//                           ordinary partial derivatives, epsVec should have one element
//                           equal to unity, and all others vanishing.
//     double epsilon = the step length to use in the finite difference approximation
//     double tol = the tolerance required for Newton's method to converge
//     int maxIters = the maximum number of allowed iterations in Newton's method
//     int & sing = an integer flag indicating if singular (in which case sing is non-zero)
DoubleMatrix pMSSM_vevNumericalDerivatives(SoftParsMssm r,double tb,int l,DoubleVector epsVec,double epsilon,
					   double tol,int maxIters,int & sing)
{
  // Extract parameter values
  double mu = r.displaySusyMu();
  double Bmu = r.displayM3Squared();
  // Calculate B
  double B = Bmu/mu;
  double m1sq = r.displayMh1Squared();
  double m2sq = r.displayMh2Squared();
  double mQlsq = r.displaySoftMassSquared(mQl, 3, 3);
  double mUrsq = r.displaySoftMassSquared(mUr, 3, 3);
  double At = r.displaySoftA(UA, 3, 3);
  double yt = r.displayYukawaElement(YU, 3, 3);

  // Compute each of the derivatives in turn
  DoubleVector vi2stpFor(2), vi1stpFor(2), vi1stpBck(2), vi2stpBck(2);

  double mutmp, Btmp, m1sqtmp, m2sqtmp, mQlsqtmp, mUrsqtmp, Attmp;

  DoubleMatrix derivs(2,7);

  // Derivatives of v1 and v2 wrt mu
  mutmp = mu+2.0*epsilon*epsVec(1);
  r.setSusyMu(mutmp);
  r.setM3Squared(B*mutmp);
  vi2stpFor = pMSSM_EWSBNewtonSolver(r, tb, l, tol, maxIters, epsilon, sing);

  mutmp = mu+epsilon*epsVec(1);
  r.setSusyMu(mutmp);
  r.setM3Squared(B*mutmp);
  vi1stpFor = pMSSM_EWSBNewtonSolver(r, tb, l, tol, maxIters, epsilon, sing);

  mutmp = mu-epsilon*epsVec(1);
  r.setSusyMu(mutmp);
  r.setM3Squared(B*mutmp);
  vi1stpBck = pMSSM_EWSBNewtonSolver(r, tb, l, tol, maxIters, epsilon, sing);

  mutmp = mu-2.0*epsilon*epsVec(1);
  r.setSusyMu(mutmp);
  r.setM3Squared(B*mutmp);
  vi2stpBck = pMSSM_EWSBNewtonSolver(r, tb, l, tol, maxIters, epsilon, sing);

  derivs(1,1) = (1.0/(12.0*epsilon))*(-1.0*vi2stpFor(1)+8.0*vi1stpFor(1)-8.0*vi1stpBck(1)+vi2stpBck(1));
  derivs(2,1) =(1.0/(12.0*epsilon))*(-1.0*vi2stpFor(2)+8.0*vi1stpFor(2)-8.0*vi1stpBck(2)+vi2stpBck(2));
  r.setSusyMu(mu);
  r.setM3Squared(B*mu);

  // Derivatives of v1 and v2 wrt B
  Btmp = B+2.0*epsilon*epsVec(2);
  r.setM3Squared(Btmp*mu);
  vi2stpFor = pMSSM_EWSBNewtonSolver(r, tb, l, tol, maxIters, epsilon, sing);

  Btmp = B+epsilon*epsVec(2);
  r.setM3Squared(Btmp*mu);
  vi1stpFor = pMSSM_EWSBNewtonSolver(r, tb, l, tol, maxIters, epsilon, sing);

  Btmp = B-epsilon*epsVec(2);
  r.setM3Squared(Btmp*mu);
  vi1stpBck = pMSSM_EWSBNewtonSolver(r, tb, l, tol, maxIters, epsilon, sing);

  Btmp = B-2.0*epsilon*epsVec(2);
  r.setM3Squared(Btmp*mu);
  vi2stpBck = pMSSM_EWSBNewtonSolver(r, tb, l, tol, maxIters, epsilon, sing);

  derivs(1,2) = (1.0/(12.0*epsilon))*(-1.0*vi2stpFor(1)+8.0*vi1stpFor(1)-8.0*vi1stpBck(1)+vi2stpBck(1));
  derivs(2,2) =(1.0/(12.0*epsilon))*(-1.0*vi2stpFor(2)+8.0*vi1stpFor(2)-8.0*vi1stpBck(2)+vi2stpBck(2));
  r.setM3Squared(B*mu);

  // Derivatives of v1 and v2 wrt m_1^2
  m1sqtmp = m1sq+2.0*epsilon*epsVec(3);
  r.setMh1Squared(m1sqtmp);
  vi2stpFor = pMSSM_EWSBNewtonSolver(r, tb, l, tol, maxIters, epsilon, sing);

  m1sqtmp = m1sq+epsilon*epsVec(3);
  r.setMh1Squared(m1sqtmp);
  vi1stpFor = pMSSM_EWSBNewtonSolver(r, tb, l, tol, maxIters, epsilon, sing);

  m1sqtmp = m1sq-epsilon*epsVec(3);
  r.setMh1Squared(m1sqtmp);
  vi1stpBck = pMSSM_EWSBNewtonSolver(r, tb, l, tol, maxIters, epsilon, sing);

  m1sqtmp = m1sq-2.0*epsilon*epsVec(3);
  r.setMh1Squared(m1sqtmp);
  vi2stpBck = pMSSM_EWSBNewtonSolver(r, tb, l, tol, maxIters, epsilon, sing);

  derivs(1,3) = (1.0/(12.0*epsilon))*(-1.0*vi2stpFor(1)+8.0*vi1stpFor(1)-8.0*vi1stpBck(1)+vi2stpBck(1));
  derivs(2,3) =(1.0/(12.0*epsilon))*(-1.0*vi2stpFor(2)+8.0*vi1stpFor(2)-8.0*vi1stpBck(2)+vi2stpBck(2));
  r.setMh1Squared(m1sq);
  // Derivatives of v1 and v2 wrt m_2^2
  m2sqtmp = m2sq+2.0*epsilon*epsVec(4);
  r.setMh2Squared(m2sqtmp);
  vi2stpFor = pMSSM_EWSBNewtonSolver(r, tb, l, tol, maxIters, epsilon, sing);

  m2sqtmp = m2sq+epsilon*epsVec(4);
  r.setMh2Squared(m2sqtmp);
  vi1stpFor = pMSSM_EWSBNewtonSolver(r, tb, l, tol, maxIters, epsilon, sing);

  m2sqtmp = m2sq-epsilon*epsVec(4);
  r.setMh2Squared(m2sqtmp);
  vi1stpBck = pMSSM_EWSBNewtonSolver(r, tb, l, tol, maxIters, epsilon, sing);

  m2sqtmp = m2sq-2.0*epsilon*epsVec(4);
  r.setMh2Squared(m2sqtmp);
  vi2stpBck = pMSSM_EWSBNewtonSolver(r, tb, l, tol, maxIters, epsilon, sing);

  derivs(1,4) = (1.0/(12.0*epsilon))*(-1.0*vi2stpFor(1)+8.0*vi1stpFor(1)-8.0*vi1stpBck(1)+vi2stpBck(1));
  derivs(2,4) =(1.0/(12.0*epsilon))*(-1.0*vi2stpFor(2)+8.0*vi1stpFor(2)-8.0*vi1stpBck(2)+vi2stpBck(2));
  r.setMh2Squared(m2sq);
  if (l != 0)
    {
      if (l != 1)
	{
	  cout << "Warning: only l = 0 or l = 1 loops are currently supported." << endl;
	  cout << "Using l = 1 loop expressions." << endl;
	}
      // Estimate derivatives of v1 and v2 wrt m_Ql^2
      mQlsqtmp = mQlsq+2.0*epsilon*epsVec(5);
      r.setSoftMassElement(mQl, 3, 3, mQlsqtmp);
      vi2stpFor = pMSSM_EWSBNewtonSolver(r, tb, l, tol, maxIters, epsilon, sing);
      
      mQlsqtmp = mQlsq+epsilon*epsVec(5);
      r.setSoftMassElement(mQl, 3, 3, mQlsqtmp);
      vi1stpFor = pMSSM_EWSBNewtonSolver(r, tb, l, tol, maxIters, epsilon, sing);
      
      mQlsqtmp = mQlsq-epsilon*epsVec(5);
      r.setSoftMassElement(mQl, 3, 3, mQlsqtmp);
      vi1stpBck = pMSSM_EWSBNewtonSolver(r, tb, l, tol, maxIters, epsilon, sing);
      
      mQlsqtmp = mQlsq-2.0*epsilon*epsVec(5);
      r.setSoftMassElement(mQl, 3, 3, mQlsqtmp);
      vi2stpBck = pMSSM_EWSBNewtonSolver(r, tb, l, tol, maxIters, epsilon, sing);
      
      derivs(1,5) = (1.0/(12.0*epsilon))*(-1.0*vi2stpFor(1)+8.0*vi1stpFor(1)-8.0*vi1stpBck(1)+vi2stpBck(1));
      derivs(2,5) =(1.0/(12.0*epsilon))*(-1.0*vi2stpFor(2)+8.0*vi1stpFor(2)-8.0*vi1stpBck(2)+vi2stpBck(2));
      r.setSoftMassElement(mQl, 3, 3, mQlsq);
      // Estimate derivatives of v1 and v2 wrt m_Ur^2
      mUrsqtmp = mUrsq+2.0*epsilon*epsVec(6);
      r.setSoftMassElement(mUr, 3, 3, mUrsqtmp);
      vi2stpFor = pMSSM_EWSBNewtonSolver(r, tb, l, tol, maxIters, epsilon, sing);
      
      mUrsqtmp = mUrsq+epsilon*epsVec(6);
      r.setSoftMassElement(mUr, 3, 3, mUrsqtmp);
      vi1stpFor = pMSSM_EWSBNewtonSolver(r, tb, l, tol, maxIters, epsilon, sing);
      
      mUrsqtmp = mUrsq-epsilon*epsVec(6);
      r.setSoftMassElement(mUr, 3, 3, mUrsqtmp);
      vi1stpBck = pMSSM_EWSBNewtonSolver(r, tb, l, tol, maxIters, epsilon, sing);
      
      mUrsqtmp = mUrsq-2.0*epsilon*epsVec(6);
      r.setSoftMassElement(mUr, 3, 3, mUrsqtmp);
      vi2stpBck = pMSSM_EWSBNewtonSolver(r, tb, l, tol, maxIters, epsilon, sing);
      
      derivs(1,6) = (1.0/(12.0*epsilon))*(-1.0*vi2stpFor(1)+8.0*vi1stpFor(1)-8.0*vi1stpBck(1)+vi2stpBck(1));
      derivs(2,6) =(1.0/(12.0*epsilon))*(-1.0*vi2stpFor(2)+8.0*vi1stpFor(2)-8.0*vi1stpBck(2)+vi2stpBck(2));
      r.setSoftMassElement(mUr, 3, 3, mUrsq);      

      // Estimate derivatives of v1 and v2 wrt A_t
      Attmp = At+2.0*epsilon*epsVec(7);
      r.setTrilinearElement(UA, 3, 3, yt*Attmp);
      vi2stpFor = pMSSM_EWSBNewtonSolver(r, tb, l, tol, maxIters, epsilon, sing);
      
      Attmp = At+epsilon*epsVec(7);
      r.setTrilinearElement(UA, 3, 3, yt*Attmp);
      vi1stpFor = pMSSM_EWSBNewtonSolver(r, tb, l, tol, maxIters, epsilon, sing);
      
      Attmp = At-epsilon*epsVec(7);
      r.setTrilinearElement(UA, 3, 3, yt*Attmp);
      vi1stpBck = pMSSM_EWSBNewtonSolver(r, tb, l, tol, maxIters, epsilon, sing);
      
      Attmp = At-2.0*epsilon*epsVec(7);
      r.setTrilinearElement(UA, 3, 3, yt*Attmp);
      vi2stpBck = pMSSM_EWSBNewtonSolver(r, tb, l, tol, maxIters, epsilon, sing);
      
      derivs(1,7) = (1.0/(12.0*epsilon))*(-1.0*vi2stpFor(1)+8.0*vi1stpFor(1)-8.0*vi1stpBck(1)+vi2stpBck(1));
      derivs(2,7) =(1.0/(12.0*epsilon))*(-1.0*vi2stpFor(2)+8.0*vi1stpFor(2)-8.0*vi1stpBck(2)+vi2stpBck(2));
      r.setTrilinearElement(UA, 3, 3, yt*At);
    }

  return derivs;

}

// pMSSM_NumericalFineTunings calculates the fine tunings numerically using the above methods.
// Inputs:
//     SoftParsMssm r = pMSSM object to calculate fine tuning for
//     double tb = tan(beta)
//     int l = # loops to use (0 or 1)
//     DoubleVector & deltas = vector to store fine tunings
//     int & sing = flag if singular or poor accuracy (non-zero if problem)
//     double & detF = determinant of tuning matrix
//     double epsilon = the step length to use in the finite difference approximation
//     double tol = the tolerance required for Newton's method to converge
//     int maxIters = the maximum number of allowed iterations in Newton's method
void pMSSM_NumericalFineTunings(SoftParsMssm r, double tb, int l, DoubleVector & deltas, int & sing, double & detF,
				double epsilon, double tol, int maxIters)
{
  // Extract parameter values
  double mu = r.displaySusyMu();
  double Bmu = r.displayM3Squared();
  // Calculate B
  double B = Bmu/mu;
  double m1sq = r.displayMh1Squared();
  double m2sq = r.displayMh2Squared();
  double mQlsq = r.displaySoftMassSquared(mQl, 3, 3);
  double mUrsq = r.displaySoftMassSquared(mUr, 3, 3);
  double At = r.displaySoftA(UA, 3, 3);

  double v = r.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;

  DoubleVector epsVec(7);


  DoubleMatrix derivs(2,7);

  // Get derivatives of v1 and v2
  epsVec(1) = 1.0;
  epsVec(2) = 1.0;
  epsVec(3) = 1.0;
  epsVec(4) = 1.0;
  epsVec(5) = 1.0;
  epsVec(6) = 1.0;
  epsVec(7) = 1.0;

  derivs = pMSSM_vevNumericalDerivatives(r,tb,l,epsVec,epsilon,tol,maxIters,sing);

  // Calculate fine tunings
  deltas(1) = fabs((2.0*mu/(v*v))*(v1*derivs(1,1)+v2*derivs(2,1)));
  deltas(2) = fabs((2.0*B/(v*v))*(v1*derivs(1,2)+v2*derivs(2,2)));
  deltas(3) = fabs((2.0*m1sq/(v*v))*(v1*derivs(1,3)+v2*derivs(2,3)));
  deltas(4) = fabs((2.0*m2sq/(v*v))*(v1*derivs(1,4)+v2*derivs(2,4)));
  deltas(5) = fabs((2.0*mQlsq/(v*v))*(v1*derivs(1,5)+v2*derivs(2,5)));
  deltas(6) = fabs((2.0*mUrsq/(v*v))*(v1*derivs(1,6)+v2*derivs(2,6)));
  deltas(7) = fabs((2.0*At/(v*v))*(v1*derivs(1,7)+v2*derivs(2,7)));

  // Also calculate detF numerically to return

  double vtmp, tbtmp;
  double fi2stpFor, fi1stpFor, fi1stpBck, fi2stpBck;
  double df1dv1,df2dv2, df1dv2, df2dv1;
  // Estimate df1/dv1
  vtmp = sqrt(v1*v1+v2*v2+4.0*epsilon*epsilon+4.0*epsilon*v1);
  tbtmp = v2/(v1+2.0*epsilon);
  r.setHvev(vtmp);
  r.setTanb(tbtmp);
  
  fi2stpFor = pMSSM_EWSBCondition1(r, tbtmp, l);
  
  vtmp = sqrt(v1*v1+v2*v2+epsilon*epsilon+2.0*epsilon*v1);
  tbtmp = v2/(v1+epsilon);
  r.setHvev(vtmp);
  r.setTanb(tbtmp);
  
  fi1stpFor = pMSSM_EWSBCondition1(r, tbtmp, l);
  
  vtmp = sqrt(v1*v1+v2*v2+epsilon*epsilon-2.0*epsilon*v1);
  tbtmp = v2/(v1-epsilon);
  r.setHvev(vtmp);
  r.setTanb(tbtmp);
  
  fi1stpBck = pMSSM_EWSBCondition1(r, tbtmp, l);
  
  vtmp = sqrt(v1*v1+v2*v2+4.0*epsilon*epsilon-4.0*epsilon*v1);
  tbtmp = v2/(v1-2.0*epsilon);
  r.setHvev(vtmp);
  r.setTanb(tbtmp);
  
  fi2stpBck = pMSSM_EWSBCondition1(r, tbtmp, l);
  
  df1dv1 = (1.0/(12.0*epsilon))*(-1.0*fi2stpFor+8.0*fi1stpFor-8.0*fi1stpBck+fi2stpBck);
  
  // Estimate df1/dv2
  vtmp = sqrt(v1*v1+v2*v2+4.0*epsilon*epsilon+4.0*epsilon*v2);
  tbtmp = (v2+2.0*epsilon)/v1;
  r.setHvev(vtmp);
  r.setTanb(tbtmp);
  
  fi2stpFor = pMSSM_EWSBCondition1(r, tbtmp, l);
  
  vtmp = sqrt(v1*v1+v2*v2+epsilon*epsilon+2.0*epsilon*v2);
  tbtmp = (v2+epsilon)/v1;
  r.setHvev(vtmp);
  r.setTanb(tbtmp);
  
  fi1stpFor = pMSSM_EWSBCondition1(r, tbtmp, l);
  
  vtmp = sqrt(v1*v1+v2*v2+epsilon*epsilon-2.0*epsilon*v2);
  tbtmp = (v2-epsilon)/v1;
  r.setHvev(vtmp);
  r.setTanb(tbtmp);
  
  fi1stpBck = pMSSM_EWSBCondition1(r, tbtmp, l);
  
  vtmp = sqrt(v1*v1+v2*v2+4.0*epsilon*epsilon-4.0*epsilon*v2);
  tbtmp = (v2-2.0*epsilon)/v1;
  r.setHvev(vtmp);
  r.setTanb(tbtmp);
  
  fi2stpBck = pMSSM_EWSBCondition1(r, tbtmp, l);
  
  df1dv2 = (1.0/(12.0*epsilon))*(-1.0*fi2stpFor+8.0*fi1stpFor-8.0*fi1stpBck+fi2stpBck);
  
  // Estimate df2/dv1
  vtmp = sqrt(v1*v1+v2*v2+4.0*epsilon*epsilon+4.0*epsilon*v1);
  tbtmp = v2/(v1+2.0*epsilon);
  r.setHvev(vtmp);
  r.setTanb(tbtmp);
  
  fi2stpFor = pMSSM_EWSBCondition2(r, tbtmp, l);
  
  vtmp = sqrt(v1*v1+v2*v2+epsilon*epsilon+2.0*epsilon*v1);
  tbtmp = v2/(v1+epsilon);
  r.setHvev(vtmp);
  r.setTanb(tbtmp);
  
  fi1stpFor = pMSSM_EWSBCondition2(r, tbtmp, l);
  
  vtmp = sqrt(v1*v1+v2*v2+epsilon*epsilon-2.0*epsilon*v1);
  tbtmp = v2/(v1-epsilon);
  r.setHvev(vtmp);
  r.setTanb(tbtmp);
  
  fi1stpBck = pMSSM_EWSBCondition2(r, tbtmp, l);
  
  vtmp = sqrt(v1*v1+v2*v2+4.0*epsilon*epsilon-4.0*epsilon*v1);
  tbtmp = v2/(v1-2.0*epsilon);
  r.setHvev(vtmp);
  r.setTanb(tbtmp);
  
  fi2stpBck = pMSSM_EWSBCondition2(r, tbtmp, l);
  
  df2dv1 = (1.0/(12.0*epsilon))*(-1.0*fi2stpFor+8.0*fi1stpFor-8.0*fi1stpBck+fi2stpBck);
  
  // Estimate df2/dv2
  vtmp = sqrt(v1*v1+v2*v2+4.0*epsilon*epsilon+4.0*epsilon*v2);
  tbtmp = (v2+2.0*epsilon)/v1;
  r.setHvev(vtmp);
  r.setTanb(tbtmp);
  
  fi2stpFor = pMSSM_EWSBCondition2(r, tbtmp, l);
  
  vtmp = sqrt(v1*v1+v2*v2+epsilon*epsilon+2.0*epsilon*v2);
  tbtmp = (v2+epsilon)/v1;
  r.setHvev(vtmp);
  r.setTanb(tbtmp);
  
  fi1stpFor = pMSSM_EWSBCondition2(r, tbtmp, l);
  
  vtmp = sqrt(v1*v1+v2*v2+epsilon*epsilon-2.0*epsilon*v2);
  tbtmp = (v2-epsilon)/v1;
  r.setHvev(vtmp);
  r.setTanb(tbtmp);
  
  fi1stpBck = pMSSM_EWSBCondition2(r, tbtmp, l);
  
  vtmp = sqrt(v1*v1+v2*v2+4.0*epsilon*epsilon-4.0*epsilon*v2);
  tbtmp = (v2-2.0*epsilon)/v1;
  r.setHvev(vtmp);
  r.setTanb(tbtmp);
  
  fi2stpBck = pMSSM_EWSBCondition2(r, tbtmp, l);
  
  df2dv2 = (1.0/(12.0*epsilon))*(-1.0*fi2stpFor+8.0*fi1stpFor-8.0*fi1stpBck+fi2stpBck);
  
  // Update guess at solution
  detF = df1dv1*df2dv2-df1dv2*df2dv1;
  if (detF < EPSTOL)
    {
      cout << "Warning: singularity encountered in fine tuning calculation." << endl;
    }
}

// pMSSM_Ztunings calculates the sensitivities using the expressions given in
// arXiv:1206.5800. The tunings are returned in a double vector, which has elements
//     1 = Z_Ab^LL
//     2 = Z_Atau^LL
//     3 = Z_M1^LL
//     4 = Z_M2^LL
//     5 = Z_ML3^LL
//     6 = Z_Me3^LL
//     7 = Z_mu^LL
//     8 = Z_mAsq^TL
//     9 = Z_tb^TL
//    10 = Z_At^NLL
//    11 = Z_MQ3^NLL
//    12 = Z_Mu3^NLL
//    13 = Z_Md3^LL
//    14 = Z_M3^NLL
// Note that in the assumptions made in the pMSSM, 5 of the 19 parameters have
// zero or negligible contribution to fine tuning.
// Inputs:
//     SoftParsMssm r = pMSSM model
//     double tb = tan(beta)
//     double X = cut-off scale
DoubleVector pMSSM_Ztunings(SoftParsMssm r, double tb, double X)
{
  // Get parameter values
  double mu = r.displaySusyMu();
  double m3sq = r.displayM3Squared();
  double M1 = r.displayGaugino(1);
  double M2 = r.displayGaugino(2);
  double M3 = r.displayGaugino(3);

  double MQ3sq = r.displaySoftMassSquared(mQl, 3, 3);
  double Mu3sq = r.displaySoftMassSquared(mUr, 3, 3);
  double Md3sq = r.displaySoftMassSquared(mDr, 3, 3);

  double ML3sq = r.displaySoftMassSquared(mLl, 3, 3);
  double Me3sq = r.displaySoftMassSquared(mEr, 3, 3);

  double At = r.displaySoftA(UA, 3, 3);
  double Ab = r.displaySoftA(DA, 3, 3);
  double Atau = r.displaySoftA(EA, 3, 3);

  double yt = r.displayYukawaElement(YU, 3, 3);
  double yb = r.displayYukawaElement(YD, 3, 3);
  double ytau = r.displayYukawaElement(YE, 3, 3);

  double g1 = r.displayGaugeCoupling(1);
  double g2 = r.displayGaugeCoupling(2);
  double g3 = r.displayGaugeCoupling(3);
  double gd = sqrt(3.0/5.0)*g1;//sqrt(0.6)*g1;
  double gbar = sqrt(g2*g2+gd*gd);
  double v = r.displayHvev();
  double Mz = gbar*v/2.0;

  cout << "g' = " << gd << endl;
  cout << "g2 = " << g2 << endl;
  cout << "g3 = " << g3 << endl;
  cout << "tb = " << tb << endl;
  cout << "yt = " << yt << endl;
  cout << "yb = " << yb << endl;
  cout << "ytau = " << ytau << endl;
  cout << "M_Z = " << Mz << endl;


  // Evaluate the various Z_i and return as a vector
  DoubleVector Z(16);
  //     1 = Z_Ab^LL
  Z(1) = (3.0*X/(2.0*PI*PI))*yb*yb*Ab*Ab/(Mz*Mz*(tb*tb-1.0));

  //     2 = Z_Atau^LL
  Z(2) = (X/(2.0*PI*PI))*ytau*ytau*Atau*Atau/(Mz*Mz*(tb*tb-1.0));

  //     3 = Z_M1^LL
  Z(3) = (X/(2.0*PI*PI))*gd*gd*M1*M1/(Mz*Mz);

  //     4 = Z_M2^LL
  Z(4) = (X/(2.0*PI*PI))*3.0*g2*g2*M2*M2/(Mz*Mz);

  //     5 = Z_ML3^LL
  Z(5) = (X/(4.0*PI*PI*Mz*Mz))*(ML3sq/(tb*tb-1.0))*(2.0*ytau*ytau+gd*gd*(1.0+tb*tb));

  //     6 = Z_Me3^LL
  Z(6) = (X/(4.0*PI*PI*Mz*Mz))*(Me3sq/(tb*tb-1.0))*(2.0*ytau*ytau-gd*gd*(1.0+tb*tb));

  //     7 = Z_mu^LL
  double s2b = 2.0*tb/(1.0+tb*tb);
  double mAsq = 2.0*m3sq/s2b;
  double t2b = 2.0*tb/(1.0-tb*tb);
  double ZmuTL = (4.0*mu*mu/(Mz*Mz))*(1.0+((mAsq+Mz*Mz)/mAsq)*t2b*t2b);
  Z(7) = ZmuTL+ZmuTL*(1.0+(X/(16.0*PI*PI))*(3.0*yt*yt+3.0*yb*yb+ytau*ytau-3.0*g2*g2-gd*gd));

  // At the moment I have pretty much no idea what tree level expressions
  // are being used to calculate the tuning in m_A and tan(beta) at tree level,
  // so I am using m_3^2 = b and including tree level tunings for the soft Higgs
  // masses m_H_u^2 and m_H_d^2 presented in hep-ph/0812.0980v3.

  //     8 = Z_b^TL
  Z(8) = (1.0+mAsq/(Mz*Mz))*sqr(2.0*tb/(1.0-tb*tb));

  //     9 = Z_tb^TL
  double c2b = (1.0-tb*tb)/(1.0+tb*tb);
  Z(9) = (4.0*tb*tb*(Mz*Mz+mAsq))/(Mz*Mz*(1.0+tb*tb*tb*tb));

  //    10 = Z_At^NLL
  double ZAtLL = (3.0*X*yt*yt*At*At*(-1.0*tb*tb))/(2.0*PI*PI*Mz*Mz*(tb*tb-1.0));
  double T1 = yb*yb*(At+Ab);
  double T2 = 12.0*yt*yt*At+yb*yb*(At+2.0*Ab);
  double T3 = (4.0/3.0)*(4.0*g3*g3*(At-M3)+(gd*gd/3.0)*(At-M1));
  double ZAtNLL = (24.0*X*X/(sqr(16.0*PI*PI)))*(At/(Mz*Mz))*(-yt*yt/(tb*tb-1))*(T1+tb*tb*(T3-T2));

  // Maybe this will work?
  // double dmHd2dAt = ((X*X)/(sqr(16*PI*PI)))*12.0*yt*yt*yb*yb*(At+Ab);
  // double dmHu2dAt = -((X/(16*PI*PI))*12.0*yt*yt*At+((X*X)/(sqr(16*PI*PI)))*(2.0*(-72.0*sqr(yt*yt)-6.0*yb*yb*yt*yt+32.0*yt*yt*g3*g3
  // 										 +(8.0/3.0)*yt*yt*gd*gd)*At-12.0*yt*yt*yb*yb*Ab
  // 									    -64.0*yt*yt*g3*g3*M3-(16.0/3.0)*yt*yt*gd*gd*M1));
  // double extraTerm = (At/(Mz*Mz))*t2b*t2b*((mAsq+Mz*Mz)/mAsq)*(dmHd2dAt+dmHu2dAt);

  // At this stage it doesn't look like it...

  Z(10) = ZAtLL + ZAtNLL;



  //    11 = Z_MQ3^NLL
  double ZmQ3LL = X*MQ3sq*(12.0*yb*yb-2.0*gd*gd-tb*tb*(12.0*yt*yt+2.0*gd*gd))/(4.0*PI*PI*Mz*Mz*(tb*tb-1.0));
  double C1 = 2.0*yb*yb*(32.0*g3*g3-4.0*gd*gd/3.0-36.0*yb*yb-12.0*yt*yt);
  double C2 = 2.0*yt*yt*(32.0*g3*g3+8.0*gd*gd/3.0-36.0*yt*yt-12.0*yb*yb);
  Z(11) = ZmQ3LL+(X*X/(128.0*PI*PI*PI*PI*Mz*Mz))*(MQ3sq/(tb*tb-1.0))*(C1-C2*tb*tb);


  //    12 = Z_Mu3^NLL
  double Zmu3LL = X*Mu3sq*(4.0*gd*gd-tb*tb*(12.0*yt*yt-4.0*gd*gd))/(8.0*PI*PI*Mz*Mz*(tb*tb-1.0));
  Z(12) = Zmu3LL+(X*X/(128.0*PI*PI*PI*PI*Mz*Mz))*(Mu3sq/(tb*tb-1.0))*(-12.0*yb*yb*yt*yt
								      -2.0*yt*yt*tb*tb*(32.0*g3*g3+8.0*gd*gd/3.0-36.0*yt*yt-6.0*yb*yb));

  //    13 = Z_Md3^LL
  Z(13) = (X/(4.0*PI*PI*Mz*Mz))*(Md3sq/(tb*tb-1.0))*(6.0*yb*yb-gd*gd);

  //    14 = Z_M3^NLL
  double alphas = g3*g3/(4.0*PI);
  Z(14) = (2.0*alphas*X*X/((3.0*PI*PI*PI)*(tb*tb-1.0)))*(M3/(Mz*Mz))*(-1.0*yb*yb*(2.0*M3-Ab)+tb*tb*yt*yt*(2.0*M3-At));

  //    15 = m_Hd^2^TL
  Z(15) = fabs(-0.5*(1.0-tb*tb)/(1.0+tb*tb)+mAsq*tb*tb/(Mz*Mz*(1.0+tb*tb))-mu*mu/(Mz*Mz))
    *fabs(1.0+(1.0+tb*tb)/(1.0-tb*tb)+((mAsq+Mz*Mz)/mAsq)*sqr(2.0*tb/(1.0-tb*tb)));

  //    16 = m_Hu^2^TL
  Z(16) = fabs(0.5*(1.0-tb*tb)/(1.0+tb*tb)+mAsq/(Mz*Mz*(1.0+tb*tb))-mu*mu/(Mz*Mz))
    *fabs(1.0-(1.0+tb*tb)/(1.0-tb*tb)+((mAsq+Mz*Mz)/mAsq)*sqr(2.0*tb/(1.0-tb*tb)));

  return Z;
}
