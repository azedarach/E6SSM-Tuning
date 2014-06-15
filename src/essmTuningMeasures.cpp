/*
  higgsNumericalTuningMeasures.cpp provides the function
  implementations for functions computing the values of
  the EWSB conditions (with or without loop corrections)
  and the numerical tuning measures.
  Author: Dylan Harries
  Date created: 30/9/2013
 */

#include "essmTuningMeasures.h"

// EWSBConditioni_TreeLevel, i = 1, 2, 3, calculates the value
// of the EWSB condition for m_H_i^2 (i=3 is m_s^2) at tree level (no tadpole
// corrections). Provided to allow for comparison with previous codes done
// at tree level.
// Inputs:
//     SoftParsEssm essmSusy = object containing the ESSM parameters
//     double s = value of the VEV s to use
//     double tb = value of tan(beta) to use
double EWSBCondition1_TreeLevel(SoftParsEssm essmSusy, double s, double tb)
{
  DoubleVector mH1Sq_vec = essmSusy.displayMh1Squared();
  double m1Sq = mH1Sq_vec.display(3);
  DoubleVector Alambda_vec = essmSusy.displayA_lambda();
  double Alambda = Alambda_vec.display(3);
  DoubleVector lambda_vec = essmSusy.displaylambda();
  double lambda = lambda_vec.display(3);

  double gdash = sqrt(3.0/5.0)*essmSusy.displayGaugeCoupling(1);
  double g_2 = essmSusy.displayGaugeCoupling(2);
  double gbar = sqrt(gdash*gdash+g_2*g_2);
  double gdash_1 = essmSusy.displaygdash_1();

  double Qtilde_1 = essmSusy.displayQH1tilde();
  double Qtilde_2 = essmSusy.displayQH2tilde();
  double Qtilde_s = essmSusy.displayQStilde();

  double cSqb = 1.0/(1.0+tb*tb);
  double sSqb = (tb*tb)/(1.0+tb*tb);
  double c2b = (1.0-tb*tb)/(1.0+tb*tb);

  double v = essmSusy.displayHvev();
  double v2Sq = v*v*sSqb;
  double v1 = v/sqrt(1.0+tb*tb);

  double fnVal = m1Sq + 0.5*lambda*lambda*(v2Sq+s*s)-lambda*Alambda*s*tb/sqrt(2.0)
    +gbar*gbar*v*v*c2b/8.0 + (0.5*Qtilde_1*gdash_1*gdash_1*
			      (Qtilde_1*v*v*cSqb + Qtilde_2*v*v*sSqb + Qtilde_s*s*s));
  fnVal = v1*fnVal;

  return fnVal;
}

double EWSBCondition2_TreeLevel(SoftParsEssm essmSusy, double s, double tb)
{
  DoubleVector mH2Sq_vec = essmSusy.displayMh2Squared();
  double m2Sq = mH2Sq_vec.display(3);
  DoubleVector Alambda_vec = essmSusy.displayA_lambda();
  double Alambda = Alambda_vec.display(3);
  DoubleVector lambda_vec = essmSusy.displaylambda();
  double lambda = lambda_vec.display(3);

  double gdash = sqrt(3.0/5.0)*essmSusy.displayGaugeCoupling(1);
  double g_2 = essmSusy.displayGaugeCoupling(2);
  double gbar = sqrt(gdash*gdash+g_2*g_2);
  double gdash_1 = essmSusy.displaygdash_1();

  double Qtilde_1 = essmSusy.displayQH1tilde();
  double Qtilde_2 = essmSusy.displayQH2tilde();
  double Qtilde_s = essmSusy.displayQStilde();

  double cSqb = 1.0/(1.0+tb*tb);
  double sSqb = (tb*tb)/(1.0+tb*tb);
  double c2b = (1.0-tb*tb)/(1.0+tb*tb);

  double v = essmSusy.displayHvev();
  double v1Sq = v*v*cSqb;
  double v2 = v*tb/sqrt(1.0+tb*tb);

  double fnVal = m2Sq + 0.5*lambda*lambda*(v1Sq+s*s)-lambda*Alambda*s/(sqrt(2.0)*tb)
    -gbar*gbar*v*v*c2b/8.0 + (0.5*Qtilde_2*gdash_1*gdash_1*
			      (Qtilde_1*v*v*cSqb + Qtilde_2*v*v*sSqb + Qtilde_s*s*s));
  fnVal = v2*fnVal;

  return fnVal;
}

double EWSBCondition3_TreeLevel(SoftParsEssm essmSusy, double s, double tb)
{
  DoubleVector msSq_vec = essmSusy.displaymSsq();
  double msSq = msSq_vec.display(3);
  DoubleVector Alambda_vec = essmSusy.displayA_lambda();
  double Alambda = Alambda_vec.display(3);
  DoubleVector lambda_vec = essmSusy.displaylambda();
  double lambda = lambda_vec.display(3);

  double gdash = sqrt(3.0/5.0)*essmSusy.displayGaugeCoupling(1);
  double g_2 = essmSusy.displayGaugeCoupling(2);
  double gbar = sqrt(gdash*gdash+g_2*g_2);
  double gdash_1 = essmSusy.displaygdash_1();

  double Qtilde_1 = essmSusy.displayQH1tilde();
  double Qtilde_2 = essmSusy.displayQH2tilde();
  double Qtilde_s = essmSusy.displayQStilde();

  double cSqb = 1.0/(1.0+tb*tb);
  double sSqb = (tb*tb)/(1.0+tb*tb);
  double c2b = (1.0-tb*tb)/(1.0+tb*tb);
  double s2b = 2.0*tb/(1.0+tb*tb);

  double v = essmSusy.displayHvev();

  double fnVal = s*msSq + 0.5*lambda*lambda*v*v*s-lambda*Alambda*v*v*s2b/(2.0*sqrt(2.0))
    +0.5*Qtilde_s*s*gdash_1*gdash_1*(Qtilde_1*v*v*cSqb + Qtilde_2*v*v*sSqb + Qtilde_s*s*s);

  

  return fnVal;
}


// EWSBConditioni_Loops, i = 1, 2, 3, calculates the value
// of the EWSB condition for m_H_i^2 (i=3 is m_s^2), including the 
// contributions from the tadpole corrections. Evaluates the tadpoles
// at the renormalisation scale chosen in the model, and includes 
// contributions from D fields.
// Inputs:
//     SoftParsEssm essmSusy = object containing the ESSM parameters
//     double s = value of the VEV s
//     double tb = value of tan(beta)
double EWSBCondition1_OneLoop(SoftParsEssm essmSusy, double s, double tb)
{
  double v = essmSusy.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);

  double fnVal = EWSBCondition1_TreeLevel(essmSusy, s, tb);

  // Add in appropriate loop correction
  fnVal = fnVal-v1*doCalcTadpoleESSMH1(essmSusy, s, tb);

  return fnVal;
}

double EWSBCondition2_OneLoop(SoftParsEssm essmSusy, double s, double tb)
{
  double v = essmSusy.displayHvev();
  double v2 = v*tb/sqrt(1.0+tb*tb);

  double fnVal = EWSBCondition2_TreeLevel(essmSusy, s, tb);

  // Add in appropriate loop correction
  fnVal = fnVal-v2*doCalcTadpoleESSMH2(essmSusy, s, tb);

  return fnVal;
}

double EWSBCondition3_OneLoop(SoftParsEssm essmSusy, double s, double tb)
{
  double fnVal = EWSBCondition3_TreeLevel(essmSusy, s, tb);

  // Add in appropriate loop correction
  fnVal = fnVal-s*doCalcTadpolesESSMS(essmSusy, s, tb);

  return fnVal;
}


// As above, but evaluates the tadpoles at the scale of M_t, and includes 
// contributions from D fields.
// Inputs:
//     SoftParsEssm essmSusy = object containing the ESSM parameters
//     double s = value of the VEV s
//     double tb = value of tan(beta)
double EWSBCondition1_OneLoop_atMt(SoftParsEssm essmSusy, double s, double tb)
{
  double v = essmSusy.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);

  double fnVal = EWSBCondition1_TreeLevel(essmSusy, s, tb);

  // Add in appropriate loop correction
  fnVal = fnVal-v1*doCalcTadpoleESSMH1_atMt(essmSusy, s, tb);

  return fnVal;
}

double EWSBCondition2_OneLoop_atMt(SoftParsEssm essmSusy, double s, double tb)
{
  double v = essmSusy.displayHvev();
  double v2 = v*tb/sqrt(1.0+tb*tb);

  double fnVal = EWSBCondition2_TreeLevel(essmSusy, s, tb);

  // Add in appropriate loop correction
  fnVal = fnVal-v2*doCalcTadpoleESSMH2_atMt(essmSusy, s, tb);

  return fnVal;
}

double EWSBCondition3_OneLoop_atMt(SoftParsEssm essmSusy, double s, double tb)
{
  double fnVal = EWSBCondition3_TreeLevel(essmSusy, s, tb);

  // Add in appropriate loop correction
  fnVal = fnVal-s*doCalcTadpolesESSMS_atMt(essmSusy, s, tb);

  return fnVal;
}

// As above, but evaluates the tadpoles at the scale of M_t, and neglects
// contributions from D fields.
// Inputs:
//     SoftParsEssm essmSusy = object containing the ESSM parameters
//     double s = value of the VEV s
//     double tb = value of tan(beta)
double EWSBCondition1_OneLoop_Roman(SoftParsEssm essmSusy, double s, double tb)
{
  double v = essmSusy.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);

  double fnVal = EWSBCondition1_TreeLevel(essmSusy, s, tb);

  // Add in appropriate loop correction
  fnVal = fnVal-v1*doCalcTadpoleESSMH1_Roman(essmSusy, s, tb);

  return fnVal;
}

double EWSBCondition2_OneLoop_Roman(SoftParsEssm essmSusy, double s ,double tb)
{
  double v = essmSusy.displayHvev();
  double v2 = v*tb/sqrt(1.0+tb*tb);

  double fnVal = EWSBCondition2_TreeLevel(essmSusy, s, tb);

  // Add in appropriate loop correction
  fnVal = fnVal-v2*doCalcTadpoleESSMH2_Roman(essmSusy, s, tb);

  return fnVal;
}

double EWSBCondition3_OneLoop_Roman(SoftParsEssm essmSusy, double s, double tb)
{
  double fnVal = EWSBCondition3_TreeLevel(essmSusy, s, tb);

  // Add in appropriate loop correction
  fnVal = fnVal-s*doCalcTadpolesESSMS_Roman(essmSusy, s, tb);

  return fnVal;
}

// As above, but evaluates the tadpoles at the scale Q defined in the model, and 
// neglects contributions from D fields.
// Inputs:
//     SoftParsEssm essmSusy = object containing the ESSM parameters
//     double s = value of the VEV s
//     double tb = value of tan(beta)
double EWSBCondition1_OneLoop_Roman_atQ(SoftParsEssm essmSusy, double s, double tb)
{
  double v = essmSusy.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);

  double fnVal = EWSBCondition1_TreeLevel(essmSusy, s, tb);

  // Add in appropriate loop correction
  fnVal = fnVal-v1*doCalcTadpoleESSMH1_Roman_atQ(essmSusy, s, tb);

  return fnVal;
}

double EWSBCondition2_OneLoop_Roman_atQ(SoftParsEssm essmSusy, double s, double tb)
{
  double v = essmSusy.displayHvev();
  double v2 = v*tb/sqrt(1.0+tb*tb);

  double fnVal = EWSBCondition2_TreeLevel(essmSusy, s, tb);

  // Add in appropriate loop correction
  fnVal = fnVal-v2*doCalcTadpoleESSMH2_Roman_atQ(essmSusy, s, tb);

  return fnVal;
}

double EWSBCondition3_OneLoop_Roman_atQ(SoftParsEssm essmSusy, double s, double tb)
{
  double fnVal = EWSBCondition3_TreeLevel(essmSusy, s, tb);

  // Add in appropriate loop correction
  fnVal = fnVal-s*doCalcTadpolesESSMS_Roman_atQ(essmSusy, s, tb);

  return fnVal;
}

// higgsImplementEWSBConstraints_TreeLevel imposes the EWSB constraints
// on the given ESSM model, setting the value of the soft masses
// m_H_1^2, m_H_2^2 and m_s^2 to satisfy the three EWSB conditions
// at tree level. Allows for comparison with previous codes done at
// tree level.
// Inputs:
//     SoftParsEssm & essmSusy = the object to impose the EWSB constraints on.
//     double s = the value of the VEV s
//     double tb = the value of tan(beta)
void ImplementEWSBConstraints_TreeLevel(SoftParsEssm & essmSusy, double s, double tb)
{
  DoubleVector Alambda_vec = essmSusy.displayA_lambda();
  double Alambda = Alambda_vec.display(3);
  DoubleVector lambda_vec = essmSusy.displaylambda();
  double lambda = lambda_vec.display(3);

  DoubleVector mH1Sq_vec = essmSusy.displayMh1Squared();
  DoubleVector mH2Sq_vec = essmSusy.displayMh2Squared();
  DoubleVector mSSq_vec = essmSusy.displaymSsq();

  double gdash = sqrt(3.0/5.0)*essmSusy.displayGaugeCoupling(1);
  double g_2 = essmSusy.displayGaugeCoupling(2);
  double gbar = sqrt(gdash*gdash+g_2*g_2);
  double gdash_1 = essmSusy.displaygdash_1();

  double Qtilde_1 = essmSusy.displayQH1tilde();
  double Qtilde_2 = essmSusy.displayQH2tilde();
  double Qtilde_s = essmSusy.displayQStilde();

  double cSqb = 1.0/(1.0+tb*tb);
  double sSqb = (tb*tb)/(1.0+tb*tb);
  double c2b = (1.0-tb*tb)/(1.0+tb*tb);
  double s2b = 2.0*tb/(1.0+tb*tb);

  double v = essmSusy.displayHvev();

  // Calculate the values of the soft squared masses using the three EWSB
  // conditions.
  double msSq = (1.0/s)*(lambda*Alambda*v*v*s2b/(2.0*sqrt(2.0))-0.5*lambda*lambda*v*v*s
			 -0.5*gdash_1*gdash_1*Qtilde_s*s*(Qtilde_1*v*v*cSqb+Qtilde_2*v*v*sSqb+Qtilde_s*s*s));
  double m1Sq = lambda*Alambda*s*tb/sqrt(2.0)-0.5*lambda*lambda*(s*s+v*v*sSqb)-gbar*gbar*v*v*c2b/8.0
    -0.5*Qtilde_1*gdash_1*gdash_1*(Qtilde_1*v*v*cSqb+Qtilde_2*v*v*sSqb+Qtilde_s*s*s);
  double m2Sq = lambda*Alambda*s/(sqrt(2.0)*tb)-0.5*lambda*lambda*(s*s+v*v*cSqb)+gbar*gbar*v*v*c2b/8.0
    -0.5*Qtilde_2*gdash_1*gdash_1*(Qtilde_1*v*v*cSqb+Qtilde_2*v*v*sSqb+Qtilde_s*s*s);

  // Update the three masses in the object, doesn't change any other parameters.
  mH1Sq_vec.set(3, m1Sq);
  mH2Sq_vec.set(3, m2Sq);
  mSSq_vec.set(3, msSq);

  essmSusy.setMh1Squared(mH1Sq_vec);
  essmSusy.setMh2Squared(mH2Sq_vec);
  essmSusy.setmSsq(mSSq_vec);
}

// higgsImplementEWSBConstraints_Loops imposes the EWSB constraints
// on the given ESSM model, setting the value of the soft masses
// m_H_1^2, m_H_2^2 and m_s^2 to satisfy the three EWSB conditions
// after including contributions from tadpole diagrams. The tadpole
// corrections used include contributions from the D fields and are 
// evaluated at the renormalisation scale Q.
// Inputs:
//     SoftParsEssm & essmSusy = the object to impose the EWSB constraints on.
//     double s = the value of the VEV s
//     double tb = the value of tan(beta)
void ImplementEWSBConstraints_OneLoop(SoftParsEssm & essmSusy, double s, double tb)
{
  DoubleVector Alambda_vec = essmSusy.displayA_lambda();
  double Alambda = Alambda_vec.display(3);
  DoubleVector lambda_vec = essmSusy.displaylambda();
  double lambda = lambda_vec.display(3);

  DoubleVector mH1Sq_vec = essmSusy.displayMh1Squared();
  DoubleVector mH2Sq_vec = essmSusy.displayMh2Squared();
  DoubleVector mSSq_vec = essmSusy.displaymSsq();

  double gdash = sqrt(3.0/5.0)*essmSusy.displayGaugeCoupling(1);
  double g_2 = essmSusy.displayGaugeCoupling(2);
  double gbar = sqrt(gdash*gdash+g_2*g_2);
  double gdash_1 = essmSusy.displaygdash_1();

  double Qtilde_1 = essmSusy.displayQH1tilde();
  double Qtilde_2 = essmSusy.displayQH2tilde();
  double Qtilde_s = essmSusy.displayQStilde();

  double cSqb = 1.0/(1.0+tb*tb);
  double sSqb = (tb*tb)/(1.0+tb*tb);
  double c2b = (1.0-tb*tb)/(1.0+tb*tb);
  double s2b = 2.0*tb/(1.0+tb*tb);

  double v = essmSusy.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;

  // Calculate the values of the soft squared masses using the three EWSB
  // conditions at tree level.
  double msSq = (1.0/s)*(lambda*Alambda*v*v*s2b/(2.0*sqrt(2.0))-0.5*lambda*lambda*v*v*s
			 -0.5*gdash_1*gdash_1*Qtilde_s*s*(Qtilde_1*v*v*cSqb+Qtilde_2*v*v*sSqb+Qtilde_s*s*s));
  double m1Sq = lambda*Alambda*s*tb/sqrt(2.0)-0.5*lambda*lambda*(s*s+v*v*sSqb)-gbar*gbar*v*v*c2b/8.0
    -0.5*Qtilde_1*gdash_1*gdash_1*(Qtilde_1*v*v*cSqb+Qtilde_2*v*v*sSqb+Qtilde_s*s*s);
  double m2Sq = lambda*Alambda*s/(sqrt(2.0)*tb)-0.5*lambda*lambda*(s*s+v*v*cSqb)+gbar*gbar*v*v*c2b/8.0
    -0.5*Qtilde_2*gdash_1*gdash_1*(Qtilde_1*v*v*cSqb+Qtilde_2*v*v*sSqb+Qtilde_s*s*s);

  // Add in appropriate loop corrections (note sign!).
  msSq = msSq + doCalcTadpolesESSMS(essmSusy, s, tb);
  m1Sq = m1Sq + doCalcTadpoleESSMH1(essmSusy, s, tb);
  m2Sq = m2Sq + doCalcTadpoleESSMH2(essmSusy, s, tb);

  // Update the three masses in the object, doesn't change any other parameters.
  mH1Sq_vec.set(3, m1Sq);
  mH2Sq_vec.set(3, m2Sq);
  mSSq_vec.set(3, msSq);

  essmSusy.setMh1Squared(mH1Sq_vec);
  essmSusy.setMh2Squared(mH2Sq_vec);
  essmSusy.setmSsq(mSSq_vec);
}

// As above, except the tadpole corrections used include contributions from the D fields 
// and are evaluated at the renormalisation scale M_t.
// Inputs:
//     SoftParsEssm & essmSusy = the object to impose the EWSB constraints on.
//     double s = the value of the VEV s
//     double tb = the value of tan(beta)
void ImplementEWSBConstraints_OneLoop_atMt(SoftParsEssm & essmSusy, double s, double tb)
{
  DoubleVector Alambda_vec = essmSusy.displayA_lambda();
  double Alambda = Alambda_vec.display(3);
  DoubleVector lambda_vec = essmSusy.displaylambda();
  double lambda = lambda_vec.display(3);

  DoubleVector mH1Sq_vec = essmSusy.displayMh1Squared();
  DoubleVector mH2Sq_vec = essmSusy.displayMh2Squared();
  DoubleVector mSSq_vec = essmSusy.displaymSsq();

  double gdash = sqrt(3.0/5.0)*essmSusy.displayGaugeCoupling(1);
  double g_2 = essmSusy.displayGaugeCoupling(2);
  double gbar = sqrt(gdash*gdash+g_2*g_2);
  double gdash_1 = essmSusy.displaygdash_1();

  double Qtilde_1 = essmSusy.displayQH1tilde();
  double Qtilde_2 = essmSusy.displayQH2tilde();
  double Qtilde_s = essmSusy.displayQStilde();

  double cSqb = 1.0/(1.0+tb*tb);
  double sSqb = (tb*tb)/(1.0+tb*tb);
  double c2b = (1.0-tb*tb)/(1.0+tb*tb);
  double s2b = 2.0*tb/(1.0+tb*tb);

  double v = essmSusy.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;

  // Calculate the values of the soft squared masses using the three EWSB
  // conditions at tree level.
  double msSq = (1.0/s)*(lambda*Alambda*v*v*s2b/(2.0*sqrt(2.0))-0.5*lambda*lambda*v*v*s
			 -0.5*gdash_1*gdash_1*Qtilde_s*s*(Qtilde_1*v*v*cSqb+Qtilde_2*v*v*sSqb+Qtilde_s*s*s));
  double m1Sq = lambda*Alambda*s*tb/sqrt(2.0)-0.5*lambda*lambda*(s*s+v*v*sSqb)-gbar*gbar*v*v*c2b/8.0
    -0.5*Qtilde_1*gdash_1*gdash_1*(Qtilde_1*v*v*cSqb+Qtilde_2*v*v*sSqb+Qtilde_s*s*s);
  double m2Sq = lambda*Alambda*s/(sqrt(2.0)*tb)-0.5*lambda*lambda*(s*s+v*v*cSqb)+gbar*gbar*v*v*c2b/8.0
    -0.5*Qtilde_2*gdash_1*gdash_1*(Qtilde_1*v*v*cSqb+Qtilde_2*v*v*sSqb+Qtilde_s*s*s);

  // Add in appropriate loop corrections (note sign!).
  msSq = msSq + doCalcTadpolesESSMS_atMt(essmSusy, s, tb);
  m1Sq = m1Sq + doCalcTadpoleESSMH1_atMt(essmSusy, s, tb);
  m2Sq = m2Sq + doCalcTadpoleESSMH2_atMt(essmSusy, s, tb);

  // Update the three masses in the object, doesn't change any other parameters.
  mH1Sq_vec.set(3, m1Sq);
  mH2Sq_vec.set(3, m2Sq);
  mSSq_vec.set(3, msSq);

  essmSusy.setMh1Squared(mH1Sq_vec);
  essmSusy.setMh2Squared(mH2Sq_vec);
  essmSusy.setmSsq(mSSq_vec);
}

// As above, except the tadpole corrections used doesn't include contributions from the D fields 
// and are evaluated at the renormalisation scale M_t.
// Inputs:
//     SoftParsEssm & essmSusy = the object to impose the EWSB constraints on.
//     double s = the value of the VEV s
//     double tb = the value of tan(beta)
void ImplementEWSBConstraints_OneLoop_Roman(SoftParsEssm & essmSusy, double s, double tb)
{
  DoubleVector Alambda_vec = essmSusy.displayA_lambda();
  double Alambda = Alambda_vec.display(3);
  DoubleVector lambda_vec = essmSusy.displaylambda();
  double lambda = lambda_vec.display(3);

  DoubleVector mH1Sq_vec = essmSusy.displayMh1Squared();
  DoubleVector mH2Sq_vec = essmSusy.displayMh2Squared();
  DoubleVector mSSq_vec = essmSusy.displaymSsq();

  double gdash = sqrt(3.0/5.0)*essmSusy.displayGaugeCoupling(1);
  double g_2 = essmSusy.displayGaugeCoupling(2);
  double gbar = sqrt(gdash*gdash+g_2*g_2);
  double gdash_1 = essmSusy.displaygdash_1();

  double Qtilde_1 = essmSusy.displayQH1tilde();
  double Qtilde_2 = essmSusy.displayQH2tilde();
  double Qtilde_s = essmSusy.displayQStilde();

  double cSqb = 1.0/(1.0+tb*tb);
  double sSqb = (tb*tb)/(1.0+tb*tb);
  double c2b = (1.0-tb*tb)/(1.0+tb*tb);
  double s2b = 2.0*tb/(1.0+tb*tb);

  double v = essmSusy.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;

  // Calculate the values of the soft squared masses using the three EWSB
  // conditions at tree level.
  double msSq = (1.0/s)*(lambda*Alambda*v*v*s2b/(2.0*sqrt(2.0))-0.5*lambda*lambda*v*v*s
			 -0.5*gdash_1*gdash_1*Qtilde_s*s*(Qtilde_1*v*v*cSqb+Qtilde_2*v*v*sSqb+Qtilde_s*s*s));
  double m1Sq = lambda*Alambda*s*tb/sqrt(2.0)-0.5*lambda*lambda*(s*s+v*v*sSqb)-gbar*gbar*v*v*c2b/8.0
    -0.5*Qtilde_1*gdash_1*gdash_1*(Qtilde_1*v*v*cSqb+Qtilde_2*v*v*sSqb+Qtilde_s*s*s);
  double m2Sq = lambda*Alambda*s/(sqrt(2.0)*tb)-0.5*lambda*lambda*(s*s+v*v*cSqb)+gbar*gbar*v*v*c2b/8.0
    -0.5*Qtilde_2*gdash_1*gdash_1*(Qtilde_1*v*v*cSqb+Qtilde_2*v*v*sSqb+Qtilde_s*s*s);

  // Add in appropriate loop corrections (note sign!).
  msSq = msSq + doCalcTadpolesESSMS_Roman(essmSusy, s, tb);
  m1Sq = m1Sq + doCalcTadpoleESSMH1_Roman(essmSusy, s, tb);
  m2Sq = m2Sq + doCalcTadpoleESSMH2_Roman(essmSusy, s, tb);

  // Update the three masses in the object, doesn't change any other parameters.
  mH1Sq_vec.set(3, m1Sq);
  mH2Sq_vec.set(3, m2Sq);
  mSSq_vec.set(3, msSq);

  essmSusy.setMh1Squared(mH1Sq_vec);
  essmSusy.setMh2Squared(mH2Sq_vec);
  essmSusy.setmSsq(mSSq_vec);
}

// As above, except the tadpole corrections used don't include contributions from the D fields 
// and are evaluated at the renormalisation scale Q.
// Inputs:
//     SoftParsEssm & essmSusy = the object to impose the EWSB constraints on.
//     double s = the value of the VEV s
//     double tb = the value of tan(beta)
void ImplementEWSBConstraints_OneLoop_Roman_atQ(SoftParsEssm & essmSusy, double s, double tb)
{
  DoubleVector Alambda_vec = essmSusy.displayA_lambda();
  double Alambda = Alambda_vec.display(3);
  DoubleVector lambda_vec = essmSusy.displaylambda();
  double lambda = lambda_vec.display(3);

  DoubleVector mH1Sq_vec = essmSusy.displayMh1Squared();
  DoubleVector mH2Sq_vec = essmSusy.displayMh2Squared();
  DoubleVector mSSq_vec = essmSusy.displaymSsq();

  double gdash = sqrt(3.0/5.0)*essmSusy.displayGaugeCoupling(1);
  double g_2 = essmSusy.displayGaugeCoupling(2);
  double gbar = sqrt(gdash*gdash+g_2*g_2);
  double gdash_1 = essmSusy.displaygdash_1();

  double Qtilde_1 = essmSusy.displayQH1tilde();
  double Qtilde_2 = essmSusy.displayQH2tilde();
  double Qtilde_s = essmSusy.displayQStilde();

  double cSqb = 1.0/(1.0+tb*tb);
  double sSqb = (tb*tb)/(1.0+tb*tb);
  double c2b = (1.0-tb*tb)/(1.0+tb*tb);
  double s2b = 2.0*tb/(1.0+tb*tb);

  double v = essmSusy.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;

  // Calculate the values of the soft squared masses using the three EWSB
  // conditions at tree level.
  double msSq = (1.0/s)*(lambda*Alambda*v*v*s2b/(2.0*sqrt(2.0))-0.5*lambda*lambda*v*v*s
			 -0.5*gdash_1*gdash_1*Qtilde_s*s*(Qtilde_1*v*v*cSqb+Qtilde_2*v*v*sSqb+Qtilde_s*s*s));
  double m1Sq = lambda*Alambda*s*tb/sqrt(2.0)-0.5*lambda*lambda*(s*s+v*v*sSqb)-gbar*gbar*v*v*c2b/8.0
    -0.5*Qtilde_1*gdash_1*gdash_1*(Qtilde_1*v*v*cSqb+Qtilde_2*v*v*sSqb+Qtilde_s*s*s);
  double m2Sq = lambda*Alambda*s/(sqrt(2.0)*tb)-0.5*lambda*lambda*(s*s+v*v*cSqb)+gbar*gbar*v*v*c2b/8.0
    -0.5*Qtilde_2*gdash_1*gdash_1*(Qtilde_1*v*v*cSqb+Qtilde_2*v*v*sSqb+Qtilde_s*s*s);

  // Add in appropriate loop corrections (note sign!).
  msSq = msSq + doCalcTadpolesESSMS_Roman_atQ(essmSusy, s, tb);
  m1Sq = m1Sq + doCalcTadpoleESSMH1_Roman_atQ(essmSusy, s, tb);
  m2Sq = m2Sq + doCalcTadpoleESSMH2_Roman_atQ(essmSusy, s, tb);

  // Update the three masses in the object, doesn't change any other parameters.
  mH1Sq_vec.set(3, m1Sq);
  mH2Sq_vec.set(3, m2Sq);
  mSSq_vec.set(3, msSq);

  essmSusy.setMh1Squared(mH1Sq_vec);
  essmSusy.setMh2Squared(mH2Sq_vec);
  essmSusy.setmSsq(mSSq_vec);
}

// FineTunings_TreeLevel calculates the fine-tuning sensitivities \Delta_a for each 
// a = \lambda, A_\lambda, m_1^2, m_2^2, m_3^2, m_Q^2, m_U^2, A_t at tree-level,
// that is, neglecting the tadpole corrections to the EWSB conditions,
// at the given point in parameter space. 
// Note that at tree-level m_Q^2, m_U^2 and A_t do not appear, so the
// associated sensitivities vanish. 
// Inputs:
//     SoftParsEssm const & essmSusy = the object to evaluate the fine tuning at
//     double s = the value of the VEV s
//     double tb = the value of tan(beta) to use
//     DoubleVector & deltas = a vector (of length 8) to store the calculated
//     sensitivities, in the order \lambda, A_\lambda, m_1^2, m_2^2, m_s^2,
//     m_Q^2, m_U^2, A_t.
//     int & sing = integer flag to indicate if QR solution is singular
void FineTunings_TreeLevel(SoftParsEssm const & essmSusy, double s, double tb, DoubleVector & deltas, int & sing, double & detF)
{
  DoubleVector Alambda_vec = essmSusy.displayA_lambda();
  double Alambda = Alambda_vec.display(3);
  DoubleVector lambda_vec = essmSusy.displaylambda();
  double lambda = lambda_vec.display(3);

  DoubleVector mH1Sq_vec = essmSusy.displayMh1Squared();
  double m1Sq = mH1Sq_vec.display(3);
  DoubleVector mH2Sq_vec = essmSusy.displayMh2Squared();
  double m2Sq = mH2Sq_vec.display(3);
  DoubleVector mSSq_vec = essmSusy.displaymSsq();
  double msSq = mSSq_vec.display(3);

  double gdash = sqrt(3.0/5.0)*essmSusy.displayGaugeCoupling(1);
  double g_2 = essmSusy.displayGaugeCoupling(2);
  double gbar = sqrt(gdash*gdash+g_2*g_2);

  double v = essmSusy.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;
  double Mz = gbar*v/2.0;

  // Construct matrix of coefficients. This is the same for all
  // 8 tunings (only the right hand side changes),
  // so we can calculate the QR decomposition once and then store
  // it. Rows 1, 2, 3 correspond to EWSB conditions 1, 2, 3,
  // and columns 1, 2, 3 correspond to dM_Z/da, ds2b/d_a, ds/da.
  DoubleMatrix F_mat = TuningsCoefficientsMatrix_TreeLevel(essmSusy, s, tb);

  // Calculate determinant of LHS matrix
  detF = F_mat(1,1)*(F_mat(2,2)*F_mat(3,3)-F_mat(2,3)*F_mat(3,2))
    - F_mat(1,2)*(F_mat(2,1)*F_mat(3,3)-F_mat(2,3)*F_mat(3,1))
    + F_mat(1,3)*(F_mat(2,1)*F_mat(3,2)-F_mat(3,1)*F_mat(2,2));

  //  cout << "Tree level tunings matrix: " << endl;
  //cout << F_mat;

  // Get the QR decomposition of the left hand side.
  DoubleMatrix Q_mat(3,3);
  QR_decomp(F_mat, Q_mat, 3, sing);

  // For each input parameter, calculate the solution to the linear system F(dO/da) = b. 
  // The result dOda(1) is dM_Z/d_a for each a. Then store the corresponding sensitivity.
  DoubleVector dOda_vec(3), b_vec(3);
  double tuning;

  // NOTE changed to calculate everything in terms of derivatives of VEVs
  // instead of derivative of M_Z, hence extra step in calculating the sensitivity.

  // Also changed to calculate tunings wrt M_Z^2 instead of M_Z, hence the
  // extra factor of 2.

  // a = \lambda
  b_vec = TuningsRHS_Lambda_TreeLevel(essmSusy, s, tb);

  QR_solve(Q_mat, F_mat, 3, b_vec, dOda_vec, sing);

  //  cout << "dv_1/dlambda (tree level) = " << dOda_vec.display(1) << endl;
  //  cout << "dv_2/dlambda (tree level) = " << dOda_vec.display(2) << endl;
  //  cout << "ds/dlambda (tree level) = " << dOda_vec.display(3) << endl;

  tuning = fabs(2.0*(lambda/Mz)*(gbar/(2.0*v))*(v1*dOda_vec.display(1)+v2*dOda_vec.display(2)));//fabs(lambda*dOda_vec.display(1)/Mz);

  deltas.set(1, tuning);

  // a = A_\lambda
  b_vec = TuningsRHS_Alambda_TreeLevel(essmSusy, s, tb);

  QR_solve(Q_mat, F_mat, 3, b_vec, dOda_vec, sing);

  //  cout << "dv_1/dA_lambda (tree level) = " << dOda_vec.display(1) << endl;
  //  cout << "dv_2/dA_lambda (tree level) = " << dOda_vec.display(2) << endl;
  //  cout << "ds/dA_lambda (tree level) = " << dOda_vec.display(3) << endl;

  tuning = fabs(2.0*(Alambda/Mz)*(gbar/(2.0*v))*(v1*dOda_vec.display(1)+v2*dOda_vec.display(2)));//fabs(Alambda*dOda_vec.display(1)/Mz);

  deltas.set(2, tuning);

  // a = m_1^2
  b_vec = TuningsRHS_Mh1Sq_TreeLevel(essmSusy, s, tb);

  QR_solve(Q_mat, F_mat, 3, b_vec, dOda_vec, sing);

  tuning = fabs(2.0*(m1Sq/Mz)*(gbar/(2.0*v))*(v1*dOda_vec.display(1)+v2*dOda_vec.display(2)));//fabs(m1Sq*dOda_vec.display(1)/Mz);

  deltas.set(3, tuning);

  // a = m_2^2
  b_vec = TuningsRHS_Mh2Sq_TreeLevel(essmSusy, s, tb);

  QR_solve(Q_mat, F_mat, 3, b_vec, dOda_vec, sing);

  tuning = fabs(2.0*(m2Sq/Mz)*(gbar/(2.0*v))*(v1*dOda_vec.display(1)+v2*dOda_vec.display(2)));//fabs(m2Sq*dOda_vec.display(1)/Mz);

  deltas.set(4, tuning);

  // a = m_s^2
  b_vec = TuningsRHS_MsSq_TreeLevel(essmSusy, s, tb);

  QR_solve(Q_mat, F_mat, 3, b_vec, dOda_vec, sing);

  tuning = fabs(2.0*(msSq/Mz)*(gbar/(2.0*v))*(v1*dOda_vec.display(1)+v2*dOda_vec.display(2)));//fabs(msSq*dOda_vec.display(1)/Mz);

  deltas.set(5, tuning);

  // a = m_Q^2 (\Delta_m_Q^2 = 0 at tree level)
  deltas.set(6, 0.0);

  // a = m_U^2 (\Delta_m_U^2 = 0 at tree level)
  deltas.set(7,0.0);

  // a = A_t (\Delta_A_t = 0 at tree level)
  deltas.set(8, 0.0);

  return;
} 

// As above, but use Cramer's rule to obtain the solution. Note returns the determinant
// of the tunings coefficients matrix.
void FineTunings_TreeLevel_CR(SoftParsEssm const & essmSusy, double s, double tb, DoubleVector & deltas, int & sing, double & detF)
{
  DoubleVector Alambda_vec = essmSusy.displayA_lambda();
  double Alambda = Alambda_vec.display(3);
  DoubleVector lambda_vec = essmSusy.displaylambda();
  double lambda = lambda_vec.display(3);

  DoubleVector mH1Sq_vec = essmSusy.displayMh1Squared();
  double m1Sq = mH1Sq_vec.display(3);
  DoubleVector mH2Sq_vec = essmSusy.displayMh2Squared();
  double m2Sq = mH2Sq_vec.display(3);
  DoubleVector mSSq_vec = essmSusy.displaymSsq();
  double msSq = mSSq_vec.display(3);

  double gdash = sqrt(3.0/5.0)*essmSusy.displayGaugeCoupling(1);
  double g_2 = essmSusy.displayGaugeCoupling(2);
  double gbar = sqrt(gdash*gdash+g_2*g_2);

  double v = essmSusy.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;
  double Mz = gbar*v/2.0;

  // Construct matrix of coefficients
  DoubleMatrix F_mat = TuningsCoefficientsMatrix_TreeLevel(essmSusy, s, tb);

  // Compute sensitivities for each parameter
  DoubleVector dOda_vec(3), b_vec(3);
  double tuning;

  b_vec = TuningsRHS_Lambda_TreeLevel(essmSusy, s, tb);
  CramersRule(F_mat, b_vec, detF, sing, dOda_vec); 
  tuning = fabs(2.0*(lambda/Mz)*(gbar/(2.0*v))*(v1*dOda_vec.display(1)+v2*dOda_vec.display(2)));
  deltas.set(1, tuning);

  b_vec = TuningsRHS_Alambda_TreeLevel(essmSusy, s, tb);
  CramersRule(F_mat, b_vec, detF, sing, dOda_vec); 
  tuning = fabs(2.0*(Alambda/Mz)*(gbar/(2.0*v))*(v1*dOda_vec.display(1)+v2*dOda_vec.display(2)));
  deltas.set(2, tuning);

  b_vec = TuningsRHS_Mh1Sq_TreeLevel(essmSusy, s, tb);
  CramersRule(F_mat, b_vec, detF, sing, dOda_vec); 
  tuning = fabs(2.0*(m1Sq/Mz)*(gbar/(2.0*v))*(v1*dOda_vec.display(1)+v2*dOda_vec.display(2)));
  deltas.set(3, tuning);

  b_vec = TuningsRHS_Mh2Sq_TreeLevel(essmSusy, s, tb);
  CramersRule(F_mat, b_vec, detF, sing, dOda_vec); 
  tuning = fabs(2.0*(m2Sq/Mz)*(gbar/(2.0*v))*(v1*dOda_vec.display(1)+v2*dOda_vec.display(2)));
  deltas.set(4, tuning);

  b_vec = TuningsRHS_MsSq_TreeLevel(essmSusy, s, tb);
  CramersRule(F_mat, b_vec, detF, sing, dOda_vec); 
  tuning = fabs(2.0*(msSq/Mz)*(gbar/(2.0*v))*(v1*dOda_vec.display(1)+v2*dOda_vec.display(2)));
  deltas.set(5, tuning);

  deltas.set(6, 0.0);

  deltas.set(7,0.0);

  deltas.set(8, 0.0);

  return;
}


// A method to generate the matrix of coefficients of the 
// derivatives necessary to evaluate the tuning measures.
// Note changed to be in terms of derivatives wrt VEVs, 
// instead of wrt M_Z, sin (2\beta) and s.
DoubleMatrix TuningsCoefficientsMatrix_TreeLevel(SoftParsEssm const & essmSusy, double s, double tb)
{
  DoubleVector Alambda_vec = essmSusy.displayA_lambda();
  double Alambda = Alambda_vec.display(3);
  DoubleVector lambda_vec = essmSusy.displaylambda();
  double lambda = lambda_vec.display(3);

  double gdash = sqrt(3.0/5.0)*essmSusy.displayGaugeCoupling(1);
  double g_2 = essmSusy.displayGaugeCoupling(2);
  double gbar = sqrt(gdash*gdash+g_2*g_2);
  double gdash_1 = essmSusy.displaygdash_1();

  double Qtilde_1 = essmSusy.displayQH1tilde();
  double Qtilde_2 = essmSusy.displayQH2tilde();
  double Qtilde_s = essmSusy.displayQStilde();

  double v = essmSusy.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;
  double Mz = gbar*v/2.0;

  double cb = 1.0/sqrt(1.0+tb*tb);
  double sb = tb/sqrt(1.0+tb*tb);
  double s2b = 2.0*tb/(1.0+tb*tb);
  double c2b = (1.0-tb*tb)/(1.0+tb*tb);
  double t2b = s2b/c2b;

  // Construct matrix of coefficients. This is the same for all
  // 6 tunings (only the right hand side changes),
  // so we can calculate the QR decomposition once and then store
  // it. Rows 1, 2, 3 correspond to EWSB conditions 1, 2, 3,
  // and columns 1, 2, 3 correspond to dM_Z/da, ds2b/d_a, ds/da.
  DoubleMatrix F_mat(3,3);

  double F11 = lambda*Alambda*s*v2/(sqrt(2.0)*v1*v1)+gbar*gbar*v1/4.0+gdash_1*gdash_1*Qtilde_1*Qtilde_1*v1;//Mz*(c2b+(4.0/(gbar*gbar))*(lambda*lambda*sb*sb+gdash_1*gdash_1*Qtilde_1*
    //			  (Qtilde_1*cb*cb+Qtilde_2*sb*sb)));
  double F12 = lambda*lambda*v2-lambda*Alambda*s/(sqrt(2.0)*v1)-gbar*gbar*v2/4.0+gdash_1*gdash_1*Qtilde_1*Qtilde_2*v2;//lambda*lambda*Mz*Mz*t2b/(gbar*gbar)-lambda*Alambda*s/(2.0*sqrt(2.0)*cb*cb*c2b)
    //-Mz*Mz*t2b/2.0+gdash_1*gdash_1*Qtilde_1*Mz*Mz*(Qtilde_2-Qtilde_1)*t2b/(gbar*gbar);
  double F13 = lambda*lambda*s-lambda*Alambda*v2/(sqrt(2.0)*v1)+gdash_1*gdash_1*Qtilde_1*Qtilde_s*s;//lambda*lambda*s-lambda*Alambda*tb/sqrt(2.0)+gdash_1*gdash_1*Qtilde_1*Qtilde_s*s;

  double F21 = lambda*lambda*v1-lambda*Alambda*s/(sqrt(2.0)*v2)-gbar*gbar*v1/4.0+gdash_1*gdash_1*Qtilde_2*Qtilde_1*v1;//Mz*(-c2b+(4.0/(gbar*gbar))*(lambda*lambda*cb*cb+gdash_1*gdash_1*Qtilde_2*
    //			   (Qtilde_1*cb*cb+Qtilde_2*sb*sb)));
  double F22 = lambda*Alambda*s*v1/(sqrt(2.0)*v2*v2)+gbar*gbar*v2/4.0+gdash_1*gdash_1*Qtilde_2*Qtilde_2*v2;//-lambda*lambda*Mz*Mz*t2b/(gbar*gbar)+lambda*Alambda*s/(2.0*sqrt(2.0)*sb*sb*c2b)
    //+Mz*Mz*t2b/2.0+gdash_1*gdash_1*Qtilde_2*Mz*Mz*(Qtilde_2-Qtilde_1)*t2b/(gbar*gbar);
  double F23 = lambda*lambda*s-lambda*Alambda*v1/(sqrt(2.0)*v2)+gdash_1*gdash_1*Qtilde_2*Qtilde_s*s;//lambda*lambda*s-lambda*Alambda/(sqrt(2.0)*tb)+gdash_1*gdash_1*Qtilde_2*Qtilde_s*s;

  double F31 = lambda*lambda*v1-lambda*Alambda*v2/(sqrt(2.0)*s)+gdash_1*gdash_1*Qtilde_s*Qtilde_1*v1;//(4.0*Mz/(gbar*gbar))*(lambda*lambda-lambda*Alambda*s2b/(sqrt(2.0)*s)
				     //	+gdash_1*gdash_1*Qtilde_s*(Qtilde_1*cb*cb+Qtilde_2*sb*sb));
  double F32 = lambda*lambda*v2-lambda*Alambda*v1/(sqrt(2.0)*s)+gdash_1*gdash_1*Qtilde_s*Qtilde_2*v2;//(Mz*Mz/(gbar*gbar))*(gdash_1*gdash_1*Qtilde_s*(Qtilde_2-Qtilde_1)*t2b
				    //-sqrt(2.0)*lambda*Alambda/s);
  double F33 = lambda*Alambda*v1*v2/(sqrt(2.0)*s*s)+gdash_1*gdash_1*Qtilde_s*Qtilde_s*s;//sqrt(2.0)*lambda*Alambda*Mz*Mz*s2b/(s*s*gbar*gbar)+gdash_1*gdash_1*Qtilde_s*Qtilde_s*s;

  F_mat(1,1) = F11;
  F_mat(1,2) = F12;
  F_mat(1,3) = F13;
  F_mat(2,1) = F21;
  F_mat(2,2) = F22;
  F_mat(2,3) = F23;
  F_mat(3,1) = F31;
  F_mat(3,2) = F32;
  F_mat(3,3) = F33;

  return F_mat;

}

DoubleVector TuningsRHS_Lambda_TreeLevel(SoftParsEssm const & essmSusy, double s, double tb)
{
  DoubleVector Alambda_vec = essmSusy.displayA_lambda();
  double Alambda = Alambda_vec.display(3);
  DoubleVector lambda_vec = essmSusy.displaylambda();
  double lambda = lambda_vec.display(3);

  double gdash = sqrt(3.0/5.0)*essmSusy.displayGaugeCoupling(1);
  double g_2 = essmSusy.displayGaugeCoupling(2);
  double gbar = sqrt(gdash*gdash+g_2*g_2);
  double gdash_1 = essmSusy.displaygdash_1();

  double Qtilde_1 = essmSusy.displayQH1tilde();
  double Qtilde_2 = essmSusy.displayQH2tilde();
  double Qtilde_s = essmSusy.displayQStilde();

  double v = essmSusy.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;
  double Mz = gbar*v/2.0;

  double cb = 1.0/sqrt(1.0+tb*tb);
  double sb = tb/sqrt(1.0+tb*tb);
  double s2b = 2.0*tb/(1.0+tb*tb);
  double c2b = (1.0-tb*tb)/(1.0+tb*tb);
  double t2b = s2b/c2b;

  DoubleVector b_vec(3);
  double b1, b2, b3;

  b1 = Alambda*s*v2/(sqrt(2.0)*v1)-lambda*(v2*v2+s*s);//Alambda*s*tb/sqrt(2.0)-lambda*(4.0*Mz*Mz*sb*sb/(gbar*gbar)+s*s);
  b2 = Alambda*s*v1/(sqrt(2.0)*v2)-lambda*(v1*v1+s*s);//Alambda*s/(sqrt(2.0)*tb)-lambda*(4.0*Mz*Mz*cb*cb/(gbar*gbar)+s*s);
  b3 = Alambda*v1*v2/(sqrt(2.0)*s)-lambda*v*v;//(2.0*Mz*Mz/(gbar*gbar))*(Alambda*s2b/(sqrt(2.0)*s)-2.0*lambda);

  b_vec.set(1, b1);
  b_vec.set(2, b2);
  b_vec.set(3, b3);

  return b_vec;
}

DoubleVector TuningsRHS_Alambda_TreeLevel(SoftParsEssm const & essmSusy, double s, double tb)
{
  DoubleVector Alambda_vec = essmSusy.displayA_lambda();
  double Alambda = Alambda_vec.display(3);
  DoubleVector lambda_vec = essmSusy.displaylambda();
  double lambda = lambda_vec.display(3);

  double gdash = sqrt(3.0/5.0)*essmSusy.displayGaugeCoupling(1);
  double g_2 = essmSusy.displayGaugeCoupling(2);
  double gbar = sqrt(gdash*gdash+g_2*g_2);
  double gdash_1 = essmSusy.displaygdash_1();

  double Qtilde_1 = essmSusy.displayQH1tilde();
  double Qtilde_2 = essmSusy.displayQH2tilde();
  double Qtilde_s = essmSusy.displayQStilde();

  double v = essmSusy.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;
  double Mz = gbar*v/2.0;

  double cb = 1.0/sqrt(1.0+tb*tb);
  double sb = tb/sqrt(1.0+tb*tb);
  double s2b = 2.0*tb/(1.0+tb*tb);
  double c2b = (1.0-tb*tb)/(1.0+tb*tb);
  double t2b = s2b/c2b;

  DoubleVector b_vec(3);
  double b1, b2, b3;

  b1 = lambda*s*v2/(sqrt(2.0)*v1);//lambda*s*tb/sqrt(2.0);//
  b2 = lambda*s*v1/(sqrt(2.0)*v2);//lambda*s/(sqrt(2.0)*tb);//
  b3 = lambda*v1*v2/(sqrt(2.0)*s);//sqrt(2.0)*lambda*Mz*Mz*s2b/(s*gbar*gbar);//

  b_vec.set(1, b1);
  b_vec.set(2, b2);
  b_vec.set(3, b3);

  return b_vec;
}

DoubleVector TuningsRHS_Mh1Sq_TreeLevel(SoftParsEssm const & essmSusy, double s, double tb)
{

  DoubleVector b_vec(3);
  double b1, b2, b3;

  b1 = -1.0;
  b2 = 0.0;
  b3 = 0.0;

  b_vec.set(1, b1);
  b_vec.set(2, b2);
  b_vec.set(3, b3);

  return b_vec;
}

DoubleVector TuningsRHS_Mh2Sq_TreeLevel(SoftParsEssm const & essmSusy, double s, double tb)
{

  DoubleVector b_vec(3);
  double b1, b2, b3;

  b1 = 0.0;
  b2 = -1.0;
  b3 = 0.0;

  b_vec.set(1, b1);
  b_vec.set(2, b2);
  b_vec.set(3, b3);

  return b_vec;
}

DoubleVector TuningsRHS_MsSq_TreeLevel(SoftParsEssm const & essmSusy, double s, double tb)
{

  DoubleVector b_vec(3);
  double b1, b2, b3;

  b1 = 0.0;
  b2 = 0.0;
  b3 = -1.0;

  b_vec.set(1, b1);
  b_vec.set(2, b2);
  b_vec.set(3, b3);

  return b_vec;
}

DoubleVector TuningsRHS_MQlSq_TreeLevel(SoftParsEssm const & essmSusy, double s, double tb)
{
  DoubleVector b_vec(3);
  b_vec.set(1, 0.0);
  b_vec.set(2, 0.0);
  b_vec.set(3, 0.0);
  return b_vec;
}

DoubleVector TuningsRHS_MUrSq_TreeLevel(SoftParsEssm const & essmSusy, double s, double tb)
{
  DoubleVector b_vec(3);
  b_vec.set(1, 0.0);
  b_vec.set(2, 0.0);
  b_vec.set(3, 0.0);
  return b_vec;
}

DoubleVector TuningsRHS_At_TreeLevel(SoftParsEssm const & essmSusy, double s, double tb)
{
  DoubleVector b_vec(3);
  b_vec.set(1, 0.0);
  b_vec.set(2, 0.0);
  b_vec.set(3, 0.0);
  return b_vec;
}

// As above, but including the leading one-loop contributions from
// top and stop loops.
void FineTunings_OneLoop(SoftParsEssm const & essmSusy, double s, double tb, DoubleVector & deltas, int & sing, double & detF)
{
  DoubleVector Alambda_vec = essmSusy.displayA_lambda();
  double Alambda = Alambda_vec.display(3);
  DoubleVector lambda_vec = essmSusy.displaylambda();
  double lambda = lambda_vec.display(3);

  DoubleVector mH1Sq_vec = essmSusy.displayMh1Squared();
  double m1Sq = mH1Sq_vec.display(3);
  DoubleVector mH2Sq_vec = essmSusy.displayMh2Squared();
  double m2Sq = mH2Sq_vec.display(3);
  DoubleVector mSSq_vec = essmSusy.displaymSsq();
  double msSq = mSSq_vec.display(3);
  double At = essmSusy.displayA_top();
  double mQlsq = essmSusy.displaySoftMassSquared(mQl,3,3);
  double mUrsq =  essmSusy.displaySoftMassSquared(mUr,3,3);


  double gdash = sqrt(3.0/5.0)*essmSusy.displayGaugeCoupling(1);
  double g_2 = essmSusy.displayGaugeCoupling(2);
  double gbar = sqrt(gdash*gdash+g_2*g_2);

  double v = essmSusy.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;
  double Mz = gbar*v/2.0;

  // Construct matrix of coefficients. This is the same for all
  // 8 tunings (only the right hand side changes),
  // so we can calculate the QR decomposition once and then store
  // it. Rows 1, 2, 3 correspond to EWSB conditions 1, 2, 3,
  // and columns 1, 2, 3 correspond to dM_Z/da, ds2b/d_a, ds/da.
  DoubleMatrix F_mat = TuningsCoefficientsMatrix_OneLoop(essmSusy, s, tb);

  
  // Calculate determinant of LHS matrix
  detF = F_mat(1,1)*(F_mat(2,2)*F_mat(3,3)-F_mat(2,3)*F_mat(3,2))
    - F_mat(1,2)*(F_mat(2,1)*F_mat(3,3)-F_mat(2,3)*F_mat(3,1))
    + F_mat(1,3)*(F_mat(2,1)*F_mat(3,2)-F_mat(3,1)*F_mat(2,2));

  //cout << "One loop tunings matrix: " << endl;
  //cout << F_mat;

  // Get the QR decomposition of the left hand side.
  DoubleMatrix Q_mat(3,3);
  QR_decomp(F_mat, Q_mat, 3, sing);

  // For each input parameter, calculate the solution to the linear system F(dO/da) = b. 
  // The result dOda(1) is dM_Z/d_a for each a. Then store the corresponding sensitivity.
  DoubleVector dOda_vec(3), b_vec(3);
  double tuning;

  // a = \lambda
  b_vec = TuningsRHS_Lambda_OneLoop(essmSusy, s, tb);



  QR_solve(Q_mat, F_mat, 3, b_vec, dOda_vec, sing);

  //  cout << "dv_1/dlambda (one loop) = " << dOda_vec.display(1) << endl;
  //  cout << "dv_2/dlambda (one loop) = " << dOda_vec.display(2) << endl;
  //  cout << "ds/dlambda (one loop) = " << dOda_vec.display(3) << endl;


  tuning = fabs(2.0*(lambda/Mz)*(gbar/(2.0*v))*(v1*dOda_vec.display(1)+v2*dOda_vec.display(2)));//fabs(lambda*dOda_vec.display(1)/Mz);

  deltas.set(1, tuning);

  // a = A_\lambda
  b_vec = TuningsRHS_Alambda_OneLoop(essmSusy, s, tb);

  QR_solve(Q_mat, F_mat, 3, b_vec, dOda_vec, sing);

  //  cout << "dv_1/dA_lambda (one loop) = " << dOda_vec.display(1) << endl;
  //  cout << "dv_2/dA_lambda (one loop) = " << dOda_vec.display(2) << endl;
  //  cout << "ds/dA_lambda (one loop) = " << dOda_vec.display(3) << endl;


  tuning = fabs(2.0*(Alambda/Mz)*(gbar/(2.0*v))*(v1*dOda_vec.display(1)+v2*dOda_vec.display(2)));//fabs(Alambda*dOda_vec.display(1)/Mz);

  deltas.set(2, tuning);

  // a = m_1^2
  b_vec = TuningsRHS_Mh1Sq_OneLoop(essmSusy, s, tb);
  QR_solve(Q_mat, F_mat, 3, b_vec, dOda_vec, sing);

  tuning = fabs(2.0*(m1Sq/Mz)*(gbar/(2.0*v))*(v1*dOda_vec.display(1)+v2*dOda_vec.display(2)));//fabs(m1Sq*dOda_vec.display(1)/Mz);

  deltas.set(3, tuning);

  // a = m_2^2
  b_vec = TuningsRHS_Mh2Sq_OneLoop(essmSusy, s, tb);

  QR_solve(Q_mat, F_mat, 3, b_vec, dOda_vec, sing);

  tuning = fabs(2.0*(m2Sq/Mz)*(gbar/(2.0*v))*(v1*dOda_vec.display(1)+v2*dOda_vec.display(2)));//fabs(m2Sq*dOda_vec.display(1)/Mz);

  deltas.set(4, tuning);

  // a = m_s^2
  b_vec = TuningsRHS_MsSq_OneLoop(essmSusy, s, tb);

  QR_solve(Q_mat, F_mat, 3, b_vec, dOda_vec, sing);

  tuning = fabs(2.0*(msSq/Mz)*(gbar/(2.0*v))*(v1*dOda_vec.display(1)+v2*dOda_vec.display(2)));//fabs(msSq*dOda_vec.display(1)/Mz);

  deltas.set(5, tuning);

  // a = m_Q^2 (\Delta_m_Q^2 = 0 at tree level)
  b_vec = TuningsRHS_MQlSq_OneLoop(essmSusy, s, tb);

  QR_solve(Q_mat, F_mat, 3, b_vec, dOda_vec, sing);

  tuning = fabs(2.0*(mQlsq/Mz)*(gbar/(2.0*v))*(v1*dOda_vec.display(1)+v2*dOda_vec.display(2)));//fabs(mQlsq*dOda_vec.display(1)/Mz);

  deltas.set(6, tuning);

  // a = m_U^2 (\Delta_m_U^2 = 0 at tree level)
  b_vec = TuningsRHS_MUrSq_OneLoop(essmSusy, s, tb);

  QR_solve(Q_mat, F_mat, 3, b_vec, dOda_vec, sing);

  tuning = fabs(2.0*(mUrsq/Mz)*(gbar/(2.0*v))*(v1*dOda_vec.display(1)+v2*dOda_vec.display(2)));//fabs(mUrsq*dOda_vec.display(1)/Mz);

  deltas.set(7, tuning);

  // a = A_t (\Delta_A_t = 0 at tree level)
  b_vec = TuningsRHS_At_OneLoop(essmSusy, s, tb);

  QR_solve(Q_mat, F_mat, 3, b_vec, dOda_vec, sing);

  tuning = fabs(2.0*(At/Mz)*(gbar/(2.0*v))*(v1*dOda_vec.display(1)+v2*dOda_vec.display(2)));//fabs(At*dOda_vec.display(1)/Mz);

  deltas.set(8, tuning);

  return;
}

void FineTunings_OneLoop_CR(SoftParsEssm const & essmSusy, double s, double tb, DoubleVector & deltas, int & sing, double & detF)
{
  DoubleVector Alambda_vec = essmSusy.displayA_lambda();
  double Alambda = Alambda_vec.display(3);
  DoubleVector lambda_vec = essmSusy.displaylambda();
  double lambda = lambda_vec.display(3);

  DoubleVector mH1Sq_vec = essmSusy.displayMh1Squared();
  double m1Sq = mH1Sq_vec.display(3);
  DoubleVector mH2Sq_vec = essmSusy.displayMh2Squared();
  double m2Sq = mH2Sq_vec.display(3);
  DoubleVector mSSq_vec = essmSusy.displaymSsq();
  double msSq = mSSq_vec.display(3);
  double At = essmSusy.displayA_top();
  double mQlsq = essmSusy.displaySoftMassSquared(mQl,3,3);
  double mUrsq =  essmSusy.displaySoftMassSquared(mUr,3,3);


  double gdash = sqrt(3.0/5.0)*essmSusy.displayGaugeCoupling(1);
  double g_2 = essmSusy.displayGaugeCoupling(2);
  double gbar = sqrt(gdash*gdash+g_2*g_2);

  double v = essmSusy.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;
  double Mz = gbar*v/2.0;

  // Construct matrix of coefficients
  DoubleMatrix F_mat = TuningsCoefficientsMatrix_OneLoop(essmSusy, s, tb);

  // Compute sensitivities for each parameter
  DoubleVector dOda_vec(3), b_vec(3);
  double tuning;

  b_vec = TuningsRHS_Lambda_OneLoop(essmSusy, s, tb);
  CramersRule(F_mat, b_vec, detF, sing, dOda_vec); 
  tuning = fabs(2.0*(lambda/Mz)*(gbar/(2.0*v))*(v1*dOda_vec.display(1)+v2*dOda_vec.display(2)));
  deltas.set(1, tuning);

  b_vec = TuningsRHS_Alambda_OneLoop(essmSusy, s, tb);
  CramersRule(F_mat, b_vec, detF, sing, dOda_vec); 
  tuning = fabs(2.0*(Alambda/Mz)*(gbar/(2.0*v))*(v1*dOda_vec.display(1)+v2*dOda_vec.display(2)));
  deltas.set(2, tuning);

  b_vec = TuningsRHS_Mh1Sq_OneLoop(essmSusy, s, tb);
  CramersRule(F_mat, b_vec, detF, sing, dOda_vec); 
  tuning = fabs(2.0*(m1Sq/Mz)*(gbar/(2.0*v))*(v1*dOda_vec.display(1)+v2*dOda_vec.display(2)));
  deltas.set(3, tuning);

  b_vec = TuningsRHS_Mh2Sq_OneLoop(essmSusy, s, tb);
  CramersRule(F_mat, b_vec, detF, sing, dOda_vec); 
  tuning = fabs(2.0*(m2Sq/Mz)*(gbar/(2.0*v))*(v1*dOda_vec.display(1)+v2*dOda_vec.display(2)));
  deltas.set(4, tuning);

  b_vec = TuningsRHS_MsSq_OneLoop(essmSusy, s, tb);
  CramersRule(F_mat, b_vec, detF, sing, dOda_vec); 
  tuning = fabs(2.0*(msSq/Mz)*(gbar/(2.0*v))*(v1*dOda_vec.display(1)+v2*dOda_vec.display(2)));
  deltas.set(5, tuning);

  b_vec = TuningsRHS_MQlSq_OneLoop(essmSusy, s, tb);
  CramersRule(F_mat, b_vec, detF, sing, dOda_vec); 
  tuning = fabs(2.0*(mQlsq/Mz)*(gbar/(2.0*v))*(v1*dOda_vec.display(1)+v2*dOda_vec.display(2)));
  deltas.set(6, tuning);

  b_vec = TuningsRHS_MUrSq_OneLoop(essmSusy, s, tb);
  CramersRule(F_mat, b_vec, detF, sing, dOda_vec); 
  tuning = fabs(2.0*(mUrsq/Mz)*(gbar/(2.0*v))*(v1*dOda_vec.display(1)+v2*dOda_vec.display(2)));
  deltas.set(7, tuning);

  b_vec = TuningsRHS_At_OneLoop(essmSusy, s, tb);
  CramersRule(F_mat, b_vec, detF, sing, dOda_vec); 
  tuning = fabs(2.0*(At/Mz)*(gbar/(2.0*v))*(v1*dOda_vec.display(1)+v2*dOda_vec.display(2)));
  deltas.set(8, tuning);

  return;
}

DoubleMatrix TuningsCoefficientsMatrix_OneLoop(SoftParsEssm const & essmSusy, double s, double tb)
{
  double v = essmSusy.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;
  double g1 = essmSusy.displayGaugeCoupling(1);
  double g2 = essmSusy.displayGaugeCoupling(2);
  double gbar = sqrt(g2*g2+0.6*g1*g1);
  double gd1 = essmSusy.displaygdash_1();

  double cb = 1.0/sqrt(1.0+tb*tb);
  double sb = tb*cb;
  double c2b = (1.0-tb*tb)/(1.0+tb*tb);

  // Get the tree level coefficients matrix
  DoubleMatrix F_mat = TuningsCoefficientsMatrix_TreeLevel(essmSusy, s, tb);

  // Calculate the tadpoles. Note that Peter's methods have the OPPOSITE
  // sign to the one we need.
  double Delta1tp = -doCalcTadpoleESSMH1(essmSusy, s, tb);
  double Delta2tp = -doCalcTadpoleESSMH2(essmSusy, s, tb);
  double Delta3tp = -doCalcTadpolesESSMS(essmSusy, s, tb);

  //cout << "Delta_1^tp = " << Delta1tp << endl;
  //cout << "Delta_2^tp = " << Delta2tp << endl;
  //cout << "Delta_3^tp = " << Delta3tp << endl;

  // Calculate the one-loop corrections to the Higgs mass matrix.
  double Delta11prime = doCalcDeltaPrime11(essmSusy, s, tb);
  double Delta12prime = doCalcDeltaPrime12(essmSusy, s, tb);
  double Delta13prime = doCalcDeltaPrime13(essmSusy, s, tb);
  double Delta22prime = doCalcDeltaPrime22(essmSusy, s, tb);
  double Delta23prime = doCalcDeltaPrime23(essmSusy, s, tb);
  double Delta33prime = doCalcDeltaPrime33(essmSusy, s, tb);

  //cout << "Delta'_11 = " << Delta11prime << endl;
  //cout << "Delta'_12 = " << Delta12prime << endl;
  //cout << "Delta'_13 = " << Delta13prime << endl;
  //cout << "Delta'_22 = " << Delta22prime << endl;
  //cout << "Delta'_23 = " << Delta23prime << endl;
  //cout << "Delta'_33 = " << Delta33prime << endl;

  // Rows 1, 2, 3 correspond to EWSB conditions 1, 2, 3,
  // and columns 1, 2, 3 correspond to dM_Z/da, ds2b/d_a, ds/da.
  // NOTE Changed to be in terms of v_1, v_2 and s instead of M_Z, sin(2\beta) and s.
  F_mat(1,1) = F_mat(1,1) + (1.0/v1)*(Delta11prime-Delta1tp);//F_mat(1,1) + (2.0*cb/gbar)*(Delta11prime/v1-Delta1tp/v1)+2.0*sb*Delta12prime/(gbar*v1);
  F_mat(1,2) = F_mat(1,2) + Delta12prime/v1;//F_mat(1,2) + v*cb*Delta12prime/(2.0*c2b*v1)-(v*sb/(2.0*c2b))*(Delta11prime/v1-Delta1tp/v1);
  F_mat(1,3) = F_mat(1,3) + Delta13prime/v1;//F_mat(1,3) + Delta13prime/v1;

  F_mat(2,1) = F_mat(2,1) + Delta12prime/v2;//F_mat(2,1) + 2.0*cb*Delta12prime/(gbar*v2)+(2.0*sb/gbar)*(Delta22prime/v2-Delta2tp/v2);
  F_mat(2,2) = F_mat(2,2) + (1.0/v2)*(Delta22prime-Delta2tp);//F_mat(2,2) + (v*cb/(2.0*c2b))*(Delta22prime/v2-Delta2tp/v2)-v*sb*Delta12prime/(2.0*c2b*v2);
  F_mat(2,3) = F_mat(2,3) + Delta23prime/v2;//F_mat(2,3) + Delta23prime/v2;

  F_mat(3,1) = F_mat(3,1) + Delta13prime/s;//F_mat(3,1) + 2.0*cb*Delta13prime/(gbar*s)+2.0*sb*Delta23prime/(gbar*s);
  F_mat(3,2) = F_mat(3,2) + Delta23prime/s;//F_mat(3,2) + v*cb*Delta23prime/(2.0*c2b*s)-v*sb*Delta13prime/(2.0*c2b*s);
  F_mat(3,3) = F_mat(3,3) + (1.0/s)*(Delta33prime-Delta3tp);//F_mat(3,3) + Delta33prime/s-Delta3tp/s;

  return F_mat;

}

// Note that the only input parameters that receive corrections from the top and stop
// loops are \lambda, A_t, m_{Q_l}^2 and m_{U_r}^2.
DoubleVector TuningsRHS_Lambda_OneLoop(SoftParsEssm const & essmSusy, double s, double tb)
{
  // Get the tree level rhs
  DoubleVector b_vec(3);

  b_vec = TuningsRHS_Lambda_TreeLevel(essmSusy, s, tb);

  // Calculate the appropriate loop corrections (second derivatives of
  // \Delta V)
  double d2DeltaVdldv1 = doCalcd2DeltaVdLambdadv1(essmSusy, s, tb);
  double d2DeltaVdldv2 = doCalcd2DeltaVdLambdadv2(essmSusy, s, tb);
  double d2DeltaVdldv3 = doCalcd2DeltaVdLambdadv3(essmSusy, s, tb);

  // Add (subtract with our sign convention) the contributions to the 
  // rhs
  b_vec(1) = b_vec(1)-d2DeltaVdldv1;
  b_vec(2) = b_vec(2)-d2DeltaVdldv2;
  b_vec(3) = b_vec(3)-d2DeltaVdldv3;

  return b_vec;
}

DoubleVector TuningsRHS_Alambda_OneLoop(SoftParsEssm const & essmSusy, double s, double tb)
{
  // Get the tree level rhs
  DoubleVector b_vec(3);

  b_vec = TuningsRHS_Alambda_TreeLevel(essmSusy, s, tb);

  // A_\lambda has no corrections from stop and top loops
  return b_vec;
}

DoubleVector TuningsRHS_Mh1Sq_OneLoop(SoftParsEssm const & essmSusy, double s, double tb)
{
  // Get the tree level rhs
  DoubleVector b_vec(3);

  b_vec = TuningsRHS_Mh1Sq_TreeLevel(essmSusy, s, tb);

  // m_H_1^2 has no corrections from stop and top loops
  return b_vec;
}

DoubleVector TuningsRHS_Mh2Sq_OneLoop(SoftParsEssm const & essmSusy, double s, double tb)
{
  // Get the tree level rhs
  DoubleVector b_vec(3);

  b_vec = TuningsRHS_Mh2Sq_TreeLevel(essmSusy, s, tb);

  // m_H_2^2 has no corrections from stop and top loops
  return b_vec;
}

DoubleVector TuningsRHS_MsSq_OneLoop(SoftParsEssm const & essmSusy, double s, double tb)
{
  // Get the tree level rhs
  DoubleVector b_vec(3);

  b_vec = TuningsRHS_MsSq_TreeLevel(essmSusy, s, tb);

  // m_s^2 has no corrections from stop and top loops
  return b_vec;
}

DoubleVector TuningsRHS_MQlSq_OneLoop(SoftParsEssm const & essmSusy, double s, double tb)
{
  // Get the tree level rhs
  DoubleVector b_vec(3);

  b_vec = TuningsRHS_MQlSq_TreeLevel(essmSusy, s, tb);
  // Calculate the appropriate loop corrections (second derivatives of
  // \Delta V)
  double d2DeltaVdmQlsqdv1 = doCalcd2DeltaVdmQlsqdv1(essmSusy, s, tb);
  double d2DeltaVdmQlsqdv2 = doCalcd2DeltaVdmQlsqdv2(essmSusy, s, tb);
  double d2DeltaVdmQlsqdv3 = doCalcd2DeltaVdmQlsqdv3(essmSusy, s, tb);

  // Add (subtract with our sign convention) the contributions to the 
  // rhs
  b_vec(1) = b_vec(1)-d2DeltaVdmQlsqdv1;
  b_vec(2) = b_vec(2)-d2DeltaVdmQlsqdv2;
  b_vec(3) = b_vec(3)-d2DeltaVdmQlsqdv3;

  return b_vec;
}

DoubleVector TuningsRHS_MUrSq_OneLoop(SoftParsEssm const & essmSusy, double s, double tb)
{
  // Get the tree level rhs
  DoubleVector b_vec(3);

  b_vec = TuningsRHS_MUrSq_TreeLevel(essmSusy, s, tb);
  // Calculate the appropriate loop corrections (second derivatives of
  // \Delta V)
  double d2DeltaVdmUrsqdv1 = doCalcd2DeltaVdmUrsqdv1(essmSusy, s, tb);
  double d2DeltaVdmUrsqdv2 = doCalcd2DeltaVdmUrsqdv2(essmSusy, s, tb);
  double d2DeltaVdmUrsqdv3 = doCalcd2DeltaVdmUrsqdv3(essmSusy, s, tb);

  // Add (subtract with our sign convention) the contributions to the 
  // rhs
  b_vec(1) = b_vec(1)-d2DeltaVdmUrsqdv1;
  b_vec(2) = b_vec(2)-d2DeltaVdmUrsqdv2;
  b_vec(3) = b_vec(3)-d2DeltaVdmUrsqdv3;

  return b_vec;
}

DoubleVector TuningsRHS_At_OneLoop(SoftParsEssm const & essmSusy, double s, double tb)
{
  // Get the tree level rhs
  DoubleVector b_vec(3);

  b_vec = TuningsRHS_At_TreeLevel(essmSusy, s, tb);
  // Calculate the appropriate loop corrections (second derivatives of
  // \Delta V)
  double d2DeltaVdAtdv1 = doCalcd2DeltaVdAtdv1(essmSusy, s, tb);
  double d2DeltaVdAtdv2 = doCalcd2DeltaVdAtdv2(essmSusy, s, tb);
  double d2DeltaVdAtdv3 = doCalcd2DeltaVdAtdv3(essmSusy, s, tb);

  // Add (subtract with our sign convention) the contributions to the 
  // rhs
  b_vec(1) = b_vec(1)-d2DeltaVdAtdv1;
  b_vec(2) = b_vec(2)-d2DeltaVdAtdv2;
  b_vec(3) = b_vec(3)-d2DeltaVdAtdv3;

  return b_vec;
}

// The following is a final fine tuning method that should incorporate the
// functionality of all of the previously defined methods above. It uses
// the version of the EWSB conditions in which the VEVs are not divided out
// to start with (i.e. unrealistic vacua are permitted).
// Inputs:
//     SoftParsEssm const & essmSusy = E6SSM model object, assumed to contain correct
//     values of parameters appearing in the Higgs scalar potential (e.g. lambda, 
//     A_lambda, and the gauge couplings).
//     double s = the value of s to use
//     double tb = the value of tan(beta) to use
//     int l = # loops to use (0 or 1)
//     DoubleMatrix const & partials = a Nx4 matrix containing the partial derivatives of the
//     EWSB conditions wrt the N input parameters. The element (i, 1) contains the value of
//     the parameter p_i. The (i, j)th, j>1, element is taken to be
//     df_(j-1)/dp_i. These should be calculated consistently with the order requested when
//     this function is called, i.e. tadpole corrections should be included if l = 1.
//     int & sing = flag to indicate if singular (non-zero if problem)
//     double & detF = determinant of CP-even Higgs mass matrix
//     
//     Returns a DoubleVector of length N containing the fine tuning sensitivities. The ith
//     element of the vector is the tuning corresponding to the parameter in the ith row
//     of the input matrix.
DoubleVector E6SSM_FineTunings(SoftParsEssm const & essmSusy, double s, double tb, int l,
			       DoubleMatrix const & partials, int & sing, double & detF)
{
  sing = 0;

  DoubleVector lambda_vec = essmSusy.displaylambda();
  double lambda = lambda_vec(3);
  DoubleVector Alambda_vec = essmSusy.displayA_lambda();
  double Alambda = Alambda_vec(3);
  DoubleVector mH1Sq_vec = essmSusy.displayMh1Squared();
  double mH1Sq = mH1Sq_vec(3);
  DoubleVector mH2Sq_vec = essmSusy.displayMh2Squared();
  double mH2Sq = mH2Sq_vec(3);
  DoubleVector msSq_vec = essmSusy.displaymSsq();
  double msSq = msSq_vec(3);

  double v = essmSusy.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;

  double cb = 1.0/sqrt(1.0+tb*tb);
  double sb = cb*tb;
  double s2b = 2.0*sb*cb;

  double g1 = essmSusy.displayGaugeCoupling(1);
  double g2 = essmSusy.displayGaugeCoupling(2);
  double gbar = sqrt(g2*g2+0.6*g1*g1);
  double gd1 = essmSusy.displaygdash_1();

  double Q1 = essmSusy.displayQH1tilde();
  double Q2 = essmSusy.displayQH2tilde();
  double Qs = essmSusy.displayQStilde();

  // Construct tree level contributions to mass matrix and RHS vector. Note I am
  // going to use the versions in which the EWSB conditions are used to eliminate
  // the soft squared masses mH1Sq, mH2Sq and mSSq.
  double F11 = mH1Sq+0.5*lambda*lambda*(v2*v2+s*s)-0.125*gbar*gbar*(v2*v2-3.0*v1*v1)
    +0.5*gd1*gd1*Q1*(3.0*Q1*v1*v1+Q2*v2*v2+Qs*s*s);
  double F22 = mH2Sq+0.5*lambda*lambda*(v1*v1+s*s)+0.125*gbar*gbar*(3.0*v2*v2-v1*v1)
    +0.5*gd1*gd1*Q2*(Q1*v1*v1+3.0*Q2*v2*v2+Qs*s*s);
  double F33 = msSq+0.5*lambda*lambda*v*v+0.5*gd1*gd1*Qs*(Q1*v1*v1+Q2*v2*v2+3.0*Qs*s*s);
  double F12 = lambda*lambda*v2*v1-lambda*Alambda*s/sqrt(2.0)-0.25*gbar*gbar*v1*v2+gd1*gd1*Q1*Q2*v1*v2;
  double F13 = lambda*lambda*s*v1-lambda*Alambda*v2/sqrt(2.0)+gd1*gd1*Q1*Qs*v1*s;
  double F23 = lambda*lambda*s*v2-lambda*Alambda*v1/sqrt(2.0)+gd1*gd1*Q2*Qs*s*v2;

  // If loops requested, add in loop corrections
  if (l != 0)
    {
      if (l != 1)
	{
	  cout << "Warning: only l = 0 or l = 1 loops currently supported." << endl;
	  cout << "Using l = 1 loops." << endl;
	}
      F11 = F11+doCalcDeltaPrime11(essmSusy, s, tb);
      F22 = F22+doCalcDeltaPrime22(essmSusy, s, tb);
      F33 = F33+doCalcDeltaPrime33(essmSusy, s, tb);
      F12 = F12+doCalcDeltaPrime12(essmSusy, s, tb);
      F13 = F13+doCalcDeltaPrime13(essmSusy, s, tb);
      F23 = F23+doCalcDeltaPrime23(essmSusy, s, tb);
    }


  // Compute coefficients appearing in the expression for the tuning sensitivity.
  double alpha = v1*(F23*F23-F22*F33)+v2*(F12*F33-F23*F13);
  double beta = v1*(F12*F33-F13*F23)+v2*(F13*F13-F11*F33);
  double gamma = v1*(F13*F22-F12*F23)+v2*(F11*F23-F13*F12);

  detF = F11*(F22*F33-F23*F23)-F12*(F12*F33-F23*F13)+F13*(F12*F23-F22*F13);

  int nRows = partials.displayRows();
  int nCols = partials.displayCols();

  DoubleVector tunings(nRows);

  // Check if sufficient information provided to compute the tuning
  // sensitivities. If not, throw an exception.
  if (nCols < 4)
    {
    ostringstream ii;
    ii << "Insufficient information provided to calculate fine tuning: needed 4 columns, got " << nCols << ".\n";
    throw(ii.str());
    }

  if (fabs(detF) < EPSTOL)
    {
      cout << "Warning: determinant of CP-even Higgs mass matrix is zero." << endl;
      cout << "Tunings cannot be trusted at this point." << endl;
      sing = 1; // As an additional flag all sensitivities are returned as being equal to -1. 
      for (int i = 1; i <= nRows; i++)
	{
	  tunings(i) = -1.0;
	}
      return tunings;
    }
  else // Otherwise det(F) is non-zero and there is a unique solution.
    {

      // Go ahead and compute the requested sensitivities.
      for (int i = 1; i <= nRows; i++)
	{
	  tunings(i) = (2.0*partials(i,1)/(v*v*detF))*(alpha*partials(i,2)+beta*partials(i,3)+gamma*partials(i,4));
	  tunings(i) = fabs(tunings(i));
	}
      return tunings;
    }
}

// This function produces the partials matrix for the specific parameter set 
// {lambda, A_lambda, m_H1^2, m_H2^2, m_S^2, m_Q^2, m_U^2, A_t} with NO RG RUNNING.
// Inputs:
//     SoftParsEssm const & essmSusy = E6SSM model object to calculate tuning for
//     double s = the value of s to use
//     double tb = the value of tan(beta) to use
//     int l = # loops to use (0 or 1)
//     DoubleMatrix & partials = a matrix to store the resulting derivatives
void calculatePartialsMatrix(SoftParsEssm const & essmSusy, double s, double tb, int l, DoubleMatrix & partials)
{
  // First check dimensions. If the number of columns is wrong, or if there are too few rows,
  // throw an exception - this method is intended to calculate the partials for a very definite parameter
  // set and so should not be used with other parameter sets in mind.
  int nRows = partials.displayRows();
  int nCols = partials.displayCols();

  if (nCols < 4 || nRows < 8)
    {
    ostringstream ii;
    ii << "Incorrect dimensions provided for partials matrix: needed 8x4 matrix, got " << nRows <<"x"<<nCols<< ".\n";
    throw(ii.str());
    }

  // Get parameter values
  DoubleVector lambda_vec = essmSusy.displaylambda();
  double lambda = lambda_vec(3);
  DoubleVector Alambda_vec = essmSusy.displayA_lambda();
  double Alambda = Alambda_vec(3);

  DoubleVector mH1Sq_vec = essmSusy.displayMh1Squared();
  double mH1Sq = mH1Sq_vec.display(3);
  DoubleVector mH2Sq_vec = essmSusy.displayMh2Squared();
  double mH2Sq = mH2Sq_vec.display(3);
  DoubleVector mSSq_vec = essmSusy.displaymSsq();
  double msSq = mSSq_vec.display(3);
  double At = essmSusy.displayA_top();
  double mQlSq = essmSusy.displaySoftMassSquared(mQl,3,3);
  double mUrSq =  essmSusy.displaySoftMassSquared(mUr,3,3);

  double v = essmSusy.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;

  // Calculate tree level partials
  double df1dlambda = lambda*v1*(v2*v2+s*s)-Alambda*s*v2/sqrt(2.0);
  double df2dlambda = lambda*v2*(v1*v1+s*s)-Alambda*s*v1/sqrt(2.0);
  double df3dlambda = lambda*v*v*s-Alambda*v1*v2/sqrt(2.0);

  double df1dAlambda = -lambda*s*v2/sqrt(2.0);
  double df2dAlambda = -lambda*s*v1/sqrt(2.0);
  double df3dAlambda = -lambda*v1*v2/sqrt(2.0);

  double df1dmH1Sq = v1;
  double df2dmH1Sq = 0.0;
  double df3dmH1Sq = 0.0;

  double df1dmH2Sq = 0.0;
  double df2dmH2Sq = v2;
  double df3dmH2Sq = 0.0;

  double df1dmsSq = 0.0;
  double df2dmsSq = 0.0;
  double df3dmsSq = s;

  double df1dmQlSq = 0.0;
  double df2dmQlSq = 0.0;
  double df3dmQlSq = 0.0;

  double df1dmUrSq = 0.0;
  double df2dmUrSq = 0.0;
  double df3dmUrSq = 0.0;

  double df1dAt = 0.0;
  double df2dAt = 0.0;
  double df3dAt = 0.0;

  // If higher order is requested, add in loop corrections. Note need to multiply by
  // appropriate VEV because of differing normalisations in helper functions below.
  if (l != 0)
    {
      if (l != 1)
	{
	  cout << "Warning: only l = 0 or l = 1 loops currently supported." << endl;
	  cout << "Using l = 1 loops." << endl;
	}

      df1dlambda = df1dlambda+v1*doCalcd2DeltaVdLambdadv1(essmSusy, s, tb);
      df2dlambda = df2dlambda+v2*doCalcd2DeltaVdLambdadv2(essmSusy, s, tb);
      df3dlambda = df3dlambda+s*doCalcd2DeltaVdLambdadv3(essmSusy, s, tb);

      df1dmQlSq = df1dmQlSq+v1*doCalcd2DeltaVdmQlsqdv1(essmSusy, s, tb);
      df2dmQlSq = df2dmQlSq+v2*doCalcd2DeltaVdmQlsqdv2(essmSusy, s, tb);
      df3dmQlSq = df3dmQlSq+s*doCalcd2DeltaVdmQlsqdv3(essmSusy, s, tb);

      df1dmUrSq = df1dmUrSq+v1*doCalcd2DeltaVdmUrsqdv1(essmSusy, s, tb);
      df2dmUrSq = df2dmUrSq+v2*doCalcd2DeltaVdmUrsqdv2(essmSusy, s, tb);
      df3dmUrSq = df3dmUrSq+s*doCalcd2DeltaVdmUrsqdv3(essmSusy, s, tb);

      df1dAt = df1dAt+v1*doCalcd2DeltaVdAtdv1(essmSusy, s, tb);
      df2dAt = df2dAt+v2*doCalcd2DeltaVdAtdv2(essmSusy, s, tb);
      df3dAt = df3dAt+s*doCalcd2DeltaVdAtdv3(essmSusy, s, tb);

    }

  // Construct partials matrix
  partials(1,1) = lambda;
  partials(1,2) = df1dlambda;
  partials(1,3) = df2dlambda;
  partials(1,4) = df3dlambda;

  partials(2,1) = Alambda;
  partials(2,2) = df1dAlambda;
  partials(2,3) = df2dAlambda;
  partials(2,4) = df3dAlambda;

  partials(3,1) = mH1Sq;
  partials(3,2) = df1dmH1Sq;
  partials(3,3) = df2dmH1Sq;
  partials(3,4) = df3dmH1Sq;

  partials(4,1) = mH2Sq;
  partials(4,2) = df1dmH2Sq;
  partials(4,3) = df2dmH2Sq;
  partials(4,4) = df3dmH2Sq;

  partials(5,1) = msSq;
  partials(5,2) = df1dmsSq;
  partials(5,3) = df2dmsSq;
  partials(5,4) = df3dmsSq;

  partials(6,1) = mQlSq;
  partials(6,2) = df1dmQlSq;
  partials(6,3) = df2dmQlSq;
  partials(6,4) = df3dmQlSq;

  partials(7,1) = mUrSq;
  partials(7,2) = df1dmUrSq;
  partials(7,3) = df2dmUrSq;
  partials(7,4) = df3dmUrSq;

  partials(8,1) = At;
  partials(8,2) = df1dAt;
  partials(8,3) = df2dAt;
  partials(8,4) = df3dAt;


}

// The following are helper functions that are useful for constructing
// the the leading one-loop tuning measures.

// doCalcdmstop1sqdvi calculates the derivative of m_stop_1^2 with respect to
// v_i for the given ESSM model, i = 1, 2, 3 (i=3 is s). 
// NB: m_stop_1^2 is taken to be the lighter stop in this method 
// -i.e. minus sign applies.
// Inputs:
//     SoftParsEssm const & essmSusy = object to calculate derivatives on
//     double s = the value of the VEV s to use
//     double tb = the value of tan(beta) to use
double doCalcdmstop1sqdv1(SoftParsEssm essmSusy, double s, double tb)
{
  double v = essmSusy.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;
  double g1 = essmSusy.displayGaugeCoupling(1);
  double g2 = essmSusy.displayGaugeCoupling(2);
  double gd1 = essmSusy.displaygdash_1();

  double yt = essmSusy.displayYukawaElement(YU, 3, 3);
  double mtop = yt*v2/sqrt(2.0);

  // Option to fix m_t(M_t) = 165 GeV
  bool usemtofmt = false;
  if(usemtofmt){
    mtop = 165;
    yt = sqrt(2.0)*mtop/v2;
    essmSusy.setYukawaElement(YU, 3, 3, yt);
  }

  DoubleVector lambda_vec = essmSusy.displaylambda();
  double lambda = lambda_vec(3);
  double At = essmSusy.displayA_top();
  double mQlsq = essmSusy.displaySoftMassSquared(mQl,3,3);
  double mUrsq =  essmSusy.displaySoftMassSquared(mUr,3,3);

  double MQQsq = mQlsq-mUrsq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2);
  double Xt = At - lambda*s/(sqrt(2.0)*tb);
  double rt = MQQsq*MQQsq+4.0*mtop*mtop*Xt*Xt;

  double deriv = 0.25*g2*g2+0.6*0.25*g1*g1-0.15*gd1*gd1;
  deriv = deriv-(1.0/sqrt(rt))*(MQQsq*0.25*(g2*g2-g1*g1)
				+4.0*mtop*mtop*Xt*(-lambda*s/(sqrt(2.0)*v1*v2)));
  deriv = deriv*(v1/2.0);

  return deriv;

}

double doCalcdmstop1sqdv2(SoftParsEssm essmSusy, double s, double tb)
{
  double v = essmSusy.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;
  double g1 = essmSusy.displayGaugeCoupling(1);
  double g2 = essmSusy.displayGaugeCoupling(2);
  double gd1 = essmSusy.displaygdash_1();

  double yt = essmSusy.displayYukawaElement(YU, 3, 3);
  double mtop = yt*v2/sqrt(2.0);
  // Option to fix m_t(M_t) = 165 GeV
  bool usemtofmt = false;
  if(usemtofmt){
    mtop = 165;
    yt = sqrt(2.0)*mtop/v2;
    essmSusy.setYukawaElement(YU, 3, 3, yt);
  }

  DoubleVector lambda_vec = essmSusy.displaylambda();
  double lambda = lambda_vec(3);
  double At = essmSusy.displayA_top();
  double mQlsq = essmSusy.displaySoftMassSquared(mQl,3,3);
  double mUrsq =  essmSusy.displaySoftMassSquared(mUr,3,3);

  double MQQsq = mQlsq-mUrsq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2);
  double Xt = At - lambda*s/(sqrt(2.0)*tb);
  double rt = MQQsq*MQQsq+4.0*mtop*mtop*Xt*Xt;

  double deriv = 2.0*yt*yt-0.25*g2*g2-0.6*0.25*g1*g1-0.1*gd1*gd1;
  deriv = deriv-(1.0/sqrt(rt))*(0.25*MQQsq*(g1*g1-g2*g2)+2.0*yt*yt*Xt*Xt
				+4.0*mtop*mtop*Xt*(lambda*s*v1/(sqrt(2.0)*v2*v2*v2)));
  deriv = deriv*(v2/2.0);

  return deriv;
}

double doCalcdmstop1sqdv3(SoftParsEssm essmSusy, double s, double tb)
{
  double v = essmSusy.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;
  double g1 = essmSusy.displayGaugeCoupling(1);
  double g2 = essmSusy.displayGaugeCoupling(2);
  double gd1 = essmSusy.displaygdash_1();

  double yt = essmSusy.displayYukawaElement(YU, 3, 3);
  double mtop = yt*v2/sqrt(2.0);
  // Option to fix m_t(M_t) = 165 GeV
  bool usemtofmt = false;
  if(usemtofmt){
    mtop = 165;
    yt = sqrt(2.0)*mtop/v2;
    essmSusy.setYukawaElement(YU, 3, 3, yt);
  }

  DoubleVector lambda_vec = essmSusy.displaylambda();
  double lambda = lambda_vec(3);
  double At = essmSusy.displayA_top();
  double mQlsq = essmSusy.displaySoftMassSquared(mQl,3,3);
  double mUrsq =  essmSusy.displaySoftMassSquared(mUr,3,3);

  double MQQsq = mQlsq-mUrsq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2);
  double Xt = At - lambda*s/(sqrt(2.0)*tb);
  double rt = MQQsq*MQQsq+4.0*mtop*mtop*Xt*Xt;

  double deriv = 0.25*gd1*gd1-(4.0*mtop*mtop/sqrt(rt))*Xt*(-lambda*v1/(sqrt(2.0)*v2*s));
  deriv = 0.5*s*deriv;

  return deriv;
}

// doCalcdmstop2sqdvi is as above but for m_stop_2^2 (+ sign applies).
double doCalcdmstop2sqdv1(SoftParsEssm essmSusy, double s, double tb)
{
  double v = essmSusy.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;
  double g1 = essmSusy.displayGaugeCoupling(1);
  double g2 = essmSusy.displayGaugeCoupling(2);
  double gd1 = essmSusy.displaygdash_1();

  double yt = essmSusy.displayYukawaElement(YU, 3, 3);
  double mtop = 165;//yt*v2/sqrt(2.0);
  // Option to fix m_t(M_t) = 165 GeV
  bool usemtofmt = false;
  if(usemtofmt){
    mtop = 165;
    yt = sqrt(2.0)*mtop/v2;
    essmSusy.setYukawaElement(YU, 3, 3, yt);
  }

  DoubleVector lambda_vec = essmSusy.displaylambda();
  double lambda = lambda_vec(3);
  double At = essmSusy.displayA_top();
  double mQlsq = essmSusy.displaySoftMassSquared(mQl,3,3);
  double mUrsq =  essmSusy.displaySoftMassSquared(mUr,3,3);

  double MQQsq = mQlsq-mUrsq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2);
  double Xt = At - lambda*s/(sqrt(2.0)*tb);
  double rt = MQQsq*MQQsq+4.0*mtop*mtop*Xt*Xt;

  double deriv = 0.25*g2*g2+0.6*0.25*g1*g1-0.15*gd1*gd1;
  deriv = deriv+(1.0/sqrt(rt))*(MQQsq*0.25*(g2*g2-g1*g1)
				+4.0*mtop*mtop*Xt*(-lambda*s/(sqrt(2.0)*v1*v2)));
  deriv = deriv*(v1/2.0);

  return deriv;
}

double doCalcdmstop2sqdv2(SoftParsEssm essmSusy, double s, double tb)
{
  double v = essmSusy.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;
  double g1 = essmSusy.displayGaugeCoupling(1);
  double g2 = essmSusy.displayGaugeCoupling(2);
  double gd1 = essmSusy.displaygdash_1();

  double yt = essmSusy.displayYukawaElement(YU, 3, 3);
  double mtop = yt*v2/sqrt(2.0);
  // Option to fix m_t(M_t) = 165 GeV
  bool usemtofmt = false;
  if(usemtofmt){
    mtop = 165;
    yt = sqrt(2.0)*mtop/v2;
    essmSusy.setYukawaElement(YU, 3, 3, yt);
  }

  DoubleVector lambda_vec = essmSusy.displaylambda();
  double lambda = lambda_vec(3);
  double At = essmSusy.displayA_top();
  double mQlsq = essmSusy.displaySoftMassSquared(mQl,3,3);
  double mUrsq =  essmSusy.displaySoftMassSquared(mUr,3,3);

  double MQQsq = mQlsq-mUrsq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2);
  double Xt = At - lambda*s/(sqrt(2.0)*tb);
  double rt = MQQsq*MQQsq+4.0*mtop*mtop*Xt*Xt;

  double deriv = 2.0*yt*yt-0.25*g2*g2-0.6*0.25*g1*g1-0.1*gd1*gd1;
  deriv = deriv+(1.0/sqrt(rt))*(0.25*MQQsq*(g1*g1-g2*g2)+2.0*yt*yt*Xt*Xt
				+4.0*mtop*mtop*Xt*(lambda*s*v1/(sqrt(2.0)*v2*v2*v2)));
  deriv = deriv*(v2/2.0);

  return deriv;
}

double doCalcdmstop2sqdv3(SoftParsEssm essmSusy, double s, double tb)
{
  double v = essmSusy.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;
  double g1 = essmSusy.displayGaugeCoupling(1);
  double g2 = essmSusy.displayGaugeCoupling(2);
  double gd1 = essmSusy.displaygdash_1();

  double yt = essmSusy.displayYukawaElement(YU, 3, 3);
  double mtop = yt*v2/sqrt(2.0);
  // Option to fix m_t(M_t) = 165 GeV
  bool usemtofmt = false;
  if(usemtofmt){
    mtop = 165;
    yt = sqrt(2.0)*mtop/v2;
    essmSusy.setYukawaElement(YU, 3, 3, yt);
  }

  DoubleVector lambda_vec = essmSusy.displaylambda();
  double lambda = lambda_vec(3);
  double At = essmSusy.displayA_top();
  double mQlsq = essmSusy.displaySoftMassSquared(mQl,3,3);
  double mUrsq =  essmSusy.displaySoftMassSquared(mUr,3,3);

  double MQQsq = mQlsq-mUrsq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2);
  double Xt = At - lambda*s/(sqrt(2.0)*tb);
  double rt = MQQsq*MQQsq+4.0*mtop*mtop*Xt*Xt;

  double deriv = 0.25*gd1*gd1+(4.0*mtop*mtop/sqrt(rt))*Xt*(-lambda*v1/(sqrt(2.0)*v2*s));
  deriv = 0.5*s*deriv;

  return deriv;
}

// The three functions below are purely used for testing the first
// derivatives of the stop squared masses. They calculate the 
// stop/top tadpole correction to the effective potential, in the
// form (1/v_i)\frac{\partial \Delta V}{\partial v_i}. The results
// should be compared against those obtained using the functions
// doCalcTadpolesESSMH<1,2,S> which are actually used in practical
// calculations. In the notation used in my notes, these functions
// calculate \Delta^{tadpole}_i for i = 1, 2, 3.
// Inputs:
//     SoftParsEssm const & essmSusy = object to calculate derivatives on
//     double s = the value of the VEV s to use
//     double tb = the value of tan(beta) to use
double calculateDeltaTadpole1(SoftParsEssm essmSusy, double s, double tb)
{
  double v = essmSusy.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;
  double g1 = essmSusy.displayGaugeCoupling(1);
  double g2 = essmSusy.displayGaugeCoupling(2);
  double gd1 = essmSusy.displaygdash_1();

  double q = essmSusy.displayMu();

  double yt = essmSusy.displayYukawaElement(YU, 3, 3);
  double mtop = yt*v2/sqrt(2.0);
  // Option to fix m_t(M_t) = 165 GeV
  bool usemtofmt = false;
  if(usemtofmt){
    mtop = 165;
    yt = sqrt(2.0)*mtop/v2;
    essmSusy.setYukawaElement(YU, 3, 3, yt);
  }

  DoubleVector lambda_vec = essmSusy.displaylambda();
  double lambda = lambda_vec(3);
  double At = essmSusy.displayA_top();
  double mQlsq = essmSusy.displaySoftMassSquared(mQl,3,3);
  double mUrsq =  essmSusy.displaySoftMassSquared(mUr,3,3);

  double MQQsq = mQlsq-mUrsq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2);
  double Xt = At - lambda*s/(sqrt(2.0)*tb);
  double rt = MQQsq*MQQsq+4.0*mtop*mtop*Xt*Xt;

  // Contributions from auxiliary D-terms. Might in future
  // want to use effective charges instead of U(1)_N charges.
  double DeltaQ3 = (gd1*gd1/80.0)*(-3.0*v1*v1-2.0*v2*v2+5.0*s*s);
  double Deltau3 = (gd1*gd1/80.0)*(-3.0*v1*v1-2.0*v2*v2+5.0*s*s);

  // The stop squared masses, in the notation m_stop_1^2 has -ve
  // sign contribution (is lighter).
  double mstop1sq = 0.5*(mQlsq+mUrsq+0.125*(g2*g2+0.6*g1*g1)*(v1*v1-v2*v2)
		     +DeltaQ3+Deltau3+2.0*mtop*mtop-sqrt(rt));
  double mstop2sq = 0.5*(mQlsq+mUrsq+0.125*(g2*g2+0.6*g1*g1)*(v1*v1-v2*v2)
		     +DeltaQ3+Deltau3+2.0*mtop*mtop+sqrt(rt));

  bool speak = false;
  if (speak)
    {
      cout << "Calculated stop masses: " << endl;
      cout << "m_stop_1^2 = " << mstop1sq << ", ";
      cout << "m_stop_1 = " << sqrt(mstop1sq) << ", ";
      cout << "m_stop_2^2 = " << mstop2sq << ", ";
      cout << "m_stop_2 = " << sqrt(mstop2sq) << endl;
    }


  // The derivatives of the stop masses wrt the VEVs.
  double dmstop1sqdv1 = doCalcdmstop1sqdv1(essmSusy, s, tb);
  double dmstop2sqdv1 = doCalcdmstop2sqdv1(essmSusy, s, tb);

  // The derivative of the top mass wrt the VEVs. This is non-zero
  // only for dm_top/dv2.
  double dmtopsqdv1 = 0.0;

  if (speak)
    {
      cout << "Calculated derivatives of stop/top masses: " << endl;
      cout << "dm_stop_1^2/dv1 = " << dmstop1sqdv1 << ", ";
      cout << "dm_stop_2^2/dv1 = " << dmstop2sqdv1 << ", ";
      cout << "dm_top/dv1 = " << dmtopsqdv1 << endl;
    }

  // The tadpole contribution in the form (1/v_i)\frac{\partial \Delta V}{\partial v_i}.
  double dDeltaVdv1 = (3.0/(32.0*PI*PI))*(2.0*a0Peter(mstop1sq, q)*dmstop1sqdv1
					  +2.0*a0Peter(mstop2sq, q)*dmstop2sqdv1
					  -4.0*a0Peter(mtop*mtop, q)*dmtopsqdv1);

  double deltaTadpole1 = (1.0/v1)*dDeltaVdv1;

  if (speak)
    {
      cout << "Calculated tadpole contributions: " << endl;
      cout << "dDeltaV/dv1 = " << dDeltaVdv1 << ", ";
      cout << "Delta^{tadpole}_1 = " << deltaTadpole1 << endl;
    }

  return deltaTadpole1;

}

double calculateDeltaTadpole2(SoftParsEssm essmSusy, double s, double tb)
{
  double v = essmSusy.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;
  double g1 = essmSusy.displayGaugeCoupling(1);
  double g2 = essmSusy.displayGaugeCoupling(2);
  double gd1 = essmSusy.displaygdash_1();

  double q = essmSusy.displayMu();

  double yt = essmSusy.displayYukawaElement(YU, 3, 3);
  double mtop = yt*v2/sqrt(2.0);
  // Option to fix m_t(M_t) = 165 GeV
  bool usemtofmt = false;
  if(usemtofmt){
    mtop = 165;
    yt = sqrt(2.0)*mtop/v2;
    essmSusy.setYukawaElement(YU, 3, 3, yt);
  }

  DoubleVector lambda_vec = essmSusy.displaylambda();
  double lambda = lambda_vec(3);
  double At = essmSusy.displayA_top();
  double mQlsq = essmSusy.displaySoftMassSquared(mQl,3,3);
  double mUrsq =  essmSusy.displaySoftMassSquared(mUr,3,3);

  double MQQsq = mQlsq-mUrsq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2);
  double Xt = At - lambda*s/(sqrt(2.0)*tb);
  double rt = MQQsq*MQQsq+4.0*mtop*mtop*Xt*Xt;

  // Contributions from auxiliary D-terms. Might in future
  // want to use effective charges instead of U(1)_N charges.
  double DeltaQ3 = (gd1*gd1/80.0)*(-3.0*v1*v1-2.0*v2*v2+5.0*s*s);
  double Deltau3 = (gd1*gd1/80.0)*(-3.0*v1*v1-2.0*v2*v2+5.0*s*s);

  // The stop squared masses, in the notation m_stop_1^2 has -ve
  // sign contribution (is lighter).
  double mstop1sq = 0.5*(mQlsq+mUrsq+0.125*(g2*g2+0.6*g1*g1)*(v1*v1-v2*v2)
		     +DeltaQ3+Deltau3+2.0*mtop*mtop-sqrt(rt));
  double mstop2sq = 0.5*(mQlsq+mUrsq+0.125*(g2*g2+0.6*g1*g1)*(v1*v1-v2*v2)
		     +DeltaQ3+Deltau3+2.0*mtop*mtop+sqrt(rt));

  bool speak = false;
  if (speak)
    {
      cout << "Calculated stop masses: " << endl;
      cout << "m_stop_1^2 = " << mstop1sq << ", ";
      cout << "m_stop_1 = " << sqrt(mstop1sq) << ", ";
      cout << "m_stop_2^2 = " << mstop2sq << ", ";
      cout << "m_stop_2 = " << sqrt(mstop2sq) << endl;
    }

  // The derivatives of the stop masses wrt the VEVs.
  double dmstop1sqdv2 = doCalcdmstop1sqdv2(essmSusy, s, tb);
  double dmstop2sqdv2 = doCalcdmstop2sqdv2(essmSusy, s, tb);

  // The derivative of the top mass wrt the VEVs. This is non-zero
  // only for dm_top/dv2.
  double dmtopsqdv2 = yt*yt*v2;

  if (speak)
    {
      cout << "Calculated derivatives of stop/top masses: " << endl;
      cout << "dm_stop_1^2/dv2 = " << dmstop1sqdv2 << ", ";
      cout << "dm_stop_2^2/dv2 = " << dmstop2sqdv2 << ", ";
      cout << "dm_top/dv2 = " << dmtopsqdv2 << endl;
    }

  // The tadpole contribution in the form (1/v_i)\frac{\partial \Delta V}{\partial v_i}.
  double dDeltaVdv2 = (3.0/(32.0*PI*PI))*(2.0*a0Peter(mstop1sq, q)*dmstop1sqdv2
					  +2.0*a0Peter(mstop2sq, q)*dmstop2sqdv2
					  -4.0*a0Peter(mtop*mtop, q)*dmtopsqdv2);

  double deltaTadpole2 = (1.0/v2)*dDeltaVdv2;

  if (speak)
    {
      cout << "Calculated tadpole contributions: " << endl;
      cout << "dDeltaV/dv2 = " << dDeltaVdv2 << ", ";
      cout << "Delta^{tadpole}_2 = " << deltaTadpole2 << endl;
    }

  return deltaTadpole2;

}

double calculateDeltaTadpole3(SoftParsEssm essmSusy, double s, double tb)
{
  double v = essmSusy.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;
  double g1 = essmSusy.displayGaugeCoupling(1);
  double g2 = essmSusy.displayGaugeCoupling(2);
  double gd1 = essmSusy.displaygdash_1();

  double q = essmSusy.displayMu();

  double yt = essmSusy.displayYukawaElement(YU, 3, 3);
  double mtop = yt*v2/sqrt(2.0);
  // Option to fix m_t(M_t) = 165 GeV
  bool usemtofmt = false;
  if(usemtofmt){
    mtop = 165;
    yt = sqrt(2.0)*mtop/v2;
    essmSusy.setYukawaElement(YU, 3, 3, yt);
  }

  DoubleVector lambda_vec = essmSusy.displaylambda();
  double lambda = lambda_vec(3);
  double At = essmSusy.displayA_top();
  double mQlsq = essmSusy.displaySoftMassSquared(mQl,3,3);
  double mUrsq =  essmSusy.displaySoftMassSquared(mUr,3,3);

  double MQQsq = mQlsq-mUrsq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2);
  double Xt = At - lambda*s/(sqrt(2.0)*tb);
  double rt = MQQsq*MQQsq+4.0*mtop*mtop*Xt*Xt;

  // Contributions from auxiliary D-terms. Might in future
  // want to use effective charges instead of U(1)_N charges.
  double DeltaQ3 = (gd1*gd1/80.0)*(-3.0*v1*v1-2.0*v2*v2+5.0*s*s);
  double Deltau3 = (gd1*gd1/80.0)*(-3.0*v1*v1-2.0*v2*v2+5.0*s*s);

  // The stop squared masses, in the notation m_stop_1^2 has -ve
  // sign contribution (is lighter).
  double mstop1sq = 0.5*(mQlsq+mUrsq+0.125*(g2*g2+0.6*g1*g1)*(v1*v1-v2*v2)
		     +DeltaQ3+Deltau3+2.0*mtop*mtop-sqrt(rt));
  double mstop2sq = 0.5*(mQlsq+mUrsq+0.125*(g2*g2+0.6*g1*g1)*(v1*v1-v2*v2)
		     +DeltaQ3+Deltau3+2.0*mtop*mtop+sqrt(rt));

  bool speak = false;
  if (speak)
    {
      cout << "Calculated stop masses: " << endl;
      cout << "m_stop_1^2 = " << mstop1sq << ", ";
      cout << "m_stop_1 = " << sqrt(mstop1sq) << ", ";
      cout << "m_stop_2^2 = " << mstop2sq << ", ";
      cout << "m_stop_2 = " << sqrt(mstop2sq) << endl;
    }

  // The derivatives of the stop masses wrt the VEVs.
  double dmstop1sqdv3 = doCalcdmstop1sqdv3(essmSusy, s, tb);
  double dmstop2sqdv3 = doCalcdmstop2sqdv3(essmSusy, s, tb);

  // The derivative of the top mass wrt the VEVs. This is non-zero
  // only for dm_top/dv2.
  double dmtopsqdv3 = 0.0;

  if (speak)
    {
      cout << "Calculated derivatives of stop/top masses: " << endl;
      cout << "dm_stop_1^2/ds = " << dmstop1sqdv3 << ", ";
      cout << "dm_stop_2^2/ds = " << dmstop2sqdv3 << ", ";
      cout << "dm_top/ds = " << dmtopsqdv3 << endl;
    }

  // The tadpole contribution in the form (1/v_i)\frac{\partial \Delta V}{\partial v_i}.
  double dDeltaVdv3 = (3.0/(32.0*PI*PI))*(2.0*a0Peter(mstop1sq, q)*dmstop1sqdv3
					  +2.0*a0Peter(mstop2sq, q)*dmstop2sqdv3
					  -4.0*a0Peter(mtop*mtop, q)*dmtopsqdv3);

  double deltaTadpole3 = (1.0/s)*dDeltaVdv3;

  if (speak)
    {
      cout << "Calculated tadpole contributions: " << endl;
      cout << "dDeltaV/ds = " << dDeltaVdv3 << ", ";
      cout << "Delta^{tadpole}_3 = " << deltaTadpole3 << endl;
    }

  return deltaTadpole3;
}

// These functions calculate the leading one-loop corrections
// to the CP-even Higgs mass matrix in the basis (v_1, v_2, v_3=s),
// given by \Delta_{ij}'=\frac{\partial^2\Delta V}{\partial v_i\partial v_j}.
// Inputs:
//     SoftParsEssm essmSusy = object to calculate derivatives on
//     double s = the value of the VEV s to use
//     double tb = the value of tan(beta) to use.
double doCalcDeltaPrime11(SoftParsEssm essmSusy, double s, double tb)
{
  double v = essmSusy.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;
  double g1 = essmSusy.displayGaugeCoupling(1);
  double g2 = essmSusy.displayGaugeCoupling(2);
  double gbar = sqrt(g2*g2+0.6*g1*g1);
  double gd1 = essmSusy.displaygdash_1();

  // Logarithms are calculated at the given Q
  double q = essmSusy.displayMu();

  double yt = essmSusy.displayYukawaElement(YU, 3, 3);
  double mtop = yt*v2/sqrt(2.0);
  // Option to fix m_t(M_t) = 165 GeV
  bool usemtofmt = false;
  if(usemtofmt){
    mtop = 165;
    yt = sqrt(2.0)*mtop/v2;
    essmSusy.setYukawaElement(YU, 3, 3, yt);
  }

  DoubleVector lambda_vec = essmSusy.displaylambda();
  double lambda = lambda_vec(3);
  double At = essmSusy.displayA_top();
  double mQlsq = essmSusy.displaySoftMassSquared(mQl,3,3);
  double mUrsq =  essmSusy.displaySoftMassSquared(mUr,3,3);

  double MQQsq = mQlsq-mUrsq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2);
  double Xt = At - lambda*s/(sqrt(2.0)*tb);
  double rt = MQQsq*MQQsq+4.0*mtop*mtop*Xt*Xt;
  double RQQ = MQQsq*(g2*g2-g1*g1);

  // Squared stop masses are calculated using physical_ESSM to make sure
  // auxiliary D-terms are included.
  DoubleVector mstop(2), mstopsq(2), mD1sq(3), mD2sq(3);
  physical_ESSM(essmSusy, mstop, mstopsq, mD1sq, mD2sq, s, tb);

  // Have to swap the mass ordering of the stops, because
  // physical_ESSM uses the opposite sign.
  double mstop1sq, mstop2sq;
  mstop1sq = mstopsq(2);
  mstop2sq = mstopsq(1);

  double term1, term2, term3, term4;

  term1 = v1*v1*sqr(0.125*gbar*gbar-0.075*gd1*gd1)
    +(1.0/rt)*sqr(0.125*v1*RQQ-2.0*mtop*mtop*Xt*s*lambda/(sqrt(2.0)*v2));
  term1 = term1*log(mstop1sq*mstop2sq/(q*q*q*q));

  term2 = (v1/(32.0*sqrt(rt)))*(gbar*gbar-0.6*gd1*gd1)
    *(v1*RQQ-16.0*mtop*mtop*Xt*s*lambda/(sqrt(2.0)*v2));
  term2 = term2*log(mstop2sq/mstop1sq);

  term3 = (0.125*gbar*gbar-0.075*gd1*gd1)*(a0Peter(mstop1sq, q)+a0Peter(mstop2sq,q));

  term4 = (1.0/sqrt(rt))*(4.0*RQQ+sqr(g2*g2-g1*g1)*v1*v1+16.0*yt*yt*s*s*lambda*lambda);
  term4 = term4-(1.0/(rt*sqrt(rt)))*sqr(v1*RQQ-16.0*mtop*mtop*Xt*s*lambda/(sqrt(2.0)*v2));
  term4 = term4*(a0Peter(mstop2sq,q)-a0Peter(mstop1sq,q))/32.0;

  double Del11prime = (3.0/(16.0*PI*PI))*(term1+term2+term3+term4);

  return Del11prime;
}

double doCalcDeltaPrime12(SoftParsEssm essmSusy, double s, double tb)
{
  double v = essmSusy.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;
  double g1 = essmSusy.displayGaugeCoupling(1);
  double g2 = essmSusy.displayGaugeCoupling(2);
  double gbar = sqrt(g2*g2+0.6*g1*g1);
  double gd1 = essmSusy.displaygdash_1();

  // Logarithms are calculated at the given Q
  double q = essmSusy.displayMu();

  double yt = essmSusy.displayYukawaElement(YU, 3, 3);
  double mtop = yt*v2/sqrt(2.0);
  // Option to fix m_t(M_t) = 165 GeV
  bool usemtofmt = false;
  if(usemtofmt){
    mtop = 165;
    yt = sqrt(2.0)*mtop/v2;
    essmSusy.setYukawaElement(YU, 3, 3, yt);
  }

  DoubleVector lambda_vec = essmSusy.displaylambda();
  double lambda = lambda_vec(3);
  double At = essmSusy.displayA_top();
  double mQlsq = essmSusy.displaySoftMassSquared(mQl,3,3);
  double mUrsq =  essmSusy.displaySoftMassSquared(mUr,3,3);
  double mueff = lambda*s/sqrt(2.0);

  double MQQsq = mQlsq-mUrsq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2);
  double Xt = At - lambda*s/(sqrt(2.0)*tb);
  double rt = MQQsq*MQQsq+4.0*mtop*mtop*Xt*Xt;
  double RQQ = MQQsq*(g2*g2-g1*g1);

  // Squared stop masses are calculated using physical_ESSM to make sure
  // auxiliary D-terms are included.
  DoubleVector mstop(2), mstopsq(2), mD1sq(3), mD2sq(3);
  physical_ESSM(essmSusy, mstop, mstopsq, mD1sq, mD2sq, s, tb);

  // Have to swap the mass ordering of the stops, because
  // physical_ESSM uses the opposite sign.
  double mstop1sq, mstop2sq;
  mstop1sq = mstopsq(2);
  mstop2sq = mstopsq(1);

  double term1, term2, term3;

  term1 = (0.125*gbar*gbar-0.075*gd1*gd1)*(yt*yt-0.125*gbar*gbar-0.05*gd1*gd1)
    +(1.0/rt)*(0.125*RQQ-yt*yt*Xt*mueff*tb)*(yt*yt*Xt*At-0.125*RQQ);
  term1 = term1*v1*v2*log(mstop1sq*mstop2sq/(q*q*q*q));

  term2 = (0.125*gbar*gbar-0.075*gd1*gd1)*(yt*yt*Xt*At-0.125*RQQ)
    +(0.125*RQQ-yt*yt*Xt*mueff*tb)*(yt*yt-0.125*gbar*gbar-0.05*gd1*gd1);
  term2 = term2*(v1*v2/sqrt(rt))*log(mstop2sq/mstop1sq);

  term3 = (sqr(g2*g2-g1*g1)*v1*v2/32.0+At*yt*yt*mueff)
    +(0.125*RQQ-yt*yt*Xt*mueff*tb)*(2.0*yt*yt*Xt*At-0.25*RQQ)*(v1*v2/rt);
  term3= term3*(1.0/sqrt(rt))*(a0Peter(mstop2sq,q)-a0Peter(mstop1sq,q));

  double Del12prime = (3.0/(16.0*PI*PI))*(term1+term2-term3);

  return Del12prime;
}

double doCalcDeltaPrime13(SoftParsEssm essmSusy, double s, double tb)
{
  double v = essmSusy.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;
  double g1 = essmSusy.displayGaugeCoupling(1);
  double g2 = essmSusy.displayGaugeCoupling(2);
  double gbar = sqrt(g2*g2+0.6*g1*g1);
  double gd1 = essmSusy.displaygdash_1();

  // Logarithms are calculated at the given Q
  double q = essmSusy.displayMu();

  double yt = essmSusy.displayYukawaElement(YU, 3, 3);
  double mtop = yt*v2/sqrt(2.0);
  // Option to fix m_t(M_t) = 165 GeV
  bool usemtofmt = false;
  if(usemtofmt){
    mtop = 165;
    yt = sqrt(2.0)*mtop/v2;
    essmSusy.setYukawaElement(YU, 3, 3, yt);
  }

  DoubleVector lambda_vec = essmSusy.displaylambda();
  double lambda = lambda_vec(3);
  double At = essmSusy.displayA_top();
  double mQlsq = essmSusy.displaySoftMassSquared(mQl,3,3);
  double mUrsq =  essmSusy.displaySoftMassSquared(mUr,3,3);
  double mueff = lambda*s/sqrt(2.0);

  double MQQsq = mQlsq-mUrsq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2);
  double Xt = At - lambda*s/(sqrt(2.0)*tb);
  double rt = MQQsq*MQQsq+4.0*mtop*mtop*Xt*Xt;
  double RQQ = MQQsq*(g2*g2-g1*g1);

  // Squared stop masses are calculated using physical_ESSM to make sure
  // auxiliary D-terms are included.
  DoubleVector mstop(2), mstopsq(2), mD1sq(3), mD2sq(3);
  physical_ESSM(essmSusy, mstop, mstopsq, mD1sq, mD2sq, s, tb);

  // Have to swap the mass ordering of the stops, because
  // physical_ESSM uses the opposite sign.
  double mstop1sq, mstop2sq;
  mstop1sq = mstopsq(2);
  mstop2sq = mstopsq(1);

  double term1, term2, term3;

  term1 = 0.125*s*v1*gd1*gd1*(0.125*gbar*gbar-0.075*gd1*gd1)
    -(0.125*RQQ*v1-yt*yt*Xt*mueff*v2)*(2.0*mtop*mtop*Xt*lambda)/(rt*sqrt(2.0)*tb);
  term1 = term1*log(mstop1sq*mstop2sq/(q*q*q*q));

  term2 = (0.125*gbar*gbar-0.075*gd1*gd1)*(2.0*v1*mtop*mtop*Xt*lambda)/(tb*sqrt(2.0*rt))
    -(0.125*s*gd1*gd1/sqrt(rt))*(0.125*RQQ*v1-2.0*mtop*mtop*Xt*s*lambda/(sqrt(2.0)*v2));
  term2 = term2*log(mstop2sq/mstop1sq);

  term3 = (yt*yt*v1*lambda*mueff/sqrt(2.0*rt))
    *(1.0-Xt*tb/mueff-4.0*Xt*Xt*mtop*mtop/rt+v1*v2*RQQ*Xt/(4.0*mueff*rt));
  term3 = term3*(a0Peter(mstop2sq, q)-a0Peter(mstop1sq,q));

  double Del13prime = (3.0/(16.0*PI*PI))*(term1-term2+term3);

  return Del13prime;
}

double doCalcDeltaPrime22(SoftParsEssm essmSusy, double s, double tb)
{
  double v = essmSusy.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;
  double g1 = essmSusy.displayGaugeCoupling(1);
  double g2 = essmSusy.displayGaugeCoupling(2);
  double gbar = sqrt(g2*g2+0.6*g1*g1);
  double gd1 = essmSusy.displaygdash_1();

  // Logarithms are calculated at the given Q
  double q = essmSusy.displayMu();

  double yt = essmSusy.displayYukawaElement(YU, 3, 3);
  double mtop = yt*v2/sqrt(2.0);
  // Option to fix m_t(M_t) = 165 GeV
  bool usemtofmt = false;
  if(usemtofmt){
    mtop = 165;
    yt = sqrt(2.0)*mtop/v2;
    essmSusy.setYukawaElement(YU, 3, 3, yt);
  }

  DoubleVector lambda_vec = essmSusy.displaylambda();
  double lambda = lambda_vec(3);
  double At = essmSusy.displayA_top();
  double mQlsq = essmSusy.displaySoftMassSquared(mQl,3,3);
  double mUrsq =  essmSusy.displaySoftMassSquared(mUr,3,3);
  double mueff = lambda*s/sqrt(2.0);

  double MQQsq = mQlsq-mUrsq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2);
  double Xt = At - lambda*s/(sqrt(2.0)*tb);
  double rt = MQQsq*MQQsq+4.0*mtop*mtop*Xt*Xt;
  double RQQ = MQQsq*(g2*g2-g1*g1);

  // Squared stop masses are calculated using physical_ESSM to make sure
  // auxiliary D-terms are included.
  DoubleVector mstop(2), mstopsq(2), mD1sq(3), mD2sq(3);
  physical_ESSM(essmSusy, mstop, mstopsq, mD1sq, mD2sq, s, tb);

  // Have to swap the mass ordering of the stops, because
  // physical_ESSM uses the opposite sign.
  double mstop1sq, mstop2sq;
  mstop1sq = mstopsq(2);
  mstop2sq = mstopsq(1);

  double term1, term2, term3, term4, term5,term6;

  term1 = sqr(yt*yt-0.125*gbar*gbar-0.05*gd1*gd1)
    +sqr(8.0*Xt*At*yt*yt-RQQ)/(64.0*rt);
  term1 = term1*v2*v2*log(mstop1sq*mstop2sq/(q*q*q*q));

  term2 = (v2*v2/(4.0*sqrt(rt)))*(yt*yt-0.125*gbar*gbar-0.05*gd1*gd1)
    *(8.0*yt*yt*Xt*At-RQQ)*log(mstop2sq/mstop1sq);

  term3 = (yt*yt-0.125*gbar*gbar-0.05*gd1*gd1)*(a0Peter(mstop1sq,q)+a0Peter(mstop2sq,q));

  term4 = (1.0/sqrt(rt))*(sqr(g2*g2-g1*g1)*v2*v2/32.0-0.125*RQQ+yt*yt*At*At
			  -(v2*v2*sqr(8.0*Xt*At*yt*yt-RQQ))/(32.0*rt));
  term4 = term4*(a0Peter(mstop2sq,q)-a0Peter(mstop1sq,q));

  term5 = 2.0*yt*yt*yt*yt*v2*v2*log(mtop*mtop/(q*q));

  term6 = 2.0*yt*yt*a0Peter(mtop*mtop, q);

  double Del22prime = (3.0/(16.0*PI*PI))*(term1+term2+term3+term4-term5-term6);

  return Del22prime;
}

double doCalcDeltaPrime23(SoftParsEssm essmSusy, double s, double tb)
{
  double v = essmSusy.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;
  double g1 = essmSusy.displayGaugeCoupling(1);
  double g2 = essmSusy.displayGaugeCoupling(2);
  double gbar = sqrt(g2*g2+0.6*g1*g1);
  double gd1 = essmSusy.displaygdash_1();

  // Logarithms are calculated at the given Q
  double q = essmSusy.displayMu();

  double yt = essmSusy.displayYukawaElement(YU, 3, 3);
  double mtop = yt*v2/sqrt(2.0);
  // Option to fix m_t(M_t) = 165 GeV
  bool usemtofmt = false;
  if(usemtofmt){
    mtop = 165;
    yt = sqrt(2.0)*mtop/v2;
    essmSusy.setYukawaElement(YU, 3, 3, yt);
  }

  DoubleVector lambda_vec = essmSusy.displaylambda();
  double lambda = lambda_vec(3);
  double At = essmSusy.displayA_top();
  double mQlsq = essmSusy.displaySoftMassSquared(mQl,3,3);
  double mUrsq =  essmSusy.displaySoftMassSquared(mUr,3,3);
  double mueff = lambda*s/sqrt(2.0);

  double MQQsq = mQlsq-mUrsq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2);
  double Xt = At - lambda*s/(sqrt(2.0)*tb);
  double rt = MQQsq*MQQsq+4.0*mtop*mtop*Xt*Xt;
  double RQQ = MQQsq*(g2*g2-g1*g1);

  // Squared stop masses are calculated using physical_ESSM to make sure
  // auxiliary D-terms are included.
  DoubleVector mstop(2), mstopsq(2), mD1sq(3), mD2sq(3);
  physical_ESSM(essmSusy, mstop, mstopsq, mD1sq, mD2sq, s, tb);

  // Have to swap the mass ordering of the stops, because
  // physical_ESSM uses the opposite sign.
  double mstop1sq, mstop2sq;
  mstop1sq = mstopsq(2);
  mstop2sq = mstopsq(1);

  double term1, term2, term3;

  term1 = 0.125*s*gd1*gd1*(yt*yt-0.125*gbar*gbar-0.05*gd1*gd1)
    -(2.0*mtop*mtop*Xt*lambda)*(yt*yt*Xt*At-0.125*RQQ)/(sqrt(2.0)*rt*tb);
  term1 = term1*v2*log(mstop1sq*mstop2sq/(q*q*q*q));
 
  term2 = (1.0/sqrt(rt))*(0.125*s*gd1*gd1*(yt*yt*Xt*At-0.125*RQQ)
			  -(2.0*mtop*mtop*Xt*lambda)*(yt*yt-0.125*gbar*gbar-0.05*gd1*gd1)/(sqrt(2.0)*tb));
  term2 = term2*v2*log(mstop2sq/mstop1sq);

  term3 = ((yt*yt*lambda*v1*At)/sqrt(2.0*rt))*(1.0-4.0*mtop*mtop*Xt*Xt/rt+v2*v2*Xt*RQQ/(4.0*At*rt));
  term3 = term3*(a0Peter(mstop2sq,q)-a0Peter(mstop1sq,q));
  
  double Del23prime = (3.0/(16.0*PI*PI))*(term1+term2-term3);

  return Del23prime;
}

double doCalcDeltaPrime33(SoftParsEssm essmSusy, double s, double tb)
{
  double v = essmSusy.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;
  double g1 = essmSusy.displayGaugeCoupling(1);
  double g2 = essmSusy.displayGaugeCoupling(2);
  double gbar = sqrt(g2*g2+0.6*g1*g1);
  double gd1 = essmSusy.displaygdash_1();

  // Logarithms are calculated at the given Q
  double q = essmSusy.displayMu();

  double yt = essmSusy.displayYukawaElement(YU, 3, 3);
  double mtop = yt*v2/sqrt(2.0);
  // Option to fix m_t(M_t) = 165 GeV
  bool usemtofmt = false;
  if(usemtofmt){
    mtop = 165;
    yt = sqrt(2.0)*mtop/v2;
    essmSusy.setYukawaElement(YU, 3, 3, yt);
  }

  DoubleVector lambda_vec = essmSusy.displaylambda();
  double lambda = lambda_vec(3);
  double At = essmSusy.displayA_top();
  double mQlsq = essmSusy.displaySoftMassSquared(mQl,3,3);
  double mUrsq =  essmSusy.displaySoftMassSquared(mUr,3,3);
  double mueff = lambda*s/sqrt(2.0);

  double MQQsq = mQlsq-mUrsq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2);
  double Xt = At - lambda*s/(sqrt(2.0)*tb);
  double rt = MQQsq*MQQsq+4.0*mtop*mtop*Xt*Xt;
  double RQQ = MQQsq*(g2*g2-g1*g1);

  // Squared stop masses are calculated using physical_ESSM to make sure
  // auxiliary D-terms are included.
  DoubleVector mstop(2), mstopsq(2), mD1sq(3), mD2sq(3);
  physical_ESSM(essmSusy, mstop, mstopsq, mD1sq, mD2sq, s, tb);

  // Have to swap the mass ordering of the stops, because
  // physical_ESSM uses the opposite sign.
  double mstop1sq, mstop2sq;
  mstop1sq = mstopsq(2);
  mstop2sq = mstopsq(1);

  double term1, term2, term3, term4;

  term1 = sqr(gd1*gd1*s)/64.0+(2.0*sqr(mtop*mtop)*Xt*Xt*lambda*lambda)/(rt*tb*tb);
  term1 = term1*log(mstop1sq*mstop2sq/(q*q*q*q));

  term2 = (gd1*gd1*mtop*mtop*Xt*mueff)/(2.0*sqrt(rt)*tb);
  term2 = term2*log(mstop2sq/mstop1sq);

  term3 = 0.125*gd1*gd1*(a0Peter(mstop1sq,q)+a0Peter(mstop2sq,q));

  term4 = ((lambda*lambda*mtop*mtop)/(sqrt(rt)*tb*tb))*(1.0-4.0*Xt*Xt*mtop*mtop/rt);
  term4 = term4*(a0Peter(mstop2sq,q)-a0Peter(mstop1sq,q));

  double Del33prime = (3.0/(16.0*PI*PI))*(term1-term2+term3+term4);

  return Del33prime;
}

// These helper functions calculate the explicit derivatives
// of the tadpole contributions wrt the input parameters, i.e. returns
// the partial derivatives (1/v_i)*\frac{\partial^2 \Delta V}{\partial a\partial v_i}
// where the VEVs v_i, i = 1, 2, 3, are treated as being fixed.
// Inputs:
//     SoftParsEssm essmSusy = the object to calculate the derivatives on
//     double s = the value of the VEV s to use
//     double tb = the value of tan(beta) to use
double doCalcd2DeltaVdLambdadv1(SoftParsEssm essmSusy, double s, double tb)
{
  double v = essmSusy.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;
  double g1 = essmSusy.displayGaugeCoupling(1);
  double g2 = essmSusy.displayGaugeCoupling(2);
  double gbar = sqrt(g2*g2+0.6*g1*g1);
  double gd1 = essmSusy.displaygdash_1();

  // Logarithms are calculated at the given Q
  double q = essmSusy.displayMu();

  double yt = essmSusy.displayYukawaElement(YU, 3, 3);
  double mtop = yt*v2/sqrt(2.0);
  // Option to fix m_t(M_t) = 165 GeV
  bool usemtofmt = false;
  if(usemtofmt){
    mtop = 165;
    yt = sqrt(2.0)*mtop/v2;
    essmSusy.setYukawaElement(YU, 3, 3, yt);
  }

  DoubleVector lambda_vec = essmSusy.displaylambda();
  double lambda = lambda_vec(3);
  double At = essmSusy.displayA_top();
  double mQlsq = essmSusy.displaySoftMassSquared(mQl,3,3);
  double mUrsq =  essmSusy.displaySoftMassSquared(mUr,3,3);
  double mueff = lambda*s/sqrt(2.0);

  double MQQsq = mQlsq-mUrsq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2);
  double Xt = At - lambda*s/(sqrt(2.0)*tb);
  double rt = MQQsq*MQQsq+4.0*mtop*mtop*Xt*Xt;

  // Squared stop masses are calculated using physical_ESSM to make sure
  // auxiliary D-terms are included.
  DoubleVector mstop(2), mstopsq(2), mD1sq(3), mD2sq(3);
  physical_ESSM(essmSusy, mstop, mstopsq, mD1sq, mD2sq, s, tb);

  // Have to swap the mass ordering of the stops, because
  // physical_ESSM uses the opposite labelling.
  double mstop1sq, mstop2sq;
  mstop1sq = mstopsq(2);
  mstop2sq = mstopsq(1);

  // Calculate the explicit input parameter dependence of the stop masses
  double dmstop1sqdlambda = 0.5*(2.0*sqrt(2.0)*mtop*mtop*Xt*s)/(sqrt(rt)*tb);
  double dmstop2sqdlambda = -0.5*(2.0*sqrt(2.0)*mtop*mtop*Xt*s)/(sqrt(rt)*tb);

  // Calculate the derivative of the stop masses wrt the appropriate VEV
  double dmstop1sqdv1 = doCalcdmstop1sqdv1(essmSusy, s, tb);
  double dmstop2sqdv1 = doCalcdmstop2sqdv1(essmSusy, s, tb);

  // Calculate the explicit input parameter dependence of the derivatives
  // of the stop masses wrt the appropriate VEV
  double coeff = ((4.0*mtop*mtop)/(sqrt(rt)))*(lambda*s/(sqrt(2.0)*v1*v2))*(s/(sqrt(2.0)*tb));
  coeff = coeff-(4.0*mtop*mtop*Xt*s)/(sqrt(2.0*rt)*v1*v2);
  coeff = coeff-(((4.0*mtop*mtop*Xt*s)/(sqrt(2.0*rt)*rt*tb))*
		 (0.25*MQQsq*(g1*g1-g2*g2)+4.0*mtop*mtop*Xt*lambda*s/(sqrt(2.0)*v1*v2)));

  double d2mstop1sqdldv1 = -0.5*v1*coeff;
  double d2mstop2sqdldv1 = 0.5*v1*coeff;

  // Calculate the full second derivative of the tadpole contribution
  double term1 = dmstop1sqdlambda*dmstop1sqdv1*log(mstop1sq/(q*q))+a0Peter(mstop1sq,q)*d2mstop1sqdldv1;
  double term2 = dmstop2sqdlambda*dmstop2sqdv1*log(mstop2sq/(q*q))+a0Peter(mstop2sq,q)*d2mstop2sqdldv1;

  double d2DeltaV = (3.0/(16.0*PI*PI))*(term1+term2);
  d2DeltaV = (1.0/v1)*d2DeltaV;

  return d2DeltaV;
}

double doCalcd2DeltaVdLambdadv2(SoftParsEssm essmSusy, double s, double tb)
{
  double v = essmSusy.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;
  double g1 = essmSusy.displayGaugeCoupling(1);
  double g2 = essmSusy.displayGaugeCoupling(2);
  double gbar = sqrt(g2*g2+0.6*g1*g1);
  double gd1 = essmSusy.displaygdash_1();

  // Logarithms are calculated at the given Q
  double q = essmSusy.displayMu();

  double yt = essmSusy.displayYukawaElement(YU, 3, 3);
  double mtop = yt*v2/sqrt(2.0);
  // Option to fix m_t(M_t) = 165 GeV
  bool usemtofmt = false;
  if(usemtofmt){
    mtop = 165;
    yt = sqrt(2.0)*mtop/v2;
    essmSusy.setYukawaElement(YU, 3, 3, yt);
  }

  DoubleVector lambda_vec = essmSusy.displaylambda();
  double lambda = lambda_vec(3);
  double At = essmSusy.displayA_top();
  double mQlsq = essmSusy.displaySoftMassSquared(mQl,3,3);
  double mUrsq =  essmSusy.displaySoftMassSquared(mUr,3,3);
  double mueff = lambda*s/sqrt(2.0);

  double MQQsq = mQlsq-mUrsq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2);
  double Xt = At - lambda*s/(sqrt(2.0)*tb);
  double rt = MQQsq*MQQsq+4.0*mtop*mtop*Xt*Xt;

  // Squared stop masses are calculated using physical_ESSM to make sure
  // auxiliary D-terms are included.
  DoubleVector mstop(2), mstopsq(2), mD1sq(3), mD2sq(3);
  physical_ESSM(essmSusy, mstop, mstopsq, mD1sq, mD2sq, s, tb);

  // Have to swap the mass ordering of the stops, because
  // physical_ESSM uses the opposite labelling.
  double mstop1sq, mstop2sq;
  mstop1sq = mstopsq(2);
  mstop2sq = mstopsq(1);

  // Calculate the explicit input parameter dependence of the stop masses
  double dmstop1sqdlambda = 0.5*(2.0*sqrt(2.0)*mtop*mtop*Xt*s)/(sqrt(rt)*tb);
  double dmstop2sqdlambda = -0.5*(2.0*sqrt(2.0)*mtop*mtop*Xt*s)/(sqrt(rt)*tb);

  // Calculate the derivative of the stop masses wrt the appropriate VEV
  double dmstop1sqdv2 = doCalcdmstop1sqdv2(essmSusy, s, tb);
  double dmstop2sqdv2 = doCalcdmstop2sqdv2(essmSusy, s, tb);

  // Calculate the explicit input parameter dependence of the derivatives
  // of the stop masses wrt the appropriate VEV
  double coeff = (1.0/sqrt(rt))*(4.0*mtop*mtop*Xt*s*v1/(sqrt(2.0)*v2*v2*v2)
				 -4.0*mtop*mtop*s*lambda*s*v1/(2.0*tb*v2*v2*v2)
				 -4.0*yt*yt*Xt*s/(sqrt(2.0)*tb));
  coeff = coeff-((4.0*mtop*mtop*Xt*s)/(sqrt(2.0*rt)*rt*tb))*
    (0.25*MQQsq*(g2*g2-g1*g1)-2.0*yt*yt*Xt*Xt-4.0*mtop*mtop*Xt*lambda*s*v1/(sqrt(2.0)*v2*v2*v2));

  double d2mstop1sqdldv2 = -0.5*v2*coeff;
  double d2mstop2sqdldv2 = 0.5*v2*coeff;

  // Calculate the full second derivative of the tadpole contribution
  double term1 = dmstop1sqdlambda*dmstop1sqdv2*log(mstop1sq/(q*q))+a0Peter(mstop1sq,q)*d2mstop1sqdldv2;
  double term2 = dmstop2sqdlambda*dmstop2sqdv2*log(mstop2sq/(q*q))+a0Peter(mstop2sq,q)*d2mstop2sqdldv2;

  double d2DeltaV = (3.0/(16.0*PI*PI))*(term1+term2);
  d2DeltaV = (1.0/v2)*d2DeltaV;

  return d2DeltaV;
}

double doCalcd2DeltaVdLambdadv3(SoftParsEssm essmSusy, double s, double tb)
{
  double v = essmSusy.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;
  double g1 = essmSusy.displayGaugeCoupling(1);
  double g2 = essmSusy.displayGaugeCoupling(2);
  double gbar = sqrt(g2*g2+0.6*g1*g1);
  double gd1 = essmSusy.displaygdash_1();

  // Logarithms are calculated at the given Q
  double q = essmSusy.displayMu();

  double yt = essmSusy.displayYukawaElement(YU, 3, 3);
  double mtop = yt*v2/sqrt(2.0);
  // Option to fix m_t(M_t) = 165 GeV
  bool usemtofmt = false;
  if(usemtofmt){
    mtop = 165;
    yt = sqrt(2.0)*mtop/v2;
    essmSusy.setYukawaElement(YU, 3, 3, yt);
  }

  DoubleVector lambda_vec = essmSusy.displaylambda();
  double lambda = lambda_vec(3);
  double At = essmSusy.displayA_top();
  double mQlsq = essmSusy.displaySoftMassSquared(mQl,3,3);
  double mUrsq =  essmSusy.displaySoftMassSquared(mUr,3,3);
  double mueff = lambda*s/sqrt(2.0);

  double MQQsq = mQlsq-mUrsq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2);
  double Xt = At - lambda*s/(sqrt(2.0)*tb);
  double rt = MQQsq*MQQsq+4.0*mtop*mtop*Xt*Xt;

  // Squared stop masses are calculated using physical_ESSM to make sure
  // auxiliary D-terms are included.
  DoubleVector mstop(2), mstopsq(2), mD1sq(3), mD2sq(3);
  physical_ESSM(essmSusy, mstop, mstopsq, mD1sq, mD2sq, s, tb);

  // Have to swap the mass ordering of the stops, because
  // physical_ESSM uses the opposite labelling.
  double mstop1sq, mstop2sq;
  mstop1sq = mstopsq(2);
  mstop2sq = mstopsq(1);

  // Calculate the explicit input parameter dependence of the stop masses
  double dmstop1sqdlambda = 0.5*(2.0*sqrt(2.0)*mtop*mtop*Xt*s)/(sqrt(rt)*tb);
  double dmstop2sqdlambda = -0.5*(2.0*sqrt(2.0)*mtop*mtop*Xt*s)/(sqrt(rt)*tb);

  // Calculate the derivative of the stop masses wrt the appropriate VEV
  double dmstop1sqdv3 = doCalcdmstop1sqdv3(essmSusy, s, tb);
  double dmstop2sqdv3 = doCalcdmstop2sqdv3(essmSusy, s, tb);

  // Calculate the explicit input parameter dependence of the derivatives
  // of the stop masses wrt the appropriate VEV
  double coeff = s*lambda*v1*mtop*mtop/(2.0*sqrt(rt)*tb*v2*s);
  coeff = coeff-mtop*mtop*Xt*v1/(sqrt(2.0*rt)*v2*s);
  coeff = coeff-4.0*sqr(mtop*mtop)*Xt*Xt*s*lambda*v1/(2.0*sqrt(rt)*rt*tb*v2*s);

  double d2mstop1sqdldv3 = -2.0*s*coeff;
  double d2mstop2sqdldv3 = 2.0*s*coeff;

  // Calculate the full second derivative of the tadpole contribution
  double term1 = dmstop1sqdlambda*dmstop1sqdv3*log(mstop1sq/(q*q))+a0Peter(mstop1sq,q)*d2mstop1sqdldv3;
  double term2 = dmstop2sqdlambda*dmstop2sqdv3*log(mstop2sq/(q*q))+a0Peter(mstop2sq,q)*d2mstop2sqdldv3;

  double d2DeltaV = (3.0/(16.0*PI*PI))*(term1+term2);
  d2DeltaV = (1.0/s)*d2DeltaV;

  return d2DeltaV;
}

double doCalcd2DeltaVdAtdv1(SoftParsEssm essmSusy, double s, double tb)
{
  double v = essmSusy.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;
  double g1 = essmSusy.displayGaugeCoupling(1);
  double g2 = essmSusy.displayGaugeCoupling(2);
  double gbar = sqrt(g2*g2+0.6*g1*g1);
  double gd1 = essmSusy.displaygdash_1();

  // Logarithms are calculated at the given Q
  double q = essmSusy.displayMu();

  double yt = essmSusy.displayYukawaElement(YU, 3, 3);
  double mtop = yt*v2/sqrt(2.0);
  // Option to fix m_t(M_t) = 165 GeV
  bool usemtofmt = false;
  if(usemtofmt){
    mtop = 165;
    yt = sqrt(2.0)*mtop/v2;
    essmSusy.setYukawaElement(YU, 3, 3, yt);
  }

  DoubleVector lambda_vec = essmSusy.displaylambda();
  double lambda = lambda_vec(3);
  double At = essmSusy.displayA_top();
  double mQlsq = essmSusy.displaySoftMassSquared(mQl,3,3);
  double mUrsq =  essmSusy.displaySoftMassSquared(mUr,3,3);
  double mueff = lambda*s/sqrt(2.0);

  double MQQsq = mQlsq-mUrsq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2);
  double Xt = At - lambda*s/(sqrt(2.0)*tb);
  double rt = MQQsq*MQQsq+4.0*mtop*mtop*Xt*Xt;

  // Squared stop masses are calculated using physical_ESSM to make sure
  // auxiliary D-terms are included.
  DoubleVector mstop(2), mstopsq(2), mD1sq(3), mD2sq(3);
  physical_ESSM(essmSusy, mstop, mstopsq, mD1sq, mD2sq, s, tb);

  // Have to swap the mass ordering of the stops, because
  // physical_ESSM uses the opposite labelling.
  double mstop1sq, mstop2sq;
  mstop1sq = mstopsq(2);
  mstop2sq = mstopsq(1);

  // Calculate the explicit input parameter dependence of the stop masses
  double dmstop1sqdAt = -0.5*(4.0*mtop*mtop*Xt)/sqrt(rt);
  double dmstop2sqdAt = 0.5*(4.0*mtop*mtop*Xt)/sqrt(rt);

  // Calculate the derivative of the stop masses wrt the appropriate VEV
  double dmstop1sqdv1 = doCalcdmstop1sqdv1(essmSusy, s, tb);
  double dmstop2sqdv1 = doCalcdmstop2sqdv1(essmSusy, s, tb);

  // Calculate the explicit input parameter dependence of the derivatives
  // of the stop masses wrt the appropriate VEV
  double coeff = (4.0*mtop*mtop*Xt/(sqrt(rt)*rt))*(0.25*MQQsq*(g1*g1-g2*g2)
						   +4.0*mtop*mtop*Xt*lambda*s/(sqrt(2.0)*v1*v2));
  coeff = coeff-(4.0*mtop*mtop*lambda*s/(sqrt(2.0*rt)*v1*v2));

  double d2mstop1sqdAtdv1 = -0.5*v1*coeff;
  double d2mstop2sqdAtdv1 = 0.5*v1*coeff;

  // Calculate the full second derivative of the tadpole contribution
  double term1 = dmstop1sqdAt*dmstop1sqdv1*log(mstop1sq/(q*q))+a0Peter(mstop1sq,q)*d2mstop1sqdAtdv1;
  double term2 = dmstop2sqdAt*dmstop2sqdv1*log(mstop2sq/(q*q))+a0Peter(mstop2sq,q)*d2mstop2sqdAtdv1;

  double d2DeltaV = (3.0/(16.0*PI*PI))*(term1+term2);
  d2DeltaV = (1.0/v1)*d2DeltaV;

  return d2DeltaV;
}

double doCalcd2DeltaVdAtdv2(SoftParsEssm essmSusy, double s, double tb)
{
  double v = essmSusy.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;
  double g1 = essmSusy.displayGaugeCoupling(1);
  double g2 = essmSusy.displayGaugeCoupling(2);
  double gbar = sqrt(g2*g2+0.6*g1*g1);
  double gd1 = essmSusy.displaygdash_1();

  // Logarithms are calculated at the given Q
  double q = essmSusy.displayMu();

  double yt = essmSusy.displayYukawaElement(YU, 3, 3);
  double mtop = yt*v2/sqrt(2.0);
  // Option to fix m_t(M_t) = 165 GeV
  bool usemtofmt = false;
  if(usemtofmt){
    mtop = 165;
    yt = sqrt(2.0)*mtop/v2;
    essmSusy.setYukawaElement(YU, 3, 3, yt);
  }

  DoubleVector lambda_vec = essmSusy.displaylambda();
  double lambda = lambda_vec(3);
  double At = essmSusy.displayA_top();
  double mQlsq = essmSusy.displaySoftMassSquared(mQl,3,3);
  double mUrsq =  essmSusy.displaySoftMassSquared(mUr,3,3);
  double mueff = lambda*s/sqrt(2.0);

  double MQQsq = mQlsq-mUrsq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2);
  double Xt = At - lambda*s/(sqrt(2.0)*tb);
  double rt = MQQsq*MQQsq+4.0*mtop*mtop*Xt*Xt;

  // Squared stop masses are calculated using physical_ESSM to make sure
  // auxiliary D-terms are included.
  DoubleVector mstop(2), mstopsq(2), mD1sq(3), mD2sq(3);
  physical_ESSM(essmSusy, mstop, mstopsq, mD1sq, mD2sq, s, tb);

  // Have to swap the mass ordering of the stops, because
  // physical_ESSM uses the opposite labelling.
  double mstop1sq, mstop2sq;
  mstop1sq = mstopsq(2);
  mstop2sq = mstopsq(1);

  // Calculate the explicit input parameter dependence of the stop masses
  double dmstop1sqdAt = -0.5*(4.0*mtop*mtop*Xt)/sqrt(rt);
  double dmstop2sqdAt = 0.5*(4.0*mtop*mtop*Xt)/sqrt(rt);

  // Calculate the derivative of the stop masses wrt the appropriate VEV
  double dmstop1sqdv2 = doCalcdmstop1sqdv2(essmSusy, s, tb);
  double dmstop2sqdv2 = doCalcdmstop2sqdv2(essmSusy, s, tb);

  // Calculate the explicit input parameter dependence of the derivatives
  // of the stop masses wrt the appropriate VEV
  double coeff = (1.0/sqrt(rt))*(4.0*yt*yt*Xt+4.0*mtop*mtop*lambda*s*v1/(sqrt(2.0)*v2*v2*v2));
  coeff = coeff+((4.0*mtop*mtop*Xt)/(rt*sqrt(rt)))*(0.25*MQQsq*(g2*g2-g1*g1)-2.0*yt*yt*Xt*Xt
						    -4.0*mtop*mtop*Xt*lambda*s*v1/(sqrt(2.0)*v2*v2*v2));

  double d2mstop1sqdAtdv2 = -0.5*v2*coeff;
  double d2mstop2sqdAtdv2 = 0.5*v2*coeff;

  // Calculate the full second derivative of the tadpole contribution
  double term1 = dmstop1sqdAt*dmstop1sqdv2*log(mstop1sq/(q*q))+a0Peter(mstop1sq,q)*d2mstop1sqdAtdv2;
  double term2 = dmstop2sqdAt*dmstop2sqdv2*log(mstop2sq/(q*q))+a0Peter(mstop2sq,q)*d2mstop2sqdAtdv2;

  double d2DeltaV = (3.0/(16.0*PI*PI))*(term1+term2);
  d2DeltaV = (1.0/v2)*d2DeltaV;

  return d2DeltaV;

}

double doCalcd2DeltaVdAtdv3(SoftParsEssm essmSusy, double s, double tb)
{
  double v = essmSusy.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;
  double g1 = essmSusy.displayGaugeCoupling(1);
  double g2 = essmSusy.displayGaugeCoupling(2);
  double gbar = sqrt(g2*g2+0.6*g1*g1);
  double gd1 = essmSusy.displaygdash_1();

  // Logarithms are calculated at the given Q
  double q = essmSusy.displayMu();

  double yt = essmSusy.displayYukawaElement(YU, 3, 3);
  double mtop = yt*v2/sqrt(2.0);
  // Option to fix m_t(M_t) = 165 GeV
  bool usemtofmt = false;
  if(usemtofmt){
    mtop = 165;
    yt = sqrt(2.0)*mtop/v2;
    essmSusy.setYukawaElement(YU, 3, 3, yt);
  }

  DoubleVector lambda_vec = essmSusy.displaylambda();
  double lambda = lambda_vec(3);
  double At = essmSusy.displayA_top();
  double mQlsq = essmSusy.displaySoftMassSquared(mQl,3,3);
  double mUrsq =  essmSusy.displaySoftMassSquared(mUr,3,3);
  double mueff = lambda*s/sqrt(2.0);

  double MQQsq = mQlsq-mUrsq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2);
  double Xt = At - lambda*s/(sqrt(2.0)*tb);
  double rt = MQQsq*MQQsq+4.0*mtop*mtop*Xt*Xt;

  // Squared stop masses are calculated using physical_ESSM to make sure
  // auxiliary D-terms are included.
  DoubleVector mstop(2), mstopsq(2), mD1sq(3), mD2sq(3);
  physical_ESSM(essmSusy, mstop, mstopsq, mD1sq, mD2sq, s, tb);

  // Have to swap the mass ordering of the stops, because
  // physical_ESSM uses the opposite labelling.
  double mstop1sq, mstop2sq;
  mstop1sq = mstopsq(2);
  mstop2sq = mstopsq(1);

  // Calculate the explicit input parameter dependence of the stop masses
  double dmstop1sqdAt = -0.5*(4.0*mtop*mtop*Xt)/sqrt(rt);
  double dmstop2sqdAt = 0.5*(4.0*mtop*mtop*Xt)/sqrt(rt);

  // Calculate the derivative of the stop masses wrt the appropriate VEV
  double dmstop1sqdv3 = doCalcdmstop1sqdv3(essmSusy, s, tb);
  double dmstop2sqdv3 = doCalcdmstop2sqdv3(essmSusy, s, tb);

  // Calculate the explicit input parameter dependence of the derivatives
  // of the stop masses wrt the appropriate VEV
  double coeff = 4.0*sqr(mtop*mtop)*Xt*Xt*lambda*v1/(sqrt(2.0*rt)*rt*v2*s);
  coeff = coeff-mtop*mtop*lambda*v1/(sqrt(2.0*rt)*v2*s);

  double d2mstop1sqdAtdv3 = -2.0*s*coeff;
  double d2mstop2sqdAtdv3 = 2.0*s*coeff;

  // Calculate the full second derivative of the tadpole contribution
  double term1 = dmstop1sqdAt*dmstop1sqdv3*log(mstop1sq/(q*q))+a0Peter(mstop1sq,q)*d2mstop1sqdAtdv3;
  double term2 = dmstop2sqdAt*dmstop2sqdv3*log(mstop2sq/(q*q))+a0Peter(mstop2sq,q)*d2mstop2sqdAtdv3;

  double d2DeltaV = (3.0/(16.0*PI*PI))*(term1+term2);
  d2DeltaV = (1.0/s)*d2DeltaV;

  return d2DeltaV;

}

double doCalcd2DeltaVdmQlsqdv1(SoftParsEssm essmSusy, double s, double tb)
{
  double v = essmSusy.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;
  double g1 = essmSusy.displayGaugeCoupling(1);
  double g2 = essmSusy.displayGaugeCoupling(2);
  double gbar = sqrt(g2*g2+0.6*g1*g1);
  double gd1 = essmSusy.displaygdash_1();

  // Logarithms are calculated at the given Q
  double q = essmSusy.displayMu();

  double yt = essmSusy.displayYukawaElement(YU, 3, 3);
  double mtop = yt*v2/sqrt(2.0);
  // Option to fix m_t(M_t) = 165 GeV
  bool usemtofmt = false;
  if(usemtofmt){
    mtop = 165;
    yt = sqrt(2.0)*mtop/v2;
    essmSusy.setYukawaElement(YU, 3, 3, yt);
  }

  DoubleVector lambda_vec = essmSusy.displaylambda();
  double lambda = lambda_vec(3);
  double At = essmSusy.displayA_top();
  double mQlsq = essmSusy.displaySoftMassSquared(mQl,3,3);
  double mUrsq =  essmSusy.displaySoftMassSquared(mUr,3,3);
  double mueff = lambda*s/sqrt(2.0);

  double MQQsq = mQlsq-mUrsq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2);
  double Xt = At - lambda*s/(sqrt(2.0)*tb);
  double rt = MQQsq*MQQsq+4.0*mtop*mtop*Xt*Xt;

  // Squared stop masses are calculated using physical_ESSM to make sure
  // auxiliary D-terms are included.
  DoubleVector mstop(2), mstopsq(2), mD1sq(3), mD2sq(3);
  physical_ESSM(essmSusy, mstop, mstopsq, mD1sq, mD2sq, s, tb);

  // Have to swap the mass ordering of the stops, because
  // physical_ESSM uses the opposite labelling.
  double mstop1sq, mstop2sq;
  mstop1sq = mstopsq(2);
  mstop2sq = mstopsq(1);

  // Calculate the explicit input parameter dependence of the stop masses
  double dmstop1sqdmQlsq = 0.5*(1.0-MQQsq/sqrt(rt));
  double dmstop2sqdmQlsq =0.5*(1.0+MQQsq/sqrt(rt));

  // Calculate the derivative of the stop masses wrt the appropriate VEV
  double dmstop1sqdv1 = doCalcdmstop1sqdv1(essmSusy, s, tb);
  double dmstop2sqdv1 = doCalcdmstop2sqdv1(essmSusy, s, tb);

  // Calculate the explicit input parameter dependence of the derivatives
  // of the stop masses wrt the appropriate VEV
  double coeff = 0.25*(g2*g2-g1*g1)/sqrt(rt);
  coeff = coeff+(MQQsq/(rt*sqrt(rt)))*(0.25*MQQsq*(g1*g1-g2*g2)
				       +4.0*mtop*mtop*Xt*lambda*s/(sqrt(2.0)*v1*v2));

  double d2mstop1sqdmQlsqdv1 = -0.5*v1*coeff;
  double d2mstop2sqdmQlsqdv1 = 0.5*v1*coeff;

  // Calculate the full second derivative of the tadpole contribution
  double term1 = dmstop1sqdmQlsq*dmstop1sqdv1*log(mstop1sq/(q*q))+a0Peter(mstop1sq,q)*d2mstop1sqdmQlsqdv1;
  double term2 = dmstop2sqdmQlsq*dmstop2sqdv1*log(mstop2sq/(q*q))+a0Peter(mstop2sq,q)*d2mstop2sqdmQlsqdv1;

  double d2DeltaV = (3.0/(16.0*PI*PI))*(term1+term2);
  d2DeltaV = (1.0/v1)*d2DeltaV;

  return d2DeltaV;

}

double doCalcd2DeltaVdmQlsqdv2(SoftParsEssm essmSusy, double s, double tb)
{
  double v = essmSusy.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;
  double g1 = essmSusy.displayGaugeCoupling(1);
  double g2 = essmSusy.displayGaugeCoupling(2);
  double gbar = sqrt(g2*g2+0.6*g1*g1);
  double gd1 = essmSusy.displaygdash_1();

  // Logarithms are calculated at the given Q
  double q = essmSusy.displayMu();

  double yt = essmSusy.displayYukawaElement(YU, 3, 3);
  double mtop = yt*v2/sqrt(2.0);
  // Option to fix m_t(M_t) = 165 GeV
  bool usemtofmt = false;
  if(usemtofmt){
    mtop = 165;
    yt = sqrt(2.0)*mtop/v2;
    essmSusy.setYukawaElement(YU, 3, 3, yt);
  }

  DoubleVector lambda_vec = essmSusy.displaylambda();
  double lambda = lambda_vec(3);
  double At = essmSusy.displayA_top();
  double mQlsq = essmSusy.displaySoftMassSquared(mQl,3,3);
  double mUrsq =  essmSusy.displaySoftMassSquared(mUr,3,3);
  double mueff = lambda*s/sqrt(2.0);

  double MQQsq = mQlsq-mUrsq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2);
  double Xt = At - lambda*s/(sqrt(2.0)*tb);
  double rt = MQQsq*MQQsq+4.0*mtop*mtop*Xt*Xt;

  // Squared stop masses are calculated using physical_ESSM to make sure
  // auxiliary D-terms are included.
  DoubleVector mstop(2), mstopsq(2), mD1sq(3), mD2sq(3);
  physical_ESSM(essmSusy, mstop, mstopsq, mD1sq, mD2sq, s, tb);

  // Have to swap the mass ordering of the stops, because
  // physical_ESSM uses the opposite labelling.
  double mstop1sq, mstop2sq;
  mstop1sq = mstopsq(2);
  mstop2sq = mstopsq(1);

  // Calculate the explicit input parameter dependence of the stop masses
  double dmstop1sqdmQlsq = 0.5*(1.0-MQQsq/sqrt(rt));
  double dmstop2sqdmQlsq =0.5*(1.0+MQQsq/sqrt(rt));

  // Calculate the derivative of the stop masses wrt the appropriate VEV
  double dmstop1sqdv2 = doCalcdmstop1sqdv2(essmSusy, s, tb);
  double dmstop2sqdv2 = doCalcdmstop2sqdv2(essmSusy, s, tb);

  // Calculate the explicit input parameter dependence of the derivatives
  // of the stop masses wrt the appropriate VEV
  double coeff = 0.25*(g1*g1-g2*g2)/sqrt(rt);
  coeff = coeff+(MQQsq/(rt*sqrt(rt)))*(0.25*MQQsq*(g2*g2-g1*g1)-2.0*yt*yt*Xt*Xt
				       -4.0*mtop*mtop*Xt*lambda*s*v1/(sqrt(2.0)*v2*v2*v2));

  double d2mstop1sqdmQlsqdv2 = -0.5*v2*coeff;
  double d2mstop2sqdmQlsqdv2 = 0.5*v2*coeff;

  // Calculate the full second derivative of the tadpole contribution
  double term1 = dmstop1sqdmQlsq*dmstop1sqdv2*log(mstop1sq/(q*q))+a0Peter(mstop1sq,q)*d2mstop1sqdmQlsqdv2;
  double term2 = dmstop2sqdmQlsq*dmstop2sqdv2*log(mstop2sq/(q*q))+a0Peter(mstop2sq,q)*d2mstop2sqdmQlsqdv2;

  double d2DeltaV = (3.0/(16.0*PI*PI))*(term1+term2);
  d2DeltaV = (1.0/v2)*d2DeltaV;

  return d2DeltaV;

}

double doCalcd2DeltaVdmQlsqdv3(SoftParsEssm essmSusy, double s, double tb)
{
  double v = essmSusy.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;
  double g1 = essmSusy.displayGaugeCoupling(1);
  double g2 = essmSusy.displayGaugeCoupling(2);
  double gbar = sqrt(g2*g2+0.6*g1*g1);
  double gd1 = essmSusy.displaygdash_1();

  // Logarithms are calculated at the given Q
  double q = essmSusy.displayMu();

  double yt = essmSusy.displayYukawaElement(YU, 3, 3);
  double mtop = yt*v2/sqrt(2.0);
  // Option to fix m_t(M_t) = 165 GeV
  bool usemtofmt = false;
  if(usemtofmt){
    mtop = 165;
    yt = sqrt(2.0)*mtop/v2;
    essmSusy.setYukawaElement(YU, 3, 3, yt);
  }

  DoubleVector lambda_vec = essmSusy.displaylambda();
  double lambda = lambda_vec(3);
  double At = essmSusy.displayA_top();
  double mQlsq = essmSusy.displaySoftMassSquared(mQl,3,3);
  double mUrsq =  essmSusy.displaySoftMassSquared(mUr,3,3);
  double mueff = lambda*s/sqrt(2.0);

  double MQQsq = mQlsq-mUrsq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2);
  double Xt = At - lambda*s/(sqrt(2.0)*tb);
  double rt = MQQsq*MQQsq+4.0*mtop*mtop*Xt*Xt;

  // Squared stop masses are calculated using physical_ESSM to make sure
  // auxiliary D-terms are included.
  DoubleVector mstop(2), mstopsq(2), mD1sq(3), mD2sq(3);
  physical_ESSM(essmSusy, mstop, mstopsq, mD1sq, mD2sq, s, tb);

  // Have to swap the mass ordering of the stops, because
  // physical_ESSM uses the opposite labelling.
  double mstop1sq, mstop2sq;
  mstop1sq = mstopsq(2);
  mstop2sq = mstopsq(1);

  // Calculate the explicit input parameter dependence of the stop masses
  double dmstop1sqdmQlsq = 0.5*(1.0-MQQsq/sqrt(rt));
  double dmstop2sqdmQlsq =0.5*(1.0+MQQsq/sqrt(rt));

  // Calculate the derivative of the stop masses wrt the appropriate VEV
  double dmstop1sqdv3 = doCalcdmstop1sqdv3(essmSusy, s, tb);
  double dmstop2sqdv3 = doCalcdmstop2sqdv3(essmSusy, s, tb);

  // Calculate the explicit input parameter dependence of the derivatives
  // of the stop masses wrt the appropriate VEV
  double coeff = mtop*mtop*Xt*lambda*v1*MQQsq/(sqrt(2.0*rt)*rt*v2*s);

  double d2mstop1sqdmQlsqdv3 = -2.0*s*coeff;
  double d2mstop2sqdmQlsqdv3 = 2.0*s*coeff;

  // Calculate the full second derivative of the tadpole contribution
  double term1 = dmstop1sqdmQlsq*dmstop1sqdv3*log(mstop1sq/(q*q))+a0Peter(mstop1sq,q)*d2mstop1sqdmQlsqdv3;
  double term2 = dmstop2sqdmQlsq*dmstop2sqdv3*log(mstop2sq/(q*q))+a0Peter(mstop2sq,q)*d2mstop2sqdmQlsqdv3;

  double d2DeltaV = (3.0/(16.0*PI*PI))*(term1+term2);
  d2DeltaV = (1.0/s)*d2DeltaV;

  return d2DeltaV;

}

double doCalcd2DeltaVdmUrsqdv1(SoftParsEssm essmSusy, double s, double tb)
{
  double v = essmSusy.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;
  double g1 = essmSusy.displayGaugeCoupling(1);
  double g2 = essmSusy.displayGaugeCoupling(2);
  double gbar = sqrt(g2*g2+0.6*g1*g1);
  double gd1 = essmSusy.displaygdash_1();

  // Logarithms are calculated at the given Q
  double q = essmSusy.displayMu();

  double yt = essmSusy.displayYukawaElement(YU, 3, 3);
  double mtop = yt*v2/sqrt(2.0);
  // Option to fix m_t(M_t) = 165 GeV
  bool usemtofmt = false;
  if(usemtofmt){
    mtop = 165;
    yt = sqrt(2.0)*mtop/v2;
    essmSusy.setYukawaElement(YU, 3, 3, yt);
  }

  DoubleVector lambda_vec = essmSusy.displaylambda();
  double lambda = lambda_vec(3);
  double At = essmSusy.displayA_top();
  double mQlsq = essmSusy.displaySoftMassSquared(mQl,3,3);
  double mUrsq =  essmSusy.displaySoftMassSquared(mUr,3,3);
  double mueff = lambda*s/sqrt(2.0);

  double MQQsq = mQlsq-mUrsq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2);
  double Xt = At - lambda*s/(sqrt(2.0)*tb);
  double rt = MQQsq*MQQsq+4.0*mtop*mtop*Xt*Xt;

  // Squared stop masses are calculated using physical_ESSM to make sure
  // auxiliary D-terms are included.
  DoubleVector mstop(2), mstopsq(2), mD1sq(3), mD2sq(3);
  physical_ESSM(essmSusy, mstop, mstopsq, mD1sq, mD2sq, s, tb);

  // Have to swap the mass ordering of the stops, because
  // physical_ESSM uses the opposite labelling.
  double mstop1sq, mstop2sq;
  mstop1sq = mstopsq(2);
  mstop2sq = mstopsq(1);

  // Calculate the explicit input parameter dependence of the stop masses
  double dmstop1sqdmUrsq = 0.5*(1.0+MQQsq/sqrt(rt));
  double dmstop2sqdmUrsq =0.5*(1.0-MQQsq/sqrt(rt));

  // Calculate the derivative of the stop masses wrt the appropriate VEV
  double dmstop1sqdv1 = doCalcdmstop1sqdv1(essmSusy, s, tb);
  double dmstop2sqdv1 = doCalcdmstop2sqdv1(essmSusy, s, tb);

  // Calculate the explicit input parameter dependence of the derivatives
  // of the stop masses wrt the appropriate VEV
  double coeff = 0.25*(g2*g2-g1*g1)/sqrt(rt);
  coeff = coeff+(MQQsq/(rt*sqrt(rt)))*(0.25*MQQsq*(g1*g1-g2*g2)
				       +4.0*mtop*mtop*Xt*lambda*s/(sqrt(2.0)*v1*v2));

  double d2mstop1sqdmUrsqdv1 = 0.5*v1*coeff;
  double d2mstop2sqdmUrsqdv1 = -0.5*v1*coeff;

  // Calculate the full second derivative of the tadpole contribution
  double term1 = dmstop1sqdmUrsq*dmstop1sqdv1*log(mstop1sq/(q*q))+a0Peter(mstop1sq,q)*d2mstop1sqdmUrsqdv1;
  double term2 = dmstop2sqdmUrsq*dmstop2sqdv1*log(mstop2sq/(q*q))+a0Peter(mstop2sq,q)*d2mstop2sqdmUrsqdv1;

  double d2DeltaV = (3.0/(16.0*PI*PI))*(term1+term2);
  d2DeltaV = (1.0/v1)*d2DeltaV;

  return d2DeltaV;

}

double doCalcd2DeltaVdmUrsqdv2(SoftParsEssm essmSusy, double s, double tb)
{
  double v = essmSusy.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;
  double g1 = essmSusy.displayGaugeCoupling(1);
  double g2 = essmSusy.displayGaugeCoupling(2);
  double gbar = sqrt(g2*g2+0.6*g1*g1);
  double gd1 = essmSusy.displaygdash_1();

  // Logarithms are calculated at the given Q
  double q = essmSusy.displayMu();

  double yt = essmSusy.displayYukawaElement(YU, 3, 3);
  double mtop = yt*v2/sqrt(2.0);
  // Option to fix m_t(M_t) = 165 GeV
  bool usemtofmt = false;
  if(usemtofmt){
    mtop = 165;
    yt = sqrt(2.0)*mtop/v2;
    essmSusy.setYukawaElement(YU, 3, 3, yt);
  }

  DoubleVector lambda_vec = essmSusy.displaylambda();
  double lambda = lambda_vec(3);
  double At = essmSusy.displayA_top();
  double mQlsq = essmSusy.displaySoftMassSquared(mQl,3,3);
  double mUrsq =  essmSusy.displaySoftMassSquared(mUr,3,3);
  double mueff = lambda*s/sqrt(2.0);

  double MQQsq = mQlsq-mUrsq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2);
  double Xt = At - lambda*s/(sqrt(2.0)*tb);
  double rt = MQQsq*MQQsq+4.0*mtop*mtop*Xt*Xt;

  // Squared stop masses are calculated using physical_ESSM to make sure
  // auxiliary D-terms are included.
  DoubleVector mstop(2), mstopsq(2), mD1sq(3), mD2sq(3);
  physical_ESSM(essmSusy, mstop, mstopsq, mD1sq, mD2sq, s, tb);

  // Have to swap the mass ordering of the stops, because
  // physical_ESSM uses the opposite labelling.
  double mstop1sq, mstop2sq;
  mstop1sq = mstopsq(2);
  mstop2sq = mstopsq(1);

  // Calculate the explicit input parameter dependence of the stop masses
  double dmstop1sqdmUrsq = 0.5*(1.0+MQQsq/sqrt(rt));
  double dmstop2sqdmUrsq =0.5*(1.0-MQQsq/sqrt(rt));

  // Calculate the derivative of the stop masses wrt the appropriate VEV
  double dmstop1sqdv2 = doCalcdmstop1sqdv2(essmSusy, s, tb);
  double dmstop2sqdv2 = doCalcdmstop2sqdv2(essmSusy, s, tb);

  // Calculate the explicit input parameter dependence of the derivatives
  // of the stop masses wrt the appropriate VEV
  double coeff = 0.25*(g1*g1-g2*g2)/sqrt(rt);
  coeff = coeff+(MQQsq/(rt*sqrt(rt)))*(0.25*MQQsq*(g2*g2-g1*g1)-2.0*yt*yt*Xt*Xt
				       -4.0*mtop*mtop*Xt*lambda*s*v1/(sqrt(2.0)*v2*v2*v2));

  double d2mstop1sqdmUrsqdv2 = 0.5*v2*coeff;
  double d2mstop2sqdmUrsqdv2 = -0.5*v2*coeff;

  // Calculate the full second derivative of the tadpole contribution
  double term1 = dmstop1sqdmUrsq*dmstop1sqdv2*log(mstop1sq/(q*q))+a0Peter(mstop1sq,q)*d2mstop1sqdmUrsqdv2;
  double term2 = dmstop2sqdmUrsq*dmstop2sqdv2*log(mstop2sq/(q*q))+a0Peter(mstop2sq,q)*d2mstop2sqdmUrsqdv2;

  double d2DeltaV = (3.0/(16.0*PI*PI))*(term1+term2);
  d2DeltaV = (1.0/v2)*d2DeltaV;

  return d2DeltaV;

}

double doCalcd2DeltaVdmUrsqdv3(SoftParsEssm essmSusy, double s, double tb)
{
  double v = essmSusy.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;
  double g1 = essmSusy.displayGaugeCoupling(1);
  double g2 = essmSusy.displayGaugeCoupling(2);
  double gbar = sqrt(g2*g2+0.6*g1*g1);
  double gd1 = essmSusy.displaygdash_1();

  // Logarithms are calculated at the given Q
  double q = essmSusy.displayMu();

  double yt = essmSusy.displayYukawaElement(YU, 3, 3);
  double mtop = yt*v2/sqrt(2.0);
  // Option to fix m_t(M_t) = 165 GeV
  bool usemtofmt = false;
  if(usemtofmt){
    mtop = 165;
    yt = sqrt(2.0)*mtop/v2;
    essmSusy.setYukawaElement(YU, 3, 3, yt);
  }

  DoubleVector lambda_vec = essmSusy.displaylambda();
  double lambda = lambda_vec(3);
  double At = essmSusy.displayA_top();
  double mQlsq = essmSusy.displaySoftMassSquared(mQl,3,3);
  double mUrsq =  essmSusy.displaySoftMassSquared(mUr,3,3);
  double mueff = lambda*s/sqrt(2.0);

  double MQQsq = mQlsq-mUrsq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2);
  double Xt = At - lambda*s/(sqrt(2.0)*tb);
  double rt = MQQsq*MQQsq+4.0*mtop*mtop*Xt*Xt;

  // Squared stop masses are calculated using physical_ESSM to make sure
  // auxiliary D-terms are included.
  DoubleVector mstop(2), mstopsq(2), mD1sq(3), mD2sq(3);
  physical_ESSM(essmSusy, mstop, mstopsq, mD1sq, mD2sq, s, tb);

  // Have to swap the mass ordering of the stops, because
  // physical_ESSM uses the opposite labelling.
  double mstop1sq, mstop2sq;
  mstop1sq = mstopsq(2);
  mstop2sq = mstopsq(1);

  // Calculate the explicit input parameter dependence of the stop masses
  double dmstop1sqdmUrsq = 0.5*(1.0+MQQsq/sqrt(rt));
  double dmstop2sqdmUrsq =0.5*(1.0-MQQsq/sqrt(rt));

  // Calculate the derivative of the stop masses wrt the appropriate VEV
  double dmstop1sqdv3 = doCalcdmstop1sqdv3(essmSusy, s, tb);
  double dmstop2sqdv3 = doCalcdmstop2sqdv3(essmSusy, s, tb);

  // Calculate the explicit input parameter dependence of the derivatives
  // of the stop masses wrt the appropriate VEV
  double coeff = mtop*mtop*Xt*lambda*v1*MQQsq/(sqrt(2.0*rt)*rt*v2*s);

  double d2mstop1sqdmUrsqdv3 = 2.0*s*coeff;
  double d2mstop2sqdmUrsqdv3 = -2.0*s*coeff;

  // Calculate the full second derivative of the tadpole contribution
  double term1 = dmstop1sqdmUrsq*dmstop1sqdv3*log(mstop1sq/(q*q))+a0Peter(mstop1sq,q)*d2mstop1sqdmUrsqdv3;
  double term2 = dmstop2sqdmUrsq*dmstop2sqdv3*log(mstop2sq/(q*q))+a0Peter(mstop2sq,q)*d2mstop2sqdmUrsqdv3;

  double d2DeltaV = (3.0/(16.0*PI*PI))*(term1+term2);
  d2DeltaV = (1.0/s)*d2DeltaV;

  return d2DeltaV;

}

/*
  --------------------------------------------------------------------------------------
  The following are methods for evaluating the derivatives of the VEVs and the tuning
  measures numerically. Wherever possible DO NOT use them, as to work they require
  turning off the use m_t(M_t) option in the Higgs codes, which may lead to bugs
  and inconsistencies. If they must be used, note also that the EWSB conditions should
  be changed to those in which the VEVs are divided out (assumes realisitic vacua in which
  all VEVs are non-vanishing). This improves the convergence of Newton's method.
  --------------------------------------------------------------------------------------
 */


// EWSBJacobian_TreeLevel and _OneLoop compute the Jacobian matrix of
// the system defined by the three EWSB conditions, with elements
// J_{ij}= \partial f_i/\partial v_j (v_3 = s). The functions f_i
// are the EWSB conditions, with f_i=\partial V/\partial v_i. The
// two versions differ in the effective potential V used; one loop
// includes the leading one loop contributions from stop and top loops.
// Inputs:
//     SoftParsEssm essmSusy = the object to calculate the Jacobian for
//     double s = the value of the VEV s to use
//     double tb = the value of tan(beta) to use 
DoubleMatrix EWSBJacobian_TreeLevel(SoftParsEssm essmSusy, double s, double tb)
{
  DoubleMatrix J(3,3);

  double v = essmSusy.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;

  double g1 = essmSusy.displayGaugeCoupling(1);
  double g2 = essmSusy.displayGaugeCoupling(2);
  double gbar = sqrt(g2*g2+0.6*g1*g1);
  double gd1 = essmSusy.displaygdash_1();

  DoubleVector Alambda_vec = essmSusy.displayA_lambda();
  double Alambda = Alambda_vec.display(3);
  DoubleVector lambda_vec = essmSusy.displaylambda();
  double lambda = lambda_vec.display(3);

  DoubleVector mH1Sq_vec = essmSusy.displayMh1Squared();
  double m1Sq = mH1Sq_vec.display(3);
  DoubleVector mH2Sq_vec = essmSusy.displayMh2Squared();
  double m2Sq = mH2Sq_vec.display(3);
  DoubleVector mSSq_vec = essmSusy.displaymSsq();
  double msSq = mSSq_vec.display(3);

  double Qtilde_1 = essmSusy.displayQH1tilde();
  double Qtilde_2 = essmSusy.displayQH2tilde();
  double Qtilde_s = essmSusy.displayQStilde();

  double df1dv1,df1dv2,df1dv3,df2dv1,df2dv2,df2dv3,df3dv1,df3dv2,df3dv3;

  df1dv1 = v1*(0.25*gbar*gbar+gd1*gd1*Qtilde_1*Qtilde_1)+lambda*Alambda*s*v2/(sqrt(2.0)*v1*v1);//m1Sq+0.5*lambda*lambda*(v2*v2+s*s)+0.125*gbar*gbar*(3.0*v1*v1-v2*v2)+0.5*gd1*gd1*Qtilde_1*
    //(3.0*Qtilde_1*v1*v1+Qtilde_2*v2*v2+Qtilde_s*s*s);
  df1dv2 = v2*(lambda*lambda-0.25*gbar*gbar+gd1*gd1*Qtilde_1*Qtilde_2)-lambda*Alambda*s/(sqrt(2.0)*v1);//v1*v2*(lambda*lambda-0.25*gbar*gbar+gd1*gd1*Qtilde_1*Qtilde_2)-lambda*Alambda*s/sqrt(2.0);
  df1dv3 = s*(lambda*lambda+gd1*gd1*Qtilde_1*Qtilde_s)-lambda*Alambda*v2/(sqrt(2.0)*v1);//s*v1*(lambda*lambda+gd1*gd1*Qtilde_1*Qtilde_s)-lambda*Alambda*v2/sqrt(2.0);

  df2dv1 = v1*(lambda*lambda-0.25*gbar*gbar+gd1*gd1*Qtilde_2*Qtilde_1)-lambda*Alambda*s/(sqrt(2.0)*v2);//df1dv2;
  df2dv2 = v2*(0.25*gbar*gbar+gd1*gd1*Qtilde_2*Qtilde_2)+lambda*Alambda*s*v1/(sqrt(2.0)*v2*v2);//m2Sq+0.5*lambda*lambda*(v1*v1+s*s)+0.125*gbar*gbar*(3.0*v2*v2-v1*v1)+0.5*gd1*gd1*Qtilde_2*
    //(Qtilde_1*v1*v1+3.0*Qtilde_2*v2*v2+Qtilde_s*s*s);
  df2dv3 = s*(lambda*lambda+gd1*gd1*Qtilde_2*Qtilde_s)-lambda*Alambda*v1/(sqrt(2.0)*v2);//s*v2*(lambda*lambda+gd1*gd1*Qtilde_2*Qtilde_s)-lambda*Alambda*v1/sqrt(2.0);

  df3dv1 = v1*(lambda*lambda+gd1*gd1*Qtilde_s*Qtilde_1)-lambda*Alambda*v2/(sqrt(2.0)*s);//df1dv3;
  df3dv2 = v2*(lambda*lambda+gd1*gd1*Qtilde_s*Qtilde_2)-lambda*Alambda*v1/(sqrt(2.0)*s);//df2dv3;
  df3dv3 = gd1*gd1*Qtilde_s*Qtilde_s*s+lambda*Alambda*v1*v2/(sqrt(2.0)*s*s);//msSq+0.5*lambda*lambda*v*v+0.5*gd1*gd1*Qtilde_s*(Qtilde_1*v1*v1+Qtilde_2*v2*v2+3.0*Qtilde_s*s*s);

  J(1,1) = df1dv1;
  J(1,2) = df1dv2;
  J(1,3) = df1dv3;
  J(2,1) = df2dv1;
  J(2,2) = df2dv2;
  J(2,3) = df2dv3;
  J(3,1) = df3dv1;
  J(3,2) = df3dv2;
  J(3,3) = df3dv3;

  //cout << "Tree level Jacobian: " << endl;
  //cout << J;

  return J;
}

DoubleMatrix EWSBJacobian_OneLoop(SoftParsEssm essmSusy, double s, double tb)
{

  double v = essmSusy.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;

  //cout << "v1 = " << v1 << endl;
  //cout << "v2 = " << v2 << endl;
  //cout << "s = " << s << endl;

  DoubleMatrix J(3,3);

  // Calculate tadpole corrections (note opposite sign convention)
  double delta1tp = -doCalcTadpoleESSMH1(essmSusy, s, tb);
  double delta2tp = -doCalcTadpoleESSMH2(essmSusy, s, tb);
  double delta3tp = -doCalcTadpolesESSMS(essmSusy, s, tb);

  // Calculate second derivatives of loop corrections to effective potential
  double delta11prime = doCalcDeltaPrime11(essmSusy, s, tb);
  double delta12prime = doCalcDeltaPrime12(essmSusy, s, tb);
  double delta13prime = doCalcDeltaPrime13(essmSusy, s, tb);
  double delta22prime = doCalcDeltaPrime22(essmSusy, s, tb);
  double delta23prime = doCalcDeltaPrime23(essmSusy, s, tb);
  double delta33prime = doCalcDeltaPrime33(essmSusy, s, tb);

  // Get the tree level Jacobian
  J = EWSBJacobian_TreeLevel(essmSusy, s, tb);

  // Add in the leading one loop corrections.
  J(1,1) = J(1,1)+(1.0/v1)*(delta11prime-delta1tp);// + delta11prime;
  J(1,2) = J(1,2)+delta12prime/v1;// + delta12prime;
  J(1,3) = J(1,3)+delta13prime/v1;// + delta13prime;
  J(2,1) = J(2,1)+delta12prime/v2;// + delta12prime;
  J(2,2) = J(2,2)+(1.0/v2)*(delta22prime-delta2tp);// + delta22prime;
  J(2,3) = J(2,3)+delta23prime/v2;// + delta23prime;
  J(3,1) = J(3,1)+delta13prime/s;// + delta13prime;
  J(3,2) = J(3,2)+delta23prime/s;// + delta23prime;
  J(3,3) = J(3,3)+(1.0/s)*(delta33prime-delta3tp);// + delta33prime;

  //cout << "One loop Jacobian: " << endl;
  //  cout << J;

  return J;
}

// qrUpdateEWSBSoln_TreeLevel and _OneLoop use a QR decomposition
// to compute the update step in Newton's method, to solve for
// the VEVs at the parameter point contained in the ESSM object.
// Returns the updated estimate for the solution to the EWSB 
// conditions. Note that the provided object is assumed to 
// contain the initial (current) guess for the solutions to
// the EWSB conditions.
// Inputs:
//     SoftParsEssm essmSusy = the object to calculate the solution for
//     double s = the value of the VEV s to use
//     double tb = the value of tan(beta) to use
//     int & sing = integer flag to indicate if singular (in which case sing is non-zero)
DoubleVector qrUpdateEWSBSoln_TreeLevel(SoftParsEssm essmSusy, double s, double tb, int & sing)
{
  // Construct tree level Jacobian
  DoubleMatrix J_mat = EWSBJacobian_TreeLevel(essmSusy, s, tb);
  // Construct right hand side
  DoubleVector rk(3);
  rk.set(1, -EWSBCondition1_TreeLevel(essmSusy, s, tb));
  rk.set(2, -EWSBCondition2_TreeLevel(essmSusy, s, tb));
  rk.set(3, -EWSBCondition3_TreeLevel(essmSusy, s, tb));

  // Solve for update step
  DoubleVector deltaStep(3);
  DoubleMatrix Q_mat(3,3);

  QR_decomp(J_mat, Q_mat, 3, sing);
  QR_solve(Q_mat, J_mat, 3, rk, deltaStep, sing);

  // Update the solution and return the new guess.
  DoubleVector initGuess(3), updateGuess(3);
  double v = essmSusy.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;

  initGuess.set(1, v1);
  initGuess.set(2, v2);
  initGuess.set(3, s);

  updateGuess = initGuess + deltaStep;

  return updateGuess;
}

DoubleVector qrUpdateEWSBSoln_OneLoop(SoftParsEssm essmSusy, double s, double tb, int & sing)
{
  // Construct one loop Jacobian
  DoubleMatrix J_mat = EWSBJacobian_OneLoop(essmSusy, s, tb);

  //cout << "Jacobian: " << endl;
  //cout << J_mat;

  // Construct right hand side
  DoubleVector rk(3);
  rk.set(1, -EWSBCondition1_OneLoop(essmSusy, s, tb));
  rk.set(2, -EWSBCondition2_OneLoop(essmSusy, s, tb));
  rk.set(3, -EWSBCondition3_OneLoop(essmSusy, s, tb));

  //cout << "Residual: " << endl;
  //cout << rk;

  // Solve for update step
  DoubleVector deltaStep(3);
  DoubleMatrix Q_mat(3,3);

  QR_decomp(J_mat, Q_mat, 3, sing);
  //cout << "QR decomposition Q: " << endl;
  //cout << Q_mat;
  //cout << "QR decomposition R: " << endl;
  //cout << J_mat;

  //cout << "Check: QR = " << endl;
  //cout << Q_mat*J_mat;

  //cout << "RHS: " << endl;
  //cout << Q_mat.transpose()*rk;

  QR_solve(Q_mat, J_mat, 3, rk, deltaStep, sing);

  //cout << "Solution: " <<endl;
  //cout << deltaStep;

  // Update the solution and return the new guess.
  DoubleVector initGuess(3), updateGuess(3);
  double v = essmSusy.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;

  initGuess.set(1, v1);
  initGuess.set(2, v2);
  initGuess.set(3, s);

  updateGuess = initGuess + deltaStep;

  return updateGuess;
}

// EWSBNewtonSolver_TreeLevel and _OneLoop computes the solution to the 
// EWSB conditions at tree level or one loop order using Newton's method.
// Returns the solution for the VEVs v_1, v_2 and s as a DoubleVector.
// Note that the provided object is taken to contain both the values of
// the input parameters and the initial guess for the VEVs (with the exception
// of s, of course). 
// Inputs:
//     SoftParsEssm essmSusy = the object to calculate the solution for
//     double s = the initial guess for the value of the VEV s
//     double tb = the initial guess for the value of tan(beta)
//     double tol = the required tolerance to test for convergence
//     int maxIters = the maximum number of allowed iterations
//     int & sing = integer flag to indicate if singular (in which case sing is non-zero)
DoubleVector EWSBNewtonSolver_TreeLevel(SoftParsEssm essmSusy, double s, double tb, double tol, int maxIters, int & sing)
{
  double v = essmSusy.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;

  DoubleVector initGuess(3);
  initGuess.set(1, v1);
  initGuess.set(2, v2);
  initGuess.set(3, s);

  DoubleVector oldVevs(3), newVevs(3), stepDiff(3);
  oldVevs = initGuess;
  double old_tb = tb;
  double old_v = v;

  double f1 = EWSBCondition1_TreeLevel(essmSusy, s, tb);
  double f2 = EWSBCondition2_TreeLevel(essmSusy, s, tb);
  double f3 = EWSBCondition3_TreeLevel(essmSusy, s, tb);

  for (int i = 0; i < maxIters; i++)
    {
      newVevs = qrUpdateEWSBSoln_TreeLevel(essmSusy, oldVevs.display(3), old_tb, sing);
      stepDiff = newVevs-oldVevs;
      oldVevs = newVevs;
      old_v = sqrt(oldVevs.display(1)*oldVevs.display(1)+oldVevs.display(2)*oldVevs.display(2));
      essmSusy.setHvev(old_v);
      old_tb = newVevs.display(2)/newVevs.display(1);
      essmSusy.setTanb(old_tb);

      f1 = EWSBCondition1_TreeLevel(essmSusy, oldVevs.display(3), old_tb);
      f2 = EWSBCondition2_TreeLevel(essmSusy, oldVevs.display(3), old_tb);
      f3 = EWSBCondition3_TreeLevel(essmSusy, oldVevs.display(3), old_tb);

      //if (sqrt(stepDiff.dot(stepDiff)) < tol) break;
      // Alternative convergence condition
      if (sqrt(f1*f1+f2*f2+f3*f3)<tol) break;
    }

  // Error flag if maximum iterations reached before solution converged
  bool speak = true;
  if (sqrt(f1*f1+f2*f2+f3*f3)>tol)
    {
      sing = 3; // arbitrary value of 3 indicates warning
      if (speak)
	{
	  cout << "Warning: solution has not converged to requested tolerance." << endl;
	  cout << "Obtained tolerance = " << sqrt(f1*f1+f2*f2+f3*f3) << endl;
	}
    }

  return newVevs;
}

DoubleVector EWSBNewtonSolver_OneLoop(SoftParsEssm essmSusy, double s, double tb, double tol, int maxIters, int & sing)
{
 double v = essmSusy.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;

  DoubleVector initGuess(3);
  initGuess.set(1, v1);
  initGuess.set(2, v2);
  initGuess.set(3, s);

  DoubleVector oldVevs(3), newVevs(3), stepDiff(3);
  oldVevs = initGuess;
  double old_tb = tb;
  double old_v = v;

  double f1 = EWSBCondition1_OneLoop(essmSusy, s, tb);
  double f2 = EWSBCondition2_OneLoop(essmSusy, s, tb);
  double f3 = EWSBCondition3_OneLoop(essmSusy, s, tb);

  //cout << "f1 = " << f1 << endl;
  //cout << "f2 = " << f2 << endl;
  //cout << "f3 = " << f3 << endl;

  for (int i = 0; i < maxIters; i++)
    {
      newVevs = qrUpdateEWSBSoln_OneLoop(essmSusy, oldVevs.display(3), old_tb, sing);
      stepDiff = newVevs-oldVevs;
      //      cout << "New vevs: " << endl;
      //     cout << newVevs;
      oldVevs = newVevs;
      old_v = sqrt(oldVevs.display(1)*oldVevs.display(1)+oldVevs.display(2)*oldVevs.display(2));
      essmSusy.setHvev(old_v);
      old_tb = newVevs.display(2)/newVevs.display(1);
      essmSusy.setTanb(old_tb);

      f1 = EWSBCondition1_OneLoop(essmSusy, oldVevs.display(3), old_tb);
      f2 = EWSBCondition2_OneLoop(essmSusy, oldVevs.display(3), old_tb);
      f3 = EWSBCondition3_OneLoop(essmSusy, oldVevs.display(3), old_tb);

      //cout << "f1 = " << f1 << endl;
      //cout << "f2 = " << f2 << endl;
      //cout << "f3 = " << f3 << endl;


      //if (sqrt(stepDiff.dot(stepDiff)) < tol) break;
      // Alternative convergence condition
      if (sqrt(f1*f1+f2*f2+f3*f3)<tol) break;
    }

  // Error flag if maximum iterations reached before solution converged
  bool speak = true;
  if (sqrt(f1*f1+f2*f2+f3*f3)>tol)
    {
      sing = 3; // arbitrary value of 3 indicates warning
      if (speak)
	{
	  cout << "Warning: solution has not converged to requested tolerance." << endl;
	  cout << "Obtained tolerance = " << sqrt(f1*f1+f2*f2+f3*f3) << endl;
	}
    }

  return newVevs;
}

// vevNumericalDerivatives_TreeLevel and _OneLoop compute estimates for the partial
// derivatives of the Higgs and singlet VEVs wrt the low energy input parameters
// \lambda, A_\lambda,... at either tree level or one loop order. This is done
// using Newton's method and a 5-pt stencil approximation to the derivative. Note that
// the result calculated is the directional derivative in the direction specified by
// epsVec. The parameter point at which the derivative is to be evaluated, together
// with the initial estimate for the VEVs, is assumed to be contained in the provided
// ESSM model object.
// Inputs:
//     SoftParsEssm essmSusy = the object to calculate the derivatives for
//     double s = the value of the VEV s to use
//     double tb = the value of tan(beta) to use
//     DoubleVector epsVec = a vector indicating the direction in the low-energy
//                           parameter space in which to calculate the directional
//                           derivative. For the tree level method, epsVec must be
//                           of length 5, corresponding to the vector in parameter
//                           space (\lambda, A_\lambda, m_1^2, m_2^2, m_S^2). For
//                           the one loop method epsVec must be of length 8, containing
//                           the above 5 values followed by (m_Q^2, m_U^2, A_t). For
//                           ordinary partial derivatives, epsVec should have one element
//                           equal to unity, and all others vanishing.
//     double epsilon = the step length to use in the finite difference approximation
//     double tol = the tolerance required for Newton's method to converge
//     int maxIters = the maximum number of allowed iterations in Newton's method
//     int & sing = an integer flag indicating if singular (in which case sing is non-zero)
DoubleVector vevNumericalDerivatives_TreeLevel(SoftParsEssm essmSusy,double s,double tb,DoubleVector epsVec,double epsilon,
					       double tol,int maxIters,int & sing)
{
  // Extract parameter point to evaluate at.
  DoubleVector lambda_vec = essmSusy.displaylambda();
  double lambda = lambda_vec.display(3);
  DoubleVector Alambda_vec = essmSusy.displayA_lambda();
  double Alambda = Alambda_vec.display(3);
  DoubleVector m1Sq_vec = essmSusy.displayMh1Squared();
  double m1Sq = m1Sq_vec.display(3);
  DoubleVector m2Sq_vec = essmSusy.displayMh2Squared();
  double m2Sq = m2Sq_vec.display(3);
  DoubleVector msSq_vec = essmSusy.displaymSsq();
  double msSq = msSq_vec.display(3);

  double v = essmSusy.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;

  DoubleVector initVevs(3), twoStepForwardVevs(3), oneStepForwardVevs(3), oneStepBackwardVevs(3), twoStepBackwardVevs(3);
  initVevs.set(1, v1);
  initVevs.set(2, v2);
  initVevs.set(3, s);

  double newLambda, newAlambda, newm1Sq, newm2Sq, newmsSq;
  int sing1 = 0, sing2 = 0, sing3 = 0, sing4 = 0;

  // Obtain VEVs two steps ahead in given direction
  newLambda = lambda+2.0*epsilon*epsVec.display(1);
  newAlambda = Alambda+2.0*epsilon*epsVec.display(2);
  newm1Sq = m1Sq+2.0*epsilon*epsVec.display(3);
  newm2Sq = m2Sq+2.0*epsilon*epsVec.display(4);
  newmsSq = msSq+2.0*epsilon*epsVec.display(5);

  lambda_vec.set(3, newLambda);
  Alambda_vec.set(3, newAlambda);
  m1Sq_vec.set(3, newm1Sq);
  m2Sq_vec.set(3, newm2Sq);
  msSq_vec.set(3, newmsSq);

  essmSusy.setlambda(lambda_vec);
  essmSusy.setA_lambda(Alambda_vec);
  essmSusy.setMh1Squared(m1Sq_vec);
  essmSusy.setMh2Squared(m2Sq_vec);
  essmSusy.setmSsq(msSq_vec);

  twoStepForwardVevs = EWSBNewtonSolver_TreeLevel(essmSusy, s, tb, tol, maxIters, sing1);

  // Obtain VEVs one step ahead in given direction
  newLambda = lambda+epsilon*epsVec.display(1);
  newAlambda = Alambda+epsilon*epsVec.display(2);
  newm1Sq = m1Sq+epsilon*epsVec.display(3);
  newm2Sq = m2Sq+epsilon*epsVec.display(4);
  newmsSq = msSq+epsilon*epsVec.display(5);

  lambda_vec.set(3, newLambda);
  Alambda_vec.set(3, newAlambda);
  m1Sq_vec.set(3, newm1Sq);
  m2Sq_vec.set(3, newm2Sq);
  msSq_vec.set(3, newmsSq);

  essmSusy.setlambda(lambda_vec);
  essmSusy.setA_lambda(Alambda_vec);
  essmSusy.setMh1Squared(m1Sq_vec);
  essmSusy.setMh2Squared(m2Sq_vec);
  essmSusy.setmSsq(msSq_vec);

  oneStepForwardVevs = EWSBNewtonSolver_TreeLevel(essmSusy, s, tb, tol, maxIters, sing2);

  // Obtain VEVs one step behind in given direction
  newLambda = lambda-epsilon*epsVec.display(1);
  newAlambda = Alambda-epsilon*epsVec.display(2);
  newm1Sq = m1Sq-epsilon*epsVec.display(3);
  newm2Sq = m2Sq-epsilon*epsVec.display(4);
  newmsSq = msSq-epsilon*epsVec.display(5);

  lambda_vec.set(3, newLambda);
  Alambda_vec.set(3, newAlambda);
  m1Sq_vec.set(3, newm1Sq);
  m2Sq_vec.set(3, newm2Sq);
  msSq_vec.set(3, newmsSq);

  essmSusy.setlambda(lambda_vec);
  essmSusy.setA_lambda(Alambda_vec);
  essmSusy.setMh1Squared(m1Sq_vec);
  essmSusy.setMh2Squared(m2Sq_vec);
  essmSusy.setmSsq(msSq_vec);

  oneStepBackwardVevs = EWSBNewtonSolver_TreeLevel(essmSusy, s, tb, tol, maxIters, sing3);

  // Obtain VEVs two steps behind in given direction
  newLambda = lambda-2.0*epsilon*epsVec.display(1);
  newAlambda = Alambda-2.0*epsilon*epsVec.display(2);
  newm1Sq = m1Sq-2.0*epsilon*epsVec.display(3);
  newm2Sq = m2Sq-2.0*epsilon*epsVec.display(4);
  newmsSq = msSq-2.0*epsilon*epsVec.display(5);

  lambda_vec.set(3, newLambda);
  Alambda_vec.set(3, newAlambda);
  m1Sq_vec.set(3, newm1Sq);
  m2Sq_vec.set(3, newm2Sq);
  msSq_vec.set(3, newmsSq);

  essmSusy.setlambda(lambda_vec);
  essmSusy.setA_lambda(Alambda_vec);
  essmSusy.setMh1Squared(m1Sq_vec);
  essmSusy.setMh2Squared(m2Sq_vec);
  essmSusy.setmSsq(msSq_vec);

  twoStepBackwardVevs = EWSBNewtonSolver_TreeLevel(essmSusy, s, tb, tol, maxIters, sing4);


  // Estimate the directional derivative using a 5-pt stencil
  DoubleVector dirDeriv(3);
  dirDeriv = (1.0/(12.0*epsilon))*(-1.0*twoStepForwardVevs+8.0*oneStepForwardVevs-8.0*oneStepBackwardVevs+twoStepBackwardVevs);

  if (sing1 == 0 && sing2 == 0 && sing3 == 0 && sing4 == 0)
    {
      sing = 0;
    }
  else
    {
      sing = 1;
    }

  return dirDeriv;

}

DoubleVector vevNumericalDerivatives_OneLoop(SoftParsEssm essmSusy,double s,double tb,DoubleVector epsVec,double epsilon,
					     double tol,int maxIters,int & sing)
{
  // Extract parameter point to evaluate at.
  DoubleVector lambda_vec = essmSusy.displaylambda();
  double lambda = lambda_vec.display(3);
  DoubleVector Alambda_vec = essmSusy.displayA_lambda();
  double Alambda = Alambda_vec.display(3);
  DoubleVector m1Sq_vec = essmSusy.displayMh1Squared();
  double m1Sq = m1Sq_vec.display(3);
  DoubleVector m2Sq_vec = essmSusy.displayMh2Squared();
  double m2Sq = m2Sq_vec.display(3);
  DoubleVector msSq_vec = essmSusy.displaymSsq();
  double msSq = msSq_vec.display(3);
  double mQlSq = essmSusy.displaySoftMassSquared(mQl, 3, 3);
  double mUrSq = essmSusy.displaySoftMassSquared(mUr, 3, 3);
  double At = essmSusy.displayA_top();

  double v = essmSusy.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;

  DoubleVector initVevs(3), twoStepForwardVevs(3), oneStepForwardVevs(3), oneStepBackwardVevs(3), twoStepBackwardVevs(3);
  initVevs.set(1, v1);
  initVevs.set(2, v2);
  initVevs.set(3, s);

  double newLambda, newAlambda, newm1Sq, newm2Sq, newmsSq, newmQlSq, newmUrSq, newAt;
  int sing1 = 0, sing2 = 0, sing3 = 0, sing4 = 0;

  // Obtain VEVs two steps ahead in given direction
  newLambda = lambda+2.0*epsilon*epsVec.display(1);
  newAlambda = Alambda+2.0*epsilon*epsVec.display(2);
  newm1Sq = m1Sq+2.0*epsilon*epsVec.display(3);
  newm2Sq = m2Sq+2.0*epsilon*epsVec.display(4);
  newmsSq = msSq+2.0*epsilon*epsVec.display(5);
  newmQlSq = mQlSq+2.0*epsilon*epsVec.display(6);
  newmUrSq = mUrSq+2.0*epsilon*epsVec.display(7);
  newAt = At+2.0*epsilon*epsVec.display(8);

  lambda_vec.set(3, newLambda);
  Alambda_vec.set(3, newAlambda);
  m1Sq_vec.set(3, newm1Sq);
  m2Sq_vec.set(3, newm2Sq);
  msSq_vec.set(3, newmsSq);

  essmSusy.setlambda(lambda_vec);
  essmSusy.setA_lambda(Alambda_vec);
  essmSusy.setMh1Squared(m1Sq_vec);
  essmSusy.setMh2Squared(m2Sq_vec);
  essmSusy.setmSsq(msSq_vec);
  essmSusy.setSoftMassElement(mQl, 3, 3, newmQlSq);
  essmSusy.setSoftMassElement(mUr, 3, 3, newmUrSq);
  essmSusy.setA_top(newAt);

  twoStepForwardVevs = EWSBNewtonSolver_OneLoop(essmSusy, s, tb, tol, maxIters, sing1);

  // Obtain VEVs one step ahead in given direction
  newLambda = lambda+epsilon*epsVec.display(1);
  newAlambda = Alambda+epsilon*epsVec.display(2);
  newm1Sq = m1Sq+epsilon*epsVec.display(3);
  newm2Sq = m2Sq+epsilon*epsVec.display(4);
  newmsSq = msSq+epsilon*epsVec.display(5);
  newmQlSq = mQlSq+epsilon*epsVec.display(6);
  newmUrSq = mUrSq+epsilon*epsVec.display(7);
  newAt = At+epsilon*epsVec.display(8);

  lambda_vec.set(3, newLambda);
  Alambda_vec.set(3, newAlambda);
  m1Sq_vec.set(3, newm1Sq);
  m2Sq_vec.set(3, newm2Sq);
  msSq_vec.set(3, newmsSq);

  essmSusy.setlambda(lambda_vec);
  essmSusy.setA_lambda(Alambda_vec);
  essmSusy.setMh1Squared(m1Sq_vec);
  essmSusy.setMh2Squared(m2Sq_vec);
  essmSusy.setmSsq(msSq_vec);
  essmSusy.setSoftMassElement(mQl, 3, 3, newmQlSq);
  essmSusy.setSoftMassElement(mUr, 3, 3, newmUrSq);
  essmSusy.setA_top(newAt);

  oneStepForwardVevs = EWSBNewtonSolver_OneLoop(essmSusy, s, tb, tol, maxIters, sing2);

  // Obtain VEVs one step behind in given direction
  newLambda = lambda-epsilon*epsVec.display(1);
  newAlambda = Alambda-epsilon*epsVec.display(2);
  newm1Sq = m1Sq-epsilon*epsVec.display(3);
  newm2Sq = m2Sq-epsilon*epsVec.display(4);
  newmsSq = msSq-epsilon*epsVec.display(5);
  newmQlSq = mQlSq-epsilon*epsVec.display(6);
  newmUrSq = mUrSq-epsilon*epsVec.display(7);
  newAt = At-epsilon*epsVec.display(8);

  lambda_vec.set(3, newLambda);
  Alambda_vec.set(3, newAlambda);
  m1Sq_vec.set(3, newm1Sq);
  m2Sq_vec.set(3, newm2Sq);
  msSq_vec.set(3, newmsSq);

  essmSusy.setlambda(lambda_vec);
  essmSusy.setA_lambda(Alambda_vec);
  essmSusy.setMh1Squared(m1Sq_vec);
  essmSusy.setMh2Squared(m2Sq_vec);
  essmSusy.setmSsq(msSq_vec);
  essmSusy.setSoftMassElement(mQl, 3, 3, newmQlSq);
  essmSusy.setSoftMassElement(mUr, 3, 3, newmUrSq);
  essmSusy.setA_top(newAt);

  oneStepBackwardVevs = EWSBNewtonSolver_OneLoop(essmSusy, s, tb, tol, maxIters, sing3);


  // Obtain VEVs two steps behind in given direction
  newLambda = lambda-2.0*epsilon*epsVec.display(1);
  newAlambda = Alambda-2.0*epsilon*epsVec.display(2);
  newm1Sq = m1Sq-2.0*epsilon*epsVec.display(3);
  newm2Sq = m2Sq-2.0*epsilon*epsVec.display(4);
  newmsSq = msSq-2.0*epsilon*epsVec.display(5);
  newmQlSq = mQlSq-2.0*epsilon*epsVec.display(6);
  newmUrSq = mUrSq-2.0*epsilon*epsVec.display(7);
  newAt = At-2.0*epsilon*epsVec.display(8);

  lambda_vec.set(3, newLambda);
  Alambda_vec.set(3, newAlambda);
  m1Sq_vec.set(3, newm1Sq);
  m2Sq_vec.set(3, newm2Sq);
  msSq_vec.set(3, newmsSq);

  essmSusy.setlambda(lambda_vec);
  essmSusy.setA_lambda(Alambda_vec);
  essmSusy.setMh1Squared(m1Sq_vec);
  essmSusy.setMh2Squared(m2Sq_vec);
  essmSusy.setmSsq(msSq_vec);
  essmSusy.setSoftMassElement(mQl, 3, 3, newmQlSq);
  essmSusy.setSoftMassElement(mUr, 3, 3, newmUrSq);
  essmSusy.setA_top(newAt);

  twoStepBackwardVevs = EWSBNewtonSolver_OneLoop(essmSusy, s, tb, tol, maxIters, sing4);

  // Estimate the directional derivative using a 5-pt stencil
  DoubleVector dirDeriv(3);
  dirDeriv = (1.0/(12.0*epsilon))*(-1.0*twoStepForwardVevs+8.0*oneStepForwardVevs-8.0*oneStepBackwardVevs+twoStepBackwardVevs);

  if (sing1 == 0 && sing2 == 0 && sing3 == 0 && sing4 == 0)
    {
      sing = 0;
    }
  else
    {
      sing = 1;
    }

  return dirDeriv;

}

/*
--------------------------------------------------------------------------------------
 */

// A method to generate the matrices Q and R in the QR decomposition
// A = QR of the n x n square matrix A. On output the original input matrix A is
// replaced by the upper triangular matrix R (the original matrix
// can always be recovered as A = QR if necessary). Not exactly optimised,
// but should be sufficient for now when only working with 3 x 3 matrices.
void QR_decomp(DoubleMatrix & A, DoubleMatrix & Q, int n, int & sing)
{
  // Construct the n-1 Householder vectors u. Because of the 
  // scaling process used in the SOFTSUSY numerics routines
  // the definition differs slightly from the usual one.
  DoubleMatrix u(n,n-1);
  // Get the QR decomposition
  DoubleVector c(n-1), d(n), u_temp(n);
  qrdcmp(A, n, c, d, sing);

  for (int k = 1; k < n; k++)
    {
      for (int i = 1; i <= n; i++)
	{
	  if (i < k)
	    {
	      u(i,k) = 0.0;
	    }
	  else
	    {
	      u(i,k) = A(i,k);
	    }
	}
    }

  // Reconstruct the matrix Q  
  DoubleMatrix eye(n,n);
  for (int j = 1; j <= n; j++) eye(j,j) = 1.0;

  Q = eye;
  DoubleMatrix Qi(n,n);

  for (int j = 1; j < n; j++)
    {
      for (int i = 1; i <= n; i++) u_temp.set(i, u(i,j));
      if (c(j) == 0.0)
	{
	  Q = Q*eye;
	}
      else
	{
	  Qi = eye - (1.0/c(j))*outerProduct(u_temp,u_temp);
	  Q = Q*Qi;
	}
    }

  // Put the output matrix A in upper triangular form.
  for (int j = 1; j <= n; j++) A(j,j) = d(j);

  for (int k = 1; k <= n; k++)
    {
      for (int j = k+1; j <= n; j++)
	{
	  A(j,k) = 0.0;
	}
    }

  return;
}

// A method to solve a linear system by using the QR decomposition to
// do back-substitution. Requires the matrices Q and R to be supplied.
void QR_solve(DoubleMatrix const & Q, DoubleMatrix const & R, int n, DoubleVector const & b, DoubleVector & x, int & sing)
{
  DoubleVector rhs = Q.transpose() * b;
  double sum;
  // Solve by back-substitution, storing result in x.
  for (int k = n; k > 0; k--)
    {
      sum = 0.0;
      for (int i = k+1; i <= n; i++)
	{
	  sum += R(k,i)*x.display(i);
	}
      if (R(k,k) != 0.0)
	{
	  x.set(k, (rhs.display(k)-sum)/R(k,k));
	}
      else if (R(k,k) == 0.0 & sum == 0.0 && rhs.display(k) == 0.0)
	{
	  sing = 1;
	  x.set(k,1.0);
	}
      else
	{
	  sing = 2;
	  x.set(k,0.0);
	}
    }
  return;
}

// A method to solve a 3x3 linear system using Cramer's rule.
// Inputs:
//     DoubleMatrix const & A = the matrix of coefficients
//     DoubleVector const & = the vector on the RHS
//     double & detA = the calculated value of the determinant of A
//     int & sing = integer flag to indicate if no unique soln (non-zero if singular)
//     DoubleVector & soln = a vector of length 3 containing the solution
void CramersRule(DoubleMatrix const & A, DoubleVector const & b, double & detA, int & sing, DoubleVector & soln)
{
  detA = A(1,1)*(A(2,2)*A(3,3)-A(3,2)*A(2,3))-A(1,2)*(A(2,1)*A(3,3)-A(3,1)*A(2,3))+A(1,3)*(A(2,1)*A(3,2)-A(3,1)*A(2,2));

  if (detA == 0.0)
  {
    sing = 1;
    soln(1) = 0.0;
    soln(2) = 0.0;
    soln(3) = 0.0;
  }
  else
    {
      // If det(A) is non-zero, get each component of the solution using Cramer's rule.
      DoubleMatrix Asub(3,3);
      double detAsub;
      for (int i = 1; i <= 3; i++)
	{

	  Asub = A;
	  Asub(1,i) = b(1);
	  Asub(2,i) = b(2);
	  Asub(3,i) = b(3);

	  detAsub = Asub(1,1)*(Asub(2,2)*Asub(3,3)-Asub(3,2)*Asub(2,3))
	    -Asub(1,2)*(Asub(2,1)*Asub(3,3)-Asub(3,1)*Asub(2,3))+Asub(1,3)*(Asub(2,1)*Asub(3,2)-Asub(3,1)*Asub(2,2));

	  soln.set(i, detAsub/detA);

	}
    }
 
  return;
}

