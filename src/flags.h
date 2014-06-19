/*
  flags.h contains the global flags used by our fine tuning code.
 */

#ifndef FLAGS_H
#define FLAGS_H

// Flag to include 1-loop top and stop tadpole contributions
// to the effective potential and fine tuning measures.
// Default is true.
extern bool INCLUDE1LPTADPOLES;

// Use m_t(m_t) as scale in calculations
extern bool USEMTOFMT;

// Number of parameters in the pMSSM as defined by Rizzo et al,
// see hep-ph/1206.5800v2
const int NUMPMSSMPARS = 20;


int const LOOPS = 2; //< number of loops to use in RG running
int const THRESH = 0; //< threshold accuracy
int const NUMLOOPSHIGGS = 2; //< number of loops in Higgs pole mass calculation
int const NUMLOOPSREWSB = 2; //< number of loops in FlexibleSUSY EWSB conditions

double const TOLEWSB = 50.0; //< tolerance to impose on minimisation conditions, GeV^2

int const UFBPROBLEM = 1; //< unbounded from below flag
int const CCBPROBLEM = 2; //< charge/colour breaking minimum flag

int const EWSBPROBLEM = 3; //< problem with iteration

int const WRONGVACUUM = 4; //< problem with unphysical vacuum

int const SQUARKTACHYON = 5; //< stop/sbottom tachyon
int const VECTORBOSONTACHYON = 6; //< vector boson tachyon
int const TADPOLESPROBLEM = 15; //< problem calculating 1-loop tadpoles

int const HIGGSPROBLEM = 10; //< problems with calculating the Higgs mass
int const NOTEXPVALID = 30; //< point not experimentally valid
int const HIGGSTACHYON = 31; //< ruled out because tachyonic
int const POLEHIGGSTACHYON = 32; //< ruled out because physical Higgs is tachyon

int const NUMERICALPROBLEM = 666; //< for serious numerical problems

int const TUNINGERROR = 35;

double const HIGGSCENT = 125.0; //< rough central value for Higgs mass, GeV
double const HIGGSERROR = 7.0; //< theory error in Higgs calculation, GeV

#endif
