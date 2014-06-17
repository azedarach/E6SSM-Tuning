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


#endif
