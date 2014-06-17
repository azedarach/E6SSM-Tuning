/*
  Defines default values of any global flags used in the fine tuning code.
 */

#include "flags.h"

// Do we include the 1-loop top and stop tadpole
// contributions to the effective potential and the tuning measure?
// Default is true, but this can be changed at run time.
bool INCLUDE1LPTADPOLES = true;

// Use m_t(m_t) as scale in calculations. Default is false
// (except in 2-loop Higgs corrections), but can be changed at run
// time.
bool USEMTOFMT = false;
