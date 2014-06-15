/*
  Defines default values of any global flags used in the fine tuning code.
 */

#include "flags.h"

// Do we include the 1-loop top and stop tadpole
// contributions to the effective potential and the tuning measure?
// Default is true, but this can be changed at run time.
bool INCLUDE1LPTADPOLES = true;
