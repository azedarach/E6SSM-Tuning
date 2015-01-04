
/** \file softpoint.h
   - Project:     SOFTSUSY 
   - Author:      Ben Allanach, Markus Bernhardt 
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/
   - Description: main calling program: command line interface
   \brief Main program: reads input in SLHA or command line format
*/

#include <iostream>
#include <sstream>
#include <cstring>
#include "softsusy_mycomplex.h"
#include "softsusy_def.h"
#include "softsusy_linalg.h"
#include "softsusy_lowe.h"
#include "softsusy_rge.h"
#include "softsusy_softsusy.h"
#include "softsusy_flavoursoft.h"
#include "softsusy_softpars.h"
#include "softsusy_physpars.h"
#include "softsusy_susy.h"
#include "softsusy_utils.h"
#include "softsusy_numerics.h"
#include "softsusy_twoloophiggs.h"
#include "softsusy_rpvneut.h"
using namespace essmsoftsusy;

/// Requested by CMS
void splitGmsb(MssmSoftsusy & m, const DoubleVector & inputParameters);

/// Does the user require gauge unification or not -- gaugeUnification changed
/// to be correct value
inline double mgutCheck(char * a, bool & gaugeUnification, 
			bool & ewsbBCscale) { 
  gaugeUnification = false; ewsbBCscale = false;
  if (!strcmp(a, "?") || !strcmp(a,"unified")) {
    gaugeUnification = true; 
    return 2.0e16;
  }
  if (!strcmp(a, "msusy")) {
    ewsbBCscale = true;
    return 1.0e3;
  }
  else return atof(a);
}

/// Incorrect input: gives advice on how to supply it
void errorCall();

