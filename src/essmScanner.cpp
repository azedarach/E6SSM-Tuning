/*
  Scan over the E6SSM parameter space and calculate the fine tuning 
  at each point.
 */

#include "essmsusy.h"
#include "flags.h"
#include "tuningnumerics.h"
#include "tuningutils.h"

using namespace softsusy;

int main(int argc, char *argv[])
{
  EssmSusy test = EssmSusy();

  DoubleVector vals(numEssmSusyPars);
  vals(vals.displayEnd()) = atan(sqrt(15.0)); //< E6SSM with U(1)_N

  test.set(vals);

  cout << "U(1)_N charges: " << endl;
  cout << test.displayU1Prime();

  return 0;
}
