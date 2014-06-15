/*
  mssmTuningAtPoint reads in an SLHA input file and then uses
  SOFTSUSY routines to calculate the spectrum for that file.
  Also calculates the fine tuning, either using our analytic
  expressions or numerically using the RGE running in SOFTSUSY.
  Then prints output if the point is allowed by any
  constraints that are imposed on it. For saving to file
  should pipe output, e.g. syntax is
  "./mssmTuningAtPoint < input_file > output_file".
  Partially written because I can't get SOFTSUSY 3.3.9 to work
  with SLHA input in which tan(beta) is specified at MX and not
  MZ, which is a feature I need for comparing with the E6SSM.
  Likewise for m_A^2(MX), which I need instead of m_A(pole).
  As a result this is essentially the code I need pulled out from
  softpoint.cpp, with the fine tuning calculation attached to
  the end and a modified output.
  Author: Dylan Harries
  Date created: 2/4/2014
 */

