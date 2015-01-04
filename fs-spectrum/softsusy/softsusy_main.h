
/** \file main.h
   - Project:     SOFTSUSY 
   - Author:      Ben Allanach 
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/

   \brief a main C++ program to calculate Higgs masses as a function of tan
   beta 
*/

/** \mainpage Detailed SOFTSUSY Documentation 

    \section install Installation or downloads
    For installation instructions or a download, please go to the 
    <a href="http://projects.hepforge.org/softsusy/">
    SOFTSUSY Homepage</a>

    \section manual Official manual
    If you use SOFTSUSY to write a paper, please cite 
    <a href="http://xxx.arxiv.org/abs/hep-ph/0104145">
    B.C. Allanach, Comput. Phys. Commun. 143 (2002) 305-331, hep-ph/0104145, 
    </a> which is the SOFTSUSY manual.
    If you use the R-parity violating aspects, please cite
    <a href="http://arxiv.org/abs/0903.1805">B.C. Allanach and M.A. Bernhardt, 
    Comput. Phys. Commun. 181 (2010) 232-245,
    arXiv:0903.1805</a>.
    For RPV neutrino oscillation output, please cite 
    <a href="http://xxx.soton.ac.uk/abs/1109.3735">B.C. Allanach,
    M. Hanussek and C.H. Kom, Comput. Phys. Commun. 183 (2012) 785,
    arXiv:1109.3735</a>. 

    \section documentation Documentation
    These web-pages contain the documentation of the latest SOFTSUSY code.
    There are class diagrams and cross-referenced links a la doxygen to help 
    you navigate.

    \section updates Official Updates
    Updates will be posted on the    
    <a href="http://projects.hepforge.org/softsusy/">
    SOFTSUSY Homepage</a>, and the manuals in the doc subdirectory
    will also be updated and released with the distribution.
 */

#include <iostream>
#include "softsusy_mycomplex.h"
#include "softsusy_def.h"
#include "softsusy_linalg.h"
#include "softsusy_lowe.h"
#include "softsusy_rge.h"
#include "softsusy_softsusy.h"
#include "softsusy_softpars.h"
#include "softsusy_susy.h"
#include "softsusy_utils.h"
#include "softsusy_numerics.h"
#include "softsusy_rpvsoft.h"
using namespace essmsoftsusy;

namespace NR {
  extern int nn;
  extern DoubleVector fvec;
}

  extern void (*nrfuncv)(int n, DoubleVector v, DoubleVector & f);

