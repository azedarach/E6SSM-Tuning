/*
  mssmtuningutils contains all of the functions needed for calculating the
  fine tuning in the (p)MSSM.
  NOTE THIS WILL ONLY RUN FOR PMSSM MODELS
 */

#include "mssmtuningutils.h"

// A function to calculate the 2x2 matrix appearing on the LHS of the calculation
// of the derivatives of the VEVs wrt the parameters in the MSSM. Note that this
// assumes that the SoftParsMssm object is already set to have its renormalisation
// scale set at M_SUSY.
DoubleMatrix doCalcLHSTuningMatrix(SoftParsMssm r, DoubleVector const & vevs)
{
  // bool INCLUDE1LPTADPOLES = true;

  DoubleMatrix lhsMat(2,2);

  if (vevs.displayEnd()-vevs.displayStart()+1 < 2)
    {
      cerr << "WARNING: incorrect number of VEVs supplied to function: skipping calculating LHS matrix." << endl;
    }
  else
    {
      
      double mu = r.displaySusyMu();
      double Bmu = r.displayM3Squared();
      // Calculate B
      double B = Bmu/mu;

      double mHdSq = r.displayMh1Squared();
      double mHuSq = r.displayMh2Squared();

      double v1 = vevs(vevs.displayStart());
      double v2 = vevs(vevs.displayEnd());
      double v = sqrt(v1*v1+v2*v2);
      double tb = v2/v1;
      
      double g1 = r.displayGaugeCoupling(1);
      double g2 = r.displayGaugeCoupling(2);
      double gbar = sqrt(g2*g2+0.6*g1*g1);
      
      
      // Construct relevant partial derivatives
      double df1dv1 = B*mu*v2/(v1*v1)+0.25*gbar*gbar*v1;//mu*mu+mHdSq+0.125*gbar*gbar*(3.0*v1*v1-v2*v2);
      double df1dv2 = -B*mu/v1-0.25*gbar*gbar*v2;//-0.25*gbar*gbar*v1*v2-Bmu;
      double df2dv1 = -B*mu/v2-0.25*gbar*gbar*v1;//df1dv2;
      double df2dv2 = B*mu*v1/(v2*v2)+0.25*gbar*gbar*v2;//mu*mu+mHuSq+0.125*gbar*gbar*(3.0*v2*v2-v1*v1);
      
      // Add in loop corrections if requested
      if (INCLUDE1LPTADPOLES)
	{
	  df1dv1 = df1dv1+(1.0/v1)*(doCalcMSSMDeltaPrime11(r, tb)-doCalcTadpoleMSSMH1(r,tb));//df1dv1+doCalcMSSMDeltaPrime11(r,tb);
	  df1dv2 = df1dv2+doCalcMSSMDeltaPrime12(r, tb)/v1;//df1dv2+doCalcMSSMDeltaPrime12(r,tb);
	  df2dv1 = df2dv1+doCalcMSSMDeltaPrime12(r, tb)/v2;//df2dv1+doCalcMSSMDeltaPrime12(r,tb);
	  df2dv2 = df2dv2+(1.0/v2)*(doCalcMSSMDeltaPrime22(r, tb)-doCalcTadpoleMSSMH2(r,tb));//df2dv2+doCalcMSSMDeltaPrime22(r,tb);
	}

      
      lhsMat(1,1) = df1dv1;
      lhsMat(1,2) = df1dv2;
      lhsMat(2,1) = df2dv1;
      lhsMat(2,2) = df2dv2;
    }
  
  return lhsMat;

}

// A function used in numerically estimating the derivatives of the EW scale parameters
// w.r.t the high scale input parameters. Passed to SOFTSUSY's calcDerivative method
// in the tuning calculation.

// Variables used for getting information into the function
static SoftParsMssm *tempSoftRunner; //< assumed to be at the input scale MX
static double currentMsusy;
static int lowParNum;
static int numTuningPars;
static DoubleVector currentPars(NUMPMSSMPARS); // 5/5/2014 current parameters
static int parChoice;
static void (*currentFTBCs)(SoftParsMssm &, DoubleVector &); // Modified 5/5/2014

double rgeDerivCalc(double x)
{
  // Store original object
  DoubleVector storedObject(tempSoftRunner->display());
  double mx = tempSoftRunner->displayMu();
  DoubleVector storedPars = currentPars;

  double lowVal;

  if (PRINTOUT > 1)
    {
      cout << *tempSoftRunner;
    }

  // 5/5/2014 Update parameter value
  currentPars(parChoice) = x;

  currentFTBCs(*tempSoftRunner, currentPars);

  if (PRINTOUT > 1)
    {
      cout << "Applied boundary conditions at MX." << endl;
      cout << *tempSoftRunner;
    }

  tempSoftRunner->runto(currentMsusy);

  if (PRINTOUT > 1)
    {
      cout << *tempSoftRunner;
    }

  switch(lowParNum)
    {
    case 1:
      {
	lowVal = tempSoftRunner->displaySusyMu();
	if (PRINTOUT > 1) cout << x << " mu(M_{SUSY}) = " << lowVal << endl;
	break;
      }
    case 2:
      {
	lowVal = tempSoftRunner->displayM3Squared()/tempSoftRunner->displaySusyMu();
	if (PRINTOUT > 1) cout << x << " B(M_{SUSY}) = " << lowVal << endl;
	break;
      }
    case 3:
      {
	lowVal = tempSoftRunner->displayMh1Squared();
	if (PRINTOUT > 1) cout << x << " mHdsq(M_{SUSY}) = " << lowVal << endl;
	break;
      }
    case 4:
      {
	lowVal = tempSoftRunner->displayMh2Squared();
	if (PRINTOUT > 1) cout << x << " mHusq(M_{SUSY}) = " << lowVal << endl;
	break;
      }
    case 5:
      {
	lowVal = tempSoftRunner->displaySoftMassSquared(mQl, 3, 3);
	if (PRINTOUT > 1) cout << x << " mqL3sq(M_{SUSY}) = " << lowVal << endl;
	break;
      }
    case 6:
      {
	lowVal = tempSoftRunner->displaySoftMassSquared(mUr, 3, 3);
	if (PRINTOUT > 1) cout << x << " mtRsq(M_{SUSY}) = " << lowVal << endl;
	break;
      }
    case 7:
      {
	lowVal = tempSoftRunner->displaySoftA(UA, 3, 3);
	if (PRINTOUT > 1) cout << x << " At(M_{SUSY}) = " << lowVal << endl;
	break;
      }
    default:
      {
      ostringstream ii;
      ii << "rgeDerivCalc called with unrecognised EW scale parameter choice " << lowParNum << "." << endl;
      throw ii.str();
      }
    }

  // Reset object
  tempSoftRunner->setMu(mx);
  tempSoftRunner->set(storedObject);
  currentPars = storedPars; // Reset parameters

  return lowVal;

}

// A function to calculate the vectors appearing on the RHS of the calculation of the
// derivatives d v_i/ dp_j. Takes as arguments a SoftParsMssm object, assumed to be evaluated
// at the scale M_SUSY, vectors of the parameters to calculate the derivatives wrt. and the vevs,
// an integer labelling which particular parameter to calculate the derivative for, and the scale
// at which that parameter is defined.
// Modified 5/5/2014 - now takes a function ftBCatMX that is used to set the values of the
// parameters, to allow for more general parameter sets.
DoubleVector doCalcRHSTuningVector(SoftParsMssm r, void (*ftBCatMX)(SoftParsMssm &, DoubleVector &), DoubleVector pars, 
				   DoubleVector const & vevs, int i, double mx, bool & hasError)
{

  DoubleVector rhsVec(2);

  if (vevs.displayEnd()-vevs.displayStart()+1 < 2)
    {
      cerr << "WARNING: incorrect number of VEVs supplied to function: skipping calculating RHS vector." << endl;
    }
  else
    {
      double v1 = vevs(vevs.displayStart());
      double v2 = vevs(vevs.displayEnd());
      double tb = v2/v1;

      double mu = r.displaySusyMu();
      double Bmu = r.displayM3Squared();
      // Calculate B
      double B = Bmu/mu;     
 
      // First calculate the elements of the matrix that appears on the RHS in general, being derivatives of the EWSB
      // conditions wrt the EW scale parameters. For the EWSB conditions used in our study, this matrix has the form:
      // 
      //     [ df1/dmu df1/dB df1/dm_Hd^2 df1/dm_Hu^2 df1/dm_Ql^2 df1/dm_uR^2 df1/dA_t ]
      //     [ df2/dmu df2/dB df2/dm_Hd^2 df2/dm_Hu^2 df2/dm_Ql^2 df2/dm_uR^2 df2/dA_t ]
      DoubleMatrix EWderivs(2,7);
      
      double df1dmu = 2*mu-B*tb;//2.0*mu*v1-B*v2;
      double df2dmu = 2*mu-B/tb;//2.0*mu*v2-B*v1;
      
      double df1dmH1Sq = 1.0;//v1;
      double df2dmH1Sq = 0.0;
      
      double df1dmH2Sq = 0.0;
      double df2dmH2Sq = 1.0;//v2;
      
      double df1dB = -mu*tb;//-mu*v2;
      double df2dB = -mu/tb;//-mu*v1;
      
      double df1dmQlSq = 0.0;
      double df2dmQlSq = 0.0;
      
      double df1dmUrSq = 0.0;
      double df2dmUrSq = 0.0;
      
      double df1dAt = 0.0;
      double df2dAt = 0.0;
      
      if (INCLUDE1LPTADPOLES)
	{
	  df1dmu = df1dmu+doCalcMSSMd2DeltaVdMudv1(r, tb);
	  df2dmu = df2dmu+doCalcMSSMd2DeltaVdMudv2(r, tb);
	  
	  df1dmQlSq = df1dmQlSq+doCalcMSSMd2DeltaVdmQlSqdv1(r, tb);
	  df2dmQlSq = df2dmQlSq+doCalcMSSMd2DeltaVdmQlSqdv2(r, tb);
	  
	  df1dmUrSq = df1dmUrSq+doCalcMSSMd2DeltaVdmUrSqdv1(r, tb);
	  df2dmUrSq = df2dmUrSq+doCalcMSSMd2DeltaVdmUrSqdv2(r, tb);
	  
	  df1dAt = df1dAt+doCalcMSSMd2DeltaVdAtdv1(r, tb);
	  df2dAt = df2dAt+doCalcMSSMd2DeltaVdAtdv2(r, tb);
	  
	}    

      //     [ df1/dmu df1/dB df1/dm_Hd^2 df1/dm_Hu^2 df1/dm_Ql^2 df1/dm_uR^2 df1/dA_t ]
      //     [ df2/dmu df2/dB df2/dm_Hd^2 df2/dm_Hu^2 df2/dm_Ql^2 df2/dm_uR^2 df2/dA_t ]      

      EWderivs(1,1) = df1dmu;
      EWderivs(2,1) = df2dmu;
      EWderivs(1,2) = df1dB;
      EWderivs(2,2) = df2dB;
      EWderivs(1,3) = df1dmH1Sq;
      EWderivs(2,3) = df2dmH1Sq;
      EWderivs(1,4) = df1dmH2Sq;
      EWderivs(2,4) = df2dmH2Sq;
      EWderivs(1,5) = df1dmQlSq;
      EWderivs(2,5) = df2dmQlSq;
      EWderivs(1,6) = df1dmUrSq;
      EWderivs(2,6) = df2dmUrSq;
      EWderivs(1,7) = df1dAt;
      EWderivs(2,7) = df2dAt;

      // Then calculate the vector that stores the derivatives of the EW scale parameters wrt the high scale parameter
      // of interest. In our study this has the form
      // [ dmu/dp dB/dp dm_Hd^2/dp dm_Hu^2/dp dm_Ql^2/dp dm_uR^2/dp dA_t/dp]^T.
      // For this version of the function, this is done numerically by running between M_SUSY and MX to estimate the
      // derivatives.

      // Get the scale M_SUSY - note that the model is assumed to be provided at this scale
      double ms = r.displayMu();      

      // Finite difference used to estimate derivative
      double epsilon = 1.0e-5;
      double h;
      double temp;

      DoubleVector paramDerivs(7);

      // Save initial values of all EW parameters.
      double mH1Sq = r.displayMh1Squared();
      double mH2Sq = r.displayMh2Squared();
      double mQlSq = r.displaySoftMassSquared(mQl, 3, 3);
      double mUrSq = r.displaySoftMassSquared(mUr, 3, 3);
      double At = r.displaySoftA(UA, 3, 3);

      double mu_final, B_final, mH1Sq_final, mH2Sq_final, mQlSq_final, mUrSq_final, At_final;

      // Run to MX
      r.runto(mx);

      currentPars = pars;

      if (PRINTOUT > 1) cout << currentPars;

      tempSoftRunner = &r;
      parChoice = i;
      currentMsusy = ms;
      currentFTBCs = ftBCatMX; // Modified 5/5/2014

      h = epsilon*fabs(pars(i));
      if (h == 0.0) h = epsilon;
      temp = pars(i);
      pars(i) = temp + h;
      h = pars(i) - temp;

 
      // The particular method used will depend on the number of parameters provided. If
      // the number of parameters = 5, we have a SUGRA model. Otherwise if the number of 
      // parameters = 25 we have a pMSSM model.
      if (pars.displayEnd()-pars.displayStart()+1 == 5) // sugra case
	{
	  numTuningPars = 5;
	}
      else if (pars.displayEnd()-pars.displayStart()+1 == 25) // pMSSM case
	{
	  // In the pMSSM, first and second generation trilinear elements must be zero.
	  r.setTrilinearElement(UA, 1, 1, 0.0);
	  r.setTrilinearElement(UA, 2, 2, 0.0);
	  r.setTrilinearElement(DA, 1, 1, 0.0);
	  r.setTrilinearElement(DA, 2, 2, 0.0);
	  r.setTrilinearElement(EA, 1, 1, 0.0);
	  r.setTrilinearElement(EA, 2, 2, 0.0);
	  
	  // Likewise first and second generation Yukawas.
	  r.setYukawaElement(YU, 1, 1, 0.0);
	  r.setYukawaElement(YU, 2, 2, 0.0);
	  r.setYukawaElement(YD, 1, 1, 0.0);
	  r.setYukawaElement(YD, 2, 2, 0.0);
	  r.setYukawaElement(YE, 1, 1, 0.0);
	  r.setYukawaElement(YE, 2, 2, 0.0);

	  numTuningPars = 20;
	}
      /* Commented out 5/5/2014
      else
	{
	  ostringstream ii;
	  ii << "doCalcRHSTuningVector called with unrecognised MSSM model." << endl;
	  throw ii.str();
	}
	End commented out 5/5/2014 */
 

      DoubleVector errVec(7);
      double err;

      if (fabs(ms-mx) < TOLERANCE)
	{
	  // No need to waste time with running
	  if (numTuningPars == 5)
	    {
	      // [ dmu/dp dB/dp dm_Hd^2/dp dm_Hu^2/dp dm_Ql^2/dp dm_uR^2/dp dA_t/dp]^T.
	      switch(i)
		{
		case 1: // m0sq
		  {
		    paramDerivs(1) = 0.0;
		    paramDerivs(2) = 0.0;
		    paramDerivs(3) = -1.0;
		    paramDerivs(4) = -1.0;
		    paramDerivs(5) = -1.0;
		    paramDerivs(6) = -1.0;
		    paramDerivs(7) = 0.0;
		    break;
		  }
		case 2: // M1/2
		  {
		    paramDerivs(1) = 0.0;
		    paramDerivs(2) = 0.0;
		    paramDerivs(3) = 0.0;
		    paramDerivs(4) = 0.0;
		    paramDerivs(5) = 0.0;
		    paramDerivs(6) = 0.0;
		    paramDerivs(7) = 0.0;
		    break;
		  }
		case 3: // mu0
		  {
		    paramDerivs(1) = -1.0;
		    paramDerivs(2) = 0.0;
		    paramDerivs(3) = 0.0;
		    paramDerivs(4) = 0.0;
		    paramDerivs(5) = 0.0;
		    paramDerivs(6) = 0.0;
		    paramDerivs(7) = 0.0;
		    break;
		  }
		case 4: // B0
		  {
		    paramDerivs(1) = 0.0;
		    paramDerivs(2) = -1.0;
		    paramDerivs(3) = 0.0;
		    paramDerivs(4) = 0.0;
		    paramDerivs(5) = 0.0;
		    paramDerivs(6) = 0.0;
		    paramDerivs(7) = 0.0;
		    break;
		  }
		case 5: // A0
		  {
		    paramDerivs(1) = 0.0;
		    paramDerivs(2) = 0.0;
		    paramDerivs(3) = 0.0;
		    paramDerivs(4) = 0.0;
		    paramDerivs(5) = 0.0;
		    paramDerivs(6) = 0.0;
		    paramDerivs(7) = -1.0;
		    break;
		  }
		default:
		  {
		    paramDerivs(1) = 0.0;
		    paramDerivs(2) = 0.0;
		    paramDerivs(3) = 0.0;
		    paramDerivs(4) = 0.0;
		    paramDerivs(5) = 0.0;
		    paramDerivs(6) = 0.0;
		    paramDerivs(7) = 0.0;
		    break;
		  }
		}
	    }
	  else
	    {
	      switch(i)
		{
		case 4: // At
		  {
		    paramDerivs(1) = 0.0;
		    paramDerivs(2) = 0.0;
		    paramDerivs(3) = 0.0;
		    paramDerivs(4) = 0.0;
		    paramDerivs(5) = 0.0;
		    paramDerivs(6) = 0.0;
		    paramDerivs(7) = -1.0;
		    break;
		  }
		case 7: // mHdsq
		  {
		    paramDerivs(1) = 0.0;
		    paramDerivs(2) = 0.0;
		    paramDerivs(3) = -1.0;
		    paramDerivs(4) = 0.0;
		    paramDerivs(5) = 0.0;
		    paramDerivs(6) = 0.0;
		    paramDerivs(7) = 0.0;
		    break;
		  }
		case 8: // mHusq
		  {
		    paramDerivs(1) = 0.0;
		    paramDerivs(2) = 0.0;
		    paramDerivs(3) = 0.0;
		    paramDerivs(4) = -1.0;
		    paramDerivs(5) = 0.0;
		    paramDerivs(6) = 0.0;
		    paramDerivs(7) = 0.0;
		    break;
		  }
		case 9: // mu
		  {
		    paramDerivs(1) = -1.0;
		    paramDerivs(2) = 0.0;
		    paramDerivs(3) = 0.0;
		    paramDerivs(4) = 0.0;
		    paramDerivs(5) = 0.0;
		    paramDerivs(6) = 0.0;
		    paramDerivs(7) = 0.0;
		    break;
		  }
		case 10: // B
		  {
		    paramDerivs(1) = 0.0;
		    paramDerivs(2) = -1.0;
		    paramDerivs(3) = 0.0;
		    paramDerivs(4) = 0.0;
		    paramDerivs(5) = 0.0;
		    paramDerivs(6) = 0.0;
		    paramDerivs(7) = 0.0;
		    break;
		  }
		case 18: // mQlsq
		  {
		    paramDerivs(1) = 0.0;
		    paramDerivs(2) = 0.0;
		    paramDerivs(3) = 0.0;
		    paramDerivs(4) = 0.0;
		    paramDerivs(5) = -1.0;
		    paramDerivs(6) = 0.0;
		    paramDerivs(7) = 0.0;
		    break;
		  }
		case 19: // muRsq
		  {
		    paramDerivs(1) = 0.0;
		    paramDerivs(2) = 0.0;
		    paramDerivs(3) = 0.0;
		    paramDerivs(4) = 0.0;
		    paramDerivs(5) = 0.0;
		    paramDerivs(6) = -1.0;
		    paramDerivs(7) = 0.0;
		    break;
		  }
		default:
		  {
		    paramDerivs(1) = 0.0;
		    paramDerivs(2) = 0.0;
		    paramDerivs(3) = 0.0;
		    paramDerivs(4) = 0.0;
		    paramDerivs(5) = 0.0;
		    paramDerivs(6) = 0.0;
		    paramDerivs(7) = 0.0;
		    break;
		  }
		}
	    }
	}
      else
	{
	  for (int j = 1; j <= 7; j++)
	    {
	      lowParNum = j;
	      paramDerivs(j) = -calcDerivative(rgeDerivCalc, temp, h, &err); // < NB changed pars(i) to temp
	      errVec(j) = err;
	      if (PRINTOUT > 1)
		{
		  cout << "# j = " << j << ", derivative = " << -paramDerivs(j) << ", error = " << errVec(j) << "." << endl;
		}
	    }
	}
      //     cout << paramDerivs;
      // Multiply the two together and return the result (note the additional minus sign already applied!).
      rhsVec = EWderivs*paramDerivs;

      // Check for errors - if one derivative failed then all of the elements of the vector will be unreliable
      for (int j = 1; j <= 7; j++)
	{
	  // Possible bug here if both the derivative and the error are very close
	  // to zero. If both are very small, for now we will just assume that
	  // the value is correct and is essentially zero.
	  if (fabs(errVec(j)) > 1.0e-8 && fabs(errVec(j)/paramDerivs(j)) > 1.0)
	    {
	      hasError = true;
	    }
	}
      if (hasError)
	{
	  for (int i = rhsVec.displayStart(); i <= rhsVec.displayEnd(); i++)
	    {
	      rhsVec(i) = numberOfTheBeast;
	    }
	}

    }

  return rhsVec;
}

// A function to calculate d log(M_Z^2)/d log(p) using the value of the parameter p
// and the already calculated derivatives d v_i/ d p_j. 
double doCalcdLogMzSqdLogParam(SoftParsMssm r, double p, DoubleVector const & vevs, DoubleVector const & dVevsdp)
{

  double deriv = 0.0;

  if ((vevs.displayEnd()-vevs.displayStart()+1 < 2) || (dVevsdp.displayEnd()-dVevsdp.displayStart()+1 < 2))
    {
      cerr << "WARNING: incorrect number of VEVs supplied to function: skipping calculating derivative." << endl;
    }
  else
    {
      double g1 = r.displayGaugeCoupling(1);
      double g2 = r.displayGaugeCoupling(2);
      double gbar = sqrt(g2*g2+0.6*g1*g1);
      
      double v1 = vevs(vevs.displayStart());
      double v2 = vevs(vevs.displayEnd());
      double v = sqrt(v1*v1+v2*v2);
      
      deriv = (2.0*p/(v*v))*(v1*dVevsdp(dVevsdp.displayStart())+v2*dVevsdp(dVevsdp.displayEnd()));
    }

  return deriv;

}

// MSSM_EWSBConditioni, i = 1, 2, calculates the value of the 
// EWSB condition for m_H_i^2.
// Inputs:
//     SoftParsMssm const & r = MSSM model object
//     DoubleVector const & vevs = the values of v1 and v2
double MSSM_EWSBCondition1(SoftParsMssm const & r)
{
  // bool INCLUDE1LPTADPOLES = true;

  double mu = r.displaySusyMu();
  // NB SOFTSUSY convention B\mu = m_3^2
  double Bmu = r.displayM3Squared();
  double mH1Sq = r.displayMh1Squared();
  double mH2Sq = r.displayMh2Squared();

  double g1 = r.displayGaugeCoupling(1);
  double g2 = r.displayGaugeCoupling(2);
  double gbar = sqrt(g2*g2+0.6*g1*g1);

  double f1;

  double tb = r.displayTanb();
  double v = r.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;
  
  if (!INCLUDE1LPTADPOLES) // Tree level
    {
      f1 = (mu*mu+mH1Sq)-Bmu*v2/v1+0.125*gbar*gbar*(v1*v1-v2*v2);
    }
  else // One-loop tadpoles included
    {
      f1 = (mu*mu+mH1Sq)-Bmu*v2/v1+0.125*gbar*gbar*(v1*v1-v2*v2)+doCalcTadpoleMSSMH1(r, tb);
    }
  
  
  return f1;
}

// Tree level checked against SOFTSUSY with model 2403883, seems 
// correct.
double MSSM_EWSBCondition2(SoftParsMssm const & r)
{
  // bool INCLUDE1LPTADPOLES = true;

  double mu = r.displaySusyMu();
  // NB SOFTSUSY convention B\mu = m_3^2
  double Bmu = r.displayM3Squared();
  double mH1Sq = r.displayMh1Squared();
  double mH2Sq = r.displayMh2Squared();

  double g1 = r.displayGaugeCoupling(1);
  double g2 = r.displayGaugeCoupling(2);
  double gbar = sqrt(g2*g2+0.6*g1*g1);

  double f2;

  double tb = r.displayTanb();
  double v = r.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;
  
  
  if (!INCLUDE1LPTADPOLES) // Tree level
    {
      f2 = (mu*mu+mH2Sq)-Bmu*v1/v2-0.125*gbar*gbar*(v1*v1-v2*v2);
    }
  else // One-loop tadpoles included
    {
      f2 = (mu*mu+mH2Sq)-Bmu*v1/v2-0.125*gbar*gbar*(v1*v1-v2*v2)+doCalcTadpoleMSSMH2(r, tb);
    }
  
  
  return f2;
    
}


static MssmSoftsusy *tempsusy;
void MSSM_EWSBCondition_NewtonFunc(DoubleVector const & guess, DoubleVector & fVals)
{
  //cout << "Calling the Newton function..." << endl;

  cout << "Current guesses for VEVs: v1 = " << guess(1) << ", v2 = " << guess(2) << endl;

  // Guess has elements v1 and v2
  double v1 = guess(1);
  double v2 = guess(2);

  tempsusy->setHvev(sqrt(sqr(v1)+sqr(v2)));
  tempsusy->setTanb(v2/v1);

  fVals(1) = MSSM_EWSBCondition1(*tempsusy);
  fVals(2) = MSSM_EWSBCondition2(*tempsusy);

  cout << "Current values of EWSB conditions: f1 = " << fVals(1) << ", f2 = " << fVals(2) << endl;
}


// Implement the EWSB conditions by varying the values of m_Hu^2 and m_Hd^2 at MX.

// Scale factors needed for Newton's method implementation in SOFTSUSY
static double mHuSqScaleFactor = 1.0e6;
static double mHdSqScaleFactor = 1.0e6;
static double MsusyScaleFactor = 1.0e3;

bool MSSM_ImplementEWSBConstraints_SoftMasses(SoftParsMssm r, double mx, double ms, bool useMxEqualsMs,
				   DoubleVector & updatedSoln, double tol)
{
  // bool INCLUDE1LPTADPOLES = true;

  // We have to use a shooting method to do this in general if MX is different from M_{SUSY}. 
  // We should start with an object given at the scale MX.
  double q = r.displayMu();

  bool hasProblem = false;

  if (fabs(q-mx) > TOLERANCE)
    {
      ostringstream ii;
      ii << "WARNING: trying to calculate required soft masses at inappropriate scale q = " << q << " GeV,\n"
	 << "         instead of at MX = " << mx << " GeV: exiting." << endl;
      throw ii.str();
    }
  else
    {
      // If MX = M_{SUSY} there is no need to use the shooting approach; just
      // rearrange the EWSB conditions to solve for m_Hu^2 and m_Hd^2.
      if (useMxEqualsMs)
	{
	  if (fabs(mx-ms) > TOLERANCE || fabs(q-ms) > TOLERANCE)
	    {
	      cerr << "WARNING: MX = M_{SUSY} requested but given values do not agree. Assuming M_{SUSY} = " << ms << "." << endl;
	      r.runto(ms);
	    }

	  // Solve for the updated soft Higgs masses

	  double mu = r.displaySusyMu();
	  // NB SOFTSUSY convention B\mu = m_3^2
	  double Bmu = r.displayM3Squared();
	  double mH1Sq, mH2Sq;
	  double g1 = r.displayGaugeCoupling(1);
	  double g2 = r.displayGaugeCoupling(2);
	  double gbar = sqrt(g2*g2+0.6*g1*g1);
	  
	  double v = r.displayHvev();
	  double tb = r.displayTanb();
	  double v1 = v/sqrt(1.0+tb*tb);
	  double v2 = v1*tb;  
	  
	  if (!INCLUDE1LPTADPOLES)
	    {
	      mH1Sq = -(mu*mu-Bmu*tb+0.125*gbar*gbar*(v1*v1-v2*v2));
	      mH2Sq = -(mu*mu-Bmu/tb-0.125*gbar*gbar*(v1*v1-v2*v2));
	    }
	  else
	    {
	      mH1Sq = -(mu*mu-Bmu*tb+0.125*gbar*gbar*(v1*v1-v2*v2)+doCalcTadpoleMSSMH1(r, tb));
	      mH2Sq = -(mu*mu-Bmu/tb-0.125*gbar*gbar*(v1*v1-v2*v2)+doCalcTadpoleMSSMH2(r, tb));
	    }
	  updatedSoln(1) = mH1Sq;
	  updatedSoln(2) = mH2Sq;
	  updatedSoln(3) = ms;
	}
      else
	{
	  // Otherwise we have to use a shooting method to determine the values
	  // of m_Hu^2 and m_Hd^2 at MX that satisfy the EWSB conditions at M_{SUSY}.

	  // Guess initial values of m_Hu^2 and m_Hd^2 at MX that solve the EWSB conditions,
	  // and guess the initial scale M_{SUSY}.

	  // The initial guess for M_{SUSY} is pretty obvious - it is the provided value ms.
	  double mSusy_guess = ms;

	  MsusyScaleFactor = ms;

	  // To guess the appropriate initial values of m_Hu^2(MX) and m_Hd^2(MX), we use 
	  // our approximate Taylor series solution to the RGEs. This should be reasonable
	  // for small values of log(MX/M_{SUSY}). Note that this assumes only non-zero
	  // couplings are third generation.
	  int nLps = r.displayLoops();
	  if (nLps != 1 && nLps != 2)
	    {
	      cerr << "WARNING: requested RG running at " << nLps << " loop order: can only do 1 or 2 loop running:" << endl;
	      cerr << "         assuming 2 loop running." << endl;
	      nLps = 2;
	    }

	  double t_run = log(mx/ms);

	  double mHuSq_guess, mHdSq_guess;

	  // If t is too large, this is an unreliable estimate, so just use a default value.
	  int TLIMIT = 10;
	  if (t_run >= TLIMIT)
	    {
	      mHdSqScaleFactor = r.displayMh1Squared();
	      mHuSqScaleFactor = r.displayMh2Squared();

	      mHuSq_guess = mHuSqScaleFactor;
	      mHdSq_guess = mHdSqScaleFactor;

	    }
	  else
	    {
	      // The Taylor series expansion is about the solution at M_{SUSY} initially.
	      SoftParsMssm s = r; // copy object to avoid numerical errors
	      s.runto(ms);
	      
	      double mHuSqInit, mHdSqInit;
	      
	      double mu = s.displaySusyMu();
	      // NB SOFTSUSY convention B\mu = m_3^2
	      double Bmu = s.displayM3Squared();
	      double mH1Sq, mH2Sq;
	      double g1 = s.displayGaugeCoupling(1);
	      double g2 = s.displayGaugeCoupling(2);
	      double gbar = sqrt(g2*g2+0.6*g1*g1);
	      
	      double v = s.displayHvev();
	      double tb = s.displayTanb();
	      double v1 = v/sqrt(1.0+tb*tb);
	      double v2 = v1*tb;  
	      
	      if (!INCLUDE1LPTADPOLES)
		{
		  mHdSqInit = -(mu*mu-Bmu*tb+0.125*gbar*gbar*(v1*v1-v2*v2));
		  mHuSqInit = -(mu*mu-Bmu/tb-0.125*gbar*gbar*(v1*v1-v2*v2));
		}
	      else
		{
		  mHdSqInit = -(mu*mu-Bmu*tb+0.125*gbar*gbar*(v1*v1-v2*v2)+doCalcTadpoleMSSMH1(s, tb));
		  mHuSqInit = -(mu*mu-Bmu/tb-0.125*gbar*gbar*(v1*v1-v2*v2)+doCalcTadpoleMSSMH2(s, tb));
		}
	      
	      // Construct required coefficients in Taylor series.
	      double mHdSqLogCoeff = doCalcMh1SquaredLogCoeff(s, nLps);
	      double mHdSqLogSqCoeff = doCalcMh1SquaredLogSqCoeff(s, 1);
	      double mHuSqLogCoeff = doCalcMh2SquaredLogCoeff(s, nLps);
	      double mHuSqLogSqCoeff = doCalcMh2SquaredLogSqCoeff(s, 1);
	      
	      mHuSq_guess = mHuSqInit + t_run*mHuSqLogCoeff + sqr(t_run)*mHuSqLogSqCoeff;
	      mHdSq_guess = mHdSqInit + t_run*mHdSqLogCoeff + sqr(t_run)*mHdSqLogSqCoeff;
	      
	      mHuSqScaleFactor = mHuSq_guess;
	      mHdSqScaleFactor = mHdSq_guess;

	    }

	  // Then shoot using Newton's method to try to get the actual solution. The
	  // BCs to satisfy are: f_1(m_Hu^2(M_{SUSY}), m_Hd^2(M_{SUSY}) = 0, 
	  // f_2(m_Hu^2(M_{SUSY}), m_Hd^2(M_{SUSY}) = 0, and 
	  // M_{SUSY} = \sqrt{m_{\tilde{t}_1}(M_{SUSY})m_{\tilde{t}_2}(M_{SUSY})}, 
	  // and all other values given at MX are fixed. Additionally all of
	  // the VEVs v_1 and v_2 are fixed to their values at M_{SUSY}. We want
	  // to try and implement a globally convergent Newton's method for this,
	  // because the approximations that we have used to obtain the initial
	  // guesses are not guaranteed to be good guesses, especially if log(MX/M_{SUSY}) is large.
	  DoubleVector solEstimate(3);

	  solEstimate(1) = mHdSq_guess/mHdSqScaleFactor;
	  solEstimate(2) = mHuSq_guess/mHuSqScaleFactor;
	  solEstimate(3) = mSusy_guess/MsusyScaleFactor;

	  hasProblem = MSSM_EWSB_NewtonShooter(r, solEstimate, tol);

	  updatedSoln(1) = solEstimate(1)*mHdSqScaleFactor;
	  updatedSoln(2) = solEstimate(2)*mHuSqScaleFactor;
	  updatedSoln(3) = solEstimate(3)*MsusyScaleFactor;

	}

    }
  return hasProblem;
}

static double muScaleFactor = 1.0e2;
static double B0ScaleFactor = 2.5e3;
static double f1ScaleFactor = 1.0e6;
static double f2ScaleFactor = 1.0e6;
bool MSSM_ImplementEWSBConstraints_SUGRA(SoftParsMssm r, double mx, double ms,
				   DoubleVector & updatedSoln, double tol)
{
  // We have to use a shooting method to do this in general if MX is different from M_{SUSY}. 
  // We should start with an object given at the scale MX.
  double q = r.displayMu();

  bool hasProblem = false;

  if (fabs(q-mx) > TOLERANCE)
    {
      ostringstream ii;
      ii << "WARNING: trying to calculate required soft masses at inappropriate scale q = " << q << " GeV,\n"
	 << "         instead of at MX = " << mx << " GeV: exiting." << endl;
      throw ii.str();
    }
  else
    {
      // We need to use a shooting method to determine the values of mu and B0 at MX
      // that will solve the EWSB conditions at M_{SUSY}.

      // Get initial guess. This is done using the current values of mu and B0.
      // Likewise for the guess for M_{SUSY}.
      muScaleFactor = r.displaySusyMu();
      B0ScaleFactor = r.displayM3Squared()/r.displaySusyMu();
      MsusyScaleFactor = ms;

      double mu_guess = muScaleFactor;
      double B0_guess = B0ScaleFactor;
      double mSusy_guess = MsusyScaleFactor;

      // Get initial guess for EWSB condition scale factors.
      SoftParsMssm s = r;
      s.runto(ms);
      f1ScaleFactor = MSSM_EWSBCondition1(s);
      f2ScaleFactor = MSSM_EWSBCondition2(s);

      DoubleVector solEstimate(3);
      
      solEstimate(1) = mu_guess/muScaleFactor;
      solEstimate(2) = B0_guess/B0ScaleFactor;
      solEstimate(3) = mSusy_guess/MsusyScaleFactor;
      
      hasProblem = SUGRA_EWSB_NewtonShooter(r, solEstimate, tol);

      updatedSoln(1) = solEstimate(1)*muScaleFactor;
      updatedSoln(2) = solEstimate(2)*B0ScaleFactor;
      updatedSoln(3) = solEstimate(3)*MsusyScaleFactor;

    }

  return hasProblem;
}

double MSSM_Msusy_Cond(SoftParsMssm r, double ms)
{
  DoubleVector mstop(2), mstopsq(2);

  double tb = r.displayTanb();

  physical_MSSM(r, mstop, mstopsq, tb);

  double f = (ms/sqrt(mstop(1)*mstop(2)))-1.0;

  return f;

}

static SoftParsMssm *tempsoftmssm;
void SUGRA_EWSB_Shooter_Functions(DoubleVector const & parVals, DoubleVector & f)
{

  SoftParsMssm s = *tempsoftmssm;

  s.setSusyMu(parVals(1)*muScaleFactor);
  s.setM3Squared(parVals(2)*B0ScaleFactor*parVals(1)*muScaleFactor);

  s.runto(parVals(3)*MsusyScaleFactor);

  f(1) = MSSM_EWSBCondition1(s)/f1ScaleFactor;
  f(2) = MSSM_EWSBCondition2(s)/f2ScaleFactor;
  f(3) = MSSM_Msusy_Cond(s, parVals(3)*MsusyScaleFactor);

}

void MSSM_EWSB_Shooter_Functions(DoubleVector const & parVals, DoubleVector & f)
{

  SoftParsMssm s = *tempsoftmssm;

  s.setMh1Squared(parVals(1)*mHdSqScaleFactor);
  s.setMh2Squared(parVals(2)*mHuSqScaleFactor);

  s.runto(parVals(3)*MsusyScaleFactor);

  f(1) = MSSM_EWSBCondition1(s)/mHdSqScaleFactor;
  f(2) = MSSM_EWSBCondition2(s)/mHuSqScaleFactor;
  f(3) = MSSM_Msusy_Cond(s, parVals(3)*MsusyScaleFactor);

}

bool SUGRA_EWSB_NewtonShooter(SoftParsMssm const & r, DoubleVector & estimate, double tol)
{
  // Copy object to actually do the running on
  SoftParsMssm s = r;

  // Set initial guess values
  double ms_guess = estimate(3)*MsusyScaleFactor;

  s.setSusyMu(estimate(1)*muScaleFactor);
  s.setM3Squared(estimate(2)*B0ScaleFactor*estimate(1)*muScaleFactor);

  tempsoftmssm = &s;

  bool hasProblem = newt(estimate, SUGRA_EWSB_Shooter_Functions);

  s.setSusyMu(estimate(1)*muScaleFactor);
  s.setM3Squared(estimate(1)*muScaleFactor*estimate(2)*B0ScaleFactor);

  s.runto(estimate(3)*MsusyScaleFactor);

  DoubleVector mstop(2), mstopsq(2);
  double tb = s.displayTanb();
  physical_MSSM(s, mstop, mstopsq, tb);

  // Have we actually converged or is newt just being ridiculous?
  double f1 = MSSM_EWSBCondition1(s);
  double f2 = MSSM_EWSBCondition2(s);
  double f3 = estimate(3)*MsusyScaleFactor/sqrt(mstop(1)*mstop(2))-1.0;

  if (fabs(f1) > tol || fabs(f2) > tol || fabs(f3) > tol)
    {
      hasProblem = true;
    }

  return hasProblem;


}

// Function for doing shooting in implementing the EWSB conditions. The SoftParsMssm object
// contains the values of the parameters and VEVs that are requested at the scale MX. The
// DoubleVector estimate initially contains the guesses for m_Hd^2, m_Hu^2 and M_{SUSY} (in that
// order), and on return contains the estimated values for the solutions.
bool MSSM_EWSB_NewtonShooter(SoftParsMssm const & r, DoubleVector & estimate, double tol)
{
  // Copy object to actually do the running on
  SoftParsMssm s = r;

  // Set initial guess values
  double ms_guess = estimate(3)*MsusyScaleFactor;

  s.setMh1Squared(estimate(1)*mHdSqScaleFactor);
  s.setMh2Squared(estimate(2)*mHuSqScaleFactor);

  DoubleVector fVals(3);

  tempsoftmssm = &s;

  // Now shoot to try to get solutions. Use globally convergent
  // Newton's method as the root finder.
  bool hasProblem = newt(estimate, MSSM_EWSB_Shooter_Functions);

  // cout << "Solution estimate: " << endl;
  // cout << "m_Hd^2 = " << estimate(1)*mHdSqScaleFactor << endl;
  // cout << "m_Hu^2 = " << estimate(2)*mHuSqScaleFactor << endl;
  // cout << "M_{SUSY} = " << estimate(3)*MsusyScaleFactor << endl;
  // cout << "m_Hd^2 scale = " << mHdSqScaleFactor << endl;
  // cout << "m_Hu^2 scale = " << mHuSqScaleFactor << endl;
  // cout << "M_{SUSY} scale = " << MsusyScaleFactor << endl;

  s.setMh1Squared(estimate(1)*mHdSqScaleFactor);
  s.setMh2Squared(estimate(2)*mHuSqScaleFactor);

  s.runto(estimate(3)*MsusyScaleFactor);

  // cout << "f1 = " << MSSM_EWSBCondition1(s) << endl;
  // cout << "f2 = " << MSSM_EWSBCondition2(s) << endl;

  DoubleVector mstop(2), mstopsq(2);
  double tb = s.displayTanb();
  physical_MSSM(s, mstop, mstopsq, tb);

  // cout << "M_{SUSY}/sqrt(mstop(1)*mstop(2)) = " << estimate(3)*MsusyScaleFactor/sqrt(mstop(1)*mstop(2)) << endl;

  // Have we actually converged or is newt just being ridiculous?
  double f1 = MSSM_EWSBCondition1(s);
  double f2 = MSSM_EWSBCondition2(s);
  double f3 = estimate(3)*MsusyScaleFactor/sqrt(mstop(1)*mstop(2))-1.0;

  // cout << "f1/mHdSq Scale Factor = " << f1/mHdSqScaleFactor << endl;
  // cout << "f2/mHdSq Scale Factor = " << f1/mHuSqScaleFactor << endl;
  // cout << "f3 = " << f3 << endl;

  if (fabs(f1) > tol || fabs(f2) > tol || fabs(f3) > tol)
    {
      hasProblem = true;
    }

  return hasProblem;

}

DoubleVector MSSM_EWSBNewtonSolver(SoftParsMssm r, double tol, int maxIters, int & sing)
{

  // Extract initial guess for soln
  double v = r.displayHvev();
  double tb = r.displayTanb();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;

  // Iterate using Newton's method to obtain new solution
  double v1old, v2old;
  double v1new, v2new;
  double f1old, f2old;
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


  for (int i = 0; i < maxIters; i++)
    {
      r.setHvev(vold);
      r.setTanb(tbold);

      // Evaluate elements of Jacobian and residual
      f1old = MSSM_EWSBCondition1(r);
      f2old = MSSM_EWSBCondition2(r);

      // Jacobian elements have been checked numerically elsewhere;
      // here we use the analytic formulae.
      df1dv1 = Bmu*v2old/(v1old*v1old)+0.25*gbar*gbar*v1old;
      df1dv2 = -Bmu/v1old-0.25*gbar*gbar*v2old;
      df2dv1 = -Bmu/v2old-0.25*gbar*gbar*v1old;
      df2dv2 = Bmu*v1old/(v2old*v2old)+0.25*gbar*gbar*v2old;
      if (INCLUDE1LPTADPOLES)
	{
	  df1dv1 = df1dv1+(1.0/v1old)*(doCalcMSSMDeltaPrime11(r, tbold)-doCalcTadpoleMSSMH1(r,tbold));
	  df1dv2 = df1dv2+doCalcMSSMDeltaPrime12(r, tbold)/v1old;
	  df2dv1 = df2dv1+doCalcMSSMDeltaPrime12(r, tbold)/v2old;
	  df2dv2 = df2dv2+(1.0/v2old)*(doCalcMSSMDeltaPrime22(r, tbold)-doCalcTadpoleMSSMH2(r,tbold));
	}
      detF = df1dv1*df2dv2-df1dv2*df2dv1;

      if (fabs(detF) < EPSTOL)
	{
	  cerr << "Warning: singularity encountered in estimating EWSB conditions solution." << endl;
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
      cerr << "Warning: max. iterations reached before solution converged to required accuracy." << endl;
      soln(1) = v1new;
      soln(2) = v2new;
    }

  return soln;

}


// Calculates the DRbar' masses of the stops in the MSSM.
// Inputs:
//     SoftParsMssm r = MSSM model
//     DoubleVector & mstop = masses of stops
//     DoubleVector & mstop sq = squared masses of stops
//     double tb = value of tan(beta)
// Checked against SOFTSUSY results for model 2403883, seem correct.
void physical_MSSM(SoftParsMssm r, DoubleVector & mstop, DoubleVector & mstopsq, double tb)
{
  bool speak = false;

  double mu = r.displaySusyMu();
  double mQlSq = r.displaySoftMassSquared(mQl, 3, 3);
  double mUrSq = r.displaySoftMassSquared(mUr, 3, 3);

  //cout << "mu = " << mu << endl;
  //cout << "mQlSq = " << mQlSq << endl;
  //cout << "mUrSq = " << mUrSq << endl;

  double v = r.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;

  double g1 = r.displayGaugeCoupling(1);
  double g2 = r.displayGaugeCoupling(2);
  double gbar = sqrt(g2*g2+0.6*g1*g1);

  double yt = r.displayYukawaElement(YU, 3, 3);
  double mt = yt*v2/sqrt(2.0);

  double At = r.displaySoftA(UA, 3, 3);

  //cout << "tb = " << tb << endl;
  // cout << "At = " << At << endl;

  double Xt = At-mu/tb;

  //  cout << "Xt = " << Xt << endl;

  mstopsq(1) = 0.5*(mQlSq+mUrSq+0.125*gbar*gbar*(v1*v1-v2*v2)+2.0*mt*mt
		    -sqrt(sqr(mQlSq-mUrSq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2))
			  +4.0*mt*mt*Xt*Xt));
  mstopsq(2) = 0.5*(mQlSq+mUrSq+0.125*gbar*gbar*(v1*v1-v2*v2)+2.0*mt*mt
		    +sqrt(sqr(mQlSq-mUrSq+0.125*(g2*g2-g1*g1)*(v1*v1-v2*v2))
			  +4.0*mt*mt*Xt*Xt));


  if (mstopsq(1) < 0.0 || mstopsq(2) < 0.0)
    {
      cerr << "Warning: tachyonic stop masses." << endl;
      cerr << "m_stop_1^2 = " << mstopsq(1) << " GeV^2." << endl;
      cerr << "m_stop_2^2 = " << mstopsq(2) << " GeV^2." << endl;
    }

  mstop(1) = sqrt(fabs(mstopsq(1)));
  mstop(2) = sqrt(fabs(mstopsq(2)));

  if (speak)
    {
      cout << "# m_stop_1 = " << mstop(1) << " GeV." << endl;
      cout << "# m_stop_2 = " << mstop(2) << " GeV." << endl;
    }
}

// Function needed to evaluate tadpoles.
// Inputs:
//     double mSq = squared mass m^2
//     double Q = renormalisation scale Q
double a0MSSM(double mSq, double Q)
{
  return mSq*(log(mSq/(Q*Q))-1.0);
}

// Calculates the tadpole corrections to the EWSB conditions. Note
// returns (1/v_i)\partial \Delta V/\partial v_i (i.e. opposite sign
// to equivalent E_6SSM methods).
// Inputs:
//     SoftParsMssm r = MSSM model
//     double tb = tan(beta)
double doCalcTadpoleMSSMH1(SoftParsMssm r, double tb)
{
  bool speak = false;

  double v = r.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);

  DoubleVector mstop(2), mstopsq(2);

  physical_MSSM(r, mstop, mstopsq, tb);

  if (mstopsq(1) < 0.0 || mstopsq(2) < 0.0)
    {
      mstopsq(1) = fabs(mstopsq(1));
      mstopsq(2) = fabs(mstopsq(2));
    }

  double q = r.displayMu();

  double a0mstop1 = a0MSSM(mstopsq(1), q);
  double dmstop1sqdv1 = doCalcMSSMdmstop1sqdv1(r, tb);
  double a0mstop2 = a0MSSM(mstopsq(2), q);
  double dmstop2sqdv1 = doCalcMSSMdmstop2sqdv1(r, tb);

  double delta1tp = (1.0/v1)*(3.0/(16.0*PI*PI))*(a0mstop1*dmstop1sqdv1+a0mstop2*dmstop2sqdv1);

  if (speak)
    {
      cout << "# t1Ov1 = " << delta1tp << endl;
    }

  return delta1tp;
}

double doCalcTadpoleMSSMH2(SoftParsMssm r, double tb)
{
  bool speak = false;

  double v = r.displayHvev();
  double v1 = v/sqrt(1.0+tb*tb);
  double v2 = v1*tb;

  DoubleVector mstop(2), mstopsq(2);

  physical_MSSM(r, mstop, mstopsq, tb);

  if (mstopsq(1) < 0.0 || mstopsq(2) < 0.0)
    {
      mstopsq(1) = fabs(mstopsq(1));
      mstopsq(2) = fabs(mstopsq(2));
    }

  double q = r.displayMu();

  double a0mstop1 = a0MSSM(mstopsq(1), q);
  double dmstop1sqdv2 = doCalcMSSMdmstop1sqdv2(r, tb);
  double a0mstop2 = a0MSSM(mstopsq(2), q);
  double dmstop2sqdv2 = doCalcMSSMdmstop2sqdv2(r, tb);

  double yt = r.displayYukawaElement(YU, 3, 3);
  double mt = yt*v2/sqrt(2.0);
  double a0mtop = a0MSSM(mt*mt, q);
  double dmtopsqdv2 = yt*yt*v2;

  double delta2tp = (1.0/v2)*(3.0/(16.0*PI*PI))*(a0mstop1*dmstop1sqdv2+a0mstop2*dmstop2sqdv2
						 -2.0*a0mtop*dmtopsqdv2);

  if (speak)
    {
      cout << "# t2Ov2 = " << delta2tp << endl;
    }

  return delta2tp;
}

// The following are helper functions useful for constructing the 
// tuning measures at one-loop order.

// MSSMdmstop1sqdvi calculates the derivative m_stop_1^2 wrt v_i.
// Note our convention is m_stop_1 is the lighter stop.
// Inputs:
//     SoftParsMssm r = MSSM model
//     double tb = tan(beta)
double doCalcMSSMdmstop1sqdv1(SoftParsMssm r, double tb)
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

double doCalcMSSMdmstop1sqdv2(SoftParsMssm r, double tb)
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
double doCalcMSSMdmstop2sqdv1(SoftParsMssm r, double tb)
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

double doCalcMSSMdmstop2sqdv2(SoftParsMssm r, double tb)
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
//     SoftParsMssm r = MSSM model
//     double tb = tan(beta)
double doCalcMSSMDeltaPrime11(SoftParsMssm r, double tb)
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
  physical_MSSM(r, mstop, mstopsq, tb);

  if (mstopsq(1) < 0.0 || mstopsq(2) < 0.0)
    {
      mstopsq(1) = fabs(mstopsq(1)); // in case tachyonic
      mstopsq(2) = fabs(mstopsq(2));
    }

  double term1 = sqr(0.125*gbar*gbar*v1)+(1.0/rt)*sqr(0.125*v1*RQQ-2.0*mt*mt*Xt*mu/v2);
  term1 = term1*log(mstopsq(1)*mstopsq(2)/(q*q*q*q));

  double term2 = (gbar*gbar*v1/(32.0*sqrt(rt)))*(v1*RQQ-16.0*mt*mt*Xt*mu/v2);
  term2 = term2*log(mstopsq(2)/mstopsq(1));

  double term3 = 0.125*gbar*gbar*(a0MSSM(mstopsq(1), q)+a0MSSM(mstopsq(2), q));

  double term4 = (1/32.0)*((1.0/sqrt(rt))*(4.0*RQQ+sqr(g2*g2-g1*g1)*v1*v1+32.0*yt*yt*mu*mu)
			   -(1.0/(rt*sqrt(rt)))*sqr(v1*RQQ-16.0*mt*mt*Xt*mu/v2));
  term4 = term4*(a0MSSM(mstopsq(2), q)-a0MSSM(mstopsq(1), q));

  double delta11p = (3.0/(16.0*PI*PI))*(term1+term2+term3+term4);

  return delta11p;
}

double doCalcMSSMDeltaPrime12(SoftParsMssm r, double tb)
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
  physical_MSSM(r, mstop, mstopsq, tb);

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
  term3 = term3*((a0MSSM(mstopsq(2), q)-a0MSSM(mstopsq(1), q)))/sqrt(rt);

  double delta12p = (3.0/(16.0*PI*PI))*(term1+term2-term3);

  return delta12p;
}

double doCalcMSSMDeltaPrime22(SoftParsMssm r, double tb)
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
  physical_MSSM(r, mstop, mstopsq, tb);

  if (mstopsq(1) < 0.0 || mstopsq(2) < 0.0)
    {
      mstopsq(1) = fabs(mstopsq(1)); // in case tachyonic
      mstopsq(2) = fabs(mstopsq(2));
    }

  double term1 = sqr(yt*yt-0.125*gbar*gbar)+sqr(8.0*yt*yt*At*Xt-RQQ)/(64.0*rt);
  term1 = term1*v2*v2*log(mstopsq(1)*mstopsq(2)/(q*q*q*q));

  double term2 =(v2*v2/(4.0*sqrt(rt)))*(yt*yt-0.125*gbar*gbar)*(8.0*yt*yt*At*Xt-RQQ);
  term2 = term2*log(mstopsq(2)/mstopsq(1));

  double term3 = (yt*yt-0.125*gbar*gbar)*(a0MSSM(mstopsq(1), q)+a0MSSM(mstopsq(2), q));

  double term4 = (1.0/sqrt(rt))*(sqr(g2*g2-g1*g1)*v2*v2/32.0-0.125*RQQ+yt*yt*At*At-(v2*v2/(32.0*rt))*sqr(8.0*Xt*At*yt*yt-RQQ));
  term4 = term4*(a0MSSM(mstopsq(2), q)-a0MSSM(mstopsq(1), q));

  double delta22p = (3.0/(16.0*PI*PI))*(term1+term2+term3+term4-2.0*yt*yt*yt*yt*v2*v2*log(mt*mt/(q*q))-2.0*yt*yt*a0MSSM(mt*mt,q));

  return delta22p;
}

// These helper functions calculate the second derivatives of 
// the one-loop corrections wrt the input parameters. They 
// return (1/v_i)\partial^2\Delta V/\partial p_j\partial v_i.
// Inputs:
//     SoftParsMssm r = MSSM model
//     double tb = tan(beta)
double doCalcMSSMd2DeltaVdMudv1(SoftParsMssm r, double tb)
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
  physical_MSSM(r, mstop, mstopsq, tb);

  if (mstopsq(1) < 0.0 || mstopsq(2) < 0.0)
    {
      mstopsq(1) = fabs(mstopsq(1)); // in case tachyonic
      mstopsq(2) = fabs(mstopsq(2));
    }

  // Calculate required first derivatives of stops
  double dmstop1sqdv1 =  doCalcMSSMdmstop1sqdv1(r, tb);
  double dmstop2sqdv1 =  doCalcMSSMdmstop2sqdv1(r, tb);
  double dmstop1sqdmu = 2.0*mt*mt*Xt/(sqrt(rt)*tb);
  double dmstop2sqdmu = -2.0*mt*mt*Xt/(sqrt(rt)*tb);

  // Calculate required second derivatives of stops
  double d2mstop1sqdmudv1 = -(v1/(2.0*sqrt(rt)))*((4.0*mt*mt*Xt/(rt*tb))*(0.25*RQQ-4.0*mt*mt*Xt*mu/(v1*v2))
						  -4.0*mt*mt*Xt/(v1*v2)+4.0*mt*mt*mu/(v2*v2));
  double d2mstop2sqdmudv1 = (v1/(2.0*sqrt(rt)))*((4.0*mt*mt*Xt/(rt*tb))*(0.25*RQQ-4.0*mt*mt*Xt*mu/(v1*v2))
						 -4.0*mt*mt*Xt/(v1*v2)+4.0*mt*mt*mu/(v2*v2));

  // Calculate total
  double d2DeltaVdmudv1 = (3.0/(16.0*PI*PI))*(1.0/v1)*(dmstop1sqdmu*dmstop1sqdv1*log(mstopsq(1)/(q*q))
					      +a0MSSM(mstopsq(1), q)*d2mstop1sqdmudv1
					      +dmstop2sqdmu*dmstop2sqdv1*log(mstopsq(2)/(q*q))
					      +a0MSSM(mstopsq(2), q)*d2mstop2sqdmudv1);
  return d2DeltaVdmudv1;
}

double doCalcMSSMd2DeltaVdMudv2(SoftParsMssm r, double tb)
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
  physical_MSSM(r, mstop, mstopsq, tb);

  if (mstopsq(1) < 0.0 || mstopsq(2) < 0.0)
    {
      mstopsq(1) = fabs(mstopsq(1)); // in case tachyonic
      mstopsq(2) = fabs(mstopsq(2));
    }

  // Calculate required first derivatives of stops
  double dmstop1sqdv2 =  doCalcMSSMdmstop1sqdv2(r, tb);
  double dmstop2sqdv2 =  doCalcMSSMdmstop2sqdv2(r, tb);
  double dmstop1sqdmu = 2.0*mt*mt*Xt/(sqrt(rt)*tb);
  double dmstop2sqdmu = -2.0*mt*mt*Xt/(sqrt(rt)*tb);

  // Calculate required second derivatives of stops
  double d2mstop1sqdmudv2 = -(v2/(2.0*sqrt(rt)))*((4.0*mt*mt*Xt/(rt*tb))*(2.0*yt*yt*At*Xt-0.25*RQQ)-2.0*yt*yt*At/tb);
  double d2mstop2sqdmudv2 = (v2/(2.0*sqrt(rt)))*((4.0*mt*mt*Xt/(rt*tb))*(2.0*yt*yt*At*Xt-0.25*RQQ)-2.0*yt*yt*At/tb);
			

  // Calculate total
  double d2DeltaVdmudv2 = (3.0/(16.0*PI*PI))*(1.0/v2)*(dmstop1sqdmu*dmstop1sqdv2*log(mstopsq(1)/(q*q))
					      +a0MSSM(mstopsq(1), q)*d2mstop1sqdmudv2
					      +dmstop2sqdmu*dmstop2sqdv2*log(mstopsq(2)/(q*q))
					      +a0MSSM(mstopsq(2), q)*d2mstop2sqdmudv2);
  return d2DeltaVdmudv2;
}

double doCalcMSSMd2DeltaVdAtdv1(SoftParsMssm r, double tb)
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
  physical_MSSM(r, mstop, mstopsq, tb);

  if (mstopsq(1) < 0.0 || mstopsq(2) < 0.0)
    {
      mstopsq(1) = fabs(mstopsq(1)); // in case tachyonic
      mstopsq(2) = fabs(mstopsq(2));
    }

  // Calculate required first derivatives of stops
  double dmstop1sqdv1 =  doCalcMSSMdmstop1sqdv1(r, tb);
  double dmstop2sqdv1 =  doCalcMSSMdmstop2sqdv1(r, tb);
  double dmstop1sqdAt = -2.0*mt*mt*Xt/sqrt(rt);
  double dmstop2sqdAt = 2.0*mt*mt*Xt/sqrt(rt);

  // Calculate required second derivatives of stops
  double d2mstop1sqdAtdv1 = (v1/(2.0*sqrt(rt)))*((4.0*mt*mt*Xt/rt)*(0.25*RQQ-4.0*mt*mt*Xt*mu/(v1*v2))+4.0*mt*mt*mu/(v1*v2));
  double d2mstop2sqdAtdv1 = -(v1/(2.0*sqrt(rt)))*((4.0*mt*mt*Xt/rt)*(0.25*RQQ-4.0*mt*mt*Xt*mu/(v1*v2))+4.0*mt*mt*mu/(v1*v2));
			

  // Calculate total
  double d2DeltaVdAtdv1 = (3.0/(16.0*PI*PI))*(1.0/v1)*(dmstop1sqdAt*dmstop1sqdv1*log(mstopsq(1)/(q*q))
					      +a0MSSM(mstopsq(1), q)*d2mstop1sqdAtdv1
					      +dmstop2sqdAt*dmstop2sqdv1*log(mstopsq(2)/(q*q))
					      +a0MSSM(mstopsq(2), q)*d2mstop2sqdAtdv1);
  return d2DeltaVdAtdv1;
}

double doCalcMSSMd2DeltaVdAtdv2(SoftParsMssm r, double tb)
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
  physical_MSSM(r, mstop, mstopsq, tb);

  if (mstopsq(1) < 0.0 || mstopsq(2) < 0.0)
    {
      mstopsq(1) = fabs(mstopsq(1)); // in case tachyonic
      mstopsq(2) = fabs(mstopsq(2));
    }

  // Calculate required first derivatives of stops
  double dmstop1sqdv2 =  doCalcMSSMdmstop1sqdv2(r, tb);
  double dmstop2sqdv2 =  doCalcMSSMdmstop2sqdv2(r, tb);
  double dmstop1sqdAt = -2.0*mt*mt*Xt/sqrt(rt);
  double dmstop2sqdAt = 2.0*mt*mt*Xt/sqrt(rt);

  // Calculate required second derivatives of stops
  double d2mstop1sqdAtdv2 = (v2/(2.0*sqrt(rt)))*((4.0*mt*mt*Xt/rt)*(2.0*yt*yt*At*Xt-0.25*RQQ)-2.0*yt*yt*(Xt+At));
  double d2mstop2sqdAtdv2 = -(v2/(2.0*sqrt(rt)))*((4.0*mt*mt*Xt/rt)*(2.0*yt*yt*At*Xt-0.25*RQQ)-2.0*yt*yt*(Xt+At));
			

  // Calculate total
  double d2DeltaVdAtdv2 = (3.0/(16.0*PI*PI))*(1.0/v2)*(dmstop1sqdAt*dmstop1sqdv2*log(mstopsq(1)/(q*q))
					      +a0MSSM(mstopsq(1), q)*d2mstop1sqdAtdv2
					      +dmstop2sqdAt*dmstop2sqdv2*log(mstopsq(2)/(q*q))
					      +a0MSSM(mstopsq(2), q)*d2mstop2sqdAtdv2);
  return d2DeltaVdAtdv2;
}

double doCalcMSSMd2DeltaVdmQlSqdv1(SoftParsMssm r, double tb)
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
  physical_MSSM(r, mstop, mstopsq, tb);

  if (mstopsq(1) < 0.0 || mstopsq(2) < 0.0)
    {
      mstopsq(1) = fabs(mstopsq(1)); // in case tachyonic
      mstopsq(2) = fabs(mstopsq(2));
    }

  // Calculate required first derivatives of stops
  double dmstop1sqdv1 =  doCalcMSSMdmstop1sqdv1(r, tb);
  double dmstop2sqdv1 =  doCalcMSSMdmstop2sqdv1(r, tb);
  double dmstop1sqdmQlSq = 0.5*(1.0-MQQSq/sqrt(rt));
  double dmstop2sqdmQlSq = 0.5*(1.0+MQQSq/sqrt(rt));

  // Calculate required second derivatives of stops
  double d2mstop1sqdmQlSqdv1 = (v1/(2.0*sqrt(rt)))*((MQQSq/rt)*(0.25*RQQ-4.0*mt*mt*Xt*mu/(v1*v2))-0.25*(g2*g2-g1*g1));
  double d2mstop2sqdmQlSqdv1 = -(v1/(2.0*sqrt(rt)))*((MQQSq/rt)*(0.25*RQQ-4.0*mt*mt*Xt*mu/(v1*v2))-0.25*(g2*g2-g1*g1));
			

  // Calculate total
  double d2DeltaVdmQlSqdv1 = (3.0/(16.0*PI*PI))*(1.0/v1)*(dmstop1sqdmQlSq*dmstop1sqdv1*log(mstopsq(1)/(q*q))
					      +a0MSSM(mstopsq(1), q)*d2mstop1sqdmQlSqdv1
					      +dmstop2sqdmQlSq*dmstop2sqdv1*log(mstopsq(2)/(q*q))
					      +a0MSSM(mstopsq(2), q)*d2mstop2sqdmQlSqdv1);
  return d2DeltaVdmQlSqdv1;
}

double doCalcMSSMd2DeltaVdmQlSqdv2(SoftParsMssm r, double tb)
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
  physical_MSSM(r, mstop, mstopsq, tb);

  if (mstopsq(1) < 0.0 || mstopsq(2) < 0.0)
    {
      mstopsq(1) = fabs(mstopsq(1)); // in case tachyonic
      mstopsq(2) = fabs(mstopsq(2));
    } 
  // Calculate required first derivatives of stops
  double dmstop1sqdv2 =  doCalcMSSMdmstop1sqdv2(r, tb);
  double dmstop2sqdv2 =  doCalcMSSMdmstop2sqdv2(r, tb);
  double dmstop1sqdmQlSq = 0.5*(1.0-MQQSq/sqrt(rt));
  double dmstop2sqdmQlSq = 0.5*(1.0+MQQSq/sqrt(rt));

  // Calculate required second derivatives of stops
  double d2mstop1sqdmQlSqdv2 = (v2/(2.0*sqrt(rt)))*((MQQSq/rt)*(2.0*yt*yt*At*Xt-0.25*RQQ)+0.25*(g2*g2-g1*g1));
  double d2mstop2sqdmQlSqdv2 = -(v2/(2.0*sqrt(rt)))*((MQQSq/rt)*(2.0*yt*yt*At*Xt-0.25*RQQ)+0.25*(g2*g2-g1*g1));
			

  // Calculate total
  double d2DeltaVdmQlSqdv2 = (3.0/(16.0*PI*PI))*(1.0/v2)*(dmstop1sqdmQlSq*dmstop1sqdv2*log(mstopsq(1)/(q*q))
					      +a0MSSM(mstopsq(1), q)*d2mstop1sqdmQlSqdv2
					      +dmstop2sqdmQlSq*dmstop2sqdv2*log(mstopsq(2)/(q*q))
					      +a0MSSM(mstopsq(2), q)*d2mstop2sqdmQlSqdv2);
  return d2DeltaVdmQlSqdv2;
}

double doCalcMSSMd2DeltaVdmUrSqdv1(SoftParsMssm r, double tb)
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
  physical_MSSM(r, mstop, mstopsq, tb);

  if (mstopsq(1) < 0.0 || mstopsq(2) < 0.0)
    {
      mstopsq(1) = fabs(mstopsq(1)); // in case tachyonic
      mstopsq(2) = fabs(mstopsq(2));
    } 
  // Calculate required first derivatives of stops
  double dmstop1sqdv1 =  doCalcMSSMdmstop1sqdv1(r, tb);
  double dmstop2sqdv1 =  doCalcMSSMdmstop2sqdv1(r, tb);
  double dmstop1sqdmUrSq = 0.5*(1.0+MQQSq/sqrt(rt));
  double dmstop2sqdmUrSq = 0.5*(1.0-MQQSq/sqrt(rt));

  // Calculate required second derivatives of stops
  double d2mstop1sqdmUrSqdv1 = -(v1/(2.0*sqrt(rt)))*((MQQSq/rt)*(0.25*RQQ-4.0*mt*mt*Xt*mu/(v1*v2))-0.25*(g2*g2-g1*g1));
  double d2mstop2sqdmUrSqdv1 = (v1/(2.0*sqrt(rt)))*((MQQSq/rt)*(0.25*RQQ-4.0*mt*mt*Xt*mu/(v1*v2))-0.25*(g2*g2-g1*g1));
			

  // Calculate total
  double d2DeltaVdmUrSqdv1 = (3.0/(16.0*PI*PI))*(1.0/v1)*(dmstop1sqdmUrSq*dmstop1sqdv1*log(mstopsq(1)/(q*q))
					      +a0MSSM(mstopsq(1), q)*d2mstop1sqdmUrSqdv1
					      +dmstop2sqdmUrSq*dmstop2sqdv1*log(mstopsq(2)/(q*q))
					      +a0MSSM(mstopsq(2), q)*d2mstop2sqdmUrSqdv1);
  return d2DeltaVdmUrSqdv1;
}

double doCalcMSSMd2DeltaVdmUrSqdv2(SoftParsMssm r, double tb)
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
  physical_MSSM(r, mstop, mstopsq, tb);

  if (mstopsq(1) < 0.0 || mstopsq(2) < 0.0)
    {
      mstopsq(1) = fabs(mstopsq(1)); // in case tachyonic
      mstopsq(2) = fabs(mstopsq(2));
    } 
  // Calculate required first derivatives of stops
  double dmstop1sqdv2 =  doCalcMSSMdmstop1sqdv2(r, tb);
  double dmstop2sqdv2 =  doCalcMSSMdmstop2sqdv2(r, tb);
  double dmstop1sqdmUrSq = 0.5*(1.0+MQQSq/sqrt(rt));
  double dmstop2sqdmUrSq = 0.5*(1.0-MQQSq/sqrt(rt));

  // Calculate required second derivatives of stops
  double d2mstop1sqdmUrSqdv2 = -(v2/(2.0*sqrt(rt)))*((MQQSq/rt)*(2.0*yt*yt*At*Xt-0.25*RQQ)+0.25*(g2*g2-g1*g1));
  double d2mstop2sqdmUrSqdv2 = (v2/(2.0*sqrt(rt)))*((MQQSq/rt)*(2.0*yt*yt*At*Xt-0.25*RQQ)+0.25*(g2*g2-g1*g1));
			

  // Calculate total
  double d2DeltaVdmUrSqdv2 = (3.0/(16.0*PI*PI))*(1.0/v2)*(dmstop1sqdmUrSq*dmstop1sqdv2*log(mstopsq(1)/(q*q))
					      +a0MSSM(mstopsq(1), q)*d2mstop1sqdmUrSqdv2
					      +dmstop2sqdmUrSq*dmstop2sqdv2*log(mstopsq(2)/(q*q))
					      +a0MSSM(mstopsq(2), q)*d2mstop2sqdmUrSqdv2);
  return d2DeltaVdmUrSqdv2;
}


// Variables used for getting information to the functions used in numerically calculating the
// fine tuning.
static SoftParsMssm *tempsoftTuning; // < a SoftParsMssm object given at the input scale MX
static double ftMsusy; // < the value of M_{SUSY} for the above object
static int ftParChoice; // < index labelling the current parameter we are varying
static DoubleVector ftPars(NUMPMSSMPARS); // < vector containing the parameters varied as part of the fine tuning
static void (*currentftBC)(SoftParsMssm & , DoubleVector &);


// A function for getting the predicted running value of M_Z^2, including 1-loop top and
// stop tadpole contributions. Used for numerically calculating the fine tuning in the
// MSSM. Passed to SOFTSUSY's calcDerivative routine as part of the fine tuning calculation.
double predpMSSMMzSqRun(double parVal)
{

  // Save copies of all initial values
  SoftParsMssm savedObject(tempsoftTuning->displaySoftPars());
  DoubleVector savedPars = ftPars;
  double savedMu = tempsoftTuning->displayMu();

  // Update parameter values
  ftPars(ftParChoice) = parVal;
  currentftBC(*tempsoftTuning, ftPars);

  // Recalculate VEVs at M_{SUSY} using Newton's method. Note that
  // this also in general changes M_{SUSY} as well. If M_{SUSY}
  // is allowed to vary then we need to iterate until we converge
  // to a new estimate for the right scale.
  bool ALLOWVARYINGMSUSY = false; // < this needs to be false to get exact agreement with our analytics

  int maxIters = 100;
  double tol = 1.0e-5;

  double ms = ftMsusy; // < guess for M_{SUSY}.
  DoubleVector vevs(2);

  int sing = 0;

  if (ALLOWVARYINGMSUSY)
    {
      // Here we need to iterate to get the new M_{SUSY}
      // when we recalculate the VEVs.
      DoubleVector mstop(2), mstopsq(2);

      for (int i = 1; i <= maxIters; i++)
	{
	  tempsoftTuning->runto(ms);

	  vevs = MSSM_EWSBNewtonSolver(*tempsoftTuning, tol, maxIters, sing);
	  if (sing != 0)
	    {
	      cerr << "WARNING: error encountered in solving MSSM EWSB conditions." << endl;
	    }

	  tempsoftTuning->setHvev(sqrt(vevs(1)*vevs(1)+vevs(2)*vevs(2)));
	  tempsoftTuning->setTanb(vevs(2)/vevs(1));

	  // Recalculate M_{SUSY} using the updated VEVs
	  physical_MSSM(*tempsoftTuning, mstop, mstopsq, vevs(2)/vevs(1));

	  if (fabs(ms-sqrt(mstop(1)*mstop(2))) < tol)
	    {
	      // Since we have converged, running to ms again is probably unnecessary.
	      ms = sqrt(mstop(1)*mstop(2));
	      break;
	    }
	  ms = sqrt(mstop(1)*mstop(2));

	  if (i == maxIters)
	    {
	      cerr << "WARNING: maximum iterations reached before estimate for M_{SUSY} converged." << endl;
	      sing = 1;
	    }

	}
    }
  else
    {

      tempsoftTuning->runto(ms);

      // If we don't vary M_{SUSY}, we just need to run down
      // and recalculate the VEVs.
      vevs = MSSM_EWSBNewtonSolver(*tempsoftTuning, tol, maxIters, sing);
      if (sing != 0)
	{
	  cerr << "WARNING: error encountered in solving MSSM EWSB conditions." << endl;
	}
    }

  tempsoftTuning->setHvev(sqrt(vevs(1)*vevs(1)+vevs(2)*vevs(2)));
  tempsoftTuning->setTanb(vevs(2)/vevs(1));

  // Calculate the new M_Z^2
  double t1Ov1 = doCalcTadpoleMSSMH1(*tempsoftTuning, vevs(2)/vevs(1));
  double t2Ov2 = doCalcTadpoleMSSMH2(*tempsoftTuning, vevs(2)/vevs(1));

  if (!INCLUDE1LPTADPOLES)
    {
      t1Ov1 = 0.0;
      t2Ov2 = 0.0;
    }

  double tb = vevs(2)/vevs(1);

  double mHdSq = tempsoftTuning->displayMh1Squared();
  double mHuSq = tempsoftTuning->displayMh2Squared();
  double mu = tempsoftTuning->displaySusyMu();
      
  double g1 = tempsoftTuning->displayGaugeCoupling(1);
  double g2 = tempsoftTuning->displayGaugeCoupling(2);
  double gbar = sqrt(g2*g2+0.6*g1*g1);
  
  double MzSq = -2.0*mu*mu+(2.0*(mHdSq+t1Ov1-tb*tb*(mHuSq+t2Ov2)))/(tb*tb-1.0);

  // Reset objects
  tempsoftTuning->setSoftPars(savedObject);
  tempsoftTuning->setMu(savedMu);
  ftPars = savedPars;

  return MzSq;

}

// Because the EWSB conditions in the MSSM are so simple, it is straightforward
// to calculate the fine tuning numerically. This will be the check on our results.
// The object r is assumed to be provided at MX. The fine tuning is calculated
// as p / M_Z^2 \frac{\partial M_Z^2}{\partial p} where p is a set of parameters. M_Z^2
// is calculated using the tree level EWSB conditions + one-loop stop and top tadpole contributions
// at the scale M_SUSY. The function BCatMX is used to set the parameters at the high input scale MX.
 DoubleVector doCalcMSSMTuningNumerically(MssmSoftsusy r, double ms, double mx, 
 					 DoubleVector pars,
 					 void (*BCatMX)(SoftParsMssm & , DoubleVector &))
 {

   int nPars = pars.displayEnd()-pars.displayStart()+1;

   DoubleVector tunings(nPars);

   // Check that we are at the right scale
   if (fabs(r.displayMu()-mx) > TOLERANCE)
    {
      cerr << "WARNING: object provided to doCalcMSSMTuningNumerically should be at scale MX = " << mx
	   << ", not Q = " << r.displayMu() << ": skipping calculation." << endl;
      for (int i = 1; i <= nPars; i++)
	{
	  tunings(i) = -numberOfTheBeast; // negative values indicate error, since tunings are positive by definition
	}
    }
  else
    {
      // Now for each parameter, estimate the numerical derivative by varying its
      // value at the input scale MX, running down to M_SUSY and recalculating M_Z.

      ftMsusy = ms;
      currentftBC = BCatMX;
      ftPars = pars;
      tempsoftTuning = &r;

      double epsilon = 1.0e-5;
      double temp;
      double h;
      double deriv;
      double err;

      // Initial value of M_Z^2
      SoftParsMssm s = r;
      s.runto(ms);

      double tb = s.displayTanb();

      double t1Ov1 = doCalcTadpoleMSSMH1(s, tb);
      double t2Ov2 = doCalcTadpoleMSSMH2(s, tb);

      if (!INCLUDE1LPTADPOLES)
	{
	  t1Ov1 = 0.0;
	  t2Ov2 = 0.0;
	}
      
      double mHdSq = s.displayMh1Squared();
      double mHuSq = s.displayMh2Squared();
      double mu = s.displaySusyMu();
      
      double g1 = s.displayGaugeCoupling(1);
      double g2 = s.displayGaugeCoupling(2);
      double gbar = sqrt(g2*g2+0.6*g1*g1);
      
      double refMzSq = -2.0*mu*mu+(2.0*(mHdSq+t1Ov1-tb*tb*(mHuSq+t2Ov2)))/(tb*tb-1.0);

      // Loop over the provided parameters and calculate the fine tuning for each
      for (int i = 1; i <= nPars; i++)
	{
	  ftParChoice = i;

	  // Initial estimate for step size h.
	  h = epsilon*fabs(pars(i));
	  if (h == 0.0) h = epsilon;
	  temp = pars(i);
	  pars(i) = temp + h;
	  h = pars(i) - temp;

	  deriv = calcDerivative(predpMSSMMzSqRun, temp, h, &err);

	  if (pars(i) > TOLERANCE && fabs(err / deriv) > 1.0) 
	    {
	      tunings(i) = -numberOfTheBeast; // tuning is inaccurate, so flag using negative value
	    }
	  else
	    {
	      tunings(i) = fabs((temp/refMzSq)*deriv);
	    }

	}

    }
  return tunings;
}

// Get the coefficient of the O(t) term in the approximate
// solution for m_H_u^2. The Softpars object rb should contain
// the values of the parameters at the starting scale Q_0, and
// nLps should be the number of loops at which the beta function
// is calculated (at the moment can only do nLps = 1). Checked against
// what I think is the one-loop beta function in SOFTSUSY, and they
// agree.
double doCalcMh2SquaredBeta1(SoftParsMssm rb)
{
  // Get parameter values
  double g1 = rb.displayGaugeCoupling(1);
  double g2 = rb.displayGaugeCoupling(2);
  double g3 = rb.displayGaugeCoupling(3);
  double yt = rb.displayYukawaElement(YU, 3, 3);
  double yb = rb.displayYukawaElement(YD, 3, 3);
  double ytau = rb.displayYukawaElement(YE, 3, 3);
  double mHuSq = rb.displayMh2Squared();
  double mHdSq = rb.displayMh1Squared();
  double M1 = rb.displayGaugino(1);
  double M2 = rb.displayGaugino(2);
  double MQ3Sq = rb.displaySoftMassSquared(mQl, 3, 3);
  double Mu3Sq = rb.displaySoftMassSquared(mUr, 3, 3);
  double Md3Sq = rb.displaySoftMassSquared(mDr, 3, 3);
  double ML3Sq = rb.displaySoftMassSquared(mLl, 3, 3);
  double Me3Sq = rb.displaySoftMassSquared(mEr, 3, 3);
  double MQ1Sq = rb.displaySoftMassSquared(mQl, 1, 1);
  double Mu1Sq = rb.displaySoftMassSquared(mUr, 1, 1);
  double Md1Sq = rb.displaySoftMassSquared(mDr, 1, 1);
  double ML1Sq = rb.displaySoftMassSquared(mLl, 1, 1);
  double Me1Sq = rb.displaySoftMassSquared(mEr, 1, 1);
  double MQ2Sq = rb.displaySoftMassSquared(mQl, 2, 2);
  double Mu2Sq = rb.displaySoftMassSquared(mUr, 2, 2);
  double Md2Sq = rb.displaySoftMassSquared(mDr, 2, 2);
  double ML2Sq = rb.displaySoftMassSquared(mLl, 2, 2);
  double Me2Sq = rb.displaySoftMassSquared(mEr, 2, 2);
  double S = mHuSq-mHdSq+(MQ3Sq-ML3Sq-2.0*Mu3Sq+Md3Sq+Me3Sq)+(MQ2Sq-ML2Sq-2.0*Mu2Sq+Md2Sq+Me2Sq)
    +(MQ1Sq-ML1Sq-2.0*Mu1Sq+Md1Sq+Me1Sq);
  double At = rb.displaySoftA(UA, 3, 3);
  double Ab = rb.displaySoftA(DA, 3, 3);
  double Atau = rb.displaySoftA(EA, 3, 3);

  double coeff;

  coeff = 6.0*yt*yt*(mHuSq+MQ3Sq+Mu3Sq+At*At)-6.0*g2*g2*M2*M2+0.6*g1*g1*(S-2.0*M1*M1);
  coeff = coeff/(16.0*PI*PI);

  return coeff;

}

// 2-loop contribution to the beta function for m_Hu^2. Checked against what I think
// is the 2-loop beta function in SOFTSUSY, and they agree.
double doCalcMh2SquaredBeta2(SoftParsMssm rb)
{
  // Get parameter values
  double g1 = rb.displayGaugeCoupling(1);
  double g2 = rb.displayGaugeCoupling(2);
  double g3 = rb.displayGaugeCoupling(3);
  double yt = rb.displayYukawaElement(YU, 3, 3);
  double yb = rb.displayYukawaElement(YD, 3, 3);
  double ytau = rb.displayYukawaElement(YE, 3, 3);
  double mHuSq = rb.displayMh2Squared();
  double mHdSq = rb.displayMh1Squared();
  double M1 = rb.displayGaugino(1);
  double M2 = rb.displayGaugino(2);
  double M3 = rb.displayGaugino(3);
  double MQ3Sq = rb.displaySoftMassSquared(mQl, 3, 3);
  double Mu3Sq = rb.displaySoftMassSquared(mUr, 3, 3);
  double Md3Sq = rb.displaySoftMassSquared(mDr, 3, 3);
  double ML3Sq = rb.displaySoftMassSquared(mLl, 3, 3);
  double Me3Sq = rb.displaySoftMassSquared(mEr, 3, 3);
  double MQ1Sq = rb.displaySoftMassSquared(mQl, 1, 1);
  double Mu1Sq = rb.displaySoftMassSquared(mUr, 1, 1);
  double Md1Sq = rb.displaySoftMassSquared(mDr, 1, 1);
  double ML1Sq = rb.displaySoftMassSquared(mLl, 1, 1);
  double Me1Sq = rb.displaySoftMassSquared(mEr, 1, 1);
  double MQ2Sq = rb.displaySoftMassSquared(mQl, 2, 2);
  double Mu2Sq = rb.displaySoftMassSquared(mUr, 2, 2);
  double Md2Sq = rb.displaySoftMassSquared(mDr, 2, 2);
  double ML2Sq = rb.displaySoftMassSquared(mLl, 2, 2);
  double Me2Sq = rb.displaySoftMassSquared(mEr, 2, 2);
  double S = mHuSq-mHdSq+(MQ3Sq-ML3Sq-2.0*Mu3Sq+Md3Sq+Me3Sq)+(MQ2Sq-ML2Sq-2.0*Mu2Sq+Md2Sq+Me2Sq)
    +(MQ1Sq-ML1Sq-2.0*Mu1Sq+Md1Sq+Me1Sq);
  double At = rb.displaySoftA(UA, 3, 3);
  double Ab = rb.displaySoftA(DA, 3, 3);
  double Atau = rb.displaySoftA(EA, 3, 3);

  double coeff;

  coeff = -36.0*sqr(yt*yt)*(mHuSq+MQ3Sq+Mu3Sq+2.0*At*At)-6.0*yt*yt*yb*yb*(mHuSq+mHdSq+2.0*MQ3Sq+Mu3Sq+Md3Sq+sqr(At+Ab))
    +32.0*g3*g3*yt*yt*(mHuSq+MQ3Sq+Mu3Sq+At*At-2.0*At*M3+2.0*M3*M3)
    +0.4*g1*g1*yt*yt*(-5.0*mHuSq+MQ3Sq+16.0*Mu3Sq+4.0*At*At-8.0*At*M1+8.0*M1*M1)
    +1.2*g1*g1*yb*yb*(3.0*mHdSq-MQ3Sq-2.0*Md3Sq)+1.2*g1*g1*ytau*ytau*(mHdSq+ML3Sq-2.0*Me3Sq)
    +3.2*g3*g3*g1*g1*(MQ3Sq-2.0*Mu3Sq+Md3Sq+MQ2Sq-2.0*Mu2Sq+Md2Sq+MQ1Sq-2.0*Mu1Sq+Md1Sq)
    +1.8*g1*g1*g2*g2*(2.0*M2*M2+2.0*M1*M1+2.0*M1*M2+mHuSq-mHdSq+MQ3Sq-ML3Sq+MQ2Sq-ML2Sq+MQ1Sq-ML1Sq)
    +3.0*sqr(g2*g2)*(11.0*M2*M2+mHuSq+mHdSq+3.0*MQ3Sq+ML3Sq+3.0*MQ2Sq+ML2Sq+3.0*MQ1Sq+ML1Sq)
    +0.04*sqr(g1*g1)*(621.0*M1*M1+18.0*mHuSq+4.0*MQ3Sq-8.0*Mu3Sq+10.0*Md3Sq+54.0*Me3Sq
		      +4.0*MQ2Sq-8.0*Mu2Sq+10.0*Md2Sq+54.0*Me2Sq+4.0*MQ1Sq-8.0*Mu1Sq+10.0*Md1Sq+54.0*Me1Sq);

  coeff = coeff/sqr(16.0*PI*PI);

  return coeff;

}

// Get the coefficient of the O(t^2) term in the approximate
// solution for m_H_u^2. The Softpars object rb should contain
// the values of the parameters at the desired scale Q_0, and
// nLogs should be the number of log coefficients to include (e.g.
// if nLogs = 1, returns the LL coefficient a^(2)/(4*pi)^4, 
// if nLogs = 2, returns the LL+NLL coefficients a^(2)/(4*pi)^4+b^(2)/(4*pi)^6, etc.)
// Currently can do nLogs = 1.
// The nLogs = 1 case has been checked numerically. It seems to be correct,
// with relative errors between this result and the numerical one using SOFTSUSY
// being between 10^-6% and 10^-4%.
double doCalcMh2SquaredLogSqCoeff(SoftParsMssm rb, int nLogs)
{
  // Get parameter values
  double g1 = rb.displayGaugeCoupling(1);
  double g2 = rb.displayGaugeCoupling(2);
  double g3 = rb.displayGaugeCoupling(3);
  double yt = rb.displayYukawaElement(YU, 3, 3);
  double yb = rb.displayYukawaElement(YD, 3, 3);
  double ytau = rb.displayYukawaElement(YE, 3, 3);
  double mHuSq = rb.displayMh2Squared();
  double mHdSq = rb.displayMh1Squared();
  double M1 = rb.displayGaugino(1);
  double M2 = rb.displayGaugino(2);
  double M3 = rb.displayGaugino(3);
  double MQ3Sq = rb.displaySoftMassSquared(mQl, 3, 3);
  double Mu3Sq = rb.displaySoftMassSquared(mUr, 3, 3);
  double Md3Sq = rb.displaySoftMassSquared(mDr, 3, 3);
  double ML3Sq = rb.displaySoftMassSquared(mLl, 3, 3);
  double Me3Sq = rb.displaySoftMassSquared(mEr, 3, 3);
  double MQ1Sq = rb.displaySoftMassSquared(mQl, 1, 1);
  double Mu1Sq = rb.displaySoftMassSquared(mUr, 1, 1);
  double Md1Sq = rb.displaySoftMassSquared(mDr, 1, 1);
  double ML1Sq = rb.displaySoftMassSquared(mLl, 1, 1);
  double Me1Sq = rb.displaySoftMassSquared(mEr, 1, 1);
  double MQ2Sq = rb.displaySoftMassSquared(mQl, 2, 2);
  double Mu2Sq = rb.displaySoftMassSquared(mUr, 2, 2);
  double Md2Sq = rb.displaySoftMassSquared(mDr, 2, 2);
  double ML2Sq = rb.displaySoftMassSquared(mLl, 2, 2);
  double Me2Sq = rb.displaySoftMassSquared(mEr, 2, 2);
  double S = mHuSq-mHdSq+(MQ3Sq-ML3Sq-2.0*Mu3Sq+Md3Sq+Me3Sq)+(MQ2Sq-ML2Sq-2.0*Mu2Sq+Md2Sq+Me2Sq)
    +(MQ1Sq-ML1Sq-2.0*Mu1Sq+Md1Sq+Me1Sq);
  double At = rb.displaySoftA(UA, 3, 3);
  double Ab = rb.displaySoftA(DA, 3, 3);
  double Atau = rb.displaySoftA(EA, 3, 3);

  double coeff;

  if (nLogs != 1)
    {
      cerr << "WARNING: can only do LL calculation at the moment: using nLogs = 1." << endl;
    }
    
  coeff = 72.0*yt*yt*yt*yt*(mHuSq+MQ3Sq+Mu3Sq+2.0*At*At)+6.0*yt*yt*yb*yb*(mHuSq+mHdSq+2.0*MQ3Sq+Mu3Sq+Md3Sq+sqr(At+Ab))
    -32.0*g3*g3*yt*yt*(mHuSq+MQ3Sq+Mu3Sq+At*At-2.0*At*M3+2.0*M3*M3)
    -18.0*g2*g2*yt*yt*(mHuSq+MQ3Sq+Mu3Sq+At*At-2.0*At*M2+2.0*M2*M2)
    -(26.0/5.0)*g1*g1*yt*yt*(mHuSq+MQ3Sq+Mu3Sq+At*At-2.0*At*M1+2.0*M1*M1)+(198.0/25.0)*g1*g1*g1*g1*(S-3*M1*M1)
    -18.0*g2*g2*g2*g2*M2*M2;
  coeff = coeff/sqr(16.0*PI*PI);
    

  return coeff;

}

// Likewise for m_Hd^2. Checked against what I think is the one-loop beta function
// in SOFTSUSY, and they agree.
double doCalcMh1SquaredBeta1(SoftParsMssm rb)
{
  // Get parameter values
  double g1 = rb.displayGaugeCoupling(1);
  double g2 = rb.displayGaugeCoupling(2);
  double g3 = rb.displayGaugeCoupling(3);
  double yt = rb.displayYukawaElement(YU, 3, 3);
  double yb = rb.displayYukawaElement(YD, 3, 3);
  double ytau = rb.displayYukawaElement(YE, 3, 3);
  double mHuSq = rb.displayMh2Squared();
  double mHdSq = rb.displayMh1Squared();
  double M1 = rb.displayGaugino(1);
  double M2 = rb.displayGaugino(2);
  double MQ3Sq = rb.displaySoftMassSquared(mQl, 3, 3);
  double Mu3Sq = rb.displaySoftMassSquared(mUr, 3, 3);
  double Md3Sq = rb.displaySoftMassSquared(mDr, 3, 3);
  double ML3Sq = rb.displaySoftMassSquared(mLl, 3, 3);
  double Me3Sq = rb.displaySoftMassSquared(mEr, 3, 3);
  double MQ1Sq = rb.displaySoftMassSquared(mQl, 1, 1);
  double Mu1Sq = rb.displaySoftMassSquared(mUr, 1, 1);
  double Md1Sq = rb.displaySoftMassSquared(mDr, 1, 1);
  double ML1Sq = rb.displaySoftMassSquared(mLl, 1, 1);
  double Me1Sq = rb.displaySoftMassSquared(mEr, 1, 1);
  double MQ2Sq = rb.displaySoftMassSquared(mQl, 2, 2);
  double Mu2Sq = rb.displaySoftMassSquared(mUr, 2, 2);
  double Md2Sq = rb.displaySoftMassSquared(mDr, 2, 2);
  double ML2Sq = rb.displaySoftMassSquared(mLl, 2, 2);
  double Me2Sq = rb.displaySoftMassSquared(mEr, 2, 2);
  double S = mHuSq-mHdSq+(MQ3Sq-ML3Sq-2.0*Mu3Sq+Md3Sq+Me3Sq)+(MQ2Sq-ML2Sq-2.0*Mu2Sq+Md2Sq+Me2Sq)
    +(MQ1Sq-ML1Sq-2.0*Mu1Sq+Md1Sq+Me1Sq);
  double At = rb.displaySoftA(UA, 3, 3);
  double Ab = rb.displaySoftA(DA, 3, 3);
  double Atau = rb.displaySoftA(EA, 3, 3);

  double coeff;

  coeff = 6.0*yb*yb*(mHdSq+MQ3Sq+Md3Sq+Ab*Ab)+2.0*ytau*ytau*(mHdSq+ML3Sq+Me3Sq+Atau*Atau)-6.0*g2*g2*M2*M2-0.6*g1*g1*(S+2.0*M1*M1);
  coeff = coeff/(16.0*PI*PI);

  return coeff;

}

// 2-loop contribution to beta. Checked against what I think is the two-loop beta
// in SOFTSUSY, and they agree.
double doCalcMh1SquaredBeta2(SoftParsMssm rb)
{
  // Get parameter values
  double g1 = rb.displayGaugeCoupling(1);
  double g2 = rb.displayGaugeCoupling(2);
  double g3 = rb.displayGaugeCoupling(3);
  double yt = rb.displayYukawaElement(YU, 3, 3);
  double yb = rb.displayYukawaElement(YD, 3, 3);
  double ytau = rb.displayYukawaElement(YE, 3, 3);
  double mHuSq = rb.displayMh2Squared();
  double mHdSq = rb.displayMh1Squared();
  double M1 = rb.displayGaugino(1);
  double M2 = rb.displayGaugino(2);
  double M3 = rb.displayGaugino(3);
  double MQ3Sq = rb.displaySoftMassSquared(mQl, 3, 3);
  double Mu3Sq = rb.displaySoftMassSquared(mUr, 3, 3);
  double Md3Sq = rb.displaySoftMassSquared(mDr, 3, 3);
  double ML3Sq = rb.displaySoftMassSquared(mLl, 3, 3);
  double Me3Sq = rb.displaySoftMassSquared(mEr, 3, 3);
  double MQ1Sq = rb.displaySoftMassSquared(mQl, 1, 1);
  double Mu1Sq = rb.displaySoftMassSquared(mUr, 1, 1);
  double Md1Sq = rb.displaySoftMassSquared(mDr, 1, 1);
  double ML1Sq = rb.displaySoftMassSquared(mLl, 1, 1);
  double Me1Sq = rb.displaySoftMassSquared(mEr, 1, 1);
  double MQ2Sq = rb.displaySoftMassSquared(mQl, 2, 2);
  double Mu2Sq = rb.displaySoftMassSquared(mUr, 2, 2);
  double Md2Sq = rb.displaySoftMassSquared(mDr, 2, 2);
  double ML2Sq = rb.displaySoftMassSquared(mLl, 2, 2);
  double Me2Sq = rb.displaySoftMassSquared(mEr, 2, 2);
  double S = mHuSq-mHdSq+(MQ3Sq-ML3Sq-2.0*Mu3Sq+Md3Sq+Me3Sq)+(MQ2Sq-ML2Sq-2.0*Mu2Sq+Md2Sq+Me2Sq)
    +(MQ1Sq-ML1Sq-2.0*Mu1Sq+Md1Sq+Me1Sq);
  double At = rb.displaySoftA(UA, 3, 3);
  double Ab = rb.displaySoftA(DA, 3, 3);
  double Atau = rb.displaySoftA(EA, 3, 3);

  double coeff;

  coeff = -36.0*sqr(yb*yb)*(mHdSq+MQ3Sq+Md3Sq+2.0*Ab*Ab)-6.0*yt*yt*yb*yb*(mHuSq+mHdSq+2.0*MQ3Sq+Mu3Sq+Md3Sq+sqr(At+Ab))
    -12.0*sqr(ytau*ytau)*(mHdSq+ML3Sq+Me3Sq+2.0*Atau*Atau)+32.0*g3*g3*yb*yb*(mHdSq+MQ3Sq+Md3Sq+Ab*Ab-2.0*M3*Ab+2.0*M3*M3)
    +1.2*g1*g1*yt*yt*(3.0*mHuSq+MQ3Sq-4.0*Mu3Sq)-0.4*g1*g1*yb*yb*(11.0*mHdSq-MQ3Sq-4.0*Md3Sq+2.0*Ab*Ab-4.0*M1*Ab+4.0*M1*M1)
    +1.2*g1*g1*ytau*ytau*(mHdSq+ML3Sq+4.0*Me3Sq+2.0*Atau*Atau-4.0*M1*Atau+4.0*M1*M1)
    -3.2*g3*g3*g1*g1*(MQ3Sq-2.0*Mu3Sq+Md3Sq+MQ2Sq-2.0*Mu2Sq+Md2Sq+MQ1Sq-2.0*Mu1Sq+Md1Sq)
    +3.0*sqr(g2*g2)*(11.0*M2*M2+mHuSq+mHdSq+3.0*MQ3Sq+ML3Sq+3.0*MQ2Sq+ML2Sq+3.0*MQ1Sq+ML1Sq)
    +1.8*g1*g1*g2*g2*(2.0*M2*M2+2.0*M1*M1+2.0*M1*M2-mHuSq+mHdSq+ML3Sq-MQ3Sq+ML2Sq-MQ2Sq+ML1Sq-MQ1Sq)
    +0.04*sqr(g1*g1)*(621.0*M1*M1+18.0*mHdSq+18.0*ML3Sq+2.0*MQ3Sq+56.0*Mu3Sq+2.0*Md3Sq-18.0*Me3Sq
		      +18.0*ML2Sq+2.0*MQ2Sq+56.0*Mu2Sq+2.0*Md2Sq-18.0*Me2Sq
		      +18.0*ML1Sq+2.0*MQ1Sq+56.0*Mu1Sq+2.0*Md1Sq-18.0*Me1Sq);

  coeff = coeff/sqr(16.0*PI*PI);

  return coeff;

}

// The nLogs = 1 case has been checked numerically. It seems to be correct,
// with relative errors between this result and the numerical one using SOFTSUSY
// being between 10^-6% and 10^-4%.
double doCalcMh1SquaredLogSqCoeff(SoftParsMssm rb, int nLogs)
{
  // Get parameter values
  double g1 = rb.displayGaugeCoupling(1);
  double g2 = rb.displayGaugeCoupling(2);
  double g3 = rb.displayGaugeCoupling(3);
  double yt = rb.displayYukawaElement(YU, 3, 3);
  double yb = rb.displayYukawaElement(YD, 3, 3);
  double ytau = rb.displayYukawaElement(YE, 3, 3);
  double mHuSq = rb.displayMh2Squared();
  double mHdSq = rb.displayMh1Squared();
  double M1 = rb.displayGaugino(1);
  double M2 = rb.displayGaugino(2);
  double M3 = rb.displayGaugino(3);
  double MQ3Sq = rb.displaySoftMassSquared(mQl, 3, 3);
  double Mu3Sq = rb.displaySoftMassSquared(mUr, 3, 3);
  double Md3Sq = rb.displaySoftMassSquared(mDr, 3, 3);
  double ML3Sq = rb.displaySoftMassSquared(mLl, 3, 3);
  double Me3Sq = rb.displaySoftMassSquared(mEr, 3, 3);
  double MQ1Sq = rb.displaySoftMassSquared(mQl, 1, 1);
  double Mu1Sq = rb.displaySoftMassSquared(mUr, 1, 1);
  double Md1Sq = rb.displaySoftMassSquared(mDr, 1, 1);
  double ML1Sq = rb.displaySoftMassSquared(mLl, 1, 1);
  double Me1Sq = rb.displaySoftMassSquared(mEr, 1, 1);
  double MQ2Sq = rb.displaySoftMassSquared(mQl, 2, 2);
  double Mu2Sq = rb.displaySoftMassSquared(mUr, 2, 2);
  double Md2Sq = rb.displaySoftMassSquared(mDr, 2, 2);
  double ML2Sq = rb.displaySoftMassSquared(mLl, 2, 2);
  double Me2Sq = rb.displaySoftMassSquared(mEr, 2, 2);
  double S = mHuSq-mHdSq+(MQ3Sq-ML3Sq-2.0*Mu3Sq+Md3Sq+Me3Sq)+(MQ2Sq-ML2Sq-2.0*Mu2Sq+Md2Sq+Me2Sq)
    +(MQ1Sq-ML1Sq-2.0*Mu1Sq+Md1Sq+Me1Sq);
  double At = rb.displaySoftA(UA, 3, 3);
  double Ab = rb.displaySoftA(DA, 3, 3);
  double Atau = rb.displaySoftA(EA, 3, 3);

  double coeff;

  if (nLogs != 1)
    {
      cerr << "WARNING: cannot yet do more than LL calculation: using nLogs = 1." << endl;
    }
  coeff = 72.0*sqr(yb*yb)*(mHdSq+MQ3Sq+Md3Sq+2.0*Ab*Ab)+16.0*sqr(ytau*ytau)*(mHdSq+ML3Sq+Me3Sq+2.0*Atau*Atau)
    +6.0*yt*yt*yb*yb*(mHuSq+mHdSq+2.0*MQ3Sq+Mu3Sq+Md3Sq+sqr(At+Ab))
    +12.0*ytau*ytau*yb*yb*(2.0*mHdSq+MQ3Sq+Md3Sq+ML3Sq+Me3Sq+sqr(Atau+Ab))
    -32.0*g3*g3*yb*yb*(mHdSq+MQ3Sq+Md3Sq+Ab*Ab-2.0*M3*Ab+2.0*M3*M3)
    -18.0*g2*g2*yb*yb*(mHdSq+MQ3Sq+Md3Sq+Ab*Ab-2.0*M2*Ab+2.0*M2*M2)
    -6.0*g2*g2*ytau*ytau*(mHdSq+ML3Sq+Me3Sq+Atau*Atau-2.0*M2*Atau+2.0*M2*M2)
    -(14.0/5.0)*g1*g1*yb*yb*(mHdSq+MQ3Sq+Md3Sq+Ab*Ab-2.0*M1*Ab+2.0*M1*M1)
    -(18.0/5.0)*g1*g1*ytau*ytau*(mHdSq+ML3Sq+Me3Sq+Atau*Atau-2.0*M1*Atau+2.0*M1*M1)
    -18.0*sqr(g2*g2)*M2*M2-(198.0/25.0)*sqr(g1*g1)*(S+3.0*M1*M1);
  coeff = coeff/sqr(16.0*PI*PI);

  return coeff;

}

// The coefficient of the term linear log term in the expansion
// of m_Hu^2. The object rb should contain the parameters evaluated
// at the desired initial scale Q_0, while nLps is the number of contributions
// to include (e.g. if nLps = 2, the 1- and 2-loop contributions are included)
double doCalcMh2SquaredLogCoeff(SoftParsMssm rb, int nLps)
{

  if (nLps == 1)
    {
      return doCalcMh2SquaredBeta1(rb);      
    }
  else if (nLps == 2)
    {
      return (doCalcMh2SquaredBeta1(rb)+doCalcMh2SquaredBeta2(rb));
    }
  else
    {
      cerr << "WARNING: can only do 1 or 2 loops: returning nLps = 1 loop result." << endl;
      return doCalcMh2SquaredBeta1(rb);      
    }
}

// The coefficient of the term linear log term in the expansion
// of m_Hd^2. The object rb should contain the parameters evaluated
// at the desired initial scale Q_0, while nLps is the number of contributions
// to include (e.g. if nLps = 2, the 1- and 2-loop contributions are included)
double doCalcMh1SquaredLogCoeff(SoftParsMssm rb, int nLps)
{

  if (nLps == 1)
    {
      return doCalcMh1SquaredBeta1(rb);      
    }
  else if (nLps == 2)
    {
      return (doCalcMh1SquaredBeta1(rb)+doCalcMh1SquaredBeta2(rb));
    }
  else
    {
      cerr << "WARNING: can only do 1 or 2 loops: returning nLps = 1 loop result." << endl;
      return doCalcMh1SquaredBeta1(rb);      
    }

}

double doCalcMuBeta1(SoftParsMssm rb)
{
  // Get parameter values
  double mu = rb.displaySusyMu();
  double g1 = rb.displayGaugeCoupling(1);
  double g2 = rb.displayGaugeCoupling(2);
  double yt = rb.displayYukawaElement(YU, 3, 3);
  double yb = rb.displayYukawaElement(YD, 3, 3);
  double ytau = rb.displayYukawaElement(YE, 3, 3);

  double coeff;

  coeff = mu*(3.0*yt*yt+3.0*yb*yb+ytau*ytau-3.0*g2*g2-0.6*g1*g1);
  coeff = coeff/(16.0*PI*PI);

  return coeff;

}

double doCalcMuBeta2(SoftParsMssm rb)
{
  // Get parameter values
  double mu = rb.displaySusyMu();
  double g1 = rb.displayGaugeCoupling(1);
  double g2 = rb.displayGaugeCoupling(2);
  double g3 = rb.displayGaugeCoupling(3);
  double yt = rb.displayYukawaElement(YU, 3, 3);
  double yb = rb.displayYukawaElement(YD, 3, 3);
  double ytau = rb.displayYukawaElement(YE, 3, 3);

  double coeff;

  coeff = mu*(-9.0*sqr(yt*yt)-9.0*sqr(yb*yb)-6.0*sqr(yt*yb)-3.0*sqr(ytau*ytau)
	      +16.0*sqr(g3*yt)+0.8*sqr(g1*yt)+16.0*sqr(g3*yb)-0.4*sqr(g1*yb)+1.2*sqr(g1*ytau)
	      +7.5*sqr(g2*g2)+1.8*sqr(g1*g2)+(207.0/50.0)*sqr(g1*g1));

  coeff = coeff/sqr(16.0*PI*PI);

  return coeff;
}

double doCalcMuLogCoeff(SoftParsMssm rb, int nLps)
{
  if (nLps == 1)
    {
      return doCalcMuBeta1(rb);      
    }
  else if (nLps == 2)
    {
      return (doCalcMuBeta1(rb)+doCalcMuBeta2(rb));
    }
  else
    {
      cerr << "WARNING: can only do 1 or 2 loops: returning nLps = 1 loop result." << endl;
      return doCalcMuBeta1(rb);      
    }
}

double doCalcMuLogSqCoeff(SoftParsMssm rb, int nLogs)
{
  // Get parameter values
  double mu = rb.displaySusyMu();
  double g1 = rb.displayGaugeCoupling(1);
  double g2 = rb.displayGaugeCoupling(2);
  double g3 = rb.displayGaugeCoupling(3);
  double yt = rb.displayYukawaElement(YU, 3, 3);
  double yb = rb.displayYukawaElement(YD, 3, 3);
  double ytau = rb.displayYukawaElement(YE, 3, 3);

  double coeff;

  coeff = 0.5*mu*(45.0*(sqr(yt*yt)+sqr(yb*yb))+9.0*sqr(ytau*ytau)+30.0*sqr(yt*yb)+6.0*sqr(yt*ytau)
		  +18.0*sqr(yb*ytau)-32.0*g3*g3*(yt*yt+yb*yb)-12.0*g2*g2*(3.0*yt*yt+3.0*yb*yb+ytau*ytau)
		  -0.8*g1*g1*(11.0*yt*yt+8.0*yb*yb+6.0*ytau*ytau)+3.0*sqr(g2*g2)-(189.0/25.0)*sqr(g1*g1)
		  +(18.0/5.0)*sqr(g1*g2));

  coeff = coeff/sqr(16.0*PI*PI);

  return coeff;
}

double doCalcBBeta1(SoftParsMssm rb)
{
 // Get parameter values
  double g1 = rb.displayGaugeCoupling(1);
  double g2 = rb.displayGaugeCoupling(2);
  double yt = rb.displayYukawaElement(YU, 3, 3);
  double yb = rb.displayYukawaElement(YD, 3, 3);
  double ytau = rb.displayYukawaElement(YE, 3, 3);
  double M1 = rb.displayGaugino(1);
  double M2 = rb.displayGaugino(2);
  double At = rb.displaySoftA(UA, 3, 3);
  double Ab = rb.displaySoftA(DA, 3, 3);
  double Atau = rb.displaySoftA(EA, 3, 3);

  double coeff;

  coeff = 6.0*yt*yt*At+6.0*yb*yb*Ab+2.0*ytau*ytau*Atau+6.0*g2*g2*M2+1.2*g1*g1*M1;
  coeff = coeff/(16.0*PI*PI);

  return coeff;

}

double doCalcBBeta2(SoftParsMssm rb)
{
  // Get parameter values
  double g1 = rb.displayGaugeCoupling(1);
  double g2 = rb.displayGaugeCoupling(2);
  double g3 = rb.displayGaugeCoupling(3);
  double yt = rb.displayYukawaElement(YU, 3, 3);
  double yb = rb.displayYukawaElement(YD, 3, 3);
  double ytau = rb.displayYukawaElement(YE, 3, 3);
  double M1 = rb.displayGaugino(1);
  double M2 = rb.displayGaugino(2);
  double M3 = rb.displayGaugino(3);
  double At = rb.displaySoftA(UA, 3, 3);
  double Ab = rb.displaySoftA(DA, 3, 3);
  double Atau = rb.displaySoftA(EA, 3, 3);

  double coeff;

  coeff = -36.0*sqr(yt*yt)*At-36.0*sqr(yb*yb)*Ab-12.0*sqr(ytau*ytau)*Atau-12.0*sqr(yt*yb)*(At+Ab)
    +32.0*sqr(g3*yt)*(At-M3)+32.0*sqr(g3*yb)*(Ab-M3)+1.6*sqr(g1*yt)*(At-M1)-0.8*sqr(g1*yb)*(Ab-M1)
    +2.4*sqr(g1*ytau)*(Atau-M1)-30.0*sqr(g2*g2)*M2-3.6*sqr(g1*g2)*(M1+M2)-(414.0/25.0)*sqr(g1*g1)*M1;

  coeff = coeff/sqr(16.0*PI*PI);

  return coeff;
}

double doCalcBLogCoeff(SoftParsMssm rb, int nLps)
{
  if (nLps == 1)
    {
      return doCalcBBeta1(rb);      
    }
  else if (nLps == 2)
    {
      return (doCalcBBeta1(rb)+doCalcBBeta2(rb));
    }
  else
    {
      cerr << "WARNING: can only do 1 or 2 loops: returning nLps = 1 loop result." << endl;
      return doCalcBBeta1(rb);      
    }
}

double doCalcBLogSqCoeff(SoftParsMssm rb, int nLogs)
{

  // Get parameter values
  double g1 = rb.displayGaugeCoupling(1);
  double g2 = rb.displayGaugeCoupling(2);
  double g3 = rb.displayGaugeCoupling(3);
  double yt = rb.displayYukawaElement(YU, 3, 3);
  double yb = rb.displayYukawaElement(YD, 3, 3);
  double ytau = rb.displayYukawaElement(YE, 3, 3);
  double M1 = rb.displayGaugino(1);
  double M2 = rb.displayGaugino(2);
  double M3 = rb.displayGaugino(3);
  double At = rb.displaySoftA(UA, 3, 3);
  double Ab = rb.displaySoftA(DA, 3, 3);
  double Atau = rb.displaySoftA(EA, 3, 3);

  double coeff;

  if (nLogs != 1)
    {
      cerr << "WARNING: cannot yet do more than LL calculation: using nLogs = 1." << endl;
    }
  coeff = 72.0*sqr(yt*yt)*At+72.0*sqr(yb*yb)*Ab+16.0*sqr(ytau*ytau)*Atau+12.0*sqr(yt*yb)*(At+Ab)
    +12.0*sqr(ytau*yb)*(Atau+Ab)-32.0*sqr(g3*yt)*(At-M3)-32.0*sqr(g3*yb)*(Ab-M3)
    -18.0*sqr(g2*yt)*(At-M2)-18.0*sqr(g2*yb)*(Ab-M2)-6.0*sqr(g2*ytau)*(Atau-M2)
    -5.2*sqr(g1*yt)*(At-M1)-2.8*sqr(g1*yb)*(Ab-M1)-3.6*sqr(g1*ytau)*(Atau-M1)+12.0*sqr(g2*g2)*M2
    +(396.0/25.0)*sqr(g1*g1)*M1;

  coeff = coeff/sqr(16.0*PI*PI);

  return coeff;
}


double doCalcAtLogCoeff(SoftParsMssm rb, int nLps)
{

  if (nLps == 1)
    {
      return doCalcAtBeta1(rb);      
    }
  else if (nLps == 2)
    {
      return (doCalcAtBeta1(rb)+doCalcAtBeta2(rb));
    }
  else
    {
      cerr << "WARNING: can only do 1 or 2 loops: returning nLps = 1 loop result." << endl;
      return doCalcAtBeta1(rb);      
    }

}

double doCalcAtBeta1(SoftParsMssm rb)
{
  // Get parameter values
  double g1 = rb.displayGaugeCoupling(1);
  double g2 = rb.displayGaugeCoupling(2);
  double g3 = rb.displayGaugeCoupling(3);
  double yt = rb.displayYukawaElement(YU, 3, 3);
  double yb = rb.displayYukawaElement(YD, 3, 3);
  double M1 = rb.displayGaugino(1);
  double M2 = rb.displayGaugino(2);
  double M3 = rb.displayGaugino(3);
  double At = rb.displaySoftA(UA, 3, 3);
  double Ab = rb.displaySoftA(DA, 3, 3);

  double coeff;

  coeff = 12.0*yt*yt*At+2.0*yb*yb*Ab+(32.0/3.0)*g3*g3*M3+6.0*g2*g2*M2+(26.0/15.0)*g1*g1*M1;
  coeff = coeff/(16.0*PI*PI);

  return coeff;

}

double doCalcAtBeta2(SoftParsMssm rb)
{
  // Get parameter values
  double g1 = rb.displayGaugeCoupling(1);
  double g2 = rb.displayGaugeCoupling(2);
  double g3 = rb.displayGaugeCoupling(3);
  double yt = rb.displayYukawaElement(YU, 3, 3);
  double yb = rb.displayYukawaElement(YD, 3, 3);
  double ytau = rb.displayYukawaElement(YE, 3, 3);
  double M1 = rb.displayGaugino(1);
  double M2 = rb.displayGaugino(2);
  double M3 = rb.displayGaugino(3);
  double At = rb.displaySoftA(UA, 3, 3);
  double Ab = rb.displaySoftA(DA, 3, 3);
  double Atau = rb.displaySoftA(EA, 3, 3);

  double coeff;

  coeff = -88.0*sqr(yt*yt)*At-20.0*sqr(yb*yb)*Ab-10.0*sqr(yt*yb)*(At+Ab)-2.0*sqr(yb*ytau)*(Ab+Atau)
    +32.0*sqr(g3*yt)*(At-M3)+12.0*sqr(g2*yt)*(At-M2)+(12.0/5.0)*sqr(g1*yt)*(At-M1)
    +(4.0/5.0)*sqr(g1*yb)*(Ab-M1)+(64.0/9.0)*sqr(g3*g3)*M3-16.0*sqr(g3*g2)*(M2+M3)-(272.0/45.0)*sqr(g3*g1)*(M1+M3)
    -30.0*sqr(g2*g2)*M2-2.0*sqr(g2*g1)*(M1+M2)-(5486.0/225.0)*sqr(g1*g1)*M1;

  coeff = coeff/sqr(16.0*PI*PI);

  return coeff;
}

double doCalcAtLogSqCoeff(SoftParsMssm rb, int nLogs)
{
  // Get parameter values
  double g1 = rb.displayGaugeCoupling(1);
  double g2 = rb.displayGaugeCoupling(2);
  double g3 = rb.displayGaugeCoupling(3);
  double yt = rb.displayYukawaElement(YU, 3, 3);
  double yb = rb.displayYukawaElement(YD, 3, 3);
  double ytau = rb.displayYukawaElement(YE, 3, 3);
  double M1 = rb.displayGaugino(1);
  double M2 = rb.displayGaugino(2);
  double M3 = rb.displayGaugino(3);
  double At = rb.displaySoftA(UA, 3, 3);
  double Ab = rb.displaySoftA(DA, 3, 3);
  double Atau = rb.displaySoftA(EA, 3, 3);

  double coeff;

  if (nLogs != 1)
    {
      cerr << "WARNING: cannot yet do more than LL calculation: using nLogs = 1." << endl;
    }
  coeff = 144.0*sqr(yt*yt)*At+24.0*sqr(yb*yb)*Ab+14.0*sqr(yt*yb)*(At+Ab)+2.0*sqr(yb*ytau)*(Atau+Ab)
    -64.0*sqr(g3*yt)*(At-M3)-36.0*sqr(g2*yt)*(At-M2)-10.4*sqr(g1*yt)*(At-M1)
    -(32.0/3.0)*sqr(g3*yb)*(Ab-M3)-6.0*sqr(g2*yb)*(Ab-M2)-(14.0/15.0)*sqr(g1*yb)*(Ab-M1)
    -64.0*sqr(g3*g3)*M3+12.0*sqr(g2*g2)*M2+(572.0/25.0)*sqr(g1*g1)*M1;

  coeff = coeff/sqr(16.0*PI*PI);

  return coeff;

}

double doCalcmqL3SquaredBeta1(SoftParsMssm rb)
{
  double g1 = rb.displayGaugeCoupling(1);
  double g2 = rb.displayGaugeCoupling(2);
  double g3 = rb.displayGaugeCoupling(3);
  double yt = rb.displayYukawaElement(YU, 3, 3);
  double yb = rb.displayYukawaElement(YD, 3, 3);
  double mHuSq = rb.displayMh2Squared();
  double mHdSq = rb.displayMh1Squared();
  double M1 = rb.displayGaugino(1);
  double M2 = rb.displayGaugino(2);
  double M3 = rb.displayGaugino(3);
  double MQ3Sq = rb.displaySoftMassSquared(mQl, 3, 3);
  double Mu3Sq = rb.displaySoftMassSquared(mUr, 3, 3);
  double Md3Sq = rb.displaySoftMassSquared(mDr, 3, 3);
  double ML3Sq = rb.displaySoftMassSquared(mLl, 3, 3);
  double Me3Sq = rb.displaySoftMassSquared(mEr, 3, 3);
  double MQ1Sq = rb.displaySoftMassSquared(mQl, 1, 1);
  double Mu1Sq = rb.displaySoftMassSquared(mUr, 1, 1);
  double Md1Sq = rb.displaySoftMassSquared(mDr, 1, 1);
  double ML1Sq = rb.displaySoftMassSquared(mLl, 1, 1);
  double Me1Sq = rb.displaySoftMassSquared(mEr, 1, 1);
  double MQ2Sq = rb.displaySoftMassSquared(mQl, 2, 2);
  double Mu2Sq = rb.displaySoftMassSquared(mUr, 2, 2);
  double Md2Sq = rb.displaySoftMassSquared(mDr, 2, 2);
  double ML2Sq = rb.displaySoftMassSquared(mLl, 2, 2);
  double Me2Sq = rb.displaySoftMassSquared(mEr, 2, 2);
  double S = mHuSq-mHdSq+(MQ3Sq-ML3Sq-2.0*Mu3Sq+Md3Sq+Me3Sq)+(MQ2Sq-ML2Sq-2.0*Mu2Sq+Md2Sq+Me2Sq)
    +(MQ1Sq-ML1Sq-2.0*Mu1Sq+Md1Sq+Me1Sq);
  double At = rb.displaySoftA(UA, 3, 3);
  double Ab = rb.displaySoftA(DA, 3, 3);

  double coeff;

  coeff = 2.0*yt*yt*(mHuSq+MQ3Sq+Mu3Sq+At*At)+2.0*yb*yb*(mHdSq+MQ3Sq+Md3Sq+Ab*Ab)
    -(32.0/3.0)*g3*g3*M3*M3-6.0*g2*g2*M2*M2+(1.0/5.0)*g1*g1*(S-(2.0/3.0)*M1*M1);
  coeff = coeff/(16.0*PI*PI);

  return coeff;
}

double doCalcmqL3SquaredBeta2(SoftParsMssm rb)
{
  // Get parameter values
  double g1 = rb.displayGaugeCoupling(1);
  double g2 = rb.displayGaugeCoupling(2);
  double g3 = rb.displayGaugeCoupling(3);
  double yt = rb.displayYukawaElement(YU, 3, 3);
  double yb = rb.displayYukawaElement(YD, 3, 3);
  double ytau = rb.displayYukawaElement(YE, 3, 3);
  double mHuSq = rb.displayMh2Squared();
  double mHdSq = rb.displayMh1Squared();
  double M1 = rb.displayGaugino(1);
  double M2 = rb.displayGaugino(2);
  double M3 = rb.displayGaugino(3);
  double MQ3Sq = rb.displaySoftMassSquared(mQl, 3, 3);
  double Mu3Sq = rb.displaySoftMassSquared(mUr, 3, 3);
  double Md3Sq = rb.displaySoftMassSquared(mDr, 3, 3);
  double ML3Sq = rb.displaySoftMassSquared(mLl, 3, 3);
  double Me3Sq = rb.displaySoftMassSquared(mEr, 3, 3);
  double MQ1Sq = rb.displaySoftMassSquared(mQl, 1, 1);
  double Mu1Sq = rb.displaySoftMassSquared(mUr, 1, 1);
  double Md1Sq = rb.displaySoftMassSquared(mDr, 1, 1);
  double ML1Sq = rb.displaySoftMassSquared(mLl, 1, 1);
  double Me1Sq = rb.displaySoftMassSquared(mEr, 1, 1);
  double MQ2Sq = rb.displaySoftMassSquared(mQl, 2, 2);
  double Mu2Sq = rb.displaySoftMassSquared(mUr, 2, 2);
  double Md2Sq = rb.displaySoftMassSquared(mDr, 2, 2);
  double ML2Sq = rb.displaySoftMassSquared(mLl, 2, 2);
  double Me2Sq = rb.displaySoftMassSquared(mEr, 2, 2);
  double At = rb.displaySoftA(UA, 3, 3);
  double Ab = rb.displaySoftA(DA, 3, 3);
  double Atau = rb.displaySoftA(EA, 3, 3);

  double coeff;

  coeff = -20.0*sqr(yt*yt)*(mHuSq+MQ3Sq+Mu3Sq+2.0*At*At)-20.0*sqr(yb*yb)*(mHdSq+MQ3Sq+Md3Sq+2.0*Ab*Ab)
    -2.0*sqr(ytau*yb)*(2.0*mHdSq+MQ3Sq+ML3Sq+Md3Sq+Me3Sq+sqr(Ab+Atau))
    +0.4*sqr(g1*yt)*(mHuSq+3.0*MQ3Sq+8.0*Mu3Sq+4.0*At*At-8.0*At*M1+8.0*M1*M1)
    +0.4*sqr(g1*yb)*(5.0*mHdSq+MQ3Sq+2.0*Ab*Ab-4.0*Ab*M1+4.0*M1*M1)
    +0.4*sqr(g1*ytau)*(mHdSq+ML3Sq-2.0*Me3Sq)
    +(16.0/3.0)*sqr(g3*g3)*(2.0*(MQ1Sq+MQ2Sq+MQ3Sq)+Mu1Sq+Mu2Sq+Mu3Sq+Md1Sq+Md2Sq+Md3Sq-8.0*M3*M3)
    +32.0*sqr(g3*g2)*(M3*M3+M2*M2+M2*M3)
    +(16.0/45.0)*sqr(g3*g1)*(3.0*(MQ1Sq+MQ2Sq+MQ3Sq-2.0*Mu1Sq-2.0*Mu2Sq-2.0*Mu3Sq+Md1Sq+Md2Sq+Md3Sq)
			     +2.0*(M3*M3+M1*M1+M1*M3))
    +3.0*sqr(g2*g2)*(mHuSq+mHdSq+3.0*(MQ1Sq+MQ2Sq+MQ3Sq)+ML1Sq+ML2Sq+ML3Sq+11.0*M2*M2)
    +0.2*sqr(g1*g2)*(2.0*M1*M1+2.0*M2*M2+2.0*M1*M2+3.0*mHuSq-3.0*mHdSq+3.0*(MQ1Sq+MQ2Sq+MQ3Sq-ML1Sq-ML2Sq-ML3Sq))
    +(1.0/75.0)*sqr(g1*g1)*(199.0*M1*M1+12.0*mHuSq-6.0*mHdSq+2.0*(MQ1Sq+MQ2Sq+MQ3Sq)-6.0*(ML1Sq+ML2Sq+ML3Sq)
			    -24.0*(Mu1Sq+Mu2Sq+Mu3Sq)+6.0*(Md1Sq+Md2Sq+Md3Sq)+42.0*(Me1Sq+Me2Sq+Me3Sq));

  coeff = coeff/sqr(16.0*PI*PI);

  return coeff;
}

double doCalcmqL3SquaredLogCoeff(SoftParsMssm rb, int nLps)
{

  if (nLps == 1)
    {
      return doCalcmqL3SquaredBeta1(rb);      
    }
  else if (nLps == 2)
    {
      return (doCalcmqL3SquaredBeta1(rb)+doCalcmqL3SquaredBeta2(rb));
    }
  else
    {
      cerr << "WARNING: can only do 1 or 2 loops: returning nLps = 1 loop result." << endl;
      return doCalcmqL3SquaredBeta1(rb);      
    }
}

double doCalcmqL3SquaredLogSqCoeff(SoftParsMssm rb, int nLogs)
{
  // Get parameter values
  double g1 = rb.displayGaugeCoupling(1);
  double g2 = rb.displayGaugeCoupling(2);
  double g3 = rb.displayGaugeCoupling(3);
  double yt = rb.displayYukawaElement(YU, 3, 3);
  double yb = rb.displayYukawaElement(YD, 3, 3);
  double ytau = rb.displayYukawaElement(YE, 3, 3);
  double mHuSq = rb.displayMh2Squared();
  double mHdSq = rb.displayMh1Squared();
  double M1 = rb.displayGaugino(1);
  double M2 = rb.displayGaugino(2);
  double M3 = rb.displayGaugino(3);
  double MQ3Sq = rb.displaySoftMassSquared(mQl, 3, 3);
  double Mu3Sq = rb.displaySoftMassSquared(mUr, 3, 3);
  double Md3Sq = rb.displaySoftMassSquared(mDr, 3, 3);
  double ML3Sq = rb.displaySoftMassSquared(mLl, 3, 3);
  double Me3Sq = rb.displaySoftMassSquared(mEr, 3, 3);
  double MQ1Sq = rb.displaySoftMassSquared(mQl, 1, 1);
  double Mu1Sq = rb.displaySoftMassSquared(mUr, 1, 1);
  double Md1Sq = rb.displaySoftMassSquared(mDr, 1, 1);
  double ML1Sq = rb.displaySoftMassSquared(mLl, 1, 1);
  double Me1Sq = rb.displaySoftMassSquared(mEr, 1, 1);
  double MQ2Sq = rb.displaySoftMassSquared(mQl, 2, 2);
  double Mu2Sq = rb.displaySoftMassSquared(mUr, 2, 2);
  double Md2Sq = rb.displaySoftMassSquared(mDr, 2, 2);
  double ML2Sq = rb.displaySoftMassSquared(mLl, 2, 2);
  double Me2Sq = rb.displaySoftMassSquared(mEr, 2, 2);
  double S = mHuSq-mHdSq+(MQ3Sq-ML3Sq-2.0*Mu3Sq+Md3Sq+Me3Sq)+(MQ2Sq-ML2Sq-2.0*Mu2Sq+Md2Sq+Me2Sq)
    +(MQ1Sq-ML1Sq-2.0*Mu1Sq+Md1Sq+Me1Sq);
  double At = rb.displaySoftA(UA, 3, 3);
  double Ab = rb.displaySoftA(DA, 3, 3);
  double Atau = rb.displaySoftA(EA, 3, 3);

  double coeff;

  if (nLogs != 1)
    {
      cerr << "WARNING: cannot yet do more than LL calculation: using nLogs = 1." << endl;
    }
  coeff = 24.0*sqr(yt*yt)*(mHuSq+MQ3Sq+Mu3Sq+2.0*At*At)+4.0*sqr(yt*yb)*(mHuSq+mHdSq+2.0*MQ3Sq+Mu3Sq+Md3Sq+sqr(At+Ab))
    -(32.0/3.0)*sqr(g3*yt)*(mHuSq+MQ3Sq+Mu3Sq+At*At-2.0*At*M3+2.0*M3*M3)
    -6.0*sqr(g2*yt)*(mHuSq+MQ3Sq+Mu3Sq+At*At-2.0*At*M2+2.0*M2*M2)
    -(26.0/15.0)*sqr(g1*yt)*(mHuSq+MQ3Sq+Mu3Sq+At*At-2.0*At*M1+2.0*M1*M1)
    +24.0*sqr(yb*yb)*(mHdSq+MQ3Sq+Md3Sq+2.0*Ab*Ab)+2.0*sqr(ytau*yb)*(2.0*mHdSq+MQ3Sq+ML3Sq+Md3Sq+Me3Sq+sqr(Atau+Ab))
    -(32.0/3.0)*sqr(g3*yb)*(mHdSq+MQ3Sq+Md3Sq+Ab*Ab-2.0*Ab*M3+2.0*M3*M3)
    -6.0*sqr(g2*yb)*(mHdSq+MQ3Sq+Md3Sq+Ab*Ab-2.0*Ab*M2+2.0*M2*M2)
    -(14.0/15.0)*sqr(g1*yb)*(mHdSq+MQ3Sq+Md3Sq+Ab*Ab-2.0*Ab*M1+2.0*M1*M1)
    +96.0*sqr(g3*g3)*M3*M3-18.0*sqr(g2*g2)*M2*M2+(66.0/25.0)*sqr(g1*g1)*(S-M1*M1);
  coeff = coeff/sqr(16.0*PI*PI);

  return coeff;
}

double doCalcmtRSquaredBeta1(SoftParsMssm rb)
{
  double g1 = rb.displayGaugeCoupling(1);
  double g3 = rb.displayGaugeCoupling(3);
  double yt = rb.displayYukawaElement(YU, 3, 3);
  double mHuSq = rb.displayMh2Squared();
  double mHdSq = rb.displayMh1Squared();
  double M1 = rb.displayGaugino(1);
  double M2 = rb.displayGaugino(2);
  double M3 = rb.displayGaugino(3);
  double MQ3Sq = rb.displaySoftMassSquared(mQl, 3, 3);
  double Mu3Sq = rb.displaySoftMassSquared(mUr, 3, 3);
  double Md3Sq = rb.displaySoftMassSquared(mDr, 3, 3);
  double ML3Sq = rb.displaySoftMassSquared(mLl, 3, 3);
  double Me3Sq = rb.displaySoftMassSquared(mEr, 3, 3);
  double MQ1Sq = rb.displaySoftMassSquared(mQl, 1, 1);
  double Mu1Sq = rb.displaySoftMassSquared(mUr, 1, 1);
  double Md1Sq = rb.displaySoftMassSquared(mDr, 1, 1);
  double ML1Sq = rb.displaySoftMassSquared(mLl, 1, 1);
  double Me1Sq = rb.displaySoftMassSquared(mEr, 1, 1);
  double MQ2Sq = rb.displaySoftMassSquared(mQl, 2, 2);
  double Mu2Sq = rb.displaySoftMassSquared(mUr, 2, 2);
  double Md2Sq = rb.displaySoftMassSquared(mDr, 2, 2);
  double ML2Sq = rb.displaySoftMassSquared(mLl, 2, 2);
  double Me2Sq = rb.displaySoftMassSquared(mEr, 2, 2);
  double S = mHuSq-mHdSq+(MQ3Sq-ML3Sq-2.0*Mu3Sq+Md3Sq+Me3Sq)+(MQ2Sq-ML2Sq-2.0*Mu2Sq+Md2Sq+Me2Sq)
    +(MQ1Sq-ML1Sq-2.0*Mu1Sq+Md1Sq+Me1Sq);
  double At = rb.displaySoftA(UA, 3, 3);


  double coeff;

  coeff = 4.0*yt*yt*(mHuSq+MQ3Sq+Mu3Sq+At*At)-(32.0/3.0)*g3*g3*M3*M3-0.8*g1*g1*(S+(8.0/3.0)*M1*M1);
  coeff = coeff/(16.0*PI*PI);

  return coeff;

}

double doCalcmtRSquaredBeta2(SoftParsMssm rb)
{
  // Get parameter values
  double g1 = rb.displayGaugeCoupling(1);
  double g2 = rb.displayGaugeCoupling(2);
  double g3 = rb.displayGaugeCoupling(3);
  double yt = rb.displayYukawaElement(YU, 3, 3);
  double yb = rb.displayYukawaElement(YD, 3, 3);
  double ytau = rb.displayYukawaElement(YE, 3, 3);
  double mHuSq = rb.displayMh2Squared();
  double mHdSq = rb.displayMh1Squared();
  double M1 = rb.displayGaugino(1);
  double M2 = rb.displayGaugino(2);
  double M3 = rb.displayGaugino(3);
  double MQ3Sq = rb.displaySoftMassSquared(mQl, 3, 3);
  double Mu3Sq = rb.displaySoftMassSquared(mUr, 3, 3);
  double Md3Sq = rb.displaySoftMassSquared(mDr, 3, 3);
  double ML3Sq = rb.displaySoftMassSquared(mLl, 3, 3);
  double Me3Sq = rb.displaySoftMassSquared(mEr, 3, 3);
  double MQ1Sq = rb.displaySoftMassSquared(mQl, 1, 1);
  double Mu1Sq = rb.displaySoftMassSquared(mUr, 1, 1);
  double Md1Sq = rb.displaySoftMassSquared(mDr, 1, 1);
  double ML1Sq = rb.displaySoftMassSquared(mLl, 1, 1);
  double Me1Sq = rb.displaySoftMassSquared(mEr, 1, 1);
  double MQ2Sq = rb.displaySoftMassSquared(mQl, 2, 2);
  double Mu2Sq = rb.displaySoftMassSquared(mUr, 2, 2);
  double Md2Sq = rb.displaySoftMassSquared(mDr, 2, 2);
  double ML2Sq = rb.displaySoftMassSquared(mLl, 2, 2);
  double Me2Sq = rb.displaySoftMassSquared(mEr, 2, 2);
  double At = rb.displaySoftA(UA, 3, 3);
  double Ab = rb.displaySoftA(DA, 3, 3);

  double coeff;

  coeff = -32.0*sqr(yt*yt)*(mHuSq+MQ3Sq+Mu3Sq+2.0*At*At)-4.0*sqr(yt*yb)*(mHuSq+mHdSq+2.0*MQ3Sq+Mu3Sq+Md3Sq+sqr(At+Ab))
    +12.0*sqr(g2*yt)*(mHuSq+MQ3Sq+Mu3Sq+At*At-2.0*M2*At+2.0*M2*M2)
    -1.6*sqr(g1*yb)*(3.0*mHdSq-MQ3Sq-2.0*Md3Sq)-1.6*sqr(g1*ytau)*(mHdSq+ML3Sq-2.0*Me3Sq)
    -0.8*sqr(g1*yt)*(-5.0*mHuSq-MQ3Sq+9.0*Mu3Sq+At*At-2.0*At*M1+2.0*M1*M1)
    +(16.0/3.0)*sqr(g3*g3)*(2.0*(MQ1Sq+MQ2Sq+MQ3Sq)+Mu1Sq+Mu2Sq+Mu3Sq+Md1Sq+Md2Sq+Md3Sq-8.0*M3*M3)
    +(64.0/45.0)*sqr(g3*g1)*(8.0*(M3*M3+M1*M1+M1*M3)-3.0*(MQ1Sq+MQ2Sq+MQ3Sq-2.0*Mu1Sq-2.0*Mu2Sq-2.0*Mu3Sq+Md1Sq+Md2Sq+Md3Sq))
    -2.4*sqr(g1*g2)*(mHuSq-mHdSq+MQ1Sq+MQ2Sq+MQ3Sq-ML1Sq-ML2Sq-ML3Sq)
    +(4.0/75.0)*sqr(g1*g1)*(856.0*M1*M1+3.0*mHuSq+21.0*mHdSq+3.0*(MQ1Sq+MQ2Sq+MQ3Sq)+21.0*(ML1Sq+ML2Sq+ML3Sq)
			    +4.0*(Md1Sq+Md2Sq+Md3Sq)-12.0*(Me1Sq+Me2Sq+Me3Sq)+64.0*(Mu1Sq+Mu2Sq+Mu3Sq));
  coeff = coeff/sqr(16.0*PI*PI);

  return coeff;
}

double doCalcmtRSquaredLogCoeff(SoftParsMssm rb, int nLps)
{
  if (nLps == 1)
    {
      return doCalcmtRSquaredBeta1(rb);      
    }
  else if (nLps == 2)
    {
      return (doCalcmtRSquaredBeta1(rb)+doCalcmtRSquaredBeta2(rb));
    }
  else
    {
      cerr << "WARNING: can only do 1 or 2 loops: returning nLps = 1 loop result." << endl;
      return doCalcmtRSquaredBeta1(rb);      
    }
}

double doCalcmtRSquaredLogSqCoeff(SoftParsMssm rb, int nLogs)
{
  // Get parameter values
  double g1 = rb.displayGaugeCoupling(1);
  double g2 = rb.displayGaugeCoupling(2);
  double g3 = rb.displayGaugeCoupling(3);
  double yt = rb.displayYukawaElement(YU, 3, 3);
  double yb = rb.displayYukawaElement(YD, 3, 3);
  double ytau = rb.displayYukawaElement(YE, 3, 3);
  double mHuSq = rb.displayMh2Squared();
  double mHdSq = rb.displayMh1Squared();
  double M1 = rb.displayGaugino(1);
  double M2 = rb.displayGaugino(2);
  double M3 = rb.displayGaugino(3);
  double MQ3Sq = rb.displaySoftMassSquared(mQl, 3, 3);
  double Mu3Sq = rb.displaySoftMassSquared(mUr, 3, 3);
  double Md3Sq = rb.displaySoftMassSquared(mDr, 3, 3);
  double ML3Sq = rb.displaySoftMassSquared(mLl, 3, 3);
  double Me3Sq = rb.displaySoftMassSquared(mEr, 3, 3);
  double MQ1Sq = rb.displaySoftMassSquared(mQl, 1, 1);
  double Mu1Sq = rb.displaySoftMassSquared(mUr, 1, 1);
  double Md1Sq = rb.displaySoftMassSquared(mDr, 1, 1);
  double ML1Sq = rb.displaySoftMassSquared(mLl, 1, 1);
  double Me1Sq = rb.displaySoftMassSquared(mEr, 1, 1);
  double MQ2Sq = rb.displaySoftMassSquared(mQl, 2, 2);
  double Mu2Sq = rb.displaySoftMassSquared(mUr, 2, 2);
  double Md2Sq = rb.displaySoftMassSquared(mDr, 2, 2);
  double ML2Sq = rb.displaySoftMassSquared(mLl, 2, 2);
  double Me2Sq = rb.displaySoftMassSquared(mEr, 2, 2);
  double S = mHuSq-mHdSq+(MQ3Sq-ML3Sq-2.0*Mu3Sq+Md3Sq+Me3Sq)+(MQ2Sq-ML2Sq-2.0*Mu2Sq+Md2Sq+Me2Sq)
    +(MQ1Sq-ML1Sq-2.0*Mu1Sq+Md1Sq+Me1Sq);
  double At = rb.displaySoftA(UA, 3, 3);
  double Ab = rb.displaySoftA(DA, 3, 3);
  double Atau = rb.displaySoftA(EA, 3, 3);

  double coeff;

  if (nLogs != 1)
    {
      cerr << "WARNING: cannot yet do more than LL calculation: using nLogs = 1." << endl;
    }
  coeff = 48.0*sqr(yt*yt)*(mHuSq+MQ3Sq+Mu3Sq+2.0*At*At)+4.0*sqr(yt*yb)*(mHuSq+mHdSq+2.0*MQ3Sq+Mu3Sq+Md3Sq+sqr(At+Ab))
    -(64.0/3.0)*sqr(g3*yt)*(mHuSq+MQ3Sq+Mu3Sq+At*At-2.0*At*M3+2.0*M3*M3)
    -12.0*sqr(g2*yt)*(mHuSq+MQ3Sq+Mu3Sq+At*At-2.0*At*M2+2.0*M2*M2)
    -(52.0/15.0)*sqr(g1*yt)*(mHuSq+MQ3Sq+Mu3Sq+At*At-2.0*At*M1+2.0*M1*M1)
    +96.0*sqr(g3*g3)*M3*M3-(264.0/25.0)*sqr(g1*g1)*(S+4.0*M1*M1);
  coeff = coeff/sqr(16.0*PI*PI);

  return coeff;
}

// Functions that are used in the routine for testing the 1-loop squared
// parts of the approximate RGE solutions
static SoftParsMssm *tempRunner;
static double q0;
double getMuBeta1(double t)
{
  // Save the initial object
  SoftParsMssm savedObject(tempRunner->displaySoftPars());
  double initialScale = tempRunner->displayMu();
  int initialLoops = tempRunner->displayLoops();

  // Get the 1-loop beta function
  tempRunner->setLoops(1);

  // Run to the new scale q
  tempRunner->runto(q0*exp(t));

  // Make sure trilinears and Yukawas are zero (the only reason
  // they do not stay zero when running is due to numerical error).
  tempRunner->setYukawaElement(YU, 1, 1, 0.0);
  tempRunner->setYukawaElement(YU, 2, 2, 0.0);
  tempRunner->setYukawaElement(YD, 1, 1, 0.0);
  tempRunner->setYukawaElement(YD, 2, 2, 0.0);
  tempRunner->setYukawaElement(YE, 1, 1, 0.0);
  tempRunner->setYukawaElement(YE, 2, 2, 0.0);

  tempRunner->setTrilinearElement(UA, 1, 1, 0.0);
  tempRunner->setTrilinearElement(UA, 2, 2, 0.0);
  tempRunner->setTrilinearElement(DA, 1, 1, 0.0);
  tempRunner->setTrilinearElement(DA, 2, 2, 0.0);
  tempRunner->setTrilinearElement(EA, 1, 1, 0.0);
  tempRunner->setTrilinearElement(EA, 2, 2, 0.0);

  SoftParsMssm beta1 = tempRunner->beta2();
  double muBeta1 = beta1.displaySusyMu();

  // Reset object
  tempRunner->setSoftPars(savedObject);
  tempRunner->setMu(initialScale);
  tempRunner->setLoops(initialLoops);

  return muBeta1;
}

double getBBeta1(double t)
{
  // Save the initial object
  SoftParsMssm savedObject(tempRunner->displaySoftPars());
  double initialScale = tempRunner->displayMu();
  int initialLoops = tempRunner->displayLoops();

  // Get the 1-loop beta function
  tempRunner->setLoops(1);

  // Run to the new scale q
  tempRunner->runto(q0*exp(t));

  // Make sure trilinears and Yukawas are zero (the only reason
  // they do not stay zero when running is due to numerical error).
  tempRunner->setYukawaElement(YU, 1, 1, 0.0);
  tempRunner->setYukawaElement(YU, 2, 2, 0.0);
  tempRunner->setYukawaElement(YD, 1, 1, 0.0);
  tempRunner->setYukawaElement(YD, 2, 2, 0.0);
  tempRunner->setYukawaElement(YE, 1, 1, 0.0);
  tempRunner->setYukawaElement(YE, 2, 2, 0.0);

  tempRunner->setTrilinearElement(UA, 1, 1, 0.0);
  tempRunner->setTrilinearElement(UA, 2, 2, 0.0);
  tempRunner->setTrilinearElement(DA, 1, 1, 0.0);
  tempRunner->setTrilinearElement(DA, 2, 2, 0.0);
  tempRunner->setTrilinearElement(EA, 1, 1, 0.0);
  tempRunner->setTrilinearElement(EA, 2, 2, 0.0);

  SoftParsMssm beta1 = tempRunner->beta2();
  double BBeta1 = (beta1.displayM3Squared()-(tempRunner->displayM3Squared()/tempRunner->displaySusyMu())*beta1.displaySusyMu())/tempRunner->displaySusyMu();

  // Reset object
  tempRunner->setSoftPars(savedObject);
  tempRunner->setMu(initialScale);
  tempRunner->setLoops(initialLoops);

  return BBeta1;
}

double getmHdSqBeta1(double t)
{
  // Save the initial object
  SoftParsMssm savedObject(tempRunner->displaySoftPars());
  double initialScale = tempRunner->displayMu();
  int initialLoops = tempRunner->displayLoops();

  // Get the 1-loop beta function
  tempRunner->setLoops(1);

  // Run to the new scale q
  tempRunner->runto(q0*exp(t));

  // Make sure trilinears and Yukawas are zero (the only reason
  // they do not stay zero when running is due to numerical error).
  tempRunner->setYukawaElement(YU, 1, 1, 0.0);
  tempRunner->setYukawaElement(YU, 2, 2, 0.0);
  tempRunner->setYukawaElement(YD, 1, 1, 0.0);
  tempRunner->setYukawaElement(YD, 2, 2, 0.0);
  tempRunner->setYukawaElement(YE, 1, 1, 0.0);
  tempRunner->setYukawaElement(YE, 2, 2, 0.0);

  tempRunner->setTrilinearElement(UA, 1, 1, 0.0);
  tempRunner->setTrilinearElement(UA, 2, 2, 0.0);
  tempRunner->setTrilinearElement(DA, 1, 1, 0.0);
  tempRunner->setTrilinearElement(DA, 2, 2, 0.0);
  tempRunner->setTrilinearElement(EA, 1, 1, 0.0);
  tempRunner->setTrilinearElement(EA, 2, 2, 0.0);

  SoftParsMssm beta1 = tempRunner->beta2();
  double mHdSqBeta1 = beta1.displayMh1Squared();

  // Reset object
  tempRunner->setSoftPars(savedObject);
  tempRunner->setMu(initialScale);
  tempRunner->setLoops(initialLoops);

  return mHdSqBeta1;
}

double getmHuSqBeta1(double t)
{
  // Save the initial object
  SoftParsMssm savedObject(tempRunner->displaySoftPars());
  double initialScale = tempRunner->displayMu();
  int initialLoops = tempRunner->displayLoops();

  // Get the 1-loop beta function
  tempRunner->setLoops(1);

  // Run to the new scale q
  tempRunner->runto(q0*exp(t));

  // Make sure trilinears and Yukawas are zero (the only reason
  // they do not stay zero when running is due to numerical error).
  tempRunner->setYukawaElement(YU, 1, 1, 0.0);
  tempRunner->setYukawaElement(YU, 2, 2, 0.0);
  tempRunner->setYukawaElement(YD, 1, 1, 0.0);
  tempRunner->setYukawaElement(YD, 2, 2, 0.0);
  tempRunner->setYukawaElement(YE, 1, 1, 0.0);
  tempRunner->setYukawaElement(YE, 2, 2, 0.0);

  tempRunner->setTrilinearElement(UA, 1, 1, 0.0);
  tempRunner->setTrilinearElement(UA, 2, 2, 0.0);
  tempRunner->setTrilinearElement(DA, 1, 1, 0.0);
  tempRunner->setTrilinearElement(DA, 2, 2, 0.0);
  tempRunner->setTrilinearElement(EA, 1, 1, 0.0);
  tempRunner->setTrilinearElement(EA, 2, 2, 0.0);


  SoftParsMssm beta1 = tempRunner->beta2();
  double mHuSqBeta1 = beta1.displayMh2Squared();

  // Reset object
  tempRunner->setSoftPars(savedObject);
  tempRunner->setMu(initialScale);
  tempRunner->setLoops(initialLoops);

  return mHuSqBeta1;
}

double getmQlSqBeta1(double t)
{
  // Save the initial object
  SoftParsMssm savedObject(tempRunner->displaySoftPars());
  double initialScale = tempRunner->displayMu();
  int initialLoops = tempRunner->displayLoops();

  // Get the 1-loop beta function
  tempRunner->setLoops(1);

  // Run to the new scale q
  tempRunner->runto(q0*exp(t));

  // Make sure trilinears and Yukawas are zero (the only reason
  // they do not stay zero when running is due to numerical error).
  tempRunner->setYukawaElement(YU, 1, 1, 0.0);
  tempRunner->setYukawaElement(YU, 2, 2, 0.0);
  tempRunner->setYukawaElement(YD, 1, 1, 0.0);
  tempRunner->setYukawaElement(YD, 2, 2, 0.0);
  tempRunner->setYukawaElement(YE, 1, 1, 0.0);
  tempRunner->setYukawaElement(YE, 2, 2, 0.0);

  tempRunner->setTrilinearElement(UA, 1, 1, 0.0);
  tempRunner->setTrilinearElement(UA, 2, 2, 0.0);
  tempRunner->setTrilinearElement(DA, 1, 1, 0.0);
  tempRunner->setTrilinearElement(DA, 2, 2, 0.0);
  tempRunner->setTrilinearElement(EA, 1, 1, 0.0);
  tempRunner->setTrilinearElement(EA, 2, 2, 0.0);

  SoftParsMssm beta1 = tempRunner->beta2();
  double mQlSqBeta1 = beta1.displaySoftMassSquared(mQl, 3, 3);

  // Reset object
  tempRunner->setSoftPars(savedObject);
  tempRunner->setMu(initialScale);
  tempRunner->setLoops(initialLoops);

  return mQlSqBeta1;
}

double getmUrSqBeta1(double t)
{
  // Save the initial object
  SoftParsMssm savedObject(tempRunner->displaySoftPars());
  double initialScale = tempRunner->displayMu();
  int initialLoops = tempRunner->displayLoops();

  // Get the 1-loop beta function
  tempRunner->setLoops(1);

  // Run to the new scale q
  tempRunner->runto(q0*exp(t));

  // Make sure trilinears and Yukawas are zero (the only reason
  // they do not stay zero when running is due to numerical error).
  tempRunner->setYukawaElement(YU, 1, 1, 0.0);
  tempRunner->setYukawaElement(YU, 2, 2, 0.0);
  tempRunner->setYukawaElement(YD, 1, 1, 0.0);
  tempRunner->setYukawaElement(YD, 2, 2, 0.0);
  tempRunner->setYukawaElement(YE, 1, 1, 0.0);
  tempRunner->setYukawaElement(YE, 2, 2, 0.0);

  tempRunner->setTrilinearElement(UA, 1, 1, 0.0);
  tempRunner->setTrilinearElement(UA, 2, 2, 0.0);
  tempRunner->setTrilinearElement(DA, 1, 1, 0.0);
  tempRunner->setTrilinearElement(DA, 2, 2, 0.0);
  tempRunner->setTrilinearElement(EA, 1, 1, 0.0);
  tempRunner->setTrilinearElement(EA, 2, 2, 0.0);

  SoftParsMssm beta1 = tempRunner->beta2();
  double mUrSqBeta1 = beta1.displaySoftMassSquared(mUr, 3, 3);

  // Reset object
  tempRunner->setSoftPars(savedObject);
  tempRunner->setMu(initialScale);
  tempRunner->setLoops(initialLoops);

  return mUrSqBeta1;
}

double getAtBeta1(double t)
{

  // Save the initial object
  SoftParsMssm savedObject(tempRunner->displaySoftPars());
  double initialScale = tempRunner->displayMu();
  int initialLoops = tempRunner->displayLoops();

  // Get the 1-loop beta function
  tempRunner->setLoops(1);

  // Run to the new scale q
  tempRunner->runto(q0*exp(t));

  // Make sure trilinears and Yukawas are zero (the only reason
  // they do not stay zero when running is due to numerical error).
  tempRunner->setYukawaElement(YU, 1, 1, 0.0);
  tempRunner->setYukawaElement(YU, 2, 2, 0.0);
  tempRunner->setYukawaElement(YD, 1, 1, 0.0);
  tempRunner->setYukawaElement(YD, 2, 2, 0.0);
  tempRunner->setYukawaElement(YE, 1, 1, 0.0);
  tempRunner->setYukawaElement(YE, 2, 2, 0.0);

  tempRunner->setTrilinearElement(UA, 1, 1, 0.0);
  tempRunner->setTrilinearElement(UA, 2, 2, 0.0);
  tempRunner->setTrilinearElement(DA, 1, 1, 0.0);
  tempRunner->setTrilinearElement(DA, 2, 2, 0.0);
  tempRunner->setTrilinearElement(EA, 1, 1, 0.0);
  tempRunner->setTrilinearElement(EA, 2, 2, 0.0);

  SoftParsMssm beta1 = tempRunner->beta2();
  double AtBeta1 = (beta1.displayTrilinear(UA, 3, 3)-
		    tempRunner->displaySoftA(UA, 3, 3)*beta1.displayYukawaElement(YU, 3, 3))
    /tempRunner->displayYukawaElement(YU, 3, 3);

  // Reset object
  tempRunner->setSoftPars(savedObject);
  tempRunner->setMu(initialScale);
  tempRunner->setLoops(initialLoops);

  return AtBeta1;
}

// This routine tests the approximate solutions to the MSSM RGEs by comparing
// our results to those calculated in SOFTSUSY. Then also numerically differentiates
// the 1-loop beta functions found in SOFTSUSY to test our results for the 1-loop
// squared contributions.
void doTestRGEApproxSolns(SoftParsMssm r)
{
  // First make sure the first and second generation Yukawas
  // and trilinears vanish - this has been assumed in our derivation
  // of the approximate solutions.
  r.setYukawaElement(YU, 1, 1, 0.0);
  r.setYukawaElement(YU, 2, 2, 0.0);
  r.setYukawaElement(YD, 1, 1, 0.0);
  r.setYukawaElement(YD, 2, 2, 0.0);
  r.setYukawaElement(YE, 1, 1, 0.0);
  r.setYukawaElement(YE, 2, 2, 0.0);
  r.setTrilinearElement(UA, 1, 1, 0.0);
  r.setTrilinearElement(UA, 2, 2, 0.0);
  r.setTrilinearElement(DA, 1, 1, 0.0);
  r.setTrilinearElement(DA, 2, 2, 0.0);
  r.setTrilinearElement(EA, 1, 1, 0.0);
  r.setTrilinearElement(EA, 2, 2, 0.0);

  // Now get the beta functions in SOFTSUSY
  r.setLoops(1);
  SoftParsMssm rbetas1lp = r.beta2();
  r.setLoops(2);
  SoftParsMssm rbetas2lp = r.beta2();

  double betamu1lp = rbetas1lp.displaySusyMu();
  double betamu2lp = rbetas2lp.displaySusyMu();
  double betamu2lppiece = betamu2lp-betamu1lp;

  double betaB1lp = (rbetas1lp.displayM3Squared()-(r.displayM3Squared()/r.displaySusyMu())*betamu1lp)/r.displaySusyMu();
  double betaB2lp = (rbetas2lp.displayM3Squared()-(r.displayM3Squared()/r.displaySusyMu())*betamu2lp)/r.displaySusyMu();
  double betaB2lppiece = betaB2lp-betaB1lp;

  double betamHdSq1lp = rbetas1lp.displayMh1Squared();
  double betamHdSq2lp = rbetas2lp.displayMh1Squared();
  double betamHdSq2lppiece = betamHdSq2lp-betamHdSq1lp;

  double betamHuSq1lp = rbetas1lp.displayMh2Squared();
  double betamHuSq2lp = rbetas2lp.displayMh2Squared();
  double betamHuSq2lppiece = betamHuSq2lp-betamHuSq1lp;

  double betamQlSq1lp = rbetas1lp.displaySoftMassSquared(mQl, 3, 3);
  double betamQlSq2lp = rbetas2lp.displaySoftMassSquared(mQl, 3, 3);
  double betamQlSq2lppiece = betamQlSq2lp-betamQlSq1lp;

  double betamUrSq1lp = rbetas1lp.displaySoftMassSquared(mUr, 3, 3);
  double betamUrSq2lp = rbetas2lp.displaySoftMassSquared(mUr, 3, 3);
  double betamUrSq2lppiece = betamUrSq2lp-betamUrSq1lp;

  double betaAt1lp = (rbetas1lp.displayTrilinear(UA, 3, 3) - 
		      r.displaySoftA(UA, 3, 3)*rbetas1lp.displayYukawaElement(YU, 3, 3))/r.displayYukawaElement(YU, 3, 3);
  double betaAt2lp = (rbetas2lp.displayTrilinear(UA, 3, 3) - 
		      r.displaySoftA(UA, 3, 3)*rbetas2lp.displayYukawaElement(YU, 3, 3))/r.displayYukawaElement(YU, 3, 3);
  double betaAt2lppiece = betaAt2lp-betaAt1lp;

  // Now get our versions
  double approxmubeta1lp = doCalcMuBeta1(r);
  double approxmubeta2lppiece = doCalcMuBeta2(r);
  double approxmubeta2lp = doCalcMuLogCoeff(r, 2);

  double approxBbeta1lp = doCalcBBeta1(r);
  double approxBbeta2lppiece = doCalcBBeta2(r);
  double approxBbeta2lp = doCalcBLogCoeff(r, 2);

  double approxmHdSqbeta1lp = doCalcMh1SquaredBeta1(r);
  double approxmHdSqbeta2lppiece = doCalcMh1SquaredBeta2(r);
  double approxmHdSqbeta2lp = doCalcMh1SquaredLogCoeff(r, 2);

  double approxmHuSqbeta1lp = doCalcMh2SquaredBeta1(r);
  double approxmHuSqbeta2lppiece = doCalcMh2SquaredBeta2(r);
  double approxmHuSqbeta2lp = doCalcMh2SquaredLogCoeff(r, 2);

  double approxmQlSqbeta1lp = doCalcmqL3SquaredBeta1(r);
  double approxmQlSqbeta2lppiece = doCalcmqL3SquaredBeta2(r);
  double approxmQlSqbeta2lp = doCalcmqL3SquaredLogCoeff(r, 2);

  double approxmUrSqbeta1lp = doCalcmtRSquaredBeta1(r);
  double approxmUrSqbeta2lppiece = doCalcmtRSquaredBeta2(r);
  double approxmUrSqbeta2lp = doCalcmtRSquaredLogCoeff(r, 2);

  double approxAtbeta1lp = doCalcAtBeta1(r);
  double approxAtbeta2lppiece = doCalcAtBeta2(r);
  double approxAtbeta2lp = doCalcAtLogCoeff(r, 2);

  // Differences
  double mu1lpdiff = betamu1lp-approxmubeta1lp;
  double mu2lppiecediff = betamu2lppiece-approxmubeta2lppiece;
  double mu2lpdiff = betamu2lp-approxmubeta2lp;

  double B1lpdiff = betaB1lp-approxBbeta1lp;
  double B2lppiecediff = betaB2lppiece-approxBbeta2lppiece;
  double B2lpdiff = betaB2lp-approxBbeta2lp;

  double mHdSq1lpdiff = betamHdSq1lp-approxmHdSqbeta1lp;
  double mHdSq2lppiecediff = betamHdSq2lppiece-approxmHdSqbeta2lppiece;
  double mHdSq2lpdiff = betamHdSq2lp-approxmHdSqbeta2lp;

  double mHuSq1lpdiff = betamHuSq1lp-approxmHuSqbeta1lp;
  double mHuSq2lppiecediff = betamHuSq2lppiece-approxmHuSqbeta2lppiece;
  double mHuSq2lpdiff = betamHuSq2lp-approxmHuSqbeta2lp;

  double mQlSq1lpdiff = betamQlSq1lp-approxmQlSqbeta1lp;
  double mQlSq2lppiecediff = betamQlSq2lppiece-approxmQlSqbeta2lppiece;
  double mQlSq2lpdiff = betamQlSq2lp-approxmQlSqbeta2lp;

  double mUrSq1lpdiff = betamUrSq1lp-approxmUrSqbeta1lp;
  double mUrSq2lppiecediff = betamUrSq2lppiece-approxmUrSqbeta2lppiece;
  double mUrSq2lpdiff = betamUrSq2lp-approxmUrSqbeta2lp;

  double At1lpdiff = betaAt1lp-approxAtbeta1lp;
  double At2lppiecediff = betaAt2lppiece-approxAtbeta2lppiece;
  double At2lpdiff = betaAt2lp-approxAtbeta2lp;

  // Print out results for comparison
  cout << "# MU: " << endl;
  cout << "# SOFTSUSY (approx) beta_mu^(1)   = "; printRow(cout, betamu1lp); cout << " ("; printRow(cout, approxmubeta1lp); 
  cout << " ), diff. = "; printRow(cout, mu1lpdiff); cout << "." << endl;
  cout << "# SOFTSUSY (approx) beta_mu^(2)   = "; printRow(cout, betamu2lppiece); cout << " ("; 
  printRow(cout, approxmubeta2lppiece);
  cout << " ), diff. = "; printRow(cout, mu2lppiecediff); cout << "." << endl;
  cout << "# SOFTSUSY (approx) beta_mu^(1+2) = "; printRow(cout, betamu2lp); cout << " ("; printRow(cout, approxmubeta2lp); 
  cout << " ), diff. = "; printRow(cout, mu2lpdiff); cout << "." << endl;
  cout << "# B: " << endl;
  cout << "# SOFTSUSY (approx) beta_B^(1)   = "; printRow(cout, betaB1lp); cout << " ("; printRow(cout, approxBbeta1lp); 
  cout << " ), diff. = "; printRow(cout, B1lpdiff); cout << "." << endl;
  cout << "# SOFTSUSY (approx) beta_B^(2)   = "; printRow(cout, betaB2lppiece); cout << " ("; 
  printRow(cout, approxBbeta2lppiece);
  cout << " ), diff. = "; printRow(cout, B2lppiecediff); cout << "." << endl;
  cout << "# SOFTSUSY (approx) beta_B^(1+2) = "; printRow(cout, betaB2lp); cout << " ("; printRow(cout, approxBbeta2lp); 
  cout << " ), diff. = "; printRow(cout, B2lpdiff); cout << "." << endl;
  cout << "# MH1SQ: " << endl;
  cout << "# SOFTSUSY (approx) beta_mHdSq^(1)   = "; printRow(cout, betamHdSq1lp); cout << " ("; 
  printRow(cout, approxmHdSqbeta1lp); 
  cout << " ), diff. = "; printRow(cout, mHdSq1lpdiff); cout << "." << endl;
  cout << "# SOFTSUSY (approx) beta_mHdSq^(2)   = "; printRow(cout, betamHdSq2lppiece); cout << " ("; 
  printRow(cout, approxmHdSqbeta2lppiece);
  cout << " ), diff. = "; printRow(cout, mHdSq2lppiecediff); cout << "." << endl;
  cout << "# SOFTSUSY (approx) beta_mHdSq^(1+2) = "; printRow(cout, betamHdSq2lp); cout << " ("; 
  printRow(cout, approxmHdSqbeta2lp); 
  cout << " ), diff. = "; printRow(cout, mHdSq2lpdiff); cout << "." << endl;
  cout << "# MH2SQ: " << endl;
  cout << "# SOFTSUSY (approx) beta_mHuSq^(1)   = "; printRow(cout, betamHuSq1lp); cout << " ("; 
  printRow(cout, approxmHuSqbeta1lp); 
  cout << " ), diff. = "; printRow(cout, mHuSq1lpdiff); cout << "." << endl;
  cout << "# SOFTSUSY (approx) beta_mHuSq^(2)   = "; printRow(cout, betamHuSq2lppiece); cout << " ("; 
  printRow(cout, approxmHuSqbeta2lppiece);
  cout << " ), diff. = "; printRow(cout, mHuSq2lppiecediff); cout << "." << endl;
  cout << "# SOFTSUSY (approx) beta_mHuSq^(1+2) = "; printRow(cout, betamHuSq2lp); cout << " ("; 
  printRow(cout, approxmHuSqbeta2lp); 
  cout << " ), diff. = "; printRow(cout, mHuSq2lpdiff); cout << "." << endl;
  cout << "# MQL3SQ: " << endl;
  cout << "# SOFTSUSY (approx) beta_mQlSq^(1)   = "; printRow(cout, betamQlSq1lp); cout << " ("; 
  printRow(cout, approxmQlSqbeta1lp); 
  cout << " ), diff. = "; printRow(cout, mQlSq1lpdiff); cout << "." << endl;
  cout << "# SOFTSUSY (approx) beta_mQlSq^(2)   = "; printRow(cout, betamQlSq2lppiece); cout << " ("; 
  printRow(cout, approxmQlSqbeta2lppiece);
  cout << " ), diff. = "; printRow(cout, mQlSq2lppiecediff); cout << "." << endl;
  cout << "# SOFTSUSY (approx) beta_mQlSq^(1+2) = "; printRow(cout, betamQlSq2lp); cout << " ("; 
  printRow(cout, approxmQlSqbeta2lp); 
  cout << " ), diff. = "; printRow(cout, mQlSq2lpdiff); cout << "." << endl;
  cout << "# MTRSQ: " << endl;
  cout << "# SOFTSUSY (approx) beta_mUrSq^(1)   = "; printRow(cout, betamUrSq1lp); cout << " ("; 
  printRow(cout, approxmUrSqbeta1lp); 
  cout << " ), diff. = "; printRow(cout, mUrSq1lpdiff); cout << "." << endl;
  cout << "# SOFTSUSY (approx) beta_mUrSq^(2)   = "; printRow(cout, betamUrSq2lppiece); cout << " ("; 
  printRow(cout, approxmUrSqbeta2lppiece);
  cout << " ), diff. = "; printRow(cout, mUrSq2lppiecediff); cout << "." << endl;
  cout << "# SOFTSUSY (approx) beta_mUrSq^(1+2) = "; printRow(cout, betamUrSq2lp); cout << " ("; 
  printRow(cout, approxmUrSqbeta2lp); 
  cout << " ), diff. = "; printRow(cout, mUrSq2lpdiff); cout << "." << endl;
  cout << "# AT: " << endl;
  cout << "# SOFTSUSY (approx) beta_At^(1)   = "; printRow(cout, betaAt1lp); cout << " ("; 
  printRow(cout, approxAtbeta1lp); 
  cout << " ), diff. = "; printRow(cout, At1lpdiff); cout << "." << endl;
  cout << "# SOFTSUSY (approx) beta_At^(2)   = "; printRow(cout, betaAt2lppiece); cout << " ("; 
  printRow(cout, approxAtbeta2lppiece);
  cout << " ), diff. = "; printRow(cout, At2lppiecediff); cout << "." << endl;
  cout << "# SOFTSUSY (approx) beta_At^(1+2) = "; printRow(cout, betaAt2lp); cout << " ("; 
  printRow(cout, approxAtbeta2lp); 
  cout << " ), diff. = "; printRow(cout, At2lpdiff); cout << "." << endl;
  cout << "#" << endl;


  // Error checking - if differences exceed particular threshold, 
  // print warning
  double tol = 1.0e-9;
  if (fabs(mu1lpdiff) > tol)
    {
      cout << "# WARNING: 1-loop beta functions for MU disagree, diff. = " << mu1lpdiff << "." << endl;
    }
  if (fabs(mu2lppiecediff) > tol)
    {
      cout << "# WARNING: 2-loop pieces of beta functions for MU disagree, diff. = " << mu2lppiecediff << "." << endl;
    }
  if (fabs(mu2lpdiff) > tol)
    {
      cout << "# WARNING: 2-loop beta functions for MU disagree, diff. = " << mu2lpdiff << "." << endl;
    }
  if (fabs(B1lpdiff) > tol)
    {
      cout << "# WARNING: 1-loop beta functions for B disagree, diff. = " << B1lpdiff << "." << endl;
    }
  if (fabs(B2lppiecediff) > tol)
    {
      cout << "# WARNING: 2-loop pieces of beta functions for B disagree, diff. = " << B2lppiecediff << "." << endl;
    }
  if (fabs(B2lpdiff) > tol)
    {
      cout << "# WARNING: 2-loop beta functions for B disagree, diff. = " << B2lpdiff << "." << endl;
    }
  if (fabs(mHdSq1lpdiff) > tol)
    {
      cout << "# WARNING: 1-loop beta functions for MH1SQ disagree, diff. = " << mHdSq1lpdiff << "." << endl;
    }
  if (fabs(mHdSq2lppiecediff) > tol)
    {
      cout << "# WARNING: 2-loop pieces of beta functions for MH1SQ disagree, diff. = " << mHdSq2lppiecediff << "." << endl;
    }
  if (fabs(mHdSq2lpdiff) > tol)
    {
      cout << "# WARNING: 2-loop beta functions for MH1SQ disagree, diff. = " << mHdSq2lpdiff << "." << endl;
    }
  if (fabs(mHuSq1lpdiff) > tol)
    {
      cout << "# WARNING: 1-loop beta functions for MH2SQ disagree, diff. = " << mHuSq1lpdiff << "." << endl;
    }
  if (fabs(mHuSq2lppiecediff) > tol)
    {
      cout << "# WARNING: 2-loop pieces of beta functions for MH2SQ disagree, diff. = " << mHuSq2lppiecediff << "." << endl;
    }
  if (fabs(mHuSq2lpdiff) > tol)
    {
      cout << "# WARNING: 2-loop beta functions for MH2SQ disagree, diff. = " << mHuSq2lpdiff << "." << endl;
    }
  if (fabs(mQlSq1lpdiff) > tol)
    {
      cout << "# WARNING: 1-loop beta functions for MQL3SQ disagree, diff. = " << mQlSq1lpdiff << "." << endl;
    }
  if (fabs(mQlSq2lppiecediff) > tol)
    {
      cout << "# WARNING: 2-loop pieces of beta functions for MQL3SQ disagree, diff. = " << mQlSq2lppiecediff << "." << endl;
    }
  if (fabs(mQlSq2lpdiff) > tol)
    {
      cout << "# WARNING: 2-loop beta functions for MQL3SQ disagree, diff. = " << mQlSq2lpdiff << "." << endl;
    }
  if (fabs(mUrSq1lpdiff) > tol)
    {
      cout << "# WARNING: 1-loop beta functions for MTRSQ disagree, diff. = " << mUrSq1lpdiff << "." << endl;
    }
  if (fabs(mUrSq2lppiecediff) > tol)
    {
      cout << "# WARNING: 2-loop pieces of beta functions for MTRSQ disagree, diff. = " << mUrSq2lppiecediff << "." << endl;
    }
  if (fabs(mUrSq2lpdiff) > tol)
    {
      cout << "# WARNING: 2-loop beta functions for MTRSQ disagree, diff. = " << mUrSq2lpdiff << "." << endl;
    }
  if (fabs(At1lpdiff) > tol)
    {
      cout << "# WARNING: 1-loop beta functions for AT disagree, diff. = " << At1lpdiff << "." << endl;
    }
  if (fabs(At2lppiecediff) > tol)
    {
      cout << "# WARNING: 2-loop pieces of beta functions for AT disagree, diff. = " << At2lppiecediff << "." << endl;
    }
  if (fabs(At2lpdiff) > tol)
    {
      cout << "# WARNING: 2-loop beta functions for AT disagree, diff. = " << At2lpdiff << "." << endl;
    }

  // Now we need to check the 1-loop squared pieces by differentiating the 1-loop beta
  // functions provided by SOFTSUSY numerically.

  // First get the values for our approximate expressions.
  double muLogSqCoeff = doCalcMuLogSqCoeff(r, 1);
  double BLogSqCoeff = doCalcBLogSqCoeff(r, 1);
  double mHdSqLogSqCoeff = doCalcMh1SquaredLogSqCoeff(r, 1);
  double mHuSqLogSqCoeff = doCalcMh2SquaredLogSqCoeff(r, 1);
  double mQlSqLogSqCoeff = doCalcmqL3SquaredLogSqCoeff(r, 1);
  double mUrSqLogSqCoeff = doCalcmtRSquaredLogSqCoeff(r, 1);
  double AtLogSqCoeff = doCalcAtLogSqCoeff(r, 1);

  // Use SOFTSUSY's calcDerivative routine to compute the 
  // derivative w.r.t. t of the 1-loop beta functions.
  double h = 1.0e-3;
  double err;
  q0 = r.displayMu();

  r.setLoops(1);

  tempRunner = &r;

  double numericalmuLogSqCoeff = 0.5*calcDerivative(getMuBeta1, 0.0, h, &err);
  double numericalBLogSqCoeff = 0.5*calcDerivative(getBBeta1, 0.0, h, &err);
  double numericalmHdSqLogSqCoeff = 0.5*calcDerivative(getmHdSqBeta1, 0.0, h, &err);
  double numericalmHuSqLogSqCoeff = 0.5*calcDerivative(getmHuSqBeta1, 0.0, h, &err);
  double numericalmQlSqLogSqCoeff = 0.5*calcDerivative(getmQlSqBeta1, 0.0, h, &err);
  double numericalmUrSqLogSqCoeff = 0.5*calcDerivative(getmUrSqBeta1, 0.0, h, &err);
  double numericalAtLogSqCoeff = 0.5*calcDerivative(getAtBeta1, 0.0, h, &err);

  // Work out the percent relative errors
  double muRelErr = 100.0*fabs((numericalmuLogSqCoeff-muLogSqCoeff)/(0.5*(numericalmuLogSqCoeff+muLogSqCoeff)));
  double BRelErr = 100.0*fabs((numericalBLogSqCoeff-BLogSqCoeff)/(0.5*(numericalBLogSqCoeff+BLogSqCoeff)));
  double mHdSqRelErr = 100.0*fabs((numericalmHdSqLogSqCoeff-mHdSqLogSqCoeff)/(0.5*(numericalmHdSqLogSqCoeff+mHdSqLogSqCoeff)));
  double mHuSqRelErr = 100.0*fabs((numericalmHuSqLogSqCoeff-mHuSqLogSqCoeff)/(0.5*(numericalmHuSqLogSqCoeff+mHuSqLogSqCoeff)));
  double mQlSqRelErr = 100.0*fabs((numericalmQlSqLogSqCoeff-mQlSqLogSqCoeff)/(0.5*(numericalmQlSqLogSqCoeff+mQlSqLogSqCoeff)));
  double mUrSqRelErr = 100.0*fabs((numericalmUrSqLogSqCoeff-mUrSqLogSqCoeff)/(0.5*(numericalmUrSqLogSqCoeff+mUrSqLogSqCoeff)));
  double AtRelErr = 100.0*fabs((numericalAtLogSqCoeff-AtLogSqCoeff)/(0.5*(numericalAtLogSqCoeff+AtLogSqCoeff)));

  cout << "# SOFTSUSY (approx) MU log^2 coeff = "; printRow(cout, numericalmuLogSqCoeff); cout << " ("; 
  printRow(cout,muLogSqCoeff); cout << " ), % error = "; printRow(cout, muRelErr); cout << " %." << endl;
  cout << "# SOFTSUSY (approx) B log^2 coeff = "; printRow(cout, numericalBLogSqCoeff); cout << " (";
  printRow(cout, BLogSqCoeff); cout << " ), % error = "; printRow(cout, BRelErr); cout << " %."  << endl;
  cout << "# SOFTSUSY (approx) MH1SQ log^2 coeff = "; printRow(cout, numericalmHdSqLogSqCoeff); cout << " (";
  printRow(cout, mHdSqLogSqCoeff); cout << " ), % error = "; printRow(cout, mHdSqRelErr); cout << " %." << endl;
  cout << "# SOFTSUSY (approx) MH2SQ log^2 coeff = "; printRow(cout, numericalmHuSqLogSqCoeff); cout << " (";
  printRow(cout, mHuSqLogSqCoeff); cout << " ), % error = "; printRow(cout, mHuSqRelErr); cout << " %." << endl;
  cout << "# SOFTSUSY (approx) MQL3SQ log^2 coeff = "; printRow(cout, numericalmQlSqLogSqCoeff); cout << " (";
  printRow(cout, mQlSqLogSqCoeff); cout << " ), % error = "; printRow(cout, mQlSqRelErr); cout << " %." << endl;
  cout << "# SOFTSUSY (approx) MTRSQ log^2 coeff = "; printRow(cout, numericalmUrSqLogSqCoeff); cout << " (";
  printRow(cout, mUrSqLogSqCoeff); cout << " ), % error = "; printRow(cout, mUrSqRelErr); cout << " %." << endl;
  cout << "# SOFTSUSY (approx) AT log^2 coeff = "; printRow(cout, numericalAtLogSqCoeff); cout << " (";
  printRow(cout, AtLogSqCoeff); cout << " ), % error = "; printRow(cout, AtRelErr); cout << " %." << endl;
  cout << "#" << endl;

  double logSqTol = 1.0e-8;
  if (muRelErr > logSqTol)
    {
      cout << "# WARNING: log^2 coefficients for MU disagree, % error = " << muRelErr << " %." << endl;
    }
  if (BRelErr > logSqTol)
    {
      cout << "# WARNING: log^2 coefficients for B disagree, % error = " << BRelErr << " %." << endl;
    }
  if (mHdSqRelErr > logSqTol)
    {
      cout << "# WARNING: log^2 coefficients for MH1SQ disagree, % error = " << mHdSqRelErr << " %." << endl;
    }
  if (mHuSqRelErr > logSqTol)
    {
      cout << "# WARNING: log^2 coefficients for MH2SQ disagree, % error = " << mHuSqRelErr << " %." << endl;
    }
  if (mQlSqRelErr > logSqTol)
    {
      cout << "# WARNING: log^2 coefficients for MQL3SQ disagree, % error = " << mQlSqRelErr << " %." << endl;
    }
  if (mUrSqRelErr > logSqTol)
    {
      cout << "# WARNING: log^2 coefficients for MTRSQ disagree, % error = " << mUrSqRelErr << " %." << endl;
    }
  if (AtRelErr > logSqTol)
    {
      cout << "# WARNING: log^2 coefficients for AT disagree, % error = " << AtRelErr << " %." << endl;
    }


}

// Functions for getting approximate derivatives of low scale parameters w.r.t high
// scale parameters. Returns the vector
// [ dmu/dp dB/dp dm_Hd^2/dp dm_Hu^2/dp dm_Ql^2/dp dm_uR^2/dp dA_t/dp]^T. The vector auxPars
// is assumed to contain the values of the gauge and Yukawa couplings at MX, in the order
// [ g1 g2 g3 yt yb ytau]^T.
DoubleVector doCalcMuDerivs(SoftParsMssm r, DoubleVector pars, double mx, bool & hasError,
			    DoubleVector const & auxPars)
{
  // Get M_{SUSY}
  double ms = r.displayMu();

  // Taylor approximation to the RGE solution is assumed to be about the values at MX.
  double t = log(ms/mx);

  DoubleVector derivs(7);

  double yt, yb, ytau, g1, g2, g3;
  // Check that all of the necessary couplings have been provided. If not, we can
  // still get them by running to MX, but we want to avoid this if possible.
  if (auxPars.displayEnd()-auxPars.displayStart()+1 != 6)
    {
      cerr << "WARNING: incorrect number of couplings provided to doCalcMuDerivs: running to get values." << endl;

      // Copy the original object to do the running.
      SoftParsMssm tempObject = r;

      tempObject.runto(mx);

      yt = tempObject.displayYukawaElement(YU, 3, 3);
      yb = tempObject.displayYukawaElement(YD, 3, 3);
      ytau = tempObject.displayYukawaElement(YE, 3, 3);
      g1 = tempObject.displayGaugeCoupling(1);
      g2 = tempObject.displayGaugeCoupling(2);
      g3 = tempObject.displayGaugeCoupling(3);
    }
  else
    {
      g1 = auxPars(1);
      g2 = auxPars(2);
      g3 = auxPars(3);
      yt = auxPars(4);
      yb = auxPars(5);
      ytau = auxPars(6);
    }

  // The vector pars contains the values of the high scale parameters. Extract
  // the needed values.

  double dmudp, dBdp, dmHdSqdp, dmHuSqdp, dmQlSqdp, dmUrSqdp, dAtdp;

  dBdp = 0.0;
  dmHdSqdp = 0.0;
  dmHuSqdp = 0.0;
  dmQlSqdp = 0.0;
  dmUrSqdp = 0.0;
  dAtdp = 0.0;

  dmudp = 1.0+(t/(16.0*PI*PI))*(3.0*yt*yt+3.0*yb*yb+ytau*ytau-3.0*g2*g2-0.6*g1*g1
				+(1.0/(16.0*PI*PI))*(-9.0*sqr(yt*yt)-9.0*sqr(yb*yb)-6.0*sqr(yt*yb)
						     -3.0*sqr(ytau*ytau)+16.0*sqr(g3*yt)+0.8*sqr(g1*yt)
						     +16.0*sqr(g3*yb)-0.4*sqr(g1*yb)+1.2*sqr(g1*ytau)
						     +7.5*sqr(g2*g2)+1.8*sqr(g1*g2)+(207.0/50.0)*sqr(g1*g1)))
    +0.5*sqr(t/(16.0*PI*PI))*(45.0*sqr(yt*yt)+45.0*sqr(yb*yb)+9.0*sqr(ytau*ytau)+30.0*sqr(yt*yb)+6.0*sqr(yt*ytau)
			      +18.0*sqr(yb*ytau)-32.0*g3*g3*(yt*yt+yb*yb)-12.0*g2*g2*(3.0*yt*yt+3.0*yb*yb+ytau*ytau)
			      -0.8*g1*g1*(11.0*yt*yt+8.0*yb*yb+6.0*ytau*ytau)+3.0*sqr(g2*g2)-(189.0/25.0)*sqr(g1*g1)
			      +3.6*sqr(g1*g2));

  derivs(1) = dmudp;
  derivs(2) = dBdp;
  derivs(3) = dmHdSqdp;
  derivs(4) = dmHuSqdp;
  derivs(5) = dmQlSqdp;
  derivs(6) = dmUrSqdp;
  derivs(7) = dAtdp;

  return derivs;
}

DoubleVector doCalcBDerivs(SoftParsMssm r, DoubleVector pars, double mx, bool & hasError,
			    DoubleVector const & auxPars)
{
  // Get M_{SUSY}
  double ms = r.displayMu();

  // Taylor approximation to the RGE solution is assumed to be about the values at MX.
  double t = log(ms/mx);

  DoubleVector derivs(7);

  double yt, yb, ytau, g1, g2, g3;
  // Check that all of the necessary couplings have been provided. If not, we can
  // still get them by running to MX, but we want to avoid this if possible.
  if (auxPars.displayEnd()-auxPars.displayStart()+1 != 6)
    {
      cerr << "WARNING: incorrect number of couplings provided to doCalcMuDerivs: running to get values." << endl;

      // Copy the original object to do the running.
      SoftParsMssm tempObject = r;

      tempObject.runto(mx);

      yt = tempObject.displayYukawaElement(YU, 3, 3);
      yb = tempObject.displayYukawaElement(YD, 3, 3);
      ytau = tempObject.displayYukawaElement(YE, 3, 3);
      g1 = tempObject.displayGaugeCoupling(1);
      g2 = tempObject.displayGaugeCoupling(2);
      g3 = tempObject.displayGaugeCoupling(3);
    }
  else
    {
      g1 = auxPars(1);
      g2 = auxPars(2);
      g3 = auxPars(3);
      yt = auxPars(4);
      yb = auxPars(5);
      ytau = auxPars(6);
    }

  double dmudp, dBdp, dmHdSqdp, dmHuSqdp, dmQlSqdp, dmUrSqdp, dAtdp;

  dBdp = 1.0;
  dmHdSqdp = 0.0;
  dmHuSqdp = 0.0;
  dmQlSqdp = 0.0;
  dmUrSqdp = 0.0;
  dAtdp = 0.0;
  dmudp = 0.0;

  derivs(1) = dmudp;
  derivs(2) = dBdp;
  derivs(3) = dmHdSqdp;
  derivs(4) = dmHuSqdp;
  derivs(5) = dmQlSqdp;
  derivs(6) = dmUrSqdp;
  derivs(7) = dAtdp;

  return derivs;
}

DoubleVector doCalcMh1SquaredDerivs(SoftParsMssm r, DoubleVector pars, double mx, 
				    bool & hasError, DoubleVector const & auxPars)
{
  // Get M_{SUSY}
  double ms = r.displayMu();

  // Taylor approximation to the RGE solution is assumed to be about the values at MX.
  double t = log(ms/mx);

  DoubleVector derivs(7);

  double yt, yb, ytau, g1, g2, g3;
  // Check that all of the necessary couplings have been provided. If not, we can
  // still get them by running to MX, but we want to avoid this if possible.
  if (auxPars.displayEnd()-auxPars.displayStart()+1 != 6)
    {
      cerr << "WARNING: incorrect number of couplings provided to doCalcMuDerivs: running to get values." << endl;

      // Copy the original object to do the running.
      SoftParsMssm tempObject = r;

      tempObject.runto(mx);

      yt = tempObject.displayYukawaElement(YU, 3, 3);
      yb = tempObject.displayYukawaElement(YD, 3, 3);
      ytau = tempObject.displayYukawaElement(YE, 3, 3);
      g1 = tempObject.displayGaugeCoupling(1);
      g2 = tempObject.displayGaugeCoupling(2);
      g3 = tempObject.displayGaugeCoupling(3);
    }
  else
    {
      g1 = auxPars(1);
      g2 = auxPars(2);
      g3 = auxPars(3);
      yt = auxPars(4);
      yb = auxPars(5);
      ytau = auxPars(6);
    }

  // The vector pars contains the values of the high scale parameters. Extract
  // the needed values.

  double dmudp, dBdp, dmHdSqdp, dmHuSqdp, dmQlSqdp, dmUrSqdp, dAtdp;

  dBdp = 0.0;
  dAtdp = 0.0;
  dmudp = 0.0;

  dmHdSqdp = 1.0+(t/(16.0*PI*PI))*(6.0*yb*yb+2.0*ytau*ytau+0.6*g1*g1
				   +(1.0/(16.0*PI*PI))*(-36.0*sqr(yb*yb)-6.0*sqr(yt*yb)-12.0*sqr(ytau*ytau)
							+32.0*sqr(g3*yb)-4.4*sqr(g1*yb)+1.2*sqr(g1*ytau)+3.0*sqr(g2*g2)
							+1.8*sqr(g1*g2)+(18.0/25.0)*sqr(g1*g1)))
    +sqr(t/(16.0*PI*PI))*(72.0*sqr(yb*yb)+6.0*sqr(yt*yb)+24.0*sqr(ytau*yb)+16.0*sqr(ytau*ytau)-32.0*sqr(g3*yb)
			  -18.0*sqr(g2*yb)-6.0*sqr(g2*ytau)-2.8*sqr(g1*yb)-3.6*sqr(g1*ytau)+(198.0/25.0)*sqr(g1*g1));

  dmHuSqdp = t/(16.0*PI*PI)*(-0.6*g1*g1+(1.0/(16.0*PI*PI))*(-6.0*sqr(yt*yb)+3.6*sqr(g1*yb)+1.2*sqr(g1*ytau)
							    -1.8*sqr(g1*g2)+3.0*sqr(g2*g2)))
    +sqr(t/(16.0*PI*PI))*(6.0*sqr(yt*yb)-(198.0/25.0)*sqr(g1*g1));

  dmQlSqdp = (t/(16.0*PI*PI))*(2.0*yb*yb-0.2*g1*g1+(1.0/(16.0*PI*PI))*(-20.0*sqr(yb*yb)-4.0*sqr(ytau*yb)+2.0*sqr(g1*yb)
								       +0.4*sqr(g1*ytau)+3.0*sqr(g2*g2)-0.6*sqr(g1*g2)
								       -(6.0/75.0)*sqr(g1*g1)))
    +sqr(t/(16.0*PI*PI))*(24.0*sqr(yb*yb)+4.0*sqr(yt*yb)+4.0*sqr(yb*ytau)-(32.0/3.0)*sqr(g3*yb)-6.0*sqr(g2*yb)
			  -(14.0/15.0)*sqr(g1*yb)-(66.0/25.0)*sqr(g1*g1));

  dmUrSqdp = (t/(16.0*PI*PI))*(0.8*g1*g1+(1.0/(16.0*PI*PI))*(-4.0*sqr(yt*yb)-4.8*sqr(g1*yb)-1.6*sqr(g1*ytau)
							     +2.4*sqr(g1*g2)+(84.0/75.0)*sqr(g1*g1)))
    +sqr(t/(16.0*PI*PI))*(4.0*sqr(yt*yb)+(264.0/25.0)*sqr(g1*g1));

  derivs(1) = dmudp;
  derivs(2) = dBdp;
  derivs(3) = dmHdSqdp;
  derivs(4) = dmHuSqdp;
  derivs(5) = dmQlSqdp;
  derivs(6) = dmUrSqdp;
  derivs(7) = dAtdp;

  return derivs;
}

DoubleVector doCalcMh2SquaredDerivs(SoftParsMssm r, DoubleVector pars, double mx, 
				    bool & hasError, DoubleVector const & auxPars)
{
  // Get M_{SUSY}
  double ms = r.displayMu();

  // Taylor approximation to the RGE solution is assumed to be about the values at MX.
  double t = log(ms/mx);

  DoubleVector derivs(7);

  double yt, yb, ytau, g1, g2, g3;
  // Check that all of the necessary couplings have been provided. If not, we can
  // still get them by running to MX, but we want to avoid this if possible.
  if (auxPars.displayEnd()-auxPars.displayStart()+1 != 6)
    {
      cerr << "WARNING: incorrect number of couplings provided to doCalcMuDerivs: running to get values." << endl;

      // Copy the original object to do the running.
      SoftParsMssm tempObject = r;

      tempObject.runto(mx);

      yt = tempObject.displayYukawaElement(YU, 3, 3);
      yb = tempObject.displayYukawaElement(YD, 3, 3);
      ytau = tempObject.displayYukawaElement(YE, 3, 3);
      g1 = tempObject.displayGaugeCoupling(1);
      g2 = tempObject.displayGaugeCoupling(2);
      g3 = tempObject.displayGaugeCoupling(3);
    }
  else
    {
      g1 = auxPars(1);
      g2 = auxPars(2);
      g3 = auxPars(3);
      yt = auxPars(4);
      yb = auxPars(5);
      ytau = auxPars(6);
    }

  // The vector pars contains the values of the high scale parameters. Extract
  // the needed values.

  double dmudp, dBdp, dmHdSqdp, dmHuSqdp, dmQlSqdp, dmUrSqdp, dAtdp;

  dBdp = 0.0;
  dAtdp = 0.0;
  dmudp = 0.0;

  dmHdSqdp = (t/(16.0*PI*PI))*(-0.6*g1*g1+(1.0/(16.0*PI*PI))*(-6.0*sqr(yt*yb)+3.6*sqr(g1*yt)+3.0*sqr(g2*g2)-1.8*sqr(g1*g2)))
    +sqr(t/(16.0*PI*PI))*(6.0*sqr(yt*yb)-(198.0/25.0)*sqr(g1*g1));

  dmHuSqdp = 1.0+(t/(16.0*PI*PI))*(6.0*yt*yt+0.6*g1*g1+(1.0/(16.0*PI*PI))*(-36.0*sqr(yt*yt)-6.0*sqr(yt*yb)+32.0*sqr(g3*yt)
									   -2.0*sqr(g1*yt)+1.8*sqr(g1*g2)+3.0*sqr(g2*g2)
									   +(18.0/25.0)*sqr(g1*g1)))
    +sqr(t/(16.0*PI*PI))*(72.0*sqr(yt*yt)+6.0*sqr(yt*yb)-32.0*sqr(g3*yt)-18.0*sqr(g2*yt)-5.2*sqr(g1*yt)+(198.0/25.0)*sqr(g1*g1));

  dmQlSqdp = (t/(16.0*PI*PI))*(2.0*yt*yt+0.2*g1*g1+(1.0/(16.0*PI*PI))*(-20.0*sqr(yt*yt)+0.4*sqr(g1*yt)+3.0*sqr(g2*g2)
								       +0.6*sqr(g1*g2)+(12.0/75.0)*sqr(g1*g1)))
    +sqr(t/(16.0*PI*PI))*(24.0*sqr(yt*yt)+4.0*sqr(yt*yb)-(32.0/3.0)*sqr(g3*yt)-6.0*sqr(g2*yt)-(26.0/15.0)*sqr(g1*yt)
			  +(66.0/25.0)*sqr(g1*g1));


  dmUrSqdp = (t/(16.0*PI*PI))*(4.0*yt*yt-0.8*g1*g1+(1.0/(16.0*PI*PI))*(-32.0*sqr(yt*yt)-4.0*sqr(yt*yb)+12.0*sqr(g2*yt)
								       +4.0*sqr(g1*yt)-2.4*sqr(g1*g2)+(12.0/75.0)*sqr(g1*g1)))
    +sqr(t/(16.0*PI*PI))*(48.0*sqr(yt*yt)+4.0*sqr(yt*yb)-(64.0/3.0)*sqr(g3*yt)-12.0*sqr(g2*yt)-(52.0/15.0)*sqr(g1*yt)
			  -(264.0/25.0)*sqr(g1*g1));

  derivs(1) = dmudp;
  derivs(2) = dBdp;
  derivs(3) = dmHdSqdp;
  derivs(4) = dmHuSqdp;
  derivs(5) = dmQlSqdp;
  derivs(6) = dmUrSqdp;
  derivs(7) = dAtdp;

  return derivs;
}

DoubleVector doCalcGauginoDerivs(SoftParsMssm r, DoubleVector pars, double mx, 
				 bool & hasError, DoubleVector const & auxPars, int whichGaugino)
{
  DoubleVector derivs(7);

  // Check that the provided index is valid
  if (whichGaugino < 1 | whichGaugino > 3)
    {
      ostringstream ii;
      ii << "# WARNING: invalid index i = " << whichGaugino << " for gaugino soft mass M_i: exiting." << endl;
      throw ii.str();
    }
  else
    {
      // Get M_{SUSY}
      double ms = r.displayMu();
      
      // Taylor approximation to the RGE solution is assumed to be about the values at MX.
      double t = log(ms/mx);
      
      double yt, yb, ytau, g1, g2, g3;
      // Check that all of the necessary couplings have been provided. If not, we can
      // still get them by running to MX, but we want to avoid this if possible.
      if (auxPars.displayEnd()-auxPars.displayStart()+1 != 6)
	{
	  cerr << "WARNING: incorrect number of couplings provided to doCalcMuDerivs: running to get values." << endl;
	  
	  // Copy the original object to do the running.
	  SoftParsMssm tempObject = r;
	  
	  tempObject.runto(mx);
	  
	  yt = tempObject.displayYukawaElement(YU, 3, 3);
	  yb = tempObject.displayYukawaElement(YD, 3, 3);
	  ytau = tempObject.displayYukawaElement(YE, 3, 3);
	  g1 = tempObject.displayGaugeCoupling(1);
	  g2 = tempObject.displayGaugeCoupling(2);
	  g3 = tempObject.displayGaugeCoupling(3);
	}
      else
	{
	  g1 = auxPars(1);
	  g2 = auxPars(2);
	  g3 = auxPars(3);
	  yt = auxPars(4);
	  yb = auxPars(5);
	  ytau = auxPars(6);
	}
      
      // The vector pars contains the values of the high scale parameters. Extract
      // the needed values.
      double M1 = pars(1);
      double M2 = pars(2);
      double M3 = pars(3);
      double At = pars(4);
      double Ab = pars(5);
      double Atau = pars(6);

      double dmudp, dBdp, dmHdSqdp, dmHuSqdp, dmQlSqdp, dmUrSqdp, dAtdp;
      
      // This is zero for all gaugino soft masses
      dmudp = 0.0;

      switch (whichGaugino)
	{
	case 1:
	  {
	    dBdp = (t/(16.0*PI*PI))*(1.2*g1*g1+(1.0/(16.0*PI*PI))*(-1.6*sqr(g1*yt)+0.8*sqr(g1*yb)-2.4*sqr(g1*ytau)
								   -3.6*sqr(g1*g2)-(414.0/25.0)*sqr(g1*g1)))
	      +sqr(t/(16.0*PI*PI))*(5.2*sqr(g1*yt)+2.8*sqr(g1*yb)+3.6*sqr(g1*ytau)+(396.0/25.0)*sqr(g1*g1));

	    dAtdp = (t/(16.0*PI*PI))*((26.0/15.0)*g1*g1+(1.0/(16.0*PI*PI))*(-2.4*sqr(g1*yt)-0.8*sqr(g1*yb)-(272.0/45.0)*sqr(g3*g1)
									    -2.0*sqr(g2*g1)-(5486.0/225.0)*sqr(g1*g1)))
	      +sqr(t/(16.0*PI*PI))*(10.4*sqr(g1*yt)+(14.0/15.0)*sqr(g1*yb)+(572.0/25.0)*sqr(g1*g1));
	    
	    dmHdSqdp = (t/(16.0*PI*PI))*(-2.4*g1*g1*M1+(1.0/(16.0*PI*PI))*(-1.6*sqr(g1*yb)*(2.0*M1-Ab)
									   +4.8*sqr(g1*ytau)*(2.0*M1-Atau)
									   +3.6*sqr(g1*g2)*(2.0*M1+M2)
									   +(1242.0/25.0)*sqr(g1*g1)*M1))
	      +sqr(t/(16.0*PI*PI))*(-5.6*sqr(g1*yb)*(2.0*M1-Ab)-7.2*sqr(g1*ytau)*(2.0*M1-Atau)-(1188.0/25.0)*sqr(g1*g1)*M1);
	    
	    dmHuSqdp = (t/(16.0*PI*PI))*(-2.4*g1*g1*M1+(1.0/(16.0*PI*PI))*(3.2*sqr(g1*yt)*(2.0*M1-At)+3.6*sqr(g1*g2)*(2.0*M1+M2)
									   +(1242.0/25.0)*sqr(g1*g1)*M1))
	      +sqr(t/(16.0*PI*PI))*(-10.4*sqr(g1*yt)*(2.0*M1-At)-(1188.0/25.0)*sqr(g1*g1)*M1);
	    
	    dmQlSqdp = (t/(16.0*PI*PI))*(-(4.0/15.0)*g1*g1*M1+(1.0/(16.0*PI*PI))*(3.2*sqr(g1*yt)*(2.0*M1-At)
										  +1.6*sqr(g1*yb)*(2.0*M1-Ab)
										  +(32.0/45.0)*sqr(g3*g1)*(2.0*M1+M3)
										  +0.4*sqr(g1*g2)*(2.0*M1+M2)
										  +(398.0/75.0)*sqr(g1*g1)*M1))
	      +sqr(t/(16.0*PI*PI))*(-(52.0/15.0)*sqr(g1*yt)*(2.0*M1-At)-(28.0/15.0)*sqr(g1*yb)*(2.0*M1-Ab)
				    -(132.0/25.0)*sqr(g1*g1)*M1);
	    
	    dmUrSqdp = (t/(16.0*PI*PI))*(-(64.0/15.0)*g1*g1*M1+(1.0/(16.0*PI*PI))*(-1.6*sqr(g1*yt)*(2.0*M1-At)
										   +(512.0/45.0)*sqr(g3*g1)*(2.0*M1+M3)
										   +(6848.0/75.0)*sqr(g1*g1)*M1))
	      +sqr(t/(16.0*PI*PI))*(-(104.0/15.0)*sqr(g1*yt)*(2.0*M1-At)-(2112.0/25.0)*sqr(g1*g1)*M1);

	    break;
	  }
	case 2:
	  {
	    dBdp = (t/(16.0*PI*PI))*(6.0*g2*g2+(1.0/(16.0*PI*PI))*(-30.0*sqr(g2*g2)-3.6*sqr(g1*g2)))
	      +sqr(t/(16.0*PI*PI))*(18.0*sqr(g2*yt)+18.0*sqr(g2*yb)+6.0*sqr(g2*ytau)+12.0*sqr(g2*g2));

	    dAtdp = (t/(16.0*PI*PI))*(6.0*g2*g2+(1.0/(16.0*PI*PI))*(-12.0*sqr(g2*yt)-16.0*sqr(g3*g2)-30.0*sqr(g2*g2)
								    -2.0*sqr(g1*g2)))
	      +sqr(t/(16.0*PI*PI))*(36.0*sqr(g2*yt)+6.0*sqr(g2*yb)+12.0*sqr(g2*g2));
	    
	    dmHdSqdp = (t/(16.0*PI*PI))*(-12.0*g2*g2*M2+(1.0/(16.0*PI*PI))*(66.0*sqr(g2*g2)*M2+3.6*sqr(g1*g2)*(2.0*M2+M1)))
	      +sqr(t/(16.0*PI*PI))*(-36.0*sqr(g2*yb)*(2.0*M2-Ab)-12.0*sqr(g2*ytau)*(2.0*M2-Atau)-36.0*sqr(g2*g2)*M2);
	    
	    dmHuSqdp = (t/(16.0*PI*PI))*(-12.0*g2*g2*M2+(1.0/(16.0*PI*PI))*(3.6*sqr(g1*g2)*(2.0*M2+M1)+66.0*sqr(g2*g2)*M2))
	      +sqr(t/(16.0*PI*PI))*(-36.0*sqr(g2*yt)*(2.0*M2-At)-36.0*sqr(g2*g2)*M2);
	    
	    dmQlSqdp = (t/(16.0*PI*PI))*(-12.0*g2*g2*M2+(1.0/(16.0*PI*PI))*(32.0*sqr(g3*g2)*(2.0*M2+M3)+66.0*sqr(g2*g2)*M2
									    +0.4*sqr(g1*g2)*(2.0*M2+M1)))
	      +sqr(t/(16.0*PI*PI))*(-12.0*sqr(g2*yt)*(2.0*M2-At)-12.0*sqr(g2*yb)*(2.0*M2-Ab)-36.0*sqr(g2*g2)*M2);
	    
	    dmUrSqdp = 24.0*sqr(g2*yt)*(2.0*M2-At)*t*(1-t)/sqr(16.0*PI*PI);

	    break;
	  }
	case 3:
	  {
	    dBdp = 32.0*sqr(g3*g3)*(yt*yt+yb*yb)*t*(t-1)/sqr(16.0*PI*PI);

	    dAtdp = (t/(16.0*PI*PI))*((32.0/3.0)*g3*g3+(1.0/(16.0*PI*PI))*(-32.0*sqr(g3*yt)+(64.0/9.0)*sqr(g3*g3)
									   -16.0*sqr(g3*g2)-(272.0/45.0)*sqr(g3*g1)))
	      +sqr(t/(16.0*PI*PI))*(64.0*sqr(g3*yt)+(32.0/3.0)*sqr(g3*yb)-64.0*sqr(g3*g3));
	    
	    dmHdSqdp = 64.0*sqr(g3*yb)*(2.0*M3-Ab)*t*(1-t)/sqr(16.0*PI*PI);
	    
	    dmHuSqdp = 64.0*sqr(g3*yt)*(2.0*M3-At)*t*(1-t)/sqr(16.0*PI*PI);
	    
	    dmQlSqdp = (t/(16.0*PI*PI))*(-(64.0/3.0)*g3*g3*M3+(1.0/(16.0*PI*PI))*(-(256.0/3.0)*sqr(g3*g3)*M3
										  +32.0*sqr(g3*g2)*(2.0*M3+M2)
										  +(32.0/45.0)*sqr(g3*g1)*(2.0*M3+M1)))
	      +sqr(t/(16.0*PI*PI))*(-(64.0/3.0)*sqr(g3*yt)*(2.0*M3-At)-(64.0/3.0)*sqr(g3*yb)*(2.0*M3-Ab)+192.0*sqr(g3*g3)*M3);
	    
	    dmUrSqdp = (t/(16.0*PI*PI))*(-(64.0/3.0)*g3*g3*M3+(1.0/(16.0*PI*PI))*(-(256.0/3.0)*sqr(g3*g3)*M3
										  +(512.0/45.0)*sqr(g1*g3)*(2.0*M3+M1)))
	      +sqr(t/(16.0*PI*PI))*(-(128.0/3.0)*sqr(g3*yt)*(2.0*M3-At)+192.0*sqr(g3*g3)*M3);

	    break;
	  }
	default:
	  {
	    ostringstream jj;
	    jj << "# WARNING: illegal index for gaugino soft mass: exiting." << endl;
	    throw jj.str();
	  }
	}
      
      derivs(1) = dmudp;
      derivs(2) = dBdp;
      derivs(3) = dmHdSqdp;
      derivs(4) = dmHuSqdp;
      derivs(5) = dmQlSqdp;
      derivs(6) = dmUrSqdp;
      derivs(7) = dAtdp;
    }
  return derivs;
}

DoubleVector doCalcSoftADerivs(SoftParsMssm r, DoubleVector pars, double mx, 
			       bool & hasError, DoubleVector const & auxPars, trilinears j, int m, int n)
{
  DoubleVector derivs(7);

  // Check that a diagonal trilinear element has been requested
  if (m != n)
    {
      ostringstream ii;
      ii << "# WARNING: can only calculate derivatives for diagonal soft trilinear elements: exiting." << endl;
      throw ii.str();
    }
  // Also check that the requested indices are valid
  else if ((m < 1 || m > 3) || (n < 1 || n > 3))
    {
 ostringstream jj;
      jj << "# WARNING: indices for soft trilinear matrices must be between 1 and 3: exiting." << endl;
      throw jj.str();
    }
  else
    {
      // Get M_{SUSY}
      double ms = r.displayMu();
      
      // Taylor approximation to the RGE solution is assumed to be about the values at MX.
      double t = log(ms/mx);
      
      double yt, yb, ytau, g1, g2, g3;
      // Check that all of the necessary couplings have been provided. If not, we can
      // still get them by running to MX, but we want to avoid this if possible.
      if (auxPars.displayEnd()-auxPars.displayStart()+1 != 6)
	{
	  cerr << "WARNING: incorrect number of couplings provided to doCalcMuDerivs: running to get values." << endl;
	  
	  // Copy the original object to do the running.
	  SoftParsMssm tempObject = r;
	  
	  tempObject.runto(mx);
	  
	  yt = tempObject.displayYukawaElement(YU, 3, 3);
	  yb = tempObject.displayYukawaElement(YD, 3, 3);
	  ytau = tempObject.displayYukawaElement(YE, 3, 3);
	  g1 = tempObject.displayGaugeCoupling(1);
	  g2 = tempObject.displayGaugeCoupling(2);
	  g3 = tempObject.displayGaugeCoupling(3);
	}
      else
	{
	  g1 = auxPars(1);
	  g2 = auxPars(2);
	  g3 = auxPars(3);
	  yt = auxPars(4);
	  yb = auxPars(5);
	  ytau = auxPars(6);
	}
      
      // The vector pars contains the values of the high scale parameters. Extract
      // the needed values.
      double M1 = pars(1);
      double M2 = pars(2);
      double M3 = pars(3);
      double At = pars(4);
      double Ab = pars(5);
      double Atau = pars(6);

      double dmudp, dBdp, dmHdSqdp, dmHuSqdp, dmQlSqdp, dmUrSqdp, dAtdp;
      
      // This is zero for all soft trilinears
      dmudp = 0.0;

      switch (j)
	{
	case UA:
	  {
	    if (m == 1 || m == 2)
	      {
		dBdp = 0.0;
		
		dAtdp = 0.0;
		
		dmHdSqdp = 0.0;
		
		dmHuSqdp = 0.0;
		
		dmQlSqdp = 0.0;
		
		dmUrSqdp = 0.0;
	      }
	    else
	      {
		dBdp = (t/(16.0*PI*PI))*(6.0*yt*yt+(1.0/(16.0*PI*PI))*(-36.0*sqr(yt*yt)-12.0*sqr(yt*yb)+32.0*sqr(g3*yt)
								       +1.6*sqr(g1*yt)))
		  +sqr(t/(16.0*PI*PI))*(72.0*sqr(yt*yt)+12.0*sqr(yt*yb)-32.0*sqr(g3*yt)-18.0*sqr(g2*yt)-5.2*sqr(g1*yt));

		dAtdp = 1.0+(t/(16.0*PI*PI))*(12.0*yt*yt+(1.0/(16.0*PI*PI))*(-88.0*sqr(yt*yt)-10.0*sqr(yt*yb)+32.0*sqr(g3*yt)
									     +12.0*sqr(g2*yt)+2.4*sqr(g1*yt)))
		  +sqr(t/(16.0*PI*PI))*(144.0*sqr(yt*yt)+14.0*sqr(yt*yb)-64.0*sqr(g3*yt)-36.0*sqr(g2*yt)-10.4*sqr(g1*yt));
		
		dmHdSqdp = 12.0*sqr(yt*yb)*(At+Ab)*t*(t-1)/sqr(16.0*PI*PI);
		
		dmHuSqdp = (t/(16.0*PI*PI))*(12.0*yt*yt*At+(1.0/(16.0*PI*PI))*(-144.0*sqr(yt*yt)*At-12.0*sqr(yt*yb)*(At+Ab)
									       +64.0*sqr(g3*yt)*(At-M3)+3.2*sqr(g1*yt)*(At-M1)))
		  +sqr(t/(16.0*PI*PI))*(288.0*sqr(yt*yt)*At+12.0*sqr(yt*yb)*(At+Ab)-64.0*sqr(g3*yt)*(At-M3)
					-36.0*sqr(g2*yt)*(At-M2)-10.4*sqr(g1*yt)*(At-M1));
		
		dmQlSqdp = (t/(16.0*PI*PI))*(4.0*yt*yt*At+(1.0/(16.0*PI*PI))*(-80.0*sqr(yt*yt)*At+3.2*sqr(g1*yt)*(At-M1)))
		  +sqr(t/(16.0*PI*PI))*(88.0*sqr(yt*yt)*At+8.0*sqr(yt*yb)*(At+Ab)-(64.0/3.0)*sqr(g3*yt)*(At-M3)
					-12.0*sqr(g2*yt)*(At-M2)-(52.0/15.0)*sqr(g1*yt)*(At-M1));
	    
		dmUrSqdp = (t/(16.0*PI*PI))*(8.0*yt*yt*At+(1.0/(16.0*PI*PI))*(-128.0*sqr(yt*yt)*At-8.0*sqr(yt*yb)*(At+Ab)
									      +24.0*sqr(g2*yt)*(At-M2)-1.6*sqr(g1*yt)*(At-M1)))
		  +sqr(t/(16.0*PI*PI))*(192.0*sqr(yt*yt)*At+8.0*sqr(yt*yb)*(At+Ab)-(128.0/3.0)*sqr(g3*yt)*(At-M3)
					-24.0*sqr(g2*yt)*(At-M2)-(104.0/15.0)*sqr(g1*yt)*(At-M1));
	      }
	    break;
	  }
	case DA:
	  {
	    if (m == 1 || m == 2)
	      {
		dBdp = 0.0;
		
		dAtdp = 0.0;
		
		dmHdSqdp = 0.0;
		
		dmHuSqdp = 0.0;
		
		dmQlSqdp = 0.0;
		
		dmUrSqdp = 0.0;
	      }
	    else
	      {
		dBdp = (t/(16.0*PI*PI))*(6.0*yb*yb+(1.0/(16.0*PI*PI))*(-36.0*sqr(yb*yb)-12.0*sqr(yt*yb)+32.0*sqr(g3*yb)
								       -0.8*sqr(g1*yb)))
		  +sqr(t/(16.0*PI*PI))*(72.0*sqr(yb*yb)+12.0*sqr(yt*yb)+12.0*sqr(ytau*yb)-32.0*sqr(g3*yb)-18.0*sqr(g2*yb)
					-2.8*sqr(g1*yb));

		dAtdp = (t/(16.0*PI*PI))*(2.0*yb*yb+(1.0/(16.0*PI*PI))*(-20.0*sqr(yb*yb)-10.0*sqr(yt*yb)-2.0*sqr(yb*ytau)
									+0.8*sqr(g1*yb)))
		  +sqr(t/(16.0*PI*PI))*(24.0*sqr(yb*yb)+14.0*sqr(yt*yb)+2.0*sqr(yb*ytau)-(32.0/3.0)*sqr(g3*yb)-6.0*sqr(g2*yb)
					-(14.0/15.0)*sqr(g1*yb));
		
		dmHdSqdp = (t/(16.0*PI*PI))*(12.0*yb*yb*Ab+(1.0/(16.0*PI*PI))*(-144.0*sqr(yb*yb)*Ab-12.0*sqr(yt*yb)*(At+Ab)
									       +64.0*sqr(g3*yb)*(Ab-M3)-1.6*sqr(g1*yb)*(Ab-M1)))
		  +sqr(t/(16.0*PI*PI))*(288.0*sqr(yb*yb)*Ab+12.0*sqr(yt*yb)*(At+Ab)+24.0*sqr(ytau*yb)*(Atau+Ab)
					-64.0*sqr(g3*yb)*(Ab-M3)-36.0*sqr(g2*yb)*(Ab-M2)-5.6*sqr(g1*yb)*(Ab-M1));
		
		dmHuSqdp = 12.0*sqr(yt*yb)*(At+Ab)*t*(t-1)/sqr(16.0*PI*PI);
		
		dmQlSqdp = (t/(16.0*PI*PI))*(4.0*yb*yb*Ab+(1.0/(16.0*PI*PI))*(-80.0*sqr(yb*yb)*Ab-4.0*sqr(ytau*yb)*(Ab+Atau)
									      -1.6*sqr(g1*yb)*(Ab-M1)))
		  +sqr(t/(16.0*PI*PI))*(88.0*sqr(yb*yb)*Ab+8.0*sqr(yt*yb)*(At+Ab)+4.0*sqr(yb*ytau)*(Ab+Atau)
					-(64.0/3.0)*sqr(g3*yb)*(Ab-M3)-12.0*sqr(g2*yb)*(Ab-M2)-(28.0/15.0)*sqr(g1*yb)*(Ab-M1));
	    
		dmUrSqdp = 8.0*sqr(yt*yb)*(At+Ab)*t*(t-1)/sqr(16.0*PI*PI);
	      }
	    break;
	  }
	case EA:
	  {
	    if (m == 1 || m == 2)
	      {
		dBdp = 0.0;
		
		dAtdp = 0.0;
		
		dmHdSqdp = 0.0;
		
		dmHuSqdp = 0.0;
		
		dmQlSqdp = 0.0;
		
		dmUrSqdp = 0.0;
	      }
	    else
	      {
		dBdp = (t/(16.0*PI*PI))*(2.0*ytau*ytau+(1.0/(16.0*PI*PI))*(-12.0*sqr(ytau*ytau)+2.4*sqr(g1*ytau)))
		  +sqr(t/(16.0*PI*PI))*(16.0*sqr(ytau*ytau)+12.0*sqr(ytau*yb)-6.0*sqr(g2*ytau)-3.6*sqr(g1*ytau));

		dAtdp = 2.0*sqr(yb*ytau)*t*(t-1)/sqr(16.0*PI*PI);
		
		dmHdSqdp = (t/(16.0*PI*PI))*(4.0*ytau*ytau*Atau+(1.0/(16.0*PI*PI))*(-48.0*sqr(ytau*ytau)*Atau
										    +4.8*sqr(g1*ytau)*(Atau-M1)))
		  +sqr(t/(16.0*PI*PI))*(24.0*sqr(ytau*yb)*(Atau+Ab)+64.0*sqr(ytau*ytau)*Atau-12.0*sqr(g2*ytau)*(Atau-M2)
					-7.2*sqr(g1*ytau)*(Atau-M1));
		
		dmHuSqdp = 0.0;
		
		dmQlSqdp = 4.0*sqr(ytau*yb)*(Ab+Atau)*t*(t-1)/sqr(16.0*PI*PI);
	    
		dmUrSqdp = 0.0;
	      }
	    break;
	  }
	default:
	  {
	    ostringstream kk;
	    kk << "# WARNING: illegal soft trilinear requested: exiting." << endl;
	    throw kk.str();
	  }
	}
      
      derivs(1) = dmudp;
      derivs(2) = dBdp;
      derivs(3) = dmHdSqdp;
      derivs(4) = dmHuSqdp;
      derivs(5) = dmQlSqdp;
      derivs(6) = dmUrSqdp;
      derivs(7) = dAtdp;
    }
  return derivs;

}

DoubleVector doCalcSoftMassSquaredDerivs(SoftParsMssm r, DoubleVector pars, double mx, 
					 bool & hasError, DoubleVector const & auxPars, softMasses j, int m, int n)
{
  DoubleVector derivs(7);

  // Has the user requested a diagonal soft mass squared (which are the only ones we can calculate for
  // at the moment)? If not, throw an error.
  if (m != n)
    {
      ostringstream ii;
      ii << "# WARNING: can only calculate derivatives for diagonal soft masses: exiting." << endl;
      throw ii.str();
    }
  // Additionally check if the indices m and n are valid (i.e. between 1 and 3)
  else if ((m < 1 || m > 3) || (n < 1 || n > 3))
    {
      ostringstream jj;
      jj << "# WARNING: indices for soft mass matrices must be between 1 and 3: exiting." << endl;
      throw jj.str();
    }
  else
    {
      // Get M_{SUSY}
      double ms = r.displayMu();
      
      // Taylor approximation to the RGE solution is assumed to be about the values at MX.
      double t = log(ms/mx);
            
      double yt, yb, ytau, g1, g2, g3;
      // Check that all of the necessary couplings have been provided. If not, we can
      // still get them by running to MX, but we want to avoid this if possible.
      if (auxPars.displayEnd()-auxPars.displayStart()+1 != 6)
	{
	  cerr << "WARNING: incorrect number of couplings provided to doCalcMuDerivs: running to get values." << endl;
	  
	  // Copy the original object to do the running.
	  SoftParsMssm tempObject = r;
	  
	  tempObject.runto(mx);
	  
	  yt = tempObject.displayYukawaElement(YU, 3, 3);
	  yb = tempObject.displayYukawaElement(YD, 3, 3);
	  ytau = tempObject.displayYukawaElement(YE, 3, 3);
	  g1 = tempObject.displayGaugeCoupling(1);
	  g2 = tempObject.displayGaugeCoupling(2);
	  g3 = tempObject.displayGaugeCoupling(3);
	}
      else
	{
	  g1 = auxPars(1);
	  g2 = auxPars(2);
	  g3 = auxPars(3);
	  yt = auxPars(4);
	  yb = auxPars(5);
	  ytau = auxPars(6);
	}
      
      // The vector pars contains the values of the high scale parameters. Extract
      // the needed values.
      
      double dmudp, dBdp, dmHdSqdp, dmHuSqdp, dmQlSqdp, dmUrSqdp, dAtdp;
      
      // Irrespective of which soft mass is requested, these three are zero.
      dBdp = 0.0;
      dAtdp = 0.0;
      dmudp = 0.0;
      
      // Determine which soft mass has ben requested and use the appropriate expressions
      // to calculate the derivatives.
      switch(j)
	{
	case mQl:
	  {
	    if (m == 1 || m == 2)
	      {
		dmHdSqdp = (t/(16.0*PI*PI))*(-0.6*g1*g1+(1.0/(16.0*PI*PI))*(-3.2*sqr(g3*g1)
								  +9.0*sqr(g2*g2)-1.8*sqr(g1*g2)+(2.0/25.0)*sqr(g1*g1)))
		  +sqr(t/(16.0*PI*PI))*(-(198.0/25.0)*sqr(g1*g1));
		
		dmHuSqdp = (t/(16.0*PI*PI))*(0.6*g1*g1+(1.0/(16.0*PI*PI))*(3.2*sqr(g1*g3)+1.8*sqr(g1*g2)+9.0*sqr(g2*g2)
									   +(4.0/25.0)*sqr(g1*g1)))
		  +sqr(t/(16.0*PI*PI))*((198.0/25.0)*sqr(g1*g1));
		
		dmQlSqdp = (t/(16.0*PI*PI))*(0.2*g1*g1+(1.0/(16.0*PI*PI))*((32.0/3.0)*sqr(g3*g3)
								      +(48.0/45.0)*sqr(g3*g1)+9.0*sqr(g2*g2)+0.6*sqr(g1*g2)
								      +(2.0/75.0)*sqr(g1*g1)))
		  +sqr(t/(16.0*PI*PI))*((66.0/25.0)*sqr(g1*g1));
		
		dmUrSqdp = (t/(16.0*PI*PI))*(-0.8*g1*g1+(1.0/(16.0*PI*PI))*((32.0/3.0)*sqr(g3*g3)
								  -(192.0/45.0)*sqr(g3*g1)-2.4*sqr(g1*g2)+(12.0/75.0)*sqr(g1*g1)))
		  +sqr(t/(16.0*PI*PI))*(-(264.0/25.0)*sqr(g1*g1));
      
	      }
	    else
	      {
		dmHdSqdp = (t/(16.0*PI*PI))*(6.0*yb*yb-0.6*g1*g1
					     +(1.0/(16.0*PI*PI))*(-36.0*sqr(yb*yb)-12.0*sqr(yt*yb)+32.0*sqr(g3*yb)
								  +1.2*sqr(g1*yt)+0.4*sqr(g1*yb)-3.2*sqr(g3*g1)
								  +9.0*sqr(g2*g2)-1.8*sqr(g1*g2)+(2.0/25.0)*sqr(g1*g1)))
		  +sqr(t/(16.0*PI*PI))*(72.0*sqr(yb*yb)+12.0*sqr(yt*yb)+12.0*sqr(ytau*yb)-32.0*sqr(g3*yb)
					-18.0*sqr(g2*yb)-2.8*sqr(g1*yb)-(198.0/25.0)*sqr(g1*g1));
		
		dmHuSqdp = (t/(16.0*PI*PI))*(6.0*yt*yt+0.6*g1*g1+(1.0/(16.0*PI*PI))*(-36.0*sqr(yt*yt)-12.0*sqr(yt*yb)
										     +32.0*sqr(g3*yt)+0.4*sqr(g1*yt)
										     -1.2*sqr(g1*yb)+3.2*sqr(g1*g3)
										     +1.8*sqr(g1*g2)+9.0*sqr(g2*g2)
										     +(4.0/25.0)*sqr(g1*g1)))
		  +sqr(t/(16.0*PI*PI))*(72.0*sqr(yt*yt)+12.0*sqr(yt*yb)-32.0*sqr(g3*yt)-18.0*sqr(g2*yt)-5.2*sqr(g1*yt)
					+(198.0/25.0)*sqr(g1*g1));
		
		dmQlSqdp = 1.0+(t/(16.0*PI*PI))*(2.0*yt*yt+2.0*yb*yb+0.2*g1*g1
						 +(1.0/(16.0*PI*PI))*(-20.0*sqr(yt*yt)-20.0*sqr(yb*yb)-2.0*sqr(ytau*yb)
								      +1.2*sqr(g1*yt)+0.4*sqr(g1*yb)+(32.0/3.0)*sqr(g3*g3)
								      +(48.0/45.0)*sqr(g3*g1)+9.0*sqr(g2*g2)+0.6*sqr(g1*g2)
								      +(2.0/75.0)*sqr(g1*g1)))
		  +sqr(t/(16.0*PI*PI))*(24.0*sqr(yt*yt)+24.0*sqr(yb*yb)+8.0*sqr(yt*yb)+2.0*sqr(yb*ytau)-(32.0/3.0)*sqr(g3*yt)
					-(32.0/3.0)*sqr(g3*yb)-6.0*sqr(g2*yt)-6.0*sqr(g2*yb)-(26.0/15.0)*sqr(g1*yt)
					-(14.0/15.0)*sqr(g1*yb)+(66.0/25.0)*sqr(g1*g1));
		
		dmUrSqdp = (t/(16.0*PI*PI))*(4.0*yt*yt-0.8*g1*g1
					     +(1.0/(16.0*PI*PI))*(-32.0*sqr(yt*yt)-8.0*sqr(yt*yb)+12.0*sqr(g2*yt)
								  +1.6*sqr(g1*yt)+1.6*sqr(g1*yb)+(32.0/3.0)*sqr(g3*g3)
								  -(192.0/45.0)*sqr(g3*g1)-2.4*sqr(g1*g2)+(12.0/75.0)*sqr(g1*g1)))
		  +sqr(t/(16.0*PI*PI))*(48.0*sqr(yt*yt)+8.0*sqr(yt*yb)-(64.0/3.0)*sqr(g3*yt)-12.0*sqr(g2*yt)-(52.0/15.0)*sqr(g1*yt)
					-(264.0/25.0)*sqr(g1*g1));
      
	      }
	    break;
	  }
	case mUr:
	  {
	    if (m == 1 || m == 2)
	      {
		dmHdSqdp = (t/(16.0*PI*PI))*(1.2*g1*g1+(1.0/(16.0*PI*PI))*(6.4*sqr(g1*g3)+(56.0/25.0)*sqr(g1*g1)))
		  +sqr(t/(16.0*PI*PI))*((396.0/25.0)*sqr(g1*g1));
		
		dmHuSqdp = (t/(16.0*PI*PI))*(-1.2*g1*g1+(1.0/(16.0*PI*PI))*(-6.4*sqr(g1*g3)-(8.0/25.0)*sqr(g1*g1)))
		  +sqr(t/(16.0*PI*PI))*(-(396.0/25.0)*sqr(g1*g1));
		
		dmQlSqdp = (t/(16.0*PI*PI))*(-0.4*g1*g1+(1.0/(16.0*PI*PI))*((16.0/3.0)*sqr(g3*g3)
								  -(96.0/45.0)*sqr(g3*g1)-(24.0/75.0)*sqr(g1*g1)))
		  +sqr(t/(16.0*PI*PI))*(-(132.0/25.0)*sqr(g1*g1));
		
		dmUrSqdp = (t/(16.0*PI*PI))*(1.6*g1*g1+(1.0/(16.0*PI*PI))*((16.0/3.0)*sqr(g3*g3)
								      +(384.0/45.0)*sqr(g3*g1)+(256.0/75.0)*sqr(g1*g1)))
		  +sqr(t/(16.0*PI*PI))*((528.0/25.0)*sqr(g1*g1));
      
	      }
	    else
	      {
		dmHdSqdp = (t/(16.0*PI*PI))*(1.2*g1*g1+(1.0/(16.0*PI*PI))*(-6.0*sqr(yt*yb)-4.8*sqr(g1*yt)
									   +6.4*sqr(g3*g1)+(56.0/25.0)*sqr(g1*g1)))
		  +sqr(t/(16.0*PI*PI))*(6.0*sqr(yt*yb)+(396.0/25.0)*sqr(g1*g1));
		
		dmHuSqdp = (t/(16.0*PI*PI))*(6.0*yt*yt-1.2*g1*g1+(1.0/(16.0*PI*PI))*(-36.0*sqr(yt*yt)-6.0*sqr(yt*yb)
										     +32.0*sqr(g3*yt)+6.4*sqr(g1*yt)
										     -6.4*sqr(g1*g3)-(8.0/25.0)*sqr(g1*g1)))
		  +sqr(t/(16.0*PI*PI))*(72.0*sqr(yt*yt)+6.0*sqr(yt*yb)-32.0*sqr(g3*yt)-18.0*sqr(g2*yt)-5.2*sqr(g1*yt)
					-(396.0/25.0)*sqr(g1*g1));
		
		dmQlSqdp = (t/(16.0*PI*PI))*(2.0*yt*yt-0.4*g1*g1
					     +(1.0/(16.0*PI*PI))*(-20.0*sqr(yt*yt)+3.2*sqr(g1*yt)+(16.0/3.0)*sqr(g3*g3)
								  -(96.0/45.0)*sqr(g3*g1)-(24.0/75.0)*sqr(g1*g1)))
		  +sqr(t/(16.0*PI*PI))*(24.0*sqr(yt*yt)+4.0*sqr(yt*yb)-(32.0/3.0)*sqr(g3*yt)-6.0*sqr(g2*yt)
					-(26.0/15.0)*sqr(g1*yt)-(132.0/25.0)*sqr(g1*g1));
		
		dmUrSqdp = 1.0+(t/(16.0*PI*PI))*(4.0*yt*yt+1.6*g1*g1
						 +(1.0/(16.0*PI*PI))*(-32.0*sqr(yt*yt)-4.0*sqr(yt*yb)+12.0*sqr(g2*yt)
								      -7.2*sqr(g1*yt)+(16.0/3.0)*sqr(g3*g3)
								      +(384.0/45.0)*sqr(g3*g1)+(256.0/75.0)*sqr(g1*g1)))
		  +sqr(t/(16.0*PI*PI))*(48.0*sqr(yt*yt)+4.0*sqr(yt*yb)-(64.0/3.0)*sqr(g3*yt)-12.0*sqr(g2*yt)
					-(52.0/15.0)*sqr(g1*yt)+(528.0/25)*sqr(g1*g1));
      
	      }
	    break;
	  }
	case mDr:
	  {
	    if (m == 1 || m == 2)
	      {
		dmHdSqdp = (t/(16.0*PI*PI))*(-0.6*g1*g1+(1.0/(16.0*PI*PI))*(-3.2*sqr(g1*g3)+(2.0/25.0)*sqr(g1*g1)))
		  +sqr(t/(16.0*PI*PI))*(-(198.0/25.0)*sqr(g1*g1));
		
		dmHuSqdp = (t/(16.0*PI*PI))*(0.6*g1*g1+(1.0/(16.0*PI*PI))*(3.2*sqr(g1*g3)+0.4*sqr(g1*g1)))
		  +sqr(t/(16.0*PI*PI))*((198.0/25.0)*sqr(g1*g1));
		
		dmQlSqdp = (t/(16.0*PI*PI))*(0.2*g1*g1+(1.0/(16.0*PI*PI))*((16.0/3.0)*sqr(g3*g3)
								  +(48.0/45.0)*sqr(g3*g1)+(6.0/75.0)*sqr(g1*g1)))
		  +sqr(t/(16.0*PI*PI))*((66.0/25.0)*sqr(g1*g1));
		
		dmUrSqdp = (t/(16.0*PI*PI))*(-0.8*g1*g1+(1.0/(16.0*PI*PI))*((16.0/3.0)*sqr(g3*g3)
									    -(192.0/45.0)*sqr(g3*g1)+(16.0/75.0)*sqr(g1*g1)))
		  +sqr(t/(16.0*PI*PI))*(-(264.0/25.0)*sqr(g1*g1));
      
	      }
	    else
	      {
		dmHdSqdp = (t/(16.0*PI*PI))*(6.0*yb*yb-0.6*g1*g1
					     +(1.0/(16.0*PI*PI))*(-36.0*sqr(yb*yb)-6.0*sqr(yt*yb)+32.0*sqr(g3*yb)
								  +1.6*sqr(g1*yb)-3.2*sqr(g1*g3)+(2.0/25.0)*sqr(g1*g1)))
		  +sqr(t/(16.0*PI*PI))*(72.0*sqr(yb*yb)+6.0*sqr(yt*yb)+12.0*sqr(ytau*yb)-32.0*sqr(g3*yb)-18.0*sqr(g2*yb)
					-2.8*sqr(g1*yb)-(198.0/25.0)*sqr(g1*g1));
		
		dmHuSqdp = (t/(16.0*PI*PI))*(0.6*g1*g1+(1.0/(16.0*PI*PI))*(-6.0*sqr(yt*yb)-2.4*sqr(g1*yb)+3.2*sqr(g1*g3)
									   +0.4*sqr(g1*g1)))
		  +sqr(t/(16.0*PI*PI))*(6.0*sqr(yt*yb)+(198.0/25.0)*sqr(g1*g1));
		
		dmQlSqdp = (t/(16.0*PI*PI))*(2.0*yb*yb+0.2*g1*g1
					     +(1.0/(16.0*PI*PI))*(-20.0*sqr(yb*yb)-2.0*sqr(ytau*yb)+(16.0/3.0)*sqr(g3*g3)
								  +(48.0/45.0)*sqr(g3*g1)+(6.0/75.0)*sqr(g1*g1)))
		  +sqr(t/(16.0*PI*PI))*(24.0*sqr(yb*yb)+4.0*sqr(yt*yb)+2.0*sqr(yb*ytau)-(32.0/3.0)*sqr(g3*yb)
					-6.0*sqr(g2*yb)-(14.0/15.0)*sqr(g1*yb)+(66.0/25.0)*sqr(g1*g1));
		
		dmUrSqdp = (t/(16.0*PI*PI))*(-0.8*g1*g1+(1.0/(16.0*PI*PI))*(-4.0*sqr(yt*yb)+3.2*sqr(g1*yb)+(16.0/3.0)*sqr(g3*g3)
									    -(192.0/45.0)*sqr(g3*g1)+(16.0/75.0)*sqr(g1*g1)))
		  +sqr(t/(16.0*PI*PI))*(4.0*sqr(yt*yb)-(264.0/25.0)*sqr(g1*g1));
      
	      }
	    break;
	  }
	case mLl:
	  {
	    if (m == 1 || m == 2)
	      {
		dmHdSqdp = (t/(16.0*PI*PI))*(0.6*g1*g1+(1.0/(16.0*PI*PI))*(3.0*sqr(g2*g2)+1.8*sqr(g1*g2)+(18.0/25.0)*sqr(g1*g1)))
		  +sqr(t/(16.0*PI*PI))*((198.0/25.0)*sqr(g1*g1));
		
		dmHuSqdp = (t/(16.0*PI*PI))*(-0.6*g1*g1+(1.0/(16.0*PI*PI))*(-1.8*sqr(g1*g2)+3.0*sqr(g2*g2)))
		  +sqr(t/(16.0*PI*PI))*(-(198.0/25.0)*sqr(g1*g1));
		
		dmQlSqdp = (t/(16.0*PI*PI))*(-0.2*g1*g1+(1.0/(16.0*PI*PI))*(3.0*sqr(g2*g2)-0.6*sqr(g1*g2)-(6.0/75.0)*sqr(g1*g1)))
		  +sqr(t/(16.0*PI*PI))*(-(66.0/25.0)*sqr(g1*g1));
		
		dmUrSqdp = (t/(16.0*PI*PI))*(0.8*g1*g1+(1.0/(16.0*PI*PI))*(2.4*sqr(g1*g2)+(84.0/75.0)*sqr(g1*g1)))
		  +sqr(t/(16.0*PI*PI))*((264.0/25.0)*sqr(g1*g1));
      
	      }
	    else
	      {
		dmHdSqdp = (t/(16.0*PI*PI))*(2.0*ytau*ytau+0.6*g1*g1
					     +(1.0/(16.0*PI*PI))*(-12.0*sqr(ytau*ytau)+1.2*sqr(g1*ytau)+3.0*sqr(g2*g2)
								  +1.8*sqr(g1*g2)+(18.0/25.0)*sqr(g1*g1)))
		  +sqr(t/(16.0*PI*PI))*(12.0*sqr(ytau*yb)+16.0*sqr(ytau*ytau)-6.0*sqr(g2*ytau)-3.6*sqr(g1*ytau)
					+(198.0/25.0)*sqr(g1*g1));
		
		dmHuSqdp = (t/(16.0*PI*PI))*(-0.6*g1*g1+(1.0/(16.0*PI*PI))*(1.2*sqr(g1*ytau)-1.8*sqr(g1*g2)+3.0*sqr(g2*g2)))
		  +sqr(t/(16.0*PI*PI))*(-(198.0/25.0)*sqr(g1*g1));
		
		dmQlSqdp = (t/(16.0*PI*PI))*(-0.2*g1*g1+(1.0/(16.0*PI*PI))*(-2.0*sqr(ytau*yb)+0.4*sqr(g1*ytau)
									    +3.0*sqr(g2*g2)-0.6*sqr(g1*g2)-(6.0/75.0)*sqr(g1*g1)))
		  +sqr(t/(16.0*PI*PI))*(2.0*sqr(yb*ytau)-(66.0/25.0)*sqr(g1*g1));
		
		dmUrSqdp = (t/(16.0*PI*PI))*(0.8*g1*g1+(1.0/(16.0*PI*PI))*(-1.6*sqr(g1*ytau)+2.4*sqr(g1*g2)
									   +(84.0/75.0)*sqr(g1*g1)))
		  +sqr(t/(16.0*PI*PI))*((264.0/25.0)*sqr(g1*g1));
      
	      }
	    break;
	  }
	case mEr:
	  {
	    if (m == 1 || m == 2)
	      {
		dmHdSqdp = (t/(16.0*PI*PI))*(-0.6*g1*g1+(1.0/(16.0*PI*PI))*(-(18.0/25.0)*sqr(g1*g1)))
		  +sqr(t/(16.0*PI*PI))*(-(198.0/25.0)*sqr(g1*g1));
		
		dmHuSqdp = (t/(16.0*PI*PI))*(0.6*g1*g1+(1.0/(16.0*PI*PI))*((54.0/25.0)*sqr(g1*g1)))
		  +sqr(t/(16.0*PI*PI))*((198.0/25.0)*sqr(g1*g1));
		
		dmQlSqdp = (t/(16.0*PI*PI))*(0.2*g1*g1+(1.0/(16.0*PI*PI))*((42.0/75.0)*sqr(g1*g1)))
		  +sqr(t/(16.0*PI*PI))*((66.0/25.0)*sqr(g1*g1));
		
		dmUrSqdp = (t/(16.0*PI*PI))*(-0.8*g1*g1+(1.0/(16.0*PI*PI))*(-(48.0/75.0)*sqr(g1*g1)))
		  +sqr(t/(16.0*PI*PI))*(-(264.0/25.0)*sqr(g1*g1));
      
	      }
	    else
	      {
		dmHdSqdp = (t/(16.0*PI*PI))*(2.0*ytau*ytau-0.6*g1*g1
					     +(1.0/(16.0*PI*PI))*(-12.0*sqr(ytau*ytau)+4.8*sqr(g1*ytau)-(18.0/25.0)*sqr(g1*g1)))
		  +sqr(t/(16.0*PI*PI))*(12.0*sqr(ytau*yb)+16.0*sqr(ytau*ytau)-6.0*sqr(g2*ytau)-3.6*sqr(g1*ytau)
					-(198.0/25.0)*sqr(g1*g1));
		
		dmHuSqdp = (t/(16.0*PI*PI))*(0.6*g1*g1+(1.0/(16.0*PI*PI))*(-2.4*sqr(g1*ytau)+(54.0/25.0)*sqr(g1*g1)))
		  +sqr(t/(16.0*PI*PI))*((198.0/25.0)*sqr(g1*g1));
		
		dmQlSqdp = (t/(16.0*PI*PI))*(0.2*g1*g1+(1.0/(16.0*PI*PI))*(-2.0*sqr(ytau*yb)-0.8*sqr(g1*ytau)
									   +(42.0/75.0)*sqr(g1*g1)))
		  +sqr(t/(16.0*PI*PI))*(2.0*sqr(yb*ytau)+(66.0/25.0)*sqr(g1*g1));
		
		dmUrSqdp = (t/(16.0*PI*PI))*(-0.8*g1*g1+(1.0/(16.0*PI*PI))*(3.2*sqr(g1*ytau)-(48.0/75.0)*sqr(g1*g1)))
		  +sqr(t/(16.0*PI*PI))*(-(264.0/25.0)*sqr(g1*g1));
      
	      }
	    break;
	  }
	default:
	  {
	    ostringstream kk;
	    kk << "# WARNING: illegal soft mass requested: exiting." << endl;
	    throw kk.str();
	  }
	} 
      
      derivs(1) = dmudp;
      derivs(2) = dBdp;
      derivs(3) = dmHdSqdp;
      derivs(4) = dmHuSqdp;
      derivs(5) = dmQlSqdp;
      derivs(6) = dmUrSqdp;
      derivs(7) = dAtdp;
      
    }
  return derivs;
}

// A function to calculate the vectors appearing on the RHS of the calculation of the
// derivatives d v_i/ dp_j. Takes as arguments a SoftParsMssm object, assumed to be evaluated
// at the scale M_SUSY, vectors of the parameters to calculate the derivatives wrt. and the vevs,
// an integer labelling which particular parameter to calculate the derivative for, and the scale
// at which that parameter is defined. This second version uses the approximate solutions
// to the MSSM RGEs to compute the tuning. Use only for pMSSM models.
DoubleVector doCalcRHSTuningVector_pMSSM_Approx(SoftParsMssm r, void (*ftBCatMX)(SoftParsMssm &, DoubleVector &), 
						DoubleVector pars, DoubleVector const & vevs, int i, double mx, 
						bool & hasError, DoubleVector const & auxPars)
{

  DoubleVector rhsVec(2);

  if (vevs.displayEnd()-vevs.displayStart()+1 < 2)
    {
      ostringstream kk;
      kk << "# WARNING: incorrect number of VEVs provided to doCalcRHSTuningVector_Approx: exiting." << endl;
      throw kk.str();
    }
  else
    {
      double v1 = vevs(vevs.displayStart());
      double v2 = vevs(vevs.displayEnd());
      double tb = v2/v1;
      
      double mu = r.displaySusyMu();
      double Bmu = r.displayM3Squared();
      // Calculate B
      double B = Bmu/mu;     
      
      // First calculate the elements of the matrix that appears on the RHS in general, being derivatives of the EWSB
      // conditions wrt the EW scale parameters. For the EWSB conditions used in our study, this matrix has the form:
      // 
      //     [ df1/dmu df1/dB df1/dm_Hd^2 df1/dm_Hu^2 df1/dm_Ql^2 df1/dm_uR^2 df1/dA_t ]
      //     [ df2/dmu df2/dB df2/dm_Hd^2 df2/dm_Hu^2 df2/dm_Ql^2 df2/dm_uR^2 df2/dA_t ]
      DoubleMatrix EWderivs(2,7);
      
      double df1dmu = 2*mu-B*tb;//2.0*mu*v1-B*v2;
      double df2dmu = 2*mu-B/tb;//2.0*mu*v2-B*v1;
      
      double df1dmH1Sq = 1.0;//v1;
      double df2dmH1Sq = 0.0;
      
      double df1dmH2Sq = 0.0;
      double df2dmH2Sq = 1.0;//v2;
      
      double df1dB = -mu*tb;//-mu*v2;
      double df2dB = -mu/tb;//-mu*v1;
      
      double df1dmQlSq = 0.0;
      double df2dmQlSq = 0.0;
      
      double df1dmUrSq = 0.0;
      double df2dmUrSq = 0.0;
      
      double df1dAt = 0.0;
      double df2dAt = 0.0;
      
      if (INCLUDE1LPTADPOLES)
	{
	  df1dmu = df1dmu+doCalcMSSMd2DeltaVdMudv1(r, tb);
	  df2dmu = df2dmu+doCalcMSSMd2DeltaVdMudv2(r, tb);
	  
	  df1dmQlSq = df1dmQlSq+doCalcMSSMd2DeltaVdmQlSqdv1(r, tb);
	  df2dmQlSq = df2dmQlSq+doCalcMSSMd2DeltaVdmQlSqdv2(r, tb);
	  
	  df1dmUrSq = df1dmUrSq+doCalcMSSMd2DeltaVdmUrSqdv1(r, tb);
	  df2dmUrSq = df2dmUrSq+doCalcMSSMd2DeltaVdmUrSqdv2(r, tb);
	  
	  df1dAt = df1dAt+doCalcMSSMd2DeltaVdAtdv1(r, tb);
	  df2dAt = df2dAt+doCalcMSSMd2DeltaVdAtdv2(r, tb);
	  
	}    

      //     [ df1/dmu df1/dB df1/dm_Hd^2 df1/dm_Hu^2 df1/dm_Ql^2 df1/dm_uR^2 df1/dA_t ]
      //     [ df2/dmu df2/dB df2/dm_Hd^2 df2/dm_Hu^2 df2/dm_Ql^2 df2/dm_uR^2 df2/dA_t ]      

      EWderivs(1,1) = df1dmu;
      EWderivs(2,1) = df2dmu;
      EWderivs(1,2) = df1dB;
      EWderivs(2,2) = df2dB;
      EWderivs(1,3) = df1dmH1Sq;
      EWderivs(2,3) = df2dmH1Sq;
      EWderivs(1,4) = df1dmH2Sq;
      EWderivs(2,4) = df2dmH2Sq;
      EWderivs(1,5) = df1dmQlSq;
      EWderivs(2,5) = df2dmQlSq;
      EWderivs(1,6) = df1dmUrSq;
      EWderivs(2,6) = df2dmUrSq;
      EWderivs(1,7) = df1dAt;
      EWderivs(2,7) = df2dAt;

      // Then calculate the vector that stores the derivatives of the EW scale parameters wrt the high scale parameter
      // of interest. In our study this has the form
      // [ dmu/dp dB/dp dm_Hd^2/dp dm_Hu^2/dp dm_Ql^2/dp dm_uR^2/dp dA_t/dp]^T.
      // For this version of the function, this is done numerically by running between M_SUSY and MX to estimate the
      // derivatives.

      // Get the scale M_SUSY - note that the model is assumed to be provided at this scale
      double ms = r.displayMu();      

      DoubleVector paramDerivs(7);

      if (pars.displayEnd()-pars.displayStart()+1 != NUMPMSSMPARS) // sugra case
	{
	  ostringstream ii;
	  ii << "doCalcRHSTuningVector_pMSSM_Approx called with incorrect number of parameters." << endl;
	  throw ii.str();
	}

      DoubleVector errVec(7);
      double err;

      if (fabs(ms-mx) < TOLERANCE)
	{
	  switch(i)
	    {
	    case 4: // At
	      {
		paramDerivs(1) = 0.0;
		paramDerivs(2) = 0.0;
		paramDerivs(3) = 0.0;
		paramDerivs(4) = 0.0;
		paramDerivs(5) = 0.0;
		paramDerivs(6) = 0.0;
		paramDerivs(7) = -1.0;
		break;
	      }
	    case 7: // mHdsq
	      {
		paramDerivs(1) = 0.0;
		paramDerivs(2) = 0.0;
		paramDerivs(3) = -1.0;
		paramDerivs(4) = 0.0;
		paramDerivs(5) = 0.0;
		paramDerivs(6) = 0.0;
		paramDerivs(7) = 0.0;
		break;
	      }
	    case 8: // mHusq
	      {
		paramDerivs(1) = 0.0;
		paramDerivs(2) = 0.0;
		paramDerivs(3) = 0.0;
		paramDerivs(4) = -1.0;
		paramDerivs(5) = 0.0;
		paramDerivs(6) = 0.0;
		paramDerivs(7) = 0.0;
		break;
	      }
	    case 9: // mu
	      {
		paramDerivs(1) = -1.0;
		paramDerivs(2) = 0.0;
		paramDerivs(3) = 0.0;
		paramDerivs(4) = 0.0;
		paramDerivs(5) = 0.0;
		paramDerivs(6) = 0.0;
		paramDerivs(7) = 0.0;
		break;
	      }
	    case 10: // B
	      {
		paramDerivs(1) = 0.0;
		paramDerivs(2) = -1.0;
		paramDerivs(3) = 0.0;
		paramDerivs(4) = 0.0;
		paramDerivs(5) = 0.0;
		paramDerivs(6) = 0.0;
		paramDerivs(7) = 0.0;
		break;
	      }
	    case 18: // mQlsq
	      {
		paramDerivs(1) = 0.0;
		paramDerivs(2) = 0.0;
		paramDerivs(3) = 0.0;
		paramDerivs(4) = 0.0;
		paramDerivs(5) = -1.0;
		paramDerivs(6) = 0.0;
		paramDerivs(7) = 0.0;
		break;
	      }
	    case 19: // muRsq
	      {
		paramDerivs(1) = 0.0;
		paramDerivs(2) = 0.0;
		paramDerivs(3) = 0.0;
		paramDerivs(4) = 0.0;
		paramDerivs(5) = 0.0;
		paramDerivs(6) = -1.0;
		paramDerivs(7) = 0.0;
		break;
	      }
	    default:
	      {
		paramDerivs(1) = 0.0;
		paramDerivs(2) = 0.0;
		paramDerivs(3) = 0.0;
		paramDerivs(4) = 0.0;
		paramDerivs(5) = 0.0;
		paramDerivs(6) = 0.0;
		paramDerivs(7) = 0.0;
		break;
	      }
	    }
	}
      else
	{
	  // Here we use the approximate expression for each of the derivatives
	  // of the low scale parameters w.r.t the high scale parameter of interest.
	  switch(i)
	    {
	    case 1:
	      {
		paramDerivs = doCalcGauginoDerivs(r, pars, mx, hasError, auxPars, 1);
		break;
	      }
	    case 2:
	      {
		paramDerivs = doCalcGauginoDerivs(r, pars, mx, hasError, auxPars, 2);
		break;
	      }
	    case 3:
	      {
		paramDerivs = doCalcGauginoDerivs(r, pars, mx, hasError, auxPars, 3);
		break;
	      }
	    case 4:
	      {
		paramDerivs = doCalcSoftADerivs(r, pars, mx, hasError, auxPars, UA, 3, 3);
		break;
	      }
	    case 5:
	      {
		paramDerivs = doCalcSoftADerivs(r, pars, mx, hasError, auxPars, DA, 3, 3);
		break;
	      }
	    case 6:
	      {
		paramDerivs = doCalcSoftADerivs(r, pars, mx, hasError, auxPars, EA, 3, 3);
		break;
	      }
	    case 7:
	      {
		paramDerivs = doCalcMh1SquaredDerivs(r, pars, mx, hasError, auxPars);
		break;
	      }
	    case 8:
	      {
		paramDerivs = doCalcMh2SquaredDerivs(r, pars, mx, hasError, auxPars);
		break;
	      }
	    case 9:
	      {
		paramDerivs = doCalcMuDerivs(r, pars, mx, hasError, auxPars);
		break;
	      }
	    case 10:
	      {
		paramDerivs = doCalcBDerivs(r, pars, mx, hasError, auxPars);
		break;
	      }
	      // Note factors of 2 for degenerate first and second generation mass parameters
	    case 11:
	      {
		paramDerivs = 2.0*doCalcSoftMassSquaredDerivs(r, pars, mx, hasError, auxPars, mLl, 1, 1);
		break;
	      }
	    case 12:
	      {
		paramDerivs = 2.0*doCalcSoftMassSquaredDerivs(r, pars, mx, hasError, auxPars, mEr, 1, 1);
		break;
	      }
	    case 13:
	      {
		paramDerivs = 2.0*doCalcSoftMassSquaredDerivs(r, pars, mx, hasError, auxPars, mQl, 1, 1);
		break;
	      }
	    case 14:
	      {
		paramDerivs = 2.0*doCalcSoftMassSquaredDerivs(r, pars, mx, hasError, auxPars, mUr, 1, 1);
		break;
	      }
	    case 15:
	      {
		paramDerivs = 2.0*doCalcSoftMassSquaredDerivs(r, pars, mx, hasError, auxPars, mDr, 1, 1);
		break;
	      }
	    case 16:
	      {
		paramDerivs = doCalcSoftMassSquaredDerivs(r, pars, mx, hasError, auxPars, mLl, 3, 3);
		break;
	      }
	    case 17:
	      {
		paramDerivs = doCalcSoftMassSquaredDerivs(r, pars, mx, hasError, auxPars, mEr, 3, 3);
		break;
	      }
	    case 18:
	      {
		paramDerivs = doCalcSoftMassSquaredDerivs(r, pars, mx, hasError, auxPars, mQl, 3, 3);
		break;
	      }
	    case 19:
	      {
		paramDerivs = doCalcSoftMassSquaredDerivs(r, pars, mx, hasError, auxPars, mUr, 3, 3);
		break;
	      }
	    case 20:
	      {
		paramDerivs = doCalcSoftMassSquaredDerivs(r, pars, mx, hasError, auxPars, mDr, 3, 3);
		break;
	      }
	    default:
	      {
		cerr << "WARNING: unrecognised pMSSM parameter " << i << " requested: ignoring it." << endl;
		hasError = true;
		break;
	      }
	    }
	  for (int j = 1; j <= 7; j++)
	    {
	      paramDerivs(j) = -paramDerivs(j);
	    }
	}
      
   
      // Multiply the two together and return the result (note the additional minus sign!).
      rhsVec = EWderivs*paramDerivs; 
    }

  return rhsVec;
}


void pMSSMftBCs(SoftParsMssm & m, DoubleVector & tuningPars)
{
  if (PRINTOUT > 1)
    {
      cout << "Applying pMSSM fine tuning boundary condition" << endl;
      cout << "Model to apply to is:" << endl;
      cout << m;
      cout << "Setting the following parameters:" << endl;
      for (int i = 1; i <= 3; i++)
	{
	  cout << "M_" << i << " = " << tuningPars(i) << endl;
	}
      cout << "y_t = " << m.displayYukawaElement(YU, 3, 3) << endl;
      cout << "A_t = " << tuningPars(4) << endl;
      cout << "y_b = " << m.displayYukawaElement(YD, 3, 3) << endl;
      cout << "A_b = " << tuningPars(5) << endl;
      cout << "y_tau = " << m.displayYukawaElement(YE, 3, 3) << endl;
      cout << "A_tau = " << tuningPars(6) << endl;
      cout << "m_Hd^2 = " << tuningPars(7) << endl;
      cout << "m_Hu^2 = " << tuningPars(8) << endl;
      cout << "mu = " << tuningPars(9) << endl;
      cout << "B = " << tuningPars(10) << endl;
      cout << "m_Ll12^2 = " << tuningPars(11) << endl;
      cout << "m_Er12^2 = " << tuningPars(12) << endl;
      cout << "m_Ql12^2 = " << tuningPars(13) << endl;
      cout << "m_Ur12^2 = " << tuningPars(14) << endl;
      cout << "m_Dr12^2 = " << tuningPars(15) << endl;
      cout << "m_Ll3^2 = " << tuningPars(16) << endl;
      cout << "m_Er3^2 = " << tuningPars(17) << endl;
      cout << "m_Ql3^2 = " << tuningPars(18) << endl;
      cout << "m_Ur3^2 = " << tuningPars(19) << endl;
      cout << "m_Dr3^2 = " << tuningPars(20) << endl;
      cout << "Applying..." << endl;
    }

  // In the pMSSM, first and second generation Yukawas
  // and trilinears are zero.
  m.setYukawaElement(YU, 1, 1, 0.0);
  m.setYukawaElement(YU, 2, 2, 0.0);
  m.setYukawaElement(YD, 1, 1, 0.0);
  m.setYukawaElement(YD, 2, 2, 0.0);
  m.setYukawaElement(YE, 1, 1, 0.0);
  m.setYukawaElement(YE, 2, 2, 0.0);

  m.setTrilinearElement(UA, 1, 1, 0.0);
  m.setTrilinearElement(UA, 2, 2, 0.0);
  m.setTrilinearElement(DA, 1, 1, 0.0);
  m.setTrilinearElement(DA, 2, 2, 0.0);
  m.setTrilinearElement(EA, 1, 1, 0.0);
  m.setTrilinearElement(EA, 2, 2, 0.0);

  // Also make sure off-diagonal elements vanish.
  m.setYukawaElement(YU, 1, 2, 0.0);
  m.setYukawaElement(YU, 1, 3, 0.0);
  m.setYukawaElement(YU, 2, 1, 0.0);
  m.setYukawaElement(YU, 2, 3, 0.0);
  m.setYukawaElement(YU, 3, 1, 0.0);
  m.setYukawaElement(YU, 3, 2, 0.0);

  m.setYukawaElement(YD, 1, 2, 0.0);
  m.setYukawaElement(YD, 1, 3, 0.0);
  m.setYukawaElement(YD, 2, 1, 0.0);
  m.setYukawaElement(YD, 2, 3, 0.0);
  m.setYukawaElement(YD, 3, 1, 0.0);
  m.setYukawaElement(YD, 3, 2, 0.0);

  m.setYukawaElement(YE, 1, 2, 0.0);
  m.setYukawaElement(YE, 1, 3, 0.0);
  m.setYukawaElement(YE, 2, 1, 0.0);
  m.setYukawaElement(YE, 2, 3, 0.0);
  m.setYukawaElement(YE, 3, 1, 0.0);
  m.setYukawaElement(YE, 3, 2, 0.0);

  m.setTrilinearElement(UA, 1, 2, 0.0);
  m.setTrilinearElement(UA, 1, 3, 0.0);
  m.setTrilinearElement(UA, 2, 1, 0.0);
  m.setTrilinearElement(UA, 2, 3, 0.0);
  m.setTrilinearElement(UA, 3, 1, 0.0);
  m.setTrilinearElement(UA, 3, 2, 0.0);

  m.setTrilinearElement(DA, 1, 2, 0.0);
  m.setTrilinearElement(DA, 1, 3, 0.0);
  m.setTrilinearElement(DA, 2, 1, 0.0);
  m.setTrilinearElement(DA, 2, 3, 0.0);
  m.setTrilinearElement(DA, 3, 1, 0.0);
  m.setTrilinearElement(DA, 3, 2, 0.0);

  m.setTrilinearElement(EA, 1, 2, 0.0);
  m.setTrilinearElement(EA, 1, 3, 0.0);
  m.setTrilinearElement(EA, 2, 1, 0.0);
  m.setTrilinearElement(EA, 2, 3, 0.0);
  m.setTrilinearElement(EA, 3, 1, 0.0);
  m.setTrilinearElement(EA, 3, 2, 0.0);

  m.setSoftMassElement(mLl, 1, 2, 0.0);
  m.setSoftMassElement(mLl, 1, 3, 0.0);
  m.setSoftMassElement(mLl, 2, 1, 0.0);
  m.setSoftMassElement(mLl, 2, 3, 0.0);
  m.setSoftMassElement(mLl, 3, 1, 0.0);
  m.setSoftMassElement(mLl, 3, 2, 0.0);

  m.setSoftMassElement(mEr, 1, 2, 0.0);
  m.setSoftMassElement(mEr, 1, 3, 0.0);
  m.setSoftMassElement(mEr, 2, 1, 0.0);
  m.setSoftMassElement(mEr, 2, 3, 0.0);
  m.setSoftMassElement(mEr, 3, 1, 0.0);
  m.setSoftMassElement(mEr, 3, 2, 0.0);

  m.setSoftMassElement(mQl, 1, 2, 0.0);
  m.setSoftMassElement(mQl, 1, 3, 0.0);
  m.setSoftMassElement(mQl, 2, 1, 0.0);
  m.setSoftMassElement(mQl, 2, 3, 0.0);
  m.setSoftMassElement(mQl, 3, 1, 0.0);
  m.setSoftMassElement(mQl, 3, 2, 0.0);

  m.setSoftMassElement(mUr, 1, 2, 0.0);
  m.setSoftMassElement(mUr, 1, 3, 0.0);
  m.setSoftMassElement(mUr, 2, 1, 0.0);
  m.setSoftMassElement(mUr, 2, 3, 0.0);
  m.setSoftMassElement(mUr, 3, 1, 0.0);
  m.setSoftMassElement(mUr, 3, 2, 0.0);

  m.setSoftMassElement(mDr, 1, 2, 0.0);
  m.setSoftMassElement(mDr, 1, 3, 0.0);
  m.setSoftMassElement(mDr, 2, 1, 0.0);
  m.setSoftMassElement(mDr, 2, 3, 0.0);
  m.setSoftMassElement(mDr, 3, 1, 0.0);
  m.setSoftMassElement(mDr, 3, 2, 0.0);

  // Set the remaining parameters.
  for (int i = 1; i <= 3; i++)
    {
      m.setGauginoMass(i, tuningPars(i));
    }

  m.setTrilinearElement(UA, 3, 3, m.displayYukawaElement(YU, 3, 3) * 
			tuningPars(4));
  m.setTrilinearElement(DA, 3, 3, m.displayYukawaElement(YD, 3, 3) * 
			tuningPars(5));
  m.setTrilinearElement(EA, 3, 3, m.displayYukawaElement(YE, 3, 3) * 
			tuningPars(6));

  m.setMh1Squared(tuningPars(7));
  m.setMh2Squared(tuningPars(8));
  m.setSusyMu(tuningPars(9));
  m.setM3Squared(tuningPars(10)*tuningPars(9));

  m.setSoftMassElement(mLl, 1, 1, tuningPars(11));
  m.setSoftMassElement(mLl, 2, 2, tuningPars(11));
  m.setSoftMassElement(mEr, 1, 1, tuningPars(12));
  m.setSoftMassElement(mEr, 2, 2, tuningPars(12));
  m.setSoftMassElement(mQl, 1, 1, tuningPars(13));
  m.setSoftMassElement(mQl, 2, 2, tuningPars(13));
  m.setSoftMassElement(mUr, 1, 1, tuningPars(14));
  m.setSoftMassElement(mUr, 2, 2, tuningPars(14));
  m.setSoftMassElement(mDr, 1, 1, tuningPars(15));
  m.setSoftMassElement(mDr, 2, 2, tuningPars(15));
  m.setSoftMassElement(mLl, 3, 3, tuningPars(16));
  m.setSoftMassElement(mEr, 3, 3, tuningPars(17));
  m.setSoftMassElement(mQl, 3, 3, tuningPars(18));
  m.setSoftMassElement(mUr, 3, 3, tuningPars(19));
  m.setSoftMassElement(mDr, 3, 3, tuningPars(20));

}

void mSUGRAftBCs(SoftParsMssm & m, DoubleVector & tuningPars)
{
  double m0sq = tuningPars(1);
  double m12 = tuningPars(2);
  double susyMu = tuningPars(3);
  double B0 = tuningPars(4);
  double A0 = tuningPars(5);

  m.universalGauginos(m12);
  m.universalTrilinears(A0);
  m.universalScalars(sqrt(m0sq));

  m.setSusyMu(susyMu);
  m.setM3Squared(susyMu*B0);

}

// This single function call will calculate the fine tuning in a pMSSM model. The given object
// is assumed to be provided at MX, and the parameters to calculate the fine tuning with
// are assumed to be already set in that object. The value ms is an initial guess for the scale
// M_{SUSY}, which will be recalculated when the EWSB conditions are imposed.
DoubleVector doCalcpMSSMFineTuning(SoftParsMssm r, double ms, bool & ewsbProblem, bool & hasError, 
				   bool useApproxSolns, double tol)
{
  // Get the current scale
  double mx = r.displayMu();

  bool MxEqualsMs = false;
  if (fabs(mx-ms) < TOLERANCE)
    {
      MxEqualsMs = true;
    }

  // In the pMSSM as defined by Rizzo et al the first and second
  // generation Yukawa and trilinear couplings are zero. Make sure
  // this is the case.
  r.setYukawaElement(YU, 1, 1, 0.0);
  r.setYukawaElement(YU, 2, 2, 0.0);
  r.setYukawaElement(YD, 1, 1, 0.0);
  r.setYukawaElement(YD, 2, 2, 0.0);
  r.setYukawaElement(YE, 1, 1, 0.0);
  r.setYukawaElement(YE, 2, 2, 0.0);

  r.setTrilinearElement(UA, 1, 1, 0.0);
  r.setTrilinearElement(UA, 2, 2, 0.0);
  r.setTrilinearElement(DA, 1, 1, 0.0);
  r.setTrilinearElement(DA, 2, 2, 0.0);
  r.setTrilinearElement(EA, 1, 1, 0.0);
  r.setTrilinearElement(EA, 2, 2, 0.0);

  // Also make sure off-diagonal elements vanish
  r.setYukawaElement(YU, 1, 2, 0.0);
  r.setYukawaElement(YU, 1, 3, 0.0);
  r.setYukawaElement(YU, 2, 1, 0.0);
  r.setYukawaElement(YU, 2, 3, 0.0);
  r.setYukawaElement(YU, 3, 1, 0.0);
  r.setYukawaElement(YU, 3, 2, 0.0);

  r.setYukawaElement(YD, 1, 2, 0.0);
  r.setYukawaElement(YD, 1, 3, 0.0);
  r.setYukawaElement(YD, 2, 1, 0.0);
  r.setYukawaElement(YD, 2, 3, 0.0);
  r.setYukawaElement(YD, 3, 1, 0.0);
  r.setYukawaElement(YD, 3, 2, 0.0);

  r.setYukawaElement(YE, 1, 2, 0.0);
  r.setYukawaElement(YE, 1, 3, 0.0);
  r.setYukawaElement(YE, 2, 1, 0.0);
  r.setYukawaElement(YE, 2, 3, 0.0);
  r.setYukawaElement(YE, 3, 1, 0.0);
  r.setYukawaElement(YE, 3, 2, 0.0);

  r.setTrilinearElement(UA, 1, 2, 0.0);
  r.setTrilinearElement(UA, 1, 3, 0.0);
  r.setTrilinearElement(UA, 2, 1, 0.0);
  r.setTrilinearElement(UA, 2, 3, 0.0);
  r.setTrilinearElement(UA, 3, 1, 0.0);
  r.setTrilinearElement(UA, 3, 2, 0.0);

  r.setTrilinearElement(DA, 1, 2, 0.0);
  r.setTrilinearElement(DA, 1, 3, 0.0);
  r.setTrilinearElement(DA, 2, 1, 0.0);
  r.setTrilinearElement(DA, 2, 3, 0.0);
  r.setTrilinearElement(DA, 3, 1, 0.0);
  r.setTrilinearElement(DA, 3, 2, 0.0);

  r.setTrilinearElement(EA, 1, 2, 0.0);
  r.setTrilinearElement(EA, 1, 3, 0.0);
  r.setTrilinearElement(EA, 2, 1, 0.0);
  r.setTrilinearElement(EA, 2, 3, 0.0);
  r.setTrilinearElement(EA, 3, 1, 0.0);
  r.setTrilinearElement(EA, 3, 2, 0.0);

  r.setSoftMassElement(mLl, 1, 2, 0.0);
  r.setSoftMassElement(mLl, 1, 3, 0.0);
  r.setSoftMassElement(mLl, 2, 1, 0.0);
  r.setSoftMassElement(mLl, 2, 3, 0.0);
  r.setSoftMassElement(mLl, 3, 1, 0.0);
  r.setSoftMassElement(mLl, 3, 2, 0.0);

  r.setSoftMassElement(mEr, 1, 2, 0.0);
  r.setSoftMassElement(mEr, 1, 3, 0.0);
  r.setSoftMassElement(mEr, 2, 1, 0.0);
  r.setSoftMassElement(mEr, 2, 3, 0.0);
  r.setSoftMassElement(mEr, 3, 1, 0.0);
  r.setSoftMassElement(mEr, 3, 2, 0.0);

  r.setSoftMassElement(mQl, 1, 2, 0.0);
  r.setSoftMassElement(mQl, 1, 3, 0.0);
  r.setSoftMassElement(mQl, 2, 1, 0.0);
  r.setSoftMassElement(mQl, 2, 3, 0.0);
  r.setSoftMassElement(mQl, 3, 1, 0.0);
  r.setSoftMassElement(mQl, 3, 2, 0.0);

  r.setSoftMassElement(mUr, 1, 2, 0.0);
  r.setSoftMassElement(mUr, 1, 3, 0.0);
  r.setSoftMassElement(mUr, 2, 1, 0.0);
  r.setSoftMassElement(mUr, 2, 3, 0.0);
  r.setSoftMassElement(mUr, 3, 1, 0.0);
  r.setSoftMassElement(mUr, 3, 2, 0.0);

  r.setSoftMassElement(mDr, 1, 2, 0.0);
  r.setSoftMassElement(mDr, 1, 3, 0.0);
  r.setSoftMassElement(mDr, 2, 1, 0.0);
  r.setSoftMassElement(mDr, 2, 3, 0.0);
  r.setSoftMassElement(mDr, 3, 1, 0.0);
  r.setSoftMassElement(mDr, 3, 2, 0.0);

  // Check that the first and second generation soft masses are degenerate.
  // If not, take the first generation mass as the common value.
  if (fabs(r.displaySoftMassSquared(mLl, 1, 1)-r.displaySoftMassSquared(mLl, 2, 2)) > TOLERANCE)
    {
      r.setSoftMassElement(mLl, 2, 2, r.displaySoftMassSquared(mLl, 1, 1));
    }
  if (fabs(r.displaySoftMassSquared(mEr, 1, 1)-r.displaySoftMassSquared(mEr, 2, 2)) > TOLERANCE)
    {
      r.setSoftMassElement(mEr, 2, 2, r.displaySoftMassSquared(mEr, 1, 1));
    }
  if (fabs(r.displaySoftMassSquared(mQl, 1, 1)-r.displaySoftMassSquared(mQl, 2, 2)) > TOLERANCE)
    {
      r.setSoftMassElement(mQl, 2, 2, r.displaySoftMassSquared(mQl, 1, 1));
    }
  if (fabs(r.displaySoftMassSquared(mUr, 1, 1)-r.displaySoftMassSquared(mUr, 2, 2)) > TOLERANCE)
    {
      r.setSoftMassElement(mUr, 2, 2, r.displaySoftMassSquared(mUr, 1, 1));
    }
  if (fabs(r.displaySoftMassSquared(mDr, 1, 1)-r.displaySoftMassSquared(mDr, 2, 2)) > TOLERANCE)
    {
      r.setSoftMassElement(mDr, 2, 2, r.displaySoftMassSquared(mDr, 1, 1));
    }

  // Check if the EWSB conditions are satisfied. If not, 
  // vary the soft masses mHd^2 and mHu^2 so that they are.
  SoftParsMssm s = r;
  s.runto(ms);
  ewsbProblem = false;
  DoubleVector updatedMasses(3);
  updatedMasses(1) = r.displayMh1Squared();
  updatedMasses(2) = r.displayMh2Squared();
  updatedMasses(3) = ms;

  if (fabs(MSSM_EWSBCondition1(s)) > tol || fabs(MSSM_EWSBCondition2(s)) > tol)
    {
      if (PRINTOUT > 1) cout << "Reapplying EWSB conditions for fine tuning calculation..." << endl;

      ewsbProblem = MSSM_ImplementEWSBConstraints_SoftMasses(r, mx, ms, MxEqualsMs, updatedMasses, tol);

      r.setMh1Squared(updatedMasses(1));
      r.setMh2Squared(updatedMasses(2));

    }

  // Get parameter values.
  DoubleVector tuningPars(NUMPMSSMPARS);
 
  tuningPars(1) = r.displayGaugino(1);
  tuningPars(2) = r.displayGaugino(2);
  tuningPars(3) = r.displayGaugino(3);
  tuningPars(4) = r.displaySoftA(UA, 3, 3);
  tuningPars(5) = r.displaySoftA(DA, 3, 3);
  tuningPars(6) = r.displaySoftA(EA, 3, 3);
  tuningPars(7) = r.displayMh1Squared();
  tuningPars(8) = r.displayMh2Squared();
  tuningPars(9) = r.displaySusyMu();
  tuningPars(10) = r.displayM3Squared()/r.displaySusyMu();
  tuningPars(11) = r.displaySoftMassSquared(mLl, 1, 1);
  tuningPars(12) = r.displaySoftMassSquared(mEr, 1, 1);
  tuningPars(13) = r.displaySoftMassSquared(mQl, 1, 1);
  tuningPars(14) = r.displaySoftMassSquared(mUr, 1, 1);
  tuningPars(15) = r.displaySoftMassSquared(mDr, 1, 1);
  tuningPars(16) = r.displaySoftMassSquared(mLl, 3, 3);
  tuningPars(17) = r.displaySoftMassSquared(mEr, 3, 3);
  tuningPars(18) = r.displaySoftMassSquared(mQl, 3, 3);
  tuningPars(19) = r.displaySoftMassSquared(mUr, 3, 3);
  tuningPars(20) = r.displaySoftMassSquared(mDr, 3, 3);

  if (PRINTOUT > 1)
    {
      cout << "Using the following high scale values: " << endl;
   for (int i = 1; i <= 3; i++)
	{
	  cout << "M_" << i << " = " << tuningPars(i) << endl;
	}
      cout << "y_t = " << r.displayYukawaElement(YU, 3, 3) << endl;
      cout << "A_t = " << tuningPars(4) << endl;
      cout << "y_b = " << r.displayYukawaElement(YD, 3, 3) << endl;
      cout << "A_b = " << tuningPars(5) << endl;
      cout << "y_tau = " << r.displayYukawaElement(YE, 3, 3) << endl;
      cout << "A_tau = " << tuningPars(6) << endl;
      cout << "m_Hd^2 = " << tuningPars(7) << endl;
      cout << "m_Hu^2 = " << tuningPars(8) << endl;
      cout << "mu = " << tuningPars(9) << endl;
      cout << "B = " << tuningPars(10) << endl;
      cout << "m_Ll12^2 = " << tuningPars(11) << endl;
      cout << "m_Er12^2 = " << tuningPars(12) << endl;
      cout << "m_Ql12^2 = " << tuningPars(13) << endl;
      cout << "m_Ur12^2 = " << tuningPars(14) << endl;
      cout << "m_Dr12^2 = " << tuningPars(15) << endl;
      cout << "m_Ll3^2 = " << tuningPars(16) << endl;
      cout << "m_Er3^2 = " << tuningPars(17) << endl;
      cout << "m_Ql3^2 = " << tuningPars(18) << endl;
      cout << "m_Ur3^2 = " << tuningPars(19) << endl;
      cout << "m_Dr3^2 = " << tuningPars(20) << endl;
    }

  DoubleVector highScaleCouplings(6);
  highScaleCouplings(1) = r.displayGaugeCoupling(1);
  highScaleCouplings(2) = r.displayGaugeCoupling(2);
  highScaleCouplings(3) = r.displayGaugeCoupling(3);
  highScaleCouplings(4) = r.displayYukawaElement(YU, 3, 3);
  highScaleCouplings(5) = r.displayYukawaElement(YD, 3, 3);
  highScaleCouplings(6) = r.displayYukawaElement(YE, 3, 3);

  r.runto(updatedMasses(3));

  if (ewsbProblem)
    {
      cerr << "WARNING: failed to implement EWSB constraints: f1 = " <<  MSSM_EWSBCondition1(r) 
	   << ", f2 = " << MSSM_EWSBCondition2(r) << endl;
    }

  DoubleVector vevs(2);
  vevs(1) = r.displayHvev()/sqrt(1.0+sqr(r.displayTanb()));
  vevs(2) = vevs(1)*r.displayTanb();

  DoubleVector fineTunings(NUMPMSSMPARS);

  int sing = 0;

  if (useApproxSolns)
    {
      fineTunings = doCalcFineTuning<SoftParsMssm>(r, pMSSMftBCs, tuningPars, vevs, highScaleCouplings, mx, sing,
							 doCalcLHSTuningMatrix, 
							 doCalcRHSTuningVector_pMSSM_Approx, doCalcdLogMzSqdLogParam);
    }
  else
    {
      fineTunings = doCalcFineTuning<SoftParsMssm>(r, pMSSMftBCs, tuningPars, vevs, mx, sing, 
						   doCalcLHSTuningMatrix, doCalcRHSTuningVector, doCalcdLogMzSqdLogParam);
    }


  if (sing != 0)
    {
      hasError = true;
    }

  return fineTunings;

}


DoubleVector doCalcSUGRAFineTuning(SoftParsMssm r, double ms, bool & hasError, double tol)
{
  // Get the current scale
  double mx = r.displayMu();

  // Check if the EWSB conditions are satisfied. If not, 
  // vary the soft masses mHd^2 and mHu^2 so that they are.
  SoftParsMssm s = r;
  s.runto(ms);
  bool ewsbProblem = false;
  DoubleVector updatedMasses(3);
  updatedMasses(1) = r.displaySusyMu();
  updatedMasses(2) = r.displayM3Squared()/r.displaySusyMu();
  updatedMasses(3) = ms;

  if (fabs(MSSM_EWSBCondition1(s)) > tol || fabs(MSSM_EWSBCondition2(s)) > tol)
    {
      ewsbProblem = MSSM_ImplementEWSBConstraints_SUGRA(r, mx, ms, updatedMasses, tol);

      r.setSusyMu(updatedMasses(1));
      r.setM3Squared(updatedMasses(2)*updatedMasses(1));

    }

  DoubleVector tuningPars(5);

  tuningPars(1) = r.displayMh1Squared();
  tuningPars(2) = r.displayGaugino(1);
  tuningPars(3) = r.displaySusyMu();
  tuningPars(4) = r.displayM3Squared()/r.displaySusyMu();
  tuningPars(5) = r.displaySoftA(UA, 3, 3);

  r.runto(updatedMasses(3));

  if (ewsbProblem)
    {
      cerr << "WARNING: failed to implement EWSB constraints: f1 = " <<  MSSM_EWSBCondition1(r) 
	   << ", f2 = " << MSSM_EWSBCondition2(r) << endl;
    }

  DoubleVector vevs(2);
  vevs(1) = r.displayHvev()/sqrt(1.0+sqr(r.displayTanb()));
  vevs(2) = vevs(1)*r.displayTanb();

  DoubleVector fineTunings(5);

  int sing = 0;

  fineTunings = doCalcFineTuning<SoftParsMssm>(r, mSUGRAftBCs, tuningPars, vevs, mx, sing, 
					       doCalcLHSTuningMatrix, doCalcRHSTuningVector, doCalcdLogMzSqdLogParam);

  if (ewsbProblem || sing != 0)
    {
      hasError = true;
    }

  return fineTunings;
}
