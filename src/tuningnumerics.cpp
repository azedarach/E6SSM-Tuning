/*
  tuningnumerics provides some useful numerical methods for calculating
  the fine tuning in a given model. Essentially just easier to use
  wrapper functions for some of the numerical utilities in SOFTSUSY.
 */

#include "tuningnumerics.h"

// A method to generate the matrices Q and R in the QR decomposition
// A = QR of the n x n square matrix A. On output the original input matrix A is
// replaced by the upper triangular matrix R (the original matrix
// can always be recovered as A = QR if necessary). Not exactly optimised,
// but should be sufficient for now when only working with small matrices.
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

// Helper function to find index of maximum element in a vector.
int findMaxIndex(DoubleVector const & vec)
{
  int maxIndex = vec.displayStart();
  double maxVal = -1e6;

  for (int i = vec.displayStart(); i <= vec.displayEnd(); i++)
    {
      if (vec.display(i) >= maxVal) 
	{
	  maxIndex = i;
	  maxVal = vec.display(i);
	}
    }

  return maxIndex;

}

// bool lnsrch(const DoubleVector & xold, const DoubleVector & auxData, double fold, const DoubleVector & g, 
// 	    DoubleVector & p, 
// 	    DoubleVector & x, double & f, double stpmax, 
// 	    void (*vecfunc)(const DoubleVector &, const DoubleVector &, DoubleVector &), 
// 	    DoubleVector & fvec) {
//   double ALF = TOLERANCE;
//   double TOLX = TOLERANCE * 1.0e-3;
  
//   int i;
//   double a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2,slope,sum,temp,
//     test,tmplam;
  
//   bool err = false;
//   sum = sqrt(p.dot(p));
//   if (sum > stpmax) p = (stpmax / sum) * p;

//   for (slope=0.0, i=1; i<=xold.displayEnd(); i++)
//     slope += g(i) * p(i);
//   test = 0.0;
//   for (i=1; i<=xold.displayEnd(); i++) {
//     temp=fabs(p(i)) / maximum(fabs(xold(i)), 1.0);
//     if (temp > test) test = temp;
//   }
//   alamin = TOLX / test;
//   alam = 1.0;
//   for (;;) {
//     x = xold + alam * p;
//     vecfunc(x, auxData, fvec); 
//     f = fvec.dot(fvec);
//     if (alam < alamin) {
//       x = xold;
//       err = true;
//       return err;
//     } else if (f <= fold + ALF * alam * slope) return err;
//     else {
//       if (alam == 1.0)
// 	tmplam = -slope / (2.0 * (f - fold - slope));
//       else {
// 	rhs1 = f - fold - alam * slope;
// 	rhs2 = f2 - fold2 - alam2 * slope;
// 	a=(rhs1 / (alam * alam) - rhs2 / (alam2 * alam2)) / (alam - alam2);
// 	b=(-alam2 * rhs1 / (alam * alam) + alam * rhs2 / (alam2 * alam2)) / 
// 	  (alam - alam2);
// 	if (a == 0.0) tmplam = -slope / (2.0 * b);
// 	else {
// 	  disc = b * b - 3.0 * a * slope;
// 	  if (disc < 0.0) throw("Roundoff problem in lnsrch.\n");
// 	  else tmplam = (-b + sqrt(disc)) / (3.0 * a);
// 	}
// 	if (tmplam > 0.5 * alam)
// 	  tmplam = 0.5 * alam;
//       }
//     }
//     alam2 = alam;
//     f2 = f;
//     fold2 = fold;
//     alam = maximum(tmplam, 0.1 * alam);
//   }
//   return err;
// }


// bool newt(DoubleVector & x, const DoubleVector & auxData,
// 	  void (*vecfunc)(const DoubleVector &, const DoubleVector &, DoubleVector &)) {
//   bool err = false; 
//   const int MAXITS = 200;    ///< max iterations
//   double TOLF   = TOLERANCE; ///< convergence on function values
//   double TOLMIN = TOLF * 1.e-2; ///< spurious convergence to min of fmin
//   double TOLX   = TOLF * 1.e-3; ///< maximum dx convergence criterion
//   const double STPMX  = 100.0; ///< maximum step length allowed in line searches
   
//   int n = x.displayEnd();
//   int i,its,j,*indx;
//   double d,den,f,fold,stpmax,sum,temp,test;
  
//   indx = ivector(1, n);
//   DoubleMatrix fjac(n, n);
//   DoubleVector g(n), p(n), xold(n);
//   DoubleVector fvec(n);

//   vecfunc(x, auxData, fvec); 
//   f = 0.5 * fvec.dot(fvec);
//   test = 0.0;
//   test = fvec.apply(fabs).max();
//   if (test < 0.01 * TOLF) {
//     err = false;
//     free_ivector(indx, 1, n); return err;
//   }
//   sum += x.dot(x);
//   stpmax = STPMX * maximum(sqrt(sum), (double) n);
//   for (its=1; its<=MAXITS; its++) {
//     //    cout << its << endl; ///< DEBUG
//     fjac = fdjac(n, x, auxData, fvec, vecfunc);
//     for (i=1;i<=n;i++) {
//       for (sum=0.0, j=1; j<=n; j++) sum += fjac(j, i) * fvec(j);
//       g(i) = sum;
//     }
//     for (i=1; i<=n; i++) xold(i) = x(i);
//     fold = f;
//     for (i=1; i<=n; i++) p(i) = -fvec(i);
//     ludcmp(fjac, n, indx, d);
//     lubksb(fjac, n, indx, p);
//     err = lnsrch(xold, auxData, fold, g, p, x, f, stpmax, vecfunc, fvec);
//     test = 0.0;
//     for (i=1; i<=n; i++)
//       if (fabs(fvec(i)) > test) test = fabs(fvec(i));
//     if (test < TOLF) {
//       err = false;
//       free_ivector(indx, 1, n);
//       return err;
//     }
//     if (err) {
//       test = 0.0;
//       den = maximum(f, 0.5 * n);
//       for (i=1; i<=n; i++) {
// 	temp = fabs(g(i)) * maximum(fabs(x(i)), 1.0) / den;
// 	if (temp > test) test = temp;
//       }
//       err = (test < TOLMIN ? true : false);
//       free_ivector(indx, 1, n);
//       return err;
//     }
//     /// Check points aren't getting too close together
//     test = 0.0;
//     for (i=1; i<=n; i++) {
//       temp = (fabs(x(i) - xold(i))) / maximum(fabs(x(i)), 1.0);
//       if (temp > test) test = temp;
//     }
//     if (test < TOLX) { 
//       free_ivector(indx, 1, n); 
//       return err; 
//     }
//   }

//   return true;
// }

// DoubleMatrix fdjac(int n, DoubleVector x, const DoubleVector & auxData, const DoubleVector & fvec,
// 		   void (*vecfunc)(const DoubleVector &, const DoubleVector &, DoubleVector &)) {
//   double EPS = maximum(TOLERANCE, 1.0e-4);
//   int i,j;
//   double h,temp;
  
//   DoubleVector f(1, n);
//   DoubleMatrix df(n, n);
//   for (j=1; j<=n; j++) {
//     temp = x(j);
//     h = EPS * fabs(temp);
//     if (h == 0.0) h = EPS;
//     x(j) = temp + h;
//     h = x(j) - temp;
//     (*vecfunc)(x, auxData, f);
//     x(j) = temp;
//     for (i=1; i<=n; i++) df(i, j) = (f(i) - fvec.display(i)) / h;
//   }
//   return df;
// }
