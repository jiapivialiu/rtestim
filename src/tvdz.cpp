#include <Rcpp.h>
#include "dptf.h"
using namespace Rcpp;

/**
 *  Solves the problem
 *  min_beta n^{-1} sum_{i=1}^n -yi * theta_i + z_i beta_i + \lambda ||D\theta||_1
 *  where beta = exp(theta) = E[y], and D is the 0-th order discrete difference
 *  operator
 *
 * @param n                    number of observations
 * @param y                    response vector
 * @param z                    covariate vector
 * @param lam                  penalty parameter
 * @param beta                 allocated space for the output
 */
void tvdz(unsigned int n, double* y, double* z, double lam, double* beta) {
  if (n == 0) return;

  double* yz = new double[n];
  for (int i = 0; i < n; i++) yz[i] = y[i] / z[i];

  if (n == 1 || lam == 0) {
    for (int i = 0; i < n; i++) beta[i] = yz[i];
    delete[] yz;
    return;
  }

  double* x = new double[2*n];
  double* a = new double[2*n];
  double* b = new double[2*n];
  double* tm = new double[n-1];
  double* tp = new double[n-1];


  int l;
  int r;
  int lo;
  int hi;
  double afirst;
  double bfirst;
  double alast;
  double blast;
  double alo;
  double blo;
  double ahi;
  double bhi;

  double* lamz = new double[n];
  for (int i = 0; i < n; i++) lamz[i] = lam / z[i];

  /* We step through the first iteration manually */
  tm[0] = -lamz[0] + yz[0];
  tp[0] = lamz[0] + yz[0];
  l = n - 1;
  r = n;
  x[l] = tm[0];
  x[r] = tp[0];
  a[l] = 1;
  b[l] = -yz[0] + lamz[1];
  a[r] = -1;
  b[r] = yz[0] + lamz[1];
  afirst = 1;
  bfirst = -yz[1] - lamz[1];
  alast = -1;
  blast = yz[1] - lamz[1];

  /* Now iterations 2 through n-1 */
  for (int k = 1; k < n - 1; k++) {
    /* Compute lo: step up from l until the derivative is greater than -lam */
    alo = afirst;
    blo = bfirst;
    for (lo = l; lo <= r; lo++) {
      if (alo * x[lo] + blo > -lamz[k + 1]) break;
      alo += a[lo];
      blo += b[lo];
    }

    /* Compute hi: step down from r until the derivative is less than lam */
    ahi = alast;
    bhi = blast;
    for (hi = r; hi >= lo; hi--) {
      if (-ahi * x[hi] - bhi < lamz[k + 1]) break;
      ahi += a[hi];
      bhi += b[hi];
    }

    /* Compute the negative knot */
    tm[k] = (-lamz[k + 1] - blo) / alo;
    l = lo - 1;
    x[l] = tm[k];

    /* Compute the positive knot */
    tp[k] = (lamz[k + 1] + bhi) / (-ahi);
    r = hi + 1;
    x[r] = tp[k];

    /* Update a and b */
    a[l] = alo;
    b[l] = blo + lamz[k + 1];
    a[r] = ahi;
    b[r] = bhi + lamz[k + 1];
    afirst = 1;
    bfirst = -yz[k + 1] - lamz[k + 1];
    alast = -1;
    blast = yz[k + 1] - lamz[k + 1];
  }

  /* Compute the last coefficient: here, the function has zero derivative */
  alo = afirst;
  blo = bfirst;
  for (lo = l; lo <= r; lo++) {
    if (alo * x[lo] + blo > 0) break;
    alo += a[lo];
    blo += b[lo];
  }
  beta[n - 1] = -blo / alo;

  /* Compute the rest of the coefficients, by the back-pointers */
  for (int k = n - 2; k >= 0; k--) {
    if (beta[k + 1] > tp[k]) beta[k] = tp[k];
    else if (beta[k + 1] < tm[k]) beta[k] = tm[k];
    else beta[k] = beta[k + 1];
  }

  delete[] x;
  delete[] a;
  delete[] b;
  delete[] tm;
  delete[] tp;
  delete[] yz;
  delete[] lamz;
}

// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::export]]
NumericVector rcpp_tvdz(NumericVector y, NumericVector z, double lambda) {
  int n = y.size();
  NumericVector beta(n);
  lambda *= n;
  tvdz(n, y.begin(), z.begin(), lambda, beta.begin());
  return beta;
}

/* The following code is from `glmgen/tvdenoising/src/tvd.cpp`. */
void tvd(int n, double *y, double lambda, double *theta) {
  /* Take care of a few trivial cases */
  if (n == 0) return;
  if (n == 1 || lambda == 0) {
    for (int i = 0; i < n; i++) theta[i] = y[i];
    return;
  }

  double *x = new double[2*n];
  double *a = new double[2*n];
  double *b = new double[2*n];
  double *tm = new double[n-1];
  double *tp = new double[n-1];

  int l;
  int r;
  int lo;
  int hi;
  double afirst;
  double bfirst;
  double alast;
  double blast;
  double alo;
  double blo;
  double ahi;
  double bhi;

  /* We step through the first iteration manually */
  tm[0] = -lambda+y[0];
  tp[0] = lambda+y[0];
  l = n-1;
  r = n;
  x[l] = tm[0];
  x[r] = tp[0];
  a[l] = 1;
  b[l] = -y[0]+lambda;
  a[r] = -1;
  b[r] = y[0]+lambda;
  afirst = 1;
  bfirst = -y[1]-lambda;
  alast = -1;
  blast = y[1]-lambda;

  /* Now step through iterations 2 through n-1 */
  for (int k = 1; k < n-1; k++) {
    /* Compute lo: step up from l until the
       derivative is greater than -lambda */
    alo = afirst;
    blo = bfirst;
    for (lo=l; lo<=r; lo++) {
      if (alo*x[lo]+blo > -lambda) break;
      alo += a[lo];
      blo += b[lo];
    }

    /* Compute hi: step down from r until the
       derivative is less than lambda */
    ahi = alast;
    bhi = blast;
    for (hi = r; hi >= lo; hi--) {
      if (-ahi*x[hi]-bhi < lambda) break;
      ahi += a[hi];
      bhi += b[hi];
    }

    /* Compute the negative knot */
    tm[k] = (-lambda-blo)/alo;
    l = lo-1;
    x[l] = tm[k];

    /* Compute the positive knot */
    tp[k] = (lambda+bhi)/(-ahi);
    r = hi+1;
    x[r] = tp[k];

    /* Update a and b */
    a[l] = alo;
    b[l] = blo+lambda;
    a[r] = ahi;
    b[r] = bhi+lambda;
    afirst = 1;
    bfirst = -y[k+1]-lambda;
    alast = -1;
    blast = y[k+1]-lambda;
  }

  /* Compute the last coefficient: this is where
     the function has zero derivative */
  alo = afirst;
  blo = bfirst;
  for (lo = l; lo <= r; lo++) {
    if (alo*x[lo]+blo > 0) break;
    alo += a[lo];
    blo += b[lo];
  }
  theta[n-1] = -blo/alo;

  /* Compute the rest of the coefficients, by the
     back-pointers */
  for (int k = n-2; k >= 0; k--) {
    if (theta[k+1] > tp[k]) theta[k] = tp[k];
    else if (theta[k+1] < tm[k]) theta[k] = tm[k];
    else theta[k] = theta[k+1];
  }

  delete[] x;
  delete[] a;
  delete[] b;
  delete[] tm;
  delete[] tp;
}


// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::export]]
NumericVector rcpp_tvd(NumericVector y, double lambda) {
  int n = y.size();
  NumericVector theta(n);
  tvd(n, y.begin(), lambda, theta.begin());
  return theta;
}



/* The following code is from `glmgen/tvdenoising/src/wtvd.cpp`. */
void wtvd(int n, double *y, double lambda, double *w, double *theta) {
  /* Take care of a few trivial cases */
  if (n == 0) return;
  if (n == 1 || lambda == 0) {
    for (int i = 0; i < n; i++) theta[i] = y[i];
    return;
  }

  double *x = new double[2*n];
  double *a = new double[2*n];
  double *b = new double[2*n];
  double *tm = new double[n-1];
  double *tp = new double[n-1];

  int l;
  int r;
  int lo;
  int hi;
  double afirst;
  double bfirst;
  double alast;
  double blast;
  double alo;
  double blo;
  double ahi;
  double bhi;

  /* We step through the first iteration manually */
  tm[0] = -lambda/w[0]+y[0];
  tp[0] = lambda/w[0]+y[0];
  l = n-1;
  r = n;
  x[l] = tm[0];
  x[r] = tp[0];
  a[l] = w[0];
  b[l] = -w[0]*y[0]+lambda;
  a[r] = -w[0];
  b[r] = w[0]*y[0]+lambda;
  afirst = w[1];
  bfirst = -w[1]*y[1]-lambda;
  alast = -w[1];
  blast = w[1]*y[1]-lambda;

  /* Now step through iterations 2 through n-1 */
  for (int k = 1; k < n-1; k++) {
    /* Compute lo: step up from l until the
       derivative is greater than -lambda */
    alo = afirst;
    blo = bfirst;
    for (lo=l; lo<=r; lo++) {
      if (alo*x[lo]+blo > -lambda) break;
      alo += a[lo];
      blo += b[lo];
    }

    /* Compute hi: step down from r until the
       derivative is less than lambda */
    ahi = alast;
    bhi = blast;
    for (hi = r; hi >= lo; hi--) {
      if (-ahi*x[hi]-bhi < lambda) break;
      ahi += a[hi];
      bhi += b[hi];
    }

    /* Compute the negative knot */
    tm[k] = (-lambda-blo)/alo;
    l = lo-1;
    x[l] = tm[k];

    /* Compute the positive knot */
    tp[k] = (lambda+bhi)/(-ahi);
    r = hi+1;
    x[r] = tp[k];

    /* Update a and b */
    a[l] = alo;
    b[l] = blo+lambda;
    a[r] = ahi;
    b[r] = bhi+lambda;
    afirst = w[k+1];
    bfirst = -w[k+1]*y[k+1]-lambda;
    alast = -w[k+1];
    blast = w[k+1]*y[k+1]-lambda;
  }

  /* Compute the last coefficient: this is where
     the function has zero derivative */
  alo = afirst;
  blo = bfirst;
  for (lo = l; lo <= r; lo++) {
    if (alo*x[lo]+blo > 0) break;
    alo += a[lo];
    blo += b[lo];
  }
  theta[n-1] = -blo/alo;

  /* Compute the rest of the coefficients, by the
     back-pointers */
  for (int k = n-2; k >= 0; k--) {
    if (theta[k+1] > tp[k]) theta[k] = tp[k];
    else if (theta[k+1] < tm[k]) theta[k] = tm[k];
    else theta[k] = theta[k+1];
  }

  delete[] x;
  delete[] a;
  delete[] b;
  delete[] tm;
  delete[] tp;
}

// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::export]]
NumericVector rcpp_wtvd(NumericVector y, double lambda, NumericVector weights) {
  int n = y.size();
  NumericVector theta(n);
  wtvd(n, y.begin(), lambda, weights.begin(), theta.begin());
  return theta;
}
