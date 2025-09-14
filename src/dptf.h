#ifndef __DPTF_H
#define __DPTF_H

Rcpp::NumericVector rcpp_tvdz(
    Rcpp::NumericVector y,
    Rcpp::NumericVector z,
    double lam
);

/* Unused, but provided */
Rcpp::NumericVector rcpp_wtvdz(
    Rcpp::NumericVector y,
    Rcpp::NumericVector z,
    double lam, Rcpp::NumericVector x
);

void tvd(int n, double *y, double lambda, double *theta);
Rcpp::NumericVector rcpp_tvd(Rcpp::NumericVector y, double lambda);

void wtvd(int n, double *y, double lambda, double *w, double *theta);
Rcpp::NumericVector rcpp_wtvd(Rcpp::NumericVector y, double lambda, Rcpp::NumericVector weights);

#endif
