#ifndef __UTILS_H
#define __UTILS_H

#include <Eigen/Sparse>
#include <Rcpp.h>
#include <RcppEigen.h>
#include <dspline.h>
using Rcpp::NumericVector;
using Eigen::SparseMatrix;
using Eigen::VectorXd;


Rcpp::NumericVector evec_to_nvec(Eigen::VectorXd evec);
Eigen::VectorXd nvec_to_evec(Rcpp::NumericVector nvec);
Eigen::SparseMatrix<double> get_Dtil(int k, Rcpp::NumericVector xd);
Eigen::SparseMatrix<double> get_D(int k, Rcpp::NumericVector xd);
Rcpp::NumericVector doDv(Rcpp::NumericVector v, int k, Rcpp::NumericVector xd);
Rcpp::NumericVector doDtv(Rcpp::NumericVector v, int k, Rcpp::NumericVector xd);
void create_lambda(Rcpp::NumericVector& lambda,
                   double& lambdamin,
                   double& lambdamax,
                   double& lambda_min_ratio,
                   int& nsol);
Rcpp::NumericVector create_lambda_test(Rcpp::NumericVector lambda,
                                       double lambdamin,
                                       double lambdamax,
                                       double lambda_min_ratio,
                                       int nsol);
double one_norm(Rcpp::NumericVector const& z);
double pois_obj(int ord,
                Rcpp::NumericVector const& y,
                Rcpp::NumericVector const& x,
                Rcpp::NumericVector const& w,
                Rcpp::NumericVector& theta,
                double lambda);

Rcpp::NumericVector centered_data(Rcpp::NumericVector const& y,
                                  Rcpp::NumericVector const& w,
                                  Rcpp::NumericVector& theta);

double line_search(double s,
                   double lambda,
                   double alpha,
                   double gamma,
                   Rcpp::NumericVector const& y,
                   Rcpp::NumericVector const& x,
                   Rcpp::NumericVector const& w,
                   int n,
                   int ord,
                   Rcpp::NumericVector& theta,
                   Rcpp::NumericVector& theta_old,
                   int M);

bool is_equal_space(const Eigen::VectorXd& x, double space_tolerance_ratio);

/* Matrix construction */
Eigen::SparseMatrix<double> identity(int n);
Eigen::SparseMatrix<double> diagonal(Eigen::ArrayXd diag);
Eigen::SparseMatrix<double> get_dk_mat(int k, NumericVector xd,
    bool tf_weighting);
Eigen::SparseMatrix<double> get_penalty_mat(int k, NumericVector xd);
Eigen::VectorXd legendre_polynomial(Eigen::VectorXd, int k, double a,
    double b);

/* Polynomial subspace projection */
Eigen::VectorXd project_polynomials(const NumericVector& x, const VectorXd& y,
    const Eigen::ArrayXd& weights, int k);

/* Tridiagonal matrix solve */
std::tuple<Eigen::VectorXd,Eigen::VectorXd,Eigen::VectorXd> extract_tridiag(
    Eigen::SparseMatrix<double> A);
Eigen::VectorXd tridiag_forward(const Eigen::VectorXd& a,
        const Eigen::VectorXd& b, const Eigen::VectorXd& c);
Eigen::VectorXd tridiag_backsolve(
        const Eigen::VectorXd& a, const::VectorXd& b,
            const Eigen::VectorXd& cp, const Eigen::VectorXd& d);

// Lambda sequence
double get_lambda_max(const NumericVector& x, const Eigen::VectorXd& y,
                      const Eigen::ArrayXd& sqrt_weights, int k);
void get_lambda_seq(Eigen::VectorXd& lambda, double lambda_max,
                    double lambda_min, double lambda_min_ratio, int n_lambda);
Eigen::VectorXd get_lambda_seq_r(Eigen::VectorXd lambda, double lambda_max,
                                 double lambda_min, double lambda_min_ratio,
                                 int n_lambda);

int calc_degrees_of_freedom(Eigen::VectorXd const &v, int k, double tol = 1e-8);



// Workarounds for interfacing with dspline / tvdenoising
Eigen::VectorXd Dkv(Eigen::VectorXd v, int k, const NumericVector& xd,
                    bool tf_weighting = false);
Eigen::VectorXd Dktv(Eigen::VectorXd v, int k, const NumericVector& xd);
Eigen::VectorXd tf_dp(Eigen::VectorXd v, double lambda);



#endif
