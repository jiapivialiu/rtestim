#include <Eigen/Sparse>
#include <RcppEigen.h>
#include "linear_system.h"
#include "kf_utils.h"
#include "utils.h"
#include "dptf.h"
#include "dptf.h"
#include "admm.h"

typedef Eigen::COLAMDOrdering<int> Ord;
Eigen::SparseQR<SparseMatrix<double>, Ord> qradmm;

using Eigen::SparseMatrix;
using Eigen::SparseQR;
using Eigen::VectorXd;
using Eigen::MatrixXd;

using namespace Rcpp;

/**
 * ADMM for Gaussian trend filtering
 * @param M maximum iteration of the algos
 * @param n signal length
 * @param korder degree of Poisson trend filtering
 * @param y observed signals
 * @param x signal locations
 * @param w signal weights
 * @param theta primal variable of length `n`
 * @param z auxiliary variable of length `n-korder`
 * @param u dual variable of length `n-korder`
 * @param rho Lagrangian parameter of ADMM
 * @param lam_z hyperparameter of the auxiliary step of ADMM
 * @param DD D^T * D
 * @param tol tolerance of stopping criteria
 * @param iter interation index
 */
void admm_gauss(int M,
                int n,
                int korder,
                Rcpp::NumericVector const& y,
                Rcpp::NumericVector const& x,
                Rcpp::NumericVector const& w,
                Rcpp::NumericVector& theta,
                Rcpp::NumericVector& z,
                Rcpp::NumericVector& u,
                double rho,
                double lam_z,
                Eigen::SparseMatrix<double> const& DD,
                double tol,
                int linear_solver,
                int& iter) {
  
  // update the weights
  NumericVector W = exp(theta) * w;         
  
  // temp var
  double r_norm = 0.0;
  double s_norm = 0.0;
  NumericVector z_old = clone(z);
  NumericVector tmp_n(n);
  NumericVector tmp_m(z.size());
  NumericVector Dth(z.size());
  VectorXd tmp_theta(n);
  //VectorXd mn(n - korder);
  
  // Eigen var
  VectorXd eW = nvec_to_evec(W);
  VectorXd ey = nvec_to_evec(y);
  VectorXd ex = nvec_to_evec(x);
  VectorXd adj_mean = nvec_to_evec(z - u);
  
  // Linear system components
  MatrixXd denseD;
  VectorXd s_seq;
  bool equal_space;
  int computation_info;
  
  // Initialize denseD and s_seq for KF solver BEFORE construct/compute
  if (linear_solver == 2) {
    equal_space = is_equal_space(ex, std::sqrt(Eigen::NumTraits<double>::epsilon()));
    denseD = MatrixXd::Zero(n - korder, korder + 1);
    s_seq = equal_space ? VectorXd::Zero(1) : VectorXd::Zero(n);
    configure_denseD(ex, denseD, s_seq, korder, equal_space);
  }
  
  //SparseMatrix<double> cDD;
  //cDD = DD * n * rho;  // a copy that doesn't change
  //for (int i = 0; i < n; i++) {
  //  cDD.diagonal()(i) += eW(i);
  //}
  //qradmm.compute(cDD);
  LinearSystem linear_system;
  linear_system.construct(ey, eW, korder, n * rho, DD, denseD, s_seq, linear_solver);
  linear_system.compute(linear_solver);

  for (iter = 0; iter < M; iter++) {
    if (iter % 1000 == 0)
      Rcpp::checkUserInterrupt();
    // solve for primal variable - theta:
    //tmp_n = doDtv(z - u, korder, x) * n * rho;
    //tmp_n += W * y;
    //tmp_theta = nvec_to_evec(tmp_n);
    //tmp_theta = qradmm.solve(tmp_theta);
    //theta = evec_to_nvec(tmp_theta);
    // theta update
    std::tie(tmp_theta, computation_info) = linear_system.solve(ey, eW,
        adj_mean, korder, x, n * rho, denseD, s_seq, linear_solver, equal_space);
    theta = evec_to_nvec(tmp_theta);

    // solve for alternating variable - z:
    Dth = doDv(theta, korder, x);
    tmp_m = Dth + u;
    z = rcpp_tvd(tmp_m, lam_z);

    // update dual variable - u:
    u += Dth - z;

    // primal residuals:
    r_norm = sqrt(mean(pow(Dth - z, 2)));
    // dual residuals:
    tmp_n = doDtv(z - z_old, korder, x);
    s_norm = rho * sqrt(mean(pow(tmp_n, 2)));
    // stopping criteria check:
    if (r_norm < tol && s_norm < tol)
      break;

    // auxiliary variables update:
    z_old = z;
  }
}

void prox_newton(int M,
                 int& Minner,
                 int Mline,
                 int n,
                 int korder,
                 Rcpp::NumericVector const& y,
                 Rcpp::NumericVector const& x,
                 Rcpp::NumericVector const& w,
                 Rcpp::NumericVector& theta,
                 Rcpp::NumericVector& z,
                 Rcpp::NumericVector& u,
                 double lambda,
                 double rho,
                 double alpha,
                 double gamma,
                 Eigen::SparseMatrix<double> const& DD,
                 double tol,
                 int linear_solver,
                 int& total_iter) {
  double s;                       // step size
  NumericVector obj_list(M + 1);  // objective list for each iterate
  double obj = 1e4;               // initialize it to be large
  NumericVector theta_old(n);     // a buffer for line search
  double lam_z = lambda / rho;
  int inner_iter;

  // initialize objective
  obj = pois_obj(korder, y, x, w, theta, lambda);
  obj_list(0) = obj;
  int iter_best = 0;
  NumericVector std_y(y);

  for (int iter = 0; iter < M; iter++) {
    if (iter % 50 == 0)
      Rcpp::checkUserInterrupt();
    theta_old = theta;

    // define new(centered) data for least squares problem
    std_y = centered_data(y, w, theta);
    // solve least squares problem (Gaussian TF)
    admm_gauss(Minner, n, korder, std_y, x, w, theta, z, u, rho, lam_z, DD, tol, 
      linear_solver, inner_iter);
    total_iter += inner_iter;
    //  line search for step size
    s = line_search(s, lambda, alpha, gamma, y, x, w, n, korder, theta, theta_old,
                    Mline);
    if (s < 0)
      break;
    // update theta
    theta = theta * s + (1 - s) * theta_old;

    // compute objective
    obj = pois_obj(korder, y, x, w, theta, lambda);
    obj_list(iter + 1) = obj;

    // check stopping criteria
    if (obj < obj_list(iter_best)) {
      if (obj_list(iter_best) - obj <= abs(obj_list(iter_best)) * tol) {
        iter_best = iter + 1;
        break;
      }
      iter_best = iter + 1;
    }

    // adjust the iterate steps
    if (iter >= iter_best + 4)
      break;
  }
}

/**
 * This is a wrapper around admm_gauss to use in test_that()
 * Tests the Gaussian ADMM step directly
 */
// [[Rcpp::export]]
Rcpp::List admm_gauss_testing(int M,
                              int korder,
                              Rcpp::NumericVector const& y,
                              Rcpp::NumericVector const& x,
                              Rcpp::NumericVector const& w,
                              Rcpp::NumericVector theta,
                              Rcpp::NumericVector z,
                              Rcpp::NumericVector u,
                              double rho,
                              double lam_z,
                              int linear_solver,
                              double tol) {
  Eigen::SparseMatrix<double> Dk;
  Eigen::SparseMatrix<double> DkDk;
  Dk = get_Dtil(korder, x);
  DkDk = Dk.transpose() * Dk;
  int n = y.size();
  int iter = 0;
  
  // Store initial values
  NumericVector theta_init = clone(theta);
  NumericVector z_init = clone(z);
  NumericVector u_init = clone(u);
  
  // Run ADMM
  admm_gauss(M, n, korder, y, x, w, theta, z, u, rho, lam_z, DkDk, tol, 
             linear_solver, iter);
  
  // Check for NAs and Infs
  int num_na_theta = 0;
  int num_inf_theta = 0;
  int num_na_z = 0;
  int num_inf_z = 0;
  
  for (int i = 0; i < n; i++) {
    if (std::isnan(theta[i])) num_na_theta++;
    if (std::isinf(theta[i])) num_inf_theta++;
  }
  
  for (int i = 0; i < z.size(); i++) {
    if (std::isnan(z[i])) num_na_z++;
    if (std::isinf(z[i])) num_inf_z++;
  }
  
  // Calculate changes
  NumericVector theta_change = theta - theta_init;
  double max_theta_change = max(abs(theta_change));
  double mean_theta_change = mean(abs(theta_change));
  
  List out = List::create(
    Named("theta") = theta,
    Named("z") = z,
    Named("u") = u,
    Named("theta_init") = theta_init,
    Named("z_init") = z_init,
    Named("u_init") = u_init,
    Named("niter") = iter,
    Named("converged") = (iter < M),
    Named("num_na_theta") = num_na_theta,
    Named("num_inf_theta") = num_inf_theta,
    Named("num_na_z") = num_na_z,
    Named("num_inf_z") = num_inf_z,
    Named("max_theta_change") = max_theta_change,
    Named("mean_theta_change") = mean_theta_change,
    Named("solver") = (linear_solver == 1 ? "QR" : "KF")
  );
  return out;
}

/**
 * This is a wrapper around the void function to use in test_that()
 */
// [[Rcpp::export]]
Rcpp::List prox_newton_testing(int M,
                               int Minner,
                               int Mline,
                               int korder,
                               Rcpp::NumericVector const& y,
                               Rcpp::NumericVector const& x,
                               Rcpp::NumericVector const& w,
                               double lambda,
                               double ls_alpha,
                               double ls_gamma,
                               int linear_solver,
                               double tol) {
  Eigen::SparseMatrix<double> Dk;
  Eigen::SparseMatrix<double> DkDk;
  Dk = get_Dtil(korder, x);
  DkDk = Dk.transpose() * Dk;
  int m = Dk.rows();
  int n = y.size();
  NumericVector beta(n);
  NumericVector z(m);
  NumericVector u(m);
  double rho = lambda;
  int iter = 0;
  prox_newton(M, Minner, Mline, n, korder, y, x, w, beta, z, u, lambda, rho,
              ls_alpha, ls_gamma, DkDk, tol, linear_solver, iter);
  List out = List::create(Named("lambda") = lambda, Named("theta") = exp(beta),
                          Named("niter") = iter);
  return out;
}
