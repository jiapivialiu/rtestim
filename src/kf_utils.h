#ifndef KF_UTILS_H
#define KF_UTILS_H

Eigen::MatrixXd computePtemp(Eigen::MatrixXd A, Eigen::MatrixXd P);
Eigen::MatrixXd smat_to_mat(const Eigen::SparseMatrix<double>& sparseMat, int k, bool equal_spaced);
void configure_denseD(Eigen::VectorXd& x, Eigen::MatrixXd& denseD, Eigen::VectorXd& s_seq, int& k, bool& equal_space);
Rcpp::List configure_denseD_test(Eigen::VectorXd x, int k);
void f1step(double y, double c, double Z, double H, const Eigen::MatrixXd& A, 
  double RQR, Eigen::VectorXd& a, Eigen::MatrixXd& P, double& vt, double& Ft,
  Eigen::VectorXd& Kt);
void df1step(double y, double Z, double H, const Eigen::MatrixXd& A, double RQR,
  Eigen::VectorXd& a, Eigen::MatrixXd& P, Eigen::MatrixXd& Pinf, int& rankp,
  double& vt, double& Ft, double& Finf, Eigen::VectorXd& Kt, Eigen::VectorXd& Kinf);

#endif
