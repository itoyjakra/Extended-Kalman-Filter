#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {

  x_ = F_ * x_;
	MatrixXd Ft = F_.transpose();
	P_ = F_ * P_ * Ft + Q_;
  std::cout << "predicted P_ = " << std::endl;
  std::cout << P_ << std::endl;
}

void KalmanFilter::Update(const VectorXd &z) {

  VectorXd z_pred = H_ * x_;
  std::cout << "z_pred" << '\n';
  std::cout << z_pred << '\n';
	VectorXd y = z - z_pred;
  std::cout << "y" << '\n';
  std::cout << y << '\n';
	MatrixXd Ht = H_.transpose();
  std::cout << "Ht" << '\n';
  std::cout << Ht << '\n';
  std::cout << "H" << '\n';
  std::cout << H_ << '\n';
  std::cout << "P_" << '\n';
  std::cout << P_ << '\n';
  std::cout << "H_ * P * Ht" << '\n';
  std::cout << H_ * P_ * Ht << '\n';
  std::cout << "R" << '\n';
  std::cout << R_ << '\n';
	MatrixXd S = H_ * P_ * Ht + R_;
  std::cout << "S" << '\n';
  std::cout << S << '\n';
	MatrixXd Si = S.inverse();
  std::cout << "Si" << '\n';
  std::cout << Si << '\n';
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;
  std::cout << "K" << '\n';
  std::cout << K << '\n';

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z, MatrixXd &Hj) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */

  std::cout << "starting UpdateEKF" << '\n';
  std::cout << "z" << '\n';
  std::cout << z << '\n';
  VectorXd z_pred = tools.CartToPolar(x_);
  std::cout << "z_pred" << '\n';
  std::cout << z_pred << '\n';
  VectorXd y = z - z_pred;
  std::cout << "y" << '\n';
  std::cout << y << '\n';
	MatrixXd Ht = Hj.transpose();
	MatrixXd S = Hj * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * Hj) * P_;

}
