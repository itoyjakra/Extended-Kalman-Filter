#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {

  VectorXd rmse(4);
  rmse << 0,0,0,0;
  VectorXd c;

	if (estimations.size() != ground_truth.size())
  {
    throw "lengths of vectors are different for RMSE calculations";
  }
  else if (estimations.size() == 0)
  {
    throw "estimate is unavailable";
  }
  else if (ground_truth.size() == 0)
  {
    throw "ground truth is unavailable";
  }
  else
  {
    for (int i=0; i<estimations.size(); i++)
    {
      c = estimations[i] - ground_truth[i];
      rmse.array() = rmse.array() + c.array() * c.array();
    }
  }
  rmse = rmse / estimations.size();
  rmse = sqrt(rmse.array());

	return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {

  const float small = 1.0e-6;
  MatrixXd Hj(3,4);

	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	float pmag = px*px + py*py;
  float pmag_root = pow(pmag, 0.5);
  float pmag_root3 = pmag * pmag_root;

	if (fabs(pmag) < small)
  {
    std::cout << "" << '\n';
    throw "division by zero encountered";
    return Hj;
  }
  else
  {
    Hj << px/pmag_root, py/pmag_root, 0, 0,
          -py/pmag, px/pmag, 0, 0,
          py*(vx*py - vy*px)/pmag_root3, px*(vy*px - vx*py)/pmag_root3, px/pmag_root, py/pmag_root;

  	return Hj;
  }
}

VectorXd Tools::PolarToCart(const VectorXd& polar)
{
  VectorXd cartesian(4);
  float rho = polar(0);
  float phi = polar(1);
  float rho_dot = polar(2);

  cartesian << rho * cos(phi), rho * sin(phi), rho_dot * cos(phi), rho_dot * sin(phi);

  return cartesian;
}

VectorXd Tools::CartToPolar(const VectorXd& cart)
{
  std::cout << "in CartToPolar at the begining" << '\n';
  float px = cart(0);
  float py = cart(1);
  float vx = cart(2);
  float vy = cart(3);
  std::cout << "px, py, vx, vy" << '\n';
  std::cout << px << '\t' << py << 't' << vx << '\t' << vy << '\n';
  float pmag = sqrt(px*px + py*py);
  VectorXd polar(3);
  

  std::cout << "in CartToPolar" << '\n';

  polar << pmag,
           atan(py/px),
           (px * vx + py * vy)/pmag;
  std::cout << "polar" << '\n';
  std::cout << polar << '\n';
  return polar;
}
