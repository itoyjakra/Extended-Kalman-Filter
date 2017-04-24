#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) 
{
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
    throw "position vector is too small";
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
  float small = 1.0e-6;
  float px = cart(0);
  float py = cart(1);
  float vx = cart(2);
  float vy = cart(3);
  float rho = sqrt(px*px + py*py);
  float phi = atan2(py, px);
  float rho_dot;
  VectorXd polar(3);

  if (rho < small)
    rho_dot = 0;
  else
    rho_dot = (px * vx + py * vy)/rho;

  polar << rho, phi, rho_dot;

  return polar;
}
