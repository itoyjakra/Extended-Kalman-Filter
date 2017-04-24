#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include "Eigen/Dense"

class Tools {
public:
  Tools();
  virtual ~Tools();
  Eigen::VectorXd CalculateRMSE(const std::vector<Eigen::VectorXd> &estimations, const std::vector<Eigen::VectorXd> &ground_truth);
  Eigen::MatrixXd CalculateJacobian(const Eigen::VectorXd& x_state);
  // Polar to Cartesian converter
  Eigen::VectorXd PolarToCart(const Eigen::VectorXd& polar);
  // Cartesian to Polar converter
  Eigen::VectorXd CartToPolar(const Eigen::VectorXd& cart);
};

#endif /* TOOLS_H_ */
