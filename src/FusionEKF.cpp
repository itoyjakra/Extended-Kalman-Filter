#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF()
{
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
  ekf_.F_ = MatrixXd(4, 4);
	ekf_.F_ << 1, 0, 1, 0,
			  0, 1, 0, 1,
			  0, 0, 1, 0,
			  0, 0, 0, 1;

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack)
{

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_)
  {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 0, 0, 0, 0;
    ekf_.P_ = MatrixXd(4, 4);
    ekf_.P_ << 1, 0, 0, 0,
               0, 1, 0, 0,
               0, 0, 1, 0,
               0, 0, 0, 1;

    std::cout << "measurement_pack" << '\n';
    std::cout << measurement_pack.sensor_type_ << '\n' << measurement_pack.raw_measurements_ << '\n';
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
    {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      ekf_.x_ = tools.PolarToCart(measurement_pack.raw_measurements_);
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER)
    {
      /**
      Initialize state.
      */
      VectorXd pos = measurement_pack.raw_measurements_;
      ekf_.x_(0) = pos(0);
      ekf_.x_(1) = pos(1);
      //ekf_.x_ = measurement_pack.raw_measurements_;
      //ekf_.Init(x, P, F, H, R, Q);
    }

    // done initializing, no need to predict or update
    previous_timestamp_ = measurement_pack.timestamp_;
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
 	previous_timestamp_ = measurement_pack.timestamp_;
  cout << "timestamp = " << previous_timestamp_ << ", dt = " << dt << endl;

  ekf_.F_(0, 2) = dt;
	ekf_.F_(1, 3) = dt;
  cout << "F = " << endl;
  cout << ekf_.F_ <<endl;

  MatrixXd G = MatrixXd(4, 2);
 	float dt2 = dt*dt;
  G << dt2/2, 0,
 			 0, dt2/2,
 			 dt, 0,
 			 0, dt;
 	cout << "G = " << endl;
 	cout << G << endl;

 	MatrixXd Gt = G.transpose();
 	MatrixXd Q_nu = MatrixXd(2, 2);
  float noise_ax = 9.0;
  float noise_ay = 9.0;
 	Q_nu << noise_ax, 0,
 				  0, noise_ay;
 	ekf_.Q_ = G * Q_nu * Gt;
  cout << "Q = " << endl;
  cout << ekf_.Q_ << endl;


  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
  {
    // Radar updates
    std::cout << "---------begin radar update------------" << '\n';
    try
    {
      Hj_ = tools.CalculateJacobian(ekf_.x_);
    }
    catch (const char* msg)
    {
      cout << msg << endl;
      Hj_ << 0, 0, 0, 0,
             0, 0, 0, 0,
             0, 0, 0, 0;
    }
    //Hj_ = tools.CalculateJacobian(ekf_.x_);
    std::cout << "Jacobian is here" << '\n';
    cout << "Hj = " << endl;
    cout << Hj_ << endl;
    
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_, Hj_);
    std::cout << "---------end radar update------------" << '\n';
  } else {
    // Laser updates
    std::cout << "---------begin laser update------------" << '\n';
    ekf_.H_ = MatrixXd(2, 4);
	  ekf_.H_ << 1, 0, 0, 0,
			        0, 1, 0, 0;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
    std::cout << "---------end laser update------------" << '\n';
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
