#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
	// Gravity (g) is 9 m/s^2. So max change maybe 1g. So STD should be about
	// half this, 4.5 m/s^2.  The lecture suggestion of 3.0 seems close enough
	// so I'll just keep it.
  std_a_ = 3.0;

  // Process noise standard deviation yaw acceleration in rad/s^2
  // Methodology:  Estimate maximum change in circle radius that can
	// occur in 1 second.  Make STD half this value.
	// Max change: From very tight circle to going straight.  We guess
	// that changing from pos to neg is not common.  Straight means psi_dot = 0
	// and a very tight circle is one round in one second, so 2*PI / sec.
	// So max change = (2*pi/sec - 0) = 2*pi / sec.  Half this is pi/sec.
	std_yawdd_ = 3.14;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

	is_initialized_ = false;

	// State dimension (checked)
	n_x_ = 5;

	// Augmented state dimension (checked)
	n_aug_ = 7;

	// Sigma point spreading parameter (checked)
	lambda_ = 3 - n_aug_;

	// the current NIS for radar
	NIS_radar_ = 0;

	// the current NIS for laser
	NIS_laser_ = 0;


	// predicted sigma points matrix
	Xsig_pred_ = MatrixXd(n_x_,2 * n_aug_ + 1);

	// time when the state is true, in microseconds
	time_us_ = 0;

	// Weights of sigma points
	weights_ = VectorXd(2*n_aug_+1);	
	weights_(0) = lambda_ / (lambda_ + n_aug_);
	for (int i = 1; i<2 * n_aug_ + 1; i++) {		
		weights_(i) = 0.5 / (n_aug_ + lambda_);
	}

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

	if (!is_initialized_) {

		// first measurement
		//ekf_.x_ = VectorXd(4);
		//ekf_.F_ = MatrixXd(4, 4);
		//ekf_.Q_ = MatrixXd(4, 4);
		//ekf_.x_ << 1, 1, 1, 1;

		if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
			/**
			Convert radar from polar to CTRV coordinates and initialize state.
			*/

			// Extract measurements 
			double rho = meas_package.raw_measurements_[0];
			double phi = meas_package.raw_measurements_[1];
			double rho_dot = meas_package.raw_measurements_[2];

			// Map to hidden state space
			double px = rho*cos(phi);
			double py = rho*sin(phi);
			double v = 0;
			double psi = 0;
			double psi_dot = 0;
					
			x_ << px, py, v, psi, psi_dot;

			P_ << std_radr_*std_radr_, 0,                   0, 0, 0,
				    0,                   std_radr_*std_radr_, 0, 0, 0,
				    0,                   0,                   1, 0, 0,
				    0,                   0,                   0, 1, 0,
				    0,                   0,                   0, 0, 1;
		}
		else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {

			// Map to hidden state space
			x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1],5,0,0;

			P_ << std_laspx_*std_laspx_, 0,                     0,    0,    0,
				    0,                     std_laspy_*std_laspy_, 0,    0,    0,
				    0,                     0,                     1, 0,    0,
				    0,                     0,                     0,    1, 0,
				    0,                     0,                     0,    0,    1;
		}

		
		
		// done initializing, no need to predict or update
		time_us_ = meas_package.timestamp_;
		is_initialized_ = true;
		return;
	}

	
	

	// Predict
	//compute the time elapsed between the current and previous measurements
	double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;	//expressed in seconds
	time_us_ = meas_package.timestamp_;
	
	Prediction(delta_t);

	//std::cout << P_ << std::endl;
	

	//Update
	if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
		// Radar updates
		UpdateRadar(meas_package);
	}
	else {
		// Laser
		UpdateLidar(meas_package);
	}

	//std::cout << x_ << std::endl;
	//std::cout << P_ << std::endl;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

	// Outline of Prediction:
	//   Generate sigma points
	//   Predict sigma points (and save!)
	//   Calculate predicted mean and covariance




	//
	// Generate Sigma Points
	//

	// Create Augmented mean and covariance
	VectorXd x_aug = VectorXd(7);
	MatrixXd P_aug = MatrixXd(7, 7);

	// Assign augmented mean.  Last 2 (augmented part) are 0 because
	// mean of noise is 0.
	x_aug.head(n_x_) = x_;
	x_aug.tail(2) << 0, 0;

	// Assign augmented covariance
	P_aug.fill(0.0);
	P_aug.topLeftCorner(n_x_, n_x_) = P_;
	P_aug.bottomRightCorner(2, 2) << std_a_*std_a_, 0,
		                               0,             std_yawdd_*std_yawdd_;
	
	//std::cout << P_aug << std::endl;
	
	// Create augmented sigma points
	MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

	// Assign augmented sigma points: First point
	Xsig_aug.col(0) = x_aug.col(0);

	// Assign augmented sigma points: Other points
	MatrixXd A = P_aug.llt().matrixL();


	//std::cout << A << std::endl;

	for (int i = 0; i<n_aug_; i++) {
		VectorXd A_col_scaled = sqrt(lambda_ + n_aug_)*A.col(i);

		Xsig_aug.col(i + 1) = x_aug.col(0) + A_col_scaled;
		Xsig_aug.col(i + 1 + n_aug_) = x_aug.col(0) - A_col_scaled;

	}

	//std::cout << Xsig_aug << std::endl;

	//
	// Predict Sigma Points
	//

	// Create change in x as mean and noise components
	VectorXd delta_x = VectorXd(n_x_);
	VectorXd noise = VectorXd(n_x_);

	// Iterate over sigma points
	for (int i = 0; i<Xsig_aug.cols(); i++) {
		
		// Extract the components of the i_th sigma point		
		double px = Xsig_aug(0,i);
		double py = Xsig_aug(1, i);
		double v = Xsig_aug(2, i);
		double psi = Xsig_aug(3, i);
		double psi_dot = Xsig_aug(4, i);
		double nu_a = Xsig_aug(5, i);
		double nu_psidd = Xsig_aug(6, i);

		// Calculate change in mean motion due to delta_t
		if (psi_dot != 0) {

			delta_x << v / psi_dot*(sin(psi + psi_dot*delta_t) - sin(psi)),
				v / psi_dot*(-cos(psi + psi_dot*delta_t) + cos(psi)),
				0,
				psi_dot*delta_t,
				0;			
		}
		else {
			delta_x << v*cos(psi)*delta_t,
				v*sin(psi)*delta_t,
				0,
				0,
				0;
		}

		// Calculate noise of motion due to delta_t
		noise << (0.5)*delta_t*delta_t*cos(psi)*nu_a,
			(0.5)*delta_t*delta_t*sin(psi)*nu_a,
			delta_t*nu_a,
			(0.5)*delta_t*delta_t*nu_psidd,
			delta_t*nu_psidd;

		// Add mean and noise components to update sigma point
		VectorXd sig_aug = Xsig_aug.col(i).head(n_x_);
		Xsig_pred_.col(i) = sig_aug + delta_x + noise;		
	}

	


	//
	// Calculate predicted mean and covariance
	//

	//predict state mean, overwrite x_ with the updated value
	x_ = Xsig_pred_*weights_;

	//predict state covariance matrix
	MatrixXd D = Xsig_pred_.colwise() - x_;

	//angle normalization
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  

		//while (D(3,i) > M_PI) D(3,i) -= 2.*M_PI;
		//while (D(3,i) < -M_PI) D(3,i) += 2.*M_PI;
	}

	// Multiply weighted matrices for covariance
	MatrixXd Dw = D*weights_.asDiagonal();
	P_ = Dw*D.transpose();

}



/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

	// Outline:
	//   Map predicted sigma points to Lidar space
	//   Update x_ and P_


	// 
	//  Map predicted sigma points to Lidar space
	//

	// Create matrix for sigma points in measurement space
	int n_z = 2;
	MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

	// Map each sigma point to radar space
	for (int i = 0; i < Xsig_pred_.cols(); i++) {

		double px = Xsig_pred_(0, i);
		double py = Xsig_pred_(1, i);
		double v = Xsig_pred_(2, i);
		double psi = Xsig_pred_(3, i);
		double psi_dot = Xsig_pred_(4, i);

		double rho = sqrt(px*px + py*py);

		Zsig.col(i) << px, py;
}

	// Calculate mean predicted measurement
	//VectorXd z_pred = VectorXd(n_z);
	//z_pred.fill(0.0);
	//for (int i = 0; i < Zsig.cols(); i++) {
	//	z_pred = z_pred + weights_(i)*Zsig.col(i);
	//}

	//predict state mean, overwrite x_ with the updated value
	VectorXd z_pred = Zsig * weights_;

	//calculate measurement covariance matrix S

	// noise matrix
	MatrixXd R = MatrixXd(n_z, n_z);

	R << std_laspx_*std_laspx_, 0,
		   0,                     std_laspy_*std_laspy_;

	// Create diff matrix, normalize angles, multiply with transpose
	MatrixXd D = Zsig.colwise() - z_pred;

	// Normalize angles
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {
		//while (Zsig(1, i) > M_PI) Zsig(1, i) -= 2.*M_PI;
		//while (Zsig(1, i) < -M_PI) Zsig(1, i) += 2.*M_PI;
	}
		
	MatrixXd Dw = D * weights_.asDiagonal();
	MatrixXd S = Dw * D.transpose() + R;
		

	//
	//   Update x_ and P_
	//
	//calculate cross correlation matrix
	//create matrix for cross correlation Tc
	MatrixXd Tc = MatrixXd(n_x_, n_z);

	//calculate cross correlation matrix
	D = Xsig_pred_.colwise() - x_;
	Dw = D * weights_.asDiagonal();
	MatrixXd Dz = Zsig.colwise() - z_pred;
	Tc = Dw*Dz.transpose();


	//calculate Kalman gain K;
	//MatrixXd K = MatrixXd(n_x_, n_z);
	MatrixXd K = Tc*S.inverse();

	// Extract lidar state
	double px = meas_package.raw_measurements_[0];
	double py = meas_package.raw_measurements_[1];

	VectorXd z = VectorXd(n_z);
	z << px,py;


	//update state mean and covariance matrix
	x_ = x_ + K*(z - z_pred);
	P_ = P_ - K*S*K.transpose();

	// Measure NIS
	NIS_laser_ = (z - z_pred).transpose()*S.inverse()*(z - z_pred);


}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

	// Outline:
	//   Map predicted sigma points to Radar space
	//   Update x_ and P_


	// 
	//  Map predicted sigma points to Radar space
	//

	// Create matrix for sigma points in measurement space
	int n_z = 3;
	MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
	
	// Map each sigma point to radar space
	for (int i = 0; i < Xsig_pred_.cols(); i++) {

		double px = Xsig_pred_(0, i);
		double py = Xsig_pred_(1, i);
		double v = Xsig_pred_(2, i);
		double psi = Xsig_pred_(3, i);
		double psi_dot = Xsig_pred_(4, i);

		double rho = sqrt(px*px + py*py);

		Zsig.col(i) << rho,
			             atan2(py, px),
			             (px*cos(psi)*v + py*sin(psi)*v) / rho;
	}

	// Calculate mean predicted measurement
	//VectorXd z_pred = VectorXd(n_z);
	//z_pred.fill(0.0);
	//for (int i = 0; i < Zsig.cols(); i++) {
	//	z_pred = z_pred + weights_(i)*Zsig.col(i);
	//}

	//predict state mean, overwrite x_ with the updated value
	VectorXd z_pred = Zsig * weights_;

	//calculate measurement covariance matrix S
	
	// noise matrix
	MatrixXd R = MatrixXd(n_z, n_z);

	R << std_radr_*std_radr_, 0, 0,
		0, std_radphi_*std_radphi_, 0,
		0, 0, std_radrd_*std_radrd_;

	// Create diff matrix, normalize angles, multiply with transpose
	MatrixXd D = Zsig.colwise() - z_pred;

	// Normalize angles
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
																						 //residual
		//angle normalization
		//while (Zsig(1, i) > M_PI) Zsig(1, i) -= 2.*M_PI;
	  //while (Zsig(1, i) < -M_PI) Zsig(1, i) += 2.*M_PI;
	}

	MatrixXd Dw = D*weights_.asDiagonal();
	MatrixXd S = Dw*D.transpose() + R;


	//
	//   Update x_ and P_
	//

	//calculate cross correlation matrix
	//create matrix for cross correlation Tc
	MatrixXd Tc = MatrixXd(n_x_, n_z);
	
	//calculate cross correlation matrix
	D = Xsig_pred_.colwise() - x_;
	Dw = D*weights_.asDiagonal();
	MatrixXd Dz = Zsig.colwise() - z_pred;
	Tc = Dw*Dz.transpose();


	//calculate Kalman gain K;
	//MatrixXd K = MatrixXd(n_x_, n_z);
	MatrixXd K = Tc*S.inverse();

	// Extract radar state
	double rho = meas_package.raw_measurements_[0];
	double phi = meas_package.raw_measurements_[1];
	double rho_dot = meas_package.raw_measurements_[2];

	VectorXd z = VectorXd(n_z);
	z << rho, phi, rho_dot;


	//update state mean and covariance matrix
	x_ = x_ + K*(z - z_pred);
	P_ = P_ - K*S*K.transpose();

	// Measure NIS
	NIS_radar_ = (z - z_pred).transpose()*S.inverse()*(z - z_pred);

	

}
