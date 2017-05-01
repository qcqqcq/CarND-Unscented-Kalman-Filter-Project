#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */

	// Initialize RMSE as vector of 4 numbers
	// x,y,vx,vy
	VectorXd rmse(4);
	rmse << 0, 0, 0, 0;
	size_t N = estimations.size();

	// Check that estimation size is not 0
	if (N == 0) {
		std::cout << "error: estimation size is 0";
		return rmse;
	}

	// Ensure estimation size is equal to ground truth size
	if (N != ground_truth.size()) {
		std::cout << "error: esimation size not equal to ground truth size";
		return rmse;
	}

	// Calculate residual, square it, then accumulate sum
	for (size_t i = 0; i < N; ++i) {
		VectorXd est = estimations[i];
		VectorXd gnd = ground_truth[i];

		VectorXd res = est - gnd;
		VectorXd res_sq = res.array()*res.array();

		rmse += res_sq;
	}

	// Get mean of sum of squared residuals
	rmse = rmse / N;

	// Take square root
	rmse = rmse.array().sqrt();

	return rmse;

}
