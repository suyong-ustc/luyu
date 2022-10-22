#pragma once
#include <armadillo>
#include "DICTypes.h"



bool isConvergent1(
	const arma::vec& dp, 
	const double& max_delta_x, 
	const double& max_delta_y,
	const double& err_thresholod = 1e-3);


bool isConvergent2(
	const arma::vec& dp, 
	const double& max_delta_x, 
	const double& max_delta_y,
	const double& err_thresholod = 1e-3);