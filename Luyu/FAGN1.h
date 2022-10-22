#pragma once
#include <armadillo>

class Interpolator;


bool ForwardAdditiveGaussNewtonMethod1(
	const arma::mat& refer_intensities, 
	const double refer_x0, 
	const double refer_y0, 
	Interpolator* deform_image_interp,
	const int half_subset_size, 
	const int max_iter_time, 
	const double zncc_threshold,
	arma::vec& deform_warp, 
	arma::uword& iter, 
	double& zncc);


static bool ForwardPseudoMatrix1(
	const arma::mat& deform_gradient_x,
	const arma::mat& deform_gradient_y,
	const arma::mat& delta_x,
	const arma::mat& delta_y,
	arma::mat& pseudo);