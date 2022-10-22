#pragma once
#include <armadillo>



/*******************************************************************************************************************/
/**********************************               形函数                ********************************************/
/******************************************************************************************************************/

bool ShapeFunction(const double& x0, const double& y0, const arma::mat& x, const arma::mat& y, const arma::vec& p, const int& order, arma::mat& deform_x, arma::mat& deform_y);

bool ZeroOrderShapeFunction(const double& x0, const double& y0, const arma::mat& x, const arma::mat& y, const arma::vec& p, arma::mat& deform_x, arma::mat& deform_y);

bool FirstOrderShapeFunction(const double& x0, const double& y0, const arma::mat& x, const arma::mat& y, const arma::vec& p, arma::mat& deform_x, arma::mat& deform_y);

bool SecondOrderShapeFunction(const double& x0, const double& y0, const arma::mat& x, const arma::mat& y, const arma::vec& p, arma::mat& deform_x, arma::mat& deform_y);

bool ThirdOrderShapeFunction(const double& x0, const double& y0, const arma::mat& x, const arma::mat& y, const arma::vec& p, arma::mat& deform_x, arma::mat& deform_y);

bool FourthOrderShapeFunction(const double& x0, const double& y0, const arma::mat& x, const arma::mat& y, const arma::vec& p, arma::mat& deform_x, arma::mat& deform_y);

bool FifthOrderShapeFunction(const double& x0, const double& y0, const arma::mat& x, const arma::mat& y, const arma::vec& p, arma::mat& deform_x, arma::mat& deform_y);



/*******************************************************************************************************************/
/**********************************            高斯牛顿迭代              ********************************************/
/******************************************************************************************************************/

class Interpolator; 

bool ZeroOrderGaussNewtonOptimization(const arma::mat& refer_intensities, const double x0, const double y0, Interpolator* deform_image_interp, arma::vec& p, double& zncc);

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


bool ForwardAdditiveGaussNewtonMethod2(
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



static bool ForwardPseudoMatrix2(
	const arma::mat& deform_gradient_x,
	const arma::mat& deform_gradient_y,
	const arma::mat& delta_x,
	const arma::mat& delta_y,
	arma::mat& pseudo);
