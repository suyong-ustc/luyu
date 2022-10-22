#include "DICConvergence.h"
using namespace arma;



bool isConvergent1(
	const vec& dp, 
	const double& max_delta_x, 
	const double& max_delta_y,
	const double& err_thresholod)
{
	// first-order shape function parameters
	const double du = dp(0);
	const double dux = dp(1);
	const double duy = dp(2);

	const double dv = dp(3);
	const double dvx = dp(4);
	const double dvy = dp(5);

	// convergence criterion
	const double err = du * du +
		max_delta_x * dux * max_delta_x * dux +
		max_delta_y * duy * max_delta_y * duy +
		dv * dv +
		max_delta_x * dvx * max_delta_x * dvx +
		max_delta_y * dvy * max_delta_y * dvy;

	bool is_convergent = false;
	if (err < err_thresholod * err_thresholod)
		is_convergent = true;

	return is_convergent;

	//// first-order shape function parameters
	//const double du = dp(0);
	//const double dv = dp(3);

	//// convergence criterion
	//const double err = 1e-6;

	//bool is_convergent = false;
	//if (abs(du) < err && abs(dv) < err)
	//	is_convergent = true;

	//return is_convergent;
}



bool isConvergent2(
	const vec& dp, 
	const double& max_delta_x, 
	const double& max_delta_y,
	const double& err_thresholod)
{
	// second-order shape function parameters
	const double du = dp(0);
	const double dux = dp(1);
	const double duy = dp(2);
	const double duxx = dp(3);
	const double duxy = dp(4);
	const double duyy = dp(5);

	const double dv = dp(6);
	const double dvx = dp(7);
	const double dvy = dp(8);
	const double dvxx = dp(9);
	const double dvxy = dp(10);
	const double dvyy = dp(11);

	// convergent condition
	const double a1 = du;
	const double a2 = dux * max_delta_x;
	const double a3 = duy * max_delta_y;
	const double a4 = 0.5 * duxx * max_delta_x * max_delta_x;
	const double a5 = duxy * max_delta_x * max_delta_y;
	const double a6 = 0.5 * duyy * max_delta_y * max_delta_y;

	const double b1 = dv;
	const double b2 = dvx * max_delta_x;
	const double b3 = dvy * max_delta_y;
	const double b4 = 0.5 * dvxx * max_delta_x * max_delta_x;
	const double b5 = dvxy * max_delta_x * max_delta_y;
	const double b6 = 0.5 * dvyy * max_delta_y * max_delta_y;

	const double err = a1 * a1 + a2 * a2 + a3 * a3 + a4 * a4 + a5 * a5 + a6 * a6 +
		b1 * b1 + b2 * b2 + b3 * b3 + b4 * b4 + b5 * b5 + b6 * b6;

	bool is_convergent = false;
	if (err < err_thresholod * err_thresholod)
		is_convergent = true;

	return is_convergent;


	//// second-order shape function parameters
	//const double du = dp(0);
	//const double dv = dp(6);

	//// convergent condition
	//const double err = 1e-4;

	//bool is_convergent = false;
	//if (abs(du) < err && abs(dv) < err)
	//	is_convergent = true;

	//return is_convergent;
}