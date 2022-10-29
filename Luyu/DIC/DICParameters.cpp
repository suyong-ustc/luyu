#include "DICParameters.h"
using namespace arma;


DICParameters::DICParameters()
{
	roi_.setRect(50, 50, 400, 400);

	grid_step_ = 1;

	subset_size_ = 19;

	zncc_threshold_ = 0.8;

	max_iteration_times_ = 10;

	error_threshold_ = 2e-4;

	seed_.setX(250);
	seed_.setY(250);

	search_radius_ = 15;

	bspline_interpolation_order_ = 5;

	shape_function_order_ = 1;
}



DICParameters::~DICParameters()
{

}



bool DICParameters::grid(mat& x, mat& y) const
{
	// 感兴趣区域
	const int xmin = roi_.left();
	const int xmax = roi_.right();
	const int ymin = roi_.top();
	const int ymax = roi_.bottom();

	if (xmin >= xmax || ymin >= ymax)
	{
		std::cerr << "Invalid ROI!" << std::endl;
		return false;
	}

	// 网格
	const int nx = 1 + (xmax - xmin) / grid_step_;
	const int ny = 1 + (ymax - ymin) / grid_step_;

	x.zeros(ny, nx);
	y.zeros(ny, nx);

	for (int r = 0; r < ny; ++r)
	{
		for (int c = 0; c < nx; ++c)
		{
			x(r, c) = xmin + c * grid_step_;
			y(r, c) = ymin + r * grid_step_;
		}

	}

	return true;
}
