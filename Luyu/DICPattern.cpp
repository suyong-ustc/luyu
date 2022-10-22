#include "DICPattern.h"
#include <iostream>
#include <fstream>
#include <vector>
using namespace arma;



DICPattern::DICPattern(const DICPattern& pattern, const double& u0, const double& v0)
	: radius_(pattern.radius_)
{
	x_ = pattern.x_ + u0;
	y_ = pattern.y_ + v0;
}



DICPattern::DICPattern(const double& radius, const vec& x, const vec& y)
	: radius_(radius), x_(x), y_(y)
{

}



DICPattern::DICPattern(const double& radius, const int& n, const double& a, const double& b)
	: radius_(radius)
{
	const int h = 2 * n + 1;
	const int w = 2 * n + 1;

	// firstly, put the speckle at the grid location
	mat speckle_position_x = zeros<mat>(h, w);
	mat speckle_position_y = zeros<mat>(h, w);

	for (int r = -n; r <= n; ++r)
	{
		for (int c = -n; c <= n; ++c)
		{
			speckle_position_x(r + n, c + n) = c * a;
			speckle_position_y(r + n, c + n) = r * a;
		}
	}

	// add random offset
	const mat x_offset = 2 * b * (randu(h, w) - 0.5);
	const mat y_offset = 2 * b * (randu(h, w) - 0.5);

	speckle_position_x = speckle_position_x + x_offset;
	speckle_position_y = speckle_position_y + y_offset;

	x_ = vectorise(speckle_position_x);
	y_ = vectorise(speckle_position_y);
}



DICPattern::DICPattern(const std::string & file_name)
{
	Read(file_name);
}



DICPattern::~DICPattern()
{
}



/***********************************************************************************************/
/*******************          Image intensity and gradients         ****************************/
/***********************************************************************************************/


double DICPattern::Value(const double& x, const double& y) const
{
	const vec dx = x_ - x;
	const vec dy = y_ - y;
	const vec g = exp(-(dx%dx + dy % dy) / (radius_*radius_));
	return accu(g);
}



double DICPattern::GradientX(const double& x, const double& y) const
{
	const vec dx = x_ - x;
	const vec dy = y_ - y;
	const double c = 2 / (radius_*radius_);
	const vec gx = c * dx % exp(-(dx%dx + dy % dy) / (radius_*radius_));
	return accu(gx);
}



double DICPattern::GradientY(const double& x, const double& y) const
{
	const vec dx = x_ - x;
	const vec dy = y_ - y;
	const double c = 2  / (radius_*radius_);
	const vec gy = c * dy % exp(-(dx%dx + dy % dy) / (radius_*radius_));
	return accu(gy);
}



mat DICPattern::Values(const mat& x, const mat& y) const
{
	const int rows = x.n_rows;
	const int cols = x.n_cols;

	mat values(rows, cols);
	for (int c = 0; c < cols; ++c)
		for (int r = 0; r < rows; ++r)
			values(r, c) = Value(x(r, c), y(r, c));

	return values;
}



mat DICPattern::GradientsX(const mat& x, const mat& y) const
{
	const int rows = x.n_rows;
	const int cols = x.n_cols;

	mat gx(rows, cols);
	for (int c = 0; c < cols; ++c)
		for (int r = 0; r < rows; ++r)
			gx(r, c) = GradientX(x(r, c), y(r, c));

	return gx;
}



mat DICPattern::GradientsY(const mat& x, const mat& y) const
{
	const int rows = x.n_rows;
	const int cols = x.n_cols;

	mat gy(rows, cols);
	for (int c = 0; c < cols; ++c)
		for (int r = 0; r < rows; ++r)
			gy(r, c) = GradientY(x(r, c), y(r, c));

	return gy;
}



mat DICPattern::Sample(const double& x0, const double& y0, const int& half_image_size) const
{
	const int h = 2 * half_image_size + 1;
	const int w = 2 * half_image_size + 1;

	mat g(h, w);
	for (int c = 0; c < w; ++c)
	{
		const double x = x0 + c - half_image_size;

		for (int r = 0; r < h; ++r)
		{
			const double y = y0 + r - half_image_size;

			g(r, c) = Value(x, y);
		}
			
	}

	return g;
}


/***********************************************************************************************/
/*********************          Read and write pattern information       ***********************/
/***********************************************************************************************/


bool DICPattern::Write(const std::string& file_name) const
{
	std::ofstream outfile;
	outfile.open(file_name.c_str(), std::ofstream::out);
	if (!outfile) {
		std::cerr << "error : unable to open output file : " << file_name << std::endl;
		return false;
	}

	outfile.precision(20);
	outfile << radius_ << '\n' << std::endl;
	x_.save(outfile, arma_ascii);
	outfile << std::endl;
	y_.save(outfile, arma_ascii);
	outfile.close();

	return true;
}


bool DICPattern::Read(const std::string file_name)
{
	std::ifstream infile;
	infile.open(file_name.c_str());
	if (!infile) {
		std::cerr << "error : unable to open input file : " << file_name << std::endl;
		return false;
	}

	infile >> radius_;
	x_.load(infile, arma_ascii);
	y_.load(infile, arma_ascii);
	infile.close();

	return true;
}