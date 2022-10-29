#include "BicubicMOMSInterpolator.h"
using namespace arma;


void BicubicMOMSInterpolator::IterateParameter(double& z, double& a)
/* Calculate z */
{
	const double b = 4 * (1 - 3 * lamda_) / (1 + 6 * lamda_);
	if (abs(b) < 2) {
		std::cerr << "Invalid MOMS coefficients!" << std::endl;
		return;
	}
	const double z1 = 0.5 * (-b + sqrt(b*b - 4));
	const double z2 = 0.5 * (-b - sqrt(b*b - 4));
	z = abs(z1) < 1 ? z1 : z2;
	a = z / (z * z - 1);
}


void BicubicMOMSInterpolator::ConstructCoefficients(const arma::mat& data)
/* Calculate interpolation coefficients */
{
	// Initialize
	double z, a;
	IterateParameter(z, a);
	const double err = 1e-7;
	const int trunc = static_cast<int>(log(err) / log(abs(-z))) + 1;

	const int h = data.n_rows;
	const int w = data.n_cols;
	lut_.set_size(h, w);
	mat q(h, w);

	// Unser's iteration for each column
	rowvec bcr(w, fill::zeros);
	for (int r = 0; r <= trunc; ++r)
		bcr += data.row(r) * pow(z, r);
	q.row(0) = bcr;
	for (int r = 1; r != h; ++r)
		q.row(r) = data.row(r) + z * q.row(r - 1);          // c+(k) = s(k) + z * c+(k-1)

	lut_.row(h - 1) = a * (q.row(h - 1) + z * q.row(h - 2));
	for (int r = h - 2; r >= 0; --r)
		lut_.row(r) = z * (lut_.row(r + 1) - q.row(r));     // c-(k) = z * (c-(k+1) - c+(k))

	// Unser's iteration for each row
	vec bcc(h, fill::zeros);
	for (int c = 0; c <= trunc; ++c)
		bcc += lut_.col(c) * pow(z, c);
	q.col(0) = bcc;
	for (int c = 1; c != w; ++c)
		q.col(c) = lut_.col(c) + z * q.col(c - 1);          // c+(k) = s(k) + z * c+(k-1)

	lut_.col(w - 1) = a * (q.col(w - 1) + z * q.col(w - 2));
	for (int c = w - 2; c >= 0; --c)
		lut_.col(c) = z * (lut_.col(c + 1) - q.col(c));     // c-(k) = z * (c-(k+1) - c+(k))
}


void BicubicMOMSInterpolator::MOMSMatrix()
/* Coefficient matrix for interpolating intensity and gradient */
{
	// Coefficient matrix for interpolating intensity
	bsv_ << 1 + 6 * lamda_  << -3 - 6 * lamda_ <<  3 << -1 << endr
		 << 4 - 12 * lamda_ << 18 * lamda_     << -6 <<  3 << endr
		 << 1 + 6 * lamda_  << 3 - 18 * lamda_ <<  3 << -3 << endr
		 << 0               << 6 * lamda_      <<  0 <<  1 << endr;
	bsv_ /= (1 + 6 * lamda_);

	// Coefficient matrix for interpolating gradient
	bsg_ << -1 - 2 * lamda_ << 2  << -1 << endr
		<< 6 * lamda_       << -4 <<  3 << endr
		<< 1 - 6 * lamda_   << 2  << -3 << endr
		<< 2 * lamda_       << 0  <<  1 << endr;
	bsg_ *= 3 / (1 + 6 * lamda_);
}


BicubicMOMSInterpolator::BicubicMOMSInterpolator(const arma::mat& data, const double& lamda) 
	: lamda_(lamda)
{
	MOMSMatrix();
	ConstructCoefficients(data);
}


BicubicMOMSInterpolator::~BicubicMOMSInterpolator()
{
}


double BicubicMOMSInterpolator::Value(const double& x, const double& y) const
{
	const int c = static_cast<int>(floor(x));
	const int r = static_cast<int>(floor(y));
	const double ex = x - c;
	const double ey = y - r;

	double dx[] = { 1, ex, ex*ex, ex*ex*ex };
	double dy[] = { 1, ey, ey*ey, ey*ey*ey };
	vec::fixed<4> deltax(dx);
	rowvec::fixed<4> deltay(dy);
	const mat coeff = lut_.submat(r - 1, c - 1, r + 2, c + 2);
	const double value = as_scalar(deltay * bsv_.t() * coeff * bsv_ * deltax);

	return value;
}


mat BicubicMOMSInterpolator::Values(const mat& x, const mat& y) const
{
	const int h = x.n_rows;
	const int w = x.n_cols;

	mat value(h, w);
	for (int c = 0; c < w; ++c)
		for (int r = 0; r < h; ++r)
			value(r, c) = Value(x(r, c), y(r, c));

	return value;
}


double BicubicMOMSInterpolator::GradientX(const double& x, const double& y) const
{
	const int c = static_cast<int>(floor(x));
	const int r = static_cast<int>(floor(y));
	const double ex = x - c;
	const double ey = y - r;

	double dx[] = { 1, ex, ex*ex };
	double dy[] = { 1, ey, ey*ey, ey*ey*ey };
	vec::fixed<3> deltax(dx);
	rowvec::fixed<4> deltay(dy);
	const mat coeff = lut_.submat(r - 1, c - 1, r + 2, c + 2);
	const double value = as_scalar(deltay * bsv_.t() * coeff * bsg_ * deltax);
	return value;
}


double BicubicMOMSInterpolator::GradientY(const double& x, const double& y) const
{
	const int c = static_cast<int>(floor(x));
	const int r = static_cast<int>(floor(y));
	const double ex = x - c;
	const double ey = y - r;

	double dx[] = { 1, ex, ex*ex, ex*ex*ex };
	double dy[] = { 1, ey, ey*ey };
	vec::fixed<4> deltax(dx);
	rowvec::fixed<3> deltay(dy);
	const mat coeff = lut_.submat(r - 1, c - 1, r + 2, c + 2);
	const double value = as_scalar(deltay * bsg_.t() * coeff * bsv_ * deltax);
	return value;
}


mat BicubicMOMSInterpolator::GradientsX(const mat& x, const mat& y) const
{
	const int h = x.n_rows;
	const int w = x.n_cols;

	mat gx(h, w);
	for (int c = 0; c < w; ++c)
		for (int r = 0; r < h; ++r)
			gx(r, c) = GradientX(x(r, c), y(r, c));

	return gx;
}


mat BicubicMOMSInterpolator::GradientsY(const mat& x, const mat& y) const
{
	const int h = x.n_rows;
	const int w = x.n_cols;

	mat gy(h, w);
	for (int c = 0; c < w; ++c)
		for (int r = 0; r < h; ++r)
			gy(r, c) = GradientY(x(r, c), y(r, c));

	return gy;
}