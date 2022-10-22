#include "BiquinticBSplineInterpolatror.h"
using namespace arma;


void BiquinticBSplineInterpolatror::IterateParameter(double& z1, double& z2, double& a1, double& a2)
/* Calculate z1, z2 */
{
	z1 = -0.4305753470999737919;
	z2 = -0.04309628820326465382;
	a1 = z1 / (z1 * z1 - 1);
	a2 = z2 / (z2 * z2 - 1);
}


void BiquinticBSplineInterpolatror::ConstructCoefficients(const mat& data)
/* Calculate interpolation coefficients for bi-quintic B-spline */
{
	// Initialize
	double z1, z2, a1, a2;
	IterateParameter(z1, z2, a1, a2);
	const double err = 1e-7;
	const int trunc1 = static_cast<int>(log(err) / log(abs(-z1))) + 1;
	const int trunc2 = static_cast<int>(log(err) / log(abs(-z2))) + 1;

	const int h = data.n_rows;
	const int w = data.n_cols;
	lut_.set_size(h, w);
	mat q(h, w);

	// Unser's iteration for each column
	rowvec bcr(w, fill::zeros);
	for (int r = 0; r <= trunc1; ++r)
		bcr += data.row(r) * pow(z1, r);
	q.row(0) = bcr;
	for (int r = 1; r != h; ++r)
		q.row(r) = data.row(r) + z1 * q.row(r - 1);          // c+(k) = s(k) + z * c+(k-1)

	lut_.row(h - 1) = a1 * (q.row(h - 1) + z1 * q.row(h - 2));
	for (int r = h - 2; r >= 0; --r)
		lut_.row(r) = z1 * (lut_.row(r + 1) - q.row(r));     // c-(k) = z * (c-(k+1) - c+(k))

	bcr.zeros();
	for (int r = 0; r <= trunc2; ++r)
		bcr += lut_.row(r) * pow(z2, r);
	q.row(0) = bcr;
	for (int r = 1; r != h; ++r)
		q.row(r) = lut_.row(r) + z2 * q.row(r - 1);          // c+(k) = s(k) + z * c+(k-1)

	lut_.row(h - 1) = a2 * (q.row(h - 1) + z2 * q.row(h - 2));
	for (int r = h - 2; r >= 0; --r)
		lut_.row(r) = z2 * (lut_.row(r + 1) - q.row(r));     // c-(k) = z * (c-(k+1) - c+(k))


	// Unser's iteration for each row
	vec bcc(h, fill::zeros);
	for (int c = 0; c <= trunc1; ++c)
		bcc += lut_.col(c) * pow(z1, c);
	q.col(0) = bcc;
	for (int c = 1; c != w; ++c)
		q.col(c) = lut_.col(c) + z1 * q.col(c - 1);          // c+(k) = s(k) + z * c+(k-1)

	lut_.col(w - 1) = a1 * (q.col(w - 1) + z1 * q.col(w - 2));
	for (int c = w - 2; c >= 0; --c)
		lut_.col(c) = z1 * (lut_.col(c + 1) - q.col(c));     // c-(k) = z * (c-(k+1) - c+(k))

	bcc.zeros();
	for (int c = 0; c <= trunc2; ++c)
		bcc += lut_.col(c) * pow(z2, c);
	q.col(0) = bcc;
	for (int c = 1; c != w; ++c)
		q.col(c) = lut_.col(c) + z2 * q.col(c - 1);          // c+(k) = s(k) + z * c+(k-1)

	lut_.col(w - 1) = a2 * (q.col(w - 1) + z2 * q.col(w - 2));
	for (int c = w - 2; c >= 0; --c)
		lut_.col(c) = z2 * (lut_.col(c + 1) - q.col(c));     // c-(k) = z * (c-(k+1) - c+(k))

}


void BiquinticBSplineInterpolatror::MOMSMatrix()
/* Coefficient matrix for interpolating intensity and gradient */
{
	// Coefficient matrix for interpolating intensity
	bsv_ << 1 << -5 << 10 << -10 << 5 << -1 << endr
		<< 26 << -50 << 20 << 20 << -20 << 5 << endr
		<< 66 << 0 << -60 << 0 << 30 << -10 << endr
		<< 26 << 50 << 20 << -20 << -20 << 10 << endr
		<< 1 << 5 << 10 << 10 << 5 << -5 << endr
		<< 0 << 0 << 0 << 0 << 0 << 1 << endr;

	// Coefficient matrix for interpolating gradient
	bsg_ << -1 << 4 << -6 << 4 << -1 << endr
		<< -10 << 8 << 12 << -16 << 5 << endr
		<< 0 << -24 << 0 << 24 << -10 << endr
		<< 10 << 8 << -12 << -16 << 10 << endr
		<< 1 << 4 << 6 << 4 << -5 << endr
		<< 0 << 0 << 0 << 0 << 1 << endr;
	bsg_ *= 5;
}



BiquinticBSplineInterpolatror::BiquinticBSplineInterpolatror(const mat& data)
{
	MOMSMatrix();
	ConstructCoefficients(data);
}


BiquinticBSplineInterpolatror::~BiquinticBSplineInterpolatror()
{
}


double BiquinticBSplineInterpolatror::Value(const double& x, const double& y) const
{
	const int c = static_cast<int>(floor(x));
	const int r = static_cast<int>(floor(y));

	const double ex = x - c;
	const double ex2 = ex * ex;
	const double ex3 = ex2 * ex;
	const double ex4 = ex3 * ex;
	const double ex5 = ex4 * ex;

	const double ey = y - r;
	const double ey2 = ey * ey;
	const double ey3 = ey2 * ey;
	const double ey4 = ey3 * ey;
	const double ey5 = ey4 * ey;

	double dx[] = { 1, ex, ex2, ex3, ex4, ex5 };
	double dy[] = { 1, ey, ey2, ey3, ey4, ey5 };
	vec::fixed<6> deltax(dx);
	rowvec::fixed<6> deltay(dy);
	const mat coeff = lut_.submat(r - 2, c - 2, r + 3, c + 3);
	const double value = as_scalar(deltay * bsv_.t() * coeff * bsv_ * deltax);

	return value;
}


mat BiquinticBSplineInterpolatror::Values(const mat& x, const mat& y) const
{
	const int h = x.n_rows;
	const int w = x.n_cols;

	mat value(h, w);
	for (int c = 0; c < w; ++c)
		for (int r = 0; r < h; ++r)
			value(r, c) = Value(x(r, c), y(r, c));

	return value;
}


double BiquinticBSplineInterpolatror::GradientX(const double& x, const double& y) const
{
	const int c = static_cast<int>(floor(x));
	const int r = static_cast<int>(floor(y));

	const double ex = x - c;
	const double ex2 = ex * ex;
	const double ex3 = ex2 * ex;
	const double ex4 = ex3 * ex;

	const double ey = y - r;
	const double ey2 = ey * ey;
	const double ey3 = ey2 * ey;
	const double ey4 = ey3 * ey;
	const double ey5 = ey4 * ey;

	double dx[] = { 1, ex, ex2, ex3, ex4 };
	double dy[] = { 1, ey, ey2, ey3, ey4, ey5 };
	vec::fixed<5> deltax(dx);
	rowvec::fixed<6> deltay(dy);
	const mat coeff = lut_.submat(r - 2, c - 2, r + 3, c + 3);
	const double value = as_scalar(deltay * bsv_.t() * coeff * bsg_ * deltax);
	return value;
}


double BiquinticBSplineInterpolatror::GradientY(const double& x, const double& y) const
{
	const int c = static_cast<int>(floor(x));
	const int r = static_cast<int>(floor(y));

	const double ex = x - c;
	const double ex2 = ex * ex;
	const double ex3 = ex2 * ex;
	const double ex4 = ex3 * ex;
	const double ex5 = ex4 * ex;

	const double ey = y - r;
	const double ey2 = ey * ey;
	const double ey3 = ey2 * ey;
	const double ey4 = ey3 * ey;

	double dx[] = { 1, ex, ex2, ex3, ex4, ex5 };
	double dy[] = { 1, ey, ey2, ey3, ey4 };
	vec::fixed<6> deltax(dx);
	rowvec::fixed<5> deltay(dy);
	const mat coeff = lut_.submat(r - 2, c - 2, r + 3, c + 3);
	const double value = as_scalar(deltay * bsg_.t() * coeff * bsv_ * deltax);
	return value;
}


mat BiquinticBSplineInterpolatror::GradientsX(const mat& x, const mat& y) const
{
	const int h = x.n_rows;
	const int w = x.n_cols;

	mat gx(h, w);
	for (int c = 0; c < w; ++c)
		for (int r = 0; r < h; ++r)
			gx(r, c) = GradientX(x(r, c), y(r, c));

	return gx;
}


mat BiquinticBSplineInterpolatror::GradientsY(const mat& x, const mat& y) const
{
	const int h = x.n_rows;
	const int w = x.n_cols;

	mat gy(h, w);
	for (int c = 0; c < w; ++c)
		for (int r = 0; r < h; ++r)
			gy(r, c) = GradientY(x(r, c), y(r, c));

	return gy;
}
