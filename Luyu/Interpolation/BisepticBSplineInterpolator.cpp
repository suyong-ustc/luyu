#include "BisepticBSplineInterpolator.h"
#define PI arma::datum::pi
using namespace arma;


void BisepticBSplineInterpolator::IterateParameter(double& z1, double& z2, double& z3, double& a1, double& a2, double& a3)
/* Calculate z1, z2, z3 */
{
	z1 = -0.5352804307964381655;
	z2 = -0.12255461519232669052;
	z3 = -0.009148694809608276929;
	a1 = z1 / (z1 * z1 - 1);
	a2 = z2 / (z2 * z2 - 1);
	a3 = z3 / (z3 * z3 - 1);
}


void BisepticBSplineInterpolator::ConstructCoefficients(const mat& data)
/* Calculate interpolation coefficients for bi-septic B-spline */
{
	// Initialize
	double z1, z2, z3, a1, a2, a3;
	IterateParameter(z1, z2, z3, a1, a2, a3);
	const double err = 1e-7;
	const int trunc1 = static_cast<int>(log(err) / log(abs(-z1))) + 1;
	const int trunc2 = static_cast<int>(log(err) / log(abs(-z2))) + 1;
	const int trunc3 = static_cast<int>(log(err) / log(abs(-z3))) + 1;

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

	bcr.zeros();
	for (int r = 0; r <= trunc3; ++r)
		bcr += lut_.row(r) * pow(z3, r);
	q.row(0) = bcr;
	for (int r = 1; r != h; ++r)
		q.row(r) = lut_.row(r) + z3 * q.row(r - 1);          // c+(k) = s(k) + z * c+(k-1)

	lut_.row(h - 1) = a3 * (q.row(h - 1) + z3 * q.row(h - 2));
	for (int r = h - 2; r >= 0; --r)
		lut_.row(r) = z3 * (lut_.row(r + 1) - q.row(r));     // c-(k) = z * (c-(k+1) - c+(k))

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

	bcc.zeros();
	for (int c = 0; c <= trunc3; ++c)
		bcc += lut_.col(c) * pow(z3, c);
	q.col(0) = bcc;
	for (int c = 1; c != w; ++c)
		q.col(c) = lut_.col(c) + z3 * q.col(c - 1);          // c+(k) = s(k) + z * c+(k-1)

	lut_.col(w - 1) = a3 * (q.col(w - 1) + z3 * q.col(w - 2));
	for (int c = w - 2; c >= 0; --c)
		lut_.col(c) = z3 * (lut_.col(c + 1) - q.col(c));     // c-(k) = z * (c-(k+1) - c+(k))

}


void BisepticBSplineInterpolator::MOMSMatrix()
/* Coefficient matrix for interpolating intensity and gradient */
{
	// Coefficient matrix for interpolating intensity
	bsv_ << 1 << -7 << 21 << -35 << 35 << -21 << 7 << -1 << endr
		<< 120 << -392 << 504 << -280 << 0 << 84 << -42 << 7 << endr
		<< 1191 << -1715 << 315 << 665 << -315 << -105 << 105 << -21 << endr
		<< 2416 << 0 << -1680 << 0 << 560 << 0 << -140 << 35 << endr
		<< 1191 << 1715 << 315 << -665 << -315 << 105 << 105 << -35 << endr
		<< 120 << 392 << 504 << 280 << 0 << -84 << -42 << 21 << endr
		<< 1 << 7 << 21 << 35 << 35 << 21 << 7 << -7 << endr
		<< 0 << 0 << 0 << 0 << 0 << 0 << 0 << 1 << endr;

	// Coefficient matrix for interpolating gradient
	bsg_ << -1 << 6 << -15 << 20 << -15 << 6 << -1 << endr
		<< -56 << 144 << -120 << 0 << 60 << -36 << 7 << endr
		<< -245 << 90 << 285 << -180 << -75 << 90 << -21 << endr
		<< 0 << -480 << 0 << 320 << 0 << -120 << 35 << endr
		<< 245 << 90 << -285 << -180 << 75 << 90 << -35 << endr
		<< 56 << 144 << 120 << 0 << -60 << -36 << 21 << endr
		<< 1 << 6 << 15 << 20 << 15 << 6 << -7 << endr
		<< 0 << 0 << 0 << 0 << 0 << 0 << 1;
	bsg_ *= 7;
}


BisepticBSplineInterpolator::BisepticBSplineInterpolator(const mat& data)
{
	MOMSMatrix();
	ConstructCoefficients(data);
}


BisepticBSplineInterpolator::~BisepticBSplineInterpolator()
{
}


double BisepticBSplineInterpolator::Value(const double& x, const double& y) const
{
	const int c = static_cast<int>(floor(x));
	const int r = static_cast<int>(floor(y));

	const double ex = x - c;
	const double ex2 = ex * ex;
	const double ex3 = ex2 * ex;
	const double ex4 = ex3 * ex;
	const double ex5 = ex4 * ex;
	const double ex6 = ex5 * ex;
	const double ex7 = ex6 * ex;

	const double ey = y - r;
	const double ey2 = ey * ey;
	const double ey3 = ey2 * ey;
	const double ey4 = ey3 * ey;
	const double ey5 = ey4 * ey;
	const double ey6 = ey5 * ey;
	const double ey7 = ey6 * ey;

	double dx[] = { 1, ex, ex2, ex3, ex4, ex5, ex6, ex7 };
	double dy[] = { 1, ey, ey2, ey3, ey4, ey5, ey6, ey7 };
	vec::fixed<8> deltax(dx);
	rowvec::fixed<8> deltay(dy);
	const mat coeff = lut_.submat(r - 3, c - 3, r + 4, c + 4);
	const double value = as_scalar(deltay * bsv_.t() * coeff * bsv_ * deltax);

	return value;
}


mat BisepticBSplineInterpolator::Values(const mat& x, const mat& y) const
{
	const int h = x.n_rows;
	const int w = x.n_cols;

	mat value(h, w);
	for (int c = 0; c < w; ++c)
		for (int r = 0; r < h; ++r)
			value(r, c) = Value(x(r, c), y(r, c));

	return value;
}


double BisepticBSplineInterpolator::GradientX(const double& x, const double& y) const
{
	const int c = static_cast<int>(floor(x));
	const int r = static_cast<int>(floor(y));

	const double ex = x - c;
	const double ex2 = ex * ex;
	const double ex3 = ex2 * ex;
	const double ex4 = ex3 * ex;
	const double ex5 = ex4 * ex;
	const double ex6 = ex5 * ex;

	const double ey = y - r;
	const double ey2 = ey * ey;
	const double ey3 = ey2 * ey;
	const double ey4 = ey3 * ey;
	const double ey5 = ey4 * ey;
	const double ey6 = ey5 * ey;
	const double ey7 = ey6 * ey;

	double dx[] = { 1, ex, ex2, ex3, ex4, ex5, ex6 };
	double dy[] = { 1, ey, ey2, ey3, ey4, ey5, ey6, ey7 };
	vec::fixed<7> deltax(dx);
	rowvec::fixed<8> deltay(dy);
	const mat coeff = lut_.submat(r - 3, c - 3, r + 4, c + 4);
	const double value = as_scalar(deltay * bsv_.t() * coeff * bsg_ * deltax);
	return value;
}


double BisepticBSplineInterpolator::GradientY(const double& x, const double& y) const
{
	const int c = static_cast<int>(floor(x));
	const int r = static_cast<int>(floor(y));

	const double ex = x - c;
	const double ex2 = ex * ex;
	const double ex3 = ex2 * ex;
	const double ex4 = ex3 * ex;
	const double ex5 = ex4 * ex;
	const double ex6 = ex5 * ex;
	const double ex7 = ex6 * ex;

	const double ey = y - r;
	const double ey2 = ey * ey;
	const double ey3 = ey2 * ey;
	const double ey4 = ey3 * ey;
	const double ey5 = ey4 * ey;
	const double ey6 = ey5 * ey;

	double dx[] = { 1, ex, ex2, ex3, ex4, ex5, ex6, ex7 };
	double dy[] = { 1, ey, ey2, ey3, ey4, ey5, ey6 };
	vec::fixed<8> deltax(dx);
	rowvec::fixed<7> deltay(dy);
	const mat coeff = lut_.submat(r - 3, c - 3, r + 4, c + 4);
	const double value = as_scalar(deltay * bsg_.t() * coeff * bsv_ * deltax);
	return value;
}


mat BisepticBSplineInterpolator::GradientsX(const mat& x, const mat& y) const
{
	const int h = x.n_rows;
	const int w = x.n_cols;

	mat gx(h, w);
	for (int c = 0; c < w; ++c)
		for (int r = 0; r < h; ++r)
			gx(r, c) = GradientX(x(r, c), y(r, c));

	return gx;
}


mat BisepticBSplineInterpolator::GradientsY(const mat& x, const mat& y) const
{
	const int h = x.n_rows;
	const int w = x.n_cols;

	mat gy(h, w);
	for (int c = 0; c < w; ++c)
		for (int r = 0; r < h; ++r)
			gy(r, c) = GradientY(x(r, c), y(r, c));

	return gy;
}