#pragma once
#include <armadillo>
#include "Interpolator.h"


class BiquinticBSplineInterpolatror : public Interpolator
{
public:
	BiquinticBSplineInterpolatror(const arma::mat& image);
	~BiquinticBSplineInterpolatror();

	double Value(const double& x, const double& y) const;
	arma::mat Values(const arma::mat& x, const arma::mat& y) const;

	double GradientX(const double& x, const double& y) const;
	double GradientY(const double& x, const double& y) const;

	arma::mat GradientsX(const arma::mat& x, const arma::mat& y) const;
	arma::mat GradientsY(const arma::mat& x, const arma::mat& y) const;

public:

	void MOMSMatrix();
	void IterateParameter(double& z1, double& z2, double& a1, double& a2);
	void ConstructCoefficients(const arma::mat& data);

	arma::mat::fixed<6, 6> bsv_;
	arma::mat::fixed<6, 5> bsg_;

	arma::mat lut_;
};

