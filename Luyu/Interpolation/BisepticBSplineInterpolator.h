#pragma once
#include <armadillo>
#include "Interpolator.h"


class BisepticBSplineInterpolator : public Interpolator
{
public:
	BisepticBSplineInterpolator(const arma::mat& data);
	~BisepticBSplineInterpolator();


	double Value(const double& x, const double& y) const;
	arma::mat Values(const arma::mat& x, const arma::mat& y) const;

	double GradientX(const double& x, const double& y) const;
	double GradientY(const double& x, const double& y) const;

	arma::mat GradientsX(const arma::mat& x, const arma::mat& y) const;
	arma::mat GradientsY(const arma::mat& x, const arma::mat& y) const;

private:

	void MOMSMatrix();
	void IterateParameter(double& z1, double& z2, double& z3, double& a1, double& a2, double& a3);
	void ConstructCoefficients(const arma::mat& data);

	arma::mat::fixed<8, 8> bsv_;
	arma::mat::fixed<8, 7> bsg_;

	arma::mat lut_;
};

