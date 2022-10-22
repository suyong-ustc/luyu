#pragma once
#include <armadillo>
#include "Interpolator.h"


class BicubicMOMSInterpolator : public Interpolator
{
public:
	BicubicMOMSInterpolator(const arma::mat& data, const double& lamda = 0);
	~BicubicMOMSInterpolator();

	double Value(const double& x, const double& y) const;
	arma::mat Values(const arma::mat& x, const arma::mat& y) const;

	double GradientX(const double& x, const double& y) const;
	double GradientY(const double& x, const double& y) const;

	arma::mat GradientsX(const arma::mat& x, const arma::mat& y) const;
	arma::mat GradientsY(const arma::mat& x, const arma::mat& y) const;

private:

	void MOMSMatrix();
	void IterateParameter(double& z, double& a);
	void ConstructCoefficients(const arma::mat& data);

	double lamda_;
	arma::mat::fixed<4,4> bsv_;
	arma::mat::fixed<4,3> bsg_;

	arma::mat lut_;
};