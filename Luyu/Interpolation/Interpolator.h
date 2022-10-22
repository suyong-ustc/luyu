#pragma once
#include <armadillo>



class Interpolator
{
public:
	Interpolator();
	virtual ~Interpolator() = 0;

	virtual double Value(const double& x, const double& y) const = 0;
	virtual arma::mat Values(const arma::mat& x, const arma::mat& y) const = 0;

	virtual double GradientX(const double& x, const double& y) const = 0;
	virtual double GradientY(const double& x, const double& y) const = 0;

	virtual arma::mat GradientsX(const arma::mat& x, const arma::mat& y) const = 0;
	virtual arma::mat GradientsY(const arma::mat& x, const arma::mat& y) const = 0;

};