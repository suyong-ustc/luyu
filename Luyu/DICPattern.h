#pragma once
#include <armadillo>


class DICPattern
{
public:
	DICPattern(const DICPattern& pattern, const double& u0, const double& v0);
	DICPattern(const double& radius, const arma::vec& x, const arma::vec& y);
	DICPattern(const double& radius, const int& n, const double& a, const double& b);
	DICPattern(const std::string& file_name);
	~DICPattern();

	double Value(const double& x, const double& y) const;
	double GradientX(const double& x, const double& y) const;
	double GradientY(const double& x, const double& y) const;

	arma::mat Values(const arma::mat& x, const arma::mat& y) const;
	arma::mat GradientsX(const arma::mat& x, const arma::mat& y) const;
	arma::mat GradientsY(const arma::mat& x, const arma::mat& y) const;

	arma::mat Sample(const double& x0, const double& y0, const int& half_image_size) const;

	bool Write(const std::string& file_name) const;
	bool Read(const std::string file_name);

private:

	double radius_;
	arma::vec x_;
	arma::vec y_;

};

