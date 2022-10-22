#pragma once
#include <armadillo>
#include "DICTypes.h"


// parameter count of different shape functions
int ShapeFunctionParameterCount(const ShapeFuncType& shape_function);


// shape function mapping
bool ShapeFunction(
	const double& x0, 
	const double& y0, 
	const arma::mat& delta_x, 
	const arma::mat& delta_y, 
	const arma::vec& p, 
	const ShapeFuncType& shape_function,
	arma::mat& deform_x, 
	arma::mat& deform_y);


bool ShapeFunction1(
	const double& x0, 
	const double& y0, 
	const arma::mat& delta_x, 
	const arma::mat& delta_y, 
	const arma::vec& p,
	arma::mat& deform_x, 
	arma::mat& deform_y);


bool ShapeFunction2(
	const double& x0, 
	const double& y0, 
	const arma::mat& delta_x, 
	const arma::mat& delta_y, 
	const arma::vec& p,
	arma::mat& deform_x, 
	arma::mat& deform_y);


bool ShapeFunction(
	const double& x0, 
	const double& y0, 
	const double& delta_x, 
	const double& delta_y, 
	const arma::vec& p, 
	const ShapeFuncType& shape_function,
	double& deform_x, 
	double& deform_y);


bool ShapeFunction1(
	const double& x0, 
	const double& y0, 
	const double& delta_x, 
	const double& delta_y, 
	const arma::vec& p,
	double& deform_x, 
	double& deform_y);

bool ShapeFunction2(
	const double& x0, 
	const double& y0, 
	const double& delta_x, 
	const double& delta_y, 
	const arma::vec& p,
	double& deform_x, 
	double& deform_y);