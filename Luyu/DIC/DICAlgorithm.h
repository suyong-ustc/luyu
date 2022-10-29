#pragma once
#include <armadillo>
#include "DICParameters.h"
#include "DICOutput.h"
#include "../Interpolation/Interpolator.h"


/*******************************************************************************************************************/
/**********************************               形函数                ********************************************/
/******************************************************************************************************************/

bool ShapeFunction(const double& x0, const double& y0, const arma::mat& x, const arma::mat& y, const arma::vec& p, const int& order, arma::mat& deform_x, arma::mat& deform_y);

bool ZeroOrderShapeFunction(const double& x0, const double& y0, const arma::mat& x, const arma::mat& y, const arma::vec& p, arma::mat& deform_x, arma::mat& deform_y);

bool FirstOrderShapeFunction(const double& x0, const double& y0, const arma::mat& x, const arma::mat& y, const arma::vec& p, arma::mat& deform_x, arma::mat& deform_y);

bool SecondOrderShapeFunction(const double& x0, const double& y0, const arma::mat& x, const arma::mat& y, const arma::vec& p, arma::mat& deform_x, arma::mat& deform_y);

bool ThirdOrderShapeFunction(const double& x0, const double& y0, const arma::mat& x, const arma::mat& y, const arma::vec& p, arma::mat& deform_x, arma::mat& deform_y);

bool FourthOrderShapeFunction(const double& x0, const double& y0, const arma::mat& x, const arma::mat& y, const arma::vec& p, arma::mat& deform_x, arma::mat& deform_y);

bool FifthOrderShapeFunction(const double& x0, const double& y0, const arma::mat& x, const arma::mat& y, const arma::vec& p, arma::mat& deform_x, arma::mat& deform_y);



/*******************************************************************************************************************/
/**********************************               伪逆矩阵                ********************************************/
/******************************************************************************************************************/

bool PseudoInverseMatrix(const arma::mat& gx, const arma::mat& gy, const arma::mat& x, const arma::mat& y, const int& order, arma::mat& pseudo);

bool ZeroOrderPseudoInverseMatrix(const arma::mat& gx, const arma::mat& gy, const arma::mat& x, const arma::mat& y, arma::mat& pseudo);

bool FirstOrderPseudoInverseMatrix(const arma::mat& gx, const arma::mat& gy, const arma::mat& x, const arma::mat& y, arma::mat& pseudo);

bool SecondOrderPseudoInverseMatrix(const arma::mat& gx, const arma::mat& gy, const arma::mat& x, const arma::mat& y, arma::mat& pseudo);

bool ThirdOrderPseudoInverseMatrix(const arma::mat& gx, const arma::mat& gy, const arma::mat& x, const arma::mat& y, arma::mat& pseudo);

bool FourthOrderPseudoInverseMatrix(const arma::mat& gx, const arma::mat& gy, const arma::mat& x, const arma::mat& y, arma::mat& pseudo);

bool FifthOrderPseudoInverseMatrix(const arma::mat& gx, const arma::mat& gy, const arma::mat& x, const arma::mat& y, arma::mat& pseudo);



bool isConverge(const arma::vec& dp);


/*******************************************************************************************************************/
/**********************************            计算全场变形             ********************************************/
/******************************************************************************************************************/


DICOutput* RigsterFullFieldDisplacement(const arma::mat& refer_image, const arma::mat& deform_image, const DICParameters& dic_parameters);

bool EstimateInitialDisplacement(const arma::mat& refer_image, const arma::mat& deform_image, const arma::mat& x, const arma::mat& y, arma::mat& u, arma::mat& v);

bool RegisterSubpixelDisplacement(const arma::mat& refer_image, const arma::mat& deform_image, const arma::mat& x, const arma::mat& y, const DICParameters& dic_parameters, const arma::mat& u, const arma::mat& v, DICOutput* dic_output);


/*******************************************************************************************************************/
/**********************************            辅助计算函数              ********************************************/
/******************************************************************************************************************/

int ParameterTotalNumber(const int& order);

double NormalizeVectorize(const arma::mat& m, arma::vec& v);

double ZNCC(const arma::mat& a, const arma::mat& b);