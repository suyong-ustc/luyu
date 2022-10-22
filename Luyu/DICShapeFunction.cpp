#include "DICShapeFunction.h"
using namespace arma;



int ShapeFunctionParameterCount(const ShapeFuncType& shape_function)
{
	int shape_func_para_count = 0;

	if (shape_function == SHAPEFUNC_FIRST)
		shape_func_para_count = FIRST_ORDER_SHAPE_FUNC_PARA_COUNT;

	if (shape_function == SHAPEFUNC_SECOND)
		shape_func_para_count = SECOND_ORDER_SHAPE_FUNC_PARA_COUNT;

	return shape_func_para_count;
}



/********************************************************************************************************/
/******************************               Shape Function           **********************************/
/********************************************************************************************************/



bool ShapeFunction(const double& x0, const double& y0, const mat& delta_x, const mat& delta_y, const vec& p, const ShapeFuncType& shape_function,
	mat& deform_x, mat& deform_y)
{
	if (shape_function == SHAPEFUNC_FIRST)
		return ShapeFunction1(x0, y0, delta_x, delta_y, p, deform_x, deform_y);

	if (shape_function == SHAPEFUNC_SECOND)
		return ShapeFunction2(x0, y0, delta_x, delta_y, p, deform_x, deform_y);

	return false;
}



bool ShapeFunction1(const double& x0, const double& y0, const mat& delta_x, const mat& delta_y, const vec&p,
	mat& deform_x, mat& deform_y)
{
	// first-order shape function parameters
	const double u = p(0);
	const double ux = p(1);
	const double uy = p(2);

	const double v = p(3);
	const double vx = p(4);
	const double vy = p(5);

	// the position estimated by the shape function
	deform_x = (1 + ux) * delta_x + uy * delta_y + u + x0;
	deform_y = vx * delta_x + (1 + vy) * delta_y + v + y0;

	return true;
}



bool ShapeFunction2(const double& x0, const double& y0, const mat& delta_x, const mat& delta_y, const vec&p,
	mat& deform_x, mat& deform_y)
{
	// second shape function parameters
	const double u = p(0);
	const double ux = p(1);
	const double uy = p(2);
	const double uxx = p(3);
	const double uxy = p(4);
	const double uyy = p(5);

	const double v = p(6);
	const double vx = p(7);
	const double vy = p(8);
	const double vxx = p(9);
	const double vxy = p(10);
	const double vyy = p(11);

	// the position estimated by the shape function
	deform_x = 0.5 * uxx * delta_x % delta_x + uxy * delta_x % delta_y + 0.5 * uyy * delta_y % delta_y +
		(1 + ux) * delta_x + uy * delta_y + u + x0;

	deform_y = 0.5 * vxx * delta_x % delta_x + vxy * delta_x % delta_y + 0.5 * vyy * delta_y % delta_y +
		vx * delta_x + (1 + vy) * delta_y + v + y0;

	return true;
}



bool ShapeFunction(const double& x0, const double& y0, const double& delta_x, const double& delta_y, const vec& p, const ShapeFuncType& shape_function,
	double& deform_x, double& deform_y)
{
	if (shape_function == SHAPEFUNC_FIRST)
		return ShapeFunction1(x0, y0, delta_x, delta_y, p, deform_x, deform_y);

	if (shape_function == SHAPEFUNC_SECOND)
		return ShapeFunction2(x0, y0, delta_x, delta_y, p, deform_x, deform_y);

	return false;

}



bool ShapeFunction1(const double& x0, const double& y0, const double& delta_x, const double& delta_y, const vec& p,
	double& deform_x, double& deform_y)
{
	// first-order shape function parameters
	const double u = p(0);
	const double ux = p(1);
	const double uy = p(2);

	const double v = p(3);
	const double vx = p(4);
	const double vy = p(5);

	// the position estimated by the shape function
	deform_x = (1 + ux) * delta_x + uy * delta_y + u + x0;
	deform_y = vx * delta_x + (1 + vy) * delta_y + v + y0;

	return true;
}



bool ShapeFunction2(const double& x0, const double& y0, const double& delta_x, const double& delta_y, const vec& p,
	double& deform_x, double& deform_y)
{
	// second shape function parameters
	const double u = p(0);
	const double ux = p(1);
	const double uy = p(2);
	const double uxx = p(3);
	const double uxy = p(4);
	const double uyy = p(5);

	const double v = p(6);
	const double vx = p(7);
	const double vy = p(8);
	const double vxx = p(9);
	const double vxy = p(10);
	const double vyy = p(11);

	// the position estimated by the shape function
	deform_x = 0.5 * uxx * delta_x * delta_x + uxy * delta_x * delta_y + 0.5 * uyy * delta_y * delta_y +
		(1 + ux) * delta_x + uy * delta_y + u + x0;

	deform_y = 0.5 * vxx * delta_x * delta_x + vxy * delta_x * delta_y + 0.5 * vyy * delta_y * delta_y +
		vx * delta_x + (1 + vy) * delta_y + v + y0;

	return true;

}