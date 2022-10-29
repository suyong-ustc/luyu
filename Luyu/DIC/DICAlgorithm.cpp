#include "DICAlgorithm.h"
#include "..\Interpolation\BicubicMOMSInterpolator.h"
#include "..\Interpolation\BiquinticBSplineInterpolatror.h"
#include "..\Interpolation\BisepticBSplineInterpolator.h"
using namespace arma;




/********************************************************************************************************/
/******************************                    �κ���              **********************************/
/********************************************************************************************************/


bool ShapeFunction(const double& x0, const double& y0, const mat& x, const mat& y, const vec& p, const int& order, mat& deform_x, mat& deform_y)
{
	bool ok = false;

	if (order == 0)
		ok = ZeroOrderShapeFunction(x0, y0, x, y, p, deform_x, deform_y);
	else if (order == 1)
		ok = FirstOrderShapeFunction(x0, y0, x, y, p, deform_x, deform_y);
	else if (order == 2)
		ok = SecondOrderShapeFunction(x0, y0, x, y, p, deform_x, deform_y);
	else if (order == 3)
		ok = ThirdOrderShapeFunction(x0, y0, x, y, p, deform_x, deform_y);
	else if (order == 4)
		ok = FourthOrderShapeFunction(x0, y0, x, y, p, deform_x, deform_y);
	else if (order == 5)
		ok = FifthOrderShapeFunction(x0, y0, x, y, p, deform_x, deform_y);

	return ok;
}



bool ZeroOrderShapeFunction(const double& x0, const double& y0, const mat& x, const mat& y, const vec& p, mat& deform_x, mat& deform_y)
{
	// ����κ�������
	const double u = p(0);
	const double v = p(1);

	// ��������κ������Ʊ��κ�����λ��
	deform_x = x + u + x0;
	deform_y = y + v + y0;

	return true;
}



bool FirstOrderShapeFunction(const double& x0, const double& y0, const mat& x, const mat& y, const vec& p, mat& deform_x, mat& deform_y)
{
	// ����κ�������
	const double u = p(0);
	const double v = p(1);

	// ��������κ������Ʊ��κ�����λ��
	deform_x = x + u + x0;
	deform_y = y + v + y0;

	// һ���κ�������
	const double ux = p(2);
	const double uy = p(3);

	const double vx = p(4);
	const double vy = p(5);

	// ����һ���κ������Ʊ��κ�����λ��
	deform_x = deform_x + ux * x + uy * y;
	deform_y = deform_y + vx * x + vy * y;

	return true;
}



bool SecondOrderShapeFunction(const double& x0, const double& y0, const mat& x, const mat& y, const vec& p, mat& deform_x, mat& deform_y)
{
	// ����κ�������
	const double u = p(0);
	const double v = p(1);

	// ��������κ������Ʊ��κ�����λ��
	deform_x = x + u + x0;
	deform_y = y + v + y0;

	// һ���κ�������
	const double ux = p(2);
	const double uy = p(3);

	const double vx = p(4);
	const double vy = p(5);

	// ����һ���κ������Ʊ��κ�����λ��
	deform_x = deform_x + ux * x + uy * y;
	deform_y = deform_y + vx * x + vy * y;

	// �����κ�������
	const double uxx = p(6);
	const double uxy = p(7);
	const double uyy = p(8);

	const double vxx = p(9);
	const double vxy = p(10);
	const double vyy = p(11);

	// �����κ������Ʊ��κ�����λ��
	const mat xx = x % x;
	const mat xy = x % y;
	const mat yy = y % y;

	deform_x = deform_x + uxx * xx + uxy * xy + uyy * yy;
	deform_y = deform_y + vxx * xx + vxy * xy + vyy * yy;

	return true;
}



bool ThirdOrderShapeFunction(const double& x0, const double& y0, const mat& x, const mat& y, const vec& p, mat& deform_x, mat& deform_y)
{
	// ����κ�������
	const double u = p(0);
	const double v = p(1);

	// ��������κ������Ʊ��κ�����λ��
	deform_x = x + u + x0;
	deform_y = y + v + y0;

	// һ���κ�������
	const double ux = p(2);
	const double uy = p(3);

	const double vx = p(4);
	const double vy = p(5);

	// ����һ���κ������Ʊ��κ�����λ��
	deform_x = deform_x + ux * x + uy * y;
	deform_y = deform_y + vx * x + vy * y;

	// �����κ�������
	const double uxx = p(6);
	const double uxy = p(7);
	const double uyy = p(8);

	const double vxx = p(9);
	const double vxy = p(10);
	const double vyy = p(11);

	// �����κ������Ʊ��κ�����λ��
	const mat xx = x % x;
	const mat xy = x % y;
	const mat yy = y % y;

	deform_x = deform_x + uxx * xx + uxy * xy + uyy * yy;
	deform_y = deform_y + vxx * xx + vxy * xy + vyy * yy;

	// �����κ�������
	const double uxxx = p(12);
	const double uxxy = p(13);
	const double uxyy = p(14);
	const double uyyy = p(15);

	const double vxxx = p(16);
	const double vxxy = p(17);
	const double vxyy = p(18);
	const double vyyy = p(19);

	// �����κ������Ʊ��κ�����λ��
	const mat xxx = xx % x;
	const mat xxy = xx % y;
	const mat xyy = x % yy;
	const mat yyy = y % yy;

	deform_x = deform_x + uxxx * xxx + uxxy * xxy + uxyy * xyy + uyyy * yyy;
	deform_y = deform_y + vxxx * xxx + vxxy * xxy + vxyy * xyy + vyyy * yyy;

	return true;
}



bool FourthOrderShapeFunction(const double& x0, const double& y0, const mat& x, const mat& y, const vec& p, mat& deform_x, mat& deform_y)
{
	// ����κ�������
	const double u = p(0);
	const double v = p(1);

	// ��������κ������Ʊ��κ�����λ��
	deform_x = x + u + x0;
	deform_y = y + v + y0;

	// һ���κ�������
	const double ux = p(2);
	const double uy = p(3);

	const double vx = p(4);
	const double vy = p(5);

	// ����һ���κ������Ʊ��κ�����λ��
	deform_x = deform_x + ux * x + uy * y;
	deform_y = deform_y + vx * x + vy * y;

	// �����κ�������
	const double uxx = p(6);
	const double uxy = p(7);
	const double uyy = p(8);

	const double vxx = p(9);
	const double vxy = p(10);
	const double vyy = p(11);

	// �����κ������Ʊ��κ�����λ��
	const mat xx = x % x;
	const mat xy = x % y;
	const mat yy = y % y;

	deform_x = deform_x + uxx * xx + uxy * xy + uyy * yy;
	deform_y = deform_y + vxx * xx + vxy * xy + vyy * yy;

	// �����κ�������
	const double uxxx = p(12);
	const double uxxy = p(13);
	const double uxyy = p(14);
	const double uyyy = p(15);

	const double vxxx = p(16);
	const double vxxy = p(17);
	const double vxyy = p(18);
	const double vyyy = p(19);

	// �����κ������Ʊ��κ�����λ��
	const mat xxx = xx % x;
	const mat xxy = xx % y;
	const mat xyy = x % yy;
	const mat yyy = y % yy;

	deform_x = deform_x + uxxx * xxx + uxxy * xxy + uxyy * xyy + uyyy * yyy;
	deform_y = deform_y + vxxx * xxx + vxxy * xxy + vxyy * xyy + vyyy * yyy;

	// �Ľ��κ�������
	const double uxxxx = p(20);
	const double uxxxy = p(21);
	const double uxxyy = p(22);
	const double uxyyy = p(23);
	const double uyyyy = p(24);

	const double vxxxx = p(25);
	const double vxxxy = p(26);
	const double vxxyy = p(27);
	const double vxyyy = p(28);
	const double vyyyy = p(29);

	// �����κ������Ʊ��κ�����λ��
	const mat xxxx = xxx % x;
	const mat xxxy = xxx % y;
	const mat xxyy = xx % yy;
	const mat xyyy = xy % yy;
	const mat yyyy = yy % yy;

	deform_x = deform_x + uxxxx * xxxx + uxxxy * xxxy + uxxyy * xxyy + uxyyy * xyyy + uyyyy * yyyy;
	deform_y = deform_y + vxxxx * xxxx + vxxxy * xxxy + vxxyy * xxyy + vxyyy * xyyy + vyyyy * yyyy;

	return true;
}



bool FifthOrderShapeFunction(const double& x0, const double& y0, const mat& x, const mat& y, const vec& p, mat& deform_x, mat& deform_y)
{
	// ����κ�������
	const double u = p(0);
	const double v = p(1);

	// ��������κ������Ʊ��κ�����λ��
	deform_x = x + u + x0;
	deform_y = y + v + y0;

	// һ���κ�������
	const double ux = p(2);
	const double uy = p(3);

	const double vx = p(4);
	const double vy = p(5);

	// ����һ���κ������Ʊ��κ�����λ��
	deform_x = deform_x + ux * x + uy * y;
	deform_y = deform_y + vx * x + vy * y;

	// �����κ�������
	const double uxx = p(6);
	const double uxy = p(7);
	const double uyy = p(8);

	const double vxx = p(9);
	const double vxy = p(10);
	const double vyy = p(11);

	// �����κ������Ʊ��κ�����λ��
	const mat xx = x % x;
	const mat xy = x % y;
	const mat yy = y % y;

	deform_x = deform_x + uxx * xx + uxy * xy + uyy * yy;
	deform_y = deform_y + vxx * xx + vxy * xy + vyy * yy;

	// �����κ�������
	const double uxxx = p(12);
	const double uxxy = p(13);
	const double uxyy = p(14);
	const double uyyy = p(15);

	const double vxxx = p(16);
	const double vxxy = p(17);
	const double vxyy = p(18);
	const double vyyy = p(19);

	// �����κ������Ʊ��κ�����λ��
	const mat xxx = xx % x;
	const mat xxy = xx % y;
	const mat xyy = x % yy;
	const mat yyy = y % yy;

	deform_x = deform_x + uxxx * xxx + uxxy * xxy + uxyy * xyy + uyyy * yyy;
	deform_y = deform_y + vxxx * xxx + vxxy * xxy + vxyy * xyy + vyyy * yyy;

	// �Ľ��κ�������
	const double uxxxx = p(20);
	const double uxxxy = p(21);
	const double uxxyy = p(22);
	const double uxyyy = p(23);
	const double uyyyy = p(24);

	const double vxxxx = p(25);
	const double vxxxy = p(26);
	const double vxxyy = p(27);
	const double vxyyy = p(28);
	const double vyyyy = p(29);

	// �����κ������Ʊ��κ�����λ��
	const mat xxxx = xxx % x;
	const mat xxxy = xxx % y;
	const mat xxyy = xx % yy;
	const mat xyyy = xy % yy;
	const mat yyyy = yy % yy;

	deform_x = deform_x + uxxxx * xxxx + uxxxy * xxxy + uxxyy * xxyy + uxyyy * xyyy + uyyyy * yyyy;
	deform_y = deform_y + vxxxx * xxxx + vxxxy * xxxy + vxxyy * xxyy + vxyyy * xyyy + vyyyy * yyyy;

	// ����κ�������
	const double uxxxxx = p(30);
	const double uxxxxy = p(31);
	const double uxxxyy = p(32);
	const double uxxyyy = p(33);
	const double uxyyyy = p(34);
	const double uyyyyy = p(35);

	const double vxxxxx = p(36);
	const double vxxxxy = p(37);
	const double vxxxyy = p(38);
	const double vxxyyy = p(39);
	const double vxyyyy = p(40);
	const double vyyyyy = p(41);

	// �����κ������Ʊ��κ�����λ��
	const mat xxxxx = xxx % xx;
	const mat xxxxy = xxx % xy;
	const mat xxxyy = xxx % yy;
	const mat xxyyy = xx % yyy;
	const mat xyyyy = xy % yyy;
	const mat yyyyy = yy % yyy;

	deform_x = deform_x + uxxxxx * xxxxx + uxxxxy * xxxxy + uxxxyy * xxxyy + uxxyyy * xxyyy + uxyyyy * xyyyy + uyyyyy * yyyyy;
	deform_y = deform_y + vxxxxx * xxxxx + vxxxxy * xxxxy + vxxxyy * xxxyy + vxxyyy * xxyyy + vxyyyy * xyyyy + vyyyyy * yyyyy;

	return true;
}



/*******************************************************************************************************************/
/**********************************               α�����                ********************************************/
/******************************************************************************************************************/



bool PseudoInverseMatrix(const mat& gx, const mat& gy, const mat& x, const mat& y, const int& order, mat& pseudo)
{
	bool ok = false;

	if (order == 0)
		ok = ZeroOrderPseudoInverseMatrix(gx, gy, x, y, pseudo);
	else if (order == 1)
		ok = FirstOrderPseudoInverseMatrix(gx, gy, x, y, pseudo);
	else if (order == 2)
		ok = SecondOrderPseudoInverseMatrix(gx, gy, x, y, pseudo);
	else if (order == 3)
		ok = ThirdOrderPseudoInverseMatrix(gx, gy, x, y, pseudo);
	else if (order == 4)
		ok = FourthOrderPseudoInverseMatrix(gx, gy, x, y, pseudo);
	else if (order == 5)
		ok = FifthOrderPseudoInverseMatrix(gx, gy, x, y, pseudo);

	return ok;
}



bool ZeroOrderPseudoInverseMatrix(const mat& gx, const mat& gy, const mat& x, const mat& y, mat& pseudo)
{
	// �ſɱȾ���
	mat jacobian(gx.n_elem, 2);

	jacobian.col(0) = vectorise(gx);
	jacobian.col(1) = vectorise(gy);

	// α�����
	pseudo = inv_sympd(jacobian.t() * jacobian) * jacobian.t();

	return true;
}



bool FirstOrderPseudoInverseMatrix(const mat& gx, const mat& gy, const mat& x, const mat& y, mat& pseudo)
{
	// �ſɱȾ���
	mat jacobian(gx.n_elem, 6);

	jacobian.col(0) = vectorise(gx);
	jacobian.col(1) = vectorise(gy);

	jacobian.col(2) = vectorise(gx % x);
	jacobian.col(3) = vectorise(gx % y);
	jacobian.col(4) = vectorise(gy % x);
	jacobian.col(5) = vectorise(gy % y);

	// α�����
	pseudo = inv_sympd(jacobian.t() * jacobian) * jacobian.t();

	return true;
}



bool SecondOrderPseudoInverseMatrix(const mat& gx, const mat& gy, const mat& x, const mat& y, mat& pseudo)
{
	// �ſɱȾ���
	mat jacobian = zeros<mat>(gx.n_elem, 12);

	jacobian.col(0) = vectorise(gx);
	jacobian.col(1) = vectorise(gy);

	jacobian.col(2) = vectorise(gx % x);
	jacobian.col(3) = vectorise(gx % y);
	jacobian.col(4) = vectorise(gy % x);
	jacobian.col(5) = vectorise(gy % y);

	jacobian.col(6) = vectorise(gx % x % x);
	jacobian.col(7) = vectorise(gx % x % y);
	jacobian.col(8) = vectorise(gx % y % y);
	jacobian.col(9) = vectorise(gy % x % x);
	jacobian.col(10) = vectorise(gy % x % y);
	jacobian.col(11) = vectorise(gy % y % y);

	// α�����
	pseudo = inv_sympd(jacobian.t() * jacobian) * jacobian.t();

	return true;
}



bool ThirdOrderPseudoInverseMatrix(const mat& gx, const mat& gy, const mat& x, const mat& y, mat& pseudo)
{
	// �ſɱȾ���
	mat jacobian = zeros<mat>(gx.n_elem, 20);

	jacobian.col(0) = vectorise(gx);
	jacobian.col(1) = vectorise(gy);

	jacobian.col(2) = vectorise(gx % x);
	jacobian.col(3) = vectorise(gx % y);
	jacobian.col(4) = vectorise(gy % x);
	jacobian.col(5) = vectorise(gy % y);

	jacobian.col(6) = vectorise(gx % x % x);
	jacobian.col(7) = vectorise(gx % x % y);
	jacobian.col(8) = vectorise(gx % y % y);
	jacobian.col(9) = vectorise(gy % x % x);
	jacobian.col(10) = vectorise(gy % x % y);
	jacobian.col(11) = vectorise(gy % y % y);

	jacobian.col(12) = vectorise(gx % x % x % x);
	jacobian.col(13) = vectorise(gx % x % x % y);
	jacobian.col(14) = vectorise(gx % x % y % y);
	jacobian.col(15) = vectorise(gx % y % y % y);
	jacobian.col(16) = vectorise(gy % x % x % x);
	jacobian.col(17) = vectorise(gy % x % x % y);
	jacobian.col(18) = vectorise(gy % x % y % y);
	jacobian.col(19) = vectorise(gy % y % y % y);

	// α�����
	pseudo = inv_sympd(jacobian.t() * jacobian) * jacobian.t();

	return true;
}



bool FourthOrderPseudoInverseMatrix(const mat& gx, const mat& gy, const mat& x, const mat& y, mat& pseudo)
{
	// �ſɱȾ���
	mat jacobian = zeros<mat>(gx.n_elem, 30);

	jacobian.col(0) = vectorise(gx);
	jacobian.col(1) = vectorise(gy);

	jacobian.col(2) = vectorise(gx % x);
	jacobian.col(3) = vectorise(gx % y);
	jacobian.col(4) = vectorise(gy % x);
	jacobian.col(5) = vectorise(gy % y);

	jacobian.col(6) = vectorise(gx % x % x);
	jacobian.col(7) = vectorise(gx % x % y);
	jacobian.col(8) = vectorise(gx % y % y);
	jacobian.col(9) = vectorise(gy % x % x);
	jacobian.col(10) = vectorise(gy % x % y);
	jacobian.col(11) = vectorise(gy % y % y);

	jacobian.col(12) = vectorise(gx % x % x % x);
	jacobian.col(13) = vectorise(gx % x % x % y);
	jacobian.col(14) = vectorise(gx % x % y % y);
	jacobian.col(15) = vectorise(gx % y % y % y);
	jacobian.col(16) = vectorise(gy % x % x % x);
	jacobian.col(17) = vectorise(gy % x % x % y);
	jacobian.col(18) = vectorise(gy % x % y % y);
	jacobian.col(19) = vectorise(gy % y % y % y);

	jacobian.col(20) = vectorise(gx % x % x % x % x);
	jacobian.col(21) = vectorise(gx % x % x % x % y);
	jacobian.col(22) = vectorise(gx % x % x % y % y);
	jacobian.col(23) = vectorise(gx % x % y % y % y);
	jacobian.col(24) = vectorise(gx % y % y % y % y);
	jacobian.col(25) = vectorise(gy % x % x % x % x);
	jacobian.col(26) = vectorise(gy % x % x % x % y);
	jacobian.col(27) = vectorise(gy % x % x % y % y);
	jacobian.col(28) = vectorise(gy % x % y % y % y);
	jacobian.col(29) = vectorise(gy % y % y % y % y);

	// α�����
	pseudo = inv_sympd(jacobian.t() * jacobian) * jacobian.t();

	return true;
}



bool FifthOrderPseudoInverseMatrix(const mat& gx, const mat& gy, const mat& x, const mat& y, mat& pseudo)
{
	// �ſɱȾ���
	mat jacobian = zeros<mat>(gx.n_elem, 42);

	jacobian.col(0) = vectorise(gx);
	jacobian.col(1) = vectorise(gy);

	jacobian.col(2) = vectorise(gx % x);
	jacobian.col(3) = vectorise(gx % y);
	jacobian.col(4) = vectorise(gy % x);
	jacobian.col(5) = vectorise(gy % y);

	jacobian.col(6) = vectorise(gx % x % x);
	jacobian.col(7) = vectorise(gx % x % y);
	jacobian.col(8) = vectorise(gx % y % y);
	jacobian.col(9) = vectorise(gy % x % x);
	jacobian.col(10) = vectorise(gy % x % y);
	jacobian.col(11) = vectorise(gy % y % y);

	jacobian.col(12) = vectorise(gx % x % x % x);
	jacobian.col(13) = vectorise(gx % x % x % y);
	jacobian.col(14) = vectorise(gx % x % y % y);
	jacobian.col(15) = vectorise(gx % y % y % y);
	jacobian.col(16) = vectorise(gy % x % x % x);
	jacobian.col(17) = vectorise(gy % x % x % y);
	jacobian.col(18) = vectorise(gy % x % y % y);
	jacobian.col(19) = vectorise(gy % y % y % y);

	jacobian.col(20) = vectorise(gx % x % x % x % x);
	jacobian.col(21) = vectorise(gx % x % x % x % y);
	jacobian.col(22) = vectorise(gx % x % x % y % y);
	jacobian.col(23) = vectorise(gx % x % y % y % y);
	jacobian.col(24) = vectorise(gx % y % y % y % y);
	jacobian.col(25) = vectorise(gy % x % x % x % x);
	jacobian.col(26) = vectorise(gy % x % x % x % y);
	jacobian.col(27) = vectorise(gy % x % x % y % y);
	jacobian.col(28) = vectorise(gy % x % y % y % y);
	jacobian.col(29) = vectorise(gy % y % y % y % y);

	jacobian.col(30) = vectorise(gx % x % x % x % x % x);
	jacobian.col(31) = vectorise(gx % x % x % x % x % y);
	jacobian.col(32) = vectorise(gx % x % x % x % y % y);
	jacobian.col(33) = vectorise(gx % x % x % y % y % y);
	jacobian.col(34) = vectorise(gx % x % y % y % y % y);
	jacobian.col(35) = vectorise(gx % y % y % y % y % y);
	jacobian.col(36) = vectorise(gy % x % x % x % x % x);
	jacobian.col(37) = vectorise(gy % x % x % x % x % y);
	jacobian.col(38) = vectorise(gy % x % x % x % y % y);
	jacobian.col(39) = vectorise(gy % x % x % y % y % y);
	jacobian.col(40) = vectorise(gy % x % y % y % y % y);
	jacobian.col(41) = vectorise(gy % y % y % y % y % y);

	// α�����
	pseudo = inv_sympd(jacobian.t() * jacobian) * jacobian.t();

	return true;
}






/*******************************************************************************************************************/
/**********************************            ����ȫ������             ********************************************/
/******************************************************************************************************************/


DICOutput* RigsterFullFieldDisplacement(const arma::mat& refer_image, const arma::mat& deform_image, const DICParameters& dic_parameters)
{
	// ��������
	mat grid_x;
	mat grid_y;
	dic_parameters.grid(grid_x, grid_y);

	// ��ʼ��������
	DICOutput* dic_output = new DICOutput(grid_x, grid_y, ParameterTotalNumber(dic_parameters.shape_function_order()));

	// ���Ƴ�ֵ
	mat u;
	mat v;
	EstimateInitialDisplacement(refer_image, deform_image, grid_x, grid_y, u, v);

	// ��ؼ���
	RegisterSubpixelDisplacement(refer_image, deform_image, grid_x, grid_y, dic_parameters, u, v, dic_output);


	return dic_output;

}



bool EstimateInitialDisplacement(const mat& refer_image, const mat& deform_image, const mat& x, const mat& y, mat& u, mat& v)
{
	// ��ʼ��
	u.zeros(x.n_rows, x.n_cols);
	v.zeros(x.n_rows, x.n_cols);

	// ����ֵ
	bool is_harmonic = true;

	if (is_harmonic)
	{
		double a = 5;
		double T = 20;
		double b = 0;

		u = a * sin(2 * datum::pi * x / T + b);
	}
	else
	{
		double a = 5;
		double x0 = 150;
		double c = 20;

		mat t = (x - x0) / c;
		u = a * exp(-t % t);
	}

	return true;
}



bool RegisterSubpixelDisplacement(const mat& refer_image, const mat& deform_image, const mat& x, const mat& y, const DICParameters& dic_parameters, const mat& u, const mat& v, DICOutput* dic_output)
{
	/**************************************************************************************************************/
	/*                                             Ԥ����                                                         */
	/**************************************************************************************************************/

	// �����ߴ�
	const int subset_size = dic_parameters.subset_size();
	const int half_subset_size = (subset_size - 1) / 2;

	// "������"�����"��������"������
	mat delta_x = zeros(subset_size, subset_size);
	mat delta_y = zeros(subset_size, subset_size);
	for (int dx = -half_subset_size; dx <= half_subset_size; ++dx)
	{
		const int c = dx + half_subset_size;
		for (int dy = -half_subset_size; dy <= half_subset_size; ++dy)
		{
			const int r = dy + half_subset_size;
			delta_x(r, c) = dx;
			delta_y(r, c) = dy;
		}

	}


	// ����ͼ��ֵϵ��
	Interpolator* deform_image_interpolator;

	if (dic_parameters.bspline_interpolation_order() == 3)
		deform_image_interpolator = new BicubicMOMSInterpolator(deform_image);
	else if (dic_parameters.bspline_interpolation_order() == 7)
		deform_image_interpolator = new BisepticBSplineInterpolator(deform_image); 
	else
		deform_image_interpolator = new BiquinticBSplineInterpolatror(deform_image);


	/**************************************************************************************************************/
	/*                                      ��˹ţ�ٷ������ص���                                                    */
	/**************************************************************************************************************/

	
	for (int r = 0; r < x.n_rows; ++r)
	{
		for (int c = 0; c < x.n_cols; ++c)
		{
			// ȷ������
			const int subset_center_x = x(r, c);
			const int subset_center_y = y(r, c);

			const int subset_x_min = subset_center_x - half_subset_size;
			const int subset_x_max = subset_center_x + half_subset_size;
			const int subset_y_min = subset_center_y - half_subset_size;
			const int subset_y_max = subset_center_y + half_subset_size;

			const mat refer_subset_intensities = refer_image.submat(subset_y_min, subset_y_max, subset_x_min, subset_x_max);

			// ������ֵ
			const int parameter_total_number = ParameterTotalNumber(dic_parameters.shape_function_order());
			vec warp_function = zeros<vec>(parameter_total_number);

			warp_function(0) = u(r, c);
			warp_function(1) = v(r, c);

			// ���ο������Ҷȹ�һ��
			vec normalized_refer_subset_intensities;
			NormalizeVectorize(refer_subset_intensities, normalized_refer_subset_intensities);

			// ��������
			uword iter = 0;

			mat deform_x, deform_y;			// �����������ͼ������
			mat deform_subset_intensities;	// ����������ĻҶ�
			mat deform_subset_gradient_x;	// ����������ĻҶ��ݶ�
			mat deform_subset_gradient_y;
			mat pseudo;						// α�����
			vec normalized_deform_subset_intensities;	// ��һ���ı��������Ҷ�

			// ����
			vec dp = ones<vec>(parameter_total_number);	// ��������
			while (!isConverge(dp) && iter <= dic_parameters.max_iteration_times())
			{
				// ȷ�������������"ͼ������"
				ShapeFunction(subset_center_x, subset_center_y, delta_x, delta_y, warp_function, dic_parameters.shape_function_order(), deform_x, deform_y);

				// ����������ĻҶ�
				deform_subset_intensities = deform_image_interpolator->Values(deform_x, deform_y);
				const double gg = NormalizeVectorize(deform_subset_intensities, normalized_deform_subset_intensities);

				// ����������ĻҶ��ݶ�
				deform_subset_gradient_x = deform_image_interpolator->GradientsX(deform_x, deform_y);
				deform_subset_gradient_y = deform_image_interpolator->GradientsY(deform_x, deform_y);

				// α�����
				PseudoInverseMatrix(deform_subset_gradient_x, deform_subset_gradient_y, delta_x, delta_x, dic_parameters.shape_function_order(), pseudo);

				// ��һ����"�ο�����"��"��������"�ĻҶȲ�
				const vec diff = normalized_refer_subset_intensities - normalized_deform_subset_intensities;

				// ��������
				dp = gg * pseudo * diff;
				warp_function = warp_function + dp;

				// ����������1
				++iter;
			}

			// ���յ����ϵ��
			ShapeFunction(subset_center_x, subset_center_y, delta_x, delta_y, warp_function, dic_parameters.shape_function_order(), deform_x, deform_y);
			deform_subset_intensities = deform_image_interpolator->Values(deform_x, deform_y);

			const double zncc = ZNCC(refer_subset_intensities, deform_subset_intensities);

			// ��¼������
			dic_output->set_warp_function(r, c, warp_function);
			dic_output->set_iteration_times(r, c, iter);
			dic_output->set_zncc(r, c, zncc);

		}
	}


	return true;

}



bool isConverge(const vec& dp)
{
	// ��������
	const double dx = dp(0);
	const double dy = dp(1);
	const double dr2 = dx * dx + dy * dy;

	// �����о�
	bool is_converge = false;
	if (dr2 < 1e-6)
		is_converge = true;

	return is_converge;
}



/*******************************************************************************************************************/
/**********************************            �������㺯��              ********************************************/
/******************************************************************************************************************/



int ParameterTotalNumber(const int& order)
{
	int n = 0;

	if (order == 0)
		n = 2;
	else if (order == 1)
		n = 6;
	else if (order == 2)
		n = 12;
	else if (order == 3)
		n = 20;
	else if (order == 4)
		n = 30;
	else if (order == 5)
		n = 42;

	return n;
}



double NormalizeVectorize(const mat& m, vec& v)
{
	// ������
	v = vectorise(m);

	// ���ֵ
	v = v - mean(v);

	// ��һ��
	const double vv = norm(v, 2);
	v = v / vv;

	return vv;
}



double ZNCC(const mat& a, const mat& b)
{
	// ������
	vec aa = vectorise(a);
	vec bb = vectorise(b);

	// ���ֵ
	aa = aa - mean(aa);
	bb = bb - mean(bb);

	// ��һ��
	aa = normalise(aa);
	bb = normalise(bb);

	// ���ϵ��
	const double zncc = dot(aa, bb);

	return zncc;
}