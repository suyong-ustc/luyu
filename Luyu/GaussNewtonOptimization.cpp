#include "GaussNewtonOptimization.h"
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

//#include "FAGN1.h"
//#include "DICConvergence.h"
//#include "DICShapeFunction.h"
//#include "DICTypes.h"
//#include "Interpolator.h"
//#include "Util.h"
//using namespace arma;
//
//
//
//bool ForwardAdditiveGaussNewtonMethod1(
//	const mat& refer_intensities,
//	const double refer_x0,
//	const double refer_y0,
//	Interpolator* deform_image_interp,
//	const int half_subset_size,
//	const int max_iter_time,
//	const double zncc_threshold,
//	vec& deform_warp,
//	arma::uword& iter,
//	double& zncc)
//{
//	// "������"�����"��������"������
//	const int subset_size = 2 * half_subset_size + 1;
//
//	mat delta_x = zeros(subset_size, subset_size);
//	mat delta_y = zeros(subset_size, subset_size);
//	for (int dx = -half_subset_size; dx <= half_subset_size; ++dx)
//	{
//		const int c = dx + half_subset_size;
//		for (int dy = -half_subset_size; dy <= half_subset_size; ++dy)
//		{
//			const int r = dy + half_subset_size;
//			delta_x(r, c) = dx;
//			delta_y(r, c) = dy;
//		}
//
//	}
//
//	// ���ο������Ҷȹ�һ��
//	vec normalized_refer_intensities;
//	NormalizeVectoerize(refer_intensities, normalized_refer_intensities);
//
//	// ��������
//	iter = 0;
//
//	mat deform_subset_x, deform_subset_y;		// �����������ͼ������
//	mat deform_intensities;						// ����������ĻҶ�
//	mat deform_gradient_x, deform_gradient_y;	// ����������ĻҶ��ݶ�
//	mat pseudo;									// α�����
//	vec normalized_deform_intensities;			// ��һ���ı��������Ҷ�
//
//	// ����
//	vec dp = ones<vec>(FIRST_ORDER_SHAPE_FUNC_PARA_COUNT);	// ��������
//	while (!isConvergent1(dp, half_subset_size, half_subset_size) && iter <= max_iter_time)
//	{
//		// ȷ�������������"ͼ������"
//		ShapeFunction1(refer_x0, refer_y0, delta_x, delta_y, deform_warp, deform_subset_x, deform_subset_y);
//
//		// ����������ĻҶ�
//		deform_intensities = deform_image_interp->Values(deform_subset_x, deform_subset_y);
//		const double gg = NormalizeVectoerize(deform_intensities, normalized_deform_intensities);
//
//		// ����α�����
//		deform_gradient_x = deform_image_interp->GradientsX(deform_subset_x, deform_subset_y);
//		deform_gradient_y = deform_image_interp->GradientsY(deform_subset_x, deform_subset_y);
//		ForwardPseudoMatrix1(deform_gradient_x, deform_gradient_y, delta_x, delta_y, pseudo);
//
//		// ��һ����"�ο�����"��"��������"�ĻҶȲ�
//		const vec diff = normalized_refer_intensities - normalized_deform_intensities;
//
//		// ��������
//		dp = gg * pseudo * diff;
//		deform_warp = deform_warp + dp;
//
//		//double u = dp(0);
//		//double ux = dp(1);
//		//double uy = dp(2);
//		//double v = dp(3);
//		//double vx = dp(4);
//		//double vy = dp(5);
//
//		// ����������1
//		++iter;
//	}
//
//	// �жϵ����Ƿ�����
//	bool is_convergent = true;
//
//	if (iter >= max_iter_time)
//	{
//		std::cerr << "The iteration time " << iter << " reach the maximum!" << std::endl;
//		is_convergent = false;
//	}
//
//	// ���յ����ϵ��
//	ShapeFunction1(refer_x0, refer_y0, delta_x, delta_y, deform_warp, deform_subset_x, deform_subset_y);
//	deform_intensities = deform_image_interp->Values(deform_subset_x, deform_subset_y);
//
//	zncc = ZNCC(refer_intensities, deform_intensities);
//
//	// �ж����ϵ���Ƿ�������ֵ
//	if (zncc < zncc_threshold)
//	{
//		std::cerr << "The ZNCC " << zncc << "is less than the threshold!" << std::endl;
//		is_convergent = false;
//	}
//
//
//	 return is_convergent;
//
//}
//
//
//
//
//
//bool ForwardPseudoMatrix1(const mat& gx, const mat& gy, const mat& dx, const mat& dy, mat& pseudo)
//{
//	// Jacobian matrix
//	mat jacobian(gx.n_elem, FIRST_ORDER_SHAPE_FUNC_PARA_COUNT);
//
//	jacobian.col(0) = vectorise(gx);
//	jacobian.col(1) = vectorise(gx % dx);
//	jacobian.col(2) = vectorise(gx % dy);
//	jacobian.col(3) = vectorise(gy);
//	jacobian.col(4) = vectorise(gy % dx);
//	jacobian.col(5) = vectorise(gy % dy);
//
//	// pseudo inverse matrix
//	pseudo = inv_sympd(jacobian.t() * jacobian) * jacobian.t();
//
//	return true;
//}



bool ForwardAdditiveGaussNewtonMethod2(
	const mat& refer_intensities,
	const double refer_x0,
	const double refer_y0,
	Interpolator* deform_image_interp,
	const int half_subset_size,
	const int max_iter_time,
	const double zncc_threshold,
	vec& deform_warp,
	uword& iter,
	double& zncc)
{
	// "������"�����"��������"������
	const int subset_size = 2 * half_subset_size + 1;

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

	// ���ο������Ҷȹ�һ��
	vec normalized_refer_intensities;
	NormalizeVectoerize(refer_intensities, normalized_refer_intensities);

	// ��������
	iter = 0;

	mat deform_subset_x, deform_subset_y;		// �����������ͼ������
	mat deform_intensities;						// ����������ĻҶ�
	mat deform_gradient_x, deform_gradient_y;	// ����������ĻҶ��ݶ�
	mat pseudo;									// α�����
	vec normalized_deform_intensities;			// ��һ���ı��������Ҷ�

	// ����
	vec dp = ones<vec>(SECOND_ORDER_SHAPE_FUNC_PARA_COUNT);	// ��������
	while (!isConvergent2(dp, half_subset_size, half_subset_size) && iter <= max_iter_time)
	{
		// ȷ�������������"ͼ������"
		ShapeFunction2(refer_x0, refer_y0, delta_x, delta_y, deform_warp, deform_subset_x, deform_subset_y);

		// ����������ĻҶ�
		deform_intensities = deform_image_interp->Values(deform_subset_x, deform_subset_y);
		const double gg = NormalizeVectoerize(deform_intensities, normalized_deform_intensities);

		// ����α�����
		deform_gradient_x = deform_image_interp->GradientsX(deform_subset_x, deform_subset_y);
		deform_gradient_y = deform_image_interp->GradientsY(deform_subset_x, deform_subset_y);
		ForwardPseudoMatrix2(deform_gradient_x, deform_gradient_y, delta_x, delta_y, pseudo);

		// ��һ����"�ο�����"��"��������"�ĻҶȲ�
		const vec diff = normalized_refer_intensities - normalized_deform_intensities;

		// ��������
		dp = gg * pseudo * diff;
		deform_warp = deform_warp + dp;

		// ����������1
		++iter;
	}


	// �жϵ����Ƿ�����
	bool is_convergent = true;

	if (iter >= max_iter_time)
	{
		std::cerr << "The iteration time " << iter << " reach the maximum!" << std::endl;
		is_convergent = false;
	}

	// ���յ����ϵ��
	ShapeFunction2(refer_x0, refer_y0, delta_x, delta_y, deform_warp, deform_subset_x, deform_subset_y);
	deform_intensities = deform_image_interp->Values(deform_subset_x, deform_subset_y);

	zncc = ZNCC(refer_intensities, deform_intensities);

	// �ж����ϵ���Ƿ�������ֵ
	if (zncc < zncc_threshold)
	{
		std::cerr << "The ZNCC " << zncc << "is less than the threshold!" << std::endl;
		is_convergent = false;
	}


	return is_convergent;

}




bool ForwardPseudoMatrix2(
	const mat& gx,
	const mat& gy,
	const mat& delta_x,
	const mat& delta_y,
	mat& pseudo)
{
	// Jacobi matrix
	mat jacobian = zeros<mat>(gx.n_elem, SECOND_ORDER_SHAPE_FUNC_PARA_COUNT);

	jacobian.col(0) = vectorise(gx);
	jacobian.col(1) = vectorise(gx % delta_x);
	jacobian.col(2) = vectorise(gx % delta_y);
	jacobian.col(3) = vectorise(0.5 * gx % delta_x % delta_x);
	jacobian.col(4) = vectorise(gx % delta_x % delta_y);
	jacobian.col(5) = vectorise(0.5 * gx % delta_y % delta_y);

	jacobian.col(6) = vectorise(gy);
	jacobian.col(7) = vectorise(gy % delta_x);
	jacobian.col(8) = vectorise(gy % delta_y);
	jacobian.col(9) = vectorise(0.5 * gy % delta_x % delta_x);
	jacobian.col(10) = vectorise(gy % delta_x % delta_y);
	jacobian.col(11) = vectorise(0.5 * gy % delta_y % delta_y);

	// pseudo inverse matrix
	pseudo = inv_sympd(jacobian.t() * jacobian) * jacobian.t();

	return true;
}