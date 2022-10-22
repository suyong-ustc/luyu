#include "GaussNewtonOptimization.h"
using namespace arma;




/********************************************************************************************************/
/******************************                    形函数              **********************************/
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
	// 零阶形函数参数
	const double u = p(0);
	const double v = p(1);

	// 根据零阶形函数估计变形后子区位置
	deform_x = x + u + x0;
	deform_y = y + v + y0;

	return true;
}



bool FirstOrderShapeFunction(const double& x0, const double& y0, const mat& x, const mat& y, const vec& p, mat& deform_x, mat& deform_y)
{
	// 零阶形函数参数
	const double u = p(0);
	const double v = p(1);

	// 根据零阶形函数估计变形后子区位置
	deform_x = x + u + x0;
	deform_y = y + v + y0;

	// 一阶形函数参数
	const double ux = p(2);
	const double uy = p(3);

	const double vx = p(4);
	const double vy = p(5);

	// 根据一阶形函数估计变形后子区位置
	deform_x = deform_x + ux * x + uy * y;
	deform_y = deform_y + vx * x + vy * y;

	return true;
}



bool SecondOrderShapeFunction(const double& x0, const double& y0, const mat& x, const mat& y, const vec& p, mat& deform_x, mat& deform_y)
{
	// 零阶形函数参数
	const double u = p(0);
	const double v = p(1);

	// 根据零阶形函数估计变形后子区位置
	deform_x = x + u + x0;
	deform_y = y + v + y0;

	// 一阶形函数参数
	const double ux = p(2);
	const double uy = p(3);

	const double vx = p(4);
	const double vy = p(5);

	// 根据一阶形函数估计变形后子区位置
	deform_x = deform_x + ux * x + uy * y;
	deform_y = deform_y + vx * x + vy * y;

	// 二阶形函数参数
	const double uxx = p(6);
	const double uxy = p(7);
	const double uyy = p(8);

	const double vxx = p(9);
	const double vxy = p(10);
	const double vyy = p(11);

	// 根据形函数估计变形后子区位置
	const mat xx = x % x;
	const mat xy = x % y;
	const mat yy = y % y;

	deform_x = deform_x + uxx * xx + uxy * xy + uyy * yy;
	deform_y = deform_y + vxx * xx + vxy * xy + vyy * yy;

	return true;
}



bool ThirdOrderShapeFunction(const double& x0, const double& y0, const mat& x, const mat& y, const vec& p, mat& deform_x, mat& deform_y)
{
	// 零阶形函数参数
	const double u = p(0);
	const double v = p(1);

	// 根据零阶形函数估计变形后子区位置
	deform_x = x + u + x0;
	deform_y = y + v + y0;

	// 一阶形函数参数
	const double ux = p(2);
	const double uy = p(3);

	const double vx = p(4);
	const double vy = p(5);

	// 根据一阶形函数估计变形后子区位置
	deform_x = deform_x + ux * x + uy * y;
	deform_y = deform_y + vx * x + vy * y;

	// 二阶形函数参数
	const double uxx = p(6);
	const double uxy = p(7);
	const double uyy = p(8);

	const double vxx = p(9);
	const double vxy = p(10);
	const double vyy = p(11);

	// 根据形函数估计变形后子区位置
	const mat xx = x % x;
	const mat xy = x % y;
	const mat yy = y % y;

	deform_x = deform_x + uxx * xx + uxy * xy + uyy * yy;
	deform_y = deform_y + vxx * xx + vxy * xy + vyy * yy;

	// 三阶形函数参数
	const double uxxx = p(12);
	const double uxxy = p(13);
	const double uxyy = p(14);
	const double uyyy = p(15);

	const double vxxx = p(16);
	const double vxxy = p(17);
	const double vxyy = p(18);
	const double vyyy = p(19);

	// 根据形函数估计变形后子区位置
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
	// 零阶形函数参数
	const double u = p(0);
	const double v = p(1);

	// 根据零阶形函数估计变形后子区位置
	deform_x = x + u + x0;
	deform_y = y + v + y0;

	// 一阶形函数参数
	const double ux = p(2);
	const double uy = p(3);

	const double vx = p(4);
	const double vy = p(5);

	// 根据一阶形函数估计变形后子区位置
	deform_x = deform_x + ux * x + uy * y;
	deform_y = deform_y + vx * x + vy * y;

	// 二阶形函数参数
	const double uxx = p(6);
	const double uxy = p(7);
	const double uyy = p(8);

	const double vxx = p(9);
	const double vxy = p(10);
	const double vyy = p(11);

	// 根据形函数估计变形后子区位置
	const mat xx = x % x;
	const mat xy = x % y;
	const mat yy = y % y;

	deform_x = deform_x + uxx * xx + uxy * xy + uyy * yy;
	deform_y = deform_y + vxx * xx + vxy * xy + vyy * yy;

	// 三阶形函数参数
	const double uxxx = p(12);
	const double uxxy = p(13);
	const double uxyy = p(14);
	const double uyyy = p(15);

	const double vxxx = p(16);
	const double vxxy = p(17);
	const double vxyy = p(18);
	const double vyyy = p(19);

	// 根据形函数估计变形后子区位置
	const mat xxx = xx % x;
	const mat xxy = xx % y;
	const mat xyy = x % yy;
	const mat yyy = y % yy;

	deform_x = deform_x + uxxx * xxx + uxxy * xxy + uxyy * xyy + uyyy * yyy;
	deform_y = deform_y + vxxx * xxx + vxxy * xxy + vxyy * xyy + vyyy * yyy;

	// 四阶形函数参数
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

	// 根据形函数估计变形后子区位置
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
	// 零阶形函数参数
	const double u = p(0);
	const double v = p(1);

	// 根据零阶形函数估计变形后子区位置
	deform_x = x + u + x0;
	deform_y = y + v + y0;

	// 一阶形函数参数
	const double ux = p(2);
	const double uy = p(3);

	const double vx = p(4);
	const double vy = p(5);

	// 根据一阶形函数估计变形后子区位置
	deform_x = deform_x + ux * x + uy * y;
	deform_y = deform_y + vx * x + vy * y;

	// 二阶形函数参数
	const double uxx = p(6);
	const double uxy = p(7);
	const double uyy = p(8);

	const double vxx = p(9);
	const double vxy = p(10);
	const double vyy = p(11);

	// 根据形函数估计变形后子区位置
	const mat xx = x % x;
	const mat xy = x % y;
	const mat yy = y % y;

	deform_x = deform_x + uxx * xx + uxy * xy + uyy * yy;
	deform_y = deform_y + vxx * xx + vxy * xy + vyy * yy;

	// 三阶形函数参数
	const double uxxx = p(12);
	const double uxxy = p(13);
	const double uxyy = p(14);
	const double uyyy = p(15);

	const double vxxx = p(16);
	const double vxxy = p(17);
	const double vxyy = p(18);
	const double vyyy = p(19);

	// 根据形函数估计变形后子区位置
	const mat xxx = xx % x;
	const mat xxy = xx % y;
	const mat xyy = x % yy;
	const mat yyy = y % yy;

	deform_x = deform_x + uxxx * xxx + uxxy * xxy + uxyy * xyy + uyyy * yyy;
	deform_y = deform_y + vxxx * xxx + vxxy * xxy + vxyy * xyy + vyyy * yyy;

	// 四阶形函数参数
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

	// 根据形函数估计变形后子区位置
	const mat xxxx = xxx % x;
	const mat xxxy = xxx % y;
	const mat xxyy = xx % yy;
	const mat xyyy = xy % yy;
	const mat yyyy = yy % yy;

	deform_x = deform_x + uxxxx * xxxx + uxxxy * xxxy + uxxyy * xxyy + uxyyy * xyyy + uyyyy * yyyy;
	deform_y = deform_y + vxxxx * xxxx + vxxxy * xxxy + vxxyy * xxyy + vxyyy * xyyy + vyyyy * yyyy;

	// 五阶形函数参数
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

	// 根据形函数估计变形后子区位置
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
//	// "子区点"相对于"子区中心"的坐标
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
//	// 将参考子区灰度归一化
//	vec normalized_refer_intensities;
//	NormalizeVectoerize(refer_intensities, normalized_refer_intensities);
//
//	// 迭代次数
//	iter = 0;
//
//	mat deform_subset_x, deform_subset_y;		// 变形子区点的图像坐标
//	mat deform_intensities;						// 变形子区点的灰度
//	mat deform_gradient_x, deform_gradient_y;	// 变形子区点的灰度梯度
//	mat pseudo;									// 伪逆矩阵
//	vec normalized_deform_intensities;			// 归一化的变形子区灰度
//
//	// 迭代
//	vec dp = ones<vec>(FIRST_ORDER_SHAPE_FUNC_PARA_COUNT);	// 迭代增量
//	while (!isConvergent1(dp, half_subset_size, half_subset_size) && iter <= max_iter_time)
//	{
//		// 确定变形子区点的"图像坐标"
//		ShapeFunction1(refer_x0, refer_y0, delta_x, delta_y, deform_warp, deform_subset_x, deform_subset_y);
//
//		// 变形子区点的灰度
//		deform_intensities = deform_image_interp->Values(deform_subset_x, deform_subset_y);
//		const double gg = NormalizeVectoerize(deform_intensities, normalized_deform_intensities);
//
//		// 计算伪逆矩阵
//		deform_gradient_x = deform_image_interp->GradientsX(deform_subset_x, deform_subset_y);
//		deform_gradient_y = deform_image_interp->GradientsY(deform_subset_x, deform_subset_y);
//		ForwardPseudoMatrix1(deform_gradient_x, deform_gradient_y, delta_x, delta_y, pseudo);
//
//		// 归一化的"参考子区"和"变形子区"的灰度差
//		const vec diff = normalized_refer_intensities - normalized_deform_intensities;
//
//		// 迭代增量
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
//		// 迭代次数加1
//		++iter;
//	}
//
//	// 判断迭代是否收敛
//	bool is_convergent = true;
//
//	if (iter >= max_iter_time)
//	{
//		std::cerr << "The iteration time " << iter << " reach the maximum!" << std::endl;
//		is_convergent = false;
//	}
//
//	// 最终的相关系数
//	ShapeFunction1(refer_x0, refer_y0, delta_x, delta_y, deform_warp, deform_subset_x, deform_subset_y);
//	deform_intensities = deform_image_interp->Values(deform_subset_x, deform_subset_y);
//
//	zncc = ZNCC(refer_intensities, deform_intensities);
//
//	// 判断相关系数是否满足阈值
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
	// "子区点"相对于"子区中心"的坐标
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

	// 将参考子区灰度归一化
	vec normalized_refer_intensities;
	NormalizeVectoerize(refer_intensities, normalized_refer_intensities);

	// 迭代次数
	iter = 0;

	mat deform_subset_x, deform_subset_y;		// 变形子区点的图像坐标
	mat deform_intensities;						// 变形子区点的灰度
	mat deform_gradient_x, deform_gradient_y;	// 变形子区点的灰度梯度
	mat pseudo;									// 伪逆矩阵
	vec normalized_deform_intensities;			// 归一化的变形子区灰度

	// 迭代
	vec dp = ones<vec>(SECOND_ORDER_SHAPE_FUNC_PARA_COUNT);	// 迭代增量
	while (!isConvergent2(dp, half_subset_size, half_subset_size) && iter <= max_iter_time)
	{
		// 确定变形子区点的"图像坐标"
		ShapeFunction2(refer_x0, refer_y0, delta_x, delta_y, deform_warp, deform_subset_x, deform_subset_y);

		// 变形子区点的灰度
		deform_intensities = deform_image_interp->Values(deform_subset_x, deform_subset_y);
		const double gg = NormalizeVectoerize(deform_intensities, normalized_deform_intensities);

		// 计算伪逆矩阵
		deform_gradient_x = deform_image_interp->GradientsX(deform_subset_x, deform_subset_y);
		deform_gradient_y = deform_image_interp->GradientsY(deform_subset_x, deform_subset_y);
		ForwardPseudoMatrix2(deform_gradient_x, deform_gradient_y, delta_x, delta_y, pseudo);

		// 归一化的"参考子区"和"变形子区"的灰度差
		const vec diff = normalized_refer_intensities - normalized_deform_intensities;

		// 迭代增量
		dp = gg * pseudo * diff;
		deform_warp = deform_warp + dp;

		// 迭代次数加1
		++iter;
	}


	// 判断迭代是否收敛
	bool is_convergent = true;

	if (iter >= max_iter_time)
	{
		std::cerr << "The iteration time " << iter << " reach the maximum!" << std::endl;
		is_convergent = false;
	}

	// 最终的相关系数
	ShapeFunction2(refer_x0, refer_y0, delta_x, delta_y, deform_warp, deform_subset_x, deform_subset_y);
	deform_intensities = deform_image_interp->Values(deform_subset_x, deform_subset_y);

	zncc = ZNCC(refer_intensities, deform_intensities);

	// 判断相关系数是否满足阈值
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