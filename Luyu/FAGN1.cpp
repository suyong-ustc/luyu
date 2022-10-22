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