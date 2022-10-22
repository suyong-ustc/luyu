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