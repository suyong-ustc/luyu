//#include "DICRegisterFullFieldDeformation.h"
//#include <QPoint>
//#include <QQueue>
//#include "DICDeformation.h"
//#include "DICGrid.h"
//#include "DICImage.h"
//#include "DICImageGradient.h"
//#include "DICParameters.h"
//#include "DICShapeFunction.h"
//#include "FAGN1.h"
//#include "FAGN2.h"
//#include "ICGN1.h"
//#include "ICGN2.h"
//#include "Interpolator.h"
//#include "Util.h"
//using namespace arma;
//
//
//
//DICDeformation* RegisterFullFieldDeformation(const DICImage& refer_image, const DICImage& deform_image, const DICParameters& dic_parameters)
//{
//	// �ο�ͼ�Ҷ��ݶ�
//	DICImageGradient refer_image_gradient(refer_image, dic_parameters.GradientEstimator());
//
//	// �������ͼ��ֵϵ��
//	Interpolator* deform_image_interp = deform_image.ConstructInterpolation(dic_parameters.InterpolationAlgorithm());
//
//	// ��ʼ�����γ�
//	DICGrid dic_grid(dic_parameters.ROI(), dic_parameters.GridStep());
//	DICDeformation* dic_deformation = new DICDeformation(dic_grid, dic_parameters.ShapeFunction());
//
//	// ���ӵ�����
//	QPoint seed_grid_id;		// ���ӵ� ��������
//	vec seed_warp_function;		// ���ӵ� �κ�������
//	double zncc;				// ���ӵ� ���ϵ��
//	uword iter;					// ���ӵ� ��������
//	bool is_convergent =
//		SeedPointSearch(
//			refer_image, 
//			refer_image_gradient, 
//			deform_image, 
//			deform_image_interp, 
//			dic_grid, 
//			dic_parameters, 
//			seed_grid_id, 
//			seed_warp_function,
//			zncc,
//			iter);
//
//	// �������ӵ������Ƿ�ɹ�
//	if (!is_convergent)
//	{
//		std::cerr << "seed point search has failed! " << std::endl;
//		delete dic_deformation;
//		return NULL;
//	}
//
//	dic_deformation->SetWarpFunction(seed_grid_id.y(), seed_grid_id.x(), seed_warp_function);
//	dic_deformation->SetValidSign(seed_grid_id.y(), seed_grid_id.x(), DICDeformation::POI_SUCESSFULLY_CALCULATED);
//	dic_deformation->SetZNCC(seed_grid_id.y(), seed_grid_id.x(), zncc);
//	dic_deformation->SetIterationTimes(seed_grid_id.y(), seed_grid_id.x(), iter);
//
//	// ���ӵ���ɢ
//	SeedPointSpread(
//		refer_image, 
//		refer_image_gradient, 
//		deform_image_interp, 
//		dic_grid, 
//		dic_parameters, 
//		seed_grid_id, 
//		dic_deformation);
//
//	return dic_deformation;
//}
//
//
//
//
///********************************************************************************************************/
///******************************                  ����������            **********************************/
///********************************************************************************************************/
//
//
//bool IntegerPixelSearch(
//	const DICImage& refer_image, 
//	const DICImage& deform_image, 
//	const DICParameters& dic_parameters,
//	const QPoint& refer_image_pos, 
//	QPoint& deform_image_integer_pos)
//{
//	// ���������Ƿ���"��������"��
//	if (!dic_parameters.ROI().contains(refer_image_pos))
//	{
//		std::cerr << "I can not conduct integer pixel search because the seed point is not in the region of interest!" << std::endl;
//		return false;
//	}
//
//	// �ο������Ҷ�
//	const mat refer_subset_intensity = refer_image.Intensities(refer_image_pos.x(), refer_image_pos.y(), dic_parameters.HalfSubsetSize());
//
//
//	// ����������
//	double zncc_max = -1.0;
//
//	const int xmin = refer_image_pos.x() - dic_parameters.SearchRadius();
//	const int xmax = refer_image_pos.x() + dic_parameters.SearchRadius();
//	const int ymin = refer_image_pos.y() - dic_parameters.SearchRadius();
//	const int ymax = refer_image_pos.y() + dic_parameters.SearchRadius();
//
//	for (int x = xmin; x <= xmax; ++x)
//	{
//		for (int y = ymin; y <= ymax; ++y)
//		{
//			// ���������Ҷ�
//			const mat deform_subset_intensity = deform_image.Intensities(x, y, dic_parameters.HalfSubsetSize());
//
//			// ���ϵ��
//			const double zncc = ZNCC(refer_subset_intensity, deform_subset_intensity);
//
//			// ����λ��
//			if (zncc > zncc_max) {
//				zncc_max = zncc;
//				deform_image_integer_pos.setX(x);
//				deform_image_integer_pos.setY(y);
//			}
//		}
//	}
//
//	return true;
//}
//
//
//
///********************************************************************************************************/
///******************************                  �����ص���            **********************************/
///********************************************************************************************************/
//
//
//bool SubpixelInterate(
//	const DICImage& refer_image,
//	const DICImageGradient& refer_image_gradient,
//	Interpolator* deform_image_interp,
//	const DICParameters& dic_parameters,
//	const double& refer_x0,
//	const double& refer_y0,
//	const double& deform_x0_est,
//	const double& deform_y0_est,
//	vec& deform_warp,
//	double& zncc,
//	uword& iter)
//{
//
//	// ������λ��
//	const int half_subset_size = dic_parameters.HalfSubsetSize();
//
//
//	// �ο�ͼ�Ҷ�
//	const mat refer_intensities = refer_image.Intensities(refer_x0, refer_y0, half_subset_size);
//
//	// ������ֵ
//	if (dic_parameters.ShapeFunction() == SHAPEFUNC_FIRST)
//	{
//		deform_warp.zeros(FIRST_ORDER_SHAPE_FUNC_PARA_COUNT);
//		deform_warp(0) = deform_x0_est - refer_x0;
//		deform_warp(3) = deform_y0_est - refer_y0;
//	}
//	else
//	{
//		deform_warp.zeros(SECOND_ORDER_SHAPE_FUNC_PARA_COUNT);
//		deform_warp(0) = deform_x0_est - refer_x0;
//		deform_warp(6) = deform_y0_est - refer_y0;
//	}
//
//
//	// �����㷨
//	bool is_ok = false;
//
//	if (dic_parameters.RegistrationStrategy() == REGISTRATION_INVERSE)
//	{
//		const mat refer_gradient_x = refer_image_gradient.GradientsX(refer_x0, refer_y0, half_subset_size);
//		const mat refer_gradient_y = refer_image_gradient.GradientsY(refer_x0, refer_y0, half_subset_size);
//
//		if (dic_parameters.ShapeFunction() == SHAPEFUNC_FIRST)
//		{
//			is_ok = 
//				InverseCompositionalGaussNewtonMethod1(
//					refer_intensities, 
//					refer_gradient_x, 
//					refer_gradient_y, 
//					refer_x0, 
//					refer_y0, 
//					deform_image_interp, 
//					half_subset_size,
//					dic_parameters.MaxIterationTimes(), 
//					dic_parameters.ZNCCThreshold(), 
//					deform_warp,
//					iter,
//					zncc);
//		} 
//		else
//		{
//			is_ok = 
//				InverseCompositionalGaussNewtonMethod2(
//					refer_intensities,
//					refer_gradient_x,
//					refer_gradient_y,
//					refer_x0,
//					refer_y0,
//					deform_image_interp,
//					half_subset_size,
//					dic_parameters.MaxIterationTimes(),
//					dic_parameters.ZNCCThreshold(),
//					deform_warp,
//					iter,
//					zncc);
//		}
//	}
//
//	// �����㷨
//	if (dic_parameters.RegistrationStrategy() == REGISTRATION_FORWARD)
//	{
//		if (dic_parameters.ShapeFunction() == SHAPEFUNC_FIRST)
//		{
//			is_ok =
//				ForwardAdditiveGaussNewtonMethod1(
//					refer_intensities,
//					refer_x0,
//					refer_y0,
//					deform_image_interp,
//					half_subset_size,
//					dic_parameters.MaxIterationTimes(),
//					dic_parameters.ZNCCThreshold(),
//					deform_warp,
//					iter,
//					zncc);
//		}
//		else
//		{
//			is_ok =
//				ForwardAdditiveGaussNewtonMethod2(
//					refer_intensities,
//					refer_x0,
//					refer_y0,
//					deform_image_interp,
//					half_subset_size,
//					dic_parameters.MaxIterationTimes(),
//					dic_parameters.ZNCCThreshold(),
//					deform_warp,
//					iter,
//					zncc);
//		}
//	}
//
//
//	return is_ok;
//}
//
//
//
///********************************************************************************************************/
///******************************                ���ӵ�����������         **********************************/
///********************************************************************************************************/
//
//
//bool SeedPointSearch(
//	const DICImage& refer_image, 
//	const DICImageGradient& refer_image_gradient, 
//	const DICImage& deform_image, 
//	Interpolator* deform_image_interp,
//	const DICGrid& dic_grid, 
//	const DICParameters& dic_parameters, 
//	QPoint& seed_grid_id, 
//	vec& seed_deform_warp_function, 
//	double& zncc,
//	uword& iter)
//{
//	// �����ӵ��Ƶ����ڽ������
//	QPoint seed_refer_image_pos = dic_parameters.Seed();
//	seed_grid_id = dic_grid.NearestGridPoint(seed_refer_image_pos.x(), seed_refer_image_pos.y());
//
//	// ���ӵ�ͼ������
//	seed_refer_image_pos = dic_grid.ImageCoordinate(seed_grid_id);
//
//
//	// ���ӵ�����������
//	QPoint seed_deform_image_integer_pos;
//	IntegerPixelSearch(
//		refer_image, 
//		deform_image, 
//		dic_parameters, 
//		seed_refer_image_pos, 
//		seed_deform_image_integer_pos);
//
//
//	// ���ӵ������ص���
//	bool is_convergent =
//		SubpixelInterate(
//			refer_image, 
//			refer_image_gradient, 
//			deform_image_interp, 
//			dic_parameters, 
//			seed_refer_image_pos.x(), 
//			seed_refer_image_pos.y(),
//			seed_deform_image_integer_pos.x(), 
//			seed_deform_image_integer_pos.y(), 
//			seed_deform_warp_function,
//			zncc,
//			iter);
//
//	return is_convergent;
//}
//
//
//
///********************************************************************************************************/
///******************************                 ���ӵ���ɢ              **********************************/
///********************************************************************************************************/
//
//
//
//bool SeedPointSpread(
//	const DICImage& refer_image, 
//	const DICImageGradient& refer_image_gradient, 
//	Interpolator* deform_image_interp,
//	const DICGrid& dic_grid, 
//	const DICParameters& dic_parameter, 
//	const QPoint& seed_grid_id, 
//	DICDeformation* dic_deformation)
//{
//	// ���ӵ���Ϣ
//	QQueue<QPoint > seed_queue;
//	seed_queue.enqueue(seed_grid_id);
//
//
//	// ���ӵ���ɢ
//	while (!seed_queue.empty()) {
//
//		// ��ǰ�����ӵ�
//		QPoint current_seed_grid_id = seed_queue.front();
//		QPoint current_seed_refer_image_pos = dic_grid.ImageCoordinate(current_seed_grid_id);
//		const vec& current_seed_warp_function = dic_deformation->WarpFunction(current_seed_grid_id.y(), current_seed_grid_id.x());
//
//
//		// ��һ�����ӵ�
//		QPoint next_seed_grid_id;
//		for (int a = 0; a < 4; a++) {
//
//			if (a == 3)			// ��
//			{
//				next_seed_grid_id = current_seed_grid_id + QPoint(-1, 0);
//			}
//			else if (a == 1)	// ��
//			{
//				next_seed_grid_id = current_seed_grid_id + QPoint(1, 0);
//			}
//			else if (a == 0)	// ��
//			{
//				next_seed_grid_id = current_seed_grid_id + QPoint(0, -1);
//			}
//			else if (a == 2)	// ��
//			{
//				next_seed_grid_id = current_seed_grid_id + QPoint(0, 1);
//			}
//	
//
//			// ��֤������ڼ���������
//			const int r = next_seed_grid_id.ry();
//			const int c = next_seed_grid_id.rx();
//			if (r < 0 || c < 0 || r >= dic_grid.GridRows() || c >= dic_grid.GridCols())
//				continue;
//
//			// ȷ��δ�������
//			if (dic_deformation->ValidSign(r, c) != DICDeformation::POI_NEED_TO_CALCULATE)
//				continue;
//
//
//			// ����ͼ������ĳ�ֵ
//			QPoint next_seed_refer_image_pos = dic_grid.ImageCoordinate(next_seed_grid_id);
//
//			double next_seed_deform_image_pos_x_est;
//			double next_seed_deform_image_pos_y_est;
//			ShapeFunction(
//				current_seed_refer_image_pos.x(),
//				current_seed_refer_image_pos.y(),
//				next_seed_refer_image_pos.x() - current_seed_refer_image_pos.x(),
//				next_seed_refer_image_pos.y() - current_seed_refer_image_pos.y(),
//				current_seed_warp_function,
//				dic_parameter.ShapeFunction(),
//				next_seed_deform_image_pos_x_est,
//				next_seed_deform_image_pos_y_est);
//			//next_seed_deform_image_pos_x_est = next_seed_refer_image_pos.x() + 0;
//			//next_seed_deform_image_pos_y_est = next_seed_refer_image_pos.y();
//
//			// ����
//			vec next_seed_warp_function;
//			double zncc;
//			uword iter;
//			bool is_converge = 
//				SubpixelInterate(
//				refer_image,
//				refer_image_gradient,
//				deform_image_interp,
//				dic_parameter,
//				next_seed_refer_image_pos.x(),
//				next_seed_refer_image_pos.y(),
//				next_seed_deform_image_pos_x_est,
//				next_seed_deform_image_pos_y_est,
//				next_seed_warp_function,
//				zncc,
//				iter);
//
//			// �ж��Ƿ�����
//			dic_deformation->SetWarpFunction(r, c, next_seed_warp_function);
//			dic_deformation->SetZNCC(r, c, zncc);
//			dic_deformation->SetIterationTimes(r, c, iter);
//
//			if (is_converge)
//			{
//				dic_deformation->SetValidSign(r, c, DICDeformation::POI_SUCESSFULLY_CALCULATED);
//			}
//			else
//			{
//				dic_deformation->SetValidSign(r, c, DICDeformation::POI_FAIL_TO_CALCULATE);
//				continue;
//			}
//
//
//			// ��ӵ����ӵ����
//			QPoint neighbour(c, r);
//			seed_queue.enqueue(neighbour);
//		}
//
//		seed_queue.dequeue();
//	}
//
//	return true;
//}
