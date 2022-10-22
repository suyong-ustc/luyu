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
//	// 参考图灰度梯度
//	DICImageGradient refer_image_gradient(refer_image, dic_parameters.GradientEstimator());
//
//	// 计算变形图插值系数
//	Interpolator* deform_image_interp = deform_image.ConstructInterpolation(dic_parameters.InterpolationAlgorithm());
//
//	// 初始化变形场
//	DICGrid dic_grid(dic_parameters.ROI(), dic_parameters.GridStep());
//	DICDeformation* dic_deformation = new DICDeformation(dic_grid, dic_parameters.ShapeFunction());
//
//	// 种子点搜索
//	QPoint seed_grid_id;		// 种子点 网格坐标
//	vec seed_warp_function;		// 种子点 形函数参数
//	double zncc;				// 种子点 相关系数
//	uword iter;					// 种子点 迭代次数
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
//	// 检验种子点搜索是否成功
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
//	// 种子点扩散
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
///******************************                  整像素搜索            **********************************/
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
//	// 检查给定点是否在"计算区域"中
//	if (!dic_parameters.ROI().contains(refer_image_pos))
//	{
//		std::cerr << "I can not conduct integer pixel search because the seed point is not in the region of interest!" << std::endl;
//		return false;
//	}
//
//	// 参考子区灰度
//	const mat refer_subset_intensity = refer_image.Intensities(refer_image_pos.x(), refer_image_pos.y(), dic_parameters.HalfSubsetSize());
//
//
//	// 整像素搜索
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
//			// 变形子区灰度
//			const mat deform_subset_intensity = deform_image.Intensities(x, y, dic_parameters.HalfSubsetSize());
//
//			// 相关系数
//			const double zncc = ZNCC(refer_subset_intensity, deform_subset_intensity);
//
//			// 更新位置
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
///******************************                  亚像素迭代            **********************************/
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
//	// 子区点位置
//	const int half_subset_size = dic_parameters.HalfSubsetSize();
//
//
//	// 参考图灰度
//	const mat refer_intensities = refer_image.Intensities(refer_x0, refer_y0, half_subset_size);
//
//	// 迭代初值
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
//	// 反向算法
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
//	// 正向算法
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
///******************************                种子点整像素搜索         **********************************/
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
//	// 将种子点移到最邻近网格点
//	QPoint seed_refer_image_pos = dic_parameters.Seed();
//	seed_grid_id = dic_grid.NearestGridPoint(seed_refer_image_pos.x(), seed_refer_image_pos.y());
//
//	// 种子点图像坐标
//	seed_refer_image_pos = dic_grid.ImageCoordinate(seed_grid_id);
//
//
//	// 种子点整像素搜索
//	QPoint seed_deform_image_integer_pos;
//	IntegerPixelSearch(
//		refer_image, 
//		deform_image, 
//		dic_parameters, 
//		seed_refer_image_pos, 
//		seed_deform_image_integer_pos);
//
//
//	// 种子点亚像素迭代
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
///******************************                 种子点扩散              **********************************/
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
//	// 种子点信息
//	QQueue<QPoint > seed_queue;
//	seed_queue.enqueue(seed_grid_id);
//
//
//	// 种子点扩散
//	while (!seed_queue.empty()) {
//
//		// 当前的种子点
//		QPoint current_seed_grid_id = seed_queue.front();
//		QPoint current_seed_refer_image_pos = dic_grid.ImageCoordinate(current_seed_grid_id);
//		const vec& current_seed_warp_function = dic_deformation->WarpFunction(current_seed_grid_id.y(), current_seed_grid_id.x());
//
//
//		// 下一个种子点
//		QPoint next_seed_grid_id;
//		for (int a = 0; a < 4; a++) {
//
//			if (a == 3)			// 左
//			{
//				next_seed_grid_id = current_seed_grid_id + QPoint(-1, 0);
//			}
//			else if (a == 1)	// 右
//			{
//				next_seed_grid_id = current_seed_grid_id + QPoint(1, 0);
//			}
//			else if (a == 0)	// 上
//			{
//				next_seed_grid_id = current_seed_grid_id + QPoint(0, -1);
//			}
//			else if (a == 2)	// 下
//			{
//				next_seed_grid_id = current_seed_grid_id + QPoint(0, 1);
//			}
//	
//
//			// 保证网格点在计算区域内
//			const int r = next_seed_grid_id.ry();
//			const int c = next_seed_grid_id.rx();
//			if (r < 0 || c < 0 || r >= dic_grid.GridRows() || c >= dic_grid.GridCols())
//				continue;
//
//			// 确保未曾计算过
//			if (dic_deformation->ValidSign(r, c) != DICDeformation::POI_NEED_TO_CALCULATE)
//				continue;
//
//
//			// 估计图像坐标的初值
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
//			// 迭代
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
//			// 判断是否收敛
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
//			// 添加到种子点队列
//			QPoint neighbour(c, r);
//			seed_queue.enqueue(neighbour);
//		}
//
//		seed_queue.dequeue();
//	}
//
//	return true;
//}
