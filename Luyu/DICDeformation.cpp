//#include "DICDeformation.h"
//using namespace arma;
//
//
//DICDeformation::DICDeformation(const DICGrid& dic_grid, const ShapeFuncType& shape_function)
//	: shape_function_(shape_function)
//{
//	const uword grid_rows = dic_grid.GridRows();
//	const uword grid_cols = dic_grid.GridCols();
//
//	// 形函数阶数
//	const uword n =	(shape_function == SHAPEFUNC_FIRST) ? FIRST_ORDER_SHAPE_FUNC_PARA_COUNT : SECOND_ORDER_SHAPE_FUNC_PARA_COUNT;
//
//	// 变形场
//	deformation_field_.zeros(grid_rows, grid_cols, n);
//
//	// 有效位
//	valid_.set_size(grid_rows, grid_cols);
//	valid_.fill(POI_NEED_TO_CALCULATE);
//
//	// 迭代次数
//	iter_.zeros(grid_rows, grid_cols);
//
//	// 相关系数
//	zncc_.set_size(grid_rows, grid_cols);
//	zncc_.fill(-10);
//
//	// 网格点图像坐标
//	x_ = dic_grid.X();
//	y_ = dic_grid.Y();
//}
//
//
//DICDeformation::~DICDeformation()
//{
//}
//
//
//
///*******************************************************************/
///*                        存储到文件                               */
///*******************************************************************/
//
//
//
//bool DICDeformation::SaveDeformationField(const std::string& u_file_name, const std::string& v_file_name, const file_type& type) const
//{
//	bool is_ok = U().save(u_file_name, type) &&  V().save(v_file_name, type);
//
//	if (!is_ok)
//		std::cerr << "I can not save the deformation field! Maybe the output file path is incorrect!" << std::endl;
//
//	return is_ok;
//}
//
//
//bool DICDeformation::SaveImageCoordinate(const std::string& x_file_name, const std::string& y_file_name, const arma::file_type& type) const
//{
//	bool is_ok = x_.save(x_file_name, type) && y_.save(y_file_name, type);
//
//	if (!is_ok)
//		std::cerr << "I can not save the image coordinates of the grids! Maybe the output file path is incorrect!" << std::endl;
//
//	return is_ok;
//}
//
//
//
//bool DICDeformation::SaveValidSign(const std::string& file_name, const file_type& type) const
//{
//	bool is_ok = valid_.save(file_name, type);
//
//	if (!is_ok)
//		std::cerr << "I can not save the valid sign! Maybe the output file path is incorrect!" << std::endl;
//
//	return is_ok;
//}
//
//
//bool DICDeformation::SaveZNCC(const std::string& file_name, const file_type& type) const
//{
//	bool is_ok = zncc_.save(file_name, type);
//
//	if (!is_ok)
//		std::cerr << "I can not save the correlation coefficients! Maybe the output file path is incorrect!" << std::endl;
//
//	return is_ok;
//}
//
//
//bool DICDeformation::SaveIterationTimes(const std::string& file_name, const file_type& type) const
//{
//	bool is_ok = iter_.save(file_name, type);
//
//	if (!is_ok)
//		std::cerr << "I can not save the correlation coefficients! Maybe the output file path is incorrect!" << std::endl;
//
//	return true;
//}