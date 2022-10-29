#include "DICOutput.h"
using namespace arma;


DICOutput::DICOutput(const mat& x, const mat& y, const int& shape_function_parameter_total_number)
{
	// 网格点图像坐标
	x_ = x;
	y_ = y;

	// 网格点的行数和列数
	const int rows = x.n_rows;
	const int cols = x.n_cols;

	// 变形场
	displacement_field_.zeros(rows, cols, shape_function_parameter_total_number);

	// 有效位
	valid_sign_.set_size(rows, cols);
	valid_sign_.fill(POI_NEED_CALCULATE);

	// 迭代次数
	iteration_times_.zeros(rows, cols);

	// 相关系数
	zncc_.set_size(rows, cols);
	zncc_.fill(-10);

}


DICOutput::~DICOutput()
{

}