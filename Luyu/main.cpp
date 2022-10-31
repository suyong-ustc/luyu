#include <iostream>
#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <armadillo>
#include "DIC/DICAlgorithm.h"
#include "DIC/DICParameters.h"
using namespace arma;



void SetDICParameters(DICParameters& dic_parameters)
{
	// 感兴趣区域
	dic_parameters.set_roi(25, 475, 25, 475);
	
	// 网格间距
	dic_parameters.set_grid_step(5);

	// 子区尺寸
	dic_parameters.set_subset_size(19);

	// ZNCC阈值
	dic_parameters.set_zncc_threshold(0.8);

	// 最大迭代次数
	dic_parameters.set_max_iteration_times(10);

	// 误差阈值
	dic_parameters.set_error_threshold(2e-4);

	// 插值阶数
	dic_parameters.set_bspline_interpolation_order(5);

	// 形函数阶数
	dic_parameters.set_shape_function_order(5);
}



bool ReadImage(const std::string& image_path, mat& image)
{
	std::cout << "Import image with path " << image_path << std::endl;

	// 读取图像
	cv::Mat cvmat = cv::imread(image_path, cv::IMREAD_GRAYSCALE);

	if (cvmat.empty())
	{
		std::cerr << "I can not import the image!" << std::endl;
		return false;
	}

	// 将 opencv 的数据格式转化为 armadillo 矩阵
	image.zeros(cvmat.rows, cvmat.cols);

	for (int r = 0; r < cvmat.rows; ++r)
	{
		for (int c = 0; c < cvmat.cols; ++c)
		{
			image(r, c) = cvmat.at<uchar>(r, c);
		}

	}

	return true;
}



int main()
{
	// 读取图像
	const std::string prefix("..\\images\\a_");
	const std::string refer_image_path(prefix + "0.bmp");
	const std::string deform_image_path(prefix + "1.bmp");

	mat refer_image;
	mat deform_image;
	ReadImage(refer_image_path, refer_image);
	ReadImage(deform_image_path, deform_image);

	// 计算参数
	DICParameters dic_parameters;
	SetDICParameters(dic_parameters);

	// 相关计算
	DICOutput* dic_output = RigsterFullFieldDisplacement(refer_image, deform_image, dic_parameters);

	dic_output->write_grid_coordinate("x.csv", "y.csv", arma::csv_ascii);
	dic_output->write_displacement_field("u.csv", "v.csv", arma::csv_ascii);
	dic_output->write_iteration_times("iter_times.csv", arma::csv_ascii);
	dic_output->write_zncc("zncc.csv", arma::csv_ascii);
	dic_output->write_valid_sign("valid_sign.csv", arma::csv_ascii);

	delete dic_output;

	return 0;

}