#include <iostream>
#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <armadillo>
#include "DIC/DICAlgorithm.h"
using namespace arma;


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

	// 划分网格
	RigsterFullFieldDisplacement(refer_image, deform_image, dic_parameters);

	return 0;

}