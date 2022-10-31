#include <iostream>
#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <armadillo>
#include "DIC/DICAlgorithm.h"
#include "DIC/DICParameters.h"
using namespace arma;



void SetDICParameters(DICParameters& dic_parameters)
{
	// ����Ȥ����
	dic_parameters.set_roi(25, 475, 25, 475);
	
	// ������
	dic_parameters.set_grid_step(5);

	// �����ߴ�
	dic_parameters.set_subset_size(19);

	// ZNCC��ֵ
	dic_parameters.set_zncc_threshold(0.8);

	// ����������
	dic_parameters.set_max_iteration_times(10);

	// �����ֵ
	dic_parameters.set_error_threshold(2e-4);

	// ��ֵ����
	dic_parameters.set_bspline_interpolation_order(5);

	// �κ�������
	dic_parameters.set_shape_function_order(5);
}



bool ReadImage(const std::string& image_path, mat& image)
{
	std::cout << "Import image with path " << image_path << std::endl;

	// ��ȡͼ��
	cv::Mat cvmat = cv::imread(image_path, cv::IMREAD_GRAYSCALE);

	if (cvmat.empty())
	{
		std::cerr << "I can not import the image!" << std::endl;
		return false;
	}

	// �� opencv �����ݸ�ʽת��Ϊ armadillo ����
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
	// ��ȡͼ��
	const std::string prefix("..\\images\\a_");
	const std::string refer_image_path(prefix + "0.bmp");
	const std::string deform_image_path(prefix + "1.bmp");

	mat refer_image;
	mat deform_image;
	ReadImage(refer_image_path, refer_image);
	ReadImage(deform_image_path, deform_image);

	// �������
	DICParameters dic_parameters;
	SetDICParameters(dic_parameters);

	// ��ؼ���
	DICOutput* dic_output = RigsterFullFieldDisplacement(refer_image, deform_image, dic_parameters);

	dic_output->write_grid_coordinate("x.csv", "y.csv", arma::csv_ascii);
	dic_output->write_displacement_field("u.csv", "v.csv", arma::csv_ascii);
	dic_output->write_iteration_times("iter_times.csv", arma::csv_ascii);
	dic_output->write_zncc("zncc.csv", arma::csv_ascii);
	dic_output->write_valid_sign("valid_sign.csv", arma::csv_ascii);

	delete dic_output;

	return 0;

}