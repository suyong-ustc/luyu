#include <iostream>
#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <armadillo>
using namespace arma;


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



bool MeshGrid(const int& xmin, const int& xmax, const int& ymin, const int& ymax, const int& grid_step, mat& x, mat& y)
{
	// �������
	if (xmin >= xmax || ymin >= ymax)
	{
		std::cerr << "Invalid ROI!" << std::endl;
		return false;
	}

	// ����ߴ�
	const int nx = 1 + (xmax - xmin) / grid_step;
	const int ny = 1 + (ymax - ymin) / grid_step;

	x.zeros(ny, nx);
	y.zeros(ny, nx);

	for (int r = 0; r < ny; ++r)
	{
		for (int c = 0; c < nx; ++c)
		{
			x(r, c) = xmin + c * grid_step;
			y(r, c) = ymin + r * grid_step;
		}

	}

	return true;
}



bool EstimateInitialDisplacement(const mat& refer_image, const mat& deform_image, const mat& x, const mat& y, mat& u, mat& v)
{
	// ��ʼ��
	u.zeros(x.n_rows, x.n_cols);
	v.zeros(x.n_rows, x.n_cols);

	// ����ֵ
	bool is_harmonic = true;

	if (is_harmonic)
	{
		double a = 5;
		double T = 20;
		double b = 0;

		u = a * sin(2 * datum::pi * x / T + b);
	}
	else
	{
		double a = 5;
		double x0 = 150;
		double c = 20;

		mat t = (x - x0) / c;
		u = a * exp(-t % t);
	}

	return true;
}



bool RegisterSubpixelDisplacement(const mat& refer_image, const mat& deform_image, const mat& x, const mat& y, mat& u, mat& v, mat& zncc);
{





}



int main()
{
	// ��ȡͼ��
	const std::string prefix("..\\Data\\a1-2_");
	const std::string refer_image_path(prefix + "0.bmp");
	const std::string deform_image_path(prefix + "1.bmp");

	mat refer_image;
	mat deform_image;
	ReadImage(refer_image_path, refer_image);
	ReadImage(deform_image_path, deform_image);

	// ��������
	const int xmin = 50;
	const int xmax = 450;
	const int ymin = 50;
	const int ymax = 450;

	const int grid_step = 10;

	mat x;
	mat y;
	MeshGrid(xmin, xmax, ymin, ymax, grid_step, x, y);

	// ���Ƴ�ֵ
	mat u;
	mat v;
	EstimateInitialDisplacement(refer_image, deform_image, x, y, u, v);

	// ��ؼ���
	mat zncc;
	RegisterSubpixelDisplacement(refer_image, deform_image, x, y, u, v, zncc);

	std::cout << "generate sensor noise ..." << std::endl;

	const double sigma = 1;			// ͼ��������׼��
	const int nnoise_add = 100;		// ������Ӵ���

	const cube refer_sensor_noise = sigma * randn(dic_parameters.ImageHeight(), dic_parameters.ImageWidth(), nnoise_add);
	const cube deform_sensor_noise = sigma * randn(dic_parameters.ImageHeight(), dic_parameters.ImageWidth(), nnoise_add);


	//  ��ͬ������
	const ivec subset_sizes = regspace<ivec>(21, 10, 81);

	for (int i = 0; i < subset_sizes.n_elem; ++i)
	{
		// �����ߴ�
		const int subset_size = subset_sizes(i);
		dic_parameters.setSubsetSize(subset_size);

		std::cout << "analyze subset " << subset_size << std::endl;

		// ����ߴ�
		DICGrid dic_grid(dic_parameters.ROI(), dic_parameters.GridStep());

		const int grid_rows = dic_grid.GridRows();
		const int grid_cols = dic_grid.GridCols();

		// ����������µ�λ�Ƴ�
		cube u_with_noise = zeros<cube>(grid_rows, grid_cols, nnoise_add);

		const int tick = nnoise_add / 10;
		int num_calculated = 0;

#	pragma omp parallel for num_threads(15)
		for (int k = 0; k < nnoise_add; ++k)
		{
			// ��ʾ����ٷֱ�
#			pragma omp critical
			++num_calculated;
			if (num_calculated % tick == 0)
				printf("%.2f%%\t", 100. * num_calculated / nnoise_add);

			// �������Ĳο�ͼ�ͱ���ͼ
			const DICImage noisy_refer_image = refer_image.AddNoise(refer_sensor_noise.slice(k));
			const DICImage noisy_deform_image = deform_image.AddNoise(deform_sensor_noise.slice(k));

			// ������γ�
			DICDeformation* dic_deformation_filed = RegisterFullFieldDeformation(noisy_refer_image, noisy_deform_image, dic_parameters);

			// �洢λ�Ƴ�U
			u_with_noise.slice(k) = dic_deformation_filed->U();

			delete dic_deformation_filed;
		}

		// ����ÿ��"�����"λ�Ƶľ�ֵ�ͱ�׼��
		std::cout << "\ncalculate mean and standard deviation ..." << std::endl;

		mat u_mean = zeros<mat>(grid_rows, grid_cols);	// λ�ƾ�ֵ
		mat u_std = zeros<mat>(grid_rows, grid_cols);	// λ�Ʊ�׼��

		for (int r = 0; r < grid_rows; r++)
		{
			for (int c = 0; c < grid_cols; c++)
			{
				const vec t = u_with_noise.tube(r, c);

				u_mean(r, c) = mean(t);
				u_std(r, c) = stddev(t);
			}
		}

		// �洢������
		std::cout << "save the DIC results ..." << std::endl;

		std::string result_prefix = prefix + std::to_string(dic_parameters.SubsetSize()) + "_Sigma" + std::to_string(int(sigma))
			+ (dic_parameters.ShapeFunction() == SHAPEFUNC_FIRST ? "_First" : "_Second");

		u_mean.save(result_prefix + "_u_mean.csv", csv_ascii);
		u_std.save(result_prefix + "_u_std.csv", csv_ascii);

		// �����˽׶η���
		std::cout << "finish analysis of subset size " << subset_size << std::endl;
	}


	// �ܺ�ʱ
	clock_t end_clock = clock();
	int time_cost = (end_clock - start_clock) / CLOCKS_PER_SEC;
	printf("it totally costs %dh %dm %ds\n", time_cost / 3600, time_cost / 60 % 60, time_cost % 60);

	return true;

}