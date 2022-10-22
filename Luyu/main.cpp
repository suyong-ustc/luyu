#include <iostream>
#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <armadillo>
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



bool MeshGrid(const int& xmin, const int& xmax, const int& ymin, const int& ymax, const int& grid_step, mat& x, mat& y)
{
	// 检查数据
	if (xmin >= xmax || ymin >= ymax)
	{
		std::cerr << "Invalid ROI!" << std::endl;
		return false;
	}

	// 矩阵尺寸
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
	// 初始化
	u.zeros(x.n_rows, x.n_cols);
	v.zeros(x.n_rows, x.n_cols);

	// 赋初值
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
	// 读取图像
	const std::string prefix("..\\Data\\a1-2_");
	const std::string refer_image_path(prefix + "0.bmp");
	const std::string deform_image_path(prefix + "1.bmp");

	mat refer_image;
	mat deform_image;
	ReadImage(refer_image_path, refer_image);
	ReadImage(deform_image_path, deform_image);

	// 划分网格
	const int xmin = 50;
	const int xmax = 450;
	const int ymin = 50;
	const int ymax = 450;

	const int grid_step = 10;

	mat x;
	mat y;
	MeshGrid(xmin, xmax, ymin, ymax, grid_step, x, y);

	// 估计初值
	mat u;
	mat v;
	EstimateInitialDisplacement(refer_image, deform_image, x, y, u, v);

	// 相关计算
	mat zncc;
	RegisterSubpixelDisplacement(refer_image, deform_image, x, y, u, v, zncc);

	std::cout << "generate sensor noise ..." << std::endl;

	const double sigma = 1;			// 图像噪声标准差
	const int nnoise_add = 100;		// 噪声添加次数

	const cube refer_sensor_noise = sigma * randn(dic_parameters.ImageHeight(), dic_parameters.ImageWidth(), nnoise_add);
	const cube deform_sensor_noise = sigma * randn(dic_parameters.ImageHeight(), dic_parameters.ImageWidth(), nnoise_add);


	//  不同的子区
	const ivec subset_sizes = regspace<ivec>(21, 10, 81);

	for (int i = 0; i < subset_sizes.n_elem; ++i)
	{
		// 子区尺寸
		const int subset_size = subset_sizes(i);
		dic_parameters.setSubsetSize(subset_size);

		std::cout << "analyze subset " << subset_size << std::endl;

		// 网格尺寸
		DICGrid dic_grid(dic_parameters.ROI(), dic_parameters.GridStep());

		const int grid_rows = dic_grid.GridRows();
		const int grid_cols = dic_grid.GridCols();

		// 有噪声情况下的位移场
		cube u_with_noise = zeros<cube>(grid_rows, grid_cols, nnoise_add);

		const int tick = nnoise_add / 10;
		int num_calculated = 0;

#	pragma omp parallel for num_threads(15)
		for (int k = 0; k < nnoise_add; ++k)
		{
			// 显示计算百分比
#			pragma omp critical
			++num_calculated;
			if (num_calculated % tick == 0)
				printf("%.2f%%\t", 100. * num_calculated / nnoise_add);

			// 有噪声的参考图和变形图
			const DICImage noisy_refer_image = refer_image.AddNoise(refer_sensor_noise.slice(k));
			const DICImage noisy_deform_image = deform_image.AddNoise(deform_sensor_noise.slice(k));

			// 计算变形场
			DICDeformation* dic_deformation_filed = RegisterFullFieldDeformation(noisy_refer_image, noisy_deform_image, dic_parameters);

			// 存储位移场U
			u_with_noise.slice(k) = dic_deformation_filed->U();

			delete dic_deformation_filed;
		}

		// 计算每个"网格点"位移的均值和标准差
		std::cout << "\ncalculate mean and standard deviation ..." << std::endl;

		mat u_mean = zeros<mat>(grid_rows, grid_cols);	// 位移均值
		mat u_std = zeros<mat>(grid_rows, grid_cols);	// 位移标准差

		for (int r = 0; r < grid_rows; r++)
		{
			for (int c = 0; c < grid_cols; c++)
			{
				const vec t = u_with_noise.tube(r, c);

				u_mean(r, c) = mean(t);
				u_std(r, c) = stddev(t);
			}
		}

		// 存储计算结果
		std::cout << "save the DIC results ..." << std::endl;

		std::string result_prefix = prefix + std::to_string(dic_parameters.SubsetSize()) + "_Sigma" + std::to_string(int(sigma))
			+ (dic_parameters.ShapeFunction() == SHAPEFUNC_FIRST ? "_First" : "_Second");

		u_mean.save(result_prefix + "_u_mean.csv", csv_ascii);
		u_std.save(result_prefix + "_u_std.csv", csv_ascii);

		// 结束此阶段分析
		std::cout << "finish analysis of subset size " << subset_size << std::endl;
	}


	// 总耗时
	clock_t end_clock = clock();
	int time_cost = (end_clock - start_clock) / CLOCKS_PER_SEC;
	printf("it totally costs %dh %dm %ds\n", time_cost / 3600, time_cost / 60 % 60, time_cost % 60);

	return true;

}