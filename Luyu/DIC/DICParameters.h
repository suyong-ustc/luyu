#pragma once
#include <QPoint>
#include <QRect>
#include <QSize>
#include <armadillo>


class DICParameters
{
public:
	DICParameters();
	~DICParameters();


	/****************************************************************************/
	/*******************                读取参数           **********************/
	/****************************************************************************/


	QRect roi() const { return roi_; }


	int grid_step() const { return grid_step_; }


	int subset_size() const { return subset_size_; }


	int half_subset_size() const { return (subset_size_ - 1) / 2; }


	double zncc_threshold() const { return zncc_threshold_; }


	int max_iteration_times() const { return max_iteration_times_; }


	double error_threshold() const { return error_threshold_; }


	QPoint Seed() const { return seed_; }


	int search_radius() const { return search_radius_; }


	int bspline_interpolation_order() const {return bspline_interpolation_order_; }


	int shape_function_order() const { return shape_function_order_; }



	/****************************************************************************/
	/*******************               设置参数            **********************/
	/****************************************************************************/


	void set_roi(const QRect& roi) { roi_ = roi; }


	void set_roi(const int& xmin, const int& xmax, const int& ymin, const int& ymax) {
		roi_.setLeft(xmin);
		roi_.setRight(xmax);
		roi_.setTop(ymin);
		roi_.setBottom(ymax);
	}


	void set_grid_step(const int& grid_step) { grid_step_ = grid_step; }


	void set_subset_size(const int& subset_size) { subset_size_ = subset_size; }


	void set_zncc_threshold(const double& zncc_threshold) { zncc_threshold_ = zncc_threshold; }


	void set_error_threshold(const double& error_threshold) { error_threshold_ = error_threshold; }


	void set_max_iteration_times(const int& max_iteration_times) { max_iteration_times_ = max_iteration_times; }


	void set_seed(const QPoint& seed) { seed_ = seed; }


	void set_seed(const int& x, const int& y) { 
		seed_.setX(x);
		seed_.setY(y);
	}


	void set_search_radius(const int& search_radius) { search_radius_ = search_radius; }


	void set_bspline_interpolation_order(const int& bspline_interpolation_order) { bspline_interpolation_order_ = bspline_interpolation_order; }


	void set_shape_function_order(const int& shape_function_order) { shape_function_order_ = shape_function_order; }



	/****************************************************************************/
	/*******************               其它函数            **********************/
	/****************************************************************************/

	bool grid(arma::mat& x, arma::mat& y) const;


private:

	QRect roi_;

	int grid_step_;

	int subset_size_;

	double zncc_threshold_;

	int max_iteration_times_;

	double error_threshold_;

	QPoint seed_;

	int search_radius_;

	int bspline_interpolation_order_;

	int shape_function_order_;

};

