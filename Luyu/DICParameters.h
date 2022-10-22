#pragma once
#include <QPoint>
#include <QRect>
#include <QSize>
#include "DICTypes.h"


class DICParameters
{
public:
	DICParameters();
	~DICParameters();

	/****************************************************************************/
	/*******************           Read parameters         **********************/
	/****************************************************************************/

	QSize ImageSize() const { 
		return image_size_; 
	}


	int ImageHeight() const { 
		return image_size_.height(); 
	}


	int ImageWidth() const { 
		return image_size_.width(); 
	}


	QRect ROI() const { 
		return roi_; 
	}


	int GridStep() const { 
		return grid_step_; 
	}


	int SubsetSize() const { 
		return subset_size_; 
	}


	int HalfSubsetSize() const { 
		return (subset_size_ - 1) / 2; 
	}


	double ZNCCThreshold() const { 
		return zncc_threshold_; 
	}


	int MaxIterationTimes() const { 
		return max_iteration_times_; 
	}


	double ErrorThreshold() const {
		return error_threshold_;
	}


	QPoint Seed() const { 
		return seed_; 
	}


	int SearchRadius() const { 
		return search_radius_; 
	}


	InterpType InterpolationAlgorithm() const {
		return interp_type_;
	}


	GradType GradientEstimator() const {
		return grad_type_;
	}


	ShapeFuncType ShapeFunction() const {
		return shape_func_type_;
	}


	RegistrationType RegistrationStrategy() const {
		return registration_type_;
	}


	/****************************************************************************/
	/*******************            Set parameters         **********************/
	/****************************************************************************/


	void setImageSize(const QSize& image_size) { 
		image_size_ = image_size; 
	}


	void setImageSize(const int& height, const int& width) {
		image_size_.setHeight(height);
		image_size_.setWidth(width);
	}


	void setROI(const QRect& roi) {
		roi_ = roi;
	}


	void setROI(const int& xmin, const int& xmax, const int& ymin, const int& ymax) {
		roi_.setLeft(xmin);
		roi_.setRight(xmax);
		roi_.setTop(ymin);
		roi_.setBottom(ymax);
	}


	void setGridStep(const int& grid_step) {
		grid_step_ = grid_step;
	}


	void setSubsetSize(const int& subset_size) {
		subset_size_ = subset_size;
	}


	void setZNCCThreshold(const double& zncc_threshold) {
		zncc_threshold_ = zncc_threshold;
	}


	void setErrorThreshold(const double& error_threshold) {
		error_threshold_ = error_threshold;
	}


	void setMaxIterationTimes(const int& max_iteration_times) {
		max_iteration_times_ = max_iteration_times;
	}


	void setSeed(const QPoint& seed) {
		seed_ = seed;
	}


	void setSeed(const int& x, const int& y) {
		seed_.setX(x);
		seed_.setY(y);
	}


	void setSearchRadius(const int& search_radius) {
		search_radius_ = search_radius;
	}


	void setInterpolationAlgorithm(const InterpType& interp_type) {
		interp_type_ = interp_type;
	}


	void setGradientEstimator(const GradType& grad_type) {
		grad_type_ = grad_type;
	}


	void setShapeFunction(const ShapeFuncType& shape_func_type) {
		shape_func_type_ = shape_func_type;
	}


	void setRegistrationStrategy(const RegistrationType& registration_type) {
		registration_type_ = registration_type;
	}


private:

	QSize image_size_;
	
	QRect roi_;

	int grid_step_;

	int subset_size_;

	double zncc_threshold_;

	int max_iteration_times_;

	double error_threshold_;

	QPoint seed_;

	int search_radius_;

	InterpType interp_type_;

	GradType grad_type_;

	ShapeFuncType shape_func_type_;

	RegistrationType registration_type_;

};

