#pragma once
#include <QPoint>
#include <armadillo>


class DICDeformation;
class DICGrid;
class DICImage;
class DICImageGradient;
class DICParameters;
class Interpolator;


// 计算全场变形 

DICDeformation* RegisterFullFieldDeformation(
	const DICImage& refer_image, 
	const DICImage& deform_image, 
	const DICParameters& dic_parameters);



// 整像素搜索 

bool IntegerPixelSearch(
	const DICImage& refer_image, 
	const DICImage& deform_image, 
	const DICParameters& dic_parameters,
	const QPoint& refer_image_pos, 
	QPoint& deform_image_integer_pos);



// 亚像素迭代 

bool SubpixelInterate(
	const DICImage& refer_image,
	const DICImageGradient& refer_image_gradient,
	Interpolator* deform_image_interp,
	const DICParameters& dic_parameters,
	const double& refer_x0,
	const double& refer_y0,
	const double& deform_x0_est,
	const double& deform_y0_est,
	arma::vec& deform_warp,
	double& zncc,
	arma::uword& iter);



// 种子点搜索 

bool SeedPointSearch(
	const DICImage& refer_image,
	const DICImageGradient& refer_image_gradient,
	const DICImage& deform_image,
	Interpolator* deform_image_interp,
	const DICGrid& dic_grid,
	const DICParameters& dic_parameter,
	QPoint& seed_grid_id,
	arma::vec& seed_deform_warp_function,
	double& zncc,
	arma::uword& iter);



// 种子点扩散 

bool SeedPointSpread(
	const DICImage& refer_image, 
	const DICImageGradient& refer_image_gradient, 
	Interpolator* deform_image_interp,
	const DICGrid& dic_grid, 
	const DICParameters& dic_parameter, 
	const QPoint& seed_grid_id, 
	DICDeformation* dic_deformation);