//#pragma once
//#include <armadillo>
//#include <QPoint>
//#include "DICGrid.h"
//#include "DICTypes.h"
//
//
//
//class DICDeformation
//{
//public:
//	DICDeformation(const DICGrid& dic_grid, const ShapeFuncType& shape_function);
//	~DICDeformation();
//
//	// 读取成员变量
//	arma::vec WarpFunction(const int& r, const int& c) const
//	{
//		return deformation_field_.tube(r, c);
//	}
//
//	arma::mat U() const
//	{
//		return deformation_field_.slice(0);
//	}
//
//	arma::mat V() const
//	{
//		return (shape_function_ == SHAPEFUNC_FIRST) ? deformation_field_.slice(3) : deformation_field_.slice(6);
//	}
//
//	arma::uword ValidSign(const int& r, const int& c) const
//	{
//		return valid_(r, c);
//	}
//
//	arma::uword IterationTimes(const int& r, const int& c) const
//	{
//		return iter_(r, c);
//	}
//
//	double ZNCC(const int& r, const int& c) const
//	{
//		return zncc_(r, c);
//	}
//
//	ShapeFuncType ShapeFunction() const
//	{
//		return shape_function_;
//	}
//
//	arma::uword Rows() const
//	{
//		return deformation_field_.n_rows;
//	}
//
//	arma::uword Cols() const
//	{
//		return deformation_field_.n_cols;
//	}
//
//
//	// 设置成员变量
//	void SetWarpFunction(const int& r, const int& c, const arma::vec& warp)
//	{
//		deformation_field_.tube(r, c) = warp;
//	}
//
//	void SetValidSign(const int& r, const int& c, const arma::uword& valid_sign)
//	{
//		valid_(r, c) = valid_sign;
//	}
//
//	void SetIterationTimes(const int& r, const int& c, const arma::uword& iter_times)
//	{
//		iter_(r, c) = iter_times;
//	}
//
//	void SetZNCC(const int& r, const int& c, const double& zncc)
//	{
//		zncc_(r, c) = zncc;
//	}
//
//
//	// 存储到文件
//	bool SaveDeformationField(const std::string& u_file_name, const std::string& v_file_name, const arma::file_type& type) const;
//	bool SaveImageCoordinate(const std::string& x_file_name, const std::string& v_file_name, const arma::file_type& type) const;
//	bool SaveValidSign(const std::string& file_name, const arma::file_type& type) const;
//	bool SaveIterationTimes(const std::string& file_name, const arma::file_type& type) const;
//	bool SaveZNCC(const std::string& file_name, const arma::file_type& type) const;
//
//private:
//	ShapeFuncType shape_function_;	// 形函数阶数
//
//	arma::cube deformation_field_;	// 变形场
//	arma::umat valid_;				// 有效位
//	arma::umat iter_;				// 迭代次数
//	arma::mat zncc_;				// 相关系数
//
//	arma::imat x_;
//	arma::imat y_;
//
//
//public:
//	static const arma::uword POI_OUT_OF_AOI = 0;				// 计算点不在AOI内
//	static const arma::uword POI_NEED_TO_CALCULATE = 1;			// 需要计算的点
//	static const arma::uword POI_SUCESSFULLY_CALCULATED = 2;	// 成功计算的点
//	static const arma::uword POI_FAIL_TO_CALCULATE = 3;			// 计算失败的点
//};
//
