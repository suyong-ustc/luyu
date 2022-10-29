#pragma once
#include <QPoint>
#include <armadillo>


class DICOutput
{
public:
	DICOutput(const arma::mat& x, const arma::mat& y, const int& shape_function_parameter_total_number);
	~DICOutput();


	/****************************************************************************/
	/*******************                读取结果           **********************/
	/****************************************************************************/

	arma::vec warp_function(const int& r, const int& c) const { return displacement_field_.tube(r, c); }

	arma::mat u() const { return displacement_field_.slice(0); }

	arma::mat v() const { return displacement_field_.slice(1); }

	double u(const int& r, const int& c) const { return displacement_field_(r, c, 0); }

	double v(const int& r, const int& c) const { return displacement_field_(r, c, 1); }


	arma::umat valid_sign() const { return valid_sign_; }

	arma::uword valid_sign(const int& r, const int& c) const { return valid_sign_(r, c); }


	arma::umat iteration_times() { return iteration_times_; }

	arma::uword iteration_times(const int& r, const int& c) const { return iteration_times_(r, c); }


	arma::mat zncc() const { return zncc_; }

	double zncc(const int& r, const int& c) const { return zncc_(r, c); }


	/****************************************************************************/
	/*******************                写入结果           **********************/
	/****************************************************************************/

	void set_warp_function(const int& r, const int& c, const arma::vec& warp) { displacement_field_.tube(r, c) = warp; }


	void set_valid_sign(const int& r, const int& c, const arma::uword& valid_sign) { valid_sign_(r, c) = valid_sign; }


	void set_iteration_times(const int& r, const int& c, const arma::uword& iter_times) { iteration_times_(r, c) = iter_times; }


	void set_zncc(const int& r, const int& c, const double& zncc) { zncc_(r, c) = zncc; }


	/****************************************************************************/
	/*******************                存储结果           **********************/
	/****************************************************************************/

	void write_displacement_field(const std::string& u_file_path, const std::string& v_file_path, const arma::file_type& type) const { u().save(u_file_path, type) && v().save(v_file_path, type); }


	bool write_grid_coordinate(const std::string& x_file_path, const std::string& y_file_path, const arma::file_type& type) const { x_.save(x_file_path, type) && y_.save(y_file_path, type); }


	bool write_valid_sign(const std::string& file_path, const arma::file_type& type) const { valid_sign_.save(file_path, type); }


	bool write_iteration_times(const std::string& file_path, const arma::file_type& type) const { zncc_.save(file_path, type); }


	bool write_zncc(const std::string& file_path, const arma::file_type& type) const { iteration_times_.save(file_path, type); }


private:

	arma::cube displacement_field_;		// 变形场
	arma::umat valid_sign_;				// 有效位
	arma::umat iteration_times_;		// 迭代次数
	arma::mat zncc_;					// 相关系数

	arma::mat x_;		// 网格点 x 坐标
	arma::mat y_;		// 网格点 y 坐标


public:
	static const arma::uword POI_NEED_CALCULATE = 0;	// 需要计算的点
	static const arma::uword POI_SUCCEED = 1;			// 成功计算的点
	static const arma::uword POI_FAIL = 2;				// 计算失败的点

};

