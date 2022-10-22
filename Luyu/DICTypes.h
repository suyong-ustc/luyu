#pragma once


// Interpolation Algorithms
enum InterpType {
	INTERP_KEYS,
	INTERP_BSPLINE3,
	INTERP_BSPLINE5,
	INTERP_BSPLINE7
};


// Gradient Estimators
enum GradType {
	GRAD_PREWITT,
	GRAD_BARRON,
	GRAD_SOBEL,
	GRAD_BSPLINE3,
	GRAD_BSPLINE5,
	GRAD_BSPLINE7
};


// Shape Functions
enum ShapeFuncType {
	SHAPEFUNC_FIRST,
	SHAPEFUNC_SECOND
};

constexpr int FIRST_ORDER_SHAPE_FUNC_PARA_COUNT = 6;
constexpr int SECOND_ORDER_SHAPE_FUNC_PARA_COUNT = 12;


// Subpixel Registration Strategy
enum RegistrationType {
	REGISTRATION_FORWARD,
	REGISTRATION_INVERSE
};