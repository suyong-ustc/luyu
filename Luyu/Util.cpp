#include "Util.h"
using namespace arma;


double NormalizeVectoerize(const mat& m, vec& v)
{
	// vectorise
	v = vectorise(m);

	// zero mean
	v = v - mean(v);

	// normalization
	const double vv = norm(v, 2);
	v = v / vv;

	return vv;
}


double ZNCC(const mat& a, const mat& b)
{
	// vectorize
	vec aa = vectorise(a);
	vec bb = vectorise(b);

	// zero mean
	aa = aa - mean(aa);
	bb = bb - mean(bb);

	// normalization
	aa = normalise(aa);
	bb = normalise(bb);

	// correlation coefficient
	const double zncc = dot(aa, bb);

	return zncc;
}