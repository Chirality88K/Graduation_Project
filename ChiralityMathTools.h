#ifndef CHIRALITYMATHTOOLS_H
#define CHIRALITYMATHTOOLS_H
#include "thirdparty/opennurbs/opennurbs.h"
#include <vector>
namespace ChiralityMath
{
	double ArcLength(const ON_BezierCurve &, double from, double to);
	double ArcLength(const ON_NurbsCurve &, double from, double to);
	double ArcLength(const ON_BezierCurve &);
	double ArcLength(const ON_NurbsCurve &);
	double Bernstein(int n, int i, double t);
	double Torsion(const ON_BezierCurve &, double t);
	ON_NurbsCurve UniformG1(ON_3dPoint ps, ON_3dPoint pe, ON_3dVector vs, ON_3dVector ve);
	void Elevate(ON_NurbsCurve &onc);
	ON_NurbsCurve CubicBsplineInterpolate_G1(const std::vector<ON_3dPoint> &Q, const std::vector<double> &knot, ON_3dVector v0, ON_3dVector vn);
	ON_NurbsSurface Skinning(const std::vector<ON_NurbsCurve> &curve_list, const std::vector<double> &u_knots,
									 const std::vector<std::pair<ON_3dVector, ON_3dVector>> &pair_tangent);
};
#endif