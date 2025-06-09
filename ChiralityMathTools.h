#ifndef CHIRALITYMATHTOOLS_H
#define CHIRALITYMATHTOOLS_H
#include "thirdparty/opennurbs/opennurbs.h"
#include <vector>
#include <functional>
namespace ChiralityMath
{
	// 二分法求函数零点，要求初始左右点函数值异号
	double Bisection(const std::function<double(double)> &, double L, double R, double eps = 1e-8);
	double ArcLength(const ON_BezierCurve &, double from, double to);
	double ArcLength(const ON_NurbsCurve &, double from, double to);
	double ArcLength(const ON_BezierCurve &);
	double ArcLength(const ON_NurbsCurve &);
	// 返回一个参数序列，将给定的曲线的弧长等分
	std::vector<double> GenerateUniformArcLength(const ON_NurbsCurve &onc, int num_param);
	double Bernstein(int n, int i, double t);
	double Torsion(const ON_NurbsCurve &, double t);
	ON_NurbsCurve UniformG1(ON_3dPoint ps, ON_3dPoint pe, ON_3dVector vs, ON_3dVector ve);
	void Elevate(ON_NurbsCurve &onc);
	ON_NurbsCurve CubicBsplineInterpolate_G1(const std::vector<ON_3dPoint> &Q, const std::vector<double> &knot, ON_3dVector v0, ON_3dVector vn);
	ON_NurbsSurface Skinning(const std::vector<ON_NurbsCurve> &curve_list, const std::vector<double> &u_knots, const std::vector<std::pair<ON_3dVector, ON_3dVector>> &pair_tangent);
	// 生成柱面，输入xy平面上的母线，沿着方向dir生成柱面，dir方向上的范围是t0到t1
	ON_NurbsSurface GenerateCylinder(const ON_NurbsCurve &parent_curve, ON_3dVector dir, double t0, double t1);
	ON_NurbsCurve ChangeDimensionFrom2To3(const ON_NurbsCurve &onc_2d);
};
#endif