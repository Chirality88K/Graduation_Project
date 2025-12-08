#ifndef CHIRALITYMATHTOOLS_H
#define CHIRALITYMATHTOOLS_H
#include "thirdparty/opennurbs/opennurbs.h"
#include <vector>
#include <functional>
#include "Frenet.h"
#include "ParameterSurface.h"
#define PI ChiralityMath::Pi

namespace ChiralityMath
{
	const double Pi = acos(-1.0);
	// 二分法求函数零点，要求初始左右点函数值异号
	double Bisection(const std::function<double(double)>&, double L, double R, double* real_error = nullptr, double eps = 1e-8);
	//Newton法求零点
	double Newton(const std::function<double(double)>& f, const std::function<double(double)>& df, double x0, double eps = 1e-8, int max_iter_num = 10000);
	double ArcLength(const ON_BezierCurve &, double from, double to);
	double ArcLength(const ON_NurbsCurve &, double from, double to);
	double ArcLength(const ON_BezierCurve &);
	double ArcLength(const ON_NurbsCurve &);
	double ComputeSignedAngle(const ON_3dVector &v1, const ON_3dVector &v2, const ON_3dVector &T = ON_3dVector::ZAxis);
	// 返回一个参数序列，将给定的曲线的弧长等分
	std::vector<double> GenerateUniformArcLength(const ON_NurbsCurve &onc, int num_param);
	double Bernstein(int n, int i, double t);
	double Torsion(const ON_NurbsCurve &, double t);
	double DiscreteCurvature(ON_3dPoint p_before, ON_3dPoint p_mid, ON_3dPoint p_after);
	ON_NurbsCurve UniformG1(ON_3dPoint ps, ON_3dPoint pe, ON_3dVector vs, ON_3dVector ve);
	ON_BezierCurve BezierG1_xOy(ON_3dPoint ps, ON_3dPoint pe, ON_3dVector vs, ON_3dVector ve);
	void Elevate(ON_NurbsCurve &onc);
	void Elevate(ON_BezierCurve &obc);
	void Elevate(std::vector<ON_3dPoint>& vp);
	ON_NurbsCurve CubicBsplineInterpolate_G1(const std::vector<ON_3dPoint> &Q, const std::vector<double> &knot, ON_3dVector v0, ON_3dVector vn);
	ON_NurbsCurve CubicBsplineInterpolate_Period(const std::vector<ON_3dPoint> &Q, const std::vector<double> &knot);
	ON_NurbsCurve CubicHomoBsplineInterpolate_G1(const std::vector<ON_4dPoint>& Q, const std::vector<double>& knot, ON_4dPoint v0, ON_4dPoint vn);
	ON_NurbsCurve CubicHomoBsplineInterpolate_Period(const std::vector<ON_4dPoint>& Q, const std::vector<double>& knot);
	ON_NurbsSurface Skinning(const std::vector<ON_NurbsCurve> &curve_list, const std::vector<double> &u_knots, const std::vector<std::pair<ON_3dVector, ON_3dVector>> &pair_tangent);
	// 生成柱面，输入xy平面上的母线，沿着方向dir生成柱面，dir方向上的范围是t0到t1
	ON_NurbsSurface GenerateCylinder(const ON_NurbsCurve &parent_curve, ON_3dVector dir, double t0, double t1);
	FrenetFrame GetFrenet(const ON_NurbsCurve &onc, double t);
	FrenetFrame GetFrenet(const ON_NurbsSurface &ons, double u, double v);
	ON_NurbsCurve GetNurbsCircle(const ON_3dPoint &center, const ON_3dPoint &start_point, ON_3dVector normal);
	ON_4dPoint GetHomoPoint(const ON_NurbsCurve& onc, double t);
	ON_4dPoint GetHomoDerivative(const ON_NurbsCurve& onc, double t);

	ON_NurbsSurface GenerateRotating(const ON_NurbsCurve &parent_curve, const ON_Line &axis);
	ON_NurbsCurve ChangeDimensionFrom2To3(const ON_NurbsCurve &onc_2d);
	// 生成随机点，距离p在min_distance和max_distance之间
	ON_3dPoint GetRandomPoint(const ON_3dPoint &p, double min_distance, double max_distance);
	double GetRandomFloat(double min, double max);
	void KnotInsertion(ON_NurbsCurve& onc, double knot_insert);
	std::vector<ON_BezierCurve> SplitIntoBezier(const ON_NurbsCurve& onc);
	using VectorField = std::function<ON_3dVector(double)>;
	std::vector<ON_NurbsSurface> WirePlasticCurvedSurface(const std::vector<ParameterCurve>& v_pc, const std::vector<VectorField>& v_tan, int sample_cnt);
	std::vector<ON_Circle> TangentToCircle(const ON_Circle& circle1, const ON_Circle& circle2, const ON_Circle& circle3, bool in1, bool in2, bool in3);


	// 求解二阶微分方程组的边值问题，X0和XN为起点终点的Frenet标架，Length为弧长s，
	// kappa_tau_param为四个double参数，a,b,c,d，
	// 分别为k(s) = as+b;\tau(s) = cs+d,
	// n为采样数
	std::vector<ON_3dPoint> SolveBoundary(FrenetFrame X0, FrenetFrame XN, double step_length, double *kappa_tau_param, int n);
};
#endif