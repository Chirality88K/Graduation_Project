#ifndef CHIRALITYMATHTOOLS_H
#define CHIRALITYMATHTOOLS_H
#include "thirdparty/opennurbs/opennurbs.h"
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
};
#endif