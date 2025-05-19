#ifndef CHIRALITYMATHTOOLS_H
#define CHIRALITYMATHTOOLS_H
#include "thirdparty/opennurbs/opennurbs.h"
namespace ChiralityMath
{
    double ArcLength(const ON_BezierCurve &, double from, double to);
    double ArcLength(const ON_NurbsCurve &, double from, double to);
    double Bernstein(int n, int i, double t);
    double Torsion(const ON_BezierCurve &, double t);
};
#endif