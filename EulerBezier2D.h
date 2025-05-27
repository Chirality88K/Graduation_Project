#ifndef EULERBEZIER2D_H
#define EULERBEZIER2D_H
#include "thirdparty/opennurbs/opennurbs.h"
#include <vector>

extern const double PI;
namespace EulerBezier2D
{
	std::vector<double> ComputeLength(const ON_BezierCurve *OBC);
	std::vector<double> ComputeAngle(const ON_BezierCurve *OBC);
	bool EulerBezierSpiralCheck(const ON_BezierCurve *OBC);
	void SmoothingBezierControlPolygon(ON_BezierCurve *OBC);
	void Elevate(ON_BezierCurve *OBC);
	void EulerBezierSpiralInterpolation(ON_BezierCurve *OBC, int max_vtx_num = 50);
	void SmoothingCorner(ON_BezierCurve *OBC, ON_3dPoint Ps, ON_3dPoint O, double alpha);
	ON_NurbsCurve GenerateSmoothingCurve(ON_3dPoint start, ON_3dPoint corner, ON_3dPoint end);
	void GenerateSymmetry(ON_BezierCurve *result, const ON_BezierCurve *OBC, ON_3dPoint O, ON_3dVector v);
	void EulerBezier2dTest(ONX_Model *model);
	void YangMethodtest(ONX_Model *model);
	void Pentagram(ONX_Model *model);
}
#endif