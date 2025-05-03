#pragma once
#include "opennurbs.h"
#include <vector>

extern const double PI;
namespace EulerBezier2D {
	std::vector<double> ComputeLength(const ON_BezierCurve* OBC);
	std::vector<double> ComputeAngle(const ON_BezierCurve* OBC);
	bool EulerBezierSpiralCheck(const ON_BezierCurve* OBC);
	void SmoothingBezierControlPolygon(ON_BezierCurve* OBC);
	void Elevate(ON_BezierCurve* OBC);
	void EulerBezierSpiralInterpolation(ON_BezierCurve* OBC, int max_vtx_num);
	void SmoothingCorner(ON_BezierCurve* OBC, ON_3dPoint Ps, ON_3dPoint O, double alpha);
	void GenerateSymmetry(ON_BezierCurve* result, const ON_BezierCurve* OBC, ON_3dPoint O, ON_3dVector v);
	void EulerBezier2dTest(ONX_Model* model);
	void YangMethodtest(ONX_Model* model);
}