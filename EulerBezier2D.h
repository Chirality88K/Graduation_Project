#ifndef EULERBEZIER2D_H
#define EULERBEZIER2D_H
#include "thirdparty/opennurbs/opennurbs.h"
#include <vector>

namespace EulerBezier2D
{
	std::vector<double> ComputeLength(const ON_BezierCurve *OBC);
	std::vector<double> ComputeAngle(const ON_BezierCurve *OBC);
	bool EulerBezierSpiralCheck(const ON_BezierCurve *OBC);
	bool EulerBezierWeakCheck(const ON_BezierCurve* OBC);
	void SmoothingBezierControlPolygon(ON_BezierCurve *OBC);
	void Elevate(ON_BezierCurve *OBC);
	void EulerBezierSpiralInterpolation(ON_BezierCurve *OBC, int max_vtx_num = 50);
	void EulerBezierWeakInterpolation(ON_BezierCurve* OBC, int max_vtx_num = 50);
	void SmoothingCorner(ON_BezierCurve *OBC, ON_3dPoint Ps, ON_3dPoint O, double alpha);
	ON_NurbsCurve GenerateSmoothingCurve(ON_3dPoint start, ON_3dPoint corner, ON_3dPoint end);
	ON_NurbsCurve SmoothingCornerWithSlope(ON_3dPoint start, ON_3dPoint corner, ON_3dPoint end, double alpha = -1.0);
	void GenerateSymmetry(ON_BezierCurve *result, const ON_BezierCurve *OBC, ON_3dPoint O, ON_3dVector v);
	void EulerBezier2dTest(ONX_Model *model);
	void YangMethodtest(ONX_Model *model);
	void Pentagram(ONX_Model *model);

	double Compute_delta_theta(ON_2dPoint ps, ON_2dPoint pe, ON_2dVector vs, ON_2dVector ve, int cv_cnt);
	ON_BezierCurve ComputeEulerBezier2D_Directly(ON_2dPoint ps, ON_2dPoint pe, ON_2dVector vs, ON_2dVector ve, int cv_cnt, double& error);
	double Compute_L_for_fixed_cv_cnt(ON_2dPoint ps, ON_2dPoint pe, ON_2dVector vs, ON_2dVector ve, int cv_cnt, double* angles = nullptr);
}
#endif