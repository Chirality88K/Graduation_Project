#ifndef EULERBSPLINE3D_H
#define EULERBSPLINE3D_H
#include "thirdparty/opennurbs/opennurbs.h"
class EulerBspline3D : public ON_NurbsCurve
{
public:
	EulerBspline3D();																					// default 啥也不干
	EulerBspline3D(ON_3dVector vs, ON_3dVector ve, double length, int min_cv_count = 2);				// 起点是原点，终点是(length,0,0)，起始切向vs，终点切向是ve
	EulerBspline3D(ON_3dPoint ps, ON_3dPoint pe, ON_3dVector vs, ON_3dVector ve, int min_cv_count = 2); // 起点是ps，终点是pe，起始切向vs，终点切向是ve
	static void EulerBspline3DTest_MidPlaneMethod(ONX_Model *model);
};
#endif