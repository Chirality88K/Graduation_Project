#ifndef EULERBEZIER3D_H
#define EULERBEZIER3D_H
#include <vector>
#include "thirdparty/opennurbs/opennurbs.h"

class EulerBezier3D : public ON_BezierCurve
{
public:
	EulerBezier3D();															 // default 啥也不干
	EulerBezier3D(ON_3dVector vs, ON_3dVector ve, double length);				 // 起点是原点，终点是(length,0,0)，起始切向vs，终点切向是ve
	EulerBezier3D(ON_3dPoint ps, ON_3dPoint pe, ON_3dVector vs, ON_3dVector ve); // 起点是ps，终点是pe，起始切向vs，终点切向是ve
	void Smoothing();
	static void EulerBezier3DTest(ONX_Model *model);
	static void EulerBezier3DTest_MidPlaneMethod(ONX_Model *model);
	static ON_3dPoint GetBallCenter(ON_3dPoint p0, ON_3dPoint p1, ON_3dVector t0, ON_3dVector t1);
	std::vector<double> ComputeAngle();		  // 长度为控制点个数n
	std::vector<double> ComputeLength();	  // 长度为控制点个数减一
	std::vector<ON_3dVector> ComputeNormal(); // 长度为控制点个数
	void Elevate();
	static void Rise_to_3D(ON_NurbsCurve *onc, ON_3dPoint pe, ON_3dVector ts, ON_3dVector te);
};

#endif // !EULERBEIZER3D_H
