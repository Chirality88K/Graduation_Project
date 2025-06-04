#ifndef EULERPOLYGON3D_H
#define EULERPOLYGON3D_H
#include "thirdparty/opennurbs/opennurbs.h"
class PolarPoint3d
{
public:
	PolarPoint3d() : mDistance(0.0), mTheta(0.0), mPhi(0.0) {}
	PolarPoint3d(double theta, double phi, double dis);
	PolarPoint3d(const ON_3dPoint &p);
	ON_3dPoint CartesianCoordinates() const;
	double mDistance; // 到原点的距离
	double mTheta;	  // 水平角[0,2\pi)
	double mPhi;	  // 俯仰角[-\pi/2,\pi/2]
};

class EulerPolygon3D
{
public:
	// 默认起点为原点，起始切向是(1,0,0)，PE位于xOy平面上
	EulerPolygon3D(PolarPoint3d PE, ON_3dVector ve);
	EulerPolygon3D(PolarPoint3d PS, PolarPoint3d PE, ON_3dVector vs,
				   ON_3dVector ve, int num = 50);
	static void EulerPolygonTest(ONX_Model *model);

private:
	std::vector<ON_3dPoint> mDiscretePolygon;
	// 计算相邻边的水平夹角，size为顶点数n，第一个和最后一个值都是0
	std::vector<double> ComputeDeltaTheta() const;
	// 计算相邻边的竖直夹角，size为顶点数n，第一个和最后一个值都是0
	std::vector<double> ComputeDeltaPhi() const;
	std::vector<double> ComputeLength() const;
	// 升阶
	void Elevate();
	void Smoothing();
	ON_NurbsCurve ToBezier() const;
};
#endif