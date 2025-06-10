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

const int ITERATIONTIMES = 20;

class EulerPolygon3D
{
public:
	enum CurveType
	{
		Bezier,
		B_spline
	};

public:
	EulerPolygon3D(ON_3dPoint PS, ON_3dPoint PE, ON_3dVector vs, ON_3dVector ve, CurveType ct = Bezier);
	ON_NurbsCurve GetCurve() const { return mCurve; }
	static void EulerPolygonTest_ForConicSpiral(ONX_Model *model);
	static void EulerPolygonTest_ForSphereSpiral(ONX_Model *model);
	static void EulerPolygonTest_ForCircularHelix(ONX_Model *model);

private:
	// 默认起点为原点，起始切向是(1,0,0)，PE位于xOy平面上
	void BuildUpToBezier(ON_3dPoint PE, ON_3dVector ve);
	void BuildUpToB_Spline(ON_3dPoint PE, ON_3dVector ve);
	// 计算相邻边的水平夹角，size为顶点数n，第一个和最后一个值都是0
	std::vector<double> ComputeDeltaTheta() const;
	// 计算相邻边的竖直夹角，size为顶点数n，第一个和最后一个值都是0
	std::vector<double> ComputeDeltaPhi() const;
	std::vector<double> ComputeLength() const;
	std::vector<double> ComputeProjectionLength() const;
	// 升阶
	void Elevate();
	void SmoothingToBezier();
	void SmoothingToB_Spline(ON_3dPoint PE, ON_3dVector ve);

private:
	std::vector<ON_3dPoint> mDiscretePolygon;
	ON_NurbsCurve mCurve;
};
#endif