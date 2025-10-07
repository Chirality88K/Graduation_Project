#ifndef DISCRETESPIRAL_H
#define DISCRETESPIRAL_H
#include "Frenet.h"
#include "ChiralityMathTools.h"
#include "write3dm.h"
#include <vector>

class DiscreteSpiral
{
public:
	DiscreteSpiral(const std::vector<ON_3dPoint> &vp, FrenetFrame b1, FrenetFrame b2) : mPoints(vp)
	{
		mBoundary_Frame[0] = b1;
		mBoundary_Frame[1] = b2;
	}
	static void DiscreteSpiralTest(ONX_Model *model);
	static void DiscreteBodyTest(ONX_Model* model);
	static void SimpleTest(ONX_Model*, const FrenetFrame& f1, const FrenetFrame& f2, const std::string& name);
	static ON_NurbsCurve Interpolate(const FrenetFrame& f1, const FrenetFrame& f2);
	static void TestExit(ONX_Model* model);

private:
	std::vector<ON_3dPoint> mPoints;
	std::vector<FrenetFrame> mFrenetFrame;
	FrenetFrame mBoundary_Frame[2];
	void ComputeFrenetFrame();
	void SmoothFrenetFrameWithAngle();
	void SmoothFrenetFrameWithCur();
	bool ComputeAngles(std::vector<double>&)const;
	std::vector<double> ComputeDiscreteCurvature() const;
	std::vector<double> ComputeDiscreteTorsion() const;
	void Elevate();
	bool EndIteration();
	ON_BezierCurve GetBezier()const;
};

void SolveBoundaryTest(ONX_Model* model);
#endif