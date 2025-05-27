#ifndef FILLET_USING_EB3D_H
#define FILLET_USING_EB3D_H
#include <functional>
#include "thirdparty/opennurbs/opennurbs.h"
class Fillet_EB3D : public ON_NurbsSurface
{
public:
	Fillet_EB3D();
	virtual ~Fillet_EB3D();
	void SetRailCurve(const ON_NurbsCurve &, const ON_NurbsCurve &);
	void SetVectorFeild(const std::function<ON_3dVector(double)> &, const std::function<ON_3dVector(double)> &);
	ON_3dVector GetTangent(bool zero_or_one, double t);
	void GenerateBone();
	void GenerateFillet();
	static void Fillet_EB3D_Test(ONX_Model *model);
	static void TwoSurfaces_Fillet_Test(ONX_Model *model);

private:
	ON_NurbsCurve mRailCurve[2];
	std::function<ON_3dVector(double)> mVectorField[2];
	std::vector<ON_NurbsCurve *> mBoneStructure;
	// 目前的u_knots是对于railcurve的定义进行等分得到的
	std::vector<double> m_u_knots;
};
#endif