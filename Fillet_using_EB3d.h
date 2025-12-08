#ifndef FILLET_USING_EB3D_H
#define FILLET_USING_EB3D_H
#include <functional>
#include "thirdparty/opennurbs/opennurbs.h"
#include "DiscreteSpiral.h"

class Fillet_EB3D : public ON_NurbsSurface
{
public:
	Fillet_EB3D();
	virtual ~Fillet_EB3D();
	void SetRailCurve(const ON_NurbsCurve &, const ON_NurbsCurve &);
	void SetRailCurve(const std::function<FrenetFrame(double)>&, const std::function<FrenetFrame(double)>&);
	//void SetVectorFeild(const std::function<ON_3dVector(double)> &, const std::function<ON_3dVector(double)> &);
	void SetFrenetField(const std::function<ON_3dVector(double)>& f1, const std::function<ON_3dVector(double)>& f2);
	ON_3dVector GetTangent(bool zero_or_one, double t);
	void GenerateBone(bool is_set_bone_num = false, int bone_num = 10);
	void GenerateFillet();
	inline const std::vector<ON_NurbsCurve*>& GetBoneCurves()const
	{
		return mBoneStructure;
	}
	static void Fillet_EB3D_Test(ONX_Model *model);
	static void TwoSurfaces_Fillet_Test(ONX_Model *model);
	static void CircleSpiral_Test(ONX_Model* model);
	static void ThreeAngle_Test(ONX_Model* model);

private:
	ON_NurbsCurve mRailCurve[2];
	std::function<FrenetFrame(double)> mRail_Param_Curve[2];
	//std::function<ON_3dVector(double)> mVectorField[2];
	std::function<ON_3dVector(double)> mFrenetField[2];
	std::vector<ON_NurbsCurve *> mBoneStructure;
	// 目前的u_knots是对于railcurve的定义进行等分得到的
	std::vector<double> m_u_knots;
};
#endif