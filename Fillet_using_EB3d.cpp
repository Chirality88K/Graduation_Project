#include "Fillet_using_EB3d.h"
#include "ChiralityMathTools.h"
#include "EulerBspline3D.h"
#include "write3dm.h"

Fillet_EB3D::Fillet_EB3D()
{
}

Fillet_EB3D::~Fillet_EB3D()
{
	for (ON_NurbsCurve *&onc : mBoneStructure)
	{
		delete onc;
		onc = nullptr;
	}
}

void Fillet_EB3D::SetRailCurve(const ON_NurbsCurve &rail1, const ON_NurbsCurve &rail2)
{
	mRailCurve[0] = rail1;
	mRailCurve[1] = rail2;
}

void Fillet_EB3D::SetVectorFeild(const std::function<ON_3dVector(double)> &v1, const std::function<ON_3dVector(double)> &v2)
{
	mVectorField[0] = v1;
	mVectorField[1] = v2;
}

ON_3dVector Fillet_EB3D::GetTangent(bool zero_or_one, double t)
{
	std::function<ON_3dVector(double)> VF = mVectorField[zero_or_one ? 1 : 0];
	double t0, t1;
	mRailCurve[zero_or_one ? 1 : 0].GetDomain(&t0, &t1);
	ON_3dVector v = VF((t - t0) / (t1 - t0));
	v.Unitize();
	return v;
}

void Fillet_EB3D::GenerateBone()
{
	const int num_of_bone = (std::max)(mRailCurve[0].CVCount(), mRailCurve[1].CVCount()) + 1;
	// double total_arclength0 = ChiralityMath::ArcLength(mRailCurve[0]);
	// double total_arclength1 = ChiralityMath::ArcLength(mRailCurve[1]);
	double start0, end0, start1, end1;
	mRailCurve[0].GetDomain(&start0, &end0);
	mRailCurve[1].GetDomain(&start1, &end1);
	int num_cv = 0;
	for (int i = 0; i <= num_of_bone; ++i)
	{
		double u0 = start0 + (end0 - start0) / double(num_of_bone) * double(i);
		double u1 = start1 + (end1 - start1) / double(num_of_bone) * double(i);
		ON_3dVector vs = mVectorField[0](double(i) / double(num_of_bone));
		ON_3dVector ve = mVectorField[1](double(i) / double(num_of_bone));
		EulerBspline3D onc(mRailCurve[0].PointAt(u0), mRailCurve[1].PointAt(u1), vs, ve);
		num_cv = (std::max)(num_cv, onc.CVCount());
	}
	for (int i = 0; i <= num_of_bone; ++i)
	{
		double u0 = start0 + (end0 - start0) / double(num_of_bone) * double(i);
		double u1 = start1 + (end1 - start1) / double(num_of_bone) * double(i);
		ON_3dVector vs = mVectorField[0](double(i) / double(num_of_bone));
		ON_3dVector ve = mVectorField[1](double(i) / double(num_of_bone));
		ON_NurbsCurve *onc = new EulerBspline3D(mRailCurve[0].PointAt(u0), mRailCurve[1].PointAt(u1), vs, ve, num_cv);
		mBoneStructure.push_back(onc);
		ChiralityDebugInfo(*onc, "Bone Structure" + std::to_string(i));
	}
}

void Fillet_EB3D::Fillet_EB3D_Test(ONX_Model *model)
{
	Fillet_EB3D test_fillet;
	ON_NurbsCurve rail[2];
	rail[0] = ChiralityMath::UniformG1(ON_3dPoint::Origin, ON_3dPoint(10, 0, 0), ON_3dVector(1, 2, 1), ON_3dVector(3, -1, -1));
	rail[1] = ChiralityMath::UniformG1(ON_3dPoint(0, 0, 20), ON_3dPoint(10, 0, 20), ON_3dVector(1, 2, 1), ON_3dVector(3, -1, -1));
	auto lambda1 = [](double t) -> ON_3dVector
	{
		ON_3dVector v1(0, -1, 1);
		ON_3dVector v2(1, 1, 1);
		return v1 * (1 - t) + v2 * t;
	};
	auto lambda2 = [](double t) -> ON_3dVector
	{
		ON_3dVector v1(0, -1, 2);
		ON_3dVector v2(1, 1, 1);
		return v1 * (1 - t) + v2 * t;
	};
	test_fillet.SetRailCurve(rail[0], rail[1]);
	test_fillet.SetVectorFeild(lambda1, lambda2);
	test_fillet.GenerateBone();
}
