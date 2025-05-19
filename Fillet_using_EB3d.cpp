#include "Fillet_using_EB3d.h"

Fillet_EB3D::Fillet_EB3D()
{
}

Fillet_EB3D::Fillet_EB3D(const ON_NurbsSurface &nurbs1, const ON_NurbsSurface &nurbs2)
{
	mSurface[0] = nurbs1;
	mSurface[1] = nurbs2;
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

void Fillet_EB3D::Generate()
{
}
