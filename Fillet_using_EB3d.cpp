#include "Fillet_using_EB3d.h"
#include "ChiralityMathTools.h"
#include "EulerBspline3D.h"

Fillet_EB3D::Fillet_EB3D()
{
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
	const int num_of_bone = (std::max)(mRailCurve[0].CVCount(),mRailCurve[1].CVCount())+1;
	//double total_arclength0 = ChiralityMath::ArcLength(mRailCurve[0]);
	//double total_arclength1 = ChiralityMath::ArcLength(mRailCurve[1]);
	double start0,end0,start1,end1;
	mRailCurve[0].GetDomain(&start0,&end0);
	mRailCurve[1].GetDomain(&start1,&end1);
	int num_cv=0;
	for(int i=0;i<=num_of_bone;++i){
		double u0 = start0+(end0-start0)/double(num_of_bone)*double(i);
		double u1 = start1+(end1-start1)/double(num_of_bone)*double(i);
		ON_3dVector vs = mVectorField[0](u0);
		ON_3dVector ve = mVectorField[1](u1);
		ON_NurbsCurve* onc = new EulerBspline3D(mRailCurve[0].PointAt(u0),mRailCurve[1].PointAt(u1),vs,ve);
		mBoneStructure.push_back(onc);
		num_cv=(std::max)(num_cv,onc->CVCount());
	}

	
}
