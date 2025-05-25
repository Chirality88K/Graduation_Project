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
	// double start0, end0, start1, end1;
	// mRailCurve[0].GetDomain(&start0, &end0);
	mRailCurve[0].SetDomain(0.0, 1.0);
	// mRailCurve[1].GetDomain(&start1, &end1);
	mRailCurve[1].SetDomain(0.0, 1.0);
	int num_cv = 0;
	for (int i = 0; i <= num_of_bone; ++i)
	{
		double u = 1.0 / double(num_of_bone) * double(i);
		ON_3dVector vs = mVectorField[0](double(i) / double(num_of_bone));
		ON_3dVector ve = mVectorField[1](double(i) / double(num_of_bone));
		EulerBspline3D onc(mRailCurve[0].PointAt(u), mRailCurve[1].PointAt(u), vs, ve);
		num_cv = (std::max)(num_cv, onc.CVCount());
	}
	for (int i = 0; i <= num_of_bone; ++i)
	{
		double u = 1.0 / double(num_of_bone) * double(i);
		m_u_knots.push_back(u);
		ON_3dVector vs = mVectorField[0](double(i) / double(num_of_bone));
		ON_3dVector ve = mVectorField[1](double(i) / double(num_of_bone));
		ON_NurbsCurve *onc = new EulerBspline3D(mRailCurve[0].PointAt(u), mRailCurve[1].PointAt(u), vs, ve, num_cv);
		mBoneStructure.push_back(onc);
		ChiralityDebugInfo(*onc, "Bone Structure" + std::to_string(i));
	}
}

void Fillet_EB3D::GenerateFillet()
{
	int n = mBoneStructure[0]->CVCount();
	std::vector<std::pair<ON_3dVector, ON_3dVector>> pair_of_tan;
	ON_3dVector v0_s = mRailCurve[0].TangentAt(0.0);
	ON_3dVector v0_e = mRailCurve[0].TangentAt(1.0);
	ON_3dVector v1_s = mRailCurve[1].TangentAt(0.0);
	ON_3dVector v1_e = mRailCurve[1].TangentAt(1.0);
	double t0, t1;
	mBoneStructure[0]->GetDomain(&t0, &t1);
	for (int i = 0; i < n; ++i)
	{
		double u = 1.0 / double(n - 1) * double(i);
		pair_of_tan.push_back({v0_s * (1 - u) + v1_s * u, v0_e * (1 - u) + v1_e * u});
	}
	std::vector<ON_NurbsCurve> temp_curve_list;
	for (auto curve : mBoneStructure)
	{
		temp_curve_list.push_back(*curve);
	}
	ON_NurbsSurface ons = ChiralityMath::Skinning(temp_curve_list, m_u_knots, pair_of_tan);
	this->Create(3, false, ons.Order(0), ons.Order(1), ons.CVCount(0), ons.CVCount(1));
	ON_3dPoint p;
	for (int i = 0; i < ons.CVCount(0); ++i)
	{
		for (int j = 0; j < ons.CVCount(1); ++j)
		{
			ons.GetCV(i, j, p);
			this->SetCV(i, j, p);
		}
	}
	for (int i = 0; i < ons.KnotCount(0); ++i)
	{
		this->SetKnot(0, i, ons.Knot(0, i));
	}
	for (int j = 0; j < ons.KnotCount(1); ++j)
	{
		this->SetKnot(1, j, ons.Knot(1, j));
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
	test_fillet.GenerateFillet();
	const int layer_index1 = model->AddLayer(L"RailCurve", ON_Color::SaturatedMagenta);
	ChiralityAddNurbsCurve(model, rail[0], L"rail curve_0", layer_index1);
	ChiralityAddNurbsCurve(model, rail[1], L"rail_curve_1", layer_index1);
	const int layer_index2 = model->AddLayer(L"BoneStructure", ON_Color::SaturatedBlue);
	for (int i = 0; i < test_fillet.mBoneStructure.size(); ++i)
	{
		ChiralityAddNurbsCurve(model, *(test_fillet.mBoneStructure[i]), L"bone curve" + std::to_wstring(i), layer_index2);
	}
	const int layer_index3 = model->AddLayer(L"Fillet Surface", ON_Color::SaturatedGold);
	ChiralityAddNurbsSurface(model, test_fillet, L"Surface", layer_index3);
}
