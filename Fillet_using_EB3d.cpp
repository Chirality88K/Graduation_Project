#include "Fillet_using_EB3d.h"
#include "ChiralityMathTools.h"
#include "EulerBspline3D.h"
#include "EulerBezier2D.h"
#include "write3dm.h"
#include "FixAxisBezier3D.h"
#include "EulerBspline2D.h"

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
	mRail_Param_Curve[0] = [rail1](double t)->FrenetFrame {
		return ChiralityMath::GetFrenet(rail1, t);
	};
	mRail_Param_Curve[1] = [rail2](double t)->FrenetFrame {
		return ChiralityMath::GetFrenet(rail2, t);
	};
}

void Fillet_EB3D::SetRailCurve(const std::function<FrenetFrame(double)>& pc1, const std::function<FrenetFrame(double)>& pc2)
{
	mRail_Param_Curve[0] = pc1;
	mRail_Param_Curve[1] = pc2;
}

void Fillet_EB3D::SetFrenetField(const std::function<ON_3dVector(double)>& f1, const std::function<ON_3dVector(double)>& f2)
{
	mFrenetField[0] = f1;
	mFrenetField[1] = f2;
}

ON_3dVector Fillet_EB3D::GetTangent(bool zero_or_one, double t)
{

	return ON_3dVector();
}

void Fillet_EB3D::GenerateBone(bool is_set_bone_num, int bone_num)
{
	int num_of_bone = bone_num;
	if (!is_set_bone_num)
	{
		num_of_bone = (std::max)(mRailCurve[0].CVCount(), mRailCurve[1].CVCount());
		mRailCurve[0].SetDomain(0.0, 1.0);
		mRailCurve[1].SetDomain(0.0, 1.0);
	}
	int num_cv = 0;
	mBoneStructure.clear();
	m_u_knots.clear();
	std::vector<ON_BezierCurve> temp_bezier;
	for (int i = 0; i <= num_of_bone; ++i)
	{
		double u = 1.0 / double(num_of_bone) * double(i);
		temp_bezier.push_back(FixAxisBezier3D::Interpolate(mRail_Param_Curve[0](u).GetPos(), mRail_Param_Curve[1](u).GetPos(),
			mFrenetField[0](u), mFrenetField[1](u)));
		num_cv = (std::max)(num_cv, temp_bezier.back().CVCount());
		m_u_knots.push_back(u);
	}
	for (ON_BezierCurve& obc : temp_bezier)
	{
		while (obc.CVCount() < num_cv)
		{
			ChiralityMath::Elevate(obc);
		}
		mBoneStructure.push_back(new ON_NurbsCurve(obc));
	}
}

void Fillet_EB3D::GenerateFillet()
{
	int n = mBoneStructure[0]->CVCount();
	std::vector<std::pair<ON_3dVector, ON_3dVector>> pair_of_tan;
	ON_3dVector v0_s = mRail_Param_Curve[0](0.0).GetAlpha();
	ON_3dVector v0_e = mRail_Param_Curve[0](1.0).GetAlpha();
	ON_3dVector v1_s = mRail_Param_Curve[1](0.0).GetAlpha();
	ON_3dVector v1_e = mRail_Param_Curve[1](1.0).GetAlpha();
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
	rail[0].SetDomain(0, 1);
	rail[1] = ChiralityMath::UniformG1(ON_3dPoint(0, 0, 20), ON_3dPoint(10, 0, 20), ON_3dVector(1, 2, 1), ON_3dVector(3, -1, -1));
	rail[1].SetDomain(0, 1);
	auto lambda1 = [rail](double t) -> ON_3dVector
	{
		ON_3dVector v1(0, -1, 1);
		ON_3dVector v2(1, 1, 1);
		ON_3dVector T = v1 * (1 - t) + v2 * t;
		return T;
	};
	auto lambda2 = [rail](double t) -> ON_3dVector
	{
		ON_3dVector v1(0, -1, 2);
		ON_3dVector v2(1, 1, 1);
		ON_3dVector T = v1 * (1 - t) + v2 * t;
		return T;
	};
	test_fillet.SetRailCurve(rail[0], rail[1]);
	test_fillet.SetFrenetField(lambda1, lambda2);
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

void Fillet_EB3D::TwoSurfaces_Fillet_Test(ONX_Model *model)
{
	// Compute 10 points
	double Acute_vertices_distance_to_origin = 10.0;
	double Blunt_vertices_distance_to_origin = Acute_vertices_distance_to_origin * sin(PI / 10.0) / sin(3 * PI / 10.0);
	double Acute_vertices_polar_angle[5] = {PI / 10.0, PI / 2.0, 162 * PI / 180.0, 234 * PI / 180.0, 306 * PI / 180.0};
	double Blunt_vertices_polar_angle[5] = {54 * PI / 180.0, 126 * PI / 180.0, 198 * PI / 180.0, 270 * PI / 180.0, 342 * PI / 180.0};
	ON_3dPoint Acute_vertices[5];
	ON_3dPoint Blunt_vertices[5];
	for (int i = 0; i < 5; ++i)
	{
		Acute_vertices[i] = ON_3dPoint(cos(Acute_vertices_polar_angle[i]), sin(Acute_vertices_polar_angle[i]), 0) * Acute_vertices_distance_to_origin;
	}
	for (int i = 0; i < 5; ++i)
	{
		Blunt_vertices[i] = ON_3dPoint(cos(Blunt_vertices_polar_angle[i]), sin(Blunt_vertices_polar_angle[i]), 0) * Blunt_vertices_distance_to_origin;
	}
	ON_3dPointArray Parray;
	for (int i = 0; i < 5; ++i)
	{
		Parray.Append(Acute_vertices[i]);
		Parray.Append(Blunt_vertices[i]);
	}
	Parray.Append(Acute_vertices[0]);
	// Compute Smoothing corner curves
	ON_NurbsCurve onc;
	for (int i = 0; i < 10; ++i)
	{
		ON_3dPoint Start = (i == 0) ? Parray[9] : Parray[i - 1];
		ON_3dPoint End = (i == 9) ? Parray[0] : Parray[i + 1];
		ON_3dPoint Corner = Parray[i];
		onc.Append(EulerBezier2D::GenerateSmoothingCurve(Start * 0.5 + Corner * 0.5, Corner, End * 0.5 + Corner * 0.5));
	}
	onc.SetDomain(0, 1);
	// Generate surface
	const int layer_index_cy = model->AddLayer(L"cylinder", ON_Color::SaturatedRed);
	const int layer_index_pl = model->AddLayer(L"plane", ON_Color::SaturatedCyan);
	ChiralityAddNurbsSurface(model, ChiralityMath::GenerateCylinder(onc, ON_3dVector(0, 0, 1), -5, 5), L"cylinder_surface", layer_index_cy);
	ON_PlaneSurface plane(ON_Plane({0, 0, 1, 0}));
	plane.Translate(ON_3dVector(-0.5, -0.5, 0));
	plane.Scale(24.0);
	plane.Translate(ON_3dVector(0, 0, -2));
	ChiralityAddPlane(model, plane, L"plane_surface", layer_index_pl);
	// Generate fillet
	Fillet_EB3D fillet;
	ON_NurbsCurve rail1 = ChiralityMath::ChangeDimensionFrom2To3(onc);
	ON_NurbsCurve rail2 = rail1;
	rail1.Translate(ON_3dVector(0, 0, 1));
	rail2.Scale(1.3);
	rail2.Translate(ON_3dVector(0, 0, -2));
	fillet.SetRailCurve(rail1, rail2);
	auto lambda1 = [&rail1](double t) -> ON_3dVector
	{
		return ON_3dVector(0, 0, -1);
	};
	auto lambda2 = [&rail2](double t) -> ON_3dVector
	{
		FrenetFrame f = ChiralityMath::GetFrenet(rail2, t);
		ON_3dPoint p = f.GetPos();
		ON_3dVector T = f.GetBeta();
		if (p.x * T.x + p.y * T.y <= 0)
		{
			T = -T;
		}
		return T;
	};
	fillet.SetFrenetField(lambda1, lambda2);
	fillet.GenerateBone();
	//fillet.GenerateFillet();
	// Add all about fillet
	const int layer_index1 = model->AddLayer(L"RailCurve", ON_Color::SaturatedMagenta);
	ChiralityAddNurbsCurve(model, rail1, L"rail curve_0", layer_index1);
	ChiralityAddNurbsCurve(model, rail2, L"rail_curve_1", layer_index1);
	const int layer_index2 = model->AddLayer(L"BoneStructure", ON_Color::SaturatedBlue);
	for (int i = 0; i < fillet.mBoneStructure.size(); ++i)
	{
		ChiralityAddNurbsCurve(model, *(fillet.mBoneStructure[i]), L"bone curve" + std::to_wstring(i), layer_index2);
	}
	//const int layer_index3 = model->AddLayer(L"Fillet Surface", ON_Color::SaturatedGold);
	//ChiralityAddNurbsSurface(model, fillet, L"Surface", layer_index3);
}

static void ElevateToSameOrder(ON_BezierCurve* obc1, ON_BezierCurve* obc2)
{
	const int n = abs(obc1->Order() - obc2->Order());
	if (obc1->Order() < obc2->Order())
	{
		for (int i = 0; i < n; ++i)
		{
			ChiralityMath::Elevate(*obc1);
		}
		return;
	}
	if (obc1->Order() > obc2->Order())
	{
		for (int i = 0; i < n; ++i)
		{
			ChiralityMath::Elevate(*obc2);
		}
		return;
	}
}

void Fillet_EB3D::CircleSpiral_Test(ONX_Model* model)
{
	double R = 5.0;
	double r = 0.6;
	int circle = 15;
	auto Curve = [R, r, circle](double theta, double offset = 0)->std::pair<ON_3dPoint, ON_3dVector> {
		double phi = circle * theta + offset;
		ON_3dPoint p = ON_3dPoint(cos(theta), sin(theta), 0) * R +
			r * ON_3dPoint(-cos(phi) * cos(theta), -cos(phi) * sin(theta), sin(phi));
		ON_3dVector v = R * ON_3dVector(-sin(theta), cos(theta), 0) +
			r * circle * ON_3dVector(sin(phi) * cos(theta) + cos(phi) * sin(theta), sin(phi) * sin(theta), cos(phi));
		return std::make_pair(p, v);
	};

	const int sample_cnt = circle;
	for (int kk = 0; kk < sample_cnt; ++kk)
	{
		double theta1 = double(kk) * PI / double(sample_cnt);
		double theta3 = double(kk + 1) * PI / double(sample_cnt);
		double theta2 = (theta1 + theta3) / 2;
		ON_3dPoint p1 = Curve(theta1).first;
		ON_3dPoint p2 = Curve(theta2).first;
		ON_3dPoint p3 = Curve(theta3).first;
		ON_3dVector v1 = Curve(theta1).second;
		ON_3dVector v2 = Curve(theta2).second;
		ON_3dVector v3 = Curve(theta3).second;
		ON_BezierCurve obc1_red = FixAxisBezier3D::Interpolate(p1, p2, v1, v2);
		ON_BezierCurve obc2_red = FixAxisBezier3D::Interpolate(p2, p3, v2, v3);

		p1 = Curve(theta1, PI).first;
		p2 = Curve(theta2, PI).first;
		p3 = Curve(theta3, PI).first;
		v1 = Curve(theta1, PI).second;
		v2 = Curve(theta2, PI).second;
		v3 = Curve(theta3, PI).second;
		ON_BezierCurve obc1_blue = FixAxisBezier3D::Interpolate(p1, p2, v1, v2);
		ON_BezierCurve obc2_blue = FixAxisBezier3D::Interpolate(p2, p3, v2, v3);

		ElevateToSameOrder(&obc1_red, &obc1_blue);
		ElevateToSameOrder(&obc2_red, &obc2_blue);

		ON_NurbsCurve blue = obc1_blue;
		blue.Append(obc2_blue);
		ON_NurbsCurve red = obc1_red;
		red.Append(obc2_red);

		ON_NurbsSurface ons(3, false, red.Order(), 2, red.CVCount(), 2);
		for (int i = 0; i < red.CVCount(); ++i)
		{
			ON_3dPoint p;
			red.GetCV(i, p);
			ons.SetCV(i, 0, p);
			blue.GetCV(i, p);
			ons.SetCV(i, 1, p);
		}
		for (int i = 0; i < red.KnotCount(); ++i)
		{
			ons.SetKnot(0, i, *(red.Knot() + i));
		}
		ons.SetKnot(1, 0, 0);
		ons.SetKnot(1, 1, 1);

		const int layer_index = model->AddLayer(L"test_surface", ON_Color::SaturatedGold);
		ChiralityAddNurbsSurface(model, ons, std::to_wstring(kk), layer_index);
		//for (int i = 0; i < circle; ++i)
		{
			//ons.Rotate(2 * PI / double(circle), ON_3dVector::ZAxis, ON_3dPoint::Origin);
			//ChiralityAddNurbsSurface(model, ons, L"test_surface", layer_index);
		}
	}
	
	
	
}

void Fillet_EB3D::ThreeAngle_Test(ONX_Model* model)
{
	double R = 5.0;
	std::vector <ON_3dPoint> vp;
	ON_3dPoint begin_p(R * cos(PI / 2 - PI / 6), R * sin(PI / 2 - PI / 6), 0.0);
	for (int i = 0; i < 3; ++i)
	{
		double mid_theta = PI / 2 + PI * 2 / 3 * i;
		ON_3dVector v(cos(mid_theta), sin(mid_theta), 0.0);
		ON_3dPoint P1(R * cos(mid_theta - PI / 6), R * sin(mid_theta - PI / 6), 0.0);
		ON_3dPoint P2(R * cos(mid_theta + PI / 6), R * sin(mid_theta + PI / 6), 0.0);
		ON_3dPoint Q1 = P1 + 1.5 * R * v;
		ON_3dPoint Q2 = P2 + 1.5 * R * v;
		vp.push_back(P1);
		vp.push_back(Q1);
		vp.push_back(Q2);
		vp.push_back(P2);
	}
	vp.push_back(begin_p);
	const int layer_index = model->AddLayer(L"test", ON_Color::SaturatedBlue);
	ChiralityAddLines(model, vp, L"Polygons", layer_index);
	double radius = 2.0;
	std::vector<ON_NurbsCurve> v_onc;
	for (int i = 0; i < vp.size() - 1; ++i)
	{
		int pre = i - 1;
		int next = i + 1;
		if (i == 0)
		{
			pre = vp.size() - 2;
		}
		v_onc.push_back(EulerBspline2D::GenerateSmoothingCorner(vp[i] + (vp[pre] - vp[i]) / (vp[pre] - vp[i]).Length() * radius,
			vp[i], vp[i] + (vp[next] - vp[i]) / (vp[next] - vp[i]).Length() * radius));
	}
	const int index = model->AddLayer(L"Curve", ON_Color::SaturatedMagenta);
	for (int i = 1; i < v_onc.size(); ++i)
	{
		ON_3dPoint Start = v_onc[0].PointAtEnd();
		ON_3dPoint End = v_onc[i].PointAtStart();
		ON_NurbsCurve line(2, false, 2, 2);
		line.SetKnot(0, 0); line.SetKnot(1, 1);
		line.SetCV(0, Start); line.SetCV(1, End);
		v_onc[0].Append(line);
		v_onc[0].Append(v_onc[i]);	
	}
	ON_3dPoint Start = v_onc[0].PointAtEnd();
	ON_3dPoint End = v_onc[0].PointAtStart();
	ON_NurbsCurve line(2, false, 2, 2);
	line.SetKnot(0, 0); line.SetKnot(1, 1);
	line.SetCV(0, Start); line.SetCV(1, End);
	v_onc[0].Append(line);
	ChiralityAddNurbsCurve(model, v_onc[0], L"Curve", index);
	ON_NurbsSurface ons = ChiralityMath::GenerateCylinder(v_onc[0], ON_3dVector::ZAxis, 2, 10);
	const int sur_index = model->AddLayer(L"Surface", ON_Color::SaturatedGold);
	ChiralityAddNurbsSurface(model, ons, L"Surface", sur_index);
	ON_Plane plane(ON_3dPoint::Origin, ON_3dVector::ZAxis);
	ON_PlaneSurface ops(plane);
	ops.Translate(ON_3dVector(-0.5, -0.5, 0.0));
	ops.Scale(40.0);
	const int plane_index = model->AddLayer(L"plane", ON_Color::SaturatedCyan);
	ChiralityAddPlane(model, ops, L"plane", plane_index);
	ON_NurbsCurve rail1 = ChiralityMath::ChangeDimensionFrom2To3(v_onc[0]);
	rail1.Translate(ON_3dVector(0, 0, 2));
	ON_NurbsCurve rail2 = ChiralityMath::ChangeDimensionFrom2To3(v_onc[0]);
	rail2.Scale(1.3);
	Fillet_EB3D fillet;
	rail1.SetDomain(0.0, 1.0);
	rail2.SetDomain(0.0, 1.0);
	fillet.SetRailCurve(rail1, rail2);
	auto lambda1 = [&rail1](double t) -> ON_3dVector
	{
		return ON_3dVector(0, 0, -1);
	};
	auto lambda2 = [&rail2](double t) -> ON_3dVector
	{
		ON_3dPoint p = rail2.PointAt(t);
		ON_3dVector T = ON_3dVector(p);
		return T;
	};
	fillet.SetFrenetField(lambda1, lambda2);
	fillet.GenerateBone();
	fillet.GenerateFillet();
	const int rail_index = model->AddLayer(L"RailCurve", ON_Color::SaturatedMagenta);
	ChiralityAddNurbsCurve(model, rail1, L"rail curve_0", rail_index);
	ChiralityAddNurbsCurve(model, rail2, L"rail_curve_1", rail_index);
	const int bone_index = model->AddLayer(L"BoneStructure", ON_Color::SaturatedBlue);
	for (int i = 0; i < fillet.mBoneStructure.size(); ++i)
	{
		ChiralityAddNurbsCurve(model, *(fillet.mBoneStructure[i]), L"bone curve" + std::to_wstring(i), bone_index);
	}
	const int fillet_index = model->AddLayer(L"Fillet Surface", ON_Color::SaturatedGold);
	ChiralityAddNurbsSurface(model, fillet, L"Surface", fillet_index);
	ChiralityDebugInfo(fillet);
}
