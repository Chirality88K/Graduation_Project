#include "thirdparty/opennurbs/opennurbs.h"
#include "ChiralityMathTools.h"
#include "write3dm.h"
#include "FixAxisBezier3D.h"
#include "Fillet_using_EB3d.h"
#include "OBJGenerator.h"
#include <iostream>

static const ON_3dPoint A(4.0, 0.0, 0.0);
static const ON_3dPoint B(5.0, 0.0, 9.0);
static const ON_Line rotate_axis(ON_3dPoint::Origin, ON_3dPoint(0, 0, 10));
static const double M_param = 0.8;
static const double N_param = 0.1;
static const ON_3dPoint M = B * M_param + A * (1 - M_param);
static const ON_3dPoint N = B * N_param + A * (1 - N_param);
static ON_BezierCurve handle_line;
static const double handle_radius = 0.4;
static const double up_handle_param = 0.07;
static const double down_handle_param = 0.9;
static const double up_cup_rail_radius = 0.9;
static const double down_cup_rail_radius = 1.0;
static ON_NurbsCurve up_cup_rail_nurbs;
static ON_NurbsCurve down_cup_rail_nurbs;
static ON_NurbsCurve up_handle_rail_nurbs;
static ON_NurbsCurve down_handle_rail_nurbs;
static const ON_3dPoint low_A = A - (B - A) / (B - A).Length() * 0.5;
static OBJGenerator Handle_obj;

void Cup_Body(ONX_Model* model)
{
	ON_BezierCurve line_cup(3, false, 2);
	line_cup.SetCV(0, low_A); line_cup.SetCV(1, B);
	ON_NurbsSurface Cup = ChiralityMath::GenerateRotating(line_cup, rotate_axis);
	const int cup_surface_layer_index = model->AddLayer(L"Cup_Surface_layer", ON_Color::SaturatedBlue);
	ChiralityAddNurbsSurface(model, Cup, L"Cup", cup_surface_layer_index);
	double r1 = sqrt(A.x * A.x + A.y * A.y);
	double r2 = sqrt(B.x * B.x + B.y * B.y);
	double h = abs(B.z - A.z);
	auto cup_position = [r1, r2, h](double u, double v)->ON_3dPoint {
		return (r1 * (1 - u) + r2 * u) * ON_3dPoint(cos(v), sin(v), 0.0) + ON_3dVector(0, 0, u * h);
	};

	auto cup_normal = [r1, r2, h](double u, double v)->ON_3dVector {
		ON_3dVector der_u = (r2 - r1) * ON_3dVector(cos(v), sin(v), 0.0) + ON_3dVector(0, 0, h);
		ON_3dVector der_v = (r1 * (1 - u) + r2 * u) * ON_3dVector(-sin(v), cos(v), 0.0);
		return ON_3dVector::CrossProduct(der_u, der_v);
	};
	double range[4] = { -(low_A - A).Length() / (A - B).Length(),1.0,0.0,2 * PI };
	ParameterSurface param_cup_surface(cup_position, cup_normal, range);
	Handle_obj.AddParameterSurface(param_cup_surface, "Cup_Body", 50, 50, false, true);
	
	auto up_cup_rail = [r1, r2, h](double t)->FrenetFrame
	{
		t = t * 2 * PI;
		double rM = r1 * (1 - M_param) + r2 * M_param;
		ON_2dPoint uv_p = up_cup_rail_radius * ON_2dPoint(cos(t), sin(t));
		double u = uv_p.x / A.DistanceTo(B) + M_param;
		double v = -uv_p.y / rM;
		ON_3dPoint m = (r1 * (1 - u) + r2 * u) * ON_3dPoint(cos(v), sin(v), 0.0) + ON_3dVector(0, 0, u * h);

		ON_3dVector der_u = (r2 - r1) * ON_3dVector(cos(v), sin(v), 0.0) + ON_3dVector(0, 0, h);
		ON_3dVector der_v = (r1 * (1 - u) + r2 * u) * ON_3dVector(-sin(v), cos(v), 0.0);
		double udt = -up_cup_rail_radius * sin(t) / A.DistanceTo(B);
		double vdt = -up_cup_rail_radius * cos(t) / rM;
		ON_3dVector der = der_u * udt + der_v * vdt;
		ON_3dVector n = ON_3dVector::CrossProduct(der_u, der_v);
		ON_3dVector beta = ON_3dVector::CrossProduct(n, der);
		return FrenetFrame(m, der, beta);
	};

	constexpr int N = 20;
	std::vector<ON_3dPoint> points;
	std::vector<double> arc_param;
	for (int i = 0; i <= N; ++i)
	{
		points.push_back(up_cup_rail(double(i % N) / double(N)).GetPos());
		arc_param.push_back(double(i) / double(N));
	}
	up_cup_rail_nurbs = ChiralityMath::CubicBsplineInterpolate_Period(points, arc_param);
	const int cup_rail_layer_index = model->AddLayer(L"Cup_Rail_layer", ON_Color::SaturatedRed);
	ChiralityAddNurbsCurve(model, up_cup_rail_nurbs, L"up_cup_rail", cup_rail_layer_index);

	auto down_cup_rail = [r1, r2, h](double t)->FrenetFrame
	{
		double n_param = N_param;
		t = t * 2 * PI;
		double rM = r1 * (1 - n_param) + r2 * n_param;
		ON_2dPoint uv_p = down_cup_rail_radius * ON_2dPoint(cos(t), sin(t));
		double u = -uv_p.x / A.DistanceTo(B) + n_param;
		double v = -uv_p.y / rM;
		ON_3dPoint m = (r1 * (1 - u) + r2 * u) * ON_3dPoint(cos(v), sin(v), 0.0) + ON_3dVector(0, 0, u * h);

		ON_3dVector der_u = (r2 - r1) * ON_3dVector(cos(v), sin(v), 0.0) + ON_3dVector(0, 0, h);
		ON_3dVector der_v = (r1 * (1 - u) + r2 * u) * ON_3dVector(-sin(v), cos(v), 0.0);
		double udt = -up_cup_rail_radius * sin(t) / A.DistanceTo(B);
		double vdt = -up_cup_rail_radius * cos(t) / rM;
		ON_3dVector der = der_u * udt + der_v * vdt;
		ON_3dVector n = ON_3dVector::CrossProduct(der_u, der_v);
		ON_3dVector beta = ON_3dVector::CrossProduct(n, der);
		return FrenetFrame(m, der, beta);
	};
	points.clear();
	for (int i = 0; i <= N; ++i)
	{
		points.push_back(down_cup_rail(double(i % N) / double(N)).GetPos());
	}
	down_cup_rail_nurbs = ChiralityMath::CubicBsplineInterpolate_Period(points, arc_param);
	ChiralityAddNurbsCurve(model, down_cup_rail_nurbs, L"down_cup_rail", cup_rail_layer_index);
}

void Handle_Body(ONX_Model* model)
{
	ON_3dVector V_M = ON_3dVector(1, 0, 1);
	ON_3dVector V_N = ON_3dVector(-1, 0, -sqrt(3));
	handle_line = FixAxisBezier3D::Interpolate(M, N, V_M, V_N);
	const int handle_layer_index = model->AddLayer(L"handle_layer", ON_Color::SaturatedRed);
	ChiralityAddNurbsCurve(model, handle_line, L"handle_line", handle_layer_index);

	auto pos = [](double u, double v)->ON_3dPoint {
		FrenetFrame f = ChiralityMath::GetFrenet(handle_line, u);
		return f.GetPos() + (f.GetBeta() * cos(v) + f.GetGamma() * sin(v)) * handle_radius;
	};

	auto der_u = [](double u, double v)->ON_3dVector {
		FrenetFrame f = ChiralityMath::GetFrenet(handle_line, u);
		double kappa = handle_line.CurvatureAt(u).Length();
		double tau = ChiralityMath::Torsion(handle_line, u);
		ON_3dVector der_u = handle_line.DerivativeAt(u) + handle_radius * (cos(v) * (-kappa * f.GetAlpha() + tau * f.GetGamma()) + sin(v) * (-tau * f.GetBeta()));
		return der_u;
	};

	auto normal = [](double u, double v)->ON_3dVector {
		FrenetFrame f = ChiralityMath::GetFrenet(handle_line, u);
		double kappa = handle_line.CurvatureAt(u).Length();
		double tau = ChiralityMath::Torsion(handle_line, u);
		ON_3dVector der_u = handle_line.DerivativeAt(u) + handle_radius * (cos(v) * (-kappa * f.GetAlpha() + tau * f.GetGamma()) + sin(v) * (-tau * f.GetBeta()));
		ON_3dVector der_v = -sin(v) * f.GetBeta() + cos(v) * f.GetGamma();
		ON_3dVector n = ON_3dVector::CrossProduct(der_u, der_v);
		n.Unitize();
		return n;
	};

	constexpr int N = 20;//number of bone curves
	std::vector<ON_BezierCurve> Outline(N);
	int max_cv_num = 0.0;
	{
		std::vector<ON_3dPoint> Up_handle_rail_points;
		std::vector<ON_3dPoint> Down_handle_rail_points;
		std::vector<double> arc_param(1, 0.0);
		ON_3dPoint Up_handle_center = handle_line.PointAt(up_handle_param);
		ON_3dPoint Down_handle_center = handle_line.PointAt(down_handle_param);
		ON_3dVector Up_circle_normal = handle_line.TangentAt(up_handle_param);
		ON_3dVector Down_circle_normal = handle_line.TangentAt(down_handle_param);
		ON_3dVector Up_circle_radius = Up_circle_normal; Up_circle_radius.Rotate(-1, 0, ON_3dVector::YAxis);
		ON_3dVector Down_circle_radius = Down_circle_normal; Down_circle_radius.Rotate(-1, 0, ON_3dVector::YAxis);
		for (int i = 0; i < N; ++i)
		{
			ON_3dPoint ps = Up_handle_center + handle_radius * Up_circle_radius;
			Up_handle_rail_points.push_back(ps);
			ON_3dPoint pe = Down_handle_center + handle_radius * Down_circle_radius;
			Down_handle_rail_points.push_back(pe);
			Up_circle_radius.Rotate(2.0 * PI / double(N), Up_circle_normal);
			Down_circle_radius.Rotate(2.0 * PI / double(N), Down_circle_normal);
			Outline[i] = FixAxisBezier3D::Interpolate(ps, pe, Up_circle_normal, Down_circle_normal);
			arc_param.push_back(double(i + 1) / double(N));
			max_cv_num = (std::max)(max_cv_num, Outline[i].CVCount());
			ChiralityAddNurbsCurve(model, Outline[i], L"test_outline" + std::to_wstring(i), handle_layer_index);
		}
		Up_handle_rail_points.push_back(Up_handle_rail_points.front());
		Down_handle_rail_points.push_back(Down_handle_rail_points.front());
		up_handle_rail_nurbs = ChiralityMath::CubicBsplineInterpolate_Period(Up_handle_rail_points, arc_param);
		down_handle_rail_nurbs = ChiralityMath::CubicBsplineInterpolate_Period(Down_handle_rail_points, arc_param);
		ChiralityAddNurbsCurve(model, up_handle_rail_nurbs, L"Up_handle_rail", handle_layer_index);
		ChiralityAddNurbsCurve(model, down_handle_rail_nurbs, L"Down_handle_rail", handle_layer_index);
	}
	for (ON_BezierCurve& obc : Outline)
	{
		while (obc.CVCount() < max_cv_num)
		{
			ChiralityMath::Elevate(obc);
		}
	}
	std::vector<ON_NurbsCurve> skining_result;
	std::vector<double> arc_knot(N + 1, 0.0);
	for (int i = 1; i <= N; ++i)
	{
		arc_knot[i] = double(i) / double(N);
	}
	for (int i = 0; i < max_cv_num; ++i)
	{
		std::vector<ON_3dPoint> Q(N + 1);
		for (int j = 0; j < N; ++j)
		{
			Outline[j].GetCV(i, Q[j]);
		}
		Q[N] = Q[0];
		skining_result.push_back(ChiralityMath::CubicBsplineInterpolate_Period(Q, arc_knot));
	}
	ON_NurbsSurface handle_surface(3, false, Outline[0].Order(), 4, max_cv_num, skining_result[0].CVCount());
	for (int i = 0; i < max_cv_num; ++i)
	{
		for (int j = 0; j < skining_result[0].CVCount();++j)
		{
			ON_3dPoint p;
			skining_result[i].GetCV(j, p);
			handle_surface.SetCV(i, j, p);
		}
	}
	for (int i = 0; i < handle_surface.KnotCount(0); ++i)
	{
		handle_surface.SetKnot(0, i, i < handle_surface.KnotCount(0) / 2 ? 0 : 1);
	}
	for (int j = 0; j < handle_surface.KnotCount(1); ++j)
	{
		handle_surface.SetKnot(1, j, skining_result[0].Knot(j));
	}

	const int handle_surface_layer_index = model->AddLayer(L"handle_surface_layer", ON_Color::SaturatedGold);
	ChiralityAddNurbsSurface(model, handle_surface, L"Handle_surface", handle_surface_layer_index);
	Handle_obj.AddNurbsSurface(handle_surface, "handle_body", 200, 50, false, true);
}

void Up_Fillet(ONX_Model* model)
{
	ON_3dVector vv = B - A;
	vv.Unitize(); vv.Rotate(1, 0, ON_3dVector::YAxis);
	ON_3dPoint temp_center = M + vv * 0.2;
	auto up_cup_rail_tan = [temp_center](double t)->ON_3dVector
	{
		t = t * 2 * PI;
		double r1 = sqrt(A.x * A.x + A.y * A.y);
		double r2 = sqrt(B.x * B.x + B.y * B.y);
		double h = abs(B.z - A.z);
		double rM = r1 * (1 - M_param) + r2 * M_param;
		ON_2dPoint uv_p = up_cup_rail_radius * ON_2dPoint(cos(t), sin(t));
		double u = uv_p.x / A.DistanceTo(B) + M_param;
		double v = -uv_p.y / rM;

		//ON_3dVector der_u = (r2 - r1) * ON_3dVector(cos(v), sin(v), 0.0) + ON_3dVector(0, 0, h);
		//ON_3dVector der_v = (r1 * (1 - u) + r2 * u) * ON_3dVector(-sin(v), cos(v), 0.0);
		///double udt = -up_cup_rail_radius * sin(t) / A.DistanceTo(B);
		//double vdt = -up_cup_rail_radius * cos(t) / rM;
		//ON_3dVector der = der_u * vdt - der_v * udt;
		ON_3dPoint m = (r1 * (1 - u) + r2 * u) * ON_3dPoint(cos(v), sin(v), 0.0) + ON_3dVector(0, 0, u * h);
		ON_3dVector der = temp_center - m;
		return der;
	};
	constexpr int N = 20;
	int max_cv_num = 0;
	std::vector<ON_BezierCurve> Outline(N);
	const int fillet_layer_index = model->AddLayer(L"Up_Fillet_Layer", ON_Color::SaturatedGreen);
	for (int i = 0; i < N; ++i)
	{
		double u = double(i) / double(N);
		Outline[i] = FixAxisBezier3D::Interpolate(up_cup_rail_nurbs.PointAt(u), up_handle_rail_nurbs.PointAt(u),
			up_cup_rail_tan(u), handle_line.TangentAt(up_handle_param));
		max_cv_num = (std::max)(max_cv_num, Outline[i].CVCount());
		//ChiralityAddNurbsCurve(model, Outline[i], L"Curves" + std::to_wstring(i), fillet_curves_layer_index);
	}
	for (ON_BezierCurve& obc : Outline)
	{
		while (obc.CVCount() < max_cv_num)
		{
			ChiralityMath::Elevate(obc);
		}
	}
	std::vector<ON_NurbsCurve> skinning_result;
	std::vector<double> s_param;
	ON_3dPoint p;
	for (int i = 0; i <= N; ++i)
	{
		s_param.push_back(double(i) / double(N));
	}
	for (int i = 0; i < max_cv_num; ++i)
	{
		std::vector<ON_3dPoint> Q;
		for (int j = 0; j < N; ++j)
		{
			Outline[j].GetCV(i, p);
			Q.push_back(p);
		}
		Q.push_back(Q[0]);
		skinning_result.push_back(ChiralityMath::CubicBsplineInterpolate_Period(Q, s_param));
	}
	ON_NurbsSurface fillet(3, false, Outline[0].Order(), 4, Outline[0].CVCount(), skinning_result[0].CVCount());
	for (int i = 0; i < Outline[0].CVCount(); ++i)
	{
		for (int j = 0; j < skinning_result[0].CVCount(); ++j)
		{
			skinning_result[i].GetCV(j, p);
			fillet.SetCV(i, j, p);
		}
	}
	for (int i = 0; i < fillet.KnotCount(0); ++i)
	{
		fillet.SetKnot(0, i, i < fillet.KnotCount(0) / 2 ? 0 : 1);
	}
	for (int j = 0; j < fillet.KnotCount(1); ++j)
	{
		fillet.SetKnot(1, j, skinning_result[0].Knot(j));
	}
	ChiralityAddNurbsSurface(model, fillet, L"up_fillet_surface", fillet_layer_index);
	Handle_obj.AddNurbsSurface(fillet, "up_fillet", 50, 50, false, true);
}

void Down_Fillet(ONX_Model* model)
{
	ON_3dVector vv = B - A;
	vv.Unitize(); vv.Rotate(1, 0, ON_3dVector::YAxis);
	ON_3dPoint temp_center = N + vv * 0.2;
	auto down_cup_rail_tan = [temp_center](double t)->ON_3dVector
	{
		t = t * 2 * PI;
		double n_param = N_param;
		double r1 = sqrt(A.x * A.x + A.y * A.y);
		double r2 = sqrt(B.x * B.x + B.y * B.y);
		double h = abs(B.z - A.z);
		double rN = r1 * (1 - n_param) + r2 * n_param;
		ON_2dPoint uv_p = down_cup_rail_radius * ON_2dPoint(cos(t), sin(t));
		double u = -uv_p.x / A.DistanceTo(B) + n_param;
		double v = -uv_p.y / rN;

		//ON_3dVector der_u = (r2 - r1) * ON_3dVector(cos(v), sin(v), 0.0) + ON_3dVector(0, 0, h);
		//ON_3dVector der_v = (r1 * (1 - u) + r2 * u) * ON_3dVector(-sin(v), cos(v), 0.0);
		///double udt = -up_cup_rail_radius * sin(t) / A.DistanceTo(B);
		//double vdt = -up_cup_rail_radius * cos(t) / rM;
		//ON_3dVector der = der_u * vdt - der_v * udt;
		ON_3dPoint m = (r1 * (1 - u) + r2 * u) * ON_3dPoint(cos(v), sin(v), 0.0) + ON_3dVector(0, 0, u * h);
		ON_3dVector der = temp_center - m;
		return der;
	};
	constexpr int N = 20;
	int max_cv_num = 0;
	std::vector<ON_BezierCurve> Outline(N);
	const int fillet_layer_index = model->AddLayer(L"Down_Fillet_Layer", ON_Color::SaturatedGreen);
	for (int i = 0; i < N; ++i)
	{
		double u = double(i) / double(N);
		Outline[i] = FixAxisBezier3D::Interpolate(down_cup_rail_nurbs.PointAt(u), down_handle_rail_nurbs.PointAt(u),
			down_cup_rail_tan(u), -handle_line.TangentAt(down_handle_param));
		max_cv_num = (std::max)(max_cv_num, Outline[i].CVCount());
		//ChiralityAddNurbsCurve(model, Outline[i], L"Curves" + std::to_wstring(i), fillet_layer_index);
	}
	for (ON_BezierCurve& obc : Outline)
	{
		while (obc.CVCount() < max_cv_num)
		{
			ChiralityMath::Elevate(obc);
		}
	}
	std::vector<ON_NurbsCurve> skinning_result;
	std::vector<double> s_param;
	ON_3dPoint p;
	for (int i = 0; i <= N; ++i)
	{
		s_param.push_back(double(i) / double(N));
	}
	for (int i = 0; i < max_cv_num; ++i)
	{
		std::vector<ON_3dPoint> Q;
		for (int j = 0; j < N; ++j)
		{
			Outline[j].GetCV(i, p);
			Q.push_back(p);
		}
		Q.push_back(Q[0]);
		skinning_result.push_back(ChiralityMath::CubicBsplineInterpolate_Period(Q, s_param));
	}
	ON_NurbsSurface fillet(3, false, Outline[0].Order(), 4, Outline[0].CVCount(), skinning_result[0].CVCount());
	for (int i = 0; i < Outline[0].CVCount(); ++i)
	{
		for (int j = 0; j < skinning_result[0].CVCount(); ++j)
		{
			skinning_result[i].GetCV(j, p);
			fillet.SetCV(i, j, p);
		}
	}
	for (int i = 0; i < fillet.KnotCount(0); ++i)
	{
		fillet.SetKnot(0, i, i < fillet.KnotCount(0) / 2 ? 0 : 1);
	}
	for (int j = 0; j < fillet.KnotCount(1); ++j)
	{
		fillet.SetKnot(1, j, skinning_result[0].Knot(j));
	}
	ChiralityAddNurbsSurface(model, fillet, L"down_fillet_surface", fillet_layer_index);
	Handle_obj.AddNurbsSurface(fillet, "down_fillet", 50, 50, false, true);
}

int main(){
    std::cout<<"Cup\n";
	SetOutput("Cup");
	const std::string filename = "Cup.3dm";
	ON::Begin();
	ONX_Model model_to_write;
	Internal_SetExampleModelProperties(model_to_write, OPENNURBS__FUNCTION__, filename.c_str());
	model_to_write.AddDefaultLayer(nullptr, ON_Color::UnsetColor);
	Cup_Body(&model_to_write);
	Handle_Body(&model_to_write);
	Up_Fillet(&model_to_write);
	Down_Fillet(&model_to_write);
	ChiralityWrite3dmModel(&model_to_write, filename);
	Handle_obj.Write("Cup_Handle");
	return 0;
}