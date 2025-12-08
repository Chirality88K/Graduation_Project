#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include "ChiralityMathTools.h"
#include "write3dm.h"
#include "BoundaryRecorder.h"
#include "FixAxisBezier3D.h"
#include "EulerBezier2D.h"

ParameterCurve GenerateCircle(double x1,double y1,double z1, double x2, double y2, double z2, double x3, double y3, double z3)
{
	ON_Circle circle(ON_3dPoint(x1, y1, z1), ON_3dPoint(x2, y2, z2), ON_3dPoint(x3, y3, z3));
	auto lambda = [circle](double t)->FrenetFrame
	{
		ON_3dPoint p = circle.PointAt(t);
		ON_3dVector alpha = circle.TangentAt(t);
		ON_3dVector gamma = circle.Normal();
		ON_3dVector beta = ON_3dVector::CrossProduct(gamma, alpha);
		return FrenetFrame(p, alpha, beta);
	};
	double range[2] = { 0,2 * PI };
	return ParameterCurve(lambda, range);
}

void TestWirePlasticCurvedSurface(ONX_Model* model)
{
	ParameterCurve c1 = GenerateCircle(10, 0, 0, 0, 10, 0, -10, 0, 0);
	ParameterCurve c2 = GenerateCircle(5, 0, 5, 0, 5, 5, -5, 0, 5);
	ParameterCurve c3 = GenerateCircle(3 * sqrt(3), 0, 8, 0, 6, 10, -3 * sqrt(3), 0, 12);
	ParameterCurve c4 = GenerateCircle(9, 0, 14, 3, 6 * sqrt(2), 20, -3, 0, 26);
	ChiralityMath::VectorField v1 = [](double t)->ON_3dVector {
		return ON_3dVector::ZAxis;
	};
	ChiralityMath::VectorField v2 = [c2](double t)->ON_3dVector {
		ON_3dVector cp = c2.GetFrame(t).GetAlpha();
		ON_3dVector v = c2.GetFrame(t).GetGamma();
		v.Rotate(-PI / 3, cp);
		return v;
	};
	ChiralityMath::VectorField v3 = [c3](double t)->ON_3dVector {
		return c3.GetFrame(t).GetGamma();
	};
	ChiralityMath::VectorField v4 = [c4](double t)->ON_3dVector {
		return c4.GetFrame(t).GetGamma();
	};
	std::vector<ON_NurbsSurface> v_ons = ChiralityMath::WirePlasticCurvedSurface({ c1,c2,c3,c4 }, { v1,v2,v3,v4 }, 10);
	std::wstring name[2] = { L"test1",L"test2" };
	int i = 0;
	for (const ON_NurbsSurface& ons : v_ons)
	{
		const int layer_index = model->AddLayer(name[i].c_str(), ON_Color::SaturatedGreen);
		ChiralityAddNurbsSurface(model, ons, name[i].c_str(), layer_index);
		++i;
	}
}

void SpecialBoundaryTest(ONX_Model* model)
{
	BoundaryRecorder::GetRecorder().Read();
	std::vector<Boundary> bo = BoundaryRecorder::GetRecorder().GetBoundary();
	int i = 0;
	for (const Boundary& b : bo)
	{
		ON_BezierCurve obc = FixAxisBezier3D::Interpolate(b.ps_, b.pe_, b.vs_, b.ve_);
		const int layer_index = model->AddLayer((std::wstring(L"Curve") + std::to_wstring(i)).c_str(), ON_Color::SaturatedGold);
		ChiralityAddNurbsCurve(model, obc, L"Curve" + std::to_wstring(i), layer_index);
		const int true_curve_layer_index = model->AddLayer((std::wstring(L"True_Curve") + std::to_wstring(i)).c_str(), ON_Color::SaturatedMagenta);
		ON_NurbsCurve true_onc = ChiralityMath::CubicBsplineInterpolate_G1({ b.ps_,b.pe_ }, { 0,1 }, b.vs_, b.ve_);
		ChiralityAddNurbsCurve(model, true_onc, L"True_Curve" + std::to_wstring(i), true_curve_layer_index);
	}
}

void DrawGoblet(ONX_Model* model)
{
	std::vector<ON_3dPoint> P = {
			ON_3dPoint(0,6,0),ON_3dPoint(1,6,0),
			ON_3dPoint(1,4,0),ON_3dPoint(3,4,0),
			ON_3dPoint(3,2,0),ON_3dPoint(6,2,0),
			ON_3dPoint(6,0,0),ON_3dPoint(0,0,0) };
	ON_NurbsCurve ONC[6];
	ON_3dPoint Q[14];
	Q[0] = P[0];
	Q[13] = P[7];
	const int curve_layer_index = model->AddLayer(L"Curve", ON_Color::SaturatedGold);
	for (int i = 1; i < 7; ++i)
	{
		ON_3dPoint start = P[i] + (P[i - 1] - P[i]) * 0.3 / (P[i - 1] - P[i]).Length();
		ONC[i - 1] = EulerBezier2D::GenerateSmoothingCurve(start, P[i], P[i + 1]);
		//ChiralityAddNurbsCurve(model, ONC[i - 1], L"Curve" + std::to_wstring(i), curve_layer_index);
		Q[i * 2 - 1] = ONC[i - 1].PointAtStart();
		Q[i * 2] = ONC[i - 1].PointAtEnd();
	}
	const int line_layer_index = model->AddLayer(L"Line", ON_Color::Black);
	ChiralityAddLines(model, P, L"lines", line_layer_index);
	ON_NurbsCurve result;
	for (int i = 0; i <= 6; ++i)
	{
		ON_NurbsCurve onc;
		ON_LineCurve line(Q[i * 2], Q[i * 2 + 1]);
		line.GetNurbForm(onc);
		result.Append(onc);
		if (i < 6)
		{
			result.Append(ONC[i]);
		}
	}
	ChiralityAddNurbsCurve(model, result, L"Curve", curve_layer_index);
	ON_NurbsSurface ons = ChiralityMath::GenerateRotating(result, ON_Line(P.back(), P.front()));
	const int surface_layer_index = model->AddLayer(L"Surface_Layer", ON_Color::SaturatedGreen);
	ChiralityAddNurbsSurface(model, ons, L"Surface", surface_layer_index);
}

void TestAlphaSlope(ONX_Model* model)
{
	ON_3dPoint A(0, 10, 0), B(0, 0, 0), C(10, 0, 0);
	constexpr int N = 10;
	double ALPHA[N] = { -5,-4,-3,-2,-1.5,-1,-0.8,-0.6,-0.4,-0.2 };
	ON_Color colors[N];
	for (int i = 0; i < N; ++i)
	{
		colors[i] = ON_Color::SaturatedRed * (1 - double(i) / N)
			+ ON_Color::SaturatedGreen * double(i) / N;
		double max_rgb = (std::max)(colors[i].FractionRed(), colors[i].FractionBlue());
		max_rgb = (std::max)(colors[i].FractionGreen(), max_rgb);
		colors[i].SetFractionalRGB(colors[i].FractionRed() / max_rgb,
			colors[i].FractionGreen() / max_rgb, colors[i].FractionBlue() / max_rgb);
	}
	for (int i = 0; i < N; ++i)
	{
		ON_NurbsCurve onc = EulerBezier2D::SmoothingCornerWithSlope((A + B) / 2, B, (B + C) / 2, ALPHA[i]);
		const int layer_index = model->AddLayer((L"Curve" + std::to_wstring(ALPHA[i])).c_str(), colors[i]);
		ChiralityAddNurbsCurve(model, onc, L"Curve" + std::to_wstring(ALPHA[i]), layer_index);
		ChiralityDebugforR(onc, "alpha=" + std::to_string(ALPHA[i]));
		ChiralityPrintCubicIntKnotBSplineForPython(onc, "alpha=" + std::to_string(ALPHA[i]));
	}
}

void RotatingStar(ONX_Model* model)
{
	const int N = 5;
	const double R = 1.5;
	const double XX = 10.0 / (1 + sqrt(3));
	std::vector<ON_3dPoint> P = { ON_3dPoint(0,10,0),ON_3dPoint(XX,XX,0),ON_3dPoint(10,0,0),
	ON_3dPoint(XX,-XX,0) ,ON_3dPoint(0,-10,0) };
	ON_NurbsCurve ONC[N];
	ON_BezierCurve obc_s, obc_e;
	EulerBezier2D::SmoothingCorner(&obc_s, P[0] + (P[1] - P[0]) / (P[1] - P[0]).Length() * R, P[0], PI / 3 * 2);
	EulerBezier2D::GenerateSymmetry(&obc_e, &obc_s, ON_3dPoint::Origin, ON_3dVector::XAxis);
	obc_s.Reverse();
	ONC[0] = obc_s;
	ONC[N - 1] = obc_e;
	const int curve_layer_index = model->AddLayer(L"Curve", ON_Color::SaturatedGold);
	for (int i = 1; i < N - 1; ++i)
	{
		ON_3dPoint start = P[i] + (P[i - 1] - P[i]) / (P[i - 1] - P[i]).Length() * R;
		ONC[i] = EulerBezier2D::GenerateSmoothingCurve(start, P[i], P[i + 1]);
	}
	const int line_layer_index = model->AddLayer(L"Line", ON_Color::Black);
	ChiralityAddLines(model, P, L"lines", line_layer_index);
	
	ON_NurbsCurve result;
	for (int i = 0; i < N; ++i)
	{
		result.Append(ONC[i]);
		if (i < N - 1)
		{
			ON_NurbsCurve onc;
			ON_LineCurve line(ONC[i].PointAtEnd(), ONC[i + 1].PointAtStart());
			line.GetNurbForm(onc);
			result.Append(onc);
		}
	}
	ChiralityAddNurbsCurve(model, result, L"Curve", curve_layer_index);
	ON_NurbsSurface ons = ChiralityMath::GenerateRotating(result, ON_Line(P.back(), P.front()));
	ons.Rotate(1, 0, ON_3dVector::XAxis, ON_3dPoint::Origin);
	const int surface_layer_index = model->AddLayer(L"Surface_Layer", ON_Color::SaturatedGreen);
	ChiralityAddNurbsSurface(model, ons, L"Surface", surface_layer_index);
	ChiralityPrintCubicIntKnotBSplineForPython(result, "Rotating_star");
}

int main()
{
	std::cout << "Test\n";
	SetOutput("Test");
	const std::string filename = "Test.3dm";
	ON::Begin();
	ONX_Model model_to_write;
	Internal_SetExampleModelProperties(model_to_write, OPENNURBS__FUNCTION__, filename.c_str());
	model_to_write.AddDefaultLayer(nullptr, ON_Color::UnsetColor);
	//SpecialBoundaryTest(&model_to_write);
	//DrawGoblet(&model_to_write);
	//TestAlphaSlope(&model_to_write);
	RotatingStar(&model_to_write);
	ChiralityWrite3dmModel(&model_to_write, filename);
}