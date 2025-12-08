#include "LogAestheticBezier.h"
#include "EditingCurvatureCurve.h"
#include "Cornu_Spiral.h"
#include "ChiralityMathTools.h"
#include "FixAxisBezier3D.h"
#include "write3dm.h"
#include "EulerBezier2D.h"
#include "EulerBspline2D.h"

void FixBoundaryTest(ONX_Model*model)
{
	ON_3dPoint ps(0, 0, 0), pe(10, 0, 0);
	ON_3dVector vs(0, 2, 0), ve(-3, -1, 0);
	int N = 20;
	ON_Color color1(255, 0, 0);
	ON_Color color2(0, 255, 0);
	for (int i = 0; i <= N; ++i)
	{
		if (i == 10)
		{
			continue;
		}
		double alpha = -4.0 + double(i) / double(N) * 8.0;
		auto lambda = [alpha](double s)->double {
			return (s - alpha) * (s - alpha) + 1;
		};
		ON_BezierCurve obc = EditingCurvatureCurve::BezierInterpolate(ps, pe, vs, ve, lambda);
		ChiralityDebugforR(obc, "EditingCurvatureCurve_alpha_" + std::to_string(alpha));
		ON_Color color;
		double r, g, b;
		r = color1.FractionRed() * (1 - double(i) / double(N)) + color2.FractionRed() * double(i) / double(N);
		g = color1.FractionGreen() * (1 - double(i) / double(N)) + color2.FractionGreen() * double(i) / double(N);
		b = color1.FractionBlue() * (1 - double(i) / double(N)) + color2.FractionBlue() * double(i) / double(N);
		double max_rgb = (std::max)(r, (std::max)(g, b));
		color.SetFractionalRGB(r / max_rgb, g / max_rgb, b / max_rgb);
		const int layer_index = model->AddLayer(std::to_wstring(alpha).c_str(), color);
		ChiralityAddNurbsCurve(model, obc, L"curve", layer_index);
	}
}

void CornuSpiralTest(ONX_Model* model)
{
	ON_3dPoint ps(0, 0, 0);
	ON_3dPoint pe(12, 3, 0);
	ON_2dVector ve(3, 4);
	const int N = 15;
	ON_Color color1(255, 0, 255);
	ON_Color color2(0, 255, 0);
	for (int i = 0; i < N; ++i)
	{
		double theta = 2 * PI / N * i;
		ON_3dVector vs(1, 0, 0);
		vs.Rotate(theta, ON_3dVector::ZAxis);
		ON_2dVector vs_2d(vs.x, vs.y);
		ON_NurbsCurve onc = Cornu_Spiral::GetNurbs(ps, pe, vs_2d, ve);
		ON_Color color;
		double r, g, b;
		r = color1.FractionRed() * (1 - double(i) / double(N)) + color2.FractionRed() * double(i) / double(N);
		g = color1.FractionGreen() * (1 - double(i) / double(N)) + color2.FractionGreen() * double(i) / double(N);
		b = color1.FractionBlue() * (1 - double(i) / double(N)) + color2.FractionBlue() * double(i) / double(N);
		double max_rgb = (std::max)(r, (std::max)(g, b));
		color.SetFractionalRGB(r / max_rgb, g / max_rgb, b / max_rgb);
		const int layer_index = model->AddLayer((L"Spiral " + std::to_wstring(i)).c_str(), color);
		ChiralityAddNurbsCurve(model, onc, L"Spiral " + std::to_wstring(i), layer_index);

	}
}

void TestBezierDeltaTheta(ONX_Model * model)
{
	ON_2dPoint PS(0, 0);
	ON_2dPoint PE(10, 0);
	const int N = 10;
	for (int i = 0; i < 1; ++i)
	{
		ON_2dVector VS(1, 2);
		VS.Rotate(2 * PI / N * i);
		for (int j = 0; j < N; ++j)
		{
			ON_2dVector VE(2, 3);
			VE.Rotate(2 * PI / N * j);
			double error;
			ON_BezierCurve obc = EulerBezier2D::ComputeEulerBezier2D_Directly(PS, PE, VS, VE, 10, error);
			const int layer_index = model->AddLayer((L"Curve" + std::to_wstring(i) + std::to_wstring(j)).c_str(), ON_Color::SaturatedGold);
			ChiralityAddNurbsCurve(model, obc, L"Curve" + std::to_wstring(i) + std::to_wstring(j), layer_index);
			ChiralityDebugforR(obc, "Plane_Curve" + std::to_string(i) + std::to_string(j) + doubleToScientificString(error));
			ChiralityPrintBezierForPython(obc, "Plane_Curve" + std::to_string(i) + std::to_string(j));
		}
	}
}

void TestBsplineDeltaTheta(ONX_Model* model) 
{
	ON_2dPoint PS(0, 0);
	ON_2dPoint PE(10, 0);
	const int N = 10;
	for (int i = 4; i < 5; ++i)
	{
		ON_2dVector VS(1, 2);
		VS.Rotate(2 * PI / N * i);
		for (int j = 0; j < N; ++j)
		{
			ON_2dVector VE(2, 3);
			VE.Rotate(2 * PI / N * j);
			double error;
			ON_NurbsCurve onc = EulerBspline2D::ComputeEulerBspline2D_Directly(PS, PE, VS, VE, error);
			const int layer_index = model->AddLayer((L"Curve" + std::to_wstring(i) + std::to_wstring(j)).c_str(), ON_Color::SaturatedGold);
			ChiralityAddNurbsCurve(model, onc, L"Curve" + std::to_wstring(i) + std::to_wstring(j), layer_index);
			//ChiralityDebugforR(onc, "Plane_Curve" + std::to_string(i) + std::to_string(j) + doubleToScientificString(error));
			ChiralityPrintCubicIntKnotBSplineForPython(onc, "Plane_Curve" + std::to_string(i) + std::to_string(j));
		}
	}
}

int main()
{
	const std::string Name = "EulerBspline";
	std::cout << Name << "\n";
	SetOutput(Name);
	const std::string filename = Name + ".3dm";
	ON::Begin();
	ONX_Model model_to_write;
	Internal_SetExampleModelProperties(model_to_write, OPENNURBS__FUNCTION__, filename.c_str());
	model_to_write.AddDefaultLayer(nullptr, ON_Color::UnsetColor);
	//TestBezierDeltaTheta(&model_to_write);
	TestBsplineDeltaTheta(&model_to_write);

	ChiralityWrite3dmModel(&model_to_write, filename);
}