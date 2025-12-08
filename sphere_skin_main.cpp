#include "write3dm.h"
#include "ChiralityMathTools.h"
#include "SphereSkinning.h"
#include "BoundaryRecorder.h"
#include "FixAxisBezier3D.h"

void CircleTangentTest(ONX_Model* model)
{
	ON_Circle c1( ON_3dPoint::Origin, 3.0);
	ON_Circle c2( ON_3dPoint(5, 0, 0), 2.0);
	ON_Circle c3(ON_3dPoint(7, 1, 0), sqrt(5) - 2.0);
	ON_NurbsCurve onc;
	const int org_circles_layer_index = model->AddLayer(L"Input Circles", ON_Color::Black);
	c1.GetNurbForm(onc);
	ChiralityAddNurbsCurve(model, onc, L"Input Circle1", org_circles_layer_index);
	c2.GetNurbForm(onc);
	ChiralityAddNurbsCurve(model, onc, L"Input Circle2", org_circles_layer_index);
	c3.GetNurbForm(onc);
	ChiralityAddNurbsCurve(model, onc, L"Input Circle3", org_circles_layer_index);
	for (unsigned int i = 0; i <= 0b111; ++i)
	{
		bool in1 = bool((i >> 2) & 1);
		bool in2 = bool((i >> 1) & 1);
		bool in3 = bool((i >> 0) & 1);
		ON_Color color;
		color.SetRGB(in1 ? 255 : 0, in2 ? 255 : 0, in3 ? 255 : 0);
		std::vector<ON_Circle> v_circle = ChiralityMath::TangentToCircle(c1, c2, c3, in1, in2, in3);
		if (v_circle.empty())
		{
			continue;
		}
		const int result_layer_index = model->AddLayer((L"Tangent Circles " + std::to_wstring(i)).c_str(), color);
		for (const ON_Circle& c : v_circle)
		{
			c.GetNurbForm(onc);
			ChiralityAddNurbsCurve(model, onc, L"Tangent Circle", result_layer_index);
		}
	}
}

void SphereSkinTest(ONX_Model* model)
{
	SphereSkinning ss;
	const int N = 7;
	ON_Color Red(255, 0, 0);
	ON_Color Green(0, 255, 0);
	ON_Color color[N];
	ON_Color resverse_color[N];
	ON_3dPoint Center_Start(0, 0, 0);
	ON_3dPoint Center_End(10, 10, 10);
	ON_3dVector Center_Tan_Start(1, 0, 0);
	ON_3dVector Center_Tan_End(-1, 1, 1);
	ON_BezierCurve Center_Rail = FixAxisBezier3D::Interpolate(Center_Start, Center_End, Center_Tan_Start, Center_Tan_End);

	for (int i = 0; i < N; ++i)
	{
		double t = double(i) / double(N - 1);
		ON_3dPoint Center = Center_Rail.PointAt(t);
		double R = (std::max)(abs(cos(t)), abs(sin(t))) * 3;
		ss.AddSphere(ON_Sphere(Center, R));
		double r = Red.FractionRed() * (1 - double(i) / double(N - 1)) + Green.FractionRed() * double(i) / double(N - 1);
		double g = Red.FractionGreen() * (1 - double(i) / double(N - 1)) + Green.FractionGreen() * double(i) / double(N - 1);
		double b = Red.FractionBlue() * (1 - double(i) / double(N - 1)) + Green.FractionBlue() * double(i) / double(N - 1);
		double max_rgb = (std::max)((std::max)(r, g), b);
		r = r / max_rgb;
		g = g / max_rgb;
		b = b / max_rgb;
		color[i].SetFractionalRGB(r, g, b);
		r = 1 - r; g = 1 - g; b = 1 - b;
		max_rgb = (std::max)((std::max)(r, g), b);
		r = r / max_rgb;
		g = g / max_rgb;
		b = b / max_rgb;
		resverse_color[i].SetFractionalRGB(r, g, b);
	}
	ON_NurbsSurface ons;
	for (int i = 0; i < N; ++i)
	{
		const int spheres_layer_index = model->AddLayer((L"Sphere" + std::to_wstring(i)).c_str(), color[i]);
		ss.GetSphere(i).GetNurbForm(ons);
		ChiralityAddNurbsSurface(model, ons, L"Sphere" + std::to_wstring(i), spheres_layer_index);
	}
	ss.GetCirclesToInterpolate();
	ss.GetConicVertex();
	std::vector<ON_Circle> param_circles = ss.Get_Circles_For_Debug();
	ON_NurbsCurve onc;
	int i = 0;
	for (const ON_Circle& c : param_circles)
	{
		const int circles_layer_index = model->AddLayer((L"Circle" + std::to_wstring(i)).c_str(), resverse_color[i]);
		c.GetNurbForm(onc);
		ChiralityAddNurbsCurve(model, onc, L"Circle" + std::to_wstring(i), circles_layer_index);
		++i;
	}
	std::vector<ON_NurbsSurface> v_ons = ss.Skinning();
	i = 0;
	for (const ON_NurbsSurface& ons : v_ons)
	{
		const int ons_layer_index = model->AddLayer((L"Surface" + std::to_wstring(i)).c_str(), ON_Color::SaturatedBlue);
		ChiralityAddNurbsSurface(model, ons, L"Surface" + std::to_wstring(i), ons_layer_index);
	}
}


int main()
{
    std::cout<<"SphereSkin\n";
	SetOutput("SphereSkin");
	const std::string filename = "SphereSkin.3dm";
	ON::Begin();
	ONX_Model model_to_write;
	Internal_SetExampleModelProperties(model_to_write, OPENNURBS__FUNCTION__, filename.c_str());
	model_to_write.AddDefaultLayer(nullptr, ON_Color::UnsetColor);
	SphereSkinTest(&model_to_write);
	ChiralityWrite3dmModel(&model_to_write, filename);
	BoundaryRecorder::GetRecorder().Write();
}