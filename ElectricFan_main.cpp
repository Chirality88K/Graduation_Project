#include "FixAxisBezier3D.h"
#include "write3dm.h"
#include "ChiralityMathTools.h"
#include "EulerBezier2D.h"

static constexpr double inner_radius = 3.0;
static constexpr double outer_radius = 15.0;
static const ON_Circle circle1(ON_3dPoint(0, 0, 10), inner_radius);
static const ON_Circle circle2(ON_3dPoint(0, 0, 10), outer_radius);
static const ON_Circle circle3(ON_3dPoint(0, 0, 8), outer_radius * 0.9);
static const ON_Circle circle4(ON_3dPoint(0, 0, 8), outer_radius * 1.1);


void DrawCircles(ONX_Model* model)
{
	const int layer_index = model->AddLayer(L"Circles", ON_Color::SaturatedBlue);
	
	//circle1.GetNurbForm(onc);
	//ON_NurbsSurface ons1 = ChiralityMath::GenerateCylinder(onc, ON_3dVector::ZAxis, -2, 0);
	//ChiralityAddNurbsSurface(model, ons1, L"cylinder1", layer_index);
	ON_3dPoint corner_start = circle1.PointAt(0.0) - 0.15 * ON_3dVector::ZAxis;
	ON_3dPoint corner = circle1.PointAt(0.0);
	ON_3dPoint corner_end = circle1.PointAt(0.0) - (circle1.PointAt(0.0) - circle1.Center()) / circle1.Radius() * 0.15;
	corner_start.Rotate(1, 0, ON_3dVector::XAxis, ON_3dPoint::Origin);
	corner.Rotate(1, 0, ON_3dVector::XAxis, ON_3dPoint::Origin);
	corner_end.Rotate(1, 0, ON_3dVector::XAxis, ON_3dPoint::Origin);
	ON_NurbsCurve round_corner = EulerBezier2D::GenerateSmoothingCurve(corner_start, corner, corner_end);
	round_corner = ChiralityMath::ChangeDimensionFrom2To3(round_corner);
	round_corner.Rotate(-1, 0, ON_3dVector::XAxis, ON_3dPoint::Origin);
	ON_NurbsCurve line_onc1, line_onc2;
	ON_Line line1(circle1.PointAt(0.0) - 2 * ON_3dVector::ZAxis, round_corner.PointAtStart());
	ON_LineCurve line_curve1(line1);
	line_curve1.GetNurbForm(line_onc1);
	ON_Line line2(round_corner.PointAtEnd(),circle1.Center());
	ON_LineCurve line_curve2(line2);
	line_curve2.GetNurbForm(line_onc2);
	ON_NurbsCurve connected = line_onc1;
	connected.Append(round_corner);
	connected.Append(line_onc2);
	ON_NurbsSurface ons1 = ChiralityMath::GenerateRotating(connected, ON_Line(ON_3dPoint::Origin, ON_3dPoint(0, 0, 1)));
	ChiralityAddNurbsSurface(model, ons1, L"rotating1", layer_index);

	ON_3dPointArray points;
	points.Append(circle3.PointAt(0.0));
	points.Append(circle2.PointAt(0.0));
	points.Append(circle4.PointAt(0.0));
	ON_PolylineCurve poly_lines(points);
	ON_NurbsCurve onc;
	poly_lines.GetNurbForm(onc);
	ON_NurbsSurface ons2 = ChiralityMath::GenerateRotating(onc, ON_Line(ON_3dPoint::Origin, ON_3dPoint(0, 0, 1)));
	ChiralityAddNurbsSurface(model, ons2, L"rotating2", layer_index);
}

void DrawClothoid(ONX_Model* model)
{
	const ON_3dVector offset = 0.3 * ON_3dVector::ZAxis;
	ON_3dPoint ps = circle1.PointAt(0.0);
	ON_3dPoint pe = circle2.PointAt(PI / 6);
	ON_3dVector vs = ps - circle1.Center();
	ON_3dVector ve = vs;
	ve.Rotate(PI / 3, ON_3dVector::ZAxis);
	ps -= offset;
	pe -= offset;
	ON_BezierCurve obc1 = FixAxisBezier3D::Interpolate(ps, pe, vs, ve);
	ps -= 2 * ON_3dVector::ZAxis;
	ps += 2 * offset;
	pe = circle3.PointAt(PI / 6 - 0.1);
	pe += offset;
	ve.Rotate(-0.1, ON_3dVector::ZAxis);
	ON_BezierCurve obc2 = FixAxisBezier3D::Interpolate(ps, pe, vs, ve);
	while (obc1.CVCount() < obc2.CVCount())
	{
		ChiralityMath::Elevate(obc1);
	}
	while (obc1.CVCount() > obc2.CVCount())
	{
		ChiralityMath::Elevate(obc2);
	}
	ON_NurbsSurface ons(3, false, obc1.Order(), 2, obc1.CVCount(), 2);
	for (int i = 0; i < ons.KnotCount(0); ++i)
	{
		ons.SetKnot(0, i, (i < ons.KnotCount(0) / 2) ? 0 : 1);
	}
	ons.SetKnot(1, 0, 0); ons.SetKnot(1, 1, 1);
	for (int i = 0; i < obc1.CVCount(); ++i)
	{
		ons.SetCV(i, 0, obc1.ControlPoint(i));
		ons.SetCV(i, 1, obc2.ControlPoint(i));
	}

	const int layer_index = model->AddLayer(L"Clothoid_Surface", ON_Color::SaturatedMagenta);
	const int N = 24;
	for (int i = 0; i < N; ++i)
	{
		ons.Rotate(PI / 12, ON_3dVector::ZAxis, circle1.Center());
		ChiralityAddNurbsSurface(model, ons, L"surface" + std::to_wstring(i), layer_index);
	}
}


int main()
{
    std::cout << "ElectricFan\n";
	SetOutput("ElectricFan");
	const std::string filename = "ElectricFan.3dm";
	ON::Begin();
	ONX_Model model_to_write;
	Internal_SetExampleModelProperties(model_to_write, OPENNURBS__FUNCTION__, filename.c_str());
	model_to_write.AddDefaultLayer(nullptr, ON_Color::UnsetColor);
	DrawCircles(&model_to_write);
	DrawClothoid(&model_to_write);
	ChiralityWrite3dmModel(&model_to_write, filename);
}