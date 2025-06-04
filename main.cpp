#ifdef WIN32
#pragma comment(lib, "Shlwapi.lib")
#endif
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include "FilletSurface.h"
#include "PiecesClothoid3d.h"
#include "Cone_Surface.h"
#include "Spiral.h"
#include "Cornu_Spiral.h"
#include <algorithm>
#include "EulerBezier2D.h"
#include "EulerBspline2D.h"
#include "EulerBezier3D.h"
#include "EulerBspline3D.h"
#include "read3dm.h"
#include "write3dm.h"
#include "Fillet_using_EB3d.h"
#include "EulerPolygon3D.h"
const double PI = acos(-1.0);

using namespace std;

void Add_Cylinder(ONX_Model *model, const wchar_t *name, double radius = 1, double height = 2, ON_Color color = ON_Color::Black, ON_Xform xform = ON_Xform::IdentityTransformation);
void Add_Circle(ONX_Model *model, const wchar_t *name, double radius = 1, ON_Color color = ON_Color::Black, ON_Xform xform = ON_Xform::IdentityTransformation);
void Add_Ball(ONX_Model *model, const wchar_t *name, double radius = 1, ON_Color color = ON_Color::Black, ON_Xform xform = ON_Xform::IdentityTransformation);
void Add_BezierCurve(ONX_Model *model, const wchar_t *name, const ON_BezierCurve *bc, ON_Color color = ON_Color::Black, ON_Xform xform = ON_Xform::IdentityTransformation);
void Get_Curvature_and_Torsion_of_NurbsCurve(const ON_NurbsCurve &onc, string filename);

// int main ( int argc, const char* argv[] )
int main()
{
	const std::string filename = "EulerPolygonTest-" + ChiralityPrintNowTime() + ".3dm";
	ON::Begin();
	ONX_Model model_to_write;
	Internal_SetExampleModelProperties(model_to_write, OPENNURBS__FUNCTION__, filename.c_str());
	model_to_write.AddDefaultLayer(nullptr, ON_Color::UnsetColor);

	// Cornu_Spiral::Cornu_test(&model_to_write);
	// EulerBspline2D::EulerBsplineTest(&model_to_write);
	// EulerBezier2D::EulerBezier2dTest(&model_to_write);
	//  EulerBezier2D::YangMethodtest(&model_to_write);
	//  EulerBezier3D::EulerBezier3DTest(&model_to_write);
	// EulerBezier3D::EulerBezier3DTest_MidPlaneMethod(&model_to_write);
	// EulerBspline3D::EulerBspline3DTest_MidPlaneMethod(&model_to_write);
	// Fillet_EB3D::Fillet_EB3D_Test(&model_to_write);
	// EulerBezier2D::Pentagram(&model_to_write);
	// Fillet_EB3D::TwoSurfaces_Fillet_Test(&model_to_write);
	// EulerBspline2D::SmoothCornerTest(&model_to_write);
	EulerPolygon3D::EulerPolygonTest(&model_to_write);

	ChiralityWrite3dmModel(&model_to_write, filename);
	ONX_Model model_to_read;
	Internal_SetExampleModelProperties(model_to_read, OPENNURBS__FUNCTION__, filename.c_str());
	ChiralityRead3dmModel(filename, &model_to_read);
	ChiralityModelDebugInfo(filename);
	return 0;
}

void Add_Cylinder(ONX_Model *model, const wchar_t *name, double radius, double height, ON_Color color, ON_Xform xform)
{
	const int bIsRational = true;
	const int dim = 3;
	const int u_degree = 2;
	const int v_degree = 1;
	const int u_cv_count = 7;
	const int v_cv_count = 2;
	double u_knot[u_cv_count + u_degree - 1] = {0, 0, 1.0 / 3, 1.0 / 3, 2.0 / 3, 2.0 / 3, 1, 1};
	double v_knot[v_cv_count + v_degree - 1] = {0, 1};
	ON_3dPoint CV[u_cv_count][v_cv_count];
	for (int k = -1; k < 2; k = k + 2)
	{
		int y = (k + 1) / 2;
		CV[0][y].Set(0, -radius, k * height / 2);
		CV[1][y].Set(-radius * sqrt(3), -radius, k * height / 2);
		CV[2][y].Set(-radius * sqrt(3) / 2, radius / 2, k * height / 2);
		CV[3][y].Set(0, 2 * radius, k * height / 2);
		CV[4][y].Set(radius * sqrt(3) / 2, radius / 2, k * height / 2);
		CV[5][y].Set(radius * sqrt(3), -radius, k * height / 2);
		CV[6][y].Set(0, -radius, k * height / 2);
	}
	double w[7] = {1, 0.5, 1, 0.5, 1, 0.5, 1};

	ON_NurbsSurface nurbs_surface(dim, bIsRational,
								  u_degree + 1, v_degree + 1,
								  u_cv_count, v_cv_count);
	int i, j;
	for (i = 0; i < nurbs_surface.KnotCount(0); i++)
	{
		nurbs_surface.SetKnot(0, i, u_knot[i]);
	}

	for (j = 0; j < nurbs_surface.KnotCount(1); j++)
	{
		nurbs_surface.SetKnot(1, j, v_knot[j]);
	}

	for (i = 0; i < nurbs_surface.CVCount(0); i++)
	{
		for (j = 0; j < nurbs_surface.CVCount(1); j++)
		{
			nurbs_surface.SetCV(i, j, CV[i][j] * w[i]);
			nurbs_surface.SetWeight(i, j, w[i]);
		}
	}
	nurbs_surface.Transform(xform);
	ON_NurbsSurface **nubs = new ON_NurbsSurface *();
	*nubs = new ON_NurbsSurface(nurbs_surface);
	model->AddDefaultLayer(nullptr, ON_Color::UnsetColor);
	const int layer_index = model->AddLayer(name, color);
	model->AddManagedModelGeometryComponent(
		*nubs,
		Internal_CreateManagedAttributes(layer_index, name));
	delete nubs;
}

void Add_Circle(ONX_Model *model, const wchar_t *name, double radius, ON_Color color, ON_Xform xform)
{
	ON_NurbsCurve *wiggle = new ON_NurbsCurve(
		3,	  // dimension
		true, // true if rational
		3,	  // order = degree+1
		7	  // number of control vertices
	);
	int i;

	ON_3dPoint CV[7];
	CV[0].Set(0, -radius, 0);
	CV[1].Set(-radius * sqrt(3), -radius, 0);
	CV[2].Set(-radius * sqrt(3) / 2, radius / 2, 0);
	CV[3].Set(0, 2 * radius, 0);
	CV[4].Set(radius * sqrt(3) / 2, radius / 2, 0);
	CV[5].Set(radius * sqrt(3), -radius, 0);
	CV[6].Set(0, -radius, 0);
	double weight[7] = {1, 0.5, 1, 0.5, 1, 0.5, 1};

	for (i = 0; i < 7; i++)
	{
		wiggle->SetCV(i, CV[i] * weight[i]);
		wiggle->SetWeight(i, weight[i]);
	}

	wiggle->SetKnot(0, 0.0);
	wiggle->SetKnot(1, 0.0);
	wiggle->SetKnot(2, 1.0 / 3);
	wiggle->SetKnot(3, 1.0 / 3);
	wiggle->SetKnot(4, 2.0 / 3);
	wiggle->SetKnot(5, 2.0 / 3);
	wiggle->SetKnot(6, 1.0);
	wiggle->SetKnot(7, 1.0);

	// model->AddManagedModelGeometryComponent(wiggle, nullptr);
	wiggle->Transform(xform);
	ON_NurbsCurve **nubs = new ON_NurbsCurve *();
	*nubs = new ON_NurbsCurve(*wiggle);
	model->AddDefaultLayer(nullptr, ON_Color::UnsetColor);
	const int layer_index = model->AddLayer(name, color);
	model->AddManagedModelGeometryComponent(
		*nubs,
		Internal_CreateManagedAttributes(layer_index, name));
	delete nubs;
	delete wiggle;
}

void Add_Ball(ONX_Model *model, const wchar_t *name, double radius, ON_Color color, ON_Xform xform)
{
	const int bIsRational = true;
	const int dim = 3;
	const int u_degree = 2;
	const int v_degree = 2;
	const int u_cv_count = 7;
	const int v_cv_count = 4;
	double u_knot[u_cv_count + u_degree - 1] = {0, 0, 0.25, 0.5, 0.5, 0.75, 1, 1};
	double v_knot[v_cv_count + v_degree - 1] = {0, 0, 0.5, 1, 1};
	ON_3dPoint CV[u_cv_count][v_cv_count];
	double weight[u_cv_count][v_cv_count];
	weight[0][0] = 1;
	weight[0][1] = 0.5;
	weight[0][2] = 0.5;
	weight[0][3] = 1;

	int i, j;
	for (i = 0; i < u_cv_count; i++)
	{
		if (i == 0 || i == 3 || i == 6)
		{
			for (j = 0; j < v_cv_count; j++)
			{
				weight[i][j] = weight[0][j];
			}
		}
		else
		{
			for (j = 0; j < v_cv_count; j++)
			{
				weight[i][j] = 0.5 * weight[0][j];
			}
		}
	}
	for (i = 0; i < u_cv_count; i++)
	{
		CV[i][0] = ON_3dPoint(0, 0, -radius) * weight[i][0];
		CV[i][3] = ON_3dPoint(0, 0, radius) * weight[i][3];
	}
	CV[0][1] = ON_3dPoint(radius, 0, -radius) * weight[0][1];
	CV[1][1] = ON_3dPoint(radius, radius, -radius) * weight[1][1];
	CV[2][1] = ON_3dPoint(-radius, radius, -radius) * weight[2][1];
	CV[3][1] = ON_3dPoint(-radius, 0, -radius) * weight[3][1];
	CV[4][1] = ON_3dPoint(-radius, -radius, -radius) * weight[4][1];
	CV[5][1] = ON_3dPoint(radius, -radius, -radius) * weight[5][1];
	CV[6][1] = ON_3dPoint(radius, 0, -radius) * weight[6][1];

	CV[0][2] = ON_3dPoint(radius, 0, radius) * weight[0][2];
	CV[1][2] = ON_3dPoint(radius, radius, radius) * weight[1][2];
	CV[2][2] = ON_3dPoint(-radius, radius, radius) * weight[2][2];
	CV[3][2] = ON_3dPoint(-radius, 0, radius) * weight[3][2];
	CV[4][2] = ON_3dPoint(-radius, -radius, radius) * weight[4][2];
	CV[5][2] = ON_3dPoint(radius, -radius, radius) * weight[5][2];
	CV[6][2] = ON_3dPoint(radius, 0, radius) * weight[6][2];

	ON_NurbsSurface nurbs_surface(dim, bIsRational,
								  u_degree + 1, v_degree + 1,
								  u_cv_count, v_cv_count);
	for (i = 0; i < nurbs_surface.KnotCount(0); i++)
	{
		nurbs_surface.SetKnot(0, i, u_knot[i]);
	}

	for (j = 0; j < nurbs_surface.KnotCount(1); j++)
	{
		nurbs_surface.SetKnot(1, j, v_knot[j]);
	}

	for (i = 0; i < nurbs_surface.CVCount(0); i++)
	{
		for (j = 0; j < nurbs_surface.CVCount(1); j++)
		{
			nurbs_surface.SetCV(i, j, CV[i][j]);
			nurbs_surface.SetWeight(i, j, weight[i][j]);
		}
	}
	nurbs_surface.Transform(xform);
	ON_NurbsSurface **nubs = new ON_NurbsSurface *();
	*nubs = new ON_NurbsSurface(nurbs_surface);
	model->AddDefaultLayer(nullptr, ON_Color::UnsetColor);
	const int layer_index = model->AddLayer(name, color);
	model->AddManagedModelGeometryComponent(
		*nubs,
		Internal_CreateManagedAttributes(layer_index, name));
	delete nubs;
}

void Add_BezierCurve(ONX_Model *model, const wchar_t *name, const ON_BezierCurve *bc, ON_Color color, ON_Xform xform)
{
	ON_NurbsCurve **onc = new ON_NurbsCurve *();
	(*onc) = new ON_NurbsCurve();
	bc->GetNurbForm(**onc);
	(*onc)->Transform(xform);
	const int layer_index = model->AddLayer(name, color);
	model->AddManagedModelGeometryComponent(
		*onc,
		Internal_CreateManagedAttributes(layer_index, name));
	delete onc;
}

void Get_Curvature_and_Torsion_of_NurbsCurve(const ON_NurbsCurve &onc, string filename)
{
	const double *knot = onc.Knot();
	int size = onc.KnotCount();
	double k0 = *knot;
	double kn = *(knot + size - 1);
	double t = 0;
	double kappa = 0;
	// double torsion = 0;
	ofstream ofs(filename);
	ON_3dPoint r;
	ON_3dVector r1;
	ON_3dVector r2;
	// ON_3dVector r3;
	for (int i = 0; i <= 1000; i++)
	{
		t = (kn - k0) / 1000 * i + k0;
		onc.Ev2Der(t, r, r1, r2);
		kappa = (ON_3dVector::CrossProduct(r1, r2)).Length() / (pow(r1.Length(), 3));
		ofs << t << " " << kappa << endl;
	}
	ofs.close();
}