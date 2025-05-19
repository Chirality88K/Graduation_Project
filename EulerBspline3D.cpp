#include "EulerBspline3D.h"
#include "EulerBspline2D.h"
#include "write3dm.h"
#include "ChiralityMathTools.h"

EulerBspline3D::EulerBspline3D()
{
}

EulerBspline3D::EulerBspline3D(ON_3dVector vs, ON_3dVector ve, double length)
{
	vs.Unitize();
	ve.Unitize();
	ON_3dVector normal(0, -vs.z - ve.z, ve.y + vs.y);
	normal.Unitize();
	ON_3dVector new_vs = vs - ON_3dVector::DotProduct(vs, normal) * normal;
	ON_3dVector new_ve = ve - ON_3dVector::DotProduct(ve, normal) * normal;
	new_vs.Unitize();
	new_ve.Unitize();
	ON_3dVector real_yaxis = ON_3dVector::CrossProduct(normal, ON_3dVector(1, 0, 0));
	ON_Plane realplane(ON_3dPoint(0, 0, 0), ON_3dVector(1, 0, 0), real_yaxis);
	ON_Plane xyplane(ON_3dPoint(0, 0, 0), ON_3dVector(0, 0, 1));
	ON_Xform rotation = ON_Xform::IdentityTransformation;
	rotation.Rotation(realplane, xyplane);
	vs.Transform(rotation);
	ve.Transform(rotation);
	new_vs.Transform(rotation);
	new_ve.Transform(rotation);
	double tan_phi0 = vs.z / sqrt(vs.x * vs.x + vs.y * vs.y);
	double tan_phi1 = ve.z / sqrt(ve.x * ve.x + ve.y * ve.y);
	ON_NurbsCurve onc2d = ChiralityMath::UniformG1(ON_3dPoint::Origin, ON_3dPoint(length, 0.0, 0.0), new_vs, new_ve);
	EulerBspline2D::EulerBsplineSpiralInterpolation(onc2d);
	this->Create(3, false, onc2d.Order(), onc2d.CVCount());
	ON_3dPoint p;
	for (int i = 0; i < onc2d.CVCount(); ++i)
	{
		onc2d.GetCV(i, p);
		this->SetCV(i, ON_3dPoint(p.x, p.y, 0));
	}
	for (int i = 0; i < onc2d.KnotCount(); ++i)
	{
		this->SetKnot(i, onc2d.Knot(i));
	}
	ON_3dPoint p0, p1;
	GetCV(0, p0);
	GetCV(2, p1);
	double length_s = p0.DistanceTo(p1) / 2;
	GetCV(CVCount() - 3, p0);
	GetCV(CVCount() - 1, p1);
	double length_e = p0.DistanceTo(p1) / 2;
	double h1 = tan_phi0 * length_s;
	double hn_1 = tan_phi1 * length_e;
	double n = double(CVCount() - 1);
	double A = h1 / (2.0 - n);
	double B = -1.0 / 6.0 * A * (1.5 * n * n - 6.0 * n + 8.0);
	for (int i = 0; i < CVCount(); ++i)
	{
		GetCV(i, p);
		SetCV(i, ON_3dPoint(p.x, p.y, A * (i - n / 2.0) * (i - n / 2.0) + B));
	}
	ON_Xform inv_rotation = rotation.Inverse();
	this->Transform(inv_rotation);
}

EulerBspline3D::EulerBspline3D(ON_3dPoint ps, ON_3dPoint pe, ON_3dVector vs, ON_3dVector ve)
{
	ON_3dPoint p_e_ = pe - ps;
	ON_3dVector T = p_e_;
	T.Unitize();
	ON_Xform rotation = ON_Xform::IdentityTransformation;
	rotation.Rotation(T, ON_3dVector::XAxis, ON_3dPoint::Origin);
	double PEX = p_e_.DistanceTo(ON_3dPoint::Origin);
	ON_3dVector real_vs = vs;
	ON_3dVector real_ve = ve;
	real_vs.Transform(rotation);
	real_ve.Transform(rotation);
	EulerBspline3D eb3d = EulerBspline3D(real_vs, real_ve, PEX);
	*this = eb3d;
	ON_Xform inv_rotation = rotation.Inverse();
	this->Transform(inv_rotation);
	this->Translate(ps);
}

void EulerBspline3D::EulerBspline3DTest_MidPlaneMethod(ONX_Model *model)
{
	ON_3dVector VS(3, -1, 0);
	ON_3dVector VE(1, 1, 1);
	ON_3dPoint PS(0, 1, 0);
	ON_3dPoint PE(10, -2, 4);
	EulerBspline3D testebc(PS, PE, VS, VE);
	ON_NurbsCurve *noc = new ON_NurbsCurve(testebc);
	const int layer_index = model->AddLayer(L"3D_EulerBspline3D_Test", ON_Color::SaturatedMagenta);
	ON_3dmObjectAttributes *attributes = new ON_3dmObjectAttributes();
	attributes->m_layer_index = layer_index;
	attributes->m_name = L"curve";
	model->AddManagedModelGeometryComponent(noc, attributes);
	PrintCurvature(*noc, "EulerBspline3DTest_MidPlaneMethod_curvature");
	PrintPosAndTan(*noc, "Pos&Tan");
}
