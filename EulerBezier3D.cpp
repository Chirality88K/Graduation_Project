#include "EulerBezier3D.h"
#include "EulerBezier2D.h"
#include "thirdparty/eigen/Eigen/Dense"
#include "write3dm.h"
#include <fstream>
using namespace std;

const int iteration_times = 30;

EulerBezier3D::EulerBezier3D()
{
}

EulerBezier3D::EulerBezier3D(ON_3dVector vs, ON_3dVector ve, double length)
{
	vs.Unitize();
	ve.Unitize();
	ON_3dVector normal(0, -vs.z - ve.z, ve.y + vs.y);
	normal.Unitize();
	ON_3dVector new_vs = vs - ON_3dVector::DotProduct(vs, normal) * normal;
	ON_3dVector new_ve = ve - ON_3dVector::DotProduct(ve, normal) * normal;
	new_vs.Unitize();
	new_ve.Unitize();
	// ON_Plane realplane(ON_3dPoint(0, 0, 0), normal);
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
	ON_BezierCurve obc2d(2, false, 4);
	obc2d.SetCV(0, ON_3dPoint(0, 0, 0));
	obc2d.SetCV(1, ON_3dPoint(new_vs) * length / 3.0);
	obc2d.SetCV(2, ON_3dPoint(length, 0, 0) - (new_ve * length / 3.0));
	obc2d.SetCV(3, ON_3dPoint(length, 0, 0));
	EulerBezier2D::EulerBezierSpiralInterpolation(&obc2d);
	this->Create(3, false, obc2d.Order());
	ON_3dPoint p;
	for (int i = 0; i < obc2d.CVCount(); ++i)
	{
		obc2d.GetCV(i, p);
		this->SetCV(i, ON_3dPoint(p.x, p.y, 0));
	}
	ON_3dPoint p0, p1;
	GetCV(0, p0);
	GetCV(1, p1);
	double length_s = p0.DistanceTo(p1);
	GetCV(CVCount() - 2, p0);
	GetCV(CVCount() - 1, p1);
	double length_e = p0.DistanceTo(p1);
	double h1 = tan_phi0 * length_s;
	double hn_1 = tan_phi1 * length_e;
	double n = double(CVCount() - 1);
	// double B = (h1 * n - h1 - hn_1) / (h1 - hn_1);
	// double A = h1 / (1 - n) / (1 - B);
	double A = h1 / (1 - n);
	for (int i = 0; i < CVCount(); ++i)
	{
		GetCV(i, p);
		SetCV(i, ON_3dPoint(p.x, p.y, A * i * (i - n)));
	}
	ON_Xform inv_rotation = rotation.Inverse();
	this->Transform(inv_rotation);
}

EulerBezier3D::EulerBezier3D(ON_3dPoint ps, ON_3dPoint pe, ON_3dVector vs, ON_3dVector ve)
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
	EulerBezier3D eb3d = EulerBezier3D(real_vs, real_ve, PEX);
	*this = eb3d;
	ON_Xform inv_rotation = rotation.Inverse();
	this->Transform(inv_rotation);
	this->Translate(ps);
}

static int sign(double x)
{
	if (x > 1e-6)
	{
		return 1;
	}
	if (x < -1e-6)
	{
		return -1;
	}
	return 0;
}

void EulerBezier3D::Smoothing()
{
	ON_3dPoint p0, p1, p2;
	ON_3dVector t0, t1, b0, b1;
	b0 = CurvatureAt(0);
	b1 = CurvatureAt(1);
	for (int it = 0; it < iteration_times; ++it)
	{
		Elevate();
		vector<double> Length = ComputeLength();
		double sum = 0;
		for (double l : Length)
		{
			sum += l;
		}
		double nowlength = 0;
		vector<ON_3dVector> binormal;
		for (int i = 0; i < CVCount(); ++i)
		{
			binormal.push_back((1 - nowlength / sum) * b0 + nowlength / sum * b1);
			nowlength += (i == CVCount() - 1 ? 0 : Length[i]);
		}
		for (int i = 2; i < CVCount() - 2; ++i)
		{
			double S0 = 0.0, S1 = 0.0, S2 = 0.0;
			vector<double> Angle = ComputeAngle();
			Length = ComputeLength();
			S0 = Length[i - 2] * Length[i - 1] * sin(Angle[i - 1]) * 0.5;
			S1 = Length[i - 1] * Length[i - 0] * sin(Angle[i - 0]) * 0.5;
			S2 = Length[i - 0] * Length[i + 1] * sin(Angle[i + 1]) * 0.5;
			double k0 = 0.0, k1 = 0.0, k2 = 0.0;
			GetCV(i - 2, p0);
			GetCV(i - 0, p2);
			k0 = 4 * S0 / (Length[i - 2] * Length[i - 1] * p0.DistanceTo(p2));
			GetCV(i - 1, p0);
			GetCV(i + 1, p2);
			k1 = 4 * S1 / (Length[i - 1] * Length[i - 0] * p0.DistanceTo(p2));
			GetCV(i - 0, p0);
			GetCV(i + 2, p2);
			k2 = 4 * S2 / (Length[i - 0] * Length[i + 1] * p0.DistanceTo(p2));

			double mean_k = (k0 + k1 + k2) / 3.0;
			GetCV(i - 1, p0);
			GetCV(i, p1);
			GetCV(i + 1, p2);
			double p0p2_length = p0.DistanceTo(p2);
			double cos2theta = sqrt(1.0 - (std::min)(1.0, 0.25 * p0p2_length * p0p2_length * mean_k * mean_k));
			SetCV(i, (p0 + p2) / 2.0 - binormal[i] / binormal[i].Length() * (0.25 * p0p2_length * p0p2_length * mean_k / (1 + cos2theta)));
		}
		Length = ComputeLength();
		sum = 0;
		for (double l : Length)
		{
			sum += l;
		}
		sum /= Length.size();
		this->GetCV(0, p0);
		this->GetCV(CVCount() - 1, p1);
		t0 = TangentAt(0);
		t1 = TangentAt(1);
		SetCV(1, p0 + t0 * sum);
		SetCV(CVCount() - 2, p1 - t1 * sum);
	}
}

std::vector<double> EulerBezier3D::ComputeLength()
{
	vector<double> result;
	ON_3dPoint p0, p1;
	for (int i = 1; i < this->CVCount(); ++i)
	{
		this->GetCV(i - 1, p0);
		this->GetCV(i, p1);
		result.push_back(p0.DistanceTo(p1));
	}
	return result;
}

std::vector<double> EulerBezier3D::ComputeAngle()
{
	vector<double> result;
	ON_3dPoint p0, p1, p2;
	double product = 0;
	result.push_back(0);
	for (int i = 1; i < this->CVCount() - 1; ++i)
	{
		this->GetCV(i - 1, p0);
		this->GetCV(i, p1);
		this->GetCV(i + 1, p2);
		product = ON_3dVector::DotProduct(p1 - p0, p2 - p1) / p0.DistanceTo(p1) / p1.DistanceTo(p2);
		if (product > 1.0)
		{
			product = 1.0;
		}
		if (product < -1.0)
		{
			product = -1.0;
		}
		result.push_back(acos(product));
	}
	result.push_back(0);
	return result;
}

std::vector<ON_3dVector> EulerBezier3D::ComputeNormal()
{
	vector<ON_3dVector> result;
	ON_3dPoint p0, p1, p2;
	ON_3dVector n;
	result.push_back(ON_3dVector(0, 0, 0));
	for (int i = 1; i < this->CVCount(); ++i)
	{
		this->GetCV(i - 1, p0);
		this->GetCV(i, p1);
		this->GetCV(i + 1, p2);
		n = ON_3dVector::CrossProduct(p1 - p0, p2 - p1);
		if (n.Length() < 1e-8)
		{
			n = ON_3dVector(0, 0, 0);
		}
		else
		{
			n.Unitize();
		}
		result.push_back(n);
	}
	result.push_back(ON_3dVector(0, 0, 0));
	return result;
}

ON_3dPoint EulerBezier3D::GetBallCenter(ON_3dPoint p0, ON_3dPoint p1, ON_3dVector t0, ON_3dVector t1)
{
	t0.Unitize();
	t1.Unitize();
	ON_3dVector dir = p1 - p0;
	dir.Unitize();
	Eigen::Matrix3d A;
	A(0, 0) = t0.x;
	A(0, 1) = t0.y;
	A(0, 2) = t0.z;
	A(1, 0) = t1.x;
	A(1, 1) = t1.y;
	A(1, 2) = t1.z;
	A(2, 0) = dir.x;
	A(2, 1) = dir.y;
	A(2, 2) = dir.z;
	Eigen::Matrix<double, 3, 1> B;
	B(0, 0) = ON_3dVector::DotProduct(p0, t0);
	B(1, 0) = ON_3dVector::DotProduct(p1, t1);
	B(2, 0) = ON_3dVector::DotProduct((p0 + p1) / 2.0, dir);
	Eigen::Matrix<double, 3, 1> X;
	X = A.partialPivLu().solve(B);
	return ON_3dPoint(X(0, 0), X(1, 0), X(2, 0));
}

void EulerBezier3D::Elevate()
{
	int n = this->CVCount();
	if (n < 2)
	{
		return;
	}
	std::vector<ON_3dPoint> parr;
	ON_3dPoint p, q;
	GetCV(0, p);
	parr.push_back(p);
	for (int i = 1; i < n; i++)
	{
		GetCV(i - 1, p);
		GetCV(i, q);
		parr.push_back(p * (i * 1.0 / n) + q * (1 - i * 1.0 / n));
	}
	GetCV(n - 1, q);
	parr.push_back(q);
	Create(3, false, n + 1);
	for (int i = 0; i < n + 1; i++)
	{
		SetCV(i, parr[i]);
	}
}

void EulerBezier3D::Rise_to_3D(ON_NurbsCurve *onc, ON_3dPoint pe, ON_3dVector ts, ON_3dVector te)
{
	ts.Unitize();
	te.Unitize();
	ON_3dVector L = pe;
	L.Unitize();
	double start = 0.0, end = 0.0;
	onc->GetDomain(&start, &end);
	ON_3dVector t0 = onc->TangentAt(start);
	ON_3dVector t1 = onc->TangentAt(end);
	double product0 = ON_3dVector::DotProduct(t0, ts);
	product0 = (std::max)(product0, -1.0);
	product0 = (std::min)(product0, 1.0);
	double theta0 = acos(product0);
	ON_3dPoint p0, p1;
	onc->GetCV(0, p0);
	onc->GetCV(1, p1);
	double h1 = tan(theta0) * p0.DistanceTo(p1);
	if (ON_3dVector::DotProduct(L, ts) < 0)
	{
		h1 = -h1;
	}
	double hn = pe.DistanceTo(ON_3dPoint(0, 0, 0));
	int n = onc->CVCount() - 1;

	double product1 = ON_3dVector::DotProduct(t1, te);
	product1 = (std::max)(product1, -1.0);
	product1 = (std::min)(product1, 1.0);
	double theta1 = acos(product1);
	onc->GetCV(n - 1, p0);
	onc->GetCV(n, p1);
	double hn_1 = tan(theta1) * p0.DistanceTo(p1);
	if (ON_3dVector::DotProduct(L, te) > 0)
	{
		hn_1 = -hn_1;
	}
	hn_1 += hn;
	Eigen::Matrix<double, 3, 3> A;
	A(0, 0) = (1.0 / n / n / n);
	A(0, 1) = (1.0 / n / n);
	A(0, 2) = (1.0 / n);
	A(1, 0) = (1.0 - 1.0 / n) * (1.0 - 1.0 / n) * (1.0 - 1.0 / n);
	A(1, 1) = (1.0 - 1.0 / n) * (1.0 - 1.0 / n);
	A(1, 2) = (1.0 - 1.0 / n);
	A(2, 0) = 1;
	A(2, 1) = 1;
	A(2, 2) = 1;
	Eigen::Matrix<double, 3, 1> B;
	B(0, 0) = h1;
	B(1, 0) = hn_1;
	B(2, 0) = hn;
	Eigen::Matrix<double, 3, 1> X;
	X = A.partialPivLu().solve(B);
	for (int i = 0; i <= n; ++i)
	{
		double k = double(i) / double(n);
		onc->GetCV(i, p0);
		p0 = ON_3dPoint(p0.x, p0.y, X(0, 0) * k * k * k + X(1, 0) * k * k + k * X(2, 0));
		onc->SetCV(i, p0);
	}
}

void EulerBezier3D::EulerBezier3DTest(ONX_Model *model)
{
	/*
	EulerBezier3D *obc = new EulerBezier3D();
	obc->Create(3, false, 5);
	obc->SetCV(0, ON_3dPoint(0, 0, 0));
	obc->SetCV(1, ON_3dPoint(0, -5, 4));
	obc->SetCV(2, ON_3dPoint(8, -5, -5));
	obc->SetCV(3, ON_3dPoint(10, 6, -3));
	obc->SetCV(4, ON_3dPoint(10, 6, 2));
	ON_NurbsCurve *onc_origin = new ON_NurbsCurve();
	obc->GetNurbForm(*onc_origin);
	const int layer_index = model->AddLayer(L"Layer1", ON_Color::Black);
	ON_3dmObjectAttributes *attributes = new ON_3dmObjectAttributes();
	attributes->m_layer_index = layer_index;
	attributes->m_name = L"origin";
	model->AddManagedModelGeometryComponent(onc_origin, attributes);
	obc->Smoothing();
	ON_NurbsCurve *onc_after = new ON_NurbsCurve();
	obc->GetNurbForm(*onc_after);
	const int layer_index__ = model->AddLayer(L"Layer2", ON_Color::SaturatedMagenta);
	ON_3dmObjectAttributes *attributes__ = new ON_3dmObjectAttributes();
	attributes__->m_layer_index = layer_index__;
	attributes__->m_name = L"after";
	model->AddManagedModelGeometryComponent(onc_after, attributes__);
	delete obc;
	*/

	ON_BezierCurve *obc = new ON_BezierCurve(2, false, 4);
	obc->SetCV(0, ON_3dPoint(0, 0, 0));
	obc->SetCV(1, ON_3dPoint(5, 0, 0));
	obc->SetCV(2, ON_3dPoint(25.0 / 4.0, 5.0 / 4.0 * sqrt(3.0), 0));
	obc->SetCV(3, ON_3dPoint(2.5, 2.5 * sqrt(3), 0));
	EulerBezier2D::EulerBezierSpiralInterpolation(obc);
	const int layer_index1 = model->AddLayer(L"2D_curve", ON_Color::Black);
	// const int layer_index2 = model->AddLayer(L"3D_curve", ON_Color::SaturatedMagenta);
	ON_NurbsCurve *onc1 = new ON_NurbsCurve();
	ON_NurbsCurve *onc2 = new ON_NurbsCurve();
	obc->GetNurbForm(*onc1);

	ON_BezierCurve *sys_obc = new ON_BezierCurve();
	EulerBezier2D::GenerateSymmetry(sys_obc, obc, ON_3dPoint(0, 0, 0), ON_3dVector(0.5, 0.5 * sqrt(3), 0));
	sys_obc->GetNurbForm(*onc2);
	onc2->Reverse();
	// ON_3dmObjectAttributes* attributes2 = new ON_3dmObjectAttributes();
	// attributes2->m_layer_index = layer_index;
	// attributes2->m_name = L"part2";
	// model->AddManagedModelGeometryComponent(onc2, attributes2);
	onc1->Append(*onc2);
	PrintCurvature(*onc1, "2DCurvesCurvature");
	ON_3dmObjectAttributes *attributes1 = new ON_3dmObjectAttributes();
	attributes1->m_layer_index = layer_index1;
	attributes1->m_name = L"part1";
	model->AddManagedModelGeometryComponent(onc1, attributes1);
	delete obc;
	delete sys_obc;
	delete onc2;
	// Rise

	ON_NurbsCurve *onc3 = new ON_NurbsCurve();
	onc3->Create(3, false, onc1->Order(), onc1->CVCount());
	ON_3dPoint P;
	for (int i = 0; i < onc1->CVCount(); ++i)
	{
		onc1->GetCV(i, P);
		onc3->SetCV(i, P);
	}
	const double *k0 = onc1->Knot();
	for (int i = 0; i < onc1->KnotCount(); ++i)
	{
		onc3->SetKnot(i, *(k0 + i));
	}
	Rise_to_3D(onc3, ON_3dPoint(0, 0, 10), ON_3dVector(5, 0, -1), ON_3dVector(1, -sqrt(3), 1));
	PrintCurvature(*onc3, "3DCurvesCurvature");
	ON_3dmObjectAttributes *attributes2 = new ON_3dmObjectAttributes();
	attributes2->m_layer_index = layer_index1;
	attributes2->m_name = L"part2";
	model->AddManagedModelGeometryComponent(onc3, attributes2);
}

void EulerBezier3D::EulerBezier3DTest_MidPlaneMethod(ONX_Model *model)
{
	ON_3dVector VS(3, -1, 0);
	ON_3dVector VE(1, 1, 1);
	ON_3dPoint PS(0, 1, 0);
	ON_3dPoint PE(10, -2, 4);
	EulerBezier3D testebc(PS, PE, VS, VE);
	ON_NurbsCurve *noc = new ON_NurbsCurve();
	testebc.GetNurbForm(*noc);
	const int layer_index = model->AddLayer(L"3D_EulerBezier_Test", ON_Color::SaturatedMagenta);
	ON_3dmObjectAttributes *attributes = new ON_3dmObjectAttributes();
	attributes->m_layer_index = layer_index;
	attributes->m_name = L"curve";
	model->AddManagedModelGeometryComponent(noc, attributes);
	PrintCurvature(*noc, "EulerBezier3DTest_MidPlaneMethod_curvature");
	PrintPosAndTan(*noc, "Pos&Tan");
}
