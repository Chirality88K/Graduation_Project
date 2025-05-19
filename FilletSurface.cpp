#include "FilletSurface.h"
#include "thirdparty/eigen/Eigen/Dense"

FilletSurface::~FilletSurface()
{
	delete Rail0VF;
	delete Rail1VF;
}

void FilletSurface::SetSurface(int i, const ON_NurbsSurface &nurbs)
{
	if (i == 0)
	{
		m_Surf0 = nurbs;
	}
	if (i == 1)
	{
		m_Surf1 = nurbs;
	}
}

void FilletSurface::SetParameter(InsectCurve0Knot k0, double par0, InsectCurve1Knot k1, double par1)
{
	m_par0 = k0;
	m_par1 = k1;
	InsectCurve0Par = par0;
	InsectCurve1Par = par1;
}

void FilletSurface::SetRailCurve(double offset0, double offset1)
{
	RailCurvePar0 = InsectCurve0Par + offset0;
	RailCurvePar1 = InsectCurve1Par + offset1;
	if (m_par0 == 0)
	{
		Rail0VF = new VectorField1d(&m_Surf0, true, RailCurvePar0);
	}
	else
	{
		Rail0VF = new VectorField1d(&m_Surf0, false, RailCurvePar0);
	}
	if (m_par1 == 0)
	{
		Rail1VF = new VectorField1d(&m_Surf1, true, RailCurvePar1);
	}
	else
	{
		Rail1VF = new VectorField1d(&m_Surf1, false, RailCurvePar1);
	}
}

void FilletSurface::ReverseRailCurve(int i)
{
	if (i == 0)
	{
		Rail0VF->ReverseCurve();
		return;
	}
	Rail1VF->ReverseCurve();
}

void FilletSurface::ReverseVector(int i)
{
	if (i == 0)
	{
		Rail0VF->ReverseVector();
		return;
	}
	Rail1VF->ReverseVector();
}

void FilletSurface::ComputeFillCurve(int num_knot, double *knot_on_rail0, double *knot_on_rail1)
{
	FillCurveArray.clear();
	Knot_on_Rail0.clear();
	Knot_on_Rail1.clear();
	ON_3dPoint P0 = ON_3dPoint::Origin;
	ON_3dPoint Pn = ON_3dPoint::Origin;
	ON_3dVector V0 = ON_3dVector::XAxis;
	ON_3dVector Vn = ON_3dVector::XAxis;
	ON_3dVector d = ON_3dVector::XAxis;
	for (int i = 0; i < num_knot; i++)
	{
		Knot_on_Rail0.push_back(knot_on_rail0[i]);
		Knot_on_Rail1.push_back(knot_on_rail1[i]);
		Rail0VF->GetAll(knot_on_rail0[i], P0, V0, d);
		Rail1VF->GetAll(knot_on_rail1[i], Pn, Vn, d);
		ON_BezierCurve obc(3, false, 4);
		obc.SetCV(0, P0);
		obc.SetCV(1, P0 + V0);
		obc.SetCV(2, Pn - Vn);
		obc.SetCV(3, Pn);
		FillCurveArray.push_back(obc);
	}
}

ON_3dPoint FilletSurface::Point_at(int sur_index, double u, double v)
{
	ON_BezierCurve FillCurve0 = FillCurveArray[sur_index];
	ON_BezierCurve FillCurve1 = FillCurveArray[sur_index + 1];
	double trueu0 = (Knot_on_Rail0[sur_index + 1] - Knot_on_Rail0[sur_index]) * u + Knot_on_Rail0[sur_index];
	double trueu1 = (Knot_on_Rail1[sur_index + 1] - Knot_on_Rail1[sur_index]) * u + Knot_on_Rail1[sur_index];

	double phi0u = 2 * u * u * u - 3 * u * u + 1;
	double phid0u = u * u * u - 2 * u * u + u;
	double phi1u = 2 * (1 - u) * (1 - u) * (1 - u) - 3 * (1 - u) * (1 - u) + 1;
	double phid1u = (1 - u) * (1 - u) * (1 - u) - 2 * (1 - u) * (1 - u) + 1 - u;

	double phi0v = 2 * v * v * v - 3 * v * v + 1;
	double phid0v = v * v * v - 2 * v * v + v;
	double phi1v = 2 * (1 - v) * (1 - v) * (1 - v) - 3 * (1 - v) * (1 - v) + 1;
	double phid1v = (1 - v) * (1 - v) * (1 - v) - 2 * (1 - v) * (1 - v) + 1 - v;

	ON_3dPoint P0v, P1v, Pu0, Pu1;
	ON_3dVector Pu0v, Pu1v, Pvu0, Pvu1;
	ON_3dVector rubbish;
	Rail0VF->GetAll(trueu0, Pu0, Pvu0, rubbish);
	Rail1VF->GetAll(trueu1, Pu1, Pvu1, rubbish);
	P0v = FillCurve0.PointAt(v);
	P1v = FillCurve1.PointAt(v);

	ON_3dPoint P00, P01, P10, P11;
	ON_3dVector Pu00, Pv00, Pu01, Pv01, Pu10, Pv10, Pu11, Pv11;
	Rail0VF->GetAll(Knot_on_Rail0[sur_index], P00, Pv00, Pu00);
	Rail0VF->GetAll(Knot_on_Rail0[sur_index + 1], P10, Pv10, Pu10);
	Rail1VF->GetAll(Knot_on_Rail1[sur_index], P01, Pv01, Pu01);
	Rail1VF->GetAll(Knot_on_Rail1[sur_index + 1], P11, Pv11, Pu11);

	if (Pu00.IsParallelTo(Pu01) == 1)
	{
		Pu0v = Pu00 * (1 - v) + Pu01 * v;
		Pu0v.Unitize();
	}
	else
	{
		ON_3dVector n = ON_CrossProduct(Pu00, Pu01);
		n.Unitize();
		double product = ON_DotProduct(Pu00, Pu01);
		double theta = acos(product);
		ON_Xform Rotate;
		Rotate.Rotation(theta * v, n, ON_3dPoint(0, 0, 0));
		Pu0v = Pu00;
		Pu0v.Transform(Rotate);
	}
	if (Pu10.IsParallelTo(Pu11) == 1)
	{
		Pu1v = Pu10 * (1 - v) + Pu11 * v;
		Pu1v.Unitize();
	}
	else
	{
		ON_3dVector n = ON_CrossProduct(Pu10, Pu11);
		n.Unitize();
		double product = ON_DotProduct(Pu10, Pu11);
		double theta = acos(product);
		ON_Xform Rotate;
		Rotate.Rotation(theta * v, n, ON_3dPoint(0, 0, 0));
		Pu1v = Pu10;
		Pu1v.Transform(Rotate);
	}

	ON_3dPoint IuP = P0v * phi0u + Pu0v * phid0u + P1v * phi1u + Pu1v * phid1u;
	ON_3dPoint IvP = Pu0 * phi0v + Pvu0 * phid0v + Pu1 * phi1v + Pvu1 * phid1v;
	/*
	ON_3dVector Puv00, Puv01, Puv10, Puv11;
	Puv00 = Rail0VF->MixedDer(0);
	Puv10 = Rail0VF->MixedDer(1);
	Puv01 = Rail1VF->MixedDer(0);
	Puv11 = Rail1VF->MixedDer(1);

	ON_3dPoint IuIvP = P00 * phi0u * phi0v + P01 * phi0u * phi1v + P10 * phi1u * phi0v +
		P11 * phi1u * phi1v +
		Pv00 * phi0u * phid0v + Pu00 * phid0u * phi0v +
		Pv01 * phi0u * phid1v + Pu10 * phid1u * phi0v +
		Pv10 * phi1u * phid0v + Pu01 * phid0u * phi1v +
		Pv11 * phi1u * phid1v + Pu11 * phid1u * phi1v +
		Puv00 * phid0u * phid0v +
		Puv01 * phid0u * phid1v + Puv10 * phid1u * phid0v + Puv11 * phid1u * phid1v;
		*/
	double alpha = v * v * (1 - v) * (1 - v) / (v * v * (1 - v) * (1 - v) + u * u * (1 - u) * (1 - u));
	double beta = 1 - alpha;
	if (u * (1 - u) == 0 && v * (1 - v) != 0)
	{
		alpha = 1;
		beta = 0;
	}
	if (u * (1 - u) != 0 && v * (1 - v) == 0)
	{
		alpha = 0;
		beta = 1;
	}
	if (u * (1 - u) == 0 && v * (1 - v) == 0)
	{
		alpha = 0.5;
		beta = 0.5;
	}
	return alpha * IuP + beta * IvP;
}

void FilletSurface::Add_FilletSurface(ONX_Model *model, const wchar_t *name, int FillCurveSample_num, int EveryPatchSample_num, ON_Color color)
{
	model->AddLayer(L"FilletSurface_mesh", color);
	bool bHasVertexNormals = false;
	bool bHasTexCoords = false;
	const int vertex_count = FillCurveSample_num * EveryPatchSample_num * (FillCurveArray.size() - 1);
	const int face_count = (FillCurveSample_num - 1) * (EveryPatchSample_num - 1) * (FillCurveArray.size() - 1);
	ON_Mesh **mesh = new ON_Mesh *();
	*mesh = new ON_Mesh(face_count, vertex_count, bHasVertexNormals, bHasTexCoords);
	for (int i = 0; i < vertex_count; i++)
	{
		(*mesh)->SetVertex(i, ON_3dPoint(0, 0, 0));
	}
	for (int i = 0; i < FillCurveArray.size() - 1; i++)
	{
		for (int j = 0; j < EveryPatchSample_num; j++)
		{
			for (int k = 0; k < FillCurveSample_num; k++)
			{
				(*mesh)->SetVertex(i * EveryPatchSample_num * FillCurveSample_num + k * EveryPatchSample_num + j,
								   Point_at(i, j / (EveryPatchSample_num - 1.0), k / (FillCurveSample_num - 1.0)));
			}
		}
	}
	int index = 0;
	for (int i = 0; i < FillCurveArray.size() - 1; i++)
	{
		for (int j = 0; j < EveryPatchSample_num - 1; j++)
		{
			for (int k = 0; k < FillCurveSample_num - 1; k++)
			{
				(*mesh)->SetQuad(index,
								 i * EveryPatchSample_num * FillCurveSample_num + k * EveryPatchSample_num + j,
								 i * EveryPatchSample_num * FillCurveSample_num + k * EveryPatchSample_num + j + 1,
								 i * EveryPatchSample_num * FillCurveSample_num + (k + 1) * EveryPatchSample_num + j + 1,
								 i * EveryPatchSample_num * FillCurveSample_num + (k + 1) * EveryPatchSample_num + j);
				index++;
			}
		}
	}
	if (!(*mesh)->HasVertexNormals())
	{
		(*mesh)->ComputeVertexNormals();
	}
	const int layer_index = model->AddLayer(name, color);
	ON_3dmObjectAttributes *attributes = new ON_3dmObjectAttributes();
	attributes->m_layer_index = layer_index;
	attributes->m_name = name;
	model->AddManagedModelGeometryComponent(*mesh, attributes);
	delete mesh;
}

void FilletSurface::Add_FillCurves(ONX_Model *model, const wchar_t *name, ON_Color color)
{
	const int layer_index = model->AddLayer(name, color);
	ON_3dmObjectAttributes *attributes = new ON_3dmObjectAttributes();
	attributes->m_layer_index = layer_index;
	attributes->m_name = name;
	for (unsigned int i = 0; i < FillCurveArray.size(); i++)
	{
		ON_NurbsCurve **onc = new ON_NurbsCurve *();
		(*onc) = new ON_NurbsCurve();
		FillCurveArray[i].GetNurbForm(**onc);
		model->AddManagedModelGeometryComponent(*onc, attributes);
		delete onc;
	}
}

VectorField1d::VectorField1d(ON_NurbsSurface *osf, bool whichone, double s)
{
	m_source_surface = new ON_NurbsSurface(*osf);
	m_source_surface->SetDomain(0, 0, 1);
	m_source_surface->SetDomain(1, 0, 1);
	is_first_const = whichone;
	par = s;
	is_curve_reverse = false;
	is_vector_reverse = false;
}

VectorField1d::~VectorField1d()
{
	delete m_source_surface;
}

void VectorField1d::ReverseCurve()
{
	is_curve_reverse = !is_curve_reverse;
}

void VectorField1d::ReverseVector()
{
	is_vector_reverse = !is_vector_reverse;
}

void VectorField1d::GetAll(double t, ON_3dPoint &p, ON_3dVector &n, ON_3dVector &v)
{
	int n_ceo = 1;
	int v_ceo = 1;
	if (is_curve_reverse)
	{
		t = 1 - t;
		v_ceo = -1;
	}
	if (is_vector_reverse)
	{
		n_ceo = -1;
	}
	ON_3dPoint pp;
	ON_3dVector du;
	ON_3dVector dv;
	if (is_first_const)
	{
		m_source_surface->Ev1Der(par, t, pp, du, dv);
		p = pp;
		n = (du * n_ceo);
		n.Unitize();
		v = (dv * v_ceo);
		v.Unitize();
		return;
	}
	else
	{
		m_source_surface->Ev1Der(t, par, pp, du, dv);
		p = pp;
		n = (dv * n_ceo);
		n.Unitize();
		v = (du * v_ceo);
		v.Unitize();
		return;
	}
}

ON_3dVector VectorField1d::MixedDer(double t)
{
	int n_ceo = 1;
	int v_ceo = 1;
	if (is_curve_reverse)
	{
		t = 1 - t;
		v_ceo = -1;
	}
	if (is_vector_reverse)
	{
		n_ceo = -1;
	}
	ON_3dVector result;
	ON_3dVector ru;
	ON_3dPoint rub;
	if (is_first_const)
	{
		m_source_surface->Ev2Der(par, t, rub, ru, ru, ru, result, ru);
		return result * v_ceo * n_ceo;
	}
	m_source_surface->Ev2Der(t, par, rub, ru, ru, ru, result, ru);
	return result * v_ceo * n_ceo;
}

void FilletSurface::SetRailUniformKnot(int num_knot, double *knot_on_rail0, double *knot_on_rail1)
{
	if (num_knot < 3)
	{
		return;
	}
	if (knot_on_rail0 == nullptr || knot_on_rail1 == nullptr)
	{
		return;
	}
	for (int i = 0; i < num_knot; i++)
	{
		*(knot_on_rail0 + i) = i * 1.0 / (num_knot - 1);
		*(knot_on_rail1 + i) = i * 1.0 / (num_knot - 1);
	}
}

void ThreeDegreeBsplineInterplate_Period(ON_NurbsCurve &onc, std::vector<ON_3dPoint> Q, double *knot)
{
	int K = Q.size() - 1;
	int L = K + 2;
	ON_NurbsCurve solve_N;
	solve_N.Create(1, false, 4, L + 1);
	solve_N.SetKnot(0, knot[0]);
	solve_N.SetKnot(1, knot[0]);
	for (int i = 0; i <= K; i++)
	{
		solve_N.SetKnot(i + 2, knot[i]);
	}
	solve_N.SetKnot(K + 3, knot[K]);
	solve_N.SetKnot(K + 4, knot[K]);
	Eigen::MatrixXd N = Eigen::MatrixXd::Zero(K + 1, K + 1);
	Eigen::MatrixXd QQ = Eigen::MatrixXd::Zero(K + 1, 3);
	ON_3dPoint One = ON_3dPoint(1, 1, 1);
	ON_3dPoint Zero = ON_3dPoint(0, 0, 0);
	for (int j = 1; j <= K - 1; j++)
	{
		solve_N.ZeroCVs();
		solve_N.SetCV(j, One);
		N(j, j - 1) = solve_N.PointAt(knot[j]).x;
		solve_N.SetCV(j, Zero);
		solve_N.SetCV(j + 1, One);
		N(j, j) = solve_N.PointAt(knot[j]).x;
		solve_N.SetCV(j + 1, Zero);
		solve_N.SetCV(j + 2, One);
		N(j, j + 1) = solve_N.PointAt(knot[j]).x;
		QQ(j, 0) = Q[j].x;
		QQ(j, 1) = Q[j].y;
		QQ(j, 2) = Q[j].z;
	}
	N(0, 0) = knot[K] - knot[K - 1];
	N(0, K) = knot[1] - knot[0];
	ON_3dPoint Q0 = Q[K] * (knot[1] - knot[0]) + Q[0] * (knot[K] - knot[K - 1]);
	QQ(0, 0) = Q0.x;
	QQ(0, 1) = Q0.y;
	QQ(0, 2) = Q0.z;
	N(K, 0) = (knot[K] - knot[K - 1]) * (-1 / (knot[1] - knot[0]) - 1 / (knot[2] - knot[0]));
	N(K, 1) = (knot[K] - knot[K - 1]) / (knot[2] - knot[0]);
	N(K, K - 1) = -(knot[1] - knot[0]) / (knot[K] - knot[K - 2]);
	N(K, K) = (knot[1] - knot[0]) * (1 / (knot[K] - knot[K - 1]) + 1 / (knot[K] - knot[K - 2]));
	ON_3dPoint QK = Q[K] * (knot[1] - knot[0]) / (knot[K] - knot[K - 1]) - Q[0] * (knot[K] - knot[K - 1]) / (knot[1] - knot[0]);
	QQ(K, 0) = QK.x;
	QQ(K, 1) = QK.y;
	QQ(K, 2) = QK.z;
	Eigen::MatrixXd P = Eigen::MatrixXd::Zero(K + 1, 3);
	P = N.partialPivLu().solve(QQ);

	onc.Create(3, false, 4, L + 1);
	for (int i = 0; i < L + 3; i++)
	{
		onc.SetKnot(i, solve_N.Knot(i));
	}
	onc.SetCV(0, Q[0]);
	for (int i = 1; i < L; i++)
	{
		onc.SetCV(i, ON_3dPoint(P(i - 1, 0), P(i - 1, 1), P(i - 1, 2)));
	}
	onc.SetCV(L, Q[K]);
}

void FilletSurface::Skinning(ON_NurbsSurface &ons, int num_fill_curve, double *u_knot)
{
	int K = num_fill_curve - 1;
	int m = K + 2;
	int p = 4;
	std::vector<ON_NurbsCurve> fill_nurbs_curve;
	std::vector<ON_NurbsCurve> InterpolateCurve;
	for (auto it : FillCurveArray)
	{
		ON_NurbsCurve onc;
		it.GetNurbForm(onc);
		fill_nurbs_curve.push_back(onc);
	}
	int n = fill_nurbs_curve[0].CVCount() - 1;
	int q = fill_nurbs_curve[0].Order();
	for (int i = 0; i <= n; i++)
	{
		std::vector<ON_3dPoint> InterpolatePoints;
		for (auto iter : fill_nurbs_curve)
		{
			ON_3dPoint p;
			iter.GetCV(i, p);
			InterpolatePoints.push_back(p);
		}
		ON_NurbsCurve onc;
		ThreeDegreeBsplineInterplate_Period(onc, InterpolatePoints, u_knot);
		InterpolateCurve.push_back(onc);
	}

	ons.Create(3, false, p, q, m + 1, n + 1);
	for (int i = 0; i < p + m - 1; i++)
	{
		ons.SetKnot(0, i, InterpolateCurve[0].Knot(i));
	}
	for (int i = 0; i < q + n - 1; i++)
	{
		ons.SetKnot(1, i, fill_nurbs_curve[0].Knot(i));
	}
	for (int i = 0; i <= m; i++)
	{
		for (int j = 0; j <= n; j++)
		{
			ON_3dPoint p;
			InterpolateCurve[j].GetCV(i, p);
			ons.SetCV(i, j, p);
		}
	}
}