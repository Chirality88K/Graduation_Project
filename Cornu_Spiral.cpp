#include "Cornu_Spiral.h"
#include <assert.h>
#include "Spiral.h"
#include "write3dm.h"
#include "ChiralityMathTools.h"

Cornu_Spiral::Cornu_Spiral(ON_2dPoint p1, ON_2dPoint p2, ON_2dVector t1, ON_2dVector t2)
{
	assert(Check_And_Init(p1, p2, t1, t2));
	ON_2dVector D = p2 - p1;
	double length = D.Length();
	t1.Unitize();
	t2.Unitize();
	double phi1 = ChiralityMath::ComputeSignedAngle(t1, D);
	double phi2 = ChiralityMath::ComputeSignedAngle(D, t2);
	assert((0 < phi1 && phi1 < phi2&& phi2 < PI) ||
		(-PI < phi1&& phi1 < 0 && phi2 > 0 && phi1 + phi2 > 0));
	double h_phi1_phi2 = Compute_S(phi1 + phi2) * cos(phi1) - Compute_C(phi1 + phi2) * sin(phi1);
	if (0 < phi1 && phi1 < phi2 && h_phi1_phi2 <= 0)
	{
		double theta = MidSection_f(phi1, phi2);
		m_a = length / ((Compute_S(theta + phi1 + phi2) - Compute_S(theta)) * sin(theta + phi1) +
						(Compute_C(theta + phi1 + phi2) - Compute_C(theta)) * cos(theta + phi1));
		m_T0 = t1;
		m_T0.Rotate(-theta);
		m_N0 = m_T0;
		m_N0.Rotate(1, 0);
		m_P0 = p1 - m_a * (Compute_C(theta) * m_T0 + Compute_S(theta) * m_N0);
		m_theta0 = theta;
		m_theta1 = theta + phi1 + phi2;
	}
	if ((0 < phi1 && phi1 < phi2 && h_phi1_phi2 > 0) || (phi1 <= 0))
	{
		double omega = MidSection_g(phi1, phi2);
		m_a = length / ((Compute_S(omega + phi1 + phi2) + Compute_S(omega)) * sin(omega + phi1) +
						(Compute_C(omega + phi1 + phi2) + Compute_C(omega)) * cos(omega + phi1));
		m_T0 = t1;
		m_T0.Rotate(-omega);
		m_N0 = m_T0;
		m_N0.Rotate(1, 0);
		m_P0 = p1 + m_a * (Compute_C(omega) * m_T0 + Compute_S(omega) * m_N0);
		m_theta0 = -omega;
		m_theta1 = omega + phi1 + phi2;
	}
}

ParameterCurve Cornu_Spiral::GetParameterCurve() const
{
	ON_3dPoint P = m_P0;
	ON_3dVector T = m_T0;
	ON_3dVector N = m_N0;
	double a = m_a;
	bool is_reverse = m_is_reverse;
	bool is_mirror = m_is_mirror;
	ON_Xform form;
	ON_3dPoint ps = m_PS;
	ON_3dPoint pe = m_PE;
	ON_3dVector normal = ON_3dVector::CrossProduct(ON_3dVector::ZAxis, pe - ps);
	form.Mirror(ps, normal);
	double startLength = GetSignedArcLength(m_theta0);
	double totalLength = GetSignedArcLength(m_theta1) - startLength;
	auto pos = [=](double s)->FrenetFrame {
		int inverse = 1;
		if (is_reverse)
		{
			s = totalLength - s;
			inverse = -1;
		}
		double t;
		if (s + startLength >= 0)
		{
			t = pow(s + startLength, 2) * PI / 2 / a / a;
		}
		else
		{
			t = -pow(s + startLength, 2) * PI / 2 / a / a;
		}
		double C = Cornu_Spiral::Compute_C(abs(t));
		double S = Cornu_Spiral::Compute_S(abs(t));
		ON_3dPoint pp;
		ON_3dVector alpha, beta;
		if (t >= 0)
		{
			pp = P + a * (C * T + S * N);
			alpha = cos(t) * T + sin(t) * N;
			beta = alpha;
			beta.Rotate(1, 0, ON_3dVector::ZAxis);
		}
		else
		{
			pp = P - a * (C * T + S * N);
			alpha = cos(t) * T - sin(t) * N;
			beta = alpha;
			beta.Rotate(1, 0, ON_3dVector::ZAxis);
		}
		if (is_mirror)
		{
			pp.Transform(form);
			alpha.Transform(form);
			beta.Transform(form);
		}
		return FrenetFrame(pp, alpha * inverse, beta * inverse);
	};

	auto der = [=](double s)->ON_3dVector {
		int inverse = 1;
		if (is_reverse)
		{
			s = totalLength - s;
			inverse = -1;
		}
		double t;
		if (s + startLength >= 0)
		{
			t = pow(s + startLength, 2) * PI / 2 / a / a;
		}
		else
		{
			t = -pow(s + startLength, 2) * PI / 2 / a / a;
		}
		ON_3dVector alpha;
		if (t >= 0)
		{
			alpha = cos(t) * T + sin(t) * N;
		}
		else
		{
			alpha = cos(t) * T - sin(t) * N;
		}
		alpha *= a / sqrt(2 * PI * abs(t));
		if (is_mirror)
		{
			alpha.Transform(form);
		}
		return alpha * inverse;
	};
	double range[2] = { 0,totalLength };
	ParameterCurve pc(pos, range);
	pc.SetDerivative(der);
	return pc;
}

ParameterCurve Cornu_Spiral::GetParamCornuSpiral(ON_2dPoint ps, ON_2dPoint pe, ON_2dVector vs, ON_2dVector ve)
{
	Cornu_Spiral c_s(ps, pe, vs, ve);
	return c_s.GetParameterCurve();
}

ON_NurbsCurve Cornu_Spiral::GetNurbs(ON_2dPoint ps, ON_2dPoint pe, ON_2dVector vs, ON_2dVector ve)
{
	Cornu_Spiral c_s(ps, pe, vs, ve);
	int n = 1000;
	double theta = 0;
	double s0 = c_s.GetSignedArcLength(c_s.m_theta0);
	std::vector<ON_3dPoint> pv;
	std::vector<double> knot;
	for (int i = 0; i <= 20; i++)
	{
		theta = (c_s.m_theta1 - c_s.m_theta0) / n * i * 50 + c_s.m_theta0;
		ON_2dPoint p = c_s.GetValue(theta);
		pv.push_back(p);
		knot.push_back(c_s.GetSignedArcLength(theta) - s0);
	}
	ON_NurbsCurve onc;
	ON_2dVector T1 = c_s.GetTangent(c_s.m_theta0);
	ON_2dVector T2 = c_s.GetTangent(c_s.m_theta1);
	if (c_s.m_is_mirror)
	{
		ON_Xform mirror;
		mirror.Mirror(c_s.m_PS, ON_3dVector::CrossProduct(ON_3dVector::ZAxis, c_s.m_PE - c_s.m_PS));
		for (ON_3dPoint& p : pv)
		{
			p.Transform(mirror);
		}
		T1.Transform(mirror);
		T2.Transform(mirror);
	}
	ThreeDegreeBsplineInterplate_Tan(onc, pv, knot, T1, T2);
	return onc;
}

double Cornu_Spiral::Compute_C(double theta)
{
	int n = 1000;
	double sum = 0;
	for (int k = 1; k <= n; k++)
	{
		sum += cos(theta * (k * k - k + 0.25) / n / n);
	}
	return sum * sqrt(2 * theta / PI) / n;
}

double Cornu_Spiral::Compute_S(double theta)
{
	int n = 1000;
	double sum = 0;
	for (int k = 1; k <= n; k++)
	{
		sum += sin(theta * (k * k - k + 0.25) / n / n);
	}
	return sum * sqrt(2 * theta / PI) / n;
}

double Cornu_Spiral::Compute_f(double theta, double phi1, double phi2)
{
	return (Compute_S(theta + phi1 + phi2) - Compute_S(theta)) * cos(theta + phi1) -
		   (Compute_C(theta + phi1 + phi2) - Compute_C(theta)) * sin(theta + phi1);
}

double Cornu_Spiral::Compute_g(double omega, double phi1, double phi2)
{
	return (Compute_S(omega + phi1 + phi2) + Compute_S(omega)) * cos(omega + phi1) -
		   (Compute_C(omega + phi1 + phi2) + Compute_C(omega)) * sin(omega + phi1);
}

double Cornu_Spiral::MidSection_f(double phi1, double phi2)
{
	double lambda = (1 - cos(phi1)) / (1 - cos(phi2));
	double theta0 = (phi1 + phi2) * lambda * lambda / (1 - lambda * lambda);
	double left = 0;
	double right = theta0;
	double leftv = Compute_f(0, phi1, phi2);
	double rightv = Compute_f(theta0, phi1, phi2);
	double mid;
	double midv;
	assert(leftv * rightv <= 0);
	while (abs(rightv) > 1e-8 && abs(leftv) > 1e-8)
	{
		mid = (left + right) / 2;
		midv = Compute_f(mid, phi1, phi2);
		if (midv * leftv > 0)
		{
			left = mid;
			leftv = midv;
			continue;
		}
		right = mid;
		rightv = midv;
	}
	if (abs(rightv) < abs(leftv))
	{
		return right;
	}
	else
	{
		return left;
	}
}

double Cornu_Spiral::MidSection_g(double phi1, double phi2)
{
	double left = 0;
	const double pi = acos(-1.0);
	double right = pi / 2 - phi1;
	double leftv = Compute_g(left, phi1, phi2);
	double rightv = Compute_g(right, phi1, phi2);
	double mid;
	double midv;
	assert(leftv * rightv <= 0);
	while (abs(rightv) > 1e-8 && abs(leftv) > 1e-8)
	{
		mid = (left + right) / 2;
		midv = Compute_g(mid, phi1, phi2);
		if (midv * leftv > 0)
		{
			left = mid;
			leftv = midv;
			continue;
		}
		right = mid;
		rightv = midv;
	}
	if (abs(rightv) < abs(leftv))
	{
		return right;
	}
	else
	{
		return left;
	}
}

bool Cornu_Spiral::Check_And_Init(ON_2dPoint& ps, ON_2dPoint& pe, ON_2dVector& vs, ON_2dVector& ve)
{
	m_PS = ps;
	m_PE = pe;
	ON_2dVector D = pe - ps;
	if (!D.Unitize() || !vs.Unitize() || !ve.Unitize())
	{
		return false;
	}
	double phi1 = ChiralityMath::ComputeSignedAngle(vs, D);
	double phi2 = ChiralityMath::ComputeSignedAngle(D, ve);
	constexpr double angle_error = 1e-3;
	if (abs(phi1) < angle_error || abs(phi2) < angle_error || 
		abs(phi1 + phi2) < angle_error || abs(phi1 - phi2) < angle_error)
	{
		return false;
	}
	if (abs(phi1 - PI) < angle_error || abs(phi1 + PI) < angle_error ||
		abs(phi2 - PI) < angle_error || abs(phi2 + PI) < angle_error)
	{
		return false;
	}
	if (phi1 > 0 && phi2 > 0)
	{
		if (phi1 < phi2)
		{
			return true;
		}
		std::swap(ps, pe);
		ON_2dVector new_vs = ve - 2 * ON_3dVector::DotProduct(ve, D) * D;
		ON_2dVector new_ve = vs - 2 * ON_3dVector::DotProduct(vs, D) * D;
		vs = new_vs;
		ve = new_ve;
		m_is_mirror = true;
		m_is_reverse = true;
		return true;
	}
	else if (phi1 < 0 && phi2 < 0)
	{
		if (phi1 < phi2)
		{
			std::swap(ps, pe);
			ON_2dVector new_vs = -ve;
			ON_2dVector new_ve = -vs;
			vs = new_vs;
			ve = new_ve;
			m_is_reverse = true;
			return true;
		}
		vs = -vs + 2 * ON_3dVector::DotProduct(vs, D) * D;
		ve = -ve + 2 * ON_3dVector::DotProduct(ve, D) * D;
		m_is_mirror = true;
		return true;
	}
	else if (phi1 < 0 && phi2 > 0)
	{
		if (phi1 + phi2 > 0)
		{
			return true;
		}
		std::swap(ps, pe);
		ON_2dVector new_vs = -ve;
		ON_2dVector new_ve = -vs;
		vs = new_vs;
		ve = new_ve;
		m_is_reverse = true;
		return true;
	}
	else if (phi1 > 0 && phi2 < 0)
	{
		if (phi1 + phi2 > 0)
		{
			std::swap(ps, pe);
			ON_2dVector new_vs = ve - 2 * ON_3dVector::DotProduct(ve, D) * D;
			ON_2dVector new_ve = vs - 2 * ON_3dVector::DotProduct(vs, D) * D;
			vs = new_vs;
			ve = new_ve;
			m_is_mirror = true;
			m_is_reverse = true;
			return true;
		}
		vs = -vs + 2 * ON_3dVector::DotProduct(vs, D) * D;
		ve = -ve + 2 * ON_3dVector::DotProduct(ve, D) * D;
		m_is_mirror = true;
		return true;
	}
	return false;
}

void Cornu_Spiral::Add_to_Model(ONX_Model *model, const wchar_t *name, ON_Color color)
{
	ON_3dPointArray list;
	int n = 500;
	double theta = 0;
	for (int i = 0; i <= n; i++)
	{
		theta = (m_theta1 - m_theta0) / n * i + m_theta0;
		ON_2dPoint p = GetValue(theta);
		list.Append(p);
	}
	ON_Polyline pl = list;
	ON_PolylineCurve *pc = new ON_PolylineCurve(pl);
	const int layer_index = model->AddLayer(name, color);
	ON_3dmObjectAttributes *attributes = new ON_3dmObjectAttributes();
	attributes->m_layer_index = layer_index;
	attributes->m_name = name;
	model->AddManagedModelGeometryComponent(pc, attributes);
}

void Cornu_Spiral::Add_Nurbs_to_Model(ONX_Model *model, const wchar_t *name, ON_Color color)
{
	int n = 1000;
	double theta = 0;
	double s0 = GetSignedArcLength(m_theta0);
	std::vector<ON_3dPoint> pv;
	std::vector<double> knot;
	for (int i = 0; i <= 20; i++)
	{
		theta = (m_theta1 - m_theta0) / n * i * 50 + m_theta0;
		ON_2dPoint p = GetValue(theta);
		pv.push_back(p);
		knot.push_back(GetSignedArcLength(theta) - s0);
	}
	ON_NurbsCurve *onc = new ON_NurbsCurve();
	ON_2dVector T1 = GetTangent(m_theta0);
	ON_2dVector T2 = GetTangent(m_theta1);
	ThreeDegreeBsplineInterplate_Tan(*onc, pv, knot, T1, T2);
	const int layer_index = model->AddLayer(name, color);
	ON_3dmObjectAttributes *attributes = new ON_3dmObjectAttributes();
	attributes->m_layer_index = layer_index;
	attributes->m_name = name;
	model->AddManagedModelGeometryComponent(onc, attributes);
}

ON_2dPoint Cornu_Spiral::GetValue(double theta)
{
	if (theta >= 0)
	{
		return m_P0 + m_a * (m_T0 * Compute_C(theta) + m_N0 * Compute_S(theta));
	}
	return m_P0 - m_a * (m_T0 * Compute_C(-theta) + m_N0 * Compute_S(-theta));
}

ON_2dVector Cornu_Spiral::GetTangent(double theta)
{
	if (theta >= 0)
	{
		return m_T0 * cos(theta) + m_N0 * sin(theta);
	}
	return m_T0 * cos(theta) - m_N0 * sin(theta);
}

double Cornu_Spiral::GetSignedArcLength(double theta) const
{
	if (theta >= 0)
	{
		return m_a * sqrt(2 * theta / PI);
	}
	return -m_a * sqrt(-2 * theta / PI);
}

double Cornu_Spiral::GetSignedCurvature(double theta) const
{
	if (theta >= 0)
	{
		return sqrt(2 * theta * PI) / m_a;
	}
	return -sqrt(-2 * theta * PI) / m_a;
}

void Cornu_Spiral::Raise_to_3D(double zheight, ONX_Model *model, const wchar_t *name, ON_Color color)
{
	double s0 = GetSignedArcLength(m_theta0);
	double L = GetSignedArcLength(m_theta1) - s0;
	int n = 1000;
	double theta = 0;
	double s = 0;
	std::vector<ON_3dPoint> pv;
	std::vector<double> knot;
	for (int i = 0; i <= 20; i++)
	{
		theta = (m_theta1 - m_theta0) / n * i * 50 + m_theta0;
		ON_3dPoint p = GetValue(theta);
		s = GetSignedArcLength(theta) - s0;
		p.z = zheight * (-2 * pow(s / L, 3) + 3 * pow(s / L, 2));
		pv.push_back(p);
		knot.push_back(s);
	}
	ON_NurbsCurve *onc = new ON_NurbsCurve();
	ON_2dVector T1 = GetTangent(m_theta0);
	ON_2dVector T2 = GetTangent(m_theta1);
	ThreeDegreeBsplineInterplate_Tan(*onc, pv, knot, T1, T2);
	const int layer_index = model->AddLayer(name, color);
	ON_3dmObjectAttributes *attributes = new ON_3dmObjectAttributes();
	attributes->m_layer_index = layer_index;
	attributes->m_name = name;
	model->AddManagedModelGeometryComponent(onc, attributes);
	PrintCurvature(*onc, "rise " + std::to_string(zheight) + " curvature");
}

void Cornu_Spiral::Cornu_test(ONX_Model *model)
{
	Cornu_Spiral *cs = new Cornu_Spiral(ON_2dPoint(0, 0), ON_2dPoint(10, 0),
										ON_2dVector(cos(10.0 / 180.0 * PI), -sin(10.0 / 180.0 * PI)),
										ON_2dVector(cos(15.0 / 180.0 * PI), sin(15.0 / 180.0 * PI)));
	// cs->Add_to_Model(model, L"test_Cornu_spiral", ON_Color::SaturatedGreen);
	cs->Add_Nurbs_to_Model(model, L"C_Shape", ON_Color::SaturatedMagenta);
	cs->Raise_to_3D(2, model, L"C_Shape_Rise2");
	cs->Raise_to_3D(-6, model, L"C_Shape_Rise-6");
	cs->Raise_to_3D(10, model, L"C_Shape_Rise10");
	delete cs;
	cs = new Cornu_Spiral(ON_2dPoint(0, 0), ON_2dPoint(10, 0),
						  ON_2dVector(-cos(75.0 / 180.0 * PI), sin(75.0 / 180.0 * PI)),
						  ON_2dVector(-cos(15.0 / 180.0 * PI), sin(15.0 / 180.0 * PI)));
	cs->Add_Nurbs_to_Model(model, L"S_Shape", ON_Color::SaturatedGold);
	cs->Raise_to_3D(2, model, L"S_Shape_Rise2");
	cs->Raise_to_3D(-6, model, L"S_Shape_Rise-6");
	cs->Raise_to_3D(10, model, L"S_Shape_Rise10");
	delete cs;
}
