#include "ChiralityMathTools.h"
#include "thirdparty/eigen/Eigen/Dense"
#include <assert.h>

double ChiralityMath::Bisection(const std::function<double(double)> &F, double L, double R, double eps)
{
	assert(L <= R);
	double mid = (L + R) / 2;
	double lv = F(L);
	double rv = F(R);
	double mv = F(mid);
	assert(lv * rv <= 0);
	int max_iter = 1000;
	while (abs(mv) >= eps && max_iter > 0)
	{
		if (lv * mv >= 0 && mv * rv <= 0)
		{
			L = mid;
		}
		else if (rv * mv >= 0 && mv * lv <= 0)
		{
			R = mid;
		}
		else
		{
			assert(0);
		}
		mid = (L + R) / 2;
		lv = F(L);
		rv = F(R);
		mv = F(mid);
		max_iter--;
	}
	return mid;
}

double ChiralityMath::ArcLength(const ON_BezierCurve &obc, double from, double to)
{
	const int N = 100;
	double sum = 0.0;
	for (int i = 1; i <= N; ++i)
	{
		sum += obc.DerivativeAt(from + (to - from) / double(N) * (double(i) - 0.5)).Length();
	}
	return sum * (to - from) / double(N);
}

double ChiralityMath::ArcLength(const ON_NurbsCurve &onc, double from, double to)
{
	const int N = 100;
	double sum = 0.0;
	for (int i = 1; i <= N; ++i)
	{
		sum += onc.DerivativeAt(from + (to - from) / double(N) * (double(i) - 0.5)).Length();
	}
	return sum * (to - from) / double(N);
}

double ChiralityMath::ArcLength(const ON_BezierCurve &obc)
{
	return ArcLength(obc, 0, 1);
}
double ChiralityMath::ArcLength(const ON_NurbsCurve &onc)
{
	double t0, t1;
	onc.GetDomain(&t0, &t1);
	return ArcLength(onc, t0, t1);
}

std::vector<double> ChiralityMath::GenerateUniformArcLength(const ON_NurbsCurve &onc, int num_param)
{
	double total_length = ArcLength(onc);
	std::vector<double> length_param(num_param, 0.0);
	for (int i = 1; i < num_param; ++i)
	{
		length_param[i] = length_param[i - 1] + total_length / (num_param - 1);
	}
	std::vector<double> param(num_param, 0.0);
	onc.GetDomain(&param[0], &param[num_param - 1]);
	for (int i = 1; i <= num_param - 2; ++i)
	{
		auto lambda = [=](double t) -> double
		{
			return ArcLength(onc, param[0], t) - length_param[i];
		};
		param[i] = Bisection(lambda, param[i - 1], param.back());
	}
	return param;
}

double ChiralityMath::Bernstein(int n, int i, double t)
{
	if (n < 1 || i < 0 || i > n)
	{
		return 0.0;
	}
	if (n == 1)
	{
		return i == 0 ? 1.0 - t : t;
	}
	ON_BezierCurve obc(2, false, n + 1);
	for (int j = 0; j < obc.CVCount(); ++j)
	{
		obc.SetCV(j, ON_3dPoint::Origin);
	}
	obc.SetCV(i, ON_3dPoint(1.0, 1.0, 1.0));
	return obc.PointAt(t).x;
}

double ChiralityMath::Torsion(const ON_NurbsCurve &onc, double t)
{
	if (onc.Degree() < 3 || onc.Dimension() < 3 || onc.IsRational())
	{
		return 0.0;
	}
	ON_NurbsCurve der_onc(3, false, onc.Order() - 1, onc.CVCount() - 1);
	int k = onc.Order();
	int n = onc.CVCount() - 1;
	for (int i = 1; i <= n; ++i)
	{
		ON_3dPoint p1, p2;
		onc.GetCV(i, p1);
		onc.GetCV(i - 1, p2);
		double t1 = onc.Knot(i + k - 1);
		double t2 = onc.Knot(i);
		der_onc.SetCV(i - 1, (k - 1) / (t1 - t2) * ON_3dPoint(p1 - p2));
	}
	for (int i = 0; i < der_onc.KnotCount(); ++i)
	{
		der_onc.SetKnot(i, onc.Knot(i + 1));
	}
	ON_3dPoint first_der_p;
	ON_3dVector second_der;
	ON_3dVector third_der;
	der_onc.Ev2Der(t, first_der_p, second_der, third_der);
	ON_3dVector first_der = ON_3dVector(first_der_p);
	ON_3dVector cross_vector = ON_3dVector::CrossProduct(first_der, second_der);
	return ON_3dVector::DotProduct(cross_vector, third_der) / cross_vector.LengthSquared();
}

ON_NurbsCurve ChiralityMath::UniformG1(ON_3dPoint ps, ON_3dPoint pe, ON_3dVector vs, ON_3dVector ve)
{
	double L = ps.DistanceTo(pe);
	vs.Unitize();
	ve.Unitize();
	Eigen::Matrix<double, 4, 4> A;
	A(0, 0) = 1.0;
	A(0, 1) = 4.0;
	A(0, 2) = 1.0;
	A(0, 3) = 0.0;
	A(1, 0) = 0.0;
	A(1, 1) = 1.0;
	A(1, 2) = 4.0;
	A(1, 3) = 1.0;
	A(2, 0) = -1.0;
	A(2, 1) = 0.0;
	A(2, 2) = 1.0;
	A(2, 3) = 0.0;
	A(3, 0) = 0.0;
	A(3, 1) = -1.0;
	A(3, 2) = 0.0;
	A(3, 3) = 1.0;
	Eigen::Matrix<double, 4, 3> B;
	B(0, 0) = 6.0 * ps.x;
	B(0, 1) = 6.0 * ps.y;
	B(0, 2) = 6.0 * ps.z;
	B(1, 0) = 6.0 * pe.x;
	B(1, 1) = 6.0 * pe.y;
	B(1, 2) = 6.0 * pe.z;
	B(2, 0) = (2.0 * L / 3.0) * vs.x;
	B(2, 1) = (2.0 * L / 3.0) * vs.y;
	B(2, 2) = (2.0 * L / 3.0) * vs.z;
	B(3, 0) = (2.0 * L / 3.0) * ve.x;
	B(3, 1) = (2.0 * L / 3.0) * ve.y;
	B(3, 2) = (2.0 * L / 3.0) * ve.z;
	Eigen::Matrix<double, 4, 3> P;
	P = A.partialPivLu().solve(B);
	ON_NurbsCurve onc;
	onc.Create(3, false, 4, 4);
	for (int i = 0; i < 6; ++i)
	{
		onc.SetKnot(i, i + 1);
	}
	for (int i = 0; i < 4; ++i)
	{
		onc.SetCV(i, ON_3dPoint(P(i, 0), P(i, 1), P(i, 2)));
	}
	return onc;
}

void ChiralityMath::Elevate(ON_NurbsCurve &onc)
{
	int v_num = onc.CVCount();
	ON_3dPoint p1, p2;
	std::vector<ON_3dPoint> vp;
	onc.GetCV(0, p1);
	vp.push_back(p1);
	for (int i = 1; i < v_num; i++)
	{
		onc.GetCV(i - 1, p1);
		onc.GetCV(i, p2);
		vp.push_back(p1 * (double(i) / double(v_num)) + p2 * (1 - double(i) / double(v_num)));
	}
	onc.GetCV(v_num - 1, p1);
	vp.push_back(p1);
	double u0, u1;
	onc.GetDomain(&u0, &u1);
	onc.EvPoint(u0, p1);
	onc.EvPoint(u1, p2);
	ON_3dVector v1 = onc.TangentAt(u0);
	ON_3dVector v2 = onc.TangentAt(u1);
	onc.Create(3, false, onc.Order(), v_num + 1);
	int i = 0;
	for (const auto &it : vp)
	{
		onc.SetCV(i, it);
		i++;
	}
}

ON_NurbsCurve ChiralityMath::CubicBsplineInterpolate_G1(const std::vector<ON_3dPoint> &Q, const std::vector<double> &knot, ON_3dVector v0, ON_3dVector vn)
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
	N(0, 0) = 3 / (knot[1] - knot[0]);
	N(K, K) = 3 / (knot[K] - knot[K - 1]);
	ON_3dPoint Q0 = 3 * Q[0] / (knot[1] - knot[0]) + v0;
	QQ(0, 0) = Q0.x;
	QQ(0, 1) = Q0.y;
	QQ(0, 2) = Q0.z;
	ON_3dPoint QK = 3 * Q[K] / (knot[K] - knot[K - 1]) - vn;
	QQ(K, 0) = QK.x;
	QQ(K, 1) = QK.y;
	QQ(K, 2) = QK.z;

	Eigen::MatrixXd P = Eigen::MatrixXd::Zero(K + 1, 3);
	P = N.partialPivLu().solve(QQ);
	ON_NurbsCurve onc;
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
	return onc;
}

ON_NurbsSurface ChiralityMath::Skinning(const std::vector<ON_NurbsCurve> &curve_list, const std::vector<double> &u_knots, const std::vector<std::pair<ON_3dVector, ON_3dVector>> &pair_tangent)
{
	int K = curve_list.size() - 1;
	int n = curve_list[0].CVCount() - 1;
	std::vector<std::vector<ON_3dPoint>> Q(n + 1, std::vector<ON_3dPoint>(K + 1, ON_3dPoint::Origin));
	for (int k = 0; k <= K; ++k)
	{
		for (int j = 0; j <= n; ++j)
		{
			curve_list[k].GetCV(j, Q[j][k]);
		}
	}
	std::vector<ON_NurbsCurve> new_curve_list;
	for (int i = 0; i <= n; ++i)
	{
		new_curve_list.push_back(ChiralityMath::CubicBsplineInterpolate_G1(Q[i], u_knots, pair_tangent[i].first, pair_tangent[i].second));
	}
	ON_NurbsSurface ons;
	ON_3dPoint p;
	ons.Create(3, false, 4, 4, new_curve_list[0].CVCount(), curve_list[0].CVCount());
	for (int i = 0; i < ons.CVCount(0); ++i)
	{
		for (int j = 0; j < ons.CVCount(1); ++j)
		{
			new_curve_list[j].GetCV(i, p);
			ons.SetCV(i, j, p);
		}
	}
	for (int i = 0; i < ons.KnotCount(0); ++i)
	{
		ons.SetKnot(0, i, new_curve_list[0].Knot(i));
	}
	for (int j = 0; j < ons.KnotCount(1); ++j)
	{
		ons.SetKnot(1, j, curve_list[0].Knot(j));
	}
	return ons;
}

ON_NurbsSurface ChiralityMath::GenerateCylinder(const ON_NurbsCurve &parent_curve, ON_3dVector dir, double t0, double t1)
{
	dir.Unitize();
	ON_NurbsSurface ons;
	ons.Create(3, parent_curve.IsRational(), parent_curve.Order(), 2, parent_curve.CVCount(), 2);
	ON_3dPoint p;
	for (int i = 0; i < parent_curve.CVCount(); ++i)
	{
		parent_curve.GetCV(i, p);
		ons.SetCV(i, 0, p + dir * t0);
		ons.SetCV(i, 1, p + dir * t1);
		if (parent_curve.IsRational())
		{
			ons.SetWeight(i, 0, parent_curve.Weight(i));
			ons.SetWeight(i, 1, parent_curve.Weight(i));
		}
	}
	for (int i = 0; i < parent_curve.KnotCount(); ++i)
	{
		ons.SetKnot(0, i, parent_curve.Knot(i));
	}
	ons.SetKnot(1, 0, 0);
	ons.SetKnot(1, 1, 1);
	return ons;
}

ON_NurbsCurve ChiralityMath::ChangeDimensionFrom2To3(const ON_NurbsCurve &onc_2d)
{
	ON_NurbsCurve onc_3d(3, onc_2d.IsRational(), onc_2d.Order(), onc_2d.CVCount());
	for (int i = 0; i < onc_2d.KnotCount(); ++i)
	{
		onc_3d.SetKnot(i, onc_2d.Knot(i));
	}
	ON_3dPoint p;
	for (int i = 0; i < onc_2d.CVCount(); ++i)
	{
		onc_2d.GetCV(i, p);
		onc_3d.SetCV(i, p);
	}
	if (onc_2d.IsRational())
	{
		for (int i = 0; i < onc_2d.CVCount(); ++i)
		{
			onc_3d.SetWeight(i, onc_2d.Weight(i));
		}
	}
	return onc_3d;
}
