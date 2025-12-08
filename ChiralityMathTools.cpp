#include "ChiralityMathTools.h"
#include "thirdparty/eigen/Eigen/Dense"
#include "FixAxisBezier3D.h"
#include <assert.h>
#include <random>
#include <map>
#include <array>

double ChiralityMath::Bisection(const std::function<double(double)> &F, double L, double R, double* real_error, double eps)
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
	if (real_error != nullptr)
	{
		(*real_error) = abs(mv);
	}
	return mid;
}

double ChiralityMath::Newton(const std::function<double(double)>& f, const std::function<double(double)>& df, double x0, double eps, int max_iter_num)
{
	double x = x0;
	int cnt = 0;
	eps = abs(eps);
	while (abs(f(x)) > eps && cnt < max_iter_num)
	{
		double dfx = df(x);
		assert(abs(dfx) > eps);
		x = x - f(x) / dfx;
		++cnt;
	}
	if (abs(f(x)) > eps)
	{
		CHIRALITY_WARN(std::string("Newton: Presicion is not enough!!!"));
	}
	return x;
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

double ChiralityMath::ComputeSignedAngle(const ON_3dVector &v1, const ON_3dVector &v2, const ON_3dVector &T)
{
	double l1 = v1.Length();
	double l2 = v2.Length();
	assert(l1 > 1e-6 && l2 > 1e-6);
	double cos_angle = ON_3dVector::DotProduct(v1, v2) / l1 / l2;
	double angle = 0.0;
	if (cos_angle > -1.0 && cos_angle < 1.0)
	{
		angle = acos(cos_angle);
	}
	if (cos_angle <= -1.0)
	{
		angle = PI;
	}
	if (cos_angle >= 1.0)
	{
		angle = 0.0;
	}
	if (ON_3dVector::DotProduct(ON_3dVector::CrossProduct(v1, v2), T) < 0)
	{
		angle = -angle;
	}
	return angle;
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
	double result[12];
	onc.Evaluate(t, 3, 3, result);
	ON_3dVector first_der(result[3], result[4], result[5]);
	ON_3dVector second_der(result[6], result[7], result[8]);
	ON_3dVector third_der(result[9], result[10], result[11]);
	ON_3dVector cross_vector = ON_3dVector::CrossProduct(first_der, second_der);
	if (cross_vector.LengthSquared()<1e-8)
	{
		return 0.0;
	}
	double torsion = ON_3dVector::DotProduct(cross_vector, third_der) / cross_vector.LengthSquared();
	return torsion;
}

double ChiralityMath::DiscreteCurvature(ON_3dPoint p_before, ON_3dPoint p_mid, ON_3dPoint p_after)
{
	double l1 = p_before.DistanceTo(p_mid);
	double l2 = p_mid.DistanceTo(p_after);
	double l3 = p_before.DistanceTo(p_after);
	assert(l1 * l2 * l3 > 1e-6);
	double A = ON_3dVector::CrossProduct(p_mid - p_before, p_after - p_mid).Length();
	return 2 * A / (l1 * l2 * l3);
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

ON_BezierCurve ChiralityMath::BezierG1_xOy(ON_3dPoint ps, ON_3dPoint pe, ON_3dVector vs, ON_3dVector ve)
{
	vs.Unitize();
	ve.Unitize();
	ON_BezierCurve obc(2, false, 4);
	obc.SetCV(0, ps);
	obc.SetCV(3, pe);
	double L = ps.DistanceTo(pe);
	obc.SetCV(1, ps + vs * L / 4);
	obc.SetCV(2, pe - ve * L / 4);
	return obc;
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

void ChiralityMath::Elevate(ON_BezierCurve &obc)
{
	int n = obc.CVCount();
	if (n < 2)
	{
		return;
	}
	std::vector<ON_3dPoint> parr;
	ON_3dPoint p, q;
	obc.GetCV(0, p);
	parr.push_back(p);
	for (int i = 1; i < n; i++)
	{
		obc.GetCV(i - 1, p);
		obc.GetCV(i, q);
		parr.push_back(p * (i * 1.0 / n) + q * (1 - i * 1.0 / n));
	}
	obc.GetCV(n - 1, q);
	parr.push_back(q);
	obc.Create(3, false, n + 1);
	for (int i = 0; i < n + 1; i++)
	{
		obc.SetCV(i, parr[i]);
	}
}

void ChiralityMath::Elevate(std::vector<ON_3dPoint>& vp)
{
	int n = vp.size();
	if (n < 2)
	{
		return;
	}
	std::vector<ON_3dPoint> parr;
	ON_3dPoint p, q;
	parr.push_back(vp[0]);
	for (int i = 1; i < n; i++)
	{
		p = vp[i - 1];
		q = vp[i];
		parr.push_back(p * (i * 1.0 / n) + q * (1 - i * 1.0 / n));
	}
	parr.push_back(q);
	vp = parr;
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

ON_NurbsCurve ChiralityMath::CubicBsplineInterpolate_Period(const std::vector<ON_3dPoint> &Q, const std::vector<double> &knot)
{
	assert(Q.front() == Q.back());
	assert(Q.size() == knot.size());
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
	N(0, 0) = 1 / (knot[1] - knot[0]);
	N(0, K) = 1 / (knot[K] - knot[K - 1]);
	ON_3dPoint Q0 = Q[0] / (knot[1] - knot[0]) + Q[K] / (knot[K] - knot[K - 1]);
	QQ(0, 0) = Q0.x;
	QQ(0, 1) = Q0.y;
	QQ(0, 2) = Q0.z;

	double rk1 = knot[K] - knot[K - 1];
	double rk2 = knot[K - 1] - knot[K - 2];
	double r0 = knot[1] - knot[0];
	double r1 = knot[2] - knot[1];
	N(K, 0) = -rk1 / (r1 + r0) - rk1 / r0;
	N(K, 1) = rk1 / (r1 + r0);
	N(K, K - 1) = -r0 / (rk1 + rk2);
	N(K, K) = r0 / rk1 + r0 / (rk1 + rk2);
	ON_3dPoint QK = r0 / rk1 * Q[K] - rk1 / r0 * Q[0];
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

ON_NurbsCurve ChiralityMath::CubicHomoBsplineInterpolate_G1(const std::vector<ON_4dPoint>& Q, const std::vector<double>& knot, ON_4dPoint v0, ON_4dPoint vn)
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
	Eigen::MatrixXd QQ = Eigen::MatrixXd::Zero(K + 1, 4);
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
		QQ(j, 3) = Q[j].w;
	}
	N(0, 0) = 3 / (knot[1] - knot[0]);
	N(K, K) = 3 / (knot[K] - knot[K - 1]);
	QQ(0, 0) = 3 * Q[0].x / (knot[1] - knot[0]) + v0.x * Q[0].w;
	QQ(0, 1) = 3 * Q[0].y / (knot[1] - knot[0]) + v0.y * Q[0].w;
	QQ(0, 2) = 3 * Q[0].z / (knot[1] - knot[0]) + v0.z * Q[0].w;
	QQ(0, 3) = 3 * Q[0].w / (knot[1] - knot[0]) + v0.w * Q[0].w;

	QQ(K, 0) = 3 * Q[K].x / (knot[K] - knot[K - 1]) - vn.x * Q[K].w;
	QQ(K, 1) = 3 * Q[K].y / (knot[K] - knot[K - 1]) - vn.y * Q[K].w;
	QQ(K, 2) = 3 * Q[K].z / (knot[K] - knot[K - 1]) - vn.z * Q[K].w;
	QQ(K, 3) = 3 * Q[K].w / (knot[K] - knot[K - 1]) - vn.w * Q[K].w;

	Eigen::MatrixXd P = Eigen::MatrixXd::Zero(K + 1, 4);
	P = N.partialPivLu().solve(QQ);
	ON_NurbsCurve onc;
	onc.Create(3, true, 4, L + 1);
	for (int i = 0; i < L + 3; i++)
	{
		onc.SetKnot(i, solve_N.Knot(i));
	}
	onc.SetCV(0, Q[0]);
	for (int i = 1; i < L; i++)
	{
		onc.SetCV(i, ON_4dPoint(P(i - 1, 0), P(i - 1, 1), P(i - 1, 2), P(i - 1, 3)));
	}
	onc.SetCV(L, Q[K]);
	return onc;
}

ON_NurbsCurve ChiralityMath::CubicHomoBsplineInterpolate_Period(const std::vector<ON_4dPoint>& Q, const std::vector<double>& knot)
{
	assert(Q.front() == Q.back());
	assert(Q.size() == knot.size());
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
	Eigen::MatrixXd QQ = Eigen::MatrixXd::Zero(K + 1, 4);
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
		QQ(j, 3) = Q[j].w;
	}
	N(0, 0) = 1 / (knot[1] - knot[0]);
	N(0, K) = 1 / (knot[K] - knot[K - 1]);
	QQ(0, 0) = Q[0].x / (knot[1] - knot[0]) + Q[K].x / (knot[K] - knot[K - 1]);
	QQ(0, 1) = Q[0].y / (knot[1] - knot[0]) + Q[K].y / (knot[K] - knot[K - 1]);
	QQ(0, 2) = Q[0].z / (knot[1] - knot[0]) + Q[K].z / (knot[K] - knot[K - 1]);
	QQ(0, 3) = Q[0].w / (knot[1] - knot[0]) + Q[K].w / (knot[K] - knot[K - 1]);

	double rk1 = knot[K] - knot[K - 1];
	double rk2 = knot[K - 1] - knot[K - 2];
	double r0 = knot[1] - knot[0];
	double r1 = knot[2] - knot[1];
	N(K, 0) = -rk1 / (r1 + r0) - rk1 / r0;
	N(K, 1) = rk1 / (r1 + r0);
	N(K, K - 1) = -r0 / (rk1 + rk2);
	N(K, K) = r0 / rk1 + r0 / (rk1 + rk2);
	QQ(K, 0) = r0 / rk1 * Q[K].x - rk1 / r0 * Q[0].x;
	QQ(K, 1) = r0 / rk1 * Q[K].y - rk1 / r0 * Q[0].y;
	QQ(K, 2) = r0 / rk1 * Q[K].z - rk1 / r0 * Q[0].z;
	QQ(K, 3) = r0 / rk1 * Q[K].w - rk1 / r0 * Q[0].w;
	
	Eigen::MatrixXd P = Eigen::MatrixXd::Zero(K + 1, 4);
	P = N.partialPivLu().solve(QQ);
	ON_NurbsCurve onc;
	onc.Create(3, true, 4, L + 1);
	for (int i = 0; i < L + 3; i++)
	{
		onc.SetKnot(i, solve_N.Knot(i));
	}
	onc.SetCV(0, Q[0]);
	for (int i = 1; i < L; i++)
	{
		onc.SetCV(i, ON_4dPoint(P(i - 1, 0), P(i - 1, 1), P(i - 1, 2), P(i - 1, 3)));
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
	ons.Create(3, false, 4, curve_list[0].Order(), new_curve_list[0].CVCount(), curve_list[0].CVCount());
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
		if (!parent_curve.IsRational())
		{
			parent_curve.GetCV(i, p);
			ons.SetCV(i, 0, p + dir * t0);
			ons.SetCV(i, 1, p + dir * t1);
		}
		else
		{
			parent_curve.GetCV(i, p);
			double w = parent_curve.Weight(i);
			ons.SetCV(i, 0, (p + dir * t0) * w);
			ons.SetCV(i, 1, (p + dir * t1) * w);
			ons.SetWeight(i, 0, w);
			ons.SetWeight(i, 1, w);
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

FrenetFrame ChiralityMath::GetFrenet(const ON_NurbsCurve &onc, double t)
{
	ON_3dPoint p;
	ON_3dVector der;
	ON_3dVector derder;
	onc.Ev2Der(t, p, der, derder);
	der.Unitize();
	ON_3dVector N = ON_3dVector::CrossProduct(der, derder);
	ON_3dVector B = ON_3dVector::CrossProduct(N, der);
	B.Unitize();
	return FrenetFrame(p, der, B);
}

FrenetFrame ChiralityMath::GetFrenet(const ON_NurbsSurface &ons, double u, double v)
{
	ON_3dPoint p;
	ON_3dVector der_u;
	ON_3dVector der_v;
	ons.Ev1Der(u, v, p, der_u, der_v);
	return FrenetFrame(p, der_u, der_v);
}

ON_NurbsCurve ChiralityMath::GetNurbsCircle(const ON_3dPoint &center, const ON_3dPoint &start_point, ON_3dVector normal)
{
	double r = center.DistanceTo(start_point);
	ON_3dVector X = start_point - center;
	ON_3dVector Y = ON_3dVector::CrossProduct(normal, X);
	Y = Y / Y.Length() * r;

	ON_NurbsCurve wiggle(3, true, 3, 7);

	ON_3dPoint CV[7];
	CV[0] = start_point;
	CV[1] = start_point + Y;
	CV[2] = start_point + Y - 2 * X;
	CV[3] = start_point - 2 * X;
	CV[4] = start_point - 2 * X - Y;
	CV[5] = start_point - Y;
	CV[6] = start_point;
	double weight[7] = {1, 0.5, 0.5, 1, 0.5, 0.5, 1};

	for (int i = 0; i < 7; i++)
	{
		wiggle.SetCV(i, CV[i] * weight[i]);
		wiggle.SetWeight(i, weight[i]);
	}

	wiggle.SetKnot(0, 0.0);
	wiggle.SetKnot(1, 0.0);
	wiggle.SetKnot(2, 0.25);
	wiggle.SetKnot(3, 0.5);
	wiggle.SetKnot(4, 0.5);
	wiggle.SetKnot(5, 0.75);
	wiggle.SetKnot(6, 1.0);
	wiggle.SetKnot(7, 1.0);
	return wiggle;
}

ON_4dPoint ChiralityMath::GetHomoPoint(const ON_NurbsCurve& onc, double t)
{
	ON_3dPoint P = onc.PointAt(t);
	if (!onc.IsRational())
	{
		return ON_4dPoint(P);
	}
	ON_NurbsCurve w_onc(1, false, onc.Order(), onc.CVCount());
	for (int i = 0; i < onc.KnotCount(); ++i)
	{
		w_onc.SetKnot(i, onc.Knot(i));
	}
	for (int i = 0; i < onc.CVCount(); ++i)
	{
		w_onc.SetCV(i, ON_3dPoint(onc.Weight(i), 0, 0));
	}
	double W = w_onc.PointAt(t).x;
	return ON_4dPoint(P.x * W, P.y * W, P.z * W, W);
}

ON_4dPoint ChiralityMath::GetHomoDerivative(const ON_NurbsCurve& onc, double t)
{
	if (!onc.IsRational())
	{
		ON_3dVector v = onc.DerivativeAt(t);
		return ON_4dPoint(v.x, v.y, v.z, 0.0);
	}
	ON_NurbsCurve w_onc(1, false, onc.Order(), onc.CVCount());
	ON_NurbsCurve p_onc(3, false, onc.Order(), onc.CVCount());
	for (int i = 0; i < onc.KnotCount(); ++i)
	{
		w_onc.SetKnot(i, onc.Knot(i));
		p_onc.SetKnot(i, onc.Knot(i));
	}
	for (int i = 0; i < onc.CVCount(); ++i)
	{
		w_onc.SetCV(i, ON_3dPoint(onc.Weight(i), 0, 0));
		ON_4dPoint p = onc.ControlPoint(i);
		p_onc.SetCV(i, ON_3dPoint(p.x, p.y, p.z));
	}
	ON_3dVector V = p_onc.DerivativeAt(t);
	double W = w_onc.DerivativeAt(t).x;
	return ON_4dPoint(V.x, V.y, V.z, W);
}

static ON_3dPoint ComputePedal(const ON_Line &line, const ON_3dPoint &p)
{
	ON_3dVector dir = line.Direction();
	dir.Unitize();
	ON_3dPoint o = line.from;
	double t = ON_3dVector::DotProduct(dir, p - o);
	return o + t * dir;
}

ON_NurbsSurface ChiralityMath::GenerateRotating(const ON_NurbsCurve &parent_curve, const ON_Line &axis)
{
	ON_NurbsSurface ons(3, true, parent_curve.Order(), 3, parent_curve.CVCount(), 7);
	ON_3dVector Z = axis.Direction();
	Z.Unitize();
	for (int i = 0; i < ons.KnotCount(0); ++i)
	{
		ons.SetKnot(0, i, parent_curve.Knot(i));
	}
	double w[8] = {0, 0, 0.25, 0.5, 0.5, 0.75, 1, 1};
	for (int j = 0; j < 8; ++j)
	{
		ons.SetKnot(1, j, w[j]);
	}
	for (int i = 0; i < parent_curve.CVCount(); ++i)
	{
		ON_3dPoint p;
		parent_curve.GetCV(i, p);
		ON_3dPoint O = ComputePedal(axis, p);
		ON_3dVector X = p - O;
		ON_3dVector Y = X;
		Y.Rotate(1, 0, Z);
		double weight[7];
		for (int j = 0; j < 7; ++j)
		{
			weight[j] = parent_curve.Weight(i);
			if (j % 3 != 0)
			{
				weight[j] *= 0.5;
			}
		}
		ons.SetCV(i, 0, p * weight[0]);
		ons.SetCV(i, 1, (p + Y) * weight[1]);
		ons.SetCV(i, 2, (p + Y - 2 * X) * weight[2]);
		ons.SetCV(i, 3, (p - 2 * X) * weight[3]);
		ons.SetCV(i, 4, (p - Y - 2 * X) * weight[4]);
		ons.SetCV(i, 5, (p - Y) * weight[5]);
		ons.SetCV(i, 6, p * weight[6]);
		for (int j = 0; j < 7; ++j)
		{
			ons.SetWeight(i, j, weight[j]);
		}
	}
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

ON_3dPoint ChiralityMath::GetRandomPoint(const ON_3dPoint &p, double min_distance, double max_distance)
{
	assert(min_distance > 0 && min_distance <= max_distance);
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dis(0.0, 1.0);
	double rand_double[3];
	for (int i = 0; i < 3; ++i)
	{
		rand_double[i] = dis(gen);
	}
	double rand_distance = rand_double[0] * max_distance + (1 - rand_double[0]) * min_distance;
	double rand_theta = rand_double[1] * 2 * PI;
	double rand_phi = rand_double[2] * (PI / 2) + (1 - rand_double[2]) * (-PI / 2);
	return p + ON_3dVector(cos(rand_theta) * cos(rand_phi),
						   sin(rand_theta) * cos(rand_phi), sin(rand_phi)) *
				   rand_distance;
}

double ChiralityMath::GetRandomFloat(double min_, double max_)
{
	assert(min_ <= max_);
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dis(0.0, 1.0);
	double rand_double = dis(gen);
	double rand_distance = rand_double * max_ + (1 - rand_double) * min_;
	return rand_distance;
}

void ChiralityMath::KnotInsertion(ON_NurbsCurve& onc, double knot_insert)
{
	if (knot_insert < onc.Knot(0) || knot_insert >= onc.Knot(onc.KnotCount() - 1))
	{
		return;
	}
	std::vector<double> Knot;
	Knot.push_back(onc.Knot(0));
	for (int i = 0; i < onc.KnotCount(); ++i)
	{
		Knot.push_back(onc.Knot(i));
	}
	Knot.push_back(Knot.back());
	int l = 0;
	while (l < Knot.size() - 1)
	{
		if (knot_insert >= Knot[l] && knot_insert < Knot[l + 1])
		{
			break;
		}
		l++;
	}
	int k = onc.Order();
	std::vector<ON_4dPoint> vp;
	for (int i = 0; i < onc.CVCount(); ++i)
	{
		vp.push_back(onc.ControlPoint(i));
	}
	std::vector<ON_4dPoint> vq(vp.size() + 1);
	for (int i = 0; i <= l - k + 1; ++i)
	{
		vq[i] = vp[i];
	}
	for (int i = l - k + 2; i <= l; ++i)
	{
		vq[i].x = ((knot_insert - Knot[i]) / (Knot[i + k - 1] - Knot[i])) * vp[i].x
			+ ((Knot[i + k - 1] - knot_insert) / (Knot[i + k - 1] - Knot[i])) * vp[i - 1].x;

		vq[i].y = ((knot_insert - Knot[i]) / (Knot[i + k - 1] - Knot[i])) * vp[i].y
			+ ((Knot[i + k - 1] - knot_insert) / (Knot[i + k - 1] - Knot[i])) * vp[i - 1].y;

		vq[i].z = ((knot_insert - Knot[i]) / (Knot[i + k - 1] - Knot[i])) * vp[i].z
			+ ((Knot[i + k - 1] - knot_insert) / (Knot[i + k - 1] - Knot[i])) * vp[i - 1].z;

		vq[i].w = ((knot_insert - Knot[i]) / (Knot[i + k - 1] - Knot[i])) * vp[i].w
			+ ((Knot[i + k - 1] - knot_insert) / (Knot[i + k - 1] - Knot[i])) * vp[i - 1].w;
	}
	for (int i = l + 1; i < vq.size(); ++i)
	{
		vq[i] = vp[i - 1];
	}
	Knot.insert(Knot.begin() + l + 1, knot_insert);
	onc.Create(onc.Dimension(), onc.IsRational(), onc.Order(), onc.CVCount() + 1);
	for (int i = 0; i < onc.KnotCount(); ++i)
	{
		onc.SetKnot(i, Knot[i + 1]);
	}
	for (int i = 0; i < onc.CVCount(); ++i)
	{
		onc.SetCV(i, vq[i]);
	}
}

std::vector<ON_BezierCurve> ChiralityMath::SplitIntoBezier(const ON_NurbsCurve& onc)
{
	ON_NurbsCurve new_onc = onc;
	std::map<double, int> knot_multi;
	const double* knot = onc.Knot();
	for (int i = 0; i < onc.KnotCount(); ++i)
	{
		knot_multi[*(knot + i)]++;
	}
	knot_multi[*(knot + 0)]++;
	knot_multi[*(knot + (onc.KnotCount() - 1))]++;
	int k = onc.Order();
	for (auto& iter : knot_multi)
	{
		while (iter.second < k - 1)
		{
			KnotInsertion(new_onc, iter.first);
			iter.second++;
		}
	}
	std::vector<ON_BezierCurve> v_obc;
	int j = 0;
	for (int i = 0; i < knot_multi.size() - 1; ++i)
	{
		ON_BezierCurve obc(onc.Dimension(), onc.IsRational(), k);
		for (int ww = 0; ww < k; ++ww)
		{
			obc.SetCV(ww, new_onc.ControlPoint(j));
			++j;
		}
		v_obc.push_back(obc);
		--j;
	}
	return v_obc;
}

std::vector<ON_NurbsSurface> ChiralityMath::WirePlasticCurvedSurface(const std::vector<ParameterCurve>& v_pc, const std::vector<VectorField>& v_tan, int sample_cnt)
{
	assert(v_pc.size() == v_tan.size());
	std::vector<ON_NurbsSurface> v_ons;
	std::vector<std::vector<ON_3dPoint>> Sample_Points;
	std::vector<std::vector<ON_3dVector>> Sample_Vectors;
	for (int i = 0; i < v_pc.size(); ++i)
	{
		std::vector<ON_3dPoint> vp;
		std::vector<ON_3dVector> vv;
		double t0, t1;
		v_pc[i].GetDomain(&t0, &t1);
		for (int j = 0; j < sample_cnt; ++j)
		{
			double t = t0 * (1 - double(j) / double(sample_cnt)) + t1 * double(j) / double(sample_cnt);
			vp.push_back(v_pc[i].GetFrame(t).GetPos());
			vv.push_back(v_tan[i](t));
		}
		Sample_Points.push_back(vp);
		Sample_Vectors.push_back(vv);
	}

	std::vector<double> param;
	for (int i = 0; i <= sample_cnt; ++i)
	{
		param.push_back(i / double(sample_cnt));
	}

	for (int i = 0; i < v_pc.size() - 1; ++i)
	{
		std::vector<ON_BezierCurve> v_clothoid;
		int max_cv = 0;
		for (int j = 0; j < sample_cnt; ++j)
		{
			ON_BezierCurve clothoid = FixAxisBezier3D::Interpolate(Sample_Points[i][j], Sample_Points[i + 1][j],
				Sample_Vectors[i][j], Sample_Vectors[i + 1][j]);
			max_cv = (std::max)(max_cv, clothoid.CVCount());
			v_clothoid.push_back(clothoid);
		}
		for (ON_BezierCurve& obc : v_clothoid)
		{
			while (obc.CVCount() < max_cv)
			{
				ChiralityMath::Elevate(obc);
			}
		}
		std::vector<ON_NurbsCurve> Inter_onc;
		for (int j = 0; j < max_cv; ++j)
		{
			std::vector<ON_3dPoint> Q;
			for (const ON_BezierCurve& obc : v_clothoid)
			{
				Q.push_back(ON_3dPoint(obc.ControlPoint(j)));
			}
			Q.push_back(Q[0]);
			Inter_onc.push_back(ChiralityMath::CubicBsplineInterpolate_Period(Q, param));
		}
		ON_NurbsSurface ons(3, false, v_clothoid[0].Order(), 4, v_clothoid[0].CVCount(), Inter_onc[0].CVCount());
		for (int j = 0; j < ons.KnotCount(0); ++j)
		{
			ons.SetKnot(0, j, (j < ons.KnotCount(0) / 2) ? 0 : 1);
		}
		for (int j = 0; j < ons.KnotCount(1); ++j)
		{
			ons.SetKnot(1, j, Inter_onc[0].Knot(j));
		}
		for (int r = 0; r < ons.CVCount(0); ++r)
		{
			for (int s = 0; s < ons.CVCount(1); ++s)
			{
				ons.SetCV(r, s, Inter_onc[r].ControlPoint(s));
			}
		}
		v_ons.push_back(ons);
	}
	return v_ons;
}

std::vector<ON_Circle> ChiralityMath::TangentToCircle(const ON_Circle& circle1, const ON_Circle& circle2, const ON_Circle& circle3, bool in1, bool in2, bool in3)
{
	std::vector<ON_Circle> result_circle;
	ON_3dPoint center[3] = { circle1.Center(),circle2.Center() ,circle3.Center() };
	assert(abs(ON_3dVector::DotProduct(center[1] - center[0], center[2] - center[1])) > 1e-6);
	ON_Plane center_plane(center[0], center[1], center[2]);
	ON_3dVector Normal = center_plane.Normal();
	assert(circle1.Normal().IsParallelTo(Normal) || circle1.Normal().IsParallelTo(-Normal));
	assert(circle2.Normal().IsParallelTo(Normal) || circle2.Normal().IsParallelTo(-Normal));
	assert(circle3.Normal().IsParallelTo(Normal) || circle3.Normal().IsParallelTo(-Normal));
	double r[3] = { in1 ? -circle1.Radius() : circle1.Radius() ,
					in2 ? -circle2.Radius() : circle2.Radius() ,
					in3 ? -circle3.Radius() : circle3.Radius() };
	double result_R;
	ON_3dPoint result_Center;
	ON_Xform Rotate;
	Rotate.Rotation(Normal, ON_3dVector::ZAxis, ON_3dPoint::Origin);
	double x[3], y[3];
	for (int i = 0; i < 3; ++i)
	{
		center[i].Transform(Rotate);
		x[i] = center[i].x;
		y[i] = center[i].y;
	}
	double ZZ = center[0].z;
	Eigen::Matrix2d M1;
	Eigen::Vector2<double> M2, M3;
	M1(0, 0) = x[0] - x[1]; M1(0, 1) = y[0] - y[1]; M1(1, 0) = x[0] - x[2]; M1(1, 1) = y[0] - y[2];
	M2(0, 0) = 0.5 * (x[0] * x[0] - x[1] * x[1] + y[0] * y[0] - y[1] * y[1] + r[1] * r[1] - r[0] * r[0]);
	M2(1, 0) = 0.5 * (x[0] * x[0] - x[2] * x[2] + y[0] * y[0] - y[2] * y[2] + r[2] * r[2] - r[0] * r[0]);
	M3(0, 0) = r[1] - r[0];
	M3(1, 0) = r[2] - r[0];
	Eigen::Matrix2d M1_inv = M1.inverse();
	Eigen::Vector2<double> M4 = M1_inv * M2;
	Eigen::Vector2<double> M5 = M1_inv * M3;
	double A = M5(0, 0); double C = M5(1, 0);
	double B = M4(0, 0); double D = M4(1, 0);
	double a = A * A + C * C - 1;
	double b = 2 * A * (B - x[0]) + 2 * C * (D - y[0]) - 2 * r[0];
	double c = (B - x[0]) * (B - x[0]) + (D - y[0]) * (D - y[0]) - r[0] * r[0];
	if (abs(a) > 1e-6)
	{
		double delta = b * b - 4 * a * c;
		if (delta >= 0)
		{
			double r1 = (-b + sqrt(delta)) / (2 * a);
			double r2 = (-b - sqrt(delta)) / (2 * a);
			if (abs(r1 - r2) < 1e-8 && r1 > 0 && r2 > 0)
			{
				result_R = (r1 + r2) / 2;
				result_Center.Set(A * result_R + B, C * result_R + D, ZZ);
				result_Center.Transform(Rotate.Inverse());
				result_circle.emplace_back(ON_Circle(center_plane, result_Center, result_R));
			}
			else
			{
				if (r1 > 0)
				{
					result_R = r1;
					result_Center.Set(A * result_R + B, C * result_R + D, ZZ);
					result_Center.Transform(Rotate.Inverse());
					result_circle.emplace_back(ON_Circle(center_plane, result_Center, result_R));
				}
				if (r2 > 0)
				{
					result_R = r2;
					result_Center.Set(A * result_R + B, C * result_R + D, ZZ);
					result_Center.Transform(Rotate.Inverse());
					result_circle.emplace_back(ON_Circle(center_plane, result_Center, result_R));
				}
			}
		}
	}
	else
	{
		result_R = -c / b;
		if (result_R > 0)
		{
			result_Center.Set(A * result_R + B, C * result_R + D, ZZ);
			result_Center.Transform(Rotate.Inverse());
			result_circle.emplace_back(ON_Circle(center_plane, result_Center, result_R));
		}
	}
	return result_circle;
}

std::vector<ON_3dPoint> ChiralityMath::SolveBoundary(FrenetFrame X0, FrenetFrame XN, double step_length, double *kappa_tau_param, int n)
{
	const Eigen::Matrix4d Identity = Eigen::Matrix4d::Identity();
	const int dim = 4 * (n + 1);
	Eigen::MatrixXd Left = Eigen::MatrixXd::Zero(dim, dim);
	Eigen::MatrixXd Right = Eigen::MatrixXd::Zero(dim, 3);
	Eigen::Matrix<double, 4, 3> M_X_0;
	M_X_0 << X0.GetPos().x, X0.GetPos().y, X0.GetPos().z,
		X0.GetAlpha().x, X0.GetAlpha().y, X0.GetAlpha().z,
		X0.GetBeta().x, X0.GetBeta().y, X0.GetBeta().z,
		X0.GetGamma().x, X0.GetGamma().y, X0.GetGamma().z;

	Eigen::Matrix<double, 4, 3> M_X_N;
	M_X_N << XN.GetPos().x, XN.GetPos().y, XN.GetPos().z,
		XN.GetAlpha().x, XN.GetAlpha().y, XN.GetAlpha().z,
		XN.GetBeta().x, XN.GetBeta().y, XN.GetBeta().z,
		XN.GetGamma().x, XN.GetGamma().y, XN.GetGamma().z;

	Right.block<4, 3>(0, 0) = M_X_0;
	Right.block<4, 3>(4 * n, 0) = M_X_N;

	for (int i = 1; i < n; ++i)
	{
		Left.block<4, 4>(i * 4, i * 4 - 4) = Identity;
		Left.block<4, 4>(i * 4, i * 4 + 4) = Identity;
	}
	Left.block<4, 4>(0, 0) = Identity;
	Left.block<4, 4>(n * 4, n * 4) = Identity;

	double a = kappa_tau_param[0];
	double b = kappa_tau_param[1];
	double c = kappa_tau_param[2];
	double d = kappa_tau_param[3];

	std::function<Eigen::Matrix4d(double)> A_s = [a, b, c, d](double s) -> Eigen::Matrix4d
	{
		double k = a * s + b;
		double r = c * s + d;
		Eigen::Matrix4d M;
		M << 0, 0, k, 0,
			0, -k * k, a, k * r,
			0, -a, -k * k - r * r, c,
			0, k * r, -c, -r * r;
		return M;
	};

	for (int i = 1; i < n; ++i)
	{
		double s = step_length * double(i);
		Left.block<4, 4>(i * 4, i * 4) = -2 * Identity - A_s(s) * step_length * step_length;
	}

	Eigen::MatrixXd X = Eigen::MatrixXd::Zero(dim, 3);
	X = Left.partialPivLu().solve(Right);
	std::vector<ON_3dPoint> re(n + 1);
	for (int i = 0; i < n + 1; ++i)
	{
		re[i] = ON_3dPoint(X(i * 4, 0), X(i * 4, 1), X(i * 4, 2));
	}
	return re;
}
