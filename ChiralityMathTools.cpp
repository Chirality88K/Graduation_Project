#include "ChiralityMathTools.h"
#include "thirdparty/eigen/Eigen/Dense"

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

double ChiralityMath::Torsion(const ON_BezierCurve &obc, double t)
{
	if (obc.Degree() < 3)
	{
		return 0.0;
	}
	ON_3dVector v1;
	ON_3dVector v2;
	ON_3dVector v3;
	ON_3dPoint p;
	obc.Ev2Der(t, p, v1, v2);
	int n = obc.Degree();
	ON_3dPoint *cp = new ON_3dPoint[n + 1];
	for (int i = 0; i < n + 1; ++i)
	{
		obc.GetCV(i, *(cp + i));
	}
	ON_BezierCurve obc3der(3, false, obc.Order() - 3);
	for (int i = 0; i < obc3der.CVCount(); ++i)
	{
		obc3der.SetCV(i, ON_3dPoint(n * (n - 1) * (n - 2) * (cp[i + 3] - 3 * cp[i + 2] + 3 * cp[i + 1] - cp[i])));
	}
	v3 = obc3der.PointAt(t);
	delete[] cp;
	return ON_3dVector::DotProduct(ON_3dVector::CrossProduct(v1, v2), v3) / (ON_3dVector::CrossProduct(v1, v2).LengthSquared());
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
