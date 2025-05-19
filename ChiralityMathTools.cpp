#include "ChiralityMathTools.h"

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
