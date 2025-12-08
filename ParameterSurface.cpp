#include "ParameterSurface.h"
#include "ChiralityMathTools.h"

static constexpr double epsilon = 1e-8;

ParameterSurface::ParameterSurface(const Position &pos, const Normal &nor, const double *param_range) : pos_(pos), normal_(nor)
{
    if (param_range != nullptr)
    {
        u_min_ = *(param_range + 0);
        u_max_ = *(param_range + 1);
        v_min_ = *(param_range + 2);
        v_max_ = *(param_range + 3);
    }
}

ON_3dPoint ParameterSurface::PointAt(double u, double v) const
{
    if (u > u_min_ - epsilon && u < u_max_ + epsilon && v > v_min_ - epsilon && v < v_max_ + epsilon)
    {
        return pos_(u, v);
    }
    return ON_3dPoint::Origin;
}

ON_3dVector ParameterSurface::NormalWithoutCheckUnit(double u, double v) const
{
    if (u > u_min_ - epsilon && u < u_max_ + epsilon && v > v_min_ - epsilon && v < v_max_ + epsilon)
    {
        return normal_(u, v);
    }
    return ON_3dVector::XAxis;
}

ON_3dVector ParameterSurface::NormalUnit(double u, double v) const
{
    if (u > u_min_ - epsilon && u < u_max_ + epsilon && v > v_min_ - epsilon && v < v_max_ + epsilon)
    {
        ON_3dVector n = normal_(u, v);
        n.Unitize();
        return n;
    }
    return ON_3dVector::ZeroVector;
}

void ParameterSurface::GetDomain(int dir, double* t0, double* t1) const
{
    if (dir == 0)
    {
        *t0 = u_min_;
        *t1 = u_max_;
        return;
    }
    if (dir == 1)
    {
        *t0 = v_min_;
        *t1 = v_max_;
        return;
    }
}

void ParameterSurface::Discretize(std::vector<std::vector<ON_3dPoint>> &result_points, std::vector<std::vector<ON_3dVector>> &result_normals, int u_sample_num, int v_sample_num) const
{
    result_points.clear();
    result_normals.clear();
    result_points.resize(u_sample_num + 1, std::vector<ON_3dPoint>(v_sample_num + 1, ON_3dPoint::Origin));
    result_normals.resize(u_sample_num + 1, std::vector<ON_3dVector>(v_sample_num + 1, ON_3dVector::ZeroVector));
    for (int i = 0; i <= u_sample_num; ++i)
    {
        for (int j = 0; j <= v_sample_num; ++j)
        {
            double u = u_min_ * (1 - double(i) / double(u_sample_num)) + u_max_ * double(i) / double(u_sample_num);
            double v = v_min_ * (1 - double(j) / double(v_sample_num)) + v_max_ * double(j) / double(v_sample_num);
            result_points[i][j] = pos_(u, v);
            result_normals[i][j] = normal_(u, v);
            result_normals[i][j].Unitize();
        }
    }
}

ON_NurbsCurve ParameterCurve::CubicNurbsApproximate(int N) const
{
    std::vector<ON_3dPoint> vp;
    std::vector<double> s_param;
    for (int i = 0; i <= N; ++i)
    {
        double t = t_min_ * (1 - double(i) / double(N)) + t_max_ * double(i) / double(N);
        s_param.push_back(t);
        vp.push_back(GetFrame(t).GetPos());
    }
    ON_3dVector vs = GetDerivative(t_min_);
    ON_3dVector ve = GetDerivative(t_max_);
    return ChiralityMath::CubicBsplineInterpolate_G1(vp, s_param, vs, ve);
}
