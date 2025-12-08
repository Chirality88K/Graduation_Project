#include "EditingCurvatureCurve.h"
#include "ChiralityMathTools.h"

ON_BezierCurve EditingCurvatureCurve::BezierInterpolate(ON_3dPoint ps, ON_3dPoint pe, ON_3dVector vs, ON_3dVector ve, const CurvatureFunction &cf)
{
    double length = ps.DistanceTo(pe) / 4;
    vs.Unitize();
    ve.Unitize();
    std::vector<ON_3dPoint> vp = { ps,ps + vs * length,pe - ve * length,pe };
    EditingCurvatureCurve log_bezier(vp, cf);
    const int max_cv_cnt = 50;
    while (log_bezier.points_.size() < max_cv_cnt && !log_bezier.Check())
    {
        ChiralityMath::Elevate(log_bezier.points_);
        const int max_iter = 10000;
        int iter_cnt = 0;
        while (iter_cnt < max_iter && !log_bezier.Check())
        {
            ++iter_cnt;
            log_bezier.Iteration();
        }
    }
    return log_bezier.GetBezier();
}

void EditingCurvatureCurve::Iteration()
{
    std::vector<double> Angle = ComputeAngle();
    for (double angle : Angle)
    {
        if (angle >= PI / 2 || angle <= -PI / 2)
        {
            return;
        }
    }
    for (int i = 2; i < points_.size() - 2; ++i)
    {
        std::vector <double> arc = ComputeDiscreteArcLength();
        double k0 = cur_func_(arc[i - 1]);
        double k1 = cur_func_(arc[i]);
        double k2 = cur_func_(arc[i + 1]);
        assert(abs(k0 - k2) > 1e-6);
        double lambda = (k1 - k0) / (k2 - k0);
        double new_angle = Angle[i - 2] * (1 - lambda) + Angle[i] * lambda;
        Angle[i - 1] = new_angle;
        ON_3dVector V = points_[i + 1] - points_[i - 1];
        V.Rotate(1, 0, ON_3dVector::ZAxis);
        points_[i] = (points_[i - 1] + points_[i + 1]) / 2 -
            (tan(new_angle / 2) * 0.5) * V;
    }
    std::vector<double> Length = ComputeLength();
    double mean_length = 0.0;
    for (double l : Length)
    {
        mean_length += l;
    }
    mean_length /= Length.size();
    points_[1] = points_[0] + VS_ * mean_length;
    points_[points_.size() - 2] = points_.back() - VE_ * mean_length;
}

ON_BezierCurve EditingCurvatureCurve::GetBezier() const
{
    ON_BezierCurve obc(2, false, points_.size());
    int i = 0;
    for (const ON_3dPoint& p : points_)
    {
        obc.SetCV(i, p);
        ++i;
    }
    return obc;
}

std::vector<double> EditingCurvatureCurve::ComputeLength() const
{
    std::vector<double> result;
    result.reserve(points_.size() - 1);
    for (int i = 1; i < points_.size(); ++i)
    {
        result.push_back(points_[i - 1].DistanceTo(points_[i]));
    }
    return result;
}

std::vector<double> EditingCurvatureCurve::ComputeAngle() const
{
    std::vector<double> result;
    result.reserve(points_.size() - 2);
    for (int i = 2; i < points_.size(); ++i)
    {
        result.push_back(ChiralityMath::ComputeSignedAngle(points_[i - 1] - points_[i - 2], points_[i] - points_[i - 1]));
    }
    return result;
}

std::vector<double> EditingCurvatureCurve::ComputeDiscreteArcLength() const
{
    std::vector<double> arc;
    double l = 0.0;
    for (int i = 1; i < points_.size(); ++i)
    {
        arc.push_back(l);
        l += points_[i].DistanceTo(points_[i - 1]);
    }
    arc.push_back(l);
    return arc;
}

bool EditingCurvatureCurve::Check() const
{
    std::vector<double> arc = ComputeDiscreteArcLength();
    std::vector<double> angle = ComputeAngle();
    if (angle.size() < 2)
    {
        return false;
    }
    for (int i = 2; i < angle.size(); ++i)
    {
        double angle_delta1 = angle[i - 1] - angle[i - 2];
        double angle_delta2 = angle[i] - angle[i - 1];
        double k_delta1 = cur_func_(arc[i]) - cur_func_(arc[i - 1]);
        double k_delta2 = cur_func_(arc[i + 1]) - cur_func_(arc[i]);
        if (abs(angle_delta1 / angle_delta2 - k_delta1 / k_delta2) > 1e-6)
        {
            return false;
        }
    }
    std::vector<double> Length = ComputeLength();
    double sum = 0.0;
    for (double l : Length)
    {
        sum += l;
    }
    double avg = sum / Length.size();
    sum = 0.0;
    for (double l : Length)
    {
        sum += pow(l - avg, 2);
    }
    double cv = sqrt(sum / Length.size()) / avg;
    if (cv > 0.1)
    {
        return false;
    }
    return true;
}
