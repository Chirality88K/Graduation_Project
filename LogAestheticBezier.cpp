#include "LogAestheticBezier.h"
#include "ChiralityMathTools.h"

ON_BezierCurve LogAestheticBezier2D::Interpolate(ON_3dPoint ps, ON_3dPoint pe, ON_3dVector vs, ON_3dVector ve, double alpha)
{
    double length = ps.DistanceTo(pe) / 4;
    vs.Unitize();
    ve.Unitize();
    std::vector<ON_3dPoint> vp = { ps,ps + vs * length,pe - ve * length,pe };
    LogAestheticBezier2D log_bezier(vp, alpha);
    const int max_cv_cnt = 50;
    while (log_bezier.points_.size() < max_cv_cnt && !log_bezier.MonoCurvatureCheck())
    {
        ChiralityMath::Elevate(log_bezier.points_);
        const int max_iter = 10000;
        int iter_cnt = 0;
        while (!log_bezier.Check() && iter_cnt < max_iter)
        {
            ++iter_cnt;
            log_bezier.Iteration();
        }
    }
    return log_bezier.GetBezier();
}

static inline double SIGN(double d)
{
    if (d > 0)
    {
        return 1.0;
    }
    if (d < 0)
    {
        return -1.0;
    }
    return 0.0;
}

void LogAestheticBezier2D::Iteration()
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
        double new_angle = (pow(abs(Angle[i - 2]), -alpha_) * SIGN(Angle[i - 2]) +
                            pow(abs(Angle[i - 1]), -alpha_) * SIGN(Angle[i - 1]) +
                            pow(abs(Angle[i]), -alpha_) * SIGN(Angle[i])) /
                           3;
        new_angle = pow(abs(new_angle), -1.0 / alpha_) * SIGN(new_angle);
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

bool LogAestheticBezier2D::Check() const
{
    std::vector<double> Angle = ComputeAngle();
    for (int i = 2; i < Angle.size(); ++i)
    {
        double pow_error = 2 * pow(abs(Angle[i - 1]), -alpha_) * SIGN(Angle[i - 1]) -
            pow(abs(Angle[i - 2]), -alpha_) * SIGN(Angle[i - 2]) -
            pow(abs(Angle[i]), -alpha_) * SIGN(Angle[i]);
        if (abs(pow_error) > 1e-6)
        {
            return false;
        }
    }
    std::vector<double> Length = ComputeLength();
    double mean_length = 0.0;
    for (double l : Length)
    {
        mean_length += l / Length.size();
    }
    double cv = 0.0;
    for (double l : Length)
    {
        cv += (l - mean_length) * (l - mean_length);
    }
    cv = sqrt(cv / Length.size()) / mean_length;
    if (cv > 0.1)
    {
        return false;
    }
    return true;
}

bool LogAestheticBezier2D::MonoCurvatureCheck() const
{
    std::vector<double> Angle = ComputeAngle();
    double Deltatheta = Angle.back() - Angle.front();
    int n = points_.size() - 1;
    assert(int(Angle.size()) == n - 1);
    double s0 = (n + 1) * sin(Angle[0]) + (n - 2) * sin(Angle[0] + Angle[1]) - 3 * (n - 1) * sin(Angle[0]) * cos(Angle[0]);
    double s1 = -(n + 1) * sin(Angle[n - 2]) - (n - 2) * sin(Angle[n - 3] + Angle[n - 2]) + 3 * (n - 1) * sin(Angle[n - 2]) * cos(Angle[n - 2]);
    if (Deltatheta * s0 >= 0 && s0 * s1 >= 0)
    {
        return Check();
    }
    else
    {
        return false;
    }
}

ON_BezierCurve LogAestheticBezier2D::GetBezier() const
{
    ON_BezierCurve obc(3, false, points_.size());
    int i = 0;
    for (const ON_3dPoint &p : points_)
    {
        obc.SetCV(i, p);
        ++i;
    }
    return obc;
}

std::vector<double> LogAestheticBezier2D::ComputeLength() const
{
    std::vector<double> result;
    result.reserve(points_.size() - 1);
    for (int i = 1; i < points_.size(); ++i)
    {
        result.push_back(points_[i - 1].DistanceTo(points_[i]));
    }
    return result;
}

std::vector<double> LogAestheticBezier2D::ComputeAngle() const
{
    std::vector<double> result;
    result.reserve(points_.size() - 2);
    for (int i = 2; i < points_.size(); ++i)
    {
        result.push_back(ChiralityMath::ComputeSignedAngle(points_[i - 1] - points_[i - 2], points_[i] - points_[i - 1]));
    }
    return result;
}
