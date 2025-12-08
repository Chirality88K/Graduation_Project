#ifndef PARAMETERSURFACE_H
#define PARAMETERSURFACE_H
#include <functional>
#include "thirdparty/opennurbs/opennurbs.h"
#include "Frenet.h"
#include <vector>

class ParameterSurface
{
    using Position = std::function<ON_3dPoint(double, double)>;
    using Normal = std::function<ON_3dVector(double, double)>;

public:
    ParameterSurface(const Position &pos, const Normal &nor, const double *param_range);
    ON_3dPoint PointAt(double u, double v) const;
    ON_3dVector NormalWithoutCheckUnit(double u, double v) const;
    ON_3dVector NormalUnit(double u, double v) const;
    void GetDomain(int dir, double *t0, double *t1) const;
    void Discretize(std::vector<std::vector<ON_3dPoint>> &result_points, std::vector<std::vector<ON_3dVector>> &result_normals, int u_sample_num = 100, int v_sample_num = 100) const;

private:
    Position pos_;
    Normal normal_;
    double u_min_ = 0.0;
    double u_max_ = 1.0;
    double v_min_ = 0.0;
    double v_max_ = 1.0;
};

class ParameterCurve
{
    using Frame = std::function<FrenetFrame(double)>;
    using Derivative = std::function<ON_3dVector(double)>;

public:
    ParameterCurve(const Frame &frame, const double *param_range) : frame_(frame)
    {
        t_min_ = *param_range;
        t_max_ = *(param_range + 1);
    }
    FrenetFrame GetFrame(double t) const
    {
        return frame_(t);
    }
    void GetDomain(double* t0, double* t1) const
    {
        *t0 = t_min_;
        *t1 = t_max_;
    }
    void SetDerivative(const Derivative& d)
    {
        der_ = d;
    }
    ON_3dVector GetDerivative(double t) const 
    {
        return der_(t);
    }
    ON_NurbsCurve CubicNurbsApproximate(int sample_cnt) const;

private:
    Frame frame_;
    Derivative der_;
    double t_min_ = 0.0;
    double t_max_ = 1.0;
};

#endif