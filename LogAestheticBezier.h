#ifndef LOGAESTHETICBEZIER_H
#define LOGAESTHETICBEZIER_H
#include "thirdparty/opennurbs/opennurbs.h"
#include <vector>
#include <assert.h>

class LogAestheticBezier2D
{
public:
    static ON_BezierCurve Interpolate(ON_3dPoint ps, ON_3dPoint pe, ON_3dVector vs, ON_3dVector ve, double alpha = -1.0);

private:
    LogAestheticBezier2D(const std::vector<ON_3dPoint> &vp, double alpha) : points_(vp), alpha_(alpha)
    {
        assert(vp.size() > 3);
        VS_ = vp[1] - vp[0];
        VE_ = vp[vp.size() - 1] - vp[vp.size() - 2];
        VS_.Unitize();
        VE_.Unitize();
    }
    void Iteration();
    bool Check() const;
    bool MonoCurvatureCheck() const;
    ON_BezierCurve GetBezier() const;
    std::vector<double> ComputeLength() const;
    std::vector<double> ComputeAngle() const;

private:
    std::vector<ON_3dPoint> points_;
    double alpha_ = -1.0;
    ON_3dVector VS_;
    ON_3dVector VE_;
};

#endif