#ifndef EDITCURATURECURVE_H
#define EDITCURATURECURVE_H
#include "thirdparty/opennurbs/opennurbs.h"
#include <vector>
#include <functional>
#include <assert.h>

class EditingCurvatureCurve
{
    using CurvatureFunction = std::function<double(double)>;

public:
    static ON_BezierCurve BezierInterpolate(ON_3dPoint ps, ON_3dPoint pe, ON_3dVector vs, ON_3dVector ve, const CurvatureFunction &cf);

private:
    EditingCurvatureCurve(const std::vector<ON_3dPoint> &vp, const CurvatureFunction &cf)
    {
        assert(vp.size() > 3);
        VS_ = vp[1] - vp[0];
        VE_ = vp[vp.size() - 1] - vp[vp.size() - 2];
        VS_.Unitize();
        VE_.Unitize();
        points_ = vp;
        cur_func_ = cf;
    }
    void Iteration();
    ON_BezierCurve GetBezier() const;
    std::vector<double> ComputeLength() const;
    std::vector<double> ComputeAngle() const;
    std::vector<double> ComputeDiscreteArcLength() const;
    bool Check() const;

private:
    std::vector<ON_3dPoint> points_;
    CurvatureFunction cur_func_;
    ON_3dVector VS_;
    ON_3dVector VE_;
};

#endif