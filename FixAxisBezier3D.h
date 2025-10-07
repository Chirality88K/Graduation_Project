#ifndef FIXAXISBEZIER3D
#define FIXAXISBEZIER3D
#include "thirdparty/opennurbs/opennurbs.h"
#include <vector>

class FixAxisBezier3D
{
public:
    FixAxisBezier3D(const std::vector<ON_3dPoint> &vp) : mCV(vp)
    {
        mPS = vp[0];
        mPE = vp.back();
        mVS = vp[1] - vp[0];
        mVS.Unitize();
        mVE = vp[vp.size() - 1] - vp[vp.size() - 2];
        mVE.Unitize();
    }
    ~FixAxisBezier3D() = default;
    ON_BezierCurve Project_to_Plane(ON_3dVector T) const;
    ON_BezierCurve Go_Back_To_Space(const ON_BezierCurve&, ON_3dVector T, double* error) const;
    std::pair<double, double> IterationForFirstTime(double* error) const;
    ON_BezierCurve BisectionIteration(double start, double end, double* error, ON_3dVector& TT) const;
    static ON_BezierCurve Interpolate(ON_3dPoint ps, ON_3dPoint pe, ON_3dPoint vs, ON_3dPoint ve);
    static void Test(ONX_Model* model);
    static void GenerateDNA(ONX_Model* model);

private:
    static void SmoothingWithFixedT(ON_BezierCurve& obc, ON_3dVector T);
    static std::vector<double> ComputeLength(const ON_BezierCurve& obc);
    static std::vector<double> ComputeRotatingAngles(const ON_BezierCurve& obc, ON_3dVector T);
    static bool CheckCurvatureMono(const ON_BezierCurve& obc);

private:
    std::vector<ON_3dPoint> mCV;
    ON_3dVector mVS;
    ON_3dVector mVE;
    ON_3dPoint mPS;
    ON_3dPoint mPE;
};

#endif