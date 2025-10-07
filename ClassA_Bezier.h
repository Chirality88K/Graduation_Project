#ifndef CLASSA_BEZIER_H
#define CLASSA_BEZIER_H
#include "thirdparty/opennurbs/opennurbs.h"
#include <array>
#include <vector>

using Matrix3x3 = std::array<double, 9>;
class ClassA_Bezier
{
public:
    ClassA_Bezier(const std::vector<ON_3dPoint>& vp):mCV(vp){}
    static void ClassATest(ONX_Model* model);
    static void ComputeGeneratorMatrix(const ON_3dPoint &ps, const ON_3dPoint &pe, ON_3dVector vs, ON_3dVector ve, Matrix3x3 &data);
    std::vector<double> ComputeAngle()const;
    std::vector<double> ComputeLength()const;
    std::vector<ON_3dVector> ComputeNormal()const;
    void Iteration();
    void Iteration_in_Plane();
    static void TestSmooth(ONX_Model* model);
    void Elevate();
    ON_BezierCurve GetBezier()const;
private:
    std::vector<ON_3dPoint> mCV;
};
#endif