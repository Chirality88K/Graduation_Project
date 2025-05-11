#ifndef CONICSPIRAL_H
#define CONICSPIRAL_H
#include "thirdparty/opennurbs/opennurbs.h"
class PolarPoint3d
{
public:
    PolarPoint3d() : mDistance(0.0), mTheta(0.0), mPhi(0.0) {}
    PolarPoint3d(double theta, double phi, double dis);
    PolarPoint3d(const ON_3dPoint &p);
    ON_3dPoint CartesianCoordinates();
    double mDistance; // 到原点的距离
    double mTheta;    // 水平角[0,2\pi)
    double mPhi;      // 俯仰角[-\pi/2,\pi/2]
};

class ConicSpiral
{
public:
    ConicSpiral(PolarPoint3d PS, PolarPoint3d PE, ON_3dVector vs, ON_3dVector ve, int num = 50);
private:
    int mNum_Points;
    std::vector<ON_3dPoint> mDiscretePolygon;
};
#endif