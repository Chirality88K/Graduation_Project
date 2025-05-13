#include "ConicSpiral.h"
extern const double PI;
PolarPoint3d::PolarPoint3d(double theta, double phi, double dis)
{
    mDistance = dis;
    while (theta < 0.0)
    {
        theta += 2 * PI;
    }
    while (theta >= 2 * PI)
    {
        theta -= 2 * PI;
    }
    mTheta = theta;
    phi = (std::min)(PI / 2, phi);
    phi = (std::max)(-PI / 2, phi);
    mPhi = phi;
}

PolarPoint3d::PolarPoint3d(const ON_3dPoint &p)
{
    mDistance = p.DistanceTo(ON_3dPoint(0, 0, 0));
    if (mDistance < 1e-8)
    {
        mTheta = 0.0;
        mPhi = 0.0;
        return;
    }
    double horizontal_distance = sqrt(p.x * p.x + p.y * p.y);
    if (horizontal_distance < 1e-8)
    {
        mTheta = 0.0;
        mPhi = p.z > 0 ? PI / 2 : -PI / 2;
        return;
    }
    double cos_theta = p.x / horizontal_distance;
    cos_theta = (std::max)(-1.0, cos_theta);
    cos_theta = (std::min)(1.0, cos_theta);
    mTheta = p.y >= 0.0 ? acos(cos_theta) : -acos(cos_theta);
    mPhi = asin(p.z / mDistance);
}

ON_3dPoint PolarPoint3d::CartesianCoordinates()
{
    ON_3dPoint p;
    p.x = mDistance * cos(mPhi) * cos(mTheta);
    p.y = mDistance * cos(mPhi) * sin(mTheta);
    p.z = mDistance * sin(mPhi);
    return p;
}

ConicSpiral::ConicSpiral(PolarPoint3d PS, PolarPoint3d PE, ON_3dVector vs, ON_3dVector ve, int num)
{
    mNum_Points = (std::max)(10, num);
    double theta0 = PS.mTheta;
    double theta1 = PE.mTheta;
    if (theta1 < theta0)
    {
        theta1 += 2 * PI;
    }
    std::vector<PolarPoint3d> polarpoints;
    std::vector <double> all_theta(mNum_Points, 0);
    std::vector <double> all_phi(mNum_Points, 0);
    for (int i = 0; i < mNum_Points; ++i)
    {
        all_theta[i] = theta0 * (1 - double(i) / double(mNum_Points - 1)) + theta1 * double(i) / double(mNum_Points - 1);
    }
}