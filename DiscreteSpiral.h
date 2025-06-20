#ifndef DISCRETESPIRAL_H
#define DISCRETESPIRAL_H
#include "Frenet.h"
#include <vector>

class DiscreteSpiral
{
public:
    DiscreteSpiral(const std::vector<ON_3dPoint> &vp) : mPoints(vp) {}

private:
    std::vector<ON_3dPoint> mPoints;
    std::vector<FrenetFrame> mFrenetFrame;
    void ComputeFrenetFrame();
    void SmoothFrenetFrame();
    std::vector<double> ComputeDiscreteCurvature() const;
    std::vector<double> ComputeDiscreteTorsion() const;
};
#endif