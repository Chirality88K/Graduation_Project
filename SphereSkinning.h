#ifndef SPHERESKINNING_H
#define SPHERESKINNING_H
#include "thirdparty/opennurbs/opennurbs.h"
#include "ChiralityMathTools.h"
#include <vector>
#include <assert.h>

class SphereSkinning
{
public:
    SphereSkinning() {}
    void AddSphere(const ON_Sphere &sp)
    {
        spheres_.push_back(sp);
    }
    void GetCirclesToInterpolate() const;
    void GetConicVertex() const;

    ON_Sphere GetSphere(int index) const
    {
        assert(index >= 0 && index < spheres_.size());
        return spheres_[index];
    }
    std::vector<ON_Circle> Get_Circles_For_Debug() const
    {
        return circles_;
    }
    std::vector<ON_NurbsSurface> Skinning() const;

private:
    void GetParamCurveAndVectorField(std::vector<ParameterCurve>& v_pc, std::vector<ChiralityMath::VectorField>& v_vf) const;

private:
    std::vector<ON_Sphere> spheres_;
    mutable std::vector<ON_Circle> circles_;
    mutable std::vector<ON_4dPoint> conic_vertices_;
};

#endif