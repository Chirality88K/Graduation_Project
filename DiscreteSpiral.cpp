#include "DiscreteSpiral.h"

void DiscreteSpiral::ComputeFrenetFrame()
{
    int num = mPoints.size();
    mFrenetFrame.clear();
    mFrenetFrame.resize(num);
    for (int i = 1; i < num - 1; ++i)
    {
        ON_3dVector a = mPoints[i + 1] - mPoints[i - 1];
        ON_3dVector v1 = mPoints[i] - mPoints[i - 1];
        ON_3dVector v2 = mPoints[i + 1] - mPoints[i];
        if (v1.IsParallelTo(v2))
        {
            CHIRALITY_WARN(std::string("v1 and v2 may be parallel, which will lead to mistakes!!"));
        }
        ON_3dVector c = ON_3dVector::CrossProduct(v1, v2);
        ON_3dVector b = ON_3dVector::CrossProduct(c, a);
        mFrenetFrame[i] = FrenetFrame(mPoints[i], a, b);
    }
}

void DiscreteSpiral::SmoothFrenetFrame()
{
}

std::vector<double> DiscreteSpiral::ComputeDiscreteCurvature() const
{
    int num = mPoints.size();
    std::vector<double> cur(num, 0.0);
    for (int i = 1; i < num - 1; ++i)
    {
        ON_3dVector v1 = mPoints[i] - mPoints[i - 1];
        ON_3dVector v2 = mPoints[i + 1] - mPoints[i];
        double product = ON_3dVector::DotProduct(v1, v2) / v1.Length() / v2.Length();
        if (product <= 0)
        {
            CHIRALITY_ERROR(std::string("angle is too large!!!"));
        }
        product = (std::min)(product, 1.0);
        cur[i] = sqrt(1.0 / product - 1.0) / (mPoints[i + 1] - mPoints[i - 1]).Length() * 2;
    }
    return cur;
}

std::vector<double> DiscreteSpiral::ComputeDiscreteTorsion() const
{
    int num = mPoints.size();
    std::vector<double> tor(num, 0.0);
    for (int i = 1; i < num - 1; ++i)
    {
        ON_3dVector v1 = mFrenetFrame[i - 1].GetGamma();
        ON_3dVector v2 = mFrenetFrame[i].GetGamma();
        ON_3dVector v3 = mFrenetFrame[i + 1].GetGamma();
        v1 = (v1 + v2) / 2;
        v2 = (v2 + v3) / 2;

        double product = ON_3dVector::DotProduct(v1, v2) / v1.Length() / v2.Length();
        if (product <= 0)
        {
            CHIRALITY_ERROR(std::string("angle is too large!!!"));
        }
        product = (std::min)(product, 1.0);
        tor[i] = sqrt(1.0 / product - 1.0) / (mPoints[i + 1] - mPoints[i - 1]).Length() * 2;
    }
    return tor;
}
