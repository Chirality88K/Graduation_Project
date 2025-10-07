#include "ClassA_Bezier.h"
#include "write3dm.h"
#include "ChiralityMathTools.h"
extern const double PI;

void ClassA_Bezier::ClassATest(ONX_Model* model)
{
    ON_3dPoint StartPoint = ON_3dPoint::Origin;
    ON_3dVector StartVector = ON_3dVector::XAxis;
    const int N = 10;
    ON_3dVector Rotating_axis = ON_3dVector(1, 0, 1);
    double rotate_angle = PI / 18.0;
    double s = 1.0 / cos(rotate_angle) * 1.1;
    ON_Xform R;
    R.Rotation(rotate_angle, Rotating_axis, ON_3dPoint::Origin);
    ON_3dVector* vector_array = new ON_3dVector[N - 1]();
    *vector_array = StartVector;
    for (int i = 1; i < N - 1; ++i) 
    {
        *(vector_array + i) = *(vector_array + (i - 1));
        (vector_array + i)->Transform(R);
        (*(vector_array + i)) *= s;
    }
    ON_3dPoint* control_p = new ON_3dPoint[N]();
    *control_p = StartPoint;
    for (int i = 1; i < N; ++i) 
    {
        *(control_p + i) = *(control_p + (i - 1)) + *(vector_array + (i - 1));
    }
    ON_BezierCurve obc(3, false, N);
    for (int i = 0; i < N; ++i)
    {
        obc.SetCV(i, *(control_p + i));
    }
    delete[]vector_array;
    delete[]control_p;
    const int layer_index = model->AddLayer(L"ClassA", ON_Color::SaturatedGold);
    ChiralityAddNurbsCurve(model, obc, L"ClassA_Test", layer_index);
    ChiralityDebugforR(obc, "ClassA_Test");
}

void ClassA_Bezier::ComputeGeneratorMatrix(const ON_3dPoint &ps, const ON_3dPoint &pe, ON_3dVector vs, ON_3dVector ve, Matrix3x3 &data)
{
    int N = 10;
    ON_3dVector vp = pe - ps;
    vp.Unitize();
    vs.Unitize();
    ve.Unitize();
    ON_3dVector n = ve - vs;
    n.Unitize();//n是vs,ve角平分面的法向量
    ON_3dVector Rotating_axis = ON_3dVector::CrossProduct(n, vp);
    Rotating_axis.Unitize();
    double cos_vs_ve = ON_3dVector::DotProduct(vs, ve);
    double theta = acos(cos_vs_ve);//angle of vs and ve
    double alpha = acos(ON_3dVector::DotProduct(vs, Rotating_axis));//angle of axis and vs(ve)
    double phi = acos(1 - (1 - cos_vs_ve) / (sin(alpha) * sin(alpha)));//real rotate angle
    ON_Xform R;
    R.Rotation(phi / double(N - 1), Rotating_axis, ON_3dPoint::Origin);

}

std::vector<double> ClassA_Bezier::ComputeAngle() const
{
    std::vector<double> re;
    re.push_back(0.0);
    for (int i = 1; i < mCV.size() - 1; ++i) 
    {
        double product = ON_3dVector::DotProduct(mCV[i] - mCV[i - 1], mCV[i + 1] - mCV[i]) / (mCV[i] - mCV[i - 1]).Length() / (mCV[i + 1] - mCV[i]).Length();
        if (product > 1) {
            product = 1;
        }
        re.push_back(acos(product));
    }
    re.push_back(0);
    return re;
}

std::vector<double> ClassA_Bezier::ComputeLength() const
{
    std::vector<double> re;
    for (int i = 1; i < mCV.size(); ++i)
    {
        re.push_back((mCV[i] - mCV[i - 1]).Length());
    }
    return re;
}

std::vector<ON_3dVector> ClassA_Bezier::ComputeNormal() const
{
    std::vector<ON_3dVector> re;
    re.push_back(ON_3dVector(0,0,0));
    for (int i = 1; i < mCV.size() - 1; ++i)
    {
        ON_3dVector v = ON_3dVector::CrossProduct(mCV[i] - mCV[i - 1], mCV[i + 1] - mCV[i]);
        v.Unitize();
        re.push_back(v);
    }
    re.push_back(ON_3dVector(0, 0, 0));
    return re;
}

void ClassA_Bezier::Iteration()
{
    std::vector<ON_3dVector> Normal = ComputeNormal();
    std::vector<double> Angle = ComputeAngle();
    std::vector<double> Length = ComputeLength();
    std::vector<ON_3dPoint> new_points = mCV;
    for (size_t i = 2; i < mCV.size() - 2; ++i)
    {
        double theta = (Angle[i - 1] + Angle[i] + Angle[i + 1]) / 3;
        ON_3dVector normal = (Normal[i - 1] + Normal[i] + Normal[i + 1]) / 3;
        normal.Unitize();
        double s1 = Length[i - 1] / Length[i - 2];
        double s2 = Length[i] / Length[i - 1];
        double s3 = Length[i + 1] / Length[i];
        double s = (s1 + s2 + s3) / 3;
        ON_3dVector beta = ON_3dVector::CrossProduct(normal, mCV[i + 1] - mCV[i - 1]);
        beta.Unitize();
        double L = (mCV[i + 1] - mCV[i - 1]).Length();
        double x = L / sqrt(1 + s * s + 2 * s * cos(theta));
        double cos_alpha = (x * x + L * L - s * s * x * x) / (2 * x * L);
        double sin_alpha = sqrt(1 - cos_alpha * cos_alpha);
        new_points[i] = mCV[i - 1] + x * cos_alpha / L * (mCV[i + 1] - mCV[i - 1]) + x * sin_alpha * (-beta);
    }
    mCV = new_points;
    Length = ComputeLength();
    double sum = 0;
    for (size_t i = 1; i < Length.size(); ++i)
    {
        sum += Length[i] / Length[i - 1];
    }
    double ratio_avg = sum / (Length.size() - 1);
    mCV[1] = mCV[0] + Length[1] / ratio_avg * (mCV[1] - mCV[0]) / (mCV[1] - mCV[0]).Length();
    size_t N = mCV.size();
    mCV[N - 2] = mCV[N - 1] - Length[N - 3] * ratio_avg * (mCV[N - 1] - mCV[N - 2]) / (mCV[N - 1] - mCV[N - 2]).Length();
}

static void UpdatePoints_in_Plane(const ON_3dPoint& p0, ON_3dPoint& p1, ON_3dPoint& p2, const ON_3dPoint& p3)
{
    ON_3dVector v0 = p1 - p0;
    ON_3dVector v2 = p3 - p2;
    ON_3dVector v3 = p3 - p0;
    v0.Unitize(); v2.Unitize(); v3.Unitize();
    double L = p0.DistanceTo(p3);
    double cos_2theta = ON_3dVector::DotProduct(v0, v2);
    double sin_2theta = sqrt(1 - cos_2theta * cos_2theta);
    double cos_alpha = ON_3dVector::DotProduct(v0, v3);
    double sin_alpha = sqrt(1 - cos_alpha * cos_alpha);
    double cos_beta = ON_3dVector::DotProduct(v2, v3);
    double sin_beta = sqrt(1 - cos_beta * cos_beta);
    double a = L / sin_2theta * sin_alpha;
    double b = L / sin_2theta * sin_beta;
    double Delta = (a + b) * (a + b) - 4 * a * b * (2 * cos_2theta - 1);
    double x1 = (a + b + sqrt(Delta)) / (2 * (2 * cos_2theta - 1));
    double x2 = (a + b - sqrt(Delta)) / (2 * (2 * cos_2theta - 1));
    double x = (std::min)(x1, x2);
    if (x < 0) {
        CHIRALITY_ERROR(std::string("x<0!!!!!!!!"));
    }
    p1 = p0 + v0 * (b - x);
    p2 = p3 - v2 * (a - x);
}

void ClassA_Bezier::Iteration_in_Plane()
{
    size_t i = 0;
    size_t N = mCV.size();
    while (i + 3 < N)
    {
        UpdatePoints_in_Plane(mCV[i], mCV[i + 1], mCV[i + 2], mCV[i + 3]);
        ++i;
    }
}

void ClassA_Bezier::TestSmooth(ONX_Model* model)
{
    for (int times = 0; times < 1; times++)
    {
        int N = 6;
        std::vector <ON_3dPoint> vp;
        vp.push_back(ON_3dPoint(0, 0, 0));
        vp.push_back(ON_3dPoint(2.5, 0, 0));
        vp.push_back(ON_3dVector(4, 1, 0));
        vp.push_back(ON_3dVector(5.7, 2.2, 0));
        vp.push_back(ON_3dVector(7, 4, 0));
        ClassA_Bezier class_A_bz(vp);
        for (int i = 0; i < 10; ++i)
        {
            class_A_bz.Elevate();
            int max_iter = 50;
            while ( max_iter > 0)
            {
                class_A_bz.Iteration_in_Plane();
                max_iter--;
            }
        }
        ON_NurbsCurve onc = class_A_bz.GetBezier();
        const int layer_index = model->AddLayer((std::wstring(L"Random_Test") + std::to_wstring(times)).c_str(), ON_Color::SaturatedGold);
        ChiralityAddNurbsCurve(model, onc, std::wstring(L"Random_Test") + std::to_wstring(times), layer_index);
        ChiralityAddLines(model, vp, std::wstring(L"Random_Test_Polygon") + std::to_wstring(times), layer_index);
        ChiralityDebugforR(onc, std::string("Random_Test") + std::to_string(times));
    }
}

void ClassA_Bezier::Elevate()
{
    size_t num_v = mCV.size();
    std::vector<ON_3dPoint> newpolygon;
    newpolygon.push_back(mCV[0]);
    for (size_t i = 1; i < num_v; ++i)
    {
        newpolygon.push_back(mCV[i - 1] * (i * 1.0 / num_v) + mCV[i] * (1.0 - i * 1.0 / num_v));
    }
    newpolygon.push_back(mCV.back());
    mCV = newpolygon;
}

ON_BezierCurve ClassA_Bezier::GetBezier() const
{
    ON_BezierCurve obc(3, false, mCV.size());
    for (size_t i = 0; i < mCV.size(); ++i)
    {
        obc.SetCV(int(i), mCV[i]);
    }
    return obc;
}
