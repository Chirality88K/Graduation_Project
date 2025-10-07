#include "FixAxisBezier3D.h"
#include "write3dm.h"
#include "EulerBezier2D.h"
#include <algorithm>
#include "ChiralityMathTools.h"
#include <sstream>
#include <iomanip>

ON_BezierCurve FixAxisBezier3D::Project_to_Plane(ON_3dVector T) const
{
    int N = mCV.size();
    T.Unitize();
    ON_BezierCurve obc(3, false, N);
    for (int i = 0; i < N; ++i)
    {
        ON_3dVector v = mCV[i] - mCV[0];
        ON_3dVector project_v_to_plane = v - ON_3dVector::DotProduct(v, T) * T;
        ON_3dPoint new_p = mCV[0] + project_v_to_plane;
        obc.SetCV(i, new_p);
    }
    return obc;
}

ON_BezierCurve FixAxisBezier3D::Go_Back_To_Space(const ON_BezierCurve& obc_plane, ON_3dVector T, double* error) const
{
    T.Unitize();
    double sin_alpha = ON_3dVector::DotProduct(T, mVS);//alpha is the angle of mVS and plane with the normal T
    double tan_alpha = sin_alpha / sqrt(1 - sin_alpha * sin_alpha);
    ON_3dPoint p0, p1;
    obc_plane.GetCV(0, p0);
    obc_plane.GetCV(1, p1);
    double H = ON_3dVector::DotProduct(mPE - mPS, T);
    double h = H / double(obc_plane.Order() - 1);
    ON_BezierCurve obc_space(3, false, obc_plane.Order());
    for (int i = 0; i < obc_plane.Order(); ++i)
    {
        obc_plane.GetCV(i, p0);
        p0 = p0 + T * h * double(i);
        obc_space.SetCV(i, p0);
    }
    if (error != nullptr) 
    {
        int N = obc_space.Order();
        ON_3dVector vs, ve;
        obc_space.GetCV(0, p0);
        obc_space.GetCV(1, p1);
        vs = p1 - p0;
        vs.Unitize();
        obc_space.GetCV(N - 2, p0);
        obc_space.GetCV(N - 1, p1);
        ve = p1 - p0;
        ve.Unitize();
        *error = (vs - mVS).Length() + (ve - mVE).Length();
    }
    return obc_space;
}

static double FindMinValue(const std::function<double(double)>& F, double x, double step)
{
    double ratio = 0.7;
    double value0 = F(x);
    for (int i = 0; i < 70; ++i)
    {
        double value1 = F(x + step);
        if (value1 < value0)
        {
            x += step;
            if (value1 < 0.01)
            {
                step *= ratio;
            }
            value0 = value1;
        }
        else
        {
            double value2 = F(x - step);
            if (value2 < value0)
            {
                x -= step;
                if (value2 < 0.01)
                {
                    step *= -ratio;
                }
                value0 = value2;
            }
            else {
                step *= ratio;
            }
        }
        if (value0 < 1e-6) {
            return value0;
        }
    }
    CHIRALITY_WARN(std::string("precision may not be enough!!!"));
    return value0;
}

std::pair<double, double> FixAxisBezier3D::IterationForFirstTime(double* error) const
{
    ON_3dVector X = mVS + mVE;
    ON_3dVector Y = ON_3dVector::CrossProduct(mVS, mVE);
    X.Unitize();
    Y.Unitize();
    double t_theta_start = 0.0;
    double t_theta_end = PI;
    int N = 30;
    std::vector<std::pair<double, double>> error_data;
    for (int i = 1; i < N; ++i)
    {
        double t_theta = t_theta_start * (1 - double(i) / double(N)) + t_theta_end * double(i) / double(N);
        ON_3dVector T = X * cos(t_theta) + Y * sin(t_theta);
        ON_BezierCurve plane_obc = Project_to_Plane(T);
        ON_Xform rotate, inv_rotate;
        rotate.Rotation(T, ON_3dVector::ZAxis, ON_3dPoint::Origin);
        inv_rotate = rotate.Inverse();
        plane_obc.Transform(rotate);
        EulerBezier2D::EulerBezierWeakInterpolation(&plane_obc, 30);
        ON_BezierCurve plane_obc_prepare_to_space(3, false, plane_obc.Order());
        for (int i = 0; i < plane_obc.Order(); ++i)
        {
            ON_3dPoint p;
            plane_obc.GetCV(i, p);
            plane_obc_prepare_to_space.SetCV(i, ON_3dPoint(p.x, p.y, 0.0));
        }
        plane_obc_prepare_to_space.Transform(inv_rotate);
        double error = 100.0;
        ON_BezierCurve space_obc = Go_Back_To_Space(plane_obc_prepare_to_space, T, &error);
        error_data.push_back(std::make_pair(t_theta, error));
    }
    std::sort(error_data.begin(), error_data.end(),
        [](const std::pair<double, double>& pair1, const std::pair<double, double>& pair2) {
            return pair1.second < pair2.second;
        });

    double new_theta = error_data[0].first;
    if (error != nullptr)
    {
        *error = error_data[0].second;
    }
    double delta_theta = (t_theta_end - t_theta_start) / N;
    t_theta_start = new_theta - delta_theta;
    t_theta_end = new_theta + delta_theta;
    return std::make_pair(t_theta_start, t_theta_end);
}

ON_BezierCurve FixAxisBezier3D::BisectionIteration(double start, double end, double* in_error, ON_3dVector& TT) const
{
    double min_error = *in_error;
    ON_BezierCurve final_result;
    ON_3dVector X = mVS + mVE;
    ON_3dVector Y = ON_3dVector::CrossProduct(mVS, mVE);
    X.Unitize();
    Y.Unitize();
    ON_3dVector Z = ON_3dVector::CrossProduct(X, Y);
    ON_3dVector pspe = mPE - mPS;
    ON_3dVector vp = pspe / pspe.Unitize();
    ON_3dVector pro_vp = vp - ON_3dVector::DotProduct(vp, Z) * Z;
    double init_cos = ON_3dVector::DotProduct(pro_vp, X) / pro_vp.Length();
    double init_tan = -ON_3dVector::DotProduct(X, pspe) / ON_3dVector::DotProduct(Y, pspe);
    double init_angle = atan(init_tan);
    if (init_angle < 0) {
        init_angle += PI;
    }
    init_angle = acos(init_cos);
    double init_step = (std::min)(init_angle, PI - init_angle) / 2;
    auto F = [this,&final_result,&X,&Y,&min_error,&TT](double theta)->double {
        ON_3dVector T = X * cos(theta) + Y * sin(theta);
        ON_BezierCurve plane_obc = Project_to_Plane(T);
        ON_Xform rotate;
        rotate.Rotation(T, ON_3dVector::ZAxis, ON_3dPoint::Origin);
        ON_Xform translation;
        translation.Translation(-ON_3dVector(mPS));
        ON_Xform trans_and_rotate = rotate * translation;
        ON_Xform reverse = trans_and_rotate.Inverse();
        ON_3dPoint plane_ps = mPS;
        ON_3dPoint plane_pe = mPE - ON_3dVector::DotProduct((mPE - mPS), T) * T;
        ON_3dVector plane_vs = mVS - ON_3dVector::DotProduct(mVS, T) * T;
        ON_3dVector plane_ve = mVE - ON_3dVector::DotProduct(mVE, T) * T;
        plane_vs.Unitize();
        plane_ve.Unitize();
        plane_ps.Transform(trans_and_rotate);
        plane_pe.Transform(trans_and_rotate);
        plane_vs.Transform(trans_and_rotate);
        plane_ve.Transform(trans_and_rotate);
        ON_BezierCurve renew_plane_obc = ChiralityMath::BezierG1_xOy(plane_ps,plane_pe, plane_vs, plane_ve);
        for (int i = 0; i < 20; ++i)
        {
            EulerBezier2D::Elevate(&renew_plane_obc);  
        }
        EulerBezier2D::SmoothingBezierControlPolygon(&renew_plane_obc);
        ON_BezierCurve plane_obc_prepare_to_space(3, false, renew_plane_obc.Order());
        for (int i = 0; i < renew_plane_obc.Order(); ++i)
        {
            ON_3dPoint p;
            renew_plane_obc.GetCV(i, p);
            plane_obc_prepare_to_space.SetCV(i, ON_3dPoint(p.x, p.y, 0.0));
        }
        plane_obc_prepare_to_space.Transform(reverse);
        double error = 100.0;
        ON_BezierCurve space_obc = Go_Back_To_Space(plane_obc_prepare_to_space, T, &error);
        if (error < min_error)
        {
            min_error = error;
            final_result = space_obc;
        }
        TT = T;
        return error;
    };
    //auto pair = IterationForFirstTime(in_error);
    //min_error = FindMinValue(F, (pair.first + pair.second) / 2, (pair.first - pair.second) / 2);
    min_error = FindMinValue(F, PI / 2, PI / 12);
    *in_error = min_error;
    std::cout << min_error << "\n";
    return final_result;
}

ON_BezierCurve FixAxisBezier3D::Interpolate(ON_3dPoint ps, ON_3dPoint pe, ON_3dPoint vs, ON_3dPoint ve)
{
    FixAxisBezier3D fix_bezier({ ps,ps + vs,pe - ve,pe });
    double error = 100;
    ON_3dVector T;
    ON_BezierCurve obc = fix_bezier.BisectionIteration(0, 0, &error, T);
    for (int i = 0; i < 0; ++i)
    {
        ChiralityMath::Elevate(obc);
        SmoothingWithFixedT(obc, T);
    }
    ChiralityDebugforR(obc, std::string("Random_Test") + "error" + doubleToScientificString(error));
    return obc;
}

void FixAxisBezier3D::Test(ONX_Model* model)
{
    for (int times = 0; times < 5; times++)
    {
        int N = 6;
        std::vector <ON_3dPoint> vp;
        vp.push_back(ON_3dPoint(0, 0, 0));
        vp.push_back(ON_3dPoint(2.5, 0, 0));
        for (int i = 0; i < N; ++i)
        {
            ON_3dPoint p;
            while (true)
            {
                p = ChiralityMath::GetRandomPoint(vp[i + 1], 2.0, 3.0);
                double product = ON_3dVector::DotProduct((p - vp[i + 1]), vp[i + 1] - vp[i]);
                double cos_theta = product / (p - vp[i + 1]).Length() / (vp[i + 1] - vp[i]).Length();
                if (cos_theta < 0.99 && cos_theta>0.1)
                {
                    break;
                }
            }
            vp.push_back(p);
        }
        FixAxisBezier3D test_fix_axis_bez(vp);
        double error = 100;
        //auto pair_theta = test_fix_axis_bez.IterationForFirstTime(&error);
        ON_3dVector T;
        ON_BezierCurve obc = test_fix_axis_bez.BisectionIteration(0, 0, &error, T);
        const int layer_index = model->AddLayer((std::wstring(L"Random_Test") + std::to_wstring(times)).c_str(), ON_Color::SaturatedGold);
        ChiralityAddNurbsCurve(model, obc, std::wstring(L"Random_Test") + std::to_wstring(times) + L"; error: " + StringToWString(doubleToScientificString(error)), layer_index);
        ChiralityAddLines(model, vp, std::wstring(L"Random_Test_Polygon") + std::to_wstring(times), layer_index);
        ChiralityDebugforR(obc, std::string("Random_Test") + std::to_string(times) + "error" + doubleToScientificString(error));
    }
}

void FixAxisBezier3D::GenerateDNA(ONX_Model* model)
{
    const int index1 = model->AddLayer(L"test1", ON_Color::SaturatedRed);
    const int index2 = model->AddLayer(L"test2", ON_Color::SaturatedBlue);
    double R = 1.0;
    double ratio = 2.0;
    double begin = 0.0;
    double end = 4 * PI;
    auto circular_helix = [a = R, b = ratio](double theta) -> FrenetFrame
    {
        ON_3dPoint p(a * cos(theta), a * sin(theta), b * theta);
        ON_3dVector der(-a * sin(theta), a * cos(theta), b);
        ON_3dVector derder(-a * cos(theta), -a * sin(theta), 0);
        ON_3dVector N = ON_3dVector::CrossProduct(der, derder);
        ON_3dVector B = ON_3dVector::CrossProduct(N, der);
        der.Unitize();
        B.Unitize();
        return FrenetFrame(p, der, B);
    };

    const int sample_cnt = 4;
    for (int i = 0; i < sample_cnt; ++i)
    {
        double theta = begin * (1 - double(i) / double(sample_cnt)) + end * double(i) / double(sample_cnt);
        FrenetFrame f1 = circular_helix(theta);
        FrenetFrame f2 = circular_helix(theta + PI);
        std::vector<ON_3dPoint> vp = { f1.GetPos(),f1.GetPos() + f1.GetAlpha(),f2.GetPos() - f2.GetAlpha(),f2.GetPos() };
        FixAxisBezier3D cab(vp);
        double error = 100;
        //auto pair_theta = test_fix_axis_bez.IterationForFirstTime(&error);
        ON_3dVector T;
        ON_BezierCurve obc = cab.BisectionIteration(0, 0, &error, T);
        ChiralityAddNurbsCurve(model, obc, std::wstring(L"DNA") + std::to_wstring(sample_cnt) + L"; error: " + StringToWString(doubleToScientificString(error)), index1);
        ChiralityDebugforR(obc, std::string("DNA") + std::to_string(sample_cnt) + "error" + doubleToScientificString(error));
        obc.Rotate(0, -1, ON_3dVector::ZAxis, ON_3dPoint::Origin);
        ChiralityAddNurbsCurve(model, obc, std::wstring(L"DNA") + std::to_wstring(sample_cnt) + L"; error: " + StringToWString(doubleToScientificString(error)), index2);
    }
}

static double ComputeMean(std::vector<double>& v)
{
    double avg = 0.0;
    for (double l : v)
    {
        avg += l;
    }
    avg /= v.size();
    return avg;
}

static double ComputeCV(std::vector<double>& v)
{
    double mean = ComputeMean(v);
    double cv = 0.0;
    for (double l : v)
    {
        cv += (l - mean) * (l - mean);
    }
    return sqrt(cv / v.size()) / mean;
}

void FixAxisBezier3D::SmoothingWithFixedT(ON_BezierCurve& obc, ON_3dVector T)
{
    T.Unitize();
    std::vector<double> angle = ComputeRotatingAngles(obc, T);
    double angle_test = 0.0;
    for (size_t i = 1; i < angle.size() - 1; ++i)
    {
        angle_test = (std::max)(angle_test, abs(angle[i] * 2 - angle[i - 1] - angle[i + 1]));
    }
    std::vector<double> length = ComputeLength(obc);
    double length_test = ComputeCV(length);
    int s_count = 0;
    int max_count = 1000;
    ON_3dVector vs = obc.TangentAt(0.0);
    ON_3dVector ve = obc.TangentAt(1.0);
    ON_3dPoint ps = obc.PointAt(0.0);
    ON_3dPoint pe = obc.PointAt(1.0);
    while ((angle_test > 1e-4 || length_test > 0.01) && s_count < max_count)
    {
        ON_3dPoint p1, p3;
        for (int i = 2; i < obc.CVCount() - 2; ++i)
        {
            obc.GetCV(i - 1, p1);
            obc.GetCV(i + 1, p3);
            ON_Xform RT;
            RT.Rotation((angle[i - 2] + angle[i - 1] + angle[i]) / 3.0, T, ON_3dPoint::Origin);
            ON_Xform I = ON_Xform::IdentityTransformation;
            ON_3dPoint p2 = ((RT + I).Inverse()) * (RT * ON_3dVector(p1) + ON_3dVector(p3));
            obc.SetCV(i, p2);
        }
        length = ComputeLength(obc);
        double avg_length = ComputeMean(length);
        obc.SetCV(1, ps + avg_length * vs);
        obc.SetCV(obc.CVCount() - 2, pe - avg_length * ve);
        length = ComputeLength(obc);
        length_test = ComputeCV(length);
        angle = ComputeRotatingAngles(obc, T);
        angle_test = 0.0;
        for (size_t i = 1; i < angle.size() - 1; ++i)
        {
            angle_test = (std::max)(angle_test, abs(angle[i] * 2 - angle[i - 1] - angle[i + 1]));
        }
        ++s_count;
    }
}

std::vector<double> FixAxisBezier3D::ComputeLength(const ON_BezierCurve& c)
{
    std::vector<double> length;
    ON_3dPoint p1;
    ON_3dPoint p2;
    for (int i = 1; i < c.CVCount(); ++i)
    {
        c.GetCV(i - 1, p1);
        c.GetCV(i, p2);
        length.push_back(p1.DistanceTo(p2));
    }
    return length;
}

std::vector<double> FixAxisBezier3D::ComputeRotatingAngles(const ON_BezierCurve& c, ON_3dVector T)
{
    T.Unitize();
    std::vector<double> angles;
    ON_3dPoint p1, p2, p3;
    ON_3dVector v1, v2;
    for (int i = 1; i < c.CVCount() - 1; ++i)
    {
        c.GetCV(i - 1, p1); c.GetCV(i, p2); c.GetCV(i + 1, p3);
        v1 = p2 - p1; v2 = p3 - p2;
        v1.Unitize(); v2.Unitize();
        v1 = v1 - ON_3dVector::DotProduct(v1, T) * T;
        v2 = v2 - ON_3dVector::DotProduct(v2, T) * T;
        double product = ON_3dVector::DotProduct(v1, v2) / v1.Length() / v2.Length();
        product = (std::min)(1.0, product);
        product = (std::max)(-1.0, product);
        if (ON_3dVector::DotProduct(T, ON_3dVector::CrossProduct(v1, v2)) >= 0)
        {
            angles.push_back(acos(product));
        }
        else {
            angles.push_back(-acos(product));
        }
    }
    return angles;
}

bool FixAxisBezier3D::CheckCurvatureMono(const ON_BezierCurve& obc)
{
    return false;
}
