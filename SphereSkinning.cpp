#include "SphereSkinning.h"
#include "ChiralityMathTools.h"

void SphereSkinning::GetCirclesToInterpolate() const
{
    assert(spheres_.size() > 2);
    circles_.clear();
    int qwe[2] = {0, spheres_.size() - 1};
    int asd[2] = {1, spheres_.size() - 2};
    ON_Circle boundary[2];
    for (int i = 0; i < 2; ++i)
    {
        const ON_Sphere &c0 = spheres_[qwe[i]];
        const ON_Sphere &c1 = spheres_[asd[i]];
        double r0 = c0.Radius();
        double r1 = c1.Radius();
        ON_3dPoint o0 = c0.Center();
        ON_3dPoint o1 = c1.Center();
        double D = o0.DistanceTo(o1);
        assert(D > 1e-6);
        ON_3dVector o0_to_o1 = (o1 - o0) / D;
        ON_Plane plane(o0, o0_to_o1);
        ON_3dPoint new_center = o0 + o0_to_o1 * (r0 * (r0 - r1) / D);
        double new_r = sqrt(r0 * r0 - std::pow(r0 * (r0 - r1) / D, 2));
        boundary[i] = ON_Circle(plane, new_center, new_r);
    }
    circles_.push_back(boundary[0]);
    for (int i = 1; i < spheres_.size() - 1; ++i)
    {
        ON_Plane plane(spheres_[i - 1].Center(), spheres_[i].Center(), spheres_[i + 1].Center());
        ON_Circle c1(plane, spheres_[i - 1].Center(), spheres_[i - 1].Radius());
        ON_Circle c2(plane, spheres_[i].Center(), spheres_[i].Radius());
        ON_Circle c3(plane, spheres_[i + 1].Center(), spheres_[i + 1].Radius());
        ON_Circle l_c, r_c;
        ON_3dPoint l_p, r_p;
        int is_has_value = 0;
        std::vector<ON_Circle> external_circle = ChiralityMath::TangentToCircle(c1, c2, c3, 0, 0, 0);
        std::vector<ON_Circle> internal_circle = ChiralityMath::TangentToCircle(c1, c2, c3, 1, 1, 1);
        if (external_circle.size() == 2)
        {
            l_c = external_circle[0];
            r_c = external_circle[1];
            l_p = c2.Center() + (l_c.Center() - c2.Center()) / (l_c.Center() - c2.Center()).Length() * c2.Radius();
            r_p = c2.Center() + (r_c.Center() - c2.Center()) / (r_c.Center() - c2.Center()).Length() * c2.Radius();
            is_has_value++;
        }
        if (internal_circle.size() == 2)
        {
            l_c = internal_circle[0];
            r_c = internal_circle[1];
            assert((l_c.Center() - c2.Center()).Length() > 1e-6);
            assert((r_c.Center() - c2.Center()).Length() > 1e-6);
            int l_sign = l_c.Radius() > c2.Radius() ? -1 : 1;
            int r_sign = r_c.Radius() > c2.Radius() ? -1 : 1;
            l_p = c2.Center() + (l_c.Center() - c2.Center()) / (l_c.Center() - c2.Center()).Length() * c2.Radius() * l_sign;
            r_p = c2.Center() + (r_c.Center() - c2.Center()) / (r_c.Center() - c2.Center()).Length() * c2.Radius() * r_sign;
            is_has_value++;
        }
        if (internal_circle.size() == 1 && external_circle.size() == 1)
        {
            l_c = internal_circle[0];
            r_c = external_circle[0];
            assert((l_c.Center() - c2.Center()).Length() > 1e-6);
            int l_sign = l_c.Radius() > c2.Radius() ? -1 : 1;
            l_p = c2.Center() + (l_c.Center() - c2.Center()) / (l_c.Center() - c2.Center()).Length() * c2.Radius() * l_sign;
            r_p = c2.Center() + (r_c.Center() - c2.Center()) / (r_c.Center() - c2.Center()).Length() * c2.Radius();
            is_has_value++;
        }
        assert(is_has_value == 1);
        ON_3dPoint new_circle_center = (l_p + r_p) / 2;
        ON_3dVector new_circle_normal = ON_3dVector::CrossProduct(plane.Normal(), r_p - l_p);
        double R = l_p.DistanceTo(r_p) / 2;
        ON_Plane new_circle_plane(new_circle_center, new_circle_normal);
        circles_.push_back(ON_Circle(new_circle_plane, new_circle_center, R));
    }
    circles_.push_back(boundary[1]);
}

void SphereSkinning::GetConicVertex() const
{
    conic_vertices_.clear();
    int N = spheres_.size();
    for (int i = 0; i < N; ++i)
    {
        ON_3dPoint SO = spheres_[i].Center();
        ON_3dPoint CO = circles_[i].Center();
        double d = SO.DistanceTo(CO);
        if (d < 1e-6)
        {
            conic_vertices_.push_back(ON_4dPoint(0, 0, 0, 0));
        }
        else
        {
            double L = std::pow(spheres_[i].Radius(), 2) / d;
            ON_3dVector Dir = (CO - SO) / d;
            ON_3dPoint vertex = SO + Dir * L;
            conic_vertices_.push_back(ON_4dPoint(vertex.x, vertex.y, vertex.z, 1.0));
        }
    }
}

std::vector<ON_NurbsSurface> SphereSkinning::Skinning() const
{
    std::vector<ParameterCurve> v_pc;
    std::vector<ChiralityMath::VectorField> v_vf;
    GetParamCurveAndVectorField(v_pc, v_vf);
    std::vector<ON_NurbsSurface> v_ons = ChiralityMath::WirePlasticCurvedSurface(v_pc, v_vf, 10);
    return v_ons;
}

void SphereSkinning::GetParamCurveAndVectorField(std::vector<ParameterCurve>& v_pc, std::vector<ChiralityMath::VectorField>& v_vf) const
{
    std::vector<ON_3dVector> v_OW;
    int n = spheres_.size();
    for (int i = 0; i < n; ++i)
    {
        ON_3dPoint O = circles_[i].Center();
        ON_4dPoint W = conic_vertices_[i];
        if (W.w == 1.0)
        {
            ON_3dVector v = ON_3dVector(W) - O;
            v.Unitize();
            v_OW.push_back(v);
        }
        else
        {
            ON_3dVector v = i < n - 1 ? spheres_[i + 1].Center() - spheres_[i].Center() : spheres_[i].Center() - spheres_[i - 1].Center();
            v.Unitize();
            v_OW.push_back(v);
        }
    }
    int end = 1;
    ON_3dVector e;
    while (end != 0)
    {
        end = 0;
        e = ChiralityMath::GetRandomPoint(ON_3dPoint::Origin, 1.0, 2.0);
        for (const ON_3dVector& v : v_OW)
        {
            if (v.IsParallelTo(e) || (-v).IsParallelTo(e))
            {
                ++end;
            }
        }
    }
    e.Unitize();

    std::vector<int> Sign;
    for (int i = 0; i < n; ++i)
    {
        double product;
        if (i < n - 1)
        {
            product = ON_3dVector::DotProduct(v_OW[i], spheres_[i + 1].Center() - spheres_[i].Center());
        }
        else
        {
            product = ON_3dVector::DotProduct(v_OW[i], spheres_[i].Center() - spheres_[i - 1].Center());
        }
        Sign.push_back(product > 0 ? 1 : -1);
    }

    for (int i = 0; i < n; ++i)
    {
        ON_3dVector v = ON_3dVector::CrossProduct(e, v_OW[i]) * Sign[i];
        v.Unitize();
        ON_3dPoint StartPoint = circles_[i].Center() + circles_[i].Radius() * v;
        ON_3dVector Normal = Sign[i] * v_OW[i];
        ON_3dPoint center = circles_[i].Center();
        auto pos = [StartPoint,Normal,center](double t)->FrenetFrame {
            ON_3dPoint p = StartPoint;
            p.Rotate(t, Normal, center);
            ON_3dVector beta = center - p;
            ON_3dVector alpha = ON_3dVector::CrossProduct(beta, Normal);
            return FrenetFrame(p, alpha, beta);
        };
        double range[2] = { 0.0,2 * PI };
        v_pc.push_back(ParameterCurve(pos, range));

        if (conic_vertices_[i].w == 1.0)
        {
            ON_3dVector Tan = (ON_3dVector(conic_vertices_[i]) - StartPoint) * Sign[i];
            Tan.Unitize();
            auto vec = [Tan,Normal](double t)->ON_3dVector {
                ON_3dVector v = Tan;
                v.Rotate(t, Normal);
                return v;
            };
            v_vf.push_back(vec);
        }
        else
        {
            auto vec = [Normal](double t)->ON_3dVector {
                return Normal;
            };
            v_vf.push_back(vec);
        }
    }
}
