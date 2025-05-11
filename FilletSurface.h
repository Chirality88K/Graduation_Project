#ifndef FILLETSURFACE_H
#define FILLETSURFACE_H
#include "thirdparty/opennurbs/opennurbs.h"
#include <vector>

class VectorField1d;

class FilletSurface :
    public ON_Mesh
{
private:
    ON_NurbsSurface m_Surf0;
    ON_NurbsSurface m_Surf1;
    double InsectCurve0Par;
    double InsectCurve1Par;
    VectorField1d* Rail0VF;
    VectorField1d* Rail1VF;
    double RailCurvePar0;
    double RailCurvePar1;
public:
    enum InsectCurve0Knot { u0, v0 } m_par0;
    enum InsectCurve1Knot { u1, v1 } m_par1;
    //ON_BezierCurve FillCurve0;
    //ON_BezierCurve FillCurve1;
    std::vector <ON_BezierCurve> FillCurveArray;
    std::vector <double> Knot_on_Rail0;
    std::vector <double> Knot_on_Rail1;
public:
    ~FilletSurface();
    void SetSurface(int i, const ON_NurbsSurface& nurbs);
    void SetParameter(InsectCurve0Knot k0, double par0, InsectCurve1Knot k1, double par1);
    void SetRailCurve(double offset0, double offset1);
    void ReverseRailCurve(int i);
    void ReverseVector(int i);
    void ComputeFillCurve(int num_knot, double* knot_on_rail0, double* knot_on_rail1);
    ON_3dPoint Point_at(int sur_index, double u, double v);
    ON_3dVector GetNormal(double u, double v);
    void Add_FillCurves(ONX_Model* model, const wchar_t* name, ON_Color color = ON_Color::Black);
    void Add_FilletSurface(ONX_Model* model, const wchar_t* name, int FillCurveSample_num, int EveryPatchSample_num, ON_Color color = ON_Color::Black);
    static void SetRailUniformKnot(int num_knot, double* knot_on_rail0, double* knot_on_rail1);
    void Skinning(ON_NurbsSurface& ons, int num_fill_curve, double* u_knot);
};

class VectorField1d
{
private:
    ON_NurbsSurface* m_source_surface;
    bool is_first_const;//true:��һ�������ǳ���
    bool is_curve_reverse;
    bool is_vector_reverse;
    double par;
public:
    VectorField1d(ON_NurbsSurface*, bool, double);
    ~VectorField1d();
    void ReverseCurve();
    void ReverseVector();
    void GetAll(double t, ON_3dPoint& p, ON_3dVector& n, ON_3dVector& v);//�㣬��������
    ON_3dVector MixedDer(double t);
};
#endif