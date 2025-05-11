#ifndef CONE_SURFACE_H
#define CONE_SURFACE_H
#include "thirdparty/opennurbs/opennurbs.h"
class Cone_Surface :
    public ON_NurbsSurface
{
public:
    Cone_Surface(double angle, double r1, double r2, double theta);
    void Add_Cone_to_Model(ONX_Model* model, const wchar_t* name, ON_Color color = ON_Color::Black);
};
#endif