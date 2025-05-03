#ifndef EULERBEZIER3D_H
#define EULERBEZIER3D_H
#include "thirdparty/opennurbs/opennurbs.h"

class EulerBezier3D :public ON_BezierCurve
{
public:
	void Smoothing();
};

void EulerBezier3D::Smoothing()
{
	
}



#endif // !EULERBEIZER3D_H
