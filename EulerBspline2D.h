#ifndef EULERBSPLINE2D_H
#define EULERBSPLINE2D_H
#include "thirdparty/opennurbs/opennurbs.h"

namespace EulerBspline2D
{
	bool EulerBsplineSpiralCheck(const ON_NurbsCurve& onc);
	void SmoothingBsplineControlPolygon(ON_NurbsCurve& onc, ON_3dPoint Pa, ON_3dPoint Pb, ON_3dVector Ta, ON_3dVector Tb);
	void EulerBsplineSpiralInterpolation(ON_NurbsCurve& onc, int max_vtx_num = 50);
	void EulerBsplineTest(ONX_Model* model);
	void Smoothing3DControlPolygon(ON_NurbsCurve& onc, ON_3dPoint Pa, ON_3dPoint Pb, ON_3dVector Ta, ON_3dVector Tb);
}

#endif
