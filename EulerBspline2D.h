#ifndef EULERBSPLINE2D_H
#define EULERBSPLINE2D_H
#include "thirdparty/opennurbs/opennurbs.h"

namespace EulerBspline2D
{
	bool EulerBsplineSpiralCheck(const ON_NurbsCurve &onc);
	void SmoothingBsplineControlPolygon(ON_NurbsCurve &onc, ON_3dPoint Pa, ON_3dPoint Pb, ON_3dVector Ta, ON_3dVector Tb);
	void EulerBsplineSpiralInterpolation(ON_NurbsCurve &onc, int max_vtx_num = 50);
	void EulerBsplineInterpolation_to_fixed_CVCount(ON_NurbsCurve &onc, int cv_count);
	void EulerBsplineTest(ONX_Model *model);
	void Smoothing3DControlPolygon(ON_NurbsCurve &onc, ON_3dPoint Pa, ON_3dPoint Pb, ON_3dVector Ta, ON_3dVector Tb);
	ON_NurbsCurve GenerateSmoothingCorner(ON_3dPoint start,ON_3dPoint corner,ON_3dPoint end);
	void SmoothCornerTest(ONX_Model* model);
}

#endif
