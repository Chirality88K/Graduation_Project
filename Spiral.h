#ifndef SPIRAL_H
#define SPIRAL_H
#include "thirdparty/opennurbs/opennurbs.h"

class Spiral
{
private:
	std::vector<ON_3dPoint> m_all_points;
	std::vector<ON_3dVector> m_all_tan;
	std::vector<double> m_all_curvature;
	ON_NurbsCurve m_three_degree_spiral;
	ON_NurbsCurve m_five_degree_spiral;
	ON_NurbsCurve m_pieces_bezier;
	double m_Length;
	int m_N;
	ON_3dVector m_beta0;
	ON_3dVector m_betaN;

public:
	Spiral(ON_3dPoint start, ON_3dPoint end, ON_3dVector alpha0, ON_3dVector alphaL, double kappa0, double kappaL, double tau0, double tauL, double L, int n);
	void ComputeNurbs_Directly();
	void Add_spiral_to_model(ONX_Model *model, const wchar_t *name, ON_Color color = ON_Color::Black);
	void ComputePiecesBezier();
	void Add_pieces_to_model(ONX_Model *model, const wchar_t *name, ON_Color color = ON_Color::Black);
	void ComputeNurbs_Five_degrees();
	void Add_five_degree_spiral_to_model(ONX_Model *model, const wchar_t *name, ON_Color color = ON_Color::Black);
	void Get_five_degree_spline(ON_NurbsCurve &onc);
};

void ThreeDegreeBsplineInterplate_Tan(ON_NurbsCurve &onc, std::vector<ON_3dPoint> Q, double *knot, ON_3dVector v0, ON_3dVector vn);
void ThreeDegreeBsplineInterplate_Tan(ON_NurbsCurve &onc, std::vector<ON_3dPoint> Q, std::vector<double> knot, ON_3dVector v0, ON_3dVector vn);
#endif