#ifndef CORNU_SPIRAL_H
#define CORNU_SPIRAL_H
#include "thirdparty/opennurbs/opennurbs.h"
#include "ParameterSurface.h"

class Cornu_Spiral
{
private:
	ON_2dPoint m_P0;
	ON_2dVector m_T0;
	ON_2dVector m_N0;
	ON_2dPoint m_PS;
	ON_2dPoint m_PE;
	double m_a;
	double m_theta0;
	double m_theta1;
	bool m_is_mirror = false;
	bool m_is_reverse = false;

private:
	
	
	static double Compute_f(double theta, double phi1, double phi2);
	static double Compute_g(double omega, double phi1, double phi2);
	static double MidSection_f(double phi1, double phi2);
	static double MidSection_g(double phi1, double phi2);
	bool Check_And_Init(ON_2dPoint& ps, ON_2dPoint& pe, ON_2dVector& vs, ON_2dVector& ve);
	ON_2dPoint GetValue(double theta);
	ON_2dVector GetTangent(double theta);
	double GetSignedArcLength(double theta) const;
	double GetSignedCurvature(double theta) const;
	ParameterCurve GetParameterCurve() const;

public:
	
	Cornu_Spiral(ON_2dPoint, ON_2dPoint, ON_2dVector, ON_2dVector);
	double GetTotalLength() const 
	{
		return GetSignedArcLength(m_theta1) - GetSignedArcLength(m_theta0);
	}
	static ParameterCurve GetParamCornuSpiral(ON_2dPoint ps, ON_2dPoint pe, ON_2dVector vs, ON_2dVector ve);
	static ON_NurbsCurve GetNurbs(ON_2dPoint ps, ON_2dPoint pe, ON_2dVector vs, ON_2dVector ve);
	static double Compute_C(double theta);
	static double Compute_S(double theta);
	void Add_to_Model(ONX_Model *model, const wchar_t *name, ON_Color color = ON_Color::Black);
	void Add_Nurbs_to_Model(ONX_Model *model, const wchar_t *name, ON_Color color = ON_Color::Black);
	void Raise_to_3D(double zheight, ONX_Model *model, const wchar_t *name, ON_Color color = ON_Color::Black);
	static void Cornu_test(ONX_Model *model);
};
#endif