#pragma once
#include "opennurbs.h"

class Cornu_Spiral
{
private:
	ON_2dPoint m_P0;
	ON_2dVector m_T0;
	ON_2dVector m_N0;
	double m_a;
	double m_theta0;
	double m_theta1;
	static double Compute_C(double theta);
	static double Compute_S(double theta);
	static double Compute_f(double theta, double phi1, double phi2);
	static double Compute_g(double omega, double phi1, double phi2);
	static double MidSection_f(double phi1, double phi2);
	static double MidSection_g(double phi1, double phi2);
	ON_2dPoint GetValue(double theta);
	ON_2dVector GetTangent(double theta);
	double GetSignedArcLength(double theta);
	double GetSignedCurvature(double theta);
public:
	Cornu_Spiral(ON_2dPoint, ON_2dPoint, ON_2dVector, ON_2dVector);
	void Add_to_Model(ONX_Model* model, const wchar_t* name, ON_Color color = ON_Color::Black);
	void Add_Nurbs_to_Model(ONX_Model* model, const wchar_t* name, ON_Color color = ON_Color::Black);
	void Raise_to_3D(double zheight, ONX_Model* model, const wchar_t* name, ON_Color color = ON_Color::Black);
	static void Cornu_test(ONX_Model* model);
};