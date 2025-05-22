#include "EulerBspline2D.h"
#include <vector>
#include <algorithm>
#include "write3dm.h"
#include "ChiralityMathTools.h"
using namespace std;
extern const double PI;

namespace EulerBspline2D
{
	bool EulerBsplineSpiralCheck(const ON_NurbsCurve &onc)
	{
		int n = onc.CVCount();
		if (n < 4)
		{
			return false;
		}
		// Compute Length
		std::vector<double> Length;
		ON_3dPoint p1, p2;
		for (int i = 0; i < n - 1; i++)
		{
			onc.GetCV(i, p1);
			onc.GetCV(i + 1, p2);
			Length.push_back((p1 - p2).Length());
		}
		// Check length is equal or not
		double sum = 0;
		for (int i = 0; i < n - 1; i++)
		{
			sum = sum + Length[i];
		}
		double avg = sum / (n - 1);
		sum = 0;
		for (int i = 0; i < n - 1; i++)
		{
			sum = sum + (Length[i] - avg) * (Length[i] - avg);
		}
		double s2 = sum / (n - 1);
		if (s2 > 0.1)
		{
			return false;
		}
		// Compute Angles
		std::vector<double> Angle;
		Angle.push_back(0);
		for (int i = 0; i < n - 2; i++)
		{
			onc.GetCV(i, p1);
			onc.GetCV(i + 1, p2);
			ON_3dVector v1 = p2 - p1;
			onc.GetCV(i + 2, p1);
			ON_3dVector v2 = p1 - p2;
			double theta = ON_3dVector::Angle(v1, v2);
			if (theta > PI / 2)
			{
				return false;
			}
			if (v1.x * v2.y - v1.y * v2.x < 0)
			{
				theta = -theta;
			}
			Angle.push_back(theta);
		}
		Angle.push_back(0);
		// Compute theta_DD
		for (int i = 2; i < n - 2; i++)
		{
			if (abs(Angle[i] * 2 - Angle[i - 1] - Angle[i + 1]) > 1e-6)
			{
				return false;
			}
		}
		// Compute delta_theta,s0,s1
		int m = n - 1;
		double Deltatheta = (Angle[m - 1] - Angle[1]) / (m - 2);
		double s0 = (m + 1) * sin(Angle[1]) + (m - 2) * sin(Angle[1] + Angle[2]) - 3 * (m - 1) * sin(Angle[1]) * cos(Angle[1]);
		double s1 = -(m + 1) * sin(Angle[n - 2]) - (m - 2) * sin(Angle[n - 3] + Angle[n - 2]) + 3 * (m - 1) * sin(Angle[n - 2]) * cos(Angle[n - 2]);
		if (Deltatheta * s0 >= 0 && s0 * s1 >= 0)
		{
			return true;
		}
		return false;
	}

	void SmoothingBsplineControlPolygon(ON_NurbsCurve &onc, ON_3dPoint Pa, ON_3dPoint Pb, ON_3dVector Ta, ON_3dVector Tb)
	{
		int n = onc.CVCount();
		int m = onc.CVCount() - 1;
		if (m <= 3)
		{
			return;
		}

		// ComputeAngle();
		std::vector<double> Angle;
		ON_3dPoint p1, p2;
		Angle.push_back(0);
		for (int i = 0; i < n - 2; i++)
		{
			onc.GetCV(i, p1);
			onc.GetCV(i + 1, p2);
			ON_3dVector v1 = p2 - p1;
			onc.GetCV(i + 2, p1);
			ON_3dVector v2 = p1 - p2;
			double theta = ON_3dVector::Angle(v1, v2);
			if (v1.x * v2.y - v1.y * v2.x < 0)
			{
				theta = -theta;
			}
			Angle.push_back(theta);
		}
		Angle.push_back(0);

		int max_count = 10000;
		for (int i = 0; i < n; i++)
		{
			if (Angle[i] > PI / 2 || Angle[i] < -PI / 2)
			{
				max_count = 2;
			}
		}
		int s_count = 1;
		double thetaDD = 10.0;
		double lengthavg = 0.0;
		onc.GetCV(0, p1);
		onc.GetCV(m, p2);
		double lengthbound = 2 * (p1 - p2).Length();

		while (s_count < max_count && thetaDD > 1e-6 && lengthavg < lengthbound)
		{
			vector<double> new_angle = Angle;

			for (int i = 2; i < m - 1; i++)
			{
				new_angle[i] = (new_angle[i - 1] + new_angle[i] + new_angle[i + 1]) / 3;
				onc.GetCV(i - 1, p1);
				onc.GetCV(i + 1, p2);
				ON_3dPoint re = (p2 - p1) / 2;
				re.Set(-re.y, re.x, 0);
				re = (p1 + p2) / 2 -
					 re * tan(new_angle[i] / 2);
				onc.SetCV(i, re);
			}

			onc.GetCV(2, p2);
			Ta.Unitize();
			Tb.Unitize();
			double la = ON_3dVector::DotProduct(p2 - Pa, Ta);
			if (la > 0)
			{
				onc.SetCV(0, p2 - 2 * la * Ta);
			}
			else
			{
				onc.SetCV(0, Pa + la * Ta);
			}
			ON_3dVector Na = Ta;
			Na.Rotate(1, 0, ON_3dVector(0, 0, 1));
			onc.GetCV(0, p1);
			double ha = 0.5 * ON_3dVector::DotProduct(Na, Pa - 0.5 * (p1 + p2));
			onc.SetCV(1, Pa + ha * Na);

			onc.GetCV(m - 2, p2);
			double lb = ON_3dVector::DotProduct(p2 - Pb, Tb);
			if (lb < 0)
			{
				onc.SetCV(m, p2 - 2 * lb * Tb);
			}
			else
			{
				onc.SetCV(m, Pb + lb * Tb);
			}
			ON_3dVector Nb = Tb;
			Nb.Rotate(1, 0, ON_3dVector(0, 0, 1));
			onc.GetCV(m, p1);
			double hb = 0.5 * ON_3dVector::DotProduct(Nb, Pb - 0.5 * (p1 + p2));
			onc.SetCV(m - 1, Pb + hb * Nb);

			// ComputeAngle() again!
			Angle.clear();
			Angle.push_back(0);
			for (int i = 0; i < n - 2; i++)
			{
				onc.GetCV(i, p1);
				onc.GetCV(i + 1, p2);
				ON_3dVector v1 = p2 - p1;
				onc.GetCV(i + 2, p1);
				ON_3dVector v2 = p1 - p2;
				double theta = ON_3dVector::Angle(v1, v2);
				if (v1.x * v2.y - v1.y * v2.x < 0)
				{
					theta = -theta;
				}
				Angle.push_back(theta);
			}
			Angle.push_back(0);

			thetaDD = 0;
			if (m > 3)
			{
				for (int i = 2; i < m - 1; i++)
				{
					thetaDD = max(thetaDD, abs(Angle[i] * 2 - Angle[i - 1] - Angle[i + 1]));
				}
			}
			s_count += 1;
		}
	}

	void EulerBsplineSpiralInterpolation(ON_NurbsCurve &onc, int max_vtx_num)
	{
		int v_num = onc.CVCount();
		while (v_num <= max_vtx_num && !EulerBsplineSpiralCheck(onc))
		{
			v_num = onc.CVCount();
			ON_3dPoint p1, p2;
			vector<ON_3dPoint> vp;
			onc.GetCV(0, p1);
			vp.push_back(p1);
			for (int i = 1; i < v_num; i++)
			{
				onc.GetCV(i - 1, p1);
				onc.GetCV(i, p2);
				vp.push_back(p1 * (double(i) / double(v_num)) + p2 * (1 - double(i) / double(v_num)));
			}
			onc.GetCV(v_num - 1, p1);
			vp.push_back(p1);
			double u0, u1;
			onc.GetDomain(&u0, &u1);
			onc.EvPoint(u0, p1);
			onc.EvPoint(u1, p2);
			ON_3dVector v1 = onc.TangentAt(u0);
			ON_3dVector v2 = onc.TangentAt(u1);
			onc.Create(2, false, onc.Order(), v_num + 1);

			int i = 0;
			for (const auto &it : vp)
			{
				onc.SetCV(i, it);
				i++;
			}
			SmoothingBsplineControlPolygon(onc, p1, p2, v1, v2);
			// onc.MakeClampedUniformKnotVector();
			/*
			const double* kn = onc.Knot();
			vector <double> seeknot;
			for (int i = 0; i < onc.KnotCount(); i++)
			{
				seeknot.push_back(*(kn + i));
			}
			*/
			for (int i = 0; i < onc.KnotCount(); i++)
			{
				onc.SetKnot(i, i + 1);
			}
		}
	}

	void EulerBsplineInterpolation_to_fixed_CVCount(ON_NurbsCurve &onc, int cv_count)
	{
		int v_num = onc.CVCount();
		while (v_num < cv_count)
		{
			ON_3dPoint p1, p2;
			vector<ON_3dPoint> vp;
			onc.GetCV(0, p1);
			vp.push_back(p1);
			for (int i = 1; i < v_num; i++)
			{
				onc.GetCV(i - 1, p1);
				onc.GetCV(i, p2);
				vp.push_back(p1 * (double(i) / double(v_num)) + p2 * (1 - double(i) / double(v_num)));
			}
			onc.GetCV(v_num - 1, p1);
			vp.push_back(p1);
			double u0, u1;
			onc.GetDomain(&u0, &u1);
			onc.EvPoint(u0, p1);
			onc.EvPoint(u1, p2);
			ON_3dVector v1 = onc.TangentAt(u0);
			ON_3dVector v2 = onc.TangentAt(u1);
			onc.Create(2, false, onc.Order(), v_num + 1);

			int i = 0;
			for (const auto &it : vp)
			{
				onc.SetCV(i, it);
				i++;
			}
			SmoothingBsplineControlPolygon(onc, p1, p2, v1, v2);
			for (int i = 0; i < onc.KnotCount(); i++)
			{
				onc.SetKnot(i, i + 1);
			}
			v_num = onc.CVCount();
		}
	}

	void EulerBsplineTest(ONX_Model *model)
	{
		ON_NurbsCurve *onc = new ON_NurbsCurve();
		*onc = ChiralityMath::UniformG1(ON_3dPoint::Origin, ON_3dPoint(10.0, 0.0, 0.0), ON_3dVector(3.0, -1.0, 0.0), ON_3dVector(1.0, 1.0, 0.0));

		const int layer_index = model->AddLayer(L"test_input_Bspline", ON_Color::SaturatedBlue);
		model->AddManagedModelGeometryComponent(onc, Internal_CreateManagedAttributes(layer_index, L"test_input_Bspline"));
		PrintPosAndTan(*onc, "origin_Bspline_pos&tan");
		PrintCurvature(*onc, "origin_Bspline_curvature");
		ON_NurbsCurve *another_onc = new ON_NurbsCurve(*onc);
		EulerBsplineSpiralInterpolation(*another_onc);
		const int layer_index2 = model->AddLayer(L"test_spiral_Bspline", ON_Color::SaturatedMagenta);
		model->AddManagedModelGeometryComponent(another_onc, Internal_CreateManagedAttributes(layer_index2, L"test_spiral_Bspline"));
		PrintPosAndTan(*another_onc, "after_smoothing_pos&tan");
		PrintCurvature(*another_onc, "after_smoothing_curvature");
	}

	void Smoothing3DControlPolygon(ON_NurbsCurve &onc, ON_3dPoint Pa, ON_3dPoint Pb, ON_3dVector Ta, ON_3dVector Tb)
	{
		int n = onc.CVCount();
		int m = onc.CVCount() - 1;
		if (m <= 3)
		{
			return;
		}

		// ComputeAngle();
		std::vector<double> Angle;
		ON_3dPoint p1, p2, p3;
		Angle.push_back(0);
		for (int i = 0; i < n - 2; i++)
		{
			onc.GetCV(i, p1);
			onc.GetCV(i + 1, p2);
			ON_3dVector v1 = p2 - p1;
			onc.GetCV(i + 2, p1);
			ON_3dVector v2 = p1 - p2;
			double theta = ON_3dVector::Angle(v1, v2);
			Angle.push_back(theta);
		}
		Angle.push_back(0);

		int max_count = 10000;
		for (int i = 0; i < n; i++)
		{
			if (Angle[i] > PI / 2 || Angle[i] < -PI / 2)
			{
				max_count = 2;
			}
		}
		int s_count = 1;
		double thetaDD = 10.0;
		double lengthavg = 0.0;
		onc.GetCV(0, p1);
		onc.GetCV(m, p2);
		double lengthbound = 2 * (p1 - p2).Length();

		// Compute Normal Vector
		std::vector<ON_3dVector> normal;
		normal.push_back(ON_3dVector(0, 0, 0));
		for (int i = 0; i < n - 2; i++)
		{
			onc.GetCV(i, p1);
			onc.GetCV(i + 1, p2);
			onc.GetCV(i + 2, p3);
			ON_3dVector vv = ON_3dVector::CrossProduct(p2 - p1, p3 - p2);
			// vv.Unitize();
			normal.push_back(vv);
		}
		normal.push_back(ON_3dVector(0, 0, 0));

		while (s_count < max_count && thetaDD > 1e-6 && lengthavg < lengthbound)
		{
			// Compute New Vertex 2--n-2
			vector<double> new_angle = Angle;
			vector<ON_3dVector> new_normal = normal;

			for (int i = 2; i < m - 1; i++)
			{
				new_angle[i] = (new_angle[i - 1] + new_angle[i] + new_angle[i + 1]) / 3;
				onc.GetCV(i - 1, p1);
				onc.GetCV(i + 1, p2);
				new_normal[i] = (new_normal[i - 1] + new_normal[i] + new_normal[i + 1]) / 3;
				ON_3dVector beta = ON_3dVector::CrossProduct(p2 - p1, new_normal[i]);
				beta.Unitize();
				ON_3dPoint re = (p1 + p2) / 2 + beta * tan(new_angle[i] / 2);
				onc.SetCV(i, re);
			}
		}
	}
}