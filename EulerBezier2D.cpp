#include "EulerBezier2D.h"
#include "write3dm.h"

namespace EulerBezier2D
{
	std::vector<double> ComputeLength(const ON_BezierCurve *OBC)
	{
		std::vector<double> l;
		ON_3dPoint p0, p1;
		for (int i = 1; i < OBC->CVCount(); ++i)
		{
			OBC->GetCV(i - 1, p0);
			OBC->GetCV(i, p1);
			l.push_back(p0.DistanceTo(p1));
		}
		return l;
	}

	std::vector<double> ComputeAngle(const ON_BezierCurve *OBC)
	{
		std::vector<double> a;
		ON_3dPoint p0, p1, p2;
		ON_3dVector v0, v1;
		a.push_back(0.0);
		double pro;
		double angle;
		for (int i = 1; i < OBC->CVCount() - 1; ++i)
		{
			OBC->GetCV(i - 1, p0);
			OBC->GetCV(i, p1);
			OBC->GetCV(i + 1, p2);
			v0 = p1 - p0;
			v1 = p2 - p1;
			pro = ON_3dVector::DotProduct(v0, v1) / v0.Length() / v1.Length();
			pro = (std::min)(1.0, pro);
			pro = (std::max)(-1.0, pro);
			angle = acos(pro);
			if (v0.x * v1.y - v0.y * v1.x < 0)
			{
				angle = -angle;
			}
			a.push_back(angle);
		}
		a.push_back(0.0);
		return a;
	}

	bool EulerBezierSpiralCheck(const ON_BezierCurve *OBC)
	{
		int n = OBC->CVCount();
		int m = n - 1;
		if (n < 2)
		{
			return false;
		}
		std::vector<double> Length = ComputeLength(OBC);
		double sum = 0;
		for (int i = 0; i < m; i++)
		{
			sum = sum + Length[i];
		}
		double avg = sum / m;
		sum = 0;
		for (int i = 0; i < m; i++)
		{
			sum = sum + (Length[i] - avg) * (Length[i] - avg);
		}
		double s2 = sum / m;
		if (s2 > 0.1)
		{
			return false;
		}
		std::vector<double> Angle = ComputeAngle(OBC);
		for (int i = 0; i < n; i++)
		{
			if (Angle[i] > PI / 2 || Angle[i] < -PI / 2)
			{
				return false;
			}
		}
		if (m > 3)
		{
			for (int i = 2; i < m - 1; i++)
			{
				if (abs(Angle[i] * 2 - Angle[i - 1] - Angle[i + 1]) > 1e-6)
				{
					return false;
				}
			}
		}
		double Deltatheta = (Angle[n - 2] - Angle[1]) / (n - 3);
		double s0 = (m + 1) * sin(Angle[1]) + (m - 2) * sin(Angle[1] + Angle[2]) - 3 * (m - 1) * sin(Angle[1]) * cos(Angle[1]);
		double s1 = -(m + 1) * sin(Angle[n - 2]) - (m - 2) * sin(Angle[n - 3] + Angle[n - 2]) + 3 * (m - 1) * sin(Angle[n - 2]) * cos(Angle[n - 2]);
		if (Deltatheta * s0 >= 0 && s0 * s1 >= 0)
		{
			return true;
		}
		else
		{
			return false;
		}
	}

	void SmoothingBezierControlPolygon(ON_BezierCurve *OBC)
	{
		int m = OBC->CVCount() - 1;
		if (m <= 3)
		{
			return;
		}
		std::vector<double> Angle = ComputeAngle(OBC);
		int max_count = 10000;
		for (int i = 0; i < m + 1; i++)
		{
			if (Angle[i] > PI / 2 || Angle[i] < -PI / 2)
			{
				max_count = 2;
			}
		}
		int s_count = 1;
		double thetaDD = 10.0;
		double lengthavg = 0.0;
		ON_3dPoint p0, pm;
		OBC->GetCV(0, p0);
		OBC->GetCV(m, pm);
		double lengthbound = 2 * (p0 - pm).Length();
		lengthbound = (std::max)(lengthavg, 10.0);

		while (s_count < max_count && thetaDD > 1e-6 && lengthavg < lengthbound)
		{
			std::vector<double> new_angle = Angle;

			for (int i = 2; i < m - 1; i++)
			{
				new_angle[i] = (new_angle[i - 1] + new_angle[i] + new_angle[i + 1]) / 3;
				OBC->GetCV(i + 1, pm);
				OBC->GetCV(i - 1, p0);
				ON_3dPoint re = (pm - p0) / 2;
				re.Set(-re.y, re.x, 0);
				re = (p0 + pm) / 2 - re * tan(new_angle[i] / 2);
				OBC->SetCV(i, re);
			}
			std::vector<double> Length = ComputeLength(OBC);

			double sum = 0;
			for (int i = 0; i < m; i++)
			{
				sum = sum + Length[i];
			}
			lengthavg = sum / m;
			OBC->GetCV(0, p0);
			OBC->GetCV(1, pm);
			ON_3dPoint p = p0 + (pm - p0) / Length[0] * lengthavg;
			OBC->SetCV(1, p);
			OBC->GetCV(m, pm);
			OBC->GetCV(m - 1, p0);
			p = pm - (pm - p0) / Length[m - 1] * lengthavg;
			OBC->SetCV(m - 1, p);

			Angle = ComputeAngle(OBC);

			double maxtheta;
			thetaDD = 0;
			if (m > 3)
			{
				for (int i = 2; i < m - 1; i++)
				{
					maxtheta = abs(Angle[i] * 2 - Angle[i - 1] - Angle[i + 1]);
					if (maxtheta > thetaDD)
					{
						thetaDD = maxtheta;
					}
				}
			}
			s_count += 1;
		}
	}

	void Elevate(ON_BezierCurve *OBC)
	{
		int n = OBC->CVCount();
		if (n < 2)
		{
			return;
		}
		std::vector<ON_3dPoint> parr;
		ON_3dPoint p, q;
		OBC->GetCV(0, p);
		parr.push_back(p);
		for (int i = 1; i < n; i++)
		{
			OBC->GetCV(i - 1, p);
			OBC->GetCV(i, q);
			parr.push_back(p * (i * 1.0 / n) + q * (1 - i * 1.0 / n));
		}
		OBC->GetCV(n - 1, q);
		parr.push_back(q);
		OBC->Create(2, false, n + 1);
		for (int i = 0; i < n + 1; i++)
		{
			OBC->SetCV(i, parr[i]);
		}
	}

	void EulerBezierSpiralInterpolation(ON_BezierCurve *OBC, int max_vtx_num)
	{
		while (OBC->CVCount() < max_vtx_num && !EulerBezierSpiralCheck(OBC))
		{
			Elevate(OBC);
			SmoothingBezierControlPolygon(OBC);
		}
	}

	void SmoothingCorner(ON_BezierCurve *OBC, ON_3dPoint Ps, ON_3dPoint O, double alpha)
	{
		int n = 4;
		ON_3dVector Ts = (O - Ps);
		Ts.Unitize();

		while (!EulerBezierSpiralCheck(OBC) && n < 20)
		{
			double deltatheta = alpha / (n - 2) / (n - 1);
			ON_3dVector temp = Ts;
			ON_3dVector D = Ts;
			for (int i = 0; i < n - 1; ++i)
			{
				temp.Rotate(deltatheta * i, ON_3dVector(0, 0, 1));
				D += temp;
			}
			double pro = ON_3dVector::DotProduct(D, Ts) / D.Length();
			pro = (std::min)(1.0, pro);
			pro = (std::max)(-1.0, pro);
			double beta = acos(pro);
			if (Ts.x * D.y - Ts.y * D.x < 0) {
				beta = -beta;
			}
			double length = cos(alpha / 2) / cos(alpha / 2 - beta) * (O - Ps).Length() / D.Length();
			temp = Ts * length;
			OBC->Create(2, false, n + 1);
			OBC->SetCV(0, Ps);
			ON_3dPoint p = Ps;
			for (int i = 0; i < n; ++i)
			{
				p += temp;
				OBC->SetCV(i + 1, p);
				temp.Rotate(deltatheta * i, ON_3dVector(0, 0, 1));
			}
			++n;
		}
	}

	ON_NurbsCurve GenerateSmoothingCurve(ON_3dPoint start, ON_3dPoint corner, ON_3dPoint end)
	{
		ON_3dVector v0 = corner - start;
		ON_3dVector v1 = end - corner;
		double product = ON_3dVector::DotProduct(v0, v1) / v1.Length() / v0.Length();
		double alpha = acos(product);
		if (v0.x * v1.y - v0.y * v1.x < 0)
		{
			alpha = -alpha;
		}
		ON_BezierCurve part1;
		ON_BezierCurve part2;
		SmoothingCorner(&part1, start, corner, alpha);
		v0.Unitize();
		v1.Unitize();
		GenerateSymmetry(&part2, &part1, corner, v1 - v0);
		ON_NurbsCurve onc1;
		part1.GetNurbForm(onc1);
		ON_NurbsCurve onc2;
		part2.GetNurbForm(onc2);
		onc2.Reverse();
		onc1.Append(onc2);
		onc1.SetDomain(0, 1);
		return onc1;
	}

	void GenerateSymmetry(ON_BezierCurve *result, const ON_BezierCurve *OBC, ON_3dPoint O, ON_3dVector v)
	{
		v.Unitize();
		result->Create(2, false, OBC->Order());
		ON_3dPoint p;
		for (int i = 0; i < OBC->CVCount(); ++i)
		{
			OBC->GetCV(i, p);
			result->SetCV(i, O + 2 * ON_3dVector::DotProduct(p - O, v) * v - (p - O));
		}
	}

	void EulerBezier2dTest(ONX_Model *model)
	{
		const int N = 6;
		ON_3dPoint P[N];
		ON_3dPointArray Parray;
		ON_3dPoint Q[2 * N];
		const int layer_index = model->AddLayer(L"EulerBezier2dTest", ON_Color::SaturatedMagenta);
		for (int i = 0; i < N; ++i)
		{
			P[i] = ON_3dPoint(cos(PI * 2 * i / N), sin(PI * 2 * i / N), 0) * 10;
			Parray.Append(P[i]);
		}
		Parray.Append(P[0]);
		for (int i = 0; i < N; ++i)
		{
			ON_BezierCurve *obc0 = new ON_BezierCurve();
			ON_BezierCurve *obc1 = new ON_BezierCurve();
			ON_3dPoint Start = (i == 0) ? P[N - 1] : P[i - 1];
			ON_3dPoint End = (i == N - 1) ? P[0] : P[i + 1];
			ON_3dPoint Corner = P[i];
			ON_3dVector v0 = Corner - Start;
			ON_3dVector v1 = End - Corner;
			double product = ON_3dVector::DotProduct(v0, v1) / v1.Length() / v0.Length();
			double alpha = acos(product);
			if (v0.x * v1.y - v0.y * v1.x < 0)
			{
				alpha = -alpha;
			}
			SmoothingCorner(obc0, Start * (1.0 / 3.0) + Corner * (2.0 / 3.0), Corner, alpha);
			v0.Unitize();
			v1.Unitize();
			GenerateSymmetry(obc1, obc0, Corner, v1 - v0);

			ON_NurbsCurve **onc0 = new ON_NurbsCurve *();
			(*onc0) = new ON_NurbsCurve();
			obc0->GetNurbForm(**onc0);
			ON_NurbsCurve **onc1 = new ON_NurbsCurve *();
			(*onc1) = new ON_NurbsCurve();
			obc1->GetNurbForm(**onc1);

			ON_3dmObjectAttributes *attributes0 = new ON_3dmObjectAttributes();
			attributes0->m_layer_index = layer_index;
			attributes0->m_name = (L"EulerBezier_No." + std::to_wstring(i + 1)).c_str();
			model->AddManagedModelGeometryComponent(*onc0, attributes0);
			ON_3dmObjectAttributes *attributes1 = new ON_3dmObjectAttributes();
			attributes1->m_layer_index = layer_index;
			attributes1->m_name = (L"EulerBezier1_No." + std::to_wstring(i + 1)).c_str();
			model->AddManagedModelGeometryComponent(*onc1, attributes1);
			delete onc0;
			delete onc1;

			Q[i * 2] = Start * (1.0 / 3.0) + Corner * (2.0 / 3.0);
			Q[i * 2 + 1] = End * (1.0 / 3.0) + Corner * (2.0 / 3.0);
		}
		const int polygon_layer_index = model->AddLayer(L"EulerBezier2dControlPoints", ON_Color::Black);
		ON_PolylineCurve *opc = new ON_PolylineCurve(ON_Polyline(Parray));
		ON_3dmObjectAttributes *attributes = new ON_3dmObjectAttributes();
		attributes->m_layer_index = polygon_layer_index;
		attributes->m_name = L"EulerBezier2dControlPoints";
		model->AddManagedModelGeometryComponent(opc, attributes);

		const int lines_layer_index = model->AddLayer(L"Lines", ON_Color::SaturatedMagenta);
		ON_3dmObjectAttributes *attributes_lines = new ON_3dmObjectAttributes();
		attributes_lines->m_layer_index = lines_layer_index;
		attributes_lines->m_name = L"Lines";
		for (int i = 0; i < N; ++i)
		{
			model->AddManagedModelGeometryComponent(new ON_LineCurve(
														Q[(2 * i + 1) % (2 * N)], Q[(2 * i + 2) % (2 * N)]),
													attributes_lines);
		}
	}

	void YangMethodtest(ONX_Model *model)
	{
		ON_BezierCurve *obc = new ON_BezierCurve(2, false, 4);
		obc->SetCV(0, ON_3dPoint(0, 0, 0));
		obc->SetCV(1, ON_3dPoint(0, -5, 0));
		obc->SetCV(2, ON_3dPoint(8, -5, 0));
		obc->SetCV(3, ON_3dPoint(10, 6, 0));
		EulerBezierSpiralInterpolation(obc, 50);

		ON_NurbsCurve **onc0 = new ON_NurbsCurve *();
		(*onc0) = new ON_NurbsCurve();
		obc->GetNurbForm(**onc0);
		const int layer_index = model->AddLayer(L"Yang", ON_Color::SaturatedMagenta);
		ON_3dmObjectAttributes *attributes = new ON_3dmObjectAttributes();
		attributes->m_layer_index = layer_index;
		attributes->m_name = L"yangtest0";

		model->AddManagedModelGeometryComponent(*onc0, attributes);

		delete onc0;
	}

	void Pentagram(ONX_Model *model)
	{
		// Compute 10 points
		double Acute_vertices_distance_to_origin = 10.0;
		double Blunt_vertices_distance_to_origin = Acute_vertices_distance_to_origin * sin(PI / 10.0) / sin(3 * PI / 10.0);
		double Acute_vertices_polar_angle[5] = {PI / 10.0, PI / 2.0, 162 * PI / 180.0, 234 * PI / 180.0, 306 * PI / 180.0};
		double Blunt_vertices_polar_angle[5] = {54 * PI / 180.0, 126 * PI / 180.0, 198 * PI / 180.0, 270 * PI / 180.0, 342 * PI / 180.0};
		ON_3dPoint Acute_vertices[5];
		ON_3dPoint Blunt_vertices[5];
		for (int i = 0; i < 5; ++i)
		{
			Acute_vertices[i] = ON_3dPoint(cos(Acute_vertices_polar_angle[i]), sin(Acute_vertices_polar_angle[i]), 0) * Acute_vertices_distance_to_origin;
		}
		for (int i = 0; i < 5; ++i)
		{
			Blunt_vertices[i] = ON_3dPoint(cos(Blunt_vertices_polar_angle[i]), sin(Blunt_vertices_polar_angle[i]), 0) * Blunt_vertices_distance_to_origin;
		}
		// connect these 10 points with lines
		ON_3dPointArray Parray;
		for (int i = 0; i < 5; ++i)
		{
			Parray.Append(Acute_vertices[i]);
			Parray.Append(Blunt_vertices[i]);
		}
		Parray.Append(Acute_vertices[0]);
		// Add these lines to model
		const int polygon_layer_index = model->AddLayer(L"Frame", ON_Color::Black);
		ON_PolylineCurve *opc = new ON_PolylineCurve(ON_Polyline(Parray));
		ON_3dmObjectAttributes *attributes = new ON_3dmObjectAttributes();
		attributes->m_layer_index = polygon_layer_index;
		attributes->m_name = L"Pentagram_frame";
		model->AddManagedModelGeometryComponent(opc, attributes);
		// Compute Smoothing corner curves
		const int curves_layer_index = model->AddLayer(L"Smoothing curves", ON_Color::SaturatedMagenta);
		for (int i = 0; i < 10; ++i)
		{
			ON_3dPoint Start = (i == 0) ? Parray[9] : Parray[i - 1];
			ON_3dPoint End = (i == 9) ? Parray[0] : Parray[i + 1];
			ON_3dPoint Corner = Parray[i];
			ChiralityAddNurbsCurve(model, GenerateSmoothingCurve(Start * 0.5 + Corner * 0.5, Corner, End * 0.5 + Corner * 0.5), L"curve" + std::to_wstring(i + 1), curves_layer_index);
		}
	}
}