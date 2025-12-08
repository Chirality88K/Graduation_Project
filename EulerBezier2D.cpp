#include "EulerBezier2D.h"
#include "write3dm.h"
#include "ChiralityMathTools.h"
#include <assert.h>

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

	bool EulerBezierWeakCheck(const ON_BezierCurve* OBC)
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
		double s2 = sqrt(sum / m) / avg;
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
		return true;
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
		//lengthbound = (std::max)(lengthavg, 10.0);

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

	void EulerBezierWeakInterpolation(ON_BezierCurve* OBC, int max_vtx_num)
	{
		while (OBC->CVCount() < max_vtx_num && !EulerBezierWeakCheck(OBC))
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

	ON_NurbsCurve SmoothingCornerWithSlope(ON_3dPoint start, ON_3dPoint corner, ON_3dPoint end, double alpha)
	{
		ON_3dVector v0 = corner - start;
		ON_3dVector v1 = end - corner;
		double product = ON_3dVector::DotProduct(v0, v1) / v1.Length() / v0.Length();
		double PHI = acos(product);
		if (v0.x * v1.y - v0.y * v1.x < 0)
		{
			PHI = -PHI;
		}
		ON_BezierCurve part1;
		ON_BezierCurve part2;
		int n = 4;
		ON_3dVector Ts = v0;
		Ts.Unitize();
		while (n < 20)
		{
			double sum = 0.0;
			for (int k = 1; k <= n - 2; ++k)
			{
				sum += pow(k, -1 / alpha);
			}
			double deltatheta_alpha = PHI / sum / 2;
			ON_3dVector temp = Ts;
			ON_3dVector D = Ts;
			for (int i = 0; i < n - 1; ++i)
			{
				if (i < 2)
				{
					temp.Rotate(deltatheta_alpha * i, ON_3dVector(0, 0, 1));
				}
				else
				{
					temp.Rotate(deltatheta_alpha * pow(i, -1 / alpha), ON_3dVector(0, 0, 1));
				}
				D += temp;
			}
			double pro = ON_3dVector::DotProduct(D, Ts) / D.Length();
			pro = (std::min)(1.0, pro);
			pro = (std::max)(-1.0, pro);
			double beta = acos(pro);
			if (Ts.x * D.y - Ts.y * D.x < 0) {
				beta = -beta;
			}
			double length = cos(PHI / 2) / cos(PHI / 2 - beta) * v0.Length() / D.Length();
			temp = Ts * length;
			double theta_n_2 = pow(n - 3, -1 / alpha) * deltatheta_alpha;
			double theta_n_1 = pow(n - 2, -1 / alpha) * deltatheta_alpha;
			double s1 = -(n + 1) * sin(theta_n_1) - (n - 2) * sin(theta_n_1 + theta_n_2) + 3 * (n - 1) * cos(theta_n_1) * sin(theta_n_1);
			if (s1 * deltatheta_alpha > 0 || n == 19)
			{
				part1.Create(3, false, n + 1);
				part1.SetCV(0, start);
				ON_3dPoint p = start;
				for (int i = 0; i < n; ++i)
				{
					p += temp;
					part1.SetCV(i + 1, p);
					if (i < 2)
					{
						temp.Rotate(deltatheta_alpha * i, ON_3dVector(0, 0, 1));
					}
					else
					{
						temp.Rotate(deltatheta_alpha * pow(i, -1 / alpha), ON_3dVector(0, 0, 1));
					}
				}
				break;
			}
			else
			{
				++n;
			}
		}
		v0.Unitize();
		v1.Unitize();
		GenerateSymmetry(&part2, &part1, corner, v0 - v1);
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
			ON_3dPoint Start = (i == 0) ? P[N - 1] : P[i - 1];
			ON_3dPoint End = (i == N - 1) ? P[0] : P[i + 1];
			ON_3dPoint Corner = P[i];
			ON_3dVector v0 = Corner - Start;
			ON_3dVector v1 = End - Corner;

			ON_NurbsCurve onc = GenerateSmoothingCurve(Start * (1.0 / 3.0) + Corner * (2.0 / 3.0), Corner, End * (1.0 / 3.0) + Corner * (2.0 / 3.0));
			
			ChiralityAddNurbsCurve(model, onc, L"EulerBezier", layer_index);
			ChiralityDebugforR(onc, "EulerBezier Debug for R" + std::to_string(i));
			Q[i * 2] = Start * (1.0 / 3.0) + Corner * (2.0 / 3.0);
			Q[i * 2 + 1] = End * (1.0 / 3.0) + Corner * (2.0 / 3.0);
		}
		
		const int polygon_layer_index = model->AddLayer(L"EulerBezier2dControlPoints", ON_Color::Black);
		ON_PolylineCurve *opc = new ON_PolylineCurve(ON_Polyline(Parray));
		ON_3dmObjectAttributes *attributes = new ON_3dmObjectAttributes();
		attributes->m_layer_index = polygon_layer_index;
		attributes->m_name = L"EulerBezier2dControlPoints";
		model->AddManagedModelGeometryComponent(opc, attributes);
		
		const int lines_layer_index = model->AddLayer(L"MidLines", ON_Color::SaturatedMagenta);
		for (int i = 0; i < N; ++i)
		{
			ON_3dmObjectAttributes* attributes_lines = new ON_3dmObjectAttributes();
			attributes_lines->m_layer_index = lines_layer_index;
			attributes_lines->m_name = (L"Lines" + std::to_wstring(i)).c_str();
			ON_LineCurve* olc = new ON_LineCurve(Q[(2 * i + 1) % (2 * N)], Q[(2 * i + 2) % (2 * N)]);
			model->AddManagedModelGeometryComponent(olc, attributes_lines);
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
			ON_NurbsCurve onc = GenerateSmoothingCurve(Start * 0.5 + Corner * 0.5, Corner, End * 0.5 + Corner * 0.5);
			ChiralityAddNurbsCurve(model, onc, L"curve" + std::to_wstring(i + 1), curves_layer_index);
			ChiralityDebugforR(onc, "Bezier 2d debug for R " + std::to_string(i));
		}
	}

	double Compute_delta_theta(ON_2dPoint ps, ON_2dPoint pe, ON_2dVector vs, ON_2dVector ve, int cv_cnt)
	{
		const double alpha = ChiralityMath::ComputeSignedAngle(pe - ps, vs);
		const double beta = ChiralityMath::ComputeSignedAngle(pe - ps, ve);
		const int n = cv_cnt - 1;
		auto F = [alpha, beta, n](double delta_theta)->double {
			double sum = 0.0;
			for (int k = 1; k <= n; ++k)
			{
				sum += sin(alpha * double(n - k) / double(n - 1) + beta * double(k - 1) / double(n - 1) - 0.5 * (k - 1) * (n - k) * delta_theta);
			}
			return sum;
		};

		auto DF = [alpha, beta, n](double delta_theta)->double {
			double sum = 0.0;
			for (int k = 2; k <= n - 1; ++k)
			{
				sum += cos(alpha * double(n - k) / double(n - 1) + beta * double(k - 1) / double(n - 1) - 0.5 * (k - 1) * (n - k) * delta_theta)
					* (-0.5 * (k - 1) * (n - k));
			}
			return sum;
		};
		//double x0 = ChiralityMath::Newton(F, DF, 0.0);
		//double x1 = ChiralityMath::Newton(F, DF, 0.001);
		//double x2 = ChiralityMath::Newton(F, DF, -0.001);
		//double x = abs(x0) < abs(x1) ? x0 : x1;
		//x = abs(x) < abs(x2) ? x : x2;

		const double F_0 = F(0.0);
		double s_theta = -5.0 / cv_cnt;
		double e_theta = 5.0 / cv_cnt;
		constexpr int N = 100;
		double t_e = 10;

		int k = (n + 1) / 2;
		double alpha_gap = alpha > 0 ? PI - alpha : PI + alpha;
		double beta_gap = beta > 0 ? PI - beta : PI + beta;
		double t_max = (std::max)(alpha_gap, beta_gap) * 2 / (k - 1) / (n - k);

		if (alpha * beta > 0)
		{
			double step = alpha > 0 ? t_max / N : -t_max / N;
			double t = step;
			while (F(t) * F_0 > 0)
			{
				t += step;
			}
			t_e = t;
		}
		else
		{
			double step = t_max / N;
			double t = step;
			while (F(t) * F_0 > 0 && F(-t) * F_0 > 0)
			{
				t += step;
			}
			t_e = F(t) * F_0 > 0 ? -t : t;
		}
		assert(t_e < 10);
		double t_s = 0.0;
		if (t_s > t_e)
		{
			std::swap(t_s, t_e);
		}
		double x = ChiralityMath::Bisection(F, t_s, t_e);
		return x;
	}

	ON_BezierCurve ComputeEulerBezier2D_Directly(ON_2dPoint ps, ON_2dPoint pe, ON_2dVector vs, ON_2dVector ve, int cv_cnt, double& error)
	{
		bool pass = false;
		const int max_cv_cnt = 50;
		cv_cnt = (std::max)(10, cv_cnt);
		cv_cnt = (std::min)(max_cv_cnt - 1, cv_cnt);
		std::vector<double> angles;
		angles.reserve(max_cv_cnt);
		double sum = 0.0;
		double alpha = ChiralityMath::ComputeSignedAngle(pe - ps, vs);
		double beta = ChiralityMath::ComputeSignedAngle(pe - ps, ve);
		while (!pass && cv_cnt < max_cv_cnt)
		{
			angles.clear();
			int n = cv_cnt - 1;
			double delta_theta = Compute_delta_theta(ps, pe, vs, ve, cv_cnt);
			pass = true;
			for (int k = 1; k <= n; ++k)
			{
				angles.push_back(alpha * double(n - k) / double(n - 1) + beta * double(k - 1) / double(n - 1)
					- 0.5 * (k - 1) * (n - k) * delta_theta);
				if (abs(angles.back()) >= PI)
				{
					pass = false;
				}
				if (k > 1 && abs(angles[k - 1] - angles[k - 2]) >= PI / 2)
				{
					pass = false;
				}
			}
			if (pass)
			{
				double a1 = angles[1] - angles[0];
				double a2 = angles[2] - angles[1];
				double an_1 = angles[n - 1] - angles[n - 2];
				double an_2 = angles[n - 2] - angles[n - 3];
				double Deltatheta = (an_2 - a1);
				double s0 = (n + 1) * sin(a1) + (n - 2) * sin(a1 + a2) - 3 * (n - 1) * sin(a1) * cos(a1);
				double s1 = -(n + 1) * sin(an_1) - (n - 2) * sin(an_2 + an_1) + 3 * (n - 1) * sin(an_1) * cos(an_1);
				pass = Deltatheta * s0 >= 0 && s0 * s1 >= 0;
			}
			++cv_cnt;
		}
		cv_cnt = angles.size() + 1;
		//assert(sum > 1e-6);
		sum = 0.0;
		for (double t : angles)
		{
			sum += cos(t);
		}
		double L = ps.DistanceTo(pe) / sum;
		if (L < 0)
		{
			L = -L;
		}
		ON_BezierCurve obc(2, false, cv_cnt);
		vs.Unitize(); ve.Unitize();
		obc.SetCV(0, ON_3dPoint(ps.x, ps.y, 0.0));
		ON_3dPoint P = ps;
		for (int i = 1; i < cv_cnt; ++i)
		{
			ON_3dVector v = pe - ps;
			v.Unitize();
			v.Rotate(angles[i - 1], ON_3dVector::ZAxis);
			P = P + L * v;
			obc.SetCV(i, P);
		}
		error = obc.PointAt(0).DistanceTo(ps) + obc.PointAt(1).DistanceTo(pe) +
			(obc.TangentAt(0) - vs).Length() + (obc.TangentAt(1) - ve).Length();
		if (error > 1e-6)
		{
			CHIRALITY_WARN(std::string("Euler_Bezier!!"));
			std::cout << "Direct Error: " << error << "\n";
		}
		return obc;
	}

	double Compute_L_for_fixed_cv_cnt(ON_2dPoint ps, ON_2dPoint pe, ON_2dVector vs, ON_2dVector ve, int cv_cnt, double* angs)
	{
		assert(cv_cnt >= 10);
		double delta_theta = Compute_delta_theta(ps, pe, vs, ve, cv_cnt);
		double alpha = ChiralityMath::ComputeSignedAngle(pe - ps, vs);
		double beta = ChiralityMath::ComputeSignedAngle(pe - ps, ve);
		std::vector<double> angles;
		int n = cv_cnt - 1;
		angles.reserve(n);
		for (int k = 1; k <= n; ++k)
		{
			angles.push_back(alpha * double(n - k) / double(n - 1) + beta * double(k - 1) / double(n - 1)
				- 0.5 * (k - 1) * (n - k) * delta_theta);
		}
		if (angs != nullptr)
		{
			for (int i = 0; i < angles.size(); ++i)
			{
				*(angs + i) = angles[i];
			}
		}
		double sum = 0.0;
		for (double t : angles)
		{
			sum += cos(t);
		}
		double L = ps.DistanceTo(pe) / sum;
		if (L < 0)
		{
			L = -L;
		}
		return L;
	}
}