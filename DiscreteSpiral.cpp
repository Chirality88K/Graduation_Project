#include "DiscreteSpiral.h"
#include <assert.h>
#include "thirdparty/eigen/Eigen/Dense"

void DiscreteSpiral::SimpleTest(ONX_Model* model, const FrenetFrame& f1, const FrenetFrame& f2, const std::string& name)
{
	ON_NurbsCurve onc = Interpolate(f1, f2);
	const int layer_index = model->AddLayer(StringToWString(name).c_str(), ON_Color::SaturatedGold);
	ChiralityAddNurbsCurve(model, onc, StringToWString(name).c_str(), layer_index);
	ChiralityDebugforR(onc, name);
}

ON_NurbsCurve DiscreteSpiral::Interpolate(const FrenetFrame& f1, const FrenetFrame& f2)
{
	double Length = f1.GetPos().DistanceTo(f2.GetPos());
	std::vector<ON_3dPoint> vp = { f1.GetPos(), f1.GetPos() + f1.GetAlpha() * Length / 4.0,
								  f2.GetPos() - f2.GetAlpha() * Length / 4.0, f2.GetPos() };
	DiscreteSpiral ds(vp, f1, f2);
	for (int i = 0; i < 10; ++i)
	{
		ds.Elevate();
		int max_iter = 200;
		while (!ds.EndIteration()&&max_iter>0)
		{
			ds.SmoothFrenetFrameWithAngle();
			max_iter--;
		}
		if (max_iter <= 0)
		{
			CHIRALITY_WARN(std::string("��������̫���ˣ���û�����ã���"));
		}
	}
	return ds.GetBezier();
}

void DiscreteSpiral::TestExit(ONX_Model* model)
{
	for (int times = 0; times < 10; times++)
	{
		int N = 6;
		std::vector <ON_3dPoint> vp;
		vp.push_back(ON_3dPoint(0, 0, 0));
		vp.push_back(ON_3dPoint(2.5, 0, 0));
		for (int i = 0; i < N; ++i)
		{
			ON_3dPoint p;
			while (true)
			{
				p = ChiralityMath::GetRandomPoint(vp[i + 1], 2.0, 3.0);
				double product = ON_3dVector::DotProduct((p - vp[i + 1]), vp[i + 1] - vp[i]);
				double cos_theta = product / (p - vp[i + 1]).Length() / (vp[i + 1] - vp[i]).Length();
				if (cos_theta < 0.99 && cos_theta>0.1)
				{
					break;
				}
			}
			vp.push_back(p);
		}
		FrenetFrame f1(vp[0], vp[1] - vp[0], vp[2] - vp[1]);
		FrenetFrame f2(vp[N + 1], vp[N + 1] - vp[N], vp[N] - vp[N - 1]);
		DiscreteSpiral ds(vp, f1, f2);
		for (int i = 0; i < 10; ++i)
		{
			ds.Elevate();
			int max_iter = 200;
			while (!ds.EndIteration() && max_iter > 0)
			{
				ds.SmoothFrenetFrameWithCur();
				//ds.SmoothFrenetFrameWithAngle();
				max_iter--;
			}
			if (max_iter <= 0)
			{
				CHIRALITY_WARN(std::string("��������̫���ˣ���û�����ã���"));
			}
		}
		ON_NurbsCurve onc = ds.GetBezier();
		const int layer_index = model->AddLayer((std::wstring(L"Random_Test") + std::to_wstring(times)).c_str(), ON_Color::SaturatedGold);
		ChiralityAddNurbsCurve(model, onc, std::wstring(L"Random_Test")+std::to_wstring(times), layer_index);
		ChiralityAddLines(model, vp, std::wstring(L"Random_Test_Polygon") + std::to_wstring(times), layer_index);
		ChiralityDebugforR(onc, std::string("Random_Test") + std::to_string(times));
	}
}

void DiscreteSpiral::DiscreteSpiralTest(ONX_Model *model)
{
	double a = 10.0;
	double b = -0.1;
	double alpha = PI / 3;
	auto conic_spiral = [a, b, alpha](double theta) -> FrenetFrame
	{
		ON_3dPoint p = ON_3dPoint(sin(alpha) * cos(theta), sin(alpha) * sin(theta), cos(alpha)) * a * exp(b * theta);
		ON_3dVector der = ON_3dVector(b * sin(alpha) * cos(theta) - sin(alpha) * sin(theta), b * sin(alpha) * sin(theta) + sin(alpha) * cos(theta), b * cos(alpha));
		der.Unitize();
		ON_3dVector derder = (a * b * exp(b * theta)) * ON_3dVector(b * sin(alpha) * cos(theta) - sin(alpha) * sin(theta), b * sin(alpha) * sin(theta) + sin(alpha) * cos(theta), b * cos(alpha)) + (a * exp(b * theta)) * ON_3dVector(-b * sin(alpha) * sin(theta) - sin(alpha) * cos(theta),
																																																									   b * sin(alpha) * cos(theta) - sin(alpha) * sin(theta), 0);
		ON_3dVector N = ON_3dVector::CrossProduct(der, derder);
		ON_3dVector B = ON_3dVector::CrossProduct(N, der);
		B.Unitize();
		return FrenetFrame(p, der, B);
	};
	ON_Color color[7] = {ON_Color::SaturatedRed, ON_Color(255, 128, 0), ON_Color::SaturatedYellow,
						 ON_Color::SaturatedGreen, ON_Color::SaturatedCyan, ON_Color::SaturatedBlue, ON_Color(76, 0, 153)};
	std::string color_name[7] = {"Red", "Orange", "Yellow", "Green", "Cyan", "Blue", "Purple"};
	std::vector<ON_NurbsCurve> vec_bec;
	vec_bec.resize(7);
	for (int i = 0; i < 7; ++i)
	{
		FrenetFrame f1 = conic_spiral(double(i) * PI / 3);
		FrenetFrame f2 = conic_spiral(double(i + 1) * PI / 3); 
		vec_bec[i] = Interpolate(f1, f2);
		const int bezier_layer_index = model->AddLayer((L"bezier_layer_" + StringToWString(color_name[i])).c_str(), color[i]);
		ChiralityAddNurbsCurve(model, vec_bec[i], (L"conic_spiral_Bezier_" + StringToWString(color_name[i])).c_str(), bezier_layer_index);
	}
	ChiralityDebugforR(vec_bec, "conic_spiral");

	a = 10.0;
	b = 2.0;
	auto circular_helix = [a, b](double theta) -> FrenetFrame
	{
		ON_3dPoint p(a * cos(theta), a * sin(theta), b * theta);
		ON_3dVector der(-a * sin(theta), a * cos(theta), b);
		ON_3dVector derder(-a * cos(theta), -a * sin(theta), 0);
		ON_3dVector N = ON_3dVector::CrossProduct(der, derder);
		ON_3dVector B = ON_3dVector::CrossProduct(N, der);
		der.Unitize();
		B.Unitize();
		return FrenetFrame(p, der, B);
	};

	for (int i = 0; i < 7; ++i)
	{
		FrenetFrame f1 = circular_helix(double(i) * PI / 3);
		FrenetFrame f2 = circular_helix(double(i + 1) * PI / 3);
		vec_bec[i] = Interpolate(f1, f2);
		const int bezier_layer_index = model->AddLayer((L"bezier_layer_" + StringToWString(color_name[i])).c_str(), color[i]);
		ChiralityAddNurbsCurve(model, vec_bec[i], (L"circular_helix_Bezier_" + StringToWString(color_name[i])).c_str(), bezier_layer_index);
	}
	ChiralityDebugforR(vec_bec, "circular_helix");

	double R = 10.0;
	double ratio = 0.03;
	auto sphere_spiral = [R, ratio](double theta) -> FrenetFrame
	{
		double phi = ratio * theta;
		ON_3dPoint p = ON_3dPoint(cos(theta) * cos(phi), cos(theta) * sin(phi), sin(theta)) * R;
		ON_3dVector der(-sin(theta) * cos(phi) - ratio * cos(theta) * sin(phi),
						-sin(theta) * sin(phi) + ratio * cos(theta) * cos(phi),
						cos(theta));
		ON_3dVector derder(-cos(theta) * cos(phi) + 2 * ratio * sin(theta) * sin(phi) - ratio * ratio * cos(theta) * cos(phi),
						   -cos(theta) * sin(phi) - 2 * ratio * sin(theta) * cos(phi) - ratio * ratio * cos(theta) * sin(phi),
						   -sin(theta));
		ON_3dVector N = ON_3dVector::CrossProduct(der, derder);
		ON_3dVector B = ON_3dVector::CrossProduct(N, der);
		der.Unitize();
		B.Unitize();
		return FrenetFrame(p, der, B);
	};

	for (int i = 0; i < 7; ++i)
	{
		FrenetFrame f1 = sphere_spiral(double(i) * PI / 3);
		FrenetFrame f2 = sphere_spiral(double(i + 1) * PI / 3);
		vec_bec[i] = Interpolate(f1, f2);
		const int bezier_layer_index = model->AddLayer((L"bezier_layer_" + StringToWString(color_name[i])).c_str(), color[i]);
		ChiralityAddNurbsCurve(model, vec_bec[i], (L"sphere_spiral_Bezier_" + StringToWString(color_name[i])).c_str(), bezier_layer_index);
	}
	ChiralityDebugforR(vec_bec, "sphere_spiral");

	FrenetFrame f1 = FrenetFrame(ON_3dPoint(0, 0, 10), ON_3dVector(0, 2, -1), ON_3dVector(0, -1, -2));
	FrenetFrame f2 = FrenetFrame(ON_3dPoint(0, 0, 0), ON_3dVector(3, 1, -2), ON_3dVector(0, 1, 0));

	SimpleTest(model, f1, f2, "Test");
}

void DiscreteSpiral::DiscreteBodyTest(ONX_Model* model)
{
	FrenetFrame f1 = FrenetFrame(ON_3dPoint(3, 0, 10), ON_3dVector(0, 2, -1), ON_3dVector(0, -1, -2));
	FrenetFrame f2 = FrenetFrame(ON_3dPoint(3, 0, 0), ON_3dVector(3, 1, -2), ON_3dVector(0, 1, 0));
	ON_NurbsCurve onc = Interpolate(f1, f2);
	const int layer_index = model->AddLayer(L"Rotating", ON_Color::SaturatedGold);
	ON_Line axis(ON_3dPoint::Origin, ON_3dPoint(0, 0, 1));
	ON_NurbsSurface ons = ChiralityMath::GenerateRotating(onc, axis);
	ChiralityAddNurbsSurface(model, ons, L"Test_Rotating", layer_index);
}

void DiscreteSpiral::ComputeFrenetFrame()
{
	size_t num = mPoints.size();
	mFrenetFrame.clear();
	mFrenetFrame.resize(num);
	mFrenetFrame[0] = mBoundary_Frame[0];
	mFrenetFrame[num - 1] = mBoundary_Frame[1];
	for (size_t i = 1; i < num - 1; ++i)
	{
		ON_3dVector v1 = mPoints[i] - mPoints[i - 1];
		ON_3dVector v2 = mPoints[i + 1] - mPoints[i];
		ON_3dVector a = v1 / v1.Length() + v2 / v2.Length();
		a.Unitize();
		if (v1.IsParallelTo(v2))
		{
			mFrenetFrame[i].Set(mPoints[i], mFrenetFrame[i - 1].GetAlpha(), mFrenetFrame[i - 1].GetBeta());
			continue;
			//CHIRALITY_WARN(std::string("v1 and v2 may be parallel, which will lead to mistakes!!"));
		}
		ON_3dVector c = ON_3dVector::CrossProduct(v1, v2);
		ON_3dVector b = ON_3dVector::CrossProduct(c, a);
		mFrenetFrame[i] = FrenetFrame(mPoints[i], a, b);
	}
}

void DiscreteSpiral::SmoothFrenetFrameWithAngle()
{
	size_t num = mPoints.size();
	std::vector<double> angles;
	while(!ComputeAngles(angles))
	{
		Elevate();
	}
	for (size_t i = 2; i < num - 2; ++i)
	{
		ComputeFrenetFrame();
		ComputeAngles(angles);
		std::vector<double> cur =  ComputeDiscreteCurvature();
		std::vector<double> tor = ComputeDiscreteTorsion();
		double this_angle = (angles[i - 1] + angles[i] + angles[i + 1]) / 3;
		double tan_half_theta = sin(this_angle) / (1.0 + cos(this_angle));
		ON_3dVector NewBeta = (mFrenetFrame[i - 1].GetBeta() * cur[i - 1] + mFrenetFrame[i].GetBeta() * cur[i] + mFrenetFrame[i + 1].GetBeta() * cur[i + 1]) / 3;
		FrenetFrame f(ON_3dPoint::Origin, mPoints[i + 1] - mPoints[i - 1], NewBeta);
		NewBeta = f.GetBeta();
		mPoints[i] = (mPoints[i + 1] + mPoints[i - 1]) / 2 - NewBeta * (mPoints[i - 1].DistanceTo(mPoints[i + 1]) / 2 * tan_half_theta);
	}
	double L = 0;
	for (int i = 1; i < num; ++i)
	{
		L += mPoints[i - 1].DistanceTo(mPoints[i]);
	}
	L /= (num - 1);
	mPoints[1] = mPoints[0] + mFrenetFrame[0].GetAlpha() * L;
	mPoints[num - 2] = mPoints[num - 1] - mFrenetFrame[num - 1].GetAlpha() * L;
}

void DiscreteSpiral::SmoothFrenetFrameWithCur()
{
	size_t num = mPoints.size();
	for (size_t i = 2; i < num - 2; ++i)
	{
		ComputeFrenetFrame();
		 std::vector<double> cur = ComputeDiscreteCurvature();
		 std::vector<double> tor = ComputeDiscreteTorsion();
		double this_cur = (cur[i - 1] + cur[i] + cur[i + 1]) / 3;
		double this_tor = (tor[i - 1] + tor[i] + tor[i + 1]) / 3;
		cur[i] = this_cur;
		tor[i] = this_tor;
		double s = mFrenetFrame[i].GetPos().DistanceTo(mFrenetFrame[i - 1].GetPos()) * 0.5 +
			mFrenetFrame[i].GetPos().DistanceTo(mFrenetFrame[i + 1].GetPos()) * 0.5;
		double der_cur_s2 = ON_3dVector::DotProduct(mFrenetFrame[i].GetBeta(), mFrenetFrame[i - 1].GetAlpha() + mFrenetFrame[i + 1].GetAlpha());
		double der_tor_s2 = -ON_3dVector::DotProduct(mFrenetFrame[i].GetBeta(), mFrenetFrame[i - 1].GetGamma() + mFrenetFrame[i + 1].GetGamma());
		Eigen::Matrix<double, 3, 3> A;
		A(0, 0) = 2 - this_cur * this_cur * s * s;
		A(0, 1) = der_cur_s2;
		A(0, 2) = s * s * this_cur * this_tor;
		A(1, 0) = -der_cur_s2;
		A(1, 1) = 2 - s * s * (this_cur * this_cur + this_tor * this_tor);
		A(1, 2) = der_tor_s2;
		A(2, 0) = s * s * this_cur * this_tor;
		A(2, 1) = -der_tor_s2;
		A(2, 2) = 2 - this_tor * this_tor * s * s;
		Eigen::Matrix<double, 3, 3> B;
		ON_3dVector tmp = mFrenetFrame[i - 1].GetAlpha() + mFrenetFrame[i + 1].GetAlpha();
		B(0, 0) = tmp.x; B(0, 1) = tmp.y; B(0, 2) = tmp.z;
		tmp = mFrenetFrame[i - 1].GetBeta() + mFrenetFrame[i + 1].GetBeta();
		B(1, 0) = tmp.x; B(1, 1) = tmp.y; B(1, 2) = tmp.z;
		tmp = mFrenetFrame[i - 1].GetGamma() + mFrenetFrame[i + 1].GetGamma();
		B(2, 0) = tmp.x; B(2, 1) = tmp.y; B(2, 2) = tmp.z;
		Eigen::Matrix<double, 3, 3> X;
		X = A.partialPivLu().solve(B);

		ON_3dVector NewAlpha = mPoints[i + 1] - mPoints[i - 1];
		ON_3dVector NewBeta = ON_3dVector(X(1, 0), X(1, 1), X(1, 2));
		FrenetFrame f(ON_3dPoint::Origin, NewAlpha, NewBeta);
		NewBeta = f.GetBeta();

		mPoints[i] = (mPoints[i + 1] + mPoints[i - 1]) / 2 - 0.5 * s * s * this_cur * NewBeta;
	}
	double L = 0;
	for (int i = 1; i < num; ++i)
	{
		L += mPoints[i - 1].DistanceTo(mPoints[i]);
	}
	L /= (num - 1);
	mPoints[1] = mPoints[0] + mFrenetFrame[0].GetAlpha() * L;
	mPoints[num - 2] = mPoints[num - 1] - mFrenetFrame[num - 1].GetAlpha() * L;
}

bool DiscreteSpiral::ComputeAngles(std::vector<double>& angle) const
{
	int num = mPoints.size();
	angle.clear();
	angle.resize(num, 0.0);
	for (int i = 1; i < num - 1; ++i)
	{
		ON_3dVector v1 = mPoints[i] - mPoints[i - 1];
		ON_3dVector v2 = mPoints[i + 1] - mPoints[i];
		double product = ON_3dVector::DotProduct(v1, v2) / v1.Length() / v2.Length();
		if (product <= 0.1)
		{
			CHIRALITY_ERROR(std::string("curvature angle is too large!!!"));
			return false;
		}
		product = (std::min)(product, 1.0);
		angle[i] = acos(product);
	}
	return true;
}

std::vector<double> DiscreteSpiral::ComputeDiscreteCurvature() const
{
	int num = mPoints.size();
	std::vector<double> cur(num, 0.0);
	for (int i = 1; i < num - 1; ++i)
	{
		ON_3dVector Alpha_2 = mFrenetFrame[i].GetAlpha();
		ON_3dVector Alpha_1 = mFrenetFrame[i - 1].GetAlpha();
		ON_3dVector Alpha_3 = mFrenetFrame[i + 1].GetAlpha();
		cur[i] = sqrt(2 - ON_3dVector::DotProduct(Alpha_1, Alpha_2) -
			ON_3dVector::DotProduct(Alpha_3, Alpha_2))
			/ (mFrenetFrame[i].GetPos().DistanceTo(mFrenetFrame[i - 1].GetPos()) * 0.5 +
				mFrenetFrame[i].GetPos().DistanceTo(mFrenetFrame[i + 1].GetPos()) * 0.5);
	}
	return cur;
}

std::vector<double> DiscreteSpiral::ComputeDiscreteTorsion() const
{
	int num = mPoints.size();
	std::vector<double> tor(num, 0.0);
	for (int i = 1; i < num - 1; ++i)
	{
		ON_3dVector Gamma_2 = mFrenetFrame[i].GetGamma();
		ON_3dVector Gamma_1 = mFrenetFrame[i - 1].GetGamma();
		ON_3dVector Gamma_3 = mFrenetFrame[i + 1].GetGamma();
		tor[i] = sqrt(2 - ON_3dVector::DotProduct(Gamma_1, Gamma_2) -
			ON_3dVector::DotProduct(Gamma_3, Gamma_2))
			/ (mFrenetFrame[i].GetPos().DistanceTo(mFrenetFrame[i - 1].GetPos()) * 0.5 +
				mFrenetFrame[i].GetPos().DistanceTo(mFrenetFrame[i + 1].GetPos()) * 0.5);
		double test_sign = ON_3dVector::DotProduct(mFrenetFrame[i - 1].GetAlpha() + mFrenetFrame[i + 1].GetAlpha(), mFrenetFrame[i].GetGamma());
		if (abs(test_sign) < 1e-6)
		{
			tor[i] = 0.0;
		}
		else if (test_sign < -1e-6)
		{
			tor[i] = -tor[i];
		}
	}
	return tor;
}

void DiscreteSpiral::Elevate()
{
	size_t num_v = mPoints.size();
	std::vector<ON_3dPoint> newpolygon;
	newpolygon.push_back(mPoints[0]);
	for (size_t i = 1; i < num_v; ++i)
	{
		newpolygon.push_back(mPoints[i - 1] * (i * 1.0 / num_v) + mPoints[i] * (1.0 - i * 1.0 / num_v));
	}
	newpolygon.push_back(mPoints.back());
	mPoints = newpolygon;
}

bool DiscreteSpiral::EndIteration()
{
	ComputeFrenetFrame();
	std::vector <double> cur = ComputeDiscreteCurvature();	
	double eps = 0.01;
	double error = 0.0;
	for (size_t i = 3; i < cur.size() - 1; ++i)
	{
		error = (std::max)(error, abs(cur[i] + cur[i - 2] - 2 * cur[i - 1]));
	}
	if (error >= eps)
	{
		return false;
	}
	double sum = 0.0;
	for (size_t i = 0; i < mPoints.size() - 1uLL; ++i)
	{
		sum += mPoints[i].DistanceTo(mPoints[i + 1uLL]);
	}
	double avglength = sum / (mPoints.size() - 1uLL);
	sum = 0.0;
	for (size_t i = 0; i < mPoints.size() - 1uLL; ++i)
	{
		sum += pow(mPoints[i].DistanceTo(mPoints[i + 1uLL]) - avglength, 2);
	}
	sum = sqrt(sum);
	sum /= avglength;
	if (sum >= 0.05)
	{
		return false;
	}
	return true;
}

ON_BezierCurve DiscreteSpiral::GetBezier() const
{
	ON_BezierCurve obc(3, false, mPoints.size());
	for (size_t i = 0uLL; i < mPoints.size(); ++i)
	{
		obc.SetCV(i, mPoints[i]);
	}
	return obc;
}

void SolveBoundaryTest(ONX_Model *model)
{
	double a = 10.0;
	double b = -0.5;
	double alpha = PI / 3;
	auto spiral_frenet = [a, b, alpha](double theta) -> FrenetFrame
	{
		ON_3dPoint p = ON_3dPoint(sin(alpha) * cos(theta), sin(alpha) * sin(theta), cos(alpha)) * a * exp(b * theta);
		ON_3dVector der = ON_3dVector(b * sin(alpha) * cos(theta) - sin(alpha) * sin(theta), b * sin(alpha) * sin(theta) + sin(alpha) * cos(theta), b * cos(alpha));
		der.Unitize();
		ON_3dVector derder = (a * b * exp(b * theta)) * ON_3dVector(b * sin(alpha) * cos(theta) - sin(alpha) * sin(theta), b * sin(alpha) * sin(theta) + sin(alpha) * cos(theta), b * cos(alpha)) + (a * exp(b * theta)) * ON_3dVector(-b * sin(alpha) * sin(theta) - sin(alpha) * cos(theta),
																																																									   b * sin(alpha) * cos(theta) - sin(alpha) * sin(theta), 0);
		ON_3dVector N = ON_3dVector::CrossProduct(der, derder);
		ON_3dVector B = ON_3dVector::CrossProduct(N, der);
		B.Unitize();
		return FrenetFrame(p, der, B);
	};

	std::vector<ON_3dPoint> vp = {ON_3dPoint(0, 0, 0), ON_3dPoint(5, 0, 0), ON_3dPoint(8, 2, 1), ON_3dPoint(10, 6, 3)};
	FrenetFrame f1 = spiral_frenet(0.0);
	FrenetFrame f2 = spiral_frenet(0.5 * PI);
	double Length = f1.GetPos().DistanceTo(f2.GetPos());
	double param[4] = {1.0 / Length, 0, 0.1, 0.0};
	std::vector<ON_3dPoint> p = ChiralityMath::SolveBoundary(f1, f2, 0.5, param, 20);
	const int lines_layer_index = model->AddLayer(L"lines_layer", ON_Color::Black);
	ChiralityAddLines(model, p, L"points", lines_layer_index);
}