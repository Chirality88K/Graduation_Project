#include "EulerPolygon3D.h"
#include "write3dm.h"
#include "ChiralityMathTools.h"
#include <string>
#include "ChiralityLog.h"
#include <locale>
#include <codecvt>
extern const double PI;
PolarPoint3d::PolarPoint3d(double theta, double phi, double dis)
{
	mDistance = dis;
	while (theta < 0.0)
	{
		theta += 2 * PI;
	}
	while (theta >= 2 * PI)
	{
		theta -= 2 * PI;
	}
	mTheta = theta;
	phi = (std::min)(PI / 2, phi);
	phi = (std::max)(-PI / 2, phi);
	mPhi = phi;
}

PolarPoint3d::PolarPoint3d(const ON_3dPoint &p)
{
	mDistance = p.DistanceTo(ON_3dPoint(0, 0, 0));
	if (mDistance < 1e-8)
	{
		mTheta = 0.0;
		mPhi = 0.0;
		return;
	}
	double horizontal_distance = sqrt(p.x * p.x + p.y * p.y);
	if (horizontal_distance < 1e-8)
	{
		mTheta = 0.0;
		mPhi = p.z > 0 ? PI / 2 : -PI / 2;
		return;
	}
	double cos_theta = p.x / horizontal_distance;
	cos_theta = (std::max)(-1.0, cos_theta);
	cos_theta = (std::min)(1.0, cos_theta);
	mTheta = p.y >= 0.0 ? acos(cos_theta) : -acos(cos_theta);
	mPhi = asin(p.z / mDistance);
}

ON_3dPoint PolarPoint3d::CartesianCoordinates() const
{
	ON_3dPoint p;
	p.x = mDistance * cos(mPhi) * cos(mTheta);
	p.y = mDistance * cos(mPhi) * sin(mTheta);
	p.z = mDistance * sin(mPhi);
	return p;
}

EulerPolygon3D::EulerPolygon3D(ON_3dPoint PS, ON_3dPoint PE, ON_3dVector vs, ON_3dVector ve, CurveType ct)
{
	ON_Plane old_plane(ON_3dPoint::Origin, vs, PE - PS);
	ON_Plane new_plane(ON_3dPoint::Origin, ON_3dVector::XAxis, ON_3dVector::YAxis);
	ON_Xform rot;
	rot.Rotation(old_plane, new_plane);
	ON_3dVector new_ve = ve;
	ON_3dPoint new_pe = PE - PS;
	new_ve.Transform(rot);
	new_pe.Transform(rot);
	if (ct == Bezier)
	{
		BuildUpToBezier(new_pe, new_ve);
	}
	else if (ct == B_spline)
	{
		BuildUpToB_Spline(new_pe, new_ve);
	}
	else
	{
		std::string str = "wrong type!";
		CHIRALITY_ERROR(str);
	}
	ON_Xform inv_rot = rot.Inverse();
	mCurve.Transform(inv_rot);
	mCurve.Translate(ON_3dVector(PS));
	for (ON_3dPoint &p : mDiscretePolygon)
	{
		p.Transform(inv_rot);
		p += PS;
	}
}

void EulerPolygon3D::EulerPolygonTest_ForConicSpiral(ONX_Model *model)
{
	double a = 10.0;
	double b = -0.5;
	double alpha = PI / 3;
	auto spiral = [a, b, alpha](double theta) -> ON_3dPoint
	{
		return ON_3dPoint(sin(alpha) * cos(theta), sin(alpha) * sin(theta), cos(alpha)) * a * exp(b * theta);
	};
	auto spiral_tan = [b, alpha](double theta) -> ON_3dVector
	{
		ON_3dVector v(b * sin(alpha) * cos(theta) - sin(alpha) * sin(theta), b * sin(alpha) * sin(theta) + sin(alpha) * cos(theta), b * cos(alpha));
		v.Unitize();
		return v;
	};

	ON_Color color[7] = {ON_Color::SaturatedRed, ON_Color(255, 128, 0), ON_Color::SaturatedYellow,
						 ON_Color::SaturatedGreen, ON_Color::SaturatedCyan, ON_Color::SaturatedBlue, ON_Color(76, 0, 153)};
	std::string color_name[7] = {"Red", "Orange", "Yellow", "Green", "Cyan", "Blue", "Purple"};
	for (int i = 0; i < 7; ++i)
	{
		double t0 = 5.0 / 7.0 * i;
		double t1 = t0 + 5.0 / 7.0;
		EulerPolygon3D ep_bezier(spiral(t0), spiral(t1), spiral_tan(t0), spiral_tan(t1), EulerPolygon3D::CurveType::Bezier);
		ON_NurbsCurve onc = ep_bezier.GetCurve();
		std::wstring_convert<std::codecvt_utf8<wchar_t>> conv;
		const int layer_index = model->AddLayer(conv.from_bytes(color_name[i]).c_str(), color[i]);
		ChiralityAddNurbsCurve(model, onc, conv.from_bytes(color_name[i] + " bezier curve").c_str(), layer_index);
		ChiralityDebugInfo(onc, color_name[i] + " bezier_debug");
		ChiralityDebugforR(onc, color_name[i] + " bezier_for_R");
	}
}

void EulerPolygon3D::EulerPolygonTest_ForSphereSpiral(ONX_Model *model)
{
	double R = 10.0;
	double ratio = 0.03;
	auto sphere_spiral = [R, ratio](double theta) -> ON_3dPoint
	{
		double phi = ratio * theta;
		return ON_3dPoint(cos(theta) * cos(phi), cos(theta) * sin(phi), sin(theta)) * R;
	};
	auto sphere_spiral_tan = [ratio](double theta) -> ON_3dVector
	{
		double phi = ratio * theta;
		return ON_3dVector(-sin(theta) * cos(phi), -sin(theta) * sin(phi), cos(theta));
	};
	double theta0 = 0.0;
	double theta1 = 0.4 * PI;
	EulerPolygon3D ep_bezier(sphere_spiral(theta0), sphere_spiral(theta1), sphere_spiral_tan(theta0), sphere_spiral_tan(theta1), EulerPolygon3D::CurveType::Bezier);
	ON_NurbsCurve onc = ep_bezier.GetCurve();
	const int layer_index = model->AddLayer(L"EulerPolygonTest_ForSphereSpiral", ON_Color::SaturatedBlue);
	ChiralityAddNurbsCurve(model, onc, L"EulerPolygonTest_ForSphereSpiral", layer_index);
	ChiralityDebugInfo(onc, "EulerPolygonTest_ForSphereSpiral bezier_debug");
	ChiralityDebugforR(onc, "EulerPolygonTest_ForSphereSpiral bezier_for_R");
}

std::vector<double> EulerPolygon3D::ComputeDeltaTheta() const
{
	std::vector<double> result(mDiscretePolygon.size(), 0.0);
	for (size_t i = 1; i < mDiscretePolygon.size() - 1; ++i)
	{
		ON_3dVector edge1 = mDiscretePolygon[i] - mDiscretePolygon[i - 1];
		ON_3dVector edge2 = mDiscretePolygon[i + 1] - mDiscretePolygon[i];
		edge1.z = 0;
		edge2.z = 0;
		double product = ON_3dVector::DotProduct(edge1, edge2) / edge1.Length() / edge2.Length();
		product = (std::min)(1.0, product);
		product = (std::max)(product, -1.0);
		double angle = acos(product);
		if (edge1.x * edge2.y - edge1.y * edge2.x < 0)
		{
			angle = -angle;
		}
		result[i] = angle;
	}
	return result;
}

std::vector<double> EulerPolygon3D::ComputeDeltaPhi() const
{
	std::vector<double> result(mDiscretePolygon.size(), 0.0);
	for (size_t i = 1; i < mDiscretePolygon.size() - 1; ++i)
	{
		ON_3dPoint edge1 = mDiscretePolygon[i] - mDiscretePolygon[i - 1];
		ON_3dPoint edge2 = mDiscretePolygon[i + 1] - mDiscretePolygon[i];
		PolarPoint3d v1(edge1);
		PolarPoint3d v2(edge2);
		result[i] = v2.mPhi - v1.mPhi;
	}
	return result;
}

std::vector<double> EulerPolygon3D::ComputeLength() const
{
	std::vector<double> result(mDiscretePolygon.size() - 1, 0.0);
	for (size_t i = 0; i < mDiscretePolygon.size() - 1; ++i)
	{
		ON_3dVector edge2 = mDiscretePolygon[i + 1] - mDiscretePolygon[i];
		result[i] = edge2.Length();
	}
	return result;
}

std::vector<double> EulerPolygon3D::ComputeProjectionLength() const
{
	std::vector<double> result(mDiscretePolygon.size() - 1, 0.0);
	for (size_t i = 0; i < mDiscretePolygon.size() - 1; ++i)
	{
		ON_3dVector edge2 = mDiscretePolygon[i + 1] - mDiscretePolygon[i];
		edge2.z = 0;
		result[i] = edge2.Length();
	}
	return result;
}

void EulerPolygon3D::Elevate()
{
	size_t num_v = mDiscretePolygon.size();
	std::vector<ON_3dPoint> newpolygon;
	newpolygon.push_back(mDiscretePolygon[0]);
	for (size_t i = 1; i < num_v; ++i)
	{
		newpolygon.push_back(mDiscretePolygon[i - 1] * (i * 1.0 / num_v) + mDiscretePolygon[i] * (1.0 - i * 1.0 / num_v));
	}
	newpolygon.push_back(mDiscretePolygon.back());
	mDiscretePolygon = newpolygon;
}

void EulerPolygon3D::SmoothingToBezier()
{
	int num_v = mDiscretePolygon.size();
	ON_3dVector vs = mDiscretePolygon[1] - mDiscretePolygon[0];
	ON_3dVector ve = mDiscretePolygon[num_v - 1] - mDiscretePolygon[num_v - 2];
	vs.Unitize();
	ve.Unitize();
	std::vector<double> Angle_Theta = ComputeDeltaTheta();
	std::vector<double> Angle_Phi = ComputeDeltaPhi();
	int max_count = 10000;
	for (int i = 0; i < num_v; ++i)
	{
		if (Angle_Theta[i] > PI / 2 || Angle_Theta[i] < -PI / 2)
		{
			max_count = 2;
		}
	}
	int s_count = 1;
	double thetaDD = 10.0;
	double lengthavg = 0.0;
	double lengthbound = 2 * (mDiscretePolygon[0] - mDiscretePolygon.back()).Length();
	while (s_count < max_count && thetaDD > 1e-6 && lengthavg < lengthbound)
	{
		std::vector<double> new_angle_theta = Angle_Theta;
		std::vector<double> new_angle_phi = Angle_Phi;

		for (int i = 2; i < num_v - 2; i++)
		{
			new_angle_theta[i] = (new_angle_theta[i - 1] + new_angle_theta[i] + new_angle_theta[i + 1]) / 3;
			new_angle_phi[i] = (new_angle_phi[i - 1] + new_angle_phi[i] + new_angle_phi[i + 1]) / 3;
			ON_3dPoint p0 = mDiscretePolygon[i - 1];
			ON_3dPoint pm = mDiscretePolygon[i + 1];
			p0.z = 0;
			pm.z = 0;
			ON_3dPoint re = (pm - p0) / 2;
			re.Set(-re.y, re.x, 0);
			re = (p0 + pm) / 2 - re * tan(new_angle_theta[i] / 2);
			double L = re.DistanceTo(p0);
			double h = mDiscretePolygon[i + 1].z - mDiscretePolygon[i - 1].z;
			double A = tan(new_angle_phi[i]);
			if (abs(A) < 1e-8)
			{
				re.z = h / 2 + mDiscretePolygon[i - 1].z;
				mDiscretePolygon[i] = re;
				continue;
			}
			double DELTA = 4 * L * L + 4 * A * A * L * L + A * A * h * h;
			double z1 = (2 * L + A * h + sqrt(DELTA)) / (2 * A);
			double z2 = (2 * L + A * h - sqrt(DELTA)) / (2 * A);
			re.z = (abs(z1 - h / 2) > abs(z2 - h / 2)) ? z2 : z1;
			re.z += mDiscretePolygon[i - 1].z;
			mDiscretePolygon[i] = re;
		}
		std::vector<double> Length = ComputeLength();
		double sum = 0;
		for (int i = 0; i < num_v - 1; i++)
		{
			sum = sum + Length[i];
		}
		lengthavg = sum / (num_v - 1);
		mDiscretePolygon[1] = mDiscretePolygon[0] + vs * lengthavg;
		mDiscretePolygon[num_v - 2] = mDiscretePolygon[num_v - 1] - ve * lengthavg;
		thetaDD = 0;
		Angle_Theta = ComputeDeltaTheta();
		Angle_Phi = ComputeDeltaPhi();
		for (int i = 2; i < num_v - 2; ++i)
		{
			thetaDD = (std::max)(thetaDD, abs(Angle_Theta[i] * 2 - Angle_Theta[i - 1] - Angle_Theta[i + 1]));
		}
		for (int i = 2; i < num_v - 2; ++i)
		{
			thetaDD = (std::max)(thetaDD, abs(Angle_Phi[i] * 2 - Angle_Phi[i - 1] - Angle_Phi[i + 1]));
		}
		s_count++;
	}
}

void EulerPolygon3D::SmoothingToB_Spline(ON_3dPoint PE, ON_3dVector ve)
{
	int num_v = mDiscretePolygon.size();
	ON_3dVector vs(1.0, 0.0, 0.0);
	vs.Unitize();
	ve.Unitize();
	ON_3dPoint ps = ON_3dPoint::Origin;
	ON_3dPoint pe = PE;
	ON_3dVector ns(0, 1, 0);
	ON_3dVector ne(-ve.y, ve.x, 0);
	ne.Unitize();
	std::vector<double> Angle_Theta = ComputeDeltaTheta();
	std::vector<double> Angle_Phi = ComputeDeltaPhi();
	int max_count = 10000;
	for (int i = 0; i < num_v; ++i)
	{
		if (Angle_Theta[i] > PI / 2 || Angle_Theta[i] < -PI / 2)
		{
			max_count = 2;
		}
	}
	int s_count = 1;
	double thetaDD = 10.0;
	double lengthavg = 0.0;
	double lengthbound = 2 * (mDiscretePolygon[0] - mDiscretePolygon.back()).Length();
	while (s_count < max_count && thetaDD > 1e-6 && lengthavg < lengthbound)
	{
		std::vector<double> new_angle_theta = Angle_Theta;
		std::vector<double> new_angle_phi = Angle_Phi;

		for (int i = 2; i < num_v - 2; i++)
		{
			new_angle_theta[i] = (new_angle_theta[i - 1] + new_angle_theta[i] + new_angle_theta[i + 1]) / 3;
			new_angle_phi[i] = (new_angle_phi[i - 1] + new_angle_phi[i] + new_angle_phi[i + 1]) / 3;
			ON_3dPoint p0 = mDiscretePolygon[i - 1];
			ON_3dPoint pm = mDiscretePolygon[i + 1];
			p0.z = 0;
			pm.z = 0;
			ON_3dPoint re = (pm - p0) / 2;
			re.Set(-re.y, re.x, 0);
			re = (p0 + pm) / 2 - re * tan(new_angle_theta[i] / 2);
			double L = re.DistanceTo(p0);
			double h = mDiscretePolygon[i + 1].z - mDiscretePolygon[i - 1].z;
			double A = tan(new_angle_phi[i]);
			double DELTA = 4 * L * L + 4 * A * A * L * L + A * A * h * h;
			double z1 = (2 * L + A * h + sqrt(DELTA)) / (2 * A);
			double z2 = (2 * L + A * h - sqrt(DELTA)) / (2 * A);
			re.z = (abs(z1 - h / 2) > abs(z2 - h / 2)) ? z2 : z1;
			re.z += mDiscretePolygon[i - 1].z;
			mDiscretePolygon[i] = re;
		}
		std::vector<double> projectLength = ComputeProjectionLength();
		double sum = 0;
		for (int i = 0; i < num_v - 1; i++)
		{
			sum = sum + projectLength[i];
		}
		lengthavg = sum / (num_v - 1);
		double theta1 = 2 * new_angle_theta[2] - new_angle_theta[3];
		mDiscretePolygon[0] = mDiscretePolygon[2] - vs * (2 * lengthavg * cos(theta1 / 2));
		mDiscretePolygon[0].z = mDiscretePolygon[2].z - (vs.z / sqrt(vs.x * vs.x + vs.y * vs.y)) *
															sqrt(pow(mDiscretePolygon[2].x - mDiscretePolygon[0].x, 2) + pow(mDiscretePolygon[2].y - mDiscretePolygon[0].y, 2));
		mDiscretePolygon[1] = (6 * ps - mDiscretePolygon[0] - mDiscretePolygon[2]) / 4;

		double thetan_1 = 2 * new_angle_theta[num_v - 3] - new_angle_theta[num_v - 4];
		mDiscretePolygon[num_v - 1] = mDiscretePolygon[num_v - 3] + ve * (2 * lengthavg * cos(thetan_1 / 2));
		mDiscretePolygon[num_v - 1].z = mDiscretePolygon[num_v - 3].z + (ve.z / sqrt(ve.x * ve.x + ve.y * ve.y)) *
																			sqrt(pow(mDiscretePolygon[num_v - 3].x - mDiscretePolygon[num_v - 1].x, 2) + pow(mDiscretePolygon[num_v - 3].y - mDiscretePolygon[num_v - 1].y, 2));
		mDiscretePolygon[num_v - 2] = (6 * pe - mDiscretePolygon[num_v - 1] - mDiscretePolygon[num_v - 3]) / 4;

		thetaDD = 0;
		Angle_Theta = ComputeDeltaTheta();
		Angle_Phi = ComputeDeltaPhi();
		for (int i = 2; i < num_v - 2; ++i)
		{
			thetaDD = (std::max)(thetaDD, abs(Angle_Theta[i] * 2 - Angle_Theta[i - 1] - Angle_Theta[i + 1]));
		}
		for (int i = 2; i < num_v - 2; ++i)
		{
			thetaDD = (std::max)(thetaDD, abs(Angle_Phi[i] * 2 - Angle_Phi[i - 1] - Angle_Phi[i + 1]));
		}
		std::vector<double> Length = ComputeLength();
		sum = 0;
		for (int i = 0; i < num_v - 1; i++)
		{
			sum = sum + Length[i];
		}
		lengthavg = sum / (num_v - 1);
		s_count++;
	}
}

void EulerPolygon3D::BuildUpToBezier(ON_3dPoint PE, ON_3dVector ve)
{
	// Init
	ve.Unitize();
	double L = PE.DistanceTo(ON_3dPoint::Origin) / 3.0;
	mDiscretePolygon.resize(4);
	mDiscretePolygon[0] = ON_3dPoint::Origin;
	mDiscretePolygon[1] = ON_3dPoint(L, 0, 0);
	mDiscretePolygon[2] = PE - ve * L;
	mDiscretePolygon[3] = PE;
	SmoothingToBezier();
	// 迭代
	for (int i = 0; i < ITERATIONTIMES; ++i)
	{
		Elevate();
		SmoothingToBezier();
	}
	ON_BezierCurve obc;
	obc.Create(3, false, mDiscretePolygon.size());
	for (int i = 0; i < mDiscretePolygon.size(); ++i)
	{
		obc.SetCV(i, mDiscretePolygon[i]);
	}
	obc.GetNurbForm(mCurve);
}

void EulerPolygon3D::BuildUpToB_Spline(ON_3dPoint PE, ON_3dVector ve)
{
	// Init
	ve.Unitize();
	ON_NurbsCurve onc = ChiralityMath::UniformG1(ON_3dPoint::Origin, PE, ON_3dVector::XAxis, ve);
	mDiscretePolygon.resize(4);
	ON_3dPoint p;
	for (int i = 0; i < 4; ++i)
	{
		onc.GetCV(i, p);
		mDiscretePolygon[i] = p;
	}
	SmoothingToB_Spline(PE, ve);
	// 迭代
	for (int i = 0; i < ITERATIONTIMES; ++i)
	{
		Elevate();
		SmoothingToB_Spline(PE, ve);
	}
	onc.Create(3, false, 4, mDiscretePolygon.size());
	for (int i = 0; i < mDiscretePolygon.size(); ++i)
	{
		onc.SetCV(i, mDiscretePolygon[i]);
	}
	for (int i = 0; i < onc.KnotCount(); ++i)
	{
		onc.SetKnot(i, i + 1);
	}
	mCurve = onc;
}
