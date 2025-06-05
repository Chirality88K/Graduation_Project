#include "EulerPolygon3D.h"
#include "write3dm.h"
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

EulerPolygon3D::EulerPolygon3D(ON_3dPoint PE, ON_3dVector ve)
{
	BuildUp(PE, ve);
}

EulerPolygon3D::EulerPolygon3D(ON_3dPoint PS, ON_3dPoint PE, ON_3dVector vs, ON_3dVector ve)
{
	ON_Plane old_plane(ON_3dPoint::Origin, vs, PE - PS);
	ON_Plane new_plane(ON_3dPoint::Origin, ON_3dVector::XAxis, ON_3dVector::YAxis);
	ON_Xform rot;
	rot.Rotation(old_plane, new_plane);
	ON_3dVector new_ve = ve;
	ON_3dPoint new_pe = PE - PS;
	new_ve.Transform(rot);
	new_pe.Transform(rot);
	BuildUp(new_pe, new_ve);
	ON_Xform inv_rot = rot.Inverse();
	for (ON_3dPoint &p : mDiscretePolygon)
	{
		p.Transform(inv_rot);
		p += PS;
	}
}

void EulerPolygon3D::EulerPolygonTest(ONX_Model *model)
{
	double a = 10.0;
	double b = 0.5;
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

	double t0 = -0.5 * PI;
	double t1 = PI * 0.5;
	EulerPolygon3D ep3d(spiral(t0), spiral(t1), spiral_tan(t0), spiral_tan(t1));
	// EulerPolygon3D ep3d(ON_3dPoint(10, 10, 0), ON_3dVector(-1, 3, 1));
	ON_NurbsCurve onc = ep3d.ToBezier();
	const int layer_index = model->AddLayer(L"test_layer", ON_Color::SaturatedBlue);
	ChiralityAddNurbsCurve(model, onc, L"test", layer_index);
	ChiralityDebugInfo(onc);
	ChiralityDebugforR(onc);
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

void EulerPolygon3D::Smoothing()
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
			double DELTA = 4 * L * L + 4 * A * A * L * L + A * A * h * h;
			double z1 = (2 * L + A * h + sqrt(DELTA)) / (2 * A);
			double z2 = (2 * L + A * h - sqrt(DELTA)) / (2 * A);
			re.z = (abs(z1 - h / 2) > abs(z2 - h / 2)) ? z2 : z1;
			re.z += mDiscretePolygon[i - 1].z;
			mDiscretePolygon[i] = re;
		}
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
		double sum = 0;
		for (int i = 0; i < num_v - 1; i++)
		{
			sum = sum + Length[i];
		}
		lengthavg = sum / (num_v - 1);
		mDiscretePolygon[1] = mDiscretePolygon[0] + vs * lengthavg;
		mDiscretePolygon[num_v - 2] = mDiscretePolygon[num_v - 1] - ve * lengthavg;
		s_count++;
	}
}

ON_NurbsCurve EulerPolygon3D::ToBezier() const
{
	ON_BezierCurve obc;
	obc.Create(3, false, mDiscretePolygon.size());
	for (int i = 0; i < mDiscretePolygon.size(); ++i)
	{
		obc.SetCV(i, mDiscretePolygon[i]);
	}
	ON_NurbsCurve onc;
	obc.GetNurbForm(onc);
	return onc;
}

void EulerPolygon3D::BuildUp(ON_3dPoint PE, ON_3dVector ve)
{
	// Init
	ve.Unitize();
	double L = PE.DistanceTo(ON_3dPoint::Origin) / 3.0;
	mDiscretePolygon.resize(4);
	mDiscretePolygon[0] = ON_3dPoint::Origin;
	mDiscretePolygon[1] = ON_3dPoint(L, 0, 0);
	mDiscretePolygon[2] = PE - ve * L;
	mDiscretePolygon[3] = PE;
	Smoothing();
	// 迭代
	const int times = 30;
	for (int i = 0; i < times; ++i)
	{
		Elevate();
		Smoothing();
	}
}
