#include "ConicSpiral.h"
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

ON_3dPoint PolarPoint3d::CartesianCoordinates()
{
	ON_3dPoint p;
	p.x = mDistance * cos(mPhi) * cos(mTheta);
	p.y = mDistance * cos(mPhi) * sin(mTheta);
	p.z = mDistance * sin(mPhi);
	return p;
}

ConicSpiral::ConicSpiral(PolarPoint3d PS, PolarPoint3d PE, ON_3dVector vs, ON_3dVector ve, int num)
{
	mNum_Points = (std::max)(10, num);
	double theta0 = PS.mTheta;
	double theta1 = PE.mTheta;
	if (theta1 < theta0)
	{
		theta1 += 2 * PI;
	}
	std::vector<double> all_theta(mNum_Points, 0);
	std::vector<double> all_distance_to_Z(mNum_Points, 0);
	std::vector<double> all_z_coor(mNum_Points, 0);
	for (int i = 0; i < mNum_Points; ++i)
	{
		all_theta[i] = theta0 * (1 - double(i) / double(mNum_Points - 1)) + theta1 * double(i) / double(mNum_Points - 1);
	}
	double dis0 = PS.mDistance * cos(PS.mPhi);
	double dis1 = PE.mDistance * cos(PE.mPhi);
	double B = log(dis1 / dis0) / (theta1 - theta0);
	double A = (theta1 * log(dis0) - theta0 * log(dis1)) / (theta1 - theta0);
	for (int i = 0; i < mNum_Points; ++i)
	{
		all_distance_to_Z[i] = A * exp(B * all_theta[i]);
	}
	double z0 = PS.mDistance * sin(PS.mPhi);
	double z1 = PE.mDistance * sin(PE.mPhi);
	double D = log(z1 / z0) / (theta1 - theta0);
	double C = (theta1 * log(z0) - theta0 * log(z1)) / (theta1 - theta0);
	for (int i = 0; i < mNum_Points; ++i)
	{
		all_z_coor[i] = C * exp(D * all_theta[i]);
	}
	for (int i = 0; i < mNum_Points; ++i)
	{
		mDiscretePolygon.push_back(ON_3dPoint(all_distance_to_Z[i] * cos(all_theta[i]), all_distance_to_Z[i] * sin(all_theta[i]), all_z_coor[i]));
	}
}

void ConicSpiral::ConicSpiralTest(ONX_Model *model)
{
}