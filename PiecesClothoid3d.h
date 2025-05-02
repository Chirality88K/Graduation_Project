#pragma once
#include "opennurbs.h"

class PiecesClothoid3d
{
private:
	std::vector <ON_3dPoint> m_ControlPoint;
	std::vector <ON_3dPoint> m_PointArray;
	double Length(int i, int j);
	std::vector <ON_3dVector> Compute_Curvature_Normal(ON_3dVector st, ON_3dVector et);
	std::vector <ON_3dVector> Compute_Tang(ON_3dVector st, ON_3dVector et);
public:
	PiecesClothoid3d(std::vector <ON_3dPoint>& p, ON_3dVector, ON_3dVector, int subtimes = 5);
	void Add_to_Model(ONX_Model* model, const wchar_t* name, ON_Color color = ON_Color::Black);
	static double ComputeDisceretCurvature(ON_3dPoint, ON_3dPoint, ON_3dPoint);
};

