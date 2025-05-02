#include "PiecesClothoid3d.h"
#include <assert.h>

PiecesClothoid3d::PiecesClothoid3d(std::vector <ON_3dPoint>& p, ON_3dVector st, ON_3dVector et, int subtimes)
{
	m_ControlPoint = p;
	m_PointArray = p;
	std::vector <ON_3dPoint> newp;
	std::vector <ON_3dVector> cur_gamma;
	int n = 0;
	st.Unitize();
	et.Unitize();
	for (int i = 0; i < subtimes; i++)
	{
		newp.clear();
		n = m_PointArray.size();
		cur_gamma = Compute_Curvature_Normal(st, et);
		for (int j = 0; j < n - 1; j++)
		{
			ON_3dVector nor = (cur_gamma[j] + cur_gamma[j + 1]) / 2;
			double cur = nor.Length();
			nor.Unitize();
			ON_3dVector al = m_PointArray[j + 1] - m_PointArray[j];
			al.Unitize();
			ON_3dVector beta = ON_3dVector::CrossProduct(nor, al);
			beta.Unitize();
			ON_3dPoint pnew = (m_PointArray[j] + m_PointArray[j + 1]) / 2 +
				beta * ((1 - sqrt(1 - cur * cur / 4 * m_PointArray[j].DistanceToSquared(m_PointArray[j + 1]))) / cur);
			newp.push_back(pnew);
		}

		std::vector <ON_3dPoint> oldP = m_PointArray;
		m_PointArray.clear();
		for (int j = 0; j < n - 1; j++)
		{
			m_PointArray.push_back(oldP[j]);
			m_PointArray.push_back(newp[j]);
		}
		m_PointArray.push_back(oldP.back());
	}
}

double PiecesClothoid3d::Length(int i, int j)
{
	assert(0 <= i && i <= j && j < m_PointArray.size());
	double sum = 0;
	for (int k = i; k < j; k++)
	{
		sum = sum + (m_PointArray[k + 1] - m_PointArray[k]).Length();
	}
	return sum;
}

void PiecesClothoid3d::Add_to_Model(ONX_Model* model, const wchar_t* name, ON_Color color)
{
	const int layer_index = model->AddLayer(name, color);
	ON_3dmObjectAttributes* attributes = new ON_3dmObjectAttributes();
	attributes->m_layer_index = layer_index;
	attributes->m_name = name;
	ON_3dPointArray sap;
	for (const auto it : m_PointArray)
	{
		sap.Append(it);
	}
	ON_Polyline pl = sap;
	ON_PolylineCurve* pc = new ON_PolylineCurve(pl);
	model->AddManagedModelGeometryComponent(pc, attributes);
}

std::vector <ON_3dVector> PiecesClothoid3d::Compute_Curvature_Normal(ON_3dVector st, ON_3dVector et)
{
	std::vector <ON_3dVector> re;
	ON_3dVector gamma = ON_3dVector::CrossProduct(-st, m_PointArray[1] - m_PointArray[0]);
	gamma.Unitize();
	double cur = 2 / (m_PointArray[1] - m_PointArray[0]).Length() * (ON_3dVector::CrossProduct
	(st, m_PointArray[1] - m_PointArray[0])).Length() / st.Length() /
		(m_PointArray[1] - m_PointArray[0]).Length();
	re.push_back(gamma * cur);
	int n = m_PointArray.size();
	for (int j = 1; j < n - 1; j++)
	{
		gamma = ON_3dVector::CrossProduct(m_PointArray[j - 1] - m_PointArray[j],
			m_PointArray[j + 1] - m_PointArray[j]);
		gamma.Unitize();
		cur = PiecesClothoid3d::ComputeDisceretCurvature(m_PointArray[j - 1], m_PointArray[j], m_PointArray[j + 1]);
		re.push_back(gamma * cur);
	}
	gamma = ON_3dVector::CrossProduct(m_PointArray[n - 2] - m_PointArray[n - 1], et);
	gamma.Unitize();
	cur = 2 / (m_PointArray[n - 1] - m_PointArray[n - 2]).Length() * (ON_3dVector::CrossProduct
	(et, m_PointArray[n - 1] - m_PointArray[n - 2])).Length() / st.Length() /
		(m_PointArray[n - 1] - m_PointArray[n - 2]).Length();
	re.push_back(gamma * cur);
	return re;
}

double PiecesClothoid3d::ComputeDisceretCurvature(ON_3dPoint p1, ON_3dPoint p2, ON_3dPoint p3)
{
	double a = ON_3dVector::CrossProduct(p1 - p2, p1 - p3).Length();
	return 2 * a / (p1 - p2).Length() / (p1 - p3).Length() / (p2 - p3).Length();
}

std::vector <ON_3dVector> PiecesClothoid3d::Compute_Tang(ON_3dVector st, ON_3dVector et)
{
	std::vector <ON_3dVector> re;
	st.Unitize();
	re.push_back(st);
	int n = m_PointArray.size();
	ON_3dVector v;
	for (int j = 1; j < n - 1; j++)
	{
		v = m_PointArray[j + 1] - m_PointArray[j - 1];
		v.Unitize();
		re.push_back(v);
	}
	et.Unitize();
	re.push_back(et);
	return re;
}