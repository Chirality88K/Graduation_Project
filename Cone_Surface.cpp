#include "Cone_Surface.h"

Cone_Surface::Cone_Surface(double angle, double r1, double r2, double theta)
{
	ON_3dPoint P0(r1, 0, r1 / tan(angle));
	ON_3dPoint Q0(r2, 0, r2 / tan(angle));
	ON_3dPoint CP(0, 0, P0.z);
	ON_3dPoint CQ(0, 0, Q0.z);
	ON_Xform Rot1;
	Rot1.Rotation(theta, ON_3dVector(0, 0, 1), CP);
	ON_Xform Rot2;
	Rot2.Rotation(theta, ON_3dVector(0, 0, 1), CQ);
	ON_3dPoint P2 = P0; ON_3dPoint Q2 = Q0;
	P2.Transform(Rot1);
	Q2.Transform(Rot2);
	ON_3dPoint MidP = (P0 + P2) / 2;
	ON_3dPoint MidQ = (Q0 + Q2) / 2;
	ON_3dPoint P1 = CP + (MidP - CP) * (2 / (1 + cos(theta)));
	ON_3dPoint Q1 = CQ + (MidQ - CQ) * (2 / (1 + cos(theta)));

	this->Create(3, true, 3, 2, 3, 2);
	this->SetKnot(0, 0, 0);
	this->SetKnot(0, 1, 0);
	this->SetKnot(0, 2, 1);
	this->SetKnot(0, 3, 1);
	this->SetKnot(1, 0, 0);
	this->SetKnot(1, 1, (P0 - Q0).Length());
	this->SetCV(0, 0, P0);
	this->SetCV(1, 0, P1 * cos(theta / 2));
	this->SetCV(2, 0, P2);
	this->SetCV(0, 1, Q0);
	this->SetCV(1, 1, Q1 * cos(theta / 2));
	this->SetCV(2, 1, Q2);
	this->SetWeight(1, 0, cos(theta / 2));
	this->SetWeight(1, 1, cos(theta / 2));
}

void Cone_Surface::Add_Cone_to_Model(ONX_Model* model, const wchar_t* name, ON_Color color)
{
	const int layer_index = model->AddLayer(name, color);
	ON_3dmObjectAttributes* attributes = new ON_3dmObjectAttributes();
	attributes->m_layer_index = layer_index;
	attributes->m_name = name;
	model->AddManagedModelGeometryComponent(this, attributes);
}