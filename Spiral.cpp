#include "Spiral.h"
#include <Eigen/Dense>
#include <vector>

void Rational_Bezier_G2(ON_BezierCurve *obc, ON_3dPoint p0, ON_3dPoint p3, ON_3dVector v0, ON_3dVector v2, double k0, double k3)
{
	ON_3dPoint p1 = p0 + v0;
	ON_3dPoint p2 = p3 - v2;
	double c0 = (ON_3dVector::CrossProduct(p1 - p0, p2 - p0)).Length() / (pow((p1 - p0).Length(), 3));
	double c1 = (ON_3dVector::CrossProduct(p3 - p2, p3 - p1)).Length() / (pow((p3 - p2).Length(), 3));
	double w1 = 2.0 / 3 * pow((c0 * c0 * c1 / (k0 * k0 * k3)), 1.0 / 3);
	double w2 = 2.0 / 3 * pow((c0 * c1 * c1 / (k0 * k3 * k3)), 1.0 / 3);
	obc->Create(3, true, 4);
	obc->SetCV(0, p0);
	obc->SetCV(1, p1 * w1);
	obc->SetCV(2, p2 * w2);
	obc->SetCV(3, p3);
	obc->SetWeight(1, w1);
	obc->SetWeight(2, w2);
}

void ThreeDegreeBsplineInterplate_Tan(ON_NurbsCurve &onc, std::vector<ON_3dPoint> Q, double *knot, ON_3dVector v0, ON_3dVector vn)
{
	int K = Q.size() - 1;
	int L = K + 2;
	ON_NurbsCurve solve_N;
	solve_N.Create(1, false, 4, L + 1);
	solve_N.SetKnot(0, knot[0]);
	solve_N.SetKnot(1, knot[0]);
	for (int i = 0; i <= K; i++)
	{
		solve_N.SetKnot(i + 2, knot[i]);
	}
	solve_N.SetKnot(K + 3, knot[K]);
	solve_N.SetKnot(K + 4, knot[K]);
	Eigen::MatrixXd N = Eigen::MatrixXd::Zero(K + 1, K + 1);
	Eigen::MatrixXd QQ = Eigen::MatrixXd::Zero(K + 1, 3);
	ON_3dPoint One = ON_3dPoint(1, 1, 1);
	ON_3dPoint Zero = ON_3dPoint(0, 0, 0);
	for (int j = 1; j <= K - 1; j++)
	{
		solve_N.ZeroCVs();
		solve_N.SetCV(j, One);
		N(j, j - 1) = solve_N.PointAt(knot[j]).x;
		solve_N.SetCV(j, Zero);
		solve_N.SetCV(j + 1, One);
		N(j, j) = solve_N.PointAt(knot[j]).x;
		solve_N.SetCV(j + 1, Zero);
		solve_N.SetCV(j + 2, One);
		N(j, j + 1) = solve_N.PointAt(knot[j]).x;
		QQ(j, 0) = Q[j].x;
		QQ(j, 1) = Q[j].y;
		QQ(j, 2) = Q[j].z;
	}
	N(0, 0) = 3 / (knot[1] - knot[0]);
	N(K, K) = 3 / (knot[K] - knot[K - 1]);
	ON_3dPoint Q0 = 3 * Q[0] / (knot[1] - knot[0]) + v0;
	QQ(0, 0) = Q0.x;
	QQ(0, 1) = Q0.y;
	QQ(0, 2) = Q0.z;
	ON_3dPoint QK = 3 * Q[K] / (knot[K] - knot[K - 1]) - vn;
	QQ(K, 0) = QK.x;
	QQ(K, 1) = QK.y;
	QQ(K, 2) = QK.z;

	Eigen::MatrixXd P = Eigen::MatrixXd::Zero(K + 1, 3);
	P = N.partialPivLu().solve(QQ);

	onc.Create(3, false, 4, L + 1);
	for (int i = 0; i < L + 3; i++)
	{
		onc.SetKnot(i, solve_N.Knot(i));
	}
	onc.SetCV(0, Q[0]);
	for (int i = 1; i < L; i++)
	{
		onc.SetCV(i, ON_3dPoint(P(i - 1, 0), P(i - 1, 1), P(i - 1, 2)));
	}
	onc.SetCV(L, Q[K]);
}

void ThreeDegreeBsplineInterplate_Tan(ON_NurbsCurve &onc, std::vector<ON_3dPoint> Q, std::vector<double> knot, ON_3dVector v0, ON_3dVector vn)
{
	int s = knot.size();
	double *arra = new double[s];
	int i = 0;
	for (auto it : knot)
	{
		*(arra + i) = it;
		i++;
	}
	ThreeDegreeBsplineInterplate_Tan(onc, Q, arra, v0, vn);
	delete[] arra;
}

void FiveDegreeBsplineInterplate_Tan_and_Cur(ON_NurbsCurve &onc, std::vector<ON_3dPoint> Q, std::vector<double> knot, ON_3dVector der1_on_knot0, ON_3dVector der1_on_knotK, ON_3dVector der2_on_knot0, ON_3dVector der2_on_knotK)
{
	int K = Q.size() - 1;
	ON_NurbsCurve solve_N;
	solve_N.Create(1, false, 6, K + 5);
	solve_N.SetKnot(0, knot[0]);
	solve_N.SetKnot(1, knot[0]);
	solve_N.SetKnot(2, knot[0]);
	solve_N.SetKnot(3, knot[0]);
	for (int i = 0; i <= K; i++)
	{
		solve_N.SetKnot(i + 4, knot[i]);
	}
	solve_N.SetKnot(K + 5, knot[K]);
	solve_N.SetKnot(K + 6, knot[K]);
	solve_N.SetKnot(K + 7, knot[K]);
	solve_N.SetKnot(K + 8, knot[K]);
	Eigen::MatrixXd N = Eigen::MatrixXd::Zero(K + 3, K + 3);
	Eigen::MatrixXd QQ = Eigen::MatrixXd::Zero(K + 3, 3);
	ON_3dPoint One = ON_3dPoint(1, 1, 1);
	ON_3dPoint Zero = ON_3dPoint(0, 0, 0);
	for (int j = 1; j <= K - 1; j++)
	{
		for (int r = -1; r < 4; r++)
		{
			solve_N.ZeroCVs();
			solve_N.SetCV(j + r + 1, One);
			// ON_3dPoint testP = solve_N.PointAt(knot[j]);
			N(j, j + r) = solve_N.PointAt(knot[j]).x;
		}
		QQ(j, 0) = Q[j].x;
		QQ(j, 1) = Q[j].y;
		QQ(j, 2) = Q[j].z;
	}
	N(0, 0) = 1;
	N(K, K + 2) = 1;
	N(K + 1, 0) = -1 / (knot[2] - knot[0]) - 1 / (knot[1] - knot[0]);
	N(K + 1, 1) = 1 / (knot[2] - knot[0]);
	N(K + 2, K + 1) = 1 / (knot[K] - knot[K - 2]);
	N(K + 2, K + 2) = -1 / (knot[K] - knot[K - 1]) - 1 / (knot[K] - knot[K - 2]);
	ON_3dPoint q0 = Q[0] + (knot[1] - knot[0]) / 5 * der1_on_knot0;
	ON_3dPoint qK = Q[K] - (knot[K] - knot[K - 1]) / 5 * der1_on_knotK;
	ON_3dPoint qK_1 = (knot[1] - knot[0]) / 20 * der2_on_knot0 - 1 / (knot[1] - knot[0]) * Q[0];
	ON_3dPoint qK_2 = (knot[K] - knot[K - 1]) / 20 * der2_on_knotK - 1 / (knot[K] - knot[K - 1]) * Q[K];
	QQ(0, 0) = q0.x;
	QQ(0, 1) = q0.y;
	QQ(0, 2) = q0.z;
	QQ(K, 0) = qK.x;
	QQ(K, 1) = qK.y;
	QQ(K, 2) = qK.z;
	QQ(K + 1, 0) = qK_1.x;
	QQ(K + 1, 1) = qK_1.y;
	QQ(K + 1, 2) = qK_1.z;
	QQ(K + 2, 0) = qK_2.x;
	QQ(K + 2, 1) = qK_2.y;
	QQ(K + 2, 2) = qK_2.z;

	Eigen::MatrixXd P = Eigen::MatrixXd::Zero(K + 3, 3);
	P = N.partialPivLu().solve(QQ);

	onc.Create(3, false, 6, K + 5);
	for (int i = 0; i < solve_N.KnotCount(); i++)
	{
		onc.SetKnot(i, solve_N.Knot(i));
	}
	onc.SetCV(0, Q[0]);
	for (int i = 1; i < K + 4; i++)
	{
		onc.SetCV(i, ON_3dPoint(P(i - 1, 0), P(i - 1, 1), P(i - 1, 2)));
	}
	onc.SetCV(K + 4, Q[K]);
}

Spiral::Spiral(ON_3dPoint start, ON_3dPoint end, ON_3dVector alpha0, ON_3dVector alphaL, double kappa0, double kappaL, double tau0, double tauL, double L, int n)
{
	alpha0.Unitize();
	alphaL.Unitize();
	Eigen::MatrixXd A = Eigen::MatrixXd::Identity(4, 4);
	for (int k = 0; k < n; k++)
	{
		Eigen::MatrixXd B = Eigen::MatrixXd::Identity(4, 4);
		double s = L / n;
		double kk = (kappaL - kappa0) * k / n + kappa0;
		double rr = (tauL - tau0) * k / n + tau0;
		B(0, 1) = s;
		B(1, 2) = s * kk;
		B(2, 1) = -B(1, 2);
		B(2, 3) = rr * s;
		B(3, 2) = -B(2, 3);

		B(1, 1) = B(1, 1) / sqrt(1 + pow(s * kk, 2));
		B(1, 2) = B(1, 2) / sqrt(1 + pow(s * kk, 2));
		B(2, 1) = B(2, 1) / sqrt(1 + pow(s * kk, 2) + pow(s * rr, 2));
		B(2, 2) = B(2, 2) / sqrt(1 + pow(s * kk, 2) + pow(s * rr, 2));
		B(2, 3) = B(2, 3) / sqrt(1 + pow(s * kk, 2) + pow(s * rr, 2));
		B(3, 2) = B(3, 2) / sqrt(1 + pow(s * rr, 2));
		B(3, 3) = B(3, 3) / sqrt(1 + pow(s * rr, 2));

		A = B * A;
	}
	Eigen::MatrixXd A11 = A.block<2, 2>(0, 0);
	Eigen::MatrixXd A12 = A.block<2, 2>(0, 2);
	Eigen::Matrix<double, 2, 3> raL;
	raL(0, 0) = end.x;
	raL(0, 1) = end.y;
	raL(0, 2) = end.z;
	raL(1, 0) = alphaL.x;
	raL(1, 1) = alphaL.y;
	raL(1, 2) = alphaL.z;
	Eigen::Matrix<double, 2, 3> ra0;
	ra0(0, 0) = start.x;
	ra0(0, 1) = start.y;
	ra0(0, 2) = start.z;
	ra0(1, 0) = alpha0.x;
	ra0(1, 1) = alpha0.y;
	ra0(1, 2) = alpha0.z;
	Eigen::MatrixXd C = A12.inverse() * (raL - A11 * ra0);
	Eigen::Matrix<double, 4, 3> X;
	X.block(0, 0, 2, 3) = ra0;
	X.block(2, 0, 2, 3) = C;
	std::vector<ON_3dPoint> all_points;
	all_points.push_back(start);
	m_all_tan.push_back(alpha0);
	m_all_curvature.push_back(kappa0);
	m_beta0 = ON_3dVector(X(2, 0), X(2, 1), X(2, 2));
	m_beta0.Unitize();
	// ON_3dVector gamma0 = ON_3dVector(X(3, 0), X(3, 1), X(3, 2));
	// double tb = m_beta0.Length();
	// double tg = gamma0.Length();
	// double pro = ON_3dVector::DotProduct(m_beta0, gamma0);
	for (int k = 0; k < n; k++)
	{
		Eigen::MatrixXd B = Eigen::MatrixXd::Identity(4, 4);
		double s = L / n;
		double kk = (kappaL - kappa0) * k / n + kappa0;
		double rr = (tauL - tau0) * k / n + tau0;
		B(0, 1) = s;
		B(1, 2) = s * kk;
		B(2, 1) = -B(1, 2);
		B(2, 3) = s * rr;
		B(3, 2) = -B(2, 3);

		B(1, 1) = B(1, 1) / sqrt(1 + pow(s * kk, 2));
		B(1, 2) = B(1, 2) / sqrt(1 + pow(s * kk, 2));
		B(2, 1) = B(2, 1) / sqrt(1 + pow(s * kk, 2) + pow(s * rr, 2));
		B(2, 2) = B(2, 2) / sqrt(1 + pow(s * kk, 2) + pow(s * rr, 2));
		B(2, 3) = B(2, 3) / sqrt(1 + pow(s * kk, 2) + pow(s * rr, 2));
		B(3, 2) = B(3, 2) / sqrt(1 + pow(s * rr, 2));
		B(3, 3) = B(3, 3) / sqrt(1 + pow(s * rr, 2));

		X = B * X;

		for (int i = 1; i < 4; i++)
		{
			ON_3dVector a = ON_3dVector(X(i, 0), X(i, 1), X(i, 2));
			a.Unitize();
			X(i, 0) = a.x;
			X(i, 1) = a.y;
			X(i, 2) = a.z;
		}

		all_points.push_back(ON_3dPoint(X(0, 0), X(0, 1), X(0, 2)));
		ON_3dVector v = ON_3dVector(X(1, 0), X(1, 1), X(1, 2));
		v.Unitize();
		m_all_tan.push_back(v);
		m_all_curvature.push_back(((kappaL - kappa0) * (k + 1)) / n + kappa0);
	}
	m_all_points = all_points;
	m_Length = L;
	m_N = n;
	m_betaN = ON_3dVector(X(2, 0), X(2, 1), X(2, 2));
	// m_betaN.Unitize();
}

void Spiral::Add_spiral_to_model(ONX_Model *model, const wchar_t *name, ON_Color color)
{
	const int layer_index = model->AddLayer(name, color);
	ON_3dmObjectAttributes *attributes = new ON_3dmObjectAttributes();
	attributes->m_layer_index = layer_index;
	attributes->m_name = name;
	model->AddManagedModelGeometryComponent(&m_three_degree_spiral, attributes);
}

void Spiral::ComputeNurbs_Directly()
{
	int stride = m_N / 20;
	std::vector<ON_3dPoint> part_of_points;
	std::vector<double> knot;
	int index = 0;
	while (index < m_N)
	{
		part_of_points.push_back(m_all_points[index]);
		knot.push_back(m_Length / m_N * index);
		index += stride;
	}
	part_of_points.push_back(m_all_points[m_N]);
	knot.push_back(m_Length);
	unsigned int size = knot.size();
	ThreeDegreeBsplineInterplate_Tan(m_three_degree_spiral, part_of_points, knot, m_all_tan[0], m_all_tan[m_N]);
}

void Spiral::ComputePiecesBezier()
{
	std::vector<ON_NurbsCurve *> vbc;
	int stride = m_N / 20;
	double unit_arclength = m_Length / m_N * stride;
	// std::vector <ON_3dPoint> part_of_points;
	// std::vector <double> knot;
	int index = 0;
	std::vector<int> index_array;
	while (index < m_N)
	{
		index_array.push_back(index);
		index += stride;
	}
	index_array.push_back(m_N);
	auto it = index_array.cbegin();
	int j = 0;
	while (it != index_array.cend() - 1)
	{
		ON_BezierCurve *bc = new ON_BezierCurve();
		Rational_Bezier_G2(bc, m_all_points[*it], m_all_points[*(it + 1)],
						   m_all_tan[*it] * unit_arclength,
						   m_all_tan[*(it + 1)] * unit_arclength,
						   m_all_curvature[*it], m_all_curvature[*(it + 1)]);
		ON_NurbsCurve *onc = new ON_NurbsCurve();
		bc->GetNurbForm(*onc);
		onc->SetDomain(unit_arclength * j, unit_arclength * j + unit_arclength);
		j++;
		vbc.push_back(onc);
		delete bc;
		it++;
	}
	int num_piece = index_array.size() - 1;
	m_pieces_bezier.Create(3, true, 4, num_piece * 3 + 1);
	for (int i = 0; i < num_piece; i++)
	{
		m_pieces_bezier.SetKnot(i * 3, unit_arclength * i);
		m_pieces_bezier.SetKnot(i * 3 + 1, unit_arclength * i);
		m_pieces_bezier.SetKnot(i * 3 + 2, unit_arclength * i);
		ON_3dPoint p;
		vbc[i]->GetCV(0, p);
		m_pieces_bezier.SetCV(i * 3, p);
		vbc[i]->GetCV(1, p);
		double w = vbc[i]->Weight(1);
		m_pieces_bezier.SetCV(i * 3 + 1, p * w);
		m_pieces_bezier.SetWeight(i * 3 + 1, w);
		vbc[i]->GetCV(2, p);
		w = vbc[i]->Weight(2);
		m_pieces_bezier.SetCV(i * 3 + 2, p * w);
		m_pieces_bezier.SetWeight(i * 3 + 2, w);
		delete vbc[i];
		vbc[i] = nullptr;
	}
	m_pieces_bezier.SetKnot(num_piece * 3, m_Length);
	m_pieces_bezier.SetKnot(num_piece * 3 + 1, m_Length);
	m_pieces_bezier.SetKnot(num_piece * 3 + 2, m_Length);
	m_pieces_bezier.SetCV(num_piece * 3, m_all_points[m_N]);
}

void Spiral::Add_pieces_to_model(ONX_Model *model, const wchar_t *name, ON_Color color)
{
	const int layer_index = model->AddLayer(name, color);
	ON_3dmObjectAttributes *attributes = new ON_3dmObjectAttributes();
	attributes->m_layer_index = layer_index;
	attributes->m_name = name;
	model->AddManagedModelGeometryComponent(&m_pieces_bezier, attributes);
}

void Spiral::ComputeNurbs_Five_degrees()
{
	int stride = m_N / 20;
	std::vector<ON_3dPoint> part_of_points;
	std::vector<double> knot;
	int index = 0;
	while (index < m_N)
	{
		part_of_points.push_back(m_all_points[index]);
		knot.push_back(m_Length / m_N * index);
		index += stride;
	}
	part_of_points.push_back(m_all_points[m_N]);
	knot.push_back(m_Length);
	unsigned int size = knot.size();

	FiveDegreeBsplineInterplate_Tan_and_Cur(m_five_degree_spiral, part_of_points,
											knot, m_all_tan[0], m_all_tan[m_N],
											m_beta0 * m_all_curvature[0], m_betaN * m_all_curvature[m_N]);
}

void Spiral::Add_five_degree_spiral_to_model(ONX_Model *model, const wchar_t *name, ON_Color color)
{
	const int layer_index = model->AddLayer(name, color);
	ON_3dmObjectAttributes *attributes = new ON_3dmObjectAttributes();
	attributes->m_layer_index = layer_index;
	attributes->m_name = name;
	model->AddManagedModelGeometryComponent(&m_five_degree_spiral, attributes);
}

void Spiral::Get_five_degree_spline(ON_NurbsCurve &onc)
{
	onc = m_five_degree_spiral;
}

void Compute_A_and_AL(Eigen::MatrixXd &AA, Eigen::MatrixXd &AAL, double kappa0, double kappaL, double tau0, double tauL, double L, int n)
{
	std::vector<Eigen::MatrixXd> Avector;
	std::vector<Eigen::MatrixXd> ALvector;
	Eigen::MatrixXd A = Eigen::MatrixXd::Identity(4, 4);
	Eigen::MatrixXd AL = Eigen::MatrixXd::Zero(4, 4);
	for (int k = 0; k < n; k++)
	{
		Eigen::MatrixXd B = Eigen::MatrixXd::Identity(4, 4);
		double s = L / n;
		double kk = (kappaL - kappa0) * k / n + kappa0;
		double rr = (tauL - tau0) * k / n + tau0;
		B(0, 1) = s;
		B(1, 2) = s * kk;
		B(2, 1) = -B(1, 2);
		B(2, 3) = rr * s;
		B(3, 2) = -B(2, 3);

		B(1, 1) = B(1, 1) / sqrt(1 + pow(s * kk, 2));
		B(1, 2) = B(1, 2) / sqrt(1 + pow(s * kk, 2));
		B(2, 1) = B(2, 1) / sqrt(1 + pow(s * kk, 2) + pow(s * rr, 2));
		B(2, 2) = B(2, 2) / sqrt(1 + pow(s * kk, 2) + pow(s * rr, 2));
		B(2, 3) = B(2, 3) / sqrt(1 + pow(s * kk, 2) + pow(s * rr, 2));
		B(3, 2) = B(3, 2) / sqrt(1 + pow(s * rr, 2));
		B(3, 3) = B(3, 3) / sqrt(1 + pow(s * rr, 2));

		Avector.push_back(B);

		Eigen::MatrixXd C = Eigen::MatrixXd::Zero(4, 4);
		C(0, 1) = 1.0 / n;
		C(1, 1) = -(L * kk * kk) / (n * n) * pow(1 + s * kk * s * kk, -1.5);
		C(1, 2) = kk / n * pow(1 + s * kk * s * kk, -0.5) + s * kk * C(1, 1);
		C(3, 3) = -(L * rr * rr) / (n * n) * pow(1 + s * rr * s * rr, -1.5);
		C(3, 2) = -rr / n * pow(1 + s * rr * s * rr, -0.5) - s * rr * C(3, 3);
		C(2, 2) = -(L * (kk * kk + rr * rr)) / (n * n) * pow(1 + s * s * (kk * kk + rr * rr), -1.5);
		C(2, 1) = -kk / n * pow(1 + s * s * (kk * kk + rr * rr), -0.5) - s * kk * C(2, 2);
		C(2, 3) = rr / n * pow(1 + s * s * (kk * kk + rr * rr), -0.5) + s * rr * C(2, 2);
		ALvector.push_back(C);
	}
	for (int i = 0; i < n; i++)
	{
		A = Avector[i] * A;
		Eigen::MatrixXd B = Eigen::MatrixXd::Identity(4, 4);
		for (int j = 0; j < n; j++)
		{
			if (j == i)
			{
				B = ALvector[j] * B;
			}
			else
			{
				B = Avector[j] * B;
			}
		}
		AL = AL + B;
	}
	AA = A;
	AAL = AL;
}

void Compute_X0theta(ON_3dVector &beta, ON_3dVector &gamma, ON_3dVector &betatheta, ON_3dVector &gammatheta, ON_3dVector alpha0, ON_3dVector beta0, double theta)
{
	ON_Xform Rotation;
	Rotation.Rotation(theta, alpha0, ON_3dPoint(0, 0, 0));
	beta0.Transform(Rotation);
	beta = beta0;
	gamma = ON_3dVector::CrossProduct(alpha0, beta);
	gamma.Unitize();
}

void GD_Method(ON_3dPoint start, ON_3dPoint end, ON_3dVector alpha0, ON_3dVector alphaL, double kappa0, double kappaL, double tau0, double tauL, double L, int n,
			   double &returnL, double &returntheta)
{
	alpha0.Unitize();
	alphaL.Unitize();
	Eigen::MatrixXd A = Eigen::MatrixXd::Identity(4, 4);
	Eigen::MatrixXd AL = Eigen::MatrixXd::Zero(4, 4);
	Compute_A_and_AL(A, AL, kappa0, kappaL, tau0, tauL, L, n);
	ON_3dVector x(1, 0, 0);
	ON_3dVector y(0, 1, 0);
	ON_3dVector beta0;
	if (abs(ON_3dVector::DotProduct(alpha0, x)) > abs(ON_3dVector::DotProduct(alpha0, y)))
	{
		beta0 = y;
	}
	else
	{
		beta0 = x;
	}
	beta0 = ON_3dVector::CrossProduct(alpha0, beta0);
	beta0.Unitize();
}