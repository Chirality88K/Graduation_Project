#ifdef WIN32
#pragma comment(lib, "Shlwapi.lib")
#endif
#include "write3dm.h"
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <codecvt>
#include <sstream>
#include <iomanip>
#include "ChiralityLog.h"
#include "ChiralityMathTools.h"

std::string output_base_dir = "./Table/notsettime/";

void SetOutput(const std::string& name)
{
	output_base_dir = "./Table/" + name + ChiralityPrintNowTime() + "/";
	std::filesystem::create_directory(output_base_dir);
}

std::string GetOutputDir()
{
	return output_base_dir;
}

void ChiralityWrite3dmModel(const ONX_Model *model, const std::string &filename)
{
	ON_TextLog error_log;
	std::string path = output_base_dir + filename;
	wchar_t *const wc = new wchar_t[path.size() + 1];
	std::mbstowcs(wc, path.c_str(), path.size() + 1);
	bool success = model->Write(wc, 0, &error_log);
	std::string output = "OpenNURBS Archive File:\t" + filename + "----";
	if (success)
	{
		output += "Successfully written.";
		CHIRALITY_INFO(output);
	}
	else
	{
		output += "Fail to be written.";
		CHIRALITY_ERROR(output);
	}
	delete[] wc;
}

std::wstring StringToWString(const std::string &s)
{
	std::wstring_convert<std::codecvt_utf8<wchar_t>> conv;
	return conv.from_bytes(s);
}

std::string ChiralityPrintNowTime()
{
	time_t seconds;
	time(&seconds);
	struct tm *p_tm = new tm();
#ifdef WIN32
	localtime_s(p_tm, &seconds);
#else
	localtime_r(&seconds, p_tm);
#endif

	std::string time_string;
	time_string += std::to_string(1900 + p_tm->tm_year);
	if (p_tm->tm_mon + 1 < 10)
	{
		time_string += "0";
	}
	time_string += std::to_string(p_tm->tm_mon + 1);
	if (p_tm->tm_mday < 10)
	{
		time_string += "0";
	}
	time_string += std::to_string(p_tm->tm_mday);
	time_string += "_";
	if (p_tm->tm_hour < 10)
	{
		time_string += "0";
	}
	time_string += std::to_string(p_tm->tm_hour);
	if (p_tm->tm_min < 10)
	{
		time_string += "0";
	}
	time_string += std::to_string(p_tm->tm_min);
	time_string += "_";
	if (p_tm->tm_sec < 10)
	{
		time_string += "0";
	}
	time_string += std::to_string(p_tm->tm_sec);
	delete p_tm;
	return time_string;
}

ON_3dmObjectAttributes *Internal_CreateManagedAttributes(int layer_index,
														 const wchar_t *name)
{
	ON_3dmObjectAttributes *attributes = new ON_3dmObjectAttributes();
	attributes->m_layer_index = layer_index;
	attributes->m_name = name;
	return attributes;
}

void Internal_SetExampleModelProperties(ONX_Model &model,
										const char *function_name,
										const char *source_file_name)
{
	const bool bHaveFunctionName =
		(nullptr != function_name && 0 != function_name[0]);
	if (!bHaveFunctionName)
		function_name = "";

	const bool bHaveFileName =
		(nullptr != source_file_name && 0 != source_file_name[0]);
	if (!bHaveFileName)
		source_file_name = "";

	model.m_sStartSectionComments =
		"This was file created by Chirality.";

	// set application information
	const ON_wString wide_function_name(function_name);
	const ON_wString wide_source_file_name(source_file_name);
	model.m_properties.m_Application.m_application_name =
		bHaveFunctionName ? ON_wString::FormatToString(
								L"OpenNURBS toolkit Example: %ls() function",
								static_cast<const wchar_t *>(wide_function_name))
						  : ON_wString(L"OpenNURBS Examples");

	model.m_properties.m_Application.m_application_URL =
		L"http://www.opennurbs.org";
	model.m_properties.m_Application.m_application_details =
		bHaveFileName ? ON_wString::FormatToString(
							L"Opennurbs examples are in the file %ls.",
							static_cast<const wchar_t *>(wide_source_file_name))
					  : ON_wString::FormatToString(
							L"Opennurbs examples are example_*.cpp files.");

	// some notes
	if (bHaveFunctionName && bHaveFileName)
	{
		model.m_properties.m_Notes.m_notes = ON_wString::FormatToString(
			L"This .3dm file was made with the OpenNURBS toolkit example function "
			L"%s() defined in source code file %ls.",
			static_cast<const wchar_t *>(wide_function_name),
			static_cast<const wchar_t *>(wide_source_file_name));
		model.m_properties.m_Notes.m_bVisible =
			model.m_properties.m_Notes.m_notes.IsNotEmpty();
	}

	// set revision history information
	model.m_properties.m_RevisionHistory.NewRevision();
}

bool Internal_WriteExampleModel(const ONX_Model &model, const wchar_t *filename,
								ON_TextLog &error_log)
{
	int version = 0;

	// writes model to archive
	return model.Write(filename, version, &error_log);
}

void PrintCurvature(const ON_BezierCurve &onc, const std::string &filename_without_extension)
{
	double k0 = 0;
	double kn = 1;
	double t = 0;
	double kappa = 0;
	double lastkappa = 0;
	ON_3dVector v1;
	ON_3dVector v2;
	ON_3dVector v3;
	ON_3dPoint dump;
	std::string filename = output_base_dir + filename_without_extension + "-" + ChiralityPrintNowTime() + ".txt";
	if (onc.Dimension() == 2)
	{
		std::ofstream ofs(filename);
		for (int i = 0; i <= 1000; i++)
		{
			t = (kn - k0) / 1000 * i + k0;
			onc.Ev2Der(t, dump, v1, v2);
			kappa = v1.x * v2.y - v1.y * v2.x;
			kappa = kappa / pow(v1.Length(), 3);
			ofs << std::fixed << std::setprecision(6) << t << "\t" << kappa;
			if (i > 0)
			{
				ofs << "\t" << std::fixed << std::setprecision(6) << (kappa - lastkappa);
			}
			ofs << std::endl;
			lastkappa = kappa;
		}
		ofs.close();
		CHIRALITY_INFO(filename + " 2-dimension bezier curve curvatures written!");
		return;
	}
	if (onc.Dimension() == 3)
	{
		std::ofstream ofs(filename);
		for (int i = 0; i <= 1000; i++)
		{
			t = (kn - k0) / 1000 * i + k0;
			onc.Ev2Der(t, dump, v1, v2);
			kappa = onc.CurvatureAt(t).Length();
			ofs << std::fixed << std::setprecision(6) << t << "\t" << kappa;
			if (i > 0)
			{
				ofs << "\t" << std::fixed << std::setprecision(6) << (kappa - lastkappa);
			}
			ofs << std::endl;
			lastkappa = kappa;
		}
		ofs.close();
		CHIRALITY_INFO(filename + " 3-dimension bezier curve curvatures written!");
		return;
	}
	CHIRALITY_ERROR(filename + "Fail to write curvatures!");
}

void PrintCurvature(const ON_NurbsCurve &onc, const std::string &filename_without_extension)
{
	double k0 = 0;
	double kn = 1;
	onc.GetDomain(&k0, &kn);
	double t = 0;
	double kappa = 0;
	double lastkappa = 0;
	ON_3dVector v1;
	ON_3dVector v2;
	ON_3dVector v3;
	ON_3dPoint dump;
	std::string filename = output_base_dir + filename_without_extension + "-" + ChiralityPrintNowTime() + ".txt";
	if (onc.Dimension() == 2)
	{
		std::ofstream ofs(filename);
		for (int i = 0; i <= 1000; i++)
		{
			t = (kn - k0) / 1000 * i + k0;
			onc.Ev2Der(t, dump, v1, v2);
			kappa = v1.x * v2.y - v1.y * v2.x;
			kappa = kappa / pow(v1.Length(), 3);
			ofs << std::fixed << std::setprecision(6) << t << "\t" << kappa;
			if (i > 0)
			{
				ofs << "\t" << std::fixed << std::setprecision(6) << (kappa - lastkappa);
			}
			ofs << std::endl;
			lastkappa = kappa;
		}
		ofs.close();
		CHIRALITY_INFO(filename + " 2-dimension nurbs curve curvatures written!");
		return;
	}
	if (onc.Dimension() == 3)
	{
		std::ofstream ofs(filename);
		for (int i = 0; i <= 1000; i++)
		{
			t = (kn - k0) / 1000 * i + k0;
			onc.Ev2Der(t, dump, v1, v2);
			kappa = onc.CurvatureAt(t).Length();
			ofs << std::fixed << std::setprecision(6) << t << "\t" << kappa;
			if (i > 0)
			{
				ofs << "\t" << std::fixed << std::setprecision(6) << (kappa - lastkappa);
			}
			ofs << std::endl;
			lastkappa = kappa;
		}
		ofs.close();
		CHIRALITY_INFO(filename + " 3-dimension nurbs curve curvatures written!");
		return;
	}
	CHIRALITY_ERROR(filename + "Fail to write curvatures!");
}

void PrintDiscreteCurvature(const std::vector<ON_3dPoint> &vp, const std::string &filename_without_extension)
{
	std::string filename = output_base_dir + filename_without_extension + "-" + ChiralityPrintNowTime() + ".txt";
	std::ofstream ofs(filename);
	double length_total = 0.0;
	for (size_t i = 1; i < vp.size() - 1ull; ++i)
	{
		length_total += vp[i - 1ull].DistanceTo(vp[i]);
		double cur_dis = ChiralityMath::DiscreteCurvature(vp[i - 1ull], vp[i], vp[i + 1ull]);
		ofs << std::fixed << std::setprecision(6) << length_total << "\t" << cur_dis << std::endl;
	}
	ofs.close();
	CHIRALITY_INFO(filename + " discrete curvatures written!");
}

void PrintPosAndTan(const ON_NurbsCurve &onc, const std::string &filename_without_extension)
{
	std::string filename = output_base_dir + filename_without_extension + "-" + ChiralityPrintNowTime() + ".txt";
	std::ofstream ofs(filename);
	double k0 = 0;
	double kn = 1;
	onc.GetDomain(&k0, &kn);
	double t = 0;
	ON_3dPoint p;
	ON_3dVector v;
	for (int i = 0; i <= 1000; i++)
	{
		t = (kn - k0) / 1000 * i + k0;
		p = onc.PointAt(t);
		v = onc.TangentAt(t);
		ofs << std::fixed << std::setprecision(6) << t << "\t";
		ofs << "(" << p.x << "," << p.y << "," << p.z << ")" << "\t";
		ofs << "(" << v.x << "," << v.y << "," << v.z << ")";
		ofs << std::endl;
	}
	ofs.close();
	CHIRALITY_INFO(filename + " nurbs curve position and tangent written!");
}

void ChiralityDebugInfo(const ON_NurbsCurve &onc, const std::string &filename_without_extension)
{
	std::string filename = output_base_dir + filename_without_extension + "-" + ChiralityPrintNowTime() + ".txt";
	std::ofstream ofs(filename);
	ofs << "Type of curve:\t" << "B-Spline curve\n";
	ofs << "Order:\t" << onc.Order() << "\tDegree:\t" << onc.Degree() << "\tNumber of control points:\t" << onc.CVCount() << "\n";
	ofs << "Control points: \n";
	ON_3dPoint p;
	for (int i = 0; i < onc.CVCount(); ++i)
	{
		onc.GetCV(i, p);
		ofs << std::fixed << std::setprecision(6)
			<< "(" << p.x << "," << p.y << "," << p.z << ")\n";
	}
	ofs << "Sample points: \n";
	double k0 = 0;
	double kn = 1;
	onc.GetDomain(&k0, &kn);
	double t = 0;
	double kappa = 0;
	double lastkappa = 0;
	ON_3dVector v;
	ON_3dVector v1;
	ON_3dVector v2;
	if (onc.Dimension() == 2)
	{
		for (int i = 0; i <= 1000; i++)
		{
			t = (kn - k0) / 1000 * i + k0;
			v = onc.TangentAt(t);
			onc.Ev2Der(t, p, v1, v2);
			kappa = v1.x * v2.y - v1.y * v2.x;
			kappa = kappa / pow(v1.Length(), 3);
			ofs << std::fixed << std::setprecision(6) << t << "\t";
			ofs << "(" << p.x << "," << p.y << "," << p.z << ")" << "\t";
			ofs << "(" << v.x << "," << v.y << "," << v.z << ")\t";
			ofs << std::fixed << std::setprecision(6) << kappa;
			if (i > 0)
			{
				ofs << "\t" << std::fixed << std::setprecision(6) << (kappa - lastkappa);
			}
			ofs << std::endl;
			lastkappa = kappa;
		}
		ofs.close();
		CHIRALITY_INFO(filename + " 2-dimension nurbs curve debug written!");
		return;
	}
	if (onc.Dimension() == 3)
	{
		for (int i = 0; i <= 1000; i++)
		{
			t = (kn - k0) / 1000 * i + k0;
			v = onc.TangentAt(t);
			onc.Ev2Der(t, p, v1, v2);
			kappa = onc.CurvatureAt(t).Length();
			ofs << std::fixed << std::setprecision(6) << t << "\t";
			ofs << "(" << p.x << "," << p.y << "," << p.z << ")" << "\t";
			ofs << "(" << v.x << "," << v.y << "," << v.z << ")\t";
			ofs << std::fixed << std::setprecision(6) << kappa;
			if (i > 0)
			{
				ofs << "\t" << std::fixed << std::setprecision(6) << (kappa - lastkappa);
			}
			ofs << std::endl;
			lastkappa = kappa;
		}
		ofs.close();
		CHIRALITY_INFO(filename + " 3-dimension nurbs curve debug written!");
		return;
	}
	ofs.close();
	CHIRALITY_ERROR(filename + "Fail to write debug information!");
}

void ChiralityDebugInfo(const ON_NurbsSurface &ons, const std::string &filename_without_extension)
{
	std::string filename = output_base_dir + filename_without_extension + "-" + ChiralityPrintNowTime() + ".TXT";
	std::ofstream ofs(filename);
	ofs << "Type of surface:\t" << "Nurbs surface\n";
	ofs << "Order:\t" << ons.Order(0) << "\t" << ons.Order(1) << "\n";
	ofs << "Degree:\t" << ons.Degree(0) << "\t" << ons.Degree(1) << "\n";
	ofs << "Number of control points:\t" << ons.CVCount(0) << "\t" << ons.CVCount(1) << "\n";
	ofs << "Number of knots:\t" << ons.KnotCount(0) << "\t" << ons.KnotCount(1) << "\n";
	ofs << "Knots_0:\n";
	for (int i = 0; i < ons.KnotCount(0); ++i)
	{
		ofs << ons.Knot(0, i) << " ";
	}
	ofs << "\nKnots_1:\n";
	for (int i = 0; i < ons.KnotCount(1); ++i)
	{
		ofs << ons.Knot(1, i) << " ";
	}
	ofs << "\nControl points: \n";
	ON_3dPoint p;
	for (int i = 0; i < ons.CVCount(0); ++i)
	{
		for (int j = 0; j < ons.CVCount(1); ++j)
		{
			ons.GetCV(i, j, p);
			ofs << "(" << i << " , " << j << ")\t";
			ofs << std::fixed << std::setprecision(6)
				<< "(" << p.x << "," << p.y << "," << p.z << ")\n";
		}
	}
	ofs.close();
}

void ChiralityDebugforR(const ON_NurbsCurve &onc, const std::string &filename_without_extension)
{
	std::string filename = output_base_dir + filename_without_extension + "-" + ChiralityPrintNowTime() + ".txt";
	std::ofstream ofs(filename);
	ON_3dPoint p;
	double k0 = 0;
	double kn = 1;
	onc.GetDomain(&k0, &kn);
	double t = 0;
	double kappa = 0;
	ON_3dVector v;
	ON_3dVector v1;
	ON_3dVector v2;
	if (onc.Dimension() == 2)
	{
		for (int i = 0; i <= 1000; i++)
		{
			t = (kn - k0) / 1000 * i + k0;
			v = onc.TangentAt(t);
			onc.Ev2Der(t, p, v1, v2);
			kappa = v1.x * v2.y - v1.y * v2.x;
			kappa = kappa / pow(v1.Length(), 3);
			ofs << std::fixed << std::setprecision(6) << t << "\t";
			// ofs << "(" << p.x << "," << p.y << "," << p.z << ")" << "\t";
			// ofs << "(" << v.x << "," << v.y << "," << v.z << ")\t";
			ofs << std::fixed << std::setprecision(6) << kappa << "\t" << std::setprecision(6) << 0.0;
			ofs << std::endl;
		}
		ofs.close();
		CHIRALITY_INFO(filename + " 2-dimension nurbs curve debug written!");
		return;
	}
	if (onc.Dimension() == 3)
	{
		for (int i = 0; i <= 1000; i++)
		{
			t = (kn - k0) / 1000 * i + k0;
			v = onc.TangentAt(t);
			onc.Ev2Der(t, p, v1, v2);
			kappa = onc.CurvatureAt(t).Length();
			ofs << std::fixed << std::setprecision(6) << t << "\t";
			// ofs << "(" << p.x << "," << p.y << "," << p.z << ")" << "\t";
			// ofs << "(" << v.x << "," << v.y << "," << v.z << ")\t";
			ofs << std::fixed << std::setprecision(6) << kappa << "\t";
			ofs << ChiralityMath::Torsion(onc, t);
			ofs << std::endl;
		}
		ofs.close();
		CHIRALITY_INFO(filename + " 3-dimension nurbs curve debug written!");
		return;
	}
	ofs.close();
	CHIRALITY_ERROR(filename + "Fail to write debug information!");
}

void ChiralityDebugforR(const std::vector<ON_NurbsCurve> &onc_list, const std::string &filename_without_extension)
{
	if (onc_list.empty())
	{
		return;
	}

	std::string filename = output_base_dir + filename_without_extension + "-" + ChiralityPrintNowTime() + ".txt";
	std::ofstream ofs(filename);
	double t = 0;
	double kappa = 0;
	double LastT = 0;
	for (const ON_NurbsCurve &onc : onc_list)
	{
		double t0, t1;
		onc_list[0].GetDomain(&t0, &t1);
		for (int i = 0; i <= 200; i++)
		{
			t = (t1 - t0) / 200 * i + t0;
			kappa = onc.CurvatureAt(t).Length();
			ofs << std::fixed << std::setprecision(6) << t + LastT << "\t";
			ofs << std::fixed << std::setprecision(6) << kappa << "\t";
			ofs << ChiralityMath::Torsion(onc, t);
			ofs << std::endl;
		}
		LastT += t1 - t0;
	}
	ofs.close();
	CHIRALITY_INFO(filename + " 3-dimension nurbs curve debug written!");
}

void ChiralityAddNurbsCurve(ONX_Model *model, const ON_NurbsCurve &onc, const std::wstring &curve_name, int layer_index)
{
	ON_3dmObjectAttributes *att = new ON_3dmObjectAttributes();
	att->m_layer_index = layer_index;
	att->m_name = curve_name.c_str();
	ON_NurbsCurve *c = new ON_NurbsCurve(onc);
	model->AddManagedModelGeometryComponent(c, att);
}

void ChiralityAddNurbsSurface(ONX_Model *model, const ON_NurbsSurface &ons, const std::wstring &surface_name, int layer_index)
{
	ON_3dmObjectAttributes *att = new ON_3dmObjectAttributes();
	att->m_layer_index = layer_index;
	att->m_name = surface_name.c_str();
	ON_NurbsSurface *c = new ON_NurbsSurface(ons);
	model->AddManagedModelGeometryComponent(c, att);
}

void ChiralityAddPlane(ONX_Model *model, const ON_PlaneSurface &p, const std::wstring &plane_name, int layer_index)
{
	ON_3dmObjectAttributes *att = new ON_3dmObjectAttributes();
	att->m_layer_index = layer_index;
	att->m_name = plane_name.c_str();
	ON_PlaneSurface *ops = new ON_PlaneSurface(p);
	model->AddManagedModelGeometryComponent(ops, att);
}

void ChiralityAddQuadMesh(ONX_Model *model, const ParameterSurface &ps, int u_sample_num, int v_sample_num, const std::wstring &mesh_name, int layer_index)
{
	ON_3dmObjectAttributes *att = new ON_3dmObjectAttributes();
	att->m_layer_index = layer_index;
	att->m_name = mesh_name.c_str();
	ON_Mesh *mesh = new ON_Mesh(u_sample_num * v_sample_num, (u_sample_num + 1) * (v_sample_num + 1), true, false);
	std::vector<std::vector<ON_3dPoint>> points;
	std::vector<std::vector<ON_3dVector>> normals;
	ps.Discretize(points, normals, u_sample_num, v_sample_num);
	for (int i = 0; i <= u_sample_num; ++i)
	{
		for (int j = 0; j <= v_sample_num; ++j)
		{
			mesh->SetVertex(i * (v_sample_num + 1) + j, points[i][j]);
			mesh->SetVertexNormal(i * (v_sample_num + 1) + j, normals[i][j]);
		}
	}
	for (int i = 0; i < u_sample_num; ++i)
	{
		for (int j = 0; j < v_sample_num; ++j)
		{
			mesh->SetQuad(i * v_sample_num + j, i * (v_sample_num + 1) + j, i * (v_sample_num + 1) + j + 1, (i + 1) * (v_sample_num + 1) + j + 1, (i + 1) * (v_sample_num + 1) + j);
		}
	}
	model->AddManagedModelGeometryComponent(mesh, att);
}

void ChiralityAddLines(ONX_Model *model, const std::vector<ON_3dPoint> &vp, const std::wstring &lines_name, int layer_index)
{
	ON_3dPointArray parray;
	for (const ON_3dPoint &p : vp)
	{
		parray.Append(p);
	}
	ON_3dmObjectAttributes *att = new ON_3dmObjectAttributes();
	att->m_layer_index = layer_index;
	att->m_name = lines_name.c_str();
	ON_PolylineCurve *opc = new ON_PolylineCurve(ON_Polyline(parray));
	model->AddManagedModelGeometryComponent(opc, att);
}

void ChiralityAddCylindricalHelix(ONX_Model *model, double R, double ratio, double begin, double end, const std::wstring &name, int layer_index)
{
	auto CylindricalHelix = [R, ratio](double theta) -> ON_3dPoint
	{
		return ON_3dPoint(R * cos(theta), R * sin(theta), ratio * theta);
	};
	const int resolution = 200;
	ON_3dPointArray parray;
	for (int i = 0; i <= resolution; ++i)
	{
		double theta = begin * (1 - double(i) / double(resolution)) + end * double(i) / double(resolution);
		parray.Append(CylindricalHelix(theta));
	}
	ON_3dmObjectAttributes *att = new ON_3dmObjectAttributes();
	att->m_layer_index = layer_index;
	att->m_name = name.c_str();
	ON_PolylineCurve *opc = new ON_PolylineCurve(ON_Polyline(parray));
	model->AddManagedModelGeometryComponent(opc, att);
}

std::string doubleToScientificString(double value)
{
	std::ostringstream oss;
	oss << std::scientific << std::setprecision(15) << value;
	return oss.str();
}

void ChiralityDrawDNA(ONX_Model *model)
{
	const int index1 = model->AddLayer(L"test1", ON_Color::SaturatedRed);
	const int index2 = model->AddLayer(L"test2", ON_Color::SaturatedBlue);
	double R = 1.0;
	double ratio = 2.0;
	double begin = 0.0;
	double end = 4 * PI;
	auto CylindricalHelix = [R, ratio](double theta) -> ON_3dPoint
	{
		return ON_3dPoint(R * cos(theta), R * sin(theta), ratio * theta);
	};
	const int resolution = 200;
	ON_3dPointArray parray1;
	ON_3dPointArray parray2;
	for (int i = 0; i <= resolution; ++i)
	{
		double theta = begin * (1 - double(i) / double(resolution)) + end * double(i) / double(resolution);
		parray1.Append(CylindricalHelix(theta));
		parray2.Append(CylindricalHelix(theta + PI) - ON_3dVector::ZAxis * ratio * PI);
	}
	ON_3dmObjectAttributes *att1 = new ON_3dmObjectAttributes();
	att1->m_layer_index = index1;
	att1->m_name = L"DNA1";
	ON_PolylineCurve *opc1 = new ON_PolylineCurve(ON_Polyline(parray1));
	model->AddManagedModelGeometryComponent(opc1, att1);

	ON_3dmObjectAttributes *att2 = new ON_3dmObjectAttributes();
	att2->m_layer_index = index2;
	att2->m_name = L"DNA2";
	ON_PolylineCurve *opc2 = new ON_PolylineCurve(ON_Polyline(parray2));
	model->AddManagedModelGeometryComponent(opc2, att2);
}

void ChiralityPrintBezierForPython(const ON_BezierCurve& obc, const std::string& name)
{
	std::string filename = output_base_dir + name + ".bezier";
	std::ofstream ofs(filename);
	ON_3dPoint p;
	bool is_3dim = obc.Dimension() == 3;
	for (int i = 0; i < obc.CVCount(); ++i)
	{
		obc.GetCV(i, p);
		ofs << p.x << "\t" << p.y;
		if (is_3dim)
		{
			ofs << "\t" << p.z << "\n";
		}
		else
		{
			ofs << "\n";
		}
	}
	ofs.close();
}

void ChiralityPrintCubicIntKnotBSplineForPython(const ON_NurbsCurve& onc, const std::string& name)
{
	std::string filename = output_base_dir + name + ".bspline";
	std::ofstream ofs(filename);
	ON_3dPoint p;
	bool is_3dim = onc.Dimension() == 3;
	for (int i = 0; i < onc.CVCount(); ++i)
	{
		ofs << "CV\t";
		onc.GetCV(i, p);
		ofs << p.x << "\t" << p.y;
		if (is_3dim)
		{
			ofs << "\t" << p.z << "\n";
		}
		else
		{
			ofs << "\n";
		}
	}
	const int N = 1000;
	double u0, u1;
	onc.GetDomain(&u0, &u1);
	for (int i = 0; i < N; ++i)
	{
		double t = (1 - double(i) / (N - 1)) * u0 + double(i) / (N - 1) * u1;
		ofs << "P\t" << t << "\t";
		ON_3dVector dp, ddp;
		onc.Ev2Der(t, p, dp, ddp);
		ofs << p.x << "\t" << p.y;
		if (is_3dim)
		{
			double cur = ON_3dVector::CrossProduct(dp, ddp).Length() / pow(dp.Length(), 3);
			ofs << "\t" << p.z << "\t" << cur << "\n";
		}
		else
		{
			double cur = (dp.x * ddp.y - dp.y * ddp.x) / pow(dp.Length(), 3);
			ofs << "\t" << cur << "\n";
		}
	}
	ofs.close();
}
