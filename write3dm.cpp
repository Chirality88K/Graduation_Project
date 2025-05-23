#include "write3dm.h"
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <time.h>

void ChiralityWrite3dmModel(const ONX_Model *model, const std::string &filename)
{
	ON_TextLog error_log;
	wchar_t *const wc = new wchar_t[filename.size() + 1];
	std::mbstowcs(wc,filename.c_str(),filename.size() +1);
	bool success = model->Write(wc, 0, &error_log);
	if (success)
	{
		std::cout << "Successfully wrote ";
		std::cout << (filename) << "\n";
	}
	else
	{
		std::cout << "Fail to wrtie ";
		std::cout << (filename) << "\n";
	}
	delete[] wc;
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
		"This was file created by OpenNURBS toolkit example code.";

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
	std::string filename = filename_without_extension + "-" + ChiralityPrintNowTime() + ".txt";
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
		std::cout << filename + " 2-dimension bezier curve curvatures written!" << std::endl;
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
		std::cout << filename + " 3-dimension bezier curve curvatures written!" << std::endl;
		return;
	}
	std::cout << "Fail to write curvatures: " + filename << std::endl;
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
	std::string filename = filename_without_extension + "-" + ChiralityPrintNowTime() + ".txt";
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
		std::cout << filename + " 2-dimension nurbs curve curvatures written!" << std::endl;
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
		std::cout << filename + " 3-dimension nurbs curve curvatures written!" << std::endl;
		return;
	}
	std::cout << "Fail to write curvatures: " + filename << std::endl;
}

void PrintPosAndTan(const ON_NurbsCurve &onc, const std::string &filename_without_extension)
{
	std::string filename = filename_without_extension + "-" + ChiralityPrintNowTime() + ".txt";
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
	std::cout << filename + " nurbs curve position and tangent written!" << std::endl;
}

void ChiralityDebugInfo(const ON_NurbsCurve &onc, const std::string &filename_without_extension)
{
	std::string filename = filename_without_extension + "-" + ChiralityPrintNowTime() + ".txt";
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
		std::cout << filename + " 2-dimension nurbs curve debug written!" << std::endl;
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
		std::cout << filename + " 3-dimension nurbs curve debug written!" << std::endl;
		return;
	}
	ofs.close();
}

void ChiralityAddNurbsCurve(ONX_Model *model, const ON_NurbsCurve &onc, const std::wstring &curve_name, int layer_index)
{
	ON_3dmObjectAttributes *att = new ON_3dmObjectAttributes();
	att->m_layer_index = layer_index;
	att->m_name = curve_name.c_str();
	ON_NurbsCurve *c = new ON_NurbsCurve(onc);
	model->AddManagedModelGeometryComponent(c, att);
}
