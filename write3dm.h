#ifndef WRITE3DM_H
#define WRITE3DM_H

#include "thirdparty/opennurbs/opennurbs.h"
#include <iostream>
#include <string>

void ChiralityWrite3dmModel(const ONX_Model *model, const std::string &filename);
std::wstring StringToWString(const std::string &s);
std::string ChiralityPrintNowTime();
ON_3dmObjectAttributes *Internal_CreateManagedAttributes(int layer_index, const wchar_t *name);
#define INTERNAL_INITIALIZE_MODEL(model) Internal_SetExampleModelProperties(model, OPENNURBS__FUNCTION__, __FILE__)
void Internal_SetExampleModelProperties(
	ONX_Model &model,
	const char *function_name,
	const char *source_file_name);

bool Internal_WriteExampleModel(
	const ONX_Model &model,
	const wchar_t *filename,
	ON_TextLog &error_log);

void PrintCurvature(const ON_BezierCurve &onc, const std::string &filename_without_extension = "BezierCurveCurvature");
void PrintCurvature(const ON_NurbsCurve &onc, const std::string &filename_without_extension = "NurbsCurveCurvature");
void PrintDiscreteCurvature(const std::vector<ON_3dPoint> &vp, const std::string &filename_without_extension = "DiscreteCurvature");
void PrintPosAndTan(const ON_NurbsCurve &onc, const std::string &filename_without_extension = "Pos&Tan");
void ChiralityDebugInfo(const ON_NurbsCurve &onc, const std::string &filename_without_extension = "Bspline Debug");
void ChiralityDebugforR(const ON_NurbsCurve &onc, const std::string &filename_without_extension = "Debug_for_R");
void ChiralityDebugforR(const std::vector<ON_NurbsCurve> &onc, const std::string &filename_without_extension = "Debug_for_R");
// 该函数内会自动使用new为Curve分配内存
void ChiralityAddNurbsCurve(ONX_Model *model, const ON_NurbsCurve &onc, const std::wstring &curve_name, int layer_index);
// 该函数内会自动使用new为Surface分配内存
void ChiralityAddNurbsSurface(ONX_Model *model, const ON_NurbsSurface &ons, const std::wstring &surface_name, int layer_index);
// 该函数内会自动使用new为Plane分配内存
void ChiralityAddPlane(ONX_Model *model, const ON_PlaneSurface &p, const std::wstring &plane_name, int layer_index);
void ChiralityAddLines(ONX_Model *model, const std::vector<ON_3dPoint> &vp, const std::wstring &lines_name, int layer_index);
void ChiralityAddCylindricalHelix(ONX_Model* model, double R, double ratio, double begin, double end, const std::wstring& name, int layer_index);
std::string doubleToScientificString(double value);
void ChiralityDrawDNA(ONX_Model* model);
#endif