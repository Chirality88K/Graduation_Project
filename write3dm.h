#ifndef WRITE3DM_H
#define WRITE3DM_H

#include "thirdparty/opennurbs/opennurbs.h"
#include <iostream>
#include <string>

void ChiralityWrite3dmModel(const ONX_Model *model, const std::string &filename);
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
void PrintPosAndTan(const ON_NurbsCurve &onc, const std::string &filename_without_extension = "Pos&Tan");
void ChiralityDebugInfo(const ON_NurbsCurve &onc, const std::string &filename_without_extension = "Bspline Debug");
#endif