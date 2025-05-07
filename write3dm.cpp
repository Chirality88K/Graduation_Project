#include "write3dm.h"

void ChiralityWrite3dmModel(const ONX_Model *model, const wchar_t *filename) {
  ON_TextLog error_log;
  bool success = model->Write(filename, 0, &error_log);
  if (success) {
    std::cout << "Successfully wrote ";
    std::wcout << std::wstring(filename) << "\n";
  } else {
    std::cout << "Fail to wrtie ";
    std::wcout << std::wstring(filename) << "\n";
  }
}

ON_3dmObjectAttributes *Internal_CreateManagedAttributes(int layer_index,
                                                         const wchar_t *name) {
  ON_3dmObjectAttributes *attributes = new ON_3dmObjectAttributes();
  attributes->m_layer_index = layer_index;
  attributes->m_name = name;
  return attributes;
}

void Internal_SetExampleModelProperties(ONX_Model &model,
                                        const char *function_name,
                                        const char *source_file_name) {
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
  if (bHaveFunctionName && bHaveFileName) {
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
                                ON_TextLog &error_log) {
  int version = 0;

  // writes model to archive
  return model.Write(filename, version, &error_log);
}