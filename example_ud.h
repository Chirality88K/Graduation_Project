#if !defined(OPENNURBS_EXAMPLE_UD_INC_)
#define OPENNURBS_EXAMPLE_UD_INC_

class CExampleWriteUserData : public ON_UserData
{
  static int m__sn;
  ON_OBJECT_DECLARE(CExampleWriteUserData);

public:
  static ON_UUID Id();

  CExampleWriteUserData();
  virtual ~CExampleWriteUserData();

  CExampleWriteUserData( const char* s);
  CExampleWriteUserData(const CExampleWriteUserData& src);
  CExampleWriteUserData& operator=(const CExampleWriteUserData& src);

  void Dump( ON_TextLog& text_log ) const override;
  bool GetDescription( ON_wString& description ) override;
  bool Archive() const override;
  bool Write(ON_BinaryArchive& file) const override;
  bool Read(ON_BinaryArchive& file) override;

  ON_wString m_str;
  int m_sn;
};

#define INTERNAL_INITIALIZE_MODEL(model) Internal_SetExampleModelProperties(model,OPENNURBS__FUNCTION__,__FILE__)

static void Internal_SetExampleModelProperties(
    ONX_Model& model,
    const char* function_name,
    const char* source_file_name
)
{
    const bool bHaveFunctionName = (nullptr != function_name && 0 != function_name[0]);
    if (!bHaveFunctionName)
        function_name = "";

    const bool bHaveFileName = (nullptr != source_file_name && 0 != source_file_name[0]);
    if (!bHaveFileName)
        source_file_name = "";

    model.m_sStartSectionComments = "This was file created by OpenNURBS toolkit example code.";

    // set application information
    const ON_wString wide_function_name(function_name);
    const ON_wString wide_source_file_name(source_file_name);
    model.m_properties.m_Application.m_application_name
        = bHaveFunctionName
        ? ON_wString::FormatToString(L"OpenNURBS toolkit Example: %ls() function", static_cast<const wchar_t*>(wide_function_name))
        : ON_wString(L"OpenNURBS Examples");

    model.m_properties.m_Application.m_application_URL = L"http://www.opennurbs.org";
    model.m_properties.m_Application.m_application_details
        = bHaveFileName
        ? ON_wString::FormatToString(L"Opennurbs examples are in the file %ls.", static_cast<const wchar_t*>(wide_source_file_name))
        : ON_wString::FormatToString(L"Opennurbs examples are example_*.cpp files.");

    // some notes
    if (bHaveFunctionName && bHaveFileName)
    {
        model.m_properties.m_Notes.m_notes
            = ON_wString::FormatToString(
                L"This .3dm file was made with the OpenNURBS toolkit example function %s() defined in source code file %ls.",
                static_cast<const wchar_t*>(wide_function_name),
                static_cast<const wchar_t*>(wide_source_file_name));
        model.m_properties.m_Notes.m_bVisible = model.m_properties.m_Notes.m_notes.IsNotEmpty();
    }

    // set revision history information
    model.m_properties.m_RevisionHistory.NewRevision();
}

static bool Internal_WriteExampleModel(
    const ONX_Model& model,
    const wchar_t* filename,
    ON_TextLog& error_log
)
{
    int version = 0;

    // writes model to archive
    return model.Write(filename, version, &error_log);
}

#endif
