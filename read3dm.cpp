#include "read3dm.h"
#include <iostream>
#include <fstream>
#include "ChiralityLog.h"

void ChiralityRead3dmModel(const std::string &filename, ONX_Model *model)
{
	ON_TextLog dump;
	wchar_t *const wc = new wchar_t[filename.size() + 1];
	std::mbstowcs(wc, filename.c_str(), filename.size() + 1);
	FILE *archive_fp = ON::OpenFile(wc, L"rb");
	delete[] wc;
	std::string output = "OpenNURBS Archive File:\t" + filename + "----";
	if (!archive_fp)
	{
		output += "Not Found.";
		CHIRALITY_ERROR(output);
		return;
	}
	ON_BinaryFile archive(ON::archive_mode::read3dm, archive_fp);
	bool rc = model->Read(archive, &dump);
	ON::CloseFile(archive_fp);
	if (rc)
	{
		output += "Successfully read.";
		CHIRALITY_INFO(output);
	}
	else
	{
		output += "Fail to be read.";
		CHIRALITY_ERROR(output);
	}
}

void ChiralityModelDebugInfo(const std::string &filename)
{
	std::wofstream debug_info(filename + "-DEBUG.txt");
	ON_wString on_w_str;
	ON_TextLog dump(on_w_str);
	wchar_t *const wc = new wchar_t[filename.size() + 1];
	std::mbstowcs(wc, filename.c_str(), filename.size() + 1);
	FILE *archive_fp = ON::OpenFile(wc, L"rb");
	std::string output = "Chirality's Debug File:\t" + filename + "-DEBUG.txt----";
	if (!archive_fp)
	{
		dump.Print("  Fail to open file.\n");
		output += "Not Found.";
		CHIRALITY_ERROR(output);
		return;
	}
	dump.Print("OpenNURBS Archive File:  %ls\n", wc);
	delete[] wc;
	dump.PushIndent();
	ON_BinaryFile archive(ON::archive_mode::read3dm, archive_fp);
	ONX_Model model;
	bool rc = model.Read(archive, &dump);
	ON::CloseFile(archive_fp);
	if (rc)
	{
		dump.Print("Successfully read.\n");
		output += "Successfully read.";
		CHIRALITY_INFO(output);
	}
	else
	{
		dump.Print("Errors during reading.\n");
		output += "Fail to be read.";
		CHIRALITY_ERROR(output);
		return;
	}
	dump.PopIndent();
	model.Dump(dump);
	dump.PopIndent();
	dump.PopIndent();

	for (int i = 0; i < on_w_str.Length(); ++i)
	{
		debug_info << on_w_str.GetAt(i);
	}
}
