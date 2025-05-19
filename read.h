#ifndef READ_H
#define READ_H
#include "opennurbs.h"
#include "example_ud.h"

static bool Dump3dmFileHelper(
	const wchar_t *sFileName, // full name of file
	ON_TextLog &dump)
{
	dump.Print("====== FILENAME: %ls\n", sFileName);
	ON_Workspace ws;
	FILE *fp = ws.OpenFile(sFileName, L"rb"); // file automatically closed by ~ON_Workspace()
	if (!fp)
	{
		dump.Print("**ERROR** Unable to open file.\n");
		return false;
	}

	ON_BinaryFile file(ON::archive_mode::read3dm, fp);

	int version = 0;
	ON_String comment_block;
	bool rc = file.Read3dmStartSection(&version, comment_block);
	if (!rc)
	{
		dump.Print("**ERROR** Read3dmStartSection() failed\n");
		return false;
	}
	dump.Print("====== VERSION: %d\n", version);
	dump.Print("====== COMMENT BLOCK:\n", version);
	dump.PushIndent();
	dump.Print(comment_block);
	dump.PopIndent();
	dump.Print("====== CHUNKS:\n", version);
	unsigned int typecode;
	while (!file.AtEnd())
	{
		typecode = file.Dump3dmChunk(dump, 0);
		if (!typecode)
			break;
		if (typecode == TCODE_ENDOFFILE)
			break;
	}
	dump.Print("====== FINISHED: %ls\n", sFileName);

	return true;
}

/*
Returns:
  True if .3dm file was successfully read into an ONX_Model class.
*/
static bool ReadFileHelper(
	const wchar_t *sFileName,
	bool bVerboseTextDump,
	bool bChunkDump,
	ON_TextLog &dump)
{
	if (bChunkDump)
	{
		return Dump3dmFileHelper(sFileName, dump);
	}

	ONX_Model model;

	dump.Print("\nOpenNURBS Archive File:  %ls\n", sFileName);

	// open file containing opennurbs archive
	FILE *archive_fp = ON::OpenFile(sFileName, L"rb");
	if (!archive_fp)
	{
		dump.Print("  Unable to open file.\n");
		return false;
	}

	dump.PushIndent();

	// create achive object from file pointer
	ON_BinaryFile archive(ON::archive_mode::read3dm, archive_fp);

	// read the contents of the file into "model"
	bool rc = model.Read(archive, &dump);

	// close the file
	ON::CloseFile(archive_fp);

	// print diagnostic
	if (rc)
		dump.Print("Successfully read.\n");
	else
		dump.Print("Errors during reading.\n");

	// create a text dump of the model
	if (bVerboseTextDump)
	{
		dump.PushIndent();
		model.Dump(dump);
		dump.PopIndent();
	}

	dump.PopIndent();

	return rc;
}

/*
Returns:
  Number of files read.
*/
static int ReadDirectoryHelper(
	int directory_depth,
	int maximum_directory_depth,
	const wchar_t *directory_name,
	const wchar_t *file_name_filter,
	bool bVerboseTextDump,
	bool bChunkDump,
	ON_TextLog &dump)
{
	int file_count = 0;
	if (directory_depth <= maximum_directory_depth)
	{
		if (0 == file_name_filter || 0 == file_name_filter[0])
			file_name_filter = L"*.3dm";

		// read files in this directory
		ON_FileIterator file_it;
		bool bFoundDirectory = false;
		for (bool bHaveFileSystemItem = (file_it.Initialize(directory_name, file_name_filter) && file_it.FirstItem());
			 bHaveFileSystemItem;
			 bHaveFileSystemItem = file_it.NextItem())
		{
			if (file_it.CurrentItemIsDirectory())
			{
				bFoundDirectory = true;
				continue;
			}

			if (false == file_it.CurrentItemIsFile())
				continue;

			if (file_it.CurrentItemIsHidden())
				continue;

			ON_wString full_path(file_it.CurrentItemFullPathName());
			if (full_path.IsEmpty())
				continue;

			if (!ON::IsOpenNURBSFile(full_path))
				continue;

			if (ReadFileHelper(full_path, bVerboseTextDump, bChunkDump, dump))
				file_count++;
		}

		// read files in subdirectories
		if (bFoundDirectory && directory_depth < maximum_directory_depth)
		{
			ON_FileIterator dir_it;
			for (bool bHaveFileSystemItem = (dir_it.Initialize(directory_name, nullptr) && dir_it.FirstItem());
				 bHaveFileSystemItem;
				 bHaveFileSystemItem = dir_it.NextItem())
			{
				if (false == dir_it.CurrentItemIsDirectory())
					continue;

				if (dir_it.CurrentItemIsHidden())
					continue;

				ON_wString full_path(dir_it.CurrentItemFullPathName());
				if (full_path.IsEmpty())
					continue;

				file_count += ReadDirectoryHelper(
					directory_depth + 1,
					maximum_directory_depth,
					full_path,
					file_name_filter,
					bVerboseTextDump,
					bChunkDump,
					dump);
			}
		}
	}

	return file_count;
}

static void print_help(const char *example_read_exe_name)
{
	if (0 == example_read_exe_name || 0 == example_read_exe_name[0])
		example_read_exe_name = "example_read";

	printf("\n");
	printf("SYNOPSIS:\n");
	printf("  %s [-out:outputfilename.txt] [-c] [-r] <file or directory names>\n", example_read_exe_name);
	printf("\n");
	printf("DESCRIPTION:\n");
	printf("  If a file is listed, it is read as an opennurbs model file.\n");
	printf("  If a directory is listed, all .3dm files in that directory\n");
	printf("  are read as opennurbs model files.\n");
	printf("\n");
	printf("  Available options:\n");
	printf("    -out:outputfilename.txt\n");
	printf("      The output is written to the named file.\n");
	printf("    -chunkdump\n");
	printf("      Does a chunk dump instead of reading the file's contents.\n");
	printf("    -recursive\n");
	printf("      Recursivly reads files in subdirectories.\n");
	printf("\n");
	printf("EXAMPLE:\n");
	printf("  %s -out:list.txt -resursive .../example_files\n", example_read_exe_name);
	printf("  with read all the opennurbs .3dm files in the\n");
	printf("  example_files/ directory and subdirectories.\n");
}

#if defined(ON_COMPILER_MSC)

// When you run a C program, you can use either of the two wildcards
// the question mark (?) and the asterisk (*) to specify filename
// and path arguments on the command line.
// By default, wildcards are not expanded in command-line arguments.
// You can replace the normal argument vector argv loading routine with
// a version that does expand wildcards by linking with the setargv.obj or wsetargv.obj file.
// If your program uses a main function, link with setargv.obj.
// If your program uses a wmain function, link with wsetargv.obj.
// Both of these have equivalent behavior.
// To link with setargv.obj or wsetargv.obj, use the /link option.
//
// For example:
// cl example.c /link setargv.obj
// The wildcards are expanded in the same manner as operating system commands.
// (See your operating system user's guide if you are unfamiliar with wildcards.)

// example_read.vcxproj linkin options include setargv.obj.

#endif
#endif