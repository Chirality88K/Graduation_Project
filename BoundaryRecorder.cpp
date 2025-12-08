#include "BoundaryRecorder.h"
#include "write3dm.h"
#include <fstream>
#include <sstream>
#include <filesystem>

static const std::string boundary_dir_ = "./Table/BoundaryRecord/";

void BoundaryRecorder::Write() const
{
	if (boundary_.empty())
	{
		return;
	}
	const std::string filename = boundary_dir_ + ChiralityPrintNowTime() + ".txt";
	if (!std::filesystem::exists(boundary_dir_))
	{
		std::filesystem::create_directory(boundary_dir_);
	}
	std::ofstream ofs(filename);
	for (const Boundary& b : boundary_)
	{
		ofs << b.ps_.x << " " << b.ps_.y << " " << b.ps_.z << " ";
		ofs << b.pe_.x << " " << b.pe_.y << " " << b.pe_.z << " ";
		ofs << b.vs_.x << " " << b.vs_.y << " " << b.vs_.z << " ";
		ofs << b.ve_.x << " " << b.ve_.y << " " << b.ve_.z << "\n";
	}
	ofs.close();
}

void BoundaryRecorder::Read()
{
	if (!std::filesystem::exists(boundary_dir_) || !std::filesystem::is_directory(boundary_dir_)) 
	{
		std::cerr << "Cannot find: " << boundary_dir_ << std::endl;
		return;
	}

	for (const auto& entry : std::filesystem::directory_iterator(boundary_dir_)) {
		if (std::filesystem::is_regular_file(entry.status())) {
			std::ifstream infile(entry.path());
			if (!infile) {
				std::cerr << "Fail to open: " << entry.path() << std::endl;
				return;
			}
			std::string line;
			while (std::getline(infile, line)) {
				std::istringstream iss(line);
				std::vector<double> row;
				double val;
				int count = 0;
				while (iss >> val) {
					row.push_back(val);
					count++;
				}
				if (count != 12) {
					continue;
				}
				boundary_.push_back(Boundary(ON_3dPoint(row[0], row[1], row[2]),
					ON_3dPoint(row[3], row[4], row[5]),
					ON_3dVector(row[6], row[7], row[8]),
					ON_3dVector(row[9], row[10], row[11])));
			}

		}
	}


	
}
