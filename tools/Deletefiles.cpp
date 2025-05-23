#include <filesystem>
#include <iostream>

int main() {
  std::filesystem::path work_path = std::filesystem::current_path();
  for (auto &entry : std::filesystem::directory_iterator(work_path)) {
    if (entry.is_regular_file() &&
        (entry.path().extension().string() == ".3dm" ||
         entry.path().extension().string() == ".txt") &&
        entry.path().filename().string() != "CMakeCache.txt") {
      std::filesystem::remove(entry.path());
      std::cout << "Delete " << entry.path().filename().string() << "!!!!!!!\n";
    }
  }
  std::cout << "Work space clean!"
            << "\n";
}