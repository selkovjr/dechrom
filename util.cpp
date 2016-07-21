#include <fstream>

bool file_exists (const char *filename) {
  std::ifstream ifile(filename);
  return (bool)ifile;
}


