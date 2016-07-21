#include "util.h"

bool file_exists (const char *filename) {
  std::ifstream ifile(filename);
  return (bool)ifile;
}

string to_string(double v) {
  stringstream ss;
  ss << v;
  return ss.str();
}

void color_switch(ostream& s, int color) {
  switch (color) {
    case RED:
      s << lightred;
      break;
    case GREEN:
      s << lightgreen;
      break;
    case BLUE:
      s << lightblue;
      break;
    default:
      s << yellow;
  }
}




