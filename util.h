#ifndef UTIL_H
#define UTIL_H

#include <string>
#include <iostream>
#include <fstream>
#include "termcolor.h"

using namespace std;
using namespace termcolor;

#define RED 0
#define GREEN 1
#define BLUE 2
#define UNKNOWN 3

bool file_exists (const char *filename);
void color_switch(ostream& s, int color);

#endif


