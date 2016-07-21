#ifndef SUBCOMMANDS_H
#define SUBCOMMANDS_H
#include <fstream>

// Entry points to tools
void run_survey (struct argp_state* state);
void run_radial (struct argp_state* state);
void run_find (struct argp_state* state);

bool file_exists (const char *filename) {
  std::ifstream ifile(filename);
  return (bool)ifile;
}

#endif
