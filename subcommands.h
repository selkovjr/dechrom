#ifndef SUBCOMMANDS_H
#define SUBCOMMANDS_H
#include <fstream>

// Entry points to tools
void run_survey (struct argp_state* state);
void run_find (struct argp_state* state);
void run_plot (struct argp_state* state);
void run_radial (struct argp_state* state);
void run_view (struct argp_state* state);

#endif
