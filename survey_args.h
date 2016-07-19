#include <argp.h>
#include "termcolor.h"
#include <iostream>
#include <iomanip>
#include <fstream>

using namespace std;
using namespace termcolor;

bool file_exists (const char *filename) {
  std::ifstream ifile(filename);
  return (bool)ifile;
}

// ------------------------
// ## survey command parser
//
struct arg_survey {
  double amin, amax;
  double bmin, bmax;
  double cmin, cmax;
  double dmin, dmax;
  unsigned nodes;
  char* ref_file;
  char* target_file;
};

static char args_doc_survey[] = "reference.tiff target.tiff";

static char doc_survey[] =
"\n"
"Calculate the TCA metric over a regular grid\n"
"surrounding the chosen set of parameters\n"
"\n"
"Inputs:\n"
"  A grayscale histogram-equalized TIFF reference (green channel)\n"
"  A grayscale histogram-equalized target TIFF (blue or red)\n"
"\n"
"Output:\n"
"  A table of parameter values and resulting TCA metric\n"
"  piped to stdout\n"
"\n"
"\v"
"For each set of parameters on a grid, the radial distortion\n"
"is applied to the target image and the TCA metric is calculatedi\n"
"as the average of absolute per-pixel differences between the\n"
"distorted  target image and the refernce image.\n"
"\n"
;

static error_t parse_survey_command(int key, char* arg, struct argp_state* state) {
  struct arg_survey* arguments = (struct arg_survey*)state->input;
  char **nonopt;

  assert( arguments );

  switch(key) {
    case 'n':
      if (sscanf(arg, "%u", &(arguments->nodes)) != 1) {
        cerr << on_red << "Expecting an integer number of grid nodes instead of " << bold << arg << reset << endl;
        exit(EXIT_FAILURE);
      }
      break;

    case 'a':
      if (sscanf(arg, "%lf:%lf", &(arguments->amin), &(arguments->amax)) != 2) {
        if (sscanf(arg, "%lf", &(arguments->amin)) == 1) {
          arguments->amax = arguments->amin;
        }
        else {
          cerr << on_red << "Expecting a floating-point number or a range (from:to) instead of " << bold << arg << reset << endl;
          exit(EXIT_FAILURE);
        }
      }
      break;

    case 'b':
      if (sscanf(arg, "%lf:%lf", &(arguments->bmin), &(arguments->bmax)) != 2) {
        if (sscanf(arg, "%lf", &(arguments->bmin)) == 1) {
          arguments->bmax = arguments->bmin;
        }
        else {
          cerr << on_red << "Expecting a floating-point number or a range (from:to) instead of " << bold << arg << reset << endl;
          exit(EXIT_FAILURE);
        }
      }
      break;

    case 'c':
      if (sscanf(arg, "%lf:%lf", &(arguments->cmin), &(arguments->cmax)) != 2) {
        if (sscanf(arg, "%lf", &(arguments->cmin)) == 1) {
          arguments->cmax = arguments->cmin;
        }
        else {
          cerr << on_red << "Expecting a floating-point number or a range (from:to) instead of " << bold << arg << reset << endl;
          exit(EXIT_FAILURE);
        }
      }
      break;

    case 'd':
      if (sscanf(arg, "%lf:%lf", &(arguments->dmin), &(arguments->dmax)) != 2) {
        if (sscanf(arg, "%lf", &(arguments->dmin)) == 1) {
          arguments->dmax = arguments->dmin;
        }
        else {
          cerr << on_red << "Expecting a floating-point number or a range (from:to) instead of " << bold << arg << reset << endl;
          exit(EXIT_FAILURE);
        }
      }
      break;

    case ARGP_KEY_NO_ARGS:
      argp_usage (state);

    case ARGP_KEY_ARG: // non-option argument
      nonopt = &state->argv[state->next];
      state->next = state->argc; // we're done
      arguments->ref_file = arg;
      arguments->target_file = nonopt[0];

      if (not file_exists(arguments->ref_file)) {
        cerr << on_red << "File does not exist: " << bold << arguments->ref_file << reset << endl;
        exit(EXIT_FAILURE);
      }

      if (not file_exists(arguments->target_file)) {
        cerr << on_red << "File does not exist: " << bold << arguments->target_file << reset << endl;
        exit(EXIT_FAILURE);
      }

      break;

    case ARGP_KEY_END:
      if (state->arg_num < 2) {
        argp_error(state, "Not enough arguments");
      }
      if (state->arg_num > 2) {
        argp_error(state, "Extra arguments");
      }
      break;

    default:
      return ARGP_ERR_UNKNOWN;
  }

  return 0;
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmissing-field-initializers"
static struct argp_option options_survey[] = {
  {"range-a",   'a', "RANGE",   0, "the range of a to survey, or a fixed point" },
  {"range-b",   'b', "RANGE",   0, "the range of b to survey, or a fixed point" },
  {"range-c",   'c', "RANGE",   0, "the range of c to survey, or a fixed point" },
  {"range-d",   'd', "RANGE",   0, "the range of d to survey, or a fixed point" },
  {"nodes",     'n', "NUMBER",  0, "the number of grid nodes in eeach parameter range" },
  { 0 }
};

static struct argp argp_survey = {
  options_survey,
  parse_survey_command,
  args_doc_survey,
  doc_survey
};
#pragma GCC diagnostic pop

#define PARSE_ARGS_SURVEY \
  struct arg_survey args; \
  int    argc = state->argc - state->next + 1; \
  char** argv = &state->argv[state->next - 1]; \
  char*  argv0 =  argv[0]; \
  argv[0] = (char *)malloc(strlen((char *)(state->name)) + strlen(" survey") + 1); \
  if (!argv[0]) argp_failure(state, 1, ENOMEM, 0); \
  sprintf(argv[0], "%s survey", state->name); \
  args.amin = 0.995; \
  args.amax = 1.005; \
  args.bmin = -0.005; \
  args.bmax = +0.005; \
  args.cmin = 0.0; \
  args.cmax = 0.0; \
  args.dmin = 0.0; \
  args.dmax = 0.0; \
  args.nodes = 5; \
  argp_parse(&argp_survey, argc, argv, ARGP_IN_ORDER, &argc, &args); \
  free(argv[0]); \
  argv[0] = argv0; \
  state->next += argc - 1; \

