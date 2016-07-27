#include <argp.h>
#include <assert.h>
#include <string.h>
#include "termcolor.h"
#include <iostream>
#include <iomanip>
#include <fstream>

#include "util.h"

using namespace std;
using namespace termcolor;

// ------------------------
// ## plot command parser
//
struct arg_plot {
  int argx, argy;
  double z[4], f;
  char const* device;
  char const* surface_file;
  char const* trace_file;
  char const* plot_file;
  bool contour;
  bool vertices;
  bool minimizer;
  bool minimizer_short;
};

static char args_doc_plot[] = "reference.tiff work.tiff";

static char doc_plot[] =
"\n"
"Plot the TCA survey surface with optional solution points and\n"
"optimization traces\n"
"\n"
"Inputs:\n"
"  TCA surface data (-s), the output from `dechrom survey`)\n"
"  (optional) optimization trace in JSON format (-t), from `dechrom find`)\n"
"  (optional) objective-function minimizer (-z), from `dechrom find`)\n"
"\n"
"Output:\n"
"  An R program piped to stdout\n"
"\n"
"\v"
"Standard and RGL outputs can be generated.\n"
"\n"
;

static error_t parse_plot_command(int key, char* arg, struct argp_state* state) {
  struct arg_plot* arguments = (struct arg_plot*)state->input;

  assert( arguments );

  switch(key) {
    case 'd':
      arguments->device = arg;
      break;

    case 's':
      if (not file_exists(arg)) {
        cerr << on_red << "TCA surface file does not exist: " << bold << arg << reset << endl;
        exit(EXIT_FAILURE);
      }
      arguments->surface_file = arg;
      break;

    case 't':
      if (not file_exists(arg)) {
        cerr << on_red << "Simplex trace file does not exist: " << bold << arg << reset << endl;
        exit(EXIT_FAILURE);
      }
      arguments->trace_file = arg;
      break;

    case 'o':
      arguments->plot_file = arg;
      break;

    case 'a':
      char x, y;
      if (sscanf(arg, "%c,%c", &x, &y) == 2) {
        if (x == y) {
          cerr << on_red << "Lens parameter names can't be the same: " << bold << arg << reset << endl;
          exit(EXIT_FAILURE);
        }
        else {
          switch (x) {
            case 'a': arguments->argx = 0; break;
            case 'b': arguments->argx = 1; break;
            case 'c': arguments->argx = 2; break;
            case 'd': arguments->argx = 3; break;
            default:
              cerr << on_red << "Expecting a lens parameter name (a, b, c, or d). Got this instead: " <<
                bold << string(1, x) << reset << on_red << " in " << bold << arg << reset << endl;
              exit(EXIT_FAILURE);
          }
          switch (y) {
            case 'a': arguments->argy = 0; break;
            case 'b': arguments->argy = 1; break;
            case 'c': arguments->argy = 2; break;
            case 'd': arguments->argy = 3; break;
            default:
              cerr << on_red << "Expecting a lens parameter name (a, b, c, or d). Got this instead: " <<
                bold << string(1, y) << reset << on_red << " in " << bold << arg << reset << endl;
              exit(EXIT_FAILURE);
          }
        }
      }
      else {
        cerr << on_red << "Expecting two lens parameter names 'x,y' where x and y are one of"
          " (a, b, c, or d). Got this instead: " << bold << arg << reset << endl;
        exit(EXIT_FAILURE);
      }
      break;

    case 'c':
      arguments->contour = true;
      break;

    case 'v':
      arguments->vertices = true;
      break;

    case 'z':
      if (arguments->contour) {
        double x, y;
        arguments->minimizer = true;
        arguments->z[0] = arguments->z[1] = arguments->z[2] = arguments->z[3] = arguments->f = 0;
        if (sscanf(arg, "%lf,%lf", &x, &y) == 2) {
          arguments->minimizer_short = true;
          if (arguments->argx == 0) {
            if (x >= 0.9 and x <= 1.1 and y >= -0.1 and y <= 0.1) {
              arguments->z[arguments->argx] = x;
              arguments->z[arguments->argy] = y;
            }
            else if (y < -0.1 or y > 0.1) {
              cerr << on_red << "Expecting a lens distortion parameter in the range [-0.1 .. 0.1]. Got this instead: " <<
                bold << y << reset << on_red << " in " << bold << arg << reset << endl;
              exit(EXIT_FAILURE);
            }
            else {
              cerr << on_red << "Expecting lens distortion parameter (a) in the range [0.9 .. 1.1]. Got this instead: " <<
                bold << x << reset << on_red << " in " << bold << arg << reset << endl;
              exit(EXIT_FAILURE);
            }
          }
          if (arguments->argy == 0) {
            if (y >= 0.9 and y <= 1.1 and x >= -0.1 and x <= 0.1) {
              arguments->z[arguments->argx] = x;
              arguments->z[arguments->argy] = y;
            }
            else if (x < -0.1 or x > 0.1) {
              cerr << on_red << "Expecting a lens distortion parameter in the range [-0.1 .. 0.1]. Got this instead: " <<
                bold << x << reset << on_red << " in " << bold << arg << reset << endl;
              exit(EXIT_FAILURE);
            }
            else {
              cerr << on_red << "Expecting lens distortion parameter (a) in the range [0.9 .. 1.1]. Got this instead: " <<
                bold << y << reset << on_red << " in " << bold << arg << reset << endl;
              exit(EXIT_FAILURE);
            }
          }
        }
        else {
          cerr << on_red << "Expecting two lens parameter names 'x,y' where x and y are one of"
            " (a, b, c, or d). Got this instead: " << bold << arg << reset << endl;
          exit(EXIT_FAILURE);
        }
      }
      break;

    case ARGP_KEY_ARG: // non-option argument
      argp_usage (state);

      break;

    default:
      return ARGP_ERR_UNKNOWN;
  }

  return 0;
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmissing-field-initializers"
static struct argp_option options_plot[] = {
  {"surface-file", 's',  "FILE",             0, "TCA surface data in tabular form <a, b, c, d, TCA>" },
  {"trace-file",   't',  "FILE",             0, "optimization trace in JSON format" },
  {"output-file",  'o',  "FILE",             0, "the plot file to be written by the generated R program (basename)" },
  {"device",       'd',  "DEV",              0, "R graphics device for plotting" },
  {"args",         'a',  "X,Y",              0, "specify co-ordinates (a, b), (a, c), etc., for the projection to plot" },
  {"contour",      'c',  0,                  0, "do a contour plot instead of a 3d surface plot" },
  {"vertices",     'v',  0,                  0, "plot all simplex vertices" },
  {"minimizer",    'z',  "[X,Y|a,b,c,d,f]",  0, "mark the locating of the minimizer. In the case of -c, same order of X and Y as in -a" },
  { 0 }
};

static struct argp argp_plot = {
  options_plot,
  parse_plot_command,
  args_doc_plot,
  doc_plot
};
#pragma GCC diagnostic pop

#define PARSE_ARGS_PLOT \
  struct arg_plot args; \
  int    argc = state->argc - state->next + 1; \
  char** argv = &state->argv[state->next - 1]; \
  char*  argv0 =  argv[0]; \
  argv[0] = (char *)malloc(strlen((char *)(state->name)) + strlen(" plot") + 1); \
  if (!argv[0]) argp_failure(state, 1, ENOMEM, 0); \
  sprintf(argv[0], "%s plot", state->name); \
  args.argx = 0; /* lens param a */ \
  args.argy = 1; /* lens param b */ \
  args.surface_file = NULL; \
  args.trace_file = NULL; \
  args.contour = false; \
  args.vertices = false; \
  args.minimizer = false; \
  args.minimizer_short = false; \
  args.device = "png"; \
  args.plot_file = "tca-surface"; \
  argp_parse(&argp_plot, argc, argv, ARGP_IN_ORDER, &argc, &args); \
  free(argv[0]); \
  argv[0] = argv0; \
  state->next += argc - 1; \

