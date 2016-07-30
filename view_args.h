#include <argp.h>
#include <assert.h>
#include "termcolor.h"
#include <iostream>
#include <iomanip>

#include "util.h"

using namespace std;
using namespace termcolor;

// ------------------------
// ## view command parser
//
struct arg_view {
  int size;
  char* input_file;
  char* output_file;
  bool r, g, b;
  bool exr;
};

static char args_doc_view[] = "input.tiff output.tiff";

static char doc_view[] =
"\n"
"Make a 3x3 tile set of critical areas of a TIFF image to allow\n"
"visual inspection of channel convergence.\n"
"\v"
;

static error_t parse_view_command(int key, char* arg, struct argp_state* state) {
  struct arg_view* arguments = (struct arg_view*)state->input;
  char **nonopt;

  assert( arguments );

  switch(key) {
    case 's':
      if (sscanf(arg, "%d", &(arguments->size)) != 1) {
        cerr << on_red << "Expecting an integer tile size. Got this: instead" << bold << arg << reset << endl;
        exit(EXIT_FAILURE);
      }
      break;

    case 'r':
      arguments->r = true;
      break;

    case 'g':
      arguments->g = true;
      break;

    case 'b':
      arguments->b = true;
      break;

    case 'x':
      arguments->exr = true;
      break;

    case ARGP_KEY_NO_ARGS:
      argp_usage (state);

    case ARGP_KEY_ARG: // non-option argument
      nonopt = &state->argv[state->next];
      state->next = state->argc; // we're done

      arguments->input_file = arg;
      if (not file_exists(arguments->input_file)) {
        cerr << on_red << "File does not exist: " << bold << arguments->input_file << reset << endl;
        exit(EXIT_FAILURE);
      }

      arguments->output_file = nonopt[0];

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
static struct argp_option options_view[] = {
  {"red",   'r',    0,     0,   "display red channel" },
  {"green", 'g',    0,     0,   "display green channel" },
  {"blue",  'b',    0,     0,   "display blue channel" },
  {"size",  's', "NUMBER", 0,   "tile size (default 1024)" },
  {"exr",   'x',    0,     0,   "input images is a  tilted EXR array" },
  { 0 }
};

static struct argp argp_view = {
  options_view,
  parse_view_command,
  args_doc_view,
  doc_view
};
#pragma GCC diagnostic pop

#define PARSE_ARGS_VIEW \
  struct arg_view args; \
  int    argc = state->argc - state->next + 1; \
  char** argv = &state->argv[state->next - 1]; \
  char*  argv0 =  argv[0]; \
  argv[0] = (char *)malloc(strlen((char *)(state->name)) + strlen(" view") + 1); \
  if (!argv[0]) argp_failure(state, 1, ENOMEM, 0); \
  sprintf(argv[0], "%s view", state->name); \
  args.exr = false; \
  args.r = false; \
  args.g = false; \
  args.b = false; \
  args.size = 1024; \
  argp_parse(&argp_view, argc, argv, ARGP_IN_ORDER, &argc, &args); \
  free(argv[0]); \
  argv[0] = argv0; \
  state->next += argc - 1;

