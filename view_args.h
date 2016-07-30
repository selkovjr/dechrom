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
  char* input_file;
  char* output_file;
  bool r, g, b;
  bool exr;
  int size;
  int center_x, center_y;
  int n_x, n_y;
  int ne_x, ne_y;
  int e_x, e_y;
  int se_x, se_y;
  int s_x, s_y;
  int sw_x, sw_y;
  int w_x, w_y;
  int nw_x, nw_y;
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
        cerr << on_red << "Expecting an integer tile size. Got this instead: " << bold << arg << reset << endl;
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

    case '0':
      if (sscanf(arg, "%d,%d", &(arguments->center_x), &(arguments->center_y)) != 2) {
        cerr << on_red << "Expecting integer offests for central patch (X,Y). Got this instead: " << bold << arg << reset << endl;
        exit(EXIT_FAILURE);
      }
      break;

    case '1':
      if (sscanf(arg, "%d,%d", &(arguments->n_x), &(arguments->n_y)) != 2) {
        cerr << on_red << "Expecting integer offests for north patch (X,Y). Got this instead: " << bold << arg << reset << endl;
        exit(EXIT_FAILURE);
      }
      break;

    case '2':
      if (sscanf(arg, "%d,%d", &(arguments->ne_x), &(arguments->ne_y)) != 2) {
        cerr << on_red << "Expecting integer offests for northeast patch (X,Y). Got this instead: " << bold << arg << reset << endl;
        exit(EXIT_FAILURE);
      }
      break;

    case '3':
      if (sscanf(arg, "%d,%d", &(arguments->e_x), &(arguments->e_y)) != 2) {
        cerr << on_red << "Expecting integer offests for east patch (X,Y). Got this instead: " << bold << arg << reset << endl;
        exit(EXIT_FAILURE);
      }
      break;

    case '4':
      if (sscanf(arg, "%d,%d", &(arguments->se_x), &(arguments->se_y)) != 2) {
        cerr << on_red << "Expecting integer offests for southeast patch (X,Y). Got this instead: " << bold << arg << reset << endl;
        exit(EXIT_FAILURE);
      }
      break;

    case '5':
      if (sscanf(arg, "%d,%d", &(arguments->s_x), &(arguments->s_y)) != 2) {
        cerr << on_red << "Expecting integer offests for south patch (X,Y). Got this instead: " << bold << arg << reset << endl;
        exit(EXIT_FAILURE);
      }
      break;

    case '6':
      if (sscanf(arg, "%d,%d", &(arguments->sw_x), &(arguments->sw_y)) != 2) {
        cerr << on_red << "Expecting integer offests for southwest patch (X,Y). Got this instead: " << bold << arg << reset << endl;
        exit(EXIT_FAILURE);
      }
      break;

    case '7':
      if (sscanf(arg, "%d,%d", &(arguments->w_x), &(arguments->w_y)) != 2) {
        cerr << on_red << "Expecting integer offests for west patch (X,Y). Got this instead: " << bold << arg << reset << endl;
        exit(EXIT_FAILURE);
      }
      break;

    case '8':
      if (sscanf(arg, "%d,%d", &(arguments->nw_x), &(arguments->nw_y)) != 2) {
        cerr << on_red << "Expecting integer offests for northwest patch (X,Y). Got this instead: " << bold << arg << reset << endl;
        exit(EXIT_FAILURE);
      }
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
  {"red",        'r',    0,     0,   "display red channel" },
  {"green",      'g',    0,     0,   "display green channel" },
  {"blue",       'b',    0,     0,   "display blue channel" },
  {"size",       's', "NUMBER", 0,   "tile size (default 1024)" },
  {"exr",        'x',    0,     0,   "input images is a  tilted EXR array" },
  {"center",     '0',  "X,Y",   0,   "central patch offsets" },
  {"north",      '1',  "X,Y",   0,   "north patch offsets" },
  {"northeast",  '2',  "X,Y",   0,   "northeast patch offsets" },
  {"east",       '3',  "X,Y",   0,   "east patch offsets" },
  {"southeast",  '4',  "X,Y",   0,   "southeast patch offsets" },
  {"south",      '5',  "X,Y",   0,   "south patch offsets" },
  {"southwest",  '6',  "X,Y",   0,   "southwest patch offsets" },
  {"west",       '7',  "X,Y",   0,   "west patch offsets" },
  {"northwest",  '8',  "X,Y",   0,   "northwest patch offsets" },
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
  args.center_x = 0; \
  args.center_y = 0; \
  args.n_x = 0; \
  args.n_y = 0; \
  args.ne_x = 0; \
  args.ne_y = 0; \
  args.e_x = 0; \
  args.e_y = 0; \
  args.se_x = 0; \
  args.se_y = 0; \
  args.s_x = 0; \
  args.s_y = 0; \
  args.sw_x = 0; \
  args.sw_y = 0; \
  args.w_x = 0; \
  args.w_y = 0; \
  args.nw_x = 0; \
  args.nw_y = 0; \
  argp_parse(&argp_view, argc, argv, ARGP_IN_ORDER, &argc, &args); \
  free(argv[0]); \
  argv[0] = argv0; \
  state->next += argc - 1;

