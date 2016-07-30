#include <argp.h>
#include <assert.h>
#include "termcolor.h"
#include <iostream>
#include <iomanip>

#include "util.h"

using namespace std;
using namespace termcolor;

// ------------------------
// ## radial command parser
//
struct arg_radial {
  double a, b, c, d;
  double shift_x, shift_y;
  char* input_file;
  char* output_file;
  bool fast;
};

static char args_doc_radial[] = "a b c d input.tiff output.tiff";

static char doc_radial[] =
"\n"
"Apply radial lens distortion to a grayscale TIFF image\n"
"\n"
"\v"
"The default remapping method is ImageMagick's Lanczos. It runs\n"
"slower but appears to result in a slightly more pleasing image.\n"
"\n"
"The faster alternative (-f) is OpenCV Lanczos. It produces visually\n"
"similar results, but runs several times faster than the default method.\n"
"\n"
"This tool was built to distort separate color planes of an image\n"
"and it expects grayscale input. However, the default ImageMagick\n"
"method will happily distort color images. The OpenCV variant will\n"
"accept color input but will produce a grayscale image of unkownn\n"
"formulation.\n"
;

static error_t parse_radial_command(int key, char* arg, struct argp_state* state) {
  struct arg_radial* arguments = (struct arg_radial*)state->input;
  char **nonopt;

  assert( arguments );

  switch(key) {
    case 'f':
      arguments->fast = true;
      break;

    case 'x':
      if (sscanf(arg, "%lf", &(arguments->shift_x)) != 1) {
        cerr << on_red << "Expecting horizontal shift amount in pixels. Got this: "
          << bold << arg << reset << endl;
        exit(EXIT_FAILURE);
      }
      break;

    case 'y':
      if (sscanf(arg, "%lf", &(arguments->shift_y)) != 1) {
        cerr << on_red << "Expecting horizontal shift amount in pixels. Got this: "
          << bold << arg << reset << endl;
        exit(EXIT_FAILURE);
      }
      break;

    case ARGP_KEY_NO_ARGS:
      argp_usage (state);

    case ARGP_KEY_ARG: // non-option argument
      nonopt = &state->argv[state->next];
      state->next = state->argc; // we're done

      if (sscanf(arg, "%lf", &(arguments->a)) != 1) {
        cerr << on_red << "Expecting a floating-point value for (a); got this: " << bold << arg << reset << endl;
        exit(EXIT_FAILURE);
      }
      if (sscanf(nonopt[0], "%lf", &(arguments->b)) != 1) {
        cerr << on_red << "Expecting a floating-point value for (b); got this: " << bold << arg << reset << endl;
        exit(EXIT_FAILURE);
      }
      if (sscanf(nonopt[1], "%lf", &(arguments->c)) != 1) {
        cerr << on_red << "Expecting a floating-point value for (c); got this: " << bold << arg << reset << endl;
        exit(EXIT_FAILURE);
      }
      if (sscanf(nonopt[2], "%lf", &(arguments->d)) != 1) {
        cerr << on_red << "Expecting a floating-point value for (d); got this: " << bold << arg << reset << endl;
        exit(EXIT_FAILURE);
      }

      arguments->input_file = nonopt[3];
      if (not file_exists(arguments->input_file)) {
        cerr << on_red << "File does not exist: " << bold << arguments->input_file << reset << endl;
        exit(EXIT_FAILURE);
      }

      arguments->output_file = nonopt[4];

      break;

    case ARGP_KEY_END:
      if (state->arg_num < 6) {
        argp_error(state, "Not enough arguments");
      }
      if (state->arg_num > 6) {
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
static struct argp_option options_radial[] = {
  {"fast",      'f',    0,     0, "trade some quality for speed by using OpenCV remap()" },
  {"shift-x",   'x', "PIXELS", 0, "shift the final distorted image horizontally" },
  {"shift-y",   'y', "PIXELS", 0, "shift the final distorted image vertically" },
  { 0 }
};

static struct argp argp_radial = {
  options_radial,
  parse_radial_command,
  args_doc_radial,
  doc_radial
};
#pragma GCC diagnostic pop

#define PARSE_ARGS_RADIAL \
  struct arg_radial args; \
  int    argc = state->argc - state->next + 1; \
  char** argv = &state->argv[state->next - 1]; \
  char*  argv0 =  argv[0]; \
  argv[0] = (char *)malloc(strlen((char *)(state->name)) + strlen(" radial") + 1); \
  if (!argv[0]) argp_failure(state, 1, ENOMEM, 0); \
  sprintf(argv[0], "%s radial", state->name); \
  args.fast = false; \
  args.shift_x = 0.0; \
  args.shift_y = 0.0; \
  argp_parse(&argp_radial, argc, argv, ARGP_IN_ORDER, &argc, &args); \
  free(argv[0]); \
  argv[0] = argv0; \
  state->next += argc - 1;

