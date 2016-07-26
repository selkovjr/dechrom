#include <argp.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#include "subcommands.h"

const char *argp_program_version = VERSION " " __DATE__;


//-------------------------------------------------------------------
// Command-line option parsing.
// Example at https://gist.github.com/sam-github/57f49711cd9073b35d22
//

// ## Top-level parser
static char doc_toplevel[1000];

static error_t parse_toplevel (int key, char *arg, struct argp_state *state) {
  switch (key) {
    case ARGP_KEY_ARG:
      assert( arg );
      if(strcmp(arg, "survey") == 0) {
        run_survey(state);
      }
      else if(strcmp(arg, "find") == 0) {
        run_find(state);
      }
      else if(strcmp(arg, "plot") == 0) {
        run_plot(state);
      }
      else if(strcmp(arg, "radial") == 0) {
        run_radial(state);
      }
      else {
        argp_error(state, "%s is not a valid command", arg);
      }
      break;

    case ARGP_KEY_END:
      if (state->arg_num < 2)
        argp_usage (state);
      break;

    default:
      return ARGP_ERR_UNKNOWN;
  }
  return 0;
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmissing-field-initializers"
static struct argp argp_toplevel = {
  0, // no options
  parse_toplevel,
  "<COMMAND> <ARGS>",
  doc_toplevel
};
#pragma GCC diagnostic pop


int main(int argc, char **argv) {
  sprintf(
    doc_toplevel,
    "\n"
    "Utilities for correcting transverse chromatic aberration\n"
    "Version: %s\n"
    "\n"
    "Command: survey  sample TCA in a volume of parameter space\n"
    "         find    find optimal distortion parameters\n"
    "         radial  apply a radial distortion to an image\n"
    "\n",
    argp_program_version
  );

  argp_parse (&argp_toplevel, argc, argv, ARGP_IN_ORDER, NULL, NULL);

  return(EXIT_SUCCESS);
}

