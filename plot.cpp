#include <iostream>
#include <iomanip>

#include "plot_args.h"
#include "termcolor.h"


using namespace std;
using namespace termcolor;

// --------------------------------------------------------------------
void run_plot (struct argp_state* state) {
  PARSE_ARGS_PLOT;

  cerr << args.argx << ", " << args.argy << endl;

} // run_survey()
