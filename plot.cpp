#include <iostream>
#include <iomanip>

#include "plot_args.h"
#include "termcolor.h"


using namespace std;
using namespace termcolor;

static string param[4] = {"a", "b", "c", "d"};

// --------------------------------------------------------------------
void run_plot (struct argp_state* state) {
  PARSE_ARGS_PLOT;

  if (args.surface_file != NULL) {
    cerr << args.argx << ", " << args.argy << endl;
    if (args.contour) {

      // library(akima)
      // reggrid <- interp(s.data$a, s.data$b, s.data$TCA, linear=T,extrap=F)
      // plot(s.data$a, s.data$b)
      // image  (reggrid,add=TRUE)
      // contour(reggrid,add=TRUE)

      // > ggplot(s.fit, aes(a, b, fill = TCA)) + geom_tile() + ggtitle('TCA surface for the red channel') + scale_fill_gradient(limits = c(min(z), max(z)), low = 'black', high = 'white') + scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0))

      cout << "library(ggplot2)\n";
      cout << "s.data <- read.table('" << args.surface_file << "', header = TRUE)\n";
      cout << "s.loess <- loess(TCA ~ " << param[args.argx] << " * " << param[args.argy] <<
        ", data = s.data, degree = 2, span = 0.25)\n";
      cout << "xmin <- min(s.data$" << param[args.argx] << ")\n";
      cout << "xmax <- max(s.data$" << param[args.argx] << ")\n";
      cout << "ymin <- min(s.data$" << param[args.argy] << ")\n";
      cout << "ymax <- max(s.data$" << param[args.argy] << ")\n";
      cout << "s.fit <- expand.grid(\n"
        "  list(\n"
        "    " << param[args.argx] << " = seq(xmin, xmax, (xmax - xmin) / 100),\n"
        "    " << param[args.argy] << " = seq(ymin, ymax, (ymax - ymin) / 100)\n"
        "  )\n"
        ")\n";
      cout << "z <- predict(s.loess, newdata=s.fit)\n";
      cout << "s.fit$TCA = as.numeric(z)\n";

      if (args.device.compare("png") == 0) {
        cout << "png('" << args.plot_file << ".png', width = 800, height = 800, res=120)\n";
      }
      cout << "p <- ggplot(s.fit, aes(a, b, fill = TCA)) + geom_tile() +\n"
        "  xlab('a') + ylab('b') +\n"
        "  ggtitle('TCA metric for the red channel') +\n"
        "  scale_fill_gradient(limits = c(min(z), max(z)), low = 'black', high = 'white') +\n"
        "  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0))\n";
      cout << "print(p)\n";

    }
    else {
    }
  }
  else {
    cerr << on_red << "Required parameter -s for TCA surface file name is missing" << reset << endl;
    exit(EXIT_FAILURE);
  }


} // run_survey()
