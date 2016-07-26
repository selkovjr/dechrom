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

      // Prepare interpolated data
      cout <<
        "library(rjson)\n"
        "library(ggplot2)\n"
        "s.data <- read.table('" << args.surface_file << "', header = TRUE)\n"
        "s.loess <- loess(TCA ~ " << param[args.argx] << " * " << param[args.argy] << ", data = s.data, degree = 2, span = 0.25)\n"
        "xmin <- min(s.data$" << param[args.argx] << ")\n"
        "xmax <- max(s.data$" << param[args.argx] << ")\n"
        "ymin <- min(s.data$" << param[args.argy] << ")\n"
        "ymax <- max(s.data$" << param[args.argy] << ")\n"
        "s.fit <- expand.grid(\n"
        "  list(\n"
        "    " << param[args.argx] << " = seq(xmin, xmax, (xmax - xmin) / 100),\n"
        "    " << param[args.argy] << " = seq(ymin, ymax, (ymax - ymin) / 100)\n"
        "  )\n"
        ")\n"
        "z <- predict(s.loess, newdata=s.fit)\n"
        "s.fit$TCA = as.numeric(z)\n"
        "title <- 'Red channel TCA survey'\n";

      // Prepare simplex trace data
      if (args.trace_file != NULL) {
        cout << "t <- fromJSON(paste(readLines('" << args.trace_file << "'), collapse=''))\n"
          "x <- numeric(0)\n"
          "y <- numeric(0)\n"
          "xend <- numeric(0)\n"
          "yend <- numeric(0)\n"
          "for (i in 1:(length(t$simplex) - 1)) {\n"
          "  prev <- t$simplex[[i]]\n"
          "  curr <- t$simplex[[i + 1]]\n"
          " x[i] = prev$x[[prev$index + 1]][" << (args.argx + 1) << "]\n"
          " y[i] = prev$x[[prev$index + 1]][" << (args.argy + 1) << "]\n"
          " xend[i] = curr$x[[curr$index + 1]][" << (args.argx + 1) << "]\n"
          " yend[i] = curr$x[[curr$index + 1]][" << (args.argy + 1) << "]\n"
          "}\n"
          "path <- subset(data.frame(x, y, xend, yend), sqrt((x - xend)^2 + (y - yend)^2) > 0.00025)\n"
          "title <- 'Red channel convergence'\n"
          ;
      }

      // Set output device
      if (string(args.device).compare("png") == 0) {
        cout << "png('" << args.plot_file << ".png', width = 800, height = 800, res=120)\n";
      }

      // Plot the surface
      cout << "p <- ggplot() +\n"
        "  geom_raster(data = s.data, aes(x = " << param[args.argx] << ", y = " << param[args.argy] << ", fill = TCA)) +\n"
        "  scale_fill_gradient(limits = c(min(s.data$TCA), max(s.data$TCA)), low = rgb(0.2, 0.05, 0.05), high = rgb(1, 0.5, 0.5)) +\n"
        "  geom_contour(data = s.fit, aes(x = " << param[args.argx] << ", y = " << param[args.argy] << ", z = TCA), color = 'white', alpha = 0.2, bins = 20) +\n";

      // Plot the trace
      if (args.trace_file != NULL) {
        cout <<
          "  geom_segment(data = path, aes(x = x, y = y, xend = xend, yend = yend),\n"
          "    size = 0.75, color = rgb(1, 0.6, 0.6, 0.5), arrow = arrow(length = unit(2, 'mm'), angle = 20)\n"
          "  ) +\n";
      }

      // Finish the plot
      cout <<
        "  ggtitle(title) +\n"
        "  xlab('" << param[args.argx] << "') + ylab('" << param[args.argy] << "') +\n"
        "  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0))\n";

      // Render
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
