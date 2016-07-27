#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>

#include "plot_args.h"
#include "termcolor.h"


using namespace std;
using namespace termcolor;

static string param[4] = {"a", "b", "c", "d"};

// --------------------------------------------------------------------
void run_plot (struct argp_state* state) {
  PARSE_ARGS_PLOT;


  if (args.surface_file != NULL) {
    vector<double> surf_min {1e8, 1e8, 1e8, 1e8, 1e8};
    vector<double> surf_max {-1e8, -1e8, -1e8, -1e8, -1e8};
    ifstream ifs (args.surface_file);
    string buf;
    string title;

    while (getline(ifs, buf)) {
      int i = 0;
      double d;
      stringstream sbuf(buf);
      while (sbuf >> d) {
        if (surf_min[i] > d) surf_min[i] = d;
        if (surf_max[i] < d) surf_max[i] = d;
        i++;
      }
    }

    if (args.contour) {
      title = "Reg-green TCA survey";
      // Prepare interpolated data
      cout <<
        "# min a = " << surf_min[0] << "\n"
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
        ;

      // Prepare simplex trace data
      if (args.trace_file != NULL) {
        title = "Red-green convergence";

        cout << "t <- fromJSON(paste(readLines('" << args.trace_file << "'), collapse=''))\n"
          "x <- numeric(0)\n"
          "y <- numeric(0)\n"
          "xend <- numeric(0)\n"
          "yend <- numeric(0)\n"
          "sx <- numeric(0)\n"
          "sy <- numeric(0)\n"
          "sxend <- numeric(0)\n"
          "syend <- numeric(0)\n"
          "for (i in 1:(length(t$simplex) - 1)) {\n"
          "  prev <- t$simplex[[i]]\n"
          "  curr <- t$simplex[[i + 1]]\n"
          "  x[i] = prev$x[[prev$index + 1]][" << (args.argx + 1) << "]\n"
          "  y[i] = prev$x[[prev$index + 1]][" << (args.argy + 1) << "]\n"
          "  xend[i] = curr$x[[curr$index + 1]][" << (args.argx + 1) << "]\n"
          "  yend[i] = curr$x[[curr$index + 1]][" << (args.argy + 1) << "]\n"
          "  ix = 0\n"
          "  for (j in 1:5) {\n"
          "    if (j != prev$index + 1) {\n"
          "      ix <- ix + 1\n"
          "        sx[4 * (i - 1) + ix] = x[i]\n"
          "        sy[4 * (i - 1) + ix] = y[i]\n"
          "        sxend[4 * (i - 1) + ix] <- prev$x[[j]][1]\n"
          "        syend[4 * (i - 1) + ix] <- prev$x[[j]][2]\n"
          "    }\n"
          "  }\n"
          "}\n"
          "path <- subset(data.frame(x, y, xend, yend), sqrt((x - xend)^2 + (y - yend)^2) > 0.00025)\n"
          "simplex <- data.frame(sx, sy, sxend, syend)\n"
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
        if (args.vertices) {
          cout <<
            "  geom_segment(data = simplex, aes(x = sx, y = sy, xend = sxend, yend = syend),\n"
            "    size = 0.5, color = rgb(0.8, 0.4, 0.4, 0.2)\n"
            "  ) +\n"
            "  geom_point(data = simplex, aes(x = sxend, y = syend),\n"
            "    size = 3, color = rgb(1, 0.6, 0.6, 0.08), shape = 21\n"
            "  ) +\n";
        }
        cout <<
          "  geom_curve(data = path, aes(x = x, y = y, xend = xend, yend = yend),\n"
          "    size = 0.75, color = rgb(1, 0.4, 0.4, 0.6), arrow = arrow(length = unit(2, 'mm'), angle = 20),\n"
          "    curvature = 0.3\n"
          "  ) +\n";
      }

      // Plot the solution
      if (args.minimizer) {
        title = "Red channel solution";
        cout << "  geom_vline(xintercept = " << args.z[args.argx] << ", size = 0.5, color = 'white', alpha = 0.5, linetype = 3) +\n";
        cout << "  geom_hline(yintercept = " << args.z[args.argy] << ", size = 0.5, color = 'white', alpha = 0.5, linetype = 3) +\n";
        cout.setf(ios::fixed, ios::floatfield);
        for (int i = 0; i < 4; i++) {
          double offset = (i + 1) * 1.5;
          double value = args.z[i];
          if (args.minimizer_short and not (i == args.argx or i == args.argy)) {
            if (surf_min[i] == surf_max[i]) {
              value = surf_min[i];
            }
            else {
              cerr << on_red << "In " << bold << args.surface_file << reset << on_red <<
                ", min(" << param[i] << ") ≠ max(" << param[i] <<
                "). I only know how to handle 2d survey data." << reset << endl;
              exit(EXIT_FAILURE);
            }
          }
          cout << "  geom_text(aes(x = Inf, y = Inf, label = '" << param[i] << " = " << setw(10) << setprecision(6) << value << "'), \n"
            "size = 3, color = 'white', alpha = 0.8, hjust = 1.1, vjust = " << offset << ") +\n";
        }
      }
      else {
        // Determine ranges or points for the invisible dimensions
        double offset = 0.5;
        for (int i = 0; i < 4; i++) {
          if (i != args.argx and i != args.argy) {
            offset += 1.5;
            stringstream range;
            if (surf_min[i] == surf_max[i]) {
              range << setprecision(6) << surf_min[i];
            }
            else {
              range << setprecision(6) << surf_min[i] << " .. " << setprecision(6) << surf_max[i];
            }
            cout << "  geom_text(aes(x = Inf, y = Inf, label = '" << param[i] << " = " << range.str() << "'), \n"
              "size = 3, color = 'white', alpha = 0.8, hjust = 1.1, vjust = " << offset << ") +\n";
          }
        }
      }

      // Finish the plot
      cout <<
        "  ggtitle('" << title << "') +\n"
        "  xlab('" << param[args.argx] << "') + ylab('" << param[args.argy] << "') +\n"
        "  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0))\n";

      // Render
      cout << "print(p)\n";
    }
    else {
      // RGL
    }
  }
  else {
    cerr << on_red << "Required parameter -s for TCA surface file name is missing" << reset << endl;
    exit(EXIT_FAILURE);
  }


} // run_survey()
