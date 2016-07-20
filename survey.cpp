#include <iostream>
#include <sstream>
#include <iomanip>
#include <functional>
#include <string>

#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include <math.h>
#include <stdio.h>
#include <unistd.h>

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <error.h>
#include <ctime>

#include "survey_args.h"
#include "cfa_mask.h"
#include "write_tiff.h"
#include "progressbar.h"
#include "termcolor.h"


using namespace std;
using namespace cv;
using namespace termcolor;

bool verbose;
bool exr_mode;

int run = 0;

Mat ref_plane, work_plane;
Mat blank, mask;
Mat dst, map_x, map_y;

#define RED 0
#define GREEN 1
#define BLUE 2
#define UNKNOWN 3

#define CFA_WIDTH 3264
#define CFA_HEIGHT 2464

int work_plane_color = UNKNOWN;

string to_string(double v) {
  stringstream ss;
  ss << v;
  return ss.str();
}

void color_switch(ostream& s, int color) {
  switch (color) {
    case RED:
      s << lightred;
      break;
    case GREEN:
      s << lightgreen;
      break;
    case BLUE:
      s << lightblue;
      break;
    default:
      s << yellow;
  }
}

double diff (double a, double b, double c, double d) {
  // private thread variables
  long i, j, pa, pb;
  double x, y, r, rd;

  clock_t start_time = clock(), end_time;
  double elapsed;

  // Fill the CFA area with white
  if (verbose) cerr << right << setw(5) << run << reset << lightgrey <<  ": making mask ... " << reset;
  blank.create( work_plane.size(), work_plane.type() );
#pragma omp parallel shared(blank, exr_mode) private(i, j)
  {
#pragma omp for collapse(2)
    for (j = 0; j < blank.rows; j++ ) {
      for (i = 0; i < blank.cols; i++ ) {
        if (exr_mode) {
          if ( // CFA interior
            i + j >= CFA_WIDTH - 1 &&                   // NW boundary
            j > i - CFA_WIDTH - 1 &&                    // NE boundary
            i + j < CFA_WIDTH + 2 * CFA_HEIGHT - 1 &&   // SE boundary
            i > j - CFA_WIDTH                           // SW boundary
          ) {
            blank.at<ushort>(j,i) = 65535;
          }
          else {
            blank.at<ushort>(j,i) = j;
          }
        }
        else {
          blank.at<ushort>(j,i) = 65535;
        }
      }
    }
  } /* end of parallel section */


  /// Create dst, map_x and map_y with the same size as work_plane:
  dst.create( work_plane.size(), work_plane.type() );
  map_x.create( work_plane.size(), CV_32FC1 );
  map_y.create( work_plane.size(), CV_32FC1 );

  run++;
  if (verbose) cerr << lightgrey <<  "computing maps ... " << reset;
#pragma omp parallel shared(map_x, map_y, a, b, c, d) private(i, j, x, y, r, rd)
  {
#pragma omp for collapse(2)
    for (j = 0; j < work_plane.rows; j++ ) {
      for (i = 0; i < work_plane.cols; i++ ) {
        if (blank.at<ushort>(j, i) == 65535) { // don't move masked pixels (black outer triangles)
          x = (double)(2 * i - work_plane.cols) / work_plane.cols;
          y = (double)(2 * j - work_plane.rows) / work_plane.rows;
          r = sqrt(x * x + y * y);
          rd = a * r + b * pow(r, 2) + c * pow(r, 3) + d * pow(r, 4);
          map_x.at<float>(j, i) = (x * rd / r + 1) * work_plane.cols / 2;
          map_y.at<float>(j, i) = (y * rd / r + 1) * work_plane.rows / 2;
        }
      }
    }
  } /* end of parallel section */

  if (verbose) cerr << lightgrey << "remapping ... " << reset;
  remap( blank, mask, map_x, map_y, CV_INTER_LINEAR, BORDER_CONSTANT, Scalar(0,0, 0) );
  remap( work_plane, dst, map_x, map_y, CV_INTER_LINEAR, BORDER_CONSTANT, Scalar(0,0, 0) );

  /// Compute difference
  if (verbose) cerr << lightgrey << "computing difference ... " << reset;
  unsigned long num = 0;
  unsigned long quot = 0;
#pragma omp parallel shared(dst, ref_plane) private(i, j, pa, pb) reduction(+:num,quot)
  {
#pragma omp for collapse(2)
    for ( j = 0; j < work_plane.rows; j++ ) {
      for ( i = 0; i < work_plane.cols; i++ ) {
        if (dst.at<ushort>(i, j) > 0 and ref_plane.at<ushort>(i,j) > 0) { // select unmasked
          pa = dst.at<ushort>(j, i);
          pb = ref_plane.at<ushort>(j,i);
          num += abs((double)pa - (double)pb);
          quot += 1;
        }
      }
    }
  } /* end of parallel section */

  double tca = 1.0 * num / quot;
  if (verbose) {
    end_time = clock();
    elapsed = double(end_time - start_time) / CLOCKS_PER_SEC;
    cerr << "\r" << flush;
    cerr << setw(100) << left << 0 << "\r" << flush; // blank line
    cerr << lightgrey << right << setw(5) << run << reset <<
      setw(10) << setprecision(6) << a <<
      setw(10) << setprecision(6) << b <<
      setw(10) << setprecision(6) << c <<
      setw(10) << setprecision(6) << d;
    color_switch(cerr, work_plane_color);
    cerr << setw(9) << setprecision(2) << tca << reset <<
      lightgrey << setw(7) << setprecision(2) << elapsed << "s" << reset << endl;
  }

  return tca;
}

// --------------------------------------------------------------------
void run_survey (struct argp_state* state) {
  PARSE_ARGS_SURVEY;

  verbose = args.verbose;

  double astep = abs(args.amax - args.amin) / (args.nodes - 1);
  double bstep = abs(args.bmax - args.bmin) / (args.nodes - 1);
  double cstep = abs(args.cmax - args.cmin) / (args.nodes - 1);
  double dstep = abs(args.dmax - args.dmin) / (args.nodes - 1);

  unsigned n = 1;
  if (astep) n *= args.nodes;
  if (bstep) n *= args.nodes;
  if (cstep) n *= args.nodes;
  if (dstep) n *= args.nodes;


  cerr << lightgrey << "Computing " << reset << n <<
    lightgrey << " iterations of " << bold << "TCA" << reset << lightgrey << " for " << reset <<

    "a" << lightgrey << (astep ? " in " : " = ") << reset <<
    args.amin << lightgrey << (astep ? ".." : "") << reset <<
    (astep ? to_string(args.amax) : "") << reset << lightgrey << ", " << reset <<

    "b" << lightgrey << (bstep ? " in " : " = ") << reset <<
    args.bmin << lightgrey << (bstep ? ".." : "") << reset <<
    (bstep ? to_string(args.bmax) : "") << reset << lightgrey << ", " << reset <<

    "c" << lightgrey << (cstep ? " in " : " = ") << reset <<
    args.cmin << lightgrey << (cstep ? ".." : "") << reset <<
    (cstep ? to_string(args.cmax) : "") << reset << lightgrey << ", " << reset <<

    "d" << lightgrey << (dstep ? " in " : " = ") << reset <<
    args.dmin << lightgrey << (dstep ? ".." : "") << reset <<
    (dstep ? to_string(args.dmax) : "") << reset << lightgrey << reset <<

    endl;

  // Load the images
  cerr << lightgrey << "loading reference plane " << reset << args.ref_file << "\r" << flush;
  ref_plane = imread( args.ref_file, CV_LOAD_IMAGE_ANYDEPTH ); // grayscale
  cerr << setw(100) << left << 0 << "\r" << flush; // blank line
  cerr << lightgrey << "reference plane: " << reset << lightgreen << args.ref_file << reset << endl;

  cerr << lightgrey << "loading work plane " << reset << args.work_plane_file << "\r" << flush;
  work_plane = imread( args.work_plane_file, CV_LOAD_IMAGE_ANYDEPTH );
  cerr << setw(100) << left << 0 << "\r" << flush; // blank line

  if (ref_plane.cols != work_plane.cols or ref_plane.rows != work_plane.rows) {
    cerr << on_red << "Input images must have identical geometry. Got "
      << bold << ref_plane.rows << reset << on_red << "×" << bold << ref_plane.cols << reset
      << on_red << " and "
      << bold << work_plane.rows << reset << on_red << "×" << bold << work_plane.cols << reset
      << endl;
    exit(EXIT_FAILURE);
  }

  if (ref_plane.cols == CFA_WIDTH and ref_plane.rows == CFA_HEIGHT) {
    exr_mode = false;
  }
  else if (
    ref_plane.rows == CFA_WIDTH + CFA_HEIGHT and
    ref_plane.cols == CFA_WIDTH + CFA_HEIGHT
  ) {
    exr_mode = true;
  }
  else {
    cerr << on_red << "Unrecognized image geometry: "
      << bold << ref_plane.rows << reset << on_red << "×" << bold << ref_plane.cols << reset
      << endl;
    cerr << lightgrey << "Expected image formats are either " << reset
      << bold << (CFA_WIDTH + CFA_HEIGHT) << reset << "×" << bold << (CFA_WIDTH + CFA_HEIGHT) << reset
      << lightgrey << " (rotated high-resolution EXR matrix) or "
      << bold << CFA_WIDTH << reset << "×" << bold << CFA_HEIGHT << reset
      << lightgrey << " (a single EXR subframe)" << reset << endl;
    exit(EXIT_FAILURE);
  }

  string work_plane = args.work_plane_file;
  if (work_plane.find("-0.") != string::npos) work_plane_color = RED;
  if (work_plane.find("-1.") != string::npos) work_plane_color = GREEN;
  if (work_plane.find("-2.") != string::npos) work_plane_color = BLUE;
  cerr << lightgrey << "work plane: " << reset;
  color_switch(cerr, work_plane_color);
  cerr << work_plane << reset;

  if (work_plane_color != RED and work_plane_color != GREEN and work_plane_color != BLUE) {
    cerr << lightgrey << " (filename does not indicate which channel it is)";
  }
  cerr << endl;

  progressbar *pbar = progressbar_new("", n);
  cerr.setf(ios::fixed, ios::floatfield);

  if (not verbose) {
    color_switch(cerr, work_plane_color);
  }
  for (unsigned ac = 0; ac < args.nodes; ac++) {
    double a = args.amin + ac * astep;
    for (unsigned bc = 0; bc < args.nodes; bc++) {
      double b = args.bmin + bc * bstep;
      for (unsigned cc = 0; cc < args.nodes; cc++) {
        double c = args.cmin + cc * cstep;
        for (unsigned dc = 0; dc < args.nodes; dc++) {
          double d = args.dmin + dc * dstep;
          double tca = diff(a, b, c, d);
          printf("%f\t%f\t%f\t%f\t%f\n", a, b, c, d, tca);
          if (not verbose) {
            progressbar_inc(pbar);
          }
          if (args.dmin == args.dmax) break;
        }
        if (args.cmin == args.cmax) break;
      }
      if (args.bmin == args.bmax) break;
    }
    if (args.amin == args.amax) break;
  }

  if (not verbose) {
    progressbar_finish(pbar);
  }
  cerr << reset;

} // run_survey()
