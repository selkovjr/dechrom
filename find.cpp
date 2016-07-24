#include <iostream>
#include <iomanip>
#include <ctime>

#include <cppoptlib/meta.h>
#include <cppoptlib/problem.h>
#include <cppoptlib/solver/bfgssolver.h>
#include <cppoptlib/solver/conjugatedgradientdescentsolver.h>
#include <cppoptlib/solver/newtondescentsolver.h>
#include <cppoptlib/solver/neldermeadsolver.h>
#include <cppoptlib/solver/lbfgssolver.h>
#include <cppoptlib/solver/cmaessolver.h>

#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include <math.h>
#include <stdio.h>
#include <omp.h>

#include "find_args.h"
#include "cfa_mask.h"
#include "progressbar.h"
#include "termcolor.h"

using namespace std;
using namespace cv;
using namespace termcolor;

static bool verbose;
static bool exr_mode;

static int run = 0;
static int work_plane_color = UNKNOWN;

static Mat ref_plane, work_plane;
static Mat blank, mask;
static Mat dst, map_x, map_y;

static progressbar *pbar;
static int iteration = 0;
static bool first_iteration = true;
static bool trace_simplices = false;
static ofstream *trace_stream;

static double diff (double a, double b, double c, double d);


namespace cppoptlib {

  // Define the target function (Lens::value())
  template<typename T> class Lens : public Problem<T, 4> {
    public:
      using typename cppoptlib::Problem<T, 4>::Scalar;
      using typename cppoptlib::Problem<T, 4>::TVector;
      using MatrixType = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;


      // this is just the objective (NOT optional)
      T value(const TVector &x) {
        return diff(x[0], x[1], x[2], x[3]);
      }

      bool detailed_callback(const cppoptlib::Criteria<T> &state, SimplexOp op, int index, const MatrixType &x, vector<Scalar> f) {
        iteration = (int)state.iterations;
        TVector xp = x.col(index).transpose();
        if (verbose) {
          cerr << green << setw(12) << op << lightgrey <<
            setw(10) << setprecision(6) << xp[0] <<
            setw(10) << setprecision(6) << xp[1] <<
            setw(10) << setprecision(6) << xp[2] <<
            setw(10) << setprecision(6) << xp[3];
          color_switch(cerr, work_plane_color);
          cerr << bold << setw(9) << setprecision(2) << f[index] << reset << endl << endl;
        }
        else {
          color_switch(cerr, work_plane_color);
          progressbar_inc(pbar);
        }
        printf("%f\t%f\t%f\t%f\t%f\n", xp[0], xp[1], xp[2], xp[3], f[index]);
        if (trace_simplices) {
          TVector x0 = x.col(0).transpose();
          TVector x1 = x.col(1).transpose();
          TVector x2 = x.col(2).transpose();
          TVector x3 = x.col(3).transpose();
          TVector x4 = x.col(4).transpose();
          *trace_stream << (first_iteration ? "" : ",\n") <<
            "    {\n"
            "      \"iter\": " << state.iterations << ",\n"
            "      \"op\": \"" << op << "\",\n"
            "      \"index\": " << index << ",\n"
            "      \"x\": [\n"
            "        [" << x0[0] << ", " << x0[1] << ", " << x0[2] << ", " << x0[3] << "],\n"
            "        [" << x1[0] << ", " << x1[1] << ", " << x1[2] << ", " << x1[3] << "],\n"
            "        [" << x2[0] << ", " << x2[1] << ", " << x2[2] << ", " << x2[3] << "],\n"
            "        [" << x3[0] << ", " << x3[1] << ", " << x3[2] << ", " << x3[3] << "],\n"
            "        [" << x4[0] << ", " << x4[1] << ", " << x4[2] << ", " << x4[3] << "]\n"
            "      ],\n"
            "      \"f\": [" << f[0] << ", " << f[1] << ", " << f[2] << ", " << f[3] << ", " << f[4] << "],\n"
            "      \"xDelta\": " << state.xDelta << ",\n"
            "      \"fDelta\": " << state.fDelta << "\n"
            "    }";
          first_iteration = false;
        }
        return true;
      }
  };

  template <typename ProblemType> class ModifiedNelderMeadSolver: public NelderMeadSolver<ProblemType> {
    public:
      using Superclass = ISolver<ProblemType, 0>;
      using typename Superclass::Scalar;
      using typename Superclass::TVector;
      using MatrixType = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

      MatrixType makeInitialSimplex (TVector &x) {
        size_t DIM = x.rows();

        // create initial simplex
        MatrixType s = MatrixType::Zero(DIM, DIM + 1);
        for (int c = 0; c < (int)DIM + 1; ++c) {
          for (int r = 0; r < (int)DIM; ++r) {
            s(r, c) = x(r);
            if (r == c - 1) {
              if (x(r) == 0) {
                s(r, c) = 0.00025;
              }
              // s(r, c) = (1 + 0.05) * x(r);
              s(r, c) = x(r) + 0.001;
            }
          }
        }

        NelderMeadSolver<ProblemType>::initialSimplexCreated = true;
        return s;
      }
  };
}


double diff (double a, double b, double c, double d) {
  // private thread variables
  long i, j, pa, pb;
  double x, y, r, rd;
  unsigned width = work_plane.cols;
  unsigned height = work_plane.rows;

  clock_t start_time = clock(), end_time;
  double elapsed;

  // Fill the CFA area with white
  if (verbose) cerr << right << setw(8) << iteration << setw(4) << run << reset << lightgrey <<  ": making mask ... " << reset;
  blank.create(work_plane.size(), work_plane.type());
#pragma omp parallel shared(blank, exr_mode) private(i, j)
  {
#pragma omp for
    for (j = 0; j < width; j++ ) {
      for (i = 0; i < height; i++ ) {
        if (exr_mode) {
          if ( // CFA interior
            i + j >= width - 1 &&               // NW boundary
            j > i - width - 1 &&                // NE boundary
            i + j < width + 2 * height - 1 &&   // SE boundary
            i > j - width                       // SW boundary
          ) {
            blank.at<ushort>(i, j) = 65535;
          }
          else {
            blank.at<ushort>(i, j) = 0;
          }
        }
        else {
          blank.at<ushort>(i, j) = 65535;
        }
      }
    }
  } /* end of parallel section */

  /// Create dst, map_x and map_y with the same size as work_plane:
  dst.create( work_plane.size(), work_plane.type() );
  map_x.create( work_plane.size(), CV_32FC1 );
  map_y.create( work_plane.size(), CV_32FC1 );

  // Map computation in the following two sections consists of two
  // slightly different variants for portrait and landscape orientation.
  // In both cases, the shorter dimension is used for radius normalization.
  run++;
  if (verbose) cerr << lightgrey <<  "computing maps ... " << reset;
  if (width >= height) { // landscape
#pragma omp parallel shared(map_x, map_y, a, b, c, d) private(i, j, x, y, r, rd)
    {
#pragma omp for
      for (j = 0; j < height; j++ ) {
        // Normalize to half-height as in IM and PanoTools.
        // The offset of 0.5 points at pixel center.
        y = (double)(2 * (j + 0.5) - height) / height;
        for (i = 0; i < width; i++ ) {
          if (blank.at<ushort>(j, i) == 65535) { // don't move masked pixels (black outer triangles)
            x = (double)(2 * (i + 0.5) - width) / height;
            r = sqrt(x * x + y * y);
            rd = a * r + b * pow(r, 2) + c * pow(r, 3) + d * pow(r, 4);
            map_x.at<float>(j, i) = (height * x * rd / r + width) / 2;
            map_y.at<float>(j, i) = (height * y * rd / r + height) / 2;
          }
        }
      }
    } /* end of parallel section */
  }

  if (width < height) { // portrait
#pragma omp parallel shared(width, height, map_x, map_y, a, b, c, d) private(i, j, x, y, r, rd)
    {
#pragma omp for
      for (j = 0; j < height; j++ ) {
        // Normalize to half-width as in IM and PanoTools.
        // The offset of 0.5 points at pixel center.
        y = (double)(2 * (j + 0.5) - height) / width;
        for (i = 0; i < width; i++ ) {
          if (blank.at<ushort>(j, i) == 65535) { // don't move masked pixels (black outer triangles)
            x = (double)(2 * (i + 0.5) - width) / width;
            r = sqrt(x * x + y * y);
            rd = a * r + b * pow(r, 2) + c * pow(r, 3) + d * pow(r, 4);
            map_x.at<float>(j, i) = (width * x * rd / r + width) / 2;
            map_y.at<float>(j, i) = (width * y * rd / r + height) / 2;
          }
        }
      }
    } /* end of parallel section */
  }

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
        if (mask.at<ushort>(j, i) == 65535) { // select unmasked
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
    cerr << lightgrey << right << setw(8) << iteration << setw(4) << run << reset <<
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
void run_find (struct argp_state* state) {
  PARSE_ARGS_FIND;

  verbose = args.verbose;
  exr_mode = args.exr;

  typedef double T;
  typedef cppoptlib::Lens<T> Lens;

  // TVector lb(4); lb << 1 - 0.05, 0 - 0.05, 0 - 0.05, 0 - 0.05;
  // TVector ub(4); lb << 1 + 0.05, 0 + 0.05, 0 + 0.05, 0 + 0.05;

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

  cerr.setf(ios::fixed, ios::floatfield);

  // initialize the optimization problem

  Lens f;
  // f.setLowerBound(lb);
  // f.setUpperBound(ub);


  // choose a starting point
  //Eigen::VectorXd x(4); x << a, b, c, d;
  Lens::TVector x(4); x << args.a, args.b, args.c, args.d;

  // choose a solver
  //cppoptlib::BfgsSolver<Lens> solver;
  //cppoptlib::LbfgsSolver<Lens> solver;
  //cppoptlib::ConjugatedGradientDescentSolver<Lens> solver;
  //cppoptlib::NewtonDescentSolver<Lens> solver;

  // cppoptlib::NelderMeadSolver<Lens> solver;
  cppoptlib::ModifiedNelderMeadSolver<Lens> solver;

  //cppoptlib::CMAesSolver<Lens> solver;

  // Create a Criteria class to set the solver's stop conditions
  Lens::TCriteria crit = Lens::TCriteria::defaults();
  crit.iterations = 100;
  crit.fDelta = 0.5;

  solver.x0 = solver.makeInitialSimplex(x);
  solver.setStopCriteria(crit);
  if (not verbose) {
    color_switch(cerr, work_plane_color);
    pbar = progressbar_new("", crit.iterations + 1);
  }

  // Prepare output
  if (args.trace) {
    trace_simplices = args.trace;
    trace_stream = &(args.trace_file);
    *trace_stream <<
      "{\n"
      "  \"simplex\": [\n";
  }

  // minimize the function
  solver.minimize(f, x);

  // Close output
  if (trace_simplices) {
    *trace_stream <<
      "\n"
      "  ],\n"
      "  \"stop\": \"" << solver.stop_condition << "\"\n"
      "}\n";
    trace_stream->close();
  }
  if (not verbose) {
    color_switch(cerr, work_plane_color);
    progressbar_finish(pbar);
    cerr << reset;
  }

  verbose = false;
  cerr << yellow << setw(12) << "stop:" << "  " << lightgrey << solver.stop_condition << reset << endl;
  cerr << lightgrey << setw(12) << "solution:" << " " << reset << setprecision(6) << x.transpose() << endl;
  color_switch(cerr, work_plane_color);
  cerr << setw(12) << "TCA:" << "  " << setprecision(2) << f(x) << reset << endl;
}

