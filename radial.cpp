#include <iostream>
#include <iomanip>

#include <Magick++.h>

#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include <math.h>
#include <omp.h>

#include <error.h>

#include "radial_args.h"
#include "termcolor.h"

using namespace std;
using namespace Magick;
using namespace cv;


// --------------------------------------------------------------------
void run_radial (struct argp_state* state) {
  PARSE_ARGS_RADIAL;
  double a = args.a;
  double b = args.b;
  double c = args.c;
  double d = args.d;
  int nthreads = omp_get_max_threads();


  if (args.fast) {
    long i, j, pa, pb;
    double x, y, r, rd;
    unsigned width, height;

    Mat work_plane, dst, map_x, map_y;

    cerr << lightgrey << "Running " << reset << "OpenCV remap()"
      << lightgrey << " with " <<  nthreads << " threads" << reset << endl;
    cerr << lightgrey << "loading input file " << reset << args.input_file << "\r" << flush;
    work_plane = imread( args.input_file, CV_LOAD_IMAGE_ANYDEPTH );
    cerr << "\r" << string(100, ' ') << "\r" << flush; // blank line
    width = work_plane.cols;
    height = work_plane.rows;

    /// Create dst, map_x and map_y with the same size as work_plane:
    dst.create( work_plane.size(), work_plane.type() );
    map_x.create( work_plane.size(), CV_32FC1 );
    map_y.create( work_plane.size(), CV_32FC1 );

    // Map computation in the following two sections consists of two
    // slightly different variants for portrait and landscape orientation.
    // In both cases, the shorter dimension is used for radius normalization.
    cerr << lightgrey <<  "computing maps ... " << reset;
    if (width >= height) { // landscape
#pragma omp parallel shared(width, height, map_x, map_y, a, b, c, d) private(i, j, x, y, r, rd)
      {
#pragma omp for
        for (j = 0; j < height; j++ ) {
          // Normalize to half-height as in IM and PanoTools.
          // The offset of 0.5 points at pixel center.
          y = (double)(2 * (j + 0.5) - height) / height;
          for (i = 0; i < width; i++ ) {
            x = (double)(2 * (i + 0.5) - width) / height;
            r = sqrt(x * x + y * y);
            rd = a * r + b * pow(r, 2) + c * pow(r, 3) + d * pow(r, 4);
            map_x.at<float>(j, i) = (height * x * rd / r + width) / 2 + args.shift_x;
            map_y.at<float>(j, i) = (height * y * rd / r + height) / 2 + args.shift_y;
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
            x = (double)(2 * (i + 0.5) - width) / width;
            r = sqrt(x * x + y * y);
            rd = a * r + b * pow(r, 2) + c * pow(r, 3) + d * pow(r, 4);
            map_x.at<float>(j, i) = (width * x * rd / r + width) / 2 + args.shift_x;
            map_y.at<float>(j, i) = (width * y * rd / r + height) / 2 + args.shift_y;
          }
        }
      } /* end of parallel section */
    }

    cerr << lightgrey << "remapping ... " << reset;
    remap( work_plane, dst, map_x, map_y, CV_INTER_LANCZOS4, BORDER_CONSTANT, Scalar(0,0, 0) );

    cerr << lightgrey << "writing ... " << reset;
    imwrite(args.output_file, dst);
    cerr << lightgrey << "\r                                                                               \r" << reset;
  }

  else {
    // Slower, high quality
    cerr << lightgrey << "Running " << reset << "Magick::BarrelDistortion"
      << lightgrey << " with " <<  nthreads << " threads" << reset << endl;
    Magick::ResourceLimits::thread(nthreads);

    double arguments[4] = {d, c, b, a};

    Magick::FilterType filter(LanczosFilter);

    string input_file(args.input_file);
    string output_file(args.output_file);

    try {
      Image image(input_file);

      image.filterType(filter);
      image.virtualPixelMethod(BlackVirtualPixelMethod);

      image.distort(Magick::BarrelDistortion, 4, arguments);

      image.write(output_file);
    }
    catch( exception &e ) {
      cerr << "Caught exception: " << e.what() << endl;
      exit(EXIT_FAILURE);
    }
  }
}
