#include <iostream>
#include <iomanip>
#include <string>

#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include <math.h>
#include <stdio.h>

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <error.h>
#include <ctime>

#include "survey_args.h"
#include "termcolor.h"
#include "cfa_mask.h"
#include "write_tiff.h"

using namespace std;
using namespace cv;
using namespace termcolor;

/// Global variables
double a;
double b;
double c;
double d;
Mat ref_image, src;

int run = 0;
char suffix[10];

long i, j, pa, pb;

 double diff () {
   Mat dst;
   Mat map_x, map_y;

   /// Create dst, map_x and map_y with the same size as src:
   dst.create( src.size(), src.type() );
   map_x.create( src.size(), CV_32FC1 );
   map_y.create( src.size(), CV_32FC1 );

   cerr << "computing maps ... ";
 #pragma omp parallel shared(map_x, map_y, a, b, c, d) private(i, j)
   {
 #pragma omp for collapse(2)
     for (j = 0; j < src.rows; j++ ) {
       for (i = 0; i < src.cols; i++ ) {
         if (src.at<ushort>(i,j) > 0) { // don't move masked pixels (black outer triangles)
           double x = (double)(2 * i - src.rows) / src.rows;
           double y = (double)(2 * j - src.cols) / src.cols;
           double r = sqrt(x * x + y * y);
           double rd = a * r + b * pow(r, 2) + c * pow(r, 3) + d * pow(r, 4);
           map_x.at<float>(j,i) = (x * rd / r + 1) * src.rows / 2;
           map_y.at<float>(j,i) = (y * rd / r + 1) * src.cols / 2;
         }
       }
     }
   } /* end of parallel section */
   cerr << "remapping ... ";
   remap( src, dst, map_x, map_y, CV_INTER_LINEAR, BORDER_CONSTANT, Scalar(0,0, 0) );

   /// Compute difference
   cerr << "computing difference ... ";
   unsigned long num = 0;
   unsigned long quot = 0;
 #pragma omp parallel shared(dst, ref_image) private(i, j, pa, pb) reduction(+:num,quot)
   {
 #pragma omp for collapse(2)
     for ( j = 0; j < src.rows; j++ ) {
       for ( i = 0; i < src.cols; i++ ) {
         if (dst.at<ushort>(i, j) > 0 and ref_image.at<ushort>(i,j) > 0) { // select unmasked
           pa = dst.at<ushort>(j, i);
           pb = ref_image.at<ushort>(j,i);
           num += abs((double)pa - (double)pb);
           quot += 1;
         }
       }
     }
   } /* end of parallel section */
   cerr << "done.\n";

   fprintf(stderr, "%d\t%f\t%f\t%f\t%f\t%f\n", ++run, a, b, c, d, 1.0 * num / quot);
   printf("%f\t%f\t%f\t%f\t%f\n", a, b, c, d, 1.0 * num / quot);
   return 1.0 * num / quot;
 }

// --------------------------------------------------------------------
void run_survey (struct argp_state* state) {
  PARSE_ARGS_SURVEY;

  cerr << args.nodes << endl;
  cerr << "a: (" << args.amin << ", " << args.amax << ")" << endl;
  cerr << "b: (" << args.bmin << ", " << args.bmax << ")" << endl;
  cerr << "c: (" << args.cmin << ", " << args.cmax << ")" << endl;
  cerr << "d: (" << args.dmin << ", " << args.dmax << ")" << endl;
  cerr << "ref: " << args.ref_file << endl;
  cerr << "target: " << args.target_file << endl;

  clock_t start_time, end_time;
  double elapsed;

  cerr.setf(ios::fixed, ios::floatfield);

  double astep = abs(args.amax - args.amin) / (args.nodes - 1);
  double bstep = abs(args.bmax - args.bmin) / (args.nodes - 1);
  double cstep = abs(args.cmax - args.cmin) / (args.nodes - 1);
  double dstep = abs(args.dmax - args.dmin) / (args.nodes - 1);

  // Load the images
  cerr << "loading reference image " << args.ref_file << endl;
  ref_image = imread( args.ref_file, CV_LOAD_IMAGE_ANYDEPTH ); // grayscale
  cerr << "loading target image " << args.target_file << endl;
  src = imread( args.target_file, CV_LOAD_IMAGE_ANYDEPTH );

  for (unsigned ac = 0; ac < args.nodes; ac++) {
    a = args.amin + ac * astep;
    for (unsigned bc = 0; bc < args.nodes; bc++) {
      b = args.bmin + bc * bstep;
      for (unsigned cc = 0; cc < args.nodes; cc++) {
        c = args.cmin + cc * cstep;
        for (unsigned dc = 0; dc < args.nodes; dc++) {
          d = args.dmin + dc * dstep;
          diff(); // parameters are global
          if (args.dmin == args.dmax) break;
        }
        if (args.cmin == args.cmax) break;
      }
      if (args.bmin == args.bmax) break;
    }
    if (args.amin == args.amax) break;
  }

} // run_survey()

