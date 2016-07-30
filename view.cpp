#include <iostream>
#include <iomanip>

#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include <math.h>
#include <omp.h>

#include <error.h>

#include "cfa_mask.h"
#include "view_args.h"
#include "termcolor.h"

using namespace std;
using namespace cv;


// --------------------------------------------------------------------
void run_view (struct argp_state* state) {
  PARSE_ARGS_VIEW;

  int padding = args.size / 10;
  int size = args.size;

  Mat input = imread(
    args.input_file,
    CV_LOAD_IMAGE_ANYDEPTH | CV_LOAD_IMAGE_COLOR
  );

  if (args.exr) {
    // Orientation needs to be handled here
    int cfaWidth = 3264;
    int cfaHeight = 2464;

    int width = input.cols;
    int height = input.rows;
    if (height != width) {
      cerr << on_red << "A tilted EXR image must have equal width and height. The input image "
        << bold << args.input_file << reset << on_red << " is "
        << bold << width << reset << on_red << "×" << bold << height << reset << endl;
      exit(EXIT_FAILURE);
    }

    int channels = input.channels();
    if (channels != 3) {
      cerr << on_red << "Expecting a 3-channel RGB image. The input image "
        << bold << args.input_file << reset << on_red << " has "
        << bold << channels << reset << on_red << (channels == 1 ? " channel" : " channels") << reset << endl;
      exit(EXIT_FAILURE);
    }

    if (args.size > width / 2) {
      cerr << on_red << "Tiles size of "
        << bold << args.size << reset << on_red << " is too large for this input image: "
        << bold << args.input_file << reset << endl;
      exit(EXIT_FAILURE);
    }

    Mat out(3 * size + 4 * padding, 3 * size + 4 * padding, input.type());
    out.setTo(Scalar(20000, 20000, 20000));

    // Center (0)
    Mat(input, Rect(width / 2 - size / 2 + args.center_x, height / 2 - size / 2 + args.center_y, size, size))
      .copyTo(Mat(out, Rect(size + 2 * padding, size + 2 * padding, size, size)));

    // North (1)
    Mat(input, Rect(cfaWidth / 2 - size / 4 + args.n_x, cfaHeight / 2 - size / 4 + args.n_y, size, size))
      .copyTo(Mat(out, Rect(size + 2 * padding, padding, size, size)));

    // Northeast (2)
    Mat(input, Rect(cfaWidth - size / 2 + args.ne_x, 0 + args.ne_y, size, size))
      .copyTo(Mat(out, Rect(2 * size + 3 * padding, padding, size, size)));

    // East (3)
    Mat(input, Rect(cfaWidth + (cfaHeight - 1.25 * size) / 2 + args.e_x, (cfaHeight - 0.75 * size) / 2 + args.e_y, size, size))
      .copyTo(Mat(out, Rect(2 * size + 3 * padding, size + 2 * padding, size, size)));

    // Southeast (4)
    Mat(input, Rect(width - size + args.se_x, cfaHeight - size / 2 + args.se_y, size, size))
      .copyTo(Mat(out, Rect(2 * size + 3 * padding, 2 * size + 3 * padding, size, size)));

    // South (5)
    Mat(input, Rect(cfaWidth + cfaHeight / 2 - 0.75 * size + args.s_x, cfaHeight + cfaWidth / 2 - 0.75 * size + args.s_y, size, size))
      .copyTo(Mat(out, Rect(size + 2 * padding, 2 * size + 3 * padding, size, size)));

    // Southwest (6)
    Mat(input, Rect(cfaHeight - size / 2 + args.sw_x, cfaWidth + cfaHeight - size + args.sw_y, size, size))
      .copyTo(Mat(out, Rect(padding, 2 * size + 3 * padding, size, size)));

    // West (7)
    Mat(input, Rect(cfaHeight / 2 - 0.25 * size + args.w_x, cfaWidth + (cfaHeight - 1.25 * size) / 2 + args.w_y, size, size))
      .copyTo(Mat(out, Rect(padding, size + 2 * padding, size, size)));

    // Northwest (8)
    Mat(input, Rect(0 + args.nw_x, cfaWidth - size / 2 + args.nw_y, size, size))
      .copyTo(Mat(out, Rect(padding, padding, size, size)));

    imwrite(args.output_file, out);
  }

  else {
    cerr << on_red << "Non-EXR processing is not done yet" << reset << endl;
    exit(EXIT_FAILURE);
  }
}
