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
    out.setTo(Scalar(0, 0, 0));

    // Center
    Mat(input, Rect(width / 2 - size / 2, height / 2 - size / 2, size, size))
      .copyTo(Mat(out, Rect(size + 2 * padding, size + 2 * padding, size, size)));

    // North
    Mat(input, Rect((int)(cfaWidth / 2 - size / 4), (int)(cfaHeight / 2 - size / 4), size, size))
      .copyTo(Mat(out, Rect(size + 2 * padding, padding, size, size)));

    // Northeast
    Mat(input, Rect((int)(cfaWidth - size / 2), 0, size, size))
      .copyTo(Mat(out, Rect(2 * size + 3 * padding, padding, size, size)));

    // East
    Mat(input, Rect((int)(cfaWidth + (cfaHeight - 1.25 * size) / 2), (int)((cfaHeight - 0.75 * size) / 2), size, size))
      .copyTo(Mat(out, Rect(2 * size + 3 * padding, size + 2 * padding, size, size)));

    // Southeast
    Mat(input, Rect(width - size, cfaHeight - size / 2, size, size))
      .copyTo(Mat(out, Rect(2 * size + 3 * padding, 2 * size + 3 * padding, size, size)));

    // South
    Mat(input, Rect((int)(cfaWidth  + cfaHeight / 2 - 0.75 * size), (int)(cfaHeight + cfaWidth / 2 - 0.75 * size), size, size))
      .copyTo(Mat(out, Rect(size + 2 * padding, 2 * size + 3 * padding, size, size)));

    // Southwest
    Mat(input, Rect(cfaHeight - size / 2, cfaWidth + cfaHeight - size, size, size))
      .copyTo(Mat(out, Rect(padding, 2 * size + 3 * padding, size, size)));

    // West
    Mat(input, Rect((int)(cfaHeight / 2 - 0.25 * size), (int)(cfaWidth + (cfaHeight - 1.25 * size) / 2), size, size))
      .copyTo(Mat(out, Rect(padding, size + 2 * padding, size, size)));

    // Northwest
    Mat(input, Rect(0, (int)(cfaWidth - size / 2), size, size))
      .copyTo(Mat(out, Rect(padding, padding, size, size)));

    imwrite(args.output_file, out);
  }

  else {
    cerr << on_red << "Non-EXR processing is not done yet" << reset << endl;
    exit(EXIT_FAILURE);
  }
}
