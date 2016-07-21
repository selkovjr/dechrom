dcraw -d -s all -4 -T 160206_172303.RAF

fuji-exr linear -h 160206_172303_* interpolated.tiff # EXR
fuji-exr linear 151018_180303_0.tiff interpolated.tiff # subframe

convert interpolated.tiff -separate interpolated-%d.tiff
convert interpolated-0.tiff -equalize equalized-0.tiff
convert interpolated-1.tiff -equalize equalized-1.tiff
convert interpolated-2.tiff -equalize equalized-2.tiff

./radial-distort-grid 0 0 0 0 equalized-1.tiff equalized-0.tiff dummy > surface.tab
./radial-distort-optimize 0.9996 0 0 0 equalized-1.tiff equalized-0.tiff > solution.tab

./radial-distort 1.002002 -0.006203 0.008245 -0.003979 interpolated-0.tiff interpolated-distorted-0.tiff
./radial-distort 1.000725 -0.000260 -0.001201 0.000909 interpolated-2.tiff interpolated-distorted-2.tiff
convert interpolated-distorted-0.tiff interpolated-1.tiff interpolated-distorted-2.tiff -set colorspace RGB -combine -set colorspace sRGB interpolated-combined.tiff
fuji-exr ssd -m 3264x2464 interpolated-distorted-0.tiff interpolated-1.tiff interpolated-distorted-2.tiff out.tiff
