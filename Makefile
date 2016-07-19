GIT_VERSION := $(shell git describe --abbrev=4 --dirty --always --tags)

OBJ = dechrom.o survey.o cfa_mask.o write_tiff.o progressbar.o
BIN = dechrom
LIBBIN=.

CXX=g++

HDRDIR=-I/opt/local/include/ -I/usr/local/include/
LIBDIR=-L/opt/local/lib/ -L/usr/local/lib/
OPNCV_LDFLAGS=-L${exec_prefix}/lib -lopencv_contrib -lopencv_stitching -lopencv_nonfree -lopencv_superres -lopencv_ocl -lopencv_ts -lopencv_videostab -lopencv_gpu -lopencv_photo -lopencv_objdetect -lopencv_legacy -lopencv_video -lopencv_ml -lopencv_calib3d -lopencv_features2d -lopencv_highgui -lIlmImf -ltiff -lopencv_imgproc -lopencv_flann -lopencv_core -L/usr/lib/x86_64-linux-gnu -lbz2 -lswscale-ffmpeg -lavutil-ffmpeg -lavformat-ffmpeg -lavcodec-ffmpeg -ldc1394 -lgthread-2.0 -lfreetype -lfontconfig -lglib-2.0 -lgobject-2.0 -lpango-1.0 -lpangoft2-1.0 -lgio-2.0 -lgdk_pixbuf-2.0 -lcairo -latk-1.0 -lpangocairo-1.0 -lgdk-x11-2.0 -lgtk-x11-2.0 -ljasper -lpng -ljpeg -lrt -lpthread -lm -ldl -lstdc++ -lz


CFLAGS += -g -O3 -D VERSION=\"$(GIT_VERSION)\" -fdiagnostics-color=auto \
  -fopenmp -funroll-loops -fomit-frame-pointer  -fno-tree-pre -falign-loops -ffast-math -ftree-vectorize \
  -Weffc++ -pedantic -Wall -Wextra  -Wno-write-strings -Wno-deprecated  $(HDRDIR)

LDFLAGS += -g $(CFLAGS) $(LIBDIR) -lncurses -lgomp -lpthread \
  `pkg-config --static --libs opencv`


LIBMX= survey.o cfa_mask.o write_tiff.o progressbar.o

default: $(OBJ) $(BIN)

$(OBJ) : %.o : %.cpp survey_args.h
	$(CXX) -c $(CFLAGS)  $< -o $@

$(BIN) : % : %.o  $(LIBMX)
	$(CXX) -o $(LIBBIN)/$@  $^ $(LDFLAGS)


.PHONY : clean
clean:
	$(RM) $(OBJ); cd $(LIBBIN); rm -f $(BIN)
